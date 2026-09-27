#include "main.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iterator>
#include <mutex>
#include <pthread.h>
#include <sched.h>
#include <sstream>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

namespace {

constexpr uint32 kMD5IV[4] = {
    0x67452301u, 0xefcdab89u, 0x98badcfeu, 0x10325476u
};

struct Task {
    std::string input;
    std::string output1;
    std::string output2;

    std::vector<unsigned char> prefix;
    std::array<uint32, 4> iv{};

    // Protected by the global scheduler mutex.
    unsigned state = 0;  // 0=pending, 1=winner is committing, 2=done
    unsigned running = 0;
    uint64_t attempts = 0;
    std::chrono::steady_clock::time_point first_start{};
};

static inline uint64_t splitmix64(uint64_t& x)
{
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

void seed_attempt(size_t task_index, uint64_t attempt)
{
    uint64_t x =
        0x243f6a8885a308d3ULL ^
        (uint64_t(task_index) * 0x9e3779b97f4a7c15ULL) ^
        (attempt * 0xd1b54a32d192ed03ULL) ^
        uint64_t(std::chrono::high_resolution_clock::now().time_since_epoch().count()) ^
        (uint64_t(uintptr_t(pthread_self())) << 1);

    seed32_1 = uint32(splitmix64(x));
    seed32_2 = uint32(splitmix64(x));
    if ((seed32_1 | seed32_2) == 0)
        seed32_2 = 0x12345678u;
}

void find_collision_local(const uint32 IV[],
                          uint32 msg1block0[], uint32 msg1block1[],
                          uint32 msg2block0[], uint32 msg2block1[])
{
    find_block0(msg1block0, IV);

    uint32 IHV[4] = { IV[0], IV[1], IV[2], IV[3] };
    md5_compress(IHV, msg1block0);

    find_block1(msg1block1, IHV);

    for (int t = 0; t < 16; ++t) {
        msg2block0[t] = msg1block0[t];
        msg2block1[t] = msg1block1[t];
    }

    msg2block0[4] += 1u << 31;
    msg2block0[11] += 1u << 15;
    msg2block0[14] += 1u << 31;

    msg2block1[4] += 1u << 31;
    msg2block1[11] -= 1u << 15;
    msg2block1[14] += 1u << 31;
}

bool prepare_task(Task& task, std::string& error)
{
    std::ifstream in(task.input.c_str(), std::ios::binary);
    if (!in) {
        error = "cannot open input file: " + task.input;
        return false;
    }

    std::vector<unsigned char> raw(
        (std::istreambuf_iterator<char>(in)),
        std::istreambuf_iterator<char>());

    const size_t padded = (raw.size() + 63u) & ~size_t(63u);
    task.prefix.assign(padded, 0);
    std::copy(raw.begin(), raw.end(), task.prefix.begin());

    task.iv = { kMD5IV[0], kMD5IV[1], kMD5IV[2], kMD5IV[3] };

    uint32 block[16];
    for (size_t off = 0; off < padded; off += 64) {
        const unsigned char* p = task.prefix.data() + off;
        for (unsigned k = 0; k < 16; ++k) {
            const unsigned char* q = p + 4 * k;
            block[k] =
                uint32(q[0]) |
                (uint32(q[1]) << 8) |
                (uint32(q[2]) << 16) |
                (uint32(q[3]) << 24);
        }
        md5_compress(task.iv.data(), block);
    }

    return true;
}

static inline void encode_block(const uint32 block[16], unsigned char out[64])
{
    for (unsigned k = 0; k < 16; ++k) {
        const uint32 x = block[k];
        out[4*k + 0] = static_cast<unsigned char>(x);
        out[4*k + 1] = static_cast<unsigned char>(x >> 8);
        out[4*k + 2] = static_cast<unsigned char>(x >> 16);
        out[4*k + 3] = static_cast<unsigned char>(x >> 24);
    }
}

bool write_all(int fd, const void* data, size_t size)
{
    const unsigned char* p = static_cast<const unsigned char*>(data);
    while (size != 0) {
        const ssize_t n = ::write(fd, p, size);
        if (n < 0) {
            if (errno == EINTR)
                continue;
            return false;
        }
        p += n;
        size -= size_t(n);
    }
    return true;
}

bool write_one(const std::string& path,
               const std::vector<unsigned char>& prefix,
               const uint32 block0[16], const uint32 block1[16],
               std::string& error)
{
    const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (fd < 0) {
        error = "cannot create output: " + path + ": " + std::strerror(errno);
        return false;
    }

    unsigned char suffix[128];
    encode_block(block0, suffix);
    encode_block(block1, suffix + 64);

    bool ok =
        (prefix.empty() || write_all(fd, prefix.data(), prefix.size())) &&
        write_all(fd, suffix, sizeof(suffix));

    const int saved_errno = errno;
    if (::close(fd) != 0 && ok) {
        ok = false;
        error = "close failed for output: " + path + ": " + std::strerror(errno);
    } else if (!ok) {
        error = "write failed for output: " + path + ": " + std::strerror(saved_errno);
    }

    return ok;
}

bool write_result(const Task& task,
                  const uint32 msg1block0[16], const uint32 msg1block1[16],
                  const uint32 msg2block0[16], const uint32 msg2block1[16],
                  std::string& error)
{
    if (!write_one(task.output1, task.prefix, msg1block0, msg1block1, error))
        return false;
    if (!write_one(task.output2, task.prefix, msg2block0, msg2block1, error))
        return false;
    return true;
}

std::vector<int> available_cpus()
{
    cpu_set_t set;
    CPU_ZERO(&set);
    std::vector<int> cpus;

    if (::sched_getaffinity(0, sizeof(set), &set) == 0) {
        for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
            if (CPU_ISSET(cpu, &set))
                cpus.push_back(cpu);
    }

    if (cpus.empty())
        cpus.push_back(0);
    if (cpus.size() > 32)
        cpus.resize(32);
    return cpus;
}

void pin_this_thread(int cpu)
{
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    (void)::pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <tasks.txt>\n", argv[0]);
        return 2;
    }

    std::ifstream task_file(argv[1]);
    if (!task_file) {
        std::fprintf(stderr, "cannot open task file: %s\n", argv[1]);
        return 1;
    }

    std::vector<Task> tasks;
    std::string line;
    while (std::getline(task_file, line)) {
        if (line.empty())
            continue;

        std::istringstream iss(line);
        Task t;
        std::string extra;
        if (!(iss >> t.input >> t.output1 >> t.output2) || (iss >> extra)) {
            std::fprintf(stderr, "malformed task line: %s\n", line.c_str());
            return 1;
        }
        tasks.push_back(std::move(t));
    }

    if (tasks.empty())
        return 0;

    // Prefix parsing and IHV derivation are deterministic and tiny compared to
    // collision search. Do them exactly once per task so hedged attempts share
    // the result instead of re-reading/re-hashing the same prefix.
    for (Task& t : tasks) {
        std::string error;
        if (!prepare_task(t, error)) {
            std::fprintf(stderr, "%s\n", error.c_str());
            return 1;
        }
    }

    const std::vector<int> cpus = available_cpus();
    const size_t worker_count = std::min(cpus.size(), tasks.size() < cpus.size()
                                                       ? cpus.size()
                                                       : cpus.size());

    std::mutex mutex;
    std::condition_variable cv;
    size_t next_unstarted = 0;
    size_t completed = 0;
    bool failed = false;
    std::string failure_message;

    auto worker = [&](size_t worker_id) {
        pin_this_thread(cpus[worker_id]);

        for (;;) {
            size_t task_index = 0;
            uint64_t attempt = 0;

            {
                std::unique_lock<std::mutex> lock(mutex);

                for (;;) {
                    if (failed || completed == tasks.size())
                        return;

                    // Throughput phase: start every independent task exactly
                    // once before spending a core on replication.
                    if (next_unstarted < tasks.size()) {
                        task_index = next_unstarted++;
                        Task& t = tasks[task_index];
                        ++t.running;
                        attempt = ++t.attempts;
                        t.first_start = std::chrono::steady_clock::now();
                        break;
                    }

                    // Tail phase: every idle core becomes a speculative copy.
                    // Prefer the fewest-running task; tie-break by oldest start.
                    size_t best = tasks.size();
                    for (size_t i = 0; i < tasks.size(); ++i) {
                        const Task& t = tasks[i];
                        if (t.state != 0)
                            continue;
                        if (best == tasks.size() ||
                            t.running < tasks[best].running ||
                            (t.running == tasks[best].running &&
                             t.first_start < tasks[best].first_start)) {
                            best = i;
                        }
                    }

                    if (best != tasks.size()) {
                        task_index = best;
                        Task& t = tasks[task_index];
                        ++t.running;
                        attempt = ++t.attempts;
                        break;
                    }

                    // All remaining tasks are currently committing winners.
                    cv.wait(lock);
                }
            }

            seed_attempt(task_index, attempt);

            uint32 msg1block0[16], msg1block1[16];
            uint32 msg2block0[16], msg2block1[16];
            find_collision_local(tasks[task_index].iv.data(),
                                 msg1block0, msg1block1,
                                 msg2block0, msg2block1);

            bool winner = false;
            {
                std::lock_guard<std::mutex> lock(mutex);
                Task& t = tasks[task_index];
                if (t.running != 0)
                    --t.running;
                if (t.state == 0) {
                    t.state = 1;
                    winner = true;
                }
            }

            if (!winner) {
                cv.notify_all();
                continue;
            }

            std::string error;
            if (!write_result(tasks[task_index],
                              msg1block0, msg1block1,
                              msg2block0, msg2block1,
                              error)) {
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    failed = true;
                    failure_message = error;
                }
                cv.notify_all();
                return;
            }

            {
                std::lock_guard<std::mutex> lock(mutex);
                Task& t = tasks[task_index];
                t.state = 2;
                ++completed;
            }
            cv.notify_all();
        }
    };

    // Detached workers let the process terminate immediately when every task
    // has a committed winner. Losing hedge searches need not run to completion.
    for (size_t i = 0; i < worker_count; ++i)
        std::thread(worker, i).detach();

    {
        std::unique_lock<std::mutex> lock(mutex);
        cv.wait(lock, [&] { return failed || completed == tasks.size(); });
    }

    if (failed) {
        std::fprintf(stderr, "%s\n", failure_message.c_str());
        std::fflush(stderr);
        std::_Exit(1);
    }

    // All winning outputs are already closed before completed is incremented.
    // _Exit kills any still-running speculative losers without join latency.
    std::_Exit(0);
}
