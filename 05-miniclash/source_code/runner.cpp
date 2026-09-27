#include <algorithm>
#include <cerrno>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <fcntl.h>
#include <sched.h>
#include <sys/prctl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {

struct Task {
    std::string input;
    std::string out1;
    std::string out2;
    bool done = false;
    unsigned running = 0;
    unsigned attempts = 0;
    std::chrono::steady_clock::time_point first_start{};
};

struct Attempt {
    std::size_t task = 0;
    std::string tmp1;
    std::string tmp2;
};

std::uint64_t mix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

unsigned allowed_cpus() {
    cpu_set_t set;
    CPU_ZERO(&set);
    if (::sched_getaffinity(0, sizeof(set), &set) == 0) {
        const int n = CPU_COUNT(&set);
        if (n > 0) return static_cast<unsigned>(n);
    }
    const long n = ::sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? static_cast<unsigned>(n) : 1U;
}

std::string dirname_of(const std::string& p) {
    const std::string::size_type pos = p.find_last_of('/');
    if (pos == std::string::npos) return ".";
    if (pos == 0) return "/";
    return p.substr(0, pos);
}

bool read_tasks(const char* path, std::vector<Task>& tasks) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "cannot read task file: " << path << "\n";
        return false;
    }

    std::string line;
    unsigned lineno = 0;
    while (std::getline(in, line)) {
        ++lineno;
        std::istringstream iss(line);
        Task t;
        std::string extra;
        if (!(iss >> t.input)) continue;
        if (!(iss >> t.out1 >> t.out2) || (iss >> extra)) {
            std::cerr << "malformed task line " << lineno << "\n";
            return false;
        }
        tasks.push_back(std::move(t));
    }
    return true;
}

void cleanup(const Attempt& a) {
    ::unlink(a.tmp1.c_str());
    ::unlink(a.tmp2.c_str());
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: " << argv[0] << " <tasks.txt>\n";
        return 2;
    }

    std::vector<Task> tasks;
    if (!read_tasks(argv[1], tasks)) return 1;
    if (tasks.empty()) return 0;

    unsigned workers = std::min<unsigned>(allowed_cpus(), static_cast<unsigned>(tasks.size()));
    if (const char* s = std::getenv("MINICLASH_JOBS")) {
        char* end = nullptr;
        const unsigned long v = std::strtoul(s, &end, 10);
        if (end != s && *end == '\0' && v > 0)
            workers = std::min<unsigned>(workers, static_cast<unsigned>(v));
    }
    workers = std::max(1U, workers);

    const std::string bin = dirname_of(argv[0]) + "/md5fastcoll";
    if (::access(bin.c_str(), X_OK) != 0) {
        std::cerr << "md5fastcoll not executable: " << bin << "\n";
        return 1;
    }

    timespec ts{};
    ::clock_gettime(CLOCK_MONOTONIC, &ts);
    const std::uint64_t nonce =
        (static_cast<std::uint64_t>(ts.tv_sec) << 32) ^
        static_cast<std::uint64_t>(ts.tv_nsec) ^
        static_cast<std::uint64_t>(::getpid());

    std::unordered_map<pid_t, Attempt> active;
    active.reserve(workers * 3U);
    std::size_t next_original = 0;
    std::size_t completed = 0;
    bool fatal = false;

    auto launch = [&](std::size_t ti) -> bool {
        Task& t = tasks[ti];
        const unsigned attempt_no = ++t.attempts;
        if (t.running == 0)
            t.first_start = std::chrono::steady_clock::now();

        const std::string suffix =
            ".mc." + std::to_string(static_cast<long long>(::getpid())) +
            "." + std::to_string(ti) + "." + std::to_string(attempt_no);
        Attempt a;
        a.task = ti;
        a.tmp1 = t.out1 + suffix + ".1";
        a.tmp2 = t.out2 + suffix + ".2";
        ::unlink(a.tmp1.c_str());
        ::unlink(a.tmp2.c_str());

        const std::uint64_t z1 = mix64(nonce ^ (static_cast<std::uint64_t>(ti) << 32) ^ attempt_no);
        const std::uint64_t z2 = mix64(z1 ^ 0xd1b54a32d192ed03ULL);
        std::uint32_t seed1 = static_cast<std::uint32_t>(z1);
        std::uint32_t seed2 = static_cast<std::uint32_t>(z2);
        if ((seed1 | seed2) == 0) seed2 = 1;
        const std::string s1 = std::to_string(seed1);
        const std::string s2 = std::to_string(seed2);

        const pid_t pid = ::fork();
        if (pid < 0) {
            std::cerr << "fork failed: " << std::strerror(errno) << "\n";
            cleanup(a);
            return false;
        }

        if (pid == 0) {
            ::prctl(PR_SET_PDEATHSIG, SIGKILL);
            if (::getppid() == 1) _exit(125);

            const int devnull = ::open("/dev/null", O_WRONLY);
            if (devnull >= 0) {
                ::dup2(devnull, STDOUT_FILENO);
                ::dup2(devnull, STDERR_FILENO);
                if (devnull > STDERR_FILENO) ::close(devnull);
            }

            ::execl(bin.c_str(), bin.c_str(),
                    "-q",
                    "--seed1", s1.c_str(),
                    "--seed2", s2.c_str(),
                    "-p", t.input.c_str(),
                    "-o", a.tmp1.c_str(), a.tmp2.c_str(),
                    static_cast<char*>(nullptr));
            _exit(127);
        }

        ++t.running;
        active.emplace(pid, std::move(a));
        return true;
    };

    auto choose_speculative = [&]() -> std::size_t {
        std::size_t best = tasks.size();
        unsigned best_running = std::numeric_limits<unsigned>::max();
        auto best_age = std::chrono::steady_clock::duration::min();
        const auto now = std::chrono::steady_clock::now();

        for (std::size_t i = 0; i < tasks.size(); ++i) {
            const Task& t = tasks[i];
            if (t.done) continue;
            if (t.running == 0) return i;
            const auto age = now - t.first_start;
            if (t.running < best_running ||
                (t.running == best_running && age > best_age)) {
                best = i;
                best_running = t.running;
                best_age = age;
            }
        }
        return best;
    };

    auto fill_slots = [&]() -> bool {
        while (!fatal && active.size() < workers && completed < tasks.size()) {
            std::size_t ti;
            if (next_original < tasks.size()) {
                ti = next_original++;
            } else {
                ti = choose_speculative();
                if (ti == tasks.size()) break;
            }
            if (!launch(ti)) return false;
        }
        return true;
    };

    if (!fill_slots()) fatal = true;

    while (!fatal && completed < tasks.size()) {
        int status = 0;
        const pid_t pid = ::waitpid(-1, &status, 0);
        if (pid < 0) {
            if (errno == EINTR) continue;
            std::cerr << "waitpid failed: " << std::strerror(errno) << "\n";
            fatal = true;
            break;
        }

        auto it = active.find(pid);
        if (it == active.end()) continue;

        const Attempt a = std::move(it->second);
        active.erase(it);
        Task& t = tasks[a.task];
        if (t.running) --t.running;

        const bool ok = WIFEXITED(status) && WEXITSTATUS(status) == 0;
        if (ok && !t.done) {
            if (::rename(a.tmp1.c_str(), t.out1.c_str()) != 0 ||
                ::rename(a.tmp2.c_str(), t.out2.c_str()) != 0) {
                std::cerr << "failed to install outputs for " << t.input
                          << ": " << std::strerror(errno) << "\n";
                cleanup(a);
                fatal = true;
                break;
            }

            t.done = true;
            ++completed;

            for (const auto& kv : active) {
                if (kv.second.task == a.task)
                    ::kill(kv.first, SIGKILL);
            }
        } else {
            cleanup(a);
        }

        if (!fill_slots()) fatal = true;
    }

    if (fatal) {
        for (const auto& kv : active) ::kill(kv.first, SIGKILL);
    }

    for (auto& kv : active) cleanup(kv.second);
    while (::waitpid(-1, nullptr, 0) > 0) {}

    if (fatal) return 1;
    return completed == tasks.size() ? 0 : 1;
}
