#include <iostream>
#include <fstream>
#include <time.h>

#include "main.hpp"
#include "timer.hpp"

using namespace std;

const uint32 MD5IV[] = { 0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476 };

unsigned load_block(istream& i, uint32 block[]);
void save_block(ostream& o, const uint32 block[]);
void find_collision(const uint32 IV[], uint32 msg1block0[], uint32 msg1block1[], uint32 msg2block0[], uint32 msg2block1[], bool verbose = false);

#if 0

// example trivial version with md5 initial value
int main()
{
	seed32_1 = uint32(time(NULL));
	seed32_2 = 0x12345678;

	uint32 IV[4] = { MD5IV[0], MD5IV[1], MD5IV[2], MD5IV[3] };
	uint32 msg1block0[16];
	uint32 msg1block1[16];
	uint32 msg2block0[16];
	uint32 msg2block1[16];
	find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1, verbose);

	ofstream ofs1("msg1", ios::binary);
	save_block(ofs1, msg1block0);
	save_block(ofs1, msg1block1);
	ofstream ofs2("msg2", ios::binary);
	save_block(ofs2, msg2block0);
	save_block(ofs2, msg2block1);
}

#else

#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <stdint.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sched.h>
#include <signal.h>
#include <unistd.h>

typedef ::uint64_t uint64;

namespace {

const char* const options_help =
	"Allowed options:\n"
	"  -h [ --help ]           Show options.\n"
	"  -q [ --quiet ]          Be less verbose.\n"
	"  -i [ --ihv ] arg        Use specified initial value. Default is MD5 initial \n"
	"                          value.\n"
	"  --seed1 arg             Specify a SEED1 value.  Default is based on time.\n"
	"  --seed2 arg             Specify a SEED2 value.  Default is a constant value.\n"
	"  -p [ --prefixfile ] arg Calculate initial value using given prefixfile. Also \n"
	"                          copies data to output files.\n"
	"  -o [ --out ] arg        Set output filenames. This must be the last option \n"
	"                          and exactly 2 filenames must be specified. \n"
	"                          Default: -o msg1.bin msg2.bin\n";

bool file_exists(const std::string& path)
{
	struct stat st;
	return ::stat(path.c_str(), &st) == 0;
}

std::string uint_to_string(unsigned value)
{
	std::ostringstream oss;
	oss << value;
	return oss.str();
}

uint32 parse_uint32(const std::string& name, const std::string& value)
{
	// An optional sign followed by digits only; no surrounding whitespace.
	// Out-of-range values are rejected, while an in-range negative value wraps.
	std::string::size_type pos = 0;
	bool negative = false;
	if (pos < value.size() && (value[pos] == '+' || value[pos] == '-'))
		negative = (value[pos++] == '-');

	bool ok = pos < value.size();
	for (std::string::size_type i = pos; ok && i < value.size(); ++i)
		if (!isdigit((unsigned char)value[i]))
			ok = false;

	uint64 magnitude = 0;
	if (ok)
	{
		errno = 0;
		magnitude = strtoull(value.c_str() + pos, 0, 10);
		if (errno == ERANGE || magnitude > 0xFFFFFFFFull)
			ok = false;
	}

	if (!ok)
		throw std::runtime_error("the argument ('" + value + "') for option '--"
			+ name + "' is invalid");

	uint32 result = uint32(magnitude);
	return negative ? uint32(0) - result : result;
}

struct options {
	bool help;
	bool quiet;
	bool has_ihv;
	bool has_prefixfile;
	bool has_out;
	bool testmd5iv;
	bool testrndiv;
	bool testreciv;
	bool testall;
	std::string ihv;
	std::string prefixfile;
	std::vector<std::string> out;

	options()
		: help(false), quiet(false), has_ihv(false), has_prefixfile(false)
		, has_out(false), testmd5iv(false), testrndiv(false), testreciv(false)
		, testall(false)
	{}
};

// Returns the argument for an option that requires one, either the remainder of
// the current token or the next argument.
std::string take_value(const std::string& name, const std::string& rest,
	int argc, char** argv, int& i)
{
	if (!rest.empty())
		return rest;
	if (i + 1 >= argc)
		throw std::runtime_error("the required argument for option '--" + name
			+ "' is missing");
	return argv[++i];
}

void parse_command_line(int argc, char** argv, options& opt,
	uint32& seed32_1, uint32& seed32_2)
{
	bool have_positional = false;

	for (int i = 1; i < argc; ++i)
	{
		std::string arg = argv[i];

		if (arg.size() > 2 && arg.compare(0, 2, "--") == 0)
		{
			std::string name = arg.substr(2);
			std::string value;
			bool has_value = false;
			std::string::size_type eq = name.find('=');
			if (eq != std::string::npos)
			{
				value = name.substr(eq + 1);
				name = name.substr(0, eq);
				has_value = true;
			}

			if (name == "help")
				opt.help = true;
			else if (name == "quiet")
				opt.quiet = true;
			else if (name == "testmd5iv")
				opt.testmd5iv = true;
			else if (name == "testrndiv")
				opt.testrndiv = true;
			else if (name == "testreciv")
				opt.testreciv = true;
			else if (name == "testall")
				opt.testall = true;
			else if (name == "ihv")
			{
				opt.ihv = has_value ? value : take_value(name, "", argc, argv, i);
				opt.has_ihv = true;
			}
			else if (name == "prefixfile")
			{
				opt.prefixfile = has_value ? value : take_value(name, "", argc, argv, i);
				opt.has_prefixfile = true;
			}
			else if (name == "seed1")
				seed32_1 = parse_uint32(name, has_value ? value : take_value(name, "", argc, argv, i));
			else if (name == "seed2")
				seed32_2 = parse_uint32(name, has_value ? value : take_value(name, "", argc, argv, i));
			else if (name == "out")
			{
				// Multitoken: consume every following argument that is not an option.
				if (has_value)
					opt.out.push_back(value);
				else
					opt.out.push_back(take_value(name, "", argc, argv, i));
				while (i + 1 < argc && argv[i + 1][0] != '-')
					opt.out.push_back(argv[++i]);
				opt.has_out = true;
			}
			else
				throw std::runtime_error("unrecognised option '" + arg + "'");
		}
		else if (arg.size() > 1 && arg[0] == '-')
		{
			// Short options may be grouped; those taking an argument consume the
			// rest of the token, or the next argument if the token ends there.
			for (std::string::size_type c = 1; c < arg.size(); ++c)
			{
				char o = arg[c];
				std::string rest = arg.substr(c + 1);
				if (o == 'h')
					opt.help = true;
				else if (o == 'q')
					opt.quiet = true;
				else if (o == 'i')
				{
					opt.ihv = take_value("ihv", rest, argc, argv, i);
					opt.has_ihv = true;
					break;
				}
				else if (o == 'p')
				{
					opt.prefixfile = take_value("prefixfile", rest, argc, argv, i);
					opt.has_prefixfile = true;
					break;
				}
				else if (o == 'o')
				{
					opt.out.push_back(take_value("out", rest, argc, argv, i));
					while (i + 1 < argc && argv[i + 1][0] != '-')
						opt.out.push_back(argv[++i]);
					opt.has_out = true;
					break;
				}
				else
					throw std::runtime_error(std::string("unrecognised option '-")
						+ o + "'");
			}
		}
		else
		{
			// Sole positional argument is the prefixfile.
			if (have_positional || opt.has_prefixfile)
				throw std::runtime_error("too many positional options have been "
					"specified on the command line");
			opt.prefixfile = arg;
			opt.has_prefixfile = true;
			have_positional = true;
		}
	}
}


struct batch_task
{
	std::vector<std::array<uint32, 16> > prefix_blocks;
	uint32 iv[4];
	std::string out1, out2;
	std::vector<pid_t> attempts;
	bool done;

	batch_task(): done(false)
	{
		iv[0] = MD5IV[0]; iv[1] = MD5IV[1]; iv[2] = MD5IV[2]; iv[3] = MD5IV[3];
	}
};

struct collision_payload
{
	uint32 msg1block0[16], msg1block1[16], msg2block0[16], msg2block1[16];
};

struct batch_attempt
{
	pid_t pid;
	size_t task;
	int cpu;
	int result_fd;
};

bool load_batch_task(const std::string& prefix, const std::string& out1,
	const std::string& out2, batch_task& task)
{
	std::ifstream ifs(prefix.c_str(), std::ios::binary);
	if (!ifs)
		return false;

	task.out1 = out1;
	task.out2 = out2;
	uint32 block[16];
	while (load_block(ifs, block))
	{
		std::array<uint32, 16> saved;
		for (unsigned i = 0; i < 16; ++i)
			saved[i] = block[i];
		task.prefix_blocks.push_back(saved);
		md5_compress(task.iv, block);
	}
	return true;
}

bool write_all_fd(int fd, const void* data, size_t bytes)
{
	const char* p = static_cast<const char*>(data);
	while (bytes)
	{
		ssize_t n = ::write(fd, p, bytes);
		if (n < 0)
		{
			if (errno == EINTR) continue;
			return false;
		}
		p += n;
		bytes -= size_t(n);
	}
	return true;
}

bool read_all_fd(int fd, void* data, size_t bytes)
{
	char* p = static_cast<char*>(data);
	while (bytes)
	{
		ssize_t n = ::read(fd, p, bytes);
		if (n == 0) return false;
		if (n < 0)
		{
			if (errno == EINTR) continue;
			return false;
		}
		p += n;
		bytes -= size_t(n);
	}
	return true;
}

bool write_batch_output(const batch_task& task, const collision_payload& payload)
{
	std::ofstream ofs1(task.out1.c_str(), std::ios::binary);
	std::ofstream ofs2(task.out2.c_str(), std::ios::binary);
	if (!ofs1 || !ofs2)
		return false;
	for (size_t b = 0; b < task.prefix_blocks.size(); ++b)
	{
		save_block(ofs1, task.prefix_blocks[b].data());
		save_block(ofs2, task.prefix_blocks[b].data());
	}
	save_block(ofs1, payload.msg1block0);
	save_block(ofs1, payload.msg1block1);
	save_block(ofs2, payload.msg2block0);
	save_block(ofs2, payload.msg2block1);
	ofs1.close();
	ofs2.close();
	if (!ofs1 || !ofs2)
	{
		(void)unlink(task.out1.c_str());
		(void)unlink(task.out2.c_str());
		return false;
	}
	return true;
}

int run_batch_attempt(const batch_task& task, int cpu, unsigned serial, int result_fd)
{
	cpu_set_t one;
	CPU_ZERO(&one);
	CPU_SET(cpu, &one);
	(void)sched_setaffinity(0, sizeof(one), &one);

	const uint64 mix =
		(uint64(uint32(time(NULL))) << 32)
		^ uint64(uint32(getpid()))
		^ (uint64(serial) * 0x9e3779b97f4a7c15ULL);
	seed32_1 = uint32(mix ^ (mix >> 32) ^ 0xa5a5a5a5U);
	seed32_2 = uint32((mix >> 17) ^ (mix << 13) ^ 0x3c6ef372U);
	if ((seed32_1 | seed32_2) == 0)
		seed32_2 = 1;

	collision_payload payload;
	find_collision(task.iv,
		payload.msg1block0, payload.msg1block1,
		payload.msg2block0, payload.msg2block1, false);
	return write_all_fd(result_fd, &payload, sizeof(payload)) ? 0 : 2;
}

int run_batch(const std::string& tasks_file)
{
	std::ifstream in(tasks_file.c_str());
	if (!in)
	{
		std::cerr << "cannot open task file: " << tasks_file << std::endl;
		return 1;
	}

	std::vector<batch_task> tasks;
	std::string line;
	while (std::getline(in, line))
	{
		if (line.empty())
			continue;
		std::istringstream iss(line);
		std::string prefix, out1, out2, extra;
		if (!(iss >> prefix >> out1 >> out2) || (iss >> extra))
		{
			std::cerr << "malformed task line" << std::endl;
			return 1;
		}
		batch_task task;
		if (!load_batch_task(prefix, out1, out2, task))
		{
			std::cerr << "cannot read prefix: " << prefix << std::endl;
			return 1;
		}
		tasks.push_back(task);
	}
	if (tasks.empty())
		return 0;

	cpu_set_t allowed;
	CPU_ZERO(&allowed);
	if (sched_getaffinity(0, sizeof(allowed), &allowed) != 0)
	{
		perror("sched_getaffinity");
		return 1;
	}
	std::vector<int> free_cpus;
	for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
		if (CPU_ISSET(cpu, &allowed))
			free_cpus.push_back(cpu);
	if (free_cpus.empty())
		return 1;
	if (free_cpus.size() > 32)
		free_cpus.resize(32);

	std::map<pid_t, batch_attempt> active;
	size_t next_unique = 0;
	size_t completed = 0;
	unsigned serial = 0;

	auto remove_pid_from_task = [&](size_t ti, pid_t pid) {
		std::vector<pid_t>& a = tasks[ti].attempts;
		a.erase(std::remove(a.begin(), a.end(), pid), a.end());
	};

	auto start_attempt = [&](size_t ti) -> bool {
		if (free_cpus.empty())
			return false;
		const int cpu = free_cpus.back();
		free_cpus.pop_back();
		const unsigned my_serial = ++serial;

		int pipefd[2];
		if (pipe(pipefd) != 0)
		{
			free_cpus.push_back(cpu);
			return false;
		}

		pid_t pid = fork();
		if (pid < 0)
		{
			close(pipefd[0]);
			close(pipefd[1]);
			free_cpus.push_back(cpu);
			return false;
		}
		if (pid == 0)
		{
			close(pipefd[0]);
			const int rc = run_batch_attempt(tasks[ti], cpu, my_serial, pipefd[1]);
			close(pipefd[1]);
			_exit(rc);
		}
		close(pipefd[1]);

		batch_attempt a;
		a.pid = pid; a.task = ti; a.cpu = cpu; a.result_fd = pipefd[0];
		active[pid] = a;
		tasks[ti].attempts.push_back(pid);
		return true;
	};

	auto release_attempt = [&](pid_t pid, bool already_reaped, bool kill_first, bool cleanup) {
		std::map<pid_t, batch_attempt>::iterator it = active.find(pid);
		if (it == active.end())
			return;
		const batch_attempt a = it->second;
		if (kill_first)
			(void)kill(pid, SIGTERM);
		if (!already_reaped)
		{
			int st = 0;
			while (waitpid(pid, &st, 0) < 0 && errno == EINTR) {}
		}
		(void)cleanup;
		close(a.result_fd);
		remove_pid_from_task(a.task, pid);
		free_cpus.push_back(a.cpu);
		active.erase(it);
	};

	auto fill = [&]() {
		while (!free_cpus.empty() && next_unique < tasks.size())
		{
			start_attempt(next_unique);
			++next_unique;
		}
		if (next_unique < tasks.size())
			return;

		while (!free_cpus.empty() && completed < tasks.size())
		{
			size_t best = tasks.size();
			size_t best_copies = size_t(-1);
			for (size_t i = 0; i < tasks.size(); ++i)
			{
				if (tasks[i].done)
					continue;
				const size_t copies = tasks[i].attempts.size();
				if (copies < best_copies)
				{
					best = i;
					best_copies = copies;
				}
			}
			if (best == tasks.size() || !start_attempt(best))
				break;
		}
	};

	fill();
	while (completed < tasks.size())
	{
		int status = 0;
		pid_t pid;
		do { pid = waitpid(-1, &status, 0); } while (pid < 0 && errno == EINTR);
		if (pid < 0)
		{
			perror("waitpid");
			break;
		}

		std::map<pid_t, batch_attempt>::iterator it = active.find(pid);
		if (it == active.end())
			continue;
		const batch_attempt winner = it->second;
		batch_task& task = tasks[winner.task];

		if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
		{
			release_attempt(pid, true, false, true);
			fill();
			continue;
		}
		if (task.done)
		{
			release_attempt(pid, true, false, true);
			fill();
			continue;
		}

		collision_payload payload;
		if (!read_all_fd(winner.result_fd, &payload, sizeof(payload))
			|| !write_batch_output(task, payload))
		{
			release_attempt(pid, true, false, true);
			fill();
			continue;
		}

		task.done = true;
		++completed;
		const std::vector<pid_t> peers = task.attempts;
		release_attempt(pid, true, false, false);
		for (size_t i = 0; i < peers.size(); ++i)
			if (peers[i] != pid)
				release_attempt(peers[i], false, true, true);
		fill();
	}

	for (std::map<pid_t, batch_attempt>::const_iterator it = active.begin();
		it != active.end(); ++it)
		(void)kill(it->first, SIGTERM);
	while (!active.empty())
		release_attempt(active.begin()->first, false, false, true);

	return completed == tasks.size() ? 0 : 1;
}

} // namespace

void test_md5iv(bool single = false);
void test_rndiv(bool single = false);
void test_reciv(bool single = false);
void test_all();

int main(int argc, char** argv)
{
	if (argc == 3 && std::string(argv[1]) == "--batch")
		return run_batch(argv[2]);

	seed32_1 = uint32(time(NULL));
	seed32_2 = 0x12345678;

	uint32 IV[4] = { MD5IV[0], MD5IV[1], MD5IV[2], MD5IV[3] };

	string outfn1 = "msg1.bin";
	string outfn2 = "msg2.bin";
	string ihv;
	string prefixfn;
	bool verbose = true;

	cout <<
		"MD5 collision generator v1.5\n"
		"by Marc Stevens (http://www.win.tue.nl/hashclash/)\n"
		<< endl;

	try
	{
		hashclash::timer runtime(true);

		options opt;
		parse_command_line(argc, argv, opt, seed32_1, seed32_2);
		ihv = opt.ihv;
		prefixfn = opt.prefixfile;

		if (opt.quiet)
			verbose = false;

		if (opt.help || argc == 1) {
			cout << options_help << endl;
			return 1;
		}

		if (opt.testmd5iv)
			test_md5iv();
		if (opt.testrndiv)
			test_rndiv();
		if (opt.testreciv)
			test_reciv();
		if (opt.testall)
			test_all();

		if (opt.has_prefixfile)
		{
			unsigned l = prefixfn.size();
			if (l >= 4 && prefixfn[l-4]=='.' && prefixfn[l-3]!='.' && prefixfn[l-2]!='.' && prefixfn[l-1]!='.')
			{
				outfn1 = prefixfn.substr(0, l-4) + "_msg1" + prefixfn.substr(l-4);
				outfn2 = prefixfn.substr(0, l-4) + "_msg2" + prefixfn.substr(l-4);
				unsigned i = 1;
				while ( file_exists(outfn1)
					 || file_exists(outfn2))
				{
					outfn1 = prefixfn.substr(0, l-4) + "_msg1_" + uint_to_string(i) + prefixfn.substr(l-4);
					outfn2 = prefixfn.substr(0, l-4) + "_msg2_" + uint_to_string(i) + prefixfn.substr(l-4);
					++i;
				}
			}
		}

		if (opt.has_out)
		{
			vector<string>& outfns = opt.out;
			if (outfns.size() != 2)
			{
				cerr << "Error: exactly two output filenames should be specified." << endl;
				return 1;
			}
			outfn1 = outfns[0];
			outfn2 = outfns[1];
		}

		if (verbose)
			cout << "Using output filenames: '" << outfn1 << "' and '" << outfn2 << "'" << endl;
		if (opt.has_prefixfile)
		{
			if (verbose)
				cout << "Using prefixfile: '" << prefixfn << "'" << endl;
			ifstream ifs(prefixfn.c_str(), ios::binary);
			if (!ifs)
			{
				cerr << "Error: cannot open inputfile: '" << prefixfn << "'" << endl;
				return 1;
			}
			uint32 block[16];
			while (true)
			{
				unsigned len = load_block(ifs, block);
				if (len)
				{
					md5_compress(IV, block);
				} else
					break;
			}
		}
		else
		{
			if (!opt.has_ihv)
				ihv = "0123456789abcdeffedcba9876543210";
			if (ihv.size() != 32)
			{
				cerr << "Error: an initial value must be specified as a hash value of 32 hexadecimal characters." << endl;
				return 1;
			} else
			{
				uint32 c;
				for (unsigned i = 0; i < 4; ++i)
				{
					IV[i] = 0;
					for (unsigned b = 0; b < 4; ++b)
					{
						stringstream ss;
						ss << ihv.substr(i*8+b*2,2);
						ss >> hex >> c;					
						IV[i] += c << (b*8);
					}
				}

			}
		}
		if (verbose)
		{
			cout << "Using initial value: " << hex;
			unsigned oldwidth = cout.width(2);
			char oldfill = cout.fill('0');
			
			for (unsigned i = 0; i < 4; ++i)
			{
				for (unsigned b = 0; b < 4; ++b)
				{
					cout.width(2);
					cout.fill('0');
					cout << ((IV[i]>>(b*8))&0xFF);
				}
			}
			cout.width(oldwidth);
			cout.fill(oldfill);
			cout << dec << endl;
		}

		if (verbose)
			cout << endl;

		uint32 msg1block0[16];
		uint32 msg1block1[16];
		uint32 msg2block0[16];
		uint32 msg2block1[16];
		find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1, verbose);

		// Defer all output creation until a collision has actually been found.
		// Hedged searches that lose the race are normally terminated here, so
		// they no longer create and rewrite temporary prefix files.
		ofstream ofs1(outfn1.c_str(), ios::binary);
		if (!ofs1)
		{
			cerr << "Error: cannot open outfile: '" << outfn1 << "'" << endl;
			return 1;
		}
		ofstream ofs2(outfn2.c_str(), ios::binary);
		if (!ofs2)
		{
			cerr << "Error: cannot open outfile: '" << outfn2 << "'" << endl;
			return 1;
		}

		if (opt.has_prefixfile)
		{
			ifstream ifs(prefixfn.c_str(), ios::binary);
			if (!ifs)
			{
				cerr << "Error: cannot reopen inputfile: '" << prefixfn << "'" << endl;
				return 1;
			}
			uint32 block[16];
			while (load_block(ifs, block))
			{
				save_block(ofs1, block);
				save_block(ofs2, block);
			}
		}

		save_block(ofs1, msg1block0);
		save_block(ofs1, msg1block1);
		save_block(ofs2, msg2block0);
		save_block(ofs2, msg2block1);
		if (verbose)
			cout << "Running time: " << runtime.time() << " s" << endl;
		return 0;
	} catch (exception& e)
	{
		cerr << "\nException caught:\n" << e.what() << endl;
		return 1;
	} catch (...)
	{
		cerr << "\nUnknown exception caught!" << endl;
		return 1;
	}
}

void test_md5iv(bool single)
{
	uint32 IV[4] = { MD5IV[0], MD5IV[1], MD5IV[2], MD5IV[3] };
	uint32 msg1block0[16];
	uint32 msg1block1[16];
	uint32 msg2block0[16];
	uint32 msg2block1[16];

	hashclash::timer runtime(true);
	while (true)
	{
		runtime.start();
		find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1);
		double time = runtime.time();
		cout << endl << "Running time: " << time << " s" << endl;
		ofstream of_timings("timings_md5iv.txt", ios::app);
		of_timings << time << endl;
		if (single) return;
	}
}

void test_rndiv(bool single)
{
	uint32 IV[4];
	uint32 msg1block0[16];
	uint32 msg1block1[16];
	uint32 msg2block0[16];
	uint32 msg2block1[16];

	hashclash::timer runtime(true);
	while (true)
	{
		runtime.start();
		IV[0] = xrng64(); IV[1] = xrng64(); IV[2] = xrng64(); IV[3] = xrng64();
		find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1);
		double time = runtime.time();
		cout << endl << "Running time: " << time << " s" << endl;
		ofstream of_timings("timings_rndiv.txt", ios::app);
		of_timings << time << endl;
		if (single) return;
	}
}

void test_reciv(bool single)
{
	uint32 IV[4];
	uint32 msg1block0[16];
	uint32 msg1block1[16];
	uint32 msg2block0[16];
	uint32 msg2block1[16];

	hashclash::timer runtime(true);
	while (true)
	{
		runtime.start();
		IV[0] = xrng64(); IV[1] = xrng64(); IV[2] = xrng64(); IV[3] = xrng64();
		IV[2] |= 1<<25; IV[2] ^= ((IV[2] & (1<<24))<<1);
		IV[3] &= ~(1<<25); IV[3] ^= ((IV[3] & (1<<24))<<1);

		find_collision(IV, msg1block0, msg1block1, msg2block0, msg2block1);
		double time = runtime.time();
		cout << endl << "Running time: " << time << " s" << endl;
		ofstream of_timings("timings_reciv.txt", ios::app);
		of_timings << time << endl;
		if (single) return;
	}
}

void test_all()
{
	while (true)
	{
		test_md5iv(true);
		test_reciv(true);
		test_rndiv(true);
	}
}

#endif


unsigned load_block(istream& i, uint32 block[])
{
	unsigned len = 0;
	char uc;
	for (unsigned k = 0; k < 16; ++k)
	{
		block[k] = 0;
		for (unsigned c = 0; c < 4; ++c)
		{
			i.get(uc);
			if (i) 
				++len;
			else
				uc = 0;
			block[k] += uint32((unsigned char)(uc))<<(c*8);
		}
	}
	return len;
}

void save_block(ostream& o, const uint32 block[])
{
	for (unsigned k = 0; k < 16; ++k)
		for (unsigned c = 0; c < 4; ++c)
			o << (unsigned char)((block[k] >> (c*8))&0xFF);
}

void find_collision(const uint32 IV[], uint32 msg1block0[], uint32 msg1block1[], uint32 msg2block0[], uint32 msg2block1[], bool verbose)
{
	if (verbose)
		cout << "Generating first block: " << flush;
	find_block0(msg1block0, IV);

	uint32 IHV[4] = { IV[0], IV[1], IV[2], IV[3] };
	md5_compress(IHV, msg1block0);

	if (verbose)
		cout << endl << "Generating second block: " << flush;
	find_block1(msg1block1, IHV);

	for (int t = 0; t < 16; ++t)
	{
		msg2block0[t] = msg1block0[t];
		msg2block1[t] = msg1block1[t];
	}
	msg2block0[4] += 1 << 31; msg2block0[11] += 1 << 15; msg2block0[14] += 1 << 31;
	msg2block1[4] += 1 << 31; msg2block1[11] -= 1 << 15; msg2block1[14] += 1 << 31;
	if (verbose)
		cout << endl;
}
