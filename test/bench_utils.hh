/**
 * @file bench_utils.hh
 * @brief Fork-based benchmark framework: every benchmark runs in a child
 *        process, reporting both wall-clock time and heap-memory delta.
 */
#ifndef BENCH_UTILS_HH
#define BENCH_UTILS_HH

#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/utsname.h>
#include <sys/wait.h>
#include <unistd.h>
#include <signal.h>
#include <malloc.h>

struct BenchResult { std::string name; double ms; };
static std::vector<BenchResult> _bench_results;

// ---- Machine info ----

inline void print_sysinfo() {
    std::cout << "==== System Info ====" << std::endl;

    struct utsname uts;
    if (uname(&uts) == 0)
        std::cout << "OS:        " << uts.sysname << " " << uts.release
                  << " " << uts.machine << std::endl;

    std::ifstream cpuinfo("/proc/cpuinfo");
    if (cpuinfo) {
        std::string line;
        while (std::getline(cpuinfo, line))
            if (line.rfind("model name", 0) == 0) {
                auto pos = line.find(':');
                if (pos != std::string::npos)
                    std::cout << "CPU:       " << line.substr(pos + 2) << std::endl;
                break;
            }
    }

    std::ifstream meminfo("/proc/meminfo");
    if (meminfo) {
        std::string line;
        while (std::getline(meminfo, line))
            if (line.rfind("MemTotal:", 0) == 0) {
                std::cout << "RAM:       " << line.substr(9) << std::endl;
                break;
            }
    }

#ifdef __VERSION__
    std::cout << "Compiler:  g++ " << __VERSION__ << std::endl;
#endif
#ifdef NDEBUG
    std::cout << "Build:     Release (-O3 -DNDEBUG)" << std::endl;
#else
    std::cout << "Build:     Debug" << std::endl;
#endif
    std::cout << std::string(70, '=') << std::endl;
}

// ---- Heap measurement ----

inline long get_heap_bytes() {
    struct mallinfo2 mi = mallinfo2();
    return static_cast<long>(mi.uordblks) + static_cast<long>(mi.hblkhd);
}

// ---- Fork-based execution core ----

struct ForkResult {
    double ms;        // average per-iteration ms, -1.0 = failed
    long   mem_bytes; // heap delta in bytes
    int    status;    // 0=ok, 1=timeout, 2=crash
};

// Payload written by the child back through the pipe
struct _ForkPayload {
    double ms;
    long   mem_bytes;
};

template<typename Fn>
ForkResult fork_run(Fn fn, int repeats, double timeout_sec) {
    int pipefd[2];
    if (pipe(pipefd) != 0) return {-1.0, 0, 2};

    pid_t pid = fork();
    if (pid < 0) { close(pipefd[0]); close(pipefd[1]); return {-1.0, 0, 2}; }

    if (pid == 0) {
        // ---- child ----
        close(pipefd[0]);
        long mem0 = get_heap_bytes();
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int r = 0; r < repeats; ++r) fn();
        auto t1 = std::chrono::high_resolution_clock::now();
        long mem1 = get_heap_bytes();
        _ForkPayload pay;
        pay.ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / repeats;
        pay.mem_bytes = mem1 - mem0;
        ssize_t wr __attribute__((unused)) = write(pipefd[1], &pay, sizeof(pay));
        close(pipefd[1]);
        _exit(0);
    }

    // ---- parent ----
    close(pipefd[1]);

    auto deadline = std::chrono::steady_clock::now()
                  + std::chrono::milliseconds(static_cast<long long>(timeout_sec * 1000));

    while (true) {
        int wstatus;
        pid_t ret = waitpid(pid, &wstatus, WNOHANG);
        if (ret > 0) {
            _ForkPayload pay{-1.0, 0};
            if (WIFEXITED(wstatus) && WEXITSTATUS(wstatus) == 0) {
                ssize_t rd __attribute__((unused)) = read(pipefd[0], &pay, sizeof(pay));
            }
            close(pipefd[0]);
            bool ok = (WIFEXITED(wstatus) && WEXITSTATUS(wstatus) == 0 && pay.ms >= 0);
            return {pay.ms, pay.mem_bytes, ok ? 0 : 2};
        }
        if (std::chrono::steady_clock::now() >= deadline) {
            kill(pid, SIGKILL);
            waitpid(pid, nullptr, 0);
            close(pipefd[0]);
            return {-1.0, 0, 1};
        }
        usleep(5000);
    }
}

// Parent-side wait helper: polls child with timeout, reads T from pipe.
// Returns true on success (child exited 0 and wrote result), false otherwise.
template<typename T>
bool _fork_wait(pid_t pid, int read_fd, double timeout_sec, T& out) {
    auto deadline = std::chrono::steady_clock::now()
                  + std::chrono::milliseconds(static_cast<long long>(timeout_sec * 1000));
    while (true) {
        int wstatus;
        pid_t ret = waitpid(pid, &wstatus, WNOHANG);
        if (ret > 0) {
            bool ok = false;
            if (WIFEXITED(wstatus) && WEXITSTATUS(wstatus) == 0) {
                ssize_t rd __attribute__((unused)) = read(read_fd, &out, sizeof(T));
                ok = true;
            }
            close(read_fd);
            return ok;
        }
        if (std::chrono::steady_clock::now() >= deadline) {
            kill(pid, SIGKILL);
            waitpid(pid, nullptr, 0);
            close(read_fd);
            return false;
        }
        usleep(5000);
    }
}

// Measure per-object KB: create N copies in a child, return KB/obj.
// Uses BENCH_MEM_INNER macro pattern — objects stay alive during measurement.
template<typename Fn>
double fork_mem(Fn create_one, int N, double timeout_sec) {
    int pipefd[2];
    if (pipe(pipefd) != 0) return -1.0;
    pid_t pid = fork();
    if (pid < 0) { close(pipefd[0]); close(pipefd[1]); return -1.0; }
    if (pid == 0) {
        close(pipefd[0]);
        long mem0 = get_heap_bytes();
        // Allocate N copies — they stay alive until _exit
        std::vector<decltype(create_one())> copies;
        copies.reserve(N);
        for (int i = 0; i < N; ++i)
            copies.push_back(create_one());
        long mem1 = get_heap_bytes();
        long delta = mem1 - mem0;
        ssize_t wr __attribute__((unused)) = write(pipefd[1], &delta, sizeof(delta));
        close(pipefd[1]);
        _exit(0);
    }
    close(pipefd[1]);
    long delta = 0;
    if (!_fork_wait(pid, pipefd[0], timeout_sec, delta)) return -1.0;
    return delta / (1024.0 * N);
}

// Raw memory delta (bytes) for an arbitrary setup lambda.
// The setup_fn must keep allocations alive until it returns.
// Since we measure AROUND the call (mem0 before, mem1 after), and the
// child _exit()s immediately, objects ARE destroyed before mem1 is read.
// Solution: we use a macro (BENCH_MEM_CMP) that inlines the measurement.
// This function is kept for cases where the lambda blocks (e.g., stores
// into a pre-declared container captured by reference).
template<typename Fn>
double fork_mem_raw(Fn setup_fn, double timeout_sec) {
    int pipefd[2];
    if (pipe(pipefd) != 0) return -1.0;
    pid_t pid = fork();
    if (pid < 0) { close(pipefd[0]); close(pipefd[1]); return -1.0; }
    if (pid == 0) {
        close(pipefd[0]);
        long mem0 = get_heap_bytes();
        setup_fn();
        long mem1 = get_heap_bytes();
        long delta = mem1 - mem0;
        ssize_t wr __attribute__((unused)) = write(pipefd[1], &delta, sizeof(delta));
        close(pipefd[1]);
        _exit(0);
    }
    close(pipefd[1]);
    long delta = 0;
    if (!_fork_wait(pid, pipefd[0], timeout_sec, delta)) return -1.0;
    return static_cast<double>(delta);
}

// ---- Formatting helpers ----

// Format milliseconds for display (right-aligned into a fixed width)
inline std::string _fmt_ms(double ms) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << ms << "ms";
    return oss.str();
}

// Format KB for display
inline std::string _fmt_kb(long mem_bytes) {
    double kb = mem_bytes / 1024.0;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << kb << "KB";
    return oss.str();
}

// ---- Table headers ----

inline void bench_header(const std::string& title) {
    std::cout << "\n==== " << title << " ====" << std::endl;
    std::cout << std::left << std::setw(48) << "Operation"
              << std::right << std::setw(12) << "Time"
              << std::setw(12) << "Mem" << std::endl;
    std::cout << std::string(72, '-') << std::endl;
}

inline void bench_cmp_header(const std::string& title, const std::string& other) {
    std::cout << "\n==== " << title << " ====" << std::endl;
    std::cout << std::left << std::setw(28) << "Operation"
              << std::right
              << std::setw(10) << "CLPoly" << std::setw(10) << "Mem"
              << std::setw(10) << other   << std::setw(10) << "Mem"
              << std::setw(8)  << "Time"  << std::setw(8)  << "Mem"
              << std::endl;
    // Subheader row with units
    std::cout << std::left << std::setw(28) << ""
              << std::right
              << std::setw(10) << "Time" << std::setw(10) << ""
              << std::setw(10) << "Time" << std::setw(10) << ""
              << std::setw(8)  << "Ratio" << std::setw(8) << "Ratio"
              << std::endl;
    std::cout << std::string(82, '-') << std::endl;
}

// ---- BENCH macro: single-side fork, time + mem ----

#define BENCH(name, repeats, code) do { \
    auto _fr = fork_run([&]() { code; }, (repeats), 300.0); \
    if (_fr.status == 1) { \
        std::cout << std::left << std::setw(48) << (name) \
                  << std::right << std::setw(12) << "TIMEOUT" \
                  << std::setw(12) << "\xe2\x80\x94" << std::endl; \
    } else if (_fr.status != 0) { \
        std::cout << std::left << std::setw(48) << (name) \
                  << std::right << std::setw(12) << "CRASH" \
                  << std::setw(12) << "\xe2\x80\x94" << std::endl; \
    } else { \
        std::cout << std::left << std::setw(48) << (name) \
                  << std::right << std::setw(12) << _fmt_ms(_fr.ms) \
                  << std::setw(12) << _fmt_kb(_fr.mem_bytes) << std::endl; \
        _bench_results.push_back({name, _fr.ms}); \
    } \
} while(0)

// ---- BENCH_CMP macro: both sides fork, time + mem + ratio ----

#define BENCH_CMP(name, repeats, timeout, clpoly_code, other_code) do { \
    auto _fr_cl = fork_run([&]() { clpoly_code; }, (repeats), (timeout)); \
    auto _fr_ot = fork_run([&]() { other_code; }, (repeats), (timeout)); \
    std::cout << std::left << std::setw(28) << (name); \
    /* CLPoly columns */ \
    if (_fr_cl.status == 0) { \
        std::cout << std::right << std::setw(10) << _fmt_ms(_fr_cl.ms) \
                  << std::setw(10) << _fmt_kb(_fr_cl.mem_bytes); \
    } else { \
        std::cout << std::right << std::setw(10) \
                  << (_fr_cl.status == 1 ? "TIMEOUT" : "CRASH") \
                  << std::setw(10) << "\xe2\x80\x94"; \
    } \
    /* Other columns */ \
    if (_fr_ot.status == 0) { \
        std::cout << std::right << std::setw(10) << _fmt_ms(_fr_ot.ms) \
                  << std::setw(10) << _fmt_kb(_fr_ot.mem_bytes); \
    } else { \
        std::cout << std::right << std::setw(10) \
                  << (_fr_ot.status == 1 ? "TIMEOUT" : "CRASH") \
                  << std::setw(10) << "\xe2\x80\x94"; \
    } \
    /* Ratio columns */ \
    if (_fr_cl.status == 0 && _fr_ot.status == 0) { \
        double _t_ratio = _fr_cl.ms / std::max(_fr_ot.ms, 0.001); \
        double _m_cl_kb = _fr_cl.mem_bytes / 1024.0; \
        double _m_ot_kb = _fr_ot.mem_bytes / 1024.0; \
        std::ostringstream _tr; _tr << std::fixed << std::setprecision(2) << _t_ratio << "x"; \
        std::cout << std::right << std::setw(8) << _tr.str(); \
        if (_m_ot_kb > 0.01) { \
            double _m_ratio = _m_cl_kb / _m_ot_kb; \
            std::ostringstream _mr; _mr << std::fixed << std::setprecision(2) << _m_ratio << "x"; \
            std::cout << std::setw(8) << _mr.str(); \
        } else { \
            std::cout << std::setw(8) << "\xe2\x80\x94"; \
        } \
    } else { \
        std::cout << std::right << std::setw(8) << "\xe2\x80\x94" \
                  << std::setw(8) << "\xe2\x80\x94"; \
    } \
    std::cout << std::endl; \
} while(0)

// ---- BENCH_MEM macro: per-polynomial memory profiling ----

inline void _bench_mem_header_once() {
    static bool printed = false;
    if (!printed) {
        std::cout << std::left << std::setw(38) << "Polynomial"
                  << std::right << std::setw(8) << "Terms"
                  << std::setw(12) << "KB/poly"
                  << std::setw(14) << "Bytes/term" << std::endl;
        std::cout << std::string(72, '-') << std::endl;
        printed = true;
    }
}

#define BENCH_MEM(name, terms, N, create_one_expr) do { \
    double _kb = fork_mem([&]() { return (create_one_expr); }, (N), 60.0); \
    size_t _terms = (terms); \
    double _bpt = (_terms > 0 && _kb > 0) ? _kb * 1024.0 / _terms : 0; \
    std::cout << std::left << std::setw(38) << (name) \
              << std::right << std::setw(8) << _terms; \
    if (_kb >= 0) { \
        std::cout << std::setw(12) << std::fixed << std::setprecision(2) << _kb \
                  << std::setw(14) << std::setprecision(0) << _bpt; \
    } else { \
        std::cout << std::setw(12) << "TIMEOUT" << std::setw(14) << "\xe2\x80\x94"; \
    } \
    std::cout << std::endl; \
} while(0)

// ---- BENCH_MEM_CMP macro: memory comparison with other library ----

inline void _bench_mem_cmp_header_once(const std::string& other) {
    static bool printed = false;
    if (!printed) {
        std::cout << std::left << std::setw(30) << "Polynomial"
                  << std::right << std::setw(8) << "Terms"
                  << std::setw(12) << "CLPoly"
                  << std::setw(12) << other
                  << std::setw(8) << "Ratio" << std::endl;
        std::cout << std::string(70, '-') << std::endl;
        printed = true;
    }
}

#define BENCH_MEM_CMP(name, terms, N, cl_expr, ot_expr) do { \
    double _cl_kb = fork_mem([&]() { return (cl_expr); }, (N), 60.0); \
    double _ot_raw = fork_mem_raw([&]() { ot_expr; }, 60.0); \
    double _ot_kb = (_ot_raw >= 0) ? _ot_raw / (1024.0 * (N)) : -1.0; \
    double _ratio = (_cl_kb > 0 && _ot_kb > 0.001) ? _cl_kb / _ot_kb : 0; \
    std::cout << std::left << std::setw(30) << (name) \
              << std::right << std::setw(8) << (terms); \
    if (_cl_kb >= 0) { \
        std::cout << std::setw(10) << std::fixed << std::setprecision(2) << _cl_kb << "KB"; \
    } else { \
        std::cout << std::setw(12) << "TIMEOUT"; \
    } \
    if (_ot_kb >= 0) { \
        std::cout << std::setw(10) << std::fixed << std::setprecision(2) << _ot_kb << "KB"; \
    } else { \
        std::cout << std::setw(12) << "TIMEOUT"; \
    } \
    if (_ratio > 0) { \
        std::cout << std::setw(8) << std::setprecision(2) << _ratio << "x"; \
    } else { \
        std::cout << std::setw(8) << "\xe2\x80\x94"; \
    } \
    std::cout << std::endl; \
} while(0)

#endif // BENCH_UTILS_HH
