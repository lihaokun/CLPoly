/**
 * @file bench_utils.hh
 * @brief Lightweight timing macros, memory measurement, and table output for benchmarks.
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
#include <sys/resource.h>
#include <sys/utsname.h>
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

// ---- RSS helpers ----

inline long get_current_rss_kb() {
    long rss = 0;
    std::ifstream status("/proc/self/status");
    if (status) {
        std::string line;
        while (std::getline(status, line))
            if (line.rfind("VmRSS:", 0) == 0) {
                std::string val;
                for (char c : line)
                    if (c >= '0' && c <= '9') val += c;
                if (!val.empty()) rss = std::stol(val);
                break;
            }
    }
    return rss;
}

inline long get_peak_rss_kb() {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss;
}

// Release freed heap memory back to OS for accurate RSS measurement
inline void flush_heap() { malloc_trim(0); }

inline void print_process_memory() {
    long cur = get_current_rss_kb();
    long peak = get_peak_rss_kb();
    std::cout << "\n---- Process Memory ----" << std::endl;
    std::cout << "  At exit:  " << std::setw(8) << cur << " KB"
              << "  (" << std::fixed << std::setprecision(1)
              << cur / 1024.0 << " MB)" << std::endl;
    std::cout << "  Peak:     " << std::setw(8) << peak << " KB"
              << "  (" << std::fixed << std::setprecision(1)
              << peak / 1024.0 << " MB)" << std::endl;
}

// ---- Memory measurement: per-object representation size ----
// Create N copies of an object via a lambda, measure heap delta, return per-object KB.
// Uses mallinfo2 for precise measurement (not RSS which has 4KB page granularity).

inline long get_heap_bytes() {
    struct mallinfo2 mi = mallinfo2();
    return static_cast<long>(mi.uordblks) + static_cast<long>(mi.hblkhd);
}

template<typename Fn>
double measure_per_object_kb(int N, Fn create_one) {
    long mem0 = get_heap_bytes();
    std::vector<decltype(create_one())> copies;
    copies.reserve(N);
    for (int i = 0; i < N; ++i)
        copies.push_back(create_one());
    long mem1 = get_heap_bytes();
    return (mem1 - mem0) / (1024.0 * N);
}

// Print table header for memory profile
inline void mem_profile_header() {
    std::cout << std::left << std::setw(38) << "Polynomial"
              << std::right << std::setw(8) << "Terms"
              << std::setw(12) << "KB/poly"
              << std::setw(14) << "Bytes/term" << std::endl;
    std::cout << std::string(72, '-') << std::endl;
}

inline void mem_profile_row(const std::string& name, size_t terms, double kb_per_poly) {
    double bytes_per_term = (terms > 0) ? kb_per_poly * 1024.0 / terms : 0;
    std::cout << std::left << std::setw(38) << name
              << std::right << std::setw(8) << terms
              << std::setw(12) << std::fixed << std::setprecision(2) << kb_per_poly
              << std::setw(14) << std::setprecision(0) << bytes_per_term << std::endl;
}

// Print table header for memory comparison
inline void mem_cmp_header(const std::string& other) {
    std::cout << std::left << std::setw(30) << "Polynomial"
              << std::right << std::setw(8) << "Terms"
              << std::setw(12) << "CLPoly"
              << std::setw(12) << other
              << std::setw(8) << "Ratio" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
}

inline void mem_cmp_row(const std::string& name, size_t terms,
                        double cl_kb, double ot_kb)
{
    double ratio = cl_kb / std::max(ot_kb, 0.001);
    std::cout << std::left << std::setw(30) << name
              << std::right << std::setw(8) << terms
              << std::setw(12) << std::fixed << std::setprecision(2) << cl_kb << " KB"
              << std::setw(12) << ot_kb << " KB"
              << std::setw(8) << std::setprecision(2) << ratio << "x" << std::endl;
}

// ---- Timing macros ----

#define BENCH(name, repeats, code) do { \
    auto _t0 = std::chrono::high_resolution_clock::now(); \
    for (int _r = 0; _r < (repeats); ++_r) { code; } \
    auto _t1 = std::chrono::high_resolution_clock::now(); \
    double _ms = std::chrono::duration<double,std::milli>(_t1-_t0).count() / (repeats); \
    std::cout << std::left << std::setw(48) << (name) \
              << std::right << std::setw(10) << std::fixed \
              << std::setprecision(3) << _ms << " ms" << std::endl; \
    _bench_results.push_back({name, _ms}); \
} while(0)

#define BENCH_CMP(name, repeats, clpoly_code, other_code) do { \
    auto _t0 = std::chrono::high_resolution_clock::now(); \
    for (int _r = 0; _r < (repeats); ++_r) { clpoly_code; } \
    auto _t1 = std::chrono::high_resolution_clock::now(); \
    double _ms_cl = std::chrono::duration<double,std::milli>(_t1-_t0).count() / (repeats); \
    _t0 = std::chrono::high_resolution_clock::now(); \
    for (int _r = 0; _r < (repeats); ++_r) { other_code; } \
    _t1 = std::chrono::high_resolution_clock::now(); \
    double _ms_ot = std::chrono::duration<double,std::milli>(_t1-_t0).count() / (repeats); \
    double _ratio = _ms_cl / std::max(_ms_ot, 0.001); \
    std::cout << std::left << std::setw(40) << (name) \
              << std::right << std::setw(10) << std::fixed << std::setprecision(3) << _ms_cl << " ms" \
              << std::setw(10) << _ms_ot << " ms" \
              << std::setw(8) << std::setprecision(2) << _ratio << "x" << std::endl; \
} while(0)

inline void bench_header(const std::string& title) {
    std::cout << "\n==== " << title << " ====" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
}

inline void bench_cmp_header(const std::string& title, const std::string& other) {
    std::cout << "\n==== " << title << " ====" << std::endl;
    std::cout << std::left << std::setw(40) << "Operation"
              << std::right << std::setw(10) << "CLPoly"
              << std::setw(10) << other
              << std::setw(8) << "Ratio" << std::endl;
    std::cout << std::string(68, '-') << std::endl;
}

#endif // BENCH_UTILS_HH
