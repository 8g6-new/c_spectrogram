#ifndef BENCH_H
    #define BENCH_H

    #include <time.h>
    #include <stdio.h>
    #include <stddef.h>
    #include <stdlib.h>
    #include <string.h>
    #include <math.h>

    /**
    * @file bench.h
    * @brief Benchmarking and value-scaling utilities using general-purpose SI prefixes.
    * @author Devadut S Balan
    * @license MIT
    */

    #define MAX_FUNS_TO_BENCH      600     /**< Maximum number of benchmarked functions. */
    #define MAX_FUNS_NAME_LENGTH   100     /**< Maximum characters in a function label. */
    #define BAR_LENGTH             20      /**< Width of progress bars in ranking output. */

    // ─── ANSI Colors ──────────────────────────────────────────────────────────────
    #define RESET           "\x1b[0m"
    #define BRIGHT_RED      "\x1b[91m"
    #define RED             "\x1b[31m"
    #define MAGENTA         "\x1b[35m"
    #define BRIGHT_YELLOW   "\x1b[93m"
    #define YELLOW          "\x1b[33m"
    #define BRIGHT_GREEN    "\x1b[92m"
    #define GREEN           "\x1b[32m"
    #define BLUE            "\x1b[34m"
    #define BRIGHT_CYAN     "\x1b[96m"
    #define BRIGHT_BLUE     "\x1b[94m"
    #define BAR_COLOR       BRIGHT_BLUE

    /**
    * @brief Unicode bar used for visual separation in terminal.
    */
    static const char line[] =
    "\n▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰▰\n";

    /**
    * @brief SI unit prefix scale (for time, data, size, etc.)
    */
    typedef struct {
        const char *suffix;       /**< SI unit suffix ("n", "µ", "m", "", "k", etc.) */
        double scale_divisor;     /**< Numerical divisor to scale a raw value */
    } scale;

    /**
    * @brief Enum index for each supported SI scale.
    */
    typedef enum {
        scale_nano_idx = 0,   /**< nano (1e-9) */
        scale_micro_idx,      /**< micro (1e-6) */
        scale_milli_idx,      /**< milli (1e-3) */
        scale_unit_idx,       /**< base (1e0) */
        scale_kilo_idx,       /**< kilo (1e3) */
        scale_mega_idx,       /**< mega (1e6) */
        scale_giga_idx,       /**< giga (1e9) */
        scale_tera_idx,       /**< tera (1e12) */
        scale_count           /**< Total number of scales */
    } scale_index;

    /**
    * @brief General-purpose SI scaling table.
    *
    * These scales can be used for:
    * - Time (e.g., µs, ms)
    * - Memory (e.g., kB, MB)
    * - Frequency or rates (e.g., MFLOP/s, MHz)
    * - Any measurable quantity with power-of-ten representation
    */
    static const scale scales[scale_count] = {
        { "n", 1e-9 },
        { "µ", 1e-6 },
        { "m", 1e-3 },
        { "",  1    },
        { "k", 1e3  },
        { "M", 1e6  },
        { "G", 1e9  },
        { "T", 1e12 }
    };

    #define STRING_LENGTH 32

    /**
    * @brief Benchmark info for a single function.
    */
    typedef struct {
        long long time_us;                            /**< Execution time in microseconds */
        char function_name[MAX_FUNS_NAME_LENGTH];     /**< Human-readable function label */
    } time_info;

    /**
    * @brief Stores timing data for all benchmarked functions.
    */
    typedef struct {
        long long  total_time;                        /**< Total time recorded */
        long long  start_time;                        /**< Internal use only */
        time_info  timings[MAX_FUNS_TO_BENCH];        /**< Per-function timing info */
        size_t     timing_index;                      /**< Number of functions tracked */
    } benchmark_t;

    /**
    * @brief Global benchmark state.
    */
    extern benchmark_t benchmarks;

    /**
    * @brief Initialize the global benchmark state.
    */
    void benchmark_init();

    /**
    * @brief Returns current high-resolution time in microseconds.
    * @return Time in µs.
    */
    long long get_time_us();

    /**
    * @brief Records the time delta since `START_TIMING()` under a named function.
    * @param function_name Label to assign to the recorded time.
    */
    void record_timing(const char *function_name);

    /**
    * @brief Compare two time_info entries (descending).
    * @param a Pointer to first time_info.
    * @param b Pointer to second time_info.
    * @return Sorting order: negative if a < b.
    */
    int compare_times(const void *a, const void *b);

    /**
    * @brief Choose the most human-readable SI scale for a value.
    * @param v Raw value (e.g., 2500000.0).
    * @return Matching scale struct with suffix and divisor.
    */
    scale get_scale(double v);

    /**
    * @brief Format a scaled value into a human-readable string.
    * @param val       Raw value to scale.
    * @param buffer    Output string buffer.
    * @param buf_size  Maximum buffer length.
    * @param unit      Unit string (e.g., "s", "FLOP/s").
    */
    void format_scaled(double val, char *buffer, size_t buf_size, const char *unit);

    /**
    * @brief Estimate MFLOP/s for a given FFT and print timing and throughput.
    * @param mu_s      Average elapsed time per FFT in microseconds.
    * @param FFT_size  Number of FFT points (N).
    */
    void FFT_bench(double mu_s, unsigned int FFT_size);

    /**
    * @brief Choose a display color for ranking output based on percentage.
    * @param percentage Time or value percentage (0.0 to 100.0)
    * @return ANSI color string.
    */
    const char* get_gradient_color(double percentage);

    /**
    * @brief Print raw timings for all benchmarked functions.
    */
    void print_bench();

    /**
    * @brief Print all benchmark data as a JSON object.
    */
    void print_bench_json();

    /**
    * @brief Pretty-print ranked benchmark chart in table and bar format.
    */
    void print_bench_ranked();

    /**
    * @brief Start timing a block of code.
    */
    #define START_TIMING()           benchmarks.start_time = get_time_us()

    /**
    * @brief Stop timing and record elapsed duration.
    * @param FUNC_NAME Label to assign this timed block.
    */
    #define END_TIMING(FUNC_NAME)    record_timing(FUNC_NAME)

    /**
    * @brief Print an informational message with file, function, and line metadata.
    *
    * Outputs to stdout with a bright cyan [INFO] label and contextual metadata
    * including source file, function name, and line number.
    *
    * @param ... A printf-style format string followed by optional arguments.
    *
    * @example
    * LOG("Loaded %zu samples from %s", num_samples, filename);
    */
    #define LOG(...) \
        do { \
            fprintf(stdout, BRIGHT_CYAN "[INFO] " RESET "[file: %s | line: %d | func: %s] ", __FILE__, __LINE__, __func__); \
            fprintf(stdout, __VA_ARGS__); \
            fprintf(stdout, "\n"); \
        } while (0)

    /**
    * @brief Print a warning message with file, function, and line metadata.
    *
    * Outputs to stderr with a bright yellow [WARN] label and contextual metadata
    * including source file, function name, and line number.
    *
    * @param ... A printf-style format string followed by optional arguments.
    *
    * @example
    * WARN("Fallback to default window type: %s", fallback_type);
    */
    #define WARN(...) \
        do { \
            fprintf(stderr, BRIGHT_YELLOW "[WARN] " RESET "[file: %s | line: %d | func: %s] ", __FILE__, __LINE__, __func__); \
            fprintf(stderr, __VA_ARGS__); \
            fprintf(stderr, "\n"); \
        } while (0)

    /**
    * @brief Print an error message with file, function, and line metadata.
    *
    * Outputs to stderr with a bright red [ERROR] label and contextual metadata
    * including source file, function name, and line number.
    *
    * @param ... A printf-style format string followed by optional arguments.
    *
    * @example
    * ERROR("Memory allocation failed for %zu bytes", buffer_size);
    */
    #define ERROR(...) \
        do { \
            fprintf(stderr, BRIGHT_RED "[ERROR] " RESET "[file: %s | line: %d | func: %s] ", __FILE__, __LINE__, __func__); \
            fprintf(stderr, __VA_ARGS__); \
            fprintf(stderr, "\n"); \
        } while (0)


#endif  // BENCH_H
