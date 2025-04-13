#ifndef BENCH_H
#define BENCH_H
    
    #include <stdio.h>
    #include <stddef.h> 
    #include <stdlib.h>
    #include <string.h>

    #include <sys/time.h>

    #define MAX_FUNS_TO_BENCH 600
    #define MAX_FUNS_NAME_LENGTH 100
    #define BAR_LENGTH 20

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

    #define BAR_COLOR       BRIGHT_BLUE


    typedef struct {
        long long time_us;
        char function_name[MAX_FUNS_NAME_LENGTH];
    } time_info;

    typedef struct {
        long long  total_time;
        long long  start_time;
        time_info  timings[MAX_FUNS_TO_BENCH];
        size_t     timing_index;
    } benchmark_t;


    long long               get_current_time_us();
    void                    record_timing(const char *function_name);
    void                    init_benchmark();
    void                    print_bench_ranked(); 
    void                    print_bench_json();
    int                     compare_times(const void *a, const void *b);
    void                    print_bench();

    extern benchmark_t benchmarks;

    #define START_TIMING() benchmarks.start_time = get_current_time_us()
    #define END_TIMING(FUNC_NAME) record_timing(FUNC_NAME)

#endif