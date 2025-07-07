#include "../../headers/utils/bench.h"

/*
 * The MIT License (MIT)
 * 
 * Copyright © 2025 Devadut S Balan
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the “Software”), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


benchmark_t benchmarks;

void benchmark_init() {
    benchmarks.total_time = 0;   
    benchmarks.timing_index = 0; 
    benchmarks.start_time = 0;
}

long long get_current_time_us() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (long long)tv.tv_sec * 1000000 + tv.tv_usec;
}

void record_timing(const char *function_name) {
    if (benchmarks.timing_index >= MAX_FUNS_TO_BENCH) {
        fprintf(stderr, "Error: Exceeded maximum number of benchmarked functions!\n");
        return;
    }

    long long end_time = get_current_time_us();
    benchmarks.timings[benchmarks.timing_index].time_us = end_time - benchmarks.start_time;
    benchmarks.total_time += benchmarks.timings[benchmarks.timing_index].time_us;

    strncpy(benchmarks.timings[benchmarks.timing_index].function_name, 
            function_name, 
            sizeof(benchmarks.timings[benchmarks.timing_index].function_name) - 1);
    
    benchmarks.timings[benchmarks.timing_index].function_name[
        sizeof(benchmarks.timings[benchmarks.timing_index].function_name) - 1] = '\0';

    benchmarks.timing_index++;
}

void format_time_us(long long time_us, char *buffer, size_t size) {
    if (time_us < 1000) {
        snprintf(buffer, size, "%4lld µs", time_us);
    } else if (time_us < 1000000) {
        snprintf(buffer, size, "%7.3f ms", time_us / 1000.0);
    } else {
        snprintf(buffer, size, "%7.3f s", time_us / 1000000.0);
    }
}

int compare_times(const void *a, const void *b) {
    return ((time_info*)b)->time_us - ((time_info*)a)->time_us;
}
const char* get_gradient_color(double percentage) {
    if (percentage >= 80.0) return BRIGHT_RED;   
    if (percentage >= 60.0) return RED;
    if (percentage >= 40.0) return MAGENTA;
    if (percentage >= 25.0) return BRIGHT_YELLOW;
    if (percentage >= 15.0) return YELLOW;
    if (percentage >= 5.0) return BRIGHT_GREEN;
    if (percentage > 0.1) return GREEN;
    return BLUE;  
}

void print_bench_ranked() {
    if (benchmarks.timing_index == 0) {
        fprintf(stdout, "\nNo benchmark data available.\n");
        return;
    }

    qsort(benchmarks.timings, benchmarks.timing_index, sizeof(time_info), compare_times);

    fprintf(stdout, "%s---------------------------------------------------------\n", BRIGHT_CYAN);
    fprintf(stdout, "| %-20s | %-12s | %-7s |\n", "Function", "Exec Time", "% of total runtime");
    fprintf(stdout, "---------------------------------------------------------%s\n", RESET);

    long long max_time = benchmarks.timings[0].time_us;

    for (size_t i = 0; i < benchmarks.timing_index; i++) {
        double percentage = (double)benchmarks.timings[i].time_us * 100.0 / benchmarks.total_time;
        int filled_length = (int)(BAR_LENGTH * percentage / 100.0);

        char time_str[15];
        format_time_us(benchmarks.timings[i].time_us, time_str, sizeof(time_str));

        const char *func_color = get_gradient_color((double)benchmarks.timings[i].time_us / max_time * 100.0);

        fprintf(stdout, "%s| %-20s | %12s | %6.4f%% |%s\n",
                func_color, 
                benchmarks.timings[i].function_name,
                time_str,
                percentage,
                RESET);

        printf("%s[", BRIGHT_CYAN);
        for (int j = 0; j < filled_length; j++) printf("▰");
        for (int j = 0; j < BAR_LENGTH - filled_length; j++) printf(" ");
        printf("]%s\n", RESET);
    }

    fprintf(stdout, "%s---------------------------------------------------------\n%s", BRIGHT_CYAN, RESET);
}


void print_bench_json() {
    fprintf(stdout, ">>>{\n");
    for (size_t i = 0; i < benchmarks.timing_index; i++) {
        double percentage = (double)benchmarks.timings[i].time_us * 100.0 / benchmarks.total_time;
        fprintf(stdout, "  \"%s\": {\"time_μs\": %lld, \"percentage\": %.2f}%s\n",
                benchmarks.timings[i].function_name,
                (long long)benchmarks.timings[i].time_us,
                percentage,
                (i < benchmarks.timing_index - 1) ? "," : "");
    }
    fprintf(stdout, "}<<<\n");
}


void print_bench() {
    for (size_t i = 0; i < benchmarks.timing_index; i++) {
        fprintf(stdout, "%s:%lld\n",benchmarks.timings[i].function_name,benchmarks.timings[i].time_us);
    }
}

