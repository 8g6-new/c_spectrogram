#ifndef OPENCV_LIKE_SCHEME_H
#define OPENCV_LIKE_SCHEME_H

    #include <stdio.h>
    #include <stdint.h>
    #include <stddef.h>
    #include <stdbool.h>


    typedef struct {
        const unsigned char* data;
        size_t size;
        const short type;
    } color_data_t;

    #define DEF_CS(NAME, SIZE,index) \
        extern const unsigned char NAME##_data[SIZE]; \
        const color_data_t NAME = { (const unsigned char*)NAME##_data, (size_t)(SIZE/4), index }; /* Cast to unsigned char* */


    extern const color_data_t *cs[];

    extern const char *color_names[];

    typedef enum {
        AUTUMN,        BONE,          JET,            WINTER,
        RAINBOW,       OCEAN,         SUMMER,         SPRING,
        COOL,          HSV,           PINK,           HOT,
        PARULA,        MAGMA,         INFERNO,        PLASMA,
        VIRIDIS,       CIVIDIS,       TWILIGHT,       TWILIGHT_SHIFTED,
        TURBO,         DEEP_GREEN,    NUM_CS_MAX
    } cs_enum;

    
    char *fetch_color_opencv_like(cs_enum type,bool log);
    
    void print_all_cs();

#endif