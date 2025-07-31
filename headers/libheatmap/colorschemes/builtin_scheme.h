#ifndef BUILTIN_SCHEME_H
#define BUILTIN_SCHEME_H
    #include <stdio.h>
    #include <stddef.h>
    #include <stdio.h>
    #include <stdlib.h>
    #include <stdbool.h>
    #include <string.h>


    // cs = color scheme

    typedef struct {
        const unsigned char* data; 
        size_t size;                
        const short int type;        
    } color_data_t;

    typedef struct {
        color_data_t discrete;
        color_data_t soft;
        color_data_t mixed;
        color_data_t mixed_exp;
    } color_data_struct;

    #define DEF_CS_REF(NAME, D_size, S_size, M_size, ME_size) \
        color_data_struct NAME = { \
            .discrete = { .data = NAME##_discrete_data, .size = (size_t)(D_size), .type = 0 }, \
            .soft = { .data = NAME##_soft_data, .size = (size_t)(S_size), .type = 1 }, \
            .mixed = { .data = NAME##_mixed_data, .size = (size_t)(M_size), .type = 2 }, \
            .mixed_exp = { .data = NAME##_mixed_exp_data, .size = (size_t)(ME_size), .type = 3 } \
        };

    #define DEF_CS(main_color, D_size, S_size, M_size, ME_size) \
        extern const unsigned char main_color##_discrete_data[D_size]; \
        extern const unsigned char main_color##_soft_data[S_size]; \
        extern const unsigned char main_color##_mixed_data[M_size]; \
        extern const unsigned char main_color##_mixed_exp_data[ME_size]; \
        DEF_CS_REF(main_color, (size_t)(D_size / 4), (size_t)(S_size / 4), (size_t)(M_size / 4), (size_t)(ME_size / 4))
    
    /*
     macro example output
    DEF_CS_EXT(Reds,48, 4100, 4100, 4100) => extern const unsigned char Reds_discrete_data[48];
                                                extern const unsigned char Reds_soft_data[4100];
                                                extern const unsigned char Reds_mixed_data[4100];
                                                extern const unsigned char Reds_mixed_exp_data[4100];
    */


    /*
    macro example output
        DEF_CS_REF(Reds, 10,1025,1025,1025); => color_data_struct Reds = {
            .discrete = { .data = Reds_discrete_data, .size = 10, .type = 0 },
            .soft = { .data = Reds_soft_data, .size = 1025, .type = 1},
            .mixed = { .data = Reds_mixed_data, .size = 1025, .type = 2 },
            .mixed_exp = { .data = Reds_mixed_exp_data, .size = 1025, .type = 3 }
        };
    */


    typedef enum {
        PRGn_discrete, PRGn_soft, PRGn_mixed, PRGn_mixed_exp,
        Blues_discrete, Blues_soft, Blues_mixed, Blues_mixed_exp,
        BrBG_discrete, BrBG_soft, BrBG_mixed, BrBG_mixed_exp,
        BuGn_discrete, BuGn_soft, BuGn_mixed, BuGn_mixed_exp,
        BuPu_discrete, BuPu_soft, BuPu_mixed, BuPu_mixed_exp,
        GnBu_discrete, GnBu_soft, GnBu_mixed, GnBu_mixed_exp,
        Greens_discrete, Greens_soft, Greens_mixed, Greens_mixed_exp,
        Greys_discrete, Greys_soft, Greys_mixed, Greys_mixed_exp,
        Oranges_discrete, Oranges_soft, Oranges_mixed, Oranges_mixed_exp,
        OrRd_discrete, OrRd_soft, OrRd_mixed, OrRd_mixed_exp,
        PiYG_discrete, PiYG_soft, PiYG_mixed, PiYG_mixed_exp,
        PuBu_discrete, PuBu_soft, PuBu_mixed, PuBu_mixed_exp,
        PuBuGn_discrete, PuBuGn_soft, PuBuGn_mixed, PuBuGn_mixed_exp,
        PuOr_discrete, PuOr_soft, PuOr_mixed, PuOr_mixed_exp,
        PuRd_discrete, PuRd_soft, PuRd_mixed, PuRd_mixed_exp,
        Purples_discrete, Purples_soft, Purples_mixed, Purples_mixed_exp,
        RdBu_discrete, RdBu_soft, RdBu_mixed, RdBu_mixed_exp,
        RdGy_discrete, RdGy_soft, RdGy_mixed, RdGy_mixed_exp,
        RdPu_discrete, RdPu_soft, RdPu_mixed, RdPu_mixed_exp,
        RdYlBu_discrete, RdYlBu_soft, RdYlBu_mixed, RdYlBu_mixed_exp,
        RdYlGn_discrete, RdYlGn_soft, RdYlGn_mixed, RdYlGn_mixed_exp,
        Reds_discrete, Reds_soft, Reds_mixed, Reds_mixed_exp,
        Spectral_discrete, Spectral_soft, Spectral_mixed, Spectral_mixed_exp,
        YlGn_discrete, YlGn_soft, YlGn_mixed, YlGn_mixed_exp,
        YlGnBu_discrete, YlGnBu_soft, YlGnBu_mixed, YlGnBu_mixed_exp,
        YlOrBr_discrete, YlOrBr_soft, YlOrBr_mixed, YlOrBr_mixed_exp,
        YlOrRd_discrete, YlOrRd_soft, YlOrRd_mixed, YlOrRd_mixed_exp,NUM_CS_MAX,
    } cs_enum;
   

    extern const char *color_subtypes[];
    extern const char *color_names[];

    extern const color_data_t *cs[];

    char *fetch_color_builtin(cs_enum type,bool log);
    void print_all_cs();
    
#endif