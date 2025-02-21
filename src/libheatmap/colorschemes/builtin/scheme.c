#include "../../../../headers/libheatmap/colorschemes/builtin_scheme.h"

DEF_CS(PRGn, 48, 4100, 4100, 4100)
DEF_CS(Blues, 40, 4100, 4100, 4100)
DEF_CS(BrBG, 48, 4100, 4100, 4100)
DEF_CS(BuGn, 40, 4100, 4100, 4100)
DEF_CS(BuPu, 40, 4100, 4100, 4100)
DEF_CS(GnBu, 40, 4100, 4100, 4100)
DEF_CS(Greens, 40, 4100, 4100, 4100)
DEF_CS(Greys, 40, 4100, 4100, 4100)
DEF_CS(Oranges, 40, 4100, 4100, 4100)
DEF_CS(OrRd, 40, 4100, 4100, 4100)
DEF_CS(PiYG, 48, 4100, 4100, 4100)
DEF_CS(PuBu, 40, 4100, 4100, 4100)
DEF_CS(PuBuGn, 40, 4100, 4100, 4100)
DEF_CS(PuOr, 48, 4100, 4100, 4100)
DEF_CS(PuRd, 40, 4100, 4100, 4100)
DEF_CS(Purples, 40, 4100, 4100, 4100)
DEF_CS(RdBu, 48, 4100, 4100, 4100)
DEF_CS(RdGy, 48, 4100, 4100, 4100)
DEF_CS(RdPu, 40, 4100, 4100, 4100)
DEF_CS(RdYlBu, 48, 4100, 4100, 4100)
DEF_CS(RdYlGn, 44, 4100, 4100, 4100)
DEF_CS(Reds, 40, 4100, 4100, 4100)
DEF_CS(Spectral, 48, 4100, 4100, 4100)
DEF_CS(YlGn, 40, 4100, 4100, 4100)
DEF_CS(YlGnBu, 40, 4100, 4100, 4100)
DEF_CS(YlOrBr, 40, 4100, 4100, 4100)
DEF_CS(YlOrRd, 40, 4100, 4100, 4100)


const color_data_t *cs[] = {
    &PRGn.discrete, &PRGn.soft, &PRGn.mixed, &PRGn.mixed_exp,
    &Blues.discrete, &Blues.soft, &Blues.mixed, &Blues.mixed_exp,
    &BrBG.discrete, &BrBG.soft, &BrBG.mixed, &BrBG.mixed_exp,
    &BuGn.discrete, &BuGn.soft, &BuGn.mixed, &BuGn.mixed_exp,
    &BuPu.discrete, &BuPu.soft, &BuPu.mixed, &BuPu.mixed_exp,
    &GnBu.discrete, &GnBu.soft, &GnBu.mixed, &GnBu.mixed_exp,
    &Greens.discrete, &Greens.soft, &Greens.mixed, &Greens.mixed_exp,
    &Greys.discrete, &Greys.soft, &Greys.mixed, &Greys.mixed_exp,
    &Oranges.discrete, &Oranges.soft, &Oranges.mixed, &Oranges.mixed_exp,
    &OrRd.discrete, &OrRd.soft, &OrRd.mixed, &OrRd.mixed_exp,
    &PiYG.discrete, &PiYG.soft, &PiYG.mixed, &PiYG.mixed_exp,
    &PuBu.discrete, &PuBu.soft, &PuBu.mixed, &PuBu.mixed_exp,
    &PuBuGn.discrete, &PuBuGn.soft, &PuBuGn.mixed, &PuBuGn.mixed_exp,
    &PuOr.discrete, &PuOr.soft, &PuOr.mixed, &PuOr.mixed_exp,
    &PuRd.discrete, &PuRd.soft, &PuRd.mixed, &PuRd.mixed_exp,
    &Purples.discrete, &Purples.soft, &Purples.mixed, &Purples.mixed_exp,
    &RdBu.discrete, &RdBu.soft, &RdBu.mixed, &RdBu.mixed_exp,
    &RdGy.discrete, &RdGy.soft, &RdGy.mixed, &RdGy.mixed_exp,
    &RdPu.discrete, &RdPu.soft, &RdPu.mixed, &RdPu.mixed_exp,
    &RdYlBu.discrete, &RdYlBu.soft, &RdYlBu.mixed, &RdYlBu.mixed_exp,
    &RdYlGn.discrete, &RdYlGn.soft, &RdYlGn.mixed, &RdYlGn.mixed_exp,
    &Reds.discrete, &Reds.soft, &Reds.mixed, &Reds.mixed_exp,
    &Spectral.discrete, &Spectral.soft, &Spectral.mixed, &Spectral.mixed_exp,
    &YlGn.discrete, &YlGn.soft, &YlGn.mixed, &YlGn.mixed_exp,
    &YlGnBu.discrete, &YlGnBu.soft, &YlGnBu.mixed, &YlGnBu.mixed_exp,
    &YlOrBr.discrete, &YlOrBr.soft, &YlOrBr.mixed, &YlOrBr.mixed_exp,
    &YlOrRd.discrete, &YlOrRd.soft, &YlOrRd.mixed, &YlOrRd.mixed_exp
};

const char *color_subtypes[] = {
    "discrete", "soft", "mixed", "mixed_exp"
};

const char *main_cs[] = {
    "PRGn", "Blues", "BrBG", "BuGn", "BuPu", "GnBu", "Greens", "Greys",
    "Oranges", "OrRd", "PiYG", "PuBu", "PuBuGn", "PuOr", "PuRd", "Purples",
    "RdBu", "RdGy", "RdPu", "RdYlBu", "RdYlGn", "Reds", "Spectral", "YlGn",
    "YlGnBu", "YlOrBr", "YlOrRd"
};

#define NUM_SCHEMES (sizeof(main_cs) / sizeof(main_cs[0]))

char *fetch_color_builtin(cs_enum type,bool log){
  if(log)
        printf("enum %d => libheatmap builtin color scheme : main type %s , subtype %s",type,main_cs[(int) type/4],color_subtypes[cs[type]->type]);

   char *buffer   = malloc(20);
   sprintf(buffer, "%s_%s",main_cs[(int) type/4],color_subtypes[cs[type]->type]);
   return buffer; 
}

void print_all_cs() {
    for (size_t i = 0; i < NUM_CS_MAX; i++) { 
        fetch_color_builtin(i,true);
    }
}
