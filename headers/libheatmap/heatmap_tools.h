#ifndef HEATMAP_TOOLS_H
#define HEATMAP_TOOLS_H
   
    #include <stdlib.h>
    #include <string.h>

    #include "heatmap.h"
    #include "../png_tools/png_tools.h"
    

    #ifdef BUILTIN
        #include "colorschemes/builtin_scheme.h"
        #define cs_from_enum fetch_color_builtin
    #endif
    #ifdef OPENCV_LIKE
        #include "colorschemes/opencv_like_scheme.h"
        #define cs_from_enum fetch_color_opencv_like
    #endif
    
    unsigned char* heatmap_get_vars(size_t w, size_t h,heatmap_t **hm_out); 
    int save_heatmap(unsigned char *image,heatmap_t **hm,char *output_file,size_t w,size_t h,unsigned char bg_clr[4],int cs_enum);

#endif