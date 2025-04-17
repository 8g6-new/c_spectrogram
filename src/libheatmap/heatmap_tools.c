#include "../../headers/libheatmap/heatmap_tools.h"


unsigned char* heatmap_get_vars(size_t w, size_t h,heatmap_t **hm_out) {
    unsigned char *image = malloc(sizeof(unsigned char) * w * h * 4);

    if (!image) {
        perror("malloc heatmap");
        return NULL;
    }
    
    memset(image, 0, w * h * 4);

    heatmap_t *hm = heatmap_new(w,h);
    
    if (!hm) {
        fprintf(stderr, "Failed to create heatmap.\n");
        return NULL;
    }

    *hm_out = hm;

    return image;
}


int save_heatmap(unsigned char *image,heatmap_t **hm,char *output_file,size_t w,size_t h,unsigned char bg_clr[4],int cs_enum){

    if(cs_enum > NUM_CS_MAX-1)
     return -1;

    // printf("data @ 20 %d and total size : %d",cs[cs_enum]->data[20],cs[cs_enum]->size);
    heatmap_colorscheme_t *scheme =  heatmap_colorscheme_load(cs[cs_enum]->data,cs[cs_enum]->size);

    heatmap_render_to(*hm,scheme,&image[0]);

    add_bg(image, w,h, bg_clr);
    save_png(output_file,image,w,h);

    free(image);
    heatmap_free(*hm);
    heatmap_colorscheme_free(scheme);
}



