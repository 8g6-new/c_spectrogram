#include "../../headers/audio_tools/audio_visualizer.h"

inline float  brachless_db(float mag,bool db){
    return !db*mag + db*log10f(1 + mag);
}

inline void spectrogram(stft_d *result, const char *output_file,float min, float max,unsigned char bg_clr[4],bool db,int cs_enum) {
    size_t w      = result->output_size;
    size_t h      = (size_t)(result->num_frequencies * max) / (result->sample_rate / 2);
    size_t min_f  = (size_t)(result->num_frequencies * min) / (result->sample_rate / 2);

    heatmap_t *hm = NULL;
    unsigned char *image = heatmap_get_vars(w, (h - min_f), &hm);

    if (!image || !hm) {
        fprintf(stderr, "Failed to allocate heatmap.\n");
        return;
    }
    #pragma omp parallel for
    for (size_t t = 0; t < w; t++) {
        for (size_t f = min_f; f < h; f++) {
            size_t inverted_axis = (h - f - 1); 
            float mag = result->magnitudes[t * result->num_frequencies + inverted_axis];
            mag =  brachless_db(mag,db);
            heatmap_add_weighted_point(hm, t, f , mag); 
        }
    }

    save_heatmap(image,&hm,output_file,w,h-min_f,bg_clr,cs_enum);

}


inline float *mel_spectrogram(stft_d *result,const char *output_file,size_t num_filters,float *mel_filter_bank,unsigned char bg_clr[4],bool db,int cs_enum){
    
    size_t w             = result->output_size;
    size_t h             = result->num_frequencies;

    
    float *mel_values    = (float*) malloc(num_filters * w * sizeof(float));

    heatmap_t *hm        = NULL;
    
    unsigned char *image = heatmap_get_vars(w, num_filters, &hm);

    if (!image || !hm) {
        fprintf(stderr, "Failed to allocate heatmap.\n");
        return NULL;
    }

    #pragma omp parallel for
    for (size_t t = 0; t < w; t++) {
        for (size_t mel = 0; mel < num_filters; mel++){
            float sum = 0; 
            for (size_t f = 0; f < h; f++) {
                size_t inverted_axis = (h-f-1); 

                float mag   = result->magnitudes[t*h+inverted_axis];

                sum+=mel_filter_bank[(num_filters - mel - 1) * h+inverted_axis] * mag*mag;
            }

            sum =  brachless_db(sum,db);
            
            mel_values[t*num_filters+(num_filters - mel - 1)] =  sum;

            heatmap_add_weighted_point(hm, t, mel,sum); 
        }
    }

    save_heatmap(image,&hm,output_file,w,num_filters,bg_clr,cs_enum);

    return mel_values;
}

inline void mfcc(float *mel_values,const char *output_file,size_t w,size_t num_filters,size_t num_coff,unsigned char bg_clr[4],int cs_enum){

    heatmap_t *hm = NULL;
    unsigned char *image = heatmap_get_vars(w, num_coff, &hm);

    #pragma omp parallel for
    for (size_t t = 0; t < w; t++) {
        for(size_t n=0;n<num_coff;n++){
            float sum = 0;
            for (size_t mel = 0; mel < num_filters; mel++){   
                sum+=(cos((M_PI/num_filters)*(mel+0.5)*n))  * mel_values[t*num_filters + mel];
            }
            heatmap_add_weighted_point(hm, t,num_coff - n-1, sum*sum);
        }
    }

    save_heatmap(image,&hm,output_file,w,num_coff,bg_clr,cs_enum);

}

