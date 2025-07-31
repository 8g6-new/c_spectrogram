#include "../../../../headers/libheatmap/colorschemes/opencv_like_scheme.h"

DEF_CS(autumn, 1024, 0);
DEF_CS(bone, 1024, 1);
DEF_CS(jet, 1024, 2);
DEF_CS(winter, 1024, 3);
DEF_CS(rainbow, 1024, 4);
DEF_CS(ocean, 1024, 5);
DEF_CS(summer, 1024, 6);
DEF_CS(spring, 1024, 7);
DEF_CS(cool, 1024, 8);
DEF_CS(hsv, 1024, 9);
DEF_CS(pink, 1024, 10);
DEF_CS(hot, 1024, 11);
DEF_CS(parula, 1024, 12);
DEF_CS(magma, 1024, 13);
DEF_CS(inferno, 1024, 14);
DEF_CS(plasma, 1024, 15);
DEF_CS(viridis, 1024, 16);
DEF_CS(cividis, 1024, 17);
DEF_CS(twilight, 1024, 18);
DEF_CS(twilight_shifted, 1024, 19);
DEF_CS(turbo, 1024, 20);
DEF_CS(deep_green, 1024, 21);


const color_data_t *cs[] = {
    &autumn, &bone, &jet, &winter, &rainbow,
    &ocean, &summer, &spring, &cool, &hsv,
    &pink, &hot, &parula, &magma, &inferno,
    &plasma, &viridis, &cividis, &twilight,
    &twilight_shifted, &turbo, &deep_green
};

const char *color_names[] = {
    "Autumn", "Bone", "Jet", "Winter", "Rainbow",
    "Ocean", "Summer", "Spring", "Cool", "HSV",
    "Pink", "Hot", "Parula", "Magma", "Inferno",
    "Plasma", "Viridis", "Cividis", "Twilight",
    "Twilight Shifted", "Turbo", "Deep Green"
};

char *fetch_color_opencv_like(cs_enum type, bool log) {
    const char *name = color_names[type];
    if (log)
        printf("enum %d  => Open CV Like color scheme : %s\n", type, name);

    char *copy = malloc(strlen(name) + 1);
    if (copy)
        strcpy(copy, name);
    return copy;
}

void print_all_cs() {
    for (size_t i = 0; i < NUM_CS_MAX; i++) {
        printf("%s\n",color_names[i]);
    }
}