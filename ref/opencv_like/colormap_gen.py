from cv2 import imread,COLOR_BGR2RGB,cvtColor
from os import listdir,mkdir

path = "./src/libheatmap/colormaps/opencv_like/gards"

fols = listdir(path)

# out = "/kaggle/working/out"
# mkdir(out)

for name in fols:
    p = ""
    k = imread(f'{path}/{name}')
    k = cvtColor(k, COLOR_BGR2RGB)
    
    name = name.split(".")[0].replace(" ", "_")
    
    p = f"#ifndef {name.upper()}\n#define VIRIDIS\n#include <stdint.h>\n"
    p += f"const unsigned char {name}_data[{len(k[0])} * 4] =" + "{\n"
    
    
    # plt.imshow(k_rgb )
    c=0
    for i in k[0]:
        p += f"         {i[0]},{i[1]},{i[2]},{int(255 /(4-c%4))},\n"
        c+=1
    
    p+="};\n#endif"
    
    with open(f"src/libheatmap/colormaps/opencv_like/{name}.c", "w") as f:
        f.write(p)
