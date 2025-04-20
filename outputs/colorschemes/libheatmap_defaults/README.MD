# Builtin Color Schemes 

Libhaetmap offers a range of built-in color schemes to accommodate various data visualization needs. These color schemes are designed to provide clarity, contrast, and smooth transitions, depending on the nature of the data being presented. The schemes are divided into four subtypes, each tailored to different use cases:

| Subtype            | Description                                                                                                  |
|--------------------|--------------------------------------------------------------------------------------------------------------|
| **Discrete**        | Only <15 RGBA points. Best for categorical data with a small range. Lacks smooth transitions. Low dynamic range. |
| **Soft**            | Uses a smooth gradient with high dynamic range, ideal for continuous data. Emphasizes subtle variations and smooth transitions. Perfect for highlighting gradual changes in data. |
| **Mixed**           | Combines discrete colors and soft gradients, providing a balance of clarity and smoothness. Ideal for datasets with both categorical and continuous elements. Higher dynamic range than discrete. |
| **Mixed Exponential** | Similar to the mixed scheme but with an exponential gradient that emphasizes higher values or spikes. Ideal for datasets with outliers or extreme values, providing stronger contrast in those areas. Even higher dynamic range. |



## Color Scheme Descriptions

| **Color Scheme** | **Description** |
|------------------|-----------------|
| **Blues**        | A range of blue shades, typically progressing from light to dark, creating a calming and cool gradient. It often starts with light sky blue and moves towards deeper, more navy-like tones. |
| **Greys**        | A monochromatic scale with varying shades of grey, typically ranging from light silver to charcoal. This scheme focuses on neutrality and subtlety without strong color dominance. |
| **PuBuGn**       | A blend of purple, blue, and green hues, often transitioning from lavender to aqua or mint. It typically features soft and cool colors, emphasizing a smooth gradient across the three main colors. |
| **RdPu**         | A color scheme featuring shades of red and purple, transitioning from pale pinkish hues to deep, dark reds and purples. It combines warmth with richness. |
| **YlGn**         | A gradient moving from yellow to green, often progressing from light, bright yellow to deeper, more olive-like greens. It is a warm-to-cool gradient dominated by earthy greens and yellows. |
| **BrBG**         | A color scale with brown and blue hues, often transitioning from shades of brown and beige to dark blues and greens. It is typically used for diverging datasets, where both extremes of the color scale represent important values. |
| **Oranges**      | A gradient of orange hues, starting from light peach to deep burnt orange. The focus is on warmth and vibrancy, using a variety of orange shades to create a dynamic range. |
| **PuOr**         | A mix of purple and orange hues, transitioning from a deep purple to a bright orange. This scheme is used to represent diverging data, with purple and orange often symbolizing two opposing extremes. |
| **RdYlBu**       | A diverging color scheme that moves from red to yellow to blue, often starting with intense red shades, shifting to bright yellow, and ending with cool blue tones. This scheme is used for data that has a central, neutral value. |
| **YlOrBr**       | A gradient transitioning from yellow to orange and brown, with bright yellow moving to a deeper orange and eventually into earthy brown tones. It's used for a more natural or earthy representation. |
| **BuGn**         | A combination of blue and green, often progressing from light, soft blue to deeper, muted greens. It emphasizes a cool, serene feeling with harmonious color transitions. |
| **OrRd**         | An orange-to-red gradient that starts with soft orange tones and moves to more intense red hues, giving a warm and vibrant appearance. Often used for data that emphasizes intensity or heat. |
| **PuRd**         | A range of purple to red hues, typically starting with light purples and moving towards darker, more intense reds. The gradient evokes a warm, rich atmosphere. |
| **RdYlGn**       | A diverging scale that moves from red to yellow to green. This is often used in heatmaps to show extreme values in red, central values in yellow, and lower values in green, providing a clear visual representation of gradients in data. |
| **YlOrRd**       | A gradient moving from yellow to orange to red, with the scheme starting from light yellow and progressing through increasingly vibrant and deep orange and red tones. It suggests progression, intensity, and heat. |
| **BuPu**         | A combination of blue and purple hues, often starting with lighter blues and transitioning into deeper purples. It is used for cool, yet rich data visualizations. |
| **PiYG**         | A palette that combines pink, purple, and green hues, often creating a visually striking contrast. It transitions from vibrant pinks to deep purples and into fresh greens, giving a bold color representation. |
| **Purples**      | A gradient of purple tones, starting with light lavender and moving into darker, more royal purples. The scheme is primarily cool and offers a subtle and elegant progression. |
| **Reds**         | A range of red tones, from light pinkish hues to dark, intense reds. This scheme evokes warmth, energy, and intensity. |
| **GnBu**         | A gradient from green to blue, moving from soft, light greens to deeper, more intense blues. It is commonly used for visualizing natural or cool data patterns. |
| **PRGn**         | A palette that includes purple and green hues, with a range from deep purples to fresh greens. It is used for diverging data sets where the two extremes represent contrasting values. |
| **RdBu**         | A diverging color scheme moving from red to blue. It is widely used for datasets where extremes need to be represented with contrasting colors, making it suitable for visualizing opposites. |
| **Spectral**     | A range of colors that typically moves from red to orange, yellow, green, and blue, encompassing the full spectrum of colors. This scheme is used for diverse and colorful data visualizations, often in maps or datasets with multiple ranges. |
| **Greens**       | A palette of green hues, ranging from soft, light greens to deeper, more intense greens. This scheme is ideal for displaying data in natural, ecological, or growth contexts. |
| **PuBu**         | A gradient of purple to blue hues, moving from

## 1. Discrete Color Schemes

These schemes are ideal for representing categorical data with distinct color ranges.

| Enum ID | Name             | Spectrogram                                                                 |
|---------|------------------|-----------------------------------------------------------------------------|
| 0       | PRGn_discrete    | ![PRGn](discrete/black_woodpecker_stft_PRGn_discrete.png)                   |
| 4       | Blues_discrete   | ![Blues](discrete/black_woodpecker_stft_Blues_discrete.png)                 |
| 8       | BrBG_discrete    | ![BrBG](discrete/black_woodpecker_stft_BrBG_discrete.png)                   |
| 12      | BuGn_discrete    | ![BuGn](discrete/black_woodpecker_stft_BuGn_discrete.png)                   |
| 16      | BuPu_discrete    | ![BuPu](discrete/black_woodpecker_stft_BuPu_discrete.png)                   |
| 20      | GnBu_discrete    | ![GnBu](discrete/black_woodpecker_stft_GnBu_discrete.png)                   |
| 24      | Greens_discrete  | ![Greens](discrete/black_woodpecker_stft_Greens_discrete.png)               |
| 28      | Greys_discrete   | ![Greys](discrete/black_woodpecker_stft_Greys_discrete.png)                 |
| 32      | Oranges_discrete | ![Oranges](discrete/black_woodpecker_stft_Oranges_discrete.png)             |
| 36      | OrRd_discrete    | ![OrRd](discrete/black_woodpecker_stft_OrRd_discrete.png)                   |
| 40      | PiYG_discrete    | ![PiYG](discrete/black_woodpecker_stft_PiYG_discrete.png)                   |
| 44      | PuBu_discrete    | ![PuBu](discrete/black_woodpecker_stft_PuBu_discrete.png)                   |
| 48      | PuBuGn_discrete  | ![PuBuGn](discrete/black_woodpecker_stft_PuBuGn_discrete.png)               |
| 52      | PuOr_discrete    | ![PuOr](discrete/black_woodpecker_stft_PuOr_discrete.png)                   |
| 56      | PuRd_discrete    | ![PuRd](discrete/black_woodpecker_stft_PuRd_discrete.png)                   |
| 60      | Purples_discrete | ![Purples](discrete/black_woodpecker_stft_Purples_discrete.png)             |
| 64      | RdBu_discrete    | ![RdBu](discrete/black_woodpecker_stft_RdBu_discrete.png)                   |
| 68      | RdGy_discrete    | ![RdGy](discrete/black_woodpecker_stft_RdGy_discrete.png)                   |
| 72      | RdPu_discrete    | ![RdPu](discrete/black_woodpecker_stft_RdPu_discrete.png)                   |
| 76      | RdYlBu_discrete  | ![RdYlBu](discrete/black_woodpecker_stft_RdYlBu_discrete.png)               |
| 80      | RdYlGn_discrete  | ![RdYlGn](discrete/black_woodpecker_stft_RdYlGn_discrete.png)               |
| 84      | Reds_discrete    | ![Reds](discrete/black_woodpecker_stft_Reds_discrete.png)                   |
| 88      | Spectral_discrete| ![Spectral](discrete/black_woodpecker_stft_Spectral_discrete.png)           |
| 92      | YlGn_discrete    | ![YlGn](discrete/black_woodpecker_stft_YlGn_discrete.png)                   |
| 96      | YlGnBu_discrete  | ![YlGnBu](discrete/black_woodpecker_stft_YlGnBu_discrete.png)               |
| 100     | YlOrBr_discrete  | ![YlOrBr](discrete/black_woodpecker_stft_YlOrBr_discrete.png)               |
| 104     | YlOrRd_discrete  | ![YlOrRd](discrete/black_woodpecker_stft_YlOrRd_discrete.png)               |



## 2. Soft Color Schemes

These color schemes are more pastel-like and are best suited for smoother, subtler contrast in the spectrogram.

| Enum ID | Name            | Spectrogram                                                              |
|---------|------------------|--------------------------------------------------------------------------|
| 1       | PRGn_soft        | ![PRGn](soft/black_woodpecker_stft_PRGn_soft.png)                        |
| 5       | Blues_soft       | ![Blues](soft/black_woodpecker_stft_Blues_soft.png)                      |
| 9       | BrBG_soft        | ![BrBG](soft/black_woodpecker_stft_BrBG_soft.png)                        |
| 13      | BuGn_soft        | ![BuGn](soft/black_woodpecker_stft_BuGn_soft.png)                        |
| 17      | BuPu_soft        | ![BuPu](soft/black_woodpecker_stft_BuPu_soft.png)                        |
| 21      | GnBu_soft        | ![GnBu](soft/black_woodpecker_stft_GnBu_soft.png)                        |
| 25      | Greens_soft      | ![Greens](soft/black_woodpecker_stft_Greens_soft.png)                    |
| 29      | Greys_soft       | ![Greys](soft/black_woodpecker_stft_Greys_soft.png)                      |
| 33      | Oranges_soft     | ![Oranges](soft/black_woodpecker_stft_Oranges_soft.png)                  |
| 37      | OrRd_soft        | ![OrRd](soft/black_woodpecker_stft_OrRd_soft.png)                        |
| 41      | PiYG_soft        | ![PiYG](soft/black_woodpecker_stft_PiYG_soft.png)                        |
| 45      | PuBu_soft        | ![PuBu](soft/black_woodpecker_stft_PuBu_soft.png)                        |
| 49      | PuBuGn_soft      | ![PuBuGn](soft/black_woodpecker_stft_PuBuGn_soft.png)                    |
| 53      | PuOr_soft        | ![PuOr](soft/black_woodpecker_stft_PuOr_soft.png)                        |
| 57      | PuRd_soft        | ![PuRd](soft/black_woodpecker_stft_PuRd_soft.png)                        |
| 61      | Purples_soft     | ![Purples](soft/black_woodpecker_stft_Purples_soft.png)                  |
| 65      | RdBu_soft        | ![RdBu](soft/black_woodpecker_stft_RdBu_soft.png)                        |
| 69      | RdGy_soft        | ![RdGy](soft/black_woodpecker_stft_RdGy_soft.png)                        |
| 73      | RdPu_soft        | ![RdPu](soft/black_woodpecker_stft_RdPu_soft.png)                        |
| 77      | RdYlBu_soft      | ![RdYlBu](soft/black_woodpecker_stft_RdYlBu_soft.png)                    |
| 81      | RdYlGn_soft      | ![RdYlGn](soft/black_woodpecker_stft_RdYlGn_soft.png)                    |
| 85      | Reds_soft        | ![Reds](soft/black_woodpecker_stft_Reds_soft.png)                        |
| 89      | Spectral_soft    | ![Spectral](soft/black_woodpecker_stft_Spectral_soft.png)                |
| 93      | YlGn_soft        | ![YlGn](soft/black_woodpecker_stft_YlGn_soft.png)                        |
| 97      | YlGnBu_soft      | ![YlGnBu](soft/black_woodpecker_stft_YlGnBu_soft.png)                    |
| 101     | YlOrBr_soft      | ![YlOrBr](soft/black_woodpecker_stft_YlOrBr_soft.png)                    |
| 105     | YlOrRd_soft      | ![YlOrRd](soft/black_woodpecker_stft_YlOrRd_soft.png)                    |


## 3. Mixed Color Schemes

These schemes blend discrete and smooth colors, giving a balance between sharp segmentation and gradient.

| Enum ID | Name             | Spectrogram                                                               |
|---------|------------------|---------------------------------------------------------------------------|
| 2       | PRGn_mixed       | ![PRGn](mixed/black_woodpecker_stft_PRGn_mixed.png)                        |
| 6       | Blues_mixed      | ![Blues](mixed/black_woodpecker_stft_Blues_mixed.png)                      |
| 10      | BrBG_mixed       | ![BrBG](mixed/black_woodpecker_stft_BrBG_mixed.png)                        |
| 14      | BuGn_mixed       | ![BuGn](mixed/black_woodpecker_stft_BuGn_mixed.png)                        |
| 18      | BuPu_mixed       | ![BuPu](mixed/black_woodpecker_stft_BuPu_mixed.png)                        |
| 22      | GnBu_mixed       | ![GnBu](mixed/black_woodpecker_stft_GnBu_mixed.png)                        |
| 26      | Greens_mixed     | ![Greens](mixed/black_woodpecker_stft_Greens_mixed.png)                    |
| 30      | Greys_mixed      | ![Greys](mixed/black_woodpecker_stft_Greys_mixed.png)                      |
| 34      | Oranges_mixed    | ![Oranges](mixed/black_woodpecker_stft_Oranges_mixed.png)                  |
| 38      | OrRd_mixed       | ![OrRd](mixed/black_woodpecker_stft_OrRd_mixed.png)                        |
| 42      | PiYG_mixed       | ![PiYG](mixed/black_woodpecker_stft_PiYG_mixed.png)                        |
| 46      | PuBu_mixed       | ![PuBu](mixed/black_woodpecker_stft_PuBu_mixed.png)                        |
| 50      | PuBuGn_mixed     | ![PuBuGn](mixed/black_woodpecker_stft_PuBuGn_mixed.png)                    |
| 54      | PuOr_mixed       | ![PuOr](mixed/black_woodpecker_stft_PuOr_mixed.png)                        |
| 58      | PuRd_mixed       | ![PuRd](mixed/black_woodpecker_stft_PuRd_mixed.png)                        |
| 62      | Purples_mixed    | ![Purples](mixed/black_woodpecker_stft_Purples_mixed.png)                  |
| 66      | RdBu_mixed       | ![RdBu](mixed/black_woodpecker_stft_RdBu_mixed.png)                        |
| 70      | RdGy_mixed       | ![RdGy](mixed/black_woodpecker_stft_RdGy_mixed.png)                        |
| 74      | RdPu_mixed       | ![RdPu](mixed/black_woodpecker_stft_RdPu_mixed.png)                        |
| 78      | RdYlBu_mixed     | ![RdYlBu](mixed/black_woodpecker_stft_RdYlBu_mixed.png)                    |
| 82      | RdYlGn_mixed     | ![RdYlGn](mixed/black_woodpecker_stft_RdYlGn_mixed.png)                    |
| 86      | Reds_mixed       | ![Reds](mixed/black_woodpecker_stft_Reds_mixed.png)                        |
| 90      | Spectral_mixed   | ![Spectral](mixed/black_woodpecker_stft_Spectral_mixed.png)                |
| 94      | YlGn_mixed       | ![YlGn](mixed/black_woodpecker_stft_YlGn_mixed.png)                        |
| 98      | YlGnBu_mixed     | ![YlGnBu](mixed/black_woodpecker_stft_YlGnBu_mixed.png)                    |
| 102     | YlOrBr_mixed     | ![YlOrBr](mixed/black_woodpecker_stft_YlOrBr_mixed.png)                    |
| 106     | YlOrRd_mixed     | ![YlOrRd](mixed/black_woodpecker_stft_YlOrRd_mixed.png)                    |


## 4. Mixed-Exponential Color Schemes

These provide sharper gradients and contrast, emphasizing transitions for better feature visibility.

| Enum ID | Name                | Spectrogram                                                                   |
|---------|---------------------|-------------------------------------------------------------------------------|
| 3       | PRGn_mixed_exp      | ![PRGn](mixed_exp/black_woodpecker_stft_PRGn_mixed_exp.png)                   |
| 7       | Blues_mixed_exp     | ![Blues](mixed_exp/black_woodpecker_stft_Blues_mixed_exp.png)                 |
| 11      | BrBG_mixed_exp      | ![BrBG](mixed_exp/black_woodpecker_stft_BrBG_mixed_exp.png)                   |
| 15      | BuGn_mixed_exp      | ![BuGn](mixed_exp/black_woodpecker_stft_BuGn_mixed_exp.png)                   |
| 19      | BuPu_mixed_exp      | ![BuPu](mixed_exp/black_woodpecker_stft_BuPu_mixed_exp.png)                   |
| 23      | GnBu_mixed_exp      | ![GnBu](mixed_exp/black_woodpecker_stft_GnBu_mixed_exp.png)                   |
| 27      | Greens_mixed_exp    | ![Greens](mixed_exp/black_woodpecker_stft_Greens_mixed_exp.png)               |
| 31      | Greys_mixed_exp     | ![Greys](mixed_exp/black_woodpecker_stft_Greys_mixed_exp.png)                 |
| 35      | Oranges_mixed_exp   | ![Oranges](mixed_exp/black_woodpecker_stft_Oranges_mixed_exp.png)             |
| 39      | OrRd_mixed_exp      | ![OrRd](mixed_exp/black_woodpecker_stft_OrRd_mixed_exp.png)                   |
| 43      | PiYG_mixed_exp      | ![PiYG](mixed_exp/black_woodpecker_stft_PiYG_mixed_exp.png)                   |
| 47      | PuBu_mixed_exp      | ![PuBu](mixed_exp/black_woodpecker_stft_PuBu_mixed_exp.png)                   |
| 51      | PuBuGn_mixed_exp    | ![PuBuGn](mixed_exp/black_woodpecker_stft_PuBuGn_mixed_exp.png)               |
| 55      | PuOr_mixed_exp      | ![PuOr](mixed_exp/black_woodpecker_stft_PuOr_mixed_exp.png)                   |
| 59      | PuRd_mixed_exp      | ![PuRd](mixed_exp/black_woodpecker_stft_PuRd_mixed_exp.png)                   |
| 63      | Purples_mixed_exp   | ![Purples](mixed_exp/black_woodpecker_stft_Purples_mixed_exp.png)             |
| 67      | RdBu_mixed_exp      | ![RdBu](mixed_exp/black_woodpecker_stft_RdBu_mixed_exp.png)                   |
| 71      | RdGy_mixed_exp      | ![RdGy](mixed_exp/black_woodpecker_stft_RdGy_mixed_exp.png)                   |
| 75      | RdPu_mixed_exp      | ![RdPu](mixed_exp/black_woodpecker_stft_RdPu_mixed_exp.png)                   |
| 79      | RdYlBu_mixed_exp    | ![RdYlBu](mixed_exp/black_woodpecker_stft_RdYlBu_mixed_exp.png)               |
| 83      | RdYlGn_mixed_exp    | ![RdYlGn](mixed_exp/black_woodpecker_stft_RdYlGn_mixed_exp.png)               |
| 87      | Reds_mixed_exp      | ![Reds](mixed_exp/black_woodpecker_stft_Reds_mixed_exp.png)                   |
| 91      | Spectral_mixed_exp  | ![Spectral](mixed_exp/black_woodpecker_stft_Spectral_mixed_exp.png)           |
| 95      | YlGn_mixed_exp      | ![YlGn](mixed_exp/black_woodpecker_stft_YlGn_mixed_exp.png)                   |
| 99      | YlGnBu_mixed_exp    | ![YlGnBu](mixed_exp/black_woodpecker_stft_YlGnBu_mixed_exp.png)               |
| 103     | YlOrBr_mixed_exp    | ![YlOrBr](mixed_exp/black_woodpecker_stft_YlOrBr_mixed_exp.png)               |
| 107     | YlOrRd_mixed_exp    | ![YlOrRd](mixed_exp/black_woodpecker_stft_YlOrRd_mixed_exp.png)               |
