const urls = [
    'https://docs.opencv.org/4.x/colorscale_autumn.jpg',
    'https://docs.opencv.org/4.x/colorscale_bone.jpg',
    'https://docs.opencv.org/4.x/colorscale_jet.jpg',
    'https://docs.opencv.org/4.x/colorscale_winter.jpg',
    'https://docs.opencv.org/4.x/colorscale_rainbow.jpg',
    'https://docs.opencv.org/4.x/colorscale_ocean.jpg',
    'https://docs.opencv.org/4.x/colorscale_summer.jpg',
    'https://docs.opencv.org/4.x/colorscale_spring.jpg',
    'https://docs.opencv.org/4.x/colorscale_cool.jpg',
    'https://docs.opencv.org/4.x/colorscale_hsv.jpg',
    'https://docs.opencv.org/4.x/colorscale_pink.jpg',
    'https://docs.opencv.org/4.x/colorscale_hot.jpg',
    'https://docs.opencv.org/4.x/colorscale_parula.jpg',
    'https://docs.opencv.org/4.x/colorscale_magma.jpg',
    'https://docs.opencv.org/4.x/colorscale_inferno.jpg',
    'https://docs.opencv.org/4.x/colorscale_plasma.jpg',
    'https://docs.opencv.org/4.x/colorscale_viridis.jpg',
    'https://docs.opencv.org/4.x/colorscale_cividis.jpg',
    'https://docs.opencv.org/4.x/colorscale_twilight.jpg',
    'https://docs.opencv.org/4.x/colorscale_twilight_shifted.jpg',
    'https://docs.opencv.org/4.x/colorscale_turbo.jpg',
    'https://docs.opencv.org/4.x/colorscale_deepgreen.jpg',
  ];
  
  const fs = await import('fs/promises');


  
  const downloadImages = async () => {
    const promises = urls.map(async (url) => {
      try {
        const response = await fetch(url);
        if (!response.ok) {
          throw new Error(`Failed to download ${url}. Status code: ${response.status}`);
        }
        const filename = url.split('/').pop().replace('colorscale_', '');
        const buffer = await response.arrayBuffer();
        await fs.writeFile(`gards/${filename}`, Buffer.from(buffer));
        console.log(`Downloaded ${filename} successfully!`);
      } catch (error) {
        console.error(error.message);
      }
    });
    await Promise.all(promises);
  };
  
  downloadImages();