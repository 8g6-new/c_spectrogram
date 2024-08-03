let   {exec}       = require('child_process')
const {promisify}  = require('util')

let  {readdir,mkdir}   = require('fs/promises')

exec = promisify(exec)

async function try_mkdir(dir){
    try {
        await mkdir(dir, { recursive: true });
    } catch (error) {
        if (error.code!== 'EEXIST') {
            throw error;
        }
    }
}


async function read(inp,out) {

    await try_mkdir(out)


    const fols=await readdir(inp);

    for(const fol of fols){
          
        let files  = []
    
        try {
            files = await readdir(`${inp}/${fol}`)
        } catch (error) {
            console.error('Error reading files:', error)
        }
        
        let s = new Date().getTime()

        let l = files.length
        
        files = await Promise.all(files.map(async(file) =>{
            let {stdout} = await exec(`./stft_heatmap "${inp}/${fol}/${file}" 512 128 hann 30 "${out}/${fol}/${file}.png"`);
            stdout;
        }))


        let e = new Date().getTime()

        let time = (e-s)
       


        console.log(`${time}ms for ${l} 1 second wav files ie ${l*1000/time} us per file aka ${l*1000/(time)} files per second`)
    }

}

read(`./test`,`./stft_spectrogram`)