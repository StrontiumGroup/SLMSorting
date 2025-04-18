# Fast SLM sorting repository
This is a ReadME-file for the Fast SLM repository, such that one is not immediately lost when looking at it. The work is the basis of the code used for the paper `Parallel assembly of neutral atom arrays with an SLM using linear phase interpolation`. The code has been developped for research purposes and as such, is not terribly well documented. For questions, please contact i.h.a.knottnerus[at]uva.nl.

## Background
The goal of the repository is to use an Ultrafast SLM of Meadowlark to dynamically control the position of single atoms in tweezers. To do so, we employ OpenCL for GPU-accellerated computation of phase masks and have tweaked the known IFTA algorithms such that we lose the iterative aspect and gain further in speed. 

Current benchmarking gives 2.1ms as the average time for updating the spot positions, updating the GPU buffers, calculating the next phase mask, and displaying this on the fast SLM. We expect to increase further by:
* Parallelizing the calculation of the next pattern with the display on the SLM. Right now, the calculation time is about the same time as the time it takes to display on to the SLM (~1ms each). Both processes are independent of each other and we could already calculate the next pattern while waiting for the SLM to finish displaying. 

## To-do
Upcoming changes:
* Exploring multidimensional kernels
* Organizing the repository
* Weed out kernels in `kernel.cl`
* Improve documentation in code
* ...

## Repository structure
The structure of the repository is as follows:
* `bin`: The folder containing the executables. Make sure that the Meadowlark libraries are in this folder, because one cannot talk to the SLM otherwise.
* `core`: In this folder, the relevant header files and c++ files are located. In principle, one does not need to change anything here. You do link it later in the building process.
* `lib`: The relevant libraries for building. Not sure if this is exactly needed, since half of these libraries need to be in the `bin`-folder as well because of how the Blink SDK is programmed.
* `examples`: Example c++ files to use the project for doing specific jobs. There is mostly poor documentation on each of these.
* `tests`: Test files for benchmarking the performance.

## Compiling the project
Let's first see how to get the code working on a device. This is very dependent on the device that you are using. The code has been tested on Windows pcs with MingW64. Assuming you have a file called `FILE.cpp` located in a subfolder of the main folder, one can build an executable using the command:
```
gcc FILE.cpp ../core/*.cpp -DVKFFT_BACKEND=3 -lstdc++ -lOpenCl -lBlink_C_wrapper -I ../core/ -o ../bin/OUTPUT.exe -lmosquitto
```

The Meadowlark SLM libraries are written such that the .dll files should be in the same folder as the executable. For the GPU code, one needs to pre-install VkFFT and OpenCL.

For a walkthrough, please check out my personal path installing the files on a new device [here](./docs/GettingStarted.md). 

## Examples
We provide small examples in `examples/` of how to use the code. In our case, we had to communicate with the code during our experiment using MQTT. Again, this is very experiment specific and we have not strived to make an universally accessible example. The basics on how to call the classes in the `core`-folder are laid out, such that it should be possible to adapt to ones own experiment.

- `example_noslm.cpp`: An example with no SLM commands. This is useful if no SLM is connected or the SDK is not installed on the PC. It sorts a 6x6 array into a 3x3 grid and stores patterns in between in the folder `/masks`. **Note** that you still need the `-lBlink_C_wrapper` tag and the Meadowlark `.dll`-files to compile the program.
- `example_slm.cpp`: Does the same sorting but does not save the masks, but rather displays them on the SLM. 
- `example_multipattern_noslm.cpp`: Another example without SLM commands, but this time sorting the 6x6 pattern into a 16-atom circle. It saves the holograms in the `/masks` folder.
- `example_multipattern_slm.cpp`: Does not save but displays the holograms on the SLM.