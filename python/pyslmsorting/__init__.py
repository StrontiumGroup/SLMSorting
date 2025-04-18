# -*- coding: utf-8 -*-
"""
pySLMSorting

Module created for SLM sorting code as presented in:

* https://arxiv.org/abs/2501.01391

Uses OpenCl and pyVKFFT for fast computation of holograms.
Note that the code only calculates holograms and does not
take care of e.g. rearranging tweezers. The user has to do
that separately. See for example, Example.ipynb.

author: Ivo Knottnerus
lastedit: 17/04/2025
"""
import numpy as np
from numpy import abs, sum, mean, amax
from PIL import Image

import os

import pyopencl as cl
import pyopencl.array as cla
from pyopencl.clrandom import rand
from pyopencl.elementwise import ElementwiseKernel
from pyvkfft.opencl import VkFFTApp

import matplotlib.pyplot as plt

def roll(x, dim):
    """
    Moves a spot coordinate to the rim of the coordinates.

    Acts as an FFTShift in coordinate space.
    """
    return x % dim

def unroll(x, dim):
    """
    The inverse of roll().
    """
    if x > dim // 2:
        return x - dim
    return x


class GPUHandler:
    """
    Class that handles all the GPU kernel executions.

    Houses the OpenCl variables and buffers. It also contains
    the kernels and the VkFFTApp to perform FFTs.
    """
    def __init__(self, parent):
        self.parent = parent

        self.app = None
        self.ctx = None
        self.cq = None

        self.knl_get_intensity = None
        self.knl_get_target = None
        self.knl_update_target = None
        self.knl_change_amps = None
        self.knl_assemble_fields = None
        self.knl_get_phase = None
        self.knl_set_field_from_spots = None
        self.knl_update_field_shortcut = None

        self.buffer_status = False
        self.field_buf = None
        self.target_buf = None
        self.phase_buf = None
        self.ui_buf = None
        self.x_buf = None
        self.y_buf = None
        self.target_intensity_buf = None
        self.intensity_buf = None
        self.initial_phase_buf = None

        self.max_iter = 1000
    
    def initOpenCl(self, do_output=True):
        """
        Initiates the OpenCl variables and VKFFTApp.

        Creates the context, selects a device (can vary per computer) and 
        command queue for OpenCl. Also loads the pyVkFFT VKFFTApp.
        """
        # Create some context on the first available GPU
        if 'PYOPENCL_CTX' in os.environ:
            ctx = cl.create_some_context()
        else:
            ctx = None
            # Find the first OpenCL GPU available and use it, unless
            for p in cl.get_platforms():
                for d in p.get_devices():
                    if d.type & cl.device_type.GPU == 1:
                        continue
                    if do_output:
                        print("Selected device: ", d.name)
                    ctx = cl.Context(devices=(d,))
                    #break
                if ctx is not None:
                    break

        self.ctx = ctx
        self.cq = cl.CommandQueue(ctx)
        self.app = VkFFTApp(
            (self.parent.fft_dimension, self.parent.fft_dimension),
            np.complex64,
            queue=self.cq,
            ndim=2,
            inplace=True,
            norm=1
        )

        self.buildKernels()

    def initBuffers(self):
        """
        Initializes the buffers for the calculation.
        """
        # Some useful memory on the GPU.
        self.field_buf = cla.empty(self.cq, (self.parent.fft_dimension, self.parent.fft_dimension), np.complex64)
        self.target_buf = cla.zeros(self.cq, (self.parent.fft_dimension, self.parent.fft_dimension), np.float32)
        self.phase_buf = cla.zeros(self.cq, (self.parent.fft_dimension, self.parent.fft_dimension), np.float32)
        self.ui_buf = cla.to_device(self.cq, self.parent.Ui.astype(np.float32))

        self.buffer_status = True

    def buildKernels(self):
        """
        Compiles the program with the kernels.

        For a full description of the kernels, check the documentation.
        """
        prg = cl.Program(
            self.ctx,
            """
            #pragma OPENCL EXTENSION cl_khr_fp64 : enable
            #define PYOPENCL_DEFINE_CDOUBLE
            #define SLM_WIDTH """+str(self.parent.width)+"""
            #define SLM_HEIGHT """+str(self.parent.height)+"""
            #include <pyopencl-complex.h>
            __kernel void get_intensity(
                __global float *intensity,
                __global int *x,
                __global int *y,
                __global cfloat_t *field
            )
            {
                int i = get_global_id(0);

                intensity[i] = cfloat_abs_squared(field[(int) y[i]*SLM_WIDTH + x[i]]);
            }

            __kernel void get_target(
                __global int *x,
                __global int *y,
                __global float *target_intensity,
                __global float *target
            ) {
                int i = get_global_id(0);
                int size = get_global_size(0);

                target[(int) (y[i]*SLM_WIDTH + x[i])] = target_intensity[i];
            }

            __kernel void set_field_from_spots(
                __global float *target_phase, 
                __global float *intensity, 
                __global int *x, 
                __global int *y,
                __global cfloat_t *field

            ){
                size_t i = get_global_id(0);
                size_t size = get_global_size(0);

                field[(int) (y[i]*SLM_WIDTH + x[i])] = cfloat_new( 
                    sqrt(intensity[i]) * cos(target_phase[i]), 
                    sqrt(intensity[i]) * sin(target_phase[i]) 
                );

            }

            __kernel void update_target(
                __global int *x,
                __global int *y,
                __global float *target_intensity,
                __global float *intensity,
                __global float *target
            ){
                int i = get_global_id(0);

                float mean_val = 0;
                int size = get_global_size(0);
                for (int j = 0; j < size; j++)
                    mean_val += sqrt((intensity[j] / target_intensity[j])) / size;

                if (intensity[i] != 0. && target_intensity[i] != 0.)
                    target[(int) (y[i]*SLM_WIDTH + x[i])] *= (mean_val / sqrt(intensity[i] / target_intensity[i]));
            }

            __kernel void update_field_shortcut(
                __global float *target_phase, 
                __global float *intensity, 
                __global int *x, 
                __global int *y,
                __global cfloat_t *field

            ){
                size_t i = get_global_id(0);

                // the factor 1000 is a dirty trick.
                field[(int) (y[i]*SLM_WIDTH + x[i])] = cfloat_new( 
                    10000. * sqrt(intensity[i]) * cos(target_phase[i]), 
                    10000. * sqrt(intensity[i]) * sin(target_phase[i]) 
                );

            }
            """
        ).build()

        self.knl_change_amps = ElementwiseKernel(
            self.ctx,
            "cfloat_t *field, float *intensity",
            """
            float value = cfloat_abs(field[i]);
            if (value != 0.) {
                field[i] = cfloat_mulr(field[i], sqrt(intensity[i]) / value);
            }
            else {
                field[i] = cfloat_addr(field[i], sqrt(intensity[i]));
            }

            """,
            "knl_change_amps",
            preamble="#define PYOPENCL_DEFINE_CDOUBLE //#include <pyopencl-complex.h>"
        )

        self.knl_assemble_field = ElementwiseKernel(
            self.ctx,
            "cfloat_t *field, float *intensity, float *phase",
            """
            field[i] = cfloat_new(sqrt(intensity[i]) * cos(phase[i]), sqrt(intensity[i]) * sin(phase[i]));
            """,
            "knl_assemble_field",
            preamble="#define PYOPENCL_DEFINE_CDOUBLE //#include <pyopencl-complex.h>"
        )

        self.knl_get_phase = ElementwiseKernel(
            self.ctx,
            "float *phase, cfloat_t *field",
            """
            phase[i] = atan2(cfloat_imag(field[i]), cfloat_real(field[i]));
            """,
            "knl_get_phase",
            preamble="#define PYOPENCL_DEFINE_CDOUBLE //#include <pyopencl-complex.h>"
        )

        # Compile the kernels for later use.
        self.knl_get_target = prg.get_target
        self.knl_update_target = prg.update_target
        self.knl_get_intensity = prg.get_intensity
        self.knl_set_field_from_spots = prg.set_field_from_spots
        self.knl_update_field_shortcut = prg.update_field_shortcut

    def copySpotsToBuffers(self, used_spots):
        """
        Copies the coordinates, intensity and phase of given spots to the GPU.
        """
        self.x_buf = cla.to_device(self.cq, used_spots[:, 0].astype(np.int32))
        self.y_buf = cla.to_device(self.cq, used_spots[:, 1].astype(np.int32))
        self.target_intensity_buf = cla.to_device(self.cq, used_spots[:, 3].astype(np.float32) / used_spots.shape[0])
        self.initial_phase_buf = cla.to_device(self.cq, used_spots[:, 4].astype(np.float32))
        self.intensity_buf = cla.zeros(self.cq, (used_spots.shape[0],), np.float32) + 1
    
    def releaseBuffers(self):
        """
        Clears the buffers. 
        """
        self.field_buf.data.release()
        self.target_buf.data.release()
        self.phase_buf.data.release()
        self.ui_buf.data.release()

        self.x_buf.data.release()
        self.y_buf.data.release()
        self.target_intensity_buf.data.release()
        self.intensity_buf.data.release()
        self.initial_phase_buf.data.release()

        self.buffer_status = False

    def calculateGSW(self, used_spots, start_hologram=None, max_delta=1e-3,
                    dev_delta=1e-3, show_plot=False, save_path=""):
        """
        Computes a hologram to produce the tweezer pattern given in used_spots.

        Uses the weighed Gerchberg-Saxton algorithm, described in e.g.:
        
        * https://doi.org/10.1364/OE.15.001913
        * https://doi.org/10.1364/oe.27.002184
        * https://doi.org/10.1364/ol.44.003178

        The function takes a potential starting hologram in `start_hologram`, which
        can be used as the starting point of the iteration. If None (default), the 
        hologram in the parent PMG instance will be used. If this is also None, a
        hologram will be made using the presented `used_spots`. This is the starting
        point of the GSW algorithm.

        In the loop of the GSW algorithm, an FFT is performed and at the coordinates 
        of the `used_spots` the intensity is compared to the target intensity. If 
        this matches, within conditions set by `max_delta`, the loop is terminated. 
        Alternatively, is the change in subsequent iterations is small compared to 
        `dev_delta`, the loop stops. Last case to terminate is reaching the 
        `max_iter`. 

        Function can plot and save the hologram directly.

        Parameters
        ----------
        used_spots : np.ndarray
            Array of the coordinates, intensities phases of the desired pattern.
        start_hologram : np.ndarray or NoneType
            Array with the start hologram for the iteration. Defaults to None.
        max_delta : float
            Maximal discrepancy in intensity allowed to stop loop. Defaults to 
            1e-3.
        dev_delta : float
            Minimial deviation between subsequent iterations to not stop loop.
            Defaults to 1e-3.
        show_plot : bool
            Whether to show a plot of the hologram and the expected intensity.
            Defaults to False.
        save_path : string
            If presented a non-zero length string, this will be used as the path
            to save the hologram as an 8-bit bitmap. Defaults to "".
        
        Returns
        -------
        hologram : np.ndarray
            Received hologram.
        """
        # We first determine the target intensity distribution.
        self.knl_get_target(
            self.cq,
            self.target_intensity_buf.shape,
            None,
            self.x_buf.data,
            self.y_buf.data,
            self.target_intensity_buf.data,
            self.target_buf.data,
        )

        # Check if a hologram was presented or still in parent memory.
        if np.any(start_hologram) or np.any(self.parent.hologram):
            if np.any(start_hologram):
                self.phase_buf = cla.to_device(self.cq, start_hologram)
            else:
                self.phase_buf = cla.to_device(self.cq, self.parent.hologram)

            # Calculate the field at SLM.
            self.cq.finish()
            self.knl_assemble_field(self.field_buf, self.ui_buf, self.phase_buf)
            self.cq.finish()
            self.app.fft(self.field_buf)
        else:
            # Else, reset field buffer and just place the spots.
            self.field_buf = cla.empty(self.cq, (self.parent.fft_dimension, self.parent.fft_dimension), np.complex64)         
            self.knl_set_field_from_spots(
                self.cq,
                self.initial_phase_buf.shape,
                None,
                self.initial_phase_buf.data,
                self.target_intensity_buf.data,
                self.x_buf.data,
                self.y_buf.data,
                self.field_buf.data
            )

        locked = False
        it = 1
        prev_delta = -1
        delta = 1
        while it <= self.max_iter:
            # Start by getting metrics on current performance.
            self.cq.finish()
            if it > 1:
                self.knl_get_intensity(
                    self.cq,
                    self.intensity_buf.shape,
                    None,
                    self.intensity_buf.data,
                    self.x_buf.data,
                    self.y_buf.data,
                    self.field_buf.data
                )
                self.cq.finish()
                
                # Then update the target using the weights.
                self.knl_update_target(
                    self.cq,
                    self.intensity_buf.shape,
                    None,
                    self.x_buf.data,
                    self.y_buf.data,
                    self.target_intensity_buf.data,
                    self.intensity_buf.data,
                    self.target_buf.data
                )
                self.cq.finish()

                # Calculate metrics. I don't know if we can speed this up further.
                # On home-pc, takes ~1ms for me.
                intensities = self.intensity_buf.get()
                mean_intensities = mean(intensities)
                delta = amax(
                    [used_spots[i, 3] - intensities[i] / mean_intensities for i in range(used_spots.shape[0])]
                )

                # Change the amplitudes of the field at focus with our constraints.
                self.knl_change_amps(self.field_buf, self.target_buf)
                self.cq.finish()


            # Calculate using GPU to field at SLM.
            self.app.ifft(self.field_buf)
            self.cq.finish()

            # See if we should stop because we had hit the metrics
            if (it > 1 and delta < max_delta) or (it > 1 and delta < 0.01 and abs(
                    prev_delta - delta) < dev_delta) or it == self.max_iter:
                break

            # Get field at SLM and calculate back to focus.
            self.knl_change_amps(self.field_buf, self.ui_buf)
            self.cq.finish()
            self.app.fft(self.field_buf)
            self.cq.finish()

            # Save performance for the next iteration and up counter.
            prev_delta = delta
            it += 1

        # Store phase in the CPU memory
        self.knl_get_phase(self.phase_buf, self.field_buf)
        hologram = np.fft.fftshift(self.phase_buf.get())
        return hologram

    def calculateSingle(self):
        """
        Computes a single iFFT of a target electric field.

        First sets the intended electric field that is given by the buffers
        in the GPU memory. Then performs a single iFFT and returns the phase.

        Returns
        -------
        hologram : np.ndarray
            The received hologram.
        """ 
        self.knl_update_field_shortcut(
            self.cq,
            self.initial_phase_buf.shape,
            None,
            self.initial_phase_buf.data,
            self.target_intensity_buf.data,
            self.x_buf.data,
            self.y_buf.data,
            self.field_buf.data
        )
        self.cq.finish()
        self.app.ifft(self.field_buf)
        self.cq.finish()
        self.knl_get_phase(self.phase_buf, self.field_buf)

        # Perhaps the cleanest would be to do another FFT here,
        # so that our target field is almost all zeros.

        hologram = np.fft.fftshift(self.phase_buf.get())
        return hologram

class PhaseMaskGenerator:
    """ Class used for creation of phase masks and iteration of them.

    More complete explanation follows later.

    """
    def __init__(self, width, height):
        """
        Initializes the PMG with dimensions width x height.
        """
        self.width = width
        self.height = height
        self.fft_dimension = max(width, height)

        # Initialize the arrays with the correct data types.
        self.spots = np.array([], dtype=np.float32)
        self.mask = np.empty([], dtype=bool)
        self.Ui = np.zeros((self.fft_dimension, self.fft_dimension))
        self.Ui[
            int((max(height, width) - height) / 2): int((max(height, width) + height) / 2),
            int((max(height, width) - width) / 2): int((max(height, width) + width) / 2)
        ] = 1. / (height * width)
        self.Ui = self.Ui.astype(np.float32)
        self.hologram = None

        self.gpu = GPUHandler(self)
        self.gpu.initOpenCl()

    def addSpot(self, x, y, z, i=1., psi=None):
        """
        Adds a spot to the np.array containing the spots.

        Parameters
        ----------
        x : np.float32
            x-coordinate, should be integer
        y : np.float32
            y-coordinate, should be integer
        z : np.float32
            z-coordinate, unused for now
        i : np.float32
            relative target intensity. Default: 1.
        psi : NoneType or np.float32
            optical phase of the tweezer
        """
        if not psi:
            psi = np.random.uniform(-np.pi, np.pi)
        if self.spots.size == 0:
            if i > 0:
                self.spots = np.array(
                    [[x % self.width, y % self.height, z, i, psi]], 
                    dtype=np.float32
                )
        else:
            double = False
            for spot in self.spots:
                if x == spot[0] and y == spot[1] and z == spot[2]:
                    double = True
                    break
            if not double:
                if i > 0:
                    self.spots = np.vstack([
                        self.spots, 
                        [x % self.width, y % self.height, z, i, psi]
                    ])

    def deleteSpot(self, i):
        """
        Removes the spot at index i.
        """
        self.spots = np.delete(
            self.spots,
            i,
            axis=0
        )

    def getSpots(self):
        """
        Returns the spots at the locations with the optical center in the middle.
        """
        return_spots = np.copy(self.spots)
        for i in range(self.spots.shape[0]):
            return_spots[i, 0] = unroll(return_spots[i, 0], self.width)
            return_spots[i, 1] = unroll(return_spots[i, 1], self.height)
        return return_spots
    
    def sortSpots(self):
        """
        Sorts the spots from low to high coordinates.

        This implementation is quite stupid in complexity, but N is 
        usually quite small (max ~1000).
        """
        return_spots = self.getSpots()
        self.spots = self.spots[return_spots[:, 1].argsort()]
        return_spots = self.getSpots()
        self.spots = self.spots[return_spots[:, 0].argsort()]

    def calculateGSW(self, start_hologram=None, max_delta=1e-3,
                    min_eff=-1, dev_delta=1e-3, dev_eff=1e-4,
                    show_plot=False, save_path=""):
        """
        Uses the weighted GS algorithm.

        This version uses fft to quickly compute 2D patterns. 
        See https://doi.org/10.1364/OL.44.003178 for the algorithm.

        Parameters
        ----------
        start_hologram : np.ndarray
            A hologram that is used as a starting point for the calculation.
            If None, first the hologram in the memory 
            (`PhaseMaskGenerator.hologram`) will be used and if this is also
            None, a random phase is assigned to each tweezer as a starting
            point.
        max_delta : float
            If deviation is below this number, the iteration will stop
        min_eff : float
            Minimally required efficiency in the sum of the tweezers
        dev_delta : float
            Minimal change in delta that is still considered as improvment
        dev_eff : float
            Minimal change in efficiency for improvement
        show_plot : bool
            Whether to display a plot of the hologram and the expected
            tweezer pattern
        save_path : string
            If presented with a string, the hologram is stored at the
            save_path
        """

        # No spots are defined, so we better stop the calculation.
        if self.spots.size == 0:
            return

        used_spots = np.copy(self.spots[self.mask])
        if not self.gpu.buffer_status:
            self.gpu.initBuffers()
        self.gpu.copySpotsToBuffers(used_spots)
        self.hologram = self.gpu.calculateGSW(
            used_spots,
            start_hologram, 
            max_delta,
            min_eff,
            dev_delta,
            dev_eff
        )
        self.gpu.releaseBuffers()

        # Optional function to see the results.
        if show_plot:
            self.display()

        # Optional function to save hologram.
        if len(save_path):
            self.save(save_path)

        return self.hologram

    def calculateSingle(self, show_plot=False, save_path=""):
        """
        Calculates only a single FFT from the currently used spots.

        This is useful for the linear interpolation rearrangement
        when all the coordinates and optical phases of the tweezers
        are known. The function copies this information to the GPU
        buffers and then runs a single FFT back and retrieves the 
        hologram.

        Parameters
        ----------
        show_plot : bool
            Whether a plot is displayed of the hologram and the expected
            pattern when using the illumination in the memory. Defaults
            to False.
        save_path : string
            Defaults to "". If a non-zero length string is presented, 
            the hologram will be saved to that path.
        
        Returns
        -------
        self.hologram : np.ndarray
            Array of floats containing the hologram that is calculated.
        """
        # No spots are defined, so we better stop the calculation.
        if self.spots.size == 0:
            return

        used_spots = np.copy(self.spots[self.mask])
        if not self.gpu.buffer_status:
            self.gpu.initBuffers()
        self.gpu.copySpotsToBuffers(used_spots)
        self.hologram = self.gpu.calculateSingle()

        # Optional function to see the results.
        if show_plot:
            self.display()

        # Optional function to save hologram.
        if len(save_path):
            self.save(save_path)

        return self.hologram


    def display(self):
        """
        Displays the current hologram and a projected pattern.

        Uses illumination pattern in the memory to calculate a single 
        forward FFT to get the expected tweezer pattern from given hologram.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        ax1.imshow(self.hologram, cmap="twilight_shifted")
        ax1.set_title("Hologram")

        ax2.imshow(
            abs(
                np.fft.fftshift(
                    np.fft.fft2(
                        self.Ui * np.exp(1j * self.hologram)
                    )
                )
            ) ** 2 * np.sum(self.mask),
            cmap="gray"
        )
        ax2.set_title("Image")
        ax2.set_xlim([self.width // 2 - 50, self.width // 2 + 50])
        ax2.set_ylim([self.height // 2 - 50, self.height // 2 + 50])
        fig.tight_layout()
        plt.show()
    
    def save(self, savepath):
        """
        Converts the hologram to 8bit and saves it to savepath.

        Parameters
        ----------
        savepath : string
            String where to save the hologram.
        """
        hologram_8bit = ((self.hologram + np.pi) / (2*np.pi) * 255) % 256
        phi_im = Image.fromarray(hologram_8bit.astype(np.uint8))
        phi_im.save(savepath)

    def loadFromHologram(self, hologram=None, n_spots=1, intensity_threshold=1e-5):
        """
        Loads the spots from the given hologram.

        Computes using the current illumination pattern a predicted
        intensity and takes the top n_spots spots out of this. If 
        these are bigger than a set threshold, the function extracts
        all values from the calculation and stores this as a spot.

        This detection leaves the spots unsorted, so one needs to resort
        afterwards. We do that here from low to high values in the
        unrolled coordinates.

        Parameters
        ----------
        hologram : np.array of floats or None
            Phase mask to be used. If None, the mask in the
            PMG memory will be used if available
        n_spots : int
            Maximal number of spots to load.
        intensity_threshold : float
            Minimal intensity that should be in the image for
            a spot to be loaded.
        """
        if hologram is None:
            if not np.any(self.hologram):
                print("No mask available in the memory to use for computing intensities")
                return
            hologram = self.hologram

        if hologram.shape != self.Ui.shape:
            print("Given phase mask does not match stored illumination pattern")
            return

        expected_field = np.fft.fftshift(np.fft.fft2(
            self.Ui * np.exp(1j * np.fft.fftshift(hologram))
        ))
        image = abs(
            expected_field
        ) ** 2 * n_spots
        idx_highest_spots = np.c_[np.unravel_index(np.argpartition(image.ravel(),-n_spots)[-n_spots:],image.shape)]

        self.spots = np.array([], dtype=np.float32)
        for y, x in idx_highest_spots:
            if image[y, x] < intensity_threshold:
                continue
            self.addSpot(
                roll(x-self.width//2, self.width), 
                roll(y-self.height//2, self.height), 
                0, image[y, x], np.angle(expected_field[y, x])
            )
        if self.spots.size > 0:
            self.spots[:, 3] /= np.mean(self.spots[:, 3])
            self.sortSpots()

    def readHologram(self,  hologram=None, n_spots=1, intensity_threshold=1e-5):
        """
        Reads the spots from the given hologram.

        Computes using the current illumination pattern a predicted
        intensity and takes the top n_spots spots out of this. If 
        these are bigger than a set threshold, the function extracts
        all values from the calculation  and returns these.

        This detection leaves the spots unsorted, so one needs to resort
        afterwards. We do that here from low to high values in the
        unrolled coordinates.

        Parameters
        ----------
        hologram : np.array of floats or None
            Phase mask to be used. If None, the mask in the
            PMG memory will be used if available
        n_spots : int
            Maximal number of spots to load.
        intensity_threshold : float
            Minimal intensity that should be in the image for
            a spot to be loaded.
        
        Returns
        -------
        return_spots : np.ndarray
            Array like spots containing the detected spots.
        """
        if hologram is None:
            if not np.any(self.hologram):
                print("No mask available in the memory to use for computing intensities")
                return
            hologram = self.hologram

        if hologram.shape != self.Ui.shape:
            print("Given phase mask does not match stored illumination pattern")
            return

        expected_field = np.fft.fftshift(np.fft.fft2(
            self.Ui * np.exp(1j * np.fft.fftshift(hologram))
        ))
        image = abs(
            expected_field
        ) ** 2 * n_spots
        idx_highest_spots = np.c_[np.unravel_index(np.argpartition(image.ravel(),-n_spots)[-n_spots:],image.shape)]

        return_spots = np.array([], dtype=np.float32)
        for y, x in idx_highest_spots:
            if image[y, x] < intensity_threshold:
                continue
            if return_spots.size == 0:
                return_spots = np.array(
                    [[roll(x-self.width//2, self.width), roll(y-self.height//2, self.height), 0, image[y, x], np.angle(expected_field[y, x])]], 
                    dtype=np.float32
                )
            else:
                return_spots = np.vstack([
                    return_spots, 
                    [roll(x-self.width//2, self.width), roll(y-self.height//2, self.height), 0, image[y, x], np.angle(expected_field[y, x])]
                ])

        if return_spots.size > 0:
            return_spots[:, 3] /= np.mean(return_spots[:, 3])

            unrolled = np.copy(return_spots)
            for i in range(unrolled.shape[0]):
                unrolled[i, 0] = unroll(return_spots[i, 0], self.width)
                unrolled[i, 1] = unroll(return_spots[i, 1], self.height)
            return_spots = return_spots[unrolled[:, 1].argsort()]

            unrolled = np.copy(return_spots)
            for i in range(unrolled.shape[0]):
                unrolled[i, 0] = unroll(return_spots[i, 0], self.width)
                unrolled[i, 1] = unroll(return_spots[i, 1], self.height)
            return_spots = return_spots[unrolled[:, 0].argsort()]
        return return_spots
    
    def __del__(self):
        """
        Destructor to remove the buffers from the GPU memory.
        """
        self.gpu.releaseBuffers()