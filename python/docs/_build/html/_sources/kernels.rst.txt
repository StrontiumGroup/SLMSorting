Kernels
=======

Below follows a list of the kernels that are used in the GPU to produce the patterns. These are all defined in the `buildKernels()` function of the `GPUHandler`. The kernels in the python code are closely related to the ones in the C++ version. Differences are mentioned.

Extensions
----------

Let us start with the required extensions. The code uses 64-bit floats for the complex numbers and 32-bit floats for the general calculation. Therefore, one needs to load `cl_khr_fp64`, define the `PYOPENCL_DEFINE_CDOUBLE` and load the header file. In addition to this, we also set the width and height of the SLM chip from the parent memory as a kernel definition, so that we can use it to locate pixels. In the kernels, this looks like:

.. code-block:: C 
    :linenos:

    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    #define PYOPENCL_DEFINE_CDOUBLE
    #define SLM_WIDTH """+str(self.parent.width)+"""
    #define SLM_HEIGHT """+str(self.parent.height)+"""
    #include <pyopencl-complex.h>

.. note:: All kernels are 1D. It could be that programming in 2D speeds up, but so far, we have not seen this.

Now we go over the relevant kernels one by one.

get_intensity()
---------------

The simplest kernel is to get the intensity of the electric field at a specific location. This is used in the weighting. It simply looks at the `field` buffer and looks at location :math:`(y, x)`. 

.. note:: Note that the choice of :math:`x` and :math:`y` leads to contradictory indices compared to :math:`(r, c)` most of the times. Since humans often think in :math:`(x, y)`, we have adopted these coordinates, but in the calculations this then translates to :math:`(y, x)` because matrices use rows and columns.

The resulting amplitude squared is saved in a vector intensity of length :math:`N`. This can efficiently be read out.

.. code-block:: C
    :linenos:

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

get_target()
------------

Sets the target intensity of the spots at locations :math:`(y, x)` with their target intensities. This is not done in the C++ code, because it is only useful in the WGS algorithm. Otherwise, we can get by with simply defining the desired electric field, see below with the shortcut, but here we need to suppress higher order effects as well and thus have an intensity pattern that contains also the zeros where we do not want spots.

.. code-block:: C
    :linenos:

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

update_target()
---------------

As part of the WGS algorithm, the target needs to be updated. Here, we need to be careful to not let the final outcome always become 1. For practical reasons, we also want to be able to produce non-uniform pattern, e.g. when iterating experimental aberrations. The basic update conditions are defined as follows:

.. math:: T_n = T_{n-1} \sqrt{\frac{T_{n-1}}{I_{n-1}}} \frac{1}{N}\sum_{i=1}^{N}\sqrt{\frac{I_{n-1}}{T_{n-1}}},


where :math:`T_n` is the target intensity for the :math:`n`-th iteration, and :math:`I_n` is the measured intensity. The sum is a normalizing average over the :math:`N` desired tweezers, to ensure the calculation stays within bounds.

.. code-block:: C
    :linenos:

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

.. note:: To calculate the average, each worker has to loop over all :math:`N` intensities. For large number of tweezers, it may be better to first calculate the averages, e.g. by transferring the information to the CPU, and then passing only a single normalization number to the kernel. 

Again, because this kernels is part of the WGS algorithm. It has no C++ counterpart in the current implementation.

set_field_from_spots()
----------------------

Creates the electric field amplitude and phase at the desired coordinates :math:`(y, x)`. This is useful as a starting point of the iteration when no hologram is provided. A slightly different variant is the core of the interpolation kernel.

.. code-block:: C
    :linenos:

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

update_field_shortcut()
-----------------------

The update kernel for the linear interpolation algorithm, which is the same for the C++ and the python code. The idea is to keep the update conditions as simple as possible. In this case, we only update the :math:`N` coordinate pixels of the electric field in the `field` buffer. To minimize any effects from residual amplitudes that are present from earlier calculations, we add a very large factor to the amplitude, effectively drowning out all information in the other pixels. This allows the update to only loop over :math:`N` coordinates, instead of the whole :math:`1024\times1024` chip size.

.. note:: An additional benefit is that this removes the need to perform an FFT back and removing the information of the previous hologram. At the SLM plane, the field will consist of a very small amplitude due to the spread normalized illumination and phases corresponding to the hologram. By using this kernel, we neglect all that information and only obtain a field that essentially looks like the field in the tweezer focus. Avoiding the single FFT does not save much time because the VkFFT implementation is very fast in C++ (~ :math:`20\,\mathrm{\mu s}`), but it alleviates the problem of having to remove the previous :math:`N` coordinates when moving tweezers. 

.. code-block:: C
    :linenos:

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