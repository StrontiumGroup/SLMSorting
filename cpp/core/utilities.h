/*
* utilities.h
*
* This is the header file with defining constants useful for the operation of the
* PMG classes, and should be imported in every other header file. It also contains
* helper functions defined in the namespace PMG.
*/


#pragma once

// Define OpenCL version
#define CL_TARGET_OPENCL_VERSION 300
#define CL_HPP_TARGET_OPENCL_VERSION 300

// Include most important libraries
#include <iostream>
#include <vector>
#include <string.h>
#include <chrono>
#include <deque>
#include "CL/cl.h"
#include "vkFFT.h"
#include "graybmp.h"
#include <thread>

// Defines the maximum number of spots and thus length of reserved arrays.
#define MAX_NUMBER_SPOTS 30000
#define MAX_NUMBER_EVENTS 100
#define MAX_NUMBER_LAYER 200

// Indices for event pointers of OpenCL kernels.
#define INDEX_KNL_ASSEMBLE_FIELD_SLM 0
#define INDEX_KNL_ASSEMBLE_FIELD_TWEEZERS 1
#define INDEX_KNL_SET_AMPLITUDE_FIELD_SLM 2
#define INDEX_KNL_SET_AMPLITUDE_FIELD_TWEEZERS 3
#define INDEX_KNL_SET_PHASE_TWEEZERS 4
#define INDEX_KNL_GET_PHASE 5
#define INDEX_KNL_CALCULATE_AMPLITUDES 6
#define INDEX_KNL_UPDATE_TARGET 7
#define INDEX_KNL_APPLY_CONSTRAINTS_GSW_LOCKED 8
#define INDEX_KNL_APPLY_CONSTRAINTS_GS_LOCKED 9
#define INDEX_KNL_APPLY_CONSTRAINTS_GSW 10
#define INDEX_KNL_APPLY_CONSTRAINTS_GS 11
#define INDEX_KNL_ASSEMBLE_FIELD_SLM_2D 12
#define INDEX_KNL_SET_AMPLITUDE_FIELD_SLM_2D 13
#define INDEX_KNL_GET_PHASE_2D 14
#define INDEX_KNL_GET_PHASE_TWEEZERS 15
#define INDEX_KNL_RESET_FIELD 16
#define INDEX_KNL_SET_TARGET 17
#define INDEX_KNL_APPLY_FFT 18
#define INDEX_KNL_APPLY_IFFT 19

#define INDEX_WRITE_X_BUFFER 50
#define INDEX_WRITE_Y_BUFFER 51
#define INDEX_WRITE_AMPLITUDE_BUFFER 52
#define INDEX_WRITE_TWEEZER_PHASE_BUFFER 53
#define INDEX_WRITE_SHOULD_BE_USED_BUFFER 54
#define INDEX_WRITE_PHASE_BUFFER 55
#define INDEX_WRITE_TARGET_BUFFER 56
#define INDEX_WRITE_INPUT_LIGHT_BUFFER 57
#define INDEX_WRITE_CALCULATED_AMPLITUDE_BUFFER 58

#define INDEX_READ_X_BUFFER 80
#define INDEX_READ_Y_BUFFER 81
#define INDEX_READ_AMPLITUDE_BUFFER 82
#define INDEX_READ_TWEEZER_PHASE_BUFFER 83
#define INDEX_READ_SHOULD_BE_USED_BUFFER 84
#define INDEX_READ_PHASE_BUFFER 85
#define INDEX_READ_FIELD_BUFFER 86

#define INDEX_HOST_GET_POSITIONS 90
#define INDEX_GET_MASK_FROM_BUFFER 91

#define INDEX_BARRIER MAX_NUMBER_EVENTS-1

// Define SLM Board number
#define SLM_BOARD_NUMBER 1

// Define Sorting Machine parameters
#define SORTING_METHOD_AUTO 0U
#define SORTING_METHOD_HUNGARIAN 1U
#define SORTING_METHOD_COMPRESSION 2U
#define SORTING_METHOD_PARTITIONED_COMPRESSION 3U

#define PATHING_METHOD_DIRECT 0U

#define SORTING_MACHINE_SUCCESS 0
#define SORTING_MACHINE_PARTITION_FAILED 1

namespace PMG {
    struct Spot
    {
        /*
         * struct Spot
         *
         * Container for Spot data, including x-position, y-position,
         * amplitude, phase, and whether the spot is to be used in
         * calculations.
         */
        unsigned int x;
        unsigned int y;
        float a;
        float psi;
        bool used;
    };

    /*
    * PMG::roll()
    *
    * Used to map coordinates from [-512,511] to [0,1023]. This 
    * serves as a FFT_Shift for the Fourier transforms.
    */
    inline unsigned int roll(
        int x,
        unsigned int dimension
    ) {
        return x % dimension;
    };

    /*
    * PMG::roll()
    *
    * Used to map coordinates from [0,1023] to [-512,511]. Serves
    * as an inverse FFT_Shift for the Fourier transforms.
    */
    inline int unroll(
        unsigned int x,
        unsigned int dimension
    ) {
        int ret;
        if (x > dimension / 2) {
            ret = x - dimension;
        }
        else {
            ret = x;
        }
        return ret;
    };

    /*
    * PMG::TestResult
    *
    * A struct for returning all the timing data and results.
    */
    struct TestResult
    {
        double duration_sorting;
        double duration_calculation;
        double duration_total;
        double duration_display;
        double duration_transfer;
        unsigned int number_of_traps;
        unsigned int number_of_patterns;
        unsigned int sorting_method;
        unsigned int calculation_method;
        unsigned int pathing_method;
    };
}


inline void setup_opencl(
    unsigned int device_nr,
    cl_device_id *device,
    cl_context *context,
    cl_command_queue *queue
) {
    cl_uint numPlatforms;
    clGetPlatformIDs(0, 0, &numPlatforms);
    cl_platform_id* platforms = (cl_platform_id*) malloc (sizeof(cl_platform_id) * numPlatforms);
    clGetPlatformIDs(numPlatforms, platforms, 0);
    uint64_t k = 0;
    for (uint64_t j = 0; j < numPlatforms; j++) {
        cl_uint numDevices;
        clGetDeviceIDs(platforms[j], CL_DEVICE_TYPE_ALL, 0, 0, &numDevices);
        cl_device_id* deviceList = (cl_device_id*) malloc (sizeof(cl_device_id) * numDevices);
        clGetDeviceIDs(platforms[j], CL_DEVICE_TYPE_ALL, numDevices, deviceList, 0);

        for (uint64_t i = 0; i < numDevices; i++) {
            if (k == device_nr) {
                *device = deviceList[i];
                *context = clCreateContext(NULL, 1, device, NULL, NULL, NULL);
                
                cl_command_queue_properties properties[] = { CL_QUEUE_PROPERTIES, /*CL_QUEUE_PROFILING_ENABLE |*/ CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE ,0 };
                *queue = clCreateCommandQueueWithProperties(*context, *device,  properties , NULL);

                i=numDevices;
                j=numPlatforms;
            }
            else {
                k++;
            }
        }
        free(deviceList);
    }
    free(platforms);
}
