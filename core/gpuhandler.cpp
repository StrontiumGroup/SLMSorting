/*
* gpuhandler.cpp
*/

#include "gpuhandler.h"

namespace PMG {
    void GPU_Handler::build_kernels() {
        knl_assemble_field_SLM = clCreateKernel(program, "assemble_field", &result);
        clSetKernelArg(knl_assemble_field_SLM, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_assemble_field_SLM, 1, sizeof(input_light_buffer), &input_light_buffer);
        clSetKernelArg(knl_assemble_field_SLM, 2, sizeof(phase_buffer), &phase_buffer);
        knl_assemble_field_SLM_2d = clCreateKernel(program, "assemble_field_2d", &result);
        clSetKernelArg(knl_assemble_field_SLM_2d, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_assemble_field_SLM_2d, 1, sizeof(input_light_buffer), &input_light_buffer);
        clSetKernelArg(knl_assemble_field_SLM_2d, 2, sizeof(phase_buffer), &phase_buffer);
        knl_assemble_field_tweezers = clCreateKernel(program, "assemble_field", &result);
        clSetKernelArg(knl_assemble_field_tweezers, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_assemble_field_tweezers, 1, sizeof(target_buffer), &target_buffer);
        clSetKernelArg(knl_assemble_field_tweezers, 2, sizeof(phase_buffer), &phase_buffer);
        knl_set_amplitude_field_SLM = clCreateKernel(program, "set_amplitude_field", &result);
        clSetKernelArg(knl_set_amplitude_field_SLM, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_set_amplitude_field_SLM, 1, sizeof(input_light_buffer), &input_light_buffer);
        knl_set_amplitude_field_SLM_2d = clCreateKernel(program, "set_amplitude_field_2d", &result);
        clSetKernelArg(knl_set_amplitude_field_SLM_2d, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_set_amplitude_field_SLM_2d, 1, sizeof(input_light_buffer), &input_light_buffer);
        knl_set_amplitude_field_tweezers = clCreateKernel(program, "set_amplitude_field", &result);
        clSetKernelArg(knl_set_amplitude_field_tweezers, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_set_amplitude_field_tweezers, 1, sizeof(target_buffer), &target_buffer);
        knl_apply_constraints_GSW_locked = clCreateKernel(program, "apply_constraints_GSW_locked", &result);
        clSetKernelArg(knl_apply_constraints_GSW_locked, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_apply_constraints_GSW_locked, 1, sizeof(x_buffer), &x_buffer);
        clSetKernelArg(knl_apply_constraints_GSW_locked, 2, sizeof(y_buffer), &y_buffer);
        clSetKernelArg(knl_apply_constraints_GSW_locked, 3, sizeof(should_be_used_buffer), &should_be_used_buffer);
        clSetKernelArg(knl_apply_constraints_GSW_locked, 4, sizeof(target_amplitude_buffer), &target_amplitude_buffer);
        clSetKernelArg(knl_apply_constraints_GSW_locked, 5, sizeof(target_phase_buffer), &target_phase_buffer);
        knl_apply_constraints_GS_locked = clCreateKernel(program, "apply_constraints_GS_locked", &result);
        clSetKernelArg(knl_apply_constraints_GS_locked, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_apply_constraints_GS_locked, 1, sizeof(x_buffer), &x_buffer);
        clSetKernelArg(knl_apply_constraints_GS_locked, 2, sizeof(y_buffer), &y_buffer);
        clSetKernelArg(knl_apply_constraints_GS_locked, 3, sizeof(should_be_used_buffer), &should_be_used_buffer);
        clSetKernelArg(knl_apply_constraints_GS_locked, 4, sizeof(target_amplitude_buffer), &target_amplitude_buffer);
        clSetKernelArg(knl_apply_constraints_GS_locked, 5, sizeof(target_phase_buffer), &target_phase_buffer);
        knl_apply_constraints_GSW = clCreateKernel(program, "apply_constraints_GSW", &result);
        clSetKernelArg(knl_apply_constraints_GSW, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_apply_constraints_GSW, 1, sizeof(x_buffer), &x_buffer);
        clSetKernelArg(knl_apply_constraints_GSW, 2, sizeof(y_buffer), &y_buffer);
        clSetKernelArg(knl_apply_constraints_GSW, 3, sizeof(should_be_used_buffer), &should_be_used_buffer);
        clSetKernelArg(knl_apply_constraints_GSW, 4, sizeof(target_amplitude_buffer), &target_amplitude_buffer);
        clSetKernelArg(knl_apply_constraints_GSW, 5, sizeof(calculated_amplitude_buffer), &calculated_amplitude_buffer);
        knl_get_phase = clCreateKernel(program, "get_phase", &result);
        clSetKernelArg(knl_get_phase, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_get_phase, 1, sizeof(return_buffer), &return_buffer);
        knl_get_phase_2d = clCreateKernel(program, "get_phase_2d", &result);
        clSetKernelArg(knl_get_phase_2d, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_get_phase_2d, 1, sizeof(return_buffer), &return_buffer);
        knl_update_target = clCreateKernel(program, "update_target", &result);
        clSetKernelArg(knl_update_target, 0, sizeof(target_buffer), &target_buffer);
        clSetKernelArg(knl_update_target, 1, sizeof(x_buffer), &x_buffer);
        clSetKernelArg(knl_update_target, 2, sizeof(y_buffer), &y_buffer);
        clSetKernelArg(knl_update_target, 3, sizeof(should_be_used_buffer), &should_be_used_buffer);
        clSetKernelArg(knl_update_target, 4, sizeof(target_amplitude_buffer), &target_amplitude_buffer);
        clSetKernelArg(knl_update_target, 5, sizeof(calculated_amplitude_buffer), &calculated_amplitude_buffer);
        knl_calculate_amplitudes = clCreateKernel(program, "calculate_amplitudes", &result);
        clSetKernelArg(knl_calculate_amplitudes, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_calculate_amplitudes, 1, sizeof(x_buffer), &x_buffer);
        clSetKernelArg(knl_calculate_amplitudes, 2, sizeof(y_buffer), &y_buffer);
        clSetKernelArg(knl_calculate_amplitudes, 3, sizeof(should_be_used_buffer), &should_be_used_buffer);
        clSetKernelArg(knl_calculate_amplitudes, 4, sizeof(calculated_amplitude_buffer), &calculated_amplitude_buffer);
        knl_get_phase_tweezers = clCreateKernel(program, "get_phase_tweezers", &result);
        clSetKernelArg(knl_get_phase_tweezers, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_get_phase_tweezers, 1, sizeof(x_buffer), &x_buffer);
        clSetKernelArg(knl_get_phase_tweezers, 2, sizeof(y_buffer), &y_buffer);
        clSetKernelArg(knl_get_phase_tweezers, 3, sizeof(target_phase_buffer), &target_phase_buffer);
        knl_set_phase_tweezers = clCreateKernel(program, "set_phase_tweezers", &result);
        clSetKernelArg(knl_set_phase_tweezers, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_set_phase_tweezers, 1, sizeof(x_buffer), &x_buffer);
        clSetKernelArg(knl_set_phase_tweezers, 2, sizeof(y_buffer), &y_buffer);
        clSetKernelArg(knl_set_phase_tweezers, 3, sizeof(should_be_used_buffer), &should_be_used_buffer);
        clSetKernelArg(knl_set_phase_tweezers, 4, sizeof(target_phase_buffer), &target_phase_buffer);
        knl_reset_field = clCreateKernel(program, "reset_field", &result);
        clSetKernelArg(knl_reset_field, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_reset_field, 1, sizeof(x_buffer), &x_buffer);
        clSetKernelArg(knl_reset_field, 2, sizeof(y_buffer), &y_buffer);
        knl_set_target = clCreateKernel(program, "set_target", &result);
        clSetKernelArg(knl_set_target, 0, sizeof(target_buffer), &target_buffer);
        clSetKernelArg(knl_set_target, 1, sizeof(x_buffer), &x_buffer);
        clSetKernelArg(knl_set_target, 2, sizeof(y_buffer), &y_buffer);
        clSetKernelArg(knl_set_target, 3, sizeof(should_be_used_buffer), &should_be_used_buffer);
        clSetKernelArg(knl_set_target, 4, sizeof(target_amplitude_buffer), &target_amplitude_buffer);
        knl_get_phase_with_corrections = clCreateKernel(program, "get_phase_with_corrections", &result);
        clSetKernelArg(knl_get_phase_with_corrections, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_get_phase_with_corrections, 1, sizeof(corrections_buffer), &corrections_buffer);
        clSetKernelArg(knl_get_phase_with_corrections, 2, sizeof(return_buffer), &return_buffer);
        knl_get_phase_with_corrections_with_roll = clCreateKernel(program, "get_phase_with_corrections_with_roll", &result);
        clSetKernelArg(knl_get_phase_with_corrections_with_roll, 0, sizeof(field_buffer), &field_buffer);
        clSetKernelArg(knl_get_phase_with_corrections_with_roll, 1, sizeof(corrections_buffer), &corrections_buffer);
        clSetKernelArg(knl_get_phase_with_corrections_with_roll, 2, sizeof(return_buffer), &return_buffer);
        clSetKernelArg(knl_get_phase_with_corrections_with_roll, 3, sizeof(cl_int), &mask_offset_x);
        clSetKernelArg(knl_get_phase_with_corrections_with_roll, 4, sizeof(cl_int), &mask_offset_y);

        if (result != CL_SUCCESS) {
            std::cout << "Ran into problems building kernels. OpenCL return code: " << result << std::endl;
            exit(1);
        }
        //std::cout << "Built kernels" << std::endl;
    }

    /*
    * GPU_Handler::inititate_spot_buffers_of_length()
    *
    * Creates the relevant buffers of size length.
    *
    * Parameters:
    * length : unsigned int
    *      Length of the vectors, i.e. number of spots.
    */
    void GPU_Handler::initiate_spot_buffers_with_length(
        unsigned int length)
    {
        if (buffer_size_spots > 0) {
            result = clReleaseMemObject(x_buffer);
            clReleaseMemObject(y_buffer);
            clReleaseMemObject(target_amplitude_buffer);
            clReleaseMemObject(target_phase_buffer);
            clReleaseMemObject(should_be_used_buffer);
            clReleaseMemObject(calculated_amplitude_buffer);

            if (result != CL_SUCCESS) {
                std::cout << "Ran into problems removing spot buffers. OpenCL return code: " << result << std::endl;
                system("pause");
                exit(1);
            }
            //std::cout << "Reset buffers to correct length" << std::endl;
        }    

        buffer_size_spots = length;

        x_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(unsigned int)*buffer_size_spots, nullptr, &result);
        y_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(unsigned int)*buffer_size_spots, nullptr, &result);
        target_amplitude_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(float)*buffer_size_spots, nullptr, &result);
        target_phase_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(float)*buffer_size_spots, nullptr, &result);
        should_be_used_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(unsigned int)*buffer_size_spots, nullptr, &result);
        calculated_amplitude_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(float)*buffer_size_spots, nullptr, &result);
       
        if (result != CL_SUCCESS) {
            std::cout << "Ran into problems setting spot buffers. OpenCL return code: " << result << std::endl;
            system("pause");
            exit(1);
        }
        //std::cout << "Reset buffers to correct length" << std::endl;
    }

    /*
    * PMG::GPU_Handler::run_FFT()
    *
    * Runs the VkFFT application once in forward direction.
    */
    void GPU_Handler::run_FFT(bool debug_timing)
    {
        if (debug_timing) {
            vkfft_time_start = std::chrono::steady_clock::now();
        }
        clWaitForEvents(1, &event_list[INDEX_BARRIER]);
        clSetUserEventStatus(event_list[INDEX_KNL_APPLY_FFT], CL_SUBMITTED);
        VkFFTAppend(&vkapp, -1, &vkparams);
        clSetUserEventStatus(event_list[INDEX_KNL_APPLY_FFT], CL_COMPLETE);
        if (debug_timing) {
            clFinish(*queue_ptr);
            vkfft_time_end = std::chrono::steady_clock::now();

            exec_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(vkfft_time_end - vkfft_time_start).count();
            printf("run_FFT: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }
    }

    /*
    * PMG::GPU_Handler::run_iFFT()
    *
    * Runs the VkFFT application once in backward direction.
    */
    void GPU_Handler::run_iFFT(bool debug_timing)
    {
        if (debug_timing) {
            vkfft_time_start = std::chrono::steady_clock::now();
        }
        clWaitForEvents(1, &event_list[INDEX_BARRIER]);
        clSetUserEventStatus(event_list[INDEX_KNL_APPLY_IFFT], CL_SUBMITTED);
        VkFFTAppend(&vkapp, 1, &vkparams);
        clSetUserEventStatus(event_list[INDEX_KNL_APPLY_IFFT], CL_COMPLETE);
        if (debug_timing) {
            clFinish(*queue_ptr);
            vkfft_time_end = std::chrono::steady_clock::now();

            exec_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(vkfft_time_end - vkfft_time_start).count();
            printf("run_iFFT: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }
    }

    /*
     * GPU_Handler::run_knl_assemble_field_tweezers()
     * 
     * Assembles the field at the SLM frame, so uses both the 
     * phase buffer and the input light buffer.
     */
    void GPU_Handler::run_knl_assemble_field_SLM(bool debug_timing)
    {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_assemble_field_SLM,
            1,  /* 1D kernel for now */
            0,
            &buffer_size_fields,
            0,
            0,
            NULL,
            &event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_assemble_field_SLM: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }   

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_assemble_field_SLM. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     * GPU_Handler::run_knl_assemble_field_SLM_2d()
     * 
     * Assembles the field at the SLM frame, so uses both the 
     * phase buffer and the input light buffer.
     */
    void GPU_Handler::run_knl_assemble_field_SLM_2d(bool debug_timing)
    {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_assemble_field_SLM_2d,
            2,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_fields_2d,
            0,
            1,
            &event_list[INDEX_BARRIER],
            &event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM_2D]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM_2D]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM_2D], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM_2D], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_assemble_field_SLM_2d: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }   

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_assemble_field_SLM_2d. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }

    }

    /*
     * GPU_Handler::run_knl_assemble_field_tweezers()
     * 
     * Assembles the field at the tweezers frame, so uses both the 
     * phase buffer and the target buffer.
     */
    void GPU_Handler::run_knl_assemble_field_tweezers(bool debug_timing)
    {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_assemble_field_tweezers,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_fields,
            0,
            0,
            NULL,
            &event_list[INDEX_KNL_ASSEMBLE_FIELD_TWEEZERS]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_ASSEMBLE_FIELD_TWEEZERS]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_ASSEMBLE_FIELD_TWEEZERS], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_ASSEMBLE_FIELD_TWEEZERS], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_assemble_field_tweezers: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_assemble_field_tweezers. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     * GPU_Handler::run_knl_assemble_field_tweezers()
     * 
     * Assembles the field at the tweezers frame, so uses both the 
     * phase buffer and the target buffer.
     */
    void GPU_Handler::run_knl_update_target(bool debug_timing)
    {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_update_target,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_spots,
            0,
            0,
            NULL,
            &event_list[INDEX_KNL_UPDATE_TARGET]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_UPDATE_TARGET]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_UPDATE_TARGET], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_UPDATE_TARGET], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_assemble_field_tweezers: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_update_target. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     * GPU_Handler::run_knl_assemble_field_tweezers()
     * 
     * Assembles the field at the tweezers frame, so uses both the 
     * phase buffer and the target buffer.
     */
    void GPU_Handler::run_knl_set_target(bool debug_timing)
    {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_set_target,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_spots,
            0,
            0,
            NULL,
            &event_list[INDEX_KNL_SET_TARGET]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_SET_TARGET]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_SET_TARGET], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_SET_TARGET], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_set_target: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_set_target. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

     /*
     * GPU_Handler::run_knl_reset_field()
     * 
     * Resets the field at the spots to 0. Useful for sorting.
     */
    void GPU_Handler::run_knl_reset_field(bool debug_timing)
    {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_reset_field,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_spots,
            0,
            1,
            &event_list[INDEX_BARRIER],
            &event_list[INDEX_KNL_RESET_FIELD]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_RESET_FIELD]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_RESET_FIELD], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_RESET_FIELD], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_reset_field: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_reset_field. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     * GPU_Handler::run_knl_get_phase()
     *
     * Runs kernel to store the phase information of the field in the return buffer.
     */
    void GPU_Handler::run_knl_get_phase(bool debug_timing)
    {
        clWaitForEvents(1, &event_list[INDEX_GET_MASK_FROM_BUFFER]);

        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_get_phase,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_fields,
            0,
            1,
            &event_list[INDEX_BARRIER],
            &event_list[INDEX_KNL_GET_PHASE]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_GET_PHASE]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_GET_PHASE], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_GET_PHASE], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_get_phase: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_get_phase. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     * GPU_Handler::run_knl_get_phase_with_corrections()
     *
     * Runs kernel to store the phase information of the field in the return buffer.
     */
    void GPU_Handler::run_knl_get_phase_with_corrections(bool debug_timing)
    {
        clWaitForEvents(1, &event_list[INDEX_BARRIER]);

        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_get_phase_with_corrections,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_fields,
            0,
            1,
            &event_list[INDEX_BARRIER],
            &event_list[INDEX_KNL_GET_PHASE]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_GET_PHASE]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_GET_PHASE], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_GET_PHASE], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("run_knl_get_phase_with_corrections: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in run_knl_get_phase_with_corrections. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     * GPU_Handler::run_knl_get_phase_with_corrections_with_roll()
     *
     * Runs kernel to store the phase information of the field in the return buffer.
     */
    void GPU_Handler::run_knl_get_phase_with_corrections_with_roll(
        int offset_x,
        int offset_y,
        bool debug_timing
    ){
        mask_offset_x = offset_x;
        mask_offset_y = offset_y;
        clSetKernelArg(knl_get_phase_with_corrections_with_roll, 3, sizeof(int), &mask_offset_x);
        clSetKernelArg(knl_get_phase_with_corrections_with_roll, 4, sizeof(int), &mask_offset_y);


        clWaitForEvents(1, &event_list[INDEX_BARRIER]);

        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_get_phase_with_corrections_with_roll,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_fields,
            0,
            1,
            &event_list[INDEX_BARRIER],
            &event_list[INDEX_KNL_GET_PHASE]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_GET_PHASE]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_GET_PHASE], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_GET_PHASE], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("run_knl_get_phase_with_corrections: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in run_knl_get_phase_with_corrections_with_roll. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     * GPU_Handler::run_knl_get_phase_2d()
     *
     * Runs kernel to store the phase information of the field in the return buffer.
     */
    void GPU_Handler::run_knl_get_phase_2d(bool debug_timing)
    {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_get_phase_2d,
            2,  /* 1D kernel for now */
            0,
            buffer_size_fields_2d,
            0,
            0,
            NULL,
            &event_list[INDEX_KNL_GET_PHASE_2D]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_GET_PHASE_2D]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_GET_PHASE_2D], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_GET_PHASE_2D], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_get_phase_2d: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_get_phase_2d. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     *
     *
     */
    void GPU_Handler::run_knl_get_phase_tweezers(bool debug_timing) {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_get_phase_tweezers,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_spots,
            0,
            1,
            &event_list[INDEX_BARRIER],
            &event_list[INDEX_KNL_GET_PHASE_TWEEZERS]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_GET_PHASE_TWEEZERS]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_GET_PHASE_TWEEZERS], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_GET_PHASE_TWEEZERS], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_get_phase_tweezers: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_get_phase_tweezers. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     *
     *
     */
    void GPU_Handler::run_knl_set_phase_tweezers(bool debug_timing) {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_set_phase_tweezers,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_spots,
            0,
            1,
            &event_list[INDEX_BARRIER],
            &event_list[INDEX_KNL_SET_PHASE_TWEEZERS]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_SET_PHASE_TWEEZERS]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_SET_PHASE_TWEEZERS], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_SET_PHASE_TWEEZERS], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_set_phase_tweezers: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_set_phase_tweezers. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     * GPU_Handler::run_knl_set_ampltiude_field_SLM
     *
     * Runs the kernel for setting the amplitude of the field at the SLM frame 
     * (i.e. updating the amplitudes to the input light buffer).
     */
    void GPU_Handler::run_knl_set_amplitude_field_SLM(bool debug_timing) {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_set_amplitude_field_SLM,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_fields,
            0,
            0,
            NULL,
            &event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_SLM]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_SLM]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_SLM], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_SLM], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_set_amplitude_field_SLM: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_set_amplitude_field_SLM. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     * GPU_Handler::get_amplitudes_tweezers()
     *
     * Runs the kernel to get the amplitude values at teh position of 
     * the tweezers. 
     * */
    void GPU_Handler::run_knl_calculate_amplitudes(
        bool debug_timing
    ) {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_calculate_amplitudes,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_spots,
            0,
            0,
            NULL,
            &event_list[INDEX_KNL_CALCULATE_AMPLITUDES]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_CALCULATE_AMPLITUDES]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_CALCULATE_AMPLITUDES], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_CALCULATE_AMPLITUDES], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("run_knl_calculate_amplitudes: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in run_knl_calculate_amplitudes. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

        /*
     * GPU_Handler::run_knl_set_ampltiude_field_SLM_2d
     *
     * Runs the kernel for setting the amplitude of the field at the SLM frame 
     * (i.e. updating the amplitudes to the input light buffer).
     */
    void GPU_Handler::run_knl_set_amplitude_field_SLM_2d(bool debug_timing) {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_set_amplitude_field_SLM_2d,
            2,  /* 1D kernel for now */
            0,
            buffer_size_fields_2d,
            0,
            1,
            &event_list[INDEX_BARRIER],
            &event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_SLM_2D]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_SLM_2D]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_SLM_2D], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_SLM_2D], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_set_amplitude_field_SLM_2d: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_set_amplitude_field_SLM. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    /*
     * GPU_Handler::run_knl_set_ampltiude_field_tweezers
     *
     * Runs the kernel for setting the amplitude of the field at the tweezer frame 
     * (i.e. updating the amplitudes to the target buffer).
     */
    void GPU_Handler::run_knl_set_amplitude_field_tweezers(bool debug_timing) {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_set_amplitude_field_tweezers,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_fields,
            0,
            1,
            &event_list[INDEX_BARRIER],
            &event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_TWEEZERS]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_TWEEZERS]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_TWEEZERS], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_SET_AMPLITUDE_FIELD_TWEEZERS], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_set_amplitude_field_tweezers: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_set_amplitude_field_tweezers. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    void GPU_Handler::run_knl_apply_constraints_GSW(bool debug_timing) {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_apply_constraints_GSW,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_spots,
            0,
            0,
            NULL,
            &event_list[INDEX_KNL_APPLY_CONSTRAINTS_GSW]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_APPLY_CONSTRAINTS_GSW]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_APPLY_CONSTRAINTS_GSW], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_APPLY_CONSTRAINTS_GSW], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_apply_constraints_GSW: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_apply_constraints_GSW. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    void GPU_Handler::run_knl_apply_constraints_GSW_locked(bool debug_timing) {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_apply_constraints_GSW_locked,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_spots,
            0,
            0,
            NULL,
            &event_list[INDEX_KNL_APPLY_CONSTRAINTS_GSW_LOCKED]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_APPLY_CONSTRAINTS_GSW_LOCKED]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_APPLY_CONSTRAINTS_GSW_LOCKED], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_APPLY_CONSTRAINTS_GSW_LOCKED], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_apply_constraints_GSW_locked: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_apply_constraints_GSW_locked. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }

    void GPU_Handler::run_knl_apply_constraints_GS_locked(bool debug_timing) {
        result |= clEnqueueNDRangeKernel(
            *queue_ptr,
            knl_apply_constraints_GS_locked,
            1,  /* 1D kernel for now */
            0,
            (const size_t *) &buffer_size_spots,
            0,
            1,
            &event_list[INDEX_BARRIER],
            &event_list[INDEX_KNL_APPLY_CONSTRAINTS_GS_LOCKED]
        );

        if (debug_timing) {
            clWaitForEvents(1, &event_list[INDEX_KNL_APPLY_CONSTRAINTS_GS_LOCKED]);

            clGetEventProfilingInfo(event_list[INDEX_KNL_APPLY_CONSTRAINTS_GS_LOCKED], CL_PROFILING_COMMAND_START, sizeof(time_start), &time_start, NULL);
            clGetEventProfilingInfo(event_list[INDEX_KNL_APPLY_CONSTRAINTS_GS_LOCKED], CL_PROFILING_COMMAND_END, sizeof(time_end), &time_end, NULL);

            exec_time_in_ns = time_end - time_start;
            printf("knl_apply_constraints_GS_locked: OpenCl Execution time is: %0.3f milliseconds \n \n", exec_time_in_ns / 1000000.0);
        }

        if (result != CL_SUCCESS) {
            std::cout << "Error in knl_apply_constraints_GS_locked. OpenCL code: " << result << std::endl;
            system("pause");
            exit(1);
        }
    }
}