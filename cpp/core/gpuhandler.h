/*
* gpuhandler.h
*
*
*/

#pragma once

// Load core libraries
#include "utilities.h"

namespace PMG {
    // ------------------------------------------- PMG :: GPU_Handler ----------------------------------- //

    class GPU_Handler
    {
    public:
        cl_ulong time_start;
        cl_ulong time_end;
        std::chrono::steady_clock::time_point vkfft_time_start;
        std::chrono::steady_clock::time_point vkfft_time_end;
        double exec_time_in_ns;
        double total_exec_time_in_ns;

        VkFFTConfiguration vkconfiguration = {};
        VkFFTApplication vkapp = {};
        VkFFTLaunchParams vkparams = {};

        // Reserved for events of kernels. See definitions in the structs.
        cl_event event_list[MAX_NUMBER_EVENTS];

        cl_int result;
        cl_int mask_offset_x;
        cl_int mask_offset_y;

        cl_device_id *device_ptr;
        cl_context *context_ptr;
        cl_command_queue *queue_ptr;
        cl_program program;

        cl_mem x_buffer;
        cl_mem y_buffer;
        cl_mem target_amplitude_buffer;
        cl_mem target_phase_buffer;
        cl_mem should_be_used_buffer;
        cl_mem calculated_amplitude_buffer;

        cl_mem field_buffer;
        cl_mem phase_buffer;
        cl_mem input_light_buffer;
        cl_mem target_buffer;
        cl_mem return_buffer;
        cl_mem corrections_buffer;

        size_t buffer_size_spots;
        size_t buffer_size_fields;
        size_t buffer_size_fields_2d[2];
        unsigned int workgroups_per_row;
        size_t group_pattern[2]; 
        size_t local_size[2];

        GPU_Handler(
            cl_device_id *set_device_ptr,
            cl_context *set_context_ptr,
            cl_command_queue *set_queue_ptr,
            unsigned int *pmg_width,
            unsigned int *pmg_height,
            const char *kernel_str)
        {
            result = CL_SUCCESS;

            device_ptr = set_device_ptr;
            context_ptr = set_context_ptr;
            queue_ptr = set_queue_ptr;

            buffer_size_fields_2d[0] = *pmg_width;
            buffer_size_fields_2d[1] = *pmg_height;
            buffer_size_fields = *pmg_width * *pmg_height;
            buffer_size_spots = 0;
            workgroups_per_row = 16;
            group_pattern[0] = workgroups_per_row;
            group_pattern[1] = buffer_size_fields_2d[1]; //My workgroups pattern
            local_size[0] = buffer_size_fields_2d[0]/group_pattern[0];
            local_size[1] = buffer_size_fields_2d[1]/group_pattern[1];
            mask_offset_x = 0;
            mask_offset_y = 0;

            event_list[INDEX_KNL_APPLY_FFT] = clCreateUserEvent(*context_ptr, NULL);
            event_list[INDEX_KNL_APPLY_IFFT] = clCreateUserEvent(*context_ptr, NULL);
            event_list[INDEX_BARRIER]= clCreateUserEvent(*context_ptr, NULL);

            
            // Allocate the field buffers already. 
            //float *field_values = (float*)malloc((sizeof(float) * 2 * buffer_size_fields));
            field_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, 2 * sizeof(float) * buffer_size_fields,  nullptr, &result);
            target_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(float) * buffer_size_fields, nullptr, &result);
            input_light_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(float) * buffer_size_fields, nullptr, &result);
            phase_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(float) * buffer_size_fields, nullptr, &result);
            return_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(uint8_t) * buffer_size_fields, nullptr, &result);
            corrections_buffer = clCreateBuffer(*context_ptr, CL_MEM_READ_WRITE, sizeof(float) * buffer_size_fields, nullptr, &result);

            if (result != CL_SUCCESS) {
                std::cout << "Ran into problems intiating field buffers. OpenCL return code: " << result << std::endl;
                system("pause");
                exit(1);
            }
            //std::cout << "Initiated buffers" << std::endl;

            // Initialize VkFFT
            vkconfiguration = {};
            vkapp = {};
            vkparams = {};

            // Here we can choose between 2D and 1D FFTs. Not sure
            // which is faster yet.
            unsigned int dim = 2;
            if (dim == 1)
            {
                vkconfiguration.FFTdim = 1;
                vkconfiguration.size[0] = *pmg_height * *pmg_width;
            }
            else if (dim == 2)
            {
                vkconfiguration.FFTdim = 2;
                vkconfiguration.size[0] = *pmg_width;
                vkconfiguration.size[1] = *pmg_height;
            }
            // Note that Visual Studio doesn't see the DVKFFT_BACKEND=3 here.
            vkconfiguration.device = device_ptr;
            vkconfiguration.context = context_ptr;
            vkconfiguration.buffer = &field_buffer;
            vkconfiguration.numberBatches = 1;
            vkconfiguration.doublePrecision = 0;
            vkconfiguration.normalize = 1;
            vkparams.buffer = &field_buffer;
            vkparams.commandQueue = queue_ptr;
            result |= initializeVkFFT(&vkapp, vkconfiguration);
            if (result != CL_SUCCESS) {
                std::cout << "Ran into problems initializing VkFFT. OpenCL return code: " << result << std::endl;
                system("pause");
                exit(1);
            }
            //std::cout << "Initialized VkFFT" << std::endl;

            // Load the kernels.
            char *source_str;
            size_t source_size;

            FILE* fp = NULL;
	        fp = fopen(kernel_str, "rb");
            fseek(fp, 0, SEEK_END);
            source_size = ftell(fp);
            fseek(fp, 0, SEEK_SET);
            source_str = (char*)malloc(sizeof(char)*(source_size)); 
            fread(source_str, 1, source_size, fp);

            program = clCreateProgramWithSource(*context_ptr, 1, (const char**) &source_str, NULL, &result);
            clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
            if (result != CL_SUCCESS) {
                std::cout << "Ran into problems building program. OpenCL return code: " << result << std::endl;
                system("pause");
                exit(1);
            }
            //std::cout << "Built program" << std::endl;
        };

        void build_kernels();
        void run_FFT(bool debug_timing = false);
        void run_iFFT(bool debug_timing = false);
        void run_knl_get_amplitudes(bool debug_timing = false);
        void run_knl_update_target(bool debug_timing = false);
        void run_knl_calculate_amplitudes(bool debug_timing = false);
        void run_knl_get_phase(bool debug_timing = false);
        void run_knl_get_phase_2d(bool debug_timing = false);
        void run_knl_get_corrected_phase(bool debug_timing = false);
        void run_knl_assemble_field_SLM(bool debug_timing = false);
        void run_knl_assemble_field_SLM_2d(bool debug_timing = false);
        void run_knl_assemble_field_tweezers(bool debug_timing = false);
        void run_knl_set_amplitude_field_SLM(bool debug_timing = false);
        void run_knl_set_amplitude_field_SLM_2d(bool debug_timing = false);
        void run_knl_set_amplitude_field_tweezers(bool debug_timing = false);
        void run_knl_set_phase_tweezers(bool debug_timing = false);
        void run_knl_get_phase_tweezers(bool debug_timing = false);
        void run_knl_change_phase_conditional(bool debug_timing = false);
        void run_knl_apply_constraints_GSW(bool debug_timing = false);
        void run_knl_apply_constraints_GS_locked(bool debug_timing = false);
        void run_knl_apply_constraints_GSW_locked(bool debug_timing = false);
        void run_knl_reset_field(bool debug_timing=false);
        void run_knl_set_target(bool debug_timing=false);
        void run_knl_get_phase_with_corrections(bool debug_timing=false);
        void run_knl_get_phase_with_corrections_with_roll(int offset_x, int offset_y, bool debug_timing=false);
        void initiate_spot_buffers_with_length(unsigned int length);

    private:
        cl_kernel knl_assemble_field_SLM;
        cl_kernel knl_assemble_field_SLM_2d;
        cl_kernel knl_assemble_field_tweezers;
        cl_kernel knl_set_amplitude_field_SLM;
        cl_kernel knl_set_amplitude_field_SLM_2d;
        cl_kernel knl_set_amplitude_field_tweezers;
        cl_kernel knl_set_phase_tweezers;
        cl_kernel knl_get_phase_tweezers;
        cl_kernel knl_get_phase;
        cl_kernel knl_get_phase_2d;
        cl_kernel knl_calculate_amplitudes;
        cl_kernel knl_update_target;
        cl_kernel knl_apply_constraints_GSW_locked;
        cl_kernel knl_apply_constraints_GSW;
        cl_kernel knl_apply_constraints_GS_locked;
        cl_kernel knl_reset_field;
        cl_kernel knl_set_target;
        cl_kernel knl_get_phase_with_corrections;   
        cl_kernel knl_get_phase_with_corrections_with_roll;    
    };
}