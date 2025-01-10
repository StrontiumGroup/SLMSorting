/*
* pmg.h
*/

#pragma once

#include "utilities.h"
#include "sortingmachine.h"
#include "gpuhandler.h"

#include "Blink_C_wrapper.h"
#include "Blink_SDK.h"

//#include <Windows.h>

namespace PMG {
    // -----------------------------------------  PMG :: PMG -------------------------------------------- //

    class PMG
        {
        public:
            cl_device_id *device_ptr;
            cl_context *context_ptr;
            cl_command_queue *queue_ptr;
            GPU_Handler *gpu_handler;

            SortingMachine *sorting_machine;

            unsigned int width;
            unsigned int height;
            unsigned int number_of_spots;

            std::vector<unsigned int> x_values;
            std::vector<unsigned int> y_values;
            std::vector<float> amplitude_values;
            std::vector<float> psi_values;
            std::vector<unsigned int> should_be_used;

            std::vector<unsigned int> start_x_values;
            std::vector<unsigned int> start_y_values;
            std::vector<unsigned int> end_x_values;
            std::vector<unsigned int> end_y_values;
            std::vector<float> start_amplitude_values;
            std::vector<float> end_amplitude_values;
            std::vector<float> start_phase_values;
            std::vector<float> end_phase_values;

            bool SLM_connected;
            uint8_t *mask_8bit;
            uint8_t *correction_mask_8bit;
            uint8_t *mask_to_display_8bit;

            PMG(
                cl_device_id *set_device_ptr,
                cl_context *set_context_ptr,
                cl_command_queue *set_queue_ptr,
                unsigned int set_width = 1024,
                unsigned int set_height = 1024,
                const char *kernel_str="./kernel.cl"
            ) {
                device_ptr = set_device_ptr;
                context_ptr = set_context_ptr;
                queue_ptr = set_queue_ptr;

                width = set_width;
                height = set_height;
                mask_8bit = new uint8_t[width * height];
                correction_mask_8bit = new uint8_t[width * height];
                mask_to_display_8bit = new uint8_t[width * height];
                number_of_spots = 0;

                gpu_handler = new GPU_Handler(
                    set_device_ptr,
                    set_context_ptr,
                    set_queue_ptr,
                    &width,
                    &height,
                    kernel_str
                );

                sorting_machine = new SortingMachine();

                SLM_connected = false;
            };

            ~PMG()
            {
                delete mask_8bit;
                delete correction_mask_8bit;
                delete gpu_handler;
                delete sorting_machine;
                SLM_power(false);
                Delete_SDK();
            }

            void fftshift_mask();
            void create_square_pattern(
                unsigned int num_x,
                unsigned int num_y,
                unsigned int spot_distance);
            std::vector<Spot> create_spots_square(
                unsigned int num_x,
                unsigned int num_y,
                unsigned int spot_distance);
            void shift_pattern(
                int dx,
                int dy);
            void shift_pattern(
                int dx,
                int dy,
                float dpsi);
            void spread_pattern(
                bool direction 
            );
            void spread_pattern(
                bool direction,
                unsigned int step,
                unsigned int number_of_steps
            );
            void set_illumination(
                unsigned int option = 0);
            void set_target_buffer();
            bool replace_spot(
                unsigned int index,
                Spot &spot);
            bool replace_spot(
                unsigned int index,
                unsigned int x,
                unsigned int y,
                float amp,
                float psi,
                unsigned int used);
            void load_start_geometry(
                std::vector<Spot> spots);
            void load_start_geometry(
                std::vector<unsigned int> target_x_values,
                std::vector<unsigned int> end_x_values);
            void load_end_geometry(
                std::vector<Spot> spots);
            void load_end_geometry(
                std::vector<unsigned int> target_x_values,
                std::vector<unsigned int> end_x_values);
            bool initiate_SLM();
            bool load_LUT_from_file(
                std::string path_to_file);
            bool load_corrections(
                std::string path_to_file);
            bool load_corrections(
                std::vector<std::string> paths_to_files);
            void add_corrections_to_mask();
            void remove_corrections_from_mask();
            bool display_mask_on_SLM();
            bool display_mask_on_SLM_immediately();
            bool display_mask_on_SLM_ext_trigger();
            bool display_mask_on_SLM_ext_trigger_immediately();
            bool copy_spots_to_buffers(bool block=false);
            bool copy_mask_to_buffer();
            bool copy_corrections_to_buffer();
            bool copy_input_light_to_buffer();
            bool set_starting_phase_from_file(
                std::string path_to_file);
            bool set_starting_phase_and_amps_from_file(
                std::string path_to_file);
            bool set_ending_phase_from_file(
                std::string path_to_file);
            bool load_random_phase();
            bool load_mask_from_file(
                std::string path_to_file,
                bool do_shift=false);
            bool get_mask_from_buffer();
            bool get_mask_from_buffer(uint8_t *save_mask);
            bool save_mask_to_file(
                std::string path_to_file);
            void build_kernels();
            void define_start();
            void reset();
            void reset_field();
            void roll_mask_x(int dx);
            void roll_mask_y(int dy);
            void run_GSW();
            void run_GSW(int token);
            void run_GSW_with_roll(int offset_x, int offset_y, int token);
            void run_GSW_loop();
            void calculate_sorting_shortcut();
            void calculate_sorting_arbitrary();
            TestResult test(
                std::string_view loaded_str,
                std::string_view target_str,
                unsigned int sorting_method = 0,
                unsigned int calculation_method = 0,
                unsigned int pathing_method = 0
            );
            TestResult test_shift_w_slm();
            void test_GSW(
                unsigned int test_num_iterations = 5, 
                bool test_debug_timings = false, 
                bool test_total_timing = true);
            void test_sorting(
                std::string_view loaded_str,
                std::string_view target_str,
                unsigned int test_num_iterations = 1,
                bool test_debug_timings = false,
                bool test_total_timing = true,
                bool save_output = false);
            void sort(
                std::string_view target_str,
                std::string_view loaded_str,
                unsigned int sorting_method = 0,
                unsigned int calculation_method = 0,
                unsigned int pathing_method = 0
            );
            void sort_arbitrary(
                std::string_view target_str,
                std::string_view loaded_str,
                unsigned int sorting_method = 0,
                unsigned int calculation_method = 0,
                unsigned int pathing_method = 0
            );
            unsigned int sort_arbitrary_return(
                std::string_view target_str,
                std::string_view loaded_str,
                unsigned int sorting_method = 0,
                unsigned int calculation_method = 0,
                unsigned int pathing_method = 0
            );
            void test_sort(
                std::string_view target_str,
                std::string_view loaded_str
            );
            void test_sort(
                std::string_view target_str,
                std::string_view loaded_str,
                int token
            );
    };
}