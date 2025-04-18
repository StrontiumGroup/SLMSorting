#define CL_TARGET_OPENCL_VERSION 300
#define CL_HPP_TARGET_OPENCL_VERSION 300
 
#include <iostream>
#include <vector>
#include <string.h>
#include <sstream>
#include <chrono>
#include <thread>
#include <format>
#include "graybmp.h"
#include "CL/cl.h"
#include "vkFFT.h"

#include "pmg.h"

PMG::PMG* pmg;

std::vector<PMG::Spot> pattern = {};

int main() {
    cl_device_id device;
    cl_context context;
    cl_command_queue queue;

    // Parameters for using the Sorting Machine
    int n = 70;
    int pitch = 12;
    float filling_factor = 0.5;
    int image_size = 1024;

    // Standard function to load GPU 0. Could be changed for other PC.
    setup_opencl(0, &device, &context, &queue);

    char device_name[100];
    clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_name), &device_name, NULL);
    std::cout<< "Using device: "<< device_name <<"\n";
    
    // Make PMG instance.
    pmg = new PMG::PMG (&device, &context, &queue, image_size, image_size, "../core/kernel.cl");
    std::cout << "loaded PMG" << std::endl;

    // We need to specify where the spots are going to be. For this we use
    // A vector of PMG::Spot structs that contain position, amplitude and phase.
    pattern = pmg->create_spots_square(n, n, pitch);
    pmg->load_start_geometry(pattern);
    pmg->copy_spots_to_buffers();
    pmg->load_end_geometry(pattern);
    
    std::cout << "Created vectors for the geometries" << std::endl;

    // Building kernels is necessary to enable GPU usage.
    pmg->build_kernels();
    pmg->set_target_buffer();
    pmg->set_illumination();

    // Next we load the start and end holograms
    std::string path_to_start_img = "../masks/70x70i.bmp";
    std::string path_to_end_img = "../masks/70x70i.bmp";

    // First load the end image so that we can be faster.
    pmg->set_ending_phase_from_file(path_to_end_img);
    pmg->set_starting_phase_from_file(path_to_start_img);

    // Display the start hologram on SLM with corrections.
    pmg->add_corrections_to_mask();

    // Loop through the next stages to see multiple points.
    for (int i = 1; i <= 10; i++) {
        // Reset the system each run to start with start coordinates etc..
        pmg->reset();

        // Create random loading string and a target (could be out of loop)
        std::string loaded_str = pmg->sorting_machine->create_loaded_string_random(n * n, filling_factor);
        std::string target_str = pmg->sorting_machine->create_target_string_center(70, 45);

        auto start = std::chrono::high_resolution_clock::now();
        PMG::TestResult test = pmg->test(
            target_str,
            loaded_str,
            0U, 1U, 0U
        );
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end - start;
        double elapsed_time = elapsed.count();
        int number_of_steps = pmg->sorting_machine->number_of_steps;

        std::cout << "Elapsed time: " << elapsed_time << " seconds\n";

        // Open file in append mode
        std::ofstream file("timing_results.txt", std::ios::app);
        if (file.is_open()) {
            file << test.duration_total << "," << test.duration_sorting << "," << test.duration_calculation << "," << test.number_of_patterns << "\n";
            file.close();
            std::cout << "Elapsed time written to timing_results.txt\n";
        } else {
            std::cerr << "Error: Unable to open file for writing!\n";
        }
    }
    return 0;

    

    std::cout << "All iterations completed.\n";
    
    system("\n\nPause");
    return 0;
}