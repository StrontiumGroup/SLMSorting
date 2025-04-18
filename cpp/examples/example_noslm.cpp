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

std::vector<PMG::Spot> spots_6x6 = {};
std::vector<unsigned int> number_of_moves = {};

void sort(std::string loaded_str) {
    std::cout << "Received message: " << loaded_str << std::endl;

    std::string target_str = "101010000000101010000000101010000000";
    unsigned int moves = pmg->sort_arbitrary_return_noslm(
        target_str,
        loaded_str,
        1U, 1U, 0U
    );
    
    if (pmg->sorting_machine->number_loaded  >= pmg->sorting_machine->number_target) {
        number_of_moves.emplace_back(moves);
    }
}

int main(void)
{
    cl_device_id device;
    cl_context context;
    cl_command_queue queue;

    // Standard function to load GPU 0. Could be changed for other PC.
    setup_opencl(0, &device, &context, &queue);

    char device_name[100];
    clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_name), &device_name, NULL);
    std::cout<< "Using device: "<< device_name <<"\n";
    
    // Make PMG instance.
    pmg = new PMG::PMG (&device, &context, &queue, 1024, 1024, "../core/kernel.cl");
    std::cout << "loaded PMG" << std::endl;

    std::vector<std::string> paths_to_corrections = {
            "../masks/black.bmp"
    };
    pmg->load_corrections(paths_to_corrections);

    // We need to specify where the spots are going to be. For this we use
    // A vector of PMG::Spot structs that contain position, amplitude and phase.
    spots_6x6 = pmg->create_spots_square(6, 6, 12);
    pmg->load_start_geometry(spots_6x6);
    pmg->copy_spots_to_buffers();
    pmg->load_end_geometry(spots_6x6);
    
    std::cout << "Created vectors for the geometries" << std::endl;

    // Building kernels is necessary to enable GPU usage.
    pmg->build_kernels();
    pmg->set_target_buffer();
    pmg->set_illumination();

    // Next we load the start and end holograms
    std::string path_to_start_img = "../masks/6x6.bmp";
    std::string path_to_end_img = "../masks/6x6.bmp";

    // First load the end image so that we can be faster.
    pmg->set_ending_phase_from_file(path_to_end_img);
    pmg->set_starting_phase_from_file(path_to_start_img);

    // For this experiment, we do not need more. The patterns are loaded and we
    // now only need to test the sorting once.
    std::string loaded_str = pmg->sorting_machine->create_loaded_string_random(6 * 6, 0.5);
    sort(loaded_str);

    system("\n\nPause");

    return 0;
}
