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
std::vector<PMG::Spot> spots_circle = {};
std::vector<unsigned int> number_of_moves = {};

void sort(std::string loaded_str) {
    std::cout << "Received message: " << loaded_str << std::endl;

    std::string target_str = "1111111111111111";
    pmg->load_end_geometry(spots_circle);
    unsigned int moves = pmg->sort_arbitrary_return_noslm(
        target_str,
        loaded_str,
        1U, 1U, 0U
    );
    
    if (pmg->sorting_machine->number_loaded  >= pmg->sorting_machine->number_target) {
        number_of_moves.emplace_back(moves);
    }

    pmg->load_start_geometry(spots_circle);
    pmg->gpu_handler->initiate_spot_buffers_with_length(pmg->number_of_spots);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();
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

    // For the circle we should manually set it, because there is no function
    // that does this in the PMG class. Luckily, it is not hard to do.
    std::vector<int> x_values_circle = {31, 28, 22, 12, 0, -12, -22, -28, -31, -28, -22, -12, 0, 12, 22, 28};
    std::vector<int> y_values_circle = {0, 12, 22, 28, 31, 28, 22, 12, 0, -12, -22, -28, -31, -28, -22, -12};
    std::vector<float> phases_circle =  {-2.961337224963786, 0.9847102148183918, -1.1476406090187652, 0.9091234106038645, -2.863836290929989, -0.6017605743886563, -3.096439295742857, -0.719928358963168, -2.9588488306614287, 0.9530670146420372, -1.1425199370999757, 0.9434168760999236, -2.8816542186364007, -0.6105938832000518, 0.06412604840669382, -0.7527727036780225};
    for (auto i=0; i<x_values_circle.size(); i++) {
        PMG::Spot spot;
        spot.x = x_values_circle[i] % pmg->width;
        spot.y = y_values_circle[i] % pmg->height;
        spot.a = 1.;
        spot.psi = phases_circle[i];
        spot.used = 1;
        spots_circle.emplace_back(spot);
    }

    pmg->load_start_geometry(spots_6x6);
    pmg->copy_spots_to_buffers();
    pmg->load_end_geometry(spots_circle);

    std::cout << "Created vectors for the geometries" << std::endl;

    // Building kernels is necessary to enable GPU usage.
    pmg->build_kernels();
    pmg->set_target_buffer();
    pmg->set_illumination();

    // Next we load the start and end holograms
    std::string path_to_start_img = "../masks/6x6.bmp";
    std::string path_to_end_img = "../masks/circle.bmp";

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
