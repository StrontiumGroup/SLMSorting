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

#include "mosquitto.h" 

#include "pmg.h"

PMG::PMG* pmg;

struct mosquitto* MQTT_client;
char payload[2000];

std::vector<PMG::Spot> spots_6x6 = {};
std::vector<PMG::Spot> spots_circle = {};
std::vector<PMG::Spot> spots_4x4 = {};

/*
* void sort()
*
* Sorts into a circle given a loaded bitstring. Returns a black screen
* on the SLM in case there were not enough atoms. 
*
* Parameters
* ----------
* loaded_str : std::string
*   Loaded bitstring, e.g. "101010101"
*/
void sort(std::string loaded_str) {
    std::cout << "Received message: " << loaded_str << std::endl;

    std::string target_str = "1111111111111111";
    pmg->load_end_geometry(spots_circle);
    pmg->sort_arbitrary(
        target_str,
        loaded_str,
        1U, 1U, 0U
    );
    
    if (pmg->sorting_machine->number_loaded  >= pmg->sorting_machine->number_target) {
        std::string path_to_end_img = "../masks/circle.bmp";
        pmg->load_mask_from_file(path_to_end_img);
        pmg->add_corrections_to_mask();
        pmg->display_mask_on_SLM();
    }

    pmg->load_start_geometry(spots_circle);

    pmg->gpu_handler->initiate_spot_buffers_with_length(pmg->number_of_spots);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();
}

/*
* void sort_to_four()
*
* Uses the current spots and configuration to sort into a 4x4 array. Typically 
* used after a sort() call that sorts an arbitrary loaded pattern into a known
* geometry. It assumes all traps are loaded.
*/
void sort_to_four() {
    if (pmg->sorting_machine->number_loaded < pmg->sorting_machine->number_target) return;
    
    pmg->load_end_geometry(spots_4x4);

    std::string target_str = "1111111111111111";
    pmg->sort_arbitrary(
        target_str,
        target_str,
        1U, 1U, 0U
    );

    std::string path_to_end_img = "../masks/4x4.bmp";
    pmg->load_mask_from_file(path_to_end_img);
    pmg->add_corrections_to_mask();
    pmg->display_mask_on_SLM();

    pmg->load_start_geometry(spots_4x4);

    pmg->gpu_handler->initiate_spot_buffers_with_length(pmg->number_of_spots);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();
}

/*
* void sort_to_circle()
*
* Uses the current spots and configuration to sort into a circle. Typically 
* used after a sort() call that sorts an arbitrary loaded pattern into a known
* geometry.
*/
void sort_to_circle() {
    if (pmg->sorting_machine->number_loaded < pmg->sorting_machine->number_target) return;
    
    pmg->load_end_geometry(spots_circle);

    std::string target_str = "1111111111111111";
    pmg->sort_arbitrary(
        target_str,
        target_str,
        1U, 1U, 0U
    );

    std::string path_to_end_img = "../masks/circle.bmp";

    pmg->load_start_geometry(spots_circle);

    pmg->gpu_handler->initiate_spot_buffers_with_length(pmg->number_of_spots);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();
}

/*
* void arm_pmg_for_sort()
*
* Resets the geometry to the 6x6 geometry we used in the start of the example
* and puts that hologram on the SLM. This is to be called after every successful
* sorting run to reset.
*/
void arm_pmg_for_sort() {
    pmg->load_start_geometry(spots_6x6);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();


    std::string path_to_start_img = "../masks/6x6.bmp";
    pmg->set_starting_phase_from_file(path_to_start_img);

    pmg->add_corrections_to_mask();
    pmg->display_mask_on_SLM();
}

/*
* Function that is called when MQTT connects. Here we just subscribe to a channel.
*/
void on_connect(struct mosquitto* MQTT_client, void* userdata, int return_code) {
    std::cout << "Connected to MQTT Server with return code: " << return_code << std::endl;

    // Subscribe to the relevant channel. I have used "fastSLM" as an example.
    mosquitto_subscribe(MQTT_client, NULL, "fastSLM/#", 0);
}

/*
* Function called on MQTT disconnect that tries to reconnect a couple times.
*/
void on_disconnect(struct mosquitto* MQTT_client, void* userdata, int return_code) {
    std::cout << "Disconnected from MQTT Server. Reconnecting..." << std::endl;

    unsigned int number_attempts = 0;
    unsigned int max_number_attempts = 300;
    unsigned int delay_between_attempts = 5;

    while (number_attempts < max_number_attempts) {
        if (mosquitto_reconnect(MQTT_client) == MOSQ_ERR_SUCCESS) {
            std::cout << "Reconnected to MQTT Server!" << std::endl;
            break;
        }
        std::cout << "Failed to reconnect. Sleeping for " << delay_between_attempts << " seconds before retrying..." << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(delay_between_attempts));
        number_attempts += 1;
    }
}

/*
* Function that is called when on the MQTT server a message is published in a channel that
* we are subscribed to. This decides what codeblock to run.
*/
void on_message(struct mosquitto* MQTT_client, void* userdata, const struct mosquitto_message* message) {
    if (strstr((char *) message->topic, "action")) {
        //std::cout << "Received message: " << (char *) message->payload << std::endl;
        if ((char *) message->payload == "reset") {
            std::cout << "Need to reset system" << std::endl;
        }
        if (strstr((char *) message->payload, "arm")) {
            std::cout << "Rearming system" << std::endl;
            arm_pmg_for_sort();
        }
        if (strstr((char *) message->payload, "circle")) {
            std::cout << "Sorting to circle" << std::endl;
            sort_to_circle();
        }
        if (strstr((char *) message->payload, "4x4")) {
            std::cout << "Sort to 4x4" << std::endl;
            sort_to_four();
        }
    }
    if (strstr((char *) message->topic, "loadstr")) {
        std::string loaded_str = (char *) message->payload;
        sort(loaded_str);
    }
}

int main(void)
{
    // -------------------------- INITIALIZE MQTT SERVER --------------------------------

    mosquitto_lib_init();

    // Now we load the MQTT client.
    MQTT_client = mosquitto_new("FastSLM", true, NULL);
    mosquitto_connect_callback_set(MQTT_client, on_connect);
    mosquitto_disconnect_callback_set(MQTT_client, on_disconnect);
    mosquitto_message_callback_set(MQTT_client, on_message);

    // Log in using your username and password.
    mosquitto_username_pw_set(MQTT_client, "USERNAME", "PASSWORD");
    if (mosquitto_connect(MQTT_client, "192.168.0.108", 1883, 60) == MOSQ_ERR_SUCCESS) {
        snprintf(payload, sizeof(payload), "Connected");
        mosquitto_publish(MQTT_client, NULL, "fastSLM/status", strlen(payload), payload, 0, false);
    }
    else {
        std::cout << "Couldn't connect to MQTT server." << std::endl;
        exit(1);
    }

    // --------------------- INITIALIZE PHASEMASKGENERATOR CLASS ------------------------

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

    // Start talking to the SLM.
    pmg->initiate_SLM();

    std::string path_to_lut = "../lib/lut.lut";
    pmg->load_LUT_from_file(path_to_lut);

    std::vector<std::string> paths_to_corrections = {
            "../masks/black.bmp"
            /* E.g.:
            "C:/Users/FastSLM/SLM/Masks/lenses/042.bmp",
            "C:/Users/FastSLM/SLM/Masks/grating/080.bmp",//,
            "C:/Users/FastSLM/SLM/Masks/zernike/sa_new/-2.10.bmp"*/
    };
    pmg->load_corrections(paths_to_corrections);

    // ----------------------------- DEFINE SPOT GEOMETRIES ------------------------------
    std::vector<PMG::Spot> start_spots = pmg->create_spots_square(6, 6, 12);
    spots_6x6 = pmg->create_spots_square(6, 6, 12);

    std::vector<unsigned int> x_values_four = {1006, 1006, 1006, 1006, 1018, 1018, 1018, 1018, 6, 6, 6, 6, 18, 18, 18, 18};
    std::vector<unsigned int> y_values_four = {1006, 1018, 6, 18, 1006, 1018, 6, 18, 1006, 1018, 6, 18, 1006, 1018, 6, 18};
    std::vector<float> phases_four = {-1.4591803058132753, 0.894663751494231, 2.878190727951471, -1.8922992138313017, -0.7918663248779386, 0.5606194849721969, -0.39316566100870853, 1.7140012571902916, 1.1800529163191757, 0.14429277880880412, 2.372908818707769, 0.5393008686094323, 1.7172775903293909, -0.08427972468486197, -1.2034505070965047, -1.6947271630328835};
    for (auto i=0; i<x_values_four.size(); i++) {
        PMG::Spot spot;
        spot.x = x_values_four[i];
        spot.y = y_values_four[i];
        spot.a = 1.;
        spot.psi = phases_four[i];
        spot.used = 1;
        spots_4x4.emplace_back(spot);
    }

    // Circle
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

    // ----------------------------- BUILD PROGRAM AND START LOOP ------------------------------


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

    // Display the start hologram on SLM with corrections.
    pmg->add_corrections_to_mask();
    pmg->display_mask_on_SLM();

    // Now we start the MQTT loop needed in our experiment to let the control sequence talk 
    // to the c++ SLM code.
    mosquitto_loop_forever(MQTT_client, -1, 1);

    system("\n\nPause");

    return 0;
}
