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

std::vector<PMG::Spot> spots_diamond = {};
std::vector<PMG::Spot> spots_triangular = {};
std::vector<PMG::Spot> spots_kagome = {};
std::vector<PMG::Spot> spots_circle = {};
std::vector<PMG::Spot> spots_four = {};
std::vector<PMG::Spot> spots_six = {};
std::vector<PMG::Spot> spots_blockade = {};
std::vector<PMG::Spot> spots_wide = {};
std::vector<unsigned int> number_of_moves = {};

void sort(std::string loaded_str) {
    std::cout << "Received message: " << loaded_str << std::endl;

    std::string target_str = "111111111";
    pmg->load_end_geometry(spots_circle);
    unsigned int moves = pmg->sort_arbitrary_return(
        target_str,
        loaded_str,
        1U, 1U, 0U
    );
    
    if (pmg->sorting_machine->number_loaded  >= pmg->sorting_machine->number_target) {
        std::string path_to_end_img = "C:/Users/FastSLM/SLM/Masks/circle.bmp";
        pmg->load_mask_from_file(path_to_end_img);
        pmg->add_corrections_to_mask();
        pmg->display_mask_on_SLM();
        number_of_moves.emplace_back(moves);
    }

    pmg->load_start_geometry(spots_circle);

    pmg->gpu_handler->initiate_spot_buffers_with_length(pmg->number_of_spots);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();
}

void sort_to_four() {
    if (pmg->sorting_machine->number_loaded < pmg->sorting_machine->number_target) return;
    
    pmg->load_end_geometry(spots_four);

    std::string target_str = "1111111111111111";
    pmg->sort_arbitrary(
        target_str,
        target_str,
        1U, 1U, 0U
    );

    std::string path_to_end_img = "C:/Users/FastSLM/SLM/Masks/4x4_2.bmp";
    pmg->load_mask_from_file(path_to_end_img);
    pmg->add_corrections_to_mask();
    pmg->display_mask_on_SLM();

    pmg->load_start_geometry(spots_four);

    pmg->gpu_handler->initiate_spot_buffers_with_length(pmg->number_of_spots);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();
}

void sort_to_circle() {
    if (pmg->sorting_machine->number_loaded < pmg->sorting_machine->number_target) return;
    
    pmg->load_end_geometry(spots_circle);

    std::string target_str = "1111111111111111";
    pmg->sort_arbitrary(
        target_str,
        target_str,
        1U, 1U, 0U
    );

    std::string path_to_end_img = "C:/Users/FastSLM/SLM/Masks/sort_all/circle2.bmp";

    pmg->load_start_geometry(spots_circle);

    pmg->gpu_handler->initiate_spot_buffers_with_length(pmg->number_of_spots);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();
}

void sort_to_diamond() {
    if (pmg->sorting_machine->number_loaded < pmg->sorting_machine->number_target)  return;

    pmg->load_end_geometry(spots_diamond);
    std::string target_str = "1111111111111111";
    pmg->sort_arbitrary(
        target_str,
        target_str,
        1U, 1U, 0U
    );
    std::string path_to_end_img = "C:/Users/FastSLM/SLM/Masks/sort_all/diamond_1.bmp";

    pmg->load_start_geometry(spots_diamond);
    pmg->gpu_handler->initiate_spot_buffers_with_length(pmg->number_of_spots);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();
}


void sort_to_kagome() {
    if (pmg->sorting_machine->number_loaded < pmg->sorting_machine->number_target)  return;

    pmg->load_end_geometry(spots_kagome);
    std::string target_str = "1111111111111111";
    pmg->sort_arbitrary(
        target_str,
        target_str,
        1U, 1U, 0U
    );

    std::string path_to_end_img = "C:/Users/FastSLM/SLM/Masks/sort_all/kagome_1.bmp";
    pmg->load_mask_from_file(path_to_end_img);
    pmg->add_corrections_to_mask();
    pmg->display_mask_on_SLM();

    pmg->load_start_geometry(spots_kagome);

    pmg->gpu_handler->initiate_spot_buffers_with_length(pmg->number_of_spots);
    pmg->reset();
    pmg->build_kernels();
}

void sort_to_triangular() {
    if (pmg->sorting_machine->number_loaded < pmg->sorting_machine->number_target)  return;

    pmg->load_end_geometry(spots_triangular);
    std::string target_str = "1111111111111111";
    pmg->sort_arbitrary(
        target_str,
        target_str,
        1U, 1U, 0U
    );
    std::string path_to_end_img = "C:/Users/FastSLM/SLM/Masks/sort_all/triangular_1.bmp";

    pmg->load_start_geometry(spots_triangular);

    pmg->gpu_handler->initiate_spot_buffers_with_length(pmg->number_of_spots);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();
}


void arm_pmg_for_sort() {
    pmg->load_start_geometry(spots_six);
    pmg->copy_spots_to_buffers();
    pmg->build_kernels();


    std::string path_to_start_img = "C:/Users/FastSLM/SLM/Masks/6x6_14.bmp";
    pmg->set_starting_phase_from_file(path_to_start_img);

    pmg->add_corrections_to_mask();
    pmg->display_mask_on_SLM();
}

void on_connect(struct mosquitto* MQTT_client, void* userdata, int return_code) {
    std::cout << "Connected to MQTT Server with return code: " << return_code << std::endl;

    // Subscribe to the relevant channel.
    mosquitto_subscribe(MQTT_client, NULL, "fastSLM/#", 0);
}

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
        if (strstr((char *) message->payload, "diamond")) {
            std::cout << "Sorting to diamond" << std::endl;
            sort_to_diamond();
        }
        if (strstr((char *) message->payload, "triangular")) {
            std::cout << "Sort to triangular" << std::endl;
            sort_to_triangular();
        }
        if (strstr((char *) message->payload, "kagome")) {
            std::cout << "Sort to kagome" << std::endl;
            sort_to_kagome();
        }
        if (strstr((char *) message->payload, "4x4")) {
            std::cout << "Sort to 4x4" << std::endl;
            sort_to_four();
        }
        if (strstr((char *) message->payload, "shift")) {
            std::string command ((char *) message->payload);
            std::cout << "Shifting" << std::endl;
            shift_and_back(command);
        }
        if (strstr((char *) message->payload, "compress")) {
            std::string command ((char *) message->payload);
            std::cout << "Compressing" << std::endl;
            compress(command);
        }
        if (strstr((char *) message->payload, "save")) {
            auto point = std::chrono::system_clock::now();
            auto timestamp = std::chrono::duration_cast<std::chrono::seconds>(point.time_since_epoch()).count();

            char log_file[100];
            snprintf(log_file, 200, "./sorting_moves_%d.txt", timestamp);
            std::ofstream file(log_file, std::ios::out );
            file << "Size\tShifting\tCalculation\tTransfer\tDisplay\tTotal\n";
            char buf[50];
            for (auto moves : number_of_moves) {
                snprintf(
                    buf, 50, "%d\n", moves
                );
                file << buf;
            }
        } 
    }
    if (strstr((char *) message->topic, "loadstr")) {
        std::string loaded_str = (char *) message->payload;
        sort(loaded_str);
    }
}

int main(void)
{
    mosquitto_lib_init();

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

    // Now we load the MQTT client with the PMG as userdata pointer.
    MQTT_client = mosquitto_new("FastSLM", true, NULL);

    mosquitto_connect_callback_set(MQTT_client, on_connect);
    mosquitto_disconnect_callback_set(MQTT_client, on_disconnect);
    mosquitto_message_callback_set(MQTT_client, on_message);

    mosquitto_username_pw_set(MQTT_client, "srmicmqtt", "90Sr@lab");
    if (mosquitto_connect(MQTT_client, "192.168.0.108", 1883, 60) == MOSQ_ERR_SUCCESS) {
        snprintf(payload, sizeof(payload), "Connected");
        mosquitto_publish(MQTT_client, NULL, "fastSLM/status", strlen(payload), payload, 0, false);
    }
    else {
        std::cout << "Couldn't connect to MQTT server." << std::endl;
        exit(1);
    }

    // Start talking to the SLM.
    pmg->initiate_SLM();

    std::string path_to_lut = "C:/Users/FastSLM/SLM/LUT Files/slm6673_at813_at33c_04.lut";
    pmg->load_LUT_from_file(path_to_lut);

    std::vector<std::string> paths_to_corrections = {
            "C:/Users/FastSLM/SLM/Masks/lenses/042.bmp",
            "C:/Users/FastSLM/SLM/Masks/grating/080.bmp",//,
            "C:/Users/FastSLM/SLM/Masks/zernike/sa_new/-2.10.bmp"
    };
    pmg->load_corrections(paths_to_corrections);

    std::vector<PMG::Spot> start_spots = pmg->create_spots_square(6, 6, 12);
    spots_six = pmg->create_spots_square(6, 6, 12);

    std::vector<unsigned int> x_values_window = {982, 982, 982, 994, 994, 994, 994, 994, 994, 994, 1006, 1006, 1006, 1006, 1018, 1018, 1018, 1018, 6, 6, 6, 6, 18, 18, 18, 18, 30, 30, 30, 30, 30, 30, 30, 42, 42, 42};
    std::vector<unsigned int> y_values_window = {1012, 0, 12, 988, 1000, 1012, 0, 12, 24, 36, 988, 1000, 24, 36, 988, 1000, 24, 36, 988, 1000, 24, 36, 988, 1000, 24, 36, 988, 1000, 1012, 0, 12, 24, 36, 1012, 0, 12};
    std::vector<float> phases_window = {
        0.7357778414625793, 3.0170158530286817, -0.7483525901799918, -3.1242329928512795, 2.632983493102855, -3.088754460740113, 0.357422040979916, 0.9790019783166862, -2.986686108872527, -1.7603065468120505, 0.9723184438952222, -3.0033838175510454, 1.8020232393893214, -1.9025848949444, -1.9183593028716124, -3.0159590585064278, -2.837663622963073, -2.894300355174755, 1.662514936738789, 0.7640451212062975, -1.580133478534146, -0.6277757251242932, -0.49898675305684304, -0.42041568687498787, 2.9438224425204345, 1.6642926490291512, -1.408100495032946, -1.421879083972996, 2.8729350268654397, 0.45129545752887373, -2.0794923257482822, 2.49562646729544, -2.1340310249169967, -2.6790977889796737, 1.8249942020591516, -1.1427135190802007
    };
    for (auto i=0; i<y_values_window.size(); i++) {
        PMG::Spot spot;
        spot.x = x_values_window[i];
        spot.y = y_values_window[i];
        spot.a = 1.;
        spot.psi = phases_window[i];
        spot.used = 1;
        spots_window.emplace_back(spot);
    }

    std::vector<unsigned int> x_values_four = {1006, 1006, 1006, 1006, 1018, 1018, 1018, 1018, 6, 6, 6, 6, 18, 18, 18, 18};
    std::vector<unsigned int> y_values_four = {1006, 1018, 6, 18, 1006, 1018, 6, 18, 1006, 1018, 6, 18, 1006, 1018, 6, 18};
    std::vector<float> phases_four = {-1.4591803058132753, 0.894663751494231, 2.878190727951471, -1.8922992138313017, -0.7918663248779386, 0.5606194849721969, -0.39316566100870853, 1.7140012571902916, 1.1800529163191757, 0.14429277880880412, 2.372908818707769, 0.5393008686094323, 1.7172775903293909, -0.08427972468486197, -1.2034505070965047, -1.6947271630328835};
    //std::vector<float> phases_four = {2.176658438814603, -2.2398677524742463, 1.135137148767967, 3.0344888804340457, 3.121592384286797, 1.5872020891491856, 2.4751810657840045, 2.286776996222293, -0.7267836569136586, 1.5693482207476752, 2.7816613081855377, 0.12178042687605288, 0.3161492593291974, 0.44509287008202947, -0.8023115266683168, -0.600168203238994};
    
    for (auto i=0; i<x_values_four.size(); i++) {
        PMG::Spot spot;
        spot.x = x_values_four[i];
        spot.y = y_values_four[i];
        spot.a = 1.;
        spot.psi = phases_four[i];
        spot.used = 1;
        spots_four.emplace_back(spot);
    }

    // Diamond
    std::vector<int> x_values_diamond = {0, -8, 8, -17, 17, -25, 25, -34, 34, 0, -8, 8, -17, 17, -25, 25};
    std::vector<int> y_values_diamond = {-34, -25, -25, -17, -17, -8, -8, 0, 0, 34, 25, 25, 17, 17, 8, 8};
    std::vector<float> phases_diamond = {-2.30764955, 1.24303689, -0.32921938, 0.3160722, -1.24558195, -0.0622023, -1.62470936, 2.93622372, -0.18670828, 1.91556674, -1.63168532, -0.06181699, -0.71353864, 0.86610658, -0.33581903, 1.24296614};
    for (auto i=0; i<x_values_diamond.size(); i++) {
        PMG::Spot spot;
        spot.x = x_values_diamond[i] % pmg->width;
        spot.y = y_values_diamond[i] % pmg->height;
        spot.a = 1.;
        spot.psi = phases_diamond[i];
        spot.used = 1;
        spots_diamond.emplace_back(spot);
    }

    // Triangle
    std::vector<int> x_values_triangular = {-18, -6, 6, 18, -12, 0, 12, 24, -18, -6, 6, 18, -12, 0, 12, 24};
    std::vector<int> y_values_triangular = {-16, -16, -16, -16, -5, -5, -5, -5, 5, 5, 5, 5, 16, 16, 16, 16};
    std::vector<float> phases_triangular = {2.10289723, -2.01884241, -0.55973281, 1.60230062, 2.43762676, 1.4578744, -0.22410987, -1.20438202, -0.48748387, -0.50706962, 0.17629372, 0.15734097, -0.15251007, 2.96998235, 0.51135102, -2.6495236};
    phases_triangular = {2.1196313513959097, -1.0220372413919296, -0.030350224970195, 1.6105390745369512, -0.7483817120367579, -1.702015876916863, 2.9277827900979343, 1.9440422472837797, 2.658290261386526, 2.6250719576084887, -3.0072829929131526, -2.9915656792185117, -0.18151800452146613, 2.940069391915269, 0.5221778284605118, -2.6778158515117876};
    for (auto i=0; i<x_values_triangular.size(); i++) {
        PMG::Spot spot;
        spot.x = x_values_triangular[i] % pmg->width;
        spot.y = y_values_triangular[i] % pmg->height;
        spot.a = 1.;
        spot.psi = phases_triangular[i];
        spot.used = 1;
        spots_triangular.emplace_back(spot);
    }

    // Kagome 
    std::vector<int> x_values_kagome = {-18, -6, 6, 18, -24, 0, 24, -18, -6, 6, 18, -12, 12, -6, 6, 0};
    std::vector<int> y_values_kagome = {-26, -26, -26, -26, -16, -16, -16, -5, -5, -5, -5, 5, 5, 16, 16, 26};
    std::vector<float> phases_kagome = {-0.36679901476931825, 1.4839674420388818, 0.9118127420180709, 0.21900472843088578, -1.2961027606230195, 1.430525753375872, -0.5366699056623159, 1.3629518376382599, -2.5966094393470254, 0.6993477708544704, -1.4941008377471865, -0.9734516370252954, 1.980646908494023, -1.5295241982390742, -0.7691074833127176, 2.0979570916020336};
    for (auto i=0; i<x_values_kagome.size(); i++) {
        PMG::Spot spot;
        spot.x = x_values_kagome[i] % pmg->width;
        spot.y = y_values_kagome[i] % pmg->height;
        spot.a = 1.;
        spot.psi = phases_kagome[i];
        spot.used = 1;
        spots_kagome.emplace_back(spot);
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

    pmg->load_start_geometry(spots_six);
    pmg->copy_spots_to_buffers();
    pmg->load_end_geometry(spots_circle);
    
    std::cout << "Created vectors for the geometries" << std::endl;

    // Building kernels is necessary to enable GPU usage.
    pmg->build_kernels();
    pmg->set_target_buffer();
    pmg->set_illumination();

    // Next we load the start and end holograms
    std::string path_to_start_img = "C:/Users/FastSLM/SLM/Masks/6x6_14.bmp";
    std::string path_to_end_img = "C:/Users/FastSLM/SLM/Masks/3x3_pairs2.bmp";

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
