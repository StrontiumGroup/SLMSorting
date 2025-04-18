/*
* sortingmachine.h
*
*
*/

//#pragma once

// Load core libraries
#include "utilities.h"

namespace PMG {
    // ------------------------------------------ PMG :: SortingMachine --------------------------------- //

    class SortingMachine {
        public:
            bool is_target[MAX_NUMBER_SPOTS];
            bool is_loaded[MAX_NUMBER_SPOTS];
            unsigned int target_id[MAX_NUMBER_SPOTS];
            unsigned int loaded_id[MAX_NUMBER_SPOTS];
            int mapping[MAX_NUMBER_SPOTS];
            int cost_array[MAX_NUMBER_SPOTS * MAX_NUMBER_SPOTS];

            unsigned int number_of_atoms;
            unsigned int number_target;
            unsigned int number_loaded;

            unsigned int current_step;
            unsigned int number_of_steps;
            unsigned int number_turnoff_frames;
            
            unsigned int sorting_method;
            unsigned int pathing_method;


            SortingMachine() {
                current_step = 0;
                number_of_steps = 1;

                number_of_atoms = 0;
                number_target = 0;
                number_loaded = 0;

                for (unsigned int i=0; i < MAX_NUMBER_SPOTS; i++) {
                    is_target[i] = false;
                    is_loaded[i] = false;
                    mapping[i] = -1;

                    for (unsigned int j=0; j < MAX_NUMBER_SPOTS; j++) {
                        cost_array[i * MAX_NUMBER_SPOTS + j] = 0;
                    }
                }
            };

            int sort(
                std::vector<unsigned int> &x_values,
                std::vector<unsigned int> &y_values,
                unsigned int dimension
            );
            int sort(
                std::vector<unsigned int> &start_x_values, 
                std::vector<unsigned int> &start_y_values, 
                std::vector<unsigned int> &end_x_values, 
                std::vector<unsigned int> &end_y_values, 
                unsigned int dimension
            );
            void reset();
            int create_cost_array(
                std::vector<unsigned int> &x_values_start,
                std::vector<unsigned int> &y_values_start,
                std::vector<unsigned int> &x_values_end,
                std::vector<unsigned int> &y_values_end,
                unsigned int dimension
            );
            int create_cost_array_square(
                std::vector<unsigned int> &x_values,
                std::vector<unsigned int> &y_values,
                unsigned int dimension
            );
            int create_cost_array_square(
                std::vector<int> &cost_array_partition,
                std::vector<unsigned int> &loaded_partition,
                std::vector<unsigned int> &target_partition,
                std::vector<unsigned int> &x_values,
                std::vector<unsigned int> &y_values,
                unsigned int dimension
            );
            int solve_maximal_cost_assignment();
            int solve_maximal_cost_assignment(
                std::vector<int> &cost_array_partition,
                std::vector<unsigned int> &loaded_partition,
                std::vector<unsigned int> &target_partition
            );
            int solve_partitioned_compression(
                std::vector<unsigned int> &x_values,
                std::vector<unsigned int> &y_values,
                unsigned int dimension
            );
            int solve_compression(
                std::vector<unsigned int> &x_values,
                std::vector<unsigned int> &y_values,
                unsigned int dimension,
                unsigned int layer_scaling = 8U,
                unsigned int first_layer = 4U
            );
            int solve_compression(
                std::vector<unsigned int> &x_values,
                std::vector<unsigned int> &y_values,
                std::vector<unsigned int> &loaded_partition,
                std::vector<unsigned int> &target_partition,
                unsigned int dimension,
                unsigned int layer_scaling = 8U,
                unsigned int first_layer = 4U
            );
            void trim_moves();
            void parse_target_str(
                std::string_view set_target_str);
            void parse_loaded_str(
                std::string_view set_loaded_str);
            bool get_next_positions(
                unsigned int iteration,
                std::vector<unsigned int> &x_values,
                std::vector<unsigned int> &y_values,
                std::vector<float> &amplitude_values,
                std::vector<float> &phase_values,
                std::vector<unsigned int> &should_be_used,
                std::vector<unsigned int> &start_x_values,
                std::vector<unsigned int> &start_y_values,
                std::vector<float> &start_amplitudes,
                std::vector<float> &start_phases,
                std::vector<float> &end_phases,
                unsigned int dimension 
            );
            bool get_next_positions(
                unsigned int iteration,
                std::vector<unsigned int> &x_values,
                std::vector<unsigned int> &y_values,
                std::vector<float> &amplitude_values,
                std::vector<float> &phase_values,
                std::vector<unsigned int> &should_be_used,
                std::vector<unsigned int> &start_x_values,
                std::vector<unsigned int> &start_y_values,
                std::vector<float> &start_amplitudes,
                std::vector<float> &start_phases,
                std::vector<unsigned int> &end_x_values,
                std::vector<unsigned int> &end_y_values,
                std::vector<float> &end_amplitudes,
                std::vector<float> &end_phases,
                unsigned int dimension 
            );
            void calculate_maximal_movement(
                std::vector<unsigned int> &x_values,
                std::vector<unsigned int> &y_values,
                unsigned int dimension
            );
            void calculate_maximal_movement(
                std::vector<unsigned int> &start_x_values,
                std::vector<unsigned int> &start_y_values,
                std::vector<unsigned int> &end_x_values,
                std::vector<unsigned int> &end_y_values,
                unsigned int dimension
            );
            std::string create_target_string_center(
                unsigned int n,
                unsigned int m
            );
            std::string create_loaded_string_random(
                unsigned int n,
                float loading_probability = 0.55
            );
            void sort_ids_from_center(
                std::vector<unsigned int> &x_values, 
                std::vector<unsigned int> &y_values, 
                unsigned int dimension
            );

    };
}