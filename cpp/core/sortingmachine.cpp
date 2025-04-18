/*
* sortingmachine.cpp
*/

#include "sortingmachine.h"

namespace PMG {
    /*
     * SortingMachine::reset()
     *
     * Resets the mapping.
     */
    void SortingMachine::reset() {
        for (auto i=0; i<MAX_NUMBER_SPOTS; ++i) {
            mapping[i] = -1;
            is_loaded[i] = false;
        }
    }

    /*
     * SortingMachine::sort()
     *
     * Solves the sorting problem using preferred method.
     * 
     * Returns an integer based on the succes off the sorting function.
     * Per method the meaning of each return code can vary. Check the
     * documentation of each function for the meaning per method.
     */
    int SortingMachine::sort(
        std::vector<unsigned int> &x_values, 
        std::vector<unsigned int> &y_values, 
        unsigned int dimension
    ) {
        int return_code;
        switch (sorting_method) {
            case SORTING_METHOD_AUTO:  // Auto-detect
                if (number_target > 300U) {
                    sort_ids_from_center(x_values, y_values, dimension);
                    sorting_method = SORTING_METHOD_PARTITIONED_COMPRESSION;
                    return_code = solve_partitioned_compression(x_values, y_values, dimension);
                    if (return_code == SORTING_MACHINE_PARTITION_FAILED) {
                        sorting_method = SORTING_METHOD_COMPRESSION;
                        return_code = solve_compression(x_values, y_values, dimension);
                    }
                }
                else {
                    sorting_method = SORTING_METHOD_HUNGARIAN;
                    return_code = create_cost_array_square(x_values, y_values, dimension);
                    if (return_code == SORTING_MACHINE_SUCCESS) {
                        return_code = solve_maximal_cost_assignment();
                    }
                }
                break;
            case SORTING_METHOD_HUNGARIAN:  // Hungarian algorithm
                return_code = create_cost_array_square(x_values, y_values, dimension);
                if (return_code == SORTING_MACHINE_SUCCESS) {
                    return_code = solve_maximal_cost_assignment();
                }
                break;
            case SORTING_METHOD_COMPRESSION:  // Compression algorithm
                sort_ids_from_center(x_values, y_values, dimension);
                return_code = solve_compression(x_values, y_values, dimension, 8U, 4U);
                break;
            case SORTING_METHOD_PARTITIONED_COMPRESSION:  // Partitioned compression algortihm.
                sort_ids_from_center(x_values, y_values, dimension);
                return_code = solve_partitioned_compression(x_values, y_values, dimension);
                if (return_code == SORTING_MACHINE_PARTITION_FAILED) {
                    sorting_method = SORTING_METHOD_COMPRESSION;
                    return_code = solve_compression(x_values, y_values, dimension, 8U, 4U);
                }
                break;
        }
        return return_code;
    }

    /*
     * SortingMachine::sort()
     *
     * Solves the sorting problem using preferred method.
     * 
     * Returns an integer based on the succes off the sorting function.
     * Per method the meaning of each return code can vary. Check the
     * documentation of each function for the meaning per method.
     */
    int SortingMachine::sort(
        std::vector<unsigned int> &start_x_values, 
        std::vector<unsigned int> &start_y_values, 
        std::vector<unsigned int> &end_x_values, 
        std::vector<unsigned int> &end_y_values, 
        unsigned int dimension
    ) {
        int return_code;
        switch (sorting_method) {
            case SORTING_METHOD_AUTO:  // Auto-detect
                if (number_target > 300U) {
                    sort_ids_from_center(start_x_values, start_y_values, dimension);
                    sorting_method = SORTING_METHOD_PARTITIONED_COMPRESSION;
                    return_code = solve_partitioned_compression(start_x_values, start_y_values, dimension);
                    if (return_code == SORTING_MACHINE_PARTITION_FAILED) {
                        sorting_method = SORTING_METHOD_COMPRESSION;
                        return_code = solve_compression(start_x_values, start_y_values, dimension);
                    }
                }
                else {
                    sorting_method = SORTING_METHOD_HUNGARIAN;
                    return_code = create_cost_array_square(start_x_values, start_y_values, dimension);
                    if (return_code == SORTING_MACHINE_SUCCESS) {
                        return_code = solve_maximal_cost_assignment();
                    }
                }
                break;
            case SORTING_METHOD_HUNGARIAN:  // Hungarian algorithm
                return_code = create_cost_array(
                    start_x_values, 
                    start_y_values,
                    end_x_values,
                    end_y_values,
                    dimension
                );
                if (return_code == SORTING_MACHINE_SUCCESS) {
                    return_code = solve_maximal_cost_assignment();
                }
                break;
            case SORTING_METHOD_COMPRESSION:  // Compression algorithm
                sort_ids_from_center(start_x_values, start_y_values, dimension);
                return_code = solve_compression(start_x_values, start_y_values, dimension, 8U, 4U);
                break;
            case SORTING_METHOD_PARTITIONED_COMPRESSION:  // Partitioned compression algortihm.
                sort_ids_from_center(start_x_values, start_y_values, dimension);
                return_code = solve_partitioned_compression(start_x_values, start_y_values, dimension);
                if (return_code == SORTING_MACHINE_PARTITION_FAILED) {
                    sorting_method = SORTING_METHOD_COMPRESSION;
                    return_code = solve_compression(start_x_values, start_y_values, dimension, 8U, 4U);
                }
                break;
        }
        return return_code;
    }

    /*
    * int SortingMachine::create_cost_array()
    *
    * Creates the cost array based on two input geometries given by the
    * x- and y-values as arguments to this function. At the end it inverts
    * the cost to make the algorithm later minimize the cost instead of
    * maximize. 
    * 
    * Parameters:
    * x_values_start : std::vector<unsigned int>
    *   The list of x-values in the starting geometry.
    * y_values_start : std::vector<unsigned int>
    *   The list of y-values in the starting geometry.
    * x_values_end : std::vector<unsigned int>
    *   The x-values in the end geometry.
    * y_values_end : std::vector<unsigned int>
    *   The y-values in the end geometry.
    * 
    * Returns:
    * return_code : int
    *   A return code based on the return codes defined in utilities.h.
    */
    int SortingMachine::create_cost_array(
        std::vector<unsigned int> &x_values_start, 
        std::vector<unsigned int> &y_values_start,
        std::vector<unsigned int> &x_values_end, 
        std::vector<unsigned int> &y_values_end, 
        unsigned int dimension
    ) {
        int cost = 0;
        int max_cost = 0;

        int xi, yi, xj, yj;
        for (auto i=0; i < number_loaded; i++) {
            for (auto j=0; j < number_loaded; j++) {
                if (j < number_target) {
                    xi = unroll(x_values_start[loaded_id[i]], dimension);
                    yi = unroll(y_values_start[loaded_id[i]], dimension);
                    xj = unroll(x_values_end[target_id[j]], dimension);
                    yj = unroll(y_values_end[target_id[j]], dimension);
                    cost = (xj - xi) * (xj - xi) + (yj - yi) * (yj - yi);
                }
                else {
                    cost = 0;
                } 
                cost_array[i * number_loaded + j] = cost;
                max_cost = std::max(cost, max_cost);    
            }
        }

        for (auto i=0; i < number_loaded; i++) {
            for (auto j=0; j < number_loaded; j++) {
                if (j < number_target) {
                    cost_array[i * number_loaded + j] = max_cost - cost_array[i * number_loaded + j];
                }
                else {
                    cost_array[i * number_loaded + j] = 0;
                }
            }
        }
        return SORTING_MACHINE_SUCCESS;
    }

    /*
    * int SortingMachine::create_cost_array_square()
    *
    * See SortingMachine::create_cost_array() for more information. This 
    * implementation assumes the geometries of the start and the end
    * are the same.
    */
    int SortingMachine::create_cost_array_square(
        std::vector<unsigned int> &x_values, 
        std::vector<unsigned int> &y_values, 
        unsigned int dimension
    ) {
        int cost = 0;
        int max_cost = 0;

        int xi, yi, xj, yj;
        for (auto i=0; i < number_loaded; i++) {
            for (auto j=0; j < number_loaded; j++) {
                if (j < number_target) {
                    xi = unroll(x_values[loaded_id[i]], dimension);
                    yi = unroll(y_values[loaded_id[i]], dimension);
                    xj = unroll(x_values[target_id[j]], dimension);
                    yj = unroll(y_values[target_id[j]], dimension);
                    cost = (xj - xi) * (xj - xi) + (yj - yi) * (yj - yi);
                }
                else {
                    cost = 0;
                } 
                cost_array[i * number_loaded + j] = cost;
                max_cost = std::max(cost, max_cost);    
            }
        }

        for (auto i=0; i < number_loaded; i++) {
            for (auto j=0; j < number_loaded; j++) {
                if (j < number_target) {
                    cost_array[i * number_loaded + j] = max_cost - cost_array[i * number_loaded + j];
                }
                else {
                    cost_array[i * number_loaded + j] = 0;
                }
            }
        }
        return SORTING_MACHINE_SUCCESS;
    }

    int SortingMachine::create_cost_array_square(
        std::vector<int> &cost_array_partition, std::vector<unsigned int> &loaded_partition, std::vector<unsigned int> &target_partition, 
        std::vector<unsigned int> &x_values, std::vector<unsigned int> &y_values, unsigned int dimension
    ) {
        cost_array_partition.clear();
        cost_array_partition.reserve(loaded_partition.size() * loaded_partition.size());

        //unsigned int n = round(sqrt(number_of_spots));
        int cost = 0;
        int max_cost = 0;
        size_t number_in_partition = loaded_partition.size();
        size_t targets_in_partition = target_partition.size(); 

        int xi, yi, xj, yj;
        for (auto i=0; i < number_in_partition; i++) {
            for (auto j=0; j < number_in_partition; j++) {
                if (j < targets_in_partition) {
                    xi = unroll(x_values[loaded_id[loaded_partition.at(i)]], dimension);
                    yi = unroll(y_values[loaded_id[loaded_partition.at(i)]], dimension);
                    xj = unroll(x_values[target_id[target_partition.at(j)]], dimension);
                    yj = unroll(y_values[target_id[target_partition.at(j)]], dimension);
                    cost = (xj - xi) * (xj - xi) + (yj - yi) * (yj - yi);
                    //std::cout << "Cost for (" << xj << ", " << yj << ") <-> (" << xi << ", " << yi << ") is " << cost << std::endl;
                }
                else {
                    cost = 0;
                } 
                cost_array_partition[i * number_in_partition + j] = cost;
                max_cost = std::max(cost, max_cost);    
            }
        }

        // Invert since we want to solve for minimal cost.
        for (auto i=0; i < number_in_partition; i++) {
            for (auto j=0; j < number_in_partition; j++) {
                if (j < targets_in_partition) {
                    cost_array_partition[i * number_in_partition + j] = max_cost - cost_array_partition[i * number_in_partition + j];
                }
                else {
                    cost_array_partition[i * number_in_partition + j] = 0;
                }
            }
        }
        return SORTING_MACHINE_SUCCESS;
    }

    int SortingMachine::solve_maximal_cost_assignment() {

        /* Initialize the labels and matching. Basic description:
            atom_labels, trap_labels			Cost vertices.
            atom_matching, trap_matching		Vectors to save the edges to which vertex one is matched.
            slack, slack_vertex				slack (excess cost) vector
            paths							Augmented paths memory.
            atom_matched, trap_matched		Boolean vectors to check whether atoms and traps are matched.
        */
        int atom_labels[MAX_NUMBER_SPOTS]; 
        int trap_labels[MAX_NUMBER_SPOTS];
        int atom_matching[MAX_NUMBER_SPOTS];
        int trap_matching[MAX_NUMBER_SPOTS];
        int slack[MAX_NUMBER_SPOTS];
        int slack_vertex[MAX_NUMBER_SPOTS];
        int paths [MAX_NUMBER_SPOTS];
        bool atom_matched[MAX_NUMBER_SPOTS];
        bool trap_matched[MAX_NUMBER_SPOTS];

        for (auto i = 0; i < number_loaded; ++i) {
            atom_labels[i] = 0;
            trap_labels[i] = 0;
            atom_matching[i] = -1;
            trap_matching[i] = -1;
            for (auto j = 0; j < number_loaded; ++j) { 
                atom_labels[i] = std::max(atom_labels[i], cost_array[i * number_loaded + j]);
            }
            //std::cout << atom_labels[i] << " " << atom_matching[i] << std::endl;
        }

        // Use augmenting paths algorithm to add matches.
        for (auto matches = 0; matches < number_loaded; ++matches) {
            std::deque<unsigned int> queue;

            // Initialize the parameters.
            for (auto i = 0; i < number_loaded; ++i) {
                atom_matched[i] = false;
                trap_matched[i] = false;
                slack[i] = (std::numeric_limits<int>::max)();
                slack_vertex[i] = -1;
                paths[i] = -1;
            }

            // Compute root vertex and slack.
            for (auto i = 0; i < number_loaded; ++i) {
                if (atom_matching[i] == -1) {
                    queue.push_back(i);
                    atom_matched[i] = true;

                    // Compute slack.
                    for (auto k = 0; k < number_loaded; ++k) {
                        if (atom_labels[i] + trap_labels[k] - cost_array[i * number_loaded + k] < slack[k]) {
                            slack[k] = atom_labels[i] + trap_labels[k] - cost_array[i * number_loaded + k];
                            slack_vertex[k] = i;
                        }
                    }
                    break;
                }
            }

            // Loop until we have found an augmenting path.
            auto i_start = 0;
            auto j_start = 0;
            auto found_path = false;
            while (!found_path) {
                while (queue.size() > 0 && !found_path) {
                    const auto i = queue.front();
                    queue.pop_front();
                    
                    for (auto j = 0; j < number_loaded; j++) {
                        if (cost_array[i * number_loaded + j] == atom_labels[i] + trap_labels[j] && !trap_matched[j]) {
                            if (trap_matching[j] == -1) {
                                i_start = i;
                                j_start = j;
                                found_path = true;
                                break;
                            }

                            trap_matched[j] = true;
                            queue.push_back(trap_matching[j]);
                            
                            paths[trap_matching[j]] = i;
                            atom_matched[trap_matching[j]] = true;
                            // Compute slack.
                            for (auto k = 0; k < number_loaded; k++) {
                                if (atom_labels[trap_matching[j]] + trap_labels[k] - cost_array[trap_matching[j] * number_loaded + k] < slack[k]) {
                                    slack[k] = atom_labels[trap_matching[j]] + trap_labels[k] - cost_array[trap_matching[j] * number_loaded + k];
                                    slack_vertex[k] = trap_matching[j];
                                }
                            }
                        }
                    }
                }
                if (found_path)
                    break;

                // No path found, so we update labels. This is the hungarian algorithm.
                auto delta = (std::numeric_limits<int>::max)();
                for (auto k = 0; k < number_loaded; ++k) {
                    if (!trap_matched[k])
                        delta = std::min(delta, slack[k]);
                }
                for (auto k = 0; k < number_loaded; ++k) {
                    if (atom_matched[k])
                        atom_labels[k] -= delta;

                    if (trap_matched[k])
                        trap_labels[k] += delta;
                    else
                        slack[k] -= delta;
                }

                // Now that we have updated the labels, there should be another free vertex.
                // We loop until we find this vertex and start the algorithm from there.
                queue.clear();
                for (auto j = 0; j < number_loaded; j++) {
                    if (!trap_matched[j] && slack[j] == 0) {
                        if (trap_matching[j] == -1) {
                            i_start = slack_vertex[j];
                            j_start = j;
                            found_path = true;
                            break;
                        }

                        trap_matched[j] = true;
                        if (!atom_matched[trap_matching[j]]) {
                            queue.push_back(trap_matching[j]);

                            atom_matched[trap_matching[j]] = true;
                            paths[trap_matching[j]] = slack_vertex[j];
                            // Compute slack for this new path.
                            for (auto k = 0; k < number_loaded; k++) {
                                if (atom_labels[trap_matching[j]] + trap_labels[k] - cost_array[trap_matching[j] * number_loaded + k] < slack[k]) {
                                    slack[k] = atom_labels[trap_matching[j]] + trap_labels[k] - cost_array[trap_matching[j] * number_loaded + k];
                                    slack_vertex[k] = trap_matching[j];
                                }
                            }
                        }
                    }
                }
            }

            // Flip edges along the path, adding one more item. This is part of the augmenting paths.
            for (
                int ci = i_start, cj = j_start, tj;
                ci != -1;
                ci = paths[ci], cj = tj
                ) {
                tj = atom_matching[ci];
                trap_matching[cj] = ci;
                atom_matching[ci] = cj;
            }
        }

        for (auto i=0; i<number_loaded; i++) {
            // Check if we don't assign the extra values of excess loaded atoms.
            if (atom_matching[i] < number_target) {
                mapping[loaded_id[i]] = target_id[atom_matching[i]];
            }
        }
        return SORTING_MACHINE_SUCCESS;
    }

    int SortingMachine::solve_maximal_cost_assignment(
        std::vector<int> &cost_array_layer, std::vector<unsigned int> &loaded_layer, std::vector<unsigned int> &target_layer) {

        /* Initialize the labels and matching. Basic description:
            atom_labels, trap_labels			Cost vertices.
            atom_matching, trap_matching		Vectors to save the edges to which vertex one is matched.
            slack, slack_vertex				slack (excess cost) vector
            paths							Augmented paths memory.
            atom_matched, trap_matched		Boolean vectors to check whether atoms and traps are matched.
        */

        size_t number_in_layer = loaded_layer.size();
        
        int atom_labels[MAX_NUMBER_LAYER]; 
        int trap_labels[MAX_NUMBER_LAYER];
        int atom_matching[MAX_NUMBER_LAYER];
        int trap_matching[MAX_NUMBER_LAYER];
        int slack[MAX_NUMBER_LAYER];
        int slack_vertex[MAX_NUMBER_LAYER];
        int paths [MAX_NUMBER_LAYER];
        bool atom_matched[MAX_NUMBER_LAYER];
        bool trap_matched[MAX_NUMBER_LAYER];

        for (auto i = 0; i < number_in_layer; ++i) {
            atom_labels[i] = 0;
            trap_labels[i] = 0;
            atom_matching[i] = -1;
            trap_matching[i] = -1;
            for (auto j = 0; j < number_in_layer; ++j) { 
                atom_labels[i] = std::max(atom_labels[i], cost_array_layer[i * number_in_layer + j]);
            }
            //std::cout << atom_labels[i] << " " << atom_matching[i] << std::endl;
        }

        // Use augmenting paths algorithm to add matches.
        for (auto matches = 0; matches < number_in_layer; ++matches) {
            std::deque<unsigned int> queue;

            // Initialize the parameters.
            for (auto i = 0; i < number_in_layer; ++i) {
                atom_matched[i] = false;
                trap_matched[i] = false;
                slack[i] = (std::numeric_limits<int>::max)();
                slack_vertex[i] = -1;
                paths[i] = -1;
            }

            // Compute root vertex and slack.
            for (auto i = 0; i < number_in_layer; ++i) {
                if (atom_matching[i] == -1) {
                    queue.push_back(i);
                    atom_matched[i] = true;

                    // Compute slack.
                    for (auto k = 0; k < number_in_layer; ++k) {
                        if (atom_labels[i] + trap_labels[k] - cost_array_layer[i * number_in_layer + k] < slack[k]) {
                            slack[k] = atom_labels[i] + trap_labels[k] - cost_array_layer[i * number_in_layer + k];
                            slack_vertex[k] = i;
                        }
                    }
                    break;
                }
            }

            // Loop until we have found an augmenting path.
            auto i_start = 0;
            auto j_start = 0;
            auto found_path = false;
            while (!found_path) {
                while (queue.size() > 0 && !found_path) {
                    const auto i = queue.front();
                    queue.pop_front();
                    
                    for (auto j = 0; j < number_in_layer; j++) {
                        if (cost_array_layer[i * number_in_layer + j] == atom_labels[i] + trap_labels[j] && !trap_matched[j]) {
                            if (trap_matching[j] == -1) {
                                i_start = i;
                                j_start = j;
                                found_path = true;
                                break;
                            }

                            trap_matched[j] = true;
                            queue.push_back(trap_matching[j]);
                            
                            paths[trap_matching[j]] = i;
                            atom_matched[trap_matching[j]] = true;
                            // Compute slack.
                            for (auto k = 0; k < number_in_layer; k++) {
                                if (atom_labels[trap_matching[j]] + trap_labels[k] - cost_array_layer[trap_matching[j] * number_in_layer + k] < slack[k]) {
                                    slack[k] = atom_labels[trap_matching[j]] + trap_labels[k] - cost_array_layer[trap_matching[j] * number_in_layer + k];
                                    slack_vertex[k] = trap_matching[j];
                                }
                            }
                        }
                    }
                }
                if (found_path)
                    break;

                // No path found, so we update labels. This is the hungarian algorithm.
                auto delta = (std::numeric_limits<int>::max)();
                for (auto k = 0; k < number_in_layer; ++k) {
                    if (!trap_matched[k])
                        delta = std::min(delta, slack[k]);
                }
                for (auto k = 0; k < number_in_layer; ++k) {
                    if (atom_matched[k])
                        atom_labels[k] -= delta;

                    if (trap_matched[k])
                        trap_labels[k] += delta;
                    else
                        slack[k] -= delta;
                }

                // Now that we have updated the labels, there should be another free vertex.
                // We loop until we find this vertex and start the algorithm from there.
                queue.clear();
                for (auto j = 0; j < number_in_layer; j++) {
                    if (!trap_matched[j] && slack[j] == 0) {
                        if (trap_matching[j] == -1) {
                            i_start = slack_vertex[j];
                            j_start = j;
                            found_path = true;
                            break;
                        }

                        trap_matched[j] = true;
                        if (!atom_matched[trap_matching[j]]) {
                            queue.push_back(trap_matching[j]);

                            atom_matched[trap_matching[j]] = true;
                            paths[trap_matching[j]] = slack_vertex[j];
                            // Compute slack for this new path.
                            for (auto k = 0; k < number_in_layer; k++) {
                                if (atom_labels[trap_matching[j]] + trap_labels[k] - cost_array_layer[trap_matching[j] * number_in_layer + k] < slack[k]) {
                                    slack[k] = atom_labels[trap_matching[j]] + trap_labels[k] - cost_array_layer[trap_matching[j] * number_in_layer + k];
                                    slack_vertex[k] = trap_matching[j];
                                }
                            }
                        }
                    }
                }
            }

            // Flip edges along the path, adding one more item. This is part of the augmenting paths.
            for (
                int ci = i_start, cj = j_start, tj;
                ci != -1;
                ci = paths[ci], cj = tj
                ) {
                tj = atom_matching[ci];
                trap_matching[cj] = ci;
                atom_matching[ci] = cj;
            }
        }

        for (auto i=0; i<number_in_layer; i++) {
            // Check if we don't assign the extra values of excess loaded atoms.
            if (atom_matching[i] < target_layer.size()) {
                mapping[loaded_id[loaded_layer.at(i)]] = target_id[target_layer.at(atom_matching[i])];
                // std::cout << loaded_id[loaded_layer.at(i)] << " --> " << target_id[target_layer.at(atom_matching[i])] << std::endl;
            }
        }
        return SORTING_MACHINE_SUCCESS;
    }

    /*
     * SortingMachine::parse_target_str()
     * 
     * Parses a bitstring to load the bool array of target atoms.
     * 
     * Parameters:
     * set_target_str : std::string_view
     *    The target bitstring.
     */
    void SortingMachine::parse_target_str(
        std::string_view set_target_str
    ) {
        if (set_target_str.length() > MAX_NUMBER_SPOTS) {
            std::cout << "Error in SortingMachine::parse_target_str: target string too long for maximum number of spots." << std::endl;
            system("pause");
            exit(1);
        }

        for (unsigned int i=0; i < MAX_NUMBER_SPOTS; i++) {
            is_target[i] = false;
        }

        number_target = 0;
        size_t pos = set_target_str.find("1", 0);
        while(pos != std::string::npos) {
            is_target[pos] = true;
            target_id[number_target] = pos;
            number_target++;
            pos = set_target_str.find("1", pos+1);
        }
    };

    /*
     * SortingMachine::parse_loaded_str()
     * 
     * Parses a bitstring to load the bool array of loaded atoms.
     */
    void SortingMachine::parse_loaded_str(
        std::string_view set_loaded_str
    ) {
        if (set_loaded_str.length() > MAX_NUMBER_SPOTS) {
            std::cout << "Error in SortingMachine::parse_target_str: target string too long for maximum number of spots." << std::endl;
            system("pause");
            exit(1);
        }

        reset();

        number_of_atoms = set_loaded_str.length();
        number_loaded = 0;
        size_t pos = set_loaded_str.find("1", 0);
        while(pos != std::string::npos) {
            is_loaded[pos] = true;
            loaded_id[number_loaded] = pos;
            pos = set_loaded_str.find("1", pos+1);
            number_loaded++;
        }
    };

    /*
    * SortingMachine::solve_partitioned_compression()
    *
    * Solves the problem with a partitioned compression algorithm.
    * Assuming no cost matrix has been assembled, but both the loaded
    * and target strings are parsed, the function does several steps:
    *   1. It partitions the data into multiple blocks, e.g. four
    *       quadrants of a square. Both loaded and target traps
    *       are sorted in this way.
    *   2. Each partition is sorted on the distance to the center
    *       or apex of that partition.
    *   3. For each parition the Hungarian algorithm is solved,
    *       resulting in a total mapping.
    * 
    * Returns an integer based upon the success of the function. 
    * Return codes are:
    *   0. Successful execution
    *   1. Partition failed because of atom number, did normal compression
    *   2. Other error
    */
    int SortingMachine::solve_partitioned_compression(
        std::vector<unsigned int> &x_values, std::vector<unsigned int> &y_values, unsigned int dimension
    ) {
        // Parition the data into four 'pizza slices'
        std::vector<unsigned int> loaded_partition_up;
        std::vector<unsigned int> target_partition_up;
        std::vector<unsigned int> loaded_partition_right;
        std::vector<unsigned int> target_partition_right;
        std::vector<unsigned int> loaded_partition_down;
        std::vector<unsigned int> target_partition_down;
        std::vector<unsigned int> loaded_partition_left;
        std::vector<unsigned int> target_partition_left;
        
        int xi, yi;
        for (auto i=0; i<number_loaded; ++i) {
            if (i < number_target) {
                xi = unroll(x_values[target_id[i]], dimension);
                yi = unroll(y_values[target_id[i]], dimension);
                if (yi >= -1 * xi && yi > xi) {
                    target_partition_up.emplace_back(i);
                }
                else if (yi > -1 * xi && yi <= xi) {
                    target_partition_right.emplace_back(i);
                }
                else if (yi <= -1 * xi && yi < xi) {
                    target_partition_down.emplace_back(i);
                }
                else {
                    target_partition_left.emplace_back(i);
                }
            }

            xi = unroll(x_values[loaded_id[i]], dimension);
            yi = unroll(y_values[loaded_id[i]], dimension);
            if (yi >= -1 * xi && yi > xi) {
                loaded_partition_up.emplace_back(i);
            }
            else if (yi > -1 * xi && yi <= xi) {
                loaded_partition_right.emplace_back(i);
            }
            else if (yi <= -1 * xi && yi < xi) {
                loaded_partition_down.emplace_back(i);
            }
            else {
                loaded_partition_left.emplace_back(i);
            }
        }

        if (target_partition_up.size() > loaded_partition_up.size() || target_partition_left.size() > loaded_partition_left.size() || 
                target_partition_down.size() > loaded_partition_down.size() || target_partition_right.size() > loaded_partition_right.size()) {
            std::cout << "Can't partition into slices. Switching to Compression algorithm" << std::endl;
            return SORTING_MACHINE_PARTITION_FAILED;
        }

        /* 
        * Now we solve the compression algorithm for each of the partitions. The scaling of the layers
        * determines how many rows of each triangle are checked per compression step. Some standards:
        *   - 1 row:    1,3,5,...       first_layer = 1U    layer_scaling = 2U
        *   - 2 rows:   4,12,20,...     first_layer = 4U    layer_scaling = 8U
        *   - 3 rows:   9,27,45,...     first_layer = 9U    layer_scaling = 18U
        *   - 4 rows:   16,48,80,...    first_layer = 16U   layer_scaling = 32U
        * Seems like the scaling is for n rows: first_layer = n^2 and layer_scaling = 2*n^2 for 
        * four quadrants. Should be some generalization.
        * 
        * In principle, the more rows, the shorter the max distance, but the longer the computation.
        * For 30x30 -> 20x20, it seems 3 rows are already enough.
        */ 
        unsigned int n_rows_compression = 2;
        unsigned int first_layer = 9U;
        unsigned int layer_scaling = 18U;
        solve_compression(x_values, y_values, loaded_partition_up, target_partition_up, dimension, layer_scaling, first_layer);
        solve_compression(x_values, y_values, loaded_partition_right, target_partition_right, dimension, layer_scaling, first_layer);
        solve_compression(x_values, y_values, loaded_partition_down, target_partition_down, dimension, layer_scaling, first_layer);
        solve_compression(x_values, y_values, loaded_partition_left, target_partition_left, dimension, layer_scaling, first_layer);
        return SORTING_MACHINE_SUCCESS;
    }

    /*
    * SortingMachine::solve_compression()
    *
    * Solves the problem with the compression algorithm from the center. 
    * This assumes that the vectors `x_values` and `y_values` are sorted
    * in such a way that the first element is closest to the center and
    * the last element is farthest.
    * 
    * At the end of the algorithm, there can be long moves because the 
    * outermost layers can be thin and far apart. For that reason the
    * function `trim_moves()` is used to break up the long moves into
    * multiple shorter ones.
    * 
    * Parameters: 
    * x_values : std::vector<unsigned int>
    *   The x_values of the tweezer traps.
    * y_values : std::vector<unsigned int>
    *   The y_values of the tweezer traps.
    * dimension : unsigned int
    *   The pixel number of the SLM chip. Used for rolling the coordinates.
    * layer_scaling : unsigned int
    *   The amount of new traps to consider per layer.
    * first_layer : unsigned int
    *   The amount of traps in the first layer.
    * 
    * Returns:
    * return_code : int
    *   An int based on the success of the operation. See `utilities.h` for more
    *   information on the different codes.
    */
    int SortingMachine::solve_compression(
        std::vector<unsigned int> &x_values, 
        std::vector<unsigned int> &y_values, 
        unsigned int dimension, 
        unsigned int layer_scaling,
        unsigned int first_layer
     ) {
        std::vector<unsigned int> loaded_layer;
        std::vector<unsigned int> target_layer;
        std::vector<int> cost_array_layer;

        unsigned int offset = 0;
        unsigned int number_in_layer = first_layer;
        while (offset < number_target) {
            loaded_layer.clear();
            target_layer.clear();

            for (auto i=offset; i < offset+number_in_layer && i<number_loaded; i++) {
                loaded_layer.emplace_back(i);
                if (i < number_target) {
                    target_layer.emplace_back(i);
                }
            }

            create_cost_array_square(cost_array_layer, loaded_layer, target_layer, x_values, y_values, dimension);
            solve_maximal_cost_assignment(cost_array_layer, loaded_layer, target_layer);

            offset += number_in_layer;
            number_in_layer += layer_scaling;
        }

        loaded_layer.clear();
        target_layer.clear();
        cost_array_layer.clear();

        return SORTING_MACHINE_SUCCESS;
    }

    /*
    * SortingMachine::solve_compression()
    *
    * Overloaded function that works for a specified loading partition
    * and a specified target partition.
    * 
    * Solves the problem with the compression algorithm from the center. 
    * This assumes that the vectors `x_values` and `y_values` are sorted
    * in such a way that the first element is closest to the center and
    * the last element is farthest.
    */
    int SortingMachine::solve_compression(
        std::vector<unsigned int> &x_values, 
        std::vector<unsigned int> &y_values, 
        std::vector<unsigned int> &loaded_partition, 
        std::vector<unsigned int> &target_partition, 
        unsigned int dimension, 
        unsigned int layer_scaling, 
        unsigned int first_layer
    ) {
        std::vector<unsigned int> loaded_layer;
        std::vector<unsigned int> target_layer;
        std::vector<int> cost_array_layer;

        size_t target_partition_size = target_partition.size();
        size_t loaded_partition_size = loaded_partition.size();

        unsigned int offset = 0;
        unsigned int number_in_layer = first_layer;
        
        //std::cout << "Assembling layers" << std::endl;
        while (offset < target_partition_size) {
            loaded_layer.clear();
            target_layer.clear();

            for (auto i=offset; i<offset + number_in_layer && i<loaded_partition_size; i++) {
                loaded_layer.emplace_back(loaded_partition[i]);
                if (i < target_partition_size) {
                    target_layer.emplace_back(target_partition[i]);
                }
            }

            //std::cout << loaded_layer.size() << " " << target_layer.size() << " " << offset << " " << number_in_layer << std::endl;

            if (loaded_layer.size() && target_layer.size()) {
                create_cost_array_square(cost_array_layer, loaded_layer, target_layer, x_values, y_values, dimension);
                solve_maximal_cost_assignment(cost_array_layer, loaded_layer, target_layer);
            }         

            offset += number_in_layer;
            if (layer_scaling + number_in_layer <= MAX_NUMBER_LAYER) {
                number_in_layer += layer_scaling;
            }
            else {
                number_in_layer = MAX_NUMBER_LAYER;
            }
        }

        loaded_layer.clear();
        target_layer.clear();
        cost_array_layer.clear();

        return SORTING_MACHINE_SUCCESS;
    }

    /*
    * SortingMachine::trim_moves()
    *
    * Breaks up the moves of the compression algorithm into more
    * shorter moves. See LSAP2 algorithm of Browaeys group for 
    * more explanation.
    * 
    * The code starts by making a vector of paths for each of the
    * mappings that we calculated. Then it will loop sequentially
    * through these to avoid collisions.
    */
    void SortingMachine::trim_moves() {

    }

    /* 
    * SortingMachine::get_next_positions()
    *
    * Gets the next values of the spots in the sorting interpolation. 
    */
    bool SortingMachine::get_next_positions(
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
    ) {
        if (iteration > (number_of_steps + number_turnoff_frames)) {
            return false;
        }

        for (unsigned int i=0; i < number_of_atoms; i++) {
            // Fade out and remove from should_be_used if not used.
            if (!is_loaded[i] || (!is_target[i] && mapping[i] < 0)) {
                if (iteration < number_turnoff_frames) {
                    amplitude_values[i] = start_amplitudes[i] * (1 - (float) iteration / (float) number_turnoff_frames);
                    if (amplitude_values[i] <= 0.f) {
                        amplitude_values[i] = 0.f;
                        should_be_used[i] = 0;
                    }
                }
                else {
                    amplitude_values[i] = 0.f;
                    should_be_used[i] = 0;
                }
                continue;
            }
            
            // If increase the amplitude value here to at least try to normalize.
            if (iteration <= number_turnoff_frames) {
                amplitude_values[i] = (1 + (float) iteration / (float) number_turnoff_frames * ((float) number_loaded / number_target - 1)) * start_amplitudes[i];
            }

            // Only change phase if already at the right position.
            if (mapping[i] == i) {
                phase_values[i] = start_phases[i] + (float) (iteration) * (
                    end_phases[mapping[i]] - start_phases[i]
                ) / ((float) (number_of_steps + number_turnoff_frames));
                continue;
            }

            // Update positions to iteration.
            if (iteration > number_turnoff_frames) {
                switch (pathing_method) {
                    case PATHING_METHOD_DIRECT:
                        x_values[i] = roll(
                            unroll(start_x_values[i], dimension) + std::round((float)(iteration - number_turnoff_frames) * (
                                unroll(start_x_values[mapping[i]], dimension) - unroll(start_x_values[i], dimension)
                            ) / ((float) number_of_steps)),
                            dimension
                        );
                        y_values[i] = roll(
                            unroll(start_y_values[i], dimension) + std::round((float)(iteration - number_turnoff_frames) * (
                                unroll(start_y_values[mapping[i]], dimension) - unroll(start_y_values[i], dimension)
                            ) / ((float) number_of_steps)),
                            dimension
                        );

                    default:
                        x_values[i] = roll(
                            unroll(start_x_values[i], dimension) + std::round((int)(iteration - number_turnoff_frames) * (
                                unroll(start_x_values[mapping[i]], dimension) - unroll(start_x_values[i], dimension)
                            ) / ((float) number_of_steps)),
                            dimension
                        );
                        y_values[i] = roll(
                            unroll(start_y_values[i], dimension) + std::round((int)(iteration - number_turnoff_frames) * (
                                unroll(start_y_values[mapping[i]], dimension) - unroll(start_y_values[i], dimension)
                            ) / ((float) number_of_steps)),
                            dimension
                        );
                }

                phase_values[i] = start_phases[i] + (float) (iteration) * (
                    end_phases[mapping[i]] - start_phases[i]
                ) / ((float) (number_of_steps + number_turnoff_frames));
            }   
        }
        return true;
    }

    /* 
    * SortingMachine::get_next_positions()
    *
    * Gets the next values of the spots in the sorting interpolation. 
    */
    bool SortingMachine::get_next_positions(
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
    ) {
        if (iteration > (number_of_steps + number_turnoff_frames)) {
            return false;
        }

        for (unsigned int i=0; i < number_of_atoms; i++) {
            // Fade out and remove from should_be_used if not used.
            if (!is_loaded[i] || (!is_target[i] && mapping[i] < 0)) {
                if (iteration < number_turnoff_frames) {
                    amplitude_values[i] = start_amplitudes[i] * (1 - (float) iteration / (float) number_turnoff_frames);
                    if (amplitude_values[i] <= 0.f) {
                        amplitude_values[i] = 0.f;
                        should_be_used[i] = 0;
                    }
                }
                else {
                    amplitude_values[i] = 0.f;
                    should_be_used[i] = 0;
                }
                continue;
            }
            
            // If increase the amplitude value here to at least try to normalize.
            if (iteration <= number_turnoff_frames) {
                amplitude_values[i] = start_amplitudes[i] + (float) iteration / (float) number_turnoff_frames * (end_amplitudes[mapping[i]] - start_amplitudes[i]);
            }

            // Only change phase if already at the right position.
            /*if (mapping[i] == i) {
                phase_values[i] = start_phases[i] + (float) (iteration) * (
                    end_phases[mapping[i]] - start_phases[i]
                ) / ((float) (number_of_steps + number_turnoff_frames));
                continue;
            }*/

            // Update positions to iteration.
            if (iteration > number_turnoff_frames) {
                switch (pathing_method) {
                    case PATHING_METHOD_DIRECT:
                        x_values[i] = roll(
                            unroll(start_x_values[i], dimension) + std::round((float)(iteration - number_turnoff_frames) * (
                                unroll(end_x_values[mapping[i]], dimension) - unroll(start_x_values[i], dimension)
                            ) / ((float) number_of_steps)),
                            dimension
                        );
                        y_values[i] = roll(
                            unroll(start_y_values[i], dimension) + std::round((float)(iteration - number_turnoff_frames) * (
                                unroll(end_y_values[mapping[i]], dimension) - unroll(start_y_values[i], dimension)
                            ) / ((float) number_of_steps)),
                            dimension
                        );

                    default:
                        x_values[i] = roll(
                            unroll(start_x_values[i], dimension) + std::round((int)(iteration - number_turnoff_frames) * (
                                unroll(end_x_values[mapping[i]], dimension) - unroll(start_x_values[i], dimension)
                            ) / ((float) number_of_steps)),
                            dimension
                        );
                        y_values[i] = roll(
                            unroll(start_y_values[i], dimension) + std::round((int)(iteration - number_turnoff_frames) * (
                                unroll(end_y_values[mapping[i]], dimension) - unroll(start_y_values[i], dimension)
                            ) / ((float) number_of_steps)),
                            dimension
                        );
                }

                phase_values[i] = start_phases[i] + (float) (iteration) * (
                    end_phases[mapping[i]] - start_phases[i]
                ) / ((float) (number_of_steps + number_turnoff_frames));
            }   
        }
        return true;
    }

    /*
    * SortingMachine::calculate_maximal_movement()
    */
    void SortingMachine::calculate_maximal_movement(
        std::vector<unsigned int> &x_values,
        std::vector<unsigned int> &y_values,
        unsigned int dimension
    ) {
        unsigned int distance;
        number_of_steps = 0;
        for (auto i=0; i<number_of_atoms; ++i) {
            if (mapping[i] == -1)
                continue;
            distance = std::max(
                std::abs(unroll(x_values[mapping[i]], dimension) - unroll(x_values[i], dimension)),
                std::abs(unroll(y_values[mapping[i]], dimension) - unroll(y_values[i], dimension))
            );
            number_of_steps = std::max(distance, number_of_steps);
        } 
    }

    /*
    * SortingMachine::calculate_maximal_movement()
    */
    void SortingMachine::calculate_maximal_movement(
        std::vector<unsigned int> &start_x_values,
        std::vector<unsigned int> &start_y_values,
        std::vector<unsigned int> &end_x_values,
        std::vector<unsigned int> &end_y_values,
        unsigned int dimension
    ) {
        unsigned int distance;
        number_of_steps = 0;
        for (auto i=0; i<number_of_atoms; ++i) {
            if (mapping[i] == -1)
                continue;
            distance = std::max(
                std::abs(unroll(end_x_values[mapping[i]], dimension) - unroll(start_x_values[i], dimension)),
                std::abs(unroll(end_y_values[mapping[i]], dimension) - unroll(start_y_values[i], dimension))
            );
            number_of_steps = std::max(distance, number_of_steps);
        } 
    }

    /*
    * SortingMachine::create_target_string_center() 
    */
    std::string SortingMachine::create_target_string_center(
        unsigned int n, unsigned int m
    ) {
        std::string output = "";
        for (auto i=0; i<n; ++i) {
            for (auto j=0; j<n; ++j) {
                if (i < (n - m) / 2. || i >= (n - (n - m) / 2.)) {
                    output += "0";
                }
                else {
                    if (j < (n - m) / 2. || j >= (n - (n - m) / 2.)) {
                        output += "0";
                    }
                    else {
                        output += "1";
                    }
                }
            }
        }
        return output;
    }

    /*
    * SortingMachine::create_loaded_string_random()
    */
    std::string SortingMachine::create_loaded_string_random(
        unsigned int n, float loading_probability
    ) {
        std::string output = "";
        srand(std::chrono::steady_clock::now().time_since_epoch().count());
        for (auto i=0; i<n; ++i) {
            if (((float) rand()) / RAND_MAX <= loading_probability) {
                output += "1";
            }
            else {
                output += "0";
            }
        }
        return output;
    }

    /*
    * SortingMachine::sort_ids_from_center()
    *
    * Sorts the `loaded_id` and `target_id` array by distance to the 
    * center. The distance is measured by unrolling the x- and y-values
    * of the spot using given `dimension`.
    */
    void SortingMachine::sort_ids_from_center(
        std::vector<unsigned int> &x_values, std::vector<unsigned int> &y_values, unsigned int dimension
    ) {
        if (number_target <= 1 || number_loaded <= 1) {
            std::cout << "Not enough atoms..." << std::endl;
            return;
        }
        std::vector<unsigned int> sorted_loaded_ids;
        std::vector<unsigned int> sorted_target_ids;
        sorted_loaded_ids.emplace_back(loaded_id[0]);
        sorted_target_ids.emplace_back(target_id[0]);

        int xi, yi, xj, yj;
        for (unsigned int i=1; i<number_loaded; i++) {
            if (i < number_target) {
                xi = unroll(x_values[target_id[i]], dimension);
                yi = unroll(y_values[target_id[i]], dimension);

                for (unsigned int j=0; j<sorted_target_ids.size(); j++) {
                    xj = unroll(x_values[sorted_target_ids[j]], dimension);
                    yj = unroll(y_values[sorted_target_ids[j]], dimension);
                    if (xi*xi + yi*yi <= xj*xj + yj*yj) {
                        sorted_target_ids.insert(sorted_target_ids.begin() + j, target_id[i]);
                        break;
                    }
                }
            }

            xi = unroll(x_values[loaded_id[i]], dimension);
            yi = unroll(y_values[loaded_id[i]], dimension);
            for (unsigned int j=0; j<sorted_loaded_ids.size(); j++) {
                xj = unroll(x_values[sorted_loaded_ids[j]], dimension);
                yj = unroll(y_values[sorted_loaded_ids[j]], dimension);
                if (xi*xi + yi*yi <= xj*xj + yj*yj) {
                    sorted_loaded_ids.insert(sorted_loaded_ids.begin() + j, loaded_id[i]);
                    break;
                }
            }
        }

        /*std::cout << "Targets sorted from center: ";
        for (auto x: sorted_target_ids) {
            std::cout << x << " ";
        }
        std::cout << std::endl;
        std::cout << "Loaded sorted from center: ";
        for (auto x: sorted_loaded_ids) {
            std::cout << x << " ";
        }
        std::cout << std::endl;*/

        // Now we resave the loaded_ids and target_ids.
        for (unsigned int i=0; i<sorted_loaded_ids.size(); i++) {
            if (i < sorted_target_ids.size()) {
                target_id[i] = sorted_target_ids[i];
            }
            loaded_id[i] = sorted_loaded_ids[i];
        }
        return;
    }
}