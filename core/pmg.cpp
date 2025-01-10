/*
* pmg.cpp
*/

#include "pmg.h"

namespace PMG
{
    

    /*
    * PMG::initiate_SLM()
    *
    * Initiates the SDK of the SLM. We have to use the Blink_C_wrapper.h files.
    * Somehow the normal Blink_SDK.h does not work, but the manual of the SLM
    * does also recommend using the Blink_C_wrapper.h file. This has also Python
    * and Matlab bindings in it, so you have to have other .dlls than just 
    * Blink_SDK.dll in the include directory (e.g. python38.dll). 
    * 
    * Return success bool upon connecting and finding more than one board.
    * 
    */
    bool PMG::initiate_SLM() {
        unsigned int bits_per_pixel = 12U;   //12U is used for 1920x1152 SLM, 8U used for the small 512x512
        bool is_nematic_type = true;
        bool RAM_write_enable = true;
        bool use_GPU_if_available = true;
        unsigned int n_boards_found = 0U;
        int success = true;

        // Try to connect to SLM.
        Create_SDK(bits_per_pixel, &n_boards_found, &success, is_nematic_type, RAM_write_enable, use_GPU_if_available, 10U, 0);

        success &= n_boards_found > 0;
        SLM_connected = success == 1;

        if (!SLM_connected) {
            return false;
        }

        // Check dimensions of the image with the SLM
        if ((int) width != Get_image_width(SLM_BOARD_NUMBER) || (int) height != Get_image_height(SLM_BOARD_NUMBER)) {
            std::cout << "SLM image dimensions do not match PMG dimensions " << width << "x" << height << std::endl;
            return false;
        }
        SetRampDelay(SLM_BOARD_NUMBER, 10U);
        SetPreRampSlope(SLM_BOARD_NUMBER, 7U);
        
        return true;
    }

    /*
    * PMG::load_LUT_from_file()
    * 
    * Loads a look-up-table which will be used to create voltages from the pixel values. 
    * 
    * Parameters:
    * path_to_file : std::string
    *   String containing the path to the LUT.
    * 
    * Returns:
    * success : bool
    *   Wanted to do int first, but have no idea about the codes. 0 is successful. 
    *   No idea about other error codes reding Meadowlark docs.
    */    
    bool PMG::load_LUT_from_file(std::string path_to_file) {
        if (!SLM_connected) {
            std::cout << "SLM is not connected. Initiate SDK first." << std::endl;
            return false;
        }

        // Try to load the LUT file into the SLM. 
        if (Load_LUT_file(SLM_BOARD_NUMBER, (char *) path_to_file.c_str()) != 0) {
            std::cout<< "Error loading LUT: " << Get_last_error_message() << std::endl;
            return false;
        }
        return true;
    }
    
    /*
    * PMG::load_corrections()
    *
    * Loads a correction image containing WFC and Zernike polynomials.
    * Should be an 8-bit grayscale BMP. 
    * 
    * Parameters:
    * path_to_file : std::string
    *   Pointer to a path to the file location.
    * 
    * Returns:
    * success : bool
    *   Success of loading.
    */
    bool PMG::load_corrections(std::string path_to_file) {
        CBitmap mask_bmp;

        if (!mask_bmp.Load(path_to_file.c_str())) {
            std::cout << "Error loading file at " << path_to_file << std::endl;
            return false;
        }

        mask_bmp.GetBits(correction_mask_8bit);
        return copy_corrections_to_buffer();
    }

    /*
    * PMG::load_corrections()
    *
    * Overload for many correction phases.
    * 
    * Parameters:
    * paths_to_files : std::vector <std::string>
    *   Vector containing multiple paths to files.
    * 
    * Returns:
    * success : bool
    *   Success of loading.
    */
    bool PMG::load_corrections(std::vector<std::string> paths_to_files) {
        uint8_t *aggregate_mask_8bit = new uint8_t[width * height];
        for (auto i=0; i<paths_to_files.size(); i++) {
            CBitmap mask_bmp;

            if (!mask_bmp.Load(paths_to_files[i].c_str())) {
                std::cout << "Error loading file at " << paths_to_files[i] << std::endl;
                return false;
            }

            mask_bmp.GetBits(aggregate_mask_8bit);

            for (auto j=0; j<width * height; j++) {
                if (i == 0) {
                    correction_mask_8bit[j] = aggregate_mask_8bit[j];
                }
                else {
                    correction_mask_8bit[j] += aggregate_mask_8bit[j];
                }
            }
        }


        delete aggregate_mask_8bit;
        return copy_corrections_to_buffer();
    }

    /*
    * PMG::roll_mask_x()
    *
    * Shifts the mask in mask_8_bit by dx pixels in the x-direction.
    */
    void PMG::roll_mask_x(int dx) {
        uint8_t *mask = new uint8_t[width * height];
        for (auto i=0; i<width*height; i++) {
            mask[i] = mask_8bit[(unsigned int) (i + dx) % (width * height)];
        }
        for (auto i=0; i < width*height; i++) {
            mask_8bit[i] = mask[i];
        }
        delete mask;
    }

    /*
    * PMG::roll_mask_y()
    *
    * Shifts the mask in mask_8_bit by dy pixels in the y-direction.
    */
    void PMG::roll_mask_y(int dy) {
        uint8_t *mask = new uint8_t[width * height];
        for (auto i=0; i<width*height; i++) {
            mask[i] = mask_8bit[(unsigned int) (i + dy * width) % (width * height)];
        }
        for (auto i=0; i < width*height; i++) {
            mask_8bit[i] = mask[i];
        }
        delete mask;
    }
    
    /*
    * PMG::add_corrections_to_mask()
    *
    * Adds the corrections to the mask to display and stores the final
    * mask in mask_to_display_8bit.
    */
    void PMG::add_corrections_to_mask() {
        for (auto i=0; i< width * height; i++) {
            mask_to_display_8bit[i] = mask_8bit[i] + correction_mask_8bit[i];
        }
    }

    /*
    * PMG::remove_corrections_from_mask()
    *
    * Removes the corrections from mask_to_display_8bit.
    */
    void PMG::remove_corrections_from_mask() {
        for (auto i=0; i< width * height; i++) {
            mask_8bit[i] = mask_to_display_8bit[i] - correction_mask_8bit[i];
        }
    }

    /*
    * PMG::display_mask_on_SLM() 
    *
    * Displays the current 8-bit mask in the memory on the SLM. After this 
    * call it is safe to start a new DMA. No external triggers are enabled 
    * here and no outputs are generated.
    */
    bool PMG::display_mask_on_SLM() {
        int result = -1;
        result = Write_image(SLM_BOARD_NUMBER, mask_to_display_8bit, width * height, false, false, false, false, 5000);
        //std::cout << "Result code after writing image: " << result << std::endl;
        //std::cout<< "Error: " << Get_last_error_message() << std::endl;
        result = ImageWriteComplete(SLM_BOARD_NUMBER, 5000);
        //std::cout << "Result code after writing image: " << result << std::endl;
        //std::cout<< "Error: " << Get_last_error_message() << std::endl;
        return true;
    }

    /*
    * PMG::display_mask_on_SLM_immediately() 
    *
    * Displays the current 8-bit mask in the memory on the SLM immediately.
    * This might give errors since we can be too fast. 
    * 
    * Update: 241105: This crashed our SLM after a few 10s of runs. Be careful!
    */
    bool PMG::display_mask_on_SLM_immediately() {
        int result = -1;
        result = Write_image(SLM_BOARD_NUMBER, mask_to_display_8bit, width * height, false, true, false, false, 5000);
        //std::cout << "Result code after writing image: " << result << std::endl;
        //std::cout<< "Error: " << Get_last_error_message() << std::endl;
        result = ImageWriteComplete(SLM_BOARD_NUMBER, 5000);
        //std::cout << "Result code after writing image: " << result << std::endl;
        //std::cout<< "Error: " << Get_last_error_message() << std::endl;
        return true;
    }

    /*
    * PMG::display_mask_on_SLM_ext_trigger() 
    *
    * Displays the current 8-bit mask in the memory on the SLM. After this call 
    * it is safe to start a new DMA. No external triggers are enabled here and
    * no outputs are generated.
    */
    bool PMG::display_mask_on_SLM_ext_trigger() {
        int result = -1;
        result = Write_image(SLM_BOARD_NUMBER, mask_to_display_8bit, width * height, true, false, false, false, 5000);
        //std::cout << "Waiting for trigger" << std::endl;

        result = ImageWriteComplete(SLM_BOARD_NUMBER, 5000);
        //std::cout << "Result code after writing image: " << result << std::endl;
        //std::cout<< "Error: " << Get_last_error_message() << std::endl;
        return true;
    }

    /*
    * PMG::display_mask_on_SLM_ext_trigger_immediately() 
    *
    * Waits for external trigger to go low to start displaying the image. This is
    * the displayed immediately. After this, it should be safe to send a new DMA.
    */
    bool PMG::display_mask_on_SLM_ext_trigger_immediately() {
        int result = -1;
        result = Write_image(SLM_BOARD_NUMBER, mask_to_display_8bit, width * height, true, true, false, false, 5000);
        //std::cout << "Waiting for trigger" << std::endl;

        result = ImageWriteComplete(SLM_BOARD_NUMBER, 5000);
        //std::cout << "Result code after writing image: " << result << std::endl;
        //std::cout<< "Error: " << Get_last_error_message() << std::endl;
        return true;
    }

    /*
    * PMG::create_square_pattern()
    *
    * **DEPRECATED**: 
    * This way of setting the spots is not the most semantic. 
    * Use `load_start_geometry()` and `load_end_geometry()`.
    * 
    * Creates a square pattern of spots of size num_x * num_y.
    * The spots are spaced by spot_distance.
    *
    * Parameters:
    * num_x : unsigned int
    *      Number of spots in x-direction.
    * num_y : unsigned int
    *      Number of spots in y-direction.
    * spot_distance : unsigned int
    *      Distance of the spots with respect to each other.
    */
    void PMG::create_square_pattern(
        unsigned int num_x,
        unsigned int num_y,
        unsigned int spot_distance)
    {
        // Clear previous spots.
        x_values.clear();
        y_values.clear();
        amplitude_values.clear();
        psi_values.clear();
        should_be_used.clear();

        start_x_values.clear();
        start_y_values.clear();
        start_amplitude_values.clear();
        start_phase_values.clear();
        end_phase_values.clear();


        // Fill the spot_list with Spots.
        for (unsigned int i = 0; i < num_x; i++)
        {
            for (unsigned int j = 0; j < num_y; j++)
            {

                unsigned int x = (int) (-((float) (num_x - 1) / 2.) * spot_distance + i * spot_distance) % width;
                unsigned int y = (int) (-((float) (num_y - 1) / 2.) * spot_distance + j * spot_distance) % height;

                //unsigned int index = i * num_y + j;
                x_values.emplace_back(x);
                y_values.emplace_back(y);
                amplitude_values.emplace_back(1. / num_x / num_y);
                psi_values.emplace_back(rand() / RAND_MAX * CL_M_PI * 2 - CL_M_PI);
                should_be_used.emplace_back(1);

                start_x_values.emplace_back(x);
                start_y_values.emplace_back(y);
                start_amplitude_values.emplace_back(1.  / num_x / num_y);
                start_phase_values.emplace_back(0.f);
                end_phase_values.emplace_back(0.f);
            }
        }

        number_of_spots = num_x * num_y;
    }

    /*
    * PMG::create_square_pattern()
    *
    * Creates a square pattern of spots of size num_x * num_y.
    * The spots are spaced by spot_distance.
    *
    * Parameters:
    * num_x : unsigned int
    *      Number of spots in x-direction.
    * num_y : unsigned int
    *      Number of spots in y-direction.
    * spot_distance : unsigned int
    *      Distance of the spots with respect to each other.
    */
    std::vector<Spot> PMG::create_spots_square(
        unsigned int num_x,
        unsigned int num_y,
        unsigned int spot_distance)
    {
        std::vector<Spot> spots;

        // Fill the spot_list with Spots.
        for (unsigned int i = 0; i < num_x; i++)
        {
            for (unsigned int j = 0; j < num_y; j++)
            {
                Spot spot;
                spot.x = (int) (-((float) (num_x - 1) / 2.) * spot_distance + i * spot_distance) % width;
                spot.y = (int) (-((float) (num_y - 1) / 2.) * spot_distance + j * spot_distance) % height;
                spot.a = 1. / num_x / num_y;
                spot.psi = rand() / RAND_MAX * CL_M_PI * 2 - CL_M_PI;
                spot.used = 1;
                spots.emplace_back(spot);
            }
        }
        return spots;
    }

    /*
    *
    *
    */
    void PMG::define_start() {
        for (auto i=0; i < number_of_spots; i++) {
            start_x_values[i] = x_values[i];
            start_y_values[i] = y_values[i];
            start_phase_values[i] = psi_values[i];
            start_amplitude_values[i] = amplitude_values[i];
        }
    }

    /*
    * PMG::shift_pattern()
    *
    * Shifts the existing spots by dx and dy.
    *
    * Parameters:
    * dx : int
    *      Number of pixels to move in x-direction.
    * dy : int
    *      Number of pixels to move in y-direction.
    */
    void PMG::shift_pattern(
        int dx,
        int dy)
    {
        for (auto i=0; i < number_of_spots; i++) {
            x_values[i] = (int) (x_values[i] + dx) % width;
            y_values[i] = (int) (y_values[i] + dy) % height;
        }
    }

    /*
    * PMG::shift_pattern()
    *
    * Shifts the existing spots by dx and dy.
    *
    * Parameters:
    * dx : int
    *      Number of pixels to move in x-direction.
    * dy : int
    *      Number of pixels to move in y-direction.
    * dpsi : float
    *      Change in phase of the tweezer.
    */
    void PMG::shift_pattern(
        int dx,
        int dy,
        float dpsi)
    {
        for (auto i=0; i < number_of_spots; i++) {
            x_values[i] = (int) (x_values[i] + dx) % width;
            y_values[i] = (int) (y_values[i] + dy) % height;
            psi_values[i] += dpsi;
        }
    }

    /*
    * PMG::spread_pattern()
    *
    * Shifts the existing spots one step inwards or outwards.
    *
    * Parameters:
    * direction : bool
    *     Whether the direction is inwards or outwards from the center.
    *     True means inwards
    */
    void PMG::spread_pattern(
        bool direction)
    {
        for (auto i=0; i < number_of_spots; i++) {
            if(x_values[i] != 0) {
                if (x_values[i] > width / 2) {
                    x_values[i] = direction ? (x_values[i] - 1) : (int) (x_values[i] + 1) % width;
                }
                else {
                    x_values[i] = direction ? (x_values[i] + 1) : (int) (x_values[i] - 1) % width;
                }
            }
            if(y_values[i] != 0) {
                if (y_values[i] > height / 2) {
                    y_values[i] = direction ? (y_values[i] - 1) : (int) (y_values[i] + 1) % height;
                }
                else {
                    y_values[i] = direction ? (y_values[i] + 1) : (int) (y_values[i] - 1) % height;
                }
            }
            psi_values[i] = start_phase_values[i];
        }
    }

    /*
    * PMG::spread_pattern()
    *
    * Shifts the existing spots one step inwards or outwards.
    *
    * Parameters:
    * direction : bool
    *     Whether the direction is inwards or outwards from the center.
    *     True means inwards
    */
    void PMG::spread_pattern(
        bool direction, unsigned int step, unsigned int number_of_steps)
    {
        for (auto i=0; i < number_of_spots; i++) {
            if(x_values[i] != 0) {
                if (x_values[i] > width / 2) {
                    x_values[i] = direction ? (start_x_values[i] - step) : (int) (start_x_values[i] + step) % width;
                }
                else {
                    x_values[i] = direction ? (start_x_values[i] + step) : (int) (start_x_values[i] - step) % width;
                }
            }
            if(y_values[i] != 0) {
                if (y_values[i] > height / 2) {
                    y_values[i] = direction ? (start_y_values[i] - step) : (int) (start_y_values[i] + step) % height;
                }
                else {
                    y_values[i] = direction ? (start_y_values[i] + step) : (int) (start_y_values[i] - step) % height;
                }
            }
            psi_values[i] = start_phase_values[i];
        }
    }

    /*
    * PMG::set_illumination
    *
    * Sets the illumination pattern for the PMG, based on given
    * choice. Options right now are:
    *      -0  : uniform intensity
    *      -1  : gaussian beam with diameter == min(height, max)
    *
    * Parameters:
    * option : unsigned int
    *      The choice in intensity pattern.
    */
    void PMG::set_illumination(
        unsigned int option)
    {
        float *input_light = new float [width * height];

        switch (option)
        {
        case 0:
            // Uniform illumination
            for (unsigned int i = 0; i < width; i++)
            {
                for (unsigned int j = 0; j < height; j++)
                {
                    input_light[(unsigned int) (i * height + j)] = 1. / width / height;
                }
            }
            break;
        case 1:
            // Uniform illumination
            for (unsigned int i = 0; i < width; i++)
            {
                for (unsigned int j = 0; j < height; j++)
                {
                    input_light[(unsigned int) (i * height + j)] = 1. / width / height;
                }
            }
            break;
        default:
            // Uniform illumination
            for (unsigned int i = 0; i < width; i++)
            {
                for (unsigned int j = 0; j < height; j++)
                {
                    input_light[(unsigned int) (i * height + j)] = 1. / width / height;
                }
            }
        }

        // Copy to the buffer.
        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->input_light_buffer, CL_TRUE, 0, sizeof(float) * width * height, input_light, 0, NULL, &gpu_handler->event_list[INDEX_WRITE_INPUT_LIGHT_BUFFER]
        );

        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error writing input light to buffer. OpenCL error code: " << gpu_handler->result << std::endl;
            system("pause");
            exit(1);
        }        

        delete input_light;
    }

    /*
    * PMG::set_target_buffer
    *
    * Sets the target buffer for the first use, based on the tweezers
    * 
    */
    void PMG::set_target_buffer() {
        float *temp_target = new float [width * height] {0.f};

        for (unsigned int i = 0; i < number_of_spots; i++) {
            if (should_be_used[i] > 0)
                temp_target[x_values[i] + y_values[i] * height] = sqrt(amplitude_values[i]);
        }

        // Copy to the buffer.
        cl_event wait_list[5] = {
            gpu_handler->event_list[INDEX_WRITE_X_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_Y_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_AMPLITUDE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_TWEEZER_PHASE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_SHOULD_BE_USED_BUFFER]
        };
        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->target_buffer, CL_TRUE, 0, sizeof(float) * width * height, temp_target, 5, wait_list, &gpu_handler->event_list[INDEX_WRITE_TARGET_BUFFER]
        );
        
        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error loading target buffer. OpenCL error code: " << gpu_handler->result << std::endl;
            system("pause");
            exit(1);
        }

        delete temp_target;
    }

    /* 
    * PMG::copy_mask_to_buffer()
    *
    * Formats the 8-bit mask to a float in [-pi, pi] and copies it to the phase buffer.
    */
    bool PMG::copy_mask_to_buffer() {
        float *mask = new float[width * height];
        for (auto i = 0; i < width * height; i++)
        {
            mask[i] = (unsigned) mask_8bit[i] / 255. * 2. * CL_M_PI - CL_M_PI;
        }

        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->phase_buffer, CL_TRUE, 0, sizeof(float) * width * height, mask, 0, NULL, &gpu_handler->event_list[INDEX_WRITE_PHASE_BUFFER]
        );
        
        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error loading mask from buffer. OpenCL error code: " << gpu_handler->result << std::endl;
            delete mask;
            return false;
        }

        delete mask;
        return true;
    }

    /* 
    * PMG::copy_corrections_to_buffer()
    *
    * Formats the 8-bit mask to a float in [-pi, pi] and copies it to the phase buffer.
    */
    bool PMG::copy_corrections_to_buffer() {
        float *mask = new float[width * height];
        for (auto i = 0; i < width * height; i++)
        {
            mask[i] = (unsigned) correction_mask_8bit[i] / 255. * 2. * CL_M_PI - CL_M_PI;
        }

        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->corrections_buffer, CL_TRUE, 0, sizeof(float) * width * height, mask, 0, NULL, &gpu_handler->event_list[INDEX_WRITE_PHASE_BUFFER]
        );
        
        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error loading mask from buffer. OpenCL error code: " << gpu_handler->result << std::endl;
            delete mask;
            return false;
        }

        delete mask;
        return true;
    }

    /*
    * PMG::replace_spot()
    *
    * Replaces the values of a spot at given index. In this case
    * the function takes a &Spot pointer.
    *
    * Parameters:
    * index : unsigned int
    *      Index in the spot_list to be changed.
    * spot : &Spot pointer
    *      Pointer to a replacement spot in the host memory.
    */
    bool PMG::replace_spot(
        unsigned int index,
        Spot &spot)
    {
        if (number_of_spots <= index)
        {
            return false;
        }
        x_values[index] = spot.x;
        y_values[index] = spot.y;
        amplitude_values[index] = spot.a;
        psi_values[index] = spot.psi;
        should_be_used[index] = (unsigned int) spot.used;
        return true;
    }

    /*
    * PMG::replace_spot()
    *
    * Replaces the values of a spot at given index. If index is
    * too high for the current spot_list, it will return false.
    * This overloaded function takes raw values instead of a Spot.
    *
    * Parameters:
    * index : unsigned int
    *      Index in the spot_list to be changed.
    * spot : &Spot pointer
    *      Pointer to a replacement spot in the host memory.
    */
    bool PMG::replace_spot(
        unsigned int index,
        unsigned int x,
        unsigned int y,
        float amp,
        float psi,
        unsigned int used)
    {
        if (number_of_spots <= index)
        {
            return false;
        }
        x_values[index] = x;
        y_values[index] = y;
        amplitude_values[index] = amp;
        psi_values[index] = psi;
        should_be_used[index] = used;
        return true;
    }

    /*
    * PMG::build_kernels()
    *
    * Builds the kernels in the GPU Handler. This is essential to use
    * when you have updated the number of spots (i.e. dimension of the
    * buffers).
    */
    void PMG::build_kernels() {
        gpu_handler->build_kernels();
    }

    /*
    * PMG::load_start_geometry()
    *
    * Copies the provided PMG::Spot structs into the vectors in the 
    * class memory for the start geometry. Also sets the currently active
    * values.
    * 
    * Parameters:
    * spots : std::vector<PMG::Spot>
    *   The vector containing structs with the spot data. 
    */
    void PMG::load_start_geometry(
        std::vector<Spot> spots
    ) {
        number_of_spots = 0U;
        x_values.clear();
        y_values.clear();
        amplitude_values.clear();
        psi_values.clear();

        start_x_values.clear();
        start_y_values.clear();
        start_amplitude_values.clear();
        start_phase_values.clear();
        should_be_used.clear();

        for (auto spot : spots) {
            x_values.emplace_back(spot.x);
            y_values.emplace_back(spot.y);
            amplitude_values.emplace_back(spot.a);
            psi_values.emplace_back(spot.psi);

            start_x_values.emplace_back(spot.x);
            start_y_values.emplace_back(spot.y);
            start_amplitude_values.emplace_back(spot.a);
            start_phase_values.emplace_back(spot.psi);
            should_be_used.emplace_back(spot.used);
            number_of_spots++;
        }
    }

    /*
    * PMG::load_end_geometry()
    *
    * Copies the provided PMG::Spot structs into the vectors in the 
    * class memory for the end geometry.
    * 
    * Parameters:
    * spots : std::vector<PMG::Spot>
    *   The vector containing structs with the spot data. 
    */
    void PMG::load_end_geometry(
        std::vector<Spot> spots
    ) {
        end_x_values.clear();
        end_y_values.clear();
        end_amplitude_values.clear();
        end_phase_values.clear();

        for (auto spot : spots) {
            end_x_values.emplace_back(spot.x);
            end_y_values.emplace_back(spot.y);
            end_amplitude_values.emplace_back(spot.a);
            end_phase_values.emplace_back(spot.psi);
        }
    }

    /*
    * PMG::copy_spots_to_buffers()
    *
    * Copies the spots to the GPU buffers in the GPU Handler.
    *
    * Parameters:
    * block: bool
    *   Whether the call should block the command queue of OpenCL.
    */
    bool PMG::copy_spots_to_buffers(bool block)
    {

        // Check if the buffers are of equal size.
        if (gpu_handler->buffer_size_spots < number_of_spots)
        {
            gpu_handler->initiate_spot_buffers_with_length(number_of_spots);
        }

        if (block) {
            cl_event wait_list[5] = {
                    gpu_handler->event_list[INDEX_WRITE_X_BUFFER],
                    gpu_handler->event_list[INDEX_WRITE_Y_BUFFER],
                    gpu_handler->event_list[INDEX_WRITE_AMPLITUDE_BUFFER],
                    gpu_handler->event_list[INDEX_WRITE_TWEEZER_PHASE_BUFFER],
                    gpu_handler->event_list[INDEX_WRITE_SHOULD_BE_USED_BUFFER]
                };
            clEnqueueBarrierWithWaitList(*queue_ptr, 5, wait_list, &gpu_handler->event_list[INDEX_BARRIER]);
        }
        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->x_buffer, CL_FALSE, 0, sizeof(unsigned int) * x_values.size(), x_values.data(), 0, NULL, &gpu_handler->event_list[INDEX_WRITE_X_BUFFER]
        );
        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->y_buffer, CL_FALSE, 0, sizeof(unsigned int) * y_values.size(), y_values.data(), 0, NULL, &gpu_handler->event_list[INDEX_WRITE_Y_BUFFER]
        );
        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->target_amplitude_buffer, CL_FALSE, 0, sizeof(float) * amplitude_values.size(), amplitude_values.data(), 0, NULL, &gpu_handler->event_list[INDEX_WRITE_AMPLITUDE_BUFFER]
        );
        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->target_phase_buffer, CL_FALSE, 0, sizeof(float) * psi_values.size(), psi_values.data(), 0, NULL, &gpu_handler->event_list[INDEX_WRITE_TWEEZER_PHASE_BUFFER]
        );
        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->should_be_used_buffer, CL_FALSE, 0, sizeof(unsigned int) * should_be_used.size(), should_be_used.data(), 0, NULL, &gpu_handler->event_list[INDEX_WRITE_SHOULD_BE_USED_BUFFER]
        );
        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->calculated_amplitude_buffer, CL_FALSE, 0, sizeof(unsigned int) * amplitude_values.size(), amplitude_values.data(), 0, NULL, &gpu_handler->event_list[INDEX_WRITE_CALCULATED_AMPLITUDE_BUFFER]
        );
        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Caught error in writing spots buffer. Open CL code: " << gpu_handler->result << std::endl;
            system("pause");
            exit(1);
        }

        return true;
    }

    /* 
     * PMG::load_random_phase()
     *
     * Creates a random phase pattern of dimensions `width` * `height`
     * and copies this to the `phase_buffer` of the GPU Handler.
     */
    bool PMG::load_random_phase() {
        float *mask = new float[width * height];
        for (unsigned int i = 0; i < width * height; i++) {
            mask[i] = (float) rand() / RAND_MAX * 2 * CL_M_PI - CL_M_PI;
        }

        gpu_handler->result |= clEnqueueWriteBuffer(
            *queue_ptr, gpu_handler->phase_buffer, CL_TRUE, 0, sizeof(float) * width * height, mask, 0, NULL, &gpu_handler->event_list[INDEX_WRITE_PHASE_BUFFER]
        );
        
        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error loading mask from buffer. OpenCL error code: " << gpu_handler->result << std::endl;
            delete mask;
            return false;
        }

        delete mask;
        return true;
    }

    /*
    * PMG::load_mask_from_file()
    *
    * Tries to open the file at the `path_to_file` and load
    * an 8-bit grayscale BMP. This is then loaded to the
    * phase_buffer of the GPU_Handler.
    *
    * Parameters:
    * path_to_file : std::string
    *      Path to the 8-bit grayscale BMP image.
    * do_shift : boolean
    *      Whether we need to do a shift after loading.
    *  
    * Returns:
    * success : bool
    *
    */
    bool PMG::load_mask_from_file(
        std::string path_to_file,
        bool do_shift)
    {
        CBitmap mask_bmp;

        if (!mask_bmp.Load(path_to_file.c_str())) {
            std::cout << "Error loading file at " << path_to_file << std::endl;
            return false;
        }

        mask_bmp.GetBits(mask_8bit);
        if (do_shift) {
            fftshift_mask();
        }

        return copy_mask_to_buffer();
    }

    /*
    * PMG::get_mask_from_buffer()
    *
    * Reads the buffer from the GPU Handler `return_buffer` and loads
    * this as an 8-bit mask in the `mask_8bit`.    * 
    */
    bool PMG::get_mask_from_buffer() {
        // Somehow I need to make this blocking?
        gpu_handler->result |= clEnqueueReadBuffer(
            *queue_ptr, gpu_handler->return_buffer, CL_TRUE, 0, sizeof(uint8_t) * width * height, mask_8bit, 1, 
            &gpu_handler->event_list[INDEX_BARRIER], &gpu_handler->event_list[INDEX_GET_MASK_FROM_BUFFER]
        );

        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error loading mask from buffer. OpenCL error code: " << gpu_handler->result << std::endl;
            return false;
        }

        return true;
    }

    /*
    * PMG::get_mask_from_buffer()
    * 
    * Reads the buffer from the GPU Handler `return_buffer` and loads
    * this as an 8-bit mask in the given integer array 
    */
    bool PMG::get_mask_from_buffer(uint8_t *save_mask) {

        clWaitForEvents(1, &gpu_handler->event_list[INDEX_BARRIER]);

        // Somehow I need to make this blocking?
        gpu_handler->result |= clEnqueueReadBuffer(
            *queue_ptr, gpu_handler->return_buffer, CL_TRUE, 0, sizeof(uint8_t) * width * height, save_mask, 1, 
            &gpu_handler->event_list[INDEX_BARRIER], &gpu_handler->event_list[INDEX_GET_MASK_FROM_BUFFER]
        );

        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error loading mask from buffer. OpenCL error code: " << gpu_handler->result << std::endl;
            return false;
        }

        return true;
    }

    /*
    * PMG::save_mask_to_file()
    *
    * Saves the `mask_8bit` to set `path_to_file`.
    */
    bool PMG::save_mask_to_file(
        std::string path_to_file
    ){
        // Wait to get the mask from the buffer.
        clWaitForEvents(1, &gpu_handler->event_list[INDEX_GET_MASK_FROM_BUFFER]);

        CBitmap mask_bmp;
        mask_bmp.SetBits(mask_8bit, width, height);
        mask_bmp.Save(path_to_file.c_str());
        return true;
    }

    /*
    * PMG::fftshift_mask()
    *
    * FFT shifts the mask_8bit in 2D.
    */
    void PMG::fftshift_mask() {
        uint8_t *buf = new uint8_t[width * height];
        size_t index;
        size_t shifted_index;

        for (auto x=0; x<width; x++) {
            for (auto y=0; y<height; y++) {
                index = y * width + x;
                shifted_index = ((y + height / 2) % height) * width + (x + width / 2) % width;
                buf[shifted_index] = mask_8bit[index];
            }
        }

        for (auto i=0; i<width*height; i++) {
            mask_8bit[i] = buf[i];
        }
        delete buf;
    }

    /*
    * PMG::set_ending_phase_from_file()
    *
    * Loads the .bmp mask at `path_to_file` and converts it to
    * floats. This float array is fed to the GPU Handler and serves
    * as an initial pattern, together with illumination set by
    * `set_illumation()` (typically uniform). Then the GPU Handler
    * propagates this and the phases of the tweezer spots at the 
    * focal plane are recorded and stored to `end_phase_values`.
    */
    bool PMG::set_ending_phase_from_file(
        std::string path_to_file
    ) {
        if (number_of_spots == 0) {
            std::cout << "Please set spots first." << std::endl;
            return false;
        }

        cl_event wait_list_preparation[2] = {
            gpu_handler->event_list[INDEX_WRITE_INPUT_LIGHT_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_PHASE_BUFFER]
        };

        clSetUserEventStatus(gpu_handler->event_list[INDEX_BARRIER], CL_COMPLETE);


        for (auto i=0; i<end_x_values.size(); i++) {
            x_values[i] = end_x_values[i];
            y_values[i] = end_y_values[i];
        }
        copy_spots_to_buffers();

        // Try to load the phase from file
        if (!load_mask_from_file(path_to_file.c_str())) {
            return false;
        }
        set_illumination();

        clEnqueueBarrierWithWaitList(*queue_ptr, 2, wait_list_preparation, &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_knl_assemble_field_SLM_2d();
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM_2D], &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_FFT();
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_FFT], &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_knl_get_phase_tweezers();

        gpu_handler->result |= clEnqueueReadBuffer(
            *queue_ptr, gpu_handler->target_phase_buffer, CL_TRUE, 0, sizeof(float) * end_phase_values.size(), end_phase_values.data(), 1, 
            &gpu_handler->event_list[INDEX_KNL_GET_PHASE_TWEEZERS], &gpu_handler->event_list[INDEX_READ_TWEEZER_PHASE_BUFFER]
        );
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_READ_TWEEZER_PHASE_BUFFER], &gpu_handler->event_list[INDEX_BARRIER]);

        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error loading end phase values from file. OpenCL error code: " << gpu_handler->result << std::endl;
            system("pause");
            exit(1);
        }
        reset();
        return true;
    }

    /*
    * PMG::set_starting_phase_from_file()
    *
    * Loads the .bmp mask at `path_to_file` and converts it to
    * floats. This float array is fed to the GPU Handler and serves
    * as an initial pattern, together with illumination set by
    * `set_illumation()` (typically uniform). Then the GPU Handler
    * propagates this and the phases of the tweezer spots at the 
    * focal plane are recorded and stored to `start_phase_values`.
    */
    bool PMG::set_starting_phase_from_file(
        std::string path_to_file
    ) {
        if (number_of_spots == 0) {
            std::cout << "Please set spots first." << std::endl;
            return false;
        }

        cl_event wait_list_preparation[2] = {
            gpu_handler->event_list[INDEX_WRITE_INPUT_LIGHT_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_PHASE_BUFFER]
        };

        clSetUserEventStatus(gpu_handler->event_list[INDEX_BARRIER], CL_COMPLETE);

        // Try to load the phase from file
        if (!load_mask_from_file(path_to_file.c_str(), false)) {
            return false;
        }
        set_illumination();

        clEnqueueBarrierWithWaitList(*queue_ptr, 2, wait_list_preparation, &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_knl_assemble_field_SLM();
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM_2D], &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_FFT();
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_FFT], &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_knl_get_phase_tweezers();

        gpu_handler->result |= clEnqueueReadBuffer(
            *queue_ptr, gpu_handler->target_phase_buffer, CL_TRUE, 0, sizeof(float) * start_phase_values.size(), start_phase_values.data(), 1, 
            &gpu_handler->event_list[INDEX_KNL_GET_PHASE_TWEEZERS], &gpu_handler->event_list[INDEX_READ_TWEEZER_PHASE_BUFFER]
        );
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_READ_TWEEZER_PHASE_BUFFER], &gpu_handler->event_list[INDEX_BARRIER]);
        

        for (auto i=0; i<number_of_spots; i++) {
            psi_values[i] = start_phase_values[i];
        }

        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error loading start phase values from file. OpenCL error code: " << gpu_handler->result << std::endl;
            system("pause");
            exit(1);
        }
        return true;
    }

    /*
    * PMG::set_starting_phase_from_file()
    *
    * Loads the .bmp mask at `path_to_file` and converts it to
    * floats. This float array is fed to the GPU Handler and serves
    * as an initial pattern, together with illumination set by
    * `set_illumation()` (typically uniform). Then the GPU Handler
    * propagates this and the phases of the tweezer spots at the 
    * focal plane are recorded and stored to `start_phase_values`.
    */
    bool PMG::set_starting_phase_and_amps_from_file(
        std::string path_to_file
    ) {
        if (number_of_spots == 0) {
            std::cout << "Please set spots first." << std::endl;
            return false;
        }

        cl_event wait_list_preparation[2] = {
            gpu_handler->event_list[INDEX_WRITE_INPUT_LIGHT_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_PHASE_BUFFER]
        };

        clSetUserEventStatus(gpu_handler->event_list[INDEX_BARRIER], CL_COMPLETE);

        // Try to load the phase from file
        if (!load_mask_from_file(path_to_file.c_str())) {
            return false;
        }
        set_illumination();

        clEnqueueBarrierWithWaitList(*queue_ptr, 2, wait_list_preparation, &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_knl_assemble_field_SLM();
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_ASSEMBLE_FIELD_SLM_2D], &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_FFT();
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_FFT], &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_knl_get_phase_tweezers();
        clFinish(*queue_ptr);
        
        
        gpu_handler->result |= clEnqueueReadBuffer(
            *queue_ptr, gpu_handler->target_phase_buffer, CL_TRUE, 0, sizeof(float) * start_phase_values.size(), start_phase_values.data(), 0, NULL, &gpu_handler->event_list[INDEX_READ_TWEEZER_PHASE_BUFFER]
        );
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_READ_TWEEZER_PHASE_BUFFER], &gpu_handler->event_list[INDEX_BARRIER]);

        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error lsadoading start phase values from file. OpenCL error code: " << gpu_handler->result << std::endl;
            system("pause");
            exit(1);
        }
        
        gpu_handler->run_knl_calculate_amplitudes();
        clFinish(*queue_ptr);

        gpu_handler->result |= clEnqueueReadBuffer(
            *queue_ptr, gpu_handler->calculated_amplitude_buffer, CL_TRUE, 0, sizeof(float) * start_amplitude_values.size(), start_amplitude_values.data(), 0, NULL, &gpu_handler->event_list[INDEX_READ_AMPLITUDE_BUFFER]
        );
        //clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_READ_AMPLITUDE_BUFFER], &gpu_handler->event_list[INDEX_BARRIER]);
        

        for (auto i=0; i<number_of_spots; i++) {
            psi_values[i] = start_phase_values[i];
        }

        if (gpu_handler->result != CL_SUCCESS) {
            std::cout << "Error loading start phase values from file. OpenCL error code: " << gpu_handler->result << std::endl;
            system("pause");
            exit(1);
        }
        return true;
    }

    /*
    * PMG::reset()
    *
    * A reset call to put all the values back to their start values.
    */
    void PMG::reset() {
        for (auto i=0; i<number_of_spots; ++i) {
            x_values[i] = start_x_values[i];
            y_values[i] = start_y_values[i];
            amplitude_values[i] = start_amplitude_values[i];
            psi_values[i] = start_phase_values[i];
            should_be_used[i] = true;
        }
        copy_spots_to_buffers(true);
    }

    /* 
    * PMG::calculate_sorting_shortcut()
    *
    * Internal function used for quickly computing the next hologram and displaying it 
    * on the SLM. Naming should be improved...
    * 
    * Only works in the case that the start and end geometry are the same and a 
    * mapping consists only of shifting indices. For more complex moves, consider
    * `PMG::calculate_sorting_arbitrary()`.
    * */
    void PMG::calculate_sorting_shortcut() {
        cl_event wait_list_spots_transfer[5] = {
            gpu_handler->event_list[INDEX_WRITE_X_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_Y_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_AMPLITUDE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_TWEEZER_PHASE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_SHOULD_BE_USED_BUFFER]
        };

        cl_event wait_list_spots_calculation[2] = {
            gpu_handler->event_list[INDEX_KNL_RESET_FIELD],
            gpu_handler->event_list[INDEX_HOST_GET_POSITIONS]
        };

        bool updated;
        for (auto k = 1; k <= sorting_machine->number_of_steps + sorting_machine->number_turnoff_frames; k++) {
            gpu_handler->run_knl_reset_field(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_RESET_FIELD], &gpu_handler->event_list[INDEX_BARRIER]);
           
            clSetUserEventStatus(gpu_handler->event_list[INDEX_HOST_GET_POSITIONS], CL_SUBMITTED);
            updated = sorting_machine->get_next_positions(
                k, x_values, y_values, amplitude_values, psi_values, should_be_used, 
                start_x_values, start_y_values, start_amplitude_values, start_phase_values, 
                end_phase_values, std::max(width, height)
            );

            if (!updated) { break; }
            clSetUserEventStatus(gpu_handler->event_list[INDEX_HOST_GET_POSITIONS], CL_COMPLETE);

            clEnqueueBarrierWithWaitList(*queue_ptr, 2, wait_list_spots_calculation, &gpu_handler->event_list[INDEX_BARRIER]);
            copy_spots_to_buffers(true);            
            clEnqueueBarrierWithWaitList(*queue_ptr, 5, wait_list_spots_transfer, &gpu_handler->event_list[INDEX_BARRIER]);

            gpu_handler->run_knl_apply_constraints_GS_locked(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_CONSTRAINTS_GS_LOCKED], &gpu_handler->event_list[INDEX_BARRIER]);

            gpu_handler->run_iFFT(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_IFFT], &gpu_handler->event_list[INDEX_BARRIER]);

            gpu_handler->run_knl_get_phase_with_corrections(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_GET_PHASE], &gpu_handler->event_list[INDEX_BARRIER]);

            get_mask_from_buffer(mask_to_display_8bit);
            display_mask_on_SLM();

            if (k == sorting_machine->number_of_steps + sorting_machine->number_turnoff_frames) break;
        }

        clFinish(*queue_ptr);
    }

    /* 
    * PMG::calculate_sorting_shortcut()
    *
    * Internal function used for quickly computing the next hologram and displaying it 
    * on the SLM. Naming should be improved...
    * 
    * Only works in the case that the start and end geometry are the same and a 
    * mapping consists only of shifting indices. For more complex moves, consider
    * `PMG::calculate_sorting_arbitrary()`.
    * */
    void PMG::calculate_sorting_arbitrary() {
        cl_event wait_list_spots_transfer[5] = {
            gpu_handler->event_list[INDEX_WRITE_X_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_Y_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_AMPLITUDE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_TWEEZER_PHASE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_SHOULD_BE_USED_BUFFER]
        };

        cl_event wait_list_spots_calculation[2] = {
            gpu_handler->event_list[INDEX_KNL_RESET_FIELD],
            gpu_handler->event_list[INDEX_HOST_GET_POSITIONS]
        };

        bool updated;
        for (auto k = 1; k <= sorting_machine->number_of_steps + sorting_machine->number_turnoff_frames; k++) {
            gpu_handler->run_knl_reset_field(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_RESET_FIELD], &gpu_handler->event_list[INDEX_BARRIER]);
           
            clSetUserEventStatus(gpu_handler->event_list[INDEX_HOST_GET_POSITIONS], CL_SUBMITTED);
            updated = sorting_machine->get_next_positions(
                k, x_values, y_values, amplitude_values, psi_values, should_be_used, 
                start_x_values, start_y_values, start_amplitude_values, start_phase_values, 
                end_x_values, end_y_values, end_amplitude_values, end_phase_values, 
                std::max(width, height)
            );

            if (!updated) { break; }
            clSetUserEventStatus(gpu_handler->event_list[INDEX_HOST_GET_POSITIONS], CL_COMPLETE);

            clEnqueueBarrierWithWaitList(*queue_ptr, 2, wait_list_spots_calculation, &gpu_handler->event_list[INDEX_BARRIER]);
            copy_spots_to_buffers(true);            
            clEnqueueBarrierWithWaitList(*queue_ptr, 5, wait_list_spots_transfer, &gpu_handler->event_list[INDEX_BARRIER]);

            gpu_handler->run_knl_apply_constraints_GS_locked(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_CONSTRAINTS_GS_LOCKED], &gpu_handler->event_list[INDEX_BARRIER]);

            gpu_handler->run_iFFT(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_IFFT], &gpu_handler->event_list[INDEX_BARRIER]);

            gpu_handler->run_knl_get_phase_with_corrections(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_GET_PHASE], &gpu_handler->event_list[INDEX_BARRIER]);

            get_mask_from_buffer(mask_to_display_8bit);
            display_mask_on_SLM();

            if (k == sorting_machine->number_of_steps + sorting_machine->number_turnoff_frames) break;
        }

        clFinish(*queue_ptr);
    }

    /*
     * PMG::test()
     *
     * Tests the program timings with various sorting_method, calculation_method
     * and pathing_method. Currently options are:
     * 
     * - sorting_method:
     *      0.  Auto-detect
     *      1.  Hungarian Algorithm
     *      2.  Compression Algorithm
     *      3.  Partitioned Compression Algorithm
     * - calculation_method:
     *      0.  None (just check sorting calculation)
     *      1.  Sorting GS
     * - pathing_method:
     *      0.  Direct linear interpolation
     * 
     * Saves different timings and steps in the timing_list TestResult. The indices are:
     *  [0]:    (double) sorting calculation time
     *  [1]:    (double) pattern calculation loop time
     *  [2]:    (double) total run time function
     *  [3]:    (unsigned int) number of patterns
     *  [4]:    (unsigned int) sorting_method
     *  [5]:    (unsigned int) calculation_method
     *  [6]:    (unsigned int) pathing_method
     */
    TestResult PMG::test(
        std::string_view target_str,
        std::string_view loaded_str,
        unsigned int sorting_method,
        unsigned int calculation_method,
        unsigned int pathing_method
    ) {
        // Initiate time points for the test.
        std::chrono::steady_clock::time_point time_end_sorting;
        std::chrono::steady_clock::time_point time_end_calculation;
        std::chrono::steady_clock::time_point time_start_test = std::chrono::steady_clock::now();

        // Set the right methods for the SortingMachine
        sorting_machine->sorting_method = sorting_method;
        sorting_machine->pathing_method = pathing_method;

        // Parse strings from the AndorServer
        sorting_machine->parse_target_str(target_str);
        sorting_machine->parse_loaded_str(loaded_str);

        // Begin with sorting.
        sorting_machine->sort(start_x_values, start_y_values, std::max(width, height));
        //std::cout << "Done sorting" << std::endl;

        // Calculate number of steps (i.e. patterns) needed.
        sorting_machine->calculate_maximal_movement(start_x_values, start_y_values, std::max(width, height));
        sorting_machine->number_turnoff_frames = 1;

        time_end_sorting = std::chrono::steady_clock::now();

        switch (calculation_method) {
            case 0:
                //std::cout << "No calculation method selected" << std::endl;
                break;
            case 1:
                calculate_sorting_shortcut();
                break;
            default:
                calculate_sorting_shortcut();
                break;
        }

        time_end_calculation = std::chrono::steady_clock::now();

        // Calculate results.
        TestResult test_result;
        test_result.duration_sorting = (double) std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_sorting - time_start_test).count() / 1000000.;
        test_result.duration_calculation = (double) std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_calculation - time_end_sorting).count() / 1000000.;
        test_result.duration_total = (double) std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_calculation - time_start_test).count() / 1000000.;
        test_result.number_of_patterns = sorting_machine->number_turnoff_frames + sorting_machine->number_of_steps;
        test_result.sorting_method = sorting_machine->sorting_method;
        test_result.calculation_method = calculation_method;
        test_result.pathing_method = pathing_method;

        return test_result;  
    }

    /*
     * PMG::test_shift()
     *
     * Shifts the tweezer positions by 10 steps and computes a hologram for each step.
     * Then saves the average time per substep in a TestResult.
     * 
     * Saves different timings and steps in the timing_list TestResult. The indices are:
     *  [0]:    (double) sorting calculation time
     *  [1]:    (double) hologram calculation loop time
     *  [2]:    (double) total run time function
     *  [3]:    (double) duration until received OK from SLM for display
     *  [4]:    (double) duration of transfer from GPU to CPU memory
     *  [5]:    (unsigned int) number of tweezers
     *  [6]:    (unsigned int) number of patterns
     *  [7]:    (unsigned int) sorting_method
     *  [8]:    (unsigned int) calculation_method
     *  [9]:    (unsigned int) pathing_method
     */
    TestResult PMG::test_shift_w_slm() {

        cl_event wait_list_spots_transfer[6] = {
            gpu_handler->event_list[INDEX_WRITE_X_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_Y_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_AMPLITUDE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_TWEEZER_PHASE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_SHOULD_BE_USED_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_CALCULATED_AMPLITUDE_BUFFER]
        };

        // Initiate time points for the test.
        std::chrono::steady_clock::time_point time_end_sorting;
        std::chrono::steady_clock::time_point time_end_calculation;
        std::chrono::steady_clock::time_point time_end_transfer;
        std::chrono::steady_clock::time_point time_end_display;
        std::chrono::steady_clock::time_point time_start_test;      

        unsigned int number_repetitions = 10U;
        double sum_duration_sorting = 0.;
        double sum_duration_calculation = 0.;
        double sum_duration_transfer = 0.;
        double sum_duration_display = 0.;
        double sum_duration_time_total = 0.;
        for (auto i=0; i<number_repetitions; i++ ) {
            time_start_test = std::chrono::steady_clock::now();

            shift_pattern(0, i);
            time_end_sorting = std::chrono::steady_clock::now();

            reset_field();
            copy_spots_to_buffers(true);            
            clEnqueueBarrierWithWaitList(*queue_ptr, 6, wait_list_spots_transfer, &gpu_handler->event_list[INDEX_BARRIER]);
            gpu_handler->run_knl_apply_constraints_GS_locked(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_CONSTRAINTS_GS_LOCKED], &gpu_handler->event_list[INDEX_BARRIER]);
            gpu_handler->run_iFFT(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_IFFT], &gpu_handler->event_list[INDEX_BARRIER]);
            gpu_handler->run_knl_get_phase_with_corrections(false);
            clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_GET_PHASE], &gpu_handler->event_list[INDEX_BARRIER]);
            time_end_calculation = std::chrono::steady_clock::now();

            get_mask_from_buffer(mask_to_display_8bit);
            time_end_transfer = std::chrono::steady_clock::now();           

            display_mask_on_SLM();
            time_end_display = std::chrono::steady_clock::now();

            sum_duration_sorting += (double) std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_sorting - time_start_test).count() / 1000000.;
            sum_duration_calculation += (double) std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_calculation - time_end_sorting).count() / 1000000.;
            sum_duration_transfer += (double) std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_transfer - time_end_calculation).count() / 1000000.;
            sum_duration_display += (double) std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_display - time_end_transfer).count() / 1000000.;
            sum_duration_time_total += (double) std::chrono::duration_cast<std::chrono::nanoseconds>(time_end_display - time_start_test).count() / 1000000.;
            
            reset();
        }

        // Calculate results.
        TestResult test_result;
        test_result.duration_sorting = sum_duration_sorting / number_repetitions;
        test_result.duration_calculation = sum_duration_calculation / number_repetitions;
        test_result.duration_transfer = sum_duration_transfer / number_repetitions;
        test_result.duration_display = sum_duration_display / number_repetitions;
        test_result.duration_total = sum_duration_time_total / number_repetitions;
        test_result.number_of_traps = number_of_spots;
        test_result.number_of_patterns = number_repetitions;
        test_result.sorting_method = 0;
        test_result.calculation_method = 0;
        test_result.pathing_method = 0;

        return test_result;  
    }

    /*
    * void PMG::run_GSW()
    *
    * Uses the internally stored spot values to run one loop of the 
    * GSW algorithm. Stores the phase mask in the mask_to_display_8bit
    * array.
    * */
    void PMG::run_GSW(int token) {
        cl_event wait_list_spots_transfer[6] = {
            gpu_handler->event_list[INDEX_WRITE_X_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_Y_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_AMPLITUDE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_TWEEZER_PHASE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_SHOULD_BE_USED_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_CALCULATED_AMPLITUDE_BUFFER]
        };

        copy_spots_to_buffers(true);            
        clEnqueueBarrierWithWaitList(*queue_ptr, 6, wait_list_spots_transfer, &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_knl_apply_constraints_GS_locked(false);
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_CONSTRAINTS_GS_LOCKED], &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_iFFT(false);
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_IFFT], &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_knl_get_phase_with_corrections(false);
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_GET_PHASE], &gpu_handler->event_list[INDEX_BARRIER]);
        get_mask_from_buffer(mask_to_display_8bit);
    }

    /*
    * void PMG::run_GSW_with_roll()
    *
    * Uses the internally stored spot values to run one loop of the 
    * GSW algorithm, with a roll of the mask in the end. Stores the 
    * phase mask in the mask_to_display_8bit array.
    * */
    void PMG::run_GSW_with_roll(int offset_x, int offset_y, int token) {
        cl_event wait_list_spots_transfer[6] = {
            gpu_handler->event_list[INDEX_WRITE_X_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_Y_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_AMPLITUDE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_TWEEZER_PHASE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_SHOULD_BE_USED_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_CALCULATED_AMPLITUDE_BUFFER]
        };

        copy_spots_to_buffers(true);            
        clEnqueueBarrierWithWaitList(*queue_ptr, 6, wait_list_spots_transfer, &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_knl_apply_constraints_GS_locked(false);
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_CONSTRAINTS_GS_LOCKED], &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_iFFT(false);
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_APPLY_IFFT], &gpu_handler->event_list[INDEX_BARRIER]);
        gpu_handler->run_knl_get_phase_with_corrections_with_roll(offset_x, offset_y, false);
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_GET_PHASE], &gpu_handler->event_list[INDEX_BARRIER]);
        get_mask_from_buffer(mask_to_display_8bit);
    }

    /*
    * void PMG::reset_field()
    *
    * Runs the reset_field kernel to make the values of the field at the
    * current spot locations zero. A quicker ways of removing spots in the
    * calculation than reinitializing a buffer.
    */
    void PMG::reset_field() {
        cl_event wait_list_spots_transfer[6] = {
            gpu_handler->event_list[INDEX_WRITE_X_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_Y_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_AMPLITUDE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_TWEEZER_PHASE_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_SHOULD_BE_USED_BUFFER],
            gpu_handler->event_list[INDEX_WRITE_CALCULATED_AMPLITUDE_BUFFER]
        };

        gpu_handler->run_knl_reset_field(false);
        clEnqueueBarrierWithWaitList(*queue_ptr, 1, &gpu_handler->event_list[INDEX_KNL_RESET_FIELD], &gpu_handler->event_list[INDEX_BARRIER]);
    }

    /*
    * void PMG::sort()
    *
    * Given a loading_str and a target_str, this function will produce
    * and display phasemasks to sort into the target string. It assumes
    * the indices in the target_str correspond to the same traps as in 
    * the loading_str. 
    * 
    * The user has options for the different methods used: sorting_method, 
    * calculation_method and pathing_method. Currently options are:
    * 
    * - sorting_method:
    *      0.  Auto-detect
    *      1.  Hungarian Algorithm
    *      2.  Compression Algorithm
    *      3.  Partitioned Compression Algorithm
    * - calculation_method:
    *      0.  None (just check sorting calculation)
    *      1.  Sorting GS
    * - pathing_method:
    *      0.  Direct linear interpolation
    */
    void PMG::sort(
        std::string_view target_str,
        std::string_view loaded_str,
        unsigned int sorting_method,
        unsigned int calculation_method,
        unsigned int pathing_method
    ) {
        // Set the right methods for the SortingMachine
        sorting_machine->sorting_method = sorting_method;
        sorting_machine->pathing_method = pathing_method;

        // Parse strings from the AndorServer
        sorting_machine->parse_target_str(target_str);
        sorting_machine->parse_loaded_str(loaded_str);

        if (sorting_machine->number_loaded < sorting_machine->number_target) {
            //std::cout << "Not enough atoms loaded" << std::endl;
            sorting_machine->number_of_steps = 1U;
            sorting_machine->number_turnoff_frames = 0U;
            load_mask_from_file("C:\\Users\\FastSLM\\SLM\\WFC Files\\blank.bmp");
            for (auto i=0; i < width * height; i++) {
                mask_to_display_8bit[i] = mask_8bit[i];
            }
            //add_corrections_to_mask();
            display_mask_on_SLM();
            return;
        }

        // Begin with sorting.
        sorting_machine->sort(start_x_values, start_y_values, std::max(width, height));
        //std::cout << "Done sorting" << std::endl;

        

        // Calculate number of steps (i.e. patterns) needed.
        sorting_machine->calculate_maximal_movement(start_x_values, start_y_values, std::max(width, height));
        sorting_machine->number_turnoff_frames = 1U;

        switch (calculation_method) {
            case 0:
                //std::cout << "No calculation method selected" << std::endl;
                break;
            case 1:
                calculate_sorting_shortcut();
                break;
            default:
                calculate_sorting_shortcut();
                break;
        }
    }

    /*
    * void PMG::sort_arbitrary()
    *
    * Given a loading_str and a target_str, this function will produce
    * and display phasemasks to sort into the target string. It assumes
    * the indices in the target_str correspond to the same traps as in 
    * the loading_str. 
    * 
    * The user has options for the different methods used: sorting_method, 
    * calculation_method and pathing_method. Currently options are:
    * 
    * - sorting_method:
    *      0.  Auto-detect
    *      1.  Hungarian Algorithm
    *      2.  Compression Algorithm
    *      3.  Partitioned Compression Algorithm
    * - calculation_method:
    *      0.  None (just check sorting calculation)
    *      1.  Sorting GS
    * - pathing_method:
    *      0.  Direct linear interpolation
    */
    void PMG::sort_arbitrary(
        std::string_view target_str,
        std::string_view loaded_str,
        unsigned int sorting_method,
        unsigned int calculation_method,
        unsigned int pathing_method
    ) {
        // Set the right methods for the SortingMachine
        sorting_machine->sorting_method = sorting_method;
        sorting_machine->pathing_method = pathing_method;

        // Parse strings from the AndorServer
        sorting_machine->parse_target_str(target_str);
        sorting_machine->parse_loaded_str(loaded_str);

        if (sorting_machine->number_loaded < sorting_machine->number_target) {
            //std::cout << "Not enough atoms loaded" << std::endl;
            sorting_machine->number_of_steps = 1U;
            sorting_machine->number_turnoff_frames = 0U;
            load_mask_from_file("C:\\Users\\FastSLM\\SLM\\WFC Files\\blank.bmp");
            for (auto i=0; i < width * height; i++) {
                mask_to_display_8bit[i] = mask_8bit[i];
            }
            //add_corrections_to_mask();
            display_mask_on_SLM();
            return;
        }

        // Begin with sorting.
        sorting_machine->sort(
            start_x_values, 
            start_y_values,
            end_x_values,
            end_y_values,
            std::max(width, height)
        );
        //std::cout << "Done sorting" << std::endl;

        // Calculate number of steps (i.e. patterns) needed.
        sorting_machine->calculate_maximal_movement(
            start_x_values, 
            start_y_values, 
            end_x_values,
            end_y_values,
            std::max(width, height));
        if (sorting_machine->number_loaded != sorting_machine->number_target) {
            sorting_machine->number_turnoff_frames = 1U;
        }
        else {
            sorting_machine->number_turnoff_frames = 0U;
        }
        
        std::cout << "Need to sort " << sorting_machine->number_of_steps << " steps " <<std::endl;

        switch (calculation_method) {
            case 0:
                //std::cout << "No calculation method selected" << std::endl;
                break;
            case 1:
                calculate_sorting_arbitrary();
                break;
            default:
                calculate_sorting_arbitrary();
                break;
        }
    }

    /*
    * void PMG::sort_arbitrary()
    *
    * Given a loading_str and a target_str, this function will produce
    * and display phasemasks to sort into the target string. It assumes
    * the indices in the target_str correspond to the same traps as in 
    * the loading_str. 
    * 
    * The user has options for the different methods used: sorting_method, 
    * calculation_method and pathing_method. Currently options are:
    * 
    * - sorting_method:
    *      0.  Auto-detect
    *      1.  Hungarian Algorithm
    *      2.  Compression Algorithm
    *      3.  Partitioned Compression Algorithm
    * - calculation_method:
    *      0.  None (just check sorting calculation)
    *      1.  Sorting GS
    * - pathing_method:
    *      0.  Direct linear interpolation
    */
    unsigned int PMG::sort_arbitrary_return(
        std::string_view target_str,
        std::string_view loaded_str,
        unsigned int sorting_method,
        unsigned int calculation_method,
        unsigned int pathing_method
    ) {
        // Set the right methods for the SortingMachine
        sorting_machine->sorting_method = sorting_method;
        sorting_machine->pathing_method = pathing_method;

        // Parse strings from the AndorServer
        sorting_machine->parse_target_str(target_str);
        sorting_machine->parse_loaded_str(loaded_str);

        if (sorting_machine->number_loaded < sorting_machine->number_target) {
            //std::cout << "Not enough atoms loaded" << std::endl;
            sorting_machine->number_of_steps = 1U;
            sorting_machine->number_turnoff_frames = 0U;
            load_mask_from_file("C:\\Users\\FastSLM\\SLM\\WFC Files\\blank.bmp");
            for (auto i=0; i < width * height; i++) {
                mask_to_display_8bit[i] = mask_8bit[i];
            }
            //add_corrections_to_mask();
            display_mask_on_SLM();
            return 0U;
        }

        // Begin with sorting.
        sorting_machine->sort(
            start_x_values, 
            start_y_values,
            end_x_values,
            end_y_values,
            std::max(width, height)
        );
        //std::cout << "Done sorting" << std::endl;

        // Calculate number of steps (i.e. patterns) needed.
        sorting_machine->calculate_maximal_movement(
            start_x_values, 
            start_y_values, 
            end_x_values,
            end_y_values,
            std::max(width, height));
        if (sorting_machine->number_loaded != sorting_machine->number_target) {
            sorting_machine->number_turnoff_frames = 1U;
        }
        else {
            sorting_machine->number_turnoff_frames = 0U;
        }
        
        std::cout << "Need to sort " << sorting_machine->number_of_steps << " steps " <<std::endl;

        switch (calculation_method) {
            case 0:
                //std::cout << "No calculation method selected" << std::endl;
                break;
            case 1:
                calculate_sorting_arbitrary();
                break;
            default:
                calculate_sorting_arbitrary();
                break;
        }

        return sorting_machine->number_of_steps;
    }
};