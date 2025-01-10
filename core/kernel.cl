# pragma OPENCL EXTENSION cl_khr_fp64 : enable
# define HEIGHT 1024
# define WIDTH 1024

kernel void set_amplitude_field(global float2 *field, global float *amplitude) {
    size_t index = get_global_id(0);

    if (amplitude[index] > 0.f) {
        float phase = atan2(field[index].y, field[index].x);

        field[index].x = amplitude[index] * cos(phase);
        field[index].y = amplitude[index] * sin(phase);
    }
    else {
        field[index].x = 0.f;
        field[index].y = 0.f;
    }
};

kernel void set_amplitude_field_2d(global float2 *field, global float *amplitude) {
    size_t index0 = get_global_id(0);
    size_t index1 = get_global_id(1);
    size_t grid_width = get_num_groups(0) * get_local_size(0);
    __private unsigned int pos = index1 * grid_width + index0;
    //__private unsigned int pos = index0 * get_global_size(1) + index1;

    if (amplitude[pos] > 0.f) {
        float phase = atan2(field[pos].y, field[pos].x);

        field[pos].x = amplitude[pos] * cos(phase);
        field[pos].y = amplitude[pos] * sin(phase);
    }
    //else {
    //    field[pos].x = 0.f;
    //    field[pos].y = 0.f;
    //}
};

kernel void assemble_field(global float2 *field, global float *amplitude, global float *phase) {
    size_t index = get_global_id(0);
    field[index].x = amplitude[index] * cos(phase[index]);
    field[index].y = amplitude[index] * sin(phase[index]);
};

kernel void assemble_field_2d(global float2 *field, global float *amplitude, global float *phase) {
    size_t index0 = get_global_id(0);
    size_t index1 = get_global_id(1);
    size_t grid_width = get_num_groups(1) * get_local_size(1);
    __private unsigned int pos = index0 * grid_width + index1;
    //__private unsigned int pos = index0 * get_global_size(1) + index1;

    field[pos].x = amplitude[pos] * cos(phase[pos]);
    field[pos].y = amplitude[pos] * sin(phase[pos]);
};

kernel void get_phase(global float2 *field, global uchar *phase){
    size_t index = get_global_id(0);
    phase[index] = atan2(field[index].y, field[index].x) / (2 * 3.14159265359f) * 255;
};

kernel void get_phase_2d(global float2 *field, global uchar *phase){
    size_t index0 = get_global_id(0);
    size_t index1 = get_global_id(1);
    size_t grid_width = get_num_groups(0) * get_local_size(0);
    __private unsigned int pos = index1 * grid_width + index0;
    //__private unsigned int pos = index0 * get_global_size(1) + index1;

    phase[pos] = atan2(field[pos].y, field[pos].x) / (2 * 3.14159265359f) * 255;
};

kernel void set_target(global float* amplitude, global unsigned int *x, global unsigned int *y, global unsigned int *used, global float* target_amplitude) {
    size_t index = get_global_id(0);
    unsigned int pos = y[index] * WIDTH + x[index];

    if (used[index] > 0) {
        amplitude[pos] = sqrt(target_amplitude[index]);
    }
    else {
        amplitude[pos] = 0.f;
    }
};

kernel void update_target(global float* amplitude, global unsigned int *x, global unsigned int *y, global unsigned int *used, global float* target_amplitude, global float* calculated_amplitude ) {
    size_t index = get_global_id(0);
    unsigned int pos = y[index] * WIDTH + x[index];

    if (used[index] > 0) {
        if (amplitude[pos] > 0.f) {
            amplitude[pos] *= sqrt(target_amplitude[index] / calculated_amplitude[index]);
        }
        else {
            amplitude[pos] = sqrt(target_amplitude[index]);
        }
    }
    else {
        amplitude[pos] = 0.f;
    }
};

kernel void reset_target(global float* amplitude, global unsigned int *x, global unsigned int *y) {
    size_t index = get_global_id(0);
    unsigned int pos = y[index] * WIDTH + x[index];
    amplitude[pos] = 0.f;
};

kernel void reset_field(global float2 *field, global unsigned int *x, global unsigned int *y) {
    size_t index = get_global_id(0);
    unsigned int pos = y[index] * WIDTH + x[index];
    field[pos].x = 0.f;
    field[pos].y = 0.f;
};

kernel void reset_field_2d(global float2 *field, global unsigned int *x, global unsigned int *y) {
    size_t index = get_global_id(0);
    unsigned int pos = y[index] * WIDTH + x[index];
    field[pos].x = 0.f;
    field[pos].y = 0.f;
};

kernel void set_phase_tweezers(global float2 *field, global unsigned int *x, global unsigned int *y, global unsigned int *used, global float* target_phase) {
    size_t index = get_global_id(0);
    unsigned int pos = y[index] * WIDTH + x[index];

    if (used[index] > 0) {
        
        float amplitude = sqrt(field[pos].x * field[pos].x + field[pos].y * field[pos].y);

        field[pos].x = amplitude * cos(target_phase[index]);
        field[pos].y = amplitude * sin(target_phase[index]);
    }
    // Not sure why, but this else makes the pattern sort of wiggled/drunk?
    //else {
    //    field[pos].x = 0.f;
    //    field[pos].y = 0.f;
    //}
};

kernel void get_phase_tweezers(global float2 *field, global unsigned int *x, global unsigned int *y, global float* phase) {
    int index = get_global_id(0);
    unsigned int pos = y[index] * WIDTH + x[index];
    phase[index] = atan2(field[pos].y, field[pos].x);
};

kernel void calculate_amplitudes(global float2 *field, global unsigned int *x, global unsigned int *y, global unsigned int *used, global float* calculated_amplitude) {
    int index = get_global_id(0);

    if (used[index] > 0) {
        unsigned int pos = y[index] * WIDTH + x[index];
        calculated_amplitude[index] = field[pos].x * field[pos].x + field[pos].y * field[pos].y;
    }
    else {
        calculated_amplitude[index] = 0.;
    }
};

kernel void apply_constraints_GSW_locked(
    global float2 *field,
    global unsigned int *x,
    global unsigned int *y,
    global unsigned int *used,
    global float *target_amplitude,
    global float *target_phase
) {
    size_t index = get_global_id(0);
    unsigned int pos = y[index] * WIDTH + x[index];

    if (used[index] > 0) {
        // First we update the weighted amplitude of the spot.
        float amplitude = field[pos].x * field[pos].x + field[pos].y * field[pos].y;
        if (amplitude > 0.f) {
            amplitude = sqrt(target_amplitude[index] / amplitude);
        }
        else {
            amplitude = sqrt(target_amplitude[index]);
        }

        // Now we compose the field with new amplitude and phase.
        field[pos].x = amplitude * cos(target_phase[index]);
        field[pos].y = amplitude * sin(target_phase[index]);
    }
    else {
        field[pos].x = 0.f;
        field[pos].y = 0.f;
    }
};

kernel void apply_constraints_GS_locked(
    global float2 *field,
    global unsigned int *x,
    global unsigned int *y,
    global unsigned int *used,
    global float *target_amplitude,
    global float *target_phase
) {
    size_t index = get_global_id(0);
    unsigned int pos = y[index] * WIDTH + x[index];

    if (used[index] > 0) {
        // Now we compose the field with new amplitude and phase.
        // Note the factor 1000 to get rid of spurious phases from higher orders!
        field[pos].x = 1000 * fabs(sqrt(target_amplitude[index])) * cos(target_phase[index]);
        field[pos].y = 1000 * fabs(sqrt(target_amplitude[index])) * sin(target_phase[index]);
    }
};

kernel void apply_constraints_GSW(
    global float2 *field,
    global unsigned int *x,
    global unsigned int *y,
    global unsigned int *used,
    global float *target_amplitude,
    global float *calculated_amplitude,
    const int number_of_spots
) {
    size_t index = get_global_id(0);
    unsigned int pos = y[index] * WIDTH + x[index];

    if (used[index] > 0) {
        float mean_amplitude = 0;
        for (unsigned int j=0; j < number_of_spots; j++) {
            if (used[j] > 0) {
                mean_amplitude += sqrt(calculated_amplitude[j] / target_amplitude[j]) / number_of_spots;
            }
        }

        // Now we compose the field with new amplitude and phase.
        if (calculated_amplitude[index] > 0.f) {
            field[pos].x *= mean_amplitude * sqrt(target_amplitude[index] / calculated_amplitude[index]);
            field[pos].y *= mean_amplitude * sqrt(target_amplitude[index] / calculated_amplitude[index]);
        }
    }
    else {
        field[pos].x = 0.;
        field[pos].y = 0.;
    }
};

kernel void apply_constraints_GS(
    global float2 *field,
    global unsigned int *x,
    global unsigned int *y,
    global unsigned int *used,
    global float *target_amplitude
) {
    size_t index = get_global_id(0);

    if (used[index] > 0) {
        unsigned int pos = y[index] * WIDTH + x[index];

        float phase = atan2(field[index].y, field[index].x);

        field[index].x = sqrt(target_amplitude[index]) * cos(phase);
        field[index].y = sqrt(target_amplitude[index]) * sin(phase);
    }
};

kernel void get_phase_with_corrections_with_roll(
    global float2 *field, 
    global float *corrections, 
    global uchar *mask,
    const int dx,
    const int dy
) {
    size_t index = get_global_id(0);
    unsigned int y = index / WIDTH;
    unsigned int x = index % WIDTH;

    size_t new_index = ((y + HEIGHT / 2 + dx) % HEIGHT) * WIDTH + (x + WIDTH / 2 + dy) % WIDTH;

    // The additional pi is unclear to me, but necessary to match Python code.
    mask[index] = convert_uchar(
        fmod(
            (
                atan2(
                    field[new_index].y, 
                    field[new_index].x
                ) + corrections[index] + 2 * 3.14159265359f
            ) / (2 * 3.14159265359f) * 255, 
            256
        )
    );
};


kernel void get_phase_with_corrections(
    global float2 *field, 
    global float *corrections, 
    global uchar *mask
) {
    size_t index = get_global_id(0);
    unsigned int y = index / WIDTH;
    unsigned int x = index % WIDTH;
    //unsigned int yshift = (y + HEIGHT / 2) % HEIGHT;
    //unsigned int xshift = (x + WIDTH / 2) % WIDTH;
    size_t new_index = (y + HEIGHT / 2) % HEIGHT * WIDTH + (x + WIDTH / 2) % WIDTH;
    //size_t new_index = xshift * WIDTH + yshift;

    // The additional pi is unclear to me, but necessary to match Python code.
    mask[index] = convert_uchar(
        fmod(
            (
                atan2(
                    field[new_index].y, 
                    field[new_index].x
                ) + corrections[index] + 2 * 3.14159265359f
            ) / (2 * 3.14159265359f) * 255, 
            256
        )
    );
};
