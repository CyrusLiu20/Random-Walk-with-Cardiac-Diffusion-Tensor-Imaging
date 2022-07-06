#ifndef SAVING_FILE_H
#define SAVING_FILE_H

#include <matio.h>
#include <eigen3/Eigen/Dense>
#include <vector>
#include "montecarlo/particle_state.h"
class saving_file{

public:

    saving_file(std::vector<particle_state> particle_states){
        
        const char *filename = "output.mat";
        mat_t *matfp = NULL; //matfp contains pointer to MAT file or NULL on failure
        matfp = Mat_CreateVer(filename, NULL, MAT_FT_MAT5); //or MAT_FT_MAT4 / MAT_FT_MAT73
        //don't forget to close file with Mat_Close(matfp);

        unsigned int first = (unsigned int)particle_states.size();
        unsigned int second = 5;
        double position[first][second] = { 0 };

        // fill 2d array
        for (int i = 0; i < first; i++){
            for (int j = 0; j < second-2; j++){
                position[i][j] = particle_states[i].position(j);
            }
            position[i][3] = particle_states[i].myoindex;
            position[i][4] = particle_states[i].flag;
        }

        unsigned int phase_row = (unsigned int)particle_states.size();
        unsigned int phase_col = 5;
        double phase[phase_row][phase_col] = { 0 };

        // fill 2d array
        for (int i = 0; i < phase_row; i++){
            for (int j = 0; j < phase_col-2; j++){
                phase[i][j] = particle_states[i].phase(j);
            }
            phase[i][3] = particle_states[i].myoindex;
            phase[i][4] = particle_states[i].flag;
        }

        // write position
        const char* fieldname2d_position = "position";
        size_t dim2d_position[2] = { second, first };
        matvar_t *variable2d_position = Mat_VarCreate(fieldname2d_position, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d_position, &position, 0); //rank 2
        Mat_VarWrite(matfp, variable2d_position, MAT_COMPRESSION_NONE);
        Mat_VarFree(variable2d_position);

        // write phase
        const char* fieldname2d_phase = "phase";
        size_t dim2d_phase[2] = { phase_col, phase_row};
        matvar_t *variable2d_phase = Mat_VarCreate(fieldname2d_phase, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d_phase, &phase, 0); //rank 2
        Mat_VarWrite(matfp, variable2d_phase, MAT_COMPRESSION_NONE);
        Mat_VarFree(variable2d_phase);
        
        Mat_Close(matfp);
    }

};

#endif