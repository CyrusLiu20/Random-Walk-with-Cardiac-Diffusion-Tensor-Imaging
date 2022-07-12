#ifndef SAVING_FILE_H
#define SAVING_FILE_H

#include <matio.h>
#include <eigen3/Eigen/Dense>
#include <vector>
#include "montecarlo/particle_state.h"
class saving_file{

public:

    saving_file(std::vector<particle_state> &particle_states, Eigen::MatrixXd &pos0, process_signal &post_processing){
        
        const char *filename = "output.mat";
        mat_t *matfp = NULL; //matfp contains pointer to MAT file or NULL on failure
        matfp = Mat_CreateVer(filename, NULL, MAT_FT_MAT5); //or MAT_FT_MAT4 / MAT_FT_MAT73
        //don't forget to close file with Mat_Close(matfp);

        // Saving all the doubles
        const char* name_MD = "MD";
        saving_file::save_double(matfp, post_processing.MD, name_MD);
        const char* name_Dx = "Dx";
        saving_file::save_double(matfp, post_processing.Dx, name_Dx);
        const char* name_Dy = "Dy";
        saving_file::save_double(matfp, post_processing.Dy, name_Dy);
        const char* name_Dz = "Dz";
        saving_file::save_double(matfp, post_processing.Dz, name_Dz);

        // Saving all the matrices
        unsigned int tensor_row = (unsigned int)post_processing.tensor.rows();
        unsigned int tensor_col = 3;
        double tensor[tensor_row][tensor_col] = { 0 };

        // fill 2d array
        for (int i = 0; i < tensor_row; i++){
            for (int j = 0; j < tensor_col; j++){
                tensor[i][j] = post_processing.tensor(i, j);
            }
        }

        // write displacement
        const char* fieldname2d_tensor = "tensor";
        size_t dim2d_tensor[2] = { tensor_col, tensor_row };
        matvar_t *variable2d_tensor = Mat_VarCreate(fieldname2d_tensor, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d_tensor, &tensor, 0); //rank 2
        Mat_VarWrite(matfp, variable2d_tensor, MAT_COMPRESSION_NONE);
        Mat_VarFree(variable2d_tensor);

        // const char* name_d = "displacement";
        // saving_file::save_matrix(matfp, post_processing.displacement, name_d);
        unsigned int displacement_row = (unsigned int)post_processing.displacement.rows();
        unsigned int displacement_col = 3;
        double displacement[displacement_row][displacement_col] = { 0 };

        // fill 2d array
        for (int i = 0; i < displacement_row; i++){
            for (int j = 0; j < displacement_col; j++){
                displacement[i][j] = post_processing.displacement(i, j);
            }
        }

        // write displacement
        const char* fieldname2d_displacement = "displacement";
        size_t dim2d_displacement[2] = { displacement_col, displacement_row };
        matvar_t *variable2d_displacement = Mat_VarCreate(fieldname2d_displacement, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d_displacement, &displacement, 0); //rank 2
        Mat_VarWrite(matfp, variable2d_displacement, MAT_COMPRESSION_NONE);
        Mat_VarFree(variable2d_displacement);

        unsigned int displacement_ECS_row = (unsigned int)post_processing.displacement_ECS.rows();
        unsigned int displacement_ECS_col = 3;
        double displacement_ECS[displacement_ECS_row][displacement_ECS_col] = { 0 };

        // fill 2d array
        for (int i = 0; i < displacement_ECS_row; i++){
            for (int j = 0; j < displacement_ECS_col; j++){
                displacement_ECS[i][j] = post_processing.displacement_ECS(i, j);
            }
        }

        // write position
        const char* fieldname2d_displacement_ECS = "displacement_ECS";
        size_t dim2d_displacement_ECS[2] = { displacement_ECS_col, displacement_ECS_row };
        matvar_t *variable2d_displacement_ECS = Mat_VarCreate(fieldname2d_displacement_ECS, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d_displacement_ECS, &displacement_ECS, 0); //rank 2
        Mat_VarWrite(matfp, variable2d_displacement_ECS, MAT_COMPRESSION_NONE);
        Mat_VarFree(variable2d_displacement_ECS);

        unsigned int displacement_ICS_row = (unsigned int)post_processing.displacement_ICS.rows();
        unsigned int displacement_ICS_col = 3;
        double displacement_ICS[displacement_ICS_row][displacement_ICS_col] = { 0 };

        // fill 2d array
        for (int i = 0; i < displacement_ICS_row; i++){
            for (int j = 0; j < displacement_ICS_col; j++){
                displacement_ICS[i][j] = post_processing.displacement_ICS(i, j);
            }
        }

        // write position
        const char* fieldname2d_displacement_ICS = "displacement_ICS";
        size_t dim2d_displacement_ICS[2] = { displacement_ICS_col, displacement_ICS_row };
        matvar_t *variable2d_displacement_ICS = Mat_VarCreate(fieldname2d_displacement_ICS, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d_displacement_ICS, &displacement_ICS, 0); //rank 2
        Mat_VarWrite(matfp, variable2d_displacement_ICS, MAT_COMPRESSION_NONE);
        Mat_VarFree(variable2d_displacement_ICS);

        // const char* name_dECS = "displacement_ECS";
        // saving_file::save_matrix(matfp, post_processing.displacement_ECS, name_dECS);

        // const char* name_dICS = "displacement_ICS";
        // saving_file::save_matrix(matfp, post_processing.displacement_ICS, name_dICS);

        // const char* name_pECS = "phase_ECS";
        // saving_file::save_matrix(matfp, post_processing.phase_ECS, name_pECS);

        // const char* name_pICS = "phase_ECS";
        // saving_file::save_matrix(matfp, post_processing.phase_ICS, name_pICS);
        

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

        unsigned int pos0_row = (unsigned int)pos0.rows();
        unsigned int pos0_col = 3;
        double initial_position[pos0_row][pos0_col] = { 0 };

        // fill 2d array
        for (int i = 0; i < pos0_row; i++){
            for (int j = 0; j < pos0_col; j++){
                initial_position[i][j] = pos0(i, j);
            }
        }

        // write position
        const char* fieldname2d_pos0 = "pos0";
        size_t dim2d_pos0[2] = { pos0_col, pos0_row };
        matvar_t *variable2d_pos0 = Mat_VarCreate(fieldname2d_pos0, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d_pos0, &initial_position, 0); //rank 2
        Mat_VarWrite(matfp, variable2d_pos0, MAT_COMPRESSION_NONE);
        Mat_VarFree(variable2d_pos0);

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

private:

    void save_double(mat_t *matfp, double input, const char* fieldname){
        double mydouble = input;
        size_t dim[2] = { 1, 1 };
        matvar_t *variable = Mat_VarCreate(fieldname, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim, &mydouble, 0);
        Mat_VarWrite(matfp, variable, MAT_COMPRESSION_NONE); //or MAT_COMPRESSION_ZLIB
        Mat_VarFree(variable);
    }

    void save_matrix(mat_t *matfp, Eigen::MatrixXd input, const char* fieldname){
        // write matrix

        unsigned int first = (unsigned int)input.rows();
        unsigned int second = (unsigned int)input.cols();
        double matrix[first][second] = { 0 };

        // fill 2d array
        for (int i = 0; i < first; i++){
            for (int j = 0; j < second; j++){
                matrix[i][j] = input(i, j);
            }
        }

        size_t dim2d[2] = {first, second};
        matvar_t *variable2d = Mat_VarCreate(fieldname, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d, &input, 0); //rank 2
        Mat_VarWrite(matfp, variable2d, MAT_COMPRESSION_NONE);
        Mat_VarFree(variable2d);
    }

    void save_matrix(mat_t *matfp, Eigen::Matrix3d input, const char* fieldname){
        // write matrix

        unsigned int first = (unsigned int)input.rows();
        unsigned int second = (unsigned int)input.cols();
        double matrix[first][second] = { 0 };

        // fill 2d array
        for (int i = 0; i < first; i++){
            for (int j = 0; j < second; j++){
                matrix[i][j] = input(i, j);
            }
        }

        size_t dim2d[2] = {first, second};
        matvar_t *variable2d = Mat_VarCreate(fieldname, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d, &input, 0); //rank 2
        Mat_VarWrite(matfp, variable2d, MAT_COMPRESSION_NONE);
        Mat_VarFree(variable2d);
    }



};

#endif