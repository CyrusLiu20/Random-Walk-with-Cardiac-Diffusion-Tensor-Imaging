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

        // Saving all the doubles
        saving_file::save_double(matfp, post_processing.MD, "MD");
        saving_file::save_double(matfp, post_processing.Dx, "Dx");
        saving_file::save_double(matfp, post_processing.Dy, "Dy");
        saving_file::save_double(matfp, post_processing.Dz, "Dz");

        // Saving all the matrices
        saving_file::save_matrix(matfp, post_processing.tensor.rows(), post_processing.tensor.cols(), post_processing.tensor.data(), "tensor");
        saving_file::save_matrix(matfp, pos0.rows(), pos0.cols(), pos0.data(), "pos0");
        saving_file::save_matrix(matfp, post_processing.position_all.rows(), post_processing.position_all.cols(), post_processing.position_all.data(), "position");
        saving_file::save_matrix(matfp, post_processing.displacement_all.rows(), post_processing.displacement_all.cols(), post_processing.displacement_all.data(), "displacement");
        saving_file::save_matrix(matfp, post_processing.displacement_ECS.rows(), post_processing.displacement_ECS.cols(), post_processing.displacement_ECS.data(), "displacement_ECS");
        saving_file::save_matrix(matfp, post_processing.displacement_ICS.rows(), post_processing.displacement_ICS.cols(), post_processing.displacement_ICS.data(), "displacement_ICS");
        saving_file::save_matrix(matfp, post_processing.phase_all.rows(), post_processing.phase_all.cols(), post_processing.phase_all.data(), "phase");
        saving_file::save_matrix(matfp, post_processing.phase_ECS.rows(), post_processing.phase_ECS.cols(), post_processing.phase_ECS.data(), "phase_ECS");
        saving_file::save_matrix(matfp, post_processing.phase_ICS.rows(), post_processing.phase_ICS.cols(), post_processing.phase_ICS.data(), "phase_ICS");


        Mat_Close(matfp);
    }

private:

    void save_double(mat_t *matfp, double &input, const char* fieldname){
        double mydouble = input;
        size_t dim[2] = { 1, 1 };
        matvar_t *variable = Mat_VarCreate(fieldname, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim, &mydouble, 0);
        Mat_VarWrite(matfp, variable, MAT_COMPRESSION_NONE); //or MAT_COMPRESSION_ZLIB
        Mat_VarFree(variable);
    }


    void save_matrix(mat_t *matfp, unsigned int first, unsigned int second, double *input, const char* fieldname){
        // write matrix

        double *matrix = input;

        size_t dim2d[2] = {first, second};
        matvar_t *variable2d = Mat_VarCreate(fieldname, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d, matrix, 0); //rank 2
        Mat_VarWrite(matfp, variable2d, MAT_COMPRESSION_NONE);
        Mat_VarFree(variable2d);
    }



};

#endif