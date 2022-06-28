#ifndef SUBSTRATE_H
#define SUBSTRATE_H

#include "myocytes.h"
#include "../geometry/polygon.h"

class substrate{

public:

    substrate(substratefile substratefile);
    void myocytes_scan(myocytes myocytes_input);

    // Build the bounding box for the block
    void buildcache();

    // Detailed version of myocytes with bounding box, volume etc
    std::vector<polygon> myos;

    // load other geometry parameters
    Eigen::Vector3d dxdydz;

    // Voxel range, volume and surface area
    boundingbox voxel;
    // substrate block range, volume and surface area
    polygon block_bb;

    Eigen::VectorXd myocytes_bbrange;

    // boundary type
    std::string boundary;

    // Diffusion properties
    std::string transit_model;
    double kappa;
    double D_e;
    double D_i;
    std::string dim;

private:

    // Stupid iteration method
    Eigen::VectorXd iteration(double a, double b, double step);
    // Custom absolute function
    double absolute(double input);

    // myocyte geometry
    myocytes myocytes_process;

    // Maximum and minimum of voxel length in y dimension
    Eigen::Vector2d y_extent;
    // Realistic y slice
    Eigen::VectorXd y_minvals;
    Eigen::MatrixXd y_slice_minmax;

    // length in y dimension
    double Ly;
    // Number of myocytes
    int N_m;

    // substrate type
    std::string substrate_type;

    // Transform parameters
    double yextent_min;
    double yextent_max;

    int deg_rot_per_L_in_y; // degree of rotation per unit length in y
    double z_amplitude; // amplitude in z for the sinusoidal transform
    double x_frequency; // frquency in x for the sinusoidal transform
    bool shift_block;  
    
    bool isIdentity;

    Eigen::Vector3d dxdydz_bb;
};

#endif