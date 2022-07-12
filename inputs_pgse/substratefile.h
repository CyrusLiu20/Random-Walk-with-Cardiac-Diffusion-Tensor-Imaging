#ifndef SUBSTRATEFILE_H
#define SUBSTRATEFILE_H

class substratefile{

public:

    substratefile(){
        voxel_range.resize(6);
        voxel_range << voxel_xmin, voxel_xmax, voxel_ymin, voxel_ymax, voxel_zmin, voxel_zmax; 
    };

    // Geometry 
    // Transform
    // y extent
    double yextent_min = -2000;
    double yextent_max = 5000;

    int deg_rot_per_L_in_y = 20; // degree of rotation per unit length in y
    double z_amplitude = 5; // amplitude in z for the sinusoidal transform
    double x_frequency = 4; // frquency in x for the sinusoidal transform
    bool shift_block = false;  

    std::string substrate_type = ""; // This will be inputted later when scanning myocytes

    // myocyte
    std::string filename = "geometry_1.mat";
    // size of the block
    double block_Lx = 495.3992; // x dimension
    double block_Ly = 392.3432; // y dimension
    double block_Lz = 126.5612; // z dimension


    // Domain
    // Imaging voxel
    Eigen::VectorXd voxel_range;
    double voxel_xmin = 0;
    double voxel_xmax = 494.899;
    double voxel_ymin = 0;
    double voxel_ymax = 391.843;
    double voxel_zmin = 58;
    double voxel_zmax = 68;


    // Transit model and permeability for the membranes
    std::string transit_model = "HybridModel"; // Type of model
    double permeability = 0.5;
    // Compartment specific diffusivity
    // Diffusitivity
    double D_ics = 4.5; // Intra-cellular diffusivity
    double D_ecs = 5.5; // Extra-cellular diffusivity


    //Dimension
    std::string dimension = "xy";
};

#endif