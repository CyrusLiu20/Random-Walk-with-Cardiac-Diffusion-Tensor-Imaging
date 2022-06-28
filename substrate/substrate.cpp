#include <algorithm>
#include "substrate.h"

substrate::substrate(substratefile substrate){

    // Voxel size
    y_extent << substrate.yextent_min, substrate.yextent_max;
    dxdydz << substrate.block_Lx, substrate.block_Ly, substrate.block_Lz;
    Ly = substrate.block_Ly;

    // y slicing
    Eigen::VectorXd low = substrate::iteration(0, -2000, 394);
    Eigen::VectorXd high = substrate::iteration(0, 5000, 394);
    y_minvals.resize(low.rows()+high.rows()-1);
    y_minvals << low(Eigen::seq(0,Eigen::last-1)), high;

    y_slice_minmax.resize(2, y_minvals.rows()-1);
    y_slice_minmax(0, Eigen::all) = y_minvals(Eigen::seq(0,Eigen::last-1));
    y_slice_minmax(1, Eigen::all) = y_minvals(Eigen::seq(1,Eigen::last));

    // Type of substrate partial or full
    substrate_type = substrate.substrate_type;

    // Initializing voxel
    voxel.initialize(substrate.voxel_range);
    transit_model = substrate.transit_model;
    kappa = substrate.permeability;
    D_i = substrate.D_ics;
    D_e = substrate.D_ecs;
    dim = substrate.dimension;

    // Copy transformation
    yextent_min = substrate.yextent_min;
    yextent_max = substrate.yextent_max;

    deg_rot_per_L_in_y = substrate.deg_rot_per_L_in_y; // degree of rotation per unit length in y
    z_amplitude = substrate.z_amplitude; // amplitude in z for the sinusoidal transform
    x_frequency = substrate.x_frequency; // frquency in x for the sinusoidal transform
    shift_block = substrate.shift_block;  

    if (substrate_type == "block"){
        boundary = "periodic";
        isIdentity = false;
    }

    // Creating substrate block bounding box

    // Parsing in Matlab (????????)
}

// DEBUG: try cross product? 
void substrate::myocytes_scan(myocytes myocytes_input){
    N_m = myocytes_input.Vertices.size();
    myos.resize(N_m);

    for (int i = 0; i < N_m; i++){
        polygon polygon_i(myocytes_input.Vertices[i], myocytes_input.Faces[i]);
        myos[i] = polygon_i;
        myos[i].bytes = sizeof(polygon_i);
    }

}

void substrate::buildcache(){
    
    // Concatonating the boundary range of each bounding box into a huge vector
    Eigen::VectorXd tmp_bbs;
    tmp_bbs.resize(6);
    myocytes_bbrange.resize(N_m*6);
    // myocytes_bbrange = Eigen::VectorXd::Zero(N_m*6);
    for (int i = 0; i < N_m; i++){
        tmp_bbs = myos[i].boundingbox.bb_range;
        myocytes_bbrange({i*6, i*6+1, i*6+2, i*6+3, i*6+4, i*6+5}) = tmp_bbs;
    }

    // Creating the substrate block bounding box
    Eigen::MatrixXd cuboid(2, 3);
    cuboid << 0, 0, 0, dxdydz(0), dxdydz(1), dxdydz(2);
    block_bb.initialize("cuboid", cuboid);
    block_bb.bytes = sizeof(block_bb);
}

// Stupid way to implement iteration
Eigen::VectorXd substrate::iteration(double a, double b, double step){    
    Eigen::VectorXd output;
    int size = std::floor(substrate::absolute((a-b)/step))+1;
    output.resize(size);

    if (a > b){
        for (int i = 0; i < size; i++){
            output(i) = a - step*i;
        }
        return output.reverse();
    }
    else{    
        for (int i = 0; i < size; i++){
            output(i) = a + step*i;
        } 
        return output;
    }
}

double substrate::absolute(double input){
    return input > 0 ? input : -input;
}

