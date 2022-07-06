#include "transform.h"
#include <math.h>

// Eigen::Vector3d transform::global2local(Eigen::Vector3d position){
transform_parameter transform::global2local(Eigen::Vector3d &position){

    transform_parameter output;
    // Eigen::Vector3d position_local = Eigen::Vector3d::Zero(3);
    double dx = dxdydz_bb(0);
    double dy = dxdydz_bb(1);
    double dz = dxdydz_bb(2);


    // Handle identity case first
    // double angle, angle_reverse; // The rotation angle in rotate_y function
    if (isIdentity){
        output.angle = 0;
        output.angle_reverse = 0; // They're the same since no rotation has been implemented
        
        // To do : transform inverse
        output.position_local = position; // No transformation (position_local = position)
        return output;
    }

    // std::cout << "Find y slice position input : " << position(1) << std::endl;
    double y_slice = transform::find_yslice(position(1));

    if (transform::absolute(y_slice) < dy/2){ // i.e. slice [0, dy], NOT [-dy, 0] (because we check the first coordinate, i.e. 0)
        output.angle = 0;
        output.angle_reverse = 0;
    }
    else{
        output.angle = (deg_rot_per_L_in_y*M_PI/180)*y_slice;
        output.angle_reverse = -output.angle;
    }

    // Rotate position
    Eigen::Vector3d position_rotated, position_slice; 
    
    // See if it's more efficient without lamda function
    // auto fn_Rot = [](Eigen::Vector3d &pos, double &theta){Eigen::Matrix3d rotation;
    //                                                     rotation << std::cos(theta), 0, std::sin(theta),
    //                                                                 0, 1, 0,
    //                                                                 -std::sin(theta), 0, std::cos(theta);
    //                                                     Eigen::Vector3d output = rotation*pos;  
    //                                                     return output;};
    // position_rotated = fn_Rot(position, angle_reverse);

    position_rotated = transform::rotate_y(position, output.angle_reverse);
    position_slice << position_rotated(0), std::fmod(position_rotated(1), dy), position_rotated(2);

    // Apply sin(x) displacement in z'
    double xCoord = position_slice(0);
    double yCoord = position_slice(1);
    double zCoord = position_slice(2) - transform::sin(dx, dz, xCoord);

    // Compute iY and iZ
    output.iY = 1+std::floor(position_rotated(1)/dy);
    output.iZ = 1+std::floor(zCoord/dz);

    // Shift axis half block to the right if iZ is odd and compute iX
    if (shift_block){
        xCoord = xCoord - std::fmod(output.iZ, 2)*dx/2;
    }
    output.iX = 1 + std::floor(xCoord/dx); // Compute iX like before with iY & iZ

    // Compute offset from rotated frame
    double ddxx = (output.iX-1)*dx;
    double ddzz = (output.iZ-1)*dz;

    output.position_local << xCoord - ddxx, yCoord, zCoord - ddzz;
    // std::cout << output.position_local << std::endl;

    return output;
}
// Notes 
// Did not use lambda functions
// some reason dz in Matlab does not change at all

Eigen::Vector3d transform::local2global(Eigen::Vector3d &pos_local, double &iX, double &iY, double &iZ){
    double iX_new = iX - 1;
    double iY_new = iY - 1;
    double iZ_new = iZ - 1;

    // Translation offset
    int shift;
    if (shift_block){
        shift = std::fmod(iZ_new+1, 2)/2;
    }
    else{
        shift = 0;
    }
    Eigen::Vector3d shifting_factor;
    shifting_factor << iX_new + shift, iY_new, iZ_new;
    Eigen::Vector3d offset = shifting_factor.array()*dxdydz_bb.array();

    // Compute position in the rotated plane
    Eigen::Vector3d position_rotated = pos_local + offset;
    position_rotated(2) = position_rotated(2) + transform::sin(dxdydz_bb(0), dxdydz_bb(2), position_rotated(0));

    // Invert rotation
    double y = dxdydz_bb(1)*iY_new;
    double angle = deg_rot_per_L_in_y*M_PI/180;
    double theta = angle*y;

    Eigen::Vector3d pos_global = transform::rotate_y(position_rotated, theta);

    return pos_global;
}

Eigen::Vector3d transform::rotate_y(Eigen::Vector3d &position, double &theta){
    Eigen::Matrix3d rotation;
    rotation << std::cos(theta), 0, std::sin(theta),
                 0, 1, 0,
                -std::sin(theta), 0, std::cos(theta);
    
    return rotation*position;
}

// Apply sin(x) displacement in z
double transform::sin(double &dx, double &dz, double &x){
    double amplitude = z_amplitude*dz;
    double val = amplitude*std::sin(x_frequency * 2*x*M_PI/dx);
    return val;
}

double transform::find_yslice(double &y_position){
    // To do : Use a binary search (now the worst possible method)
    bool found = false;
    for (int i = 0; i < y_slice_minmax.cols(); i++){
        if (y_position > y_slice_minmax(0, i) and y_position < y_slice_minmax(1, i)){
            found =  true;
            return y_slice_minmax(0, i);
        }
    }

    if (not(found)){
        std::cout << "Transform::find_yslice::where', 'Corresponding slice not found'" << std::endl;
        return -1;
    }

    return -1; // just in case
}

// custom stupid absolute function
double transform::absolute(double &input){
    return input > 0 ? input : -input;
}

// // custom stupid vector mod function
// Eigen::Vector3d transform::vect_mod(Eigen::Vector3d &original, Eigen::Vector3d &vector_mod){
//     Eigen::Vector3d output;
//     output << original(0) % vector_mod(0), original(1) % vector_mod(1), original(2) % vector_mod(2);
//     return output;
// }
