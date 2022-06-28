#include <cstdio>
#include <iostream>
#include <stdlib.h>

// header file
#include "simulation.h"

simulation::simulation(walkers obj_input){
    obj = obj_input;
}

bool simulation::seedParticlesInBox(Eigen::MatrixXd boundingboxes_input, int particlesPerBox_input = -1){
    bool flag = true;


    // Initialize Mersenne Twister pseudo-random number generator
    // std::mt19937 gen(rd());
    std::mt19937 gen(obj.get_rng_seed());

    // Obtain the dimension
    cols = (int)boundingboxes_input.cols();
    // Obtain the number of bounding boxes
    N_b = (int)boundingboxes_input.rows();

    printf("Number of rows : %d ; Number of columns : %d\n", N_b, cols);

    // initialization of matrices and vectors
    if (not(simulation::initialize(N_b))){
        flag = false;
    }

    // Check the dimensions of the bounding boxes
    switch (cols){
        case 2:
            boundingboxes(Eigen::all, Eigen::seqN(0,2)) = boundingboxes_input;
            dimension = 1;
            break;
        case 4:
            boundingboxes(Eigen::all, Eigen::seqN(0,4)) = boundingboxes_input;
            dimension = 2;
            break;
        case 6:
            boundingboxes = boundingboxes_input;
            dimension = 3;
            break;
        default:
            printf("simulation::seedParticlesInBox::inconsistent, BBs must be [M,N] with M=2,4,6 and N=N_p\n");
            flag = false;               
    }

    std::cout << "The updated bounding boxes : " << std::endl << boundingboxes << std::endl;

    // Check bounding box inconsistentcy
    for (int i = 0; i < N_b; i ++){
        bool inconsistent = (boundingboxes(i,0) > boundingboxes(i,1) or boundingboxes(i,2) > boundingboxes(i,3) or boundingboxes(i,4) > boundingboxes(i,5));
        printf("Testing : inconsistency : %s\n", inconsistent ? "true" : "false");
        if (inconsistent){
            printf("ParticleWalker:seedParticlesInBox:inconsistent, BBs must be [[xmin;xmax;ymin;ymax;zmin;zmax]]\n");
            flag = false;
        }
    }

    // Check overlap
    // To do

    if (particlesPerBox_input == -1){
        // Calculate the side lengths of bounding boxes
        Eigen::MatrixXd max, min;
        max = boundingboxes(Eigen::all, {1,3,5});
        min = boundingboxes(Eigen::all, {0,2,4});
        sidelengths_raw = max - min;

        // Check if and dimension has zero length in specified dimensions (e.g. only x and y)
        hasZeroDim = (sidelengths_raw(Eigen::all, Eigen::seqN(0,dimension)).array() == 0.0).any();
        if (hasZeroDim){
            printf("simulation::seedParticlesInBox::inconsistent, each dimension must be either zero or non-zero across all boxes");
            flag = false;
        }
        printf("Testing : Zero dimension : %s\n", hasZeroDim ? "true" : "false");

        // Removing zero dimensions (?????)
        sidelengths = sidelengths_raw(Eigen::all, Eigen::seqN(0, dimension));

        // Calculating the box volume of all the boxes
        boxVolumes = sidelengths.rowwise().prod();
        total_boxVolume = boxVolumes.sum();

        // Calculating the particles per box
        particlesPerBox_theo = boxVolumes/total_boxVolume*obj.get_N_p();
        particlesPerBox_prel = Eigen::floor(particlesPerBox_theo.array());
        int missing_particles = obj.get_N_p() - particlesPerBox_prel.sum();

        // Refilling missing particles
        particlesPerBox = simulation::refill(particlesPerBox_prel, missing_particles);
    }
    else if(particlesPerBox_input > 0){
        particlesPerBox.resize(N_b);
        particlesPerBox = Eigen::VectorXd::Ones(N_b, 1) * particlesPerBox;
        // Validation (to do) don't know how
    }

    // Check if there are still any missing particles
    printf("Testing : difference in particle number : %d\n", particlesPerBox.sum() - obj.get_N_p());
    if (particlesPerBox.sum() != obj.get_N_p()){
        printf("simulation::seedParticlesInBox::missing_particles, please check your preliminary particles number");
        flag = false;
    }
    else{
        printf("Testing : Missing particles : Passed\n");
    }

    std::cout << "The updated sidelengths : " << std::endl << sidelengths << std::endl;
    std::cout << "The box volumes : " << boxVolumes.transpose() << std::endl;
    std::cout << "The total box volume : " << total_boxVolume << std::endl;
    std::cout << "The theoretical particles per box : " << particlesPerBox_theo.transpose() << std::endl;
    std::cout << "The particles per box : " << particlesPerBox.transpose() << std::endl;

    // Seeding particles in bounding boxes
    if (not(simulation::seeding())){
        printf("simulation::seedParticlesInBox::seeding");
    }

    return flag;
}

bool simulation::initialize(int N_b){
    bool flag = true;
    // Initializing the bounding boxes, side lengths of bounding boxes, and box volumes
    boundingboxes = Eigen::MatrixXd::Zero(N_b, 6);
    sidelengths = Eigen::MatrixXd::Zero(N_b, 3);
    boxVolumes = Eigen::VectorXd::Zero(N_b, 1);
    // Initializing the preliminary and theoretical particles per box
    particlesPerBox_theo = Eigen::VectorXd::Zero(N_b,1);
    particlesPerBox_prel = Eigen::VectorXd::Zero(N_b,1);

    return flag;
}

// A custom function for refilling missing particles
Eigen::VectorXd simulation::refill(Eigen::VectorXd particlesPerBox, int missing_particles){

    // Initialize Mersenne Twister pseudo-random number generator
    // std::mt19937 gen(rd());
    std::mt19937 gen(obj.get_rng_seed());


    int rand_index;
    int number = particlesPerBox.rows();

    std::uniform_int_distribution<> rand_index_gen(0, number-1);
    for (int i = 0; i < missing_particles; i++){
        rand_index = rand_index_gen(gen);
        particlesPerBox(rand_index) = particlesPerBox(rand_index) + 1;
    }

    return particlesPerBox;
}

bool simulation::seeding(){

    // Flag for seeding
    bool flag = true;

    // Initialize Mersenne Twister pseudo-random number generator
    std::mt19937 gen(obj.get_rng_seed());
    // std::mt19937 gen(rd());

    // Initialization
    Eigen::VectorXd boundingBox = Eigen::VectorXd::Zero(N_b, 1).transpose();
    int ip_SetLast = 0; // Index of the particles
    int index_first, index_last, nP_iB;
    
    for (int iB = 0; iB < N_b; iB++){
        boundingBox = boundingboxes(iB, Eigen::all);
        nP_iB = particlesPerBox(iB);
        index_first = ip_SetLast;
        index_last = ip_SetLast + nP_iB;
        ip_SetLast += nP_iB;

        std::uniform_real_distribution<> x_interval(boundingBox(0), boundingBox(1));
        std::uniform_real_distribution<> y_interval(boundingBox(2), boundingBox(3));
        std::uniform_real_distribution<> z_interval(boundingBox(4), boundingBox(5));


        for (int j = index_first; j < index_last; j++){

            obj.position(j, 0) = x_interval(gen);
            obj.position(j, 1) = y_interval(gen);
            obj.position(j, 2) = z_interval(gen);
        }
    }

    return flag;
}

Eigen::MatrixXd simulation::get_position(){
    return obj.position;
}

// bool simulation::performScan()
// Notes
// dimension must be sequential