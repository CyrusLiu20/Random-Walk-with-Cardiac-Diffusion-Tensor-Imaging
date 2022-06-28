#include <cstdio>
#include <chrono>
// #include <eigen3/Eigen/Dense>
#include "eigen3/Eigen/Dense"
#include <iostream>
#include <matio.h>

// Inputs
#include "inputs/montecarlofile.h"
#include "inputs/sequencefile.h"
#include "inputs/substratefile.h"

// Substrate
#include "substrate/read_myocytes.h"
#include "substrate/myocytes.h"
#include "substrate/substrate.h"

// Monte Carlo
#include "montecarlo/walkers.h"
#include "montecarlo/simulation.h"

// // MRI
#include "MRI/sequence.h"
#include "MRI/scansequence.h"

// Geometry
#include "geometry/polygon.h"
#include "geometry/boundingbox.h"

// CPP files VScode debugging use
#include "substrate/read_myocytes.cpp"
#include "montecarlo/walkers.cpp"
#include "montecarlo/simulation.cpp"
#include "MRI/scansequence.cpp"
#include "substrate/substrate.cpp"
#include "geometry/polygon.cpp"

int main(int argc, char *argv[]){

    // Monitor the computational time
    auto begin = std::chrono::steady_clock::now();

    // Configurating sequence
    ScanSequence scanner;
    sequencefile sequence_raw;
    sequence sequence;
    sequence = scanner.create(sequence_raw);



    // Configurating substrate
    substratefile substrate_raw;
    read_myocytes myocyte_scanner(substrate_raw.filename);
    myocytes myocytes;
    if (myocyte_scanner.scanned()){
        myocytes.load(myocyte_scanner);  
        substrate_raw.substrate_type = "block";
    }
    else{
        // To do: load myocytes from substratefile directly
    }
    substrate substrate(substrate_raw);
    substrate.myocytes_scan(myocytes);
    substrate.buildcache();
    std::vector<polygon> myos = substrate.myos;



    // Configurating monte carlo simulator
    montecarlofile montecarlo_raw;
    walkers particles(montecarlo_raw);
    Eigen::VectorXd buffer = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd seedbox;
    double simulation_time_total;

    int dim = substrate_raw.dimension.length();
    if (dim == 0){
        printf("run::sim::No dimension has been selected");
    }

    if (montecarlo_raw.seedbox_type == "text"){ // User used text input
        if (montecarlo_raw.seedbox_text.find("voxel") != std::string::npos) {
            seedbox = substrate_raw.voxel_range;
        }
        else if (montecarlo_raw.seedbox_text.find("origin")){
            seedbox = Eigen::VectorXd::Zero(6);
        }
        else{
            printf("run_sim::Unknown seedbox option, such input is not supported");
        }

        if (montecarlo_raw.seedbox_text.find("buffer")){
            simulation_time_total = sequence.dt.sum();
            for (int i = 0; i < dim; i++){
                buffer({i*2,i*2+1}) << -1, 1;
            }
            buffer = buffer * std::sqrt(6*substrate_raw.D_ecs*simulation_time_total);
            seedbox = seedbox + buffer;
        }
    }
    else{
        seedbox = montecarlo_raw.seedbox_vect;
    }

    simulation system(particles);
    Eigen::MatrixXd seedboxes = seedbox.transpose();
    system.seedParticlesInBox(seedboxes);
    Eigen::MatrixXd pos0 = system.get_position();

    std::cout << pos0 << std::endl;

    // Eigen::MatrixXd boundingbox(2,6); // Temporary bounding box
    // boundingbox << 1, 1.2, 3.2, 5.6, 4.3, 6.4,
    //                1, 1.2, 3.2, 5.6, 4.3, 6.4;
    // printf("Testing : seedInParticleBox testing results : %s\n", system.seedParticlesInBox(boundingbox) ? "True" : "false");

    auto finish = std::chrono::steady_clock::now();
	double run_time = std::chrono::duration_cast<std::chrono::milliseconds>(finish - begin).count()/(double)1000;///1000000.0;
    printf("Simulation time : %lf (seconds)\n", run_time);

    // Testing site

    return 0;
}

// Notes