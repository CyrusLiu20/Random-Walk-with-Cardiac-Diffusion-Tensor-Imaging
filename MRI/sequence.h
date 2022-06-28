#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <eigen3/Eigen/Dense>

class sequence{

public:

    // Construtors
    sequence() = default;
	sequence(const sequence &pSrc) = default;	// copying class
	sequence(sequence &&pSrc) = default;	// moving class
	~sequence() = default;

    // Operators overloading
    sequence &operator =(const sequence &pSrc){
        // Transfer of resources
        if(this != &pSrc){
            N = pSrc.N;
            bvalue = pSrc.bvalue;
            dt = pSrc.dt;
            gG = pSrc.gG;
        }
        return *this;
    };	// Transfer of resourses using equal sign

    // Number of dt
    int N;
    // What is this?
    double bvalue;
    // Time step
    Eigen::VectorXd dt;
    // Something gradient
    Eigen::MatrixXd gG;

};

#endif