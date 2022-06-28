#ifndef POLYGON_H
#define POLYGON_H

#include <eigen3/Eigen/Dense>
#include "boundingbox.h"

class polygon{

public:

    polygon() = default;
    polygon(Eigen::MatrixXd vertices_input, Eigen::MatrixXd faces_input);
    polygon(std::string box_type, Eigen::MatrixXd bounding_box);

    // Create substrate bounding box
    void initialize(std::string box_type, Eigen::MatrixXd bounding_box);

    // Inputs the particle position
    bool containspoint(Eigen::Vector3d point);
    // Whether particle intersects boundary
    // bool intersect();

    // Vertices of one myocyte
    Eigen::MatrixXd Vertices;
    // Faces of one myocyte
    Eigen::MatrixXd Faces;
    // Transposed vertices of one myocyte
    Eigen::MatrixXd Vertices_t;
    // Transposed faces of one myocyte
    Eigen::MatrixXd Faces_t;

    boundingbox boundingbox;
    
    
    // Number of faces
    int n_faces;
    // Number of vertices
    int n_vertices;
    // Myocyte volume
    double volume;
    // Myocyte surface area
    double surface_area;
    // Number of bytes
    int bytes;

    // Data retreival
    Eigen::Vector3d get_minXYZ();

private:

    // Compute the mesh property
    void volume_compute(Eigen::MatrixXd &vertices_input, Eigen::MatrixXd &faces_input);
    void surface_compute(Eigen::MatrixXd &vertices_input, Eigen::MatrixXd &faces_input);

    // The minimum values of the vertices
    Eigen::Vector3d minXYZ;
    // The maximum values of the vertices
    Eigen::Vector3d maxXYZ;

    // Bounding box range
    Eigen::VectorXd bb_range;

    // Mean of the vertices
    Eigen::Vector3d vertices_mean;
    // A shifted vertices
    Eigen::MatrixXd vertices_shifted;
    // Volumes of each tetrahedron
    Eigen::VectorXd volume_tetra; 
    // Areas of each tetrahedron
    Eigen::VectorXd area_surface; 

};

#endif