#include "polygon.h"

polygon::polygon(Eigen::MatrixXd vertices_input, Eigen::MatrixXd faces_input){

    // Assigning vertices and faces
    Vertices = vertices_input;
    Vertices_t = vertices_input.transpose();
    Faces = faces_input;
    Faces_t = faces_input.transpose();

    // Number of faces and vertices
    n_vertices = vertices_input.rows();
    n_faces = faces_input.rows();

    polygon::volume_compute(vertices_input, faces_input);
    polygon::surface_compute(vertices_input, faces_input);


    // Computing the bounding box range
    bb_range.resize(6);
    minXYZ << vertices_input(Eigen::all, 0).minCoeff(), vertices_input(Eigen::all, 1).minCoeff(), vertices_input(Eigen::all, 2).minCoeff();
    maxXYZ << vertices_input(Eigen::all, 0).maxCoeff(), vertices_input(Eigen::all, 1).maxCoeff(), vertices_input(Eigen::all, 2).maxCoeff();
    bb_range << minXYZ(0), maxXYZ(0), minXYZ(1), maxXYZ(1), minXYZ(2), maxXYZ(2);

    boundingbox.initialize(bb_range);
}

// Substrate block bounding box
void polygon::initialize(std::string box_type, Eigen::MatrixXd bounding_box){
    Eigen::VectorXd rangespec(6);
    Eigen::MatrixXd vertices(8, 3);
    Eigen::MatrixXd faces(12,3);
    Eigen::Vector3d origin, range;
    if (box_type == "cuboid"){
        rangespec << bounding_box(0, 0), bounding_box(1, 0), bounding_box(0, 1), bounding_box(1, 1), bounding_box(0, 2), bounding_box(1, 2);
        faces << 1,3,4, 1,4,2, 5,6,8, 5,8,7, 2,4,8, 2,8,6,
                 1,5,7, 1,7,3, 1,2,6, 1,6,5, 3,7,8, 3,8,4;
        vertices << 0,0,0, 1,0,0, 0,1,0, 1,1,0, 0,0,1, 1,0,1, 0,1,1, 1,1,1;
        vertices = vertices.array() - 0.5;
        origin = (bounding_box(1, Eigen::all) + bounding_box(0, Eigen::all))*0.5;
        range = bounding_box(1, Eigen::all) - bounding_box(0, Eigen::all);

        vertices = vertices.array().rowwise()*range.transpose().array();
        vertices = vertices.rowwise() + origin.transpose();

        // Assigning vertices and faces
        Vertices = vertices;
        Vertices_t = vertices.transpose();
        Faces = faces;
        Faces_t = faces.transpose();

        // Number of faces and vertices
        n_vertices = Vertices.rows();
        n_faces = Faces.rows();

        // Compute mesh property
        polygon::volume_compute(vertices, faces);
        polygon::surface_compute(vertices, faces);

        // Computing the bounding box range
        bb_range.resize(6);
        minXYZ << vertices(Eigen::all, 0).minCoeff(), vertices(Eigen::all, 1).minCoeff(), vertices(Eigen::all, 2).minCoeff();
        maxXYZ << vertices(Eigen::all, 0).maxCoeff(), vertices(Eigen::all, 1).maxCoeff(), vertices(Eigen::all, 2).maxCoeff();
        bb_range << minXYZ(0), maxXYZ(0), minXYZ(1), maxXYZ(1), minXYZ(2), maxXYZ(2);

        boundingbox.initialize(bb_range);

    }
}

void polygon::volume_compute(Eigen::MatrixXd &vertices_input, Eigen::MatrixXd &faces_input){

    // Compute the volume and surface areas of the tetrahedron
    vertices_mean << vertices_input(Eigen::all, 0).mean(), vertices_input(Eigen::all, 1).mean(), vertices_input(Eigen::all, 2).mean(); 
    vertices_shifted = vertices_input.rowwise() - vertices_mean.transpose();
    // Since each vertex is associated with different face
    volume_tetra.resize(n_faces);
    Eigen::Matrix3d tetra; // To compute the volume of each tetrahedron
    for (int i = 0; i < n_faces; i++){
        tetra = vertices_shifted({faces_input(i, 0)-1, faces_input(i, 1)-1, faces_input(i, 2)-1}, Eigen::all);
        volume_tetra(i) = tetra.determinant()/6;
    }

    volume = volume > 0 ? volume : -volume;
}


// This is too slow!!!!!!
void polygon::surface_compute(Eigen::MatrixXd &vertices_input, Eigen::MatrixXd &faces_input){

    // Initialize area
    area_surface.resize(n_faces);
    Eigen::Matrix3d xyz, xy, yz, zx;
    Eigen::Vector3d ones = Eigen::VectorXd::Ones(3);
    Eigen::Vector3d x, y, z;

    xy = Eigen::Matrix3d::Ones(3,3);
    yz = Eigen::Matrix3d::Ones(3,3);
    zx = Eigen::Matrix3d::Ones(3,3);


    for (int i = 0; i < n_faces; i++){
        xyz = vertices_input({faces_input(i, 0)-1, faces_input(i, 1)-1, faces_input(i, 2)-1}, {0,1,2});
        x = xyz(Eigen::all, 0);
        y = xyz(Eigen::all, 1);
        z = xyz(Eigen::all, 2);  

        // x-y principal plane
        xy(Eigen::all,0) = x;
        xy(Eigen::all,1) = y;

        // y-z principal plane
        yz(Eigen::all,0) = y;
        yz(Eigen::all,1) = z;

        // z-x principal plane
        zx(Eigen::all,0) = z;
        zx(Eigen::all,1) = x;

        area_surface(i) = 0.5*std::sqrt(std::pow(xy.determinant(), 2) + std::pow(yz.determinant(), 2) + std::pow(zx.determinant(), 2));

    }
    surface_area = area_surface.sum();
}

Eigen::Vector3d polygon::get_minXYZ(){
    return minXYZ;
}