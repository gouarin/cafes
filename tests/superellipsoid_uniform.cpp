#include <particle/geometry/super_ellipsoid.hpp>
#include <particle/geometry/triangulation.hpp>
#include <iostream>
#include <chrono>
#include "hdf5.h"

int main()
{
  {
    cafes::geometry::super_ellipsoid<2> se{ {1.,1.}, {1.,1.}, 10. };

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    auto points = se.surface(.1);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


    cafes::geometry::super_ellipsoid<3> se3d{ {1.,1.,1.}, {1.,1.,1.}, 2., 2. };

    start = std::chrono::system_clock::now();
    auto points3d = se3d.surface({{.1, .1}});
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


    auto triangles = cafes::geometry::triangulation(points3d);
    std::cout << points.size() << " " << triangles.size() << "\n";

    hid_t file_id = H5Fcreate("test_mesh.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimsc[2] = {points.size(), 2};              // dataset dimensions
    hid_t dataspace_id= H5Screate_simple(2, dimsc, NULL);
    hid_t dataset_id = H5Dcreate(file_id,"points", 
                                 H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT,
                                 H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
             H5P_DEFAULT,points.data());
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    hsize_t dims[2] = {triangles.size(), 3};              // dataset dimensions
    dataspace_id= H5Screate_simple(2,dims,NULL);
    dataset_id = H5Dcreate(file_id,"triangles", 
                                 H5T_NATIVE_INT, dataspace_id,H5P_DEFAULT,
                                 H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
             H5P_DEFAULT,triangles.data());
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);

    }
    return 0;
}