    
#ifndef MRI_SCANS_TRILINEAR_H
#define MRI_SCANS_TRILINEAR_H

#include <algorithm>
#include "mriScansNearest.hpp"


// The density function, which takes p as input, where p 
// ranges from (0, 0, 0) 
// to (256, 256, 109), and then returns the density, using trilinear
// interpolation as an approach (for a continuous sampling).
real mriDensityTrilinear(
    const std::vector<std::vector<std::vector<uint16_t>>>& density,
    const Vector3& p
){
    int nx = density.size();
    int ny = density[0].size();
    int nz = density[0][0].size();

    // The corner indexes
    int x0 = std::clamp(static_cast<int>(p.x()), 0, nx - 1);
    int y0 = std::clamp(static_cast<int>(p.y()), 0, ny - 1);
    int z0 = std::clamp(static_cast<int>(p.z()), 0, nz - 1);

    // The other corner
    int x1 = std::min(x0 + 1, nx - 1);
    int y1 = std::min(y0 + 1, ny - 1);
    int z1 = std::min(z0 + 1, nz - 1);

    // Distances (the alpha value)
    real dx = p.x() - x0;
    real dy = p.y() - y0;
    real dz = p.z() - z0;

    // The densities
    real v000 = density[x0][y0][z0];
    real v100 = density[x1][y0][z0];
    real v010 = density[x0][y1][z0];
    real v110 = density[x1][y1][z0];
    real v001 = density[x0][y0][z1];
    real v101 = density[x1][y0][z1];
    real v011 = density[x0][y1][z1];
    real v111 = density[x1][y1][z1];

    // Interpolation in x
    real v00 = v000 * (1 - dx) + v100 * dx;
    real v10 = v010 * (1 - dx) + v110 * dx;
    real v01 = v001 * (1 - dx) + v101 * dx;
    real v11 = v011 * (1 - dx) + v111 * dx;

    // Interpolation in y
    real v0 = v00 * (1 - dy) + v10 * dy;
    real v1 = v01 * (1 - dy) + v11 * dy;

    // Interpolation in z
    real v = v0 * (1 - dz) + v1 * dz;
    
    return v;
}


// Tests marching cube using a the mri density
void mriScanTrilinearTest(){

    // These inputs are the image dimensions and number of images
    std::vector<std::vector<std::vector<uint16_t>>> density = 
        importData(256, 256, 109);

    // We have to use a lambda to capture the data
    auto f = [&density](const Vector3& p) -> real {
        return mriDensityTrilinear(density, p);
    };

    // The issue with nearest neighbor interpolation is that
    // we are limited by the image dimensions and number of slices;
    // using more cubes than that won't get us a smoother result.
    marchingCubes(f, Vector3(256.0, 256.0, 109.0), 1500, 200, 200, 
        100, false, "output/mri_trilinear_test.obj");
}

#endif
