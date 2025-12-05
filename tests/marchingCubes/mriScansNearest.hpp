    
#ifndef MRI_SCANS_NEAREST_H
#define MRI_SCANS_NEAREST_H

#include <algorithm>
#include "../../marchingCubes.hpp"
   

// This is a function that returns the pixel values of MRI scans
// as a 3D array. The xy axes are the image axes, and the z axis
// varies along the different slices (images).
// We also send in the dimensions and number of images.
std::vector<std::vector<std::vector<uint16_t>>> importData(
    int nx, int ny, int nz){

    std::vector<std::vector<std::vector<uint16_t>>> density(nx, 
    std::vector<std::vector<uint16_t>>(ny, 
        std::vector<uint16_t>(nz, 0)));

    // The slices
    for (int z = 0; z < nz; z++) {

        // This is just the slice file name (MRbrain.1, MRbrain.2, etc ...)
        std::string filename = "input/mri_data/MRbrain." 
            + std::to_string(z + 1);

        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            std::cerr << "Error opening the file:" << filename << "\n";
            exit(1);
        }

        // With the file open, we can start reading pixel data.
        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {

                // We can then read 2 bytes (16 bits) at a time
                // and populate the density function.
                char c[2];
                file.read(c, 2);
                if (!file) {
                    std::cerr << "Error reading pixel at:" 
                        << x << ", " << y << ", " << z << "\n";
                    exit(1);
                }

                // Converting the bit to a uint16 
                // All this is doing is shifting the first 8 bits (most
                // significant) to the right, and then following them by
                // the next (least significant) bits from the second byte.
                density[x][y][z] = (
                    static_cast<unsigned char>(c[0]) << 8) |
                    static_cast<unsigned char>(c[1]
                );
            }
        }
        file.close();
    }

    return density;
}


// The density function, which takes p as input, where p 
// ranges from (0, 0, 0) 
// to (256, 256, 109), and then returns the density, using the nearest
// neighbor approach (for a continuous sampling).
real mriDensityNearest(
    const std::vector<std::vector<std::vector<uint16_t>>>& density,
    const Vector3& p
){
    int nx = density.size();
    int ny = density[0].size();
    int nz = density[0][0].size();

    int x = std::clamp(static_cast<int>(p.x()), 0, nx - 1);
    int y = std::clamp(static_cast<int>(p.y()), 0, ny - 1);
    int z = std::clamp(static_cast<int>(p.z()), 0, nz - 1);

    return static_cast<real>(density[x][y][z]);
}


// Tests marching cube using a the mri density
void mriScanNearestTest(){

    // These inputs are the image dimensions and number of images
    std::vector<std::vector<std::vector<uint16_t>>> density = 
        importData(256, 256, 109);

    // We have to use a lambda to capture the data
    auto f = [&density](const Vector3& p) -> real {
        return mriDensityNearest(density, p);
    };

    // The issue with nearest neighbor interpolation is that
    // we are limited by the image dimensions and number of slices;
    // using more cubes than that won't get us a smoother result.
    marchingCubes(f, Vector3(256.0, 256.0, 109.0), 1500, 256, 256, 
        109, false, "output/mri_nearest_test.obj");
}

#endif
