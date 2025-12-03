
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdint>
#include <cmath>

#include "vector3.hpp"
#include "lookupTables.hpp"


int main() {


    //================================ DATA =================================//


    // Our 3D data is described by 109 slices (2D images).
    // Each slice is 256 by 256 pixels, so the domain of x, y, and z
    // ranges from (0, 0, 0) to (255, 255, 108) (here the multiple
    // image slices give us z, the third dimension).
    const int NX = 256;
    const int NY = 256;
    const int NZ = 109;

    // This boolean indicates whether or not triangle vertices are
    // calculated by linearly interpolating the surface and intersecting
    // the cube edges, or just taking the edge middle point. 
    const bool USE_MIDDLE_POINT = true; 

    // This gives us a f(x, y, z), which is a scalar field representing
    // the measured density at a point (x, y, z). This will depend on
    // the imaging technique used (MRI, CT, XRAY...).
    // The output is a 16-bit integer. 
    std::vector<std::vector<std::vector<uint16_t>>> density(NX, 
        std::vector<std::vector<uint16_t>>(NY, 
            std::vector<uint16_t>(NZ, 0)));

    // We will also calculate the gradient of the above scalar field
    // at each of the voxels. We can only approximate it of course, as the
    // data is discrete.
    std::vector<std::vector<std::vector<Vector3>>> gradient(NX, 
        std::vector<std::vector<Vector3>>(NY, 
        std::vector<Vector3>(NZ)));

    // The surface defined by f(x, y, z) = suruface value is essentially 
    // just the surface that seperates the outside of the captured object 
    // from the inside.
    // To find a suitable surface value, we can visualize some slices while
    // masking pixels below different thresholds. In practice, this
    // value here seemed to work well.
    const uint16_t SURFACE_VALUE = 1500;


    // Now to populate our function, we start by picking a slice,
    // which is to say, loop over all z-values.
    for (int z = 0; z < NZ; z++) {
        // This is just the slice file name (MRbrain.1, MRbrain.2, etc ...)
        std::string filename = "data/MRbrain." + std::to_string(z + 1);
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            std::cerr << "Error opening the file:" << filename << "\n";
            exit(1);
        }

        // With the file open, we can start reading pixel data.
        for (int y = 0; y < NY; y++) {
            for (int x = 0; x < NX; x++) {

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


    // Now we approximate the gradient
    for (int x = 0; x < NX; x++) {
        for (int y = 0; y < NY; y++) {
            for(int z = 0; z < NZ; z++){

                int x_p1 = std::min(x+1, NX-1);  
                int x_m1 = std::max(x-1, 0);  
                int y_p1 = std::min(y+1, NY-1);  
                int y_m1 = std::max(y-1, 0);  
                int z_p1 = std::min(z+1, NZ-1);  
                int z_m1 = std::max(z-1, 0); 

                gradient[x][y][z] = Vector3(
                    (density[x_p1][y][z] - density[x_m1][y][z]) / (x_p1 - x_m1),
                    (density[x][y_p1][z] - density[x][y_m1][z]) / (y_p1 - y_m1),
                    (density[x][y][z_p1] - density[x][y][z_m1]) / (z_p1 - z_m1)
                );

                gradient[x][y][z].normalize();
            }
        }
    }


    //=========================== MARCHING CUBES ============================//


    // The marching cubes algorithm now requires us to divide the 3D
    // space into a grid of 255 * 255 * 108 cubes.
    // Based on whether their vertices are inside or outside
    // the surface (intensity at vertex below or above the surface value),
    // we can tell how (if at all) the surface intersets the cube,
    // which in turn tells how to generate a triangle mesh.

    // For a less accurate mesh, we can use less cubes.

    // Pixel offset for each cube corner
    int offset[8][3]{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

    
    struct Cube {
        // Corners are true if inside the surface, false otherwise
        bool corners[8];

        // The edges are true if the surface passes through them,
        // which we can deduce by cross referencing the corners with
        // the reference table.
        bool edges[12];

        // Each triangle's vertices are stored as the index of the
        // edge they belong to. 
        std::vector<std::vector<int>> triangles;

        // Each edge's intersection with the surface in world coordinates.
        float vertex_coordinate[12][3];
        
        // Each intersection point's normal (assumed to be normalized)
        float vertex_normal[12][3];
    };


    // Maps each edge index to the two vertices it connects
    // (as shifts).
    const int edge_shifts[12][6] = {
        {0,0,0, 1,0,0},
        {1,0,0, 1,1,0},
        {1,1,0, 0,1,0},
        {0,1,0, 0,0,0},
        {0,0,1, 1,0,1},
        {1,0,1, 1,1,1},
        {1,1,1, 0,1,1},
        {0,1,1, 0,0,1},
        {0,0,0, 0,0,1},
        {1,0,0, 1,0,1},
        {1,1,0, 1,1,1},
        {0,1,0, 0,1,1}
    };

    // Number of cubes.
    // In order to match the voxels, we can just use NX-1, NY-1, and NZ-1
    int CX = NX-1;
    int CY = NY-1;
    int CZ = NZ-1;
    
    // Sides of each cube
    int DX = (NX-1)/CX;
    int DY = (NY-1)/CY;
    int DZ = (NZ-1)/CZ;

    // Array of cubes
    std::vector<std::vector<std::vector<Cube>>> cubes(CX, 
    std::vector<std::vector<Cube>>(CY, 
        std::vector<Cube>(CZ, Cube{})));

    for (int x = 0; x < CX; x++) {
        for (int y = 0; y < CY; y++) {
            for(int z = 0; z < CZ; z++){

                // Corners (checked against the surface value)
                for(int cor = 0; cor < 8; cor++){
                    cubes[x][y][z].corners[cor] = 
                        density[(x + offset[cor][0]) * DX]
                        [(y + offset[cor][1]) * DY]
                        [(z + offset[cor][2]) * DZ] > SURFACE_VALUE;
                }

                // Edges that contain an intersection
                // We will look at the reference table, specifically at
                // the index that matches the cube vertex configuration 
                // (turning the boolean array into a uint8).
                uint8_t configuration = 0;
                for (int i = 0; i < 8; i++){
                    configuration |= ((cubes[x][y][z].corners[i] ? 1 : 0) << i);
                }
                
                int edges = edge_table[configuration];
                for(int i = 0; i < 12; i++){
                    cubes[x][y][z].edges[i] = (edges >> i) & 1;
                }

                // Next we have to calculate the edge intersection positions
                // in world coordinates.
                for(int i = 0; i < 12; i++){

                    // if the edge has an intersection at all
                    if(!cubes[x][y][z].edges[i]){
                        continue;
                    }

                    if(USE_MIDDLE_POINT){
                        
                        // The two points of each edge
                        float x0 = (x + edge_shifts[i][0]) * DX;
                        float y0 = (y + edge_shifts[i][1]) * DY;
                        float z0 = (z + edge_shifts[i][2]) * DZ;
                        float x1 = (x + edge_shifts[i][3]) * DX;
                        float y1 = (y + edge_shifts[i][4]) * DY;
                        float z1 = (z + edge_shifts[i][5]) * DZ;

                        cubes[x][y][z].vertex_coordinate[i][0] = 0.5f * (x0 + x1);
                        cubes[x][y][z].vertex_coordinate[i][1] = 0.5f * (y0 + y1); 
                        cubes[x][y][z].vertex_coordinate[i][2] = 0.5f * (z0 + z1);    
                        
                        // Likewise, to calculate normal, we just use the
                        // middle point between the gradients on the two ends
                        // of the egde.
                        Vector3 normal;
                        normal = gradient[x0][y0][z0], gradient[x1][y1][z1];
                        normal = normal * 0.5;
                        normal.normalize();

                        cubes[x][y][z].vertex_normal[i][0] = normal.x();
                        cubes[x][y][z].vertex_normal[i][1] = normal.y();
                        cubes[x][y][z].vertex_normal[i][2] = normal.z();
                    }
                    else{

                    }
                    
                }

                // We can then determine which triangles are formed, using
                // the second lookup table.
                // There are at most 5 triangles.
                cubes[x][y][z].triangles.reserve(5);

                for(int i = 0; i < 16; i+=3){
                    if(tri_table[configuration][i] == -1){
                        break;
                    }
                    
                    cubes[x][y][z].triangles.push_back(std::vector<int>{
                        tri_table[configuration][i], 
                        tri_table[configuration][i+1],
                        tri_table[configuration][i+2]
                    });
                }
            }
        }
    }


    //=======================================================================//


    std::ofstream out("mesh.obj");

    // Vertices
    for (int x = 0; x < CX; x++) {
        for (int y = 0; y < CY; y++) {
            for(int z = 0; z < CZ; z++){

                for(int tri = 0; tri < cubes[x][y][z].triangles.size(); tri++){
                    // For each triangle vertex, we save the coordinates
                    for(int i = 0; i < 3; i++){
                        out << "v "
                            << cubes[x][y][z].vertex_coordinate
                                [cubes[x][y][z].triangles[tri][i]][0] << " "
                            << cubes[x][y][z].vertex_coordinate
                                [cubes[x][y][z].triangles[tri][i]][1] << " "
                            << cubes[x][y][z].vertex_coordinate
                                [cubes[x][y][z].triangles[tri][i]][2] << "\n";
                    }
                }
            }
        }
    }


    // Normals
    for (int x = 0; x < CX; x++) {
        for (int y = 0; y < CY; y++) {
            for(int z = 0; z < CZ; z++){

                for(int tri = 0; tri < cubes[x][y][z].triangles.size(); tri++){
                    // For each triangle vertex, we save the normals
                    for(int i = 0; i < 3; i++){
                        out << "vn "
                            << cubes[x][y][z].vertex_normal
                                [cubes[x][y][z].triangles[tri][i]][0] << " "
                            << cubes[x][y][z].vertex_normal
                                [cubes[x][y][z].triangles[tri][i]][1] << " "
                            << cubes[x][y][z].vertex_normal
                                [cubes[x][y][z].triangles[tri][i]][2] << "\n";
                    }
                }
            }
        }
    }
    
    // In a .obj, the indexes are 1-based
    int v_count = 1;

    for (int x = 0; x < CX; x++) {
        for (int y = 0; y < CY; y++) {
            for(int z = 0; z < CZ; z++){
                for(int tri = 0; tri < cubes[x][y][z].triangles.size(); tri++){
                    // Now we save the triangle faces
                    out << "f " << v_count << "//" << v_count << " " 
                                << v_count+1  << "//" << v_count + 1 << " " 
                                << v_count+2  << "//" << v_count + 2 << "\n";
                    v_count += 3;
                }
            }
        }
    }

    out.close();

    std::cout << "Terminated Successfully.\n";

    return 0;
}
