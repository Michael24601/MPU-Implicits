
#ifndef MARCHING_CUBES_H
#define MARCHING_CUBES_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <functional>

#include "vector3.hpp"
#include "lookupTables.hpp"
#include "cubeOffset.hpp"


// This is the data associated with each cube in marching cubes.
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
    Vector3 vertexCoordinate[12];
    
    // Each intersection point's normal (assumed to be normalized)
    Vector3 vertexNormal[12];
};


// Instead of expecting data, it just expects us to send it a function
// object that satisfies some basic assumptions, a surface value b
// such that f(x, y, z) = b defines our surface, the range
// of points, the number of desired cubes in each direction
// (here nx, ny, and nz is the number of cubes + 1), and whether to
// use the middle point or interpolate to get the triangle vertices
// along the edges.
std::vector<std::vector<std::vector<Cube>>> marchingCubes(
    const std::function<real(const Vector3&)>& function,
    const Vector3& range, 
    real surfaceValue, 
    int nx, int ny, int nz,
    bool useMiddlePoints
) {


    //================================ DATA =================================//


    // This boolean indicates whether or not triangle vertices are
    // calculated by linearly interpolating the surface and intersecting
    // the cube edges, or just taking the edge middle point. 
    const bool USE_MIDDLE_POINT = useMiddlePoints; 

    // If we already know how many cubes we want to use,
    // as well as our support, then we can precompute some values
    // for efficiency, including f(x, y, z) and nabla f(x, y, z),
    // the gradient.

    // Number of points.
    int NX = nx;
    int NY = ny;
    int NZ = nz;

    // Number of cubes.
    int CX = NX-1;
    int CY = NY-1;
    int CZ = NZ-1;
    
    // Sides of each cube
    real DX = (range.x())/CX;
    real DY = (range.y())/CY;
    real DZ = (range.z())/CZ;

    // Our 3D data is described by a desnity function f(x, y, z).
    // The original data could be a signed distance function, a
    // collection of scans, etc... It doesn't matter. This function
    // expects the density to be in the format of a function that
    // can be evaluated at any point between (0, 0, 0)
    // and the range point, so a cuboid.
    // Instead of reevaluating the function over and over again
    // during the algorithm, we can store the values of the function
    // at the cube vertices.
    std::vector<std::vector<std::vector<real>>> density(NX, 
        std::vector<std::vector<real>>(NY, 
            std::vector<real>(NZ)));

    // We will also calculate the gradient of the above scalar field
    // at each of the voxels. This can be approximated if the data
    // is discrete (as was the case in the original paper),
    // where the density comprised of several 2D scans that moved
    // along the z axis.
    std::vector<std::vector<std::vector<Vector3>>> gradient(NX, 
        std::vector<std::vector<Vector3>>(NY, 
        std::vector<Vector3>(NZ)));

    // The surface defined by f(x, y, z) = suruface value is essentially 
    // just the surface that seperates the outside of the captured object 
    // from the inside.
    // In the case of a signed distance function, this is just 0.
    const real SURFACE_VALUE = surfaceValue;

    // Now to populate our function.
    for (int z = 0; z < NZ; z++) {
        for (int y = 0; y < NY; y++) {
            for (int x = 0; x < NX; x++) {
                density[x][y][z] = 
                    function(Vector3(x * DX, y * DY, z * DZ));
            }
        }
    }


    // Now we approximate the gradient (easier than differentiating
    // f(x, y, z), and fine as far as marching cubes goes.
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
                    (density[x_p1][y][z] - density[x_m1][y][z]) 
                        / ((x_p1 - x_m1) * DX),
                    (density[x][y_p1][z] - density[x][y_m1][z]) 
                        / ((y_p1 - y_m1) * DY),
                    (density[x][y][z_p1] - density[x][y][z_m1]) 
                        / ((z_p1 - z_m1) * DZ)
                );

                gradient[x][y][z].normalize();
            }
        }
    }


    //=========================== MARCHING CUBES ============================//


    // Based on whether their vertices are inside or outside
    // the surface (intensity at vertex below or above the surface value),
    // we can tell how (if at all) the surface intersets the cube,
    // which in turn tells how to generate a triangle mesh.

    // For a less accurate mesh, we can use less cubes.

    // Array of cubes
    std::vector<std::vector<std::vector<Cube>>> cubes(CX, 
    std::vector<std::vector<Cube>>(CY, 
        std::vector<Cube>(CZ, Cube{})));

    for (int x = 0; x < CX; x++) {
        for (int y = 0; y < CY; y++) {
            for(int z = 0; z < CZ; z++){

                // Here we populate the corners array, by checking
                // if the corners are inside or outside the surface.
                for(int cor = 0; cor < 8; cor++){
                    cubes[x][y][z].corners[cor] = 
                        density[x + cubeCornerOffset[cor][0]]
                        [y + cubeCornerOffset[cor][1]]
                        [z + cubeCornerOffset[cor][2]] > SURFACE_VALUE;
                }

                // The first table told us which edges will have points
                // on them based on the configuration of points inside
                // and outside the surface.
                // We will look at the reference table, specifically at
                // the index that matches the cube vertex configuration.
                // Since we have 8 corners, we can store this in a uint8,
                // and use that as our array index.
                uint8_t configuration = 0;
                for (int i = 0; i < 8; i++){
                    configuration |= ((cubes[x][y][z].corners[i] ? 1 : 0) << i);
                }
                
                int edges = EDGE_TABLE[configuration];
                for(int i = 0; i < 12; i++){
                    // The table itself stores the 12 booleans inside
                    // a 12 byte integer, so we use that to get the
                    // status of each edge.
                    cubes[x][y][z].edges[i] = (edges >> i) & 1;
                }

                // Next we have to calculate the edge intersection positions
                // in world coordinates.
                for(int i = 0; i < 12; i++){

                    // if the edge has an intersection at all
                    if(!cubes[x][y][z].edges[i]){
                        continue;
                    }

                    // The two points of each edge
                    int x0 = (x + edgeOffset[i][0]);
                    int y0 = (y + edgeOffset[i][1]);
                    int z0 = (z + edgeOffset[i][2]);
                    int x1 = (x + edgeOffset[i][3]);
                    int y1 = (y + edgeOffset[i][4]);
                    int z1 = (z + edgeOffset[i][5]);

                    // This is the interpolation variable along the edge
                    real alpha;

                    if(USE_MIDDLE_POINT){
                        // In this case we always pick the middle point
                        // of the edge as the triangle vertex point/
                        alpha = 0.5;
                    }
                    else{
                        // In this case, we linearly interpolate
                        // between both edge points to get a better
                        // approximation of the intersection point along
                        // the edge.
                        // We know one is positive and one is negative
                        // (at least relative to the surface value).
                        alpha = density[x1][y1][z1] - SURFACE_VALUE;
                        alpha /= (density[x1][y1][z1] - density[x0][y0][z0]);
                    }

                    // The edge corners in world coordinates
                    Vector3 p0(x0 * DX, y0 * DY, z0 * DZ);
                    Vector3 p1(x1 * DX, y1 * DY, z1 * DZ);

                    // The vertex along the edge, using the alpha
                    // we computed.
                    cubes[x][y][z].vertexCoordinate[i] = p0 * alpha 
                        + p1 * (1.0 - alpha);
                    
                    // Likewise, to calculate normal, we just use the
                    // same alpha.
                    Vector3 normal = gradient[x0][y0][z0] * alpha + 
                        gradient[x1][y1][z1] * (1.0 - alpha);
                    normal.normalize();

                    cubes[x][y][z].vertexNormal[i] = normal;
                    
                }

                // We can then determine which triangles are formed, using
                // the second lookup table.
                // The second table takes the configuration of points
                // that are inside and outside as input, and returns
                // the indexes of the points in triplets such that these
                // form triangles.
                cubes[x][y][z].triangles.reserve(5);

                // There are at most 5 triangles, so 15 total indexes
                // to loop through.
                for(int i = 0; i < 16; i+=3){

                    // Once the configuration starts returning -1,
                    // it means there are no longer triangles.
                    if(TRI_TABLE[configuration][i] == -1){
                        break;
                    }
                    
                    cubes[x][y][z].triangles.push_back(std::vector<int>{
                        TRI_TABLE[configuration][i], 
                        TRI_TABLE[configuration][i+1],
                        TRI_TABLE[configuration][i+2]
                    });
                }
            }
        }
    }

    return cubes;
}


// This function just takes the array of cubes containing some triangles
// and converts them to a .obj file.


void writeToFile(
    std::vector<std::vector<std::vector<Cube>>> cubes,
    std::string filename
){

    std::ofstream out(filename);

    if(cubes.empty() || cubes[0].empty()){
        throw std::invalid_argument("The cubes array is empty");
    }
    
    // Number of cubes.
    int CX = cubes.size();
    int CY = cubes[0].size();
    int CZ = cubes[0][0].size();

    // Vertices
    for (int x = 0; x < CX; x++) {
        for (int y = 0; y < CY; y++) {
            for(int z = 0; z < CZ; z++){

                for(int tri = 0; tri < cubes[x][y][z].triangles.size(); tri++){
                    // For each triangle vertex, we save the coordinates
                    for(int i = 0; i < 3; i++){
                        out << "v "
                            << cubes[x][y][z].vertexCoordinate
                                [cubes[x][y][z].triangles[tri][i]][0] << " "
                            << cubes[x][y][z].vertexCoordinate
                                [cubes[x][y][z].triangles[tri][i]][1] << " "
                            << cubes[x][y][z].vertexCoordinate
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
                            << cubes[x][y][z].vertexNormal
                                [cubes[x][y][z].triangles[tri][i]][0] << " "
                            << cubes[x][y][z].vertexNormal
                                [cubes[x][y][z].triangles[tri][i]][1] << " "
                            << cubes[x][y][z].vertexNormal
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
}


#endif