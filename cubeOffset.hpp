
#ifndef CUBE_OFFSET_H
#define CUBE_OFFSET_H

#include "vector3.hpp"

// Cube offset to the corners from the center
// (recorded using vectors for convenience).
Vector3 cubeCenterOffset[8]{
    Vector3(-1, -1, -1),
    Vector3(1, -1, -1),
    Vector3(1, 1, -1),
    Vector3(-1, 1, -1),
    Vector3(-1, -1, 1),
    Vector3(1, -1, 1),
    Vector3(1, 1, 1),
    Vector3(-1, 1, 1)
};


// This is like the inverse of the above map, takes the
// offssets (from the corner) and returns the index
int inverseOffset[2][2][2] = {
    { {0, 4}, {3, 7} },
    { {1, 5}, {2, 6} }
};


// Cube offset to the corners from the min corner in x, y, z
int cubeCornerOffset[8][3]{
    {0, 0, 0}, 
    {1, 0, 0}, 
    {1, 1, 0}, 
    {0, 1, 0},
    {0, 0, 1}, 
    {1, 0, 1}, 
    {1, 1, 1}, 
    {0, 1, 1}
};


// Maps each edge index to the two vertices it connects
// (as offsets).
const int edgeOffset[12][6] = {
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

#endif