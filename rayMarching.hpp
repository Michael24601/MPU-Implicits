
#ifndef RAY_MARCHING_H
#define RAY_MARCHING_H

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"
#include <functional>
#include "vector3.hpp"
#include "vector2.hpp"

class RayMarching{
  
private:

    struct Ray{
        // Assumed to be normalized
        Vector3 direction;
        Vector3 origin;

        Ray(const Vector3& direction, const Vector3& origin) :
            direction{direction}, origin{origin}{}

        // Returns point this far along the point
        Vector3 getPoint(real t){
            return origin + direction * t;
        }
    };

    struct PerspectiveCamera{
        // Stores the vectors from the center of the pixel
        // grid to the border along x and y.
        Vector3 halfWidth;
        Vector3 halfHeight;
        Vector3 origin;
        // Can give us the position of the aperture
        real focalLength;

        PerspectiveCamera(const Vector3& halfWidth, 
            const Vector3& halfHeight, const Vector3& origin, 
            real focalLength) : 
            halfWidth{halfWidth}, halfHeight{halfHeight}, origin{origin},
            focalLength{focalLength}{}

        
        // Given a point between (-1, -1) and (1, 1), shoots a ray
        // from the aperture towards said point onto the pixel
        // at that point on the pixel grid.
        Ray shootRay(const Vector2& point) const{

            assert((point.x() >= -1.0 && point.x() <= 1.0 
                && point.y() >= -1.0 && point.y() <= 1.0)
                && "The given point is not normalized\n");

            // Cross product, normalized
            Vector3 forwardDirection = (halfWidth % halfHeight).normalized();
            Vector3 aperture = origin - forwardDirection * focalLength;
            
            Vector3 pixelPoint = origin + halfWidth * point.x() 
                + halfHeight * point.y();
            Vector3 direction = (pixelPoint - aperture).normalized();

            Ray ray(direction, pixelPoint);
            return ray;
        }
    };

public:


    // Controls ray marcher step size (will be scaled by the scene)
    // This needs to be small to avoid staircase artifcats, where
    // large portions of the rendered image look like a solid color.
    static real RAY_MARCHER_STEP;
    static real FINITE_DIFFERENCE_STEP;


    // Takes as input an implicit function f(x, y, z) = 0 that
    // defines some surface, some camera parameters, and then
    // outputs an image rendered from that perspective.
    // We expect f to return a negative value inside the surface,
    // and a positive one outisde, but it doesn't have to be a
    // signed distance function.
    // We also want the main diagonal length of the cuboid containing
    // the entire surface as well as the centroid, 
    // so we get a sense for the scale.
    // For visualization purposes, this function also takes as input
    // a function that returns a color for each input (RGB, from 0 to 255).
    // If we don't want any colors, we can send a function that returns
    // pure white.
    static void run(
        const std::function<real(const Vector3&)>& function,
        const std::function<Vector3(const Vector3&)>& color,
        const Vector3& halfWidth, const Vector3& halfHeight, 
        const Vector3& origin, real focalLength,
        const Vector3& centroid, real scale,
        int imageWidth, int imageHeight,
        const Vector3& backgroundColor,
        std::string filename
    ){

        PerspectiveCamera camera(halfWidth, halfHeight, origin, focalLength);
        
        // A rough estimate of how far a ray might have to march before
        // it exits the scene entirely.
        real maxDistance = (centroid - origin).length() + scale;
        real dt = RAY_MARCHER_STEP * maxDistance;

        std::vector<unsigned char> imageBuffer(imageWidth * imageHeight * 3);

        real w = static_cast<real>(imageWidth);
        real h = static_cast<real>(imageHeight);
        for(int i = 0; i < imageWidth; i++){
            // We want to map it such that pixel i is at 
            // the center of the first pixel, and pixel W-1 is at the
            // center of the last pixel.
            // Each pixel is 1/W wide.
            real x = (i + 0.5) / w * 2.0 - 1.0;
            for(int j = 0; j < imageHeight; j++){

                // Same for y
                real y = (j + 0.5) / h * 2.0 - 1.0;

                // We can then shoot the ray
                Ray ray = camera.shootRay(Vector2(x, y));
                // How far along the ray we've marched
                real t = 0;
                bool intersected = false;
                while(t < maxDistance && !intersected){
                    real value = function(ray.getPoint(t));
                    if(value < 0){
                        intersected = true;
                    }
                    else{
                        // Having a fixed step size leads to incredibly
                        // slow rendering time, however, for the purposes
                        // of visualisation (not in real time) it will
                        // have to do.
                        t += dt;
                    }
                }
                
                // Flattened index (note that we flip the height
                // so it turns out correct in the PNG).
                int index = ((imageHeight - 1 - j) * imageWidth + i) * 3;

                // If we do intersect a point, we can color it
                // however we want.
                if(intersected){
                    // We can use finite differences to get the gradient.
                    // For simplicity,  we can use the same h for all
                    // directions.
                    // We multiply it by the diagonal to fit the scale.
                    real h = (FINITE_DIFFERENCE_STEP * scale);
                    real denom = 1.0 / (2 * h);
                    Vector3 gradient;
                    Vector3 p = ray.getPoint(t);
                    gradient.setX((function(p + Vector3(h, 0, 0)) 
                        - function(p - Vector3(h, 0, 0))) / denom);
                    gradient.setY((function(p + Vector3(0, h, 0)) 
                        - function(p - Vector3(0, h, 0))) / denom);
                    gradient.setZ((function(p + Vector3(0, 0, h)) 
                        - function(p - Vector3(0, 0, h))) / denom);
                    gradient.normalize();
                    
                    real brightness = -gradient * ray.direction;
                    brightness = std::max(0.0, brightness);
                    // We can get the color from the function,
                    // and then scale using the brightness
                    Vector3 pointColor = color(p) * brightness;

                    imageBuffer[index] = 
                        static_cast<unsigned char>(pointColor.x());
                    imageBuffer[index + 1] = 
                        static_cast<unsigned char>(pointColor.y());
                    imageBuffer[index + 2] = 
                        static_cast<unsigned char>(pointColor.z());
                }
                // Otherwise we just keep the background color
                else{
                    imageBuffer[index] = 
                        static_cast<unsigned char>(backgroundColor.x());
                    imageBuffer[index + 1] = 
                        static_cast<unsigned char>(backgroundColor.y());
                    imageBuffer[index + 2] = 
                        static_cast<unsigned char>(backgroundColor.z());
                }
            }
        }

        // Finally, we can create the PNG
        stbi_write_png(filename.c_str(), imageWidth, imageHeight, 3, 
            imageBuffer.data(), imageWidth * 3);

    }

};

#endif