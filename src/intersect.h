#ifndef _INTERSECT_H_
#define _INTERSECT_H_

#include "scene.h"

// intersection record
struct intersection3f {
    bool        hit;        // whether it hits something
    float       ray_t;      // ray parameter for the hit
    vec3f       pos;        // hit position
    vec3f       norm;       // hit normal
    vec3f       tan = zero3f;
    vec2f       texcoord;   // hit texture coordinates
    Material*   mat;        // hit material
    vec3f       in_translucency = zero3f;
    
    // constructor (defaults to no intersection)
    intersection3f() : hit(false) { }
    
    // constructor to override default intersection
    explicit intersection3f(bool hit) : hit(hit) { }
};

#define ray3f_epsilon 0.0005f
#define ray3f_rayinf 1000000.0f

// 3D Ray
struct ray3f {
    vec3f e;        // origin
    vec3f d;        // direction
    float tmin;     // min t value
    float tmax;     // max t value
    
    // Default constructor
    ray3f() : e(zero3f), d(z3f), tmin(ray3f_epsilon), tmax(ray3f_rayinf) { }
    
    // Element-wise constructor
    ray3f(const vec3f& e, const vec3f& d) :
    e(e), d(d), tmin(ray3f_epsilon), tmax(ray3f_rayinf) { }
    
    // Element-wise constructor
    ray3f(const vec3f& e, const vec3f& d, float tmin, float tmax) :
    e(e), d(d), tmin(tmin), tmax(tmax) { }
    
    // Eval ray at a specific t
    vec3f eval(float t) const { return e + d * t; }
    
    // Create a ray from a segment
    static ray3f make_segment(const vec3f& a, const vec3f& b) { return ray3f(a,normalize(b-a),ray3f_epsilon,dist(a,b)-2*ray3f_epsilon); }
};

// transform a ray by a frame
inline ray3f transform_ray(const frame3f& f, const ray3f& v) { return ray3f(transform_point(f,v.e), transform_vector(f,v.d), v.tmin, v.tmax); }
// transform a ray by a frame inverse
inline ray3f transform_ray_inverse(const frame3f& f, const ray3f& v) { return ray3f(transform_point_inverse(f,v.e),transform_vector_inverse(f,v.d),v.tmin,v.tmax); }

void triangulate(Scene* scene);

// prepare scene acceleration and triangulate meshes
void accelerate(Scene* scene);

// intersects the scene and return the first intrerseciton
intersection3f intersect(Scene* scene, ray3f ray);

// intersects the scene and return any intrerseciton
bool intersect_shadow(Scene* scene, ray3f ray);

// intersects the scene's surfaces and return the first intrerseciton (used for raytracing homework)
intersection3f intersect_surfaces(Scene* scene, ray3f ray);

#endif
