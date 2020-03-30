#ifndef _MONTECARLO_H_
#define _MONTECARLO_H_

#include "vmath.h"

#include <random>

// Random number generator
struct Rng {
    // C++11 random number engine
    std::minstd_rand                        engine;
    
    // Seed the generator
    void seed(unsigned int seed) { engine.seed(seed); }
	
    // Generate a float in [0,1)
	float next_float() { return std::uniform_real_distribution<float>(0,1)(engine); }
    // Generate a float in [v.x,v.y)
	float next_float(const vec2f& v) { return std::uniform_real_distribution<float>(v.x,v.y)(engine); }
    
	// Generate 2 floats in [0,1)^2
    vec2f next_vec2f() { return vec2f(next_float(),next_float()); }
	// Generate 3 floats in [0,1)^3
    vec3f next_vec3f() { return vec3f(next_float(),next_float(),next_float()); }
    
    // Generate an int in [v.x,v.y)
    int next_int(const vec2i& v) { return std::uniform_int_distribution<int>(v.x,v.y)(engine); }
    
    // Create and seed nrngs generators
    static std::vector<Rng> generate_seeded(int nrngs) {
        std::seed_seq sseq{0,1,2,3,4,5,6,7,8,9};
        auto seeds = std::vector<int>(nrngs);
        sseq.generate(seeds.begin(), seeds.end());
        auto rngs = std::vector<Rng>(nrngs);
        for(int i = 0; i < nrngs; i ++) {
            rngs[i].seed(seeds[i]);
        }
        return rngs;
    }
};

// A set of per-pixel randon number generators seeded automatically
struct RngImage {
    // Default constructor
    RngImage() : _w(0), _h(0) { }
    // Size Constructor (sets width and height)
	RngImage(int w, int h) : _w(w), _h(h), _d(Rng::generate_seeded(w*h)) { }
    
    // element access
	Rng& at(int i, int j) { return _d[j*_w+i]; }
    
private:
	int _w, _h;
	vector<Rng> _d;
};

// hemispherical direction with uniform distribution
inline vec3f sample_direction_hemispherical_uniform(const vec2f& ruv) {
    auto z = ruv.y;
    auto r = sqrt(1-z*z);
    auto phi = 2 * pi * ruv.x;
    return vec3f(r*cos(phi), r*sin(phi), z);
}

// pdf for hemispherical direction with uniform distribution
inline float sample_direction_hemispherical_uniform_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : 1/(2*pi);
}

// spherical direction with uniform distribution
inline vec3f sample_direction_spherical_uniform(const vec2f ruv) {
    auto z = 2*ruv.y-1;
    auto r = sqrt(1-z*z);
    auto phi = 2 * pi * ruv.x;
    return vec3f(r*cos(phi), r*sin(phi), z);
}

// pdf for spherical direction with uniform distribution
inline float sample_direction_spherical_uniform_pdf(const vec3f& w) {
    return 1/(4*pi);
}

// hemispherical direction with cosine distribution
inline vec3f sample_direction_hemispherical_cosine(const vec2f& ruv) {
    auto z = sqrt(ruv.y);
    auto r = sqrt(1-z*z);
    auto phi = 2 * pi * ruv.x;
    return vec3f(r*cos(phi), r*sin(phi), z);
}

// pdf for hemispherical direction with cosine distribution
inline float sample_direction_hemispherical_cosine_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : w.z/pi;
}

// hemispherical direction with cosine power distribution
inline vec3f sample_direction_hemispherical_cospower(const vec2f& ruv, float n) {
    auto z = pow(ruv.y,1/(n+1));
    auto r = sqrt(1-z*z);
    auto phi = 2 * pi * ruv.x;
    return vec3f(r*cos(phi), r*sin(phi), z);
}

// pdf for hemispherical direction with cosine power distribution
inline float sample_direction_hemispherical_cospower_pdf(const vec3f& w, float n) {
    return (w.z <= 0) ? 0 : pow(w.z,n) * (n+1) / (2*pi);
}

// index with uniform distribution
inline int sample_index_uniform(float r, int size) {
    return clamp((int)(r * size), 0, size-1);
}

// pdf for index with uniform distribution
inline float sample_index_uniform_pdf(int size) {
    return 1.0f / size;
}

// computes the sample number in each dimension for stratified sampling
inline int sample_stratify_samplesnumber(int samples) {
    return (int)round(sqrt(samples));
}

// computes a sample for stratified sampling
inline vec2f sample_stratify_sample(const vec2f& uv, int sample, int samples_x, int samples_y) {
    int sample_x = sample % samples_x;
    int sample_y = sample / samples_x;
    return vec2f((sample_x + uv.x / samples_x), (sample_y + uv.y / samples_y));
}

// power distribution heuristics
inline float sample_power_heuristics(float fPdf, float gPdf) {
    return (fPdf*fPdf) / (fPdf*fPdf + gPdf*gPdf);
}

#endif
