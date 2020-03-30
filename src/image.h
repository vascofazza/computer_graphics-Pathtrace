#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "common.h"
#include "vmath.h"

// A generic image
struct image3f {
    // Default Constructor (empty image)
	image3f() : _w(0), _h(0) { }
    // Size Constructor (sets width and height)
	image3f(int w, int h) : _w(w), _h(h), _d(_w*_h,vec3f()) { }
    // Size Constructor with initialization (sets width and height and initialize pixels)
	image3f(int w, int h, const vec3f& v) : _w(w), _h(h), _d(_w*_h,vec3f()) {
        for(auto &p : _d) p = v;
    }
    
    // image width
    int width() const { return _w; }
    // image height
    int height() const { return _h; }
	
    // element access
	vec3f& at(int i, int j) { return _d[j*_w+i]; }
    // element access
	const vec3f& at(int i, int j) const { return _d[j*_w+i]; }
    
    // data access
    vec3f* data() { return _d.data(); }
    // data access
    const vec3f* data() const { return _d.data(); }
    
    // flips this image along the y axis returning a new image
    image3f flipy() const {
        image3f ret(width(),height());
        for(int j = 0; j < height(); j ++) {
            for(int i = 0; i < width(); i ++) {
                ret.at(i,height()-1-j) = at(i,j);
            }
        }
        return ret;
    }
    
    // apply gamma correction
    image3f gamma(float gamma) const {
        image3f ret(width(),height());
        for(int j = 0; j < height(); j ++) {
            for(int i = 0; i < width(); i ++) {
                ret.at(i,j) = pow(at(i,j),gamma);
            }
        }
        return ret;
    }
    
    // apply a scale to the image
    image3f scale(float s) {
        image3f ret(width(),height());
        for(int j = 0; j < height(); j ++) {
            for(int i = 0; i < width(); i ++) {
                ret.at(i,j) = at(i,j)*s;
            }
        }
        return ret;
    }
    
private:
	int _w, _h;
	vector<vec3f> _d;
};

// Write an floating point color PFM image file
void write_pfm(const string& filename, const image3f& img, bool flipY = false);
// Write an 8-bit color compressed PNG file (sets PNG alpha to 1 everywhere)
void write_png(const string& filename, const image3f& img, bool flipY = false);

// Load a PFM or PPM color image and return it as a floating point color image
image3f read_pnm(const string& filename, bool flipY);
// Load a compressed PNG color image and return it as a floating point color image
image3f read_png(const string& filename, bool flipY);

#endif
