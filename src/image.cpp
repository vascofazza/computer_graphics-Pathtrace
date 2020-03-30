#include "image.h"
#include "lodepng.h"

static void _read_pnm(const string& filename, char& type,
               int& width, int& height, int& nc,
               float& scale, unsigned char*& buffer) {
	char identifier[4]; 
    bool ascii = false;
    
	FILE *f = fopen(filename.c_str(), "rb");
    if (f == 0) {
        error_if_not(f != 0, "failed to open image file %s", filename.c_str());
        nc = width = height = 0;
        buffer = nullptr;
        return;
    }
    
	error_if_not(fscanf(f, "%s", identifier) == 1, "error reading image file");
	string id = identifier;
	
	if("Pf" == id) { nc = 1; ascii = false; type = 'f'; }
	else if ("PF"  == id) { nc = 3; ascii = false; type = 'f'; }
	else if ("P2"  == id) { nc = 1; ascii = true; type = 'B'; }
	else if ("P3"  == id) { nc = 3; ascii = true; type = 'B'; }
	else if ("P5"  == id) { nc = 1; ascii = false; type = 'B'; }
	else if ("P6"  == id) { nc = 3; ascii = false; type = 'B'; }
	else {	error("unknown image format"); }
    
	error_if_not(fscanf(f, "%d%d\n", &width, &height) == 2, "error reading image file %s", filename.c_str());
    if(type == 'B') { 
        buffer = new unsigned char[width*height*nc]; 
        int maxv;
        error_if_not(fscanf(f, "%d\n", &maxv) == 1, "error reading image file %s", filename.c_str());
        error_if_not(maxv == 255, "unsupported max value");
        scale = 1.0f / maxv;
        if(!ascii) {
            int bread = (int)fread(buffer, sizeof(unsigned char), width*height*nc, f);
            error_if_not(bread == width*height*nc, "error reading image file %s", filename.c_str());
        } else {
            unsigned char* buf = (unsigned char*)buffer;
            for(int i = 0; i < width*height*nc; i ++) {
                int v;
                error_if_not(fscanf(f, "%d", &v) == 1, "error reading image file %s", filename.c_str());
                buf[i] = (unsigned char)v;
            }
        }
    } else if(type == 'f') { 
        buffer = (unsigned char*)(new float[width*height*nc]);
        char scale_string[16];
        fgets(scale_string, 16, f);
        sscanf(scale_string, "%f\n", &scale);
        error_if_not(scale < 0, "only support little endian pfm");
        scale = abs(scale);
        // flip y while reading
        for(int j = height-1; j >= 0; j --) {
            float* buf = (float*)buffer;
            int bytes = (int)fread(buf + (j*width*nc), sizeof(float), width*nc, f);
            error_if_not(bytes == width*nc, "error reading image file");
        }
    } else error("unknown type");
    
    fclose(f);
}

image3f read_pnm(const string& filename, bool flipY) {
    int width, height, nc; float scale; unsigned char* buffer; char type;
    _read_pnm(filename, type, width, height, nc, scale, buffer);
    if (not buffer) {
        error_if_not(false, "failed to load image file %s", filename.c_str());
        return image3f();
    }
    error_if_not(nc == 3 && (type == 'f' || type == 'B'), "unsupported image format in file %s", filename.c_str());
    
    image3f img(width,height);
    if(type == 'f') {
        float* buf = (float*)buffer;
        for(int i = 0; i < img.width()*img.height(); i ++) {
            img.data()[i] = vec3f(buf[i*3+0],buf[i*3+1],buf[i*3+2]) * scale;
        }
    } else if(type == 'B') {
        unsigned char* buf = (unsigned char*)buffer;
        for(int i = 0; i < img.width()*img.height(); i ++) {
            img.data()[i] = vec3f((float)buf[i*3+0],(float)buf[i*3+1],(float)buf[i*3+2]) * scale;
        }        
    }
    if (buffer) delete [] buffer;
    
    if(flipY) img = img.flipy();
    
    return img;
}

static void _write_pnm(const char *filename, char type,
                         int width, int height, int nc,
                         bool ascii, unsigned char* buffer) {
	FILE *f = fopen(filename, "wb");
	error_if_not(f != 0, "failed to create image file %s", filename);
    
    string magic; 
    int scale = -1; 
    int ds = 0;
    if(type == 'f' && nc == 1 && !ascii) { magic = "Pf"; scale = -1; 
        ds = sizeof(float); }
	else if(type == 'f' && nc == 3 && !ascii) { magic = "PF"; scale = -1; 
        ds = sizeof(float); }
	else if(type == 'B' && nc == 1 && !ascii) { magic = "P5"; scale = 255; 
        ds = sizeof(unsigned char); }
	else if(type == 'B' && nc == 3 && !ascii) { magic = "P6"; scale = 255; 
        ds = sizeof(unsigned char); }
	else if(type == 'B' && nc == 1 && ascii) { magic = "P2"; scale = 255; 
        ds = sizeof(unsigned char); }
	else if(type == 'B' && nc == 3 && ascii) { magic = "P3"; scale = 255; 
        ds = sizeof(unsigned char); }
    else error("unsupported image format");
    
    error_if_not(fprintf(f, "%s\n", magic.c_str()) > 0, "error writing file %s", filename);
    error_if_not(fprintf(f, "%d %d\n", width, height) > 0, "error writing file %s", filename);
    error_if_not(fprintf(f, "%d\n", scale) > 0, "error writing file %s", filename);
    
    if(!ascii) {
        if(type == 'f') {
            for(int j = height-1; j >= 0; j --) {
                float* buf = (float*)buffer;
                error_if_not((int)fwrite(buf + j*width*nc, ds, width*nc, f) == width*nc, "error writing file %s", filename);
            }
        } else {
            error_if_not((int)fwrite(buffer, ds, width*height*nc, f) == width*height*nc, "error writing file %s", filename);
        }
    } else {
        unsigned char* buf = (unsigned char*)buffer;
        for(int i = 0; i < width*height*nc; i ++) {
            int v = buf[i];
            error_if_not(fprintf(f, "%d \n", v) != 0, "error writing file %s", filename);
        }
    }
    
	fclose(f);
}

void write_pfm(const string& filename, const image3f& img, bool flipY) {
    _write_pnm(filename.c_str(), 'f', img.width(), img.height(), 3, false,
               (flipY)?(unsigned char*)img.flipy().data():(unsigned char*)img.data());
}

image3f read_png(const string& filename, bool flipY) {
    vector<unsigned char> pixels;
    unsigned width, height;
	
    unsigned error = lodepng::decode(pixels, width, height, filename);
    error_if_not(not error,"cannot read png image: %s", filename.c_str());
	
    error_if_not(pixels.size() == width*height*4, "bad reading");
    
	image3f img(width,height);
	for(int i = 0; i < width*height; i ++) {
        int x = i%width;
        int y = i/width;
        int ii = (flipY) ? x + (height-y-1)*width : x + y*width;
		img.data()[ii].x = pixels[i*4+0] / 255.0f;
		img.data()[ii].y = pixels[i*4+1] / 255.0f;
		img.data()[ii].z = pixels[i*4+2] / 255.0f;
	}
	
	return img;
}

void write_png(const string& filename, const image3f& img, bool flipY) {
    vector<unsigned char> img_png(img.width()*img.height()*4);
    for(int x = 0; x < img.width(); x++ ) {
        for( int y = 0; y < img.height(); y++ ) {
            int i_img = y * img.width() + x;
            int i_png = ( ( flipY ? (img.height()-1-y) : y ) * img.width() + x ) * 4;
            img_png[i_png+0] = (unsigned char)clamp(img.data()[i_img].x * 255, 0.0f, 255.0f);
            img_png[i_png+1] = (unsigned char)clamp(img.data()[i_img].y * 255, 0.0f, 255.0f);
            img_png[i_png+2] = (unsigned char)clamp(img.data()[i_img].z * 255, 0.0f, 255.0f);
            img_png[i_png+3] = 255;
        }
    }
    unsigned error = lodepng::encode(filename, img_png, img.width(), img.height());
    error_if_not(not error, "cannot write png image: %s", filename.c_str());
}
