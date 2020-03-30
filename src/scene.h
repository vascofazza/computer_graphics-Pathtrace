#ifndef _SCENE_H_
#define _SCENE_H_

#include "common.h"
#include "json.h"
#include "vmath.h"
#include "image.h"
#include <random>

// forward declarations
struct BVHAccelerator;

// blinn-phong material
// textures are scaled by the respective coefficient and may be missing
struct Material {
    vec3f       ke = zero3f;        // emission coefficient
    vec3f       kd = one3f;         // diffuse coefficient
    vec3f       ks = zero3f;        // specular coefficient
    float       n = 10;             // specular exponent
    vec3f       kr = zero3f;        // reflection coefficient
    vec3f       translucency = zero3f;
    float       hair_radius = 0.0035f;
    float       t_factor = 65.f;
    
    float       blur = 0;           //enable blur reflection
    float       blur_size = 0.1;    //blur cone size
    
    int         hair_count = 0;    // number of hair to grow
    float       hair_length = 1.f; // length of hair
    
    float       ruffle = 0.65f;     // hair ruffle o.65f
    
    image3f*    ke_txt = nullptr;   // emission texture
    image3f*    kd_txt = nullptr;   // diffuse texture
    image3f*    ks_txt = nullptr;   // specular texture
    image3f*    kr_txt = nullptr;   // reflection texture
    image3f*    norm_txt = nullptr; // normal texture
    
    bool        double_sided = false;   // double-sided material
    bool        microfacet = false; // use microfacet formulation
};

// Keyframed Animation Data
struct FrameAnimation {
    frame3f                 rest_frame = identity_frame3f; // animation rest frame
    vector<int>             keytimes;   // key frame times
    vector<vec3f>           translation;// translation key frames
    vector<vec3f>           rotation;   // rotation key frames
};

// Mesh Skinning Data
struct MeshSkinning {
    vector<vec3f>           rest_pos;      // rest pos
    vector<vec3f>           rest_norm;     // rest norm
    vector<vec4i>           bone_ids;      // skin bones
    vector<vec4f>           bone_weights;  // skin weights
    vector<vector<mat4f>>   bone_xforms;   // bone xforms (bone index is the first index)
};

// Mesh Simulation Data
struct MeshSimulation {
    // simulation init data
    vector<vec3f>           init_pos;  // initial position
    vector<vec3f>           init_vel;  // initial velocity
    vector<float>           mass;      // mass
    vector<bool>            pinned;    // whether a point is pinned
    
    // inter-particle springs
    struct Spring {
        vec2i               ids;       // springs vertex ids
        float               restlength;// springs rest length
        float               ks;        // springs static constant
        float               kd;        // springs dynamic constant
    };
    vector<Spring>              springs;    // springs
    
    // simulation compute data
    vector<vec3f>           vel;       // velocity
    vector<vec3f>           force;     // forces
};

// Mesh Collision Data
struct MeshCollision {
    float                   radius;     // collision radius
    bool                    isquad;     // whether the collision object is a sphere or quad
};

// indexed mesh data structure with vertex positions and normals,
// a list of indices for triangle and quad faces, material and frame
struct Mesh {
    frame3f         frame = identity_frame3f;   // frame
    vector<vec3f>   pos;                        // vertex position
    vector<vec3f>   norm;                       // vertex normal
    vector<vec2f>   texcoord;                   // vertex texcture coordinates
    vector<vec3i>   triangle;                   // triangle
    vector<vec4i>   quad;                       // quad
    vector<int>     point;                      // point
    vector<vec2i>   line;                       // line
    vector<vec4i>   spline;                     // cubic bezier segments
    vector<vec3f>   in_radiance;                // entrance radiance per vertex
    Material*       mat = new Material();       // material
    Material*       hair_mat = nullptr;  // hair material
    float           area = 0;                   // area
    float           invert = 1.f;
    std::discrete_distribution<> distribution;
    
    int  subdivision_catmullclark_level = 0;        // catmullclark subdiv level
    bool subdivision_catmullclark_smooth = false;   // catmullclark subdiv smooth
    int  subdivision_bezier_level = 0;              // bezier subdiv level
    
    FrameAnimation* animation = nullptr;        // animation data
    MeshSkinning*   skinning = nullptr;         // skinning data
    MeshSimulation* simulation = nullptr;       // simulation data
    MeshCollision*  collision = nullptr;        // collision data
    
    BVHAccelerator* bvh = nullptr;              // bvh accelerator for intersection
    vector<int>     bvh_surfaces;               //boundigbox for surfaces
};

// surface made of eitehr a spehre or a quad (as determined by
// isquad. the sphere is centered frame.o with radius radius.
// the quad is at frame.o with normal frame.z and axes frame.x, frame.y.
// the quad side is 2*radius.
struct Surface {
    frame3f     frame = identity_frame3f;   // frame
    float       radius = 1;                 // radius
    bool        isquad = false;             // whether it's a quad
    float       height = -1;                 //shape height
    Material*   mat = new Material();       // material
    vec3f       tan = zero3f;
    size_t      id = -1;
    FrameAnimation* animation = nullptr;    // animation data

    Mesh*       _display_mesh = nullptr;    // display mesh
};

// point light at frame.o with intensity intensity
struct Light {
    frame3f     frame = identity_frame3f;       // frame
    vec3f       intensity = one3f;              // intersntiy
};

// perspective camera at frame.o with direction (-z)
// and image plane orientation (x,y); the image plane
// is at a distance dist with size (width,height);
// the camera is focussed at a distance focus.
struct Camera {
    frame3f frame = identity_frame3f;   // frame
    float   width = 1;                  // image plane width
    float   height = 1;                 // image plane height
    float   dist = 1;                   // image plane distance
    float   focus = 1;                  // distance of focus
    float   aperture = 0;               // aperture size
    bool    circle = false;             // type of lens
};

// Scene Animation Data
struct SceneAnimation {
    int     time = 0;                       // current animation time
    int     length = 0;                     // animation length
    float   dt = 1/30.0f;                   // time in seconds for each time step
    int     simsteps = 100;                 // simulation steps for time step of animation
    vec3f   gravity = {0,-9.8f,0};          // acceleration of gravity
    vec2f   bounce_dump = {0.001f,0.5f};    // loss of velocity at bounce (parallel,ortho)
    bool    gpu_skinning = false;           // whether to skin on the gpu
    bool    loop = false;                   // whether to loop the animation
};

// scene comprised of a camera, a list of meshes,
// and a list of lights. rendering parameters are
// also included, namely the background color (color
// if a ray misses) the ambient illumination, the
// image resolution (image_width, image_height) and
// the samples per pixel (image_samples).
struct Scene {
    Camera*             camera = new Camera();  // camera
    vector<Mesh*>       meshes;                 // meshes
    vector<Surface*>    surfaces;               // surfaces
    vector<Surface*>    surfaces_aux;           // surfaces accelerated
    vector<Light*>      lights;                 // lights
    
    vec3f               background = one3f*0.2; // background color
    image3f*            background_txt = nullptr;// background texture
    vec3f               ambient = one3f*0.2;    // ambient illumination

    SceneAnimation*     animation = new SceneAnimation();    // scene animation data

    int                 image_width = 512;      // image resolution in x
    int                 image_height = 512;     // image resolution in y
    int                 image_samples = 1;      // samples per pixels in each direction
    
    bool                draw_wireframe = false; // whether to use wireframe for interactive drawing
    bool                draw_animated = false;  // whether to draw with animation
    bool                draw_gpu_skinning = false;  // whether skinning is performed on the gpu
    bool                draw_captureimage = false;  // whether to capture the image in the next frame
    bool                draw_normals = false;       // whether to draw normals for debugging
    
    int                 path_max_depth = 2;     // maximum path depth
    bool                path_sample_brdf = true;// sample brdf in path tracing
    bool                path_shadows = true;    // whether to compute shadows
};

// grab all scene textures
vector<image3f*> get_textures(Scene* scene);

// create a Camera at eye, pointing towards center with up vector up, and with specified image plane params
Camera* lookat_camera(vec3f eye, vec3f center, vec3f up, float width, float height, float dist);

// set camera view with a "turntable" modification
void set_view_turntable(Camera* camera, float rotate_phi, float rotate_theta, float dolly, float pan_x, float pan_y);

// load a scene from a json file
Scene* load_json_scene(const string& filename);

Scene* load_obj_scene(const string& filename);

#endif

