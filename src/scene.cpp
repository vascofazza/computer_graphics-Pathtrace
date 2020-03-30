#include "scene.h"
#include "tesselation.h"
#include "yocto_obj.h"

vector<image3f*> get_textures(Scene* scene) {
    auto textures = set<image3f*>();
    for(auto mesh : scene->meshes) {
        if(mesh->mat->ke_txt) textures.insert(mesh->mat->ke_txt);
        if(mesh->mat->kd_txt) textures.insert(mesh->mat->kd_txt);
        if(mesh->mat->ks_txt) textures.insert(mesh->mat->ks_txt);
        if(mesh->mat->norm_txt) textures.insert(mesh->mat->norm_txt);
    }
    for(auto surface : scene->surfaces) {
        if(surface->mat->ke_txt) textures.insert(surface->mat->ke_txt);
        if(surface->mat->kd_txt) textures.insert(surface->mat->kd_txt);
        if(surface->mat->ks_txt) textures.insert(surface->mat->ks_txt);
        if(surface->mat->norm_txt) textures.insert(surface->mat->norm_txt);
    }
    if(scene->background_txt) textures.insert(scene->background_txt);
    return vector<image3f*>(textures.begin(),textures.end());
}

Camera* lookat_camera(vec3f eye, vec3f center, vec3f up, float width, float height, float dist) {
    auto camera = new Camera();
    camera->frame = lookat_frame(eye, center, up, true);
    camera->width = width;
    camera->height = height;
    camera->dist = dist;
    camera->focus = length(eye-center);
    return camera;
}

void set_view_turntable(Camera* camera, float rotate_phi, float rotate_theta, float dolly, float pan_x, float pan_y) {
    auto phi = atan2(camera->frame.z.z,camera->frame.z.x) + rotate_phi;
    auto theta = clamp(acos(camera->frame.z.y) + rotate_theta, 0.001f,pif-0.001f);
    auto new_z = vec3f(sin(theta)*cos(phi),cos(theta),sin(theta)*sin(phi));
    auto new_center = camera->frame.o-camera->frame.z*camera->focus;
    auto new_o = new_center + new_z * camera->focus;
    camera->frame = lookat_frame(new_o,new_center,y3f,true);
    camera->focus = dist(new_o,new_center);
    
    // apply dolly
    vec3f c = camera->frame.o - camera->frame.z * camera->focus;
    camera->focus = max(camera->focus+dolly,0.00001f);
    camera->frame.o = c + camera->frame.z * camera->focus;
    
    // apply pan
    camera->frame.o += camera->frame.x * pan_x + camera->frame.y * pan_y;
}

void json_set_values(const jsonvalue& json, float* value, int n) {
    error_if_not(n == json.array_size(), "incorrect array size");
    for(auto i : range(n)) value[i] = json.array_element(i).as_float();
}
void json_set_values(const jsonvalue& json, int* value, int n) {
    error_if_not(n == json.array_size(), "incorrect array size, %d,%d", n,json.array_size());
    for(auto i : range(n)) value[i] = json.array_element(i).as_int();
}

void json_set_value(const jsonvalue& json, bool& value) { value = json.as_bool(); }
void json_set_value(const jsonvalue& json, int& value) { value = json.as_int(); }
void json_set_value(const jsonvalue& json, float& value) { value = json.as_double(); }
void json_set_value(const jsonvalue& json, vec2f& value) { json_set_values(json, &value.x, 2); }
void json_set_value(const jsonvalue& json, vec3f& value) { json_set_values(json, &value.x, 3); }
void json_set_value(const jsonvalue& json, vec4f& value) { json_set_values(json, &value.x, 4); }
void json_set_value(const jsonvalue& json, vec2i& value) { json_set_values(json, &value.x, 2); }
void json_set_value(const jsonvalue& json, vec3i& value) { json_set_values(json, &value.x, 3); }
void json_set_value(const jsonvalue& json, vec4i& value) { json_set_values(json, &value.x, 4); }
void json_set_value(const jsonvalue& json, vector<bool>& value) {
    value.resize(json.array_size());
    for(auto i : range(value.size())) value[i] = json.array_element(i).as_bool();
}
void json_set_value(const jsonvalue& json, vector<float>& value) { value.resize(json.array_size()); json_set_values(json, &value[0], value.size()); }
void json_set_value(const jsonvalue& json, vector<vec2f>& value) { value.resize(json.array_size()/2); json_set_values(json, &value[0].x, value.size()*2); }
void json_set_value(const jsonvalue& json, vector<vec3f>& value) { value.resize(json.array_size()/3); json_set_values(json, &value[0].x, value.size()*3); }
void json_set_value(const jsonvalue& json, vector<vec4f>& value) { value.resize(json.array_size()/4); json_set_values(json, &value[0].x, value.size()*4); }
void json_set_value(const jsonvalue& json, vector<int>& value) { value.resize(json.array_size()); json_set_values(json, &value[0], value.size()); }
void json_set_value(const jsonvalue& json, vector<vec2i>& value) { value.resize(json.array_size()/2); json_set_values(json, &value[0].x, value.size()*2); }
void json_set_value(const jsonvalue& json, vector<vec3i>& value) { value.resize(json.array_size()/3); json_set_values(json, &value[0].x, value.size()*3); }
void json_set_value(const jsonvalue& json, vector<vec4i>& value) { value.resize(json.array_size()/4); json_set_values(json, &value[0].x, value.size()*4); }
void json_set_value(const jsonvalue& json, vector<mat4f>& value) { value.resize(json.array_size()/16); json_set_values(json, &value[0].x.x, value.size()*16); }
void json_set_value(const jsonvalue& json, frame3f& value) {
    value = identity_frame3f;
    if(json.object_contains("from") or json.object_contains("to") or json.object_contains("up")) {
        auto from = z3f, to = zero3f, up = y3f;
        if(json.object_contains("from")) json_set_value(json.object_element("from"), from);
        if(json.object_contains("to")) json_set_value(json.object_element("to"), to);
        if(json.object_contains("up")) json_set_value(json.object_element("up"), up);
        value = lookat_frame(from, to, up);
    } else {
        if(json.object_contains("o")) json_set_value(json.object_element("o"), value.o);
        if(json.object_contains("x")) json_set_value(json.object_element("x"), value.x);
        if(json.object_contains("y")) json_set_value(json.object_element("y"), value.y);
        if(json.object_contains("z")) json_set_value(json.object_element("z"), value.z);
        //value = orthonormalize_zyx(value); //TODO DEBUG
    }
}
void json_set_value(const jsonvalue& json, vector<vector<mat4f>>& value) {
    value.resize(json.array_size());
    for(auto i : range(value.size())) json_set_value(json.array_element(i), value[i]);
}

template<typename T>
void json_set_optvalue(const jsonvalue& json, T& value, const string& name) {
    if(not json.object_contains(name)) return;
    json_set_value(json.object_element(name), value);
}

Camera* json_parse_camera(const jsonvalue& json) {
    auto camera = new Camera();
    json_set_optvalue(json, camera->frame, "frame");
    json_set_optvalue(json, camera->width, "width");
    json_set_optvalue(json, camera->height, "height");
    json_set_optvalue(json, camera->focus, "focus");
    json_set_optvalue(json, camera->dist, "dist");
    return camera;
}

Camera* json_parse_lookatcamera(const jsonvalue& json) {
    auto from = z3f, to = zero3f, up = y3f;
    auto width = 1.0f, height = 1.0f, dist = 1.0f;
    json_set_optvalue(json, from, "from");
    json_set_optvalue(json, to, "to");
    json_set_optvalue(json, up, "up");
    json_set_optvalue(json, width, "width");
    json_set_optvalue(json, height, "height");
    json_set_optvalue(json, dist, "dist");
    vec3f point = zero3f;
    json_set_optvalue(json, point, "focus_to");
    auto cam = lookat_camera(from, to, up, width, height, dist);
    json_set_optvalue(json, cam->aperture, "aperture");
    if(!(point == zero3f)) cam->focus = length(from-point);
    json_set_optvalue(json, cam->focus, "focus");
    json_set_optvalue(json, cam->circle, "circle");
    return cam;
}

vector<string>          json_texture_paths;
map<string,image3f*>    json_texture_cache;

void json_texture_path_push(string filename) {
    auto pos = filename.rfind("/");
    auto dirname = (pos == string::npos) ? string() : filename.substr(0,pos+1);
    json_texture_paths.push_back(dirname);
}
void json_texture_path_pop() { json_texture_paths.pop_back(); }

void json_parse_opttexture(jsonvalue json, image3f*& txt, string name) {
    if(not json.object_contains(name)) return;
    auto filename = json.object_element(name).as_string();
    if(filename.empty()) { txt = nullptr; return; }
    auto dirname = json_texture_paths.back();
    auto fullname = dirname + filename;
    if (json_texture_cache.find(fullname) == json_texture_cache.end()) {
        auto ext = fullname.substr(fullname.size()-3);
        if(ext == "pfm") {
            auto image = read_pnm("models/pisa_latlong.pfm", true);
            image = image.gamma(1/2.2);
            json_texture_cache[fullname] = new image3f(image);
        } else if(ext == "png") {
            auto image = read_png(fullname,true);
            json_texture_cache[fullname] = new image3f(image);
        } else error("unsupported image format %s\n", ext.c_str());
    }
    txt = json_texture_cache[fullname];
}

Material* json_parse_material(const jsonvalue& json) {
    auto material = new Material();
    json_set_optvalue(json, material->ke, "ke");
    json_set_optvalue(json, material->kd, "kd");
    json_set_optvalue(json, material->ks, "ks");
    json_set_optvalue(json, material->kr, "kr");
    json_set_optvalue(json, material->n, "n");
    json_set_optvalue(json, material->blur, "blur");
    json_set_optvalue(json, material->blur_size, "blur_size");
    json_set_optvalue(json, material->microfacet, "microfacet");
    json_parse_opttexture(json, material->ke_txt, "ke_txt");
    json_parse_opttexture(json, material->kd_txt, "kd_txt");
    json_parse_opttexture(json, material->ks_txt, "ks_txt");
    json_parse_opttexture(json, material->norm_txt, "norm_txt");
    json_set_optvalue(json, material->hair_count, "hair_count");
    json_set_optvalue(json, material->hair_length, "hair_length");
    json_set_optvalue(json, material->hair_radius, "hair_radius");
    json_set_optvalue(json, material->translucency, "translucency");
    json_set_optvalue(json, material->t_factor, "t_factor");
    json_set_optvalue(json, material->ruffle, "ruffle");
    return material;
}

FrameAnimation* json_parse_frame_animation(const jsonvalue& json) {
    auto animation = new FrameAnimation();
    json_set_optvalue(json, animation->rest_frame, "rest_frame");
    json_set_optvalue(json, animation->keytimes, "keytimes");
    json_set_optvalue(json, animation->translation, "translation");
    json_set_optvalue(json, animation->rotation, "rotation");
    return animation;
}

Surface* json_parse_surface(const jsonvalue& json) {
    auto surface = new Surface();
    json_set_optvalue(json, surface->frame, "frame");
    json_set_optvalue(json, surface->radius,"radius");
    json_set_optvalue(json, surface->isquad,"isquad");
    json_set_optvalue(json, surface->height,"height");
    if(json.object_contains("material")) surface->mat = json_parse_material(json.object_element("material"));
    if(json.object_contains("animation")) surface->animation = json_parse_frame_animation(json.object_element("animation"));
    return surface;
}

vector<Surface*> json_parse_surfaces(const jsonvalue& json) {
    auto surfaces = vector<Surface*>();
    for(auto value : json.as_array_ref())
        surfaces.push_back( json_parse_surface(value) );
    return surfaces;
}

MeshSkinning* json_parse_mesh_skinning(const jsonvalue& json) {
    auto skinning = new MeshSkinning();
    json_set_optvalue(json, skinning->rest_pos, "rest_pos");
    json_set_optvalue(json, skinning->rest_norm, "rest_norm");
    json_set_optvalue(json, skinning->bone_ids, "bone_ids");
    json_set_optvalue(json, skinning->bone_weights, "bone_weights");
    json_set_optvalue(json, skinning->bone_xforms, "bone_xforms");
    return skinning;
}

MeshSimulation* json_parse_mesh_simulation(const jsonvalue& json) {
    auto simulation = new MeshSimulation();
    json_set_optvalue(json, simulation->init_pos, "init_pos");
    json_set_optvalue(json, simulation->init_vel, "init_vel");
    json_set_optvalue(json, simulation->mass, "mass");
    json_set_optvalue(json, simulation->pinned, "pinned");
    json_set_optvalue(json, simulation->vel, "vel");
    json_set_optvalue(json, simulation->force, "force");
    if (json.object_contains("springs")) {
        for(auto& elem : json.object_element("springs").as_array_ref()) {
            auto spring = MeshSimulation::Spring();
            json_set_optvalue(elem, spring.ids, "ids");
            json_set_optvalue(elem, spring.restlength, "restlength");
            json_set_optvalue(elem, spring.ks, "ks");
            json_set_optvalue(elem, spring.kd, "kd");
            simulation->springs.push_back(spring);
        }
    }
    return simulation;
}

frame3f get_frame(yobj::float4x4 f)
{
    return frame3f(vec3f(f[3][0],f[3][1],f[3][2]), vec3f(f[0][0],f[0][1],f[0][2]), vec3f(f[1][0],f[1][1],f[1][2]), vec3f(f[2][0],f[2][1],f[2][2]));
}

vector<Mesh*> obj_load_meshes(yobj::fl_scene& scene)
{
    vector<Mesh*> res;
    for(int m = 0; m < scene.shapes.size(); m++)
    {
        Mesh* mesh = new Mesh();
        for(int i = 0; i < scene.shapes[m].pos.size(); i++)
        {
            mesh->pos.push_back(vec3f(scene.shapes[m].pos[i][0], scene.shapes[m].pos[i][1],scene.shapes[m].pos[i][2]));
            //mesh->norm.push_back(vec3f(scene.shapes[0].norm[i][0], scene.shapes[0].norm[i][1],scene.shapes[0].norm[i][2]));
        }
        for(int i = 0; i< scene.shapes[m].triangles.size(); i++)
        {
            mesh->triangle.push_back(vec3i(scene.shapes[m].triangles[i][0], scene.shapes[m].triangles[i][1],scene.shapes[m].triangles[i][2]));
        }
        for(int i = 0; i < scene.shapes[m].texcoord.size(); i++)
        {
            mesh->texcoord.push_back(vec2f(scene.shapes[m].texcoord[i][0], scene.shapes[m].texcoord[i][1]));
        }
        if(!scene.materials.empty())
        {
            if(scene.materials[0].kd_txt>0)
            {
                auto image = new image3f(read_png(scene.textures[scene.materials[0].kd_txt].path,true));
                mesh->mat->kd_txt = image;
            }
            if(scene.materials[0].ks_txt>0)
            {
                auto image = new image3f(read_png(scene.textures[scene.materials[0].ks_txt].path,true));
                mesh->mat->ks_txt = image;
            }
            if(scene.materials[0].rs_txt>0)
            {
                auto image = new image3f(read_png(scene.textures[scene.materials[0].rs_txt].path,true));
                mesh->mat->kr_txt = image;
            }
            mesh->mat->ks = vec3f(scene.materials[m].ks[0],scene.materials[m].ks[1],scene.materials[m].ks[2]);
            mesh->mat->kd = vec3f(scene.materials[m].kd[0],scene.materials[m].kd[1],scene.materials[m].kd[2]);
            mesh->mat->ke = vec3f(scene.materials[m].ke[0],scene.materials[m].ke[1],scene.materials[m].ke[2]);
            //mesh->mat->kr = vec3f(scene.materials[0].[0],scene.materials[0].kr[1],scene.materials[0].kr[2]);
            mesh->mat->n = scene.materials[0].rs;
        }
        smooth_normals(mesh);
        res.push_back(mesh);
    }
    return res;
}


Scene* load_obj_scene(const string& filename)
{
    Scene* scene = new Scene();
    yobj::obj obj;
    string err;
    auto res = yobj::load_obj(filename, obj, err);
    if(!res){
        error("%s\n", &err);
        exit(1);
    }
    auto s = yobj::flatten_obj(obj);
    if(s.cameras.size() > 0)
    {
        //scene->camera->aperture = s.cameras[0].aperture; TODO
        scene->camera->focus = s.cameras[0].focus;
        scene->camera->frame = get_frame(s.cameras[0].xform);
    }
    scene->meshes = obj_load_meshes(s);
    return scene;
}

Mesh* json_load_obj(const string& filename)
{
    Mesh* mesh = new Mesh();
    yobj::obj obj;
    string err;
    auto res = yobj::load_obj(filename, obj, err);
    if(!res){
        error("%s\n", &err);
        exit(1);
    }
    auto scene = yobj::flatten_obj(obj);
    if(res)
    {
        int last = 0;
        for(int m = 0; m < scene.shapes.size(); m++)
        {
            for(int i = 0; i < scene.shapes[m].pos.size(); i++)
            {
                mesh->pos.push_back(vec3f(scene.shapes[m].pos[i][0], scene.shapes[m].pos[i][1],scene.shapes[m].pos[i][2]));
                //mesh->norm.push_back(vec3f(scene.shapes[0].norm[i][0], scene.shapes[0].norm[i][1],scene.shapes[0].norm[i][2]));
            }
            for(int i = 0; i< scene.shapes[m].triangles.size(); i++)
            {
                mesh->triangle.push_back(vec3i(last+scene.shapes[m].triangles[i][0], last+scene.shapes[m].triangles[i][1],last+scene.shapes[m].triangles[i][2]));
            }
            for(int i = 0; i < scene.shapes[m].texcoord.size(); i++)
            {
                mesh->texcoord.push_back(vec2f(scene.shapes[m].texcoord[i][0], scene.shapes[m].texcoord[i][1]));
            }
            if(!scene.materials.empty())
            {
                if(scene.materials[0].kd_txt>0)
                {
                    auto image = new image3f(read_png(scene.textures[scene.materials[0].kd_txt].path,true));
                    mesh->mat->kd_txt = image;
                }
                if(scene.materials[0].ks_txt>0)
                {
                    auto image = new image3f(read_png(scene.textures[scene.materials[0].ks_txt].path,true));
                    mesh->mat->ks_txt = image;
                }
                if(scene.materials[0].rs_txt>0)
                {
                    auto image = new image3f(read_png(scene.textures[scene.materials[0].rs_txt].path,true));
                    mesh->mat->kr_txt = image;
                }
                mesh->mat->ks = vec3f(scene.materials[m].ks[0],scene.materials[m].ks[1],scene.materials[m].ks[2]);
                mesh->mat->kd = vec3f(scene.materials[m].kd[0],scene.materials[m].kd[1],scene.materials[m].kd[2]);
                mesh->mat->ke = vec3f(scene.materials[m].ke[0],scene.materials[m].ke[1],scene.materials[m].ke[2]);
                //mesh->mat->kr = vec3f(scene.materials[0].[0],scene.materials[0].kr[1],scene.materials[0].kr[2]);
                mesh->mat->n = scene.materials[0].rs;
            }
            last=mesh->pos.size();
        }
    }
    smooth_normals(mesh);
    return mesh;
}

Mesh* json_parse_mesh(const jsonvalue& json) {
    auto mesh = new Mesh();
    if(json.object_contains("json_mesh")) {
        json_texture_path_push(json.object_element("json_mesh").as_string());
        mesh = json_parse_mesh(load_json(json.object_element("json_mesh").as_string()));
        json_texture_path_pop();
    }
    if(json.object_contains("obj_mesh"))
    {
        mesh = json_load_obj(json.object_element("obj_mesh").as_string());
    }
    vec3f rot = zero3f;
    bool apply = false;
    json_set_optvalue(json, apply, "apply");
    json_set_optvalue(json, rot, "rotation");
    json_set_optvalue(json, mesh->frame, "frame");
    if(apply){
        auto how = length(mesh->frame.x) < 1.f;
        vec3f pos = mesh->frame.o;
        if(!how)mesh->frame.o = zero3f;
        auto rot_frame = identity_frame3f;
        if(!(rot == zero3f))
        {
            rot_frame = transform_frame(rotation_matrix(rot.x, x3f), rot_frame, false);
            rot_frame = transform_frame(rotation_matrix(rot.y, y3f), rot_frame, false);
            rot_frame = transform_frame(rotation_matrix(rot.z, z3f), rot_frame, false);
        }
        //mesh->frame.o = pos;
        for(int i : range(mesh->pos.size()))
        {
            mesh->pos[i] = transform_point(rot_frame, mesh->pos[i]);
            mesh->pos[i] = transform_point(mesh->frame, mesh->pos[i]);
            mesh->norm[i] = transform_normal(rot_frame, mesh->norm[i]);
            mesh->norm[i] = transform_normal(mesh->frame, mesh->norm[i]);
        }
        mesh->frame = identity_frame3f;
        if(!how)
            mesh->frame.o = pos;
    }
    else{
        mesh->frame = orthonormalize_zyx(mesh->frame);
    }
    json_set_optvalue(json, mesh->invert, "invert");
    json_set_optvalue(json, mesh->pos, "pos");
    json_set_optvalue(json, mesh->norm, "norm");
    json_set_optvalue(json, mesh->texcoord, "texcoord");
    json_set_optvalue(json, mesh->triangle, "triangle");
    json_set_optvalue(json, mesh->quad, "quad");
    json_set_optvalue(json, mesh->point, "point");
    json_set_optvalue(json, mesh->line, "line");
    json_set_optvalue(json, mesh->spline, "spline");
    if(json.object_contains("material")) mesh->mat = json_parse_material(json.object_element("material"));
    if(json.object_contains("hair_material")) mesh->hair_mat = json_parse_material(json.object_element("hair_material"));
    json_set_optvalue(json, mesh->subdivision_catmullclark_level, "subdivision_catmullclark_level");
    json_set_optvalue(json, mesh->subdivision_catmullclark_smooth, "subdivision_catmullclark_smooth");
    json_set_optvalue(json, mesh->subdivision_bezier_level, "subdivision_bezier_level");
    if(json.object_contains("animation")) mesh->animation = json_parse_frame_animation(json.object_element("animation"));
    if(json.object_contains("skinning")) mesh->skinning = json_parse_mesh_skinning(json.object_element("skinning"));
    if(json.object_contains("json_skinning")) mesh->skinning = json_parse_mesh_skinning(load_json(json.object_element("json_skinning").as_string()));
    if(json.object_contains("simulation")) mesh->simulation = json_parse_mesh_simulation(json.object_element("simulation"));
    if (mesh->skinning) {
        if (mesh->skinning->rest_pos.empty()) mesh->skinning->rest_pos = mesh->pos;
        if (mesh->skinning->rest_norm.empty()) mesh->skinning->rest_norm = mesh->norm;
        if (mesh->pos.empty()) mesh->pos = mesh->skinning->rest_pos;
        if (mesh->norm.empty()) mesh->norm = mesh->skinning->rest_norm;
    }
    if(mesh->norm.size() < mesh->pos.size())
        smooth_normals(mesh);
    return mesh;
}

vector<Mesh*> json_parse_meshes(const jsonvalue& json) {
    auto meshes = vector<Mesh*>();
    for(auto& value : json.as_array_ref())
        meshes.push_back( json_parse_mesh(value) );
    return meshes;
}

Light* json_parse_light(const jsonvalue& json) {
    auto light = new Light();
    json_set_optvalue(json, light->frame, "frame");
    json_set_optvalue(json, light->intensity, "intensity");
    return light;
}

vector<Light*> json_parse_lights(const jsonvalue& json) {
    auto lights = vector<Light*>();
    for(auto& value : json.as_array_ref())
        lights.push_back( json_parse_light(value) );
    return lights;
}

SceneAnimation* json_parse_scene_animation(const jsonvalue& json) {
    auto animation = new SceneAnimation();
    json_set_optvalue(json, animation->time, "time");
    json_set_optvalue(json, animation->length, "length");
    json_set_optvalue(json, animation->dt, "dt");
    json_set_optvalue(json, animation->simsteps, "simsteps");
    json_set_optvalue(json, animation->gravity, "gravity");
    json_set_optvalue(json, animation->bounce_dump, "bounce_dump");
    json_set_optvalue(json, animation->loop, "loop");
    return animation;
}

Scene* json_parse_scene(const jsonvalue& json) {
    // prepare scene
    auto scene = new Scene();
    // camera
    if (json.object_contains("camera")) scene->camera = json_parse_camera(json.object_element("camera"));
    if (json.object_contains("lookat_camera")) scene->camera = json_parse_lookatcamera(json.object_element("lookat_camera"));
    // surfaces
    if(json.object_contains("surfaces")) scene->surfaces = json_parse_surfaces(json.object_element("surfaces"));
    // meshes
    if(json.object_contains("json_meshes")) {
        json_texture_path_push(json.object_element("json_meshes").as_string());
        scene->meshes = json_parse_meshes(load_json(json.object_element("json_meshes").as_string()));
        json_texture_path_pop();
    }
    if(json.object_contains("meshes")) {
        for(Mesh* m : json_parse_meshes(json.object_element("meshes")))
            scene->meshes.push_back(m);
    }
    // lights
    if(json.object_contains("lights")) scene->lights = json_parse_lights(json.object_element("lights"));
    // animation
    if(json.object_contains("animation")) scene->animation = json_parse_scene_animation(json.object_element("animation"));
    // rendering parameters
    json_set_optvalue(json, scene->image_width, "image_width");
    json_set_optvalue(json, scene->image_height, "image_height");
    json_set_optvalue(json, scene->image_samples, "image_samples");
    json_set_optvalue(json, scene->background, "background");
    json_parse_opttexture(json, scene->background_txt, "background_txt");
    json_set_optvalue(json, scene->ambient, "ambient");
    json_set_optvalue(json, scene->path_max_depth, "path_max_depth");
    json_set_optvalue(json, scene->path_sample_brdf, "path_sample_brdf");
    json_set_optvalue(json, scene->path_shadows, "path_shadows");
    // done
    return scene;
}

Scene* load_json_scene(const string& filename) {
    json_texture_cache.clear();
    json_texture_paths = { "" };
    auto scene = json_parse_scene(load_json(filename));
    json_texture_cache.clear();
    json_texture_paths = { "" };
    return scene;
}

