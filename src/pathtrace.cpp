#include "scene.h"
#include "intersect.h"
#include "montecarlo.h"
#include "tesselation.h"
#include "animation.h"
#include "yocto_obj.h"

#include <thread>
using std::thread;

static const vec3f half3f = one3f/2.f;

static float russian_th = 0.1f;

vec3f get_random_point_triangle(Mesh* mesh, vec3i tri, std::minstd_rand& generator);
vec3f get_normal(Mesh* mesh, vec3i tri, vec3f point);
Surface* get_cylinder(vec3f a, vec3f b, float radius);

// lookup texture value
vec3f lookup_scaled_texture(vec3f value, image3f* texture, vec2f uv, bool tile = true) {
    // YOUR CODE GOES HERE ----------------------
    if(not texture)
        return value;
    int width = texture->width();
    int height = texture->height();
    int i = (int)(uv.x*width);
    int j = (int)(uv.y*height);
    int i1 = i+1;
    int j1 = j+1;
    auto s = uv.x*width-i;
    auto t = uv.y*height-j;
    if(tile)
    {
        i = i%width;
        if(i<0) i+=width;
        i1 = i1%width;
        if(i1<0) i1+=width;
        j = j%height;
        if(j<0) j+=height;
        j1 = j1%height;
        if(j1<0) j1+=height;
    }
    else
    {
        i = clamp(i, 0, width-1); //-2?
        i1 = clamp(i1,0,width-1);
        j = clamp(j, 0, height-1);
        j1 = clamp(j1,0,height-1);
    }
//    auto s = uv.x*width-i;
//    auto t = uv.y*height-j;
    return value*(texture->at(i,j)*(1-s)*(1-t)+texture->at(i,j1)*(1-s)*t+texture->at(i1,j)*s*(1-t)+texture->at(i1,j1)*s*t);
}

// compute the brdf
vec3f eval_brdf(vec3f kd, vec3f ks, float n, vec3f v, vec3f l, vec3f norm, vec3f tan, bool microfacet, bool hair) {
    // YOUR CODE GOES HERE ----------------------
    auto h = normalize(v+l);
    if(microfacet)
    {
        auto d = (n+2)/(2*pif)*pow(max(0.f,dot(h,norm)),n);
        auto f = ks+(one3f-ks)*pow((1-dot(h,l)),5); //? one3f zero3f WARNING
        auto g = min(1.f, min((2*dot(h,norm)*dot(v,norm))/dot(v,h),(2*dot(norm,h)*dot(l,norm))/dot(l,h)));
        return (d*g*f)/(4.f*dot(l,norm)*dot(v,norm));
    }
//    else if (hair)
//    {
//        auto res = zero3f;
//        if(kd == zero3f && ks == zero3f)
//            return res;
//        auto ndo = dot(norm, v), ndi = dot(norm,l), ndh = dot(h,norm);
//        auto so = sqrt(1-ndo*ndo), si = sqrt(1 - ndi * ndi), sh = sqrt(1 - ndh * ndh);
//        // diffuse term (Kajiya-Kay)
//        if (si > 0 && so > 0 && !(kd == zero3f)) {
//            auto diff = kd * si / pif;
//            res += diff;
//        }
//        
//        // specular term (Kajiya-Kay)
//        if (si > 0 && so > 0 && sh > 0 && !(ks == zero3f)) {
//            auto ns = 2 / (n * n) - 2;
//            auto d = (ns + 2) * pow(sh, ns) / (2 + pif);
//            auto spec = ks * si * d / (4 * si * so);
//            res += spec;
//        }
//        return res;
//    }
     // placeholder (non-microfacet model)
    return kd/pif + ks*(n+8)/(8*pif) * pow(max(0.0f,hair?sqrt(1- dot(tan, h)*dot(tan,h)):dot(norm,h)),n);// * max(0.f,sqrt(1- dot(norm, h)*dot(norm,h))); // placeholder (non-microfacet model)
}

// evaluate the environment map
vec3f eval_env(vec3f ke, image3f* ke_txt, vec3f dir) {
    // YOUR CODE GOES HERE ----------------------
    if(not ke_txt) return ke;
    auto u = atan2(dir.x,dir.z)/(2*pif); //WARNING
    auto v = 1 - acos(dir.y)/pif; //WARNING
    if (isnan(u))
        u = 0.f;
    if (isnan(v))
        v = 0.f;
    return lookup_scaled_texture(ke, ke_txt, vec2f(u,v), true);
}

// pick a direction according to the cosine (returns direction and its pdf)
pair<vec3f,float> sample_cosine(vec3f norm, vec2f ruv) {
    auto frame = frame_from_z(norm);
    auto l_local = sample_direction_hemispherical_cosine(ruv);
    auto pdf = sample_direction_hemispherical_cosine_pdf(l_local);
    auto l = transform_direction(frame, l_local);
    return {l,pdf};
}

// pick a direction according to the brdf (returns direction and its pdf)
pair<vec3f,float> sample_brdf(vec3f kd, vec3f ks, float n, vec3f v, vec3f norm, vec2f ruv, float rl) {
    if(ks == zero3f) return sample_cosine(norm, ruv);
    auto frame = frame_from_z(norm);
    auto dw = mean(kd) / (mean(kd) + mean(ks));
    auto v_local = transform_direction_inverse(frame, v);
    auto l_local = zero3f, h_local = zero3f;
    if(rl < dw) {
        l_local = sample_direction_hemispherical_cosine(ruv);
        h_local = normalize(l_local+v_local);
    } else {
        h_local = sample_direction_hemispherical_cospower(ruv, n);
        l_local = -v_local + h_local*2*dot(v_local,h_local);
    }
    auto l = transform_direction(frame, l_local);
    auto dpdf = sample_direction_hemispherical_cosine_pdf(l_local);
    auto spdf = sample_direction_hemispherical_cospower_pdf(h_local,n) / (4*dot(v_local,h_local));
    auto pdf = dw * dpdf + (1-dw) * spdf;
    return {l,pdf};
}

// compute the color corresponing to a ray by pathtrace
vec3f pathtrace_ray(Scene* scene, ray3f ray, Rng* rng, int depth, bool reflected = false) {
    
    //russian roulette
    if(rng->next_float() < russian_th) return zero3f;
    
    // get scene intersection
    auto intersection = intersect(scene,ray);
    
    // if not hit, return background (looking up the texture by converting the ray direction to latlong around y)
    if(not intersection.hit) {
        // YOUR CODE GOES HERE ----------------------
        //environment
        return eval_env(scene->background, scene->background_txt, ray.d);
    }
    
    // setup variables for shorter code
    auto pos = intersection.pos;
    auto norm = intersection.norm;
    auto tan = intersection.tan;
    auto v = -ray.d;
    auto hair = intersection.mat->hair_count > 0;
    
    // compute material values by looking up textures
    // YOUR CODE GOES HERE ----------------------
    vec2f uv = intersection.texcoord;
    auto ke = lookup_scaled_texture(intersection.mat->ke, intersection.mat->ke_txt, uv);
    auto kd = lookup_scaled_texture(intersection.mat->kd, intersection.mat->kd_txt, uv);
    auto ks = lookup_scaled_texture(intersection.mat->ks, intersection.mat->ks_txt, uv);
    norm = lookup_scaled_texture(norm, intersection.mat->norm_txt, uv);
    auto kr = intersection.mat->kr;
    auto n = intersection.mat->n;
    auto mf = intersection.mat->microfacet;
    
    bool radiance = depth < 0;
    
    // accumulate color starting with ambient
    auto c = scene->ambient * kd;
    
    // add emission if on the first bounce
    // YOUR CODE GOES HERE ----------------------
    if(depth <= 0 || reflected)
        c+=ke;
    
    if(radiance) depth = 0;
    
    // foreach point light
    for(auto light : scene->lights) {
        if(radiance) continue;
        // compute light response
        auto cl = light->intensity / (lengthSqr(light->frame.o - pos));
        // compute light direction
        auto l = normalize(light->frame.o - pos);
        // compute the material response (brdf*cos)
        auto brdfcos = max(dot(norm,l),0.0f) * eval_brdf(kd, ks, n, v, l, norm, tan, mf, hair);
        // multiply brdf and light
        auto shade = cl * brdfcos;
        // check for shadows and accumulate if needed
        if(shade == zero3f) continue;
        // if shadows are enabled
        if(scene->path_shadows) {
            // perform a shadow check and accumulate
            if(not intersect_shadow(scene,ray3f::make_segment(pos,light->frame.o))) c += shade;
        } else {
            // else just accumulate
            c += shade;
        }
    }
    
    // YOUR AREA LIGHT CODE GOES HERE ----------------------
    // foreach surface
    for(auto surf : scene->surfaces)
    {
        // skip if no emission from surface
        if(surf->mat->ke == zero3f || radiance)
            continue;
        // pick a point on the surface, grabbing normal, area and texcoord
        vec3f lpos = zero3f;
        vec3f lnorm = zero3f;
        float area = 0;
        vec2f texcoord = zero2f;
        // check if quad
        if(surf->isquad)
        {
            // generate a 2d random number
            vec2f rnd = rng->next_vec2f();
            // compute light position, normal, area
            lpos = transform_point(surf->frame, surf->radius*2.f * vec3f(rnd.x-0.5f,rnd.y-0.5f,0.f));
            lnorm = surf->frame.z;
            area = pow(2*surf->radius,2);
            // set tex coords as random value got before
            texcoord = rnd;
        }
        // else
        else
        {
            // generate a 2d random number
            vec2f rnd = rng->next_vec2f();
            // compute light position, normal, area
            auto theta = 2*pif*rnd.x;
            auto phi = acos(2 * rnd.y - 1); //WARNING
            if(isnan(phi))
                phi = 0;
            auto radius = surf->radius;
            auto point = vec3f(radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), radius*cos(phi));
            lpos = surf->frame.o+point;
            lnorm = normalize(lpos-surf->frame.o);
            area = 4*pif*pow(radius,2);
            // set tex coords as random value got before
            texcoord = rnd;
        }

        // get light emission from material and texture
        auto ke = lookup_scaled_texture(surf->mat->ke, surf->mat->ke_txt, texcoord);
        // compute light direction
        auto ldir = normalize(lpos-pos);
        // compute light response
        auto cl = ke *area*max(0.f,-dot(lnorm, ldir))/ (lengthSqr(lpos - pos));
        // compute the material response (brdf*cos)
        auto brdfcos = max(dot(norm,ldir),0.0f) * eval_brdf(kd, ks, n, v, ldir, norm, tan, mf, hair);
        // multiply brdf and light
        auto shade = cl * brdfcos;
        // check for shadows and accumulate if needed
        if(shade == zero3f) continue;
        // if shadows are enabled
        if(scene->path_shadows) {
            // perform a shadow check and accumulate
            if(not intersect_shadow(scene,ray3f::make_segment(pos,lpos)))
            {
                c += shade;
            }
        } else {
            // else just accumulate
            c += shade;
        }
    }
    
    // foreach mesh
    for(auto mesh : scene->meshes)
    {
        // skip if no emission from surface
        if(mesh->mat->ke == zero3f || radiance)
            continue;
        // pick a point on the surface, grabbing normal, area and texcoord
        vec3f lpos = zero3f;
        vec3f lnorm = zero3f;
        float area = mesh->area;
        vec2f texcoord = zero2f; //TODO
        
        vec3i t = mesh->triangle[mesh->distribution(rng->engine)];
        lpos = get_random_point_triangle(mesh, t, rng->engine);
        lnorm = get_normal(mesh, t, lpos);
        
        // get light emission from material and texture
        auto ke = lookup_scaled_texture(mesh->mat->ke, mesh->mat->ke_txt, texcoord);
        // compute light direction
        auto ldir = normalize(lpos-pos);
        // compute light response
        auto cl = ke *area*max(0.f,-dot(lnorm, ldir))/ (lengthSqr(lpos - pos));
        // compute the material response (brdf*cos)
        auto brdfcos = max(dot(norm,ldir),0.0f) * eval_brdf(kd, ks, n, v, ldir, tan, norm, mf, hair);
        // multiply brdf and light
        auto shade = cl * brdfcos;
        // check for shadows and accumulate if needed
        if(shade == zero3f) continue;
        // if shadows are enabled
        if(scene->path_shadows) {
            // perform a shadow check and accumulate
            if(not intersect_shadow(scene,ray3f::make_segment(pos,lpos)))
            {
                c += shade;
            }
        } else {
            // else just accumulate
            c += shade;
        }
    }
    
    // YOUR ENVIRONMENT LIGHT CODE GOES HERE ----------------------
    // sample the brdf for environment illumination if the environment is there
    if(scene->background_txt)
    {
        // pick direction and pdf
        auto ruv = rng->next_vec2f();
//        auto dir = sample_direction_hemispherical_uniform(uv);
//        auto pdf = sample_direction_hemispherical_uniform_pdf(dir);
        auto pp = sample_brdf(kd, ks, n, v, norm, ruv, rng->next_float()); //rng.nextfloat();
        // compute the material response (brdf*cos)
        auto brdfcos = max(dot(norm,pp.first),0.0f) * eval_brdf(kd, ks, n, v, pp.first, norm, tan, mf, hair);
        // accumulate recersively scaled by brdf*cos/pdf
        auto cl = eval_env(scene->background, scene->background_txt, pp.first)/pp.second;
        auto shade = cl * brdfcos;// * pathtrace_ray(scene, ray3f(pos, pp.first), rng, depth+1);
        // if shadows are enabled
        if(scene->path_shadows) {
            // perform a shadow check and accumulate
            if(not intersect_shadow(scene,ray3f(pos,pp.first)))
            {
                c += shade;
            }
        } else {
            // else just accumulate
            c += shade;
        }

    }

    // YOUR INDIRECT ILLUMINATION CODE GOES HERE ----------------------
    // sample the brdf for indirect illumination
    if(depth < scene->path_max_depth)
    {
        // pick direction and pdf
        auto ruv = rng->next_vec2f();
        auto pp = sample_brdf(kd, ks, n, v, norm, ruv, rng->next_float()); //rng.nextfloat();
        // compute the material response (brdf*cos)
        auto brdfcos = max(dot(norm,pp.first),0.0f) * eval_brdf(kd, ks, n, v, pp.first, norm, tan, mf, hair);
        // accumulate recersively scaled by brdf*cos/pdf
        c += brdfcos * pathtrace_ray(scene, ray3f(pos, pp.first), rng, depth+1, radiance)/pp.second;
    }
    
    //REFLECTION
    if(!(kr == zero3f) && depth < scene->path_max_depth)
    {
        auto dir = normalize(ray.d - 2*(dot(intersection.norm,ray.d))*intersection.norm); //mirror
        if(intersection.mat->blur > 0)
        {
            auto blurred = zero3f;
            for(int s : range(intersection.mat->blur))
            {
                dir = dir+(rng->next_vec3f()-half3f)*intersection.mat->blur_size; //cono
                ray3f ray = ray3f(pos, dir);
                blurred+=kr*pathtrace_ray(scene, ray, rng, depth+1, true);
            }
            c += blurred / intersection.mat->blur;
        }
        else
        {
            c+=kr*pathtrace_ray(scene, ray3f(pos, dir), rng, depth+1,true);
        }
    }
    
    //SINGLE SCATTERING
//    if(false)
//    {
//        vec3f To = normalize(refract(v,norm,1/1.3f));
//        
//        for(int i : range(scene->image_samples))
//        {
//            auto pdf_pair = sample_cosine(, )
//            fresnel
//        }
//    }
    
    //TRANSLUCENCY
    if(!(intersection.mat->translucency == zero3f) && !radiance)
    {
        c+= kd*intersection.mat->translucency*intersection.in_translucency;/// (lengthSqr(zero3f - pos));
    }
    
    // return the accumulated color
    return c;
}

static const float aperture_time = 3.f;

// pathtrace an image
void pathtrace(Scene* scene, image3f* image, RngImage* rngs, int offset_row, int skip_row, bool verbose) {
    if(verbose) message("\n  rendering started        ");
    auto aperture = scene->camera->aperture;
    if(verbose) message("FOCUS: %f\n", scene->camera->focus);
    auto time = scene->animation->time;
    // foreach pixel
    for(auto j = offset_row; j < scene->image_height; j += skip_row ) {
        if(verbose) message("\r  rendering %03d/%03d        ", j, scene->image_height);
        for(auto i = 0; i < scene->image_width; i ++) {
            // init accumulated color
            image->at(i,j) = zero3f;
            // grab proper random number generator
            auto rng = &rngs->at(i, j);
            // foreach sample
            for(auto jj : range(scene->image_samples)) {
                for(auto ii : range(scene->image_samples)) {
                    //DOF lens point
                    // compute ray-camera parameters (u,v) for the pixel and the sample
                    auto u = (i+(ii+0.5f)/scene->image_samples)/scene->image_width;//                    auto u = (i + (ii + rng->next_float())/scene->image_samples) /scene->image_width;
                    auto v = (j+(jj+0.5f)/scene->image_samples)/scene->image_height;//                    auto v = (j + (jj + rng->next_float())/scene->image_samples) /scene->image_height;
                            if(time >= 0) //motion blur
                            {
                                auto ap = rng->next_vec2f();
                                auto rnd = sqrt(ap.x)*cos(pif/2.f*ap.y);
                                auto t = float(rnd*aperture_time + time);//(rng->next_float() < 0.5f && time > 0)? time -1: time;
                                animate_frame(scene,t);
                                animate_skin(scene,t);
                            }
                            auto la = rng->next_vec2f();
                            vec3f e;
                            if(scene->camera->circle)
                            {
                                auto theta = pif2*la.x;
                                auto r = sqrt(la.y)*aperture/2.f;
                                e = vec3f(cos(theta)*r, sin(theta)*r,0.f);///scene->image_samples; //-0.5, 0.5
                            }
                            else
                            {
                                e = vec3f(aperture*(la.x-0.5f),aperture*(la.y-0.5f),0.f);///scene->image_samples;
                            }
                            auto q = vec3f((u-0.5f)*scene->camera->width,(v-0.5f)*scene->camera->height,-scene->camera->dist)*scene->camera->focus;
                            auto ray = transform_ray(scene->camera->frame, ray3f(e, normalize(q-e)));
                            auto trace = pathtrace_ray(scene,ray,rng,0);
                            image->at(i,j) += isnan(trace)? zero3f: trace;///(aperture == 0? 1 : (scene->image_samples*scene->image_samples));
                }
            }
            // scale by the number of samples
            image->at(i,j) /= pow(scene->image_samples, 2)*(1-russian_th);
        }
    }
    if(verbose) message("\r  rendering done        \n");
    
}

// pathtrace an image with multithreading if necessary
image3f pathtrace(Scene* scene, bool multithread) {
    //set roulette just for more than 1step
    if(scene[0].path_max_depth == 0)
        russian_th = 0.f;
    // allocate an image of the proper size
    auto image = image3f(scene[0].image_width, scene[0].image_height);
    
    // create a random number generator for each pixel
    auto rngs = RngImage(scene[0].image_width, scene[0].image_height);
    
    bool animate = scene[0].animation->time != -1;
   
    // if multitreaded
    if(multithread) {
        // get pointers
        auto image_ptr = &image;
        auto rngs_ptr = &rngs;
        // allocate threads and pathtrace in blocks
        auto threads = vector<thread>();
        auto nthreads = thread::hardware_concurrency();
        for(auto tid : range(nthreads)) threads.push_back(thread([=](){
            return pathtrace(scene+(animate? tid : 0),image_ptr,rngs_ptr,tid,nthreads,tid==0);}));
        for(auto& thread : threads) thread.join();
    } else {
        // pathtrace all rows
        pathtrace(scene, &image, &rngs, 0, 1, true);
    }
    
    // done
    return image;
}

float get_triangle_area(Mesh* m, int idx)
{
    auto tri = m->triangle[idx];
    vec3f a,b,c;
    a = m->pos[tri.x];
    b = m->pos[tri.y];
    c = m->pos[tri.z];
    auto ab = b-a;
    auto ac = c-a;
    return length(cross(ab,ac))/2.f;
}

vec3f get_random_point_triangle(Mesh* mesh, vec3i tri, std::minstd_rand& generator)
{
    auto A = mesh->pos[tri.x];
    auto B = mesh->pos[tri.y];
    auto C = mesh->pos[tri.z];
    //u = 1 - sqrt(r1) and v = r2 * sqrt(r1)
    auto r1 = sqrt(generator()/((float)generator.max()));
    auto r2 = generator()/((float)generator.max());
    auto u = 1 - r1;
    auto v = r2 * r1;
    auto p = A + u*(B-A) + v*(C-A);
    return transform_point(mesh->frame, p);
}

vec3f get_normal(Mesh* mesh, vec3i tri, vec3f point)
{
    point = transform_point_inverse(mesh->frame, point);
    return normalize((((mesh->norm[tri.x] + mesh->norm[tri.y] + mesh->norm[tri.z])/3.f)/((mesh->pos[tri.x] + mesh->pos[tri.y] + mesh->pos[tri.z])/3.f))*point);
}

void create_distribution(Mesh* mesh)
{
    //to triangle mesh;
    for(auto quad : mesh->quad){
        mesh->triangle.push_back({quad.x,quad.y,quad.z});
        mesh->triangle.push_back({quad.z,quad.w,quad.x});
    }
    mesh->quad.clear();
    smooth_normals(mesh);
    
    //init weights
    vector<float> weights = vector<float>();
    float totArea = 0;
    for(int i = 0; i < mesh->triangle.size(); i++){
        auto area = get_triangle_area(mesh, i);
        weights.push_back(area);
        totArea+=area;
    }
    std::discrete_distribution<> dist(weights.begin(), weights.end());
    mesh->distribution = dist;
    mesh->area = totArea;
}

void create_light_structures(Scene* scene)
{
    for(auto mesh : scene->meshes)
    {
        if(!(mesh->mat->ke == zero3f))
        {
            create_distribution(mesh);
        }
    }
}

void preprocess_hair(Scene* scene)
{
    Rng rng;
    rng.seed(std::random_device()());
    for(auto mesh : scene->meshes)
    {
        if(mesh->mat->hair_count)
        {
            create_distribution(mesh);
            Mesh *m = new Mesh();
            message("Creating %d hair...\t", mesh->mat->hair_count);
            for(int i : range(mesh->mat->hair_count))
            {
                int t = mesh->distribution(rng.engine);
                auto tri = mesh->triangle[t];
                auto rnd_pos = get_random_point_triangle(mesh, tri, rng.engine);
                auto normal = get_normal(mesh, tri, rnd_pos);

                int pos = m->pos.size();
//                auto dir = normalize(normal+(rng.next_vec3f()-half3f)*0.2f); //cono
                m->pos.push_back(rnd_pos);
                m->pos.push_back(rnd_pos + normalize(normal+(rng.next_vec3f()-half3f)*mesh->mat->ruffle)*0.2f*mesh->mat->hair_length); //0.4
                m->pos.push_back(rnd_pos + normalize(normal+(rng.next_vec3f()-half3f)*mesh->mat->ruffle)*0.4f*mesh->mat->hair_length+vec3f(0.0f,-0.4f,0.0f)*0.5*mesh->mat->hair_length);
                m->pos.push_back(rnd_pos + normalize(normal+(rng.next_vec3f()-half3f)*mesh->mat->ruffle)*0.55f*mesh->mat->hair_length+vec3f(0.0f,-0.8f,0.0f)*0.5*mesh->mat->hair_length);
                m->spline.push_back({pos, pos+1, pos+2, pos+3});
            }
            m->subdivision_bezier_level = 4;
            subdivide_bezier(m);
            
            //creo cilindretti
            for(auto line : m->line)
            {
                const float radius = mesh->mat->hair_radius*mesh->mat->hair_length;
                auto cyl = get_cylinder(m->pos[line.x], m->pos[line.y], radius);
                cyl->mat = !mesh->hair_mat? mesh->mat : mesh->hair_mat;
                cyl->mat->hair_count = 1;
                frame3f f1 = frame_from_z(cyl->tan);
                f1.o = m->pos[line.x];
                //boundig box per accelerazione
                vec3f p1 = vec3f(radius, radius, 0);
                vec3f p2 = vec3f(-radius, radius, 0);
                vec3f p3 = vec3f(radius, -radius, 0);
                vec3f p4 = vec3f(-radius, -radius, 0);
                p1 = transform_point(f1, p1);
                p2 = transform_point(f1, p2);
                p3 = transform_point(f1, p3);
                p4 = transform_point(f1, p4);
                int pos = m->pos.size();
                m->pos.push_back(p1);
                m->pos.push_back(p2);
                m->pos.push_back(p3);
                m->pos.push_back(p4);
                m->triangle.push_back({pos,pos+1,pos+2});
                m->triangle.push_back({pos+2,pos+3,pos+1});
                f1.o = m->pos[line.y];
                vec3f p5 = vec3f(radius, radius, 0);
                vec3f p6 = vec3f(-radius, radius, 0);
                vec3f p7 = vec3f(radius, -radius, 0);
                vec3f p8 = vec3f(-radius, -radius, 0);
                p5 = transform_point(f1, p5);
                p6 = transform_point(f1, p6);
                p7 = transform_point(f1, p7);
                p8 = transform_point(f1, p8);
                pos = m->pos.size();
                m->pos.push_back(p5);
                m->pos.push_back(p6);
                m->pos.push_back(p7);
                m->pos.push_back(p8);
                m->triangle.push_back({pos,pos+1,pos+2});
                m->triangle.push_back({pos+2,pos+3,pos+1});
                pos-=4;
                //top
                m->triangle.push_back({pos,pos+1,pos+5});
                m->triangle.push_back({pos+5,pos+4,pos});
                //right
                m->triangle.push_back({pos+1,pos+5,pos+7});
                m->triangle.push_back({pos+7,pos+1,pos+3});
                //left
                m->triangle.push_back({pos,pos+4,pos+6});
                m->triangle.push_back({pos+6,pos,pos+2});
                //bottom
                m->triangle.push_back({pos+2,pos+6,pos+7});
                m->triangle.push_back({pos+7,pos+2,pos+3});
                
                int s_idx = scene->surfaces_aux.size();
                scene->surfaces_aux.push_back(cyl);
                for(auto x : range(12))
                    m->bvh_surfaces.push_back(s_idx);
                
            }
            mesh->mat->hair_count = 0;
            //smooth_normals(m);
            scene->meshes.push_back(m);
            message("OK\n");
        }
    }
    //scene->meshes.clear();
}

//calcolo, sfoco e peso nel tracer/intersector
void compute_inradiance(Scene *s, bool multithread = true)
{
    const static int samples = 1000;
    for(Mesh *m : s->meshes)
    {
        if(!(m->mat->translucency == zero3f))
        {
            smooth_normals(m);
            m->in_radiance.resize(m->pos.size());
            message("\r  Computing input radiance...\t");
            auto nthreads = thread::hardware_concurrency();
            // create a random number generator for each pixel
            
            // if multitreaded
            if(multithread) {
                // get pointers
                // allocate threads and pathtrace in blocks
                auto threads = vector<thread>();
                int rr = m->pos.size()/float(nthreads);
                int progress = (m->pos.size()/nthreads/100)+1;
                for(auto tid : range(nthreads)) threads.push_back(thread([=](){
                    Rng rng;
                    rng.seed(tid+1);
                    for(int x : range(tid*rr, tid == nthreads -1 ? m->pos.size(): (tid+1)*rr))
                    {
                        if(tid == 0 && x%progress == 0)
                            message("\r  Computing input radiance %05d/%05d        ", x*nthreads, m->pos.size());
                        for(int k : range(samples))
                        {
                            auto rnd = rng.next_vec2f();
                            auto dir = sample_cosine(m->norm[x], rnd);
                            m->in_radiance[x] += pathtrace_ray(s, transform_ray_inverse(m->frame,ray3f((m->invert*m->pos[x]), dir.first)), &rng, -1)/(1/pif2); //(dir.second);//*(1-dir.second);
                        }
                        m->in_radiance[x] /=float(samples)*(1-russian_th);
                    }
                }));
                for(auto& thread : threads) thread.join();
            } else {
                // pathtrace all rows
                Rng rng;
                rng.seed(samples);
                for(int x : range(m->pos.size()))
                {
                    if(x%samples == 0)
                        message("\r  Computing input radiance %05d/%05d        ", x, m->pos.size());
                    for(int k : range(samples))
                    {
                        auto rnd = rng.next_vec2f();
                        auto dir = sample_cosine(m->norm[x], rnd);
                        m->in_radiance[x] += pathtrace_ray(s, ray3f(m->pos[x], dir.first), &rng, -1)/(1/pif2); //(dir.second);//*(1-dir.second);
                    }
                    m->in_radiance[x] /=float(samples)*(1-russian_th);
                }
            }
            message("OK\n");
            //sfocare
            message("\r  Blurring input radiance...        ");
            // if multitreaded
            if(multithread) {
                // get pointers
                // allocate threads and pathtrace in blocks
                auto threads = vector<thread>();
                vector<vec3f> copy(m->pos.size(), zero3f);
                auto c = &copy;
                int progress = (m->pos.size()/nthreads/100)+1;
                int rr = m->pos.size()/float(nthreads);
                for(auto tid : range(nthreads)) threads.push_back(thread([=, &copy](){
                    //auto copy = *c;
                    for(int i : range(tid*(rr), tid == nthreads -1 ? m->pos.size(): (tid+1)*rr))
                    {
                        if(tid == 0 && i%progress == 0)
                            message("\r  Blurring input radiance %05d/%05d        ", i*nthreads, m->pos.size());
                        //copy[i] = (m->in_radiance[i]);
                        float c = 0;
                        for(auto j : range(m->pos.size()))
                        {
                            auto w = exp(m->mat->t_factor*-distSqr(m->pos[i], m->pos[j]));
                            copy[i] += m->in_radiance[j]*w;
                            c+=w;
                        }
                        copy[i] /= float(c);//m->pos.size();
                    }
                }));
                for(auto& thread : threads) thread.join();
                 m->in_radiance = copy;
            } else {
                vector<vec3f> copy;
                for(auto i : range(m->pos.size()))
                {
                    if(i%samples == 0)
                        message("\r  Blurring input radiance %05d/%05d        ", i, m->pos.size());
                    copy.push_back(zero3f);//m->in_radiance[i]);
                    float c = 0;
                    for(auto j : range(m->pos.size()))
                    {

                        auto w = exp(m->mat->t_factor*-distSqr(m->pos[i], m->pos[j]));
                        copy[i] += m->in_radiance[j]*w;
                        c+=w;///0.25f;
                    }
                    copy[i] /= float(c);
                }
                m->in_radiance = copy;
            }
            message("OK\n");
        }
    }
}

template <typename... Ts>
std::string fmt (const std::string &fmt, Ts... vs)
{
    char b;
    unsigned required = std::snprintf(&b, 0, fmt.c_str(), vs...) + 1;
    // See comments: the +1 is necessary, while the first parameter
    //               can also be set to nullptr
    
    char bytes[required];
    std::snprintf(bytes, required, fmt.c_str(), vs...);
    
    return std::string(bytes);
}

Surface* get_cylinder(vec3f a, vec3f b, float radius)
{
    Surface *s = new Surface();
    s->height = length(b-a);
    vec3f mid = (a+b)/2.f;
    auto tan = normalize(b-a);
//    float x2 = 1;
//    float y2 = 1;
//    float z2 = (-tan.x*x2-tan.y*y2)/tan.z;
    auto n = normalize(cross(tan, one3f));//
    auto center = mid+n;
    frame3f f = lookat_frame(mid, center, tan);
    s->frame = f;
    s->radius = radius;
    s->tan = tan; //serve a creare i frame dei box
    return s;
}

// runs the raytrace over all tests and saves the corresponding images
int main(int argc, char** argv) {
    auto args = parse_cmdline(argc, argv,
        { "05_pathtrace", "raytrace a scene",
            {  {"resolution", "r", "image resolution", "int", true, jsonvalue() }  },
            {  {"scene_filename", "", "scene filename", "string", false, jsonvalue("scene.json")},
               {"image_filename", "", "image filename", "string", true, jsonvalue("")}  }
        });
    auto scene_filename = args.object_element("scene_filename").as_string();
    auto image_filename = (args.object_element("image_filename").as_string() != "") ?
        args.object_element("image_filename").as_string() :
        scene_filename.substr(0,scene_filename.size()-5)+".png";
    auto nthreads = thread::hardware_concurrency();
    auto scene1 = scene_filename[scene_filename.size()-1] == 'j' ? load_obj_scene(scene_filename) : load_json_scene(scene_filename);
    
    //DEBUG
//    yobj::obj obj;
//    yobj::fl_scene sce;
//    triangulate(scene1);
//    int k = 0;
//    yobj::fl_camera cam;
//    auto mat = frame_to_matrix(scene1->camera->frame);
//    cam.xform = {mat.x.x,mat.y.x,mat.z.x,mat.w.x,
//                 mat.x.y,mat.y.y,mat.z.y,mat.w.y,
//                 mat.x.z,mat.y.z,mat.z.z,mat.w.z,
//                 mat.x.w,mat.y.w,mat.z.w,mat.w.w,
//};
//    sce.cameras.push_back(cam);
//    for(Mesh* m : scene1->meshes)
//    {
//        smooth_normals(m);
//        yobj::fl_shape s;
//        for(auto i : range(m->pos.size()))
//        {
//            auto v = m->pos[i];
//            s.pos.push_back({v.x,v.y,v.z});
//            auto n = m->norm[i];
//            s.norm.push_back({n.x,n.y,n.z});
//            if(m->texcoord.size() > i){
//            auto t = m->texcoord[i];
//                s.texcoord.push_back({t.x,t.y});}
//        }
//        for(auto t : m->triangle)
//        {
//            s.triangles.push_back({t.x,t.y,t.z});
//        }
//        s.name = k++;
//        sce.shapes.push_back(s);
//    }
//    obj = yobj::unflatten_obj(sce);
//    yobj::save_obj("/Users/fscozzafava/Downloads/export.obj", obj, scene_filename);
//    message("%s",&scene_filename);
    
    scene1->animation->time = scene1->animation->length > 0? 0:-1;
    
    bool animate = scene1->animation->time != -1;
    
    Scene *scene = new Scene[nthreads];
    scene[0] = *scene1;
    if(animate)
    {
        for(int x : range(nthreads-1))
            scene[x+1] = *(scene_filename[scene_filename.size()-1] == 'j' ? load_obj_scene(scene_filename) : load_json_scene(scene_filename));
    }
    
    if(not args.object_element("resolution").is_null()) {
        if(animate)
        {
            for(int x : range(nthreads-1)){
                scene[x+1].image_height = args.object_element("resolution").as_int();
                scene[x+1].image_width = scene[x+1].camera->width * scene[x+1].image_height / scene[x+1].camera->height;
            }
        }
        scene->image_height = args.object_element("resolution").as_int();
        scene->image_width = scene->camera->width * scene->image_height / scene->camera->height;
    }
    
    do{
        preprocess_hair(scene);
        if(!animate){
            subdivide(scene);
            message("Creating acceleration structues...\t");
            accelerate(scene);
            message("OK\n");
            compute_inradiance(scene);
        }
        else
            for(int x : range(nthreads))
                triangulate(scene+x);
        create_light_structures(scene);
        message("rendering %s ... ", scene_filename.c_str());
        auto image = pathtrace(scene,true);
        if(animate)
        {
            image_filename = scene_filename.substr(0,scene_filename.size()-5)+"_KEY"+fmt("%05d",scene->animation->time)+".png";
            for(int x : range(nthreads))
                scene[x].animation->time++;
            if(scene->animation->time >= scene->animation->length)
                break; //fine animazione
//            animate_frame(scene, scene->animation->time); //aggiorniamo di un frame per aggiustare la struttura di accelerazione
//            animate_skin(scene, scene->animation->time);
        }
        write_png(image_filename, image, true);
    }while(scene->animation->length > 0);
    delete[] scene;
    message("done\n");
}
