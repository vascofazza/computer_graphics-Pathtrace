#include "animation.h"
#include "tesselation.h"
#include <tuple>
#include <unordered_set>
#include <unordered_map>

// compute the frame from an animation
frame3f animate_compute_frame(FrameAnimation* animation, float time) {
    // grab keyframe interval
    auto interval = 0;
    for(auto t : animation->keytimes) if(time < t) break; else interval++;
    interval--;
    // get translation and rotation matrices
    auto t = float(time-animation->keytimes[interval])/float(animation->keytimes[interval+1]-animation->keytimes[interval]);
    auto m_t = translation_matrix(animation->translation[interval]*(1-t)+animation->translation[interval+1]*t);
    auto m_rz = rotation_matrix(animation->rotation[interval].z*(1-t)+animation->rotation[interval+1].z*t,z3f);
    auto m_ry = rotation_matrix(animation->rotation[interval].y*(1-t)+animation->rotation[interval+1].y*t,y3f);
    auto m_rx = rotation_matrix(animation->rotation[interval].x*(1-t)+animation->rotation[interval+1].x*t,x3f);
    // compute combined xform matrix
    auto m = m_t * m_rz * m_ry * m_rx;
    // return the transformed frame
    return transform_frame(m, animation->rest_frame);
}

// update mesh frames for animation
void animate_frame(Scene* scene, float time) {
    // YOUR CODE GOES HERE ---------------------
    // foreach mesh
    for(auto mesh : scene->meshes)
    {
        // if not animation, continue
        if(not mesh->animation)
            continue;
        // update frame
        mesh->frame = animate_compute_frame(mesh->animation, time);
    }
    // foreach surface
    for(auto surf : scene->surfaces){
        // if not animation, continue
        if(not surf->animation)
            continue;
        // update frame
        surf->frame = animate_compute_frame(surf->animation, time);
        // update the _display_mesh
        surf->_display_mesh->frame = surf->frame;
    }
}

// skinning scene
void animate_skin(Scene* scene, int time) {
    // YOUR CODE GOES HERE ---------------------
    // foreach mesh
    for(auto mesh : scene->meshes)
    {
        // if no skinning, continue
        if(not mesh->skinning) continue;
        // foreach vertex index
        for(int i = 0; i < mesh->pos.size(); i++) //ogni vertice e' assegnato ad un bone o piu a seconda dell idx
        {
            // set pos/norm to zero
            mesh->pos[i] = zero3f;
            mesh->norm[i] = zero3f;
            // for each bone slot (0..3)
            for(int j = 0; j < 4; j++)
            {
                // get bone weight and index
                auto idx = mesh->skinning->bone_ids[i][j];
                auto w = mesh->skinning->bone_weights[i][j];
                // if index < 0, continue
                if(idx < 0) continue; //se il bone e' attivo
                // grab bone xform
                auto xform = mesh->skinning->bone_xforms[time][idx];
                // update position and normal
                mesh->pos[i] += transform_point(xform, mesh->skinning->rest_pos[i])*w;
                mesh->norm[i] += transform_normal(xform, mesh->skinning->rest_norm[i])*w;
            }
            // normalize normal
            mesh->norm[i] = normalize(mesh->norm[i]);
        }
    }
}

// scene reset
void animate_reset(Scene* scene) {
    scene->animation->time = 0;
    for(auto mesh : scene->meshes) {
        if(mesh->animation) {
            mesh->frame = mesh->animation->rest_frame;
        }
        if(mesh->skinning) {
            mesh->pos = mesh->skinning->rest_pos;
            mesh->norm = mesh->skinning->rest_norm;
        }
        if(mesh->simulation) {
            mesh->pos = mesh->simulation->init_pos;
            mesh->simulation->vel = mesh->simulation->init_vel;
            mesh->simulation->force.resize(mesh->simulation->init_pos.size());
        }
    }
}

// scene update
//void animate_update(Scene* scene) {
//    if(scene->animation->time >= scene->animation->length-1) {
//        if(scene->animation->loop) animate_reset(scene);
//        else return;
//    } else scene->animation->time ++;
//    animate_frame(scene);
//    if(not scene->animation->gpu_skinning) animate_skin(scene);
//    simulate(scene);
//}
