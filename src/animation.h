#ifndef _ANIMATION_H_
#define _ANIMATION_H_

#include "scene.h"

// keyframe animation
void animate_frame(Scene* scene, float time);

// skinning scene
void animate_skin(Scene* scene, int time);

// particle simulation
//void simulate(Scene* scene);

// scene reset
void animate_reset(Scene* scene);

// scene update
void animate_update(Scene* scene);

//void Tetrahedralize(Mesh* mesh);

#endif
