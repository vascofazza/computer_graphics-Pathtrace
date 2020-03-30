#ifndef _TESSELATION_H_
#define _TESSELATION_H_

#include "scene.h"

// map used to uniquify edges
struct EdgeMap {
    map<pair<int,int>,int>  _edge_map;  // internal map
    vector<vec2i>           _edge_list; // internal list to generate unique ids
    
    // create an edge map for a collection of triangles and quads
    EdgeMap(vector<vec3i> triangle, vector<vec4i> quad) {
        for(auto f : triangle) { _add_edge(f.x,f.y); _add_edge(f.y,f.z); _add_edge(f.z,f.x); }
        for(auto f : quad) { _add_edge(f.x,f.y); _add_edge(f.y,f.z); _add_edge(f.z,f.w); _add_edge(f.w,f.x); }
    }
    
    // internal function to add an edge
    void _add_edge(int i, int j) {
        if(_edge_map.find(make_pair(i,j)) == _edge_map.end()) {
            _edge_map[make_pair(i,j)] = _edge_list.size();
            _edge_map[make_pair(j,i)] = _edge_list.size();
            _edge_list.push_back(vec2i(i,j));
        }
    }
    
    // edge list
    const vector<vec2i>& edges() const { return _edge_list; }
    
    // get an edge from two vertices
    int edge_index(vec2i e) const {
        error_if_not(not (_edge_map.find(make_pair(e.x,e.y)) == _edge_map.end()), "non existing edge");
        return _edge_map.find(make_pair(e.x, e.y))->second;
    }
};

// set face normals (duplicating vertices)
void facet_normals(Mesh* mesh);

// compute smoothed normals
void smooth_normals(Mesh* mesh);

// compute smoothed line tangents
void smooth_tangents(Mesh* lines);

// subdivide the scene
void subdivide(Scene* scene);

// apply catmull-clark subdivision to the mesh recursively
void subdivide_catmullclark(Mesh* subdiv);

// apply bezier spline subdivision
void subdivide_bezier(Mesh* splines);

// make display meshes for surfaces
void subdivide_surface(Surface* surface);

// make a surface mesh
Mesh* make_surface_mesh(frame3f frame, float radius, bool isquad, Material* mat, float offset = 0);

#endif
