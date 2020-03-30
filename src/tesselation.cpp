#include "tesselation.h"

// make normals for each face - duplicates all vertex data
void facet_normals(Mesh* mesh) {
    // allocates new arrays
    auto pos = vector<vec3f>();
    auto norm = vector<vec3f>();
    auto texcoord = vector<vec2f>();
    auto triangle = vector<vec3i>();
    auto quad = vector<vec4i>();
    // froeach triangle
    for(auto f : mesh->triangle) {
        // grab current pos size
        auto nv = (int)pos.size();
        // compute face face normal
        auto fn = normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x]));
        // add triangle
        triangle.push_back({nv,nv+1,nv+2});
        // add vertex data
        for(auto i : range(3)) {
            pos.push_back(mesh->pos[f[i]]);
            norm.push_back(fn);
            if(not mesh->texcoord.empty()) texcoord.push_back(mesh->texcoord[f[i]]);
        }
    }
    // froeach quad
    for(auto f : mesh->quad) {
        // grab current pos size
        auto nv = (int)pos.size();
        // compute face normal
        auto fn = normalize(normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x])) +
                            normalize(cross(mesh->pos[f.z]-mesh->pos[f.x], mesh->pos[f.w]-mesh->pos[f.x])));
        // add quad
        quad.push_back({nv,nv+1,nv+2,nv+3});
        // add vertex data
        for(auto i : range(4)) {
            pos.push_back(mesh->pos[f[i]]);
            norm.push_back(fn);
            if(not mesh->texcoord.empty()) texcoord.push_back(mesh->texcoord[f[i]]);
        }
    }
    // set back mesh data
    mesh->pos = pos;
    mesh->norm = norm;
    mesh->texcoord = texcoord;
    mesh->triangle = triangle;
    mesh->quad = quad;
}

// smooth out normal - does not duplicate data
void smooth_normals(Mesh* mesh) {
    // set normals array to the same length as pos and init all elements to zero
    mesh->norm = vector<vec3f>(mesh->pos.size(),zero3f);
    // foreach triangle
    for(auto f : mesh->triangle) {
        // compute face normal
        auto fn = normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x]));
        // accumulate face normal to the vertex normals of each face index
        for (auto i : range(3)) mesh->norm[f[i]] += fn;
    }
    // foreach quad
    for(auto f : mesh->quad) {
        // compute face normal
        auto fn = normalize(normalize(cross(mesh->pos[f.y]-mesh->pos[f.x], mesh->pos[f.z]-mesh->pos[f.x])) +
                            normalize(cross(mesh->pos[f.z]-mesh->pos[f.x], mesh->pos[f.w]-mesh->pos[f.x])));
        // accumulate face normal to the vertex normals of each face index
        for (auto i : range(4)) mesh->norm[f[i]] += fn;
    }
    // normalize all vertex normals
    for (auto& n : mesh->norm) n = normalize(n);
}

// smooth out tangents
void smooth_tangents(Mesh* polyline) {
    // set tangent array
    polyline->norm = vector<vec3f>(polyline->pos.size(),zero3f);
    // foreach line
    for(auto l : polyline->line) {
        // compute line tangent
        auto lt = normalize(polyline->pos[l.y]-polyline->pos[l.x]);
        // accumulate segment tangent to vertex tangent on each vertex
        for (auto i : range(2)) polyline->norm[l[i]] += lt;
    }
    // normalize all vertex tangents
    for (auto& t : polyline->norm) t = normalize(t);
}

// apply Catmull-Clark mesh subdivision
// does not subdivide texcoord
void subdivide_catmullclark(Mesh* subdiv) {
    // skip is needed
    if(not subdiv->subdivision_catmullclark_level) return;
    // allocate a working Mesh copied from the subdiv
    auto mesh = new Mesh(*subdiv);
    // foreach level
    for(auto l : range(subdiv->subdivision_catmullclark_level)) {
        // make empty pos and quad arrays
        auto pos = vector<vec3f>();
        auto quad = vector<vec4i>();
        // create edge_map from current mesh
        auto edge_map = EdgeMap(mesh->triangle,mesh->quad);
        // linear subdivision - create vertices
        // copy all vertices from the current mesh
        for(auto p : mesh->pos) pos.push_back(p);
        // add vertices in the middle of each edge (use EdgeMap)
        for(auto e : edge_map.edges()) pos.push_back((mesh->pos[e.x]+mesh->pos[e.y])/2);
        // add vertices in the middle of each triangle
        for(auto f : mesh->triangle) pos.push_back((mesh->pos[f.x]+mesh->pos[f.y]+mesh->pos[f.z])/3);
        // add vertices in the middle of each quad
        for(auto f : mesh->quad) pos.push_back((mesh->pos[f.x]+mesh->pos[f.y]+mesh->pos[f.z]+mesh->pos[f.w])/4);
        // subdivision pass --------------------------------
        // compute an offset for the edge vertices
        auto evo = (int)mesh->pos.size();
        // compute an offset for the triangle vertices
        auto tvo = evo + (int)edge_map.edges().size();
        // compute an offset for the quad vertices
        auto qvo = tvo + (int)mesh->triangle.size();
        // foreach triangle
        for(auto fi : range(mesh->triangle.size())) {
            auto f = mesh->triangle[fi];
            // add three quads to the new quad array
            for(auto i : range(3))
                quad.push_back({f[i],evo+edge_map.edge_index({f[i],f[(i+1)%3]}),
                    tvo+fi,evo+edge_map.edge_index({f[i],f[(i+2)%3]})});
        }
        // foreach quad
        for(auto fi : range(mesh->quad.size())) {
            auto f = mesh->quad[fi];
            // add four quads to the new quad array
            for(auto i : range(4))
                quad.push_back({f[i],evo+edge_map.edge_index({f[i],f[(i+1)%4]}),
                    qvo+fi,evo+edge_map.edge_index({f[i],f[(i+3)%4]})});
        }
        // averaging pass ----------------------------------
        // create arrays to compute pos averages (avg_pos, avg_count)
        // arrays have the same length as the new pos array, and are init to zero
        auto avg_pos = vector<vec3f>(pos.size(),zero3f);
        auto avg_count = vector<int>(pos.size(),0);
        // for each new quad
        for(auto f : quad) {
            // compute quad center using the new pos array
            auto fc = (pos[f.x]+pos[f.y]+pos[f.z]+pos[f.w])/4;
            // foreach vertex index in the quad
            for(auto i : range(4)) {
                // accumulate face center to the avg_pos and add 1 to avg_count
                avg_pos[f[i]] += fc;
                avg_count[f[i]] += 1;
            }
        }
        // normalize avg_pos with its count avg_count
        for(auto i : range(pos.size())) avg_pos[i] /= avg_count[i];
        // correction pass ----------------------------------
        // foreach pos, compute correction p = p + (avg_p - p) * (4/avg_count)
        for(auto i : range(pos.size())) pos[i] = pos[i] + (avg_pos[i] - pos[i]) * (4.0f / avg_count[i]);
        // set new arrays pos, quad back into the working mesh; clear triangle array
        mesh->pos = pos;
        mesh->triangle = vector<vec3i>();
        mesh->quad = quad;
    }
    // clear subdivision
    mesh->subdivision_catmullclark_level = 0;
    // according to smooth, either smooth_normals or facet_normals
    if(subdiv->subdivision_catmullclark_smooth) smooth_normals(mesh);
    else facet_normals(mesh);
    // copy back
    *subdiv = *mesh;
    // clear
    delete mesh;
}

// subdivide bezier spline into line segments (assume bezier has only bezier segments and no lines)
void subdivide_bezier(Mesh* bezier) {
    // skip is needed
    if(not bezier->subdivision_bezier_level) return;
    // allocate a working polyline from bezier
    auto lines = new Mesh(*bezier);
    
//    for(auto segment : lines->spline) {
//        lines->line.push_back({segment.x,segment.y});
//        lines->line.push_back({segment.y,segment.z});
//        lines->line.push_back({segment.z,segment.w});
//    }
    // foreach level
    for(auto l : range(bezier->subdivision_bezier_level)) {
        // make new arrays of positions and bezier segments
        auto pos = vector<vec3f>();
        auto segments = vector<vec4i>();
        // copy all the vertices into the new array (this waste space but it is easier for now)
        for(auto p : lines->pos) pos.push_back(p);
        // foreach bezier segment
        for (auto s : lines->spline) {
            // apply subdivision algorithm
            // prepare indices for two new segments
            auto s0 = vec4i{s.x,0,0,(int)pos.size()};
            auto s1 = vec4i{(int)pos.size(),0,0,s.w};
            // add mid point
            pos.push_back(lines->pos[s.x]*(1.f/8)+lines->pos[s.y]*(3.f/8)+lines->pos[s.z]*(3.f/8)+lines->pos[s.w]*(1.f/8));
            // add points for first segment and fix segment indices
            s0.y = pos.size(); s0.z = pos.size()+1;
            pos.push_back(lines->pos[s.x]*(1.f/2)+lines->pos[s.y]*(1.f/2));
            pos.push_back(lines->pos[s.x]*(1.f/4)+lines->pos[s.y]*(1.f/2)+lines->pos[s.z]*(1.f/4));
            // add points for second segment and fix segment indices
            s1.y = pos.size(); s1.z = pos.size()+1;
            pos.push_back(lines->pos[s.y]*(1.f/4)+lines->pos[s.z]*(1.f/2)+lines->pos[s.w]*(1.f/4));
            pos.push_back(lines->pos[s.z]*(1.f/2)+lines->pos[s.w]*(1.f/2));
            // add indices for both segments into new segments array
            segments.push_back(s0);
            segments.push_back(s1);
        }
        // set new arrays pos, segments into the working lineset
        lines->pos = pos;
        lines->spline = segments;
    }
    // copy bezier segments into line segments
    lines->line.clear();
    for(auto segment : lines->spline) {
        lines->line.push_back({segment.x,segment.y});
        lines->line.push_back({segment.y,segment.z});
        lines->line.push_back({segment.z,segment.w});
    }
    // clear bezier array from lines
    lines->spline.clear();
    lines->subdivision_bezier_level = 0;
    // run smoothing to get proper tangents
    smooth_tangents(lines);
    // copy back
    *bezier = *lines;
    // clear
    delete lines;
}

Mesh* make_surface_mesh(frame3f frame, float radius, bool isquad, Material* mat, float offset) {
    auto mesh = new Mesh{};
    mesh->frame = frame;
    mesh->mat = mat;
    if(isquad) {
        mesh->pos = { {-radius,-radius,-offset}, {radius,-radius,-offset},
            {radius,radius,-offset}, {-radius,radius,-offset} };
        mesh->norm = {z3f,z3f,z3f,z3f};
        mesh->quad = { {0,1,2,3} };
    } else {
        map<pair<int,int>,int> vid;
        for(auto j : range(64+1)) {
            for(auto i : range(128+1)) {
                auto u = 2 * pif * i / 64.0f, v = pif * j / 32.0f;
                auto d = vec3f{cos(u)*sin(v),sin(u)*sin(v),cos(v)};
                vid[{i,j}] = mesh->pos.size();
                mesh->pos.push_back(d*radius*(1-offset));
                mesh->norm.push_back(d);
            }
        }
        for(auto j : range(64)) {
            for(auto i : range(128)) {
                mesh->quad.push_back({vid[{i,j}],vid[{i+1,j}],vid[{i+1,j+1}],vid[{i,j+1}]});
            }
        }
    }
    return mesh;
}

void subdivide_surface(Surface* surface) {
    surface->_display_mesh = make_surface_mesh(
        surface->frame, surface->radius, surface->isquad, surface->mat);
}

void subdivide(Scene* scene) {
    for(auto mesh : scene->meshes) {
        if(mesh->subdivision_catmullclark_level) subdivide_catmullclark(mesh);
        if(mesh->subdivision_bezier_level) subdivide_bezier(mesh);
    }
    for(auto surface : scene->surfaces) {
        subdivide_surface(surface);
    }
}
