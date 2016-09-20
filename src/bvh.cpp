#include "bvh.h"

#include "CMU462/CMU462.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <algorithm>
#include <stack>

using namespace std;

namespace CMU462 { namespace StaticScene {
    
    void BVHAccel::construct(BVHNode* root, size_t start, size_t range, size_t max_leaf_size){
        if (range < max_leaf_size) {
            return;
        }
        
        BBox bb = root->bb;
        double xbound_min = bb.min.x;
        double ybound_min = bb.min.y;
        double zbound_min = bb.min.z;
        
        double xdiff = bb.max.x - bb.min.x;
        double ydiff = bb.max.y - bb.min.y;
        double zdiff = bb.max.z - bb.min.z;
        
        std::vector<Primitive *> l_primitives;
        std::vector<Primitive *> r_primitives;
        std::vector<Primitive *> u_primitives;
        
        double xcost = INF_D;
        double ycost = INF_D;
        double zcost = INF_D;
        double minCost = INF_D;
        double part_pos = 0;
        BBox part_b_l = BBox();
        BBox part_b_r = BBox();
        int l_size = 0;
        int r_size = 0;
        
        int axis = 0; //x=1, y=2, z=3
        int partition = 32;
        double sr = root->bb.surface_area();
        
        for (int i=1; i<partition; i++) {
            double xc = xbound_min + i * xdiff / (double)partition;
            double yc = ybound_min + i * ydiff / (double)partition;
            double zc = zbound_min + i * zdiff / (double)partition;
            
            BBox xlbox = BBox();
            BBox xrbox = BBox();
            BBox ylbox = BBox();
            BBox yrbox = BBox();
            BBox zlbox = BBox();
            BBox zrbox = BBox();
        //double xc = xdiff/2 + xbound_min;
            
            int xlbox_size = 0;
            int xrbox_size = 0;
            int ylbox_size = 0;
            int yrbox_size = 0;
            int zlbox_size = 0;
            int zrbox_size = 0;
            
            
            for (size_t j=0; j<range; j++) {
                
                if (primitives[start+j]->get_bbox().centroid().x < xc) {
                    xlbox.expand(primitives[start+j]->get_bbox());
                    xlbox_size++;
                }else{
                    xrbox.expand(primitives[start+j]->get_bbox());
                    xrbox_size++;
                }
                if (primitives[start+j]->get_bbox().centroid().y < yc) {
                    ylbox.expand(primitives[start+j]->get_bbox());
                    ylbox_size++;
                }else{
                    yrbox.expand(primitives[start+j]->get_bbox());
                    yrbox_size++;
                }
                
                if (primitives[start+j]->get_bbox().centroid().z < zc) {
                    zlbox.expand(primitives[start+j]->get_bbox());
                    zlbox_size++;
                }else{
                    zrbox.expand(primitives[start+j]->get_bbox());
                    zrbox_size++;
                }
                 
            }
        
            if (xlbox_size >= max_leaf_size && xrbox_size >= max_leaf_size) {
            //if (xlbox_size != 0 && xrbox_size!=0) {
                xcost = (double)(xlbox_size * xlbox.surface_area() + xrbox_size * xrbox.surface_area())/(double)sr;
            }else{
                xcost = INF_D;
            }
            
            if (ylbox_size >= max_leaf_size && yrbox_size >= max_leaf_size){
            //if (xlbox_size !=0 && xrbox_size!=0)  {
                ycost = (double)(ylbox_size * ylbox.surface_area() + yrbox_size * yrbox.surface_area())/(double)sr;
    
            }else{
                ycost = INF_D;
            }
            if (zlbox_size >= max_leaf_size && zrbox_size >= max_leaf_size){
            //if (xlbox_size !=0 && xrbox_size!=0)  {
                zcost = (double)(zlbox_size * zlbox.surface_area() + zrbox_size * zrbox.surface_area())/(double)sr;
            }else{
                zcost = INF_D;
            }
            
            if (xcost<minCost) {
                minCost = xcost;
                axis = 1;
                part_pos = xc;
                part_b_l = xlbox;
                part_b_r = xrbox;
            }
            if (ycost<minCost) {
                minCost = ycost;
                axis = 2;
                part_pos = yc;
                part_b_l = ylbox;
                part_b_r = yrbox;
            }
            if (zcost<minCost) {
                minCost = zcost;
                axis = 3;
                part_pos = zc;
                part_b_l = zlbox;
                part_b_r = zrbox;
            }
        }
        
        
        for (size_t j=0; j<range; j++) {
            if (axis == 1){
                if (primitives[start+j]->get_bbox().centroid().x < part_pos){
                    l_primitives.push_back(primitives[start+j]);
                    l_size++;
                }else{
                    r_primitives.push_back(primitives[start+j]);
                    r_size++;
                }
            }else if (axis == 2){
                if (primitives[start+j]->get_bbox().centroid().y < part_pos) {
                    l_primitives.push_back(primitives[start+j]);
                    l_size++;
                }else{
                    r_primitives.push_back(primitives[start+j]);
                    r_size++;
                }
            }else if (axis == 3){
                if (primitives[start+j]->get_bbox().centroid().z < part_pos) {
                    l_primitives.push_back(primitives[start+j]);
                    l_size++;
                }else{
                    r_primitives.push_back(primitives[start+j]);
                    r_size++;
                }
            }
            else
                return;
        }
        //cout<<l_primitives.size()+r_primitives.size()<<endl;
        
        for (size_t j = 0; j<l_size; j++) {
            primitives[start+j] = l_primitives[j];
        }
        for (size_t j = 0; j<r_size; j++) {
            primitives[start+l_primitives.size()+j] = r_primitives[j];
        }
        
        //primitives = l_primitives;
        
        root->l = new BVHNode(part_b_l, start, l_size);
        root->r = new BVHNode(part_b_r, start+l_size, r_size);
        //cout<<endl<<l_size+r_size<<" "<<range<<endl;
        
        construct(root->l, start, l_primitives.size(), max_leaf_size);
        construct(root->r, start+l_primitives.size(), r_primitives.size(), max_leaf_size);
        
    }

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  this->primitives = _primitives;

  // TODO:
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  BBox bb;
  for (size_t i = 0; i < primitives.size(); ++i) {
    bb.expand(primitives[i]->get_bbox());
  }

  root = new BVHNode(bb, 0, primitives.size());
    construct(root, 0, primitives.size(), max_leaf_size);
}

BVHAccel::~BVHAccel() {

  // TODO:
  // Implement a proper destructor for your BVH accelerator aggregate
    delete root;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

bool BVHAccel::intersect(const Ray &ray) const {

  // TODO:
  // Implement ray - bvh aggregate intersection test. A ray intersects
  // with a BVH aggregate if and only if it intersects a primitive in
  // the BVH that is not an aggregate.
  bool hit = false;
  for (size_t p = 0; p < primitives.size(); ++p) {
    //if(primitives[p]->intersect(ray)) hit = true;
  }

  //return hit;
    return intersect(ray, root);

}

bool BVHAccel::intersect(const Ray &ray, Intersection *i) const {

  // TODO:
  // Implement ray - bvh aggregate intersection test. A ray intersects
  // with a BVH aggregate if and only if it intersects a primitive in
  // the BVH that is not an aggregate. When an intersection does happen.
  // You should store the non-aggregate primitive in the intersection data
  // and not the BVH aggregate itself.

  bool hit = false;
    
  for (size_t p = 0; p < primitives.size(); ++p) {
    //if(primitives[p]->intersect(ray, i)) hit = true;
  }

  //return hit;
    return intersect(ray, i, root);
}

    bool BVHAccel::intersect(const Ray& r, Intersection* i, BVHNode* root) const{
        double t0 = r.min_t;
        double t1 = r.max_t;
        BBox bb = root->bb;
        
        if (!bb.intersect(r, t0, t1)||t0>i->t) {
            return false;
        }
        
        if (root->isLeaf()) {
            bool hit = false;
            size_t start = root->start;
            size_t range = root->range;
            double min_t = INF_D;
            
            for (size_t p = start; p < start + range; ++p) {
                Intersection isect;
                if(primitives[p]->intersect(r, &isect)&&isect.t<t1&&isect.t>t0){
                    if (hit == false) {
                         hit = true;
                    }
                    if (isect.t<min_t) {
                        i->bsdf = primitives[p]->get_bsdf();
                        i->primitive = primitives[p];
                        i->t = isect.t;
                        i->n = isect.n;
                        min_t = isect.t;
                    }
                }
            }
            return hit;
        } else {
            
            bool r_hit = intersect(r, i, root->r);
            bool l_hit = intersect(r, i, root->l);
            if (!r_hit&&!l_hit) {
                return false;
            } else {
                return true;
            }
        }
    }
    bool BVHAccel::intersect(const Ray& r, BVHNode * root) const{
        double t0 = r.min_t;
        double t1 = r.max_t;
        BBox bb = root->bb;
        
        if (!bb.intersect(r, t0, t1)) {
            return false;
        }
        
        if (root->isLeaf()) {
            bool hit = false;
            size_t start = root->start;
            size_t range = root->range;
            
            for (size_t p = start; p < start + range; ++p) {
                if(primitives[p]->intersect(r)){
                    hit = true;
                    break;
                }
            }
            return hit;
        } else {
            bool r_hit = intersect(r, root->r);
            bool l_hit = intersect(r, root->l);
            if (!r_hit&&!l_hit) {
                return false;
            } else {
                return true;
            }
        }
    }
    
}  // namespace StaticScene
}  // namespace CMU462
