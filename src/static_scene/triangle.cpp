#include "triangle.h"

#include "CMU462/CMU462.h"
#include "GL/glew.h"

namespace CMU462 { namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3) :
    mesh(mesh), v1(v1), v2(v2), v3(v3) { }

BBox Triangle::get_bbox() const {
  
  // TODO: 
  // compute the bounding box of the triangle
    
    Vector3D min = Vector3D(std::min(mesh->positions[v1].x, std::min(mesh->positions[v2].x,mesh->positions[v3].x)),
                            std::min(mesh->positions[v1].y, std::min(mesh->positions[v2].y,mesh->positions[v3].y)),
                            std::min(mesh->positions[v1].z, std::min(mesh->positions[v2].z,mesh->positions[v3].z)));
    Vector3D max = Vector3D(std::max(mesh->positions[v1].x, std::max(mesh->positions[v2].x,mesh->positions[v3].x)),
                            std::max(mesh->positions[v1].y, std::max(mesh->positions[v2].y,mesh->positions[v3].y)),
                            std::max(mesh->positions[v1].z, std::max(mesh->positions[v2].z,mesh->positions[v3].z)));
  return BBox(min,max);
}

bool Triangle::intersect(const Ray& r) const {
    //cout<<"do!"<<endl;
  // TODO: implement ray-triangle intersection
    Vector3D p0 = mesh->positions[v1];
    Vector3D p1 = mesh->positions[v2];
    Vector3D p2 = mesh->positions[v3];
    Vector3D o = r.o;
    Vector3D d = r.d;
    
    Vector3D e1 = p1 - p0;
    Vector3D e2 = p2 - p0;
    Vector3D s = o - p0;
    
    if (dot(cross(e1, d), e2) == 0) {
        return false;
    }
    
    Vector3D intersection = (1/dot(cross(e1, d), e2)) * Vector3D((-1)*dot(cross(s, e2), d),dot(cross(e1, d), s),(-1)*dot(cross(s, e2), e1));
    double u = intersection.x;
    double v = intersection.y;
    double t = intersection.z;
    
    if (u>0 && v>0 && u+v<1 && t<=r.max_t && t>=r.min_t) {
        //cout<<"^^"<<endl;
        return true;
    }else{
        return false;
    }
}

bool Triangle::intersect(const Ray& r, Intersection *isect) const {
  
  // TODO: 
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
   
    Vector3D p0 = mesh->positions[v1];
    Vector3D p1 = mesh->positions[v2];
    Vector3D p2 = mesh->positions[v3];
    Vector3D o = r.o;
    Vector3D d = r.d;
    
    
    Vector3D e1 = p1 - p0;
    Vector3D e2 = p2 - p0;
    Vector3D s = o - p0;
    
    if (dot(cross(e1, d), e2) == 0) {
        return false;
    }
    
    Vector3D intersection = (1/dot(cross(e1, d), e2)) * Vector3D((-1)*dot(cross(s, e2), d),dot(cross(e1, d), s),(-1)*dot(cross(s, e2), e1));
    double u = intersection.x;
    double v = intersection.y;
    double t = intersection.z;
    
    Vector3D n;
    
    n = (1-u-v) * mesh->normals[v1] + u*mesh->normals[v2] + v*mesh->normals[v3];
    
    n.normalize();
    isect->n = n;
    isect->t = t;
    isect->primitive = this;
    isect->bsdf = mesh->get_bsdf();
    
    
    if (u>0 && v>0 && u+v<1 && t<=r.max_t && t>=r.min_t) {
        //cout<<"do!"<<endl;
        if (dot(isect->n, d) > 0) {
            isect->n *= -1;
        }
        return true;
    }else{
        return false;
    }
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}



} // namespace StaticScene
} // namespace CMU462
