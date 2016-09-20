#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CMU462 { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
    Vector3D ori = r.o;
    Vector3D dir = r.d;
    
    Vector3D o = this->o;
    
    Vector3D L = o-ori;
    double tca = dot(L , dir);
    double d2 = dot(L, L) - tca*tca;
    double thc = sqrt(r2 - d2);
    double t11 = tca - thc;
    double t22 = tca + thc;
    
    if (t11 > t22) std::swap(t11, t22);
    if (t11>t1) {
        t1 = t11;
    }
    if (t22<t2) {
        t2 = t22;
    }
    
    if (d2 <= this->r2) {
        return true;
    }else{
        return false;
    }
}

bool Sphere::intersect(const Ray& r) const {

  // TODO:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
    Vector3D ori = r.o;
    Vector3D dir = r.d;
    Vector3D o = this->o;
    
    Vector3D L = o-ori;
    double tca = dot(L , dir);
    double d2 = dot(L, L) - tca*tca;
    double thc = sqrt(r2 - d2);
    double t11 = tca - thc;
    double t22 = tca + thc;
    
    if (t11 > t22) std::swap(t11, t22);
    if (t11<0) {
        t11 = t22;
        if (t11 < 0) {
            return false;
        }
    }
    
    if (d2 > r2) {
        return false;
    }
    
    return true;
}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO:
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
    Vector3D ori = r.o;
    Vector3D dir_o = r.d;
    Vector3D dir = dir_o.unit();
    //Vector3D o = this->o;
    
    Vector3D L = o-ori;
    double tca = dot(L , dir);
    if (tca<0) return false;
    double d2 = dot(L, L) - tca*tca;
    double thc = sqrt(r2 - d2);
    double t11 = tca - thc;
    double t22 = tca + thc;
    i->primitive = this;
    i->bsdf = this->get_bsdf();
    
    if (t11 > t22) std::swap(t11, t22);
    
    if (t11<0) {
        t11 = t22;
        if (t11<0) {
            return false;
        }
    }
    
    i->t = t11;
    i->n = ((ori + t11 * dir)-o).unit();
    
    //if (dot(r.d, i->n) > 0) i->n = -i->n;
    
    if (d2 > r2) {
        return false;
    }
    
    return true;
}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CMU462
