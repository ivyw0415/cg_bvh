#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CMU462 {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();
    
    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    return albedo * (1.0 / PI) * wi.z;
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    *wi = sampler.get_sample();
    *pdf = 1.0 / 2.0 / PI;
    
  return albedo * (1.0 / PI) * wi->z;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  //return Spectrum(0,0,0);
    if (dot(wo, Vector3D(0,0,1))>=0 && wo.z == wi.z && wo.x + wi.x == 0 && wo.y + wi.y == 0) {
        return reflectance;
    }
    return Spectrum(0,0,0);
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
  // Implement MirrorBSDF
    reflect(wo, wi);
    *pdf = 1;
    return reflectance;
}

// Glossy BSDF //

/*
Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0f;
  return reflect(wo, wi, reflectance);
}
*/

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
  // Implement RefractionBSDF
    return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    float rs = 0;
    float rp = 0;
    
    if (wo.unit().z == wi.unit().z && wo.unit().x + wi.unit().x == 0 && wo.unit().y + wi.unit().y == 0) {
        //std::cout<<":<"<<std::endl;
        double cosi = fabs(wo.unit().z);
        double etai = 1;
        double etat = ior;
        double eta = etai/etat;
        double cost = sqrt(1-(eta * eta)*(1-cosi*cosi));
        rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        double fresnel = (rs*rs + rp*rp)/2.0;
        return reflectance * fresnel;
    } else if (wo.z * wi.z < 0){
        if (wi.unit().x*ior + wo.unit().x == 0 || wo.unit().y + wi.unit().y*ior == 0 ){
            double cosi = fabs(wo.unit().z);
            double etai = 1;
            double etat = ior;
            double cost = fabs(wi.unit().z);
            rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        } else if (wo.unit().x*ior + wi.unit().x == 0 || wi.unit().y + wo.unit().y*ior == 0){
            double cosi = fabs(wi.unit().z);
            double etat = ior;
            double etai = 1;
            double cost = fabs(wo.unit().z);
            rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        } else {
            //std::cout<<":<"<<std::endl;
            return Spectrum(0,0,0);
        }
        double fresnel = (rs*rs + rp*rp)/2.0;
        return reflectance * (1-fresnel);

    } else {
        return Spectrum(0,0,0);
    }
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
//  // Compute Fresnel coefficient and either reflect or refract based on it.
    float rs = 0;
    float rp = 0;
    //std::cout<<wo.z<<std::endl;
    if(!refract(wo, wi, ior)){
        reflect(wo,wi);
        *pdf = 1;
        return reflectance;
    }else if (wo.z < 0) {
            refract(wo, wi, ior);
            double cost = fabs(wo.unit().z);
            double etai = 1;
            double etat = ior;
            double eta = etat/etai;
            double cosi = sqrt( 1 - (eta*eta)*(1-cost*cost));
            rs = ((etat * cost) - (etai * cosi)) / ((etat * cost) + (etai * cosi));
            rp = ((etai * cost) - (etat * cosi)) / ((etai * cost) + (etat * cosi));
            double fresnel = (rs*rs + rp*rp)/2.0;
            double rand = (double)std::rand()/(double)RAND_MAX;
            if (rand<fresnel) {
                *pdf = fresnel;
                *wi = -wo + (2 * dot(wo, Vector3D(0,0,-1))) * Vector3D(0,0,-1);
                return reflectance * fresnel;
            }else{
                *pdf = 1-fresnel;
                return reflectance * *pdf;
            }
        
    } else{
        double cosi = wo.unit().z;
        double etai = 1;
        double etat = ior;
        double eta = etai/etat;
        double cost = sqrt(1 - (eta*eta)*(1-cosi*cosi));
        rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        double fresnel = (rs*rs + rp*rp)/2.0;
        *pdf = fresnel;
        double rand = (double)std::rand()/(double)RAND_MAX;
        if (rand<fresnel) {
            *wi = -wo + (2 * dot(wo, Vector3D(0,0,1))) * Vector3D(0,0,1);
            return reflectance * fresnel;
        } else {
            if (fabs(cosi) == 1) {
                *wi = -wo;
                *pdf = 1-fresnel;
            }else{
                refract(wo, wi, ior);
                *pdf = 1-fresnel;
            }
            return reflectance * *pdf;
        }
    }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
    Vector3D n = Vector3D(0,0,1);
    *wi = -wo + (2 * dot(wo, n)) * n;
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
    Vector3D n = Vector3D(0,0,1);
    double cosi = wo.z;
    double etai = 1.0;
    double etat = ior;
    double eta = etai/etat;
    
    
    if (cosi>0) {
        float k = 1-eta*eta*(1-cosi*cosi);
        *wi = eta * (-wo) + (eta*(cosi) - sqrt(k))*n;
        return true;
    } else {
        //std::cout<<wo.z<<" "<<wi->z<<std::endl;
        std::swap(etai, etat);
        cosi = -cosi;
        eta = etai/etat;
        double cost = sqrt(1 - (eta*eta)*(1-cosi*cosi));
        n = -n;
        float k = 1-eta*eta*(1-cosi*cosi);
        if (k<0) {
           return false;
        }
        *wi = eta * (-wo) + (eta*(cosi) - sqrt(k))*n;
        return true;
    }
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
//  *wi  = sampler.get_sample(pdf);
    *wi = sampler.get_sample();
    *pdf = 1.0/2.0/PI;
  return Spectrum();
}

} // namespace CMU462
