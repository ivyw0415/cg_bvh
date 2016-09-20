#include "environment_light.h"

namespace CMU462 { namespace StaticScene {

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {
  // TODO: initialize things here as needed
    this->envMap = envMap;
}

Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight,
                                    float* pdf) const {
  // TODO: Implement
    *pdf = 1.0/4.0/PI;
    
    double Xi1 = (double)(std::rand()) / RAND_MAX;
    double Xi2 = (double)(std::rand()) / RAND_MAX;
    
    double theta = acos(2*Xi2-1);
    double phi = 2.0 * PI * Xi1;
    
    double xs = sinf(theta)*cosf(phi);
    double ys = sinf(theta) * sinf(phi);
    double zs = cosf(phi);
    
    *wi = Vector3D(xs, ys, zs);
    *distToLight = INF_D;
    
    double x = std::fabs(phi/(2.0*PI) * envMap->w);
    double y = std::fabs((theta)/PI * envMap->h);
    
    int xx = floor(x);
    int yy = floor(y);
    
    double x_ratio = x - xx;
    double y_ratio = y - yy;
    double x_oppo = 1 - x_ratio;
    double y_oppo = 1 - y_ratio;
    
    Spectrum result;
    
    if (xx>=0 && xx<envMap->w-1 && y>=0 && y<envMap->h-1) {
        result = (envMap->data[xx + yy * envMap->w]*x_oppo + envMap->data[(xx+1) + yy * envMap->w]*x_ratio)*y_oppo +
        (envMap->data[xx + (yy+1) * envMap->w]*x_oppo + envMap->data[(xx+1) + (yy+1) * envMap->w]*x_ratio)*y_ratio;
    }else{
        return envMap->data[xx+yy*envMap->w];
    }
    return result;
}

Spectrum EnvironmentLight::sample_dir(const Ray& r) const {
  // TODO: Implement
    double theta = acos(r.d.y);
    double phi = acos(r.d.z);
    
    double x = std::fabs(phi/(2.0*PI) * envMap->w);
    double y = std::fabs((theta)/PI * envMap->h);
    
    int xx = floor(x);
    int yy = floor(y);
    
    double x_ratio = x - xx;
    double y_ratio = y - yy;
    double x_oppo = 1 - x_ratio;
    double y_oppo = 1 - y_ratio;
    
    Spectrum result;
    
    if (xx>=0 && xx<envMap->w-1 && y>=0 && y<envMap->h-1) {
        result = (envMap->data[xx + yy * envMap->w]*x_oppo + envMap->data[(xx+1) + yy * envMap->w]*x_ratio)*y_oppo +
        (envMap->data[xx + (yy+1) * envMap->w]*x_oppo + envMap->data[(xx+1) + (yy+1) * envMap->w]*x_ratio)*y_ratio;
    }else{
        return envMap->data[xx+yy*envMap->w];
    }
    return result;
}

} // namespace StaticScene
} // namespace CMU462
