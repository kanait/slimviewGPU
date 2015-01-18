#ifndef BASISFUNCTION
#define BASISFUNCTION

#include <cmath>
#include <iostream>
using namespace std;

#define EE16 2.1653645317858030703

class BasisFunction {
  
public:

  float centerX, centerY, centerZ;
  float support;
  
  bool leaf;
  int level;
  int visit;

  int id_;
  
  int childN;
  BasisFunction** child;
  
  virtual inline bool rayIntersectLAZ(float &z, float x, float y) = 0;
  virtual inline void normalLA(float g[3], float x, float y, float z) = 0;
  virtual inline void gradientLA(float g[3], float x, float y, float z) = 0;
  virtual inline void gradient12LA(float g[3], float gg[6], float x, float y, float z) = 0;
  virtual inline void gradient123LA(float g[3], float gg[6], float ggg[9], float x, float y, float z) = 0;
  virtual inline float normalZAtC() = 0;
  virtual inline void flipOrientation() = 0;
  
  BasisFunction(){
    leaf = true;
    child = NULL;
    childN = 0;
    id_ = -1;
    
    visit = 0;
  }
  
  virtual ~BasisFunction(){
    if(child != NULL)
      delete[] child;
  }
  
  void setChild(BasisFunction** bfs, int bfN){
    child = bfs;
    childN = bfN;
    leaf = false;
  }
  
  double weight(float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    
    return weight2(vx*vx+vy*vy+vz*vz, (double)support*support);
  }
  
  static inline double weight2(double d, double R){
    if(R < d)
      return 0;
    else{
      /*
      d = 1.5*sqrt(d/R);
      if(d < 0.5)
        return (-d*d + 0.75);
      else
        return (0.5*(1.5-d)*(1.5-d));
       */
      d /= R;
      if(d < 0.25f)
        return exp(-8.0*d);
      else{
        d = 1.0 - sqrt(d);
        d = d*d;
        return EE16*d*d;
      }
    }
  }
  
  static inline double weight(double d, float R){
    //spline patched Gaussian
    if(R < d)
      return 0;
    else{
      /*
      d = 1.5*(d/R);
      if(d < 0.5)
        return (-d*d + 0.75);
      else
        return (0.5*(1.5-d)*(1.5-d));
      */
      d /= R;
      if(d < 0.5f)
        return exp(-8*d*d);
      else{
        d = 1.0 - d;
        d = d*d;
        return EE16*d*d;
      }
    }
  }
  
  inline bool rayIntersectSuppZ(float p[], float start[]){
    float hx = start[0] - centerX;
    float hy = start[1] - centerY;
    float D = support*support - hx*hx - hy*hy;
    if(D < 0){
      return false;
    }
    float k = (float)sqrt(D);
    p[0] = start[0];
    p[1] = start[1];
    p[2] = centerZ-k;
    return true;
  }
  
  bool rayIntersectSupp(float p[], float start[], float dire[]){
    float dot = dire[0]*(centerX - start[0])
              + dire[1]*(centerY - start[1])
              + dire[2]*(centerZ - start[2]);
    float hx = dot*dire[0] + start[0] - centerX;
    float hy = dot*dire[1] + start[1] - centerY;
    float hz = dot*dire[2] + start[2] - centerZ;
    float D = support*support - hx*hx - hy*hy - hz*hz;
    if(D < 0){
      return false;
    }
    float k = (float)sqrt(D);
    p[0] = hx+centerX-k*dire[0];
    p[1] = hy+centerY-k*dire[1];
    p[2] = hz+centerZ-k*dire[2];
    return true;
  }
  
  inline void normalSupp(float n[], float x, float y, float z){
    float nx = x - centerX;
    float ny = y - centerY;
    float nz = z - centerZ;
    float len = -(float)(1.0/sqrt(nx*nx + ny*ny + nz*nz));
    n[0] = nx*len;
    n[1] = ny*len;
    n[2] = nz*len;
  }
};

class Cubic : public BasisFunction{
public:
  float cXXX, cYYY, cZZZ, cXXY, cYYZ, cZZX, cXYY, cYZZ, cZXX, cXYZ;
  float cXX, cYY, cZZ, cXY, cYZ, cZX, cX, cY, cZ, c0;
  
  inline void flipOrientation(){
    cXXX *= -1;
    cYYY *= -1;
    cZZZ *= -1;
    cXXY *= -1;
    cYYZ *= -1;
    cZZX *= -1;
    cXYY *= -1;
    cYZZ *= -1;
    cZXX *= -1;
    cXYZ *= -1;
    
    cXX *= -1;
    cYY *= -1;
    cZZ *= -1;
    cXY *= -1;
    cYZ *= -1;
    cZX *= -1;
    
    cX *= -1;
    cY *= -1;
    cZ *= -1;
    c0 *= -1;
  }
  
#define DEG120 2.0943951023931954923
  
  inline bool rayIntersectLAZ(float &z, float x, float y){
    float sx = x - centerX;
    float sy = y - centerY;
    float dxy2 = sx*sx + sy*sy;
    float sup2 = support*support;
    
    if(cZZZ != 0){
      double a = (cZZX*sx + cYZZ*sy + cZZ)/cZZZ;
      double b = ((cZXX*sx + cXYZ*sy + cZX)*sx + (cYYZ*sy + cYZ)*sy + cZ)/cZZZ;
      double c = (((cXXX*sx + cXXY*sy + cXX)*sx + cXY*sy + cX)*sx + ((cYYY*sy + cXYY*sx + cYY)*sy + cY)*sy + c0)/cZZZ;
      
      double Q = (a*a - 3.0*b)/9.0;
      double R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
      double Q3 = Q*Q*Q;
      double a3 = a/3.0;
      double S = R*R - Q3;
      if(S < 0){
        double s = acos(R/sqrt(Q3))/3.0;
        double Q2 = -2.0*sqrt(Q);
        
        double t1 = Q2*cos(s) - a3;
        double t2 = Q2*cos(s+DEG120) - a3;
        double t3 = Q2*cos(s-DEG120) - a3;
        
        double tmp;
        if(fabs(t1) > fabs(t2)){
          tmp = t1;
          t1 = t2;
          t2 = tmp;
        }
        if(fabs(t2) > fabs(t3)){
          tmp = t2;
          t2 = t3;
          t3 = tmp;
        }
        if(fabs(t1) > fabs(t2)){
          tmp = t1;
          t1 = t2;
          t2 = tmp;
        }
        
        z = (float)t1;
        if(dxy2 + z*z > sup2){
          z = (float)t2;
          if(dxy2 + z*z > sup2){
            z = (float)t3;
            if(dxy2 + z*z > sup2)
              return false;
          }
        }
      }
      else{
        double A = pow(fabs(R) + sqrt(S), 1.0/3.0);
        if(R > 0)
          A = -A;
        double t;
        if((float)A != 0)
          t = (A + Q/A) - a3;
        else
          t = A - a3;
        
        z = (float)t;
        
        if(dxy2 + z*z > sup2)
          return false;
      }
    }
    else if(cZZ != 0){
      float a = cZZ;
      float b = cZX*sx + cYZ*sy + cZ;
      float c = (cXX*sx + cXY*sy + cX)*sx + (cYY*sy + cY)*sy + c0;
      
      float D = b*b - 4.0f*a*c;
      if(D < 0)
        return false;
      
      double q;
      if(b > 0)
        q = -0.5*(b + sqrt(D));
      else
        q = -0.5*(b - sqrt(D));
      float t1 = (float)(q/a);
      float t2 = (float)(c/q);
      if(fabs(t1) > fabs(t2)){
        float tmp = t1;
        t1 = t2;
        t2 = tmp;
      }
      
      z = t1;
      if(dxy2 + z*z > sup2){
        z = t2;
        if(dxy2 + z*z > sup2)
          return false;
      }
    }
    else{
      float a = cZX*sx + cYZ*sy + cZ;
      float b = (cXX*sx + cXY*sy + cX)*sx + (cYY*sy + cY)*sy + c0;
      
      float t = -b/a;
      z = t;
      
      if(dxy2 + z*z > sup2)
        return false;
    }
    z += centerZ;
    return true;
  }
  
#undef DEG120
  
  inline void normalLA(float g[3], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    
    g[0] = (3.0f*cXXX*vx + 2.0f*(cXXY*vy + cXX))*vx + (cXYZ*vz + cXYY*vy + cXY)*vy + (2.0f*cZXX*vx + cZZX*vz +  cZX)*vz + cX;
    g[1] = (3.0f*cYYY*vy + 2.0f*(cYYZ*vz + cYY))*vy + (cXYZ*vx + cYZZ*vz + cYZ)*vz + (2.0f*cXYY*vy + cXXY*vx +  cXY)*vx + cY;
    g[2] = (3.0f*cZZZ*vz + 2.0f*(cZZX*vx + cZZ))*vz + (cXYZ*vy + cZXX*vx + cZX)*vx + (2.0f*cYZZ*vz + cYYZ*vy +  cYZ)*vy + cZ;
    
    float len = (float)(1.0/sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]));
    g[0] *= len;
    g[1] *= len;
    g[2] *= len;
  }
  
  inline void gradientLA(float g[3], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    g[0] = (3.0f*cXXX*vx + 2.0f*(cXXY*vy + cXX))*vx + (cXYZ*vz + cXYY*vy + cXY)*vy + (2.0f*cZXX*vx + cZZX*vz +  cZX)*vz + cX;
    g[1] = (3.0f*cYYY*vy + 2.0f*(cYYZ*vz + cYY))*vy + (cXYZ*vx + cYZZ*vz + cYZ)*vz + (2.0f*cXYY*vy + cXXY*vx +  cXY)*vx + cY;
    g[2] = (3.0f*cZZZ*vz + 2.0f*(cZZX*vx + cZZ))*vz + (cXYZ*vy + cZXX*vx + cZX)*vx + (2.0f*cYZZ*vz + cYYZ*vy +  cYZ)*vy + cZ;
    
  }
  
  inline void gradient12LA(float g[3], float gg[6], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    g[0] = (3.0f*cXXX*vx + 2.0f*(cXXY*vy + cXX))*vx + (cXYZ*vz + cXYY*vy + cXY)*vy + (2.0f*cZXX*vx + cZZX*vz +  cZX)*vz + cX;
    g[1] = (3.0f*cYYY*vy + 2.0f*(cYYZ*vz + cYY))*vy + (cXYZ*vx + cYZZ*vz + cYZ)*vz + (2.0f*cXYY*vy + cXXY*vx +  cXY)*vx + cY;
    g[2] = (3.0f*cZZZ*vz + 2.0f*(cZZX*vx + cZZ))*vz + (cXYZ*vy + cZXX*vx + cZX)*vx + (2.0f*cYZZ*vz + cYYZ*vy +  cYZ)*vy + cZ;
    
    gg[0] = 6*cXXX*vx + 2*(cXXY*vy + cXX + cZXX*vz);
    gg[1] = cXYZ*vz + 2*(cXXY*vx + cXYY*vy) + cXY;
    gg[2] = cXYZ*vy + 2*(cZXX*vx + cZZX*vz) + cZX;
    gg[3] = 6*cYYY*vy + 2*(cYYZ*vz + cYY + cXYY*vx);
    gg[4] = cXYZ*vx + 2*(cYYZ*vy + cYZZ*vz) + cYZ;
    gg[5] = 6*cZZZ*vz + 2*(cZZX*vx + cZZ + cYZZ*vy);
  }
  
  inline void gradient123LA(float g[3], float gg[6], float ggg[9], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    g[0] = (3.0f*cXXX*vx + 2.0f*(cXXY*vy + cXX))*vx + (cXYZ*vz + cXYY*vy + cXY)*vy + (2.0f*cZXX*vx + cZZX*vz +  cZX)*vz + cX;
    g[1] = (3.0f*cYYY*vy + 2.0f*(cYYZ*vz + cYY))*vy + (cXYZ*vx + cYZZ*vz + cYZ)*vz + (2.0f*cXYY*vy + cXXY*vx +  cXY)*vx + cY;
    g[2] = (3.0f*cZZZ*vz + 2.0f*(cZZX*vx + cZZ))*vz + (cXYZ*vy + cZXX*vx + cZX)*vx + (2.0f*cYZZ*vz + cYYZ*vy +  cYZ)*vy + cZ;
    
    gg[0] = 6*cXXX*vx + 2*(cXXY*vy + cXX + cZXX*vz);
    gg[1] = cXYZ*vz + 2*(cXXY*vx + cXYY*vy) + cXY;
    gg[2] = cXYZ*vy + 2*(cZXX*vx + cZZX*vz) + cZX;
    gg[3] = 6*cYYY*vy + 2*(cYYZ*vz + cYY + cXYY*vx);
    gg[4] = cXYZ*vx + 2*(cYYZ*vy + cYZZ*vz) + cYZ;
    gg[5] = 6*cZZZ*vz + 2*(cZZX*vx + cZZ + cYZZ*vy);
    
    ggg[0] = 6*cXXX;
    ggg[1] = 2*cXXY;
    ggg[2] = 2*cZXX;
    ggg[3] = 2*cXYY;
    ggg[4] = cXYZ;
    ggg[5] = 2*cZZX;
    ggg[6] = 6*cYYY;
    ggg[7] = 2*cYYZ;
    ggg[8] = 2*cYZZ;
    ggg[9] = 6*cZZZ;
  }
  
  inline float normalZAtC(){
    return cZ;
  }
};

class Quadratic : public BasisFunction{
public:
  float cXX, cYY, cZZ, cXY, cYZ, cZX, cX, cY, cZ, c0;
  
  inline void flipOrientation(){
    cXX *= -1;
    cYY *= -1;
    cZZ *= -1;
    cXY *= -1;
    cYZ *= -1;
    cZX *= -1;
    
    cX *= -1;
    cY *= -1;
    cZ *= -1;
    c0 *= -1;
  }
  
  inline bool rayIntersectLAZ(float &z, float x, float y){
    float sx = x - centerX;
    float sy = y - centerY;
    float dxy2 = sx*sx + sy*sy;
    float sup2 = support*support;
    
    if(cZZ != 0){
      float a = cZZ;
      float b = cZX*sx + cYZ*sy + cZ;
      float c = (cXX*sx + cXY*sy + cX)*sx + (cYY*sy + cY)*sy + c0;
      
      float D = b*b - 4.0f*a*c;
      if(D < 0)
        return false;
      
      double q;
      if(b > 0)
        q = -0.5*(b + sqrt(D));
      else
        q = -0.5*(b - sqrt(D));
      float t1 = (float)(q/a);
      float t2 = (float)(c/q);
      if(fabs(t1) > fabs(t2)){
        float tmp = t1;
        t1 = t2;
        t2 = tmp;
      }
      
      z = t1;
      if(dxy2 + z*z > sup2){
        z = t2;
        if(dxy2 + z*z > sup2)
          return false;
      }
    }
    else{
      float a = cZX*sx + cYZ*sy + cZ;
      float b = (cXX*sx + cXY*sy + cX)*sx + (cYY*sy + cY)*sy + c0;
      
      float t = -b/a;
      z = t;
      
      if(dxy2 + z*z > sup2)
        return false;
    }
    z += centerZ;
    return true;
  }
  
  inline void normalLA(float g[3], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    
    g[0] = 2.0f*cXX*vx + cXY*vy + cZX*vz + cX;
    g[1] = 2.0f*cYY*vy + cYZ*vz + cXY*vx + cY;
    g[2] = 2.0f*cZZ*vz + cZX*vx + cYZ*vy + cZ;
    
    float len = (float)(1.0/sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]));
    g[0] *= len;
    g[1] *= len;
    g[2] *= len;
  }
  
  inline void gradientLA(float g[3], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    
    g[0] = 2.0f*cXX*vx + cXY*vy + cZX*vz + cX;
    g[1] = 2.0f*cYY*vy + cYZ*vz + cXY*vx + cY;
    g[2] = 2.0f*cZZ*vz + cZX*vx + cYZ*vy + cZ;
  }
  
  inline void gradient12LA(float g[3], float gg[6], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    
    g[0] = 2.0f*cXX*vx + cXY*vy + cZX*vz + cX;
    g[1] = 2.0f*cYY*vy + cYZ*vz + cXY*vx + cY;
    g[2] = 2.0f*cZZ*vz + cZX*vx + cYZ*vy + cZ;
    
    gg[0] = 2.0f*cXX;
    gg[1] = cXY;
    gg[2] = cZX;
    gg[3] = 2.0f*cYY;
    gg[4] = cYZ;
    gg[5] = 2.0f*cZZ;
  }
  
  inline void gradient123LA(float g[3], float gg[6], float ggg[9], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    
    g[0] = 2.0f*cXX*vx + cXY*vy + cZX*vz + cX;
    g[1] = 2.0f*cYY*vy + cYZ*vz + cXY*vx + cY;
    g[2] = 2.0f*cZZ*vz + cZX*vx + cYZ*vy + cZ;
    
    gg[0] = 2.0f*cXX;
    gg[1] = cXY;
    gg[2] = cZX;
    gg[3] = 2.0f*cYY;
    gg[4] = cYZ;
    gg[5] = 2.0f*cZZ;
    
    ggg[0] = 0;
    ggg[1] = 0;
    ggg[2] = 0;
    ggg[3] = 0;
    ggg[4] = 0;
    ggg[5] = 0;
    ggg[6] = 0;
    ggg[7] = 0;
    ggg[8] = 0;
    ggg[9] = 0;
  }
  
  inline float normalZAtC(){
    return cZ;
  }
};

class Linear : public BasisFunction{
public:
  float cX, cY, cZ, c0;
  
  inline void flipOrientation(){
    cX *= -1;
    cY *= -1;
    cZ *= -1;
    c0 *= -1;
  }
  
  inline bool rayIntersectLAZ(float &z, float x, float y){
    float sx = x - centerX;
    float sy = y - centerY;
    float dxy2 = sx*sx + sy*sy;
    float sup2 = support*support;
    
    float a = cZ;
    float b = cX*sx + cY*sy + c0;
    
    float t = -b/a;
    z = t;
    
    if(dxy2 + z*z > sup2)
      return false;
    
    z += centerZ;
    return true;
  }
  
  inline void normalLA(float g[3], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    
    g[0] = cX;
    g[1] = cY;
    g[2] = cZ;
    
    float len = (float)(1.0/sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]));
    g[0] *= len;
    g[1] *= len;
    g[2] *= len;
  }
  
  inline void gradientLA(float g[3], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    
    g[0] = cX;
    g[1] = cY;
    g[2] = cZ;
  }
  
  inline void gradient12LA(float g[3], float gg[6], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    
    g[0] = cX;
    g[1] = cY;
    g[2] = cZ;
    
    gg[0] = 0;
    gg[1] = 0;
    gg[2] = 0;
    gg[3] = 0;
    gg[4] = 0;
    gg[5] = 0;
  }
  
  inline void gradient123LA(float g[3], float gg[6], float ggg[9], float x, float y, float z){
    float vx = x - centerX;
    float vy = y - centerY;
    float vz = z - centerZ;
    
    g[0] = cX;
    g[1] = cY;
    g[2] = cZ;
    
    gg[0] = 0;
    gg[1] = 0;
    gg[2] = 0;
    gg[3] = 0;
    gg[4] = 0;
    gg[5] = 0;
    
    ggg[0] = 0;
    ggg[1] = 0;
    ggg[2] = 0;
    ggg[3] = 0;
    ggg[4] = 0;
    ggg[5] = 0;
    ggg[6] = 0;
    ggg[7] = 0;
    ggg[8] = 0;
    ggg[9] = 0;
  }
  
  inline float normalZAtC(){
    return cZ;
  }
};

#endif
