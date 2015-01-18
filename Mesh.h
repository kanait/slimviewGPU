#ifndef MESH_H
#define MESH_H

#include <math.h>

class Mesh{
public:
  
  int vertexN;
  int faceN;
  float (*vertex)[3];
  float (*o_vertex)[3];
  int (*face)[3];
  float (*normal)[3];
  float (*o_normal)[3];
  
  Mesh(int vN, int fN){
    vertexN = vN;
    faceN = fN;
    vertex = new float[vN][3];
    o_vertex = new float[vN][3];
    normal = new float[vN][3];
    o_normal = new float[vN][3];
    face = new int[fN][3];
  }
  
  ~Mesh(){
    delete[] vertex;
    delete[] o_vertex;
    delete[] normal;
    delete[] o_normal;
    delete[] face;
  }
  
  inline void setVertex(int i, float x, float y, float z){
    o_vertex[i][0] = x;
    o_vertex[i][1] = y;
    o_vertex[i][2] = z;
  }
  
  inline void setNormal(int i, float x, float y, float z){
    o_normal[i][0] = x;
    o_normal[i][1] = y;
    o_normal[i][2] = z;
  }
  
  inline void setFace(int i, int f0, int f1, int f2){
    face[i][0] = f0;
    face[i][1] = f1;
    face[i][2] = f2;
  }
  
  inline void getFaceBoundX(float b[2], int i){
    int* f = face[i];
    b[0] = b[1] = vertex[f[0]][0];
    
    float p = vertex[f[1]][0];
    if(p < b[0])
      b[0] = p;
    else if(p > b[1])
      b[1] = p;
    
    p = vertex[f[2]][0];
    if(p < b[0])
      b[0] = p;
    else if(p > b[1])
      b[1] = p;
  }
  
  inline void getFaceBoundY(float b[2], int i){
    int* f = face[i];
    b[0] = b[1] = vertex[f[0]][1];
    
    float p = vertex[f[1]][1];
    if(p < b[0])
      b[0] = p;
    else if(p > b[1])
      b[1] = p;
    
    p = vertex[f[2]][1];
    if(p < b[0])
      b[0] = p;
    else if(p > b[1])
      b[1] = p;
  }
  
  inline bool rayIntersectFaceZ(float &z, float n[3],
                                int i, float x, float y){
    int* f = face[i];
    float *p0 = vertex[f[0]];
    float *p1 = vertex[f[1]];
    float *p2 = vertex[f[2]];
    
    float a0 = (p1[0] - x)*(p2[1] - y) - (p1[1] - y)*(p2[0] - x);
    float a1 = (p2[0] - x)*(p0[1] - y) - (p2[1] - y)*(p0[0] - x);
    float a2 = (p0[0] - x)*(p1[1] - y) - (p0[1] - y)*(p1[0] - x);
    
    if(a0*a1 < 0 || a1*a2 < 0 || a2*a0 < 0)
      return false;
    
    z = (a0*p0[2] + a1*p1[2] + a2*p2[2])/(a0+a1+a2);
    
    float *n0 = normal[f[0]];
    float *n1 = normal[f[1]];
    float *n2 = normal[f[2]];
    
    float nx = (a0*n0[0] + a1*n1[0] + a2*n2[0])/(a0+a1+a2);
    float ny = (a0*n0[1] + a1*n1[1] + a2*n2[1])/(a0+a1+a2);
    float nz = (a0*n0[2] + a1*n1[2] + a2*n2[2])/(a0+a1+a2);
    float l = (float)sqrt(nx*nx + ny*ny + nz*nz);
    if(l != 0){
      l = 1.0f/l;
      n[0] = l*nx;
      n[1] = l*ny;
      n[2] = l*nz;
    }
    else
      n[0] = n[1] = n[2] = 0;
	return true;
  }
  
  void computeNormal(){
    int i;
    for(i=0; i<vertexN; i++)
      o_normal[i][0] = o_normal[i][1] = o_normal[i][2] = 0;
    
    for(i=0; i<faceN; i++){
      int* f = face[i];
      float* p0 = o_vertex[f[0]];
      float* p1 = o_vertex[f[1]];
      float* p2 = o_vertex[f[2]];
      
      float nx = (p1[1]-p0[1])*(p2[2]-p0[2]) - (p1[2]-p0[2])*(p2[1]-p0[1]);
      float ny = (p1[2]-p0[2])*(p2[0]-p0[0]) - (p1[0]-p0[0])*(p2[2]-p0[2]);
      float nz = (p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0]);
      
      float* n = o_normal[f[0]];
      n[0] += nx;
      n[1] += ny;
      n[2] += nz;
      
      n = o_normal[f[1]];
      n[0] += nx;
      n[1] += ny;
      n[2] += nz;
      
      n = o_normal[f[2]];
      n[0] += nx;
      n[1] += ny;
      n[2] += nz;
    }
    
    for(i=0; i<vertexN; i++){
      float* n = o_normal[i];
      float l = (float)sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
      if(l != 0){
        n[0] /= l;
        n[1] /= l;
        n[2] /= l;
      }
    }
  }
  
  void bound(float min[3], float max[3]){
    int i, j;
    for(i=0; i<3; i++)
      min[i] = max[i] = o_vertex[0][i];
    for(j=0; j<vertexN; j++){
      float* p = o_vertex[j];
      for(i=0; i<3; i++){
        if(p[i] < min[i])
          min[i] = p[i];
        else if(p[i] > max[i])
          max[i] = p[i];
      }
    }
  }
  
  void rescale(float scale){
    float min[3], max[3];
    bound(min, max);
    float dia = (float)sqrt((max[0]-min[0])*(max[0]-min[0]) +
                            (max[1]-min[1])*(max[1]-min[1]) + 
                            (max[2]-min[2])*(max[2]-min[2]));
    float s = scale/dia;
    int i,j;
    for(i=0; i<3; i++){
      float m = 0.5f*(min[i]+max[i]);
      for(j=0; j<vertexN; j++){
        o_vertex[j][i] = s*(o_vertex[j][i] - m);
      }
    }
  }
};

#endif
