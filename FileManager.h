#ifndef FILEMANAGER
#define FILEMANAGER

#include "BasisFunction.h"
#include "Mesh.h"

#include <cstdio>
#include <fstream>
using namespace std;

#include <gl\gl.h>
#include <gl\glu.h>

class FileManager{
public:
  static Mesh* readPolyFile( char* name ){
    FILE* in = fopen(name, "r");
    
    int vN, fN;
    fscanf(in, "%d", &vN);
    fscanf(in, "%d", &fN);
    Mesh* mesh = new Mesh(vN, fN);
    
    int i;
    for(i=0; i<vN; i++){
      float x, y, z;
      fscanf(in, "%f %f %f", &x, &y, &z);
      mesh->setVertex(i, x, y, z);
    }
    
    for(i=0; i<fN; i++){
      int dummy, f0, f1, f2;
      fscanf(in, "%d %d %d %d", &dummy, &f0, &f1, &f2);
      mesh->setFace(i, f0, f1, f2);
    }
    
    mesh->computeNormal();
    
    fclose(in);
    
    return mesh;
  }
  
  static Mesh* readPlynFile( char* name ) {
    FILE* in = fopen(name, "r");
    
    int vN, fN;
    fscanf(in, "%d", &vN);
    fscanf(in, "%d", &fN);
    Mesh* mesh = new Mesh(vN, fN);
    
    int i;
    for(i=0; i<vN; i++){
      float x, y, z;
      fscanf(in, "%f %f %f", &x, &y, &z);
      mesh->setVertex(i, x, y, z);
    }
    
    for(i=0; i<fN; i++){
      int dummy, f0, f1, f2;
      fscanf(in, "%d %d %d %d", &dummy, &f0, &f1, &f2);
      mesh->setFace(i, f0, f1, f2);
    }
    
    for(i=0; i<vN; i++){
      float x, y, z;
      fscanf(in, "%f %f %f", &x, &y, &z);
      mesh->setNormal(i, x, y, z);
    }
    
    fclose(in);
    
    return mesh;
  }
  
  static int readPwnFile( char* name, float (*&point)[3], float (*&normal)[3]){
    FILE* in = fopen(name, "r");
    
    int N;
    fscanf(in, "%d", &N);
    
    point = new float[N][3];
    normal = new float[N][3];
    
    int i;
    for(i=0; i<N; i++){
      float* p = point[i];
      fscanf(in, "%f %f %f", &p[0], &p[1], &p[2]);
    }
    
    for(i=0; i<N; i++){
      float *n = normal[i];
      fscanf(in, "%f %f %f", &n[0], &n[1], &n[2]);
    }
    
    fclose(in);
    
    return N;
  }
  
  static int readBallFile( char* name, BasisFunction*** &lists, int* &Ns,
                          int &deg, bool& oriented ) {
    FILE* in = fopen(name, "r+b");
    
    int N[2];
    fread(N, sizeof(int), 2, in);
    deg = N[0];
    if(deg < 0){
      deg = -deg;
      oriented = false;
    }
    else
      oriented = true;
    int listN = N[1];
    
    lists = new BasisFunction**[2*listN];
    Ns = new int[2*listN];
    int i, j;
    float c[24];
    for(i=0; i<2*listN; i++){
      fread(N, sizeof(int), 1, in);
      int M = Ns[i] = N[0];
      if(M == 0)
        continue;
      BasisFunction** list = lists[i] = new BasisFunction*[M];
      for(j=0; j<M; j++){
        if(deg == 3){
          fread(c, sizeof(float), 24, in);
          
          Cubic* bf = new Cubic();
          list[j] = bf;
          
          bf->centerX = c[0];
          bf->centerY = c[1];
          bf->centerZ = c[2];
          
          bf->cXXX = c[3];
          bf->cYYY = c[4];
          bf->cZZZ = c[5];
          
          bf->cXXY = c[6];
          bf->cYYZ = c[7];
          bf->cZZX = c[8];
          
          bf->cXYY = c[9];
          bf->cYZZ = c[10];
          bf->cZXX = c[11];
          
          bf->cXYZ = c[12];
          
          bf->cXX = c[13];
          bf->cYY = c[14];
          bf->cZZ = c[15];
          
          bf->cXY = c[16];
          bf->cYZ = c[17];
          bf->cZX = c[18];
          
          bf->cX = c[19];
          bf->cY = c[20];
          bf->cZ = c[21];
          
          bf->c0 = c[22];
          
          bf->support = c[23];
          
          bf->leaf = (i%2 == 0);
        }
        else if(deg == 2){
          fread(c, sizeof(float), 14, in);
          
          Quadratic* bf = new Quadratic();
          list[j] = bf;
          
          bf->centerX = c[0];
          bf->centerY = c[1];
          bf->centerZ = c[2];
          
          bf->cXX = c[3];
          bf->cYY = c[4];
          bf->cZZ = c[5];
          
          bf->cXY = c[6];
          bf->cYZ = c[7];
          bf->cZX = c[8];
          
          bf->cX = c[9];
          bf->cY = c[10];
          bf->cZ = c[11];
          
          bf->c0 = c[12];
          
          bf->support = c[13];
          
          bf->leaf = (i%2 == 0);
        }
        else{
          fread(c, sizeof(float), 8, in);
          
          Linear* bf = new Linear();
          list[j] = bf;
          
          bf->centerX = c[0];
          bf->centerY = c[1];
          bf->centerZ = c[2];
          
          bf->cX = c[3];
          bf->cY = c[4];
          bf->cZ = c[5];
          
          bf->c0 = c[6];
          
          bf->support = c[7];
          
          bf->leaf = (i%2 == 0);
        }
      }
    }
    
    fclose(in);
    
    return listN;
  }
  
  static void writePpmFile( char* name, GLubyte *image, int w, int h ) {
    FILE* out = fopen(name, "w");
    fprintf(out, "P3\n");
    fprintf(out, "%d %d\n", w, h);
    fprintf(out, "255\n");
    GLubyte* im = image;
    for(int i=0; i<h; i++){
      for(int j=0; j<w; j++){
        fprintf(out, "%d ", (int)(*im));
        im++;
        fprintf(out, "%d ", (int)(*im));
        im++;
        fprintf(out, "%d ", (int)(*im));
        im++;
      }
      fprintf(out, "\n");
    }
    fclose(out);
  }
};

#endif
