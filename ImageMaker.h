#ifndef IMAGEMAKER
#define IMAGEMAKER

#include <cmath>
using namespace std;

#include <gl\gl.h>
#include <gl\glu.h>

#include "BasisFunction.h"
#include "BfAxisKdTree.h"
#include "Mesh.h"

#define EMPTY 1000000
#define RIDGE 10
#define VALLEY -10
#define SILHOUETTE -10000
#define NORMAL 0

#define FEW_PIX 4

class ImageMaker{

public:
  float DERIVATIVE;
  float DERIVATIVE_I;
  
  int width, height;
  
  bool multi;
  bool both;
  
  BasisFunction* root;
  BasisFunction* bft;
  int degree;
  
  int pointN;
  float (*point)[3];
  float (*normal)[3];
  float (*o_point)[3];
  float (*o_normal)[3];
  
  float **front;
  float **depth;
  float **radius;
  float (**normals)[3];
  float (**second)[6];
  float (**third)[10];
  float **value;
  float **totalW;
  bool spec;
  
  float view_point[3];
  float view_center[3];
  float view_up[3];
  float view_port[4];
  
  GLubyte *image;
  
  float splat;
  bool update;
  
  float sizeX, sizeY;
  float centerX, centerY;
  
  float rate;
  
  int lightN;
  float (*lightC)[3];
  float ambient[3];
  float (*lightD)[3];
  float materialD[3];
  float materialS[3];
  
  bool point_render;
  bool mesh_render;
  
  Mesh* mesh;
  
  float Rt[9], R[9];
  float is, is2, is3, s;
  float Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz;
  float tx, ty, tz;
  
  BasisFunction **stack;
  int top;
  int visitID;
  int leaf_counter;
  int total_leaf;
  
  int bfListN;
  BasisFunction*** bf_lists;
  int* bfNs;
  
  ImageMaker(int w, int h){
    both = false;
    DERIVATIVE = 1.5f;
    DERIVATIVE_I = 2.0f/3.0f;
    
    point_render = false;
    mesh_render = false;
    
    splat = 3.0f;
    width = w;
    height = h;
    
    rate = 1.0f;
    
    view_point[0] = 0;
    view_point[1] = 0;
    view_point[2] = 100;
    
    view_center[0] = 0;
    view_center[1] = 0;
    view_center[2] = 0;
    
    view_up[0] = 0;
    view_up[1] = 1;
    view_up[2] = 0;
    
    view_port[0] = -10;
    view_port[1] = -10;
    view_port[2] = 10;
    view_port[3] = 10;
    
    allocateBuffers();
    
    materialD[0] = 0.5f;//1;//255.0f/255.0f;//0.5f;//1.0; //1.0f;
    materialD[1] = 0.5f;//1;//125.0f/255.0f;//0.5f;////0.5; //1.0f;
    materialD[2] = 1.0f;//1;//180.0f/255.0f;//1.0f;//0.25f;
    
    materialS[0] = 0.5f;
    materialS[1] = 0.5f;
    materialS[2] = 0.5f;
    /*
      lightN = 3;
    
      lightD = new float[lightN][3];
    
      lightD[0][0] = 0;
      lightD[0][1] = -(float)(1.0/sqrt(2));
      lightD[0][2] = (float)(1.0/sqrt(2));
    
      lightD[1][0] = -(float)(1.0/sqrt(3));
      lightD[1][1] = -(float)(1.0/sqrt(3));
      lightD[1][2] = (float)(1.0/sqrt(3));
    
      lightD[2][0] = -(float)(1.0/sqrt(2));
      lightD[2][1] = 0;
      lightD[2][2] = (float)(1.0/sqrt(2));
    
      lightC = new float[lightN][3];
    
      lightC[0][0] = 1.0f;
      lightC[0][1] = 0;
      lightC[0][2] = 0;
    
      lightC[1][0] = 0;
      lightC[1][1] = 1.0f;
      lightC[1][2] = 0;
    
      lightC[2][0] = 0;
      lightC[2][1] = 0;
      lightC[2][2] = 1.0f;
    */
    
    lightN = 1;
    
    lightD = new float[lightN][3];
    
    lightD[0][0] = 0;
    lightD[0][1] = -(float)(1.0/std::sqrt(2.0f));
    lightD[0][2] = (float)(1.0/std::sqrt(2.0f));
    
    /*lightD[1][0] = -(float)(1.0/sqrt(3));
      lightD[1][1] = -(float)(1.0/sqrt(3));
      lightD[1][2] = (float)(1.0/sqrt(3));
    */
    lightC = new float[lightN][3];
    
    lightC[0][0] = 0.75f;
    lightC[0][1] = 0.75f;
    lightC[0][2] = 0.75f;
    /*
      lightC[1][0] = 0.5f;
      lightC[1][1] = 0.5f;
      lightC[1][2] = 0.5f;
    */
    ambient[0] = 0.2f;
    ambient[1] = 0.2f;
    ambient[2] = 0.2f;
    
    spec = true;
    
    pointN = 0;
    mesh = NULL;
    
    visitID = 0;
    
    root = NULL;
    bft = NULL;
  }
  
  ~ImageMaker(){
    clearBuffers();
    clearData();
  }
  
  void allocateBuffers(){
    image = new GLubyte[height*width*3];
    depth = new float*[height];
    front = new float*[height];
    radius = new float*[height];
    value = new float*[height];
    totalW = new float*[height];
    normals = (float(**)[3])new float*[height];
    second = (float(**)[6])new float*[height];
    third = (float(**)[10])new float*[height];
    for(int i=0; i<height; i++){
      depth[i] = new float[width];
      front[i] = new float[width];
      value[i] = new float[width];
      radius[i] = new float[width];
      totalW[i] = new float[width];
      normals[i] = new float[width][3];
      second[i] = new float[width][6];
      third[i] = new float[width][10];
    }
  }
  
  void clearBuffers(){
    for(int i=0; i<height; i++){
      delete[] depth[i];
      delete[] front[i];
      delete[] radius[i];
      delete[] totalW[i];
      delete[] value[i];
      delete[] normals[i];
      delete[] second[i];
      delete[] third[i];
    }
    delete[] image;
    delete[] depth;
    delete[] front;
    delete[] radius;
    delete[] totalW;
    delete[] value;
    delete[] normals;
    delete[] second;
    delete[] third;
  }
  
  void reseizeImageBuffer(int w, int h){
    clearBuffers();
    
    width = w;
    height = h;
    
    allocateBuffers();
  }
  
  void setBfLists(int deg, int listN, BasisFunction*** bf_lists, int* bfNs){
    clearData();
    
    total_leaf = 0;
    
    this->bfListN = listN;
    this->bf_lists = bf_lists;
    this->bfNs = bfNs;
    
    degree = deg;
    
    point_render = false;
    mesh_render = false;
    
    int bfN_leaf = bfNs[0];
    total_leaf += bfN_leaf;
    int bfN_inter = bfNs[1];
    int bfN_total = bfN_leaf + bfN_inter;
    BasisFunction** bfs_child = new BasisFunction*[bfN_total];
    int i, j, k;
    for(i=0; i<bfN_leaf; i++)
      bfs_child[i] = bf_lists[0][i];
    for(i=0; i<bfN_inter; i++)
      bfs_child[i+bfN_leaf] = bf_lists[1][i];
    
    root = allocateBf();
    root->setChild(bfs_child, bfN_total);
    root->centerX = root->centerY = root->centerZ = 0;
    root->support = 1000000;
    
    stack = new BasisFunction*[100*listN + bfN_total];
    
    for(i=1; i<listN; i++){
      bfN_leaf = bfNs[2*i];
      total_leaf += bfN_leaf;
      bfN_inter = bfNs[2*i+1];
      bfN_total = bfN_leaf + bfN_inter;
      
      if(bfN_total == 0)
        continue;
      
      BasisFunction** bfs_current = new BasisFunction*[bfN_total];
      for(j=0; j<bfN_leaf; j++)
        bfs_current[j] = bf_lists[2*i][j];
      for(j=0; j<bfN_inter; j++)
        bfs_current[j+bfN_leaf] = bf_lists[2*i+1][j];
      
      BfAxisKdTree* tree = new BfAxisKdTree(bfs_current, bfN_total);
      
      int bfN_prev = bfNs[2*i-1];
      BasisFunction** bfs_prev = bf_lists[2*i-1];
      for(j=0; j<bfN_prev; j++){
        BasisFunction* bf = bfs_prev[j];
        //bf->support *= DERIVATIVE;
        int *list, listNB;
        tree->collectBfIndexInBf(list, listNB, bf);
        //bf->support *= DERIVATIVE_I;
        if(listNB == 0)
          continue;
        
        bfs_child = new BasisFunction*[listNB];
        for(k=0; k<listNB; k++)
          bfs_child[k] = bfs_current[list[k]];
        bf->setChild(bfs_child, listNB);
      }
      delete tree;
      delete[] bfs_current;
    }
    
    both = false;
    visitID = 0;
    bft = allocateBf();
  }
  
  void setBoth(bool flag){
    both = flag;
    if(!both){
      DERIVATIVE = 1.5f;
      DERIVATIVE_I = 2.0f/3.0f;
    }
    else{
      DERIVATIVE = 2.0f;
      DERIVATIVE_I = 0.5;
    }
  }
  
  BasisFunction* allocateBf(){
    if(degree == 3)
      return new Cubic;
    else if(degree == 2)
      return new Quadratic;
    else
      return new Linear;
  }
  
  void setPoints(int N, float (*p)[3], float (*n)[3]){
    clearData();
    
    point_render = true;
    mesh_render = false;
    
    pointN = N;
    o_point = p;
    o_normal = n;
    point = new float[N][3];
    normal = new float[N][3];
    
    both = true;
  }
  
  void rescalePoints(float scale){
    int i, j;
    float min[3], max[3];
    
    for(i=0; i<3; i++){
      min[i] = max[i] = o_point[0][i];
      for(j=1; j<pointN; j++){
        float p = o_point[j][i];
        if(p < min[i])
          min[i] = p;
        else if(p > max[i])
          max[i] = p;
      }
    }
    
    float sx = max[0] - min[0];
    float sy = max[1] - min[1];
    float sz = max[2] - min[2];
    float s = scale/(float)sqrt(sx*sx + sy*sy + sz*sz);
    
    for(i=0; i<3; i++){
      float m = 0.5f*(min[i] + max[i]);
      for(j=0; j<pointN; j++){
        o_point[j][i] = s*(o_point[j][i] - m);
      }
    }
  }
  
  void setMesh(Mesh* m){
    clearData();
    
    point_render = false;
    mesh_render = true;
    
    mesh = m;
    
    both = true;
  }
  
  void clearData(){
    if(pointN != 0){
      delete[] point;
      delete[] normal;
      delete[] o_point;
      delete[] o_normal;
      pointN = 0;
    }
    
    if(mesh != NULL){
      delete mesh;
      mesh = NULL;
    }
    
    if(root != NULL){
      for(int i=0; i<2*bfListN; i++){
        if(bfNs[i] == 0)
          continue;
        for(int j=0; j<bfNs[i]; j++)
          delete bf_lists[i][j];
        delete[] bf_lists[i];
      }
      delete[] bf_lists;
      delete[] bfNs;
      
      delete root;
      delete bft;
      root = NULL;
    }
  }
  
  void createImageRef(int N){
    GLubyte *im = image;
    for(int i=0; i<height; i++){
      for(int j=0; j<width; j++){
        float* n = normals[i][width-j-1];
        if(n[0] == 0 && n[1] == 0 && n[2] == 0){
          *im++ = (GLubyte)255;
          *im++ = (GLubyte)255;
          *im++ = (GLubyte)255;
          continue;
        }
        
        float R = ambient[0];
        float G = ambient[1];
        float B = ambient[2];
        
        for(int k=0; k<lightN; k++){
          float* ld = lightD[k];
          //float ld[3] = {0,-1.0f/sqrt(2),1/sqrt(2)};
          float dot = ld[0]*n[0] + ld[1]*n[1] + ld[2]*n[2];
          if(dot < 0)
            continue;
          float m;
          if((int)(N*dot)%2 == 0)
            m = 1;
          else
            m = 0;
          float* lc = lightC[k];
          R += m*lc[0]*dot;
          G += m*lc[1]*dot;
          B += m*lc[2]*dot;
        }
        if(spec){
          for(int k=0; k<lightN; k++){
            float* ld = lightD[k];
	    //float ld[3] = {0,-1.0f/sqrt(2),1/sqrt(2)};
            
            float d2 = 2.0f*(ld[0]*n[0] + ld[1]*n[1] + ld[2]*n[2]);
            float dot = d2*n[2] - ld[2];
            if(dot > 1.0f)
              dot = 1.0f;
            else if(dot < 0)
              continue;
            
            dot = (float)pow(dot, 10);
            float* lc = lightC[k];
            R += materialS[0]*lc[0]*dot;
            G += materialS[1]*lc[1]*dot;
            B += materialS[2]*lc[2]*dot;
          }
        }
        
        if(R > 1.0f)
          R = 1.0f;
        else if(R < 0.0f)
          R = 0.0f;
        *im++ = (GLubyte)(255*R);
        
        if(G > 1.0f)
          G = 1.0f;
        else if(G < 0.0f)
          G = 0.0f;
        *im++ = (GLubyte)(255*G);
        
        if(B > 1.0f)
          B = 1.0f;
        else if(B < 0.0f)
          B = 0.0f;
        *im++ = (GLubyte)(255*B);
      }
    }
  }
  
  void createImage(){
    GLubyte *im = image;
    for(int i=0; i<height; i++){
      for(int j=0; j<width; j++){
        float* n = normals[i][width-j-1];
        if(n[0] == 0 && n[1] == 0 && n[2] == 0){
          *im++ = (GLubyte)255;
          *im++ = (GLubyte)255;
          *im++ = (GLubyte)255;
          continue;
        }
        
        float R = ambient[0];
        float G = ambient[1];
        float B = ambient[2];
        
        for(int k=0; k<lightN; k++){
          float* ld = lightD[k];
          float dot = ld[0]*n[0] + ld[1]*n[1] + ld[2]*n[2];
          if(dot < 0)
            continue;
          
          float* lc = lightC[k];
          //if(j < width/2){
          R += materialD[0]*lc[0]*dot;
          G += materialD[1]*lc[1]*dot;
          B += materialD[2]*lc[2]*dot;
          /*}
	    else{
	    float m;
	    if((int)(20*dot)%2 == 0)
            m = 1;
	    else
            m = 0;
	    R += m*lc[0]*dot;
	    G += m*lc[1]*dot;
	    B += m*lc[2]*dot;
	    }*/
        }
        if(spec){
          for(int k=0; k<lightN; k++){
            float* ld = lightD[k];
            //float ld[3] = {-1/sqrt(3), -1/sqrt(3), 1/sqrt(3)};
            float d2 = 2.0f*(ld[0]*n[0] + ld[1]*n[1] + ld[2]*n[2]);
            float dot = d2*n[2] - ld[2];
            if(dot > 1.0f)
              dot = 1.0f;
            else if(dot < 0)
              continue;
            
            dot = (float)pow(dot, 10);
            float* lc = lightC[k];
            R += materialS[0]*lc[0]*dot;
            G += materialS[1]*lc[1]*dot;
            B += materialS[2]*lc[2]*dot;
          }
        }
        
        if(R > 1.0f)
          R = 1.0f;
        else if(R < 0.0f)
          R = 0.0f;
        *im++ = (GLubyte)(255*R);
        
        if(G > 1.0f)
          G = 1.0f;
        else if(G < 0.0f)
          G = 0.0f;
        *im++ = (GLubyte)(255*G);
        
        if(B > 1.0f)
          B = 1.0f;
        else if(B < 0.0f)
          B = 0.0f;
        *im++ = (GLubyte)(255*B);
      }
    }
  }
  
  void createValueImage(){
    int i;
    GLubyte *im = image;
    
    float mid, vari;
    float S = 0;
    
    int count = 0;
    for(i=0; i<height; i++){
      for(int j=0; j<width; j++){
        if(value[i][j] != EMPTY){
          S += value[i][j];
          count++;
        }
      }
    }
    mid = S/count;
    
    S = 0;
    for(i=0; i<height; i++){
      for(int j=0; j<width; j++){
        if(value[i][j] != EMPTY){
          float v = value[i][j];
          S += (v-mid)*(v-mid);
        }
      }
    }
    vari = (float)sqrt(S/count);
    
    for(i=0; i<height; i++){
      for(int j=0; j<width; j++){
        float v = value[i][width-j-1];
        if(v == EMPTY){
          *im++ = (GLubyte)255;
          *im++ = (GLubyte)255;
          *im++ = (GLubyte)255;
          continue;
        }
        float c[3];
        convertColor(c, mid, vari, v);
        
        float mR = c[0];
        float mG = c[1];
        float mB = c[2];
        
        float* n = normals[i][width-j-1];
        
        float R = 2*ambient[0];
        float G = 2*ambient[1];
        float B = 2*ambient[2];
        
        for(int k=0; k<lightN; k++){
          float* ld = lightD[k];
          float dot = ld[0]*n[0] + ld[1]*n[1] + ld[2]*n[2];
          if(dot < 0)
            continue;
          
          float *lc = lightC[k];
          R += mR*lc[0]*dot;
          G += mG*lc[1]*dot;
          B += mB*lc[2]*dot;
        }
        if(spec){
          for(int k=0; k<lightN; k++){
            float* ld = lightD[k];
            
            float d2 = 2.0f*(ld[0]*n[0] + ld[1]*n[1] + ld[2]*n[2]);
            float dot = d2*n[2] - ld[2];
            if(dot > 1.0f)
              dot = 1.0f;
            else if(dot < 0)
              continue;
            
            dot = (float)pow(dot, 10);
            
            float *lc = lightC[k];
            R += materialS[0]*lc[0]*dot;
            G += materialS[1]*lc[1]*dot;
            B += materialS[2]*lc[2]*dot;
          }
        }
        
        if(R > 1.0f)
          R = 1.0f;
        else if(R < 0.0f)
          R = 0.0f;
        *im++ = (GLubyte)(255*R);
        
        if(G > 1.0f)
          G = 1.0f;
        else if(G < 0.0f)
          G = 0.0f;
        *im++ = (GLubyte)(255*G);
        
        if(B > 1.0f)
          B = 1.0f;
        else if(B < 0.0f)
          B = 0.0f;
        *im++ = (GLubyte)(255*B);
      }
    }
  }
  
  void createRvImage(bool shading){
    int i;
    GLubyte *im = image;
    
    for(i=0; i<height; i++){
      for(int j=0; j<width; j++){
        float v = value[i][width-j-1];
        if(v == EMPTY){
          *im++ = (GLubyte)255;
          *im++ = (GLubyte)255;
          *im++ = (GLubyte)255;
          continue;
        }
        
        float R, G, B;
        
        if(v == RIDGE){
          if(shading)
            R = 1;
          else
            R = 0;
          G = B = 0;
        }
        else if(v == VALLEY){
          if(shading)
            B = 1;
          else
            B = 0;
          G = R = 0;
        }
        else if(v == SILHOUETTE){
          R = G = B = 0;
        }
        else if(shading){
          float* n = normals[i][width-j-1];
          
          R = 2*ambient[0];
          G = 2*ambient[1];
          B = 2*ambient[2];
          
          float mR = 1;
          float mG = 1;
          float mB = 1;
          
          for(int k=0; k<lightN; k++){
            float* ld = lightD[k];
            float dot = ld[0]*n[0] + ld[1]*n[1] + ld[2]*n[2];
            if(dot < 0)
              continue;
            
            float *lc = lightC[k];
            R += mR*lc[0]*dot;
            G += mG*lc[1]*dot;
            B += mB*lc[2]*dot;
          }
          if(spec){
            for(int k=0; k<1; k++){
              float* ld = lightD[0];
              
              float d2 = 2.0f*(ld[0]*n[0] + ld[1]*n[1] + ld[2]*n[2]);
              float dot = d2*n[2] - ld[2];
              if(dot > 1.0f)
                dot = 1.0f;
              else if(dot < 0)
                continue;
              
              dot = (float)pow(dot, 10);
              
              float *lc = lightC[k];
              R += materialS[0]*lc[0]*dot;
              G += materialS[1]*lc[1]*dot;
              B += materialS[2]*lc[2]*dot;
            }
          }
        }
        else
          R = G = B = 1.0f;
        
        if(R > 1.0f)
          R = 1.0f;
        else if(R < 0.0f)
          R = 0.0f;
        *im++ = (GLubyte)(255*R);
        
        if(G > 1.0f)
          G = 1.0f;
        else if(G < 0.0f)
          G = 0.0f;
        *im++ = (GLubyte)(255*G);
        
        if(B > 1.0f)
          B = 1.0f;
        else if(B < 0.0f)
          B = 0.0f;
        *im++ = (GLubyte)(255*B);
      }
    }
  }
  
  void createBinaryImage(){
    int i;
    GLubyte *im = image;
    
    for(i=0; i<height; i++){
      for(int j=0; j<width; j++){
        float v = value[i][width-j-1];
        if(v == EMPTY){
          *im++ = (GLubyte)255;
          *im++ = (GLubyte)255;
          *im++ = (GLubyte)255;
          continue;
        }
        
        float mR, mG, mB;
        if(v == 0){
          mR = 1.0f;
          mG = 0.5f;
          mB = 1.0f;
        }
        else{
          mR = 0.0f;
          mG = 0.0f;
          mB = 0.5f;
        }
        
        float* n = normals[i][width-j-1];
        
        float R = 2*ambient[0];
        float G = 2*ambient[1];
        float B = 2*ambient[2];
        
        for(int k=0; k<lightN; k++){
          float* ld = lightD[k];
          float dot = ld[0]*n[0] + ld[1]*n[1] + ld[2]*n[2];
          if(dot < 0)
            continue;
          
          float* lc = lightC[k];
          R += mR*lc[0]*dot;
          G += mG*lc[1]*dot;
          B += mB*lc[2]*dot;
        }
        if(spec){
          for(int k=0; k<lightN; k++){
            float* ld = lightD[k];
            
            float d2 = 2.0f*(ld[0]*n[0] + ld[1]*n[1] + ld[2]*n[2]);
            float dot = d2*n[2] - ld[2];
            if(dot > 1.0f)
              dot = 1.0f;
            else if(dot < 0)
              continue;
            
            dot = (float)pow(dot, 10);
            
            float* lc = lightC[k];
            R += 2*materialS[0]*lc[0]*dot;
            G += 2*materialS[1]*lc[1]*dot;
            B += 2*materialS[2]*lc[2]*dot;
          }
        }
        
        if(R > 1.0f)
          R = 1.0f;
        else if(R < 0.0f)
          R = 0.0f;
        *im++ = (GLubyte)(255*R);
        
        if(G > 1.0f)
          G = 1.0f;
        else if(G < 0.0f)
          G = 0.0f;
        *im++ = (GLubyte)(255*G);
        
        if(B > 1.0f)
          B = 1.0f;
        else if(B < 0.0f)
          B = 0.0f;
        *im++ = (GLubyte)(255*B);
      }
    }
  }
  
  void createNormalsSupp(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        normals[i][j][0] = normals[i][j][1] = normals[i][j][2] = 0;
        value[i][j] = EMPTY;
      }
    }
    
    leaf_counter = 0;
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    rate = 1.0f;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      
      float r = bft->support*rate;
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX){
        leaf_counter++;
        
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            if((k-cu)*(k-cu) > dj)
              continue;
            
            float f = front[j][k];
            if(f < cd)
              continue;
            
            float start[3] = {(float)k, (float)j, 0};
            
            float q[3];
            if(!bft->rayIntersectSuppZ(q, start))
              continue;
            
            float qd = q[2];
            
            if(qd < f){
              front[j][k] = qd;
              float n1[3];
              bft->normalSupp(n1, q[0], q[1], q[2]);
              
              float* n = normals[j][k];
              n[0] = n1[0];
              n[1] = n1[1];
              n[2] = n1[2];
              
              value[j][k] = -r; 
            }
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfs[i];
          }
        }
      }
    }
  }
  
  void createNormalsSuppB(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        normals[i][j][0] = normals[i][j][1] = normals[i][j][2] = 0;
        value[i][j] = EMPTY;
      }
    }
    
    leaf_counter = 0;
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      
      float r = bft->support*rate;
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX){
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            if((k-cu)*(k-cu) > dj)
              continue;
            
            float f = front[j][k];
            if(f < cd)
              continue;
            
            float start[3] = {(float)k, (float)j, 0};
            
            float q[3];
            if(!bft->rayIntersectSuppZ(q, start))
              continue;
            
            float qd = q[2];
            
            if(qd < f){
              front[j][k] = qd;
              float n1[3];
              bft->normalSupp(n1, q[0], q[1], q[2]);
              
              float* n = normals[j][k];
              n[0] = n1[0];
              n[1] = n1[1];
              n[2] = n1[2];
              
              if(bft->leaf)
                value[j][k] = 0;
              else
                value[j][k] = 1;
            }
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfs[i];
          }
        }
      }
    }
  }
  
  void createNormalsFirstInt(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 1000000;
	value[i][j] = EMPTY;
        normals[i][j][0] = normals[i][j][1] = normals[i][j][2] = 0;
      }
    }

    leaf_counter = 0;
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      
      float r = bft->support*rate;
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX){
        transPoly(bft, bf);
        
        if(!both && bft->normalZAtC() > 0)
          continue;
        
        leaf_counter++;
        
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            if((k-cu)*(k-cu) > dj)
              continue;
            
            float f = front[j][k];
            if(f < cd)
              continue;
            
            float qd;
            if(!bft->rayIntersectLAZ(qd, (float)k, (float)j))
              continue;
            
            if(qd < f){
              front[j][k] = qd;
              
              float n1[3];
              bft->normalLA(n1, (float)k, (float)j, qd);
              if(n1[2] > 0){
                if(!both)
                  continue;
                n1[0] = -n1[0];
                n1[1] = -n1[1];
                n1[2] = -n1[2];
              }
              
              float* n = normals[j][k];
              n[0] = -n1[0];
              n[1] = -n1[1];
              n[2] = -n1[2];

	      value[j][k] = -r;
            }
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfs[i];
          }
        }
      }
    }
  }
  
  void createNormalsP(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        normals[i][j][0] = normals[i][j][1] = normals[i][j][2] = 0;
      }
    }
    
    float r = splat;
    float r2 = r*r;
    for(i=0; i<pointN; i++){
      float* p = point[i];
      float* n = normal[i];
      float cu = p[0];
      float cv = p[1];
      float cd = p[2];
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
	continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      for(int j = sj; j < ej; j++){
        for(int k = sk; k < ek; k++){
          if((k-cu)*(k-cu) + (j-cv)*(j-cv) > r2)
            continue;
          
          if(front[j][k] > cd){
            float *m = normals[j][k];
            if(n[2] > 0){
              if(both){
                n[0] = -n[0];
                n[1] = -n[1];
                n[2] = -n[2];
              }
              else
                continue;
            }
            
            m[0] = -n[0];
            m[1] = -n[1];
            m[2] = -n[2];
            
            front[j][k] = cd;
          }
        }
      }
    }
  }
  
  void createNormalsFullV(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        normals[i][j][0] = normals[i][j][1] = normals[i][j][2] = 0;
        value[i][j] = 0;
        totalW[i][j] = 0;
      }
    }
    
    computeFrontDepth();
    
    leaf_counter = 0;
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      bft->support *= DERIVATIVE;
      
      float r = bft->support;
      
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ;// - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX*DERIVATIVE){
        transPoly(bft, bf);
        
        if(!both && bft->normalZAtC() > 0)
          continue;
        
        leaf_counter++;
        r = bft->support;
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            float f = front[j][k];
            if((k-cu)*(k-cu) + (f-cd)*(f-cd) > dj)
              continue;
            
            float qd;
            if(!bft->rayIntersectLAZ(qd, (float)k, (float)j))
              continue;
            
            float n1[3];
            bft->normalLA(n1, (float)k, (float)j, qd);
            if(n1[2] > 0){
              if(!both)
                continue;
              n1[0] = -n1[0];
              n1[1] = -n1[1];
              n1[2] = -n1[2];
            }
            
            float w = -(float)bft->weight((float)k, (float)j, qd);
            float *n = normals[j][k];
            n[0] += w*n1[0];
            n[1] += w*n1[1];
            n[2] += w*n1[2];
            
            totalW[j][k] += w;
            value[j][k] += -w*r;
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfi;
          }
        }
      }
    }
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        float *n = normals[i][j];
        float len = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        if(len != 0){
          len = (float)(1.0/sqrt(len));
          n[0] *= len;
          n[1] *= len;
          n[2] *= len;
          
          value[i][j] /= totalW[i][j];
        }
        else{
          n[0] = n[1] = n[2] = 0;
          value[i][j] = EMPTY;
        }
      }
    }
  }
  
  void createNormalM(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        normals[i][j][0] = normals[i][j][1] = normals[i][j][2] = 0;
      }
    }
    
    int faceN = mesh->faceN;
    
    for(i=0; i<faceN; i++){
      float bx[2], by[2];
      mesh->getFaceBoundX(bx, i);
      mesh->getFaceBoundY(by, i);
      
      int sk = (int)bx[0];
      if(sk < 0)
        sk = 0;
      int ek = (int)bx[1]+1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)by[0];
      if(sj < 0)
        sj = 0;
      int ej = (int)by[1]+1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
	continue;
      
      for(int j = sj; j < ej; j++){
        for(int k = sk; k < ek; k++){
          float qd, n1[3];
          if(!mesh->rayIntersectFaceZ(qd, n1, i, (float)k, (float)j))
            continue;
          
          if(n1[2] > 0){
            if(!both)
              continue;
            n1[0] = -n1[0];
            n1[1] = -n1[1];
            n1[2] = -n1[2];
          }
          
          if(qd  < front[j][k]){
            front[j][k] = qd;
            float *n = normals[j][k];
            n[0] = -n1[0];
            n[1] = -n1[1];
            n[2] = -n1[2];
          }
        }
      }
    }
  }
  
  void createNormalsFull(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        normals[i][j][0] = normals[i][j][1] = normals[i][j][2] = 0;
      }
    }
    
    computeFrontDepth();
    
    leaf_counter = 0;
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      bft->support *= DERIVATIVE;
      float r = bft->support;
      
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ;// - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX*DERIVATIVE){
        transPoly(bft, bf);
        if(!both && bft->normalZAtC() > 0)
          continue;
        
        leaf_counter++;
        r = bft->support;
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            float f = front[j][k];
            if((k-cu)*(k-cu) + (f-cd)*(f-cd) > dj)
              continue;
            
            float qd;
            if(!bft->rayIntersectLAZ(qd, (float)k, (float)j))
              continue;
            
            float n1[3];
            bft->normalLA(n1, (float)k, (float)j, qd);
            if(n1[2] > 0){
              if(!both)
                continue;
              n1[0] = -n1[0];
              n1[1] = -n1[1];
              n1[2] = -n1[2];
            }
            
            float w = -(float)bft->weight((float)k, (float)j, qd);
            float *n = normals[j][k];
            n[0] += w*n1[0];
            n[1] += w*n1[1];
            n[2] += w*n1[2];
            
            totalW[j][k] += w;
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfi;
          }
        }
      }
    }
    
    normalizeNormals();
  }
  
  void createSC1(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        totalW[i][j] = radius[i][j] = depth[i][j] =
	  normals[i][j][0] = normals[i][j][1] = normals[i][j][2] =
	  second[i][j][0] = second[i][j][1] = second[i][j][2] =
	  second[i][j][3] = second[i][j][4] = second[i][j][5] = 0;
        value[i][j] = EMPTY;
      }
    }
    
    computeFrontDepth();
    
    leaf_counter = 0;
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      bft->support *= DERIVATIVE;
      
      float r = bft->support;
      
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ;// - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX*DERIVATIVE){
        transPoly(bft, bf);
        
        if(!both && bft->normalZAtC() > 0)
          continue;
        
        leaf_counter++;
        
        r = bft->support;
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            float f = front[j][k];
            if((k-cu)*(k-cu) + (f-cd)*(f-cd) > dj)
              continue;
            
            float qd;
            if(!bft->rayIntersectLAZ(qd, (float)k, (float)j))
              continue;
            
            float n1[3], gg1[6];
            bft->gradient12LA(n1, gg1, (float)k, (float)j, qd);
            if(n1[2] > 0){
              if(!both)
                continue;
              n1[0] = -n1[0];
              n1[1] = -n1[1];
              n1[2] = -n1[2];
              gg1[0] = -gg1[0];
              gg1[1] = -gg1[1];
              gg1[2] = -gg1[2];
              gg1[3] = -gg1[3];
              gg1[4] = -gg1[4];
              gg1[5] = -gg1[5];
            }
            float w = (float)bft->weight((float)k, (float)j, qd);
            float *n = normals[j][k];
            n[0] -= w*n1[0];
            n[1] -= w*n1[1];
            n[2] -= w*n1[2];
            float *gg = second[j][k];
            gg[0] -= w*gg1[0];
            gg[1] -= w*gg1[1];
            gg[2] -= w*gg1[2];
            gg[3] -= w*gg1[3];
            gg[4] -= w*gg1[4];
            gg[5] -= w*gg1[5];
            totalW[j][k] += w;
            depth[j][k] += w*qd;
            radius[j][k] += w*bft->support;
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfi;
          }
        }
      }
    }
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        float W = totalW[i][j];
        if(W == 0){
          depth[i][j] = 1000000;
          continue;
        }
        
        float *g = normals[i][j];
        float *gg = second[i][j];
        W = 1.0f/W;
        g[0] *= W;
        g[1] *= W;
        g[2] *= W;
        gg[0] *= W;
        gg[1] *= W;
        gg[2] *= W;
        gg[3] *= W;
        gg[4] *= W;
        gg[5] *= W;
        
        depth[i][j] *= W;
        radius[i][j] *= W;
        
	/* temporary derivatives*/
        float fxfx = g[0]*g[0];
        float fxfy = g[0]*g[1];
        float fxfz = g[0]*g[2];
        float fyfy = g[1]*g[1];
        float fyfz = g[1]*g[2];
        float fzfz = g[2]*g[2];
        
        float fxx = gg[0];
        float fxy = gg[1];
        float fxz = gg[2];
        float fyy = gg[3];
        float fyz = gg[4];
        float fzz = gg[5];
        
        float g2 = g[0]*g[0]+g[1]*g[1]+g[2]*g[2];
        float g1 = (float)sqrt(g2);
        if(g1 == 0)
          continue;
        float g3 = g2*g1;
        float g4 = g2*g2;
        
        g[0] /= g1;
        g[1] /= g1;
        g[2] /= g1;
        
        /* mean and gaussian curvatures */
        float H = (fxx*(fyfy+fzfz) + fyy*(fxfx+fzfz) + fzz*(fxfx+fyfy)
		   - 2*(fxy*fxfy+fxz*fxfz+fyz*fyfz))/2;
        H /= g3;
        
        float K = fxfx*(fyy*fzz-fyz*fyz)
          + fyfy*(fxx*fzz-fxz*fxz)
	  + fzfz*(fxx*fyy-fxy*fxy)
	  + 2*(fxfy*(fxz*fyz-fxy*fzz)
               + fxfz*(fxy*fyz-fxz*fyy)
	       + fyfz*(fxy*fxz-fxx*fyz));
        K /= g4;
        
        /* principal curvatures */
        float discr = (float)sqrt(fabs(H*H-K));
        
        float Kmax = H + discr;
        float Kmin = H - discr;
        
        /* matrix entries */
        float m11 = ((-1 + fxfx/g2)*fxx)/g1 + (fxfy*fxy)/g3 + (fxfz*fxz)/g3;
        float m12 = ((-1 + fxfx/g2)*fxy)/g1 + (fxfy*fyy)/g3 + (fxfz*fyz)/g3;
        float m13 = ((-1 + fxfx/g2)*fxz)/g1 + (fxfy*fyz)/g3 + (fxfz*fzz)/g3;
        float m21 = (fxfy*fxx)/g3 + ((-1 + fyfy/g2)*fxy)/g1 + (fyfz*fxz)/g3;
        float m22 = (fxfy*fxy)/g3 + ((-1 + fyfy/g2)*fyy)/g1 + (fyfz*fyz)/g3;
        float m23 = (fxfy*fxz)/g3 + ((-1 + fyfy/g2)*fyz)/g1 + (fyfz*fzz)/g3;
        float m31 = (fxfz*fxx)/g3 + (fyfz*fxy)/g3 + ((-1 + fzfz/g2)*fxz)/g1;
        float m32 = (fxfz*fxy)/g3 + (fyfz*fyy)/g3 + ((-1 + fzfz/g2)*fyz)/g1;
        float m33 = (fxfz*fxz)/g3 + (fyfz*fyz)/g3 + ((-1 + fzfz/g2)*fzz)/g1;
        
        /* solve for eigenvectors */
        float tmp1 = m11+Kmax;
        float tmp2 = m22+Kmax;
        float tmp3 = m33+Kmax;
        
        float ux[3], uy[3], uz[3], len[3];
        ux[0] = m12*m23-m13*tmp2;
        uy[0] = m13*m21-m23*tmp1;
        uz[0] = tmp1*tmp2-m12*m21;
        len[0] = (float)sqrt(ux[0]*ux[0]+uy[0]*uy[0]+uz[0]*uz[0]);
        
        ux[1] = m12*tmp3-m13*m32;
        uy[1] = m13*m31-tmp1*tmp3;
        uz[1] = tmp1*m32-m12*m31;
        len[1] = (float)sqrt(ux[1]*ux[1]+uy[1]*uy[1]+uz[1]*uz[1]);
        
        ux[2] = tmp2*tmp3-m23*m32;
        uy[2] = m23*m31-m21*tmp3;
        uz[2] = m21*m32-m31*tmp2;
        len[2] = (float)sqrt(ux[2]*ux[2]+uy[2]*uy[2]+uz[2]*uz[2]);
        
        int index = 0;
        double max = len[0];
        if ( len[1] > max ) {
          index  = 1;
          max = len[1];
        }
        if ( len[2] > max ) {
          index = 2;
          max = len[2];
        }
        
        float Tmax[3];
        Tmax[0] = ux[index]/len[index];
        Tmax[1] = uy[index]/len[index];
        Tmax[2] = uz[index]/len[index];
        
        /* second tangent is cross product of first tangent and normal */
        float Tmin[3];
        Tmin[0] = Tmax[1]*g[2]-Tmax[2]*g[1];
        Tmin[1] = Tmax[2]*g[0]-Tmax[0]*g[2];
        Tmin[2] = Tmax[0]*g[1]-Tmax[1]*g[0];
        
        value[i][j] = NORMAL;
        
        float c = (float)cos(atan2(Tmin[2], Tmax[2]));
        c *= c;
        gg[0] = Kmax*c + Kmin*(1.0f-c);
        
        gg[1] = g[0];
        gg[2] = g[1];
      }
    }
    
    for(i=1; i<height-1; i++){
      for(j=1; j<width-1; j++){
        if(value[i][j] == EMPTY && value[i][j+1] == EMPTY &&
           value[i+1][j] == EMPTY)
          continue;
        
        if(value[i][j] != EMPTY && value[i][j+1] != EMPTY &&
           value[i+1][j] != EMPTY){
          if(fabs(depth[i][j] - depth[i][j+1]) > radius[i][j] ||
             fabs(depth[i][j] - depth[i+1][j]) > radius[i][j]){
            value[i][j] = SILHOUETTE;
            continue;
          }
          
          if(normals[i][j][2] > 0.94f)
            continue;
          
          if(second[i][j][0]*second[i][j+1][0] >= 0 &&
             second[i][j][0]*second[i+1][j][0] >= 0)
            continue;
          
          float dx = (second[i][j][1] + second[i][j+1][1] + second[i+1][j][1])/3;
          float dy = (second[i][j][2] + second[i][j+1][2] + second[i+1][j][2])/3;
          if((second[i][j][0] - second[i][j+1][0])*dx +
             (second[i][j][0] - second[i+1][j][0])*dy > 0)
            value[i][j] = SILHOUETTE;
        }
        else
          value[i][j] = SILHOUETTE;
      }
    }
    
    return;
    
    for(i=1; i<height-1; i++){
      for(j=1; j<width-1; j++){
        if(value[i][j] == EMPTY)
          continue;
        
        float Kr = second[i][j][0];
        
        float tx = second[i][j][1];
        float ty = second[i][j][2];
        float t1 = ty/tx;
        float t2 = tx/ty;
        float k11, k12, k21, k22, t;
        float d11, d12, d21, d22;
        if(0 <= t1 && t1 < 1){
          if(tx > 0){
            k11 = second[i][j+1][0];
            k12 = second[i+1][j+1][0];
            k21 = second[i][j-1][0];
            k22 = second[i-1][j-1][0];
            
            d11 = depth[i][j+1];
            d12 = depth[i+1][j+1];
            d21 = depth[i][j-1];
            d22 = depth[i-1][j-1];
          }
          else{
            k21 = second[i][j+1][0];
            k22 = second[i+1][j+1][0];
            k11 = second[i][j-1][0];
            k12 = second[i-1][j-1][0];
            
            d21 = depth[i][j+1];
            d22 = depth[i+1][j+1];
            d11 = depth[i][j-1];
            d12 = depth[i-1][j-1];
          }
          t = t1;
        }
        else if(-1 <= t1 && t1 < 0){
          if(tx > 0){
            k11 = second[i][j+1][0];
            k12 = second[i-1][j+1][0];
            k21 = second[i][j-1][0];
            k22 = second[i+1][j-1][0];
            
            d11 = depth[i][j+1];
            d12 = depth[i-1][j+1];
            d21 = depth[i][j-1];
            d22 = depth[i+1][j-1];
          }
          else{
            k21 = second[i][j+1][0];
            k22 = second[i-1][j+1][0];
            k11 = second[i][j-1][0];
            k12 = second[i+1][j-1][0];
            
            d21 = depth[i][j+1];
            d22 = depth[i-1][j+1];
            d11 = depth[i][j-1];
            d12 = depth[i+1][j-1];
          }
          t = -t1;
        }
        else if(0 <= t2 && t2 <= 1){
          if(ty > 0){
            k11 = second[i+1][j][0];
            k12 = second[i+1][j+1][0];
            k21 = second[i-1][j][0];
            k22 = second[i-1][j-1][0];
            
            d11 = depth[i+1][j];
            d12 = depth[i+1][j+1];
            d21 = depth[i-1][j];
            d22 = depth[i-1][j-1];
          }
          else{
            k21 = second[i+1][j][0];
            k22 = second[i+1][j+1][0];
            k11 = second[i-1][j][0];
            k12 = second[i-1][j-1][0];
            
            d21 = depth[i+1][j];
            d22 = depth[i+1][j+1];
            d11 = depth[i-1][j];
            d12 = depth[i-1][j-1];
          }
          t = t2;
        }
        else if(-1 < t2 && t2 <= 0){
          if(ty > 0){
            k11 = second[i+1][j][0];
            k12 = second[i+1][j-1][0];
            k21 = second[i-1][j][0];
            k22 = second[i-1][j+1][0];
            
            d11 = depth[i+1][j];
            d12 = depth[i+1][j-1];
            d21 = depth[i-1][j];
            d22 = depth[i-1][j+1];
          }
          else{
            k21 = second[i+1][j][0];
            k22 = second[i+1][j-1][0];
            k11 = second[i-1][j][0];
            k12 = second[i-1][j+1][0];
            
            d21 = depth[i+1][j];
            d22 = depth[i+1][j-1];
            d11 = depth[i-1][j];
            d12 = depth[i-1][j+1];
          }
          t = -t2;
        }
        else
          continue;
        float kf = (1.0f-t)*k11 + t*k12;
        float kb = (1.0f-t)*k21 + t*k22;
        
        float r = radius[i][j];
        float dc = depth[i][j];
        float df = (1.0f-t)*d11 + t*d12;
        float db = (1.0f-t)*d21 + t*d22;
        
        if(((Kr*kf <= 0 /* || Kr*kb <= 0*/) && kf - kb <= 0) ||
           db - dc > 2*r/DERIVATIVE /*|| fabs(db - dc) > r*/){
          value[i][j] = SILHOUETTE;
        }
      }
    }
  }
  
  void createSC2(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        totalW[i][j] = 
	  normals[i][j][0] = normals[i][j][1] = normals[i][j][2] = 0;
        value[i][j] = EMPTY;
      }
    }
    
    computeFrontDepth();
    
    leaf_counter = 0;
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      bft->support *= DERIVATIVE;
      float r = bft->support;
      
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ;// - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX*DERIVATIVE){
        transPoly(bft, bf);
        
        if(!both && bft->normalZAtC() > 0)
          continue;
        
        leaf_counter++;
        
        r = bft->support;
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            float f = front[j][k];
            if((k-cu)*(k-cu) + (f-cd)*(f-cd) > dj)
              continue;
            
            float qd;
            if(!bft->rayIntersectLAZ(qd, (float)k, (float)j))
              continue;
            
            float n1[3];
            bft->normalLA(n1, (float)k, (float)j, qd);
            if(n1[2] > 0){
              if(!both)
                continue;
              n1[0] = -n1[0];
              n1[1] = -n1[1];
              n1[2] = -n1[2];
            }
            
            float w = -(float)bft->weight((float)k, (float)j, qd);
            float *n = normals[j][k];
            n[0] += w*n1[0];
            n[1] += w*n1[1];
            n[2] += w*n1[2];
            
            totalW[j][k] += w;
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfi;
          }
        }
      }
    }
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        float W = totalW[i][j];
        float *n = normals[i][j];
        if(W == 0){
          n[2] = 1000000;
          continue;
        }
        
        float len = (float)(1.0/sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]));
        n[0] *= len;
        n[1] *= len;
        n[2] *= len;
        
        value[i][j] = NORMAL;
      }
    }
    
    for(i=1; i<height-1; i++){
      for(j=1; j<width-1; j++){
        if(value[i][j] == EMPTY)
          continue;
        
        float* n = normals[i][j];
        
        if(n[2] > 0.94f)
          continue;
        
        float vn = n[2];
        float tx = n[0];
        float ty = n[1];
        
        float t1 = ty/tx;
        float t2 = tx/ty;
        float k11, k12, k21, k22, t;
        if(0 <= t1 && t1 < 1){
          k11 = normals[i][j+1][2];
          k12 = normals[i+1][j+1][2];
          k21 = normals[i][j-1][2];
          k22 = normals[i-1][j-1][2];
          t = t1;
        }
        else if(-1 <= t1 && t1 < 0){
          k11 = normals[i][j+1][2];
          k12 = normals[i-1][j+1][2];
          k21 = normals[i][j-1][2];
          k22 = normals[i+1][j-1][2];
          t = -t1;
        }
        else if(0 <= t2 && t2 <= 1){
          k11 = normals[i+1][j][2];
          k12 = normals[i+1][j+1][2];
          k21 = normals[i-1][j][2];
          k22 = normals[i-1][j-1][2];
          t = t2;
        }
        else if(-1 < t2 && t2 <= 0){
          k11 = normals[i-1][j][2];
          k12 = normals[i-1][j+1][2];
          k21 = normals[i+1][j][2];
          k22 = normals[i+1][j-1][2];
          t = -t2;
        }
        else
          continue;
        if(vn < (1.0-t)*k11 + t*k12 - 0.0001f && vn < (1.0-t)*k21 + t*k22 - 0.0001f){
          value[i][j] = SILHOUETTE;
          /*
            if(value[i][j+1] != EMPTY)
	    value[i][j+1] = RIDGE;
            if(value[i+1][j+1] != EMPTY)
	    value[i+1][j+1] = RIDGE;
            if(value[i+1][j] != EMPTY)
	    value[i+1][j] = RIDGE;
            if(value[i+1][j-1] != EMPTY)
	    value[i+1][j-1] = RIDGE;
            if(value[i][j-1] != EMPTY)
	    value[i][j-1] = RIDGE;
            if(value[i-1][j-1] != EMPTY)
	    value[i-1][j-1] = RIDGE;
            if(value[i-1][j] != EMPTY)
	    value[i-1][j] = RIDGE;
            if(value[i][j+1] != EMPTY)
	    value[i][j+1] = RIDGE;*/
        }
      }
    }
  }
  
  void createRidgeValleyByCanny(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        totalW[i][j] = radius[i][j] = depth[i][j] = 
	  normals[i][j][0] = normals[i][j][1] = normals[i][j][2] =
	  second[i][j][0] = second[i][j][1] = second[i][j][2] =
	  second[i][j][3] = second[i][j][4] = second[i][j][5] = 0;
        value[i][j] = EMPTY;
      }
    }
    
    computeFrontDepth();
    
    leaf_counter = 0;
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      bft->support *= DERIVATIVE;
      float r = bft->support;
      
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ;// - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX*DERIVATIVE){
        transPoly(bft, bf);
        
        if(!both && bft->normalZAtC() > 0)
          continue;
        
        leaf_counter++;
        
        r = bft->support;
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            float f = front[j][k];
            if((k-cu)*(k-cu) + (f-cd)*(f-cd) > dj)
              continue;
            
            float qd;
            if(!bft->rayIntersectLAZ(qd, (float)k, (float)j))
              continue;
            
            float n1[3], gg1[6];
            bft->gradient12LA(n1, gg1, (float)k, (float)j, qd);
            if(n1[2] > 0){
              if(!both)
                continue;
              n1[0] = -n1[0];
              n1[1] = -n1[1];
              n1[2] = -n1[2];
              gg1[0] = -gg1[0];
              gg1[1] = -gg1[1];
              gg1[2] = -gg1[2];
              gg1[3] = -gg1[3];
              gg1[4] = -gg1[4];
              gg1[5] = -gg1[5];
            }
            float w = (float)bft->weight((float)k, (float)j, qd);
            float *n = normals[j][k];
            n[0] -= w*n1[0];
            n[1] -= w*n1[1];
            n[2] -= w*n1[2];
            float *gg = second[j][k];
            gg[0] -= w*gg1[0];
            gg[1] -= w*gg1[1];
            gg[2] -= w*gg1[2];
            gg[3] -= w*gg1[3];
            gg[4] -= w*gg1[4];
            gg[5] -= w*gg1[5];
            depth[j][k] += w*qd;
            radius[j][k] += w*bft->support;
            totalW[j][k] += w;
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfi;
          }
        }
      }
    }
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        float W = totalW[i][j];
        if(W == 0)
          continue;
        
        float *g = normals[i][j];
        float *gg = second[i][j];
        W = 1.0f/W;
        g[0] *= W;
        g[1] *= W;
        g[2] *= W;
        gg[0] *= W;
        gg[1] *= W;
        gg[2] *= W;
        gg[3] *= W;
        gg[4] *= W;
        gg[5] *= W;
        depth[i][j] *= W;
        radius[i][j] *= W;
        
	/* temporary derivatives*/
        float fxfx = g[0]*g[0];
        float fxfy = g[0]*g[1];
        float fxfz = g[0]*g[2];
        float fyfy = g[1]*g[1];
        float fyfz = g[1]*g[2];
        float fzfz = g[2]*g[2];
        
        float fxx = gg[0];
        float fxy = gg[1];
        float fxz = gg[2];
        float fyy = gg[3];
        float fyz = gg[4];
        float fzz = gg[5];
        
        float g2 = g[0]*g[0]+g[1]*g[1]+g[2]*g[2];
        float g1 = (float)sqrt(g2);
        if(g1 == 0)
          continue;
        float g3 = g2*g1;
        float g4 = g2*g2;
        
        g[0] /= g1;
        g[1] /= g1;
        g[2] /= g1;
        
        /* mean and gaussian curvatures */
        float H = (fxx*(fyfy+fzfz) + fyy*(fxfx+fzfz) + fzz*(fxfx+fyfy)
		   - 2*(fxy*fxfy+fxz*fxfz+fyz*fyfz))/2;
        H /= g3;
        
        float K = fxfx*(fyy*fzz-fyz*fyz)
          + fyfy*(fxx*fzz-fxz*fxz)
	  + fzfz*(fxx*fyy-fxy*fxy)
	  + 2*(fxfy*(fxz*fyz-fxy*fzz)
               + fxfz*(fxy*fyz-fxz*fyy)
	       + fyfz*(fxy*fxz-fxx*fyz));
        K /= g4;
        
        /* principal curvatures */
        float discr = (float)sqrt(fabs(H*H-K));
        
        float Kmax = H + discr;
        
        float Kmin = H - discr;
        
        /* matrix entries */
        float m11 = ((-1 + fxfx/g2)*fxx)/g1 + (fxfy*fxy)/g3 + (fxfz*fxz)/g3;
        float m12 = ((-1 + fxfx/g2)*fxy)/g1 + (fxfy*fyy)/g3 + (fxfz*fyz)/g3;
        float m13 = ((-1 + fxfx/g2)*fxz)/g1 + (fxfy*fyz)/g3 + (fxfz*fzz)/g3;
        float m21 = (fxfy*fxx)/g3 + ((-1 + fyfy/g2)*fxy)/g1 + (fyfz*fxz)/g3;
        float m22 = (fxfy*fxy)/g3 + ((-1 + fyfy/g2)*fyy)/g1 + (fyfz*fyz)/g3;
        float m23 = (fxfy*fxz)/g3 + ((-1 + fyfy/g2)*fyz)/g1 + (fyfz*fzz)/g3;
        float m31 = (fxfz*fxx)/g3 + (fyfz*fxy)/g3 + ((-1 + fzfz/g2)*fxz)/g1;
        float m32 = (fxfz*fxy)/g3 + (fyfz*fyy)/g3 + ((-1 + fzfz/g2)*fyz)/g1;
        float m33 = (fxfz*fxz)/g3 + (fyfz*fyz)/g3 + ((-1 + fzfz/g2)*fzz)/g1;
        
        /* solve for eigenvectors */
        float tmp1 = m11+Kmax;
        float tmp2 = m22+Kmax;
        float tmp3 = m33+Kmax;
        
        float ux[3], uy[3], uz[3], len[3];
        ux[0] = m12*m23-m13*tmp2;
        uy[0] = m13*m21-m23*tmp1;
        uz[0] = tmp1*tmp2-m12*m21;
        len[0] = (float)sqrt(ux[0]*ux[0]+uy[0]*uy[0]+uz[0]*uz[0]);
        
        ux[1] = m12*tmp3-m13*m32;
        uy[1] = m13*m31-tmp1*tmp3;
        uz[1] = tmp1*m32-m12*m31;
        len[1] = (float)sqrt(ux[1]*ux[1]+uy[1]*uy[1]+uz[1]*uz[1]);
        
        ux[2] = tmp2*tmp3-m23*m32;
        uy[2] = m23*m31-m21*tmp3;
        uz[2] = m21*m32-m31*tmp2;
        len[2] = (float)sqrt(ux[2]*ux[2]+uy[2]*uy[2]+uz[2]*uz[2]);
        
        int index = 0;
        double max = len[0];
        if ( len[1] > max ) {
          index  = 1;
          max = len[1];
        }
        if ( len[2] > max ) {
          index = 2;
          max = len[2];
        }
        
        float Tmax[3];
        Tmax[0] = ux[index]/len[index];
        Tmax[1] = uy[index]/len[index];
        Tmax[2] = uz[index]/len[index];
        
        /* second tangent is cross product of first tangent and normal */
        float Tmin[3];
        Tmin[0] = Tmax[1]*g[2]-Tmax[2]*g[1];
        Tmin[1] = Tmax[2]*g[0]-Tmax[0]*g[2];
        Tmin[2] = Tmax[0]*g[1]-Tmax[1]*g[0];
        
        value[i][j] = NORMAL;
        
        gg[0] = Kmax;
        gg[1] = Tmax[0];
        gg[2] = Tmax[1];
        
        gg[3] = Kmin;
        gg[4] = Tmin[0];
        gg[5] = Tmin[1];
      }
    }
    
    for(i=1; i<height-1; i++){
      for(j=1; j<width-1; j++){
        if(value[i][j] == EMPTY && value[i][j+1] == EMPTY &&
           value[i+1][j] == EMPTY)
          continue;
        
        if(value[i][j] != EMPTY && value[i][j+1] != EMPTY &&
           value[i+1][j] != EMPTY){
          if(fabs(depth[i][j] - depth[i][j+1]) > radius[i][j] ||
             fabs(depth[i][j] - depth[i+1][j]) > radius[i][j]){
            value[i][j] = SILHOUETTE;
            continue;
          }
        }
        else{
          value[i][j] = SILHOUETTE;
          continue;
        }
        
        float* tensor = second[i][j];
        float Kmax = tensor[0];
        float Kmin = tensor[3];
        if(Kmax > 0 && Kmax > fabs(Kmin)){
          float tx = tensor[1];
          float ty = tensor[2];
          float t1 = ty/tx;
          float t2 = tx/ty;
          float k11, k12, k21, k22, t;
          if(0 <= t1 && t1 < 1){
            k11 = second[i][j+1][0];
            k12 = second[i+1][j+1][0];
            k21 = second[i][j-1][0];
            k22 = second[i-1][j-1][0];
            t = t1;
          }
          else if(-1 <= t1 && t1 < 0){
            k11 = second[i][j+1][0];
            k12 = second[i-1][j+1][0];
            k21 = second[i][j-1][0];
            k22 = second[i+1][j-1][0];
            t = -t1;
          }
          else if(0 <= t2 && t2 <= 1){
            k11 = second[i+1][j][0];
            k12 = second[i+1][j+1][0];
            k21 = second[i-1][j][0];
            k22 = second[i-1][j-1][0];
            t = t2;
          }
          else if(-1 < t2 && t2 <= 0){
            k11 = second[i-1][j][0];
            k12 = second[i-1][j+1][0];
            k21 = second[i+1][j][0];
            k22 = second[i+1][j-1][0];
            t = -t2;
          }
          else
            continue;
          if(Kmax > (1.0-t)*k11 + t*k12 && Kmax > (1.0-t)*k21 + t*k22){
            value[i][j] = RIDGE;
            second[i][j][1] = Kmax;// - 0.5f*((1.0-t)*(k11+k21) + t*(k12*k22));
          }
        }
        else if(Kmin < 0 && -Kmin > fabs(Kmax)){
          float tx = tensor[4];
          float ty = tensor[5];
          float t1 = ty/tx;
          float t2 = tx/ty;
          float k11, k12, k21, k22, t;
          if(0 <= t1 && t1 < 1){
            k11 = second[i][j+1][3];
            k12 = second[i+1][j+1][3];
            k21 = second[i][j-1][3];
            k22 = second[i-1][j-1][3];
            t = t1;
          }
          else if(-1 <= t1 && t1 < 0){
            k11 = second[i][j+1][3];
            k12 = second[i-1][j+1][3];
            k21 = second[i][j-1][3];
            k22 = second[i+1][j-1][3];
            t = -t1;
          }
          else if(0 <= t2 && t2 <= 1){
            k11 = second[i+1][j][3];
            k12 = second[i+1][j+1][3];
            k21 = second[i-1][j][3];
            k22 = second[i-1][j-1][3];
            t = t2;
          }
          else if(-1 < t2 && t2 <= 0){
            k11 = second[i-1][j][3];
            k12 = second[i-1][j+1][3];
            k21 = second[i+1][j][3];
            k22 = second[i+1][j-1][3];
            t = -t2;
          }
          else
            continue;
          if(Kmin < (1.0-t)*k11 + t*k12 && Kmin < (1.0-t)*k21 + t*k22){
            value[i][j] = VALLEY;
            second[i][j][4] = Kmin;// - 0.5f*((1.0-t)*(k11+k21) + t*(k12*k22));
          }
        }
      }
    }
  }
  
  void createRidgeValley(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        totalW[i][j] = depth[i][j] = 
	  normals[i][j][0] = normals[i][j][1] = normals[i][j][2] =
	  second[i][j][0] = second[i][j][1] = second[i][j][2] =
	  second[i][j][3] = second[i][j][4] = second[i][j][5] =
	  third[i][j][0] = third[i][j][1] = third[i][j][2] =
	  third[i][j][3] = third[i][j][4] = third[i][j][5] =
	  third[i][j][6] = third[i][j][7] = third[i][j][8] =
	  third[i][j][9] = 0;
        value[i][j] = EMPTY;
      }
    }
    
    computeFrontDepth();
    
    leaf_counter = 0;
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      bft->support *= DERIVATIVE;
      float r = bft->support;
      
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ;// - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX*DERIVATIVE){
        transPoly(bft, bf);
        
        if(!both && bft->normalZAtC() > 0)
          continue;
        
        leaf_counter++;
        
        r = bft->support;
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            float f = front[j][k];
            if((k-cu)*(k-cu) + (f-cd)*(f-cd) > dj)
              continue;
            
            float qd;
            if(!bft->rayIntersectLAZ(qd, (float)k, (float)j))
              continue;
            
            float n1[3], gg1[6], ggg1[10];
            bft->gradient123LA(n1, gg1, ggg1, (float)k, (float)j, qd);
            if(n1[2] > 0){
              if(!both)
                continue;
              n1[0] = -n1[0];
              n1[1] = -n1[1];
              n1[2] = -n1[2];
              gg1[0] = -gg1[0];
              gg1[1] = -gg1[1];
              gg1[2] = -gg1[2];
              gg1[3] = -gg1[3];
              gg1[4] = -gg1[4];
              gg1[5] = -gg1[5];
              ggg1[0] = -ggg1[0];
              ggg1[1] = -ggg1[1];
              ggg1[2] = -ggg1[2];
              ggg1[3] = -ggg1[3];
              ggg1[4] = -ggg1[4];
              ggg1[5] = -ggg1[5];
              ggg1[6] = -ggg1[6];
              ggg1[7] = -ggg1[7];
              ggg1[8] = -ggg1[8];
              ggg1[9] = -ggg1[9];
            }
            float w = (float)bft->weight((float)k, (float)j, qd);
            float *n = normals[j][k];
            n[0] -= w*n1[0];
            n[1] -= w*n1[1];
            n[2] -= w*n1[2];
            float *gg = second[j][k];
            gg[0] -= w*gg1[0];
            gg[1] -= w*gg1[1];
            gg[2] -= w*gg1[2];
            gg[3] -= w*gg1[3];
            gg[4] -= w*gg1[4];
            gg[5] -= w*gg1[5];
            float *ggg = third[j][k];
            ggg[0] -= w*ggg1[0];
            ggg[1] -= w*ggg1[1];
            ggg[2] -= w*ggg1[2];
            ggg[3] -= w*ggg1[3];
            ggg[4] -= w*ggg1[4];
            ggg[5] -= w*ggg1[5];
            ggg[6] -= w*ggg1[6];
            ggg[7] -= w*ggg1[7];
            ggg[8] -= w*ggg1[8];
            ggg[9] -= w*ggg1[9];
            depth[j][k] += w*qd;
            totalW[j][k] += w;
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfi;
          }
        }
      }
    }
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        float W = totalW[i][j];
        if(W == 0)
          continue;
        
        value[i][j] = NORMAL;
        
        W = 1.0f/W;
        float *g = normals[i][j];
        g[0] *= W;
        g[1] *= W;
        g[2] *= W;
        float *gg = second[i][j];
        gg[0] *= W;
        gg[1] *= W;
        gg[2] *= W;
        gg[3] *= W;
        gg[4] *= W;
        gg[5] *= W;
        float *ggg = third[i][j];
        ggg[0] *= W;
        ggg[1] *= W;
        ggg[2] *= W;
        ggg[3] *= W;
        ggg[4] *= W;
        ggg[5] *= W;
        ggg[6] *= W;
        ggg[7] *= W;
        ggg[8] *= W;
        ggg[9] *= W;
        depth[i][j] *= W;
        
        float Kmax, Kmin, Rmax, Rmin, Tmax[3], Tmin[3];
        curvatureDerivative(Kmax, Kmin, Rmax, Rmin, Tmax, Tmin,
                            g, gg, ggg);
        gg[0] = Kmax;
        gg[1] = Tmax[0];
        gg[2] = Tmax[1];
        gg[3] = Kmin;
        gg[4] = Tmin[0];
        gg[5] = Tmin[1];
        ggg[0] = Rmax;
        ggg[1] = Rmin;
      }
    }
    
    normalizeNormals();
    
    for(i=1; i<height-1; i++){
      for(j=1; j<width-1; j++){
        if(value[i][j] == EMPTY)
          continue;
        
        if((value[i][j+1] != EMPTY &&
            ridgeTest(second[i][j], second[i][j+1],
                      third[i][j][0], third[i][j+1][0], 0)) ||
           (value[i+1][j] != EMPTY &&
            ridgeTest(second[i][j], second[i+1][j],
                      third[i][j][0], third[i+1][j][0], 1)))
          value[i][j] = RIDGE;
        
        if((value[i][j+1] != EMPTY &&
            valleyTest(second[i][j], second[i][j+1],
                       third[i][j][1], third[i][j+1][1], 0)) ||
           (value[i+1][j] != EMPTY &&
            valleyTest(second[i][j], second[i+1][j],
                       third[i][j][1], third[i+1][j][1], 1)))
          value[i][j] = VALLEY;
      }
    }
  }
  
  //mode = 1 mean 
  //mode = 2 gauss
  //mode = 3 k_max
  //mode = 4 k_min
  //mode = 5 total
  //mode = 6 radial curvature
  void createCurvature(int mode){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        front[i][j] = 100000;
        totalW[i][j] = 
	  normals[i][j][0] = normals[i][j][1] = normals[i][j][2] =
	  second[i][j][0] = second[i][j][1] = second[i][j][2] =
	  second[i][j][3] = second[i][j][4] = second[i][j][5] = 0;
        value[i][j] = EMPTY;
      }
    }
    
    computeFrontDepth();
    
    leaf_counter = 0;
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      bft->support *= DERIVATIVE;
      float r = bft->support;
      
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ;// - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX*DERIVATIVE){
        transPoly(bft, bf);
        
        if(!both && bft->normalZAtC() > 0)
          continue;
        
        leaf_counter++;
        
        r = bft->support;
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            float f = front[j][k];
            if((k-cu)*(k-cu) + (f-cd)*(f-cd) > dj)
              continue;
            
            float qd;
            if(!bft->rayIntersectLAZ(qd, (float)k, (float)j))
              continue;
            
            float n1[3], gg1[6];
            bft->gradient12LA(n1, gg1, (float)k, (float)j, qd);
            if(n1[2] > 0){
              if(!both)
                continue;
              n1[0] = -n1[0];
              n1[1] = -n1[1];
              n1[2] = -n1[2];
              gg1[0] = -gg1[0];
              gg1[1] = -gg1[1];
              gg1[2] = -gg1[2];
              gg1[3] = -gg1[3];
              gg1[4] = -gg1[4];
              gg1[5] = -gg1[5];
            }
            float w = (float)bft->weight((float)k, (float)j, qd);
            float *n = normals[j][k];
            n[0] -= w*n1[0];
            n[1] -= w*n1[1];
            n[2] -= w*n1[2];
            float *gg = second[j][k];
            gg[0] -= w*gg1[0];
            gg[1] -= w*gg1[1];
            gg[2] -= w*gg1[2];
            gg[3] -= w*gg1[3];
            gg[4] -= w*gg1[4];
            gg[5] -= w*gg1[5];
            totalW[j][k] += w;
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfi;
          }
        }
      }
    }
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        float W = totalW[i][j];
        if(W == 0)
          continue;
        
        float *g = normals[i][j];
        float *gg = second[i][j];
        W = 1.0f/W;
        g[0] *= W;
        g[1] *= W;
        g[2] *= W;
        gg[0] *= W;
        gg[1] *= W;
        gg[2] *= W;
        gg[3] *= W;
        gg[4] *= W;
        gg[5] *= W;
        
	/* temporary derivatives*/
        float fxfx = g[0]*g[0];
        float fxfy = g[0]*g[1];
        float fxfz = g[0]*g[2];
        float fyfy = g[1]*g[1];
        float fyfz = g[1]*g[2];
        float fzfz = g[2]*g[2];
        
        float fxx = gg[0];
        float fxy = gg[1];
        float fxz = gg[2];
        float fyy = gg[3];
        float fyz = gg[4];
        float fzz = gg[5];
        
        float g2 = g[0]*g[0]+g[1]*g[1]+g[2]*g[2];
        float g1 = (float)sqrt(g2);
        if(g1 == 0)
          continue;
        float g3 = g2*g1;
        float g4 = g2*g2;
        
        g[0] /= g1;
        g[1] /= g1;
        g[2] /= g1;
        
        /* mean and gaussian curvatures */
        float H = (fxx*(fyfy+fzfz) + fyy*(fxfx+fzfz) + fzz*(fxfx+fyfy)
		   - 2*(fxy*fxfy+fxz*fxfz+fyz*fyfz))/2;
        H /= g3;
        
        if(mode == 1){
          value[i][j] = -H;
          continue;
        }
        
        float K = fxfx*(fyy*fzz-fyz*fyz)
          + fyfy*(fxx*fzz-fxz*fxz)
	  + fzfz*(fxx*fyy-fxy*fxy)
	  + 2*(fxfy*(fxz*fyz-fxy*fzz)
               + fxfz*(fxy*fyz-fxz*fyy)
	       + fyfz*(fxy*fxz-fxx*fyz));
        K /= g4;
        
        if(mode == 2){
          if(K > 0)
            value[i][j] = (float)sqrt(fabs(K));
          else
            value[i][j] = -(float)sqrt(fabs(K));
          continue;
        }
        
        /* principal curvatures */
        float discr = (float)sqrt(fabs(H*H-K));
        
        float Kmax = H + discr;
        
        if(mode == 3){
          value[i][j] = -Kmax;
          continue;
        }
        
        float Kmin = H - discr;
        
        if(mode == 4){
          value[i][j] = -Kmin;
          continue;
        }
        
        if(mode == 5){
          value[i][j] = (float)sqrt(Kmax*Kmax + Kmin*Kmin);
          continue;
        }
        
        /* matrix entries */
        float m11 = ((-1 + fxfx/g2)*fxx)/g1 + (fxfy*fxy)/g3 + (fxfz*fxz)/g3;
        float m12 = ((-1 + fxfx/g2)*fxy)/g1 + (fxfy*fyy)/g3 + (fxfz*fyz)/g3;
        float m13 = ((-1 + fxfx/g2)*fxz)/g1 + (fxfy*fyz)/g3 + (fxfz*fzz)/g3;
        float m21 = (fxfy*fxx)/g3 + ((-1 + fyfy/g2)*fxy)/g1 + (fyfz*fxz)/g3;
        float m22 = (fxfy*fxy)/g3 + ((-1 + fyfy/g2)*fyy)/g1 + (fyfz*fyz)/g3;
        float m23 = (fxfy*fxz)/g3 + ((-1 + fyfy/g2)*fyz)/g1 + (fyfz*fzz)/g3;
        float m31 = (fxfz*fxx)/g3 + (fyfz*fxy)/g3 + ((-1 + fzfz/g2)*fxz)/g1;
        float m32 = (fxfz*fxy)/g3 + (fyfz*fyy)/g3 + ((-1 + fzfz/g2)*fyz)/g1;
        float m33 = (fxfz*fxz)/g3 + (fyfz*fyz)/g3 + ((-1 + fzfz/g2)*fzz)/g1;
        
        /* solve for eigenvectors */
        float tmp1 = m11+Kmax;
        float tmp2 = m22+Kmax;
        float tmp3 = m33+Kmax;
        
        float ux[3], uy[3], uz[3], len[3];
        ux[0] = m12*m23-m13*tmp2;
        uy[0] = m13*m21-m23*tmp1;
        uz[0] = tmp1*tmp2-m12*m21;
        len[0] = (float)sqrt(ux[0]*ux[0]+uy[0]*uy[0]+uz[0]*uz[0]);
        
        ux[1] = m12*tmp3-m13*m32;
        uy[1] = m13*m31-tmp1*tmp3;
        uz[1] = tmp1*m32-m12*m31;
        len[1] = (float)sqrt(ux[1]*ux[1]+uy[1]*uy[1]+uz[1]*uz[1]);
        
        ux[2] = tmp2*tmp3-m23*m32;
        uy[2] = m23*m31-m21*tmp3;
        uz[2] = m21*m32-m31*tmp2;
        len[2] = (float)sqrt(ux[2]*ux[2]+uy[2]*uy[2]+uz[2]*uz[2]);
        
        int index = 0;
        double max = len[0];
        if ( len[1] > max ) {
          index  = 1;
          max = len[1];
        }
        if ( len[2] > max ) {
          index = 2;
          max = len[2];
        }
        
        float Tmax[3];
        Tmax[0] = ux[index]/len[index];
        Tmax[1] = uy[index]/len[index];
        Tmax[2] = uz[index]/len[index];
        
        /* second tangent is cross product of first tangent and normal */
        float Tmin[3];
        Tmin[0] = Tmax[1]*g[2]-Tmax[2]*g[1];
        Tmin[1] = Tmax[2]*g[0]-Tmax[0]*g[2];
        Tmin[2] = Tmax[0]*g[1]-Tmax[1]*g[0];
        
        value[i][j] = NORMAL;
        
        if(mode == 6){
          float c = (float)cos(atan2(Tmin[2], Tmax[2]));
          c *= c;
          value[i][j] = -(Kmax*c + Kmin*(1.0f-c));
        }
      }
    }
  }
  
  void transformM(float R[9], float s, float t[3]){
    for(int i=0; i<mesh->vertexN; i++){
      float *p = mesh->vertex[i];
      MAT_BY_VEC(p, R, mesh->o_vertex[i]);
      p[0] = s*p[0] + t[0];
      p[1] = s*p[1] + t[1];
      p[2] = s*p[2] + t[2];
      MAT_BY_VEC(mesh->normal[i], R, mesh->o_normal[i]);
    }
  }
  
  void transformP(float R[9], float s, float t[3]){
    for(int i=0; i<pointN; i++){
      float *p = point[i];
      MAT_BY_VEC(p, R, o_point[i]);
      p[0] = s*p[0] + t[0];
      p[1] = s*p[1] + t[1];
      p[2] = s*p[2] + t[2];
      MAT_BY_VEC(normal[i], R, o_normal[i]);
    }
  }
  
  void flipOrientation(){
    for(int i=0; i<2*bfListN; i++){
      for(int j=0; j<bfNs[i]; j++){
        bf_lists[i][j]->flipOrientation();
      }
    }
  }
  
  inline void transPoly(BasisFunction* bf_out, BasisFunction* bf_in){
    if(degree == 3){
      Cubic* b1 = (Cubic*)bf_in;
      Cubic* b2 = (Cubic*)bf_out;
      
      b2->cXXX = is3*(Rxx*Rxx*Rxx*b1->cXXX + Ryx*Ryx*Ryx*b1->cYYY + Rzx*Rzx*Rzx*b1->cZZZ + Rxx*Rxx*Ryx*b1->cXXY + Ryx*Ryx*Rzx*b1->cYYZ + Rzx*Rzx*Rxx*b1->cZZX + Rxx*Ryx*Ryx*b1->cXYY + Ryx*Rzx*Rzx*b1->cYZZ + Rzx*Rxx*Rxx*b1->cZXX + Rxx*Ryx*Rzx*b1->cXYZ);
      b2->cYYY = is3*(Rxy*Rxy*Rxy*b1->cXXX + Ryy*Ryy*Ryy*b1->cYYY + Rzy*Rzy*Rzy*b1->cZZZ + Rxy*Rxy*Ryy*b1->cXXY + Ryy*Ryy*Rzy*b1->cYYZ + Rzy*Rzy*Rxy*b1->cZZX + Rxy*Ryy*Ryy*b1->cXYY + Ryy*Rzy*Rzy*b1->cYZZ + Rzy*Rxy*Rxy*b1->cZXX + Rxy*Ryy*Rzy*b1->cXYZ);
      b2->cZZZ = is3*(Rxz*Rxz*Rxz*b1->cXXX + Ryz*Ryz*Ryz*b1->cYYY + Rzz*Rzz*Rzz*b1->cZZZ + Rxz*Rxz*Ryz*b1->cXXY + Ryz*Ryz*Rzz*b1->cYYZ + Rzz*Rzz*Rxz*b1->cZZX + Rxz*Ryz*Ryz*b1->cXYY + Ryz*Rzz*Rzz*b1->cYZZ + Rzz*Rxz*Rxz*b1->cZXX + Rxz*Ryz*Rzz*b1->cXYZ);
      
      b2->cXXY = is3*(3.0f*(Rxx*Rxx*Rxy*b1->cXXX + Ryx*Ryx*Ryy*b1->cYYY + Rzx*Rzx*Rzy*b1->cZZZ) + (2.0f*Rxx*Rxy*Ryx + Rxx*Rxx*Ryy)*b1->cXXY + (2.0f*Ryx*Ryy*Rzx + Ryx*Ryx*Rzy)*b1->cYYZ + (2.0f*Rzx*Rzy*Rxx + Rzx*Rzx*Rxy)*b1->cZZX + (Rxx*Ryy*Ryx + 2.0f*Rxx*Ryx*Ryy)*b1->cXYY + (Ryx*Rzy*Rzx + 2.0f*Ryx*Rzx*Rzy)*b1->cYZZ + (Rzx*Rxy*Rxx + 2.0f*Rzx*Rxx*Rxy)*b1->cZXX + (Rxy*Ryx*Rzx + Rxx*Ryy*Rzx + Rxx*Ryx*Rzy)*b1->cXYZ);
      b2->cYYZ = is3*(3.0f*(Rxy*Rxy*Rxz*b1->cXXX + Ryy*Ryy*Ryz*b1->cYYY + Rzy*Rzy*Rzz*b1->cZZZ) + (2.0f*Rxy*Rxz*Ryy + Rxy*Rxy*Ryz)*b1->cXXY + (2.0f*Ryy*Ryz*Rzy + Ryy*Ryy*Rzz)*b1->cYYZ + (2.0f*Rzy*Rzz*Rxy + Rzy*Rzy*Rxz)*b1->cZZX + (Rxy*Ryz*Ryy + 2.0f*Rxy*Ryy*Ryz)*b1->cXYY + (Ryy*Rzz*Rzy + 2.0f*Ryy*Rzy*Rzz)*b1->cYZZ + (Rzy*Rxz*Rxy + 2.0f*Rzy*Rxy*Rxz)*b1->cZXX + (Rxz*Ryy*Rzy + Rxy*Ryz*Rzy + Rxy*Ryy*Rzz)*b1->cXYZ);
      b2->cZZX = is3*(3.0f*(Rxz*Rxz*Rxx*b1->cXXX + Ryz*Ryz*Ryx*b1->cYYY + Rzz*Rzz*Rzx*b1->cZZZ) + (2.0f*Rxz*Rxx*Ryz + Rxz*Rxz*Ryx)*b1->cXXY + (2.0f*Ryz*Ryx*Rzz + Ryz*Ryz*Rzx)*b1->cYYZ + (2.0f*Rzz*Rzx*Rxz + Rzz*Rzz*Rxx)*b1->cZZX + (Rxz*Ryx*Ryz + 2.0f*Rxz*Ryz*Ryx)*b1->cXYY + (Ryz*Rzx*Rzz + 2.0f*Ryz*Rzz*Rzx)*b1->cYZZ + (Rzz*Rxx*Rxz + 2.0f*Rzz*Rxz*Rxx)*b1->cZXX + (Rxx*Ryz*Rzz + Rxz*Ryx*Rzz + Rxz*Ryz*Rzx)*b1->cXYZ);
      
      b2->cXYY = is3*(3.0f*(Rxx*Rxy*Rxy*b1->cXXX + Ryx*Ryy*Ryy*b1->cYYY + Rzx*Rzy*Rzy*b1->cZZZ) + (2.0f*Rxx*Rxy*Ryy + Rxy*Rxy*Ryx)*b1->cXXY + (2.0f*Ryx*Ryy*Rzy + Ryy*Ryy*Rzx)*b1->cYYZ + (2.0f*Rzx*Rzy*Rxy + Rzy*Rzy*Rxx)*b1->cZZX + (Rxx*Ryy*Ryy + 2.0f*Rxy*Ryy*Ryx)*b1->cXYY + (Ryx*Rzy*Rzy + 2.0f*Ryy*Rzy*Rzx)*b1->cYZZ + (Rzx*Rxy*Rxy + 2.0f*Rzy*Rxy*Rxx)*b1->cZXX + (Rxy*Ryy*Rzx + Rxy*Ryx*Rzy + Rxx*Ryy*Rzy)*b1->cXYZ);
      b2->cYZZ = is3*(3.0f*(Rxy*Rxz*Rxz*b1->cXXX + Ryy*Ryz*Ryz*b1->cYYY + Rzy*Rzz*Rzz*b1->cZZZ) + (2.0f*Rxy*Rxz*Ryz + Rxz*Rxz*Ryy)*b1->cXXY + (2.0f*Ryy*Ryz*Rzz + Ryz*Ryz*Rzy)*b1->cYYZ + (2.0f*Rzy*Rzz*Rxz + Rzz*Rzz*Rxy)*b1->cZZX + (Rxy*Ryz*Ryz + 2.0f*Rxz*Ryz*Ryy)*b1->cXYY + (Ryy*Rzz*Rzz + 2.0f*Ryz*Rzz*Rzy)*b1->cYZZ + (Rzy*Rxz*Rxz + 2.0f*Rzz*Rxz*Rxy)*b1->cZXX + (Rxz*Ryz*Rzy + Rxz*Ryy*Rzz + Rxy*Ryz*Rzz)*b1->cXYZ);
      b2->cZXX = is3*(3.0f*(Rxz*Rxx*Rxx*b1->cXXX + Ryz*Ryx*Ryx*b1->cYYY + Rzz*Rzx*Rzx*b1->cZZZ) + (2.0f*Rxz*Rxx*Ryx + Rxx*Rxx*Ryz)*b1->cXXY + (2.0f*Ryz*Ryx*Rzx + Ryx*Ryx*Rzz)*b1->cYYZ + (2.0f*Rzz*Rzx*Rxx + Rzx*Rzx*Rxz)*b1->cZZX + (Rxz*Ryx*Ryx + 2.0f*Rxx*Ryx*Ryz)*b1->cXYY + (Ryz*Rzx*Rzx + 2.0f*Ryx*Rzx*Rzz)*b1->cYZZ + (Rzz*Rxx*Rxx + 2.0f*Rzx*Rxx*Rxz)*b1->cZXX + (Rxx*Ryx*Rzz + Rxx*Ryz*Rzx + Rxz*Ryx*Rzx)*b1->cXYZ);
      
      b2->cXYZ = is3*(6.0f*(Rxx*Rxy*Rxz*b1->cXXX + Ryx*Ryy*Ryz*b1->cYYY + Rzx*Rzy*Rzz*b1->cZZZ) + 2.0f*(Rxy*Rxz*Ryx + Rxx*Rxz*Ryy + Rxx*Rxy*Ryz)*b1->cXXY + 2.0f*(Ryy*Ryz*Rzx + Ryx*Ryz*Rzy + Ryx*Ryy*Rzz)*b1->cYYZ + 2.0f*(Rzy*Rzz*Rxx + Rzx*Rzz*Rxy + Rzx*Rzy*Rxz)*b1->cZZX + 2.0f*(Rxz*Ryx*Ryy + Rxy*Ryx*Ryz + Rxx*Ryy*Ryz)*b1->cXYY + 2.0f*(Ryz*Rzx*Rzy + Ryy*Rzx*Rzz + Ryx*Rzy*Rzz)*b1->cYZZ + 2.0f*(Rzz*Rxy*Rxx + Rzy*Rxx*Rxz + Rzx*Rxy*Rxz)*b1->cZXX + (Rxz*Ryy*Rzx + Rxy*Ryz*Rzx + Rxz*Ryx*Rzy + Rxx*Ryz*Rzy + Rxy*Ryx*Rzz + Rxx*Ryy*Rzz)*b1->cXYZ);
      
      float Q[9];
      Q[0] = b1->cXX;
      Q[1] = Q[3] = 0.5f*b1->cXY;
      Q[2] = Q[6] = 0.5f*b1->cZX;
      Q[4] = b1->cYY;
      Q[5] = Q[7] = 0.5f*b1->cYZ;
      Q[8] = b1->cZZ;
      
      float QR[9];
      MAT_BY_MAT(QR, Q, Rt);
      MAT_BY_MAT(Q, R, QR);
      
      b2->cXX = Q[0]*is2;
      b2->cXY = 2.0f*Q[1]*is2;
      b2->cZX = 2.0f*Q[2]*is2;
      b2->cYY = Q[4]*is2;
      b2->cYZ = 2.0f*Q[5]*is2;
      b2->cZZ = Q[8]*is2;
      
      float b[3] = {b1->cX, b1->cY, b1->cZ};
      float Rb[3];
      MAT_BY_VEC(Rb, R, b);
      
      b2->cX = Rb[0]*is;
      b2->cY = Rb[1]*is;
      b2->cZ = Rb[2]*is;
      
      b2->c0 = b1->c0;
    }
    else if(degree == 2){
      Quadratic* b1 = (Quadratic*)bf_in;
      Quadratic* b2 = (Quadratic*)bf_out;
      
      float Q[9];
      Q[0] = b1->cXX;
      Q[1] = Q[3] = 0.5f*b1->cXY;
      Q[2] = Q[6] = 0.5f*b1->cZX;
      Q[4] = b1->cYY;
      Q[5] = Q[7] = 0.5f*b1->cYZ;
      Q[8] = b1->cZZ;
      
      float QR[9];
      MAT_BY_MAT(QR, Q, Rt);
      MAT_BY_MAT(Q, R, QR);
      
      b2->cXX = Q[0]*is2;
      b2->cXY = 2.0f*Q[1]*is2;
      b2->cZX = 2.0f*Q[2]*is2;
      b2->cYY = Q[4]*is2;
      b2->cYZ = 2.0f*Q[5]*is2;
      b2->cZZ = Q[8]*is2;
      
      float b[3] = {b1->cX, b1->cY, b1->cZ};
      float Rb[3];
      MAT_BY_VEC(Rb, R, b);
      
      b2->cX = Rb[0]*is;
      b2->cY = Rb[1]*is;
      b2->cZ = Rb[2]*is;
      
      b2->c0 = b1->c0;
    }
    else{
      Linear* b1 = (Linear*)bf_in;
      Linear* b2 = (Linear*)bf_out;
      
      float b[3] = {b1->cX, b1->cY, b1->cZ};
      float Rb[3];
      MAT_BY_VEC(Rb, R, b);
      
      b2->cX = Rb[0]*is;
      b2->cY = Rb[1]*is;
      b2->cZ = Rb[2]*is;
      
      b2->c0 = b1->c0;
    }
  }
  
  inline void MAT_BY_VEC(float y[3], float A[9], float x[3]){
    y[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2];
    y[1] = A[3]*x[0] + A[4]*x[1] + A[5]*x[2];
    y[2] = A[6]*x[0] + A[7]*x[1] + A[8]*x[2];
  }
  
  inline void MAT_BY_MAT(float C[9], float A[9], float B[9]){
    C[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
    C[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
    C[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];
    
    C[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
    C[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
    C[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];
    
    C[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
    C[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
    C[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
  }
  
  void convertColor(float color[3], float mid, float vari, float value){
    float s = (value-mid)/vari*2.0f;
    //float s = (value-mid)/vari*1.0f;
    /*if(s < -2.1f){
      color[0] = 0;
      color[1] = 0;
      color[2] = 0;
      }*/
    if(s < -2.1f){
      color[0] = 0;
      color[1] = 0;
      color[2] = 0.25f;
    }
    else if(s < -1.5f){
      color[0] = 0;
      color[1] = 0;
      color[2] = 0.25f + (s + 2.1f)/0.8f;
    }
    else if(s < -0.9f){
      color[0] = 0;
      color[1] = (s + 1.5f)/0.6f;
      color[2] = 1;
    }
    else if(s < -0.3f){
      color[0] = 0;
      color[1] = 1;
      color[2] = -(s + 0.3f)/0.6f;
    }
    else if(s < 0.3f){
      color[0] = (s + 0.3f)/0.6f;
      color[1] = 1;
      color[2] = 0;
    }
    else if(s < 0.9f){
      color[0] = 1;
      color[1] = -(s - 0.9f)/0.6f;
      color[2] = 0;
    }
    else if(s < 1.5f){
      color[0] = 1;
      color[1] = 0;
      color[2] = (s - 0.9f)/0.6f;
    }
    else if(s < 2.1f){
      color[0] = 1;
      color[1] = (s - 1.5f)/0.8f;
      color[2] = 1;
    }
    else{
      color[0] = 1;
      color[1] = 0.75;
      color[2] = 1;
    }
    /*
      else if(s < 2.1f){
      color[0] = 1;
      color[1] = (s - 1.5f)/0.6f;
      color[2] = 1;
      }
      else{
      color[0] = 1;
      color[1] = 1;
      color[2] = 1;
      }*/
  }
  
  void createSC2FromNormals(){
    int i,j;
    
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        float *n = normals[i][j];
        if(n[0] == 0 && n[1] == 0 && n[2] == 0){
          n[2] = 1000000;
          value[i][j] = EMPTY;
        }
        else
          value[i][j] = NORMAL;
      }
    }
    
    for(i=1; i<height-1; i++){
      for(j=1; j<width-1; j++){
        if(value[i][j] == EMPTY)
          continue;
        
        float* n = normals[i][j];
        
        if(n[2] > 0.94f)
          continue;
        
        float vn = n[2];
        float tx = n[0];
        float ty = n[1];
        
        float t1 = ty/tx;
        float t2 = tx/ty;
        float k11, k12, k21, k22, t;
        if(0 <= t1 && t1 < 1){
          k11 = normals[i][j+1][2];
          k12 = normals[i+1][j+1][2];
          k21 = normals[i][j-1][2];
          k22 = normals[i-1][j-1][2];
          t = t1;
        }
        else if(-1 <= t1 && t1 < 0){
          k11 = normals[i][j+1][2];
          k12 = normals[i-1][j+1][2];
          k21 = normals[i][j-1][2];
          k22 = normals[i+1][j-1][2];
          t = -t1;
        }
        else if(0 <= t2 && t2 <= 1){
          k11 = normals[i+1][j][2];
          k12 = normals[i+1][j+1][2];
          k21 = normals[i-1][j][2];
          k22 = normals[i-1][j-1][2];
          t = t2;
        }
        else if(-1 < t2 && t2 <= 0){
          k11 = normals[i-1][j][2];
          k12 = normals[i-1][j+1][2];
          k21 = normals[i+1][j][2];
          k22 = normals[i+1][j-1][2];
          t = -t2;
        }
        else
          continue;
        if(vn < (1.0-t)*k11 + t*k12 - 0.0001f && vn < (1.0-t)*k21 + t*k22 - 0.0001f){
          value[i][j] = SILHOUETTE;
        }
      }
    }
  }
  
  void curvatureDerivative(float &Kmax, float &Kmin,
                           float &Rmax, float &Rmin,
                           float Tmax[3], float Tmin[3],
                           float g[3], float gg[6], float ggg[10]){
    
    double pd1[3], pd2[3][3], pd3[3][3][3];
    
    pd1[0] = g[0];
    pd1[1] = g[1];
    pd1[2] = g[2];
    
    pd2[0][0] = gg[0];
    pd2[0][1] = gg[1];
    pd2[0][2] = gg[2];
    pd2[1][0] = gg[1];
    pd2[1][1] = gg[3];
    pd2[1][2] = gg[4];
    pd2[2][0] = gg[2];
    pd2[2][1] = gg[4];
    pd2[2][2] = gg[5];
    
    pd3[0][0][0] = ggg[0];
    pd3[0][0][1] = ggg[1];
    pd3[0][0][2] = ggg[2];
    pd3[0][1][0] = ggg[1];
    pd3[0][1][1] = ggg[3];
    pd3[0][1][2] = ggg[4];
    pd3[0][2][0] = ggg[2];
    pd3[0][2][1] = ggg[4];
    pd3[0][2][2] = ggg[5];
    
    pd3[1][0][0] = ggg[1];
    pd3[1][0][1] = ggg[3];
    pd3[1][0][2] = ggg[4];
    pd3[1][1][0] = ggg[3];
    pd3[1][1][1] = ggg[6];
    pd3[1][1][2] = ggg[7];
    pd3[1][2][0] = ggg[4];
    pd3[1][2][1] = ggg[7];
    pd3[1][2][2] = ggg[8];
    
    pd3[2][0][0] = ggg[2];
    pd3[2][0][1] = ggg[4];
    pd3[2][0][2] = ggg[5];
    pd3[2][1][0] = ggg[4];
    pd3[2][1][1] = ggg[7];
    pd3[2][1][2] = ggg[8];
    pd3[2][2][0] = ggg[5];
    pd3[2][2][1] = ggg[8];
    pd3[2][2][2] = ggg[9];
    
    /* temporary derivatives*/
    double fxfx = pd1[0]*pd1[0];
    double fxfy = pd1[0]*pd1[1];
    double fxfz = pd1[0]*pd1[2];
    double fyfy = pd1[1]*pd1[1];
    double fyfz = pd1[1]*pd1[2];
    double fzfz = pd1[2]*pd1[2];
    
    double fxx = pd2[0][0];
    double fxy = pd2[0][1];
    double fxz = pd2[0][2];
    double fyy = pd2[1][1];
    double fyz = pd2[1][2];
    double fzz = pd2[2][2];
    
    double g2 = pd1[0]*pd1[0]+pd1[1]*pd1[1]+pd1[2]*pd1[2];
    if(g2 == 0){
      Kmax = 0;
      Kmin = 0;
      
      Rmax=0;
      Rmin=0;
      
      Tmax[0] = Tmax[1] = Tmax[2] = 0;
      Tmin[0] = Tmin[1] = Tmin[2] = 0;
      
      return;
    }
    double g1 = sqrt(g2);
    double g3 = g2*g1;
    double g4 = g2*g2;
    
    /* mean and gaussian curvatures */
    double H = (fxx*(fyfy+fzfz) + fyy*(fxfx+fzfz) + fzz*(fxfx+fyfy)
                - 2*(fxy*fxfy+fxz*fxfz+fyz*fyfz))/2;
    H /= g3;
    
    double K = fxfx*(fyy*fzz-fyz*fyz)
      + fyfy*(fxx*fzz-fxz*fxz)
      + fzfz*(fxx*fyy-fxy*fxy)
      + 2*(fxfy*(fxz*fyz-fxy*fzz)
	   + fxfz*(fxy*fyz-fxz*fyy)
	   + fyfz*(fxy*fxz-fxx*fyz));
    
    K /= g4;
    
    /* principal curvatures */
    double discr = sqrt(fabs(H*H-K));
    
    Kmax = (float)(H + discr);
    Kmin = (float)(H - discr);
    
    Rmax=0;
    Rmin=0;
    
    Tmax[0] = Tmax[1] = Tmax[2] = 0;
    Tmin[0] = Tmin[1] = Tmin[2] = 0;
    
    double EPS = 0.0000001;
    if(discr > EPS){/* it is not an umbilic */
      
      /* matrix entries */
      double m11 = ((-1 + fxfx/g2)*fxx)/g1 + (fxfy*fxy)/g3 + (fxfz*fxz)/g3;
      double m12 = ((-1 + fxfx/g2)*fxy)/g1 + (fxfy*fyy)/g3 + (fxfz*fyz)/g3;
      double m13 = ((-1 + fxfx/g2)*fxz)/g1 + (fxfy*fyz)/g3 + (fxfz*fzz)/g3;
      double m21 = (fxfy*fxx)/g3 + ((-1 + fyfy/g2)*fxy)/g1 + (fyfz*fxz)/g3;
      double m22 = (fxfy*fxy)/g3 + ((-1 + fyfy/g2)*fyy)/g1 + (fyfz*fyz)/g3;
      double m23 = (fxfy*fxz)/g3 + ((-1 + fyfy/g2)*fyz)/g1 + (fyfz*fzz)/g3;
      double m31 = (fxfz*fxx)/g3 + (fyfz*fxy)/g3 + ((-1 + fzfz/g2)*fxz)/g1;
      double m32 = (fxfz*fxy)/g3 + (fyfz*fyy)/g3 + ((-1 + fzfz/g2)*fyz)/g1;
      double m33 = (fxfz*fxz)/g3 + (fyfz*fyz)/g3 + ((-1 + fzfz/g2)*fzz)/g1;
      
      /* solve for eigenvectors */
      double tmp1 = m11+Kmax;
      double tmp2 = m22+Kmax;
      double tmp3 = m33+Kmax;
      
      double ux[3], uy[3], uz[3], len[3];
      ux[0] = m12*m23-m13*tmp2;
      uy[0] = m13*m21-m23*tmp1;
      uz[0] = tmp1*tmp2-m12*m21;
      len[0] = sqrt(ux[0]*ux[0]+uy[0]*uy[0]+uz[0]*uz[0]);
      
      ux[1] = m12*tmp3-m13*m32;
      uy[1] = m13*m31-tmp1*tmp3;
      uz[1] = tmp1*m32-m12*m31;
      len[1] = sqrt(ux[1]*ux[1]+uy[1]*uy[1]+uz[1]*uz[1]);
      
      ux[2] = tmp2*tmp3-m23*m32;
      uy[2] = m23*m31-m21*tmp3;
      uz[2] = m21*m32-m31*tmp2;
      len[2] = sqrt(ux[2]*ux[2]+uy[2]*uy[2]+uz[2]*uz[2]);
      
      int index = 0;
      double max = len[0];
      if ( len[1] > max ) {
        index  = 1;
        max = len[1];
      }
      if ( len[2] > max ) {
        index = 2;
        max = len[2];
      }
      
      Tmax[0] = (float)(ux[index]/len[index]);
      Tmax[1] = (float)(uy[index]/len[index]);
      Tmax[2] = (float)(uz[index]/len[index]);
      
      /* second tangent is cross product of first tangent and normal */
      double N[3];
      int i;
      for(i=0; i<3; i++){
        N[i] = - pd1[i]/g1;
      }
      
      Tmin[0] = (float)(Tmax[1]*N[2]-Tmax[2]*N[1]);
      Tmin[1] = (float)(Tmax[2]*N[0]-Tmax[0]*N[2]);
      Tmin[2] = (float)(Tmax[0]*N[1]-Tmax[1]*N[0]);
      
      Rmax = 0;
      Rmin = 0;
      
      for(i=0; i<3; i++){
        for(int j=0; j<3; j++){
          Rmax += (float)(pd2[i][j]*Tmax[i]*N[j]);
          Rmin += (float)(pd2[i][j]*Tmin[i]*N[j]);
        }
      }
      
      Rmax *= (3*Kmax);
      Rmin *= (3*Kmin);
      
      for(i=0; i<3; i++)
        for(int j=0; j<3; j++)
          for(int k=0; k<3; k++){
            Rmax += (float)(pd3[i][j][k]*Tmax[i]*Tmax[j]*Tmax[k]);
            Rmin += (float)(pd3[i][j][k]*Tmin[i]*Tmin[j]*Tmin[k]);
          }
      
      Rmax /= (float)g1;
      Rmin /= (float)g1;
    }
  }
  
  inline bool ridgeTest(float t1[6], float t2[6],
                        float Rmax1, float Rmax2,
                        int axis){
    float kmax1 = t1[0];
    float kmax2 = t2[0];
    float kmin1 = t1[3];
    float kmin2 = t2[3];
    float tx1 = t1[1];
    float ty1 = t1[2];
    float tx2 = t2[1];
    float ty2 = t2[2];
    
    if(kmax1 < fabs(kmin1) || kmax2 < fabs(kmin2))
      return false;
    if(tx1*tx2 + ty1*ty2 < 0){
      tx2 = -tx2;
      ty2 = -ty2;
      Rmax2 = -Rmax2;
    }
    
    if(Rmax1*Rmax2 > 0)
      return false;
    
    if(axis == 0)
      return (Rmax1*tx1 > 0) || (Rmax2*tx2 < 0);
    else
      return (Rmax1*ty1 > 0) || (Rmax2*ty2 < 0);
  }
  
  inline bool valleyTest(float t1[6], float t2[6],
                         float Rmin1, float Rmin2,
                         int axis){
    float kmax1 = t1[0];
    float kmax2 = t2[0];
    float kmin1 = t1[3];
    float kmin2 = t2[3];
    float tx1 = t1[4];
    float ty1 = t1[5];
    float tx2 = t2[4];
    float ty2 = t2[5];
    
    if(kmin1 > -fabs(kmax1) || kmin2 > -fabs(kmax2))
      return false;
    if(tx1*tx2 + ty1*ty2 < 0){
      tx2 = -tx2;
      ty2 = -ty2;
      Rmin2 = -Rmin2;
    }
    
    if(Rmin1*Rmin2 > 0)
      return false;
    
    if(axis == 0)
      return (Rmin1*tx1 < 0) || (Rmin2*tx2 > 0);
    else
      return (Rmin1*ty1 < 0) || (Rmin2*ty2 > 0);
  }
  
  void thresholdingRV(float TR, float TV, bool bold){
    int i, j;
    for(i=0; i<height; i++)
      for(j=0; j<width; j++)
        totalW[i][j] = 0;
    
    int *stack = new int[height*width];
    int *list = new int[height*width];
    int top, listN;
    for(i=0; i<height; i++)
      for(j=0; j<width; j++){
        if(totalW[i][j] != 0)
          continue;
        
        if(value[i][j] == RIDGE){
          top = 1;
          listN = 0;
          stack[0] = i*width + j;
          totalW[i][j] = -1;
          float stren = 0;
          while(top > 0){
            int s = stack[--top];
            list[listN++] = s;
            int u = (int)(s/width);
            int v = s%width;
            int x, y;
            for(y=-1; y<2; y++)
              for(x=-1; x<2; x++){
                int x1 = v+x;
                int y1 = u+y;
                if((x == 0 && y == 0) ||
                   x1 < 0 || x1 >= width ||
                   y1 < 0 || y1 >= height)
                  continue;
                
                if(value[y1][x1] != RIDGE ||
                   totalW[y1][x1] != 0)
                  continue;
                
                stack[top++] = y1*width + x1;
                totalW[y1][x1] = -1;
                float d2;
                if(x == 0 || y == 0)
                  d2 = 1;
                else
                  d2 = 1.41421f;
                d2 += (depth[i][j]-depth[y1][x1])*(depth[i][j]-depth[y1][x1]);
                
                stren += (float)(sqrt(d2)*(second[i][j][1] + second[y1][x1][1]));
              }
          }
          if(0.5f*stren > TR){
            if(bold)
              for(int k=0; k<listN; k++){
                int l = list[k];
                int x = l%width;
                int y = (int)(l/width);
                value[y][x+1] = RIDGE;
                totalW[y][x+1] = -1;
                value[y+1][x] = RIDGE;
                totalW[y+1][x] = -1;
                value[y][x-1] = RIDGE;
                totalW[y][x-1] = -1;
                value[y-1][x] = RIDGE;
                totalW[y-1][x] = -1;
              }
          }
          else //if(0.5f*stren < TR)
            for(int k=0; k<listN; k++){
              int l = list[k];
              value[(int)(l/width)][l%width] = NORMAL;
            }
        }
        else if(value[i][j] == VALLEY){
          top = 1;
          listN = 0;
          stack[0] = i*width + j;
          totalW[i][j] = -1;
          float stren = 0;
          while(top > 0){
            int s = stack[--top];
            list[listN++] = s;
            int u = (int)(s/width);
            int v = s%width;
            int x, y;
            for(y=-1; y<2; y++)
              for(x=-1; x<2; x++){
                int x1 = v+x;
                int y1 = u+y;
                if((x == 0 && y == 0) ||
                   x1 < 0 || x1 >= width ||
                   y1 < 0 || y1 >= height)
                  continue;
                
                if(value[y1][x1] != VALLEY ||
                   totalW[y1][x1] != 0)
                  continue;
                
                stack[top++] = y1*width + x1;
                totalW[y1][x1] = -1;
                float d2;
                if(x == 0 || y == 0)
                  d2 = 1;
                else
                  d2 = 1.41421f;
                d2 += (depth[i][j]-depth[y1][x1])*(depth[i][j]-depth[y1][x1]);
                
                stren -= (float)(sqrt(d2)*(second[i][j][4] + second[y1][x1][4]));
              }
          }
          if(0.5f*stren > TV){
            if(bold)
              for(int k=0; k<listN; k++){
                int l = list[k];
                int x = l%width;
                int y = (int)(l/width);
                value[y][x+1] = VALLEY;
                totalW[y][x+1] = -1;
                value[y+1][x] = VALLEY;
                totalW[y+1][x] = -1;
                value[y][x-1] = VALLEY;
                totalW[y][x-1] = -1;
                value[y-1][x] = VALLEY;
                totalW[y-1][x] = -1;
              }
          }
          else //if(0.5f*stren < TV)
            for(int k=0; k<listN; k++){
              int l = list[k];
              value[(int)(l/width)][l%width] = NORMAL;
            }
        }
      }
    delete[] stack;
    delete[] list;
  }
  
  void thresholdingRV2(float strong, float week){
    int i, j;
    for(i=0; i<height; i++)
      for(j=0; j<width; j++)
        totalW[i][j] = 0;
    
    int *stack = new int[height*width];
    int *list = new int[height*width];
    int top, listN;
    for(i=0; i<height; i++)
      for(j=0; j<width; j++){
        if(totalW[i][j] != 0)
          continue;
        
        if(value[i][j] == RIDGE){
          if(second[i][j][1] < week)
            value[i][j] = NORMAL;
          else{
            top = 1;
            listN = 0;
            stack[0] = i*width + j;
            totalW[i][j] = -1;
            bool flag = false;
            while(top > 0){
              int s = stack[--top];
              list[listN++] = s;
              int u = (int)(s/width);
              int v = s%width;
              if(second[u][v][1] > strong)
                flag = true;
              int x, y;
              for(y=-1; y<2; y++)
                for(x=-1; x<2; x++){
                  int x1 = v+x;
                  int y1 = u+y;
                  if((x == 0 && y == 0) ||
                     x1 < 0 || x1 >= width ||
                     y1 < 0 || y1 >= height)
                    continue;
                  
                  if(value[y1][x1] != RIDGE ||
                     second[y1][x1][1] < week || 
                     totalW[y1][x1] != 0)
                    continue;
                  
                  stack[top++] = y1*width + x1;
                  totalW[y1][x1] = -1;
                }
            }
            if(flag)
              for(int k=0; k<listN; k++){
                int l = list[k];
                int x = l%width;
                int y = (int)(l/width);
                if(second[y][x][1] < strong)
                  continue;
                value[y][x+1] = RIDGE;
                totalW[y][x+1] = -1;
                value[y+1][x] = RIDGE;
                totalW[y+1][x] = -1;
                value[y][x-1] = RIDGE;
                totalW[y][x-1] = -1;
                value[y-1][x] = RIDGE;
                totalW[y-1][x] = -1;
              }
            else //if(0.5f*stren < TR)
              for(int k=0; k<listN; k++){
                int l = list[k];
                value[(int)(l/width)][l%width] = NORMAL;
              }
          }
        }
        else if(value[i][j] == VALLEY){
          if(-second[i][j][3] < week)
            value[i][j] = NORMAL;
          else{
            top = 1;
            listN = 0;
            stack[0] = i*width + j;
            totalW[i][j] = -1;
            bool flag = false;
            while(top > 0){
              int s = stack[--top];
              list[listN++] = s;
              int u = (int)(s/width);
              int v = s%width;
              if(-second[u][v][3] > strong)
                flag = true;
              int x, y;
              for(y=-1; y<2; y++)
                for(x=-1; x<2; x++){
                  int x1 = v+x;
                  int y1 = u+y;
                  if((x == 0 && y == 0) ||
                     x1 < 0 || x1 >= width ||
                     y1 < 0 || y1 >= height)
                    continue;
                  
                  if(value[y1][x1] != VALLEY ||
                     -second[y1][x1][3] < week || 
                     totalW[y1][x1] != 0)
                    continue;
                  
                  stack[top++] = y1*width + x1;
                  totalW[y1][x1] = -1;
                }
            }
            if(flag)
              for(int k=0; k<listN; k++){
                int l = list[k];
                int x = l%width;
                int y = (int)(l/width);
                if(-second[y][x][3] < strong)
                  continue;
                value[y][x+1] = VALLEY;
                totalW[y][x+1] = -1;
                value[y+1][x] = VALLEY;
                totalW[y+1][x] = -1;
                value[y][x-1] = VALLEY;
                totalW[y][x-1] = -1;
                value[y-1][x] = VALLEY;
                totalW[y-1][x] = -1;
              }
            else //if(0.5f*stren < TR)
              for(int k=0; k<listN; k++){
                int l = list[k];
                value[(int)(l/width)][l%width] = NORMAL;
              }
          }
        }
      }
    delete[] stack;
    delete[] list;
  }
  
  void inline computeFrontDepth(){
    visitID++;
    stack[0] = root;
    top = 0;
    root->visit = visitID;
    while(top >= 0){
      BasisFunction* bf = stack[top--];
      transCenterAndSupp(bft, bf);
      
      float r = bft->support;
      float cu = bft->centerX;
      float cv = bft->centerY;
      float cd = bft->centerZ - r;
      
      int sk = (int)(cu-r) + 1;
      if(sk < 0)
        sk = 0;
      int ek = (int)(cu+r) + 1;
      if(ek > width)
        ek = width;
      if(sk >= width || ek <= 0)
        continue;
      
      int sj = (int)(cv-r) + 1;
      if(sj < 0)
        sj = 0;
      int ej = (int)(cv+r) + 1;
      if(ej > height)
        ej = height;
      if(sj >= height || ej <= 0)
        continue;
      
      if(bf->leaf || r < FEW_PIX){
        transPoly(bft, bf);
        
        if(!both && bft->normalZAtC() > 0)
          continue;
        
        float r2 = r*r;
        for(int j = sj; j < ej; j++){
          float dj = r2 - (j-cv)*(j-cv);
          for(int k = sk; k < ek; k++){
            if((k-cu)*(k-cu) > dj)
              continue;
            
            float f = front[j][k];
            if(f < cd)
              continue;
            
            float qd;
            if(!bft->rayIntersectLAZ(qd, (float)k, (float)j))
              continue;
            
            if(!both){
              float n1[3];
              bft->gradientLA(n1, (float)k, (float)j, qd);
              if(n1[2] > 0)
                continue;
            }
            
            if(qd < f)
              front[j][k] = qd;
          }
        }
      }
      else{
        int N = bf->childN;
        BasisFunction** bfs = bf->child;
        for(int i=0; i<N; i++){
          BasisFunction* bfi = bfs[i];
          if(bfi->visit != visitID){
            bfi->visit = visitID;
            stack[++top] = bfs[i];
          }
        }
      }
    }
  }
  
  inline void normalizeNormals(){
    int i, j;
    for(i=0; i<height; i++){
      for(j=0; j<width; j++){
        float *n = normals[i][j];
        float len = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
        if(len != 0){
          len = (float)(1.0/sqrt(len));
          n[0] *= len;
          n[1] *= len;
          n[2] *= len;
        }
        else{
          n[0] = n[1] = n[2] = 0;
        }
      }
    }
  }
  
  void setTransform(float R1[9], float s, float t[3]){
    Rxx = Rt[0] = R[0] = R1[0];
    Rxy = Rt[1] = R[3] = R1[3];
    Rxz = Rt[2] = R[6] = R1[6];
    Ryx = Rt[3] = R[1] = R1[1];
    Ryy = Rt[4] = R[4] = R1[4];
    Ryz = Rt[5] = R[7] = R1[7];
    Rzx = Rt[6] = R[2] = R1[2];
    Rzy = Rt[7] = R[5] = R1[5];
    Rzz = Rt[8] = R[8] = R1[8];
    
    tx = t[0];
    ty = t[1];
    tz = t[2];
    
    this->s = s;
    is = 1.0f/s;
    is2 = 1.0f/(s*s);
    is3 = is*is2;
  }
  
  inline void transCenterAndSupp(BasisFunction* bf_out, BasisFunction* bf_in){
    float c[3] = {bf_in->centerX, bf_in->centerY, bf_in->centerZ};
    float Rc[3];
    MAT_BY_VEC(Rc, R, c);
    
    bf_out->centerX = s*Rc[0] + tx;
    bf_out->centerY = s*Rc[1] + ty;
    bf_out->centerZ = s*Rc[2] + tz;
    
    bf_out->support = s*bf_in->support;
    
    bf_out->leaf = bf_in->leaf;
  }
};

#endif
