#ifndef BFAXISKDTREE
#define BFAXISKDTREE 1

#include "BasisFunction.h"
#include <cstdio>
#include <iostream>
using namespace std;

#define SEARCH_MAX 20

class BfAxisKdTree;

class BfAxisKdCell {

public:

  BfAxisKdTree* tree;
  int start, end;

  BfAxisKdCell *fore, *back;

  char s_axis;
  float middle;

  bool leaf;

  BfAxisKdCell(BfAxisKdTree* tree, int s, int e){
    this->tree = tree;
    start = s;
    end = e;

    leaf = true;

    fore = NULL;
    back = NULL;

    split();
  }

  ~BfAxisKdCell(){
    if(!leaf){
      delete fore;
      delete back;
    }
  }

#if 1
  void split();
  void collectBfIndexInBf(BasisFunction* bf);
#endif

};

class BfAxisKdTree{
public:
  int bfN;
  BasisFunction** bfs;
  BfAxisKdCell* root;
  int *table;
  int listN, *list;

  BfAxisKdTree(BasisFunction** bfs, int bfN){
    this->bfs = bfs;
    int N = bfN;

    table = new int[N];
    for(int i=0; i<N; i++)
      table[i] = i;

    root = new BfAxisKdCell(this, 0, N);

    list = new int[N];
  }

  ~BfAxisKdTree(){
    delete root;
    delete[] list;
    delete[] table;
  }

  void collectBfIndexInBf(int *&list, int &listN, BasisFunction* bf){
    this->listN = 0;
    root->collectBfIndexInBf(bf);
    listN = this->listN;
    list = this->list;
  }

  void swapIndex(int i, int j){
    int tmp = table[i];
    table[i] = table[j];
    table[j] = tmp;
  }

  inline int getIndex(int i){
    return table[i];
  }
};

void BfAxisKdCell::split(){
  if(end-start < SEARCH_MAX)
    return;

  int i,j;
  BasisFunction** bfs = tree->bfs;

  float min[3], max[3];
  BasisFunction* bf = bfs[tree->getIndex(start)];
  max[0] = min[0] = bf->centerX;
  max[1] = min[1] = bf->centerY;
  max[2] = min[2] = bf->centerZ;
  for(i=start+1; i<end; i++){
    bf = bfs[tree->getIndex(i)];

    if(bf->centerX < min[0])      min[0] = bf->centerX;
    else if(bf->centerX > max[0]) max[0] = bf->centerX;

    if(bf->centerY < min[1])      min[1] = bf->centerY;
    else if(bf->centerY > max[1]) max[1] = bf->centerY;

    if(bf->centerZ < min[2])      min[2] = bf->centerZ;
    else if(bf->centerZ > max[2]) max[2] = bf->centerZ;
  }
  if(max[0] - min[0] > max[1] - min[1]){
    s_axis = 0;
    middle = 0.5f*(max[0] + min[0]);
  }
  else{
    s_axis = 1;
    middle = 0.5f*(max[1] + min[1]);
  }
  if(max[2] - min[2] > max[s_axis] - min[s_axis]){
    s_axis = 2;
    middle = 0.5f*(max[2] + min[2]);
  }

  i = start;
  j = end-1;
  while(i <= j){
    bf = bfs[tree->getIndex(i)];
    float pi;
    if(s_axis == 0)      pi = bf->centerX;
    else if(s_axis == 1) pi = bf->centerY;
    else                 pi = bf->centerZ;
    if(pi > middle)
      i++;
    else{
      tree->swapIndex(i, j);
      j--;
    }
  }

  if(start == i || i == end)
    return;

  leaf = false;
  fore = new BfAxisKdCell(tree, start, i);
  back = new BfAxisKdCell(tree, i, end);

  return;
}

void BfAxisKdCell::collectBfIndexInBf(BasisFunction* bf){
  float r = bf->support;
  float c[3] = {bf->centerX, bf->centerY, bf->centerZ};
  if(leaf){
    BasisFunction** bfs = tree->bfs;
    double r2 = r*r;
    for(int i=start; i<end; i++){
      int j = tree->getIndex(i);
      BasisFunction* bfi = bfs[j];
      float vx = c[0] - bfi->centerX;
      float vy = c[1] - bfi->centerY;
      float vz = c[2] - bfi->centerZ;
      if(vx*vx+vy*vy+vz*vz < r2)
        tree->list[tree->listN++] = j;
    }
  }
  else{
    float d = c[s_axis] - middle;
    if(d+r >= 0)
      fore->collectBfIndexInBf(bf);
    if(d-r <= 0)
      back->collectBfIndexInBf(bf);
  }
}

#endif // BFAXISKDTREE
