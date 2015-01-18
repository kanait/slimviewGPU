////////////////////////////////////////////////////////////////////
//
// $Id: $
//
// Copyright (c) 2005 by Takashi Kanai
// All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _VIEWFRUSTUM_HXX
#define _VIEWFRUSTUM_HXX 1

class ViewFrustum {

public:

  ViewFrustum(){};
  ~ViewFrustum(){};
  
  void update( float modl[16], float proj[16] ) {
    
    float   clip[16];
    float   t;

    // Combine The Two Matrices (Multiply Projection By Modelview) 
    // But Keep In Mind This Function Will Only Work If You Do NOT
    // Rotate Or Translate Your Projection Matrix
    clip[ 0] = modl[ 0] * proj[ 0];
    clip[ 1] = modl[ 1] * proj[ 5];
    clip[ 2] = modl[ 2] * proj[10] + modl[ 3] * proj[14];
    clip[ 3] = modl[ 2] * proj[11];

    clip[ 4] = modl[ 4] * proj[ 0];
    clip[ 5] = modl[ 5] * proj[ 5];
    clip[ 6] = modl[ 6] * proj[10] + modl[ 7] * proj[14];
    clip[ 7] = modl[ 6] * proj[11];

    clip[ 8] = modl[ 8] * proj[ 0];
    clip[ 9] = modl[ 9] * proj[ 5];
    clip[10] = modl[10] * proj[10] + modl[11] * proj[14];
    clip[11] = modl[10] * proj[11];

    clip[12] = modl[12] * proj[ 0];
    clip[13] = modl[13] * proj[ 5];
    clip[14] = modl[14] * proj[10] + modl[15] * proj[14];
    clip[15] = modl[14] * proj[11];

    // Extract The Numbers For The RIGHT Plane
    frustum_[0][0] = clip[ 3] - clip[ 0];
    frustum_[0][1] = clip[ 7] - clip[ 4];
    frustum_[0][2] = clip[11] - clip[ 8];
    frustum_[0][3] = clip[15] - clip[12];
  
    // Normalize The Result
    t = sqrt( frustum_[0][0] * frustum_[0][0] + frustum_[0][1] * frustum_[0][1] + frustum_[0][2] * frustum_[0][2] );
    frustum_[0][0] /= t;
    frustum_[0][1] /= t;
    frustum_[0][2] /= t;
    frustum_[0][3] /= t;

    // Extract The Numbers For The LEFT Plane
    frustum_[1][0] = clip[ 3] + clip[ 0];
    frustum_[1][1] = clip[ 7] + clip[ 4];
    frustum_[1][2] = clip[11] + clip[ 8];
    frustum_[1][3] = clip[15] + clip[12];

    // Normalize The Result
    t = sqrt( frustum_[1][0] * frustum_[1][0] + frustum_[1][1] * frustum_[1][1] + frustum_[1][2] * frustum_[1][2] );
    frustum_[1][0] /= t;
    frustum_[1][1] /= t;
    frustum_[1][2] /= t;
    frustum_[1][3] /= t;

    // Extract The BOTTOM Plane
    frustum_[2][0] = clip[ 3] + clip[ 1];
    frustum_[2][1] = clip[ 7] + clip[ 5];
    frustum_[2][2] = clip[11] + clip[ 9];
    frustum_[2][3] = clip[15] + clip[13];

    // Normalize The Result
    t = sqrt( frustum_[2][0] * frustum_[2][0] + frustum_[2][1] * frustum_[2][1] + frustum_[2][2] * frustum_[2][2] );
    frustum_[2][0] /= t;
    frustum_[2][1] /= t;
    frustum_[2][2] /= t;
    frustum_[2][3] /= t;

    // Extract The TOP Plane
    frustum_[3][0] = clip[ 3] - clip[ 1];
    frustum_[3][1] = clip[ 7] - clip[ 5];
    frustum_[3][2] = clip[11] - clip[ 9];
    frustum_[3][3] = clip[15] - clip[13];

    // Normalize The Result
    t = sqrt( frustum_[3][0] * frustum_[3][0] + frustum_[3][1] * frustum_[3][1] + frustum_[3][2] * frustum_[3][2] );
    frustum_[3][0] /= t;
    frustum_[3][1] /= t;
    frustum_[3][2] /= t;
    frustum_[3][3] /= t;

    // Extract The FAR Plane
    frustum_[4][0] = clip[ 3] - clip[ 2];
    frustum_[4][1] = clip[ 7] - clip[ 6];
    frustum_[4][2] = clip[11] - clip[10];
    frustum_[4][3] = clip[15] - clip[14];

    // Normalize The Result
    t = sqrt( frustum_[4][0] * frustum_[4][0] + frustum_[4][1] * frustum_[4][1] + frustum_[4][2] * frustum_[4][2] );
    frustum_[4][0] /= t;
    frustum_[4][1] /= t;
    frustum_[4][2] /= t;
    frustum_[4][3] /= t;

    // Extract The NEAR Plane
    frustum_[5][0] = clip[ 3] + clip[ 2];
    frustum_[5][1] = clip[ 7] + clip[ 6];
    frustum_[5][2] = clip[11] + clip[10];
    frustum_[5][3] = clip[15] + clip[14];

    // Normalize The Result
    t = sqrt( frustum_[5][0] * frustum_[5][0] + frustum_[5][1] * frustum_[5][1] + frustum_[5][2] * frustum_[5][2] );
    frustum_[5][0] /= t;
    frustum_[5][1] /= t;
    frustum_[5][2] /= t;
    frustum_[5][3] /= t;
  }

  bool sphereInFrustum( float x, float y, float z, float radius )// Bude koule vid?t na scen??
  {
    for ( int i = 0; i < 6; ++i )
      {
	if( frustum_[i][0] * x + frustum_[i][1] * y + frustum_[i][2] * z + frustum_[i][3] <= -radius )
	  {
	    return false;
	  }
      }
    return true;
  };

  bool pointInFrustum( float x, float y, float z )// Bude koule vid?t na scen??
  {
    for ( int i = 0; i < 6; ++i )
      {
	if( frustum_[i][0] * x + frustum_[i][1] * y + frustum_[i][2] * z + frustum_[i][3] <= 0.0f )
	  {
	    return false;
	  }
      }
    return true;
  };

private:

  float frustum_[6][4];
  
};

#endif // _VIEWFRUSTUM_HXX
