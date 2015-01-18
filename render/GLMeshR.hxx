////////////////////////////////////////////////////////////////////
//
// $Id: $
//
//   OpenGL MeshR draw class
//
// Copyright (c) 2004-2010 by Takashi Kanai. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _GLMESHR_HXX
#define _GLMESHR_HXX 1

#include "envDep.h"
#if defined(_WINDOWS)
#include "stdafx.h"
#endif
#include "mydef.h"

#include <GL/gl.h>

#include "MeshR.hxx"
#include "GLMesh.hxx"

class GLMeshR : public GLMesh {

public:

  GLMeshR() : mesh_(NULL) {};
  ~GLMeshR() {};

  void clear() { if ( mesh_ ) deleteMesh(); };
  void setMesh( MeshR& mesh ) { mesh_ = &mesh; };
  void setMeshP( MeshR* mesh ) { mesh_ = mesh; };
  MeshR& mesh() { return *mesh_; };
  MeshR* meshP() { return mesh_; };

  void deleteMesh() { delete mesh_; mesh_ = NULL; };
  bool empty() const { 
    if ( mesh_ == NULL ) return true;
    return ( mesh_->numPoints() ) ? false : true; 
  };

  void draw() {
    if ( empty() ) return;

    ::glPushMatrix();
    if ( isDrawShading_ ) { drawShading(); }
    if ( isDrawWireframe_ ) { drawWireframe(); };
    if ( isDrawPoint_ ) { drawPoint(); };
    ::glPopMatrix();
  };
  //void draw() { draw( mesh_ ); };

  void drawShading() {
    if ( empty() ) return;

    if ( isSmoothShading_ )
      ::glShadeModel( GL_SMOOTH );
    else
      ::glShadeModel( GL_FLAT );

    ::glEnable( GL_LIGHTING );

    if ( isSmoothShading_ )
      {
	::glNormalPointer( GL_FLOAT, 0, 
			   (GLfloat *) &(mesh().normals()[0]) );
	::glEnableClientState( GL_NORMAL_ARRAY );
      }
    else
      {
	::glDisableClientState( GL_NORMAL_ARRAY );
      }

    if ( isDrawTexture_ )
      {
	//       ::glEnable( GL_TEXTURE_2D );
	::glTexCoordPointer( mesh().n_tex(), GL_FLOAT, 0, 
			     (GLfloat *) &(mesh().texcoords()[0]) );
	::glEnableClientState( GL_TEXTURE_COORD_ARRAY );
      }
    else
      {
	::glDisableClientState( GL_TEXTURE_COORD_ARRAY );
      }

    if ( isSmoothShading_ )
      {
	::glVertexPointer( 3, GL_FLOAT, 0, 
			   (GLfloat*) &(mesh().points()[0]) );
	::glEnableClientState( GL_VERTEX_ARRAY );
      }
    else
      {
	::glDisableClientState( GL_VERTEX_ARRAY );
      }

    mtl_.bind();

    std::vector<float>& points = mesh().points();
    std::vector<float>& normals = mesh().normals();
    std::vector<float>& fnormals = mesh().fnormals();
    std::vector<unsigned int>& indices = mesh().indices();
    std::vector<unsigned int>& nindices = mesh().nindices();

    if ( isSmoothShading_ )
      {
	::glDrawElements( GL_TRIANGLES, mesh().indices_size(), 
			  GL_UNSIGNED_INT, &(mesh().indices()[0]) );
      }
    else
      {
#if 0 // for large model
	::glBegin( GL_TRIANGLES );
	unsigned int n_id = 0;
	for ( int i = 0; i < mesh().numFaces(); ++i )
	  {
	    unsigned int id0 = indices[TRIANGLE * i];
	    unsigned int id1 = indices[TRIANGLE * i + 1];
	    unsigned int id2 = indices[TRIANGLE * i + 2];
	    ::glNormal3fv( &(fnormals[ nXYZ * n_id++ ]) );
	    ::glVertex3fv( &(points[ nXYZ * id0 ]) );
	    ::glVertex3fv( &(points[ nXYZ * id1 ]) );
	    ::glVertex3fv( &(points[ nXYZ * id2 ]) );
	  }
	::glEnd();
#endif
#if 1 // for fandisk
	::glBegin( GL_TRIANGLES );
	unsigned int n_id = 0;
	for ( int i = 0; i < mesh().numFaces(); ++i )
	  {
	    unsigned int id0 = indices[TRIANGLE * i];
	    unsigned int id1 = indices[TRIANGLE * i + 1];
	    unsigned int id2 = indices[TRIANGLE * i + 2];
	    ::glNormal3fv( &(normals[ nXYZ * n_id++ ]) );
	    ::glVertex3fv( &(points[ nXYZ * id0 ]) );
	    ::glNormal3fv( &(normals[ nXYZ * n_id++ ]) );
	    ::glVertex3fv( &(points[ nXYZ * id1 ]) );
	    ::glNormal3fv( &(normals[ nXYZ * n_id++ ]) );
	    ::glVertex3fv( &(points[ nXYZ * id2 ]) );
	  }
	::glEnd();
#endif
      }

    if ( isDrawTexture_ )
      {
	::glDisableClientState( GL_TEXTURE_COORD_ARRAY );
      }
  
    if ( isSmoothShading_ )
      {
	::glDisableClientState( GL_NORMAL_ARRAY );
	::glDisableClientState( GL_VERTEX_ARRAY );
      }
  

    ::glDisable( GL_LIGHTING );
  };

  void drawWireframe() {
    if ( empty() ) return;

    ::glColor3fv( wireColor() );
    ::glLineWidth( wireSize() );

    ::glEnableClientState( GL_VERTEX_ARRAY );
    ::glVertexPointer( 3, GL_FLOAT, 0, 
		       (GLfloat*) &(mesh().points()[0]) );

    for ( int i = 0; i < mesh().numFaces(); ++i )
      {
	::glDrawElements( GL_LINE_LOOP, 3, GL_UNSIGNED_INT, 
			  &(mesh().indices()[3*i]) );
      }

    ::glDisableClientState( GL_VERTEX_ARRAY );
  };

  void drawPoint() {
    if ( empty() ) return;

    ::glDisable( GL_LIGHTING );

    ::glColor3f( .0f, .0f, .0f );

    ::glEnableClientState( GL_VERTEX_ARRAY );
    ::glVertexPointer( 3, GL_FLOAT, 0, 
		       (GLfloat*) &(mesh().points()[0]) );

    ::glDrawArrays( GL_POINTS, 0, mesh().numPoints() );

    ::glDisableClientState( GL_VERTEX_ARRAY );
  };
  //   void drawColor();

private:

  MeshR* mesh_;
  
};

#endif // _GLMESHR_HXX
