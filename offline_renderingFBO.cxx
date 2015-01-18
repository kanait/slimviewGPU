////////////////////////////////////////////////////////////////////
//
// $Id: $
//
// Copyright (c) 2004-2005 by Keio Research Institute at SFC
// All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#include "envDep.h"

#include "GL/glew.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <GL/glext.h>

#include "GLPanel.hxx"
#include "GLMaterial.hxx"
#include "PNGImage.hxx"

static bool isCapture = false;
static bool isCaptureDepth = false;

////////////////////////////////////////////////////////////////////////////////

GLPanel pane;
int width = 512;
int height = 512;
int capture_width = 4096;
int capture_height = 4096;
GLMaterial glm;

////////////////////////////////////////////////////////////////////////////////

#include "RenderTextureFBO.h"

RenderTextureFBO* rt;
GLenum texTarget = GL_TEXTURE_2D;
GLenum texInternalFormat = GL_RGB;
GLenum texFormat = GL_RGB;
GLenum texType = GL_UNSIGNED_BYTE;

void initExtensions()
{
  GLenum err = glewInit();
  if( err != GLEW_OK ) cout << "Error: %s" << glewGetErrorString(err) << endl;
}

void initFrameBufferObjects()
{
  GLsizei sizeW;
  glGetIntegerv( GL_MAX_RENDERBUFFER_SIZE_EXT, &sizeW );

  // cout << "max size: " << sizeW << endl;
  // capture_width = 4096;
  // capture_height = 4096;

//   rt = new RenderTextureFBO( capture_width, capture_height, texTarget, texInternalFormat, texFormat, texType, 1 );
//   rt->Activate();
//   pane.initGL();
//   rt->Deactivate();
}

////////////////////////////////////////////////////////////////////////////////

#include "mydefine.h"

void display()
{
  if ( isCapture == true ) rt->Activate();

  {
    if ( isCapture == true ) pane.clear( capture_width, capture_height );
    else pane.clear( width, height );

    pane.setView();
    pane.setLight();
    
    ::glShadeModel( GL_SMOOTH );
    ::glEnable( GL_LIGHTING );

    glm.bind();

    glutSolidTeapot( 1.0f );
    
    ::glDisable( GL_LIGHTING );

    pane.finish();
  }

  if ( isCapture == true ) 
    {
      cout << "aaa" << endl;
      rt->Read();
      cout << "aaa" << endl;
      PNGImage pi( capture_width, capture_height, isCaptureDepth ); 
      pi.capture_and_write( "screen.png" );
      cout << "aaa" << endl;
      rt->Deactivate();

      isCapture = false;
      isCaptureDepth = false;
      ::glutPostRedisplay();

      delete rt;
    }

  ::glutSwapBuffers();
  ::glutReportErrors();
}

////////////////////////////////////////////////////////////////////////////////////

void mouse( int button, int state, int x, int y )
{
  pane.setScreenXY( x, y );
  
  // Rotate
  if ( state == GLUT_DOWN )
    {
      if ( button == GLUT_LEFT_BUTTON )        pane.startRotate();
      else if ( button == GLUT_MIDDLE_BUTTON ) pane.startZoom();
      else if ( button == GLUT_RIGHT_BUTTON )  pane.startMove();
      ::glutPostRedisplay();
    }
  else if ( state == GLUT_UP )
    {
      pane.finishRMZ();
      ::glutPostRedisplay();
    }
}

void motion( int x, int y )
{
  if ( pane.isRotate() )
    {
      pane.updateRotate( x, y );
      ::glutPostRedisplay();
    }
  if ( pane.isMove() )
    {
      pane.updateMove( x, y );
      ::glutPostRedisplay();
    }
  if ( pane.isZoom() )
    {
      pane.updateZoom( x, y );
      ::glutPostRedisplay();
    }
}

void keyboard( unsigned char c, int x, int y )
{
  switch ( c )
    {
    case 'q':
      
      exit(0);
      break;

    case 'o':

      isCapture = true;
      isCaptureDepth = true;
      texInternalFormat = GL_RGBA;
      texFormat = GL_RGBA;
      rt = new RenderTextureFBO( capture_width, capture_height, texTarget, texInternalFormat, texFormat, texType, 1 );

      break;

    case 'p':

      cout << "bbb" << endl;
      isCapture = true;
      texInternalFormat = GL_RGB;
      texFormat = GL_RGB;
      rt = new RenderTextureFBO( capture_width, capture_height, texTarget, texInternalFormat, texFormat, texType, 1 );
      cout << "bbb" << endl;

      break;

    default:
      break;
      
    }
  
  ::glutPostRedisplay();
}

int main( int argc, char* argv[] )
{
  ::glutInitWindowSize( width, height );
  ::glutInit( &argc, argv );
  ::glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
  ::glutCreateWindow( argv[0] );

  ::glutDisplayFunc( display );
  ::glutKeyboardFunc( keyboard );
  ::glutMouseFunc( mouse );
  ::glutMotionFunc( motion );
  
  pane.init( width, height );
  pane.initGL();
  pane.initGLEW();
  
  // pane.setIsLightOn( 1, false );
  pane.setViewPoint( Point3f( .0f, .0f, 8.0f ) );
  pane.setViewVector( Vector3f( .0f, .0f, -8.0f ) ); 

  cout << "aa" << endl;
  // initExtensions();
  cout << "aa" << endl;
  initFrameBufferObjects();

  // pane.initGL();

  ::glutMainLoop();
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////

