////////////////////////////////////////////////////////////////////
//
// $Id: bspsgpu.cxx,v 1.10 2005/01/05 14:07:53 kanai Exp $
//
// Copyright (c) 2004-2005 by Keio Research Institute at SFC
// All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#include "envDep.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
using namespace std;

#include <Point3.h>
#include <Vector3.h>
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

#include "MeshL.hxx"
#include "SMFLIO.hxx"

#include "GLPanel.hxx"
#include "PNGImage.hxx"
#include "GLMaterial.hxx"
#include "GLMeshL.hxx"

GLPanel pane;
int width = 512;
int height = 512;

std::vector<float> point_data;

void display()
{
  pane.setBackgroundColor( 1.0f, 1.0f, 1.0f );
  pane.setIsGradientBackground( true );
  pane.clear( width, height );
  pane.setView();
//   pane.setLight();
  
  glColor3f( 1.0f, .0f, .0f );
  glBegin( GL_POINTS );
  for ( int i = 0; i < point_data.size() / 3; ++i )
    {
      glVertex3f( point_data[3*i], point_data[3*i+1], point_data[3*i+2] );
    }
  glEnd();

  pane.finish();

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

    case 'p':

      // screen capture
      pane.setRecordPNG( true );
      //pane.setIsCaptureDepth( true );
      pane.setPNGFile( "save_screen.png" );
      
      break;
      
    default:
      break;
      
    }
  
  ::glutPostRedisplay();
}

int main( int argc, char **argv )
{
  char filename[BUFSIZ];
  if ( argc < 2 )
    {
      strcpy( filename, "tmp000.pnt" );
    }
  else
    {
      strcpy( filename, argv[1] );
    }

  std::ifstream ifs; ifs.open( filename );
  if ( !ifs.is_open() )
    {
      std::cerr << "Cannot open " << filename << std::endl;
      return false;
    }

  std::string cline;
  while ( getline(ifs, cline, '\n') )
    {
      std::istringstream isstr( cline );
      float x;
      isstr >> x; point_data.push_back( x );
      isstr >> x; point_data.push_back( x );
      isstr >> x; point_data.push_back( x );
    }
  ifs.close();
  
  ::glutInitWindowSize( width, height );
  ::glutInit( &argc, argv );
  //::glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_DEPTH | GLUT_STENCIL );
  ::glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_DEPTH );
  ::glutCreateWindow( argv[0] );

  ::glutDisplayFunc( display );
  ::glutKeyboardFunc( keyboard );
  ::glutMouseFunc( mouse );
  ::glutMotionFunc( motion );
  
  pane.init( width, height );
  //pane.setIsLightOn( 1, false );
  pane.setViewPoint( Point3f( .0f, .0f, 50.0f ) );
  pane.setViewVector( Vector3f( .0f, .0f, -50.0f ) ); 
  pane.initGL();

  ::glutMainLoop();
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////

