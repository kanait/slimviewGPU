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

#include "GLPanel.hxx"
#include "PNGImage.hxx"
#include "GLMaterial.hxx"

GLPanel pane;
int width = 512;
int height = 512;
bool isCapture = false;

//std::vector<float> point_data;

#include "FileManager.h"
#include "BasisFunction.h"
#include "ImageMaker.h"

FileManager filemanager;
BasisFunction*** lists;
ImageMaker imagemaker( width, height );
int* Ns;
int deg;
bool oriented;
int listN;

void display()
{
  pane.setBackgroundColor( 1.0f, 1.0f, 1.0f );
  pane.setIsGradientBackground( false );
  pane.clear( width, height );
  pane.setView();
//   pane.setLight();
  
  float mat[16];
  glGetFloatv( GL_MODELVIEW_MATRIX, mat );

  Point3f vRight( mat[0], mat[4], mat[8] );
  Point3f vUp( mat[1], mat[5], mat[9] );

  //
  // draw quads
  //

#if 0
  imagemaker.visitID++;
  imagemaker.stack[0] = imagemaker.root;
  imagemaker.top = 0;
  imagemaker.root->visit = imagemaker.visitID;
  
  while( imagemaker.top >= 0 )
    {
      BasisFunction* bf = imagemaker.stack[imagemaker.top--];

      if ( bf->leaf ) 
	{
	  
//  	  if(!(imagemaker.both) && (bf->normalZAtC() > 0) )
//  	    continue;

	  Point3f vCenter( bf->centerX, bf->centerY, bf->centerZ );
	  float pSize = bf->support;


// 	  if ( pSize > 1.0 ) continue;

	  Point3f vPoint0( vCenter + ((-1.0f*vRight - vUp) * pSize) );
	  Point3f vPoint1( vCenter + (( vRight - vUp) * pSize) );
	  Point3f vPoint2( vCenter + (( vRight + vUp) * pSize) );
	  Point3f vPoint3( vCenter + ((-1.0f*vRight + vUp) * pSize) );

	  // 2D coordinate of a parameter texture
 	  //glColor4f( texcoordParam[k].x, texcoordParam[k].y, 1.0f, pSize / FID_DIV );

	  glColor3f( .6f, .6f, .6f );
	  glBegin( GL_QUADS );
	  glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
	  glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
	  glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
	  glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );
	  glEnd();

	  glColor3f( .0f, .0f, .0f );
	  glBegin( GL_LINE_LOOP );
	  glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
	  glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
	  glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
	  glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );
	  glEnd();
	}
      else
	{
	  int N = bf->childN;
	  BasisFunction** bfs = bf->child;
	  for( int i=0; i<N; i++ )
	    {
	      BasisFunction* bfi = bfs[i];
	      if( bfi->visit != imagemaker.visitID )
		{
		  bfi->visit = imagemaker.visitID;
		  imagemaker.stack[++(imagemaker.top)] = bfs[i];
		}
	    }
	}
    }
#endif

#if 1
  for( int i = 0; i < 2*listN; i++ )
    {
      if( Ns[i] == 0 ) continue;
      if ( i%2 != 0 ) continue;

      BasisFunction** list = lists[i];
      for( int j = 0; j < Ns[i]; j++ )
	{
	  
	  if ( !(list[j]->leaf) ) continue;

	  Point3f vCenter( list[j]->centerX, list[j]->centerY, list[j]->centerZ );
	  float pSize = list[j]->support;

// 	  if ( pSize > 1.0 ) continue;

	  Point3f vPoint0( vCenter + ((-1.0f*vRight - vUp) * pSize) );
	  Point3f vPoint1( vCenter + (( vRight - vUp) * pSize) );
	  Point3f vPoint2( vCenter + (( vRight + vUp) * pSize) );
	  Point3f vPoint3( vCenter + ((-1.0f*vRight + vUp) * pSize) );

	  // 2D coordinate of a parameter texture
 	  //glColor4f( texcoordParam[k].x, texcoordParam[k].y, 1.0f, pSize / FID_DIV );

	  glColor3f( .85f, .85f, .85f );
	  glBegin( GL_QUADS );
	  glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
	  glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
	  glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
	  glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );
	  glEnd();

	  glColor3f( .0f, .0f, .0f );
	  glBegin( GL_LINE_LOOP );
	  glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
	  glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
	  glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
	  glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );
	  glEnd();
	}
    }
#endif
  

  pane.finish();

  if ( isCapture )
    {
      PNGImage pi( width, height, false ); 
      pi.capture_and_write( "save_screen.png" );
      isCapture = false;
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

    case 'p':

      isCapture = true;
//       // screen capture
//       pane.setRecordPNG( true );
//       //pane.setIsCaptureDepth( true );
//       pane.setPNGFile( "save_screen.png" );
      
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
      strcpy( filename, "moai_0.02.slim2" );
    }
  else
    {
      strcpy( filename, argv[1] );
    }

  // initialize SLIM
  listN = filemanager.readBallFile( filename, lists, Ns, deg, oriented );
  imagemaker.setBfLists( deg, listN, lists, Ns );
  imagemaker.setBoth( !oriented );
  
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
  pane.setMagObject( 10.0f );

  ::glutMainLoop();
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////

