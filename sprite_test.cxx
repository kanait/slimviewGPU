////////////////////////////////////////////////////////////////////
//
// $Id: bspsgpu.cxx,v 1.10 2005/01/05 14:07:53 kanai Exp $
//
// Copyright (c) 2004-2005 by Keio Research Institute at SFC
// All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#include "envDep.h"

#include "GL/glew.h"
#include "GL/wglew.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <GL/glext.h>

#include <Cg/cg.h>
#include <Cg/cgGL.h>

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

////////////////////////////////////////////////////////////////////////////////

GLPanel pane;
int width = 512;
int height = 512;

////////////////////////////////////////////////////////////////////////////////

#include "DebugGPU.hxx"

DebugGPU debugGPU;

////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////

#if 0
#define GLH_EXT_SINGLE_FILE
#include "glh/glh_extensions.h"

#define EXTENSIONS_REQUIRED "GL_NV_texture_rectangle " \
                            "GL_NV_vertex_program " \
                            "GL_NV_vertex_program2 " \
                            "GL_NV_vertex_program3 " \
                            "GL_NV_fragment_program " \
                            "GL_NV_fragment_program2 " \
                            "GL_NV_float_buffer " \
                            "GL_ARB_fragment_program " \
                            "GL_ARB_multitexture " \
                            "GL_ARB_point_sprite " \
                            "WGL_ARB_render_texture "\
                            "WGL_ARB_pixel_format " \
                            "WGL_ARB_pbuffer " \
                            "WGL_ATI_pixel_format_float "

#endif

void initExtensions()
{
#if 0
  if ( !glh_init_extensions(EXTENSIONS_REQUIRED) ) 
    {
      cerr << "Unable to load the following extension(s): "
	   << glh_get_unsupported_extensions() << "\n\nExiting...\n" << endl;
      exit(-1);
    }
#endif

  GLenum err = glewInit();
  if( err != GLEW_OK ) cout << "Error: %s" << glewGetErrorString(err) << endl;
}

////////////////////////////////////////////////////////////////////////////////

CGcontext context;
CGprogram vertProgram;
CGprogram fragProgram;
CGprofile vertProfile;
CGprofile fragProfile;

static void handleCgError() 
{
  fprintf( stderr, "Cg error: %s\n", cgGetErrorString(cgGetError()) );
}

void initCg()
{
  cgSetErrorCallback( handleCgError );
  context = cgCreateContext();
  
  // profiles
  vertProfile = cgGLGetLatestProfile( CG_GL_VERTEX );
  cgGLSetOptimalOptions( vertProfile );
  fragProfile  = cgGLGetLatestProfile( CG_GL_FRAGMENT );
  cgGLSetOptimalOptions( fragProfile );

  // vertex program (1st-pass)
  vertProgram  = cgCreateProgramFromFile( context, CG_SOURCE, "sprite_vp.cg",
 					  vertProfile, NULL, NULL );
  cgGLLoadProgram( vertProgram );

  // fragment program (1st-pass)
  fragProgram = cgCreateProgramFromFile( context, CG_SOURCE, "sprite_fp.cg",
					 fragProfile, NULL, NULL );
  cgGLLoadProgram( fragProgram );
}

void destroyCg()
{
  cgDestroyProgram( vertProgram );
  cgDestroyProgram( fragProgram );
  cgDestroyContext( context );
}

////////////////////////////////////////////////////////////////////////////////

#include "mydefine.h"

#define DRAW_BILLBOARD 0
#define DRAW_SPRITE 1
int render_mode = DRAW_SPRITE;

void display()
{
  debugGPU.begin();

  pane.setBackgroundColor( .0f, .0f, .0f );
  pane.setIsGradientBackground( false );
  pane.clear( width, height );
  pane.setView();
  //   pane.setLight();

  if ( render_mode == DRAW_SPRITE )
    {

      cgGLEnableProfile( vertProfile );
      cgGLBindProgram( vertProgram );

      cgGLSetStateMatrixParameter( cgGetNamedParameter(vertProgram, "ModelViewProj"),
				   CG_GL_MODELVIEW_PROJECTION_MATRIX,
				   CG_GL_MATRIX_IDENTITY );

      GLint viewport[4];
      GLdouble projMatrix[16];
      ::glGetDoublev( GL_PROJECTION_MATRIX, projMatrix );
      ::glGetIntegerv( GL_VIEWPORT, viewport );

      //
      // see eq.5 in [Botsch04]
      //
      // projMatrix[5]: 2.0 * n / (t - b) ... 
      // viewport[3]:   h_vp
      float size_param[4];
      size_param[0] = (float) projMatrix[5];
      size_param[1] = (float) viewport[3];
      size_param[2] = .0f;
      size_param[3] = .0f;
      cgGLSetParameter4fv( cgGetNamedParameter(vertProgram, "size_param"), size_param );


      cgGLEnableProfile( fragProfile );
      cgGLBindProgram( fragProgram );

      cgGLSetStateMatrixParameter(cgGetNamedParameter(fragProgram, "ModelViewProjInv"),
				  CG_GL_MODELVIEW_PROJECTION_MATRIX,
				  CG_GL_MATRIX_INVERSE );
      

      glEnable( GL_POINT_SPRITE_ARB );
      glTexEnvi( GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE );
      glEnable( GL_VERTEX_PROGRAM_POINT_SIZE_ARB );
  
#if 0
      glBegin( GL_POINTS );

      glColor4f( 1.0f, .0f, .0f, 3.0f / FID_DIV );
      glVertex3f( 1.0f, .0f, .0f );
  
      glEnd();
#endif

#if 1
      glBegin( GL_POINTS );
      //   int k = 0;
      for( int i = 0; i < 2*listN; i++ )
	{
	  if( Ns[i] == 0 ) continue;
	  if ( i%2 != 0 ) continue;

	  BasisFunction** list = lists[i];
	  for( int j = 0; j < Ns[i]; j++ )
	    {
	      if ( !(list[j]->leaf) ) continue;

	      // 2D coordinate of a parameter texture
	      glColor4f( 1.0f, .0f, .0f, list[j]->support / FID_DIV );
	      // 	  ++k;

	      //
	      // (center_x, center_y, center_z, radius of support)
	      //
	      glVertex3f( list[j]->centerX, list[j]->centerY, list[j]->centerZ );
	    }
	}
      glEnd();
#endif

      glDisable( GL_VERTEX_PROGRAM_POINT_SIZE_ARB );
      glDisable( GL_POINT_SPRITE_ARB );
  
      cgGLDisableProfile( fragProfile );
      cgGLDisableProfile( vertProfile );

      debugGPU.capture( true );

      pane.finish();

      debugGPU.end();

//       debugGPU.renderTexture();

    }
  else // #else
    {

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
  
      glColor3f( 1.0f, .0f, .0f );
      glBegin( GL_QUADS );

      while( imagemaker.top >= 0 )
	{
	  BasisFunction* bf = imagemaker.stack[imagemaker.top--];
      
	  //       r = bf->support * size_param.x / hpos.w + size_param.y;

	  //       if ( bf->leaf || (r < 4) )
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

	      glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
	      glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
	      glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
	      glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );

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

      glEnd();

#endif

#if 1
      glColor3f( 1.0f, .0f, .0f );
      glBegin( GL_QUADS );
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

	      glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
	      glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
	      glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
	      glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );

#if 0
	      glColor3f( .0f, .0f, .0f );
	      glBegin( GL_LINE_LOOP );
	      glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
	      glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
	      glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
	      glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );
	      glEnd();
#endif
	    }
	}

      glEnd();

#endif
  
      pane.finish();

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

      // screen capture
      pane.setRecordPNG( true );
      //pane.setIsCaptureDepth( true );
      pane.setPNGFile( "save_screen.png" );
      
      break;

    case 'n':

      if ( render_mode == DRAW_BILLBOARD )     render_mode = DRAW_SPRITE;
      else if ( render_mode == DRAW_SPRITE ) render_mode = DRAW_BILLBOARD;

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
  ::glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
  ::glutCreateWindow( argv[0] );

  ::glutDisplayFunc( display );
  ::glutKeyboardFunc( keyboard );
  ::glutMouseFunc( mouse );
  ::glutMotionFunc( motion );
  
  pane.init( width, height );
  pane.setViewPoint( Point3f( .0f, .0f, 50.0f ) );
  pane.setViewVector( Vector3f( .0f, .0f, -50.0f ) ); 

  initExtensions();

  debugGPU.init( width, height, DEBUG_POINTS );
//   debugGPU.init( width, height, DEBUG_TEXTURE );
  debugGPU.Activate();
  pane.initGL();
  debugGPU.Deactivate();
  debugGPU.setIsValid( false );
  
  initCg();

  pane.initGL();

  ::glutMainLoop();
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////

