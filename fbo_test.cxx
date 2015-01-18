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

#if 0
#include "DebugGPU.hxx"

DebugGPU debugGPU;
#endif

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
                            "GL_NV_float_buffer " \
                            "GL_ARB_vertex_program " \
                            "GL_ARB_fragment_program " \
                            "GL_ARB_multitexture " \
                            "GL_ARB_point_sprite " \
                            "GL_EXT_framebuffer_object "
#endif

#include "RenderTextureFBO.h"

RenderTextureFBO* rt;
// GLenum texTarget = GL_TEXTURE_2D;
// GLenum texInternalFormat = GL_FLOAT_RGBA_NV;
GLenum texTarget = GL_TEXTURE_RECTANGLE_NV;
GLenum texInternalFormat = GL_RGBA;
GLenum texFormat = GL_RGBA;
GLenum texType = GL_FLOAT;

GLuint textureProgram = 0, renderProgram = 0;

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

void initFrameBufferObjects()
{
  rt = new RenderTextureFBO( width, height, texTarget, texInternalFormat, texFormat, texType, 1 );
  rt->Activate();
  pane.initGL();
  rt->Deactivate();

  // load fragment programs
  const char* strTextureProgram2D = 
    "!!ARBfp1.0\n"
    "TEX result.color, fragment.texcoord[0], texture[0], 2D;\n"
    "END\n";

  const char* strTextureProgramRECT = 
    "!!ARBfp1.0\n"
    "TEX result.color, fragment.texcoord[0], texture[0], RECT;\n"
    "END\n";

  glGenProgramsARB( 1, &textureProgram );
  glBindProgramARB( GL_FRAGMENT_PROGRAM_ARB, textureProgram );
  // load correct program based on texture target
  if ( texTarget == GL_TEXTURE_RECTANGLE_NV ) 
    {
      glProgramStringARB(GL_FRAGMENT_PROGRAM_ARB, GL_PROGRAM_FORMAT_ASCII_ARB,
			 (GLsizei)strlen(strTextureProgramRECT), strTextureProgramRECT);
    }
  else 
    {
      glProgramStringARB(GL_FRAGMENT_PROGRAM_ARB, GL_PROGRAM_FORMAT_ASCII_ARB,
			 (GLsizei)strlen(strTextureProgram2D), strTextureProgram2D);
    }
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
//   debugGPU.begin();

  rt->Activate();
  {
    pane.setBackgroundColor( .0f, .0f, .0f );
    pane.setIsGradientBackground( false );
    pane.clear( width, height );
    pane.setView();
    //   pane.setLight();

#if 0
    glBindProgramARB(GL_FRAGMENT_PROGRAM_ARB, renderProgram);
    glEnable(GL_FRAGMENT_PROGRAM_ARB);
#endif

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

    if ( render_mode == DRAW_SPRITE )
      {
	glEnable( GL_POINT_SPRITE_ARB );
	glTexEnvi( GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE );
  	glEnable( GL_VERTEX_PROGRAM_POINT_SIZE_ARB );
  
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

 	glDisable( GL_VERTEX_PROGRAM_POINT_SIZE_ARB );
	glDisable( GL_POINT_SPRITE_ARB );

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

// 	glColor3f( 1.0f, .0f, .0f );
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

		glColor4f( 1.0f, .0f, .0f, list[j]->support / FID_DIV );
		// 2D coordinate of a parameter texture
		//glColor4f( texcoordParam[k].x, texcoordParam[k].y, 1.0f, pSize / FID_DIV );

		glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
		glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
		glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
		glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );
	      }
	  }

	glEnd();

      }

    cgGLDisableProfile( fragProfile );
    cgGLDisableProfile( vertProfile );

    //       debugGPU.capture( true );

    pane.finish();

    //       debugGPU.end();

    //       debugGPU.renderTexture();
  }
  
  rt->Deactivate();

  pane.clear( width, height );
  pane.initView();

  glBindProgramARB( GL_FRAGMENT_PROGRAM_ARB, textureProgram );
  glEnable( GL_FRAGMENT_PROGRAM_ARB );

  rt->Bind();
  rt->drawQuad();

  glDisable(GL_FRAGMENT_PROGRAM_ARB);

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
      
      delete rt;
      glDeleteProgramsARB(1, &textureProgram);
#if 0
      glDeleteProgramsARB(1, &renderProgram);
#endif
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
  pane.setMagObject( 10.0f );

  initExtensions();
  initFrameBufferObjects();

//   debugGPU.init( width, height, DEBUG_POINTS );
//   debugGPU.init( width, height, DEBUG_TEXTURE );
//   debugGPU.Activate();
//   pane.initGL();
//   debugGPU.Deactivate();
//   debugGPU.setIsValid( false );
  
  initCg();

  pane.initGL();

  ::glutMainLoop();
  return 0;
}


////////////////////////////////////////////////////////////////////////////////////

