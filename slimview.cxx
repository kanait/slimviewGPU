////////////////////////////////////////////////////////////////////
//
// $Id$
// $Date$
// $Author$
// $Rev$
// $URL$
//
// Copyright (c) 2005-2012 by Takashi Kanai. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#include "envDep.h"

#include "GL/glew.h"
#include "GL/wglew.h"

#include <GL/glut.h>


// #include <GL/gl.h>
// #include <GL/glu.h>
// #include <GL/glext.h>

#include <Cg/cg.h>
#include <Cg/cgGL.h>

#include <vector>
#include <list>
using namespace std;

#include <Point2.h>
#include <Point3.h>
#include <Vector3.h>
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

#include "mydef.h"
#include "nvtimer.h"

#ifndef RENDER_FBO

#include "RenderTexture.h"

#else // RENDER_FBO

#include "RenderTextureFBO.h"

#endif // RENDER_FBO

#include "GLPanel.hxx"
#include "VWIO.hxx"
#include "PNGImage.hxx"
#include "GLMaterial.hxx"

#include "MeshR.hxx"
#include "SMFRIO.hxx"
#include "GLMeshR.hxx"

#include "mydefine.h"

// old slim format
#ifdef OLD_SLIM
#include "FileManager.h"
#include "BasisFunction.h"
#include "ImageMaker.h"
#endif // OLD_SLIM

// new slim format
#ifdef NEW_SLIM
#include "SlimTree.hxx"
#include "SlimBall.hxx"
#include "BasisFunction.hxx"
#include "SlimTreeIO.hxx"
#endif // NEW_SLIM

////////////////////////////////////////////////////////////////////////////////////

int width = 512;
int height = 512;
//int width = 800;
//int height = 800;
// int width = 1024;
// int height = 1024;
// int width = 1280;
// int height = 1024;

float bgColor[4] = { 1.0f, 1.0f, 1.0f, 1.0f };

static float myLight[] = {
  0.0f, 0.0f, 100.0f, 1.0f,
  0.2f, 0.2f, 0.2f, 1.0f,
  1.0f, 1.0f, 1.0f, 1.0f,
  0.8f, 0.8f, 0.8f, 1.0f,
  1.0f, // 0.1f,
  0.0f,  // 0.05f
  0.0f  // 0.05f
};

static float myMatl[] = {
  0.2f, 0.2f, 0.2f, 1.0f, 
  0.8f, 0.8f, 0.8f, 1.0f, 
  0.0f, 0.0f, 0.0f, 1.0f, 
  1.0f, 1.0f, 1.0f, 1.0f, 
  80.0f
};

static float myBallLight[] = {
  0.0f, 0.0f, 100.0f, 1.0f,
  0.2f, 0.2f, 0.2f, 1.0f,
  1.0f, 1.0f, 1.0f, 1.0f,
  0.8f, 0.8f, 0.8f, 1.0f,
  1.0f, // 0.1f,
  0.0f,  // 0.05f
  0.0f  // 0.05f
};

static float myBallMatl[] = {
  0.2f, 0.2f, 0.2f, 1.0f, 
  0.9f, 0.5f, 0.5f, 1.0f, 
  0.0f, 0.0f, 0.0f, 1.0f, 
  1.0f, 1.0f, 1.0f, 1.0f, 
  80.0f
};

GLPanel pane;

MeshR meshR;
GLMeshR glmeshR;
GLMaterial ballMtl( myBallMatl );

////////////////////////////////////////////////////////////////////////////////////

timer fps(10);
float max_fps = .0f;
char buf[BUFSIZ];
char txt[BUFSIZ];
int slimball_count = 0;
int slimball_level = 10;
double slimball_error = 1.0e-05;

bool b[256];

#define RENDER_GPU  0
#define RENDER_BALL 1
#define RENDER_NONE 2
int render_mode = RENDER_GPU;

#define LOD_NONE  0
#define LOD_LEVEL 1
#define LOD_ERROR 2
int lod_mode = LOD_NONE;
#ifdef OLD_SLIM
bool isLOD = false;
#endif

//int draw_mode = DRAW_BILLBOARD;
int draw_mode = DRAW_POINTSPRITE;

// TEMP
PNGImage pi;

////////////////////////////////////////////////////////////////////////////////////

// const int wgl_buffers[] = {
//   //WGL_FRONT_LEFT_ARB,
//   WGL_AUX0_ARB,
//   WGL_AUX1_ARB
// };

// const GLenum colorAttachment_[] = {
//   GL_COLOR_ATTACHMENT0_EXT,
//   GL_COLOR_ATTACHMENT1_EXT,
// };

////////////////////////////////////////////////////////////////////////////////////
//
// pbuffer
//
int pbuffer_width  = width;
int pbuffer_height = height;

#ifndef RENDER_FBO

const GLenum buffers[] = {
  GL_AUX0,
  GL_AUX1,
};

RenderTexture* pbuffer = NULL;
RenderTexture* pbuffer16 = NULL;

#else // RENDER_FBO

const GLenum buffers[] = {
  GL_COLOR_ATTACHMENT0_EXT,
  GL_COLOR_ATTACHMENT1_EXT
};

RenderTextureFBO* pbuffer = NULL;

GLenum texTarget0 = GL_TEXTURE_RECTANGLE_EXT;
GLenum texInternalFormat0 = GL_RGBA32F_ARB;
GLenum texFormat0 = GL_RGBA;
GLenum texType0 = GL_FLOAT;

RenderTextureFBO* pbuffer16 = NULL;

GLenum texTarget1 = GL_TEXTURE_RECTANGLE_EXT;
//GLenum texInternalFormat1 = GL_RGBA32F_ARB;
GLenum texInternalFormat1 = GL_RGBA16F_ARB;
GLenum texFormat1 = GL_RGBA;
GLenum texType1 = GL_FLOAT;

#endif // RENDER_FBO

//
// for screen capture
//
int cap_width  = width;
int cap_height = height;

#ifndef RENDER_FBO

RenderTexture* cap_buffer = NULL;

#else //RENDER_FBO

RenderTextureFBO* cap_buffer = NULL;

GLenum texTarget2 = GL_TEXTURE_2D;
GLenum texInternalFormat2 = GL_RGB;
GLenum texFormat2 = GL_RGB;
GLenum texType2 = GL_UNSIGNED_BYTE;

#endif //RENDER_FBO

bool isCapture = false;

////////////////////////////////////////////////////////////////////////////////////
//
// Ball Color
//

#include "ColorConv.hxx"

static float max_dis = 1.0;
static float min_dis = 1.0e-10;
static Point3f max_col( 1.0f, 0.0f, 0.0f );
static Point3f min_col( 0.0f, 0.0f, 1.0f );

void convColor( float dis, float* color )
{
  ColorConv cc;
  Point3f min_col_hsv;
  Point3f max_col_hsv;
  cc.rgbtohsv( max_col, max_col_hsv );
  cc.rgbtohsv( min_col, min_col_hsv );

  float max_dis_log = (float) std::log10( max_dis );
  float min_dis_log = (float) std::log10( min_dis );
  
  float dis_log = (float) std::log10( dis );

  float t;
  if ( dis > max_dis ) 
    {
      t = 1.0f;
    }
  else if ( dis < min_dis )
    {
      t = 0.0f;
    }
  else
    {
      t = (dis_log - min_dis_log) / (max_dis_log - min_dis_log);
    }

  Point3f col_hsv; col_hsv.interpolate( min_col_hsv, max_col_hsv, t );
  Point3f col; 
  cc.hsvtorgb( col_hsv, col );
  color[0] = col.x;
  color[1] = col.y;
  color[2] = col.z;
}


////////////////////////////////////////////////////////////////////////////////////
//
// Cg resources
//

CGcontext context;

//
// vertex program
//
CGprogram vert1Program;
CGprogram vert2Program;
CGprofile vertProfile;

// fragment program
CGprogram frag1Program; // 1st
CGprogram frag2Program; // 2nd
CGprogram frag3Program; // 3rd
// CGprogram debugProgram; // for debug
CGprofile fragProfile;

// Material properties
CGparameter g_CGparam_DiffMat;
CGparameter g_CGparam_SpecMat;
CGparameter g_CGparam_AmbMat;
CGparameter g_CGparam_ShineMat;
CGparameter g_CGparam_EmisMat;

// Lighting properties
CGparameter g_CGparam_LightVec;
CGparameter g_CGparam_LightPos;
CGparameter g_CGparam_DiffLight;
CGparameter g_CGparam_SpecLight;
CGparameter g_CGparam_AmbLight;
CGparameter g_CGparam_AttenLight;
CGparameter g_CGparam_SpotLight;

////////////////////////////////////////////////////////////////////////////////////

// texture coordinate
std::vector<Point2f> texcoordParam;


#ifndef RENDER_FBO
// floating-point texture for pbuffer
GLuint qhatT;     // 1st pass
GLuint positionT; // 2nd pass
GLuint normalT;   // 2nd pass
#endif // RENDER_FBO

// tmp
static int count_frame = 0;

// #define OGL_QUERY 1

#ifdef OGL_QUERY
// occlusion query
GLuint oq_slim1;
GLuint oq_slim2;
#endif

////////////////////////////////////////////////////////////////////////////////////

// old slim
#ifdef OLD_SLIM
FileManager filemanager;
BasisFunction*** lists;
ImageMaker imagemaker( pbuffer_width, pbuffer_height );
int* Ns;
int deg;
bool oriented;
int listN;
#endif // OLD_SLIM

// new slim
#ifdef NEW_SLIM
SlimTreef slimtree;
#endif // NEW_SLIM

////////////////////////////////////////////////////////////////////////////////////

void checkGLErrors(char *s)
{
  GLenum error;
  while ((error = glGetError()) != GL_NO_ERROR) {
    fprintf(stderr, "%s: error - %s\n", s, (char *) gluErrorString(error));
  }
}

static void handleCgError() 
{
  fprintf(stderr, "Cg error: %s\n", cgGetErrorString(cgGetError()));
}

////////////////////////////////////////////////////////////////////////////////////

void initFPBuffer( int pW, int pH, int cW, int cH )
{
#if 0
  if ( !glh_init_extensions(EXTENSIONS_REQUIRED) ) 
    {
      cerr << "Unable to load the following extension(s): "
	   << glh_get_unsupported_extensions() << "\n\nExiting...\n" << endl;
      exit(-1);
    }
#endif

#if 0  
  GLenum err = glewInit();
  if( err != GLEW_OK ) cout << "Error: %s" << glewGetErrorString(err) << endl;

#endif

  wglSwapIntervalEXT(0);

#if 0
  GLint max_width;
  GLint max_height;
  ::glGetIntegerv( WGL_MAX_PBUFFER_WIDTH_ARB,  &max_width );
  ::glGetIntegerv( WGL_MAX_PBUFFER_HEIGHT_ARB, &max_height );
  cout << "max: " << max_width << " x " << max_height << " pbuffer size will be accepted. " << endl;
#endif

  //
  // 32bit floating-point buffer used in the first pass
  //

  if ( pbuffer ) delete pbuffer;

#ifndef RENDER_FBO

  pbuffer = new RenderTexture( "float=32 rgba depth textureRECT",
			       pW, pH, GL_TEXTURE_RECTANGLE_EXT );

#else // RENDER_FBO
  
  pbuffer = new RenderTextureFBO( pW, pH, 
				  texTarget0, texInternalFormat0, texFormat0, texType0, 
				  1 );

#endif // RENDER_FBO

  pbuffer->Activate();
  pane.initGL();
  pbuffer->Deactivate();

  //
  // 16bit floating-point buffer used in the second pass (including blending operations)
  //

  if ( pbuffer16 ) delete pbuffer16;

#ifndef RENDER_FBO

  pbuffer16 = new RenderTexture( "float=16 aux=2 rgba textureRECT",
				 pW, pH, GL_TEXTURE_RECTANGLE_EXT );

#else // RENDER_FBO

  pbuffer16 = new RenderTextureFBO( pW, pH, 
				    texTarget1, texInternalFormat1, texFormat1, texType1, 
				    2 );

#endif // RENDER_FBO

  pbuffer16->Activate();
  pane.initGL();
  pbuffer16->Deactivate();

  if (cap_buffer) delete cap_buffer;
  
#ifndef RENDER_FBO

  cap_buffer = new RenderTexture( "rgba textureRECT", cW, cH, GL_TEXTURE_RECTANGLE_EXT );

#else // RENDER_FBO

  cap_buffer = new RenderTextureFBO( pW, pH, 
				    texTarget2, texInternalFormat2, texFormat2, texType2, 
				     1 );

#endif // RENDER_FBO

  cap_buffer->Activate();
  pane.initGL();
  cap_buffer->Deactivate();
}

void destroyFPBuffer()
{
  delete pbuffer;
  pbuffer = NULL;
  delete pbuffer16;
  pbuffer16 = NULL;
  delete cap_buffer;
  cap_buffer = NULL;
}

//
// initialize textures to store the contents of floating-point buffers
//
void initRenderTexture( int pW, int pH )
{
#ifndef RENDER_FBO
  // q_hat 
  glGenTextures( 1, &qhatT );
  glBindTexture( GL_TEXTURE_RECTANGLE_EXT, qhatT );
  glTexImage2D( GL_TEXTURE_RECTANGLE_EXT, 0, GL_FLOAT_RGBA_NV, pW, pH, 
		0, GL_RGBA, GL_FLOAT, NULL );

  // position
  glGenTextures( 1, &positionT );
  glBindTexture( GL_TEXTURE_RECTANGLE_EXT, positionT );
  glTexImage2D( GL_TEXTURE_RECTANGLE_EXT, 0, GL_FLOAT_RGBA_NV, pW, pH, 
		0, GL_RGBA, GL_FLOAT, NULL );

  // normal vector
  glGenTextures( 1, &normalT );
  glBindTexture( GL_TEXTURE_RECTANGLE_EXT, normalT );
  glTexImage2D( GL_TEXTURE_RECTANGLE_EXT, 0, GL_FLOAT_RGBA_NV, pW, pH, 
		0, GL_RGBA, GL_FLOAT, NULL );
#endif // RENDER_FBO

  // initialize GL_ARB_point_sprite
  float maxSize;
  glGetFloatv( GL_POINT_SIZE_MAX_ARB, &maxSize );
  cout << "max point sprite size: " << maxSize << endl;

#ifdef OGL_QUERY
  GLint bitsSupported;
  glGetQueryivARB( GL_SAMPLES_PASSED_ARB, GL_QUERY_COUNTER_BITS_ARB, &bitsSupported );
  cout << "Number of counter bits = " << bitsSupported << endl;

  glGenQueriesARB( 1, &oq_slim1 );
  glGenQueriesARB( 1, &oq_slim2 );
#endif

}

void destroyRenderTextures()
{
#ifndef RENDER_FBO
  glDeleteTextures( 1, &qhatT );
  glDeleteTextures( 1, &positionT );
  glDeleteTextures( 1, &normalT );
#endif // RENDER_FBO
}

////////////////////////////////////////////////////////////////////////////////////

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
  vert1Program  = cgCreateProgramFromFile( context, CG_SOURCE, "shader_vp_1st.cg",
 					  vertProfile, NULL, NULL );
  cgGLLoadProgram( vert1Program );

  // fragment program (1st-pass)
  frag1Program = cgCreateProgramFromFile( context, CG_SOURCE, "shader_fp_1st.cg",
					  fragProfile, NULL, NULL );
  cgGLLoadProgram( frag1Program );

  // vertex program (2nd-pass)
  vert2Program  = cgCreateProgramFromFile( context, CG_SOURCE, "shader_vp_2nd.cg",
 					  vertProfile, NULL, NULL );
  cgGLLoadProgram( vert2Program );

  // fragment program (2nd-pass)
  frag2Program = cgCreateProgramFromFile( context, CG_SOURCE, "shader_fp_2nd.cg",
					  fragProfile, NULL, NULL );
  cgGLLoadProgram( frag2Program );

  // fragment program (3rd-pass)
  frag3Program = cgCreateProgramFromFile( context, CG_SOURCE, "shader_fp_3rd.cg",
					  fragProfile, NULL, NULL );
  cgGLLoadProgram( frag3Program );

  // Material properties
  g_CGparam_DiffMat = cgGetNamedParameter(frag3Program, "mat.DiffMat");
  g_CGparam_SpecMat = cgGetNamedParameter(frag3Program, "mat.SpecMat");
  g_CGparam_AmbMat = cgGetNamedParameter(frag3Program, "mat.AmbMat");
  g_CGparam_ShineMat = cgGetNamedParameter(frag3Program, "mat.ShineMat");
  g_CGparam_EmisMat = cgGetNamedParameter(frag3Program, "mat.EmisMat");

  // Lighting properties
  g_CGparam_LightVec = cgGetNamedParameter(frag3Program, "lt.LightVec");
  g_CGparam_LightPos = cgGetNamedParameter(frag3Program, "lt.LightPos");
  g_CGparam_DiffLight = cgGetNamedParameter(frag3Program, "lt.DiffLight");
  g_CGparam_SpecLight = cgGetNamedParameter(frag3Program, "lt.SpecLight");
  g_CGparam_AmbLight = cgGetNamedParameter(frag3Program, "lt.AmbLight");
  g_CGparam_AttenLight = cgGetNamedParameter(frag3Program, "lt.AttenLight");
  g_CGparam_SpotLight = cgGetNamedParameter(frag3Program, "lt.SpotLight");
  
#if 0
  // for debug
  debugProgram = cgCreateProgramFromFile( context, CG_SOURCE, "debug.cg",
					  fragProfile, NULL, NULL );
  cgGLLoadProgram( debugProgram );
#endif
}

void destroyCg()
{
  cgDestroyProgram( vert1Program );
  cgDestroyProgram( vert2Program );
  cgDestroyProgram( frag1Program );
  cgDestroyProgram( frag2Program );
  cgDestroyProgram( frag3Program );
//   cgDestroyProgram( debugProgram ); // for debug
  cgDestroyContext( context );
}

////////////////////////////////////////////////////////////////////////////////////

#ifdef OLD_SLIM
void setPolyTexture()
{
  GLuint polyParam;
  glGenTextures( 1, &polyParam );
  glBindTexture( GL_TEXTURE_RECTANGLE_EXT, polyParam );
  glTexParameteri( GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  glTexParameteri( GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

  GLsizei sizeW;
  glGetIntegerv( GL_MAX_RECTANGLE_TEXTURE_SIZE_EXT, &sizeW );

  // Store coefficients of polynomial to a texture
  // Each polynomial uses 4 pixels to store coefficients and a center point
  std::vector<float> pp;
#if 0
  int max_width = 512;
  int tex_width  = max_width * 4; // pixel width
#endif
  int tex_width = sizeW;
  int max_width = (int) (tex_width / 4.0);
  cout << "max slimball size: " << max_width * sizeW << endl;

  int flo_width = tex_width * 4;
  for ( int k = 0; k < flo_width; ++k ) pp.push_back( 0.0f );
  int tex_height = 1;

  int count = 0;
  int pp_count = 0;
  int px = 0;
  int py = 0;
  int id_count = 0;
  int leaf_count = 0;
  for( int i = 0; i < 2*listN; i++ )
    {
      if( Ns[i] == 0 ) continue;
//       if ( i%2 != 0 ) continue;

      Quadratic** list = (Quadratic**) lists[i];
      for( int j = 0; j < Ns[i]; j++ )
	{
// 	  if ( !(list[j]->leaf) ) continue;
	  if ( list[j]->leaf ) ++(leaf_count);

	  list[j]->id_ = id_count; ++id_count;
	  
	  texcoordParam.push_back( Point2f( (float) px / FID_DIV, 
					    (float) py / FID_DIV ) );
	  
	  // px, py
	  pp[ pp_count ] = list[j]->centerX; ++( pp_count );
	  pp[ pp_count ] = list[j]->centerY; ++( pp_count );
	  pp[ pp_count ] = list[j]->centerZ; ++( pp_count );
	  pp[ pp_count ] = list[j]->support; ++( pp_count );

	  // px+1, py
	  pp[ pp_count ] = list[j]->cXX; ++( pp_count );
          pp[ pp_count ] = list[j]->cYY; ++( pp_count );
	  pp[ pp_count ] = list[j]->cZZ; ++( pp_count );
	  pp[ pp_count ] = .0f; ++( pp_count );
          
	  // px+2, py
	  pp[ pp_count ] = list[j]->cXY; ++( pp_count );
	  pp[ pp_count ] = list[j]->cYZ; ++( pp_count );
	  pp[ pp_count ] = list[j]->cZX; ++( pp_count );
	  pp[ pp_count ] = .0f; ++( pp_count );
          
	  // px+3, py
	  pp[ pp_count ] = list[j]->cX; ++( pp_count );
	  pp[ pp_count ] = list[j]->cY; ++( pp_count );
	  pp[ pp_count ] = list[j]->cZ; ++( pp_count );
	  pp[ pp_count ] = list[j]->c0; ++( pp_count );

// 	  // px, py
// 	  pp.push_back( list[j]->centerX );
// 	  pp.push_back( list[j]->centerY );
// 	  pp.push_back( list[j]->centerZ );
// 	  pp.push_back( list[j]->support );

// 	  // px+1, py
// 	  pp.push_back( list[j]->cXX );
//           pp.push_back( list[j]->cYY );
// 	  pp.push_back( list[j]->cZZ );
// 	  pp.push_back( .0f );
          
// 	  // px+2, py
// 	  pp.push_back( list[j]->cXY );
// 	  pp.push_back( list[j]->cYZ );
// 	  pp.push_back( list[j]->cZX );
// 	  pp.push_back( .0f );
          
// 	  // px+3, py
// 	  pp.push_back( list[j]->cX );
// 	  pp.push_back( list[j]->cY );
// 	  pp.push_back( list[j]->cZ );
// 	  pp.push_back( list[j]->c0 );

	  px += 4;
	  count++;
	  if ( count == max_width )
	    {
	      count = 0;
	      ++tex_height;

	      for ( int k = 0; k < flo_width; ++k ) pp.push_back( 0.0f );
	      
	      px = 0;
	      ++py;
	    }

	}
    }
  
  cout << "Slim points ... all: " << id_count << " leaf " << leaf_count << endl;
  cout << "Param texture: width " << tex_width << " height " << tex_height 
       << " pp size " << pp.size() << endl;
  
  //GL_FLOAT_RGBA_NV
  glTexImage2D( GL_TEXTURE_RECTANGLE_EXT, 0, GL_RGBA32F_ARB, 
		tex_width, tex_height, 0, GL_RGBA, GL_FLOAT, &pp[0] );
  cgGLSetTextureParameter( cgGetNamedParameter(frag1Program, "polyParam"), polyParam );
  cgGLSetTextureParameter( cgGetNamedParameter(frag2Program, "polyParam"), polyParam );
}
#endif // OLD_SLIM

////////////////////////////////////////////////////////////////////////////////////

#ifdef NEW_SLIM
void setSlimTreePolyTexture()
{
  GLuint polyParam;
  glGenTextures( 1, &polyParam );
  glBindTexture( GL_TEXTURE_RECTANGLE_EXT, polyParam );
  glTexParameteri( GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  glTexParameteri( GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

  GLsizei sizeW;
  glGetIntegerv( GL_MAX_RENDERBUFFER_SIZE_EXT, &sizeW );

  // Store coefficients of polynomial to a texture
  // Each polynomial uses 4 pixels to store coefficients and a center point
#if 0
  int max_width = 512;
  int tex_width  = max_width * 4; // pixel width
#endif
  int tex_width = sizeW;
  int max_width = (int) (tex_width / 4.0);
  cout << "max slimball size: " << max_width * sizeW << endl;

  int flo_width = tex_width * 4;
  std::vector<float> pp( flo_width );
//   for ( int k = 0; k < flo_width; ++k ) pp.push_back( 0.0f );

  int tex_height = 1;

  int count = 0;
  int pp_count = 0;
  int px = 0;
  int py = 0;

  std::vector<SlimBallf*>& slimballs = slimtree.slimballs();
  for ( int i = 0; i < slimballs.size(); ++i )
    {
      SlimBallf* sb    = slimballs[i];

      texcoordParam.push_back( Point2f( (float) px / FID_DIV, 
					(float) py / FID_DIV ) );

      // px, py
      Point3f& center = sb->center();
      pp[ pp_count++ ] = center.x;
      pp[ pp_count++ ] = center.y;
      pp[ pp_count++ ] = center.z;
      pp[ pp_count++ ] = sb->support();

      // px+1, py
      Quadraticf& bf   = (Quadraticf&) sb->bf();
//       cout << bf << endl;
      std::vector<float>& c = bf.coeffs();
      pp[ pp_count++ ] = c[0]; // XX
      pp[ pp_count++ ] = c[1]; // YY
      pp[ pp_count++ ] = c[2]; // ZZ
      pp[ pp_count++ ] = .0f;
          
      // px+2, py
      pp[ pp_count++ ] = c[3];
      pp[ pp_count++ ] = c[4];
      pp[ pp_count++ ] = c[5];
      pp[ pp_count++ ] = .0f;
          
      // px+3, py
      pp[ pp_count++ ] = c[6];
      pp[ pp_count++ ] = c[7];
      pp[ pp_count++ ] = c[8];
      pp[ pp_count++ ] = c[9];

      px += 4;
      count++;
      if ( count == max_width )
	{
	  count = 0;
	  ++tex_height;

	  for ( int k = 0; k < flo_width; ++k ) pp.push_back( 0.0f );
	  
	  px = 0;
	  ++py;
	}
    }

  cout << "Param texture: width " << tex_width << " height " << tex_height 
       << " pp size " << pp.size() << endl;
  
  // GL_FLOAT_RGBA_NV
  glTexImage2D( GL_TEXTURE_RECTANGLE_EXT, 0, GL_RGBA32F_ARB, 
		tex_width, tex_height, 0, GL_RGBA, GL_FLOAT, &pp[0] );
  cgGLSetTextureParameter( cgGetNamedParameter(frag1Program, "polyParam"), polyParam );
  cgGLSetTextureParameter( cgGetNamedParameter(frag2Program, "polyParam"), polyParam );
}
#endif // NEW_SLIM

////////////////////////////////////////////////////////////////////////////////////

void setGaussKernelTexture()
{
  GLuint gaussKernel;
  glGenTextures( 1, &gaussKernel );
  glBindTexture( GL_TEXTURE_RECTANGLE_EXT, gaussKernel );
  glTexParameteri( GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  glTexParameteri( GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

  std::vector<float> gk;
  int pix_width = GAUSS_WIDTH;
  int flo_width = pix_width * 4;
  int pix_height = 1;

  for ( int k = 0; k < flo_width; ++k ) gk.push_back( 0.0f );

  int pix_count = 0;
  for ( int i = 0; i < pix_width; ++i )
    {
      float r = (float) i / (float) (pix_width-1);
#if 0
      gk[ pix_count ] = (float) std::exp( -r*r ); 
#endif

#if 1
      if ( r < .5f )
	gk[ pix_count ] = - r * r + .75f;
      else if ( (r >= .5f) && (r <= 1.0f) )
	gk[ pix_count ] = 0.5f * (1.5f - r) * (1.5f - r);
#endif
      ++( pix_count );

      gk[ pix_count ] = .0f; ++( pix_count );
      gk[ pix_count ] = .0f; ++( pix_count );
      gk[ pix_count ] = .0f; ++( pix_count );
    }
  
  // GL_FLOAT_RGBA_NV
  glTexImage2D( GL_TEXTURE_RECTANGLE_EXT, 0, GL_RGBA32F_ARB,
		pix_width, pix_height, 0, GL_RGBA, GL_FLOAT, &gk[0] );
  cgGLSetTextureParameter( cgGetNamedParameter(frag2Program, "gaussKernel"), gaussKernel );
}

////////////////////////////////////////////////////////////////////////////////////

void MakeGlutWindowCurrent()
{
  static int glutWinId = glutGetWindow();

  glutSetWindow( glutWinId );
}

////////////////////////////////////////////////////////////////////////////////////

void drawQuad( int w, int h )
{
  glBegin( GL_QUADS );
  glTexCoord2f(0, 0); glVertex2f(-1, -1);
  glTexCoord2f(0, h); glVertex2f(-1,  1);
  glTexCoord2f(w, h); glVertex2f( 1,  1);
  glTexCoord2f(w, 0); glVertex2f( 1, -1);
  glEnd();
}

float zEye( float m[16], float p[16], float x, float y, float z )
{
  float a = m[0]*x+m[4]*y+m[ 8]*z+m[12];
  float b = m[1]*x+m[5]*y+m[ 9]*z+m[13];
  float c = m[2]*x+m[6]*y+m[10]*z+m[14];
  float d = m[3]*x+m[7]*y+m[11]*z+m[15];
  
  return (p[2]*a+p[6]*b+p[10]*c+p[14]*d);
}

#include "ViewFrustum.hxx"

ViewFrustum vf;

#define SLIMVIEW_FEW_PIX 4

////////////////////////////////////////////////////////////////////////////////////
//
// old slim
//
////////////////////////////////////////////////////////////////////////////////////

#ifdef OLD_SLIM
void drawBillboardsLOD( float size_param )
{
  float mat[16];
  glGetFloatv( GL_MODELVIEW_MATRIX, mat );

  Point3f vRight( mat[0], mat[4], mat[8] );
  Point3f vUp( mat[1], mat[5], mat[9] );

  float prj[16];
  glGetFloatv( GL_PROJECTION_MATRIX, prj );
  
  vf.update( mat, prj );

  imagemaker.visitID++;
  imagemaker.stack[0] = imagemaker.root;
  imagemaker.top = 0;
  imagemaker.root->visit = imagemaker.visitID;

  slimball_count = 0;

  glBegin( GL_QUADS );

  while( imagemaker.top >= 0 )
    {
      BasisFunction* bf = imagemaker.stack[imagemaker.top--];

      // view frustum culling
      bool f = vf.sphereInFrustum( bf->centerX, bf->centerY, bf->centerZ, bf->support );
      if ( f != true ) continue;
      
      // radius in screen coordinate
      float z = zEye( mat, prj, bf->centerX, bf->centerY, bf->centerZ );
      float rad_screen = bf->support * size_param / z;
      
      if ( bf->leaf | (rad_screen < SLIMVIEW_FEW_PIX) )
	{
	  //  	  if(!(imagemaker.both) && (bf->normalZAtC() > 0) )
	  //  	    continue;

	  Point3f vCenter( bf->centerX, bf->centerY, bf->centerZ );
	  float pSize = bf->support;

	  Point3f vPoint0( vCenter + ((-1.0f*vRight - vUp) * pSize) );
	  Point3f vPoint1( vCenter + (( vRight - vUp) * pSize) );
	  Point3f vPoint2( vCenter + (( vRight + vUp) * pSize) );
	  Point3f vPoint3( vCenter + ((-1.0f*vRight + vUp) * pSize) );
	  
	  // 2D coordinate of a parameter texture
	  //glColor4f( texcoordParam[k].x, texcoordParam[k].y, 1.0f, pSize / FID_DIV );
	  glColor4f( texcoordParam[bf->id_].x, texcoordParam[bf->id_].y, 1.0f, pSize / FID_DIV );

	  glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
	  glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
	  glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
	  glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );
	  
	  ++( slimball_count );
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
}

void drawBillboards()
{
  float mat[16];
  glGetFloatv( GL_MODELVIEW_MATRIX, mat );

  Point3f vRight( mat[0], mat[4], mat[8] );
  Point3f vUp( mat[1], mat[5], mat[9] );

  //
  // draw quads
  //

  slimball_count = 0;

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

	  Point3f vPoint0( vCenter + ((-1.0f*vRight - vUp) * pSize) );
	  Point3f vPoint1( vCenter + (( vRight - vUp) * pSize) );
	  Point3f vPoint2( vCenter + (( vRight + vUp) * pSize) );
	  Point3f vPoint3( vCenter + ((-1.0f*vRight + vUp) * pSize) );

	  // 2D coordinate of a parameter texture, (none), support
	  glColor4f( texcoordParam[list[j]->id_].x, texcoordParam[list[j]->id_].y, 1.0f, pSize / FID_DIV );

	  glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
	  glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
	  glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
	  glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );

	  ++slimball_count;


	}
    }
  glEnd();
}

void drawPointSpritesLOD( float size_param )
{
  float mat[16];
  glGetFloatv( GL_MODELVIEW_MATRIX, mat );

  float prj[16];
  glGetFloatv( GL_PROJECTION_MATRIX, prj );
  
  vf.update( mat, prj );

  glEnable( GL_POINT_SPRITE_ARB );
  glTexEnvi( GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE );
  glEnable( GL_VERTEX_PROGRAM_POINT_SIZE_ARB );

  slimball_count = 0;
  
  glBegin( GL_POINTS );

  imagemaker.visitID++;
  imagemaker.stack[0] = imagemaker.root;
  imagemaker.top = 0;
  imagemaker.root->visit = imagemaker.visitID;
  while( imagemaker.top >= 0 )
    {
      BasisFunction* bf = imagemaker.stack[imagemaker.top--];

      // view frustum culling
      bool f = vf.sphereInFrustum( bf->centerX, bf->centerY, bf->centerZ, bf->support );
      if ( f != true ) continue;
      
      // radius in screen coordinate
      float z = zEye( mat, prj, bf->centerX, bf->centerY, bf->centerZ );
      float rad_screen = bf->support * size_param / z;
      
      if ( bf->leaf | (rad_screen < SLIMVIEW_FEW_PIX) )
	{
	  //  	  if(!(imagemaker.both) && (bf->normalZAtC() > 0) )
	  //  	    continue;

	  // 2D coordinate of a parameter texture
	  glColor4f( texcoordParam[bf->id_].x, texcoordParam[bf->id_].y, 1.0f, bf->support / FID_DIV );

	  glVertex3f( bf->centerX, bf->centerY, bf->centerZ );

	  ++( slimball_count );


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

  glDisable( GL_VERTEX_PROGRAM_POINT_SIZE_ARB );
  glDisable( GL_POINT_SPRITE_ARB );
}

void drawPointSprites()
{
  glEnable( GL_POINT_SPRITE_ARB );
  glTexEnvi( GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE );
  glEnable( GL_VERTEX_PROGRAM_POINT_SIZE_ARB );

  slimball_count = 0;
  
  glBegin( GL_POINTS );
  for( int i = 0; i < 2*listN; i++ )
    {
      if( Ns[i] == 0 ) continue;
      if ( i%2 != 0 ) continue;
      BasisFunction** list = lists[i];
      for( int j = 0; j < Ns[i]; j++ )
	{
	  if ( !(list[j]->leaf) ) continue;
	  
	  // 2D coordinate of a parameter texture, ID, support
	  glColor4f( texcoordParam[list[j]->id_].x, texcoordParam[list[j]->id_].y, 1.0f, list[j]->support / FID_DIV );
	  //
	  // (center_x, center_y, center_z, radius of support)
	  //
	  glVertex3f( list[j]->centerX, list[j]->centerY, list[j]->centerZ );

	  ++( slimball_count );

	}
    }
  glEnd();

  glDisable( GL_VERTEX_PROGRAM_POINT_SIZE_ARB );
  glDisable( GL_POINT_SPRITE_ARB );
}

void displayBall()
{
  if ( isCapture ) cap_buffer->Activate();

  pane.setBackgroundColor( bgColor[0], bgColor[1], bgColor[2] );
  pane.setIsGradientBackground( false );
  pane.clear( pbuffer_width, pbuffer_height );
  pane.setView();
  pane.setLight();

  ::glShadeModel ( GL_SMOOTH );
  ::glEnable( GL_LIGHTING );
  
  ballMtl.bind();

  slimball_count = 0;

  for( int i = 0; i < 2*listN; i++ )
    {
      if( Ns[i] == 0 ) continue;
      if ( i%2 != 0 ) continue;

      BasisFunction** list = lists[i];
      for( int j = 0; j < Ns[i]; j++ )
	{
	  if ( !(list[j]->leaf) ) continue;
	  
	  glPushMatrix();
	  glTranslatef( list[j]->centerX, list[j]->centerY, list[j]->centerZ );
	  glutSolidSphere( list[j]->support, 10, 10 );
	  glPopMatrix();

	  ++( slimball_count );

	}
    }

  pane.finish();
}
  
#endif // OLD_SLIM

////////////////////////////////////////////////////////////////////////////////////
//
// new slim
//
////////////////////////////////////////////////////////////////////////////////////

#ifdef NEW_SLIM

void drawSlimTreeBillboards()
{
  float mat[16];
  glGetFloatv( GL_MODELVIEW_MATRIX, mat );

  Point3f vRight( mat[0], mat[4], mat[8] );
  Point3f vUp( mat[1], mat[5], mat[9] );

  //
  // draw quads
  //

  slimball_count = 0;

  glBegin( GL_QUADS );

  std::vector<SlimBallf*>& slimballs = slimtree.slimballs();
  for ( int i = 0; i < slimballs.size(); ++i )
    {
      SlimBallf* slimball = slimballs[i];
      
      if ( slimball->childs().size() ) continue;

      Point3f vCenter( slimball->center() );
      float pSize = slimball->support();

      Point3f vPoint0( vCenter + ((-1.0f*vRight - vUp) * pSize) );
      Point3f vPoint1( vCenter + (( vRight - vUp) * pSize) );
      Point3f vPoint2( vCenter + (( vRight + vUp) * pSize) );
      Point3f vPoint3( vCenter + ((-1.0f*vRight + vUp) * pSize) );
      
      // 2D coordinate of a parameter texture, (none), support
      int id = slimball->id();
      glColor4f( texcoordParam[id].x, texcoordParam[id].y, 1.0f, pSize / FID_DIV );

      glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
      glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
      glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
      glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );

      ++( slimball_count );
    }

  glEnd();
}

int drawSlimBallBillboardsLOD( SlimBallf* slimball, int level, Point3f& vRight, Point3f& vUp )
{
  if ( !slimball ) return 0;

  int count = 0;
  if ( level == slimball_level )
    {
      glBegin( GL_QUADS );
      
      Point3f vCenter( slimball->center() );
      float pSize = slimball->support();

      Point3f vPoint0( vCenter + ((-1.0f*vRight - vUp) * pSize) );
      Point3f vPoint1( vCenter + (( vRight - vUp) * pSize) );
      Point3f vPoint2( vCenter + (( vRight + vUp) * pSize) );
      Point3f vPoint3( vCenter + ((-1.0f*vRight + vUp) * pSize) );
      
      // 2D coordinate of a parameter texture, (none), support
      int id = slimball->id();
      glColor4f( texcoordParam[id].x, texcoordParam[id].y, 1.0f, pSize / FID_DIV );

      glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
      glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
      glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
      glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );
      
      glEnd();

      ++( count );

      return count;
    }

  if ( slimball->childs().size() )
    {
      for ( int i = 0; i < slimball->childs().size(); ++i )
	count += drawSlimBallBillboardsLOD( slimball->child(i), level+1, vRight, vUp );
    }

  return count;
}

void drawSlimTreeBillboardsLOD()
{
  float mat[16];
  glGetFloatv( GL_MODELVIEW_MATRIX, mat );

  Point3f vRight( mat[0], mat[4], mat[8] );
  Point3f vUp( mat[1], mat[5], mat[9] );

  //
  // draw quads
  //

  SlimBallf* slimball = slimtree.root();
  slimball_count = drawSlimBallBillboardsLOD( slimball, 0, vRight, vUp );
}

int drawSlimBallBillboardsLODError( SlimBallf* slimball, int level, Point3f& vRight, Point3f& vUp )
{
  if ( !slimball ) return 0;

//   cout << "error " << slimball->userDefined() << " threshold " << slimball_error << endl;
  int count = 0;
  if ( slimball->userDefined() < slimball_error )
    {
      glBegin( GL_QUADS );
      
      Point3f vCenter( slimball->center() );
      float pSize = slimball->support();

      Point3f vPoint0( vCenter + ((-1.0f*vRight - vUp) * pSize) );
      Point3f vPoint1( vCenter + (( vRight - vUp) * pSize) );
      Point3f vPoint2( vCenter + (( vRight + vUp) * pSize) );
      Point3f vPoint3( vCenter + ((-1.0f*vRight + vUp) * pSize) );
      
      // 2D coordinate of a parameter texture, (none), support
      int id = slimball->id();
      glColor4f( texcoordParam[id].x, texcoordParam[id].y, 1.0f, pSize / FID_DIV );

      glVertex3f( vPoint0.x, vPoint0.y, vPoint0.z );
      glVertex3f( vPoint1.x, vPoint1.y, vPoint1.z );
      glVertex3f( vPoint2.x, vPoint2.y, vPoint2.z );
      glVertex3f( vPoint3.x, vPoint3.y, vPoint3.z );
      
      glEnd();

      ++( count );

      return count;
    }

  if ( slimball->childs().size() )
    {
      for ( int i = 0; i < slimball->childs().size(); ++i )
	count += drawSlimBallBillboardsLODError( slimball->child(i), level+1, vRight, vUp );
    }

  return count;
}

void drawSlimTreeBillboardsLODError()
{
  float mat[16];
  glGetFloatv( GL_MODELVIEW_MATRIX, mat );

  Point3f vRight( mat[0], mat[4], mat[8] );
  Point3f vUp( mat[1], mat[5], mat[9] );

  //
  // draw quads
  //

  SlimBallf* slimball = slimtree.root();
  slimball_count = drawSlimBallBillboardsLODError( slimball, 0, vRight, vUp );
}

void drawSlimTreePointSprites()
{
  glEnable( GL_POINT_SPRITE_ARB );
  glTexEnvi( GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE );
  glEnable( GL_VERTEX_PROGRAM_POINT_SIZE_ARB );
  
  slimball_count = 0;

  glBegin( GL_POINTS );

  std::vector<SlimBallf*>& slimballs = slimtree.slimballs();
  for ( int i = 0; i < slimballs.size(); ++i )
    {
      SlimBallf* slimball = slimballs[i];
      
      if ( slimball->childs().size() ) continue;

      // 2D coordinate of a parameter texture, ID, support
      float pSize = slimball->support();
      int id = slimball->id();
      glColor4f( texcoordParam[id].x, texcoordParam[id].y, 1.0f, pSize / FID_DIV );
      //
      // (center_x, center_y, center_z, radius of support)
      //
      Point3f vCenter( slimball->center() );
      glVertex3f( vCenter.x, vCenter.y, vCenter.z );

      ++( slimball_count );
    }

  glEnd();

  glDisable( GL_VERTEX_PROGRAM_POINT_SIZE_ARB );
  glDisable( GL_POINT_SPRITE_ARB );
}

void drawSlimTreeBall()
{
  if ( isCapture ) cap_buffer->Activate();

  pane.setBackgroundColor( bgColor[0], bgColor[1], bgColor[2] );
  pane.setIsGradientBackground( false );
  pane.clear( pbuffer_width, pbuffer_height );
  pane.setView();
  pane.setIsLightOn( 1, false );
  pane.setLight();

  ::glShadeModel ( GL_SMOOTH );
  ::glEnable( GL_LIGHTING );
  
//   ballMtl.bind();

  slimball_count = 0;

  std::vector<SlimBallf*>& slimballs = slimtree.slimballs();
  for ( int i = 0; i < slimballs.size(); ++i )
    {
      if ( slimballs[i]->childs().size() ) continue;

      glPushMatrix();
      Point3f& center = slimballs[i]->center();
      glTranslatef( center.x, center.y, center.z );

      float color[3];
      convColor( slimballs[i]->userDefined(), color );
      ballMtl.setDiffuseColor( color[0], color[1], color[2], 1.0f );
      ballMtl.bind();

      glutSolidSphere( slimballs[i]->support(), 50, 50 );
      glPopMatrix();

      ++( slimball_count );
    }

  pane.finish();
}

int drawSlimBallBallLODError( SlimBallf* slimball, int level )
{
  if ( !slimball ) return 0;

  int count = 0;
  if ( slimball->userDefined() < slimball_error )
    {
      glPushMatrix();
      Point3f& center = slimball->center();
      glTranslatef( center.x, center.y, center.z );
      
      float color[3];
      convColor( slimball->userDefined(), color );
      ballMtl.setDiffuseColor( color[0], color[1], color[2], 1.0f );
      ballMtl.bind();
      
      glutSolidSphere( slimball->support(), 10, 10 );
      glPopMatrix();

      ++( count );

      return count;
    }

  if ( slimball->childs().size() )
    {
      for ( int i = 0; i < slimball->childs().size(); ++i )
	count += drawSlimBallBallLODError( slimball->child(i), level+1 );
    }

  return count;
}

void drawSlimTreeBallLODError()
{
  if ( isCapture ) cap_buffer->Activate();

  pane.setBackgroundColor( bgColor[0], bgColor[1], bgColor[2] );
  pane.setIsGradientBackground( false );
  pane.clear( pbuffer_width, pbuffer_height );
  pane.setView();
  pane.setLight();

  ::glShadeModel ( GL_SMOOTH );
  ::glEnable( GL_LIGHTING );

//   ballMtl.bind();

  slimball_count = 0;
  SlimBallf* slimball = slimtree.root();
  slimball_count += drawSlimBallBallLODError( slimball, 0 );

  pane.finish();
}


  
#endif // NEW_SLIM

////////////////////////////////////////////////////////////////////////////////////

void displayPoly()
{
  if ( isCapture ) cap_buffer->Activate();

  pane.setBackgroundColor( bgColor[0], bgColor[1], bgColor[2] );
  pane.setIsGradientBackground( false );
  pane.clear( pbuffer_width, pbuffer_height );
  pane.setView();
  pane.setLight();

  if ( !(glmeshR.empty()) ) glmeshR.draw();

  pane.finish();
}
  
////////////////////////////////////////////////////////////////////////////////////

// #define SLIM_DEBUG 1
// #define SLIM_DEBUG_II 1

#ifdef OLD_SLIM
void drawSlimTree( float* size_param )
{
  if ( draw_mode == DRAW_BILLBOARD )
    {
      if ( lod_mode == LOD_LEVEL )
	drawBillboardsLOD( size_param[0] * size_param[1] );
      else
	drawBillboards();
    }
  else
    {
      if ( lod_mode == LOD_LEVEL )
	drawPointSpritesLOD( size_param[0] * size_param[1] );
      else
	drawPointSprites();
    }
}

#endif // OLD_SLIM

#ifdef NEW_SLIM
void drawSlimTree( float* size_param )
{
  if ( draw_mode == DRAW_BILLBOARD )
    {
      if ( lod_mode == LOD_LEVEL )
	{
	  drawSlimTreeBillboardsLOD();
	}
      else if ( lod_mode == LOD_ERROR )
	{
	  drawSlimTreeBillboardsLODError();
	}
      else
	drawSlimTreeBillboards();
    }
  else
    {
      drawSlimTreePointSprites();
    }
}

void displaySlimTreeBall()
{
  if ( lod_mode == LOD_ERROR )
    drawSlimTreeBallLODError();
  else
    drawSlimTreeBall();
}

#endif // NEW_SLIM

void displaySLIMGPU()
{
  ////////////////////////////////////////////////////////////////////////////////////
  //
  // 1st pass
  //
  ////////////////////////////////////////////////////////////////////////////////////

  pbuffer->Activate();

  pane.setBackgroundColor( .0f, .0f, .0f );
  pane.setIsGradientBackground( false );
  pane.clear( pbuffer_width, pbuffer_height );
  pane.setView();

#ifdef OGL_QUERY
  glBeginQueryARB( GL_SAMPLES_PASSED_ARB, oq_slim1 );
#endif

  //
  // vertex program
  //
  cgGLEnableProfile( vertProfile );
  cgGLBindProgram( vert1Program );

  cgGLSetStateMatrixParameter( cgGetNamedParameter(vert1Program, "ModelViewProj"),
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
  float size_param[2];
  size_param[0] = (float) projMatrix[5];
  size_param[1] = (float) viewport[3];
  cgGLSetParameter2fv( cgGetNamedParameter(vert1Program, "size_param"), size_param );

  //
  // fragment program
  //
  cgGLEnableProfile( fragProfile );
  cgGLBindProgram( frag1Program );

  cgGLEnableTextureParameter( cgGetNamedParameter(frag1Program, "polyParam") );

  cgGLSetStateMatrixParameter( cgGetNamedParameter(frag1Program, "ModelView"),
 			       CG_GL_MODELVIEW_MATRIX,
			       CG_GL_MATRIX_IDENTITY );

  cgGLSetStateMatrixParameter( cgGetNamedParameter(frag1Program, "ModelViewInv"),
 			       CG_GL_MODELVIEW_MATRIX,
			       CG_GL_MATRIX_INVERSE );

  cgGLSetStateMatrixParameter( cgGetNamedParameter(frag1Program, "ModelViewProj"),
			       CG_GL_MODELVIEW_PROJECTION_MATRIX,
 			       CG_GL_MATRIX_IDENTITY );

  cgGLSetStateMatrixParameter(cgGetNamedParameter(frag1Program, "ModelViewProjInv"),
			      CG_GL_MODELVIEW_PROJECTION_MATRIX,
			      CG_GL_MATRIX_INVERSE );

  //
  // get width and height of near plane
  //
  float left1, right1, bottom1, top1, near1, far1;
  pane.getProjectionParametersFromPerspective( &left1, &right1, &bottom1, &top1, &near1, &far1 );
  float viewParam[4];
  viewParam[0] = near1;
  viewParam[1] = far1;
  cgGLSetParameter4fv( cgGetNamedParameter(frag1Program, "viewParam"), viewParam );

  // set draw mode ... 0: Billboard, 1; Point Sprite
  cgGLSetParameter1f( cgGetNamedParameter(frag1Program, "drawMode"), draw_mode );

  drawSlimTree( size_param );

  cgGLDisableTextureParameter( cgGetNamedParameter(frag1Program, "polyParam") );
  cgGLDisableProfile( fragProfile );

  cgGLDisableProfile( vertProfile );

  pane.finish();

#ifdef OGL_QUERY
  glEndQueryARB( GL_SAMPLES_PASSED_ARB );
#endif

  // TEMP
  //
  // for capturing depth
  //
#if 0

  std::vector<float> pd( pbuffer_width * pbuffer_height );
  glReadPixels(0, 0, pbuffer_width, pbuffer_height, GL_DEPTH_COMPONENT, GL_FLOAT, &pd[0] );
  float amax = .0f;
  float amin = .0f;
//   unsigned int amax = 0;
//   unsigned int amin = 0;
  bool famax = false;
   for ( int i = 0; i < pbuffer_width * pbuffer_height; ++i )
     {
       if ( famax )
	 {
	   if ( (pd[i] != 1.0) && (amax < pd[i]) ) amax = pd[i];
	 }
       else
	 {
	   if ( pd[i] != 1.0 )
	     {
	       amax = pd[i];
	       famax = true;
	     }
	 }
       if ( i )
	 {
	   if ( amin > pd[i] ) amin = pd[i];
	 }
       else
	 {
	   amin = pd[i];
	 }
     }
//    cout << " amax " << amax << " amin " << amin << endl;

   std::vector<unsigned char> pdb( pbuffer_width * pbuffer_height );
   float sub = amax - amin;
   for ( int i = 0; i < pbuffer_width * pbuffer_height; ++i )
     {
       pd[i] = (pd[i] - amin) / sub;
       if ( pd[i] > 1.0 ) pd[i] = 1.0;
       pdb[i] = (unsigned char) ((1.0-pd[i]) * 255.0);
//        if ( pd[i] > 1.0 ) pdb[i] = (unsigned char) 255;
     }

#if 1
  char filename[BUFSIZ];
  sprintf( filename, "tmp%03d.png", count_frame++ );
  pi.writeGrayImage( pbuffer_width, pbuffer_height, &pdb[0], filename );
#endif
  
  // TEMP
#endif

#if !defined(SLIM_DEBUG)

#ifndef RENDER_FBO
  // copy fpbuffer to a texture
  glBindTexture( GL_TEXTURE_RECTANGLE_EXT, qhatT );
  glCopyTexSubImage2D( GL_TEXTURE_RECTANGLE_EXT, 0, 0, 0, 0, 0, pbuffer_width, pbuffer_height );
  cgGLSetTextureParameter( cgGetNamedParameter(frag2Program, "qhatT"), qhatT );
#endif

#endif // #ifndef SLIM_DEBUG

  
#if defined(SLIM_DEBUG)

  std::vector<float> point_data( pbuffer_width*pbuffer_height*4 );

  //glReadBuffer( WGL_FRONT_LEFT_ARB );
  glReadPixels(0, 0, pbuffer_width, pbuffer_height, GL_RGBA, GL_FLOAT, &point_data[0] );

  char filename[BUFSIZ];
  sprintf( filename, "tmp%03d.pnt", count_frame );
  std::ofstream ofs( filename );
  int pixel_count = 0;
  float max_val;

  for ( int i = 0; i < pbuffer_width*pbuffer_height; ++i )
    {
      if ( (point_data[4*i] == 0.0) && (point_data[4*i+1] == 0.0) && (point_data[4*i+2] == 0.0) )
	{
	}
      else
	{
	  ofs << point_data[4*i] << " " << point_data[4*i+1] << " "
		<< point_data[4*i+2] << " " << point_data[4*i+3] << endl;

	  if ( pixel_count )
	    {
	      if ( max_val < point_data[4*i] )
		max_val = point_data[4*i];
	      if ( max_val < point_data[4*i+1] )
		max_val = point_data[4*i+1];
	    }
	  else
	    {
	      max_val = point_data[4*i];
	      if ( max_val < point_data[4*i+1] )
		max_val = point_data[4*i+1];
	    }

	  ++pixel_count;
	}
    }

  cout << "max val " << max_val << endl;
  cout << "pixel count " << pixel_count << endl;

  ++count_frame;

  ofs.close();

#endif // #ifdef SLIM_DEBUG

  pbuffer->Deactivate();

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // 2nd pass
  //
  ////////////////////////////////////////////////////////////////////////////////////

//   cout << "bb" << endl;

#if !defined(SLIM_DEBUG)

//   glPointParameterfARB( GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT );
  
  pbuffer16->Activate();

  for ( int i = 0; i < 2; ++i )
    {

#ifndef RENDER_FBO

      glDrawBuffer( buffers[i] );

#else // RENDER_FBO

      pbuffer16->Draw( i );

#endif // RENDER_FBO

      pane.setBackgroundColor( .0f, .0f, .0f );
      pane.setIsGradientBackground( false );
      pane.clear( pbuffer_width, pbuffer_height );
      pane.setView();
    }

#ifdef OGL_QUERY
  glBeginQueryARB( GL_SAMPLES_PASSED_ARB, oq_slim2 );
#endif

  //
  // vertex program
  //
  cgGLEnableProfile( vertProfile );
  cgGLBindProgram( vert2Program );

  cgGLSetStateMatrixParameter( cgGetNamedParameter(vert2Program, "ModelViewProj"),
			       CG_GL_MODELVIEW_PROJECTION_MATRIX,
			       CG_GL_MATRIX_IDENTITY );

  // assign the same parameter as a 1st pass
  cgGLSetParameter2fv( cgGetNamedParameter(vert2Program, "size_param"), size_param );

  //
  // fragment program
  //
  cgGLEnableProfile( fragProfile );
  cgGLBindProgram( frag2Program );

#ifdef RENDER_FBO
  pbuffer->Bind(0);
  cgGLSetTextureParameter( cgGetNamedParameter(frag2Program, "qhatT"), pbuffer->tex(0) );
#endif

  cgGLEnableTextureParameter( cgGetNamedParameter(frag2Program, "qhatT") );

  cgGLEnableTextureParameter( cgGetNamedParameter(frag2Program, "polyParam") );
  cgGLEnableTextureParameter( cgGetNamedParameter(frag2Program, "gaussKernel") );


  cgGLSetStateMatrixParameter( cgGetNamedParameter(frag2Program, "ModelViewInv"),
			       CG_GL_MODELVIEW_MATRIX,
			       CG_GL_MATRIX_INVERSE );

  cgGLSetStateMatrixParameter(cgGetNamedParameter(frag2Program, "ModelViewProjInv"),
			      CG_GL_MODELVIEW_PROJECTION_MATRIX,
			      CG_GL_MATRIX_INVERSE );

//   GLint viewport[4];
//   glGetIntegerv( GL_VIEWPORT, viewport );

//   float viewParam[4];
  viewParam[0] = (float) viewport[0];
  viewParam[1] = (float) viewport[1];
  viewParam[2] = (float) viewport[2];
  viewParam[3] = (float) viewport[3];
  cgGLSetParameter4fv( cgGetNamedParameter(frag2Program, "viewParam"), viewParam );

  // set draw mode ... 0: Billboard, 1; Point Sprite
  cgGLSetParameter1f( cgGetNamedParameter(frag2Program, "drawMode"), draw_mode );

  // Use Multiple Draw Buffers
  // set destination buffers
  glDrawBuffersATI( sizeof(buffers) / sizeof(GLenum), buffers );


#define DRAW_BLEND 1

  //
  // draw quads
  //

#ifdef DRAW_BLEND
  glDisable ( GL_DEPTH_TEST );
  glEnable( GL_BLEND );
  glBlendFunc ( GL_SRC_ALPHA, GL_ONE );
  glDepthMask( GL_FALSE );
#endif

  drawSlimTree( size_param );

#ifdef DRAW_BLEND
  glEnable ( GL_DEPTH_TEST );
  glDepthMask( GL_TRUE );
  glDisable( GL_BLEND );
#endif

  cgGLDisableTextureParameter( cgGetNamedParameter(frag2Program, "qhatT") );
  pbuffer->Release(0);

  cgGLDisableTextureParameter( cgGetNamedParameter(frag2Program, "gaussKernel") );
  cgGLDisableTextureParameter( cgGetNamedParameter(frag2Program, "polyParam") );

  cgGLDisableProfile( fragProfile );

  cgGLDisableProfile( vertProfile );

  pane.finish();

#ifdef OGL_QUERY
  glEndQueryARB( GL_SAMPLES_PASSED_ARB );
#endif
  
#if !defined(SLIM_DEBUG_II)
  
  //
  // copy fpbuffer (buffer #0) to a position texture
  //

#ifndef RENDER_FBO

  glReadBuffer( buffers[0] );

  glBindTexture( GL_TEXTURE_RECTANGLE_EXT, positionT );
  glCopyTexSubImage2D( GL_TEXTURE_RECTANGLE_EXT, 0, 0, 0, 0, 0, pbuffer_width, pbuffer_height );
  cgGLSetTextureParameter( cgGetNamedParameter(frag3Program, "positionT"), positionT );

#else // RENDER_FBO

#endif // RENDER_FBO

  //
  // copy fpbuffer (buffer #1) to a normal texture
  //

#ifndef RENDER_FBO

  glReadBuffer( buffers[1] );

  glBindTexture( GL_TEXTURE_RECTANGLE_EXT, normalT );
  glCopyTexSubImage2D( GL_TEXTURE_RECTANGLE_EXT, 0, 0, 0, 0, 0, pbuffer_width, pbuffer_height );
  cgGLSetTextureParameter( cgGetNamedParameter(frag3Program, "normalT"), normalT );

#else // RENDER_FBO

#endif // RENDER_FBO

#endif // !defined(SLIM_DEBUG_II)

#if defined(SLIM_DEBUG_II)

  point_data.clear(); point_data.resize( pbuffer_width*pbuffer_height*4 );
  glReadBuffer( buffers[0] );
  glReadPixels(0, 0, pbuffer_width, pbuffer_height, GL_RGBA, GL_FLOAT, &point_data[0] );

  char filename[BUFSIZ];
  sprintf( filename, "tmpII%03d.pnt", count_frame );
  std::ofstream ofs( filename );
  int pixel_count = 0;
  float max_val;

  for ( int i = 0; i < pbuffer_width*pbuffer_height; ++i )
    {
      if ( (point_data[4*i] == 0.0) && (point_data[4*i+1] == 0.0) && (point_data[4*i+2] == 0.0) )
	{
	}
      else
	{
	  ofs << point_data[4*i] /point_data[4*i+3] << " " 
	      << point_data[4*i+1] / point_data[4*i+3] << " "
	      << point_data[4*i+2] / point_data[4*i+3] << endl;

	  if ( pixel_count )
	    {
	      if ( max_val < point_data[4*i] )
		max_val = point_data[4*i];
	      if ( max_val < point_data[4*i+1] )
		max_val = point_data[4*i+1];
	    }
	  else
	    {
	      max_val = point_data[4*i];
	      if ( max_val < point_data[4*i+1] )
		max_val = point_data[4*i+1];
	    }

	  ++pixel_count;
	}
    }

  cout << "max val " << max_val << endl;
  cout << "pixel count " << pixel_count << endl;

  ++count_frame;

  ofs.close();

#endif // defined(SLIM_DEBUG_II)

  pbuffer16->Deactivate();

//   glPointParameterfARB( GL_POINT_SPRITE_COORD_ORIGIN, GL_UPPER_LEFT ); 
  
#endif // #ifndef SLIM_DEBUG

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // 3rd pass
  //
  ////////////////////////////////////////////////////////////////////////////////////

#if !defined(SLIM_DEBUG) && !defined(SLIM_DEBUG_II)

  if ( isCapture ) cap_buffer->Activate();
    
  pane.clear( pbuffer_width, pbuffer_height );

#if !defined(RENDER_FBO)

#else // RENDER_FBO

  pane.setView();

#endif // RENDER_FBO

  cgGLEnableProfile( fragProfile );
  cgGLBindProgram( frag3Program );

  // matrices
  cgGLSetStateMatrixParameter( cgGetNamedParameter(frag3Program, "ModelViewProj"),
			       CG_GL_MODELVIEW_PROJECTION_MATRIX,
			       CG_GL_MATRIX_IDENTITY );
  
  cgGLSetStateMatrixParameter( cgGetNamedParameter(frag3Program, "ModelView"),
			       CG_GL_MODELVIEW_MATRIX,
			       CG_GL_MATRIX_IDENTITY );
  
  cgGLSetStateMatrixParameter( cgGetNamedParameter(frag3Program, "ModelViewIT"),
			       CG_GL_MODELVIEW_MATRIX,
			       CG_GL_MATRIX_INVERSE_TRANSPOSE );
  
#if !defined(RENDER_FBO)

#else // RENDER_FBO

  pane.initView();

#endif // RENDER_FBO

  // background color
  cgGLSetParameter4fv( cgGetNamedParameter(frag3Program, "bgColor"), bgColor );

  // materials
  cgGLSetParameter4f( g_CGparam_AmbMat,   myMatl[0],  myMatl[1],  myMatl[2],  myMatl[3] );
  cgGLSetParameter4f( g_CGparam_DiffMat,  myMatl[4],  myMatl[5],  myMatl[6],  myMatl[7] );
  cgGLSetParameter4f( g_CGparam_EmisMat,  myMatl[8],  myMatl[9],  myMatl[10], myMatl[11] );
  cgGLSetParameter4f( g_CGparam_SpecMat,  myMatl[12], myMatl[13], myMatl[14], myMatl[15] );
  cgGLSetParameter4f( g_CGparam_ShineMat, myMatl[16], 0.0f,       0.0f,       0.0f );
  
  // light
  float current_light[3];
  pane.getRealLightPosition( current_light );
  cgGLSetParameter4f( g_CGparam_LightPos,  current_light[0], current_light[1],  
		      current_light[2], 1.0f  );
  cgGLSetParameter4f( g_CGparam_LightVec,  -current_light[0], -current_light[1], 
		      -current_light[2], 1.0f );
  cgGLSetParameter4f( g_CGparam_AmbLight,  myLight[4],  myLight[5],  myLight[6],  myLight[7] );
  cgGLSetParameter4f( g_CGparam_DiffLight, myLight[8],  myLight[9],  myLight[10], myLight[11] );
  cgGLSetParameter4f( g_CGparam_SpecLight, myLight[12], myLight[13], myLight[14], myLight[15] );
  cgGLSetParameter3f( g_CGparam_AttenLight, myLight[16], myLight[17], myLight[18] );  // constant, linear, quadratic
  cgGLSetParameter2f( g_CGparam_SpotLight, cosf(180.0f), 0.0f ); //cosf(spot_cuttoff_angle), spot power

#ifdef RENDER_FBO
  pbuffer16->Bind( 0 );
  cgGLSetTextureParameter( cgGetNamedParameter(frag3Program, "positionT"), pbuffer16->tex(0) );
#endif

  cgGLEnableTextureParameter( cgGetNamedParameter(frag3Program, "positionT") );
  
#ifdef RENDER_FBO
  pbuffer16->Bind( 1 );
  cgGLSetTextureParameter( cgGetNamedParameter(frag3Program, "normalT"), pbuffer16->tex(1) );
#endif

  cgGLEnableTextureParameter( cgGetNamedParameter(frag3Program, "normalT") );
  
  {
    drawQuad( pbuffer_width, pbuffer_height );
  }

  cgGLDisableTextureParameter( cgGetNamedParameter(frag3Program, "normalT") );
  cgGLDisableTextureParameter( cgGetNamedParameter(frag3Program, "positionT") );
  pbuffer16->Release(0);
  pbuffer16->Release(1);
  
  cgGLDisableProfile( fragProfile );

#endif  // #ifndef SLIM_DEBUG

  ////////////////////////////////////////////////////////////////////////////////////

#if defined(SLIM_DEBUG) || defined(SLIM_DEBUG_II)

  // for debug
  pane.setBackgroundColor( 1.0f, 1.0f, 1.0f );
  pane.setIsGradientBackground( true );
  pane.clear( pbuffer_width, pbuffer_height );
  pane.setView();
  
  glColor3f( 1.0f, .0f, .0f );
  glBegin( GL_POINTS );
  for ( int i = 0; i < pbuffer_width*pbuffer_height; ++i )
    {
      glVertex3f( point_data[4*i], point_data[4*i+1], point_data[4*i+2] );
    }
  glEnd();

#endif // #ifdef SLIM_DEBUG

  ////////////////////////////////////////////////////////////////////////////////////

  pane.finish();

#ifdef OGL_QUERY
  GLuint sl1, sl2;
  glGetQueryObjectuivARB( oq_slim1, GL_QUERY_RESULT_ARB, &sl1 );
  glGetQueryObjectuivARB( oq_slim2, GL_QUERY_RESULT_ARB, &sl2 );

  cout << "occlusion query ... 1st pass: " << sl1 << " 2nd pass: " << sl2 << endl;
#endif
}

void display()
{
  MakeGlutWindowCurrent();

  pane.setBackgroundColor( bgColor[0], bgColor[1], bgColor[2] );
  pane.clear( pbuffer_width, pbuffer_height );
  pane.initView();

  switch ( render_mode )
    {
    case RENDER_GPU:
      displaySLIMGPU();
      break;

    case RENDER_BALL:
#ifdef OLD_SLIM
      displayBall();
#endif // OLD_SLIM
#ifdef NEW_SLIM
      displaySlimTreeBall();
#endif // NEW_SLIM
      break;
      
    case RENDER_NONE:
      displayPoly();
      break;
      
    default:
      break;
    }

  checkGLErrors("display");

  fps.frame();
  if ( fps.timing_updated() )
    {
      float f = fps.get_fps();
      if ( max_fps < f ) max_fps = f;
      sprintf( buf,"%.1f fps max: %.1f #balls: %d",
               f, max_fps, slimball_count );
    }

  char str0[BUFSIZ], str1[BUFSIZ];
  if ( draw_mode == DRAW_BILLBOARD )
    strcpy( str0, "Billboard" );
  else if ( draw_mode == DRAW_POINTSPRITE )
    strcpy( str0, "Point sprite" );
  if ( lod_mode == LOD_NONE )
    strcpy( str1, "Leaf" );
  else 
    strcpy( str1, "LOD" );
      
  sprintf( txt, "%s - %s - %s", str0, str1, buf );
  ::glutSetWindowTitle( txt );  

  if ( isCapture )
    {
      PNGImage pi( cap_width, cap_height, false ); 
      pi.capture_and_write( "save_screen.png" );
      isCapture = false;
      cap_buffer->Deactivate();
      ::glutPostRedisplay();
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
      int mod = glutGetModifiers();
      if ( button == GLUT_LEFT_BUTTON && (mod != GLUT_ACTIVE_SHIFT) && (mod != GLUT_ACTIVE_CTRL) ) 
	{
	  pane.startRotate();
	}
      else if ( (button == GLUT_LEFT_BUTTON) && (mod == GLUT_ACTIVE_SHIFT) ) 
	{
	  pane.startZoom();
	}
      else if ( (button == GLUT_LEFT_BUTTON) && (mod == GLUT_ACTIVE_CTRL) ) 
	{
	  pane.startMove();
	}
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

void reshapeFunc( int w, int h )
{
  glViewport( 0, 0, w, h );

  pane.setW( w );
  pane.setH( h );

  destroyRenderTextures();
  destroyFPBuffer();

  width = w; height = h;
  pbuffer_width = w; pbuffer_height = h;
  cap_width = w; cap_height = h;

  initFPBuffer( pbuffer_width, pbuffer_height, cap_width, cap_height );
  initRenderTexture( pbuffer_width, pbuffer_height );

  MakeGlutWindowCurrent();
  pane.initGL();
}

void keyboard( unsigned char c, int x, int y )
{
  VWIO vw_out;

  b[c] = ! b[c];

  switch ( c )
    {
    case 'q':

      cout << "max fps: " << max_fps << endl;
      destroyRenderTextures();
      destroyFPBuffer();
      destroyCg();
      exit(0);

      break;

    case 'v':

      cout << "output to .vw file ... " << endl;
      vw_out.outputToFile( "tmp.vw", pane.manip() );
      cout << "done." << endl;
      
      break;

    case 'p':

      // screen capture
      isCapture = true;
      
      break;

    case 'm':

      if ( lod_mode == LOD_ERROR )
	{
	  slimball_error *= 0.1;
	  if ( slimball_error < 1.0e-15 ) slimball_error = 1.0e-15;
	  cout << "slimball error threshold: " << slimball_error << endl;
	}
      else if ( lod_mode == LOD_LEVEL )
	{
	  ++(slimball_level);
	  cout << "slimball level: " << slimball_level << endl;
	}
      break;

    case 'n':

      if ( lod_mode == LOD_ERROR )
	{
	  slimball_error *= 10.0;
	  if ( slimball_error > 1.0 ) slimball_error = 1.0;
	  cout << "slimball error threshold: " << slimball_error << endl;
	}
      else if ( lod_mode == LOD_LEVEL )
	{
	  --(slimball_level);
	  if ( slimball_level < 0 ) slimball_level = 0;
	  cout << "slimball level: " << slimball_level << endl;
	}
      
      break;

    case '1':
      render_mode = RENDER_GPU;
      draw_mode   = DRAW_BILLBOARD;
      lod_mode = LOD_NONE;
      break;
      
    case '2':
      render_mode = RENDER_GPU;
      draw_mode   = DRAW_POINTSPRITE;
      lod_mode = LOD_NONE;
      break;

    case '3':
      render_mode = RENDER_GPU;
      draw_mode   = DRAW_BILLBOARD;
      lod_mode = LOD_LEVEL;
      break;

    case '4':
      render_mode = RENDER_GPU;
      draw_mode   = DRAW_POINTSPRITE;
      lod_mode = LOD_LEVEL;
      break;

    case '5':
      render_mode = RENDER_GPU;
      draw_mode   = DRAW_BILLBOARD;
      lod_mode = LOD_ERROR;
      break;

    case '6':
      render_mode = RENDER_BALL;
      lod_mode = LOD_NONE;
      break;

    case '7':
      render_mode = RENDER_BALL;
      lod_mode = LOD_ERROR;
      break;

    case '8':
      render_mode = RENDER_NONE;
      lod_mode = LOD_NONE;
      break;
    }
  
  ::glutPostRedisplay();
}

void main_menu( int i )
{
  keyboard( (unsigned char) i, 0, 0 );
}

void init_menu()
{
  int draw_menu = glutCreateMenu( main_menu );
  glutAddMenuEntry( "SLIM with Billboard [1]", '1' );
  glutAddMenuEntry( "SLIM with Point Sprite [2]", '2' );
  glutAddMenuEntry( "SLIM with Billboard LOD [3]", '3' );
  glutAddMenuEntry( "SLIM with Point Sprite LOD [4]", '4' );
  glutAddMenuEntry( "Support Ball [5]", '5' );
 
  glutCreateMenu( main_menu );
  glutAddSubMenu("Draw...", draw_menu );
  glutAddMenuEntry( "Quit [q]", 'q' );
  glutAttachMenu( GLUT_RIGHT_BUTTON );
}

void idle()
{
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

#ifdef OLD_SLIM
  // initialize old SLIM
  listN = filemanager.readBallFile( filename, lists, Ns, deg, oriented );
  imagemaker.setBfLists( deg, listN, lists, Ns );
  imagemaker.setBoth( !oriented );
#endif // OLD_SLIM

#ifdef NEW_SLIM
  SlimTreefIO slimtreeIO( slimtree );
  slimtreeIO.inputFromFile( filename );
#endif // NEW_SLIM
  
  ::glutInitWindowSize( width, height );
  ::glutInit( &argc, argv );
  //::glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_DEPTH | GLUT_STENCIL );
  ::glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
  ::glutCreateWindow( argv[0] );

  ::glutDisplayFunc( display );
  ::glutKeyboardFunc( keyboard );
  ::glutIdleFunc( idle );
  ::glutMouseFunc( mouse );
//   ::glutReshapeFunc( reshapeFunc );
  ::glutMotionFunc( motion );
//   init_menu();
  
  pane.init( width, height );
  //pane.setIsLightOn( 1, false );

#ifdef OLD_SLIM
  pane.setViewPoint( Point3f( .0f, .0f, 55.0f ) );
  pane.setViewVector( Vector3f( .0f, .0f, -55.0f ) ); 
  pane.setMagObject( 10.0f );
#endif // OLD_SLIM

#ifdef NEW_SLIM
  pane.setViewPoint( Point3f( .0f, .0f, 30.0f ) );
  pane.setViewVector( Vector3f( .0f, .0f, -30.0f ) ); 
  pane.setMagObject( 10.0f );
#endif // NEW_SLIM

  // Initialize float pbuffer, Cg, Texture
  pane.initGLEW();
  initFPBuffer( pbuffer_width, pbuffer_height, cap_width, cap_height );
  initRenderTexture( pbuffer_width, pbuffer_height );
  initCg();
  setGaussKernelTexture();

#ifdef OLD_SLIM
  setPolyTexture();
#endif // OLD_SLIM

#ifdef NEW_SLIM
  setSlimTreePolyTexture();
#endif // NEW_SLIM

  // Initialize some state for the GLUT window's rendering context.
  MakeGlutWindowCurrent();
  pane.initGL();
  pane.setIsLightOn( 1, false );
  pane.setLightParameters( 0, myLight );

#if 0  
  if ( argc >= 3 )
    {
      SMFRIO rio_input( meshR );
      rio_input.inputFromFile( argv[2] );
      meshR.normalize();
      meshR.scale( 10.0 );
      if ( meshR.normals().empty() ) meshR.createVertexNormals();
      
      glmeshR.setMesh( meshR );
      glmeshR.setIsSmoothShading( true );
      glmeshR.setMaterial( myMatl );
    }
#endif
  
  if ( argc == 3 )
    {
      VWIO vw_in;
      vw_in.inputFromFile( argv[2], pane.manip() );
    }

  ::glutMainLoop();
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////

