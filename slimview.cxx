////////////////////////////////////////////////////////////////////
//
// $Id: slimview.cxx 2025/09/06 16:20:08 kanai Exp 
//
// Copyright (c) 2021-2025 Takashi Kanai
// Released under the MIT license
//
////////////////////////////////////////////////////////////////////

#include "envDep.h"

#include "GL/glew.h"
#include "GL/wglew.h"

#include <GLFW/glfw3.h>

// #include <GL/gl.h>
// #include <GL/glu.h>
// #include <GL/glext.h>

#include <Cg/cg.h>
#include <Cg/cgGL.h>

#include <vector>
#include <iomanip>
using namespace std;

#include <Point2.h>
#include <Point3.h>
#include <Vector3.h>
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif  // VM_INCLUDE_NAMESPACE

#include "mydef.h"

#ifndef RENDER_FBO

#include "RenderTexture.h"

#else  // RENDER_FBO

#include "RenderTextureFBO.h"

#endif  // RENDER_FBO

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
#endif  // OLD_SLIM

// new slim format
#ifdef NEW_SLIM
#include "SlimTree.hxx"
#include "SlimBall.hxx"
#include "BasisFunction.hxx"
#include "SlimTreeIO.hxx"
#endif  // NEW_SLIM

// #include "shader_init.hxx"
// ShaderPipelines shaders;

////////////////////////////////////////////////////////////////////////////////////

// keyboard
bool shift_key_pressed = false;
bool control_key_pressed = false;
// mouse
bool left_button_pressed = false;
bool right_button_pressed = false;

// int width = 512;
// int height = 512;
int width = 800;
int height = 800;
//int width = 1024;
//int height = 1024;
// int width = 1280;
// int height = 1024;

float bgColor[4] = {1.0f, 1.0f, 1.0f, 1.0f};

static float myLight[] = {
    0.0f, 0.0f, 100.0f, 1.0f, 0.2f, 0.2f, 0.2f, 1.0f,
    1.0f, 1.0f, 1.0f,   1.0f, 0.8f, 0.8f, 0.8f, 1.0f,
    1.0f,  // 0.1f,
    0.0f,  // 0.05f
    0.0f   // 0.05f
};

static float myMatl[] = {0.2f, 0.2f, 0.2f, 1.0f, 0.8f, 0.8f, 0.8f, 1.0f, 0.0f,
                         0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 80.0f};

static float myBallLight[] = {
    0.0f, 0.0f, 100.0f, 1.0f, 0.2f, 0.2f, 0.2f, 1.0f,
    1.0f, 1.0f, 1.0f,   1.0f, 0.8f, 0.8f, 0.8f, 1.0f,
    1.0f,  // 0.1f,
    0.0f,  // 0.05f
    0.0f   // 0.05f
};

static float myBallMatl[] = {0.2f, 0.2f, 0.2f, 1.0f, 0.9f, 0.5f,
                             0.5f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                             1.0f, 1.0f, 1.0f, 1.0f, 80.0f};

GLPanel pane;

MeshR meshR;
GLMeshR glmeshR;
GLMaterial ballMtl(myBallMatl);

////////////////////////////////////////////////////////////////////////////////////

#include "c11timer.hxx"
C11Timer c11fps;
double max_c11fps = 0.0;

// timer fps(10);
// float max_fps = .0f;
char buf[BUFSIZ];
char txt[BUFSIZ];
int slimball_count = 0;
int slimball_level = 10;
double slimball_error = 1.0e-05;

bool b[256];

#define RENDER_GPU 0
#define RENDER_BALL 1
#define RENDER_NONE 2
int render_mode = RENDER_GPU;

#define LOD_NONE 0
#define LOD_LEVEL 1
#define LOD_ERROR 2
int lod_mode = LOD_NONE;
#ifdef OLD_SLIM
bool isLOD = false;
#endif

// int draw_mode = DRAW_BILLBOARD;
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
int pbuffer_width = width;
int pbuffer_height = height;

#ifndef RENDER_FBO

const GLenum buffers[] = {
    GL_AUX0,
    GL_AUX1,
};

RenderTexture* pbuffer = NULL;
RenderTexture* pbuffer16 = NULL;

#else  // RENDER_FBO

const GLenum buffers[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT};

RenderTextureFBO* pbuffer = NULL;

GLenum texTarget0 = GL_TEXTURE_RECTANGLE_EXT;
GLenum texInternalFormat0 = GL_RGBA32F_ARB;
GLenum texFormat0 = GL_RGBA;
GLenum texType0 = GL_FLOAT;

RenderTextureFBO* pbuffer16 = NULL;

GLenum texTarget1 = GL_TEXTURE_RECTANGLE_EXT;
// GLenum texInternalFormat1 = GL_RGBA32F_ARB;
GLenum texInternalFormat1 = GL_RGBA16F_ARB;
GLenum texFormat1 = GL_RGBA;
GLenum texType1 = GL_FLOAT;

#endif  // RENDER_FBO

//
// for screen capture
//
int cap_width = width;
int cap_height = height;

#ifndef RENDER_FBO

RenderTexture* cap_buffer = NULL;

#else  // RENDER_FBO

RenderTextureFBO* cap_buffer = NULL;

GLenum texTarget2 = GL_TEXTURE_2D;
GLenum texInternalFormat2 = GL_RGB;
GLenum texFormat2 = GL_RGB;
GLenum texType2 = GL_UNSIGNED_BYTE;

#endif  // RENDER_FBO

bool isCapture = false;

////////////////////////////////////////////////////////////////////////////////////
//
// Ball Color
//

#include "ColorConv.hxx"

static float max_dis = 1.0;
static float min_dis = 1.0e-10;
static Point3f max_col(1.0f, 0.0f, 0.0f);
static Point3f min_col(0.0f, 0.0f, 1.0f);

void convColor(float dis, float* color) {
  ColorConv cc;
  Point3f min_col_hsv;
  Point3f max_col_hsv;
  cc.rgbtohsv(max_col, max_col_hsv);
  cc.rgbtohsv(min_col, min_col_hsv);

  float max_dis_log = (float)std::log10(max_dis);
  float min_dis_log = (float)std::log10(min_dis);

  float dis_log = (float)std::log10(dis);

  float t;
  if (dis > max_dis) {
    t = 1.0f;
  } else if (dis < min_dis) {
    t = 0.0f;
  } else {
    t = (dis_log - min_dis_log) / (max_dis_log - min_dis_log);
  }

  Point3f col_hsv;
  col_hsv.interpolate(min_col_hsv, max_col_hsv, t);
  Point3f col;
  cc.hsvtorgb(col_hsv, col);
  color[0] = col.x;
  color[1] = col.y;
  color[2] = col.z;
}

////////////////////////////////////////////////////////////////////////////////////
//
// Draw sphere (replacement for glutSolidSphere)
//

// Create the unit sphere mesh once (client arrays, no VBO/IBO)
struct SphereMesh {
  std::vector<float>    interleaved;  // [px py pz nx ny nz u v]*
  std::vector<uint32_t> indices;      // GL_UNSIGNED_INT
  int triCount = 0;
};
SphereMesh sm;

static SphereMesh makeUnitSphereClient(int n1, int n2)
{
  SphereMesh m;
  if (n1 < 2 || n2 < 3) return m;

  const float PI = 3.14159265358979323846f;
  const int rows = n1 + 1;
  const int cols = n2 + 1;

  m.interleaved.reserve(rows * cols * 8);
  m.indices.reserve(n1 * n2 * 6);

  for (int i = 0; i < rows; ++i) {
    float v   = static_cast<float>(i) / n1;   // [0,1]
    float phi = PI * v;                        // [0,\pi]
    float sphi = std::sin(phi), cphi = std::cos(phi);
    for (int j = 0; j < cols; ++j) {
      float u     = static_cast<float>(j) / n2; // [0,1]
      float theta = 2.0f * PI * u;              // [0,2\pi]
      float st = std::sin(theta), ct = std::cos(theta);

      // 単位球の法線=位置
      float nx = sphi * ct, ny = cphi, nz = sphi * st;
      float px = nx,        py = ny,        pz = nz; // Location on the sphere at radius r = 1

      m.interleaved.push_back(px);
      m.interleaved.push_back(py);
      m.interleaved.push_back(pz);
      m.interleaved.push_back(nx);
      m.interleaved.push_back(ny);
      m.interleaved.push_back(nz);
      m.interleaved.push_back(u);
      m.interleaved.push_back(v);
    }
  }

  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      uint32_t i0 = i * cols + j;
      uint32_t i1 = i * cols + (j + 1);
      uint32_t i2 = (i + 1) * cols + j;
      uint32_t i3 = (i + 1) * cols + (j + 1);
      // CCW
      m.indices.push_back(i0); m.indices.push_back(i2); m.indices.push_back(i1);
      m.indices.push_back(i1); m.indices.push_back(i2); m.indices.push_back(i3);
    }
  }
  m.triCount = static_cast<int>(m.indices.size() / 3);
  return m;
}

// Render with variable radius r (using matrix scaling)
static void drawSphereClientWithRadius(const SphereMesh& m, float r)
{
  if (m.triCount == 0 || r <= 0.0f) return;

  const GLsizei stride = 8 * sizeof(float);
  const GLvoid* base   = static_cast<const GLvoid*>(m.interleaved.data());

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);

  glVertexPointer  (3, GL_FLOAT, stride, reinterpret_cast<const GLvoid*>(
                      reinterpret_cast<uintptr_t>(base) + 0));
  glNormalPointer  (   GL_FLOAT, stride, reinterpret_cast<const GLvoid*>(
                         reinterpret_cast<uintptr_t>(base) + 3 * sizeof(float)));
  glTexCoordPointer(2, GL_FLOAT, stride, reinterpret_cast<const GLvoid*>(
                      reinterpret_cast<uintptr_t>(base) + 6 * sizeof(float)));

  // ---- Apply radius r here ----
  glPushMatrix();
  glScalef(r, r, r);   // Multiply position by r (normals require separate handling)

  // Handling normals:
  // 1) For uniform scaling, use GL_RESCALE_NORMAL (cheaper)
  // 2) If compatibility is a concern, use GL_NORMALIZE (slightly more expensive)
  // glEnable(GL_RESCALE_NORMAL);   // OpenGL 1.2+ (assumes uniform scaling)
  glEnable(GL_NORMALIZE);           // Use this if you want guaranteed normalization

  glDrawElements(GL_TRIANGLES,
                 static_cast<GLsizei>(m.indices.size()),
                 GL_UNSIGNED_INT,
                 m.indices.data());

  glDisable(GL_NORMALIZE);
  // glDisable(GL_RESCALE_NORMAL);
  glPopMatrix();
  // --------------------------------

  glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);
}

#if 0
void DrawSphere(float r, int n1, int n2) {
  if (r <= 0.0f || n1 < 2 || n2 < 3) return;

  const float PI = 3.14159265358979323846f;
  const float dPhi = PI / static_cast<float>(n1);           // Latitude angle (0..\pi)
  const float dTheta = 2.0f * PI / static_cast<float>(n2);  // Longitude angle (0..2\pi)

  for (int i = 0; i < n1; ++i) {
    float phi0 = i * dPhi;
    float phi1 = (i + 1) * dPhi;

    // Render each latitude band as a triangle strip
    glBegin(GL_TRIANGLE_STRIP);
    for (int j = 0; j <= n2; ++j)
    {
      float theta = j * dTheta;

      // Points of the upper ring (phi0)
      float x0 = std::sin(phi0) * std::cos(theta);
      float y0 = std::cos(phi0);
      float z0 = std::sin(phi0) * std::sin(theta);

      // Points of the lower ring (phi1)
      float x1 = std::sin(phi1) * std::cos(theta);
      float y1 = std::cos(phi1);
      float z1 = std::sin(phi1) * std::sin(theta);

      // Normal = position on the unit sphere (for smooth shading)
      glNormal3f(x1, y1, z1);
      glVertex3f(r * x1, r * y1, r * z1);

      glNormal3f(x0, y0, z0);
      glVertex3f(r * x0, r * y0, r * z0);
    }
    glEnd();
  }
}
#endif

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
CGprogram frag1Program;  // 1st
CGprogram frag2Program;  // 2nd
CGprogram frag3Program;  // 3rd
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
GLuint qhatT;      // 1st pass
GLuint positionT;  // 2nd pass
GLuint normalT;    // 2nd pass
#endif             // RENDER_FBO

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
ImageMaker imagemaker(pbuffer_width, pbuffer_height);
int* Ns;
int deg;
bool oriented;
int listN;
#endif  // OLD_SLIM

// new slim
#ifdef NEW_SLIM
SlimTreef slimtree;
#endif  // NEW_SLIM

////////////////////////////////////////////////////////////////////////////////////

void checkGLErrors(char* s) {
  GLenum error;
  while ((error = glGetError()) != GL_NO_ERROR) {
    fprintf(stderr, "%s: error - %s\n", s, (char*)gluErrorString(error));
  }
}

static void handleCgError() {
  fprintf(stderr, "Cg error: %s\n", cgGetErrorString(cgGetError()));
}

////////////////////////////////////////////////////////////////////////////////////

void initFPBuffer(int pW, int pH, int cW, int cH) {
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

  if (pbuffer) delete pbuffer;

#ifndef RENDER_FBO

  pbuffer = new RenderTexture("float=32 rgba depth textureRECT", pW, pH,
                              GL_TEXTURE_RECTANGLE_EXT);

#else  // RENDER_FBO

  pbuffer = new RenderTextureFBO(pW, pH, texTarget0, texInternalFormat0,
                                 texFormat0, texType0, 1);

#endif  // RENDER_FBO

  pbuffer->Activate();
  pane.initGL();
  pbuffer->Deactivate();

  //
  // 16bit floating-point buffer used in the second pass (including blending
  // operations)
  //

  if (pbuffer16) delete pbuffer16;

#ifndef RENDER_FBO

  pbuffer16 = new RenderTexture("float=16 aux=2 rgba textureRECT", pW, pH,
                                GL_TEXTURE_RECTANGLE_EXT);

#else  // RENDER_FBO

  pbuffer16 = new RenderTextureFBO(pW, pH, texTarget1, texInternalFormat1,
                                   texFormat1, texType1, 2);

#endif  // RENDER_FBO

  pbuffer16->Activate();
  pane.initGL();
  pbuffer16->Deactivate();

  if (cap_buffer) delete cap_buffer;

#ifndef RENDER_FBO

  cap_buffer =
      new RenderTexture("rgba textureRECT", cW, cH, GL_TEXTURE_RECTANGLE_EXT);

#else  // RENDER_FBO

  cap_buffer = new RenderTextureFBO(pW, pH, texTarget2, texInternalFormat2,
                                    texFormat2, texType2, 1);

#endif  // RENDER_FBO

  cap_buffer->Activate();
  pane.initGL();
  cap_buffer->Deactivate();
}

void destroyFPBuffer() {
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
void initRenderTexture(int pW, int pH) {
#ifndef RENDER_FBO
  // q_hat
  glGenTextures(1, &qhatT);
  glBindTexture(GL_TEXTURE_RECTANGLE_EXT, qhatT);
  glTexImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, GL_FLOAT_RGBA_NV, pW, pH, 0,
               GL_RGBA, GL_FLOAT, NULL);

  // position
  glGenTextures(1, &positionT);
  glBindTexture(GL_TEXTURE_RECTANGLE_EXT, positionT);
  glTexImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, GL_FLOAT_RGBA_NV, pW, pH, 0,
               GL_RGBA, GL_FLOAT, NULL);

  // normal vector
  glGenTextures(1, &normalT);
  glBindTexture(GL_TEXTURE_RECTANGLE_EXT, normalT);
  glTexImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, GL_FLOAT_RGBA_NV, pW, pH, 0,
               GL_RGBA, GL_FLOAT, NULL);
#endif  // RENDER_FBO

  // initialize GL_ARB_point_sprite
  float maxSize;
  glGetFloatv(GL_POINT_SIZE_MAX_ARB, &maxSize);
  cout << "max point sprite size: " << maxSize << endl;

#ifdef OGL_QUERY
  GLint bitsSupported;
  glGetQueryivARB(GL_SAMPLES_PASSED_ARB, GL_QUERY_COUNTER_BITS_ARB,
                  &bitsSupported);
  cout << "Number of counter bits = " << bitsSupported << endl;

  glGenQueriesARB(1, &oq_slim1);
  glGenQueriesARB(1, &oq_slim2);
#endif
}

void destroyRenderTextures() {
#ifndef RENDER_FBO
  glDeleteTextures(1, &qhatT);
  glDeleteTextures(1, &positionT);
  glDeleteTextures(1, &normalT);
#endif  // RENDER_FBO
}

////////////////////////////////////////////////////////////////////////////////////

void initCg() {
  cgSetErrorCallback(handleCgError);
  context = cgCreateContext();

  // profiles
  vertProfile = cgGLGetLatestProfile(CG_GL_VERTEX);
  cgGLSetOptimalOptions(vertProfile);
  fragProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
  cgGLSetOptimalOptions(fragProfile);

  // vertex program (1st-pass)
  vert1Program = cgCreateProgramFromFile(context, CG_SOURCE, "shader_vp_1st.cg",
                                         vertProfile, NULL, NULL);
  cgGLLoadProgram(vert1Program);

  // fragment program (1st-pass)
  frag1Program = cgCreateProgramFromFile(context, CG_SOURCE, "shader_fp_1st.cg",
                                         fragProfile, NULL, NULL);
  cgGLLoadProgram(frag1Program);

  // vertex program (2nd-pass)
  vert2Program = cgCreateProgramFromFile(context, CG_SOURCE, "shader_vp_2nd.cg",
                                         vertProfile, NULL, NULL);
  cgGLLoadProgram(vert2Program);

  // fragment program (2nd-pass)
  frag2Program = cgCreateProgramFromFile(context, CG_SOURCE, "shader_fp_2nd.cg",
                                         fragProfile, NULL, NULL);
  cgGLLoadProgram(frag2Program);

  // fragment program (3rd-pass)
  frag3Program = cgCreateProgramFromFile(context, CG_SOURCE, "shader_fp_3rd.cg",
                                         fragProfile, NULL, NULL);
  cgGLLoadProgram(frag3Program);

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

void destroyCg() {
  cgDestroyProgram(vert1Program);
  cgDestroyProgram(vert2Program);
  cgDestroyProgram(frag1Program);
  cgDestroyProgram(frag2Program);
  cgDestroyProgram(frag3Program);
  //   cgDestroyProgram( debugProgram ); // for debug
  cgDestroyContext(context);
}

////////////////////////////////////////////////////////////////////////////////////

#ifdef OLD_SLIM
void setPolyTexture() {
  GLuint polyParam;
  glGenTextures(1, &polyParam);
  glBindTexture(GL_TEXTURE_RECTANGLE_EXT, polyParam);
  glTexParameteri(GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

  GLsizei sizeW;
  glGetIntegerv(GL_MAX_RECTANGLE_TEXTURE_SIZE_EXT, &sizeW);

  // Store coefficients of polynomial to a texture
  // Each polynomial uses 4 pixels to store coefficients and a center point
  std::vector<float> pp;
#if 0
  int max_width = 512;
  int tex_width  = max_width * 4; // pixel width
#endif
  int tex_width = sizeW;
  int max_width = (int)(tex_width / 4.0);
  cout << "max slimball size: " << max_width * sizeW << endl;

  int flo_width = tex_width * 4;
  for (int k = 0; k < flo_width; ++k) pp.push_back(0.0f);
  int tex_height = 1;

  int count = 0;
  int pp_count = 0;
  int px = 0;
  int py = 0;
  int id_count = 0;
  int leaf_count = 0;
  for (int i = 0; i < 2 * listN; i++) {
    if (Ns[i] == 0) continue;
    //       if ( i%2 != 0 ) continue;

    Quadratic** list = (Quadratic**)lists[i];
    for (int j = 0; j < Ns[i]; j++) {
      // 	  if ( !(list[j]->leaf) ) continue;
      if (list[j]->leaf) ++(leaf_count);

      list[j]->id_ = id_count;
      ++id_count;

      texcoordParam.push_back(
          Point2f((float)px / FID_DIV, (float)py / FID_DIV));

      // px, py
      pp[pp_count] = list[j]->centerX;
      ++(pp_count);
      pp[pp_count] = list[j]->centerY;
      ++(pp_count);
      pp[pp_count] = list[j]->centerZ;
      ++(pp_count);
      pp[pp_count] = list[j]->support;
      ++(pp_count);

      // px+1, py
      pp[pp_count] = list[j]->cXX;
      ++(pp_count);
      pp[pp_count] = list[j]->cYY;
      ++(pp_count);
      pp[pp_count] = list[j]->cZZ;
      ++(pp_count);
      pp[pp_count] = .0f;
      ++(pp_count);

      // px+2, py
      pp[pp_count] = list[j]->cXY;
      ++(pp_count);
      pp[pp_count] = list[j]->cYZ;
      ++(pp_count);
      pp[pp_count] = list[j]->cZX;
      ++(pp_count);
      pp[pp_count] = .0f;
      ++(pp_count);

      // px+3, py
      pp[pp_count] = list[j]->cX;
      ++(pp_count);
      pp[pp_count] = list[j]->cY;
      ++(pp_count);
      pp[pp_count] = list[j]->cZ;
      ++(pp_count);
      pp[pp_count] = list[j]->c0;
      ++(pp_count);

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
      if (count == max_width) {
        count = 0;
        ++tex_height;

        for (int k = 0; k < flo_width; ++k) pp.push_back(0.0f);

        px = 0;
        ++py;
      }
    }
  }

  cout << "Slim points ... all: " << id_count << " leaf " << leaf_count << endl;
  cout << "Param texture: width " << tex_width << " height " << tex_height
       << " pp size " << pp.size() << endl;

  // GL_FLOAT_RGBA_NV
  glTexImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, GL_RGBA32F_ARB, tex_width,
               tex_height, 0, GL_RGBA, GL_FLOAT, &pp[0]);
  cgGLSetTextureParameter(cgGetNamedParameter(frag1Program, "polyParam"),
                          polyParam);
  cgGLSetTextureParameter(cgGetNamedParameter(frag2Program, "polyParam"),
                          polyParam);
}
#endif  // OLD_SLIM

////////////////////////////////////////////////////////////////////////////////////

#ifdef NEW_SLIM
void setSlimTreePolyTexture() {
  GLuint polyParam;
  glGenTextures(1, &polyParam);
  glBindTexture(GL_TEXTURE_RECTANGLE_EXT, polyParam);
  glTexParameteri(GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

  GLsizei sizeW;
  glGetIntegerv(GL_MAX_RENDERBUFFER_SIZE_EXT, &sizeW);

  // Store coefficients of polynomial to a texture
  // Each polynomial uses 4 pixels to store coefficients and a center point
#if 0
  int max_width = 512;
  int tex_width  = max_width * 4; // pixel width
#endif
  int tex_width = sizeW;
  int max_width = (int)(tex_width / 4.0);
  cout << "max slimball size: " << max_width * sizeW << endl;

  int flo_width = tex_width * 4;
  std::vector<float> pp(flo_width);
  //   for ( int k = 0; k < flo_width; ++k ) pp.push_back( 0.0f );

  int tex_height = 1;

  int count = 0;
  int pp_count = 0;
  int px = 0;
  int py = 0;

  std::vector<SlimBallf*>& slimballs = slimtree.slimballs();
  for (int i = 0; i < slimballs.size(); ++i) {
    SlimBallf* sb = slimballs[i];

    texcoordParam.push_back(Point2f((float)px / FID_DIV, (float)py / FID_DIV));

    // px, py
    Point3f& center = sb->center();
    pp[pp_count++] = center.x;
    pp[pp_count++] = center.y;
    pp[pp_count++] = center.z;
    pp[pp_count++] = sb->support();

    // px+1, py
    Quadraticf& bf = (Quadraticf&)sb->bf();
    //       cout << bf << endl;
    std::vector<float>& c = bf.coeffs();
    pp[pp_count++] = c[0];  // XX
    pp[pp_count++] = c[1];  // YY
    pp[pp_count++] = c[2];  // ZZ
    pp[pp_count++] = .0f;

    // px+2, py
    pp[pp_count++] = c[3];
    pp[pp_count++] = c[4];
    pp[pp_count++] = c[5];
    pp[pp_count++] = .0f;

    // px+3, py
    pp[pp_count++] = c[6];
    pp[pp_count++] = c[7];
    pp[pp_count++] = c[8];
    pp[pp_count++] = c[9];

    px += 4;
    count++;
    if (count == max_width) {
      count = 0;
      ++tex_height;

      for (int k = 0; k < flo_width; ++k) pp.push_back(0.0f);

      px = 0;
      ++py;
    }
  }

  cout << "Param texture: width " << tex_width << " height " << tex_height
       << " pp size " << pp.size() << endl;

  // GL_FLOAT_RGBA_NV
  glTexImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, GL_RGBA32F_ARB, tex_width,
               tex_height, 0, GL_RGBA, GL_FLOAT, &pp[0]);
  cgGLSetTextureParameter(cgGetNamedParameter(frag1Program, "polyParam"),
                          polyParam);
  cgGLSetTextureParameter(cgGetNamedParameter(frag2Program, "polyParam"),
                          polyParam);
}
#endif  // NEW_SLIM

////////////////////////////////////////////////////////////////////////////////////

void setGaussKernelTexture() {
  GLuint gaussKernel;
  glGenTextures(1, &gaussKernel);
  glBindTexture(GL_TEXTURE_RECTANGLE_EXT, gaussKernel);
  glTexParameteri(GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

  std::vector<float> gk;
  int pix_width = GAUSS_WIDTH;
  int flo_width = pix_width * 4;
  int pix_height = 1;

  for (int k = 0; k < flo_width; ++k) gk.push_back(0.0f);

  int pix_count = 0;
  for (int i = 0; i < pix_width; ++i) {
    float r = (float)i / (float)(pix_width - 1);
#if 0
      gk[ pix_count ] = (float) std::exp( -r*r );
#endif

#if 1
    if (r < .5f)
      gk[pix_count] = -r * r + .75f;
    else if ((r >= .5f) && (r <= 1.0f))
      gk[pix_count] = 0.5f * (1.5f - r) * (1.5f - r);
#endif
    ++(pix_count);

    gk[pix_count] = .0f;
    ++(pix_count);
    gk[pix_count] = .0f;
    ++(pix_count);
    gk[pix_count] = .0f;
    ++(pix_count);
  }

  // GL_FLOAT_RGBA_NV
  glTexImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, GL_RGBA32F_ARB, pix_width,
               pix_height, 0, GL_RGBA, GL_FLOAT, &gk[0]);
  cgGLSetTextureParameter(cgGetNamedParameter(frag2Program, "gaussKernel"),
                          gaussKernel);
}

////////////////////////////////////////////////////////////////////////////////////

void drawQuad(int w, int h) {
  glBegin(GL_QUADS);
  glTexCoord2f(0, 0);
  glVertex2f(-1, -1);
  glTexCoord2f(0, h);
  glVertex2f(-1, 1);
  glTexCoord2f(w, h);
  glVertex2f(1, 1);
  glTexCoord2f(w, 0);
  glVertex2f(1, -1);
  glEnd();
}

float zEye(float m[16], float p[16], float x, float y, float z) {
  float a = m[0] * x + m[4] * y + m[8] * z + m[12];
  float b = m[1] * x + m[5] * y + m[9] * z + m[13];
  float c = m[2] * x + m[6] * y + m[10] * z + m[14];
  float d = m[3] * x + m[7] * y + m[11] * z + m[15];

  return (p[2] * a + p[6] * b + p[10] * c + p[14] * d);
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
void drawBillboardsLOD(float size_param) {
  float mat[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, mat);

  Point3f vRight(mat[0], mat[4], mat[8]);
  Point3f vUp(mat[1], mat[5], mat[9]);

  float prj[16];
  glGetFloatv(GL_PROJECTION_MATRIX, prj);

  vf.update(mat, prj);

  imagemaker.visitID++;
  imagemaker.stack[0] = imagemaker.root;
  imagemaker.top = 0;
  imagemaker.root->visit = imagemaker.visitID;

  slimball_count = 0;

  glBegin(GL_QUADS);

  while (imagemaker.top >= 0) {
    BasisFunction* bf = imagemaker.stack[imagemaker.top--];

    // view frustum culling
    bool f =
        vf.sphereInFrustum(bf->centerX, bf->centerY, bf->centerZ, bf->support);
    if (f != true) continue;

    // radius in screen coordinate
    float z = zEye(mat, prj, bf->centerX, bf->centerY, bf->centerZ);
    float rad_screen = bf->support * size_param / z;

    if (bf->leaf | (rad_screen < SLIMVIEW_FEW_PIX)) {
      //  	  if(!(imagemaker.both) && (bf->normalZAtC() > 0) )
      //  	    continue;

      Point3f vCenter(bf->centerX, bf->centerY, bf->centerZ);
      float pSize = bf->support;

      Point3f vPoint0(vCenter + ((-1.0f * vRight - vUp) * pSize));
      Point3f vPoint1(vCenter + ((vRight - vUp) * pSize));
      Point3f vPoint2(vCenter + ((vRight + vUp) * pSize));
      Point3f vPoint3(vCenter + ((-1.0f * vRight + vUp) * pSize));

      // 2D coordinate of a parameter texture
      // glColor4f( texcoordParam[k].x, texcoordParam[k].y, 1.0f, pSize /
      // FID_DIV );
      glColor4f(texcoordParam[bf->id_].x, texcoordParam[bf->id_].y, 1.0f,
                pSize / FID_DIV);

      glVertex3f(vPoint0.x, vPoint0.y, vPoint0.z);
      glVertex3f(vPoint1.x, vPoint1.y, vPoint1.z);
      glVertex3f(vPoint2.x, vPoint2.y, vPoint2.z);
      glVertex3f(vPoint3.x, vPoint3.y, vPoint3.z);

      ++(slimball_count);
    } else {
      int N = bf->childN;
      BasisFunction** bfs = bf->child;
      for (int i = 0; i < N; i++) {
        BasisFunction* bfi = bfs[i];
        if (bfi->visit != imagemaker.visitID) {
          bfi->visit = imagemaker.visitID;
          imagemaker.stack[++(imagemaker.top)] = bfs[i];
        }
      }
    }
  }

  glEnd();
}

void drawBillboards() {
  float mat[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, mat);

  Point3f vRight(mat[0], mat[4], mat[8]);
  Point3f vUp(mat[1], mat[5], mat[9]);

  //
  // draw quads
  //

  slimball_count = 0;

  glBegin(GL_QUADS);

  for (int i = 0; i < 2 * listN; i++) {
    if (Ns[i] == 0) continue;
    if (i % 2 != 0) continue;
    BasisFunction** list = lists[i];
    for (int j = 0; j < Ns[i]; j++) {
      if (!(list[j]->leaf)) continue;

      Point3f vCenter(list[j]->centerX, list[j]->centerY, list[j]->centerZ);
      float pSize = list[j]->support;

      Point3f vPoint0(vCenter + ((-1.0f * vRight - vUp) * pSize));
      Point3f vPoint1(vCenter + ((vRight - vUp) * pSize));
      Point3f vPoint2(vCenter + ((vRight + vUp) * pSize));
      Point3f vPoint3(vCenter + ((-1.0f * vRight + vUp) * pSize));

      // 2D coordinate of a parameter texture, (none), support
      glColor4f(texcoordParam[list[j]->id_].x, texcoordParam[list[j]->id_].y,
                1.0f, pSize / FID_DIV);

      glVertex3f(vPoint0.x, vPoint0.y, vPoint0.z);
      glVertex3f(vPoint1.x, vPoint1.y, vPoint1.z);
      glVertex3f(vPoint2.x, vPoint2.y, vPoint2.z);
      glVertex3f(vPoint3.x, vPoint3.y, vPoint3.z);

      ++slimball_count;
    }
  }
  glEnd();
}

void drawPointSpritesLOD(float size_param) {
  float mat[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, mat);

  float prj[16];
  glGetFloatv(GL_PROJECTION_MATRIX, prj);

  vf.update(mat, prj);

  glEnable(GL_POINT_SPRITE_ARB);
  glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);

  slimball_count = 0;

  glBegin(GL_POINTS);

  imagemaker.visitID++;
  imagemaker.stack[0] = imagemaker.root;
  imagemaker.top = 0;
  imagemaker.root->visit = imagemaker.visitID;
  while (imagemaker.top >= 0) {
    BasisFunction* bf = imagemaker.stack[imagemaker.top--];

    // view frustum culling
    bool f =
        vf.sphereInFrustum(bf->centerX, bf->centerY, bf->centerZ, bf->support);
    if (f != true) continue;

    // radius in screen coordinate
    float z = zEye(mat, prj, bf->centerX, bf->centerY, bf->centerZ);
    float rad_screen = bf->support * size_param / z;

    if (bf->leaf | (rad_screen < SLIMVIEW_FEW_PIX)) {
      //  	  if(!(imagemaker.both) && (bf->normalZAtC() > 0) )
      //  	    continue;

      // 2D coordinate of a parameter texture
      glColor4f(texcoordParam[bf->id_].x, texcoordParam[bf->id_].y, 1.0f,
                bf->support / FID_DIV);

      glVertex3f(bf->centerX, bf->centerY, bf->centerZ);

      ++(slimball_count);

    } else {
      int N = bf->childN;
      BasisFunction** bfs = bf->child;
      for (int i = 0; i < N; i++) {
        BasisFunction* bfi = bfs[i];
        if (bfi->visit != imagemaker.visitID) {
          bfi->visit = imagemaker.visitID;
          imagemaker.stack[++(imagemaker.top)] = bfs[i];
        }
      }
    }
  }

  glEnd();

  glDisable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
  glDisable(GL_POINT_SPRITE_ARB);
}

void drawPointSprites() {
  glEnable(GL_POINT_SPRITE_ARB);
  glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);

  slimball_count = 0;

  glBegin(GL_POINTS);
  for (int i = 0; i < 2 * listN; i++) {
    if (Ns[i] == 0) continue;
    if (i % 2 != 0) continue;
    BasisFunction** list = lists[i];
    for (int j = 0; j < Ns[i]; j++) {
      if (!(list[j]->leaf)) continue;

      // 2D coordinate of a parameter texture, ID, support
      glColor4f(texcoordParam[list[j]->id_].x, texcoordParam[list[j]->id_].y,
                1.0f, list[j]->support / FID_DIV);
      //
      // (center_x, center_y, center_z, radius of support)
      //
      glVertex3f(list[j]->centerX, list[j]->centerY, list[j]->centerZ);

      ++(slimball_count);
    }
  }
  glEnd();

  glDisable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
  glDisable(GL_POINT_SPRITE_ARB);
}

void displayBall() {
  if (isCapture) cap_buffer->Activate();

  pane.setBackgroundColor(bgColor[0], bgColor[1], bgColor[2]);
  pane.setIsGradientBackground(false);
  pane.clear(pbuffer_width, pbuffer_height);
  pane.setView();
  pane.setLight();

  ::glShadeModel(GL_SMOOTH);
  ::glEnable(GL_LIGHTING);

  ballMtl.bind();

  slimball_count = 0;

  for (int i = 0; i < 2 * listN; i++) {
    if (Ns[i] == 0) continue;
    if (i % 2 != 0) continue;

    BasisFunction** list = lists[i];
    for (int j = 0; j < Ns[i]; j++) {
      if (!(list[j]->leaf)) continue;

      glPushMatrix();
      glTranslatef(list[j]->centerX, list[j]->centerY, list[j]->centerZ);

      //DrawSphere( list[j]->support, 10, 10 );
      drawSphereClientWithRadius(sm, list[j]->support());
      
      glPopMatrix();

      ++(slimball_count);
    }
  }

  pane.finish();
}

#endif  // OLD_SLIM

////////////////////////////////////////////////////////////////////////////////////
//
// new slim
//
////////////////////////////////////////////////////////////////////////////////////

#ifdef NEW_SLIM

void drawSlimTreeBillboards() {
  float mat[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, mat);

  Point3f vRight(mat[0], mat[4], mat[8]);
  Point3f vUp(mat[1], mat[5], mat[9]);

  //
  // draw quads
  //

  slimball_count = 0;

  glBegin(GL_QUADS);

  std::vector<SlimBallf*>& slimballs = slimtree.slimballs();
  for (int i = 0; i < slimballs.size(); ++i) {
    SlimBallf* slimball = slimballs[i];

    if (slimball->childs().size()) continue;

    Point3f vCenter(slimball->center());
    float pSize = slimball->support();

    Point3f vPoint0(vCenter + ((-1.0f * vRight - vUp) * pSize));
    Point3f vPoint1(vCenter + ((vRight - vUp) * pSize));
    Point3f vPoint2(vCenter + ((vRight + vUp) * pSize));
    Point3f vPoint3(vCenter + ((-1.0f * vRight + vUp) * pSize));

    // 2D coordinate of a parameter texture, (none), support
    int id = slimball->id();
    glColor4f(texcoordParam[id].x, texcoordParam[id].y, 1.0f, pSize / FID_DIV);

    glVertex3f(vPoint0.x, vPoint0.y, vPoint0.z);
    glVertex3f(vPoint1.x, vPoint1.y, vPoint1.z);
    glVertex3f(vPoint2.x, vPoint2.y, vPoint2.z);
    glVertex3f(vPoint3.x, vPoint3.y, vPoint3.z);

    ++(slimball_count);
  }

  glEnd();
}

int drawSlimBallBillboardsLOD(SlimBallf* slimball, int level, Point3f& vRight,
                              Point3f& vUp) {
  if (!slimball) return 0;

  int count = 0;
  if (level == slimball_level) {
    glBegin(GL_QUADS);

    Point3f vCenter(slimball->center());
    float pSize = slimball->support();

    Point3f vPoint0(vCenter + ((-1.0f * vRight - vUp) * pSize));
    Point3f vPoint1(vCenter + ((vRight - vUp) * pSize));
    Point3f vPoint2(vCenter + ((vRight + vUp) * pSize));
    Point3f vPoint3(vCenter + ((-1.0f * vRight + vUp) * pSize));

    // 2D coordinate of a parameter texture, (none), support
    int id = slimball->id();
    glColor4f(texcoordParam[id].x, texcoordParam[id].y, 1.0f, pSize / FID_DIV);

    glVertex3f(vPoint0.x, vPoint0.y, vPoint0.z);
    glVertex3f(vPoint1.x, vPoint1.y, vPoint1.z);
    glVertex3f(vPoint2.x, vPoint2.y, vPoint2.z);
    glVertex3f(vPoint3.x, vPoint3.y, vPoint3.z);

    glEnd();

    ++(count);

    return count;
  }

  if (slimball->childs().size()) {
    for (int i = 0; i < slimball->childs().size(); ++i)
      count +=
          drawSlimBallBillboardsLOD(slimball->child(i), level + 1, vRight, vUp);
  }

  return count;
}

void drawSlimTreeBillboardsLOD() {
  float mat[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, mat);

  Point3f vRight(mat[0], mat[4], mat[8]);
  Point3f vUp(mat[1], mat[5], mat[9]);

  //
  // draw quads
  //

  SlimBallf* slimball = slimtree.root();
  slimball_count = drawSlimBallBillboardsLOD(slimball, 0, vRight, vUp);
}

int drawSlimBallBillboardsLODError(SlimBallf* slimball, int level,
                                   Point3f& vRight, Point3f& vUp) {
  if (!slimball) return 0;

  int count = 0;
  if (slimball->userDefined() < slimball_error) {
    glBegin(GL_QUADS);

    Point3f vCenter(slimball->center());
    float pSize = slimball->support();

    Point3f vPoint0(vCenter + ((-1.0f * vRight - vUp) * pSize));
    Point3f vPoint1(vCenter + ((vRight - vUp) * pSize));
    Point3f vPoint2(vCenter + ((vRight + vUp) * pSize));
    Point3f vPoint3(vCenter + ((-1.0f * vRight + vUp) * pSize));

    // 2D coordinate of a parameter texture, (none), support
    int id = slimball->id();
    glColor4f(texcoordParam[id].x, texcoordParam[id].y, 1.0f, pSize / FID_DIV);

    glVertex3f(vPoint0.x, vPoint0.y, vPoint0.z);
    glVertex3f(vPoint1.x, vPoint1.y, vPoint1.z);
    glVertex3f(vPoint2.x, vPoint2.y, vPoint2.z);
    glVertex3f(vPoint3.x, vPoint3.y, vPoint3.z);

    glEnd();

    ++(count);

    return count;
  }

  if (slimball->childs().size()) {
    for (int i = 0; i < slimball->childs().size(); ++i)
      count += drawSlimBallBillboardsLODError(slimball->child(i), level + 1,
                                              vRight, vUp);
  }

  return count;
}

void drawSlimTreeBillboardsLODError() {
  float mat[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, mat);

  Point3f vRight(mat[0], mat[4], mat[8]);
  Point3f vUp(mat[1], mat[5], mat[9]);

  //
  // draw quads
  //

  SlimBallf* slimball = slimtree.root();
  slimball_count = drawSlimBallBillboardsLODError(slimball, 0, vRight, vUp);
}

void drawSlimTreePointSprites() {
  glEnable(GL_POINT_SPRITE_ARB);
  glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);

  slimball_count = 0;

  glBegin(GL_POINTS);

  std::vector<SlimBallf*>& slimballs = slimtree.slimballs();
  for (int i = 0; i < slimballs.size(); ++i) {
    SlimBallf* slimball = slimballs[i];

    if (slimball->childs().size()) continue;

    // 2D coordinate of a parameter texture, ID, support
    float pSize = slimball->support();
    int id = slimball->id();
    glColor4f(texcoordParam[id].x, texcoordParam[id].y, 1.0f, pSize / FID_DIV);
    //
    // (center_x, center_y, center_z, radius of support)
    //
    Point3f vCenter(slimball->center());
    glVertex3f(vCenter.x, vCenter.y, vCenter.z);

    ++(slimball_count);
  }

  glEnd();

  glDisable(GL_VERTEX_PROGRAM_POINT_SIZE_ARB);
  glDisable(GL_POINT_SPRITE_ARB);
}

void drawSlimTreeBall() {
  if (isCapture) cap_buffer->Activate();

  pane.setBackgroundColor(bgColor[0], bgColor[1], bgColor[2]);
  pane.setIsGradientBackground(false);
  pane.clear(pbuffer_width, pbuffer_height);
  pane.setView();
  pane.setIsLightOn(1, false);
  pane.setLight();

  ::glShadeModel(GL_SMOOTH);
  ::glEnable(GL_LIGHTING);

  //   ballMtl.bind();

  slimball_count = 0;

  std::vector<SlimBallf*>& slimballs = slimtree.slimballs();
  for (int i = 0; i < slimballs.size(); ++i) {
    if (slimballs[i]->childs().size()) continue;

    glPushMatrix();
    Point3f& center = slimballs[i]->center();
    glTranslatef(center.x, center.y, center.z);

    float color[3];
    convColor(slimballs[i]->userDefined(), color);
    ballMtl.setDiffuseColor(color[0], color[1], color[2], 1.0f);
    ballMtl.bind();

    //DrawSphere( slimballs[i]->support(), 10, 10 );
    drawSphereClientWithRadius(sm, slimballs[i]->support());

    glPopMatrix();

    ++(slimball_count);
  }

  pane.finish();
}

int drawSlimBallBallLODError(SlimBallf* slimball, int level) {
  if (!slimball) return 0;

  int count = 0;
  if (slimball->userDefined() < slimball_error) {
    glPushMatrix();
    Point3f& center = slimball->center();
    glTranslatef(center.x, center.y, center.z);

    float color[3];
    convColor(slimball->userDefined(), color);
    ballMtl.setDiffuseColor(color[0], color[1], color[2], 1.0f);
    ballMtl.bind();

    //DrawSphere( slimball->support(), 10, 10 );
    drawSphereClientWithRadius(sm, slimball->support());

    glPopMatrix();

    ++(count);

    return count;
  }

  if (slimball->childs().size()) {
    for (int i = 0; i < slimball->childs().size(); ++i)
      count += drawSlimBallBallLODError(slimball->child(i), level + 1);
  }

  return count;
}

void drawSlimTreeBallLODError() {
  if (isCapture) cap_buffer->Activate();

  pane.setBackgroundColor(bgColor[0], bgColor[1], bgColor[2]);
  pane.setIsGradientBackground(false);
  pane.clear(pbuffer_width, pbuffer_height);
  pane.setView();
  pane.setLight();

  ::glShadeModel(GL_SMOOTH);
  ::glEnable(GL_LIGHTING);

  //   ballMtl.bind();

  slimball_count = 0;
  SlimBallf* slimball = slimtree.root();
  slimball_count += drawSlimBallBallLODError(slimball, 0);

  pane.finish();
}

#endif  // NEW_SLIM

////////////////////////////////////////////////////////////////////////////////////

void displayPoly() {
  if (isCapture) cap_buffer->Activate();

  pane.setBackgroundColor(bgColor[0], bgColor[1], bgColor[2]);
  pane.setIsGradientBackground(false);
  pane.clear(pbuffer_width, pbuffer_height);
  pane.setView();
  pane.setLight();

  if (!(glmeshR.empty())) glmeshR.draw();

  pane.finish();
}

////////////////////////////////////////////////////////////////////////////////////

// #define SLIM_DEBUG 1
// #define SLIM_DEBUG_II 1

#ifdef OLD_SLIM
void drawSlimTree(float* size_param) {
  if (draw_mode == DRAW_BILLBOARD) {
    if (lod_mode == LOD_LEVEL)
      drawBillboardsLOD(size_param[0] * size_param[1]);
    else
      drawBillboards();
  } else {
    if (lod_mode == LOD_LEVEL)
      drawPointSpritesLOD(size_param[0] * size_param[1]);
    else
      drawPointSprites();
  }
}

#endif  // OLD_SLIM

#ifdef NEW_SLIM
void drawSlimTree(float* size_param) {
  if (draw_mode == DRAW_BILLBOARD) {
    if (lod_mode == LOD_LEVEL) {
      drawSlimTreeBillboardsLOD();
    } else if (lod_mode == LOD_ERROR) {
      drawSlimTreeBillboardsLODError();
    } else
      drawSlimTreeBillboards();
  } else {
    drawSlimTreePointSprites();
  }
}

void displaySlimTreeBall() {
  if (lod_mode == LOD_ERROR)
    drawSlimTreeBallLODError();
  else
    drawSlimTreeBall();
}

#endif  // NEW_SLIM

void displaySLIMGPU() {
  ////////////////////////////////////////////////////////////////////////////////////
  //
  // 1st pass
  //
  ////////////////////////////////////////////////////////////////////////////////////

  pbuffer->Activate();

  pane.setBackgroundColor(.0f, .0f, .0f);
  pane.setIsGradientBackground(false);
  pane.clear(pbuffer_width, pbuffer_height);
  pane.setView();

#ifdef OGL_QUERY
  glBeginQueryARB(GL_SAMPLES_PASSED_ARB, oq_slim1);
#endif

  //
  // vertex program
  //
  cgGLEnableProfile(vertProfile);
  cgGLBindProgram(vert1Program);

  cgGLSetStateMatrixParameter(
      cgGetNamedParameter(vert1Program, "ModelViewProj"),
      CG_GL_MODELVIEW_PROJECTION_MATRIX, CG_GL_MATRIX_IDENTITY);

  GLint viewport[4];
  GLdouble projMatrix[16];
  ::glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  ::glGetIntegerv(GL_VIEWPORT, viewport);

  //
  // see eq.5 in [Botsch04]
  //
  // projMatrix[5]: 2.0 * n / (t - b) ...
  // viewport[3]:   h_vp
  float size_param[2];
  size_param[0] = (float)projMatrix[5];
  size_param[1] = (float)viewport[3];
  cgGLSetParameter2fv(cgGetNamedParameter(vert1Program, "size_param"),
                      size_param);

  //
  // fragment program
  //
  cgGLEnableProfile(fragProfile);
  cgGLBindProgram(frag1Program);

  cgGLEnableTextureParameter(cgGetNamedParameter(frag1Program, "polyParam"));

  cgGLSetStateMatrixParameter(cgGetNamedParameter(frag1Program, "ModelView"),
                              CG_GL_MODELVIEW_MATRIX, CG_GL_MATRIX_IDENTITY);

  cgGLSetStateMatrixParameter(cgGetNamedParameter(frag1Program, "ModelViewInv"),
                              CG_GL_MODELVIEW_MATRIX, CG_GL_MATRIX_INVERSE);

  cgGLSetStateMatrixParameter(
      cgGetNamedParameter(frag1Program, "ModelViewProj"),
      CG_GL_MODELVIEW_PROJECTION_MATRIX, CG_GL_MATRIX_IDENTITY);

  cgGLSetStateMatrixParameter(
      cgGetNamedParameter(frag1Program, "ModelViewProjInv"),
      CG_GL_MODELVIEW_PROJECTION_MATRIX, CG_GL_MATRIX_INVERSE);

  //
  // get width and height of near plane
  //
  float left1, right1, bottom1, top1, near1, far1;
  pane.getProjectionParametersFromPerspective(&left1, &right1, &bottom1, &top1,
                                              &near1, &far1);
  float viewParam[4];
  viewParam[0] = near1;
  viewParam[1] = far1;
  cgGLSetParameter4fv(cgGetNamedParameter(frag1Program, "viewParam"),
                      viewParam);

  // set draw mode ... 0: Billboard, 1; Point Sprite
  cgGLSetParameter1f(cgGetNamedParameter(frag1Program, "drawMode"), draw_mode);

  drawSlimTree(size_param);

  cgGLDisableTextureParameter(cgGetNamedParameter(frag1Program, "polyParam"));
  cgGLDisableProfile(fragProfile);

  cgGLDisableProfile(vertProfile);

  pane.finish();

#ifdef OGL_QUERY
  glEndQueryARB(GL_SAMPLES_PASSED_ARB);
#endif

#ifndef RENDER_FBO
  // copy fpbuffer to a texture
  glBindTexture(GL_TEXTURE_RECTANGLE_EXT, qhatT);
  glCopyTexSubImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, 0, 0, 0, 0, pbuffer_width,
                      pbuffer_height);
  cgGLSetTextureParameter(cgGetNamedParameter(frag2Program, "qhatT"), qhatT);
#endif

  pbuffer->Deactivate();

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // 2nd pass
  //
  ////////////////////////////////////////////////////////////////////////////////////

  //   glPointParameterfARB( GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT );

  pbuffer16->Activate();

  for (int i = 0; i < 2; ++i) {
#ifndef RENDER_FBO

    glDrawBuffer(buffers[i]);

#else  // RENDER_FBO

    pbuffer16->Draw(i);

#endif  // RENDER_FBO

    pane.setBackgroundColor(.0f, .0f, .0f);
    pane.setIsGradientBackground(false);
    pane.clear(pbuffer_width, pbuffer_height);
    pane.setView();
  }

#ifdef OGL_QUERY
  glBeginQueryARB(GL_SAMPLES_PASSED_ARB, oq_slim2);
#endif

  //
  // vertex program
  //
  cgGLEnableProfile(vertProfile);
  cgGLBindProgram(vert2Program);

  cgGLSetStateMatrixParameter(
      cgGetNamedParameter(vert2Program, "ModelViewProj"),
      CG_GL_MODELVIEW_PROJECTION_MATRIX, CG_GL_MATRIX_IDENTITY);

  // assign the same parameter as a 1st pass
  cgGLSetParameter2fv(cgGetNamedParameter(vert2Program, "size_param"),
                      size_param);

  //
  // fragment program
  //
  cgGLEnableProfile(fragProfile);
  cgGLBindProgram(frag2Program);

#ifdef RENDER_FBO
  pbuffer->Bind(0);
  cgGLSetTextureParameter(cgGetNamedParameter(frag2Program, "qhatT"),
                          pbuffer->tex(0));
#endif

  cgGLEnableTextureParameter(cgGetNamedParameter(frag2Program, "qhatT"));

  cgGLEnableTextureParameter(cgGetNamedParameter(frag2Program, "polyParam"));
  cgGLEnableTextureParameter(cgGetNamedParameter(frag2Program, "gaussKernel"));

  cgGLSetStateMatrixParameter(cgGetNamedParameter(frag2Program, "ModelViewInv"),
                              CG_GL_MODELVIEW_MATRIX, CG_GL_MATRIX_INVERSE);

  cgGLSetStateMatrixParameter(
      cgGetNamedParameter(frag2Program, "ModelViewProjInv"),
      CG_GL_MODELVIEW_PROJECTION_MATRIX, CG_GL_MATRIX_INVERSE);

  //   GLint viewport[4];
  //   glGetIntegerv( GL_VIEWPORT, viewport );

  //   float viewParam[4];
  viewParam[0] = (float)viewport[0];
  viewParam[1] = (float)viewport[1];
  viewParam[2] = (float)viewport[2];
  viewParam[3] = (float)viewport[3];
  cgGLSetParameter4fv(cgGetNamedParameter(frag2Program, "viewParam"),
                      viewParam);

  // set draw mode ... 0: Billboard, 1; Point Sprite
  cgGLSetParameter1f(cgGetNamedParameter(frag2Program, "drawMode"), draw_mode);

  // Use Multiple Draw Buffers
  // set destination buffers
  glDrawBuffersATI(sizeof(buffers) / sizeof(GLenum), buffers);

#define DRAW_BLEND 1

  //
  // draw quads
  //

#ifdef DRAW_BLEND
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  glDepthMask(GL_FALSE);
#endif

  drawSlimTree(size_param);

#ifdef DRAW_BLEND
  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
#endif

  cgGLDisableTextureParameter(cgGetNamedParameter(frag2Program, "qhatT"));
  pbuffer->Release(0);

  cgGLDisableTextureParameter(cgGetNamedParameter(frag2Program, "gaussKernel"));
  cgGLDisableTextureParameter(cgGetNamedParameter(frag2Program, "polyParam"));

  cgGLDisableProfile(fragProfile);

  cgGLDisableProfile(vertProfile);

  pane.finish();

#ifdef OGL_QUERY
  glEndQueryARB(GL_SAMPLES_PASSED_ARB);
#endif

  //
  // copy fpbuffer (buffer #0) to a position texture
  //

#ifndef RENDER_FBO

  glReadBuffer(buffers[0]);

  glBindTexture(GL_TEXTURE_RECTANGLE_EXT, positionT);
  glCopyTexSubImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, 0, 0, 0, 0, pbuffer_width,
                      pbuffer_height);
  cgGLSetTextureParameter(cgGetNamedParameter(frag3Program, "positionT"),
                          positionT);

#else  // RENDER_FBO

#endif  // RENDER_FBO

  //
  // copy fpbuffer (buffer #1) to a normal texture
  //

#ifndef RENDER_FBO

  glReadBuffer(buffers[1]);

  glBindTexture(GL_TEXTURE_RECTANGLE_EXT, normalT);
  glCopyTexSubImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, 0, 0, 0, 0, pbuffer_width,
                      pbuffer_height);
  cgGLSetTextureParameter(cgGetNamedParameter(frag3Program, "normalT"),
                          normalT);

#else  // RENDER_FBO

#endif  // RENDER_FBO

  pbuffer16->Deactivate();

  //   glPointParameterfARB( GL_POINT_SPRITE_COORD_ORIGIN, GL_UPPER_LEFT );

  ////////////////////////////////////////////////////////////////////////////////////
  //
  // 3rd pass
  //
  ////////////////////////////////////////////////////////////////////////////////////

  if (isCapture) cap_buffer->Activate();

  pane.clear(pbuffer_width, pbuffer_height);

#if !defined(RENDER_FBO)

#else  // RENDER_FBO

  pane.setView();

#endif  // RENDER_FBO

  cgGLEnableProfile(fragProfile);
  cgGLBindProgram(frag3Program);

  // matrices
  cgGLSetStateMatrixParameter(
      cgGetNamedParameter(frag3Program, "ModelViewProj"),
      CG_GL_MODELVIEW_PROJECTION_MATRIX, CG_GL_MATRIX_IDENTITY);

  cgGLSetStateMatrixParameter(cgGetNamedParameter(frag3Program, "ModelView"),
                              CG_GL_MODELVIEW_MATRIX, CG_GL_MATRIX_IDENTITY);

  cgGLSetStateMatrixParameter(cgGetNamedParameter(frag3Program, "ModelViewIT"),
                              CG_GL_MODELVIEW_MATRIX,
                              CG_GL_MATRIX_INVERSE_TRANSPOSE);

#if !defined(RENDER_FBO)

#else  // RENDER_FBO

  pane.initView();

#endif  // RENDER_FBO

  // background color
  cgGLSetParameter4fv(cgGetNamedParameter(frag3Program, "bgColor"), bgColor);

  // materials
  cgGLSetParameter4f(g_CGparam_AmbMat, myMatl[0], myMatl[1], myMatl[2],
                     myMatl[3]);
  cgGLSetParameter4f(g_CGparam_DiffMat, myMatl[4], myMatl[5], myMatl[6],
                     myMatl[7]);
  cgGLSetParameter4f(g_CGparam_EmisMat, myMatl[8], myMatl[9], myMatl[10],
                     myMatl[11]);
  cgGLSetParameter4f(g_CGparam_SpecMat, myMatl[12], myMatl[13], myMatl[14],
                     myMatl[15]);
  cgGLSetParameter4f(g_CGparam_ShineMat, myMatl[16], 0.0f, 0.0f, 0.0f);

  // light
  float current_light[3];
  pane.getRealLightPosition(current_light);
  cgGLSetParameter4f(g_CGparam_LightPos, current_light[0], current_light[1],
                     current_light[2], 1.0f);
  cgGLSetParameter4f(g_CGparam_LightVec, -current_light[0], -current_light[1],
                     -current_light[2], 1.0f);
  cgGLSetParameter4f(g_CGparam_AmbLight, myLight[4], myLight[5], myLight[6],
                     myLight[7]);
  cgGLSetParameter4f(g_CGparam_DiffLight, myLight[8], myLight[9], myLight[10],
                     myLight[11]);
  cgGLSetParameter4f(g_CGparam_SpecLight, myLight[12], myLight[13], myLight[14],
                     myLight[15]);
  cgGLSetParameter3f(g_CGparam_AttenLight, myLight[16], myLight[17],
                     myLight[18]);  // constant, linear, quadratic
  cgGLSetParameter2f(g_CGparam_SpotLight, cosf(180.0f),
                     0.0f);  // cosf(spot_cuttoff_angle), spot power

#ifdef RENDER_FBO
  pbuffer16->Bind(0);
  cgGLSetTextureParameter(cgGetNamedParameter(frag3Program, "positionT"),
                          pbuffer16->tex(0));
#endif

  cgGLEnableTextureParameter(cgGetNamedParameter(frag3Program, "positionT"));

#ifdef RENDER_FBO
  pbuffer16->Bind(1);
  cgGLSetTextureParameter(cgGetNamedParameter(frag3Program, "normalT"),
                          pbuffer16->tex(1));
#endif

  cgGLEnableTextureParameter(cgGetNamedParameter(frag3Program, "normalT"));

  {
    drawQuad(pbuffer_width, pbuffer_height);
  }

  cgGLDisableTextureParameter(cgGetNamedParameter(frag3Program, "normalT"));
  cgGLDisableTextureParameter(cgGetNamedParameter(frag3Program, "positionT"));
  pbuffer16->Release(0);
  pbuffer16->Release(1);

  cgGLDisableProfile(fragProfile);

  ////////////////////////////////////////////////////////////////////////////////////


  pane.finish();

#ifdef OGL_QUERY
  GLuint sl1, sl2;
  glGetQueryObjectuivARB(oq_slim1, GL_QUERY_RESULT_ARB, &sl1);
  glGetQueryObjectuivARB(oq_slim2, GL_QUERY_RESULT_ARB, &sl2);

  cout << "occlusion query ... 1st pass: " << sl1 << " 2nd pass: " << sl2
       << endl;
#endif
}

void display() {

#if 0
  pane.setBackgroundColor( bgColor[0], bgColor[1], bgColor[2] );
  pane.clear( pbuffer_width, pbuffer_height );
  pane.initView();
#endif

  switch (render_mode) {
    case RENDER_GPU:
      displaySLIMGPU();
      break;

    case RENDER_BALL:
#ifdef OLD_SLIM
      displayBall();
#endif  // OLD_SLIM
#ifdef NEW_SLIM
      displaySlimTreeBall();
#endif  // NEW_SLIM
      break;

    case RENDER_NONE:
      displayPoly();
      break;

    default:
      break;
  }

  checkGLErrors((char*)"display");
}

////////////////////////////////////////////////////////////////////////////////////

static void error_callback(int error, const char* description) {
  fputs(description, stderr);
}

static void framebuffer_size_callback(GLFWwindow* window, int width,
                                      int height) {
  glViewport(0, 0, width, height);
}

// keyboard event processing function
static void key_callback(GLFWwindow* window, int key, int scancode, int action,
                         int mods) {
  VWIO vw_out;

  // ESC
  if ((key == GLFW_KEY_ESCAPE) && (action == GLFW_PRESS)) {
    cout << "max fps: " << max_c11fps << endl;
    destroyRenderTextures();
    destroyFPBuffer();
    destroyCg();
    glfwSetWindowShouldClose(window, GL_TRUE);
    return;
  }

  // q
  else if ((key == GLFW_KEY_Q) && (action == GLFW_PRESS)) {
    cout << "max fps: " << max_c11fps << endl;
    destroyRenderTextures();
    destroyFPBuffer();
    destroyCg();
    glfwSetWindowShouldClose(window, GL_TRUE);
    return;
  }

  // v (output to .vw file)
  else if ((key == GLFW_KEY_V) && (action == GLFW_PRESS)) {
    cout << "output to .vw file ... " << endl;
    vw_out.outputToFile("tmp.vw", pane.manip());
    cout << "done." << endl;
    return;
  }

  // p (screen capture)
  else if ((key == GLFW_KEY_P) && (action == GLFW_PRESS)) {
    //
    isCapture = true;
    return;
  }

  // m
  else if ((key == GLFW_KEY_M) && (action == GLFW_PRESS)) {
    if (lod_mode == LOD_ERROR) {
      slimball_error *= 0.1;
      if (slimball_error < 1.0e-15) slimball_error = 1.0e-15;
      cout << "slimball error threshold: " << slimball_error << endl;
    } else if (lod_mode == LOD_LEVEL) {
      ++(slimball_level);
      cout << "slimball level: " << slimball_level << endl;
    }
    return;
  }

  // 1 (Billboard, LOD None)
  else if ((key == GLFW_KEY_1) && (action == GLFW_PRESS)) {
    render_mode = RENDER_GPU;
    draw_mode = DRAW_BILLBOARD;
    lod_mode = LOD_NONE;
    max_c11fps = 0.0;
    return;
  }

  // 2 (Point Sprite, LOD None)
  else if ((key == GLFW_KEY_2) && (action == GLFW_PRESS)) {
    render_mode = RENDER_GPU;
    draw_mode = DRAW_POINTSPRITE;
    lod_mode = LOD_NONE;
    max_c11fps = 0.0;
    return;
  }

  // 3 (Billboard, LOD Level)
  else if ((key == GLFW_KEY_3) && (action == GLFW_PRESS)) {
    render_mode = RENDER_GPU;
    draw_mode = DRAW_BILLBOARD;
    lod_mode = LOD_LEVEL;
    max_c11fps = 0.0;
    return;
  }

  // 4 (Point Sprite, LOD Level)
  else if ((key == GLFW_KEY_4) && (action == GLFW_PRESS)) {
    render_mode = RENDER_GPU;
    draw_mode = DRAW_POINTSPRITE;
    lod_mode = LOD_LEVEL;
    max_c11fps = 0.0;
    return;
  }

  // 5 (Billboard, LOD Error)
  else if ((key == GLFW_KEY_5) && (action == GLFW_PRESS)) {
    render_mode = RENDER_GPU;
    draw_mode = DRAW_BILLBOARD;
    lod_mode = LOD_ERROR;
    max_c11fps = 0.0;
    return;
  }

  // 6 (Ball, LOD None)
  else if ((key == GLFW_KEY_6) && (action == GLFW_PRESS)) {
    render_mode = RENDER_BALL;
    lod_mode = LOD_NONE;
    max_c11fps = 0.0;
    return;
  }

  // 7 (Ball, LOD Error)
  else if ((key == GLFW_KEY_7) && (action == GLFW_PRESS)) {
    render_mode = RENDER_BALL;
    lod_mode = LOD_ERROR;
    max_c11fps = 0.0;
    return;
  }

  // 8 (None, LOD None)
  else if ((key == GLFW_KEY_8) && (action == GLFW_PRESS)) {
    render_mode = RENDER_NONE;
    lod_mode = LOD_NONE;
    max_c11fps = 0.0;
    return;
  }

  // shift
  else if ((key == GLFW_KEY_LEFT_SHIFT) && (action == GLFW_PRESS)) {
    shift_key_pressed = true;
    return;
  } else if ((key == GLFW_KEY_LEFT_SHIFT) && (action == GLFW_RELEASE)) {
    shift_key_pressed = false;
    return;
  }

  // control
  else if ((key == GLFW_KEY_LEFT_CONTROL) && (action == GLFW_PRESS)) {
    control_key_pressed = true;
    return;
  } else if ((key == GLFW_KEY_LEFT_CONTROL) && (action == GLFW_RELEASE)) {
    control_key_pressed = false;
    return;
  }
}

// mouse event processing function
static void mousebutton_callback(GLFWwindow* window, int button, int action,
                                 int mods) {
  double xd, yd;
  glfwGetCursorPos(window, &xd, &yd);
  pane.setScreenXY((int)xd, (int)yd);

  if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_PRESS)) {
    left_button_pressed = true;
    pane.startRotate();
    pane.startZoom();
    pane.startMove();
  } else if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_RELEASE)) {
    left_button_pressed = false;
    pane.finishRMZ();
  } else if ((button == GLFW_MOUSE_BUTTON_2) && (action == GLFW_PRESS)) {
    right_button_pressed = true;
  } else if ((button == GLFW_MOUSE_BUTTON_2) && (action == GLFW_RELEASE)) {
    right_button_pressed = false;
  }
}

// cursor event processing function
static void cursorpos_callback(GLFWwindow* window, double xd, double yd) {
  int x = (int)xd;
  int y = (int)yd;

  if (left_button_pressed && !shift_key_pressed && !control_key_pressed) {
    pane.updateRotate(x, y);
  } else if (left_button_pressed && shift_key_pressed && !control_key_pressed) {
    pane.updateZoom(x, y);
  } else if (left_button_pressed && !shift_key_pressed && control_key_pressed) {
    pane.updateMove(x, y);
  }
}

// mouse wheel
static void scroll_callback(GLFWwindow* window, double xoffset,
                            double yoffset) {
  pane.updateWheelZoom(yoffset);
}

// window resize
static void windowsize_callback(GLFWwindow* window, int w, int h) {
  width = w;
  height = h;
  pane.changeSize(w, h);
  pane.setW(w);
  pane.setH(h);

  destroyRenderTextures();
  destroyFPBuffer();

  width = w;
  height = h;
  pbuffer_width = w;
  pbuffer_height = h;
  cap_width = w;
  cap_height = h;

  initFPBuffer(pbuffer_width, pbuffer_height, cap_width, cap_height);
  initRenderTexture(pbuffer_width, pbuffer_height);
}

////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " in.slim2t" << std::endl;
    return EXIT_FAILURE;
  }

  char filename[BUFSIZ];
  strcpy(filename, argv[1]);

#ifdef OLD_SLIM
  // initialize old SLIM
  listN = filemanager.readBallFile(filename, lists, Ns, deg, oriented);
  imagemaker.setBfLists(deg, listN, lists, Ns);
  imagemaker.setBoth(!oriented);
#endif  // OLD_SLIM

#ifdef NEW_SLIM
  SlimTreefIO slimtreeIO(slimtree);
  slimtreeIO.inputFromFile(filename);
#endif  // NEW_SLIM

  // GLFW
  // Window initialization starts here
  glfwSetErrorCallback(error_callback);
  if (!glfwInit()) return EXIT_FAILURE;

  //  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  //  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  //  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  //  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  //  glfwWindowHint(GLFW_SAMPLES, 4);  // MSAA 4x sampling

  GLFWwindow* window =
      glfwCreateWindow(width, height, "GLFW Window", NULL, NULL);
  if (!window) {
    std::cerr << "Failed to create GLFW window" << std::endl;
    const char* error;
    glfwGetError(&error);
    if (error) {
      std::cerr << "GLFW Error: " << error << std::endl;
    }
    glfwTerminate();
    return EXIT_FAILURE;
  }

  glfwMakeContextCurrent(window);

  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
  glfwSetKeyCallback(window, key_callback);
  glfwSetMouseButtonCallback(window, mousebutton_callback);
  glfwSetCursorPosCallback(window, cursorpos_callback);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetWindowSizeCallback(window, windowsize_callback);

  glfwShowWindow(window);
  glfwFocusWindow(window);

  pane.init(width, height);
  // pane.setIsLightOn( 1, false );

#ifdef OLD_SLIM
  Point3f p(.0f, .0f, 55.0f);
  pane.setViewPoint(p);
  Vector3f v(.0f, .0f, -55.0f);
  pane.setViewVector(v);
  pane.setMagObject(10.0f);
#endif  // OLD_SLIM

#ifdef NEW_SLIM
  Point3f p(.0f, .0f, 30.0f);
  pane.setViewPoint(p);
  Vector3f v(.0f, .0f, -30.0f);
  pane.setViewVector(v);
  pane.setMagObject(10.0f);
#endif  // NEW_SLIM

  // Initialize float pbuffer, Cg, Texture
  pane.initGLEW();
  initFPBuffer(pbuffer_width, pbuffer_height, cap_width, cap_height);
  initRenderTexture(pbuffer_width, pbuffer_height);
  initCg();
  setGaussKernelTexture();

#ifdef OLD_SLIM
  setPolyTexture();
#endif  // OLD_SLIM

#ifdef NEW_SLIM
  setSlimTreePolyTexture();
#endif  // NEW_SLIM

  // Initialize some state for the GLUT window's rendering context.
  glfwMakeContextCurrent(window);

  pane.initGL();
  pane.setIsLightOn(1, false);
  pane.setLightParameters(0, myLight);

  // GLSL shaders
  // shaders = initShaders("shaders");

  // Initial viewport setup (required to avoid small window size on Mac)
  int fbWidth, fbHeight;
  glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
  glViewport(0, 0, fbWidth, fbHeight);

  glfwSwapInterval(0);
  // Window initialization ends here

  // SphereMesh initialization
  sm = makeUnitSphereClient(/*n1=*/20, /*n2=*/20);

  while (!glfwWindowShouldClose(window)) {
    // clear and initialization
    pane.setBackgroundColor(bgColor[0], bgColor[1], bgColor[2]);
    glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
    pane.clear(fbWidth, fbHeight);
    pane.initView();
    // pane.update(glmeshl.material());

    display();

    pane.finish();

    // fps measurement
    double f = c11fps.CheckGetFPS();
    if (max_c11fps < f) max_c11fps = f;

    std::stringstream ss;
    ss << std::fixed << std::setprecision(3) << std::setw(8) << f
       << " fps - max " << std::setw(8) << max_c11fps << " fps";
    std::string buf = ss.str();

    string str0, str1;
    if ( draw_mode == DRAW_BILLBOARD )
      str0 = "Billboard";
    else if ( draw_mode == DRAW_POINTSPRITE )
      str0 = "Point Sprite";
    if ( lod_mode == LOD_NONE )
      str1 = "Leaf";
    else
      str1 = "LOD";
    if (render_mode == RENDER_BALL) {
      str0 = "Ball";
    }
    
    std::string txt = "GLFW Window - " + buf + " " + str0 + " " + str1;
    glfwSetWindowTitle(window, txt.c_str());

    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  glfwDestroyWindow(window);
  glfwTerminate();

  return EXIT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////////
