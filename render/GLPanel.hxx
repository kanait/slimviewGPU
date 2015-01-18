////////////////////////////////////////////////////////////////////
//
// $Id: GLPanel.hxx $
//
// Copyright (c) 2002-2010 by Takashi Kanai. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _GLPANEL_HXX
#define _GLPANEL_HXX 1

#include "envDep.h"
#include "mydef.h"

#include <string>
#include <vector>
using namespace std;

#define GLEW_STATIC 1
#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glut.h>

#include "Arcball.hxx"
#include "GLLight.hxx"
#include "GLMaterial.hxx"

class GLPanel {

private:
  
  int width_;
  int height_;

  // background color
  float bgrgb_[3];

  // projection parameters
  double fov_;
  double aspect_;
  double nearP_;
  double farP_;

  //
  // Light
  //

  static const int num_lights_ = 4;
  std::vector<GLLight> light_;
  
  //SelList sel_list_;

  // Transformation Flag
  bool rotateFlag_;
  bool moveFlag_; 
  bool zoomFlag_;

  // Display Flag
  bool drawWireframe_;
  bool drawShading_;

  // Gradient Background flag
  bool gradientBackground_;

  // Arcball
  Arcball manip_;

  // Transformation in 2D
  float move2d_x_;
  float move2d_y_;
  float scale2d_;
  int x0_;
  int y0_;
  int s0_;

  // View Pont, View Vector
  Point3f view_point_;
  Vector3f view_vector_;

public:
  
  GLPanel() {};
  virtual ~GLPanel() {};

  void init( int w, int h ) {
    bgrgb_[0] = 1.0f;
    bgrgb_[1] = 1.0f;
    bgrgb_[2] = 1.0f;

    setW( w );
    setH( h );

    initViewParameters( w, h );

    rotateFlag_ = false;
    moveFlag_ = false;
    zoomFlag_ = false;
    gradientBackground_ = true;
  
    // 2D transformation
    resetView2d();
  };

  void initViewParameters( int w, int h ) {
    // set projection parameters
    fov_ = 30.0f;
    aspect_ = (double) w / (double) h;
    //   nearP_ = 1.0f;
    //   farP_ = 10.0f;
    nearP_ = .01f;
    farP_ = 100000.0f;

    // view point, view vector
    view_point_.set( 0.0f, 0.0f, 3.0f );
    view_vector_.set( 0.0f, 0.0f, -3.0f );

    // for Arcball
    manip_.init();
    manip_.setHalfWHL( (int) (w/2.0), (int) (h/2.0) );
  };

  void initGLEW( bool flag = false ) {
    GLenum err = glewInit();
    if( err != GLEW_OK ) cout << "Error: %s" << glewGetErrorString(err) << endl;
    
    if ( flag ) wglSwapIntervalEXT(0);
  };

  //
  // 2D Functions
  //

  void initGL2d() {
    ::glDisable( GL_ALPHA_TEST );
    ::glDisable( GL_BLEND );
    ::glDisable( GL_DEPTH_TEST );
    ::glDisable( GL_LIGHTING );
    ::glDisable( GL_TEXTURE_1D );
    ::glDisable( GL_TEXTURE_2D );
    ::glDisable( GL_POLYGON_OFFSET_FILL );
    //   ::glShadeModel( GL_FLAT );

    setIsGradientBackground( false );
  };

  void clear2d() {
    ::glViewport( 0, 0, w(), h() );
    ::glClearColor( bgrgb_[0], bgrgb_[1], bgrgb_[2], 0.0f );
    ::glClear( GL_COLOR_BUFFER_BIT );
    ::glDisable( GL_DEPTH_TEST );
  };

  void setView2d() {
    ::glMatrixMode( GL_PROJECTION );
    ::glLoadIdentity();
    ::gluOrtho2D( .0f, w(), .0f, h() );
    ::glMatrixMode( GL_MODELVIEW );
    ::glLoadIdentity();

    ::glTranslatef( move2d_x_, move2d_y_, .0f );
    ::glScalef( scale2d_, scale2d_, 1.0f );
  };

  void finish2d() {
    ::glFinish();
  };

  //
  // 3D Functions
  //

  void initGL() { initGL( false, false ); };
  void initGL( bool isTransparency, bool isLineSmooth ) {
    // initialize lights
    // number of lights initialize
    light_.resize( num_lights_ );

    initLight();

    light_[0].setIsOn( true );
    light_[1].setIsOn( false );
    light_[2].setIsOn( false );
    light_[3].setIsOn( false );
  
    ::glPolygonMode( GL_BACK, GL_FILL );
    //   ::glCullFace( GL_BACK );
    //   ::glEnable( GL_CULL_FACE );

    // enable depth testing
    ::glEnable( GL_DEPTH_TEST );
    //   ::glDepthFunc( GL_LEQUAL );
    ::glDepthFunc( GL_LESS );
    //   ::glClearDepth( 1.0f );
  
    ::glEnable( GL_NORMALIZE );

    // transparency settings
    if ( isTransparency )
      {
	::glEnable( GL_ALPHA_TEST );
	::glEnable( GL_BLEND );
	::glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
      }

    ::glEnable( GL_POLYGON_OFFSET_FILL );
    ::glEnable( GL_POLYGON_OFFSET_LINE );
    ::glPolygonOffset( (float) 1.0, (float) 1e-5 );

    // line anti-aliasing
    if ( isLineSmooth )
      {
	::glEnable( GL_LINE_SMOOTH );
	::glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
      }

    initView();
  };

  void initLight() {
    //    ::glLightModelfv(GL_LIGHT_MODEL_AMBIENT, GL_FALSE);
    //    ::glLightModeli( GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE );
    //    ::glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );
    for ( int i = 0; i < num_lights_; ++i )
      {
	light_[i].setID( i );
	light_[i].init();
	light_[i].setPos( &(light_position_[4*i]) );
	light_[i].bind();
      }
    //    ::glEnable( GL_LIGHTING );
  };

  void initView() {
    ::glMatrixMode( GL_PROJECTION );
    ::glLoadIdentity();
    ::glMatrixMode( GL_MODELVIEW );
    ::glLoadIdentity();
  };

  void changeSize( int w, int h ) {
    setW( w );
    setH( h );

    aspect_ = (double) w / (double) h;

    manip_.setHalfWHL( (int) ((float) w / 2.0), (int) ((float) h / 2.0) );
  };

  //
  // draw functions
  //
  void clear() { clear( w(), h() ); };
  void clear( int w, int h ) {
    ::glViewport( 0, 0, w, h );
    ::glClearColor( bgrgb_[0], bgrgb_[1], bgrgb_[2], 0.0f );
    ::glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    // draw gradient background
    if ( isGradientBackground() )
      drawGradientBackground();

    ::glEnable( GL_DEPTH_TEST );
  };

  void drawGradientBackground() {
    ::glDisable(GL_DEPTH_TEST);
    ::glPushMatrix();
    ::glLoadIdentity();
    ::glShadeModel( GL_SMOOTH );
    ::glMatrixMode(GL_PROJECTION);
    ::glPushMatrix();
    ::glLoadIdentity();

    ::glBegin(GL_QUADS);
    ::glColor3f(0.0F, 0.0F, 0.1F); ::glVertex2f(-1.0F, -1.0F);
    ::glColor3f(0.0F, 0.0F, 0.1F); ::glVertex2f( 1.0F, -1.0F);
    ::glColor3f(0.4F, 0.4F, 1.0F); ::glVertex2f( 1.0F,  1.0F);
    ::glColor3f(0.4F, 0.4F, 1.0F); ::glVertex2f(-1.0F,  1.0F);
    //   ::glColor3f(0.2F, 0.0F, 1.0F); ::glVertex2f(-1.0F, -1.0F);
    //   ::glColor3f(0.2F, 0.0F, 1.0F); ::glVertex2f( 1.0F, -1.0F);
    //   ::glColor3f(0.0F, 0.0F, 0.1F); ::glVertex2f( 1.0F,  1.0F);
    //   ::glColor3f(0.0F, 0.0F, 0.1F); ::glVertex2f(-1.0F,  1.0F);
    ::glEnd();
    ::glEnable(GL_DEPTH_TEST);
    ::glShadeModel( GL_FLAT );

    ::glPopMatrix();
    ::glMatrixMode(GL_MODELVIEW);
    ::glPopMatrix();
  };

  void drawAxis() {
    ::glDisable( GL_LIGHTING );

    ::glColor3f( .0f, .0f, .0f );

    ::glBegin( GL_LINES );
    ::glVertex3f(  .0f, .0f, .0f );
    ::glVertex3f( 1.0f, .0f, .0f );
    ::glEnd();

    ::glBegin( GL_LINES );
    ::glVertex3f(  .0f, .0f, .0f );
    ::glVertex3f( .0f, 1.0f, .0f );
    ::glEnd();

    ::glBegin( GL_LINES );
    ::glVertex3f(  .0f, .0f, .0f );
    ::glVertex3f( .0f, .0f, 1.0f );
    ::glEnd();

#if 0
    glQuickText glqt;
    glqt.printfAt( 1.1f, .0f, .0f, .05f, "x");
    glqt.printfAt( .0f, 1.1f, .0f, .05f, "y");
    glqt.printfAt( .0f, .0f, 1.1f, .05f, "z");
#endif
  };

  void setProjectionView() {
    // set Perspective View
    ::glMatrixMode( GL_PROJECTION );
    ::glLoadIdentity();

    setPerspective();
  };

  void setPerspective() {
    aspect_ = (double) w() / (double) h();
    ::gluPerspective( fov_, aspect_, nearP_, farP_ );
  }

  void setModelView() {
    ::glMatrixMode( GL_MODELVIEW );
    ::glLoadIdentity();

    Point3f look = view_point_ + view_vector_;
    ::gluLookAt( view_point_.x, view_point_.y, view_point_.z,
		 look.x, look.y, look.z, 0.0, 1.0, 0.0 );

    // from Arcball
    ::glTranslatef( 0, 0, manip_.seezo() );
    ::glMultMatrixf( (GLfloat *)(&(manip_.mNow().m00)) );
    Point3f o = manip_.offset();
    ::glTranslatef( -o.x, -o.y, -o.z );
  };
  
  void setView() {
    setProjectionView();
    setModelView();
  }

  bool getProjectionParameters( float* leftOut, float* rightOut,
				float* botOut,  float* topOut,
				float* nearOut, float* farOut ) {
    float m[16];
    ::glGetFloatv( GL_PROJECTION_MATRIX, m );
    bool isPerspective;

    if (m[15] == 0.0) 
      {
	// perspective
	//       float p[16];
	const float x = m[0];  // 2N / (R-L)
	const float y = m[5];  // 2N / (T-B)
	const float a = m[8];  // (R+L) / (R-L)
	const float b = m[9];  // (T+B) / (T-B)
	const float c = m[10]; // -(F+N) / (F-N)
	const float d = m[14]; // -2FN / (F-N)
	//
	// These equations found with simple algebra, knowing the arithmetic
	// use to set up a typical perspective projection matrix in OpenGL.
	//
	const float nearZ = -d / (1.0 - c);
	const float farZ = (c - 1.0) * nearZ / (c + 1.0);
	const float left = nearZ * (a - 1.0) / x;
	const float right = 2.0 * nearZ / x + left;
	const float bottom = nearZ * (b - 1.0) / y;
	const float top = 2.0 * nearZ / y + bottom;

	isPerspective = true;
	*leftOut = left;
	*rightOut = right;
	*botOut = bottom;
	*topOut = top;
	*nearOut = nearZ;
	*farOut = farZ;
      }
    else 
      {
	// orthographic
	const float x = m[0];  //  2 / (R-L)
	const float y = m[5];  //  2 / (T-B)
	const float z = m[10]; // -2 / (F-N)
	const float a = m[12]; // -(R+L) / (R-L)
	const float b = m[13]; // -(T+B) / (T-B)
	const float c = m[14]; // -(F+N) / (F-N)
	//
	// again, simple algebra
	//
	const float right  = -(a - 1.0) / x;
	const float left   = right - 2.0 / x;
	const float top    = -(b - 1.0) / y;
	const float bottom = top - 2.0 / y;
	const float farZ   = (c - 1.0) / z;
	const float nearZ  = farZ + 2.0 / z;

	isPerspective = false;
	*leftOut = left;
	*rightOut = right;
	*botOut = bottom;
	*topOut = top;
	*nearOut = nearZ;
	*farOut = farZ;
      }
  
    return isPerspective;
  };

  void getProjectionParametersFromPerspective( float* leftOut, float* rightOut,
					       float* botOut,  float* topOut,
					       float* nearOut, float* farOut ) {
    *topOut   = nearP_ * std::tan(fov_ * M_PI / 360.0);
    *botOut   = -(*topOut);
    *leftOut  = (*botOut) * aspect_;
    *rightOut = (*topOut) * aspect_;
    *nearOut  = nearP_;
    *farOut   = farP_;
  };

  void setLight() {
    for ( int i = 0; i < num_lights_; ++i )
      {
	if ( isLightOn(i) ) light_[i].on(); 
	else light_[i].off();
      }
  };

  void finish() {
    ::glFlush();

    ::glFinish();
  }

  void draw() {
    clear();
    setView();
    setLight();

    ::glEnable( GL_LIGHTING );
    GLMaterial glm;
    glm.bind();

    Point3f vec( 0.0, 0.0, 0.0 );
    double radius = .5;

    GLUquadricObj   *qobj;
    GLint slices = 50,staks = 50;
	
    if ((qobj = gluNewQuadric()) != NULL) {
      ::glPushMatrix();
      ::glShadeModel( GL_SMOOTH );
      ::glTranslatef( vec.x, vec.y, vec.z );
      ::gluSphere(qobj, radius, slices, staks);
      ::glPopMatrix();
      ::gluDeleteQuadric(qobj);
    }

    ::glDisable( GL_LIGHTING );

    finish();
  };

  Point3f& view_point() { return view_point_; };
  void setViewPoint( Point3f& p ) { view_point_.set(p); };
  void setViewPoint( float x, float y, float z ) { view_point_.set(Point3f(x, y, z)); };
  Vector3f& view_vector() { return view_vector_; };
  void setViewVector( Vector3f& p ) { view_vector_.set(p); };
  void setViewVector( float x, float y, float z ) { view_vector_.set(Vector3f(x, y, z)); };
  
  void setViewParameters( double width, double height, double fov, double nearP, double farP ) {
    fov_ = fov;
    nearP_ = nearP;
    farP_ = farP;
    aspect_ = (double) width / (double) height;
  };
  void setNearFarPlanes( double nearP, double farP ) {
    nearP_ = nearP;
    farP_ = farP;
  };
  void setFOV( double fov ) {
    fov_ = fov;
  };

  void setMagObject( float f ) { manip_.setMagObject( f ); };
    
  // get "real" view point
  void getRealViewPoint( Point3f& p ) {
    Point3f look( view_point_ + view_vector_ );
    Vector3f sub( view_point_ - look );
    Vector3f sub_tran;
    Matrix4f inv( manip_.mNow() );
    inv.transform( sub, &sub_tran );
    p.set( sub_tran + look );
  };
    
  void getRealViewPoint( float* pos ) {
    Point3f p;
    getRealViewPoint( p );
    pos[0] = p.x; pos[1] = p.y; pos[2] = p.z; 
  };

  // lights

  const int num_lights() const { return num_lights_; };

  void setLightPos( int i, float lpos[] ) {
    setLightPos( i, lpos[0], lpos[1], lpos[2], lpos[3] );
  };

  void initLightPos( int i ) {
    light_[i].setPos( light_position_[4*i],
		      light_position_[4*i+1],
		      light_position_[4*i+2],
		      light_position_[4*i+3] );
  };

  void setLightPos( int i, Point3f& p ) {
    light_[i].setPos( p.x, p.y, p.z, 0.0f );
  };
    
  void setLightPos( int i, float l0, float l1, float l2, float l3 ) {
//     light_position_[4*i]   = l0;
//     light_position_[4*i+1] = l1;
//     light_position_[4*i+2] = l2;
//     light_position_[4*i+3] = l3;
    light_[i].setPos( l0, l1, l2, l3 );
  };

  void getLightPos( int i, Point3f& p ) { 
    float w;
    light_[i].getPos( &(p.x), &(p.y), &(p.z), &w );
  };

  void getInitLightPos( int i, Point3f& p ) {
    p.x = light_position_[4*i];
    p.y = light_position_[4*i+1];
    p.z = light_position_[4*i+2];
  };

  void getLightPos( int i, float lpos[] ) { 
    lpos[0] = light_position_[4*i];
    lpos[1] = light_position_[4*i+1];
    lpos[2] = light_position_[4*i+2];
    lpos[3] = light_position_[4*i+3];
  };
  void getLightVec( int i, float lvec[] ) { 
    lvec[0] = -light_position_[4*i];
    lvec[1] = -light_position_[4*i+1];
    lvec[2] = -light_position_[4*i+2];
    lvec[3] = 1.0f;
  };

  void getRealLightPosition( Point3f& p ) {
    getRealLightPosition( 0, p );
  };
    
  // get "real" light position
  void getRealLightPosition( int i, Point3f& p ) {
    Point3f light_pos( light_position_[4*i], light_position_[4*i+1], light_position_[4*i+2] );
    Point3f light_vec( -light_position_[4*i], -light_position_[4*i+1], -light_position_[4*i+2] );
    Point3f look( light_pos + light_vec );
    Vector3f sub( light_pos - look );
    Vector3f sub_tran;
    Matrix4f inv( manip_.mNow() );
    inv.transform( sub, &sub_tran );
    p.set( sub_tran + look );
  };

  void getRealLightPosition( int i, float* pos ) {
    Point3f p;
    getRealLightPosition( i, p );
    pos[0] = p.x; pos[1] = p.y; pos[2] = p.z; 
  };
    
  void getRealLightPosition( float* pos ) {
    getRealLightPosition( 0, pos );
  };

  void setLightParameters( int i, float* light ) {
    light_[i].set( light );
  };
    
#if 0  
  // wireframe color functions
  float* wfColor() { return &(wfrgb_[0]); };
  void setWireframeColor( unsigned char r, unsigned char g,
			   unsigned char b ) {
    wfrgb_[0] = (float) r / 255.0;
    wfrgb_[1] = (float) g / 255.0;
    wfrgb_[2] = (float) b / 255.0;
  };
#endif
  
  // background color functions
  float* bgColor() { return &(bgrgb_[0]); };
  void setBackgroundColor( unsigned char r, unsigned char g,
			   unsigned char b ) {
    bgrgb_[0] = (float) r / 255.0;
    bgrgb_[1] = (float) g / 255.0;
    bgrgb_[2] = (float) b / 255.0;
  };
  void setBackgroundColor( float r, float g, float b ) {
    bgrgb_[0] = r;
    bgrgb_[1] = g;
    bgrgb_[2] = b;
  };
  
  void setSize( int w, int h ) {
    setW( w );
    setH( h );
  };
  void setW( int w ) { width_ = w; };
  void setH( int h ) { height_ = h; };
  int w() const { return width_; };
  int h() const { return height_; };

  void reshape( int w, int h ) {
    setW( w );
    setH( h );
    manip_.setHalfWHL( (int) (w/2.0), (int) (h/2.0) );
  };

  // rotate, translation and zoom function
  
  void setScreenXY( int x, int y ) {
    manip_.setScrnXY( x, y );
    x0_ = x;
    y0_ = y;
    s0_ = y;
  };

  void startRotate() {
    setIsRotate( true );
    manip_.vFrom( manip_.mouse_on_sphere( manip_.scrn_x(),
					    manip_.scrn_y(),
					    manip_.halfW(),
					    manip_.halfH() ) );  
  };

  void startMove() { setIsMove( true ); };
  void startZoom() { setIsZoom( true ); };

  void updateRotate( int x, int y ) {
    manip_.vTo( manip_.mouse_on_sphere(x, y, manip_.halfW(),
					 manip_.halfH()) );
    manip_.updateRotate( x, y );
    manip_.setScrnXY( x, y );  
  };

  void updateMove( int x, int y ) {
    manip_.updateMove( x, y, manip_.scrn_x(), manip_.scrn_y() );
    manip_.setScrnXY( x, y );
  };
  
  void updateZoom( int x, int y ) {
    manip_.updateZoom( x, y, manip_.scrn_x(), manip_.scrn_y() );
    manip_.setScrnXY( x, y );
  };

  
  void updateMove2d( int x, int y ) {
    move2d_x_ += (float) (x - x0_); x0_ = x;
    move2d_y_ -= (float) (y - y0_); y0_ = y;
  };
  
  void updateZoom2d( int x, int y ) {
    scale2d_ -= .1f * (float) (y - s0_); s0_ = y;
    if ( scale2d_ < .01f ) scale2d_ = .01f;
  };

  
  void finishRMZ() {
    setIsRotate( false );
    setIsZoom( false );
    setIsMove( false );
    manip_.qDown( manip_.qNow() );
    manip_.mDown( manip_.mNow() );
    x0_ = 0;
    y0_ = 0;
    s0_ = 0;
  };

  void resetView2d() {
    move2d_x_ = .0f;
    move2d_y_ = .0f;
    scale2d_  = 1.0f;
    x0_ = 0;
    y0_ = 0;
    s0_ = 0;
  };
  
  // 
  bool isRotate() const { return rotateFlag_; };
  bool isMove() const { return moveFlag_; };
  bool isZoom() const { return zoomFlag_; };
  void setIsRotate( bool f ) { rotateFlag_ = f; };
  void setIsMove( bool f ) { moveFlag_ = f; };
  void setIsZoom( bool f ) { zoomFlag_ = f; };

  //
  // flag functions
  //

  // light functions
  void setIsLightOn( unsigned short id , bool f ) { light_[id].setIsOn( f ); };
  bool isLightOn( unsigned short id ) const { return light_[id].isOn(); };

  // gradient background functions
  void setIsGradientBackground( bool f ) { gradientBackground_ = f; };
  bool isGradientBackground() const { return gradientBackground_; };
  
  // arcball
  Arcball& manip() { return manip_; };

};

#endif // _GLPANEL_HXX
