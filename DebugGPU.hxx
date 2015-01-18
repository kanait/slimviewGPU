////////////////////////////////////////////////////////////////////
//
// $Id: $
//
// Copyright (c) 2005 by Takashi Kanai
// All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _DEBUGGPU_HXX
#define _DEBUGGPU_HXX 1

#include "envDep.h"

#include <vector>
using namespace std;

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>

#if 0
#define GLH_EXT_SINGLE_FILE
#include <glh/glh_extensions.h>
#endif

#include "RenderTexture.h"

#define nRGBA_ 4

// mode
#define DEBUG_POINTS   0
#define DEBUG_TEXTURE  1

class DebugGPU {

public:

  void setIsValid( bool f ) { isValid_ = f; };
  bool isValid() const { return isValid_; };
  
  DebugGPU() {
    rt_ = NULL;
    n_frame_ = -1;
    isCapture_ = false;
    mode_ = DEBUG_POINTS;
    isValid_ = true;
  };
  
  DebugGPU( int width, int height ) {
    rt_ = NULL;
    width_ = width;
    height_ = height;
    n_frame_ = -1;
    isCapture_ = false;
  };
  
  ~DebugGPU() {
    if ( rt_ ) delete rt_;
  };

  void init( int width, int height, int mode ) {
    width_ = width;
    height_ = height;
    mode_ = mode;
    if ( mode_ == DEBUG_POINTS )
      {
	rt_ = new RenderTexture( "float=32 rgba textureRECT", width_, height_, GL_TEXTURE_RECTANGLE_NV );
      }
    else
      {
	rt_ = new RenderTexture( "rgba textureRECT", width_, height_, GL_TEXTURE_RECTANGLE_NV );
      }
  };

  void begin() {
    
    if ( !isValid_ ) return;

    Activate();
    points_.clear();
    isCapture_ = false;
    ++( n_frame_ );
  };
  
  void end() {

    if ( !isValid_ ) return;

    Deactivate();
  };

  void Activate() { rt_->Activate(); };
  void Deactivate() { rt_->Deactivate(); };

  void capture( bool isWriteToFile ) {

    if ( !isValid_ ) return;

    if ( mode_ != DEBUG_POINTS ) return;
    
    isCapture_ = true;
    points_.resize( width_ * height_ * nRGBA_ );
    glReadPixels( 0, 0, width_, height_, GL_RGBA, GL_FLOAT, &points_[0] );
    
    if ( isWriteToFile )
      {
	char filename[BUFSIZ];
	sprintf( filename, "tmp%03d.pnt", n_frame_ );
	writeToFile( filename );
      }
  };

  void writeToFile( char* filename ) {

    std::ofstream ofs( filename );
    for ( int i = 0; i < width_ * height_; ++i )
      {
	if ( (points_[nRGBA_*i] == .0f) && (points_[nRGBA_*i+1] == .0f) && (points_[nRGBA_*i+2] == .0f) )
	  {
	  }
	else
	  {
	    ofs << points_[nRGBA_*i] << " " << points_[nRGBA_*i+1] << " "
		<< points_[nRGBA_*i+2] << " " << points_[nRGBA_*i+3] << endl;
	  }
      }
    ofs.close();
  };

  void renderPoints() {

    if ( !isValid_ ) return;

    if ( mode_ != DEBUG_POINTS ) return;
    if ( !isCapture_ ) return;
      
    glBegin( GL_POINTS );
    for ( int i = 0; i < width_ * height_; ++i )
      {
	glVertex3f( points_[nRGBA_*i], points_[nRGBA_*i+1], points_[nRGBA_*i+2] );
      }
    glEnd();
  };

  void drawQuad() {

    glBegin( GL_QUADS );
    glTexCoord2f(0, 0); glMultiTexCoord2fARB(GL_TEXTURE1_ARB, 0, 0); glVertex2f(-1, -1);
    glTexCoord2f(width_, 0); glMultiTexCoord2fARB(GL_TEXTURE1_ARB, width_, 0); glVertex2f(1, -1);
    glTexCoord2f(width_, height_); glMultiTexCoord2fARB(GL_TEXTURE1_ARB, width_, height_); glVertex2f(1, 1);
    glTexCoord2f(0, height_); glMultiTexCoord2fARB(GL_TEXTURE1_ARB, 0, height_); glVertex2f(-1, 1);
    glEnd();
  };

  void renderTexture() {

    if ( !isValid_ ) return;

    if ( mode_ != DEBUG_TEXTURE ) return;

    glDisable(GL_DEPTH_TEST);
    glDrawBuffer( GL_BACK );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity();

    glEnable( GL_TEXTURE_RECTANGLE_NV );
    glActiveTextureARB( GL_TEXTURE0_ARB );
    glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );    

    rt_->Bind();
    drawQuad();
    rt_->Release();
    
    glDisable( GL_TEXTURE_RECTANGLE_NV );

    glPopMatrix();
  }
    
private:

  RenderTexture* rt_;

  int mode_;

  int width_;
  int height_;

  int n_frame_;

  bool isValid_;

  bool isCapture_;
  std::vector<float> points_;

};

#endif // _DEBUGGPU_HXX



