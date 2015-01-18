#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

int width = 512;
int height = 512;

void drawQuad( int w, int h )
{
  glBegin( GL_QUADS );
  glTexCoord2f(0, 0); glVertex2f(-.5, -.5);
  glTexCoord2f(0, h); glVertex2f(-.5,  .5);
  glTexCoord2f(w, h); glVertex2f( .5,  .5);
  glTexCoord2f(w, 0); glVertex2f( .5, -.5);
  glEnd();
}

void display()
{
  glViewport( 0, 0, width, height );
  glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

/*   glEnable (GL_DEPTH); */
/*   glEnable (GL_CULL_FACE); */

  // 深度テストを無効にしてα合成有効
  glDisable (GL_DEPTH_TEST);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE);
/*   glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA); */
/*   glBlendFunc( GL_ONE, GL_ONE ); */
  
  glColor4f( 1.0, 1.0, 0.0, 0.25 );
  drawQuad( width, height );
  drawQuad( width, height );

  glFinish();

  glutSwapBuffers();
}

int main( int argc, char*argv[] )
{
  glutInitWindowSize( width, height );
  glutInit( &argc, argv );

  //glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_DEPTH | GLUT_STENCIL );
  glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_DEPTH );
  glutCreateWindow("test");

  glutDisplayFunc( display );

  glutMainLoop();
}
