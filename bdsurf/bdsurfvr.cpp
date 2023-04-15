/* Manifold-based Surface Construction
   Copyright (C) 2004 Lexing Ying, Denis Zorin, New York University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */
#include "bdsurfvr.hpp"
#include "trackball.hpp"
#include "IImage.hpp"

using std::ofstream;

// ---------------------------------------------------------------------- 
BdSurfVr::BdSurfVr(int* argc, char** argv, BdSurf* bdsurf, int lvl)
{
  //set + precopute
  _bdsurf = bdsurf;
  _cnt = 0;
  int TTL = pow2(lvl);
  int NS = TTL/4*3-1;	 //int NS = pow2(_lvl-1);
  //int numf = _bdsurf->numf();
  int numv = _bdsurf->numv();
  _Vfxy.resize(numv);	 _Vfpos.resize(numv);	 _Vfnor.resize(numv);  //Point3 bbmax[3];
  for(int V=0; V<numv; V++) {
	 int K = _bdsurf->valence(V);
	 if(K==0) continue;
	 _Vfxy[V].resize(K);		_Vfpos[V].resize(K);		_Vfnor[V].resize(K);
	 for(int f=0; f<K; f++) {
		NumMat<Point2>& cmxy  = _Vfxy[V][f];  cmxy.resize(NS+1,NS+1);
		NumMat<Point3>& cmpos = _Vfpos[V][f]; cmpos.resize(NS+1,NS+1);
		NumMat<Point3>& cmnor = _Vfnor[V][f]; cmnor.resize(NS+1,NS+1);
		double step = 1.0/double(TTL);
		for(int j=0; j<=NS; j++)
		  for(int i=0; i<=NS; i++) {
			 double cd[2];			 cd[0] = i*step;			 cd[1] = j*step;
			 double xy[2];			 iC( _bdsurf->Vfcd2xy(BdSurf::EVAL_VL, V, f, cd, xy) ); 
			 Point3 ret[6];
			 iC( _bdsurf->eval(  BdSurf::EVAL_VL|BdSurf::EVAL_FD, V, xy, ret) );
			 //-----------
			 cmxy(i,j) = Point2(xy[0], xy[1]);				  cmpos(i,j) = ret[0];				  cmnor(i,j) = cross(ret[1],ret[2]).dir();
		  }
	 }
  }
  
  //window size and other stuff
  _bdsurf->bbx(_bbmin, _bbmax);  //cerr<<_bbmin<<"  "<<_bbmax<<endl;
  _ws = Index2(640,480);
  _eyepos = (_bbmax[1]-_bbmin[1])/tan(M_PI/6.0)*1.5;  //_eyepos = 10.0;
  trackball(_quat, 0, 0, 0, 0);
  _proj = 1;
  //control
  for(int i=0; i<3; i++)		_mouse_down[i] = 0;
  _vrop = RENDER_SURF;//_vrop[0] = Fl3d::VR_D; //fluid option  _vrop[0] = 0; //fluid option  //_vrop[1] = Sh3d::VR_FRAME; //shell option
  for(int i=0; i<4; i++) _startquat[i] = _quat[i];
  _startdrag = Index2(0,0);
  //constants
  _self = this;  //_rendered = false;
  //---------------------------
  //create
  glutInit(argc, argv);
  glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(_ws[0], _ws[1]);
  glutCreateWindow("BdSurfVr");  glCheck();  //glutReshapeWindow(_ws[0], _ws[1]);
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClearDepth(1.0);
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glutSwapBuffers ();
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glutSwapBuffers ();
  //---------------------------
  //callbacks
  glutReshapeFunc(BdSurfVr::reshapeWrapper);
  glutDisplayFunc(BdSurfVr::displayWrapper);
  glutMouseFunc(  BdSurfVr::mouseWrapper);
  glutMotionFunc( BdSurfVr::motionWrapper);
  glutKeyboardFunc(BdSurfVr::keyboardWrapper);
  glutIdleFunc(   BdSurfVr::idleWrapper);  //char* p = (char*)glGetString(GL_EXTENSIONS); cerr<<p<<endl;
  glutPostRedisplay();
  //---------------------------
  //lighting
  GLfloat lightPos0[4] = { 1, 1, 2, 0 };
  GLfloat lightPos1[4] = { 0, -1, 0, 0 };
  GLfloat lightCol0[4] = { 1, 1, 1, 1 };
  GLfloat lightCol1[4] = { 0.1, 0.1, 0.1, 1 };
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightCol0);
  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);
  glLightfv(GL_LIGHT1, GL_DIFFUSE,  lightCol1);
  glEnable(GL_LIGHT1);
  
  glShadeModel(GL_SMOOTH);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  
  glEnable(GL_POINT_SMOOTH);  //glEnable(GL_LINE_SMOOTH);  //glEnable(GL_POLYGON_SMOOTH);
  glCheck();
}
//-------------------------------------------
void BdSurfVr::reshape(int w, int h)
{
  _ws = Index2(w,h);
  //VIEW 
  glViewport(0,0,w,h);
  //PROJECTION
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if(_proj==1) {
	 gluPerspective(30, double(w)/double(h), 0.1, 100);
  } else {
	 Point3 r = (_bbmax-_bbmin)/2.0 * 1.5;  double q = double(w)/double(h);  glOrtho(-r[1]*q, r[1]*q, -r[1], r[1], 0, 10);
  }
  glMatrixMode(GL_MODELVIEW);
}
//-------------------------------------------
void BdSurfVr::display()
{
  //clear
  glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );  
  //glClearColor(0.3, 0.5, 0.9, 0.0);
  glClearDepth(1.0f);
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glEnable(GL_DEPTH_TEST);
  //modelview
  glLoadIdentity();
  gluLookAt(0,0,_eyepos, //eye
				0,0,0, //look at
				0,1,0); //up dir
  GLfloat m[4][4];
  build_rotmatrix(m, _quat);
  glMultMatrixf(&m[0][0]);
  Point3 t = -(_bbmin+_bbmax)/2.0;
  glTranslatef(t[0],t[1],t[2]);
  glCheck();
  //render
  double zOffset = 1e-4;
  if(_vrop & RENDER_SURF) {
	 //--------------------------
	 glDepthRange(zOffset,1.0);//	 	 if(_mapOn==1) cubeMapOn();	 //glCullFace(GL_BACK);	 glEnable(GL_CULL_FACE);
	 glEnable(GL_DEPTH_TEST);
	 glColor3f(0.5,0.5,1);
	 glEnable(GL_LIGHTING);
	 for(int V=0; V<_Vfpos.size(); V++)
		for(int f=0; f<_Vfpos[V].size(); f++) {
		  NumMat<Point3>& cmpos = _Vfpos[V][f];
		  NumMat<Point3>& cmnor = _Vfnor[V][f];
		  int HLF = cmpos.m()-1;
		  for(int j=0; j<HLF; j++) {
			 glBegin(GL_QUAD_STRIP);
			 for(int i=0; i<HLF+1; i++) {
				glNormal3dv(cmnor(i,j+1));				  				glVertex3dv(cmpos(i,j+1));
				glNormal3dv(cmnor(i,j));				  				glVertex3dv(cmpos(i,j));
			 }
			 glEnd();
		  }
		}
	 glFlush();
  }
  glCheck();
  if(_vrop & RENDER_FRAME) {
	 //--------------------------
	 glDepthRange(0,1.0-zOffset);
	 glDisable(GL_LIGHTING);
	 glColor3f(0.0, 0.0, 0.0);
	 glLineWidth(1.0);
	 for(int V=0; V<_Vfpos.size(); V++) {
		for(int f=0; f<_Vfpos[V].size(); f++) {
		  NumMat<Point3>& cmpos = _Vfpos[V][f];		  //NumMat<Point3>& cmnor = _Vfnor[V][f];
		  for(int j=0; j<cmpos.n(); j++) {
			 assert(!glCheck());
			 glBegin(GL_LINE_STRIP);
			 for(int i=0; i<cmpos.m(); i++) {
				glVertex3dv(cmpos(i,j));
			 }
			 glEnd();
		  }
		  for(int i=0; i<cmpos.m(); i++) {
			 assert(!glCheck());
			 glBegin(GL_LINE_STRIP);
			 for(int j=0; j<cmpos.n(); j++) {
				glVertex3dv(cmpos(i,j));
			 }
			 glEnd();
		  }
		}
		glFlush();
	 }
  }
  if(_vrop & RENDER_CTRL) {
	 //--------------------------
	 glDisable(GL_LIGHTING);
	 //face
	 glDepthRange(zOffset,1.0);
	 glColor3f(1.0, 1.0, 1.0);
	 //glColor3f(0.7, 0.7, 0.7);
	 glBegin(GL_QUADS);
	 GpMesh& gpmesh = _bdsurf->gpmesh();
	 for(int F=0; F<gpmesh.numf(); F++) {
		vector<int> Vs(4);		gpmesh.F2Vs(F, Vs);
		for(int v=0; v<4; v++) {
		  int V = Vs[v];
		  glVertex3dv(gpmesh.points()[V].array());
		}
	 }
	 glEnd();
	 //line
	 glDepthRange(0,1.0-zOffset);
	 glColor3f(0.0, 0.0, 0.0);
	 glLineWidth(4.0);
	 for(int F=0; F<gpmesh.numf(); F++) {
		vector<int> Vs(4);		gpmesh.F2Vs(F, Vs);
		glBegin(GL_LINE_LOOP);
		for(int v=0; v<4; v++) {
		  int V = Vs[v];
		  glVertex3dv(gpmesh.points()[V].array());
		}
		glEnd();
	 }
	 glFlush();	 
  }
  glCheck();
  //swap buffer
  glutSwapBuffers();  //_rendered = true;
}
//-------------------------------------------
void BdSurfVr::mouse(int button, int state, int x, int y)
{
  _startdrag = Index2(x,y);
  for(int i=0; i<4; i++) _startquat[i] = _quat[i];   //_startquat = _quat; //copy rotation first
  _mouse_down[button] = (state==GLUT_DOWN);
}
//-------------------------------------------
void BdSurfVr::motion(int x, int y)
{
  int w=_ws[0], h=_ws[1];
  int sx = _startdrag[0];  int sy = _startdrag[1];
  double x1 = 2.0*   sx/(w-1.0)  - 1.0;  double y1 = 2.0*(1-sy/(h-1.0)) - 1.0;
  double x2 = 2.0*    x/(w-1.0)  - 1.0;  double y2 = 2.0*(1- y/(h-1.0)) - 1.0;
  if(       _mouse_down[0]) { //left down, rotate, non-accumulative
	 float curspin[4];
	 trackball(curspin, x1, y1, x2, y2);
	 add_quats(curspin, _startquat, _quat);
  } else if(_mouse_down[2]) { //right down, zoom, accumulative
	 _eyepos += 1.0*(y2-y1);
	 _startdrag = Index2(x,y);
  }
}
//-------------------------------------------
void BdSurfVr::keyboard(unsigned char key, int , int )
{
  switch(key) {
  case 'q':	 exit(0);	 break;
  case 's': _vrop = _vrop ^ RENDER_SURF; break;
  case 'f': _vrop = _vrop ^ RENDER_FRAME; break;
  case 'c': _vrop = _vrop ^ RENDER_CTRL; break;
  case 'o': trackball(_quat, 0, 0, 0, 0); _eyepos = (_bbmax[1]-_bbmin[1])/tan(M_PI/6.0)*1.5; break;
  case 'p':
	 _proj = 1-_proj; glMatrixMode(GL_PROJECTION);
	 glLoadIdentity();
	 if(_proj==1) {		gluPerspective(30, double(_ws[0])/double(_ws[1]), 0.1, 100);
	 } else {	 Point3 r = (_bbmax-_bbmin)/2.0 * 1.5;  double q = double(_ws[0])/double(_ws[1]);  glOrtho(-r[1]*q, r[1]*q, -r[1], r[1], 0, 10);
	 }
	 glMatrixMode(GL_MODELVIEW);	 break;
  case 'w':
	 {
		char tmp[200]; sprintf(tmp, "IMAGE_%d.png", _cnt++);
		//dump image
		GLint viewport[4];  glGetIntegerv(GL_VIEWPORT, viewport);
		int w = viewport[2];		  int h = viewport[3];
		IImage image; image.resize(w,h);
		glReadBuffer(GL_FRONT);		glReadPixels(0, 0, w, h, GL_BGRA, GL_UNSIGNED_BYTE, image.data());
		unsigned char* p = (unsigned char*) image.data();
		unsigned char* buf = new unsigned char[4 * w];
		for (unsigned j = 0; j < h/2; ++j) {
		  unsigned char* src = p + j * 4 * w;		  unsigned char* dst = p + (h - 1 - j) * 4 * w;
		  memcpy(buf, src, 4 * w);		  memcpy(src, dst, 4 * w);		  memcpy(dst, buf, 4 * w);
		}
		image.save(tmp); cerr<<"    save image file"<<tmp<<endl;
	 } break;
  case 'h':
	 {
		cerr<<"control keys:\n";
		cerr<<"q: quit\n";
		cerr<<"s: surface rendering on/off\n";
		cerr<<"f: frame rendering on/off\n";
		cerr<<"c: control mesh rendering on/off\n";
		cerr<<"o: back to the default viewing position\n";
		cerr<<"p: ortho/perspective projection\n";
		cerr<<"w: dump image\n";
		cerr<<"left  mouse key: rotate\n";
		cerr<<"right mouse key: scale\n";
		cerr<<"h: help\n";
	 } break;
  }
  glutPostRedisplay();
}
//-------------------------------------------
void BdSurfVr::idle()
{
  glutPostRedisplay();
}

//-----------------------------------------------------------
//static stuff
BdSurfVr* BdSurfVr::_self=NULL;
void BdSurfVr::reshapeWrapper(int w,int h) {  _self->reshape(w,h);}
void BdSurfVr::displayWrapper() {  _self->display();}
void BdSurfVr::mouseWrapper(int b,int s,int x,int y) {  _self->mouse(b,s,x,y);}
void BdSurfVr::motionWrapper(int x, int y) {  _self->motion(x,y);}
void BdSurfVr::keyboardWrapper(unsigned char c, int x, int y) {  _self->keyboard(c, x, y);}
void BdSurfVr::idleWrapper() {  _self->idle();}


/*
//MINE
GLfloat lightPos0[4] = { 1, 1, 1, 0 };
GLfloat lightSpe0[4] = { 1.5, 1.5, 1.5, 1 };
GLfloat lightDif0[4] = { 0.5, 0.5, 0.5, 1 };
GLfloat lightAmb0[4] = { 0.1, 0.1, 0.1, 1 };
glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpe0);
glLightfv(GL_LIGHT0, GL_DIFFUSE,  lightDif0);
glLightfv(GL_LIGHT0, GL_AMBIENT,  lightAmb0);
glEnable(GL_LIGHT0);
//GLfloat lightPos1[4] = { 1, 1, 0, 0 };
//GLfloat lightCol1[4] = { 0.375, 0.375, 0.375, 1 };
//glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);
//glLightfv(GL_LIGHT1, GL_SPECULAR, lightCol1);
//glLightfv(GL_LIGHT1, GL_DIFFUSE, lightCol1);
//glLightfv(GL_LIGHT1, GL_AMBIENT, lightCol1);
//glEnable(GL_LIGHT1);
*/
//glutInitWindowPosition(0, 0);  //glutInitWindowSize(_ws[0], _ws[1]);  //glutCreateWindow("BdSurfVr");  glCheck();
