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
#ifndef _BDSURFVR_HPP_
#define _BDSURFVR_HPP_

#include "bdsurf.hpp"
#include "common/vec2t.hpp"
#include "vr.hpp"

//------------------------------------
class BdSurfVr {
public:
  enum { RENDER_SURF = 1, RENDER_FRAME = 2, RENDER_CTRL = 4 };
protected:
  //object
  BdSurf* _bdsurf;
  vector< vector< NumMat<Point2> > > _Vfxy;
  vector< vector< NumMat<Point3> > > _Vfpos;
  vector< vector< NumMat<Point3> > > _Vfnor;
  int _cnt;
  //size
  Point3 _bbmin, _bbmax; //bounding box info
  Index2 _ws; //window size
  float _eyepos;
  float _quat[4];
  int   _proj;
  //control, mouse, quat, ...
  int _mouse_down[3];
  int _vrop;
  float _startquat[4];
  Index2 _startdrag;
  //some constant
  static BdSurfVr* _self;  //bool _rendered;
public:
  BdSurfVr(int* argc, char** argv, BdSurf* bdsurf, int lvl);
  void reshape(int,int);
  void display();
  void mouse(int,int,int,int);
  void motion(int,int);
  void keyboard(unsigned char, int, int);
  void idle();
  static void reshapeWrapper(int,int);
  static void displayWrapper();
  static void mouseWrapper(int,int,int,int);
  static void motionWrapper(int,int);
  static void keyboardWrapper(unsigned char,int,int);
  static void idleWrapper();
};

#endif
