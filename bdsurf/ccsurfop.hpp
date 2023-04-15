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
#ifndef _CCSURFOP_HPP_
#define _CCSURFOP_HPP_

#include "ccsurf.hpp"



class CCSurfOp {
public:
  typedef pair<int,int> intpair;
  typedef CCSurf::DataLevel DataLevel;
  
  static int share(    GpMesh&, int lvl, DataLevel& wk);
  
  static int midsub(   GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
  static int subdivide(GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
  static int limit(    GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
  static int normal(   GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
  static int du(       GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
  static int dv(       GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
  static int duu(      GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
  static int duv(      GpMesh&, int lvl, DataLevel& fm, DataLevel& to);
  static int dvv(      GpMesh&, int lvl, DataLevel& fm, DataLevel& to);

  
  static int eno2index(int lvl, int eno, int t, int& i, int& j);
  //static int vno2index(int lvl, int vno, int &i, int& j);
  
  static int getcenter( GpMesh&, int lvl, DataLevel& wk, int V, Point3&);
  static int putcenter( GpMesh&, int lvl, DataLevel& wk, int V, Point3&);
  static int getonering(GpMesh&, int lvl, DataLevel& wk, int V, vector<Point3>&);
  static int putonering(GpMesh&, int lvl, DataLevel& wk, int V, vector<Point3>&);
  //static int gettworing(GpMesh&, int lvl, DataLevel& wk, int V, vector<Point3>&);
  //static int puttworing(GpMesh&, int lvl, DataLevel& wk, int V, vector<Point3>&);
};

//static int share( GpMesh&, int lvl, DataLevel& wk);
//static int buffer(GpMesh&, int lvl, DataLevel& fm, int F, Matrix<Point3>& to);
//todo: addto, scale, assign
//static int rotate(   GpMesh&, int lvl, DataLevel& wk, Quaternion r);
//static int shift(    GpMesh&, int lvl, DataLevel& wk, Point3 s);



#endif
