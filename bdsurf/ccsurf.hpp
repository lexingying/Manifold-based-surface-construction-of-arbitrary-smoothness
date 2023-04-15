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
#ifndef _CCSURF_HPP_
#define _CCSURF_HPP_

#include "gpmesh.hpp"
#include "common/nummat.hpp"
#include "ccsubmatlib.hpp"

//-------------------------------------------------
template <class F>
class CCRect {
public:
  int _m, _n;
  F* _data;
public:
  CCRect(int m=0, int n=0): _m(m), _n(n) {
	 _data = new F[(_m+2)*(_n+2)]; assert( _data!=NULL );
  }
  CCRect(const CCRect& C): _m(C._m), _n(C._n) {
	 _data = new F[(_m+2)*(_n+2)]; assert( _data!=NULL );
	 memcpy( _data, C._data, (_m+2)*(_n+2)*sizeof(F) );
  }
  ~CCRect() {
	 delete[] _data; _data = NULL;
  }
  CCRect& operator=(const CCRect& C) {
	 delete[] _data; _data = NULL;
	 _m = C._m; _n=C._n;
	 _data = new F[(_m+2)*(_n+2)]; assert( _data!=NULL );
	 memcpy( _data, C._data, (_m+2)*(_n+2)*sizeof(F) );
	 return *this;
  }
  void resize(int m, int n)  {
	 if(_m!=m || _n!=n) {
		delete[] _data; _data = NULL;
		_m = m; _n = n;
		_data = new F[(_m+2)*(_n+2)]; assert( _data!=NULL );
	 }
  }
  const F& operator()(int i, int j) const  { 
	 assert( i>=-1 && i<_m+1 && j>=-1 && j<_n+1 );
	 return _data[(i+1)+(j+1)*(_m+2)];
  }
  F& operator()(int i, int j)  { 
	 assert( i>=-1 && i<_m+1 && j>=-1 && j<_n+1 );
	 return _data[(i+1)+(j+1)*(_m+2)];
  }
  int m() const { return _m; }
  int n() const { return _n; }
  F* data() const { return _data; }
};

//-------------------------------------------------
//subdivision surface
class CCSurf: public ComObject {
public:
  enum { EVAL_VL=1, EVAL_FD=2, EVAL_SD=4 };  //enum { MAXLEVELNUM=7 };
  typedef pair<int,int> intpair;
  typedef vector< CCRect<Point3> > DataLevel; //vector of faces
  
protected:
  GpMesh _gpmesh;  //int _inlvl; // input level
  vector<DataLevel> _posvec; //vector of levels

public:
  CCSurf(const string& p): ComObject(p), _gpmesh(p) {;}
  ~CCSurf() {;}
  
  int setup(map<string,string>& opts);  //int dump( ostream& ot, int lvl);
  int eval(int flags, int fid, double* uv, Point3*); //evaluation routine
  Point3 ctr() { return _gpmesh.ctr(); }
  void bbx(Point3& bbmin, Point3& bbmax) { _gpmesh.bbx(bbmin,bbmax); }
  
  //ops
  int subdivide(int lvl);
  
  //access
  int numf() { return _gpmesh.numf(); }
  int numv() { return _gpmesh.numv(); }
  GpMesh& gpmesh() { return _gpmesh; }  //CCSubMatLib* sm() { return _sm; }
  vector<DataLevel>& posvec() { return _posvec; }
};



#endif
