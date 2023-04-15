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
#ifndef _BDSURF_HPP_
#define _BDSURF_HPP_

#include "common/nummat.hpp"
#include "ccsurf.hpp"

class Vxy;
class Vfcd;
class Fcd;

//-------------------------------------------------
//blended surface
class BdSurf: public ComObject {
public:
  //------------------
  enum {	 EVAL_VL=1,	 EVAL_FD=2,	 EVAL_SD=4  };
  typedef NumMat<Point3> Pt3NumMat;
  //------------------
  class Vxy {
  public:
	 int V;
	 double xy[2];
  public:
	 Vxy(): V(-1) {;}
	 Vxy(int tV, double* txy): V(tV) { 		xy[0]=txy[0]; xy[1]=txy[1];   }
	 Vxy(const Vxy& c): V(c.V) { 		xy[0]=c.xy[0]; xy[1]=c.xy[1];   }
	 Vxy& operator=(const Vxy& c) { 		V=c.V; xy[0]=c.xy[0]; xy[1]=c.xy[1]; return *this; 	 }
  };
  class Vfcd {
  public:
	 int V;
	 int f;
	 double cd[2];
  public:
	 Vfcd(): V(-1), f(-1) {;}
	 Vfcd(int tV, int tf, double* tcd): V(tV), f(tf) {	 cd[0]=tcd[0]; cd[1]=tcd[1];  }
	 Vfcd(const Vfcd& c): V(c.V), f(c.f) {	 cd[0]=c.cd[0]; cd[1]=c.cd[1];  }
	 Vfcd& operator=(const Vfcd& c) { 	 V=c.V; f=c.f; cd[0]=c.cd[0]; cd[1]=c.cd[1]; return *this;   }
  };
  /*
  //------------------
  class SData { //sampling data
  public:
	 double _init; //starting position
	 double _step; //step
	 int _nspl; //number of samples
	 BolNumMat _bvlmat; //for square
	 Pt3NumMat _posmat; //position info
	 };*/
  //------------------
  typedef pair<int,int> intpair;
  typedef vector< CCRect<Point3> > DataLevel;
protected:
  CCSurf _ccsurf; //a CC subdivision surface  //GpMesh _gpmesh;
  vector< NumVec<Point3> > _cgr;
  double _EVAL_UB; //THE RANGE WHERE THE SURFACE IS BLENDED
  double _INTE_UB; //THE BOUND FOR INTERNAL EVALUATION //the maximum range allows for internal evaluation
  double _CONS_UB; //THE BOUND USED FOR BLEND THE SURFACE
  //NOT USED ANYMORE
  int _ctrllvl; //control level, where points are used to approximate around extraordinary vertex
  int _pouctrl; //what partition of unity to use, C^\infty or Spline based
  int _fourgen; //control for vertex valence 4
  int _chttyp; //chart type: FULLYCOMPLEX = 0, CHARACTERISTIC = 1, EQUALDISTANCE = 2
  int _bsstyp; //basis type: POLY = 0, TENSOR_IN_POLAR = 1, ...
  int _stpwts;  //setup weight scheme
  //---------fast eval data
  //for each chart
  //vector<SData>    _sdatavec;
  Point3 _bbmin, _bbmax;
public:
  //------------------
  BdSurf(const string& p);
  ~BdSurf();
  int setup(map<string,string>& opts);  //int setup(istream&, double h);
  Point3 ctr() { return gpmesh().ctr(); }
  void bbx(Point3& bbmin, Point3& bbmax) { gpmesh().bbx(bbmin,bbmax); }
  //evaluation function
  int eval(  int flags, int V, double* xy, Point3* res);
  int sgeval(int flags, int V, double* xy, Point3* res);
  
  CCSurf& ccsurf() { return _ccsurf; }
  GpMesh& gpmesh() { return _ccsurf.gpmesh(); }
  int numf() { return gpmesh().numf(); }
  int numv() { return gpmesh().numv(); } //number of charts
  int valence(int V) { return gpmesh().Vvalence(V); }
  
  int compose(int flags, int udof, double* U, double* C, double* R);
  int Vxy2fcd(int flags, int, double*, int&,  double*); //x,y should be away from (0,0)
  int Vfcd2xy(int flags, int, int, double*, double*); //cd, should be away from (0,0)
  int Vfcd2Fcd(int flags, int, int, double*, int&, double*); //no requirement
  int Fcdv2Vfcd(int flags, int, double*, int, int&, int&, double*); //no requirement
  int Vfcd2Pou(int flags, int, int, double*, double, double, double*); //no requirement
  
  int chartbound(int V, double range, double&);
  
  double EVAL_UB() { return _EVAL_UB; }
  double INTE_UB() { return _INTE_UB; }
  double CONS_UB() { return _CONS_UB; }  
  
protected:
  int pou1d(int flags, double, double, double, double*);
  int charteval(int flags, int V, double*, Point3*);
  bool needgen(int K) { return (K!=4 || _fourgen==1); }
};

//GENERATION OF SURFACE, USING BLENDIeNG [1/8, 7/8] 
//CHART RANGE, USE [1/4,3/4]
//internal evaluation function, 
//int internal_eval(  int flags, int V, double* xy, Point3* res);
//int internal_sgeval(int flags, int V, double* xy, Point3* res);



#endif
