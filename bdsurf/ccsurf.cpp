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
#include "ccsurf.hpp"
#include "ccsurfop.hpp"
#include "common/vecmatop.hpp"

using std::min;
using std::max;

// ---------------------------------------------------------------------- 
int CCSurf::setup(map<string,string>& opts)
{
  //1. gpmesh initialize
  iC( _gpmesh.setup(opts) );
  
  //2. add gpmesh points
  _posvec.resize(1);
  DataLevel& tp = _posvec[0];
  tp.resize(numf());
  vector<Point3>& gpp = _gpmesh.points();
  for(int F=0; F<numf(); F++) {
	 tp[F].resize(2,2);
	 vector<int> Vs; _gpmesh.F2Vs(F, Vs);
	 tp[F](0,0) = gpp[Vs[0]];
	 tp[F](1,0) = gpp[Vs[1]];
	 tp[F](1,1) = gpp[Vs[2]];
	 tp[F](0,1) = gpp[Vs[3]];
  }
  iC( CCSurfOp::share(_gpmesh, 0, tp) ); //LEXING: VERY IMPORTANT
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurf::eval(int flags, int F, double* cd, Point3* ret)
{
  double c=cd[0];  double d=cd[1];
  //if(!(c>=0. && c<=1. && d>=0. && d<=1.))	 cerr<<"BAD "<<cd[0]<<" "<<cd[1]<<endl;
  //int lastlvl = lastposlvl(); 
  int lastlvl = _posvec.size()-1;
  int sz = pow2(lastlvl);
  double ss = 1.0/sz;
  int v = 0;
  bool adj = false;
  if(       c<ss && d<ss) {
	 adj = true; v = 0;
  } else if(c>1-ss && d<ss) {
	 adj = true; v = 1;
  } else if(c>1-ss && d>1-ss) {
	 adj = true; v = 2;
  } else if(c<ss && d>1-ss) {
	 adj = true; v = 3;
  } else{
	 adj = false;
  }
  bool bad;
  if(adj==false) {
	 bad = false;
  } else {
	 vector<int> tmp; _gpmesh.F2Vs(F, tmp);
	 int V = tmp[v];
	 if(_gpmesh.Vvalence(V) == 4) {
		bad = false;
	 } else {
		bad = true;
	 }
  }
  //LEXING: FOR THE TIME BEING, this is good enough
  iA(bad==false);
  
  CCRect<Point3> PD = _posvec[lastlvl][F];  //Matrix<Point3> mp(sz+3, sz+3, false, PD.data()); //matrix
  int    mn[2];  mn[0] = sz+3;  mn[1] = sz+3;
  double ef[2];  ef[0] = double(sz+3)/double(sz);  ef[1] = double(sz+3)/double(sz);
  double ts[2];  ts[0] = cd[0]*sz;  ts[1] = cd[1]*sz;
  int    ij[2];  ij[0] = min(max(int(floor(ts[0])),0),sz-1);  ij[1] = min(max(int(floor(ts[1])),0),sz-1);
  double lf[2];  lf[0] = ts[0]-ij[0];	 lf[1] = ts[1]-ij[1];
  ij[0]++;  ij[1]++; 
  //CALL spline evaluation function
  iC( spev2d(flags, DMFLAG_CLOSED, 3, (double*)(PD.data()), mn, ef, ij, lf, (double*)ret) );
  
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurf::subdivide(int lvl)
{
  if(lvl==_posvec.size()-1)
	 _posvec.push_back( DataLevel() );
  iC( CCSurfOp::subdivide(_gpmesh, lvl, _posvec[lvl], _posvec[lvl+1]) );
  return (0);
}

/*// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "CCSurf::dump"
int CCSurf::dump(ostream& fout, int lvl)
{
  ebiFunctionBegin;
  ebiFunctionReturn(0);
}*/


