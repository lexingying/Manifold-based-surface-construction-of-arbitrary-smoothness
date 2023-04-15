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
#include "bdsurf.hpp"
#include "ccsubmatlib.hpp"
#include "ccsurf.hpp"
#include "ccsurfop.hpp"
#include "common/vecmatop.hpp"
#include "common/vec2t.hpp"

using std::min;
using std::max;
using std::abs;
using std::set;
using std::cerr;
using std::endl;

// ---------------------------------------------------------------------- 
BdSurf::BdSurf(const string& p): ComObject(p), _ccsurf(p)
{
  _EVAL_UB = 3.0/4.0;
  _INTE_UB = 63.0/64.0;
  _CONS_UB = 7.0/8.0;
  
  //DEFAULT VALUES
  _ctrllvl = 2;
  _pouctrl = 1;
  _fourgen = 0;
  _chttyp  = 1;
  _bsstyp  = 0;
  _stpwts  = 0;
}
// ---------------------------------------------------------------------- 
BdSurf::~BdSurf()
{
}

// ---------------------------------------------------------------------- 
int BdSurf::setup(map<string,string>& opts)
{
  //1. setup ccsurf
  iC( _ccsurf.setup(opts) );  iA(_ccsurf.posvec().size()==1); //LEXING: no multires version
  for(int i=0; i<3; i++) {
	 iC( _ccsurf.subdivide(i) );
  }
  //2. construction of _cgr
  int numv = gpmesh().numv();
  int ctrlsize = pow2(_ctrllvl);  iA(_ctrllvl==2);
  double ctrlstep = 1.0/ctrlsize;
  DataLevel& ctrlpos = _ccsurf.posvec()[_ctrllvl]; //control pos
  DataLevel  ctrllmt;  iC( CCSurfOp::limit(gpmesh(), _ctrllvl, ctrlpos, ctrllmt) );
  DataLevel  ctrldu;   iC( CCSurfOp::du(   gpmesh(), _ctrllvl, ctrlpos, ctrldu) );
  DataLevel  ctrldv;   iC( CCSurfOp::dv(   gpmesh(), _ctrllvl, ctrlpos, ctrldv) );
  
  _cgr.resize(numv);
  for(int V=0; V<numv; V++) {
	 int K = valence(V);	 //LEXING: FOUR
	 if(K==0) continue; //unused vertices
	 if(needgen(K)==false) continue;	//unused or regular vertices without 4gen
	 //2.1. gather	 //vector<Vfcd> tmpVfcd;
	 vector<Point2> tmpcd;  //weight matrix
	 vector<Point2> tmpxy;	//xy position	 //vector<int>    tmptyp; //0--value, 1--du, 2--dv
	 vector<Point3> tmpval;	//	 vector<Point3> tmpdu;	 vector<Poitn3> tmpdv;
	 set<intpair> Vijset;
	 int ss = (int)ceil( _CONS_UB * (ctrlsize) ); iA(ss==ctrlsize);
	 for(int f=0; f<K; f++) {
		for(int j=0; j<ss; j++)
		  for(int i=0; i<ss; i++) {
			 //----------remove redundancy
			 if(i==0 && j==0 && f!=0) continue; 
			 if(i==0 && j>0) continue;
			 //----------get position in different parameterization
			 double cd[2];			 cd[0] = i*ctrlstep;			 cd[1] = j*ctrlstep; //V position
			 double xy[2];			 iC( Vfcd2xy(EVAL_VL, V, f, cd, xy) );     //C position
			 int F;
			 double st[2];			 iC( Vfcd2Fcd(EVAL_VL, V, f, cd, F, st) ); //F position
			 int gh[2];			 gh[0] = (int)round(st[0]*ctrlsize);			 gh[1] = (int)round(st[1]*ctrlsize);
			 //
			 int MUL = pow2(20);
			 intpair tVij( (int)round(xy[0]*MUL), (int)round(xy[1]*MUL) );
			 if(Vijset.find(tVij)==Vijset.end()) { //not there yet
				Vijset.insert(tVij); 				//mark this point as processed
				//1. pos
				tmpcd.push_back( Point2(cd) );
				tmpxy.push_back( Point2(xy) );				//tmptyp.push_back(0);
				tmpval.push_back(ctrllmt[F](gh[0],gh[1]));
			 } else {				//REDUDANCY IS ALREADY REMOVED
				iA(0);
			 }
		  }
	 }
	 //LEX cerr<<" NUMBER OF CONSTRAINTS "<<tmpxy.size()<<endl;
	 //2.2. get associated weights
	 int WTS = _stpwts; iA(WTS==0); //WEIGHT CHANGE DISABLED
	 vector<double> tmpsqD(tmpcd.size(), 1);	 iA(tmpsqD.size()==tmpcd.size());
	 //2.3. solve least square system
	 int cm = tmpxy.size();
	 assert(_bsstyp==0);
	 //--------------------------------
	 int DEG;	 DEG = K+2;
	 int cn = (DEG+1)*DEG/2;		//LEX cerr<<"DEG "<<DEG<<" cm "<<cm<<" cn "<<cn<<endl;
	 DblNumMat M(cm,cn);
	 for(int i=0; i<cm; i++) {
		double x = tmpxy[i](0);		  double y = tmpxy[i](1);
		//construct prep
		double xs[40], ys[40];
		xs[0] = 1;		ys[0] = 1;
		for(int d=1; d<DEG; d++) {
		  xs[d] = xs[d-1]*x;		  ys[d] = ys[d-1]*y;
		}
		//----------------
		int j = 0;
		for(int d=0; d<DEG; d++)
		  for(int t=0; t<=d; t++)
			 M(i,j++) = xs[d-t] * ys[t];
		iA(j==cn);
		for(int j=0; j<cn; j++) M(i,j) *= tmpsqD[i];
	 }
	 DblNumMat piM(cn,cm);		iC( pinv(M, 1e-10, piM) ); //pinv( sqrt(D)*M )
	 _cgr[V].resize(cn);
	 for(int i=0; i<cn; i++) {
		Point3& cur = _cgr[V](i);
		cur = Point3(0.0);
		for(int j=0; j<cm; j++)
		  cur += piM(i,j) * tmpsqD[j] * tmpval[j];
	 }
  }
  /*
  //3. provide a dense sampling
  _sdatavec.resize(numv);
  _bbmin = Point3(SCL_MAX);  _bbmax = Point3(-SCL_MAX);
  cerr.precision(16);
  for(int V=0; V<numv; V++) {
	 double& init = _sdatavec[V]._init;	 double& step = _sdatavec[V]._step;	 int& nspl = _sdatavec[V]._nspl;
	 BolNumMat& bvlmat = _sdatavec[V]._bvlmat;	 Pt3NumMat& posmat = _sdatavec[V]._posmat;
	 //1. decide size
	 double bnd;	 iC( chartbound(V, EVAL_UB(), bnd) );
	 double xy[2];	 xy[0] = 0;	 xy[1] = 0;
	 Point3 ret[3];	 iC( internal_eval(EVAL_VL|EVAL_FD, V, xy, ret) );
	 double tmp = spacing / max(ret[1].l2(), ret[2].l2());
	 int tnum = int(ceil(bnd/tmp));	 //cerr<<bnd/tmp-1e-9<<" "<<tnum<<" | ";
	 step = bnd/tnum;
	 init = -(tnum+1)*step;
	 nspl = 2*tnum+3;
	 posmat.resize(nspl,nspl);
	 BolNumMat pvlmat(nspl, nspl); clear(pvlmat); 
	 for(int i=0; i<nspl; i++)
		for(int j=0; j<nspl; j++) {
		  double xy[2];		  xy[0] = init + i*step;		  xy[1] = init + j*step;
		  int f;		  double cd[2];		  iC( Vxy2fcd(EVAL_VL, V, xy, f, cd) );
		  pvlmat(i,j) = (cd[0]<INTE_UB() && cd[1]<INTE_UB());
		  if(pvlmat(i,j)==true) { //if nice, evaluation
			 iC( internal_eval(EVAL_VL, V, xy, &(posmat(i,j))) );
			 _bbmin = min(_bbmin, posmat(i,j));			 _bbmax = max(_bbmax, posmat(i,j)); //LEXING: compute bounding box
		  }
		}
	 bvlmat.resize(nspl-1,nspl-1); clear(bvlmat);
	 for(int i=1; i<nspl-2; i++)
		for(int j=1; j<nspl-2; j++) {
		  bool good = true;
		  for(int ni=i-1; ni<i+3; ni++)
			 for(int nj=j-1; nj<j+3; nj++)
				good = good && pvlmat(ni,nj);
		  bvlmat(i,j) = good; //good if all 4by4 grid is good
		}
	 cerr<<nspl-3<<" ";
  }
  cerr<<endl;
  */
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::chartbound(int V, double range, double& rad)
{
  double bbd = 0;
  int K = valence(V);  //cerr<<"valence "<<K<<endl;
  DblNumVec& D = (CCSubMatLib::getlibptr())->chooseMatrix(CCSubMatLib::INTERIOR, K).D();
  iA(_chttyp==1);
  // constants
  double p, q;
  if(       _chttyp==0) { //FULLY COMPLEX
	 p = 4.0/K;	 q = 4.0/K;
  } else if(_chttyp==1) { //CHARACTERISTIC MAP
	 p = log(1.0/D(1))/log(2.0);	 q = 4.0/K;
  } else if(_chttyp==2) { //UNIFORM SPACING
	 p = 1.0;	 q = 4.0/K;
  }
  //evaluate
  double r = sqrt(2.0)*range;	 double a = M_PI/4.0;
  double rp = pow(r,p);	 double aq = a*q;
  double xn = rp*cos(aq);	 double yn = rp*sin(aq);
  for(int f=0; f<K; f++) {
	 double ang = (2*M_PI*f)/K;
	 double xt = xn*cos(ang)-yn*sin(ang);		double yt = xn*sin(ang)+yn*cos(ang);
	 bbd = max(max(bbd,abs(xt)),abs(yt));
  }
  rad = bbd;
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::eval(int flags, int V, double* xy, Point3* ret)
{
  iA( (flags==EVAL_VL) || (flags==EVAL_VL|EVAL_FD) || (flags==EVAL_VL|EVAL_FD|EVAL_SD) );
  
  int Vi = V;  double* xyi = xy;
  int fi;  double cdi[12];  iC( Vxy2fcd(flags, Vi, xyi, fi, cdi) );
  iA(abs(cdi[0])<=_INTE_UB && abs(cdi[1])<=_INTE_UB); //IMPORTANT: check valid range for internal evaluation
  
  if(abs(cdi[0])<=1.0-_CONS_UB && abs(cdi[1])<=1.0-_CONS_UB) {
	 //--------------------one chart
	 int K = valence(Vi);
	 if(needgen(K)==true) { //extraordinary vertex or valence 4 with generation
		iC( charteval(flags, Vi, xyi, ret) );
	 } else { //regular
		int F;
		double st[12];		iC( Vfcd2Fcd(flags, Vi, fi, cdi, F, st) );
		//eval wrt to st
		Point3 buf[6];		iC( _ccsurf.eval(flags, F, st, buf) );
		//wrt to cdi
		Point3 tmp[6];
		iC( compose(flags, 3, (double*)buf, st, (double*)tmp) );
		//wrt to xyi
		iC( compose(flags, 3, (double*)tmp, cdi, (double*)ret) );
	 }
  } else {
	 //--------------------multiple chart
	 int F;
	 double st[12];	 iC( Vfcd2Fcd(flags, Vi, fi, cdi, F, st) );
	 //left weight
	 double lb[6];	 lb[0] = 1;	 for(int g=1; g<6; g++) lb[g] = 0;
	 //accumulation of results
	 Point3 ap[6];	 for(int i=0; i<6; i++)		ap[i]=Point3(0,0,0);
	 //extraordinary vertex
	 for(int v=0; v<4; v++) {
		//get cdj
		int Vj;		int fj;		double cdj[12];		iC( Fcdv2Vfcd(flags, F, st, v, Vj, fj, cdj) );
		int K = valence(Vj); //cerr<<"K "<<K<<endl;
		if( abs(cdj[0])<_CONS_UB && abs(cdj[1])<_CONS_UB && needgen(K) ) {
		  //within evaluation range and vertex is extraordinary
		  //get xyj
		  double xyj[12];		  iC( Vfcd2xy(flags, Vj, fj, cdj, xyj) );
		  //evaluate bj wrt cdj
		  double bj[6];		  iC( Vfcd2Pou(flags, Vj, fj, cdj, 1-_CONS_UB, _CONS_UB, bj) );
		  //evaluate pj wrt xyj
		  Point3 pj[6];		  iC( charteval(flags, Vj, xyj, pj) );
		  //evaluate pj wrt cdj
		  Point3 pjn[6]; 
		  iC( compose(flags, 3, (double*)pj, xyj, (double*)pjn) );
		  //evaluate pj wrt st
		  Point3 pjs[6];
		  iC( compose(flags, 3, (double*)pjn, cdj, (double*)pjs) );
		  //evaluate bj wrt st
		  double bjs[6];
		  iC( compose(flags, 1, (double*)bj,  cdj, (double*)bjs) );
		  //combine bjs+pjs into ap and lb
		  if(flags & EVAL_VL) {
			 ap[0] += bjs[0] * pjs[0];
			 lb[0] -= bjs[0];
		  }
		  if(flags & EVAL_FD) {
			 ap[1] += bjs[1] * pjs[0] + bjs[0] * pjs[1];
			 ap[2] += bjs[2] * pjs[0] + bjs[0] * pjs[2];
			 lb[1] -= bjs[1];
			 lb[2] -= bjs[2];
		  }
		  if(flags & EVAL_SD) {
			 ap[3] += bjs[3] * pjs[0] + 2*bjs[1]*pjs[1] + bjs[0] * pjs[3];
			 ap[4] += bjs[4] * pjs[0] + bjs[2]*pjs[1] + bjs[1]*pjs[2] + bjs[0] * pjs[4];
			 ap[5] += bjs[5] * pjs[0] + 2*bjs[2]*pjs[2] + bjs[0] * pjs[5];
			 lb[3] -= bjs[3];
			 lb[4] -= bjs[4];
			 lb[5] -= bjs[5];
		  }
		}
	 }
	 //do regular part if necessary
	 if(lb[0]>1e-8) {		//cerr<<"CCSURF EVAL"<<endl;
		Point3 lp[6];		iC( _ccsurf.eval(flags, F, st, lp) );
		if(flags & EVAL_VL) {
		  ap[0] += lb[0] * lp[0];
		}
		if(flags & EVAL_FD) {
		  ap[1] += lb[1] * lp[0] + lb[0] * lp[1];
		  ap[2] += lb[2] * lp[0] + lb[0] * lp[2];
		}
		if(flags & EVAL_SD) {
		  ap[3] += lb[3] * lp[0] + 2*lb[1]*lp[1] + lb[0] * lp[3];
		  ap[4] += lb[4] * lp[0] + lb[2]*lp[1] + lb[1]*lp[2] + lb[0] * lp[4];
		  ap[5] += lb[5] * lp[0] + 2*lb[2]*lp[2] + lb[0] * lp[5];
		}
	 }
	 Point3 apn[6];
	 iC( compose(flags, 3, (double*)ap, st, (double*)apn) );
	 iC( compose(flags, 3, (double*)apn, cdi, (double*)ret) );
  }
  
  return 0;
}
// ---------------------------------------------------------------------- 
int BdSurf::sgeval(int flags, int V, double* xy, Point3* ret)
{
  iA( (flags==EVAL_VL) || (flags==EVAL_VL|EVAL_FD) || (flags==EVAL_VL|EVAL_FD|EVAL_SD) );
  int Vi = V;
  double* xyi = xy;
  int fi;  double cdi[12];  iC( Vxy2fcd(flags, Vi, xyi, fi, cdi) );
  iA(abs(cdi[0])<=_INTE_UB && abs(cdi[1])<=_INTE_UB);
  int K = valence(Vi);
  if(needgen(K)==true) { //extraordinary vertex or valence 4 with generation
	 iC( charteval(flags, Vi, xyi, ret) );
  } else { //regular
	 int F;
	 double st[12];		iC( Vfcd2Fcd(flags, Vi, fi, cdi, F, st) );
	 //eval wrt to st
	 Point3 buf[6];		iC( _ccsurf.eval(flags, F, st, buf) );
	 //wrt to cdi
	 Point3 tmp[6];
	 iC( compose(flags, 3, (double*)buf, st, (double*)tmp) );
	 //wrt to xyi
	 iC( compose(flags, 3, (double*)tmp, cdi, (double*)ret) );
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::compose(int flags, int dof, double* U, double* C, double* R)
{
  //u--u(c,d), c--c(x,y), r--u(x,y)
  double* u  =&(U[0]);  double* uc =&(U[dof]);  double* ud =&(U[2*dof]);  double* ucc=&(U[3*dof]);  double* ucd=&(U[4*dof]);  double* udd=&(U[5*dof]);
  double& c  = C[0];  double& cx = C[2];  double& cy = C[4];  double& cxx= C[6];  double& cxy= C[8];  double& cyy= C[10];
  double& d  = C[1];  double& dx = C[3];  double& dy = C[5];  double& dxx= C[7];  double& dxy= C[9];  double& dyy= C[11];
  double* r  =&(R[0]);  double* ux =&(R[dof]);  double* uy =&(R[2*dof]);  double* uxx=&(R[3*dof]);  double* uxy=&(R[4*dof]);  double* uyy=&(R[5*dof]);
  //0 order
  if(flags & EVAL_VL) {
	 for(int d=0; d<dof; d++) r[d] = u[d];
  }
  if(flags & EVAL_FD) {
	 for(int i=0; i<dof; i++) {
		ux[i] = uc[i]*cx + ud[i]*dx;
		uy[i] = uc[i]*cy + ud[i]*dy;
	 }
  }
  if(flags & EVAL_SD) {
	 for(int i=0; i<dof; i++) {
		uxx[i] = ucc[i]*cx*cx + 2*ucd[i]*cx*dx + udd[i]*dx*dx       + uc[i]*cxx + ud[i]*dxx;
		uxy[i] = ucc[i]*cx*cy + ucd[i]*(cx*dy+cy*dx) + udd[i]*dx*dy + uc[i]*cxy + ud[i]*dxy;
		uyy[i] = ucc[i]*cy*cy + 2*ucd[i]*cy*dy + udd[i]*dy*dy       + uc[i]*cyy + ud[i]*dyy;
	 }
  }
  return (0);
}
// ---------------------------------------------------------------------- 
int BdSurf::Vxy2fcd(int flags, int V, double* xy, int& f, double* cd)
{
  iA( (flags==EVAL_VL) || (flags==EVAL_VL|EVAL_FD) || (flags==EVAL_VL|EVAL_FD|EVAL_SD) );
  int K = valence(V);
  DblNumVec& D = (CCSubMatLib::getlibptr())->chooseMatrix(CCSubMatLib::INTERIOR, K).D();
  // constants
  double p, q;
  if(       _chttyp==0) { //COMPLEX
	 p = K/4.0;	 q = K/4.0;
  } else if(_chttyp==1) { //CHARACTERISTIC MAP
	 p = log(2.0)/log(1.0/D(1));	 q = K/4.0;
  } else if(_chttyp==2) { //ISODISTANCE
	 p = 1.0;	 q = K/4.0;
  }
  
  //-------------------------------------------------------------------------
  if(abs(xy[0])<=SCL_EPS && abs(xy[1])<=SCL_EPS) {
	 //----------------------------	 //very close to zero point, map to face zero (defaultly)
	 double* val = cd;	 double* ext = cd+2;	 double* fnl = cd+6;
	 if(flags & EVAL_VL) {
		val[0] = 0;		val[1] = 0;
	 }
	 if(flags & EVAL_FD) { //ONLY VALID FOR VALENCE 4 
		if(K==4) {
		  ext[0] = 1;		ext[1] = 0;		ext[2] = 0;		ext[3] = 1;
		} else {
		  ext[0] = 0;		ext[1] = 0;		ext[2] = 0;		ext[3] = 0; //CASE INVALID
		}
	 }
	 if(flags & EVAL_SD) { //ONLY VALID FOR VALENCE 4
		if(K==4) {
		  for(int g=0; g<6; g++)		  fnl[g] = 0;
		} else {
		  for(int g=0; g<6; g++)		  fnl[g] = 0; //CASE INVALID
		}
	 }
	 f = 0; //VERY IMPORTANT
  } else {
	 //----------------------------
	 //compose
	 double x = xy[0];  double y = xy[1];
	 double ang = atan2(y,x);  if(ang<0) ang=ang+2*M_PI;
	 int i = min(int(floor(ang*K/(2*M_PI))), K-1);  //double rst = ang - 2*M_PI/K * double(i);
	 double R = -2*M_PI*double(i)/double(K);  //Matrix2 M2;	 M2(0,0) = cos(R);	 M2(0,1) = -sin(R);	 M2(1,0) = sin(R);	 M2(1,1) = cos(R);
	 //------------------
	 double st2xy[12];
	 if(flags & EVAL_VL) {
		double* tmp = st2xy;  tmp[0] = cos(R)*x - sin(R)*y;  tmp[1] = sin(R)*x + cos(R)*y;
	 }
	 if(flags & EVAL_FD) {
		double* tmp = st2xy+2;  tmp[0] = cos(R);  tmp[1] = sin(R);  tmp[2] = -sin(R);  tmp[3] = cos(R);
	 }
	 if(flags & EVAL_SD) {
		double* tmp = st2xy+6; for(int g=0; g<6; g++) tmp[g]=0;
	 }
	 //------------------
	 double ra2st[12];
	 double s = st2xy[0];  double t = st2xy[1];
	 double r2 = s*s+t*t;  double r = sqrt(r2);  double a = atan2(t,s);  double r3 = r2*r;  double r4 = r2*r2;
	 if(flags & EVAL_VL) {
		double* tmp = ra2st;  tmp[0] = r;  tmp[1] = a;
	 }
	 if(flags & EVAL_FD) {
		double* tmp = ra2st+2;  tmp[0] = s/r;  tmp[1] = -t/(r*r);  tmp[2] = t/r;  tmp[3] = s/(r*r); 
	 }
	 if(flags & EVAL_SD) {
		double* tmp = ra2st+6;
		tmp[0] = 1/r-s*s/r3;  tmp[1] = 2*t*s/r4;  tmp[2] = -s*t/r3;
		tmp[3] = (-s*s+t*t)/r4;  tmp[4] = 1/r-t*t/r3;  tmp[5] = -2*t*s/r4;
	 }
	 //------------------
	 double cd2ra[12];
	 double rp = pow(r,p); double rpm1 = pow(r,p-1);  double rpm2 = pow(r,p-2);  double aq = a*q;
	 if(flags & EVAL_VL) {
		double* tmp = cd2ra;  tmp[0] = rp*cos(aq);  tmp[1] = rp*sin(aq);
	 }
	 if(flags & EVAL_FD) {
		double* tmp = cd2ra+2;  tmp[0] = rpm1*p*cos(aq);  tmp[1] = rpm1*p*sin(aq);  tmp[2] = -rp*sin(aq)*q;  tmp[3] = rp*cos(aq)*q;
	 }
	 if(flags & EVAL_SD) {
		double* tmp = cd2ra+6;
		tmp[0] = rpm2*p*cos(aq)*(p-1);  tmp[1] = rpm2*p*sin(aq)*(p-1);  tmp[2] = -rpm1*p*sin(aq)*q;
		tmp[3] = rpm1*p*cos(aq)*q;  tmp[4] = -rp*cos(aq)*q*q;  tmp[5] = -rp*sin(aq)*q*q;
	 }
	 double ra2xy[12];
	 iC( compose(flags, 2, ra2st, st2xy, ra2xy) );
	 iC( compose(flags, 2, cd2ra, ra2xy, cd) );   //double cd2xy[12];
	 f = i; //VERY IMPORTANT
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int BdSurf::Vfcd2xy(int flags, int V, int f, double* cd, double* xy)
{
    iA( (flags==EVAL_VL) || (flags==EVAL_VL|EVAL_FD) || (flags==EVAL_VL|EVAL_FD|EVAL_SD) );
  int K = valence(V);//valence
  DblNumVec& D = (CCSubMatLib::getlibptr())->chooseMatrix(CCSubMatLib::INTERIOR, K).D();
  //constants
  double p, q;
  if(       _chttyp==0) { //FULLY COMPLEX
	 p = 4.0/K;	 q = 4.0/K;
  } else if(_chttyp==1) { //CHARACTERISTIC MAP
	 p = log(1.0/D(1))/log(2.0);	 q = 4.0/K;
  } else if(_chttyp==2) {
	 p = 1.0;	 q = 4.0/K;
  }
  //---------------------
  if( abs(cd[0])<=SCL_EPS && abs(cd[1])<=SCL_EPS ) {
	 //-------------------	 //CANNOT TRANSFROM AT C==0 && D==0	 //xy[0] = 0;	 xy[1] = 0;
	 double* val = xy;	 double* ext = xy+2;	 double* fnl = xy+6;
	 if(flags & EVAL_VL) {
		val[0] = 0;		val[1] = 0;
	 }
	 if(flags & EVAL_FD) {
		//TODO: make it complete
		assert(0);
	 }
	 if(flags & EVAL_SD) {
		//TODO: make it complete
		assert(0);
	 }
  } else {
	 double c = cd[0]; double d = cd[1];
	 //-------------------
	 double ra2cd[12];
	 double r2 = c*c+d*d;  double r = sqrt(r2);  double a = atan2(d,c);  double r3 = r2*r;  double r4 = r2*r2;
	 if(flags & EVAL_VL) {
		double* tmp = ra2cd;  tmp[0] = r;  tmp[1] = a;
	 }
	 if(flags & EVAL_FD) {
		double* tmp = ra2cd+2;  tmp[0] = c/r;  tmp[1] = -d/(r*r);  tmp[2] = d/r;  tmp[3] = c/(r*r); 
	 }
	 if(flags & EVAL_SD) {
		double* tmp = ra2cd+6;
		tmp[0] = 1/r-c*c/r3;  tmp[1] = 2*d*c/r4;  tmp[2] = -c*d/r3;
		tmp[3] = (-c*c+d*d)/r4;  tmp[4] = 1/r-d*d/r3;  tmp[5] = -2*d*c/r4;
	 }
	 //-------------------
	 double st2ra[12];
	 double rp = pow(r,p); double rpm1 = pow(r,p-1);  double rpm2 = pow(r,p-2);  double aq = a*q;
	 if(flags & EVAL_VL) {
		double* tmp = st2ra;  tmp[0] = rp*cos(aq);  tmp[1] = rp*sin(aq);
	 }
	 if(flags & EVAL_FD) {
		double* tmp = st2ra+2;  tmp[0] = rpm1*p*cos(aq);  tmp[1] = rpm1*p*sin(aq);  tmp[2] = -rp*sin(aq)*q;  tmp[3] = rp*cos(aq)*q;
	 }
	 if(flags & EVAL_SD) {
		double* tmp = st2ra+6;
		tmp[0] = rpm2*p*cos(aq)*(p-1);  tmp[1] = rpm2*p*sin(aq)*(p-1);  tmp[2] = -rpm1*p*sin(aq)*q;
		tmp[3] = rpm1*p*cos(aq)*q;  tmp[4] = -rp*cos(aq)*q*q;  tmp[5] = -rp*sin(aq)*q*q;
	 }
	 //-------------------
	 double xy2st[12];
	 double R = 2*M_PI*double(f)/double(K);
	 double s = st2ra[0];  double t = st2ra[1];
	 if(flags & EVAL_VL) {
		double* tmp = xy2st;  tmp[0] = cos(R)*s - sin(R)*t;  tmp[1] = sin(R)*s + cos(R)*t;
	 }
	 if(flags & EVAL_FD) {
		double* tmp = xy2st+2;  tmp[0] = cos(R);  tmp[1] = sin(R);  tmp[2] = -sin(R);  tmp[3] = cos(R);
	 }
	 if(flags & EVAL_SD) {
		double* tmp = xy2st+6; for(int g=0; g<6; g++) tmp[g]=0;
	 }
	 double st2cd[12];
	 iC( compose(flags, 2, st2ra, ra2cd, st2cd) );
	 iC( compose(flags, 2, xy2st, st2cd, xy) );
  }
  return (0);
}
// ---------------------------------------------------------------------- 
int BdSurf::Vfcd2Fcd(int flags, int V, int f, double* cd, int& F, double* st)
{
  iA( (flags==EVAL_VL) || (flags==EVAL_VL|EVAL_FD) || (flags==EVAL_VL|EVAL_FD|EVAL_SD) ); //LEXING
  
  //return val
  double* val = st;
  double* ext = st+2; //first order
  double* fnl = st+6; //second order
  
  int v; gpmesh().Vf2Fv(V,f,F,v);
  switch(v) {
  case 0:
	 if(flags & EVAL_VL) {
		val[0]=cd[0];		val[1]=cd[1];
	 }
	 if(flags & EVAL_FD) {
		ext[0] = 1;		ext[1] = 0;		ext[2] = 0;		ext[3] = 1;
	 }
	 if(flags & EVAL_SD) {
		for(int g=0; g<6; g++) fnl[g] = 0;
	 }
	 break;
  case 1:
	 if(flags & EVAL_VL) {
		val[0]=1-cd[1];		val[1]=cd[0];
	 }
	 if(flags & EVAL_FD) {
		ext[0] = 0;		ext[1] = 1;		ext[2] = -1;		ext[3] = 0;
	 }
	 if(flags & EVAL_SD) {
		for(int g=0; g<6; g++) fnl[g] = 0;
	 }
	 break;
  case 2:
	 if(flags & EVAL_VL) {
		val[0]=1-cd[0];		val[1]=1-cd[1];
	 }
	 if(flags & EVAL_FD) {
		ext[0] = -1;		ext[1] = 0;		ext[2] = 0;		ext[3] = -1;
	 }
	 if(flags & EVAL_SD) {
		for(int g=0; g<6; g++) fnl[g] = 0;
	 }
	 break;
  case 3:
	 if(flags & EVAL_VL) {
		val[0]=cd[1];		val[1]=1-cd[0];
	 }
	 if(flags & EVAL_FD) {
		ext[0] = 0;		ext[1] = -1;		ext[2] = 1;		ext[3] = 0;
	 }
	 if(flags & EVAL_SD) {
		for(int g=0; g<6; g++) fnl[g] = 0;
	 }
	 break;
  default:
	 iA(0);
  }
  return 0;
}

// ---------------------------------------------------------------------- 
#undef __FUNCT__
#define __FUNCT__ "BdSurf::Fcdv2Vfcd"
int BdSurf::Fcdv2Vfcd(int flags, int F, double* cd, int v, int& V, int& f, double* st)
{
  iA( (flags==EVAL_VL) || (flags==EVAL_VL|EVAL_FD) || (flags==EVAL_VL|EVAL_FD|EVAL_SD) ); //LEXING
    
  //return val
  double* val = st;
  double* ext = st+2;
  double* fnl = st+6;
  
  gpmesh().Fv2Vf(F, v, V, f);//intpair Vf = gpmesh().Fv2Vf(F,v);  V = Vf.first;  f = Vf.second;
  switch(v) {
  case 0:
	 if(flags & EVAL_VL) {
		val[0] = cd[0];		val[1] = cd[1];
	 }
	 if(flags & EVAL_FD) {
		ext[0] = 1;		ext[1] = 0;		ext[2] = 0;		ext[3] = 1;
	 }
	 if(flags & EVAL_SD) {
		for(int g=0; g<6; g++) fnl[g] = 0;
	 }
	 break;
  case 1:
	 if(flags & EVAL_VL) {
		val[0] = cd[1];		val[1] = 1-cd[0];
	 }
	 if(flags & EVAL_FD) {
		ext[0] = 0;		ext[1] = -1;		ext[2] = 1;		ext[3] = 0;
	 }
	 if(flags & EVAL_SD) {
		for(int g=0; g<6; g++) fnl[g] = 0;
	 }
	 break;
  case 2:
	 if(flags & EVAL_VL) {
		val[0] = 1-cd[0];		val[1] = 1-cd[1];
	 }
	 if(flags & EVAL_FD) {
		ext[0] = -1;		ext[1] = 0;		ext[2] = 0;		ext[3] = -1;
	 }
	 if(flags & EVAL_SD) {
		for(int g=0; g<6; g++) fnl[g] = 0;
	 }
	 break;
  case 3:
	 if(flags & EVAL_VL) {
		val[0] = 1-cd[1];		val[1] = cd[0];
	 }
	 if(flags & EVAL_FD) {
		ext[0] = 0;		ext[1] = 1;		ext[2] = -1;		ext[3] = 0;
	 }
	 if(flags & EVAL_SD) {
		for(int g=0; g<6; g++) fnl[g] = 0;
	 }
	 break;
  default:
	 iA(0);
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::Vfcd2Pou(int flags, int V, int f, double* cd, double LB, double UB, double* pou)
{
  iA( (flags==EVAL_VL) || (flags==EVAL_VL|EVAL_FD) || (flags==EVAL_VL|EVAL_FD|EVAL_SD) );
  double resu[3];  iC( pou1d(flags, cd[0], LB, UB, resu) );
  double resv[3];  iC( pou1d(flags, cd[1], LB, UB, resv) );
  double* val = pou;
  double* ext = pou+1;
  double* fnl = pou+3;
  if(flags & EVAL_VL) {
	 val[0] = resu[0]*resv[0];
  }
  if(flags & EVAL_FD) {
	 ext[0] = resu[1]*resv[0];
	 ext[1] = resu[0]*resv[1]; 
  }
  if(flags & EVAL_FD) {
	 fnl[0] = resu[2]*resv[0];
	 fnl[1] = resu[1]*resv[1];
	 fnl[2] = resu[0]*resv[2];
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int BdSurf::pou1d(int flags, double u, double LB, double UB, double* res)
{
  iA( (flags==EVAL_VL) || (flags==EVAL_VL|EVAL_FD) || (flags==EVAL_VL|EVAL_FD|EVAL_SD) );
  double* val = res;
  double* ext = res+1;
  double* fnl = res+2;
  double t = (u-LB)/(UB-LB);
  double ERR = 1e-7;
  if(       t<=0+ERR) {
	 //--------------------------------------------------
	 if(flags & EVAL_VL) val[0]=1;
	 if(flags & EVAL_FD) ext[0]=0;
	 if(flags & EVAL_SD) fnl[0]=0;
  } else if(t>=1-ERR) {
	 //--------------------------------------------------
	 if(flags & EVAL_VL) val[0]=0;
	 if(flags & EVAL_FD) ext[0]=0;
	 if(flags & EVAL_SD) fnl[0]=0;
  } else {
	 //--------------------------------------------------
	 if(       _pouctrl==0 ) {
		//--------------------------------
		double s = 1-t;
		double t2 = t*t;		double t3 = t2*t;		double t4 = t2*t2;
		double s2 = s*s;		double s3 = s2*s;		double s4 = s2*s2;
		double a = 2*exp(-1/t)/(t-1);
		double b = 2*exp(-1/s)/(s-1);
		double da =  a*(1/(t*t) - 1/(t-1));
		double db = -b*(1/(s*s) - 1/(s-1));
		double dda = a*(-4*t3+7*t2-4*t+1+2*t4)/((t-1)*(t-1))/t4;
		double ddb = b*(-4*s3+7*s2-4*s+1+2*s4)/((s-1)*(s-1))/s4;
		double ea = exp(a);
		double eb = exp(b);
		double f = ea;
		double g = ea+eb;
		double df = ea*da;
		double dg = ea*da + eb*db;
		double ddf = ea*da*da + ea*dda;
		double ddg = ea*da*da + ea*dda + eb*db*db + eb*ddb;
		if(flags & EVAL_VL) {
		  val[0] = f/g;
		}
		if(flags & EVAL_FD) {
		  ext[0] = (df/g - f/(g*g)*dg) / (UB-LB); //ext[0] = (dea*eab-ea*deab) / (eab*eab) * 1.0/(UB-LB); //scaling
		}
		if(flags & EVAL_SD) {
		  fnl[0] = (ddf/g - 2*df/(g*g)*dg + 2*f/(g*g*g)*dg*dg - f/(g*g)*ddg) / ((UB-LB)*(UB-LB));
		}
	 } else if(_pouctrl==1) {
		//--------------------------------
		double poucp[6] = {1,1,1,0,0,0}; //ctrl points for pou
		double x = t*3 + 1.0;
		int i = min(max(int(floor(x)),1),3);
		double f = x - i;		//double tmp[2]; //iC( spev1d(flags, _poucp, i, f, tmp) );//res[0] = tmp[0];res[1] = tmp[1] * 3 / (UB-LB);
		iC( spev1d(flags, DMFLAG_CLOSED, 1, poucp, 6, 2*(UB-LB), i, f, res) ); //cerr<<res[0]<<endl;
	 }
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int BdSurf::charteval(int flags, int V, double* xy, Point3* ret)
{
  iA( (flags==EVAL_VL) || (flags==EVAL_VL|EVAL_FD) || (flags==EVAL_VL|EVAL_FD|EVAL_SD) );
  int K = valence(V);  //iA(K!=4);
  NumVec<Point3>& vp = _cgr[V];
  double x = xy[0];	 double y = xy[1];
  //------------------------
  if(       _bsstyp==0) { //POLYNOMIAL
	 //---------------------------
	 int DEG;	 DEG = K+2;
	 int cn = (DEG+1)*DEG/2;	 iA(cn==vp.m());
	 double xs[40], ys[40];
	 xs[0] = 1;		ys[0] = 1;
	 for(int d=1; d<DEG; d++) {
		xs[d] = xs[d-1]*x;
		ys[d] = ys[d-1]*y;
	 }
	 double tmp[200];
	 Point3* cur = ret;
	 //value
	 if(flags & EVAL_VL) {
		{
		  int j = 0;
		  for(int d=0; d<DEG; d++) {
			 for(int t=0; t<=d; t++)
				tmp[j++] = xs[d-t] * ys[t];
		  }
		  iA(j==cn);
		  *cur = Point3( 0.0); //clear
		  for(int i=0; i<cn; i++)			 (*cur) += tmp[i] * vp(i);
		  cur = cur+1;
		}
	 }
	 //fd
	 if(flags & EVAL_FD) {
		{
		  int j = 0;
		  for(int d=0; d<DEG; d++) {
			 for(int t=0; t<d; t++)
				tmp[j++] = (d-t) * xs[d-t-1] * ys[t];
			 tmp[j++] = 0;
		  }
		  iA(j==cn);
		  *cur = Point3( 0.0); //clear
		  for(int i=0; i<cn; i++)			 (*cur) += tmp[i] * vp(i);
		  cur = cur+1;
		}
		{
		  int j = 0;
		  for(int d=0; d<DEG; d++) {
			 tmp[j++] = 0;
			 for(int t=1; t<=d; t++) {
				tmp[j++] = t * xs[d-t] * ys[t-1];
			 }
		  }
		  iA(j==cn);
		  *cur = Point3( 0.0); //clear
		  for(int i=0; i<cn; i++)			 (*cur) += tmp[i] * vp(i);
		  cur = cur+1;
		}
	 }
	 //sd
	 if(flags & EVAL_SD) {
		{
		  int j=0;
		  for(int d=0; d<DEG; d++) {
			 for(int t=0; t<=d; t++)
				if(d-t-2>=0) tmp[j++] = (d-t)*(d-t-1)*xs[d-t-2] * ys[t];
				else tmp[j++] = 0;
		  }
		  iA(j==cn);
		  *cur = Point3(0.0);
		  for(int i=0; i<cn; i++)			 (*cur) += tmp[i] * vp(i);
		  cur = cur+1;
		}
		{
		  int j = 0;
		  for(int d=0; d<DEG; d++) {
			 for(int t=0; t<=d; t++)
				if(d-t-1>=0 && t-1>=0) tmp[j++] = (d-t)*t*xs[d-t-1]*ys[t-1];
				else tmp[j++] = 0;
		  }
		  iA(j==cn);
		  *cur = Point3(0.0);
		  for(int i=0; i<cn; i++)			 (*cur) += tmp[i] * vp(i);
		  cur = cur+1;
		}
		{
		  int j=0;
		  for(int d=0; d<DEG; d++) {
			 for(int t=0; t<=d; t++)
				if(t-2>=0) tmp[j++] = t*(t-1)*xs[d-t]*ys[t-2];
				else tmp[j++] = 0;
		  }
		  iA(j==cn);
		  *cur = Point3(0.0);
		  for(int i=0; i<cn; i++)			 (*cur) += tmp[i] * vp(i);
		  cur = cur+1;
		}
	 }
  } else {
	 iA(0);
  }
  return 0;
}

