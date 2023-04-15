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
#include "ccsurfop.hpp"

// ---------------------------------------------------------------------- 
int CCSurfOp::share(GpMesh& gm, int lvl, DataLevel& wk)
{
  int numf = gm.numf();
  int size = pow2(lvl)+1;
  //--------------------
  //sharing edges
  for(int F=0; F<numf; F++) {
	 for(int e=0; e<4; e++) {
		int ci, cj; iC( eno2index(lvl,e,0,ci,cj) );
		int cdi,cdj;iC( eno2index(lvl,e,1,cdi,cdj) );
		cdi=cdi-ci; cdj=cdj-cj;
		ci=ci+cdj; cj=cj-cdi;
		iA( (ci==-1)||(cj==-1)||(ci==size)||(cj==size));
		CCRect<Point3>& CD = wk[F];
		
		int nF, ne;		gm.Fe2Fe(F, e, nF, ne); //intpair nFe = gm.Fe2Fe(F,e);		//int nF = nFe.first;		//int ne = nFe.second;
		
		int ni, nj; iC( eno2index(lvl,ne,size-1,ni,nj) );
		int ndi,ndj;iC( eno2index(lvl,ne,size-2,ndi,ndj) );
		ndi=ndi-ni; ndj=ndj-nj;
		ni=ni+ndj; nj=nj-ndi;
		iA( (ni==1)||(nj==1)||(ni==size-2)||(nj==size-2));
		CCRect<Point3>& ND = wk[nF];

		for(int k=0; k<size; k++)
		  CD(ci+k*cdi,cj+k*cdj) = ND(ni+k*ndi,nj+k*ndj);
	 }
  }
  //--------------------
  //sharing vertices
  for(int F=0; F<numf; F++) {
	 for(int e=0; e<4; e++) {
		int ci, cj; iC( eno2index(lvl,e,0,ci,cj) );
		int cdi,cdj;iC( eno2index(lvl,e,1,cdi,cdj) );
		cdi=cdi-ci; cdj=cdj-cj;
		ci=ci+cdj-cdi; cj=cj-cdi-cdj;
		CCRect<Point3>& CD = wk[F];
		
		int nF, ne;		gm.Fe2Fe(F, e, nF, ne);		//intpair nFe = gm.Fe2Fe(F,e);		int nF = nFe.first;		int ne = nFe.second;
		
		int ni, nj; iC( eno2index(lvl,ne,size-1,ni,nj) );
		int ndi,ndj;iC( eno2index(lvl,ne,size-2,ndi,ndj) );
		ndi=ndi-ni; ndj=ndj-nj;
		ni=ni+ndj-ndi; nj=nj-ndi-ndj;
		CCRect<Point3>& ND = wk[nF];
		CD(ci,cj) = ND(ni,nj);
	 }
  }
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::midsub(GpMesh&, int, DataLevel&, DataLevel&)
{
  //TODO
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::subdivide(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{
  CCSubMatLib& submat = *(CCSubMatLib::getlibptr());
  
  int numf = gm.numf();
  int numv = gm.numv();
  int size = pow2(lvl)+1;
  
  iA(fm.size()==numf);
  to.resize(numf);
  for(int F=0; F<numf; F++) {
	 iA(fm[F].m()==size && fm[F].n()==size);
	 to[F].resize(2*size-1,2*size-1); 	 //iA(to[F].size()==2*size-1);
  }
  
  //subdivide face  //Matrix<Point3> BD(size+2,size+2);
  for(int F=0; F<numf; F++) {
	 CCRect<Point3>& CD = fm[F];
	 CCRect<Point3>& ND = to[F];
	 {
		double cwt = 1.0/64.0;
		double ewt = 6.0/64.0;
		double swt = 36.0/64.0;
		for(int i=0; i<size; i++)
		  for(int j=0; j<size; j++)
			 ND(2*i,2*j) = 
				cwt*CD(i-1,j-1) + ewt*CD(i,j-1) + cwt*CD(i+1,j-1) +
				ewt*CD(i-1,j  ) + swt*CD(i,j  ) + ewt*CD(i+1,j  ) +
				cwt*CD(i-1,j+1) + ewt*CD(i,j+1) + cwt*CD(i+1,j+1);
	 }
	 {
		double cwt = 1.0/16.0;
		double ewt = 6.0/16.0;
		for(int i=0; i<size-1; i++)
		  for(int j=0; j<size; j++)
			 ND(2*i+1,2*j) = 
				cwt*CD(i,j-1) + cwt*CD(i+1,j-1) +
				ewt*CD(i,j  ) + ewt*CD(i+1,j  ) +
				cwt*CD(i,j+1) + cwt*CD(i+1,j+1);
		for(int i=0; i<size; i++)
		  for(int j=0; j<size-1; j++)
			 ND(2*i,2*j+1) = 
				cwt*CD(i-1,j  ) + ewt*CD(i,j  ) + cwt*CD(i+1,j  ) +
				cwt*CD(i-1,j+1) + ewt*CD(i,j+1) + cwt*CD(i+1,j+1);
	 }
	 {
		double cwt = 1.0/4.0;
		  for(int i=0; i<size-1; i++)
			 for(int j=0; j<size-1; j++)
				ND(2*i+1,2*j+1)= 
				  cwt*( CD(i,j  ) + CD(i+1,j  ) +
						  CD(i,j+1) + CD(i+1,j+1) );
	 }
  }
  
  //subdivide vertex
  for(int V=0; V<numv; V++) {
	 int K = gm.Vvalence(V);	 iA(K>=3 || K==0);
	 if(K==0) continue; //unused vertices
	 if(K!=4) {
		vector<Point3> onering;
		iC( getonering(gm, lvl, fm, V, onering) );
		Point3 center(0.0);
		DblNumMat& SM = submat.chooseMatrix(CCSubMatLib::INTERIOR, K).SM();
		for(int j=0; j<onering.size(); j++) {
		  center += SM(0,j) * onering[j];
		}
		iC( putcenter(gm, lvl+1, to, V, center ) );
	 }
  }
  
  //LEXING: VERY IMPORTANT
  iC( share(gm, lvl+1, to) );
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::limit(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{
  CCSubMatLib& submat = *(CCSubMatLib::getlibptr());
  
  int numf = gm.numf();
  int numv = gm.numv();
  int size = pow2(lvl)+1;

  iA(fm.size()==numf);
  to.resize(numf);
  for(int F=0; F<numf; F++) {
	 iA(fm[F].m()==size && fm[F].n()==size);
	 to[F].resize(size,size); //iA(to[F].size()==size);
  }
  
  for(int F=0; F<numf; F++) {
	 CCRect<Point3>& CD = fm[F];
	 CCRect<Point3>& ND = to[F];
	 
	 double cwt = 1.0/36.0;
	 double ewt = 4.0/36.0;
	 double swt = 16.0/36.0;
	 for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
		  ND(i,j)=
			 cwt*CD(i-1,j-1) + ewt*CD(i,j-1) + cwt*CD(i+1,j-1) +
			 ewt*CD(i-1,j  ) + swt*CD(i,j  ) + ewt*CD(i+1,j  ) +
			 cwt*CD(i-1,j+1) + ewt*CD(i,j+1) + cwt*CD(i+1,j+1);
  }
  
  for(int V=0; V<numv; V++) {
	 int K = gm.Vvalence(V);	 iA(K>=3 || K==0);
	 if(K==0) continue; //unused vertices
	 if(K!=4) {
		vector<Point3> onering;
		iC( getonering(gm, lvl, fm, V, onering) );
		Point3 limit( 0.0);
		DblNumMat& IV = submat.chooseMatrix(CCSubMatLib::INTERIOR, K).IV();
		for(int j=0; j<onering.size(); j++) {
		  limit += IV(0,j) * onering[j];
		}
		iC( putcenter(gm, lvl, to, V, limit) );
	 }
  }
  
  //LEXING: VERY IMPORTANT
  iC( share(gm, lvl, to) );
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::normal(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{
  CCSubMatLib& submat = *(CCSubMatLib::getlibptr());

  int numf = gm.numf();
  int numv = gm.numv();
  int size = pow2(lvl)+1;
  
  iA(fm.size()==numf);
  to.resize(numf);
  for(int F=0; F<numf; F++) {
	 iA(fm[F].m()==size && fm[F].n()==size);
	 to[F].resize(size,size); //iA(to[F].size()==size);
  }

  for(int F=0; F<numf; F++) {
	 CCRect<Point3>& CD = fm[F];
	 CCRect<Point3>& ND = to[F];
	 
	 for(int i=0; i<size; i++)
		for(int j=0; j<size; j++) {
		  Point3 du = 
			 ( (CD(i+1,j  )-CD(i-1,j  ))*4.0 + 
				(CD(i+1,j+1)-CD(i-1,j+1)) +
				(CD(i+1,j-1)-CD(i-1,j-1)) )  / 12.0;
		  Point3 dv = 
			 ( (CD(i,  j+1)-CD(i,  j-1)  )*4.0 + 
				(CD(i+1,j+1)-CD(i+1,j-1)) +
				(CD(i-1,j+1)-CD(i-1,j-1)) ) / 12.0;
		  Point3 nor( cross(du,dv) ); nor=nor.dir();
		  ND(i,j) = nor;
		}
  }
  //corners
  for(int V=0; V<numv; V++) {
	 int K = gm.Vvalence(V);	 iA(K>=3 || K==0);
	 if(K==0) continue;
	 if(K!=4) {
		vector<Point3> onering;
		iC( getonering(gm, lvl, fm, V, onering) );
		DblNumMat& IV = submat.chooseMatrix(CCSubMatLib::INTERIOR, K).IV();
		Point3 du( 0.0);
		for(int j=0; j<onering.size(); j++) {
		  du += IV(1,j) * onering[j];
		}
		Point3 dv( 0.0);
		for(int j=0; j<onering.size(); j++) {
		  dv += IV(2,j) * onering[j];
		}
		Point3 nor( cross(du,dv) ); nor=nor.dir();
		iC( putcenter(gm, lvl, to, V, nor) );
	 }
  }
  
  //LEXING: VERY IMPORTANT
  iC( share(gm, lvl, to) );
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::du(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{
  CCSubMatLib& submat = *(CCSubMatLib::getlibptr());

  int numf = gm.numf();
  int numv = gm.numv();
  int size = pow2(lvl)+1;
  double h = 1.0/pow2(lvl);
  iA(fm.size()==numf);
  to.resize(numf);
  for(int F=0; F<numf; F++) {
	 iA(fm[F].m()==size && fm[F].n()==size);
	 to[F].resize(size,size); //iA(to[F].size()==size);
  }
  for(int F=0; F<numf; F++) {
	 CCRect<Point3>& CD = fm[F];
	 CCRect<Point3>& ND = to[F];
	 
	 for(int i=0; i<size; i++)
		for(int j=0; j<size; j++) {
		  ND(i,j) = 
			 ( (CD(i+1,j  )-CD(i-1,j  ))*4.0 + 
				(CD(i+1,j+1)-CD(i-1,j+1)) +
				(CD(i+1,j-1)-CD(i-1,j-1)) )  / (12.0*h);
		}
  }
  //LEXING: IMPORTANT, DO NOT SHARE
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::dv(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{
  CCSubMatLib& submat = *(CCSubMatLib::getlibptr());

  int numf = gm.numf();
  int numv = gm.numv();
  int size = pow2(lvl)+1;
  double h = 1.0/pow2(lvl);
  iA(fm.size()==numf);
  to.resize(numf);
  for(int F=0; F<numf; F++) {
	 iA(fm[F].m()==size && fm[F].n()==size);
	 to[F].resize(size,size); //iA(to[F].size()==size);
  }
  for(int F=0; F<numf; F++) {
	 CCRect<Point3>& CD = fm[F];
	 CCRect<Point3>& ND = to[F];
	 
	 for(int i=0; i<size; i++)
		for(int j=0; j<size; j++) {
		  ND(i,j) = 
			 ( (CD(i,  j+1)-CD(i,  j-1)  )*4.0 + 
				(CD(i+1,j+1)-CD(i+1,j-1)) +
				(CD(i-1,j+1)-CD(i-1,j-1)) ) / (12.0*h);
		}
  }
  //LEXING: IMPORTANT, DO NOT SHARE
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::duu(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{
  CCSubMatLib& submat = *(CCSubMatLib::getlibptr());

  int numf = gm.numf();  int numv = gm.numv();
  int size = pow2(lvl)+1;  double h = 1.0/pow2(lvl);
  iA(fm.size()==numf);
  to.resize(numf);
  for(int F=0; F<numf; F++) {
	 iA(fm[F].m()==size && fm[F].n()==size);
	 to[F].resize(size,size); //iA(to[F].size()==size);
  }
  for(int F=0; F<numf; F++) {
	 CCRect<Point3>& CD = fm[F];
	 CCRect<Point3>& ND = to[F];
	 for(int i=0; i<size; i++)
		for(int j=0; j<size; j++) {
		  ND(i,j) = 
			 ( (CD(i+1,j+1)-2.*CD(i,j+1)+CD(i-1,j+1)) +
				(CD(i+1,j  )-2.*CD(i,j  )+CD(i-1,j  ))*4.0 + 
				(CD(i+1,j-1)-2.*CD(i,j-1)+CD(i-1,j-1)) )  / (6.0*h*h);
		}
  }
  //LEXING: IMPORTANT, DO NOT SHARE
  return (0);
}
// ---------------------------------------------------------------------- 
int CCSurfOp::dvv(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{
  CCSubMatLib& submat = *(CCSubMatLib::getlibptr());

  int numf = gm.numf();  int numv = gm.numv();
  int size = pow2(lvl)+1;  double h = 1.0/pow2(lvl);
  iA(fm.size()==numf);
  to.resize(numf);
  for(int F=0; F<numf; F++) {
	 iA(fm[F].m()==size && fm[F].n()==size);
	 to[F].resize(size,size); //iA(to[F].size()==size);
  }
  for(int F=0; F<numf; F++) {
	 CCRect<Point3>& CD = fm[F];
	 CCRect<Point3>& ND = to[F];
	 for(int i=0; i<size; i++)
		for(int j=0; j<size; j++) {
		  ND(i,j) = 
			 ( (CD(i+1,j+1)-2.*CD(i+1,j)+CD(i+1,j-1)) + 
				(CD(i,  j+1)-2.*CD(i  ,j)+CD(i,  j-1)  )*4.0 + 
				(CD(i-1,j+1)-2.*CD(i-1,j)+CD(i-1,j-1)) ) / (6.0*h*h);
		}
  }
  //LEXING: IMPORTANT, DO NOT SHARE
  return (0);
}
// ---------------------------------------------------------------------- 
int CCSurfOp::duv(GpMesh& gm, int lvl, DataLevel& fm, DataLevel& to)
{
  CCSubMatLib& submat = *(CCSubMatLib::getlibptr());

  int numf = gm.numf();  int numv = gm.numv();
  int size = pow2(lvl)+1;  double h = 1.0/pow2(lvl);
  iA(fm.size()==numf);
  to.resize(numf);
  for(int F=0; F<numf; F++) {
	 iA(fm[F].m()==size && fm[F].n()==size);
	 to[F].resize(size,size); //iA(to[F].size()==size);
  }
  for(int F=0; F<numf; F++) {
	 CCRect<Point3>& CD = fm[F];
	 CCRect<Point3>& ND = to[F];
	 for(int i=0; i<size; i++)
		for(int j=0; j<size; j++) {
		  ND(i,j) = 
			 ( (CD(i+1,j+1)-CD(i+1,j-1)) -
				(CD(i-1,j+1)-CD(i-1,j-1)) ) / (4.0*h*h);
		}
  }
  //LEXING: IMPORTANT, DO NOT SHARE
  return (0);
}

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// ---------------------------------------------------------------------- 
int CCSurfOp::eno2index(int lvl, int e, int t, int& i, int& j)
{
  int size = pow2(lvl)+1;
  switch(e) {
  case 0: i = t; j = 0; break;
  case 1: i = size-1; j = t; break;
  case 2: i = size-1-t; j = size-1; break;
  case 3: i = 0; j = size-1-t; break;
  default: iA(0);
  }
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::getcenter(GpMesh& gm, int lvl, DataLevel& wk, int V, Point3& val)
{
  //just get one
  int F, v; gm.Vf2Fv(V, 0, F, v);   //intpair Fv = gm.Vf2Fv(V,0);
  
  int i, j;
  iC( eno2index(lvl, v, 0, i, j) );
  val = wk[F](i,j);
  
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::putcenter(GpMesh& gm, int lvl, DataLevel& wk, int V, Point3& val)
{
  int K = gm.Vvalence(V);
  for(int f=0; f<K; f++) {
	 int F, v;	 gm.Vf2Fv(V, f, F, v);	 //intpair Fv = gm.Vf2Fv(V,f);
	 int i, j;
	 iC( eno2index(lvl, v, 0, i, j) );
	 wk[F](i,j) = val;
  }
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::getonering(GpMesh& gm, int lvl, DataLevel& wk, int V, vector<Point3>& val)
{
  int K = gm.Vvalence(V);
  val.resize(2*K+1);
  
  for(int f=0; f<K; f++) {
	 int F, v;	 gm.Vf2Fv(V, f, F, v);	 //intpair Fv = gm.Vf2Fv(V,f);
	 int i,j;
	 iC( eno2index(lvl, v, 0, i, j) );
	 int di,dj;
	 iC( eno2index(lvl, v, 1, di, dj) );
	 val[0] = wk[F](i,j);
	 val[1  +f] = wk[F](di,dj);
	 val[1+K+f] = wk[F](di-(dj-j),dj+(di-i));
  }
  return (0);
}

// ---------------------------------------------------------------------- 
int CCSurfOp::putonering(GpMesh& gm, int lvl, DataLevel& wk, int V, vector<Point3>& val)
{
  int K = gm.Vvalence(V);
  iA(val.size()==2*K+1);
  for(int f=0; f<K; f++) {
	 int F, v;	 gm.Vf2Fv(V, f, F, v);	 //intpair Fv = gm.Vf2Fv(V,f);
	 int i,j;
	 iC( eno2index(lvl, v, 0, i, j) );
	 int di,dj;
	 iC( eno2index(lvl, v, 1, di, dj) );
	 wk[F](i,j) = val[0];
	 wk[F](di,dj) = val[1  +f];
	 wk[F](di-(dj-j),dj+(di-i)) = val[1+K+f];
  }
  return (0);
}
