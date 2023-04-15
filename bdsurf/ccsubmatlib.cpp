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
#include "ccsubmatlib.hpp"

using std::ifstream;

CCSubMatLib* CCSubMatLib::_libptr = NULL;

//-----------------------------------------------------------------------
CCSubMatLib* CCSubMatLib::getlibptr()
{
  if(_libptr==NULL) {
	 _libptr = new CCSubMatLib("ccsubmatlib_");	 //string data;	 //file2string(PETSC_COMM_WORLD, "submat.dat", data);
	 ifstream tmp("ccsubmat.dat");	 //istringstream tmp(data);
	 _libptr->setup(tmp);
  }
  return _libptr;
}

// ---------------------------------------------------------------------- 
int CCSubMatLib::setup(istream& fin)
{
  iA(fin.good());
  _interior.resize(MAXVALENCE);
  _boundary.resize(MAXVALENCE);
  _corner.resize(MAXVALENCE);
  
  char tmp[100];
  fin>>tmp;
  while(fin.eof()==false) {
	 iA(strcmp(tmp, "interior")==0);
	 vector<Entry>& cur = _interior;
	 //get valence
	 int k;
	 fin>>k;	 iA(k>=3 && k<MAXVALENCE);
	 Entry& entry = cur[k];
	 //get size
	 int size;
	 fin>>size;	 iA(size == 6*k+1);
	 //read in SM
	 DblNumMat& sm = entry.SM();
	 sm.resize(size,size);
	 for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
		  fin>>sm(i,j); //read row by row
	 //read in IV
	 DblNumMat& iv = entry.IV();
	 iv.resize(size,size);
	 for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
		  fin>>iv(i,j); //read row by row
	 //read in V
	 DblNumMat& v = entry.V();
	 v.resize(size,size);
	 for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
		  fin>>v(i,j);  //read row by row
	 //read in D
	 DblNumVec& d = entry.D();
	 d.resize(size);
	 for(int i=0; i<size; i++)
		fin>>d(i);
	 
	 fin>>tmp;
  }
  return (0);
}


CCSubMatLib::Entry& CCSubMatLib::chooseMatrix(int flag, int k)
{
  if(flag==INTERIOR) {
	 iA(k<_interior.size());
	 return _interior[k];
  } else if(flag==BOUNDARY) {
	 iA(k<_boundary.size());
	 return _boundary[k];
  } else { //(flag==CORNER)
	 iA(k<_corner.size());
	 return _corner[k];
  }
}

