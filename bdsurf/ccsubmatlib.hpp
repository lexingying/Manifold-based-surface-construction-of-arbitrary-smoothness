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
#ifndef _CCSUBMATLIB_HPP_
#define _CCSUBMATLIB_HPP_

#include "comobject.hpp"
#include "common/nummat.hpp"

using std::vector;
using std::istream;

class CCSubMatLib: ComObject {
public:
  class Entry {
  public:
	 DblNumMat& SM() { return _SM; }
	 DblNumMat& IV() { return _IV; }
	 DblNumMat& V()  { return _V;  }
	 DblNumVec& D()  { return _D;  }
  protected:
	 DblNumMat _SM;
	 DblNumMat _IV;
	 DblNumMat _V;
	 DblNumVec _D;
  };
  enum {
	 MAXVALENCE = 13
  };
  enum {
	 INTERIOR=1,
	 BOUNDARY=2,
	 CORNER  =3
  };
public:
  CCSubMatLib(const string& p): ComObject(p) {;}
  ~CCSubMatLib() {;}
  //ops
  int setup(istream& sin);
  Entry& chooseMatrix(int flag, int k);
protected:
  vector<Entry> _interior;
  vector<Entry> _boundary;
  vector<Entry> _corner;
  
  //------------------
public:
  static CCSubMatLib* _libptr;
  static CCSubMatLib* getlibptr();
};

#endif
