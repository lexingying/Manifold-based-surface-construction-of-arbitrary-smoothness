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
#ifndef _GPMESH_HPP_
#define _GPMESH_HPP_

#include "common/vec3t.hpp"
#include "comobject.hpp"

using std::vector;
using std::pair;

//----------------------
//quadrilateral mesh
//----------------------
class GpMesh: public ComObject {
public:
  typedef pair<int,int> intpair;
  //-----------
  //Topology related data types
  struct Edge {
	 int SV; //starting vert
	 int CE; //conjugate(neighbor) edge
	 int F;  //face
	 int PE; //prev edge
	 int NE; //next edge	 //int GI; //group id
  };
  struct Vert {
	 int SE; //starting edge	 //int GI; //group id
  };
  struct Face {
	 int SE; //starting edge	 //int GI; //group id
  };
  //-----------
  /*
  //Group data structure
  class Group { //Group == Object
  protected:
	 Point3 _intpt; //interior point
	 int _orient; //orientation
	 int _vcnt, _ecnt, _fcnt;
  public:
	 Point3& intpt() { return _intpt; }
	 int& orient() { return _orient; }
	 int& vcnt() { return _vcnt; }
	 int& ecnt() { return _ecnt; }
	 int& fcnt() { return _fcnt; }
  };
  */
protected:
  //  mesh
  vector<Edge> _edges;
  vector<Vert> _verts;
  vector<Face> _faces;  //vector<Group> _groups; //group information vector
  // points
  vector<Point3> _points;
public:
  GpMesh(const string& p): ComObject(p) {;}
  ~GpMesh() {;}
  int setup(map<string,string>& opts);  //int dump( ostream& ot);
  Point3 ctr();
  void bbx(Point3& bbmin, Point3& bbmax);
  
  vector<Edge>& edges() { return _edges; }
  vector<Vert>& verts() { return _verts; }
  vector<Face>& faces() { return _faces; }  //vector<Group>& groups() { return _groups; }
  vector<Point3>& points() { return _points; }
  
  int numv() { return _verts.size(); }
  int nume() { return _edges.size(); }
  int numf() { return _faces.size(); }  //int numgrp() { return _groups.size(); }
    //  edge function
  bool Emajor(int E) {	 int CE = _edges[E].CE;	 return (E > CE);  }
  bool Einterior(int E) {	 int CE = _edges[E].CE;	 return (CE!=-1);  }
  int Econj(int E) {	 return _edges[E].CE;  }
  void E2Vs(int E, vector<int>& Vs) { //from E to get all Vs of E's face
	 Vs.resize(4);
	 int ME = E;          Vs[0] = _edges[ME].SV;
	 ME = _edges[ME].NE;  Vs[1] = _edges[ME].SV;
	 ME = _edges[ME].NE;  Vs[2] = _edges[ME].SV;
	 ME = _edges[ME].NE;  Vs[3] = _edges[ME].SV;
  }
  //  vert function
  bool Vclosed(int V) {	 int SE = _verts[V].SE;	 int CE = _edges[SE].CE;	 return (CE!=-1); } //for open ones, SE is the first one
  int Vvalence(int V) {
	 int cnt = 0;	 int SE = _verts[V].SE;	 int ME = SE;	 bool start=true;
	 while( (start==true || ME!=SE) && _edges[_edges[ME].PE].CE!=-1) {
		start = false;
		ME = _edges[_edges[ME].PE].CE;
		cnt ++;
	 }
	 return cnt;
  }
  //  face function
  void F2Vs(int F, vector<int>& Vs) {	 int E = _faces[F].SE;	 E2Vs(E, Vs);  }
  void F2Es(int F, vector<int>& Es) {
	 Es.resize(4);  //assert(Es.size()==3);
	 int ME = _faces[F].SE;  Es[0] = ME;
	 ME = _edges[ME].NE;     Es[1] = ME;
	 ME = _edges[ME].NE;     Es[2] = ME;
	 ME = _edges[ME].NE;     Es[3] = ME;
  }
  //others
  void Fe2Fe(int F, int e, int& cF, int& ce) {
	 vector<int> Es(4);	 F2Es(F, Es); //get all edges
	 int E = Es[e]; //get eth edge
	 int cE = _edges[E].CE; //find neighbor edge
	 cF = _edges[cE].F; //find neighbor face
	 F2Es(cF, Es); //get all edges
	 for(int i=0; i<4; i++)		if(Es[i]==cE)		  ce = i; //compute ce
	 return;
  }
  void Vf2Fv(int V, int f, int& F, int& v) {
	 int ME = _verts[V].SE;
	 for(int i=0; i<f; i++)
		ME = _edges[_edges[ME].PE].CE;
	 F = _edges[ME].F; //face
	 vector<int> Es(4);	 F2Es(F, Es);
	 for(int i=0; i<4; i++)		if(Es[i]==ME)		  v = i;
	 return;
  }
  void Fv2Vf(int F, int v, int& V, int& f) {
	 vector<int> Vs(4);	 F2Vs(F, Vs);
	 V = Vs[v];
	 int ME = _verts[V].SE;
	 int cnt = 0;
	 while( _edges[ME].F != F) {
		ME = _edges[_edges[ME].PE].CE;
		cnt ++;
	 }
	 f = cnt;
	 return;
  }
  //  aux functions
  int nextvno(int vno) { return (vno+1)%4; }
  int prevvno(int vno) { return (vno+3)%4; }
  int nexteno(int eno) { return (eno+1)%4; }
  int preveno(int eno) { return (eno+3)%4; }
  int fromvno(int eno) { return eno; }
  int gotovno(int eno) { return (eno+1)%4; }

private:
  int init_readmesh(istream&, vector<Point3>&, vector< vector<int> >&);
  int init_createmesh(        vector<Point3>&, vector< vector<int> >&);
};



#endif
