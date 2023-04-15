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
#include "gpmesh.hpp"

using std::map;
using std::ifstream;
using std::istringstream;

// ---------------------------------------------------------------------- 
int GpMesh::setup(map<string,string>& opts)
{
  char meshfile[100];
  map<string,string>::iterator mi;
  mi = opts.find("-" + prefix() + "meshfile"); iA(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>meshfile; }
  ifstream fin(meshfile);

  vector<Point3> points;
  vector< vector<int> > coords;  //vector<int> gids;
  init_readmesh(fin, points, coords); //group info written there
  init_createmesh(points, coords);
  
  return (0);
}

// ---------------------------------------------------------------------- 
int GpMesh::init_readmesh(istream& fin, vector<Point3>& points, vector< vector<int> >& coords)
{
  points.clear();
  coords.clear();
  //gids.clear();
  string tmp;
  int cnt;
  fin>>tmp; iA(tmp=="points");
  fin>>cnt;
  for(int i=0; i<cnt; i++) {
	 double x,y,z;
	 fin>>x>>y>>z;
	 points.push_back( Point3(x,y,z) );
  }
  fin>>tmp; iA(tmp=="coords");
  fin>>cnt;
  for(int i=0; i<cnt; i++) {
	 vector<int> is(4); //4 vertices 
	 fin>>is[0]>>is[1]>>is[2]>>is[3];
	 coords.push_back(is);
  }
  /*
  fin>>tmp; iA(tmp=="gids");
  fin>>cnt;
  for(int i=0; i<cnt; i++) {
	 int g; fin>>g;
	 gids.push_back(g);
  }
  iA(gids.size()==coords.size()); //both==number of faces
  fin>>tmp; iA(tmp=="groups");
  fin>>cnt; _groups.resize(cnt);
  double x,y,z; 
  for(int i=0; i<cnt; i++) {
	 fin>>tmp; iA(tmp=="intpt");
	 fin>>x>>y>>z; _groups[i].intpt()=Point3(x,y,z); //interior point
	 fin>>tmp; iA(tmp=="orient");
	 fin>>_groups[i].orient(); //orientation
  }
  */
  return (0);
}

// ---------------------------------------------------------------------- 
int GpMesh::init_createmesh(vector<Point3>& points, vector< vector<int> >& coords)
{
  //check
  for(int F=0; F<coords.size(); F++)	 iA(coords[F].size()==4);
  //nv nf ne
  int nv = points.size();
  int nf = coords.size();
  int ne=0;
  for(int F=0; F<coords.size(); F++) 	 ne += coords[F].size();
  _edges.resize(ne);  _verts.resize(nv);  _faces.resize(nf);
  //_faces, the easiest
  for(int F=0; F<nf; F++) {
	 _faces[F].SE = 4*F;	 //_faces[F].GI = gids[F]; //group id
  }
  //_edges,
  map<intpair,int> VV2E;
  int ec=0; //edge count
  for(int F=0; F<nf; F++) {
	 for(int e=0; e<4; e++) {
		int fm = coords[F][fromvno(e)];
		int to = coords[F][gotovno(e)];
		_edges[ec+e].SV = fm;
		_edges[ec+e].CE = -1; //tmp
		_edges[ec+e].F = F;
		_edges[ec+e].PE = ec+preveno(e);
		_edges[ec+e].NE = ec+nexteno(e);		//_edges[ec+e].GI = gids[F]; //group id
		VV2E[intpair(fm,to)] = ec+e;
	 }
	 ec += 4;
  }
  iA(ec==ne);
  ec = 0;
  for(int F=0; F<nf; F++) {
	 for(int e=0; e<4; e++) {
		int fm = coords[F][fromvno(e)];
		int to = coords[F][gotovno(e)];
		map<intpair,int>::iterator mi = VV2E.find(intpair(to,fm));
		if(mi!=VV2E.end())
		  _edges[ec+e].CE = (*mi).second;
	 }
	 ec += 4;
  }
  iA(ec==ne);
  //_verts,
  //for(int F=0; F<nf; F++)	 for(int v=0; v<4; v++)		_verts[ coords[F][v] ].GI = gids[F]; //group id
  for(int V=0; V<nv; V++)
	 _verts[V].SE = -1; //set invalid
  for(int E=0; E<ne; E++) {
	 int V = _edges[E].SV; //start v
	 if(_verts[V].SE==-1) { //do only once
		bool start = true;
		int ME = E;
		while( _edges[ME].CE !=-1 && (start==true||ME!=E)  ) {
		  start = false;
		  ME = _edges[_edges[ME].CE].NE;
		}
		_verts[V].SE = ME;
	 }
  }
  /*
  //_groups, finishing the job
  for(int i=0; i<_groups.size(); i++) {
	 _groups[i].vcnt()=0;
	 _groups[i].ecnt()=0;
	 _groups[i].fcnt()=0;
  }
  for(int i=0; i<_verts.size(); i++)	 _groups[ _verts[i].GI ].vcnt()++;
  for(int i=0; i<_edges.size(); i++)	 _groups[ _edges[i].GI ].ecnt()++;
  for(int i=0; i<_faces.size(); i++)	 _groups[ _faces[i].GI ].fcnt()++;
  */
  //points
  _points = points;
  
  return (0);
}


// ---------------------------------------------------------------------- 
Point3 GpMesh::ctr()
{
  Point3 sum(0,0,0);
  for(int V=0; V<_points.size(); V++)
	 sum += _points[V];
  return sum/double(_points.size());
}

// ---------------------------------------------------------------------- 
void GpMesh::bbx(Point3& bbmin, Point3& bbmax)
{
  bbmin = Point3( SCL_MAX);
  bbmax = Point3(-SCL_MAX);
  for(int V=0; V<_points.size(); V++) {
	 bbmin = min(bbmin, _points[V]);
	 bbmax = max(bbmax, _points[V]);
  }
}



