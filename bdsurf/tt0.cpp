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
#include "ccsurfvr.hpp"

using namespace std;

int optionsCreate(const char* optfile, map<string,string>& options)
{
  options.clear();
  ifstream fin(optfile);
  if(fin.good()==false) {
	 cerr<<"wrong option file"<<endl;	 exit(1);
  }
  string name;  fin>>name;
  while(fin.good()) {
	 char cont[100];	 fin.getline(cont, 99);
	 options[name] = string(cont);
	 fin>>name;
  }
  fin.close();
  return 0;
}

int main(int argc, char** argv)
{
  srand48( (long)time(NULL) );
  iA(argc==2);
  map<string,string> opts;
  optionsCreate(argv[1], opts);

  //1. ccsurf
  CCSurf ccsurf("ccsurf_");
  iC( ccsurf.setup(opts) );
  
  map<string,string>::iterator mi;
  mi = opts.find("-lvl"); iA(mi!=opts.end());
  int lvl;  { istringstream ss((*mi).second);  ss>>lvl; }
  for(int i=0; i<lvl; i++) {	 iC( ccsurf.subdivide(i) );  }
  
  CCSurfVr ccsurfvr(&argc, argv, &ccsurf, lvl);
  glutMainLoop();
  
  return 0;
}
