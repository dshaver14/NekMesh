#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "core/meshgen.h"

int make_type1_subchannel(int,int,double,int);
int make_type2_subchannel(int,int,int,double,double,int);
int make_type3_subchannel(int,int,int,double,double,int,int);
int set_E_bcs(boundary_condition*);

int nquad=0,npts=0,nelem=0,nvert=0,ncide=0,nfld=1;

//block blocks[max_block];
quad elems[max_elem];
point verts[4*max_elem];
edge cides[4*max_elem];

point *points;
static point origin={0.0, 0.0, 0.0};

#include "core/spaces.c"
#include "core/triangles.c"
#include "core/ptops.c"
#include "core/conops.c"
#include "core/output.c"
#include "core/read_inp.c"

int set_E_bcs(boundary_condition *BC){

  for(int i=0;i<4;i++){
    BC->ID[i]=0;
    for(int j=0;j<2;j++) sprintf(BC->cbc[i][j],"E  ");
  }

return 0; 
}

int main(int argc, char *argv[]){

  char inname[64],reaname[64];

  point corners[4];
  point base[4];
  point translate;
  boundary_condition BC;

  int i,j,k,n;

  static double pio6=M_PI/6., pio3=M_PI/3.;
  translate.z=0.0;

  for(i=1;i<argc;i++){
    if(strncmp(argv[i],"-fi",3)==0){
      i++;
      strcpy(inname,argv[i]);
    }else if(strncmp(argv[i],"-rea",4)==0){
      i++;
      strcpy(reaname,argv[i]);
    }
  }

//read the points from the input file
//read_inp(inname);
  
  int N_pin_rings = 2;
  double pitch = 1.15676;
  double deltag=0.10*(pitch-1.);
  double delta=deltag*1.5;
  int nx=6,ny=4;
//double apoth =        ((double)N_pin_rings - 1.0)*pitch*cos(pio6)+0.5+0.12;
//double gap = apoth - (((double)N_pin_rings - 1.0)*pitch*cos(pio6)+0.5);
  double apoth = pitch/2.0;
  points=malloc(4*sizeof(point));


  base[0].x=0.5*cos(pio3);
  base[0].y=0.5*sin(pio3);
  base[1].x=0.5;
  base[1].y=0.0;
  base[2].x=apoth/cos(pio6);
  base[2].y=0.0;
  base[3]=rotate_point(base[2],pio3,origin);

  set_E_bcs(&BC);

  sprintf(BC.cbc[0][0],"pin");
  sprintf(BC.cbc[0][1],"pin");

  sprintf(BC.cbc[2][0],"can");
  sprintf(BC.cbc[2][1],"can");

  for(int isd=0;isd<6;isd++){
    double theta=(double)isd*pio3;

    for(int i=0;i<4;i++) corners[i]=rotate_point(base[i],theta,origin);
    
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,BC.cbc);

  }

//sprintf(reaname,"type3.dat");
//output_pts(points,npts3,reaname);

  write_rea(reaname);

  sprintf(reaname,"onepin.vtk");
  write_vtk(reaname);
    
//  sprintf(reaname,"pts.dat");
//  output_pts(points,npts,reaname);
  sprintf(reaname,"vts.dat");
  if(nvert>0) output_pts(verts,nvert,reaname);
  
return 0;
}

