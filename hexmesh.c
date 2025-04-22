#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "core/meshgen.h"

int make_type1_subchannel(int,int,double,int);

int nquad=0,npts=0,nelem=0,nvert=0,ncide=0,nfld=2;

block blocks[max_block];
quad elems[max_elem];
point verts[4*max_elem];
edge cides[4*max_elem];

point *points;
point origin={0.0, 0.0};

#include "core/spaces.c"
#include "core/triangles.c"
#include "core/ptops.c"
#include "core/conops.c"
#include "core/output.c"
#include "core/read_inp.c"

int main(int argc, char *argv[]){

  char inname[64],reaname[64];
  char bcs[4][2][4];

  point corners[4];
  point translate;

  int i,j,k,n;

  double pio6=M_PI/6., pio3=2.*pio6;
  int a=2,b=3,c=3,d=3,e=3,f=2,g=3;

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
  
  int N_pin_rings = 3;
  double pitch = 1.10,delta=0.1*(pitch-1.);
  double apoth = pitch*((double)N_pin_rings - 0.5)*sin(pio6);
  points=malloc(32*sizeof(point));

//Layout the canonical type 1 subchannel
  translate.x=-0.5;translate.y=0.0;
  point type1[16];
  point base[16];

  type1[0]=origin;
  type1[1].x=pitch/2.;type1[1].y=-0.5*tan(pio6)*pitch;
  type1[2]=rotate_point(type1[1],2.*pio3,type1[0]);
  type1[3]=rotate_point(type1[1],4.*pio3,type1[0]);
  type1[4]=midpoint(type1[1],type1[3]);
  type1[5]=line_circle_intercept(type1[4],type1[1],type1[1],0.5);
  type1[6]=line_circle_intercept(type1[0],type1[1],type1[1],0.5);
  type1[7]=line_circle_intercept(type1[2],type1[1],type1[1],0.5);
  for(int iang=1;iang<3;iang++){
    double theta = 2.*((double)iang)*pio3;
    for(int ipt=0;ipt<4;ipt++) type1[4+ipt+4*iang]=rotate_point(type1[4+ipt],theta,type1[0]);
  }

//sprintf(reaname,"type1.dat");
//output_pts(points,16,reaname);

//int iring=1;
//int nchans=6;

  translate.x=-type1[3].x;translate.y=-type1[3].y;
  for(int ipt=0;ipt<16;ipt++) points[ipt]=translate_point(type1[ipt],translate);

//ring 1
  for(int iang=0;iang<6;iang++){
    for(int ipt=0;ipt<16;ipt++) points[ipt]=rotate_point(points[ipt],pio3,origin);
    make_type1_subchannel(2,3,delta,0);
  }

//ring 2
  int iring=2;
  
  translate.x+=pitch;
  for(int ipt=0;ipt<16;ipt++) points[ipt]=translate_point(type1[ipt],translate);

  for(int iang=0;iang<6;iang++){
    for(int ipt=0;ipt<16;ipt++) points[ipt]=rotate_point(points[ipt],pio3,origin);
    make_type1_subchannel(2,3,delta,0);
  }

  translate.x-=pitch/2.0;
  translate.y+=pitch*sin(pio3);
  for(int ipt=0;ipt<16;ipt++) points[ipt]=translate_point(type1[ipt],translate);

  for(int iang=0;iang<6;iang++){
    for(int ipt=0;ipt<16;ipt++) points[ipt]=rotate_point(points[ipt],pio3,origin);
    make_type1_subchannel(2,3,delta,0);
  }

  translate.y-=pitch/(2.*cos(pio6));
  for(int ipt=0;ipt<16;ipt++) points[ipt]=rotate_point(type1[ipt],M_PI,origin);
  for(int ipt=0;ipt<16;ipt++) points[ipt]=translate_point(points[ipt],translate);

  for(int iang=0;iang<6;iang++){
    for(int ipt=0;ipt<16;ipt++) points[ipt]=rotate_point(points[ipt],pio3,origin);
    make_type1_subchannel(2,3,delta,0);
  }

/*
  for(k=0;k<6;k++){

    for(i=0;i<4;i++) for(j=0;j<2;j++) sprintf(bcs[i][j],"W  ");

    corners[0]=points[4+4*k];
    corners[1]=points[3+4*k];
    corners[2]=points[5+4*k];
    corners[3]=points[6+4*k];
    make_cgquad_space(4,3,-R,delta,0.0,corners,bcs);

    if(k==5){
      corners[0]=points[3];
      corners[3]=points[5];
    }else{
      corners[0]=points[7+4*k];
      corners[3]=points[9+4*k];
    }
    corners[1]=points[4+4*k];
    corners[2]=points[6+4*k];
    make_cgquad_space(4,3,-R,delta,0.0,corners,bcs);
    
  }

  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[3][0],"SYM");
  sprintf(bcs[0][1],"I  ");
  sprintf(bcs[3][1],"I  ");
  make_gquad_space(2*a,b,delta,corners,bcs);

  corners[0]=points[9];
  corners[1]=points[4];
  corners[2]=points[3];
  corners[3]=points[6];
  for(i=0;i<4;i++) for(j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[2][0],"SYM");
  sprintf(bcs[2][1],"I  ");
  make_quad_space(g,2*a,corners,bcs);

  corners[0]=points[2];
  corners[1]=points[7];
  corners[2]=points[9];
  corners[3]=points[1];
  for(i=0;i<4;i++) for(j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[3][0],"W  ");
  sprintf(bcs[0][1],"I  ");
  sprintf(bcs[3][1],"I  ");
  make_cquad_space(b,c,R,delta,delta,corners,bcs);

  corners[0]=points[7];
  corners[1]=points[8];
  corners[2]=points[4];
  corners[3]=points[9];
  for(i=0;i<4;i++) for(j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[0][1],"I  ");
  make_cquad_space(g,c,R,delta,0.0,corners,bcs);

  corners[0]=points[8];
  corners[1]=points[5];
  corners[2]=points[3];
  corners[3]=points[4];
  for(i=0;i<4;i++) for(j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[1][0],"SYM");
  sprintf(bcs[0][1],"I  ");
  sprintf(bcs[1][1],"I  ");
  make_cquad_space(2*a,c,R,delta,0.0,corners,bcs);
*/

//write_rea(reaname);
    
//  sprintf(reaname,"pts.dat");
//  output_pts(points,npts,reaname);
  sprintf(reaname,"vts.dat");
  if(nvert>0) output_pts(verts,nvert,reaname);
  
return 0;
}

int make_type1_subchannel(int nx,int ny,double delta,int pid0){

  char bcs[4][2][4];
  point corners[4];

  point *shift=&points[pid0];

  for(int i=0;i<4;i++) for(int j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[0][1],"f  ");
  
  for(int ibox=0;ibox<3;ibox++){
    corners[0]=shift[5+ibox*4];
    corners[1]=shift[6+ibox*4];
    corners[2]=shift[0];
    corners[3]=shift[4+ibox*4];
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);
  }

  for(int ibox=0;ibox<3;ibox++){
    corners[3]=shift[0];
    corners[0]=shift[6+ibox*4];
    corners[1]=shift[7+ibox*4];
    if(ibox==2) {corners[2]=shift[4];}
    else {corners[2]=shift[8+ibox*4];}
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);
  }

return 0;
}

