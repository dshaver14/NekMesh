#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "core/meshgen.h"

int make_type1_subchannel(int,int,double,int);
int make_type2_subchannel(int,int,int,double,double,int);
int make_type3_subchannel(int,int,int,double,double,int);

int nquad=0,npts=0,nelem=0,nvert=0,ncide=0,nfld=2;

block blocks[max_block];
quad elems[max_elem];
point verts[4*max_elem];
edge cides[4*max_elem];

point *points;
static point origin={0.0, 0.0};

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

  static double pio6=M_PI/6., pio3=M_PI/3.;

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
  double pitch = 1.10;
  double delta=0.05*(pitch-1.);
  double deltag=delta*2.5;
  int nx=4,ny=4,ng=5;
  double apoth =        ((double)N_pin_rings - 1.0)*pitch*cos(pio6)+0.5+0.01;
  double gap = apoth - (((double)N_pin_rings - 1.0)*pitch*cos(pio6)+0.5);
  double gap_min = ((double)ny*delta+(double)ng*deltag)*1.01;
  if(gap<gap_min) {
    printf("Warning: gap between pin and wall too small for requested BL spacing. %f > %f\n",gap,gap_min);
    gap=gap_min;
  }
  points=malloc(16*sizeof(point));

//Layout the canonical type 1 subchannel (center)
  int npts1=16;
  point type1[npts1];
  
  type1[0]=origin;
  type1[1].x=0.5*pitch;type1[1].y=-0.5*tan(pio6)*pitch;
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
  sprintf(inname,"type1.dat");
  output_pts(type1,npts1,inname);

//Layout the canonical type 2 subchannel (edge)
//incorporate apothem and gap size later 2025-04-22
  int npts2=16;
  point type2[npts2];

  type2[0]=origin;
  type2[1].x=0.0;type2[1].y=-pitch/(2.*cos(pio6));
  type2[2]=rotate_point(type2[1],2.*pio3,type2[0]);
  type2[3]=rotate_point(type2[1],4.*pio3,type2[0]);
  type2[4]=midpoint(type2[1],type2[3]);
  type2[5]=line_circle_intercept(type2[3],type2[1],type2[1],0.5);
  type2[6]=line_circle_intercept(type2[0],type2[1],type2[1],0.5);
  type2[7]=line_circle_intercept(type2[2],type2[1],type2[1],0.5);
  translate.x=(gap+0.5)*cos(pio6);translate.y=(gap+0.5)*sin(pio6);
  for(int iang=1;iang<3;iang++){
    double theta = 2.*((double)iang)*pio3;
    for(int ipt=0;ipt<4;ipt++) type2[4+ipt+4*iang]=rotate_point(type2[4+ipt],theta,type2[0]);
  }
  //adjust 9, 10, 11 for gap size
  double s=0.5+gap-(pitch-0.5)*cos(pio6);
  translate.x=s*cos(pio6);translate.y=s*sin(pio6);
  type2[9]=translate_point(type2[9],translate);
  type2[11]=translate_point(type2[11],translate);
  type2[10]=midpoint(type2[11],type2[9]);
/*
  type2[8]=line_circle_intercept(type2[1],type2[2],type2[2],0.5);
  type2[10]=line_circle_intercept(type2[3],type2[2],type2[2],0.5);
  type2[9]=midpoint(type2[8],type2[10]);
  type2[11]=line_circle_intercept(type2[2],type2[3],type2[3],0.5);
  type2[12]=line_circle_intercept(type2[0],type2[3],type2[3],0.5);
  type2[13]=line_circle_intercept(type2[1],type2[3],type2[3],0.5);
  //adjust 8, 9, 10 for gap size
  type2[9]=rotate_point(type2[5],2.*pio3,type2[0]);
  type2[11]=rotate_point(type2[7],2.*pio3,type2[0]);
  type2[10]=midpoint(type2[9],type2[11]);
  type2[12]=midpoint(type2[2],type2[3]);
  type2[13]=rotate_point(type2[5],4.*pio3,type2[0]);
  type2[14]=rotate_point(type2[6],4.*pio3,type2[0]);
  type2[15]=rotate_point(type2[7],4.*pio3,type2[0]);
  //adjust 9, 10, 11 for gap size
*/
  sprintf(inname,"type2.dat");
  output_pts(type2,npts2,inname);

//Layout the canonical type 3 subchannel (corner)
//incorporate apothem and gap size later 2025-04-22
//int npts3=8;
  int npts3=15;
  point type3[npts3];

  type3[0]=origin;
  type3[1].x=-0.5*pitch;type3[1].y=-0.5*pitch*tan(pio6);
  type3[2]=rotate_point(type3[1],2.*pio3,type3[0]);
  type3[3]=rotate_point(type3[1],4.*pio3,type3[0]);
  type3[4]=midpoint(type3[1],type3[3]);
  type3[5]=line_circle_intercept(type3[1],type3[3],type3[1],0.5);
  type3[6]=line_circle_intercept(type3[1],type3[2],type3[1],0.5);
  type3[7]=midpoint(type3[1],type3[2]);
  //adjust 8, 9 for gap size
  type3[8]=type3[1];
  type3[8].x+=(0.5+gap)/cos(pio6);//type3[8].y=0.0;
  type3[9]=rotate_point(type3[8],pio3,type3[1]);
  type3[10]=rotate_point(type3[6],pio6,type3[1]);
  type3[11]=line_line_intercept(type3[8],type3[9],type3[10],type3[0]);
  type3[12]=rotate_point(type3[7],pio6,type3[1]);
  type3[13]=line_circle_intercept(type3[1],type3[2],type3[2],0.5);
  type3[13]=translate_point(type3[13],translate);
  type3[14]=reflect_point(type3[13],type3[1],type3[0]);

/*
  type3[5]=line_circle_intercept(type3[0],type3[1],type3[1],0.5);
  type3[6]=line_circle_intercept(type3[2],type3[1],type3[1],0.5);
  type3[7]=midpoint(type3[1],type3[2]);
  type3[11]=midpoint(type3[1],type3[3]);
  //adjust 8, 10, 12 for gap size
  type3[8]=line_circle_intercept(type3[1],type3[2],type3[2],0.5);
  type3[10]=line_circle_intercept(type3[1],type3[3],type3[3],0.5);
  type3[9]=midpoint(type3[8],type3[10]);
  type3[12]=midpoint(type3[5],type3[9]);
*/
  sprintf(inname,"type3.dat");
  output_pts(type3,npts3,inname);

 
//make the mesh 
  for(int ipring=1;ipring<=N_pin_rings;ipring++){
    //regular channels and corner channels
    translate.x=((double)ipring-0.5)*pitch; translate.y=0.5*pitch*tan(pio6);
    for(int ichan=0;ichan<ipring;ichan++){
      if(ipring<N_pin_rings){
        for(int ipt=0;ipt<npts1;ipt++) points[ipt]=translate_point(type1[ipt],translate);

        for(int iang=0;iang<6;iang++){
          for(int ipt=0;ipt<npts1;ipt++) points[ipt]=rotate_point(points[ipt],pio3,origin);
          make_type1_subchannel(nx,ny,delta,0);
        }
      }else{
        for(int ipt=0;ipt<npts3;ipt++) points[ipt]=translate_point(type3[ipt],translate);

        for(int iang=0;iang<6;iang++){
          for(int ipt=0;ipt<npts3;ipt++) points[ipt]=rotate_point(points[ipt],pio3,origin);
          make_type3_subchannel(nx,ny,ng,delta,deltag,0);
        }
      }
      translate.x-=0.5*pitch;
      translate.y+=pitch*sin(pio3);
    }
    //flipped channels and edge channels
    translate.x=((double)ipring-1.0)*pitch; translate.y=pitch*tan(pio6);
    for(int ichan=0;ichan<(ipring-1);ichan++){
      if(ipring<N_pin_rings){
        for(int ipt=0;ipt<npts1;ipt++) points[ipt]=rotate_point(type1[ipt],M_PI,origin);
        for(int ipt=0;ipt<npts1;ipt++) points[ipt]=translate_point(points[ipt],translate);

        for(int iang=0;iang<6;iang++){
          for(int ipt=0;ipt<npts1;ipt++) points[ipt]=rotate_point(points[ipt],pio3,origin);
          make_type1_subchannel(nx,ny,delta,0);
        }
      }else{
        for(int ipt=0;ipt<npts2;ipt++) points[ipt]=translate_point(type2[ipt],translate);

        for(int iang=0;iang<6;iang++){
          for(int ipt=0;ipt<npts2;ipt++) points[ipt]=rotate_point(points[ipt],pio3,origin);
          make_type2_subchannel(nx,ny,ng,delta,deltag,0);
        }
      }
      translate.x-=0.5*pitch;
      translate.y+=pitch*sin(pio3);
    }
  }

//sprintf(reaname,"type3.dat");
//output_pts(points,npts3,reaname);

  write_rea(reaname);
    
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

int make_type2_subchannel(int nx,int ny,int ng,double delta,double deltag,int pid0){

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
    if(ibox==1) make_gquad_space(nx,ng,deltag,corners,bcs);
    else make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);
  }

  for(int ibox=0;ibox<3;ibox++){
    corners[3]=shift[0];
    corners[0]=shift[6+ibox*4];
    corners[1]=shift[7+ibox*4];
    if(ibox==2) {corners[2]=shift[4];}
    else {corners[2]=shift[8+ibox*4];}
    if(ibox==1) make_gquad_space(nx,ng,deltag,corners,bcs);
    else make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);
  }


/*
  if(nx%2==1){

    corners[0]=shift[5];
    corners[1]=shift[6];
    corners[2]=shift[9];
    corners[3]=shift[4];
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

    corners[0]=shift[12];
    corners[1]=shift[13];
    corners[2]=shift[4];
    corners[3]=shift[9];
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

    sprintf(bcs[2][0],"W  ");
    sprintf(bcs[2][1],"f  ");

    corners[0]=shift[6];
    corners[1]=shift[7];
    corners[2]=shift[8];
    corners[3]=shift[9];
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

    corners[0]=shift[11];
    corners[1]=shift[12];
    corners[2]=shift[9];
    corners[3]=shift[10];
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

  }else{

    corners[0]=shift[5];
    corners[1]=shift[6];
    corners[2]=shift[0];
    corners[3]=shift[4];
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

    corners[0]=shift[12];
    corners[1]=shift[13];
    corners[2]=shift[4];
    corners[3]=shift[0];
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

    corners[0]=shift[6];
    corners[1]=shift[7];
    corners[2]=shift[8];
    corners[3]=shift[0];
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

    corners[0]=shift[11];
    corners[1]=shift[12];
    corners[2]=shift[0];
    corners[3]=shift[10];
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

    corners[0]=shift[8];
    corners[1]=shift[10];
    corners[2]=shift[0];
    make_tri_space(nx/2,corners,bcs);

    corners[0]=shift[9];
    corners[1]=shift[10];
    corners[2]=shift[0];
    make_tri_space(nx/2,corners,bcs);
   }
*/

return 0;
}

int make_type3_subchannel(int nx,int ny,int ng,double delta,double deltag,int pid0){

  char bcs[4][2][4];
  point corners[4];

  point *shift=&points[pid0];

  for(int i=0;i<4;i++) for(int j=0;j<2;j++) sprintf(bcs[i][j],"E  ");
  sprintf(bcs[0][0],"W  ");
  sprintf(bcs[0][1],"f  ");
//sprintf(bcs[2][0],"W  ");
//sprintf(bcs[2][1],"f  ");

  corners[0]=shift[5];
  corners[1]=shift[10];
  corners[2]=shift[12];
  corners[3]=shift[4];
  make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

  corners[0]=shift[10];
  corners[1]=shift[6];
  corners[2]=shift[7];
  corners[3]=shift[12];
  make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

  corners[0]=shift[8];
  corners[1]=shift[9];
  corners[2]=shift[4];
  corners[3]=shift[7];
//make_gquad_space(2*nx,ng,deltag,corners,bcs);
  
/*
  corners[0]=shift[5];
  corners[1]=shift[6];
  corners[2]=shift[7];
  corners[3]=shift[12];
  make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,bcs);

  corners[0]=shift[8];
  corners[1]=shift[9];
  corners[2]=shift[12];
  corners[3]=shift[7];
  make_gquad_space(nx,ny,delta,corners,bcs);

  corners[0]=shift[9];
  corners[1]=shift[10];
  corners[2]=shift[11];
  corners[3]=shift[12];
  make_gquad_space(nx,ny,delta,corners,bcs);
*/
return 0;
}
  
