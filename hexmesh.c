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
static point origin={0.0, 0.0};

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
  point translate;

  int i,j,k,n;

  static double pio6=M_PI/6., pio3=M_PI/3.;
  translate.z=0.0;

  for(int ivrt=0;ivrt<4*max_elem;ivrt++){
    verts[ivrt].x=0.0;
    verts[ivrt].y=0.0;
    verts[ivrt].z=0.0;
  }

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
  double pitch = 1.15676;
  double deltag=0.10*(pitch-1.);
  double delta=deltag*1.5;
  int nx=2,ny=3,ng=3;
//double apoth =        ((double)N_pin_rings - 1.0)*pitch*cos(pio6)+0.5+0.12;
//double gap = apoth - (((double)N_pin_rings - 1.0)*pitch*cos(pio6)+0.5);
  double apoth = 7.32432/2.0;
  double gap = apoth - (((double)4 - 1.0)*pitch*cos(pio6)+0.5);
  double gap_min = ((double)ny*delta+(double)ng*deltag)*1.01;
  if(gap<gap_min) {
    printf("Warning: gap between pin and hexcan wall too small for requested BL spacing. %f < %f\n",gap,gap_min);
    gap=gap_min;
  } else if(gap<(pitch-1.0)) printf("Warning: gap between pin and hexcan wall less than gap between pins! %f < %f \n",gap,pitch-1.0);

  points=malloc(16*sizeof(point));

//Layout the canonical type 1 subchannel (center)
  int npts1=16;
  point type1[npts1];
  for(int ipt=0;ipt<npts1;ipt++) type1[ipt].z=0.0;
  
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
  int npts2=16;
  point type2[npts2];
  for(int ipt=0;ipt<npts2;ipt++) type2[ipt].z=0.0;

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

  sprintf(inname,"type2.dat");
  output_pts(type2,npts2,inname);

//Layout the canonical type 3 subchannel (corner)
//incorporate apothem and gap size later 2025-04-22
  int npts3=15;
  point type3[npts3];
  for(int ipt=0;ipt<npts3;ipt++) type3[ipt].z=0.0;

  type3[0]=origin;
  type3[1].x=-0.5*pitch;type3[1].y=-0.5*pitch*tan(pio6);
  type3[2]=rotate_point(type3[1],2.*pio3,type3[0]);
  type3[3]=rotate_point(type3[1],4.*pio3,type3[0]);
  type3[4]=midpoint(type3[1],type3[3]);
  type3[5]=line_circle_intercept(type3[3],type3[1],type3[1],0.5);
  type3[6]=line_circle_intercept(type3[2],type3[1],type3[1],0.5);
  type3[7]=midpoint(type3[1],type3[2]);
  //adjust 8, 9 for gap size
  type3[8]=type3[1];
  type3[8].x+=(0.5+gap)/cos(pio6);//type3[8].y=0.0;
  type3[9]=rotate_point(type3[8],pio3,type3[1]);
  type3[10]=rotate_point(type3[6],pio6,type3[1]);
  type3[11]=line_line_intercept(type3[8],type3[9],type3[10],type3[0]);
  type3[12]=rotate_point(type3[7],pio6,type3[1]);
//type3[12]=midpoint(type3[10],type3[11]);
  type3[13]=line_circle_intercept(type3[1],type3[2],type3[2],0.5);
  type3[13]=translate_point(type3[13],translate);
  type3[14]=reflect_point(type3[13],type3[1],type3[0]);

  double del1=distance(type3[11],type3[10]);
  double del2=distance(type3[12],type3[10]);
  
  if((del1-del2)<0.0){
    printf("Error: hexcan too tight for adopted meshing strategy!\n");
    return 0;
  }else if((del1-del2)<(double)ng*deltag){
    printf("Error: gap size too small for requested hexcan BL spacing!\n");
    return 0;
  }

  sprintf(inname,"type3.dat");
  output_pts(type3,npts3,inname);
 
//make the mesh 
  for(int ipring=1;ipring<=N_pin_rings;ipring++){
    //regular channels and corner channels
    translate.x=((double)ipring-0.5)*pitch; translate.y=0.5*pitch*tan(pio6);
    for(int ichan=0;ichan<ipring;ichan++){
      int t3s=1;
      if(N_pin_rings>1){
        t3s=3;
        if(ichan==0) t3s=2;
        if(ichan==(ipring-1)) t3s=4;
      }
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
          make_type3_subchannel(nx,ny,ng,delta,deltag,t3s,0);
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

  boundary_condition BC;
  point corners[4];

  point *shift=&points[pid0];

  set_E_bcs(&BC);
  BC.ID[0]=1;
  sprintf(BC.cbc[0][0],"pin");
  sprintf(BC.cbc[0][1],"pin");
  
  for(int ibox=0;ibox<3;ibox++){
    corners[0]=shift[5+ibox*4];
    corners[1]=shift[6+ibox*4];
    corners[2]=shift[0];
    corners[3]=shift[4+ibox*4];
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,BC.cbc);
  }

  for(int ibox=0;ibox<3;ibox++){
    corners[3]=shift[0];
    corners[0]=shift[6+ibox*4];
    corners[1]=shift[7+ibox*4];
    if(ibox==2) {corners[2]=shift[4];}
    else {corners[2]=shift[8+ibox*4];}
    make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,BC.cbc);
  }

return 0;
}

int make_type2_subchannel(int nx,int ny,int ng,double delta,double deltag,int pid0){

  point corners[4];
  boundary_condition BC;

  point *shift=&points[pid0];

  set_E_bcs(&BC);

  for(int ibox=0;ibox<3;ibox++){
    BC.ID[0]=1;
    sprintf(BC.cbc[0][0],"pin");
    sprintf(BC.cbc[0][1],"pin");
    corners[0]=shift[5+ibox*4];
    corners[1]=shift[6+ibox*4];
    corners[2]=shift[0];
    corners[3]=shift[4+ibox*4];
    if(ibox==1){  
      BC.ID[0]=2;
      sprintf(BC.cbc[0][0],"can");
      sprintf(BC.cbc[0][1],"can");
      make_gquad_space(nx,ng,deltag,corners,BC.cbc);
    }
    else make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,BC.cbc);
  }

  for(int ibox=0;ibox<3;ibox++){
    sprintf(BC.cbc[0][0],"pin");
    sprintf(BC.cbc[0][1],"pin");
    corners[3]=shift[0];
    corners[0]=shift[6+ibox*4];
    corners[1]=shift[7+ibox*4];
    if(ibox==2) {corners[2]=shift[4];}
    else {corners[2]=shift[8+ibox*4];}
    if(ibox==1){
      sprintf(BC.cbc[0][0],"can");
      sprintf(BC.cbc[0][1],"can");
      make_gquad_space(nx,ng,deltag,corners,BC.cbc);
    } else make_cgquad_space(nx,ny,-0.5,delta,0.0,corners,BC.cbc);
  }

return 0;
}

int make_type3_subchannel(int nx,int ny,int ng,double delta,double deltag,int subtype,int pid0){

  point corners[4];
  boundary_condition BC;

  point *shift=&points[pid0];

  if(subtype<1||subtype>4) {
    printf("Error: unsupported subtype in make_type3_subchannel!\n");
    return 0;
  }

  set_E_bcs(&BC);
  sprintf(BC.cbc[0][0],"pin");
  sprintf(BC.cbc[0][1],"pin");

  corners[0]=shift[5];
  corners[1]=shift[10];
  corners[2]=shift[12];
  corners[3]=shift[4];
  make_garc_space(nx,ny,-0.5,1,delta,corners,BC.cbc);

  corners[0]=shift[10];
  corners[1]=shift[6];
  corners[2]=shift[7];
  corners[3]=shift[12];
  make_garc_space(nx,ny,-0.5,1,delta,corners,BC.cbc);

  set_E_bcs(&BC);
  sprintf(BC.cbc[2][0],"can");
  sprintf(BC.cbc[2][1],"can");

  double rr=0.5+distance(shift[6],shift[7]);
  if(subtype==1||subtype==2) corners[2]=shift[8];
  if(subtype==3||subtype==4) corners[2]=shift[13];
  corners[3]=shift[11];
  corners[0]=shift[12];
  corners[1]=shift[7];
  make_cgquad_space(nx,ng,-rr,-deltag,0.0,corners,BC.cbc);

  corners[2]=shift[11];
  if(subtype==1||subtype==4) corners[3]=shift[9];
  if(subtype==2||subtype==3) corners[3]=shift[14];
  corners[0]=shift[4];
  corners[1]=shift[12];
  make_cgquad_space(nx,ng,-rr,-deltag,0.0,corners,BC.cbc);

return 0;
}
