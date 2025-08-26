point get_fermat(point *p){
  static double psi = M_PI/3.0;

  point v[2],out;

  v[0]=rotate_point(p[1],-psi,p[0]);
  v[1]=rotate_point(p[1],psi,p[2]);
  out=line_line_intercept(p[0],v[1],v[0],p[2]);

  return out;
}

int make_tri_space(int ne,point *p,char bcs[4][2][4]){

  char bc[4][2][4];
  int i,j;

  point cp;
  point mp[3];
  point corners[4];

  for(i=0;i<3;i++){
    j=(i+1)%3;
    mp[i]=midpoint(p[i],p[j]);
  }
//cp=line_line_intercept(p[0],mp[1],p[2],mp[0]);
  cp=get_fermat(mp);

  corners[2]=cp;
  sprintf(bc[1][0],"E  ");
  sprintf(bc[1][1],"E  ");
  sprintf(bc[2][0],"E  ");
  sprintf(bc[2][1],"E  ");
  for(i=0;i<3;i++){
    j=(i+2)%3;
    corners[0]=p[i];
    corners[1]=mp[i];
    corners[3]=mp[j];
    strcpy(bc[0][0],bcs[i][0]);
    strcpy(bc[0][1],bcs[i][1]);
    strcpy(bc[3][0],bcs[j][0]);
    strcpy(bc[3][1],bcs[j][1]);
    make_quad_space(ne,ne,corners,bc);
  }

return 0;
}

int make_vtri_space(int nb,int nh,double *x,double *y){

  int i,j,ns,nr;

  double mx[3],my[3],gx=0.0,gy=0.0;
  double xq[4],yq[4];

/*
  for(i=0;i<3;i++){
    j=(i+1)%3;
    mx[i]=0.5*(x[i]+x[j]);
    my[i]=0.5*(y[i]+y[j]);
    gx+=x[i];
    gy+=y[i];
  }
  gx/=3.0;
  gy/=3.0;
*/
  mx[0]=0.5*(x[0]+x[1]);
  my[0]=0.5*(y[0]+y[1]);

  mx[1]=(nb*x[1]+nh*x[2])/(nh+nb);
  my[1]=(nb*y[1]+nh*y[2])/(nh+nb);

  mx[2]=(nh*x[2]+nb*x[0])/(nh+nb);
  my[2]=(nh*y[2]+nb*y[0])/(nh+nb);

  gx=(nh*x[2]+nb*nb*mx[0])/(nh+nb*nb);
  gy=(nh*y[2]+nb*nb*my[0])/(nh+nb*nb);

  for(i=0;i<3;i++){
    j=(i+2)%3;
    xq[0]=x[i];  yq[0]=y[i];
    xq[1]=mx[i]; yq[1]=my[i];
    xq[2]=gx;    yq[2]=gy;
    xq[3]=mx[j]; yq[3]=my[j];

    ns=nb; nr=nb;
    if(i==0) nr=nh;
    if(i==1) ns=nh;

//  make_quad_space(ns,nr,xq,yq);
  }

return 0;
}

int make_ctri_space(int nb,int nh,double R,point *p,char bcs[4][2][4]){

//int i,j,ns,nr;

//double mx[3],my[3],gx=0.0,gy=0.0;
//double xq[4],yq[4],xc,yc,aa,w[3];
  point m[3],cent,pc,vtx[4];

  double theta;

  pc=circle_center_2pR(p[0],p[1],R);
  theta=angle(pc,p[0],p[1])*R/fabs(R);
  m[0]=rotate_point(p[0],theta/2.0,pc);
  m[1]=midpoint(p[1],p[2]);
  m[2]=midpoint(p[2],p[0]);

  cent=centroid(p,3);

  vtx[0]=p[0];
  vtx[1]=m[0];
  vtx[2]=cent;
  vtx[3]=m[2];
  make_cquad_space(nb/2,nh/2,R,vtx,bcs);

  vtx[0]=m[0];
  vtx[1]=p[1];
  vtx[2]=m[1];
  vtx[3]=cent;
  make_cquad_space(nb/2,nh/2,R,vtx,bcs);

  vtx[0]=p[2];
  vtx[1]=m[2];
  vtx[2]=cent;
  vtx[3]=m[1];
  make_quad_space(nh/2,nb/2,vtx,bcs);

//w[0]=2.0*(double)nb;
//w[1]=w[0];
//w[2]=(double)nh;
//for(i=1;i<3;i++){
//  j=(i+1)%3;
//  (m+i)->x=((p+i)->x*w[i]+(p+j)->x*w[j])/(w[i]+w[j]);
//  (m+i)->y=((p+i)->y*w[i]+(p+j)->y*w[j])/(w[i]+w[j]);
//}

//p1.x=x[0];
//p1.y=y[0];
//p2.x=x[1];
//p2.y=y[1];

//xc=pc.x;
//yc=pc.y;

//aa=distance(m[0],pc); //sqrt((mx[0]-xc)*(mx[0]-xc)+(my[0]-yc)*(my[0]-yc));
//mx[0]=(mx[0]-xc)*fabs(R)/aa+xc;
//my[0]=(my[0]-yc)*fabs(R)/aa+yc;

//gx=(3.0*nb*mx[0]+nh*x[2])/(nh+3.0*nb);
//gy=(3.0*nb*my[0]+nh*y[2])/(nh+3.0*nb);

//xq[0]=x[0];  yq[0]=y[0];
//xq[1]=mx[0]; yq[1]=my[0];
//xq[2]=gx;    yq[2]=gy;
//xq[3]=mx[2]; yq[3]=my[2];
//make_cquad_space(nb,nh,R,dr0,0.0,xq,yq);

//xq[0]=mx[0]; yq[0]=my[0];
//xq[1]=x[1];  yq[1]=y[1];
//xq[2]=mx[1]; yq[2]=my[1];
//xq[3]=gx;    yq[3]=gy;
//make_cquad_space(nb,nh,R,dr0,0.0,xq,yq);

//xq[0]=x[2];  yq[0]=y[2];
//xq[1]=mx[2]; yq[1]=my[2];
//xq[2]=gx;    yq[2]=gy;
//xq[3]=mx[1]; yq[3]=my[1];
//make_quad_space(nb,nb,xq,yq);

return 0;
}

