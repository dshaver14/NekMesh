point line_circle_intercept(point p1,point p2,point pc,double rr){

  point pi;
  double dx2,dy2;
  double dx1=p2.x-p1.x;
  double dy1=p2.y-p1.y;
  if(fabs(dx1)>1.0e-12){
    double m=dy1/dx1;
    double c=-m*p1.x+p1.y;
    double A=(m*m+1.0);
    double B=2.0*(m*c-m*pc.y-pc.x);
    double C=pc.y*pc.y-rr*rr+pc.x*pc.x-2.0*c*pc.y+c*c;
    double des=B*B-4.0*A*C;
    if(des<0.0){
      printf("Error in line_circle_intercept: undefined intercept!\n%f, %f\t%f, %f\n%f, %f\t,%f\n",p1.x,p1.y,p2.x,p2.y,pc.x,pc.y,rr);
      pi.x=0.0;
      pi.y=0.0;
      return pi;
    }
    //Always want the intercept closest to x1,y1
    dx1=(-B-sqrt(des))/(2.0*A);
    dx2=(-B+sqrt(des))/(2.0*A);
    dy1=m*dx1+c-p1.y;
    dy2=m*dx2+c-p1.y;
    dx1-=p1.x;
    dx2-=p1.x;
    dx1=sqrt(dx1*dx1+dy1*dy1);
    dx2=sqrt(dx2*dx2+dy2*dy2);
    if(dx1<dx2){
      pi.x=(-B-sqrt(des))/(2.0*A);
    }else{
      pi.x=(-B+sqrt(des))/(2.0*A);
    }
    pi.y=m*pi.x+c;
  }else{
    pi.x=p1.x;
    dy1= sqrt(rr*rr-(p1.x-pc.x)*(p1.x-pc.x))+pc.y;
    dy2=-sqrt(rr*rr-(p1.x-pc.x)*(p1.x-pc.x))+pc.y;
    if(fabs(dy1-p1.y)<fabs(dy2-p1.y)){
     pi.y=dy1;
    }else{
     pi.y=dy2;
    }
  }

  return pi;
}

point line_line_intercept(point p1,point p2,point p3, point p4){
  point p0={0.0,0.0};
  double denom=(p1.x-p2.x)*(p3.y-p4.y)-(p1.y-p2.y)*(p3.x-p4.x);
  if(fabs(denom)<1.0e-12){
    printf("Error in line_line_intercept: lines are parallel!");
    return p0;
  }
  p0.x=(p1.x*p2.y-p1.y*p2.x)*(p3.x-p4.x)-(p1.x-p2.x)*(p3.x*p4.y-p3.y*p4.x);
  p0.x/=denom;
  p0.y=(p1.x*p2.y-p1.y*p2.x)*(p3.y-p4.y)-(p1.y-p2.y)*(p3.x*p4.y-p3.y*p4.x);
  p0.y/=denom;

  return p0;
}

point circle_center_2pR(point p0,point p1,double R){

  double delta = distance(p0,p1);
  double aa = delta/2.0;
  double bb = sqrt(R*R-aa*aa);
  double xa = 0.5*(p1.x-p0.x); double ya = 0.5*(p1.y-p0.y);
  point pc;
  if(R>0.0){ 
    pc.x = p0.x+xa-bb*ya/aa; 
    pc.y = p0.y+ya+bb*xa/aa;
  }else if(R<0.0){
    pc.x = p0.x+xa+bb*ya/aa;
    pc.y = p0.y+ya-bb*xa/aa;
  }

  return pc;
}

point rotate_point(point pin,double theta,point zero){
  point pout;

  pout.x = (pin.x-zero.x)*cos(theta)-(pin.y-zero.y)*sin(theta)+zero.x;
  pout.y = (pin.x-zero.x)*sin(theta)+(pin.y-zero.y)*cos(theta)+zero.y;

 return pout;
}

double distance(point p1,point p2){
  double dist=sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
  return dist;
}

point reflect_point(point pin,point p0,point p1){
  point pout;
  double a=p1.y-p0.y;
  double b=p0.x-p1.x;
  double c=(p1.x-p0.x)*p0.y-(p1.y-p0.y)*p0.x;
  double RHS=-2.0*(a*pin.x+b*pin.y+c)/(a*a+b*b);
  pout.x=RHS*a+pin.x;
  pout.y=RHS*b+pin.y;
  return pout;
}

point translate_point(point pin,point del){
  point pout;
  pout.x=pin.x+del.x;
  pout.y=pin.y+del.y;
  pout.z=pin.z+del.z;
  return pout;
}

point midpoint(point p1,point p2){
  point pout;
  pout.x=0.5*(p1.x+p2.x);
  pout.y=0.5*(p1.y+p2.y);
  return pout;
}

point linpoint(double del,point p1,point p2){
  point pout;
  pout.x=p1.x+(p2.x-p1.x)*del;
  pout.y=p1.y+(p2.y-p1.y)*del;
  return pout;
}

point centroid(point *p,int np){
  point pout={0.0,0.0};
  int i;
  for(i=0;i<np;i++){
    pout.x+=(p+i)->x;
    pout.y+=(p+i)->y;
  }
  pout.x/=(double)np;
  pout.y/=(double)np;

  return pout;
}

double angle(point pc,point p1,point p2){

  double del1,del2,dot,theta;

  dot=(p1.x-pc.x)*(p2.x-pc.x)+(p1.y-pc.y)*(p2.y-pc.y);
  del1=distance(p1,pc);
  del2=distance(p2,pc);
  theta=acos(dot/(del1*del2));

return theta;
}

double crossprod2d(point pc,point p1,point p2){

  point v1,v2;

  v1.x=p1.x-pc.x;
  v1.y=p1.y-pc.y;
  v2.x=p2.x-pc.x;
  v2.y=p2.y-pc.y;

  return v1.x*v2.y-v1.y*v2.x;
}
int rotate_all(double theta,point pc){

  int i;point pttmp;
  for(i=0;i<nvert;i++){
    pttmp=rotate_point(verts[i],theta,pc);
    verts[i]=pttmp;
  }

}
point set_point(double xin,double yin){
  point pout; 
  pout.x=xin;
  pout.y=yin;

  return pout;
}
