int make_quad_space(int nr,int ns,point *p,char bcs[4][2][4]){

  int i,j,k,nvert0=nvert;
  connector s2,s4,rr;

  printf("Generating %d elements in a quadrilateral space\n",ns*nr);

  initialize_con(&s2,ns);
  initialize_con(&s4,ns);
  initialize_con(&rr,nr);

//      3---s3---2
//      |        |
// j,s  s4      s2
//  ᴧ   |        |
//  |   0---s1---1
//  |
//  *---->i,r

  set_lin_side(p[1],p[2],&s2); //side 2 (1-2)
  set_lin_side(p[0],p[3],&s4); //side 4 (0-3)
  for(j=0;j<=ns;j++){
    set_lin_side(s4.p[j],s2.p[j],&rr);
    for(i=0;i<=nr;i++) verts[nvert0+i+j*(nr+1)]=rr.p[i];
    nvert+=rr.n;
  }

  dealloc_con(&s2);
  dealloc_con(&s4);
  dealloc_con(&rr);

  for(j=0;j<ns;j++){
    for(i=0;i<nr;i++){
      (elems+nelem)->vid[0]=(i+0)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[1]=(i+1)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[2]=(i+1)+(j+1)*(nr+1)+nvert0;
      (elems+nelem)->vid[3]=(i+0)+(j+1)*(nr+1)+nvert0;
      for(k=0;k<nfld;k++){
        if(nr==1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
        }else if(i==0){
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
        if(ns==1){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
        }else if(j==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(j==ns-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
      }
      nelem++;
    }
  }

  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_cquad_space(int nt,int nr,double R1,point *p,char bcs[4][2][4]){

  int i,j,k,nvert0=nvert;
  connector s2,s4,rr;
  double *rb;
  point pc;

  printf("Generating %d elements in a curved quadrilateral space\n",nt*nr);


//      3----s3----2
//      |          |
//      |          |
//      s4        s2
// j,r  |          |
//  ᴧ   | /--s1--\ |
//  |   0/        \1
//  |
//  *---->i,t  (r,theta)

  initialize_con(&s2,nr);
  initialize_con(&s4,nr);
  initialize_con(&rr,nt);
  rb=malloc(sizeof(double)*nr);

  set_lin_side(p[1],p[2],&s2);
  set_lin_side(p[0],p[3],&s4);
  for(j=0;j<nr;j++){
    rb[j]=R1/(1.-(double)j/(double)nr);
    set_arc_side(s4.p[j],s2.p[j],rb[j],&rr);
    for(i=0;i<=nt;i++) verts[nvert0+i+j*(nt+1)]=rr.p[i];
    nvert+=rr.n;
  } 
  set_lin_side(s4.p[nr],s2.p[nr],&rr);
  for(i=0;i<=nt;i++) verts[nvert0+i+nr*(nt+1)]=rr.p[i];
  nvert+=rr.n;

  dealloc_con(&s2);
  dealloc_con(&s4);
  dealloc_con(&rr);

  for(i=0;i<nr;i++){
    for(j=0;j<nt;j++){
      (elems+nelem)->vid[0]=(j+0)+(i+0)*(nt+1)+nvert0;
      (elems+nelem)->vid[1]=(j+1)+(i+0)*(nt+1)+nvert0;
      (elems+nelem)->vid[2]=(j+1)+(i+1)*(nt+1)+nvert0;
      (elems+nelem)->vid[3]=(j+0)+(i+1)*(nt+1)+nvert0;
      for(k=0;k<nfld;k++){
        if(nr==1){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
        }else if(i==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
        if(nt==1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
        }else if(j==0){ //j is build from side 4 [3][k] to 2 [1][k]
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(j==nt-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
      }
      if(i<nr-1){
        (cides+ncide)->elid=nelem;
        (cides+ncide)->esid=1;
        (cides+ncide)->curve=*(rb+i);
        (cides+ncide)->ccurve='C';
        ncide++;
        (cides+ncide)->elid=nelem;
        (cides+ncide)->esid=3;
        (cides+ncide)->curve=-*(rb+i+1);
        (cides+ncide)->ccurve='C';
        ncide++;
      }else{
        (cides+ncide)->elid=nelem;
        (cides+ncide)->esid=1;
        (cides+ncide)->curve=*(rb+i);
        (cides+ncide)->ccurve='C';
        ncide++;
      }
      nelem++;
    }
  }

  free(rb);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_gquad_space(int nr,int ns,double dn0,point *p,char bcs[4][2][4]){

  int i,j,k,nvert0=nvert;
  double dr,ds;
  connector r2,r4;

  printf("Generating %d elements in a quadrilateral space with geometric growth\n",ns*nr);

  initialize_con(&r2,ns);
  initialize_con(&r4,ns);

//side 2 (1-2)
  set_geo_side(p[1],p[2],dn0,&r2);
//side 4 (0-3)
  set_geo_side(p[0],p[3],dn0,&r4);

  for(j=0;j<=ns;j++){
    dr=(r2.p[j].x-r4.p[j].x)/nr;
    ds=(r2.p[j].y-r4.p[j].y)/nr;
    for(i=0;i<=nr;i++){
      verts[nvert0+i+j*(nr+1)].x=r4.p[j].x+(double)i*dr;
      verts[nvert0+i+j*(nr+1)].y=r4.p[j].y+(double)i*ds;
      nvert++;
    }
  }

  for(j=0;j<ns;j++){
    for(i=0;i<nr;i++){
      (elems+nelem)->vid[0]=(i+0)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[1]=(i+1)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[2]=(i+1)+(j+1)*(nr+1)+nvert0;
      (elems+nelem)->vid[3]=(i+0)+(j+1)*(nr+1)+nvert0;
      for(k=0;k<nfld;k++){
        if(nr==1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
        }else if(i==0){
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
        if(ns==1){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
        }else if(j==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(j==ns-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
      }
      nelem++;
    }
  }

  free(r2.p);free(r4.p);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_g2quad_space(int nr,int ns,double dr0,double ds0,point *p,char bcs[4][2][4]){

  int i,j,k,nvert0=nvert;
  point r0,r1;
  connector r,r2,r4;

  printf("Generating %d elements with geometric growth in two directions\n",ns*nr);

  initialize_con(&r2,ns);
  initialize_con(&r4,ns);
  initialize_con(&r,nr);

//side 2 (1-2)
  set_geo_side(p[1],p[2],ds0,&r2);
//side 4 (0-3)
  set_geo_side(p[0],p[3],ds0,&r4);

  for(j=0;j<=ns;j++){
    if(dr0<0.0){
      r0=r4.p[j];
      r1=r2.p[j];
    }else{
      r1=r4.p[j];
      r0=r2.p[j];
    }
    get_g_side(r0,r1,fabs(dr0),&r);
    if(dr0>0.0) invert_connector(&r);
    for(i=0;i<=nr;i++) verts[nvert0+i+j*(nr+1)]=r.p[i];
    nvert+=r.n;
  }

  for(j=0;j<ns;j++){
    for(i=0;i<nr;i++){
      (elems+nelem)->vid[0]=(i+0)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[1]=(i+1)+(j+0)*(nr+1)+nvert0;
      (elems+nelem)->vid[2]=(i+1)+(j+1)*(nr+1)+nvert0;
      (elems+nelem)->vid[3]=(i+0)+(j+1)*(nr+1)+nvert0;
      for(k=0;k<nfld;k++){
        if(nr==1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
        }else if(i==0){
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
        if(ns==1){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
        }else if(j==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(j==ns-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
      }
      nelem++;
    }
  }

  free(r2.p);free(r4.p),free(r.p);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}

int make_cgquad_space(int nt,int nr,double R,double dr0,double dt0,point *p,char bcs[4][2][4]){
  
  int i,j,k,nvert0=nvert;
  double theta,thet0,delta;
  double *ta,*rb;
  point pc,r0,r1;
  connector r,r2,r3,r4;

  printf("Generating %d elements in a curved quadrilateral space with geometric growth\n",nt*nr);
  
//for(k=0;k<nfld;k++) for(i=0;i<4;i++) printf("\t\"%s\"\n",bcs[i][k]);
  
  pc=circle_center_2pR(p[0],p[1],R);

  initialize_con(&r,nr);
  initialize_con(&r2,nr);
  initialize_con(&r3,nt);
  initialize_con(&r4,nr);

  ta=(double *)malloc(sizeof(double)*(nt+1));
  rb=(double *)malloc(sizeof(double)*(nr+1));

//      3-----s3-----2
//      |            |
//      |            |
//      s4          s2
// j,r  |            |
//  ᴧ   | /---s1---\ |
//  |   0/          \1
//  |
//  *---->i,t  (r,theta)

  *rb=R;
  *(rb+1)=R-fabs(dr0);

  if(dr0>0.0){
  //geometric growth from arc side (s1)
    //side 2 (1-2)
    r0=p[1];
    r1=line_circle_intercept(p[1],p[2],pc,*(rb+1));
    delta=distance(r0,r1);
    get_g_side(p[1],p[2],delta,&r2);

    //side 4 (0-3)
    r0=p[0];
    r1=line_circle_intercept(p[0],p[3],pc,*(rb+1));
    delta=distance(r0,r1);
    get_g_side(p[0],p[3],delta,&r4);
  }else if(dr0<0.0){
  //geometric growth from linear side (s3)
    //side 2 (1-2)
    get_g_side(p[2],p[1],-dr0,&r2);
    invert_connector(&r2);
    //side 4 (0-3)
    get_g_side(p[3],p[0],-dr0,&r4);
    invert_connector(&r4);
    *(rb+1)=R+(double)nr*fabs(dr0);
  }else{
    set_lin_side(p[1],p[2],&r2);
    set_lin_side(p[0],p[3],&r4);
    dr0=0.5*(distance(p[1],p[2])+distance(p[0],p[3]))/((double)nr);
  }

  delta=distance(r4.p[0],r2.p[0]);
  theta=2.0*asin(delta/(2.0*fabs(R)));
  
  for(i=2;i<=nr;i++){
    if(R>0.0) *(rb+i)=*(rb+i-1)+0.5*(distance(r2.p[i],r2.p[i-1])+distance(r4.p[i],r4.p[i-1]));
    if(R<0.0) *(rb+i)=*(rb+i-1)-pow(1.1,(double)(i-1))*fmax(distance(r2.p[i],r2.p[i-1]),distance(r4.p[i],r4.p[i-1]));
  }

//side 3 (3-2), theta1a, theta1b
  if(dt0>0.0){ //geometric growth from side 4
    get_g_side(p[3],p[2],dt0,&r3);
    thet0=get_theta_0(r4.p[0],p[3],0,dt0,pc);
    get_g1D(0.0,-theta,nt,thet0,ta);
  }else if(dt0<0.0){ //geometric growth from side 2
    get_g_side(p[2],p[3],-dt0,&r3);
    thet0=get_theta_0(r2.p[0],p[2],1,-dt0,pc);
    get_g1D(0.0,-theta,nt,thet0,ta);
    //re-orient from side 4
    invert_connector(&r3);
    for(i=0;i<=nt;i++) *(ta+i)=-(theta+*(ta+i));
    invert(ta,nt+1);
  }else{ //linear distribution
    for(i=0;i<=nt;i++){
      r3.p[i].x=p[3].x+(p[2].x-p[3].x)*((double)i/(double)(nt));
      r3.p[i].y=p[3].y+(p[2].y-p[3].y)*((double)i/(double)(nt));
      *(ta+i)=-((double)i)/((double)nt)*theta;
    }
  }

  for(i=0;i<=nt;i++){
    int j0;
    //j=0
    theta=-R/fabs(R)* *(ta+i);
    verts[nvert0+i]=rotate_point(r4.p[0],theta,pc);
    //j=1
    if(dr0>=0.0){
      verts[nvert0+i+nt+1]=line_circle_intercept(verts[nvert0+i],r3.p[i],pc,*(rb+1));
      j0=2;  
    }else j0=1;
 
    for(j=j0;j<=nr;j++){
      if(j==j0){ 
        if(dr0>0.0) { get_g_side(verts[nvert0+i],r3.p[i],dr0,&r);
        }else if (dr0<0.0){
          get_g_side(r3.p[i],verts[nvert0+i],-dr0,&r);
          invert_connector(&r);
        }else {set_lin_side(verts[nvert0+i],r3.p[i],&r);}
      }
      verts[nvert0+i+j*(nt+1)]=r.p[j];
    }
  }
  nvert+=(nr+1)*(nt+1);

  for(i=0;i<nr;i++){
    for(j=0;j<nt;j++){
      (elems+nelem)->vid[0]=(j+0)+(i+0)*(nt+1)+nvert0;
      (elems+nelem)->vid[1]=(j+1)+(i+0)*(nt+1)+nvert0;
      (elems+nelem)->vid[2]=(j+1)+(i+1)*(nt+1)+nvert0;
      (elems+nelem)->vid[3]=(j+0)+(i+1)*(nt+1)+nvert0;
      for(k=0;k<nfld;k++){
        if(nr==1){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
        }else if(i==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
        if(nt==1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
        }else if(j==0){
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(j==nt-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
      }
      if(dr0>0.0){
        if(i<nr-1){
          (cides+ncide)->elid=nelem;
          (cides+ncide)->esid=1;
          (cides+ncide)->curve=*(rb+i);
          (cides+ncide)->ccurve='C';
          ncide++;
          (cides+ncide)->elid=nelem;
          (cides+ncide)->esid=3;
          (cides+ncide)->curve=-*(rb+i+1);
          (cides+ncide)->ccurve='C';
          ncide++;
        }else{
          (cides+ncide)->elid=nelem;
          (cides+ncide)->esid=1;
          (cides+ncide)->curve=*(rb+i);
          (cides+ncide)->ccurve='C';
          ncide++;
        }
      }else{
        if(i==0){
          (cides+ncide)->elid=nelem;
          (cides+ncide)->esid=1;
          (cides+ncide)->curve=R;
          (cides+ncide)->ccurve='C';
          ncide++;
        }
      }
      nelem++;
    }
  }

  dealloc_con(&r);
  dealloc_con(&r2);
  dealloc_con(&r3);
  dealloc_con(&r4);
/*
  free(r.p);
  free(r2.p);
  free(r3.p);
  free(r4.p);
*/
  free(ta);free(rb);
  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}
//---------------------------------------------------------------------------------------
int make_arc_space(int nt,int nr,double R1,point *p,char bcs[4][2][4]){

  int i,j,k,nvert0=nvert;
  double dr,R3,R3a;
  connector s4,s2,rr;
  point pc;

//       3-----s3---2
//      /            \
//      |            |
//      s4          s2
// j,r  |            |
//  ᴧ   | /---s1---\ |
//  |   0/          \1
//  |
//  *---->i,t  (r,theta)

  printf("Generating %d elements in a concentric arc space\n",nt*nr);

//  check if sides 1 and 3 are concentric
  pc=circle_center_2pR(p[0],p[1],R1);
  R3=distance(p[2],pc);
  R3a=distance(p[3],pc);
  if(fabs(R3-R3a)>1.0e-8) {
    printf("Error in arc_space, sides 1 and 3 not concentric!\n");
    for(i=0;i<4;i++)printf("\t corner %d: %f, %f\n",i,p[i].x,p[i].y);
    return -1;
  }

  initialize_con(&s4,nr);
  initialize_con(&s2,nr);
  initialize_con(&rr,nt);

  R3*=R1/fabs(R1);
  dr=(R1-R3)/((double)nr);

  for(i=0;i<=nr;i++)s4.p[i]=line_circle_intercept(p[0],p[3],pc,R1-(double)i*dr); //side 4
  for(i=0;i<=nr;i++)s2.p[i]=line_circle_intercept(p[1],p[2],pc,R1-(double)i*dr); //side 2

//vertices
  for(j=0;j<=nr;j++){
    set_arc_side(s4.p[j],s2.p[j],R1-(double)j*dr,&rr);
    for(i=0;i<=nt;i++) verts[nvert0+i+j*(nt+1)]=rr.p[i];
    nvert+=rr.n;
  }

  dealloc_con(&s4);
  dealloc_con(&s2);
  dealloc_con(&rr);

//elements
  for(i=0;i<nr;i++){
    for(j=0;j<nt;j++){
      (elems+nelem)->vid[0]=(j+0)+(i+0)*(nt+1)+nvert0;
      (elems+nelem)->vid[1]=(j+1)+(i+0)*(nt+1)+nvert0;
      (elems+nelem)->vid[2]=(j+1)+(i+1)*(nt+1)+nvert0;
      (elems+nelem)->vid[3]=(j+0)+(i+1)*(nt+1)+nvert0;
      for(k=0;k<nfld;k++){
        if(nr==1){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
        }else if(i==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
        if(nt==1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
        }else if(j==0){
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(j==nt-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
      }
      (cides+ncide)->elid=nelem;
      (cides+ncide)->esid=1;
      (cides+ncide)->curve=R1-(double)i*dr; 
      (cides+ncide)->ccurve='C';
      ncide++; 
      (cides+ncide)->elid=nelem;
      (cides+ncide)->esid=3;
      (cides+ncide)->curve=-R1+(double)(i+1)*dr; 
      (cides+ncide)->ccurve='C';
      ncide++; 
      nelem++;
    }
  }

  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}
//---------------------------------------------------------------------------------------
int make_garc_space(int nt,int nr,double R1,int isd,double r0,point *p,char bcs[4][2][4]){

  int i,j,k,nvert0=nvert;
  double dr,R3,R3a,*ss,t0;
  connector s4,s2,rr;
  point pc,p2a,p3a;

//   geometric growth from side isd
//       3-----s3---2
//      /            \
//      |            |
//      s4          s2
// i,r  |            |
//  ᴧ   | /---s1---\ |
//  |   0/          \1
//  |
//  *---->j,t  (r,theta)

  printf("Generating %d elements in a concentric arc space with geometric growth from side %d\n",nt*nr,isd);

//  check if sides 1 and 3 are concentric
  pc=circle_center_2pR(p[0],p[1],R1);
  R3=distance(p[2],pc);
  R3a=distance(p[3],pc);
  if(fabs(R3-R3a)>1.0e-8) {
    printf("Error in arc_space, sides 1 and 3 not concentric!\n");
    for(i=0;i<4;i++)printf("\t corner %d: %f, %f\n",i,p[i].x,p[i].y);
    return -1;
  }

  p2a=linpoint(1.05,p[1],p[2]);
  p3a=linpoint(1.05,p[0],p[3]);

  initialize_con(&s4,nr);
  initialize_con(&s2,nr);
  initialize_con(&rr,nt);
  ss=malloc(sizeof(double)*(nr+1));

  R3*=R1/fabs(R1);
  if(isd==1||isd==5||isd==6){
    get_g1D(R1,R3,nr,r0,ss);
  }else if(isd==3||isd==7||isd==8){ 
    get_g1D(R3,R1,nr,r0,ss); 
    invert(ss,nr+1);
  }else{
    dr=(R3-R1)/(double)(nr);
    for(i=0;i<=nr;i++)ss[i]=R1+(double)i*dr;
  }

  for(i=0;i<=nr;i++)s4.p[i]=line_circle_intercept(p[0],p[3],pc,ss[i]); //side 4
//for(i=0;i<=nr;i++)printf("s4 %f, %f\n",s4.p[i].x,s4.p[i].y);
  for(i=0;i<=nr;i++)s2.p[i]=line_circle_intercept(p[1],p[2],pc,ss[i]); //side 2
//for(i=0;i<=nr;i++)printf("s2 %f, %f\n",s2.p[i].x,s2.p[i].y);

//vertices
  for(j=0;j<=nr;j++){
    if(isd==2||isd==5||isd==7){
      t0=get_theta_0(s2.p[j],p2a,1,r0,pc);
//    printf("t0 = %f\n",t0);
      set_garc_side(s2.p[j],s4.p[j],pc,t0,&rr);
      invert_connector(&rr);
    }else if(isd==4||isd==6||isd==8){
      t0=get_theta_0(s4.p[j],p3a,0,r0,pc);
//    printf("t0 = %f\n",t0);
      set_garc_side(s4.p[j],s2.p[j],pc,t0,&rr);
    }else{
      set_arc_side(s4.p[j],s2.p[j],ss[j],&rr);
    }
    for(i=0;i<=nt;i++) verts[nvert0+i+j*(nt+1)]=rr.p[i];
    nvert+=rr.n;
  }

dealloc_con(&s4);
dealloc_con(&s2);
dealloc_con(&rr);

//elements
  for(i=0;i<nr;i++){
    for(j=0;j<nt;j++){
      (elems+nelem)->vid[0]=(j+0)+(i+0)*(nt+1)+nvert0;
      (elems+nelem)->vid[1]=(j+1)+(i+0)*(nt+1)+nvert0;
      (elems+nelem)->vid[2]=(j+1)+(i+1)*(nt+1)+nvert0;
      (elems+nelem)->vid[3]=(j+0)+(i+1)*(nt+1)+nvert0;
      for(k=0;k<nfld;k++){
        if(nr==1){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
        }else if(i==0){
          strcpy((elems+nelem)->BC[0][k],bcs[0][k]);
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }else if(i==nr-1){
          strcpy((elems+nelem)->BC[2][k],bcs[2][k]);
          sprintf((elems+nelem)->BC[0][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[0][k],"E  ");
          sprintf((elems+nelem)->BC[2][k],"E  ");
        }
        if(nt==1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
        }else if(j==0){
          strcpy((elems+nelem)->BC[3][k],bcs[3][k]);
          sprintf((elems+nelem)->BC[1][k],"E  ");
        }else if(j==nt-1){
          strcpy((elems+nelem)->BC[1][k],bcs[1][k]);
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }else{
          sprintf((elems+nelem)->BC[1][k],"E  ");
          sprintf((elems+nelem)->BC[3][k],"E  ");
        }
      }
      (cides+ncide)->elid=nelem;
      (cides+ncide)->esid=1;
      (cides+ncide)->curve=ss[i]; 
      (cides+ncide)->ccurve='C';
      ncide++; 
      (cides+ncide)->elid=nelem;
      (cides+ncide)->esid=3;
      (cides+ncide)->curve=-ss[i+1]; 
      (cides+ncide)->ccurve='C';
      ncide++; 
      nelem++;
    }
  }

  printf("nelem = %d, nvert = %d\n",nelem,nvert);

return 0;
}
//---------------------------------------------------------------------------------------
int set_4_corners(point *cnr,int c1,int c2,int c3,int c4){

     cnr[0]=points[c1];
     cnr[1]=points[c2];
     cnr[2]=points[c3];
     cnr[3]=points[c4];

return 0;
}
int set_4_bcs(char bcs[4][2][4],int ifld,char *b1,char *b2,char *b3,char *b4){

  strncpy(bcs[0][ifld-1],b1,4);
  strncpy(bcs[1][ifld-1],b2,4);
  strncpy(bcs[2][ifld-1],b3,4);
  strncpy(bcs[3][ifld-1],b4,4);

return 0;
}
