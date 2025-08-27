int write_rea(char *fname){

  char line[256],bcout[4];

  FILE *reain,*reaout;

  int i,j,k,ivert,inc=0;
  double dx;
  double xmin=1.0e30,xmax=-1.0e30,ymin=1.0e30,ymax=-1.0e30;

  if(nelem>=1000000){
    printf("Error: rea printout for 1e6 elements or more NOT supported!\n\tNo rea file produced\n");
    return -1;
  }

  printf("\n\telement consistency check...\n");
  for(i=0;i<nelem;i++){
    if(axi_consistency_check(elems+i)!=0){
      inc++;
    }
  }
  if(inc>0) printf("\trotated %d elements to align 'A  ' BCs with face 1\n",inc);

  for(i=0;i<nvert;i++){
    xmax=fmax(xmax,(verts+i)->x);
    xmin=fmin(xmin,(verts+i)->x);
    ymax=fmax(ymax,(verts+i)->y);
    ymin=fmin(ymin,(verts+i)->y);
  }

  dx=fmax((xmax-xmin),(ymax-ymin));

  reaout=fopen(fname,"w");

  printf("\n\twriting %d elements to file \"%s\"\n\n",nelem,fname);

  reain=fopen("core/reatop.txt","r");
  while(fgets(line,256,reain)!=NULL) fprintf(reaout,"%s",line);
  fclose(reain);

  fprintf(reaout,"           13  LOGICAL SWITCHES FOLLOW\n");
  fprintf(reaout,"  T     IFFLOW\n");
  if(nfld==1){ 
    fprintf(reaout,"  F     IFHEAT\n");
  }else{
    fprintf(reaout,"  T     IFHEAT\n");
  }
  fprintf(reaout,"  T     IFTRAN\n");
  fprintf(reaout,"  T T F F F F F F F F F IFNAV & IFADVC (convection in P.S. fields)\n");
  fprintf(reaout,"  F F T T T T T T T T T T IFTMSH (IF mesh for this field is T mesh)\n");
  fprintf(reaout,"  F     IFAXIS\n");
  fprintf(reaout,"  F     IFSTRS\n");
  fprintf(reaout,"  F     IFSPLIT\n");
  fprintf(reaout,"  F     IFMGRID\n");
  fprintf(reaout,"  F     IFMODEL\n");
  fprintf(reaout,"  F     IFKEPS\n");
  fprintf(reaout,"  F     IFMVBD\n");
  fprintf(reaout,"  F     IFCHAR\n");

  fprintf(reaout,"%10.6f %10.6f %10.6f %10.6f   XFAC,YFAC,XZERO,YZERO\n",1.2*dx,1.2*dx,xmin-0.1*dx,ymin-0.1*dx);
  fprintf(reaout,"**MESH DATA** 1st line is X of corner 1,2,3,4. 2nd line is Y.\n");

  fprintf(reaout,"    %5d         2      %5d NEL,NDIM,NELV\n",nelem,nelem);
  for(i=0;i<nelem;i++){
    fprintf(reaout,"            ELEMENT%5d [    1a]    GROUP     0\n",i+1);
    for(j=0;j<4;j++){
      ivert=(elems+i)->vid[j];
      fprintf(reaout,"%14.6e",(verts+ivert)->x);
    }
    fprintf(reaout,"\n");
    for(j=0;j<4;j++){
      ivert=(elems+i)->vid[j];
      fprintf(reaout,"%14.6e",(verts+ivert)->y);
    }
    fprintf(reaout,"\n");
  }

  fprintf(reaout,"  ***** CURVED SIDE DATA *****\n");
  if(ncide==0){
    fprintf(reaout,"0 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n");
  }else{
    fprintf(reaout,"%6d Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n",ncide);
    if(nelem<1000){for(i=0;i<ncide;i++)
      fprintf(reaout,"%3d%3d%14.6f%14.6g%14.6g%14.6g%14.6g %c\n",
        (cides+i)->esid,(cides+i)->elid+1,(cides+i)->curve,0.0,0.0,0.0,0.0,(cides+i)->ccurve);
    }else{ for(i=0;i<ncide;i++)
      fprintf(reaout,"%2d%6d%14.6f%14.6g%14.6g%14.6g%14.6g %c\n",
        (cides+i)->esid,(cides+i)->elid+1,(cides+i)->curve,0.0,0.0,0.0,0.0,(cides+i)->ccurve);
    }
  }

  fprintf(reaout,"  ***** BOUNDARY CONDITIONS *****\n");
  fprintf(reaout,"  ***** FLUID   BOUNDARY CONDITIONS *****\n");
  k=0;
  for(i=0;i<nelem;i++){
    for(j=0;j<4;j++){
      strcpy(bcout,(elems+i)->BC[j][k]);
      if(nelem<999){
        fprintf(reaout," %s%3d  %1d %13.5f %13.5f %13.5f %13.5f %13.5f\n",bcout,i+1,j+1,0.0,0.0,0.0,0.0,0.0);
      }else{
        fprintf(reaout," %s%5d%1d %13.5f %13.5f %13.5f %13.5f %13.5f\n",bcout,i+1,j+1,0.0,0.0,0.0,0.0,0.0);
      }
    }
  }
  if(nfld==1){
    fprintf(reaout,"  ***** NO THERMAL BOUNDARY CONDITIONS *****\n");
  }else{
    fprintf(reaout,"  ***** THERMAL BOUNDARY CONDITIONS *****\n");
    k=1;
    for(i=0;i<nelem;i++){
      for(j=0;j<4;j++){
        strcpy(bcout,(elems+i)->BC[j][k]);
        if(nelem<999){
          fprintf(reaout," %s%3d  %1d %13.5f %13.5f %13.5f %13.5f %13.5f\n",bcout,i+1,j+1,0.0,0.0,0.0,0.0,0.0);
        }else{
          fprintf(reaout," %s%5d%1d %13.5f %13.5f %13.5f %13.5f %13.5f\n",bcout,i+1,j+1,0.0,0.0,0.0,0.0,0.0);
        }
      }
    }
  }

  reain=fopen("core/reabot.txt","r");
  while(fgets(line,256,reain)!=NULL) fprintf(reaout,"%s",line);
  fclose(reain);

  fclose(reaout);

return 0;
}
//-------------------------------------------------------------------
int axi_consistency_check(quad *elem){
  int i,j;

  for(i=1;i<4;i++){
    if(strncmp("A  ",elem->BC[i][0],3)==0){
//    printf("Error: side %d set to \"A  \"\n",i+1);
      for(j=i;j>0;j--){
//      printf("%d %d %d %d\n",elem->vid[0],elem->vid[1],elem->vid[2],elem->vid[3]);
//      printf("\"%s\" \"%s\" \"%s\" \"%s\"\n",elem->BC[0][0],elem->BC[1][0],elem->BC[2][0],elem->BC[3][0]);
        rotate_element(elem);
//      printf("%d %d %d %d\n",elem->vid[0],elem->vid[1],elem->vid[2],elem->vid[3]);
//      printf("\"%s\" \"%s\" \"%s\" \"%s\"\n\n",elem->BC[0][0],elem->BC[1][0],elem->BC[2][0],elem->BC[3][0]);
      }   
      return i;
    }   
  }

return 0;
}
//-------------------------------------------------------------------
int rotate_element(quad *elem){
  char cbuf[2][4];
  int i,ibuf;

  ibuf=elem->vid[0];
  strcpy(cbuf[0],elem->BC[0][0]);
  strcpy(cbuf[1],elem->BC[0][1]);
  for(i=0;i<3;i++){
    elem->vid[i]=elem->vid[i+1];
    strcpy(elem->BC[i][0],elem->BC[i+1][0]);
    strcpy(elem->BC[i][1],elem->BC[i+1][0]);
  }
  elem->vid[3]=ibuf;
  strcpy(elem->BC[3][0],cbuf[0]);
  strcpy(elem->BC[3][1],cbuf[1]);

return 0;
}
//-------------------------------------------------------------------
int write_vtk(char *fname){

  FILE *vtkout;

  vtkout = fopen(fname,"w");

  fprintf(vtkout,"%s\n","# vtk DataFile Version 3.0");
  fprintf(vtkout,"%s\n",fname);
  fprintf(vtkout,"%s\n","ASCII");
  fprintf(vtkout,"%s\n","DATASET UNSTRUCTURED_GRID");

//Points section
  fprintf(vtkout,"%s%d%s\n","POINTS ",nvert," float");
  for(int ivrt=0;ivrt<nvert;ivrt++) fprintf(vtkout,"%f %f %f\n",verts[ivrt].x,verts[ivrt].y,verts[ivrt].z);
  fprintf(vtkout,"\n"); //write blank line

//Elements section
  fprintf(vtkout,"%s %d %d\n","CELLS",nelem,nelem*5);
  for(int iel=0;iel<nelem;iel++) fprintf(vtkout,"%d %d %d %d %d\n",4,
    elems[iel].vid[0],
    elems[iel].vid[1],
    elems[iel].vid[2],
    elems[iel].vid[3]);
  fprintf(vtkout,"\n"); //write blank line

//Element type section
  fprintf(vtkout,"%s %d\n","CELL_TYPES",nelem);
  for(int iel=0;iel<nelem;iel++) fprintf(vtkout,"%d\n",9);

  fclose(vtkout);

return 0;
}
//-------------------------------------------------------------------
int reset(void){
  printf("\n\treseting mesh generator\n\n");
  nelem = 0;
  nvert = 0;
  ncide = 0;
}
//-------------------------------------------------------------------
int output_pts(point *pts,int npts,char *fnm){

  int i;
  FILE *ptout;

  ptout=fopen(fnm,"w");
  for(i=0;i<npts;i++) fprintf(ptout,"%13.7f %13.7f \"%d\"\n",(pts+i)->x,(pts+i)->y,i);
  fclose(ptout);

  return i;
}
//-------------------------------------------------------------------
int gmsh_pts(point *pts,int npts,char *fnm){

  int i;
  FILE *ptout;

  ptout=fopen(fnm,"w");
  for (i=0;i<npts;i++) fprintf(ptout,"Point(%d) = {%f, %f, %f};\n",i+1,(pts+i)->x,(pts+i)->y,0.0);
  fclose(ptout);

  return i;
}
//-------------------------------------------------------------------

