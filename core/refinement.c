int make_pframe_space(int Nin,int bin,int c,point *p,char bcs[4][2][4]){

  int a;
  double fact1,fact2;
  int i,j,k,nvert0=nvert;
  connector s2,s4,rr;
  point p4,p5,p6,p7,p8,p9;

  if(((N-b)%2)>0){
    printf("Error: invalid ratio in picture frame: %d to %d\n",N,b);
    return 0;
  }

  a=(N-b)/2;

  printf("Generating %d elements in a quadrilateral picture frame space\n",(2*a+b)*c+(a*b));

//      3--a--5--b--4--a--2
//      |     |     |     |
//      |     c     c     |
//      |     |     |     |
//      c     6--b--7     c 
//      |    /       \    |
//      |   /         \   | 
//      |  a           a  |
// j,s  | /             \ |
//  ᴧ   |/               \|
//  |   0--------b--------1
//  |
//  *---->i,r

  fact1=(double)a/(double)(2*a+b);
  p4=linpoint(fact1,p[2],p[3]);
  p5=linpoint(fact1,p[3],p[2]);
  fact2=(double)c/(double)(c+a)
  p6=linpoint(fact2,p[3],p[0]);
  p7=linpoint(fact2,p[2],p[1]);
  p8=linpoint(fact1,p6,p7);
  p9=linpoint(fact1,p7,p6);

//make the 4 elements
  set_4_bcs()

//take points 0 - 5, generate 6 & 7

//take points 0 - 3, generate 4 - 7

