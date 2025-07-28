#define max_elem 25000
#define max_block 25

typedef struct { double x; double y; double z;} point;
typedef struct { int ID[4]; char cbc[4][2][4];} boundary_condition;
typedef struct { int n; point *p;} connector;
typedef struct { int elid; int esid; double curve; char ccurve;} edge;
typedef struct { int vid[4]; char BC[4][2][4];} quad;
typedef struct { int vid[8];} hex;
typedef struct { point vert[8]; int btid; int elid[2]; int nbid[8]; char BC[6][2][4];} block;
typedef struct { char name[4]; double val;} var;

//point generation operations
typedef struct {char fname[256];int npts;} list;
typedef struct {int pt1,pt2;} mid;
typedef struct {int pt1,pt2; double fra;} lin;
typedef struct {int pt1,org; double ang;} rot;
typedef struct {int pt1,pt2,org; double rad;} lci;
typedef struct {int pt1,pt2,pt3,pt4;} lli;
typedef struct {int pt1,pt2,pt3;} ref;
typedef struct {int pt1; double dxt,dyt;} trn;
typedef struct {double xc,yc;} apt;
typedef struct {double rc,tc;} rpt;

int make_quad_space(int,int,point*,char[4][2][4]);
int make_cquad_space(int,int,double,point*,char[4][2][4]);
int make_gquad_space(int,int,double,point*,char[4][2][4]);
int make_g2quad_space(int,int,double,double,point*,char[4][2][4]);
int make_cgquad_space(int,int,double,double,double,point*,char[4][2][4]);
int make_arc_space(int,int,double,point*,char[4][2][4]);
int make_garc_space(int,int,double,int,double,point*,char[4][2][4]);

int make_tri_space(int,point*,char[4][2][4]);
int make_ctri_space(int,int,double,point*,char[4][2][4]);

double distance(point,point);
double angle(point,point,point);
double growth_ratio(int,double,double);
int initialize_con(connector*,int);
int dealloc_con(connector*);
int set_lin_side(point,point,connector*);
int set_arc_side(point,point,double,connector*);
int set_garc_side(point,point,point,double,connector*);
int set_geo_side(point,point,double,connector*);
int get_g_side(point,point,double,connector*);
int get_g1D(double,double,int,double,double*);
double get_theta_0(point,point,int,double,point);
point line_circle_intercept(point,point,point,double);
point line_line_intercept(point,point,point,point);
point circle_center_2pR(point,point,double);
point rotate_point(point,double,point);
point reflect_point(point,point,point);
point translate_point(point,point);
point midpoint(point,point);
point linpoint(double,point,point);
point centroid(point*,int);
point set_point(double,double);
int invert(double*,int);
int invert_connector(connector*);

int write_rea(char*);
int axi_consistency_check(quad*);
int rotate_element(quad*);
int reset(void);
int output_pts(point*,int,char*);

int read_inp(char*);
int read_vars(FILE*,var**,int*);
int read_points(FILE*,var*,int);
int read_list(list,int*);
double get_var(char*,var*,int);
int sanitize_string(char*,int);

int make_vtri_space(int,int,double*,double*);
