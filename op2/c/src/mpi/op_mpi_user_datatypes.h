 
	template <typename T>
	void gather_data_hdf5(op_dat dat, char *usr_ptr, int low, int high);

	template <typename T, void (*F)(FILE *, int, int, T *, const char *)>
	void write_file(op_dat dat, const char *file_name);

	template <typename T, const char *fmt>
	void write_txt(FILE *fp, int g_size, int elem_size, T *g_array,const char *file_name);

	template <typename T>
	void write_bin(FILE *fp, int g_size, int elem_size, T *g_array,const char *file_name);

	
#ifndef TOL
#define TOL 1e-18
#endif

#ifndef QUAD
#define QUAD 4
#endif

#ifndef ELEMTYPE
#define ELEMTYPE 4
#endif

#ifndef NAuxVar
#define NAuxVar 1
#endif

#ifndef NDIM
#define NDIM 2
#endif

#ifndef MAXDGORD
#define MAXDGORD 3
#endif

#ifndef MAXDGDOF
#define MAXDGDOF 10
#endif

#ifndef MAXFVORD
#define MAXFVORD 3
#endif

#ifndef MAXFVDOF
#define MAXFVDOF 10
#endif

#ifndef MAXIMDOF
#define MAXIMDOF 1
#endif

#ifndef MAXSEGQDPTS
#define MAXSEGQDPTS 5
#endif

#ifndef MAXQDORDTRI
#define MAXQDORDTRI 6
#endif

#ifndef MAXQDPTSTRI
#define MAXQDPTSTRI 12
#endif

#ifndef ADVECTI
#define ADVECTI 12
#endif

#ifndef BURGERS
#define BURGERS 23
#endif

#ifndef EULNSEQ
#define EULNSEQ 34
#endif

#ifndef EQUATIO
#define EQUATIO 34
#endif

#ifndef NEQ
#define NEQ 4
#endif

#ifndef CONSERV
#define CONSERV 123
#endif

#ifndef PRIMITI
#define PRIMITI 234
#endif

#ifndef VARIABLE
#define VARIABLE 123
#endif

#ifndef AUXVARIABLE
#define AUXVARIABLE 123
#endif

#ifndef NTM
#define NTM 2
#endif

#ifndef GLOBTIME
#define GLOBTIME 0
#endif

#ifndef LOCALEU
#define LOCALEU 1
#endif

#ifndef LOCALNS
#define LOCALNS 2
#endif

#ifndef UNSTEADY
#define UNSTEADY 0
#endif

#ifndef STEADY
#define STEADY 1
#endif

#ifndef LEFT__ELEM
#define LEFT__ELEM 1
#endif

#ifndef RIGHT_ELEM
#define RIGHT_ELEM 0
#endif

#ifndef M_E
#define M_E 2.718281828459045
#endif

#ifndef M_LOG2E
#define M_LOG2E 1.4426950408889634
#endif

#ifndef M_LOG10E
#define M_LOG10E 0.4342944819032518
#endif

#ifndef M_LN2
#define M_LN2 0.6931471805599453
#endif

#ifndef M_LN10
#define M_LN10 2.302585092994046
#endif

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

#ifndef M_PI_2
#define M_PI_2 1.5707963267948966
#endif

#ifndef M_PI_4
#define M_PI_4 0.7853981633974483
#endif

#ifndef M_1_PI
#define M_1_PI 0.3183098861837907
#endif

#ifndef M_2_PI
#define M_2_PI 0.6366197723675814
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.1283791670955126
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730951
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.7071067811865476
#endif

#ifndef FLUXNND
#define FLUXNND 1
#endif

#ifndef FLUXROE
#define FLUXROE 2
#endif

#ifndef FLUXLLF
#define FLUXLLF 3
#endif

#ifndef FLUXKIN
#define FLUXKIN 4
#endif

#ifndef FLUXOSH
#define FLUXOSH 5
#endif

#ifndef ISOTWALL
#define ISOTWALL 1
#endif

#ifndef ADIAWALL
#define ADIAWALL 2
#endif

#ifndef BC_BLKCONNECTION
#define BC_BLKCONNECTION -1
#endif

#ifndef BC_NONE
#define BC_NONE 0x0000
#endif

#ifndef BC_INTERBLK
#define BC_INTERBLK 0x0001     
#endif

#ifndef BC_SOLIDSURFACE
#define BC_SOLIDSURFACE 0x0002 
#endif

#ifndef BC_SYMMETRY
#define BC_SYMMETRY 0x0003     
#endif

#ifndef BC_FARFIELD
#define BC_FARFIELD 0x0004     
#endif

#ifndef BC_INFLOW
#define BC_INFLOW 0x0005       
#endif

#ifndef BC_OUTFLOW
#define BC_OUTFLOW 0x0006      
#endif

#ifndef BC_I_POLE
#define BC_I_POLE 000071       
#endif

#ifndef BC_J_POLE
#define BC_J_POLE 000072       
#endif

#ifndef BC_K_POLE
#define BC_K_POLE 000073       
#endif

#ifndef BC_GENERIC1
#define BC_GENERIC1 0x0008     
#endif

#ifndef BC_GENERIC2
#define BC_GENERIC2 0x0009     
#endif

#ifndef BC_GENERIC3
#define BC_GENERIC3 0x000A     
#endif

#ifndef BC_OTHERS
#define BC_OTHERS 0x0010       
#endif

#ifndef N_BCTYPE
#define N_BCTYPE 9
#endif

#ifndef TETREXTETRA
#define TETREXTETRA 1
#endif

#ifndef TETREXHEXA
#define TETREXHEXA 2
#endif

#ifndef TETREXPRISM
#define TETREXPRISM 3
#endif

#ifndef TETREXPYRAMID
#define TETREXPYRAMID 4
#endif

#ifndef TETREXQUAD
#define TETREXQUAD 5
#endif

#ifndef TETREXTRI
#define TETREXTRI 6
#endif

#ifndef MAXQDPTS1D
#define MAXQDPTS1D 5
#endif

#ifndef Entropy_Delta
#define Entropy_Delta 0.1
#endif

typedef struct{
	double uh[NTM][MAXDGDOF][NEQ]; 
	double rhs[MAXDGDOF][NEQ];     

	double wh[MAXFVDOF][NEQ];
	double grad[NDIM][MAXFVDOF][NEQ];

	
	double dt;

	unsigned int pdg; 
	unsigned int pfv; 

	
	double isDiscon[NEQ];
}uvar;
MPI_Datatype MPI_UVAR;

typedef struct{
	int Node[ELEMTYPE];
	unsigned int Side[ELEMTYPE];
	
	

	

	

	
	double xyc[NDIM];
	double volr; 

	
	
	
	
	

	double LocdGlb[NDIM][NDIM]; 
	double dh;
	double dxy[NDIM];
}elem;
MPI_Datatype MPI_ELEM;

typedef struct{
	unsigned int Node[2];
	int ofElem[2]; 

	unsigned int lsdidx[2]; 

	
	double dl;	   
	double norm[NDIM]; 

	double etaf;	    
	unsigned int nqdpt; 

	double flux[5][NEQ]; 
}side;
MPI_Datatype MPI_SIDE;

typedef struct{
	double xy[NDIM];

	

	unsigned int nNSide; 
	unsigned int nNElem; 

	unsigned int *ofSide; 
	unsigned int *ofElem; 
	double *ofElemWht;    
}node;
MPI_Datatype MPI_NODE;

typedef struct{
	unsigned int nele;
	unsigned int nsid;
	unsigned int nnod;
	unsigned int nbosid;
	unsigned int ninsid;
	unsigned int ncvbd;

	elem *pele;
	side *psid;
	node *pnod;
	uvar *pvar;
	unsigned int *pbosid;
	unsigned int *pinsid;
	
}grid;
MPI_Datatype MPI_GRID;

typedef struct{
	int ivis;      
	int fluxtype;  
	int pDG;       
	int pFV;       
	int isteady;   
	int itmax;     
	int wffre;     
	int wsfre;     
	int dtfre;     
	int iconti;    
	int resi;      
	int CFLIncStp; 
	
	int iRuKu;
	int bcTtype;	 
	int ImResSmoStp; 

	double ImResSmoCoe; 
	int ChgMaxStp;
	double M_TVB;	    
	double ChgMaxReduc; 
	double ChgMaxRatio; 
	double Re;	    
	double Ma;	    
	double RatioP;	    
	double Tinfd;	    
	double Twald;	    
	double AoA;	    
	double CFL;	    
	double tmax;	    

	char gfile[256]; 

	double a;
	double b;

#if EQUATIO == BURGERS
	
#endif

	double rhoinf, uinf, vinf, Preinf, Tinf, einf;
	double cosAoA, sinAoA;
	double Ttinf, Pretinf, Sinfn1, Rhoq2d2n1;

	double Twal;	
	double c4suth;	
	double R;	
	double Tcoe;	
	double KapdMiu; 

	
	
}para;
MPI_Datatype MPI_PARA;

typedef struct{

}cvbd2d;
MPI_Datatype MPI_CVBD2D;

void struct_uvar()
{
MPI_Datatype array_of_types[8] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT,MPI_DOUBLE};
uvar uvar_user;
MPI_Aint array_of_displacements[8];
MPI_Get_address(&uvar_user.uh, &array_of_displacements[0]);
MPI_Get_address(&uvar_user.rhs, &array_of_displacements[1]);
MPI_Get_address(&uvar_user.wh, &array_of_displacements[2]);
MPI_Get_address(&uvar_user.grad, &array_of_displacements[3]);
MPI_Get_address(&uvar_user.dt, &array_of_displacements[4]);
MPI_Get_address(&uvar_user.pdg, &array_of_displacements[5]);
MPI_Get_address(&uvar_user.pfv, &array_of_displacements[6]);
MPI_Get_address(&uvar_user.isDiscon, &array_of_displacements[7]);
array_of_displacements[8] = array_of_displacements[8] - array_of_displacements[0];
array_of_displacements[7] = array_of_displacements[7] - array_of_displacements[0];
array_of_displacements[6] = array_of_displacements[6] - array_of_displacements[0];
array_of_displacements[5] = array_of_displacements[5] - array_of_displacements[0];
array_of_displacements[4] = array_of_displacements[4] - array_of_displacements[0];
array_of_displacements[3] = array_of_displacements[3] - array_of_displacements[0];
array_of_displacements[2] = array_of_displacements[2] - array_of_displacements[0];
array_of_displacements[1] = array_of_displacements[1] - array_of_displacements[0];
array_of_displacements[0] = 0;
int array_of_blocklengths[8] = { 2*10*4, 10*4, 10*4, 2*10*4, 1, 1, 1, 4 };
MPI_Type_create_struct(8, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_UVAR);
MPI_Type_commit(&MPI_UVAR);
}

void _mpi_allgather(uvar*l, uvar *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_uvar();
MPI_Allgatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_UVAR, comm);
}

void _mpi_gather(uvar *l, uvar *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_uvar();
MPI_Gatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_UVAR,MPI_ROOT, comm);
}

void struct_elem()
{
MPI_Datatype array_of_types[7] = {MPI_INT,MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
elem elem_user;
MPI_Aint array_of_displacements[7];
MPI_Get_address(&elem_user.Node, &array_of_displacements[0]);
MPI_Get_address(&elem_user.Side, &array_of_displacements[1]);
MPI_Get_address(&elem_user.xyc, &array_of_displacements[2]);
MPI_Get_address(&elem_user.volr, &array_of_displacements[3]);
MPI_Get_address(&elem_user.LocdGlb, &array_of_displacements[4]);
MPI_Get_address(&elem_user.dh, &array_of_displacements[5]);
MPI_Get_address(&elem_user.dxy, &array_of_displacements[6]);
array_of_displacements[7] = array_of_displacements[7] - array_of_displacements[0];
array_of_displacements[6] = array_of_displacements[6] - array_of_displacements[0];
array_of_displacements[5] = array_of_displacements[5] - array_of_displacements[0];
array_of_displacements[4] = array_of_displacements[4] - array_of_displacements[0];
array_of_displacements[3] = array_of_displacements[3] - array_of_displacements[0];
array_of_displacements[2] = array_of_displacements[2] - array_of_displacements[0];
array_of_displacements[1] = array_of_displacements[1] - array_of_displacements[0];
array_of_displacements[0] = 0;
int array_of_blocklengths[7] = { 4, 4, 2, 1, 2*2, 1, 2 };
MPI_Type_create_struct(7, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_ELEM);
MPI_Type_commit(&MPI_ELEM);
}

void _mpi_allgather(elem*l, elem *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_elem();
MPI_Allgatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_ELEM, comm);
}

void _mpi_gather(elem *l, elem *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_elem();
MPI_Gatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_ELEM,MPI_ROOT, comm);
}

void struct_side()
{
MPI_Datatype array_of_types[8] = {MPI_INT,MPI_INT,MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_DOUBLE};
side side_user;
MPI_Aint array_of_displacements[8];
MPI_Get_address(&side_user.Node, &array_of_displacements[0]);
MPI_Get_address(&side_user.ofElem, &array_of_displacements[1]);
MPI_Get_address(&side_user.lsdidx, &array_of_displacements[2]);
MPI_Get_address(&side_user.dl, &array_of_displacements[3]);
MPI_Get_address(&side_user.norm, &array_of_displacements[4]);
MPI_Get_address(&side_user.etaf, &array_of_displacements[5]);
MPI_Get_address(&side_user.nqdpt, &array_of_displacements[6]);
MPI_Get_address(&side_user.flux, &array_of_displacements[7]);
array_of_displacements[8] = array_of_displacements[8] - array_of_displacements[0];
array_of_displacements[7] = array_of_displacements[7] - array_of_displacements[0];
array_of_displacements[6] = array_of_displacements[6] - array_of_displacements[0];
array_of_displacements[5] = array_of_displacements[5] - array_of_displacements[0];
array_of_displacements[4] = array_of_displacements[4] - array_of_displacements[0];
array_of_displacements[3] = array_of_displacements[3] - array_of_displacements[0];
array_of_displacements[2] = array_of_displacements[2] - array_of_displacements[0];
array_of_displacements[1] = array_of_displacements[1] - array_of_displacements[0];
array_of_displacements[0] = 0;
int array_of_blocklengths[8] = { 2, 2, 2, 1, 2, 1, 1, 5*4 };
MPI_Type_create_struct(8, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_SIDE);
MPI_Type_commit(&MPI_SIDE);
}

void _mpi_allgather(side*l, side *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_side();
MPI_Allgatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_SIDE, comm);
}

void _mpi_gather(side *l, side *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_side();
MPI_Gatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_SIDE,MPI_ROOT, comm);
}

void struct_node()
{
MPI_Datatype array_of_types[6] = {MPI_DOUBLE,MPI_INT,MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
node node_user;
MPI_Aint array_of_displacements[6];
MPI_Get_address(&node_user.xy, &array_of_displacements[0]);
MPI_Get_address(&node_user.nNSide, &array_of_displacements[1]);
MPI_Get_address(&node_user.nNElem, &array_of_displacements[2]);
MPI_Get_address(&node_user.ofSide, &array_of_displacements[3]);
MPI_Get_address(&node_user.ofElem, &array_of_displacements[4]);
MPI_Get_address(&node_user.ofElemWht, &array_of_displacements[5]);
array_of_displacements[6] = array_of_displacements[6] - array_of_displacements[0];
array_of_displacements[5] = array_of_displacements[5] - array_of_displacements[0];
array_of_displacements[4] = array_of_displacements[4] - array_of_displacements[0];
array_of_displacements[3] = array_of_displacements[3] - array_of_displacements[0];
array_of_displacements[2] = array_of_displacements[2] - array_of_displacements[0];
array_of_displacements[1] = array_of_displacements[1] - array_of_displacements[0];
array_of_displacements[0] = 0;
int array_of_blocklengths[6] = { 2, 1, 1, 1, 1, 1 };
MPI_Type_create_struct(6, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_NODE);
MPI_Type_commit(&MPI_NODE);
}

void _mpi_allgather(node*l, node *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_node();
MPI_Allgatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_NODE, comm);
}

void _mpi_gather(node *l, node *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_node();
MPI_Gatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_NODE,MPI_ROOT, comm);
}

void struct_grid()
{
MPI_Datatype array_of_types[12] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
grid grid_user;
MPI_Aint array_of_displacements[12];
MPI_Get_address(&grid_user.nele, &array_of_displacements[0]);
MPI_Get_address(&grid_user.nsid, &array_of_displacements[1]);
MPI_Get_address(&grid_user.nnod, &array_of_displacements[2]);
MPI_Get_address(&grid_user.nbosid, &array_of_displacements[3]);
MPI_Get_address(&grid_user.ninsid, &array_of_displacements[4]);
MPI_Get_address(&grid_user.ncvbd, &array_of_displacements[5]);
MPI_Get_address(&grid_user.pele, &array_of_displacements[6]);
MPI_Get_address(&grid_user.psid, &array_of_displacements[7]);
MPI_Get_address(&grid_user.pnod, &array_of_displacements[8]);
MPI_Get_address(&grid_user.pvar, &array_of_displacements[9]);
MPI_Get_address(&grid_user.pbosid, &array_of_displacements[10]);
MPI_Get_address(&grid_user.pinsid, &array_of_displacements[11]);
array_of_displacements[12] = array_of_displacements[12] - array_of_displacements[0];
array_of_displacements[11] = array_of_displacements[11] - array_of_displacements[0];
array_of_displacements[10] = array_of_displacements[10] - array_of_displacements[0];
array_of_displacements[9] = array_of_displacements[9] - array_of_displacements[0];
array_of_displacements[8] = array_of_displacements[8] - array_of_displacements[0];
array_of_displacements[7] = array_of_displacements[7] - array_of_displacements[0];
array_of_displacements[6] = array_of_displacements[6] - array_of_displacements[0];
array_of_displacements[5] = array_of_displacements[5] - array_of_displacements[0];
array_of_displacements[4] = array_of_displacements[4] - array_of_displacements[0];
array_of_displacements[3] = array_of_displacements[3] - array_of_displacements[0];
array_of_displacements[2] = array_of_displacements[2] - array_of_displacements[0];
array_of_displacements[1] = array_of_displacements[1] - array_of_displacements[0];
array_of_displacements[0] = 0;
int array_of_blocklengths[12] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
MPI_Type_create_struct(12, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_GRID);
MPI_Type_commit(&MPI_GRID);
}

void _mpi_allgather(grid*l, grid *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_grid();
MPI_Allgatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_GRID, comm);
}

void _mpi_gather(grid *l, grid *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_grid();
MPI_Gatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_GRID,MPI_ROOT, comm);
}

void struct_para()
{
MPI_Datatype array_of_types[39] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_DOUBLE,MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_CHAR,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
para para_user;
MPI_Aint array_of_displacements[39];
MPI_Get_address(&para_user.ivis, &array_of_displacements[0]);
MPI_Get_address(&para_user.fluxtype, &array_of_displacements[1]);
MPI_Get_address(&para_user.pDG, &array_of_displacements[2]);
MPI_Get_address(&para_user.pFV, &array_of_displacements[3]);
MPI_Get_address(&para_user.isteady, &array_of_displacements[4]);
MPI_Get_address(&para_user.itmax, &array_of_displacements[5]);
MPI_Get_address(&para_user.wffre, &array_of_displacements[6]);
MPI_Get_address(&para_user.wsfre, &array_of_displacements[7]);
MPI_Get_address(&para_user.dtfre, &array_of_displacements[8]);
MPI_Get_address(&para_user.iconti, &array_of_displacements[9]);
MPI_Get_address(&para_user.resi, &array_of_displacements[10]);
MPI_Get_address(&para_user.CFLIncStp, &array_of_displacements[11]);
MPI_Get_address(&para_user.iRuKu, &array_of_displacements[12]);
MPI_Get_address(&para_user.bcTtype, &array_of_displacements[13]);
MPI_Get_address(&para_user.ImResSmoStp, &array_of_displacements[14]);
MPI_Get_address(&para_user.ImResSmoCoe, &array_of_displacements[15]);
MPI_Get_address(&para_user.ChgMaxStp, &array_of_displacements[16]);
MPI_Get_address(&para_user.M_TVB, &array_of_displacements[17]);
MPI_Get_address(&para_user.ChgMaxReduc, &array_of_displacements[18]);
MPI_Get_address(&para_user.ChgMaxRatio, &array_of_displacements[19]);
MPI_Get_address(&para_user.Re, &array_of_displacements[20]);
MPI_Get_address(&para_user.Ma, &array_of_displacements[21]);
MPI_Get_address(&para_user.RatioP, &array_of_displacements[22]);
MPI_Get_address(&para_user.Tinfd, &array_of_displacements[23]);
MPI_Get_address(&para_user.Twald, &array_of_displacements[24]);
MPI_Get_address(&para_user.AoA, &array_of_displacements[25]);
MPI_Get_address(&para_user.CFL, &array_of_displacements[26]);
MPI_Get_address(&para_user.tmax, &array_of_displacements[27]);
MPI_Get_address(&para_user.gfile, &array_of_displacements[28]);
MPI_Get_address(&para_user.a, &array_of_displacements[29]);
MPI_Get_address(&para_user.b, &array_of_displacements[30]);
MPI_Get_address(&para_user.rhoinf, &array_of_displacements[31]);
MPI_Get_address(&para_user.cosAoA, &array_of_displacements[32]);
MPI_Get_address(&para_user.Ttinf, &array_of_displacements[33]);
MPI_Get_address(&para_user.Twal, &array_of_displacements[34]);
MPI_Get_address(&para_user.c4suth, &array_of_displacements[35]);
MPI_Get_address(&para_user.R, &array_of_displacements[36]);
MPI_Get_address(&para_user.Tcoe, &array_of_displacements[37]);
MPI_Get_address(&para_user.KapdMiu, &array_of_displacements[38]);
array_of_displacements[39] = array_of_displacements[39] - array_of_displacements[0];
array_of_displacements[38] = array_of_displacements[38] - array_of_displacements[0];
array_of_displacements[37] = array_of_displacements[37] - array_of_displacements[0];
array_of_displacements[36] = array_of_displacements[36] - array_of_displacements[0];
array_of_displacements[35] = array_of_displacements[35] - array_of_displacements[0];
array_of_displacements[34] = array_of_displacements[34] - array_of_displacements[0];
array_of_displacements[33] = array_of_displacements[33] - array_of_displacements[0];
array_of_displacements[32] = array_of_displacements[32] - array_of_displacements[0];
array_of_displacements[31] = array_of_displacements[31] - array_of_displacements[0];
array_of_displacements[30] = array_of_displacements[30] - array_of_displacements[0];
array_of_displacements[29] = array_of_displacements[29] - array_of_displacements[0];
array_of_displacements[28] = array_of_displacements[28] - array_of_displacements[0];
array_of_displacements[27] = array_of_displacements[27] - array_of_displacements[0];
array_of_displacements[26] = array_of_displacements[26] - array_of_displacements[0];
array_of_displacements[25] = array_of_displacements[25] - array_of_displacements[0];
array_of_displacements[24] = array_of_displacements[24] - array_of_displacements[0];
array_of_displacements[23] = array_of_displacements[23] - array_of_displacements[0];
array_of_displacements[22] = array_of_displacements[22] - array_of_displacements[0];
array_of_displacements[21] = array_of_displacements[21] - array_of_displacements[0];
array_of_displacements[20] = array_of_displacements[20] - array_of_displacements[0];
array_of_displacements[19] = array_of_displacements[19] - array_of_displacements[0];
array_of_displacements[18] = array_of_displacements[18] - array_of_displacements[0];
array_of_displacements[17] = array_of_displacements[17] - array_of_displacements[0];
array_of_displacements[16] = array_of_displacements[16] - array_of_displacements[0];
array_of_displacements[15] = array_of_displacements[15] - array_of_displacements[0];
array_of_displacements[14] = array_of_displacements[14] - array_of_displacements[0];
array_of_displacements[13] = array_of_displacements[13] - array_of_displacements[0];
array_of_displacements[12] = array_of_displacements[12] - array_of_displacements[0];
array_of_displacements[11] = array_of_displacements[11] - array_of_displacements[0];
array_of_displacements[10] = array_of_displacements[10] - array_of_displacements[0];
array_of_displacements[9] = array_of_displacements[9] - array_of_displacements[0];
array_of_displacements[8] = array_of_displacements[8] - array_of_displacements[0];
array_of_displacements[7] = array_of_displacements[7] - array_of_displacements[0];
array_of_displacements[6] = array_of_displacements[6] - array_of_displacements[0];
array_of_displacements[5] = array_of_displacements[5] - array_of_displacements[0];
array_of_displacements[4] = array_of_displacements[4] - array_of_displacements[0];
array_of_displacements[3] = array_of_displacements[3] - array_of_displacements[0];
array_of_displacements[2] = array_of_displacements[2] - array_of_displacements[0];
array_of_displacements[1] = array_of_displacements[1] - array_of_displacements[0];
array_of_displacements[0] = 0;
int array_of_blocklengths[39] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 256, 1, 1,,,, 1, 1, 1, 1, 1 };
MPI_Type_create_struct(39, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_PARA);
MPI_Type_commit(&MPI_PARA);
}

void _mpi_allgather(para*l, para *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_para();
MPI_Allgatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_PARA, comm);
}

void _mpi_gather(para *l, para *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
struct_para();
MPI_Gatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_PARA,MPI_ROOT, comm);
}

void fetch_data_hdf5_user_datatype(op_dat dat, char *usr_ptr, int low, int high)
{
if (strcmp(dat->type,"uvar") == 0)
gather_data_hdf5<uvar>(dat, usr_ptr, low, high);
else if (strcmp(dat->type,"elem") == 0)
gather_data_hdf5<elem>(dat, usr_ptr, low, high);
else if (strcmp(dat->type,"side") == 0)
gather_data_hdf5<side>(dat, usr_ptr, low, high);
else if (strcmp(dat->type,"node") == 0)
gather_data_hdf5<node>(dat, usr_ptr, low, high);
else if (strcmp(dat->type,"grid") == 0)
gather_data_hdf5<grid>(dat, usr_ptr, low, high);
else if (strcmp(dat->type,"para") == 0)
gather_data_hdf5<para>(dat, usr_ptr, low, high);
else 
printf("Unknown type %s, cannot error in fetch_data_hdf5() \n", dat->type);
}

extern const char fmt_struct[] = "%f ";


void print_dat_to_txtfile_mpi_user_datatypes(op_dat dat, const char *file_name)
{
if (strcmp(dat->type,"uvar") == 0)
write_file<uvar, write_txt<uvar, fmt_struct>>(dat, file_name);
else if (strcmp(dat->type,"elem") == 0)
write_file<elem, write_txt<elem, fmt_struct>>(dat, file_name);
else if (strcmp(dat->type,"side") == 0)
write_file<side, write_txt<side, fmt_struct>>(dat, file_name);
else if (strcmp(dat->type,"node") == 0)
write_file<node, write_txt<node, fmt_struct>>(dat, file_name);
else if (strcmp(dat->type,"grid") == 0)
write_file<grid, write_txt<grid, fmt_struct>>(dat, file_name);
else if (strcmp(dat->type,"para") == 0)
write_file<para, write_txt<para, fmt_struct>>(dat, file_name);
else 
printf("Unknown type %s, cannot be written to file %s\n", dat->type,file_name);
}


void print_dat_to_binfile_mpi_user_datatypes(op_dat dat, const char *file_name)
{
if (strcmp(dat->type,"uvar") == 0)
write_file<uvar, write_bin<uvar>>(dat, file_name);
else if (strcmp(dat->type,"elem") == 0)
write_file<elem, write_bin<elem>>(dat, file_name);
else if (strcmp(dat->type,"side") == 0)
write_file<side, write_bin<side>>(dat, file_name);
else if (strcmp(dat->type,"node") == 0)
write_file<node, write_bin<node>>(dat, file_name);
else if (strcmp(dat->type,"grid") == 0)
write_file<grid, write_bin<grid>>(dat, file_name);
else if (strcmp(dat->type,"para") == 0)
write_file<para, write_bin<para>>(dat, file_name);
else 
printf("Unknown type %s, cannot be written to file %s\n", dat->type,file_name);
}

