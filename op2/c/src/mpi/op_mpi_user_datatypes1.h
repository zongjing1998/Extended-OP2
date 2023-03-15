template <typename T>
void gather_data_hdf5(op_dat dat, char *usr_ptr, int low, int high);

template <typename T, void (*F)(FILE *, int, int, T *, const char *)>
void write_file(op_dat dat, const char *file_name);

template <typename T, const char *fmt>
void write_txt(FILE *fp, int g_size, int elem_size, T *g_array,
	       const char *file_name);

template <typename T>
void write_bin(FILE *fp, int g_size, int elem_size, T *g_array,
	       const char *file_name);

typedef struct
{
	double uh[2][10][4];
	double rhs[10][4];
	double wh[10][4];
	double grad[2][10][4];
	double dt;
	int pdg;
	int pfv;
	double isDiscon[4];
} uvar;
MPI_Datatype MPI_UVAR;

typedef struct
{
	int Node[2];
	int ofElem[2];
	int lsdidx[2];
	double dl;
	double norm[2];
	double etaf;
	int nqdpt;
	double flux[5][4];
} side;
MPI_Datatype MPI_SIDE;

typedef struct
{
	int Node[4];
	int Side[4];
	double xyc[2];
	double volr;
	double LocdGlb[2][2];
	double dh;
	double dxy[2];
} elem;
MPI_Datatype MPI_ELEM;

void struct_uvar()
{
	int array_of_blocklengths[8] = {80, 40, 40, 80, 1, 1, 1, 4};
	MPI_Aint array_of_displacements[8] = {0, 80 * sizeof(double), 120 * sizeof(double),
					      160 * sizeof(double), 240 * sizeof(double), 241 * sizeof(double),
					      241 * sizeof(double) + sizeof(int), 241 * sizeof(double) + 2 * sizeof(int)};
	MPI_Datatype array_of_types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
					  MPI_DOUBLE, MPI_INT, MPI_INT, MPI_DOUBLE};
	MPI_Type_create_struct(8, array_of_blocklengths, array_of_displacements, array_of_types, &MPI_UVAR);
	MPI_Type_commit(&MPI_UVAR);
}
void _mpi_allgather(uvar *l, uvar *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
	struct_uvar();
	MPI_Allgatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_UVAR, comm);
}
void _mpi_gather(uvar *l, uvar *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
	struct_uvar();
	MPI_Gatherv(l, size, MPI_UVAR, g, recevcnts, displs, MPI_UVAR, MPI_ROOT, comm);
}

void struct_side()
{
	int side_of_blocklengths[8] = {2, 2, 2, 1, 2, 1, 1, 20};
	MPI_Aint side_of_displacements[8] = {0, 2 * sizeof(int), 4 * sizeof(int),
					     6 * sizeof(int), 6 * sizeof(int) + sizeof(double), 6 * sizeof(int) + 3 * sizeof(double),
					     6 * sizeof(int) + 4 * sizeof(double), 7 * sizeof(int) + 4 * sizeof(double)};
	MPI_Datatype side_of_types[8] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE,
					 MPI_DOUBLE, MPI_INT, MPI_DOUBLE};
	MPI_Type_create_struct(8, side_of_blocklengths, side_of_displacements, side_of_types, &MPI_SIDE);
	MPI_Type_commit(&MPI_SIDE);
}
void _mpi_allgather(side *l, side *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
	struct_side();
	MPI_Allgatherv(l, size, MPI_SIDE, g, recevcnts, displs, MPI_SIDE, comm);
}
void _mpi_gather(side *l, side *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
	struct_side();
	MPI_Gatherv(l, size, MPI_SIDE, g, recevcnts, displs, MPI_SIDE, MPI_ROOT, comm);
}

void struct_elem()
{
	int elem_of_blocklengths[7] = {4, 4, 2, 1, 4, 1, 2};
	MPI_Aint elem_of_displacements[7] = {0, 4 * sizeof(int), 8 * sizeof(int), 8 * sizeof(int) + 2 * sizeof(double),
					     8 * sizeof(int) + 3 * sizeof(double), 8 * sizeof(int) + 7 * sizeof(double), 8 * sizeof(int) + 8 * sizeof(double)};
	MPI_Datatype elem_of_types[7] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Type_create_struct(7, elem_of_blocklengths, elem_of_displacements, elem_of_types, &MPI_ELEM);
	MPI_Type_commit(&MPI_ELEM);
}
void _mpi_allgather(elem *l, elem *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
	struct_elem();
	MPI_Allgatherv(l, size, MPI_ELEM, g, recevcnts, displs, MPI_ELEM, comm);
}
void _mpi_gather(elem *l, elem *g, int size, int *recevcnts, int *displs, MPI_Comm comm)
{
	struct_elem();
	MPI_Gatherv(l, size, MPI_ELEM, g, recevcnts, displs, MPI_ELEM, MPI_ROOT, comm);
}

void fetch_data_hdf5_user_datatype(op_dat dat, char *usr_ptr, int low, int high)
{
	if (strcmp(dat->type, "uvar") == 0)
		gather_data_hdf5<uvar>(dat, usr_ptr, low, high);
	else if (strcmp(dat->type, "side") == 0)
		gather_data_hdf5<side>(dat, usr_ptr, low, high);
	else if (strcmp(dat->type, "elem") == 0)
		gather_data_hdf5<elem>(dat, usr_ptr, low, high);
	else
		printf("Unknown type %s, cannot error in fetch_data_hdf5() \n", dat->type);
}

extern const char fmt_uvar[] = "%f ";
extern const char fmt_side[] = "%f ";

void print_dat_to_txtfile_mpi_user_datatypes(op_dat dat, const char *file_name)
{
	if (strcmp(dat->type, "uvar") == 0)
		write_file<uvar, write_txt<uvar, fmt_uvar>>(dat, file_name);
	else if (strcmp(dat->type, "side") == 0)
		write_file<side, write_txt<side, fmt_side>>(dat, file_name);
	else
		printf("Unknown type %s, cannot be written to file %s\n", dat->type,
		       file_name);
}

void print_dat_to_binfile_mpi_user_datatypes(op_dat dat, const char *file_name)
{
	if (strcmp(dat->type, "uvar") == 0)
		write_file<uvar, write_bin<uvar>>(dat, file_name);
	else if (strcmp(dat->type, "side") == 0)
		write_file<side, write_bin<side>>(dat, file_name);
	else
		printf("Unknown type %s, cannot be written to file %s\n", dat->type,
		       file_name);
}
