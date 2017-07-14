
//structure to contain field information
struct field {
        double *data;
        int im[3];
        int ib[3];
};

int getoutput(char *key);
void write_text_2d(struct field* fp,int* im, FILE* file, int output_every, int t, int e);
void write_text_3d(struct field* fp,int* im, FILE* file, int output_every, int t, int e);
void write_binary_2d(struct field* fp,int* im, FILE* file, int output_every, int t, int e);
void write_binary_3d(struct field* fp,int* im, FILE* file, int output_every, int t, int e);
void create_field(struct field *fp, int *im, int *ib, int nr_comps);

// Everything we need to measure wallclock-time
#define BILLION 1E9

// Call this macro once at the beginning of your code with a useful profname
#define init_profs(profname) struct timespec start_time_##profname, stop_time_##profname; \
  double elapsed_time_##profname

// The time-measuring is then done via
//
//  prof_start(profname);
//  do whatever you want to measure here
//  prof_end(profname);
//
// and is printed to standard output
//

#define profs_start(profname)                                           \
  if( clock_gettime( CLOCK_REALTIME, &start_time_##profname) == -1 )            \
    {                                                           \
      perror( "clock gettime" );                                \
      exit( EXIT_FAILURE );                                     \
    }


#define profs_end(profname) if(clock_gettime( CLOCK_REALTIME, &stop_time_##profname) == -1 ) \
    {                                                                   \
      perror( "clock gettime" );                                        \
      exit( EXIT_FAILURE );                                             \
    }                                                                   \
  elapsed_time_##profname = ( stop_time_##profname.tv_sec - start_time_##profname.tv_sec ) \
    + ( stop_time_##profname.tv_nsec - start_time_##profname.tv_nsec )/BILLION; \
  printf("%s took %f seconds.\n", #profname, elapsed_time_##profname)

static int verb;

//macro to linearize 3d space
#define F3_OFF(pf, fldnr, jx,jy,jz)                                     \
  (((((fldnr                                    \
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))                              \
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))                               \
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#define F3(pf, fldnr, jx,jy,jz)         \
  (((pf)->data)[F3_OFF(pf, fldnr, jx,jy,jz)])

//macro to perform looping over whole 3d space with correct argument order
//ghost cells may be omitted via l = r = 0
//x is the fastest changing dimension (compare to F3 macro)
#define foreach_3d(ix, iy, iz, l, r) {                                  \
  int __ilo[3] = { -l, -l, -l };                \
  int __ihi[3] = { im[0] + r,                   \
                   im[1] + r,                   \
                   im[2] + r };         \
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {                        \
    for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {                      \
      for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

//macro to reduce confusion of IDEs
#define foreach_3d_end                          \
  } } }

#define foreach_2d(ix, iy, l, r){                                       \
  int __ilo[2] = { -l, -l };            \
  int __ihi[2] = { im[0] + r,                   \
                   im[1] + r };         \
  for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {                        \
    for (int ix = __ilo[0]; ix < __ihi[0]; ix++)                        \

#define foreach_2d_end                          \
  } }

//macro to perform looping over whole 3d space with correct argument order
//ghost cells are included
//x is the fastest changing dimension (compare to F3 macro)
#define foreach_3d_g(ix, iy, iz) {                                      \
  int __ilo[3] = { -ib[0], -ib[1], -ib[2] };            \
  int __ihi[3] = { im[0] + ib[0],                       \
                   im[1] + ib[1],                       \
                   im[2] + ib[2] };             \
  for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {                        \
    for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {                      \
      for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

//return true if position [x, y, z] can be found inside the pml region of the domain
#define inside_pml(x, y, z)                                             \
  (((x) < pml[0]) || ((x) >= (im[0] - pml[0])) ||                       \
   ((y) <  pml[1]) || ((y) >= (im[1] - pml[1])) ||                      \
   ((z) <  pml[2]) || ((z) >= (im[2] - pml[2])))


//output options
enum{
        ASC2D,ASC3D,BIN2D,BIN3D
} ;

typedef struct { char *key; int val; } selec;

static selec output_selection[] = {
    { "BIN3D", BIN3D }, { "bin3d", BIN3D }, { "binary3d", BIN3D },{ "BIN2D", BIN2D }, { "bin2d", BIN2D }, { "binary2d", BIN2D }, { "ASC3D", ASC3D }, { "asc3d", ASC3D }, { "ascii3d", ASC3D },{ "ASC2D", ASC2D }, { "asc2d", ASC2D }, { "ascii2d", ASC2D }
};

#define NROUTS (sizeof(output_selection)/sizeof(selec))


//field component identifiers
enum {
        EX,EY,EZ,
        HX,HY,HZ,
        NR_COMPS,
};

//possible boundary conditions
enum{
        REFLECTING,
        PERIODIC,
        ABSORBING } ;


#define sqr(x) ((x) * (x))


