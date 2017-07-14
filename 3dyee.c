#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include "3dyee.h"


//function to create memory structure for field information
void create_field(struct field *fp, int *im, int *ib, int nr_comps)
{
        int size = nr_comps;
        for (int i = 0 ;i < 3; i++)
        {
                fp->im[i] = im[i] + 2 * ib[i];
                fp->ib[i] = -ib[i];
                size *= fp->im[i];
        }
        fp->data = (double*) calloc(size, sizeof(*fp->data));
}

//function to destroy memory structure for field information
void free_field(struct field *fp)
{
        free(fp->data);
}


//value of a 3d gaussian (with maximum at xm, width dxm and frequency f, 
//flying in x direction)
//at position [xx, yy, zz]
//needed for initialization
double gauss_3d(double xx, double yy, double zz, double *xm, double *dxm, double f)
{
	double xr = xx - xm[0];
	double yr = yy - xm[1];
	double zr = zz - xm[2];
	return (exp(-sqr(xr/dxm[0])) *
	       exp(-sqr(yr/dxm[1])) *
	       exp(-sqr(zr/dxm[2]))) * sin(f * 2. * M_PI * (xr));
}

//periodic boundaries for magnetic fields (one ghost cell)
void exchange_hfields(struct field *fp, int *im)
{
	for (int iz = 0; iz < im[2]; iz++) {
		for (int iy = 0; iy < im[1]; iy++) {
			F3(fp, HX, -1, iy, iz) = F3(fp, HX, im[0] - 1, iy, iz);
			F3(fp, HY, -1, iy, iz) = F3(fp, HY, im[0] - 1, iy, iz);
			F3(fp, HZ, -1, iy, iz) = F3(fp, HZ, im[0] - 1, iy, iz);
	
			F3(fp, HX, im[0], iy, iz) = F3(fp, HX, 0, iy, iz);
			F3(fp, HY, im[0], iy, iz) = F3(fp, HY, 0, iy, iz);
			F3(fp, HZ, im[0], iy, iz) = F3(fp, HZ, 0, iy, iz);
		}
	}
	for (int iz = 0; iz < im[2]; iz++) {
		for (int ix = 0; ix < im[0]; ix++) {
			F3(fp, HX, ix, -1, iz) = F3(fp, HX, ix, im[1] - 1, iz);
			F3(fp, HY, ix, -1, iz) = F3(fp, HY, ix, im[1] - 1, iz);
			F3(fp, HZ, ix, -1, iz) = F3(fp, HZ, ix, im[1] - 1, iz);
	
			F3(fp, HX, ix, im[1], iz) = F3(fp, HX, ix, 0, iz);
			F3(fp, HY, ix, im[1], iz) = F3(fp, HY, ix, 0, iz);
			F3(fp, HZ, ix, im[1], iz) = F3(fp, HZ, ix, 0, iz);
		}
	}
	
	for (int iy = 0; iy < im[1]; iy++) {
		for (int ix = 0; ix < im[0]; ix++) {
			F3(fp, HX, ix, iy, -1) = F3(fp, HX, ix, iy, im[2] - 1);
			F3(fp, HY, ix, iy, -1) = F3(fp, HY, ix, iy, im[2] - 1);
			F3(fp, HZ, ix, iy, -1) = F3(fp, HZ, ix, iy, im[2] - 1);
	
			F3(fp, HX, ix, iy, im[2]) = F3(fp, HX, ix, iy, 0);
			F3(fp, HY, ix, iy, im[2]) = F3(fp, HY, ix, iy, 0);
			F3(fp, HZ, ix, iy, im[2]) = F3(fp, HZ, ix, iy, 0);
		}
	}
}

//periodic boundaries for electric fields (one ghost cell)
void exchange_efields(struct field *fp, int *im)
{
	for (int iz = 0; iz < im[2]; iz++) {
		for (int iy = 0; iy < im[1]; iy++) {
			F3(fp, EX, -1, iy, iz) = F3(fp, EX, im[0] - 1, iy, iz);
			F3(fp, EY, -1, iy, iz) = F3(fp, EY, im[0] - 1, iy, iz);
			F3(fp, EZ, -1, iy, iz) = F3(fp, EZ, im[0] - 1, iy, iz);

			F3(fp, EX, im[0], iy, iz) = F3(fp, EX, 0, iy, iz);
			F3(fp, EY, im[0], iy, iz) = F3(fp, EY, 0, iy, iz);
			F3(fp, EZ, im[0], iy, iz) = F3(fp, EZ, 0, iy, iz);
		}
	}
	for (int iz = 0; iz < im[2]; iz++) {
		for (int ix = 0; ix < im[0]; ix++) {
			F3(fp, EX, ix, -1, iz) = F3(fp, EX, ix, im[1] - 1, iz);
			F3(fp, EY, ix, -1, iz) = F3(fp, EY, ix, im[1] - 1, iz);
			F3(fp, EZ, ix, -1, iz) = F3(fp, EZ, ix, im[1] - 1, iz);

			F3(fp, EX, ix, im[1], iz) = F3(fp, EX, ix, 0, iz);
			F3(fp, EY, ix, im[1], iz) = F3(fp, EY, ix, 0, iz);
			F3(fp, EZ, ix, im[1], iz) = F3(fp, EZ, ix, 0, iz);
		}
	}

	for (int iy = 0; iy < im[1]; iy++) {
		for (int ix = 0; ix < im[0]; ix++) {
			F3(fp, EX, ix, iy, -1) = F3(fp, EX, ix, iy, im[2] - 1);
			F3(fp, EY, ix, iy, -1) = F3(fp, EY, ix, iy, im[2] - 1);
			F3(fp, EZ, ix, iy, -1) = F3(fp, EZ, ix, iy, im[2] - 1);

			F3(fp, EX, ix, iy, im[2]) = F3(fp, EX, ix, iy, 0);
			F3(fp, EY, ix, iy, im[2]) = F3(fp, EY, ix, iy, 0);
			F3(fp, EZ, ix, iy, im[2]) = F3(fp, EZ, ix, iy, 0);
		}
	}
}

int main(int argc, char **argv)
{
	
        verb=0;
	double length[3] = {1., 1., 1.};	//some normalized box size
	int im[3] = {50,50,50};				//number of grid-points
	int ib[3] = {1,1,1};				//number of ghost cells (must be 1)
	int pml[3] = {7,7,7};                   //seven cells for the pmls at each boundary	
	int tmax = 100;						//number of desired timesteps
 	double sigma_j[3] = {0., 0., 0.};	//constant conductivity, here we take the vacuum
	double sigma_pml[3] = {30., 30., 30.};	//constant conductivity for PMLs
	double *sigma = sigma_j; // this is a pointer to a double array which is set to sigma_j or sigma_pml depending on whether we are in the PMLs or not
	
	int output_type=BIN2D;   
	int bound = REFLECTING;
	int output_every = 10;				//intervall between field outputs
	void (*write_output)(struct field* fp,int* im, FILE* file, int output_every, int t, int e);

	char *name="out";       //name of the output file, set by flag -n
	char *output="BIN";       //type of the output, set by flag -o
	char *grid = NULL;
        char *cond = NULL;
        char *time = NULL;
 	int n = 0;  		//number of cells to be passed on cli
 	double sig = 0;  		//conductivity to be passed on cli
        int opt;

        printf("parse parameters: \n \
 give	-g <int> to set grid cell number \n \
	-t <int> to set number of timesteps \n \
	-s <int> to set conductivity in box \n \
	-o <string> to set type of output \n \
	-n <string> to set name for outputfile \n \
	-p to set periodic boundaries \n \
	-r to set reflecting boundaries \n \
	-a to set absorbing boundaries (pmls) \n \
	-v to set to verbose (profiling output)\n");

        while ((opt = getopt(argc,argv,"arpvg:t:s::n:o:")) != -1) {
        switch (opt){
                case 'v': verb=1; if (verb==1) { printf("set to verbose \n");} break;
                case 'r': bound=REFLECTING;if (verb) { printf("set boundary to reflecting");}  break;
                case 'p': bound=PERIODIC;if (verb) { printf("set boundary to periodic");}  break;
                case 'a': bound=ABSORBING;if (verb) { printf("set boundary to absorbing");} break;
                case 'g': grid=optarg;if (verb) { printf("n = %s\n",grid);}  break;
                case 't': time=optarg;if (verb) { printf("tmax = %s\n",time);}  break;
                case 's': cond=optarg;if (verb) { printf("sigma_j = %s\n",cond);}  break;
                case 'n': name=optarg;if (verb) { printf("name = %s\n",name);}  break;
                case 'o': output=optarg;if (verb) { printf("output_type = %s\n",output);}  break;
                default:
                        fprintf(stderr, "use [-rpgts] [file...]\n",argv[0]);
                }
        }

	
        //convert to integer/double
        if(grid) n=strtol(grid,NULL,0);
        if(time) tmax=strtol(time,NULL,0) ;
        if(cond) sig=strtod(cond,NULL) ;
	
	if(n>0) {
		im[0]=n; im[1]=n; im[2]=n;
	}
	
	if (verb) {printf("grid %i \n",im[0]); }
	if(sig>0) {
		if (verb) { printf("setting sigma_j to %f \n",sig); }
		sigma_j[0]=sig; sigma_j[1]=sig; sigma_j[2]=sig;
	}

	output_type=getoutput(output);
	
	switch(output_type){
        case ASC2D: {write_output = &write_text_2d; if (verb) { printf("output number %i \n",output_type); } break;}
        case BIN2D: {write_output = &write_binary_2d; if (verb) { printf("output number %i \n",output_type); }  break;}
	case ASC3D: {write_output = &write_text_3d; if (verb) { printf("output number %i \n",output_type); } break;}
        case BIN3D: {write_output = &write_binary_3d; if (verb) { printf("output number %i \n",output_type); } break;}
        default:if (verb) { printf("no valid output type specified: options are ASC2D/ASC3D for textual output and BIN2D/BIN3D for binary\n");}
        }
	
        printf("set parameters to %d,%d,%f \n", im[0],tmax,sigma_j[0] );
	
	char filename[50];
        sprintf(filename, "%s",name);
	FILE *file;
	file = fopen (filename,"wb");
	//write timesteps and grid size to file 
	
	if ((output_type==ASC2D) || (output_type==ASC3D)){
		fprintf(file,"%i\n", tmax);
		(output_type==ASC3D)? fprintf(file,"%i \n", im[0]):fprintf(file,"%i \n", 1);
		fprintf(file,"%i \n", im[1]);
		fprintf(file,"%i \n", im[2]);
		fprintf(file,"%i\n",output_every);
	}

	if ( (output_type==BIN2D) || (output_type==BIN3D)){
		int one;
		one=1;
		fwrite(&tmax, sizeof(tmax), 1, file);
		(output_type==ASC3D)? fwrite(&im[0], sizeof(int), 1, file):fwrite(&one, sizeof(int), 1, file);
		fwrite(&im[1], sizeof(int), 1, file);
		fwrite(&im[2], sizeof(int), 1, file);
		fwrite(&output_every, sizeof(int), 1, file);
	}

	//field structure creation
	struct field f;
	struct field *fp = &f;
	create_field(fp, im, ib, NR_COMPS);
	
	//spatial and temporal step size calculation in 3d
	double dx[3];
	for (int i = 0; i< 3; i++) {
		dx[i] = length[i] / im[i];
	}
	double inv_sum = 0.;
	for (int d=0;d<3;d++) {
		if (im[d] > 1) {
			inv_sum += 1. / (dx[d] * dx[d]);
		}
	}
	assert(inv_sum);  //Simulation has 0 dimensions
	double dt = .75 * sqrt(1./inv_sum);

	double x0[3] = {.2, .8, .5};		//center point for gaussian pulse
	double sigma_gauss[3] = {.1, .1, .1};		//standard deviation for gaussian pulse
	double freq = 10.;					//frequency of gaussian pulse

	//initialization of field data
	foreach_3d_g (ix, iy, iz) {
		F3(fp, EZ, ix, iy, iz) = gauss_3d(ix * dx[0], iy * dx[1], (iz+.5) * dx[2], x0, sigma_gauss, freq);
		F3(fp, HY, ix, iy, iz) = -gauss_3d((ix + .5) * dx[0], iy * dx[1], (iz +.5)* dx[2], x0, sigma_gauss, freq);
	} foreach_3d_end;
	double cnx = dt / dx[0];
	double cny = dt / dx[1];
	double cnz = dt / dx[2];

	init_profs(fdtd_integration_loop);
        init_profs(e_field_update);
        init_profs(h_field_update);
        init_profs(output);

        profs_start(fdtd_integration_loop);

	//integration loop
	for(int t = 0; t < tmax; t++)
	{
		if (bound == REFLECTING) {exchange_hfields(fp, im);}	//periodic boundary conditions
		//e field update
		if ((verb) && ((t%output_every) == 0)) profs_start(e_field_update);
		  
		foreach_3d(ix, iy, iz, 0, 1) {
		  if(bound==ABSORBING) {
		    if (inside_pml(ix, iy, iz)) sigma = sigma_pml;
		    else sigma = sigma_j;
		  }
		  
		  F3(fp, EX, ix, iy, iz) =
		    (cny * (F3(fp, HZ, ix,iy,iz) - F3(fp, HZ, ix,iy-1,iz)) -
		     cnz * (F3(fp, HY, ix,iy,iz) - F3(fp, HY, ix,iy,iz-1)) +
		     F3(fp, EX, ix, iy, iz) * (1. - sigma[0]/2 * dt)) / (1. + sigma[0]/2 * dt);
			F3(fp, EY, ix, iy, iz) =
				(cnz * (F3(fp, HX, ix,iy,iz) - F3(fp, HX, ix,iy,iz-1)) -
				 cnx * (F3(fp, HZ, ix,iy,iz) - F3(fp, HZ, ix-1,iy,iz)) +
				 F3(fp, EY, ix, iy, iz) * (1. - sigma[1]/2 * dt)) / (1. + sigma[1]/2 * dt);
			F3(fp, EZ, ix, iy, iz) =
				(cnx * (F3(fp, HY, ix,iy,iz) - F3(fp, HY, ix-1,iy,iz)) -
				 cny * (F3(fp, HX, ix,iy,iz) - F3(fp, HX, ix,iy-1,iz)) +
				 F3(fp, EZ, ix, iy, iz) * (1. - sigma[2]/2 * dt)) / (1. + sigma[2]/2 * dt);
		} foreach_3d_end;
		
		if ((verb) && ((t%output_every) == 0)) { profs_end(e_field_update);}
		if (bound == REFLECTING) {exchange_efields(fp, im);}	//periodic boundary conditions

		//h field update

		if ((verb) && ((t%output_every) == 0)) profs_start(h_field_update);
		
		foreach_3d(ix, iy, iz, 1, 0) {
		  if(bound==ABSORBING) {
		    if (inside_pml(ix, iy, iz)) sigma = sigma_pml;
		    else sigma = sigma_j;
		  }

		  F3(fp, HX, ix,iy,iz) =
		    (cny * (F3(fp, EZ, ix,iy,iz) - F3(fp, EZ, ix,iy+1,iz)) -
		     cnz * (F3(fp, EY, ix,iy,iz) - F3(fp, EY, ix,iy,iz+1)) + 
		     F3(fp, HX, ix,iy,iz) * (1. - sigma[0]/2 * dt)) / (1. + sigma[0]/2 * dt);
		  
		  F3(fp, HY, ix,iy,iz) =
		    (cnz * (F3(fp, EX, ix,iy,iz) - F3(fp, EX, ix,iy,iz+1)) -
		     cnx * (F3(fp, EZ, ix,iy,iz) - F3(fp, EZ, ix+1,iy,iz)) +
		     F3(fp, HY, ix,iy,iz) * (1. - sigma[1]/2 * dt)) / (1. + sigma[1]/2 * dt);
		  
		  F3(fp, HZ, ix,iy,iz) =
		    (cnx * (F3(fp, EY, ix,iy,iz) - F3(fp, EY, ix+1,iy,iz)) -
		     cny * (F3(fp, EX, ix,iy,iz) - F3(fp, EX, ix,iy+1,iz)) +
		     F3(fp, HZ, ix,iy,iz) * (1. - sigma[2]/2 * dt)) / (1. + sigma[2]/2 * dt);
		} foreach_3d_end;
		
		if ((verb) && ((t%output_every) == 0)) { profs_end(h_field_update);}

		if ((verb) && ((t%output_every) == 0)) {profs_start(output);}
		(*write_output)(fp,im,file,output_every,t,EZ);		
        	if ((verb) && ((t%output_every) == 0)) { profs_end(output); }
	
	}
        profs_end(fdtd_integration_loop);


	fclose(file);
	free_field(fp);
}
