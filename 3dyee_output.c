#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include "3dyee.h"

int getoutput(char *key) 
{ 
    int i; 
    for (i=0; (i < NROUTS ); i++) { 
        selec *out = output_selection + i; 
        //printf("checking %s,%s \n",out->key,key); 
        if (strcmp(out->key, key) == 0){ 
           if (verb) { printf("setting output to: %s \n",out->key);} 
            return out->val; 
        } 
    } 
    printf("no valid output type specified: options are ASC/asc/ascii for textual output and BIN/bin/binary for binary\n Setting to default (binary)."); 
    return 0; 
}

void write_binary_3d(struct field* fp,int* im, FILE* file, int output_every, int t, int e){
                if ((t%output_every) == 0) {
                        foreach_3d(ix, iy, iz, 0, 0) {
                                fwrite(&F3(fp, e, ix, iy, iz), sizeof(double), 1, file);
                        } foreach_3d_end;
                }
}

void write_text_3d(struct field* fp,int* im, FILE* file, int output_every, int t, int e){
                if ((t%output_every) == 0) {
                        foreach_3d(ix, iy, iz, 0, 0) {
                                fprintf(file,"%.9g\n",F3(fp, e, ix, iy, iz));
                        } foreach_3d_end;
                }
}

void write_binary_2d(struct field* fp,int* im, FILE* file, int output_every, int t, int e){
                if ((t%output_every) == 0) {
                        foreach_2d(ix, iy, 0, 0){
                                fwrite(&F3(fp, e, ix, iy, im[2]/2), sizeof(double), 1, file);
                        } foreach_2d_end;
                }
}

void write_text_2d(struct field* fp,int* im, FILE* file, int output_every, int t, int e){
                if ((t%output_every) == 0) {
                        foreach_2d(ix, iy, 0, 0){
                                fprintf(file,"%.9g\n",F3(fp, e, ix, iy, im[2]/2));
                        } foreach_2d_end;
                }
}
