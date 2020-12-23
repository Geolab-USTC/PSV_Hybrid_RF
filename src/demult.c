#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*
#include "libget.h"
  @(#)demult.c	1.2 
*/

int main(int ac,char **av)
{
  int ntr,npts,i,j,k,check,size_f;
  float *data,*trace;
  int nx, nl, nt, nz;
  int nl_skip = 0, kdx = 1;
  FILE *fop, *open_file(), *fop1; char greenfile[50], greenfile1[50];
  int COMP = 2;
  int diffile = 0;
  int finish;
  long int offset, rec_len;
  int read_record();

  size_f = sizeof(float);
  
  setpar(ac,av);
  mstpar("nx",        "d", &nx);
  getpar("nlskip",    "d", &nl_skip);
  getpar("kdx",       "d", &kdx);
  mstpar("nt",        "d", &nt);
  getpar("greenfile", "s", greenfile);
  getpar("diffile",   "d", &diffile);
  getpar("COMP",      "d", &COMP);
  if(diffile != 0)
      mstpar("greenfile1","s", greenfile1);
  endpar();

  COMP *= 2;

  fop  = open_file(greenfile,"rb"); 
  fread(&nz,sizeof(int),1,fop);
//       fprintf(stderr,"nz=%d\n",nz);
       printf("nz=%d\n",nz+1);

  nl_skip = ((int) ((nl_skip+0.01)/kdx))*kdx;
  nx += nl_skip;

  ntr  = COMP*(nx+nz + 1);
  npts = nt;
  printf("ntr = %d, npts = %d\n",ntr,npts);    /* Ling: 01/31/04 */

/*
  if( (data = (float *)malloc(size_f * ntr * npts)) != NULL || diffile == 1){
*/
  if( diffile != 1){
      if( (data = (float *)malloc(size_f * ntr * npts)) == NULL){
          fprintf(stderr,"annot allocate memory for data\n");
          exit(-1);
      }
      check = fread(data,size_f,ntr*npts,fop);
      if (check != ntr*npts) {
        fprintf(stderr,"demult: too little data for specified nx and nt\n");
        exit(1);
      }
      fclose(fop);
  
      fop  = open_file(greenfile,"wb"); 
      trace = (float *)malloc(size_f * ntr); assert(trace != NULL);
      for (i=0;i<npts;++i) {
        for (j=i,k=0;j<ntr*npts;j+=npts,k++) {
          trace[k] = data[j];
        }
        check = fwrite(trace,size_f,ntr,fop);
        if (check != ntr) {
          fprintf(stderr,"demult: write failed for trace %d\n",i);
          exit(1);
        }
      }
      fclose(fop);
      free(data);
      finish = 0;
   } else {
       fop1 = open_file(greenfile1,"wb");
       fprintf(stderr,"ntr=%d\n",ntr);
       trace = (float *)malloc(size_f * ntr); assert(trace != NULL);
       rec_len = npts * size_f;
       for (i=0;i<npts;++i) {
           offset = i * size_f +sizeof(int);
           read_record(fop,offset,trace,ntr,size_f,rec_len);

           check = fwrite(trace,size_f,ntr,fop1);
           if (check != ntr) {
               fprintf(stderr,"demult: write failed for trace %d\n",i);
               exit(1);
           }
        }
        fclose(fop1);
        fclose(fop);
        finish = 1;
   }

  return finish;
}
