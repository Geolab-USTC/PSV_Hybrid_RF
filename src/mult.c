#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*
#include "libget.h"
  @(#)mult.c	1.2 
*/

int main(int ac,char **av)
{
  int ntr,npts,i,j,k,check,size_f;
  float *data,*trace;
  int nx, nl, nt, nz;
  int COMP = 2;
  FILE *fop, *open_file(); char greenfile[50];

  size_f = sizeof(float);
  
  setpar(ac,av);
  mstpar("nx",        "d", &nx);
  mstpar("nt",        "d", &nt);
  mstpar("nz",        "d", &nz);
  getpar("COMP",      "d", &COMP);
  getpar("greenfile", "s", greenfile);
  endpar();

  fop  = open_file(greenfile,"rb"); 
  ntr  = 2*(nx +nz+1)*COMP;
  npts = nt;
  printf("ntr=%d npts=%d\n",ntr,npts);
  data = (float *)malloc(size_f * ntr * npts); assert(data != NULL);
  check = fread(data,size_f,ntr*npts,fop);
  if (check != ntr*npts) {
    fprintf(stderr,"demult: too little data for specified nx and nt\n");
    exit(1);
  }
  fclose(fop);
  
  fop  = open_file(greenfile,"wb"); 
  trace = (float *)malloc(size_f * npts); assert(trace != NULL);
  for (i=0;i<ntr;++i) {
    for (j=i,k=0;j<ntr*npts;j+=ntr,k++) {
      trace[k] = data[j];
    }
    check = fwrite(trace,size_f,npts,fop);
    if (check != npts) {
      fprintf(stderr,"demult: write failed for trace %d\n",i);
      exit(1);
    }
  }
  fclose(fop);

  return(0);
}
