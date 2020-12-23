#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void writegreen(FILE *fop, int order, int ntr, int npts, int ntsample,int nt_kir, int COMP, float *data, float reduce)
{

  int i, j, k, m, check;
  float *trace;
  int nreduce;

  trace = (float *)malloc(sizeof(float) * nt_kir); assert(trace != NULL);
  for (i=0; i<ntr; i ++) {
      nreduce = (int) (i * reduce);
      for(m=0; m<2; m++){
          for(j=i+(order+m)*ntr + nreduce*(COMP*ntr), k=0; (j < COMP*ntr*npts && k<nt_kir); j += (COMP*ntr)*ntsample, k++)
	      trace[k] = data[j];
          for(j=k; j<nt_kir; j++)
	      trace[j] = 0.0;    /* paded with zero */

          check = fwrite(trace,sizeof(float),nt_kir,fop);
          if (check != nt_kir) {
              fprintf(stderr,"demult2kir: write failed for trace %d\n",i);
              exit(1);
          }
      }
  }
  free(trace);
}

int main(int ac,char **av)
{
  int ntr,npts,i,j,k,check,size_f;
  float *data;
  int nx, nl, nh, nt, nz;
  float dp, dt;
  int nt_kir, kdx=1, kdt=1, ntsample;
  FILE *fop, *open_file(); 
  char kirfile[80];
  char kirfile_x[80];
  char kirfile_z[80];
  float reduce_vel = 0.0; float h, nreduce;
  int COMP = 2;
  int rcore = 1;
  int tops = 0, topp = 0; 

  size_f = sizeof(float);
  
  setpar(ac,av);
  mstpar("nl",        "d", &nl);
  mstpar("nh",        "d", &nh);
  mstpar("nx",        "d", &nx);
  mstpar("kdx",       "d", &kdx);
  mstpar("kdt",       "d", &kdt);
  mstpar("nt",        "d", &nt);
  mstpar("nt_kir",    "d", &nt_kir);
  mstpar("dp",        "f", &dp);
  mstpar("dt",        "f", &dt);
  mstpar("h",         "f", &h);
  mstpar("ntsample",  "d", &ntsample);
  getpar("rcore",     "d", &rcore);
  getpar("tops",      "d", &tops);
  getpar("topp",      "d", &topp);
  getpar("COMP",      "d", &COMP);
  if(rcore != 0 || COMP == 1)
      mstpar("fdout",     "s", kirfile);
  if(tops != 0)
      mstpar("topsfile",     "s", kirfile);
  if(topp != 0)
      mstpar("toppfile",     "s", kirfile);
  mstpar("kirfile_x", "s", kirfile_x);
  getpar("reduce_vel","f", &reduce_vel);
  if(COMP == 2)
      mstpar("kirfile_z", "s", kirfile_z);
  endpar();

  COMP *= 2;
  npts = (int) ( (float)((nt-3)*dp) /dt);
  npts /= kdt;
  ntr  = (int)((nx -nh-3 + 0.01)/kdx);     /* note nx is different from final nx in FD run*/

  fop  = open_file(kirfile,"rb"); 
  data = (float *)malloc(size_f * COMP*ntr * npts); assert(data != NULL);
  check = fread(data,size_f,COMP*ntr*npts,fop);
  if (check !=COMP* ntr*npts) {
    fprintf(stderr,"demult2kir: too little data for specified nx and nt\n");
    exit(1);
  }
  fclose(fop);

  h *= kdx;
  if(reduce_vel > 1.e-13){
      nreduce = h/reduce_vel/(dt*kdt);
  } else {
      nreduce = 0.0;
  }
  
  fop  = open_file(kirfile_x,"wb"); 
  writegreen(fop, 0,  ntr, npts, ntsample, nt_kir, COMP, data, nreduce);
  fclose(fop);
  if(COMP > 2){
      fop  = open_file(kirfile_z,"wb"); 
      writegreen(fop, 2,  ntr, npts, ntsample, nt_kir, COMP, data, nreduce);
      fclose(fop);
  }
  free(data);

  return(1);
}
