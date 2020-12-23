
/*********************************************************************/
/*                                                                   */
/*          Hybrid code for P-SV system                              */
/*                                                                   */
/*  Written by: Lianxing Wen, Seismological Lab, Caltech             */
/*  Last modified: Tuesday, Oct. 31, 1995.                           */
/*                                                                   */
/*********************************************************************/

/* modified history:
    Lianxing Wen          8/22/1997 add output options to the top
    Lianxing Wen          11/2/1997 add options to skip the first 
			  several output points-> nl_skip
    Lianxing Wen          7/4/1999 add options for FD model inputs 
*/

#include   <stdio.h>
#include   <stdlib.h>
#include   <sys/file.h>
#include   <fcntl.h>
#include   <math.h>
#include   <malloc.h>
#include   "nrutil.h"
#include   "psvfd.h"

int   nx;            /* x-dimension of mesh */
int   nz;            /* z-dimension of mesh (half) */
int   nz2;
int   nl=15;         /* left region of mesh */
int   nf;            /* fluid meshes */
int   nd=15;         /* incident medium meshes */
int   nt;            /* number of time steps */
int   nh=10;         /* absorbing zone mesh */
int   nh1;           /* absorbing zone mesh */
int   atten=2;       /* attenuating mode */
int   afilter=4;     /* filter mode */
float h;             /* spatial mesh interval */
float dt;            /* time digitization interval */
float ass=0.92;      /* the reduction of the wavefield */ 
int   source=0;      /* output source type =0 2 P wave; =1 2 S wave */

/* filenames */
char greenfile[80];  /* file for the Green's function from GRT */
char raymodel[80];   /* file for model in GRT calculation */
char fdmodel[80];    /* file for model of heterogeneous region */

/* following for the Kirchoff interfacing */
char kirfile[80];    /* file for Kirchoff output seismograms */
int   kdx=1;         /* The grids separation of output to Kirchhoff */ 
int   kdt=1;         /* The time separation of output to Kirchhoff */ 

float *uwt34, *taper, *filter;

char pvel[60], svel[60], den[60];
int  readpars = 0;
int  reppars  = 0;
int  startiz  = 0;
int  zrandom;
int  seed;
float scalex, scalez, variation, vargrad = 0.0;

main(ac,av)
int ac; 
char **av;
{


   int   flat =1;       /* flag for earth: 1-> flat layers 0->spherical */
   int   it;            /* timestep counter */

    int   nl_skip=0;     /* left region of skip */
   /* following for table */
   /* float ntable [5*NTSIZE]; */
   float *ntable;
   int   *ptable;
   int table_size, getabc();

   FILE *infile, *fop_kir, *open_file();
   struct vel_stress *p1, *pp1;
   int len;

/* following for the I/O handling */
   long int rec_len;

   void zap_fill(), AFilter(), atten_taper(), snapshot(), tstep();
   void kir_record(); 

   int  mtotal, *ns; float dp, rdt, **sy, **matrix(); 
   int i, j; 
   void interpl();

/* for check only */
   int ix1[700], iz1[700], ix2[700], iz2[700], ix3[700], iz3[700], itrace;
   int xtrace = 10, ztrace = 20;
   FILE *fopu1, *fopu2, *fopu3, *fopw1, *fopw2, *fopw3;
   float xmin, zs, rhos, tstart;

   int ist=100, incre=50;
   int plot_trace=0; int plot_snap=0, plot_w=0;
   int iz = 20;  /*output point of kirchhoff */
   float zcore;  /* the distance of the kiroutput from the core */
   int rcore = 1;
   float ztop = 4;
   int izt, tops = 0, topp = 0;
   char topsfile[80], toppfile[80];
   FILE *fop_tops, *fop_topp;
   

   int kirrecord(struct vel_stress *p, int *ptable, float *na, int iz, 
	  float h, float dt, int nx, int nl, int nl_skip, int nh, int kdx, 
	  FILE *kirout, int source);

   setpar(ac,av);
   mstpar("greenfile",     "s", greenfile);
   mstpar("raymodel",      "s", raymodel);
   mstpar("fdmodel",       "s", fdmodel);
   getpar("NTSIZE",        "d", &NTSIZE);
   mstpar("nx",            "d", &nx);
   mstpar("nl",            "d", &nl);
   getpar("nlskip",        "d", &nl_skip);
   mstpar("nt",            "d", &nt);
   mstpar("h",             "f", &h);
   mstpar("dp",            "f", &dp);
   dt = dp;
   getpar("dt",            "f", &dt);
   getpar("nh",            "d", &nh);
   getpar("nf",            "d", &nf);
   getpar("nd",            "d", &nd);
   getpar("ass",           "f", &ass);
   getpar("atten",         "d", &atten);
   getpar("afilter",       "d", &afilter);
   getpar("source",        "d", &source);
   getpar("kdx",           "d", &kdx);
   getpar("kdt",           "d", &kdt);
   getpar("rcore",         "d", &rcore);
   if(rcore == 1){
       mstpar("zcore",         "f", &zcore);
       mstpar("fdout",         "s", kirfile);
   }
   getpar("tops",          "d", &tops);
   getpar("topp",          "d", &topp);
   getpar("ztop",          "f", &ztop);
   if(tops != 0)
       mstpar("topsfile",  "s", topsfile);
   if(topp != 0)
       mstpar("toppfile",  "s", toppfile);

   getpar("flat",          "d", &flat);

   /* for test only  */
   getpar("xmin",          "f", &xmin);
   getpar("zs",            "f", &zs);
   getpar("rhos",          "f", &rhos);
   getpar("tstart",        "f", &tstart);
   getpar("startsnap",     "d", &ist);
   getpar("incresnap",     "d", &incre);
   getpar("plot_snap",     "d", &plot_snap);
   getpar("plot_trace",    "d", &plot_trace);
   getpar("plot_w",        "d", &plot_w);
   getpar("xtrace",        "d", &xtrace);
   getpar("ztrace",        "d", &ztrace);

   getpar("readpars",      "d", &readpars);
   if(readpars > 0){
       mstpar("pvelocity", "s", pvel);
       mstpar("svelocity", "s", svel);
       mstpar("density",   "s", den);
/*
       mstpar("zrandom",   "d", &zrandom);
       getpar("reppars",   "d", &reppars);
       getpar("startiz",   "d", &startiz);
*/
       mstpar("scalez",   "f", &scalez);
       mstpar("scalex",   "f", &scalex);
       mstpar("seed",     "d", &seed);
       mstpar("variation","f", &variation);
       getpar("vargrad",  "f", &vargrad);
   }

   /* modify the left boundary */

   nl_skip = ((int) ((nl_skip+0.01)/kdx)) * kdx;
   nx     += nl_skip;
   xmin   -= nl_skip *h;

   /* following for reading the record */
   rec_len = (long) (nt*sizeof(float));

   table_size=getabc(flat);
   endpar();

}


#define D    2   /* density */
#define SV   4   /* s-wave velocity */
#define PV   8   /* p-wave velocity */
#define NOTFOUND    0
#define NUMBER      1
#define MATRIX      2

int getabc(flat)

/*
	   Subroutine GETABC evalutes the media parameters at
	each grid point and tabulates them, along with    the
	necessary combinations for the finite      difference
	stencils and boundary conditions. Comparisons eliminate 
	redundant table entries to save space.
*/
int flat;
/* REQUIRES: True                                               */
/* GLOBALS:  nx, nl, nf, nz, nd, nz2: Integer
	     dt, h,                 : Floats                   */


/* MODIFIES: *pptable, *ntable, nx, nz, nd                      */
/* ENSURES:  ta(0)=mu    *dt/h
             tb(0)=ld    *dt/h
	     tc(0)=Bij   *dt/h   (1-facu)/(1+facu)
	     td(0)=ld2mu *dt/h
	     te(0)=Bi1j1 *dt/h   (1-facw)/(1+facw)

             nz2=....
	     nx = nx+nl; nd = -nd;                              */
      

{
   int iz1, iz2, iz10;
   float *na; int *ptable;
   int key, cnt, side;

   float dval, svval, pvval, tdtest, tetest;
   int   dkey, svkey, pvkey;
   int   df,   svf,   pvf;

   int nx2, nf2, nd2;
   int len, i, j, k, l, test, y[6], iz, ix, j2;
   int ii;
   int pt[2];
   int syze;
   float *a, *b, *c; float *pa, *pb, *pc, *ppa, *ppb, *ppc;
   float *pvp, *pvs, *pden;
   float *ppvp, *ppvs, *ppden;
   float vsmax, vpmax, fac, facu, facw;
   float mu, ld, ld2mu, Bij, Bi1j1;

   int jo, nb; float *cc, *ss, *dd, *tth;
   int lfinal, *nen, *ncoun, *nna, *nray;
   int nrec, nfil; float thickness;
   float fvp, fden;
   float ivp, ivs, iden;

   int *np, nll; float *vp, *vs, *rho, **an, **bn; float **matrix();
   int firstlayer;
   void free_matrix(); float lip(); float f0, x;

   FILE *af, *bf, *cf;
   FILE *fop, *open_file();
   int *ratio;
   FILE *ftmp;
   int x0, x1, z0, z1;
    float var;

   double sqrt(), fabs();
   void   copy(), read_model();
   int messy_mod(int nx, int nz, float *pvel, float *svel, float *density, char *model);


/* reading the medium parameters from the raymodel used in the GRT 
   calculation; the left ix<nl; bottom iz<0 and fluid iz>nz regions 
   are specified. 
*/

   firstlayer = 0;

   read_model(flat,raymodel,&jo,&nb,&lfinal,&cc,&ss,&dd,&tth,
		       &nen,&ncoun,&nna,&nray);

   nrec = nna[nen[0]-2];
   nfil = nfinal(jo,ss,nrec);

   thickness=0.0;
   if(nfil >nrec){
       for(i=nrec; i<nfil-1;i++)
           thickness += tth[i];
   } else {
       for(i=nfil; i<nrec-1;i++)
           thickness += tth[i];
   }

   nz2 = (int) (tth[nrec-1]/h) - 3;
   if(nd > nz2){
       fprintf(stderr,"nd change from %d to %d\n",nd,nz2);
       nd = nz2;
   }

/*
*/ 

   nz2 = (int) ((thickness+1.e-3)/(0.5*h));

   nx2  =2*nx;
   nd2  =2*nd;
   nx2 +=2*nl;
   nf2  =2*nf;

   nd2 *= (-1);
   nd  *= (-1);

    fprintf(stderr,"nz2=%d nf2=%d \n", nz2, nf2);	
    free(nen);free(nna);free(ncoun);free(nray);
 
    if(readpars == 0){
   } else {

       scalex /= (0.5 *h); scalez /= (0.5 *h);
       thickness= 0.0;
       for(j=nrec; j<nfil-1; j++)
          thickness += tth[j];
       af = open_file(pvel,"wb");
       bf = open_file(svel,"wb");
       cf = open_file(den,"wb");
       pvp  = (float *) malloc (nx2 * (nz2 + nf2) * sizeof(float));
       pvs  = (float *) malloc (nx2 * (nz2 + nf2) * sizeof(float));
       pden = (float *) malloc (nx2 * (nz2 + nf2) * sizeof(float));
       ppvp = pvp; ppvs = pvs; ppden = pden;
       iz1 = (int) ((thickness+1.e-3)/(0.5*h)) + 2;
       iz10 = iz1;

       vargrad *= (0.5*h)/(6371.0/1221.5); /*times inner-core stretch */
       for(i=nfil-1; i<nfil+400; i++){

           thickness= 0.0;
           for(j=nrec; j<=i; j++)
	       thickness += tth[j];

           iz2 = (int) ((thickness+1.e-3)/(0.5*h)) + 2;
           var = variation - (0.5*(iz2 + iz1) - iz10)*vargrad;
           if(var < 0) var = 0.0;
           if(iz2 > nz2 + nf2 -5) iz2=nz2+nf2-5;

           if(iz2 < iz1)break;

           /* constructing for a model input */
           ivp  = cc[i]; 
           ivs  = ss[i]; 
           iden = dd[i];
           x0   = 0;
           x1   = (nx2-2*nh-2) - (2*nl + 2) -5 +1;
           z0   = -1;
           z1   = iz2 - iz1 +2 +100;

           ftmp = open_file("tmpmodel","w");

           /* for a polygon */
           fprintf(ftmp,"%f %f %f\n", ivp, ivs, iden);
           fprintf(ftmp,"%d %f %f %f %d\n", 4, ivp, ivs, iden, -1);
           fprintf(ftmp,"%f %f %f %d\n", var*ivp, scalex, scalez, seed);
           fprintf(ftmp,"%d %d\n", x0, z0);
           fprintf(ftmp,"%d %d\n", x1, z0);
           fprintf(ftmp,"%d %d\n", x1, z1);
           fprintf(ftmp,"%d %d\n", x0, z1);
           fprintf(ftmp,"%d %d\n", x0, z0);
           fclose(ftmp);

/*

*/

           messy_mod(x1,iz2-iz1+1,ppvp,ppvs,ppden,"tmpmodel");

           fwrite(ppvp,sizeof(float),x1*(iz2-iz1+1),af);
           fwrite(ppvs,sizeof(float),x1*(iz2-iz1+1),bf);
           fwrite(ppden,sizeof(float),x1*(iz2-iz1+1),cf);
	   iz1 = iz2+1;
        }
        free(pvp); free(pvs); free(pden);

   }
    free(cc); free(ss); free(dd); free(tth);
   fclose(af);
   fclose(bf);
   fclose(cf);


   return 1;
}

