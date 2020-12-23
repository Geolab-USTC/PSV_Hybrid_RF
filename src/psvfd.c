
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
float scalex, scalez, variation;

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
   }

   /* modify the left boundary */

   nl_skip = ((int) ((nl_skip+0.01)/kdx)) * kdx;
   nx     += nl_skip;
   xmin   -= nl_skip *h;

   /* following for reading the record */
   rec_len = (long) (nt*sizeof(float));

   table_size=getabc(flat,&ptable,&ntable);
   fprintf(stderr,"table_size=%d\n",table_size);
   endpar();

   /* allocate velocity stress arrays */
   p1= (struct vel_stress *) malloc (nx*(nz-nd+1)*sizeof(struct vel_stress));
   if(p1==NULL)
   { 
      fprintf(stderr,"cannot allocate p1 memory\n"); 
      exit(-1); 
   }
    fprintf(stderr,"pass 1\n");

   zap_fill(p1,0.0,nx*(nz-nd+1));
   pp1 =p1-nx*nd; 

   if(afilter){ 
       filter   = (float *) malloc (2*nh*sizeof(float));
       AFilter(nh,filter); 
   }    
   if(atten){
       taper = (float *) malloc (2*nh*sizeof(float));
       atten_taper(nh,ass,taper);
   }
  
/* following for the interpolation */
   infile=open_file(greenfile,"rb");
   mtotal = 4*(nx - nl + nz+1 - nd);
   uwt34   = (float *) malloc (mtotal*sizeof(float));
   sy      = matrix(0,(mtotal-1),0,3);
   ns      = (int *) malloc (4*sizeof(int));

   rdt = dt/dp;
   for(j=0; j<4; j++){
       ns[j]  = j;
       for(i=0; i<mtotal; i++)
	   sy[i][j] = 0.0;
   }
   nt = (int) (((float)((nt-3)*dp))/dt);

   
   if(rcore != 0)
       fop_kir=open_file(kirfile,"wb");
   if(topp != 0)
       fop_topp=open_file(toppfile,"wb");
   if(tops != 0)
       fop_tops=open_file(topsfile,"wb");


   if(plot_trace){

       /*  Output for test only */
       fopu1 = open_file("seis.u","w");
       fopw1 = open_file("seis.w","w");

       for(itrace=0; itrace<xtrace; itrace++){
           for(j=0; j<ztrace; j++){
               ix1[itrace+j*xtrace] =  nl + 100*itrace;
               iz1[itrace+j*xtrace] = 8*j;
           }
       }

       fprintf(fopu1,"%f %f %f %f %f %f %f %f %f \n",
	      xmin,xmin,nz*h,nt*dt,h,dt,rhos,zs,tstart);
       fprintf(fopu1,"%d %d\n",nt,xtrace*ztrace);
       fprintf(fopw1,"%f %f %f %f %f %f %f %f %f \n",
	      xmin,xmin,nz*h,nt*dt,h,dt,rhos,zs,tstart);
       fprintf(fopw1,"%d %d\n",nt,xtrace*ztrace);

       for(itrace=0; itrace<xtrace; itrace++)
           for(j=0; j<ztrace; j++)
               fprintf(fopu1,"%f ",ix1[itrace+j*xtrace]*h);
           fprintf(fopu1,"\n");

       for(itrace=0; itrace<xtrace; itrace++)
           for(j=0; j<ztrace; j++)
               fprintf(fopu1,"%f ",-iz1[itrace+j*xtrace]*h);
           fprintf(fopu1,"\n");

       for(itrace=0; itrace<xtrace; itrace++)
           for(j=0; j<ztrace; j++)
               fprintf(fopw1,"%f ",ix1[itrace+j*xtrace]*h);
           fprintf(fopw1,"\n");

       for(itrace=0; itrace<xtrace; itrace++)
           for(j=0; j<ztrace; j++)
               fprintf(fopw1,"%f ",-iz1[itrace+j*xtrace]*h);
           fprintf(fopw1,"\n");


       /*  End of output header for test only */
   }

   iz  = (int) (zcore*6371.0/3480.0/h);     /* zcore km away from the CMB */
   iz  = (int) (zcore/h);     /* zcore km away from the CMB */
   iz += nz2/2;
   izt = -(int)(ztop*6371.0/3480.0/h);    /* ztop km above the interface */
   izt = -(int)(ztop/h);    /* ztop km above the interface */
   nh1 = nh;
   if(nh > -nd - 7)
      nh1 = -nd - 7;
           fprintf(stderr,"nd=%d izt=%d\n",nd,izt);
   
   nl_skip = (int) ((nl_skip+0.01)/kdx);
   /* begin time iteration */
   for(it=0; it< nt; it++) {

         interpl(infile,rec_len,it,rdt,mtotal,ns,sy,uwt34,nx-nl,nz+1-nd);
         tstep(pp1,ptable,ntable);

         if(!(it % kdt)){
             if(rcore != 0)
                 kirrecord(pp1,ptable,ntable,iz,h,dt,nx,nl,nl_skip,nh,kdx,fop_kir,source);
	     if(topp != 0)
                 kirrecord(pp1,ptable,ntable,izt,h,dt,nx,nl,nl_skip,nh,kdx,fop_topp,0);
	     if(tops != 0)
                 kirrecord(pp1,ptable,ntable,izt,h,dt,nx,nl,nl_skip,nh,kdx,fop_tops,1);
         }

         if(plot_trace){ 
	     for(itrace=0; itrace<xtrace*ztrace; itrace++)
                 fprintf(fopu1,"%e ",pp1[iz1[itrace]*nx+ix1[itrace]].u);
             fprintf(fopu1,"\n");
	     for(itrace=0; itrace<xtrace*ztrace; itrace++)
                 fprintf(fopw1,"%e ",pp1[iz1[itrace]*nx+ix1[itrace]].w);
             fprintf(fopw1,"\n");

         }    
         if(plot_snap){
	     /* Output snapshots */
             if(it==ist)snapshot("out1",pp1,nx,nd,nz,plot_w);
             if(it==ist+incre)snapshot("out2",pp1,nx,nd,nz,plot_w);
             if(it==ist+2*incre)snapshot("out3",pp1,nx,nd,nz,plot_w);
             if(it==ist+3*incre)snapshot("out4",pp1,nx,nd,nz,plot_w);
             if(it==ist+4*incre)snapshot("out5",pp1,nx,nd,nz,plot_w);
             if(it==ist+5*incre)snapshot("out6",pp1,nx,nd,nz,plot_w);
             if(it==ist+6*incre)snapshot("out7",pp1,nx,nd,nz,plot_w);
             if(it==ist+7*incre)snapshot("out8",pp1,nx,nd,nz,plot_w);
             if(it==ist+8*incre)snapshot("out9",pp1,nx,nd,nz,plot_w);
             if(it==ist+9*incre)snapshot("out10",pp1,nx,nd,nz,plot_w);
             if(it==ist+10*incre)snapshot("out11",pp1,nx,nd,nz,plot_w);
             if(it==ist+11*incre)snapshot("out12",pp1,nx,nd,nz,plot_w);
             if(it==ist+12*incre)snapshot("out13",pp1,nx,nd,nz,plot_w);
             if(it==ist+13*incre)snapshot("out14",pp1,nx,nd,nz,plot_w);
             if(it==ist+14*incre)snapshot("out15",pp1,nx,nd,nz,plot_w);
             if(it==ist+15*incre)snapshot("out16",pp1,nx,nd,nz,plot_w);
             if(it==ist+16*incre)snapshot("out17",pp1,nx,nd,nz,plot_w);
             if(it==ist+17*incre)snapshot("out18",pp1,nx,nd,nz,plot_w);
             if(it==ist+18*incre)snapshot("out19",pp1,nx,nd,nz,plot_w);
             if(it==ist+19*incre)snapshot("out20",pp1,nx,nd,nz,plot_w);
             if(it==ist+20*incre)snapshot("out21",pp1,nx,nd,nz,plot_w);
             if(it==ist+21*incre)snapshot("out22",pp1,nx,nd,nz,plot_w);
             if(it==ist+22*incre)snapshot("out23",pp1,nx,nd,nz,plot_w);
             if(it==ist+23*incre)snapshot("out24",pp1,nx,nd,nz,plot_w);
	 }

   }
   fclose(infile); fclose(fop_kir);

   if(plot_trace){ fclose(fopw1); fclose(fopu1);} 
}

void snapshot(file,p1,nx,nd,nz,plot_w)
char *file; 
struct vel_stress *p1;
int nx,nd,nz; int plot_w;
{
    int ix, iz, increx=3, increz=3;
    FILE *fop, *open_file();
    double fabs();

    fop=open_file(file,"w");
    for(ix=0; ix<nx; ix += increx)
        for(iz=nd; iz<=nz; iz +=increz){
            if(plot_w)fprintf(fop,"%d %d %e \n",ix,iz,p1[iz*nx+ix].w);
            if(!plot_w)fprintf(fop,"%d %d %e \n",ix,iz,p1[iz*nx+ix].u);
        }
    fclose(fop);        
}    


#define   BEFORE              1
#define   AFTER               0
#define   Finite_Diff(ix,iz)  inner(p1,ptable,ntable,ix,iz,stress,order)

void tstep(p1,ptable,ntable)
int *ptable;
float *ntable;
struct vel_stress *p1;

/*
	   Subroutine TSTEP is responsible for generating the next
	timestep from the previous timestep p1,  for the grid and 
	material parameters represented by ptable, ntable, nx, nl 
	and nz.  It takes care of applying boundary conditions and 
	calling  the finite difference subroutines.
*/

/* REQUIRES: p1, ptable, ntable                               */
/* GLOBALS:  nx, nz, nl, nh, nd  Integers                     */
/* MODIFIES: p1[nd*nx..nz*nx+nx-1]                            */
/* ENSURES:  p1 overwrites p1 at previous step                */ 

{
   int ix, iz;
   struct vel_stress *pp1;
   int *pptable;
   int offset, cnt, stress, order;

   void absorb(), inner(), inner1(), inner0(), atten_wave();

   for(stress=1; stress>=0;stress--){
       if(!stress){
         absorb(p1,ptable,ntable,      nd*nx, 1, nx,   nx,BEFORE); 
         absorb(p1,ptable,ntable,    nz*nx+1, 1,-nx, nx-2,BEFORE); 
         absorb(p1,ptable,ntable,(nd+2)*nx-1,nx, -1,nz-nd,BEFORE); 
         absorb(p1,ptable,ntable,  (nd+1)*nx,nx,  1,nz-nd,BEFORE); 
       }

       order =4;

   /*  use fourth order for most grid points */
       for(iz=nd+2; iz <nz-1; iz++)
           for(ix=2; ix<nl-1; ix++)
              Finite_Diff(ix,iz);

       for(ix=nl+2; ix<nx-2; ix++){
           for(iz=2; iz<nz-1; iz++)
              Finite_Diff(ix,iz);
           for(iz=nd+2; iz < -1; iz++)
              Finite_Diff(ix,iz);
       }             

      /* apply second order code to edges */
      
       order =2;

       for(ix=2; ix<nl-1; ix++){
          Finite_Diff(ix,nz-1);
          Finite_Diff(ix,nd+1);
       }         

       for(ix=nl+2; ix<nx-2; ix++){
          Finite_Diff(ix,nz-1);
          Finite_Diff(ix,1);
          Finite_Diff(ix,-1);
          Finite_Diff(ix,nd+1);
       }         

       for(iz=nd+1; iz<nz; iz++){
          Finite_Diff(1,iz);
          Finite_Diff(nl-1,iz);
          if(iz==0)continue;
          Finite_Diff(nl+1,iz);
          Finite_Diff(nx-2,iz);
       }         

      offset = nx*nd + nl;
      pp1     = p1     + offset;
      pptable = ptable + offset;
      inner1(pp1,pptable,ntable,stress);


      offset = nl + 1;
      pp1     = p1     + offset;
      pptable = ptable + offset;
      inner0(pp1,pptable,ntable,stress);
      
      if(stress){
         
          for(ix=1; ix<nl; ix++){          
              Finite_Diff(ix,nz);
              Finite_Diff(ix,nd);
          }


          for(ix=nl+1; ix<nx-1; ix++){
              Finite_Diff(ix,nz);
              Finite_Diff(ix,nd);
          }

          for(iz=nd; iz<nz+1; iz++){
              Finite_Diff(0,iz);
              if(iz==0)continue;
              Finite_Diff(nx-1,iz);
          }              

      } else {
         absorb(p1,ptable,ntable,      nd*nx, 1, nx,   nx,AFTER); 
         absorb(p1,ptable,ntable,    nz*nx+1, 1,-nx, nx-2,AFTER); 
         absorb(p1,ptable,ntable,(nd+2)*nx-1,nx, -1,nz-nd,AFTER); 
         absorb(p1,ptable,ntable,  (nd+1)*nx,nx,  1,nz-nd,AFTER); 
     }
   }


   if(atten)atten_wave(p1,nx,nz,nd,nh,atten,taper);

}

void zap_fill(pp,value,n)
struct vel_stress *pp;
float value;
int n;

/* REQUIRES: pp[0..n-1]                                            */
/* MODIFIES: pp[0..n-1]                                            */
/* ENSURES:  pp[0..n-1].stress = value                             */

{
   while(n--){
      pp->u=pp->w=pp->t11=pp->t12=pp->t22=value;
      pp++;
   }
}

/* Subroutine INNER applies finite-difference stencil */

void inner(pp1,ptable,na,ix,iz,stress,order)
struct vel_stress *pp1;
int *ptable;
float *na;
int ix, iz, stress, order;

/* REQUIRES: ix:       in Integer, x-dimension grid         
             iz:       in integer, z-dimension grid
             stress:   in 0 or 1, 0 for velocities, 1 for stresses
             order:    in 2 or 4, second or fourth order stencil    */ 
/* GLOBALS:  nx, nh, nd, afilter nz:  Integers
             filter                :  Float Pointers                */       
/* MODIFIES: pp1[iz*nx+ix], NOT T12     at left and top boundaries
			    NOT T11 T22 at bottom and right boundaries
			    NOT U W     at all boundaries*           /
/* ENSURES:  pp1[iz*nx+ix].stress is the value of next time step    */
/* CALLED: by TSTEP                                                 */

{
   float uxx, uzz, wxx, wzz, t11xx, t12xx, t12zz, t22zz;
   float *pfilter;
   int offset, *pt;
   struct vel_stress *p1;

   offset = iz*nx+ix;
   pt     = ptable + offset;
   p1     = pp1    + offset;

   if(stress){
       if(!(ix==nx-1 || iz==nd)) uxx   =  p1[ 1].u-p1[  0].u;
       if(!(ix==nx-1 || iz==nd)) wzz   =  p1[ 0].w-p1[-nx].w;
       if(!(ix==0    || iz==nz)) uzz   =  p1[nx].u-p1[  0].u;
       if(!(ix==0    || iz==nz)) wxx   =  p1[ 0].w-p1[ -1].w;
       if(order==4){
           uxx  -= (p1[ 2].u-p1[   -1].u)/27.0;
	   uzz  -= (p1[2*nx].u-p1[-nx].u)/27.0;
	   wxx  -= (p1[ 1].w-p1[-2   ].w)/27.0;
	   wzz  -= (p1[nx].w-p1[-2*nx].w)/27.0;
       }
       
       if(afilter >0 && iz>nz-nh){
           pfilter = filter+2*(iz-nz+nh-1);
           uxx    *= pfilter[0];
           wxx    *= pfilter[1];
       }    
       if(afilter >1 && ix>=nx-nh){
           pfilter = filter+2*(ix-nx+nh);
           uzz    *= pfilter[0];
           wzz    *= pfilter[1];
       }    
       if(afilter >2 && iz<nd+nh1){
           pfilter = filter+2*(nd+nh1-1-iz);
           uxx    *= pfilter[1];
           wxx    *= pfilter[0];
       }    
           
       if(afilter >3 && ix<nh){
           pfilter = filter+2*(nh-1-ix);
           uzz    *= pfilter[1];
           wzz    *= pfilter[0];
       }    
       
	if(!(ix==nx-1 || iz==nd))p1->t11 +=  td(0)* uxx + tb(0)*wzz;
	if(!(ix==nx-1 || iz==nd))p1->t22 +=  tb(0)* uxx + td(0)*wzz;
	if(!(ix==0    || iz==nz))p1->t12 +=  ta(0)*(uzz +       wxx);  
      }
   else {
      
       t11xx = p1[  0].t11 - p1[ -1].t11;
       t12zz = p1[  0].t12 - p1[-nx].t12;
       t12xx = p1[  1].t12 - p1[  0].t12;
       t22zz = p1[ nx].t22 - p1[  0].t22;
       if(order==4){
           t11xx  -= (p1[   1].t11 - p1[   -2].t11)/27.0;
           t12zz  -= (p1[  nx].t12 - p1[-2*nx].t12)/27.0;
           t12xx  -= (p1[   2].t12 - p1[   -1].t12)/27.0;
           t22zz  -= (p1[2*nx].t22 - p1[  -nx].t22)/27.0;
       }
       if(afilter >0 && iz>nz-nh){
           pfilter  = filter+2*(iz-nz+nh-1);
           t12xx   *= pfilter[1];
           t11xx   *= pfilter[0];
       }    
       if(afilter >1 && ix>=nx-nh){
           pfilter  = filter+2*(ix-nx+nh);
           t22zz   *= pfilter[1];
           t12zz   *= pfilter[0];
       }    
       if(afilter >2 && iz<nd+nh1){
           pfilter  = filter+2*(nd+nh1-1-iz);
           t12xx   *= pfilter[0];
           t11xx   *= pfilter[1];
       }    
       if(afilter >3 && ix<nh){
           pfilter  = filter+2*(nh-1-ix);
           t12zz   *= pfilter[1];
           t22zz   *= pfilter[0];
       }    
       

       p1->u += tc(0)*(  t11xx +    t12zz);
       p1->w += te(0)*(  t12xx +    t22zz);
   }
}


/* Subroutine INNER0 applies 2nd finite-difference stencil at iz=0  */
/* T11 T22 AND U iz=0 are the reflected wavefield
   T12 and W          are the total wavefield                       */

void inner0(p1,pt,na,stress)
struct vel_stress *p1;
int *pt;
float *na;
int stress;

/* REQUIRES: stress:   in 0 or 1, 0 for velocities, 1 for stresses
             p1[offset], offset starts at ix=0 iz=0 if stress=1
                                       at ix=1 iz=0 if stress=0                          
             cnt:      in Interger, cnt=nx    if stress=1
                                    cnt=nx-2  if stress=0           */               
/* GLOBALS:  nx, nh, afilter, nl, nz:  Integers
             filter, uwt34          :  Float Pointers               */       
/* MODIFIES: p1[offset] T11 T22: iz=0, 0<=ix<nx-1
                        T12    : iz=0, 1<=ix<=nx-1
                        U,W    : iz=0, 1<=ix<nx-1                   */
/* ENSURES:  p1[iz*nx+ix].stress is the value of next time step     */
/* CALLED: by TSTEP                                                 */

{
   int ix;
   float *u3, *t223, *w4, *t124;
   float uxx, uzz, wxx, wzz, t11xx, t12xx, t12zz, t22zz;
   float *pfilter;

   u3     = uwt34              + 1;
   t223   = uwt34 +  (nx - nl) + 1;
   t124   = uwt34 +2*(nx - nl) + 1;
   w4     = uwt34 +3*(nx - nl) + 1;

   p1--; 
   pt--;
   u3--;  t223--; w4--; t124--;
   for(ix=nl+1; ix<nx; ix++)
   {
      p1++;
      pt++;
      u3++; t223++; w4++; t124++;
      if(stress){

          uxx   =  p1[ 1].u-p1[  0].u;
          uzz   =  p1[nx].u-(p1[  0].u+u3[0]);
          wxx   =  p1[ 0].w-p1[ -1].w;
          wzz   = (p1[ 0].w-w4[0])-p1[-nx].w;
          
          if(afilter>1 && ix>=nx-nh){
              pfilter = filter+2*(ix-nx+nh);
              uzz    *= pfilter[0];
              wzz    *= pfilter[1];
          }    
          
	  if(ix!=nx-1){
	      p1->t11 +=  td(0)* uxx + tb(0)*wzz;
	      p1->t22 +=  tb(0)* uxx + td(0)*wzz;
          }

	  p1->t12 +=  ta(0)*(uzz +       wxx);  
	
      }
      else {
      
          t11xx = p1[  0].t11 - p1[ -1].t11;
          t12zz = (p1[  0].t12-t124[0]) - p1[-nx].t12;
          t12xx = p1[  1].t12 - p1[  0].t12;
          t22zz = p1[ nx].t22 - (p1[  0].t22+t223[0]);
          if(afilter>1 && ix>=nx-nh){
              pfilter  = filter+2*(ix-nx+nh);
              t22zz   *= pfilter[1];
              t12zz   *= pfilter[0];
          }    
          
          if(ix!=nx-1){
	      p1->u += tc(0)*(  t11xx +    t12zz);
              p1->w += te(0)*(  t12xx +    t22zz);
          }

      }
   }
}

/* note the filter is not working now */
/* Subroutine INNER1 applies 2nd finite-difference stencil at ix=nl iz>0 */
/* T11 T22 AND W ix=nl are the total wavefield
   T12 and U          are the reflected wavefield                     */
   
void inner1(p1,pt,na,stress)
struct vel_stress *p1;
int *pt;
float *na;
int stress;

/* REQUIRES: stress:   in 0 or 1, 0 for velocities, 1 for stresses
             p1[offset], offset starts at ix=0 iz=1 if stress=1
             cnt:      in Interger, cnt=nz+1    if stress=1
                                    cnt=nz      if stress=0           */               
/* GLOBALS:  nx, nh, afilter, nl, nz:  Integers
             filter, uwt34          :  Float Pointers               */       
/* MODIFIES: p1[offset] T11 T22  : ix=nl, 1<=iz<=nz
                        T12 W U  : ix=nl, 1<=iz<=nz-1
/* ENSURES:  p1[iz*nx+ix].stress is the value of next time step     */
/* CALLED: by TSTEP                                                 */

{
   float *u0, *t110, *w0, *t120;
   int iz; 
   float uxx, uzz, wxx, wzz, t11xx, t12xx, t12zz, t22zz;
   float *pfilter;
   float t22nl, wnl; 

   pfilter = filter;
   u0      = uwt34  + 4*(nx -nl);
   t120    = uwt34  + 4*(nx -nl)+  (nz+1 -nd);
   t110    = uwt34  + 4*(nx -nl)+2*(nz+1 -nd);
   w0      = uwt34  + 4*(nx -nl)+3*(nz+1 -nd);
   wnl     = uwt34[3*(nx-nl)];
   t22nl   = uwt34[(nx-nl)];

   p1 -=nx; pt -=nx;
   u0--; t110--; w0--; t120--;
   for(iz=nd; iz<=nz; iz++)
   {
      p1 +=nx; pt +=nx;
      u0++; t110++; w0++; t120++;
      if(stress){
          
	  if(iz > nd){
              uxx   =  p1[ 1].u-(p1[  0].u + u0[0]);
              wzz   =  p1[ 0].w-p1[-nx].w;
	      if(iz == 0) wzz -= wnl;
              if(afilter>0 && iz>=nz-nh+1)
                  uxx *= pfilter[(iz-(nz-nh+1))*2];
              if(afilter>0 && iz<=nd+nh1)
                  uxx *= filter[2*(nd+nh1-1-iz) + 1];

	      p1->t11 +=  td(0)* uxx + tb(0)*wzz;
	      p1->t22 +=  tb(0)* uxx + td(0)*wzz;
          }

          if(iz < nz){
              uzz   =  p1[nx].u-p1[  0].u;
              wxx   = (p1[ 0].w-w0[0])-p1[ -1].w;
              if(afilter>0 && iz>=nz-nh+1)
                  uxx *= pfilter[(iz-(nz-nh+1))*2 + 1];
              if(afilter >2 && iz<nd+nh1)
                  wxx *= filter[2*(nd+nh1-1-iz) + 0];

	      p1->t12 +=  ta(0)*(uzz + wxx);  
          }


      }
      else {
          
	  if(iz > nd && iz < nz){
              t11xx = (p1[  0].t11-t110[0]) - p1[ -1].t11;
              t12zz = p1[  0].t12 - p1[-nx].t12;
              t12xx = p1[  1].t12 - (p1[  0].t12+t120[0]);
              t22zz = p1[ nx].t22 - p1[  0].t22;
	      if(iz == 0) t22zz -= t22nl;
              if(afilter>0 && iz>=nz-nh+1){
                  t11xx *= pfilter[0];
                  t12xx *= pfilter[1];
                  pfilter += 2;
              }    
              if(afilter >2 && iz<nd+nh1){
                   pfilter  = filter+2*(nd+nh1-1-iz);
                   t12xx   *= pfilter[0];
                   t11xx   *= pfilter[1];
              }    
              p1->u += tc(0)*(  t11xx +    t12zz);
              p1->w += te(0)*(  t12xx +    t22zz);
	  }
      }
   }
}


#define D    2   /* density */
#define SV   4   /* s-wave velocity */
#define PV   8   /* p-wave velocity */
#define NOTFOUND    0
#define NUMBER      1
#define MATRIX      2

int getabc(flat,pptable,ntable)

/*
	   Subroutine GETABC evalutes the media parameters at
	each grid point and tabulates them, along with    the
	necessary combinations for the finite      difference
	stencils and boundary conditions. Comparisons eliminate 
	redundant table entries to save space.
*/
int flat;
int **pptable;
float **ntable;

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
   int iz1, iz2;
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
   ivp = cc[nrec-1]; ivs = ss[nrec-1]; iden= dd[nrec-1];
   fvp = cc[nfil-1];                   fden= dd[nfil-1];


   nx2  =2*nx;
   nd2  =2*nd;
   nx2 +=2*nl;
   nf2  =2*nf;

   len= nx2*(nz2+1+nf2+nd2+1)*sizeof(float);
   a=        (float *) malloc(len);
   b=        (float *) malloc(len);
   c=        (float *) malloc(len);
   *pptable= (int   *) malloc(nx2*(nz2+nf2+1+nd2+4)/4*sizeof(int));
   ppa =a; ppb= b; ppc=c;
   a += nx2*nd2; b +=nx2*nd2; c +=nx2*nd2; *pptable +=nx2*nd2/4;
   ptable=*pptable;
   if(a == NULL || b == NULL || c == NULL || *pptable == NULL)
   { 
      fprintf(stderr,"cannot allocate a, b, c memory\n"); 
      exit(-1); 
   }
   nd2 *= (-1);
   nd  *= (-1);

   cnt=-nd2*nx2+2*nx2;
   ivp *= ivp; ivs *= ivs;
   for(pa=a+2*nx2-1,pb=b+2*nx2-1,pc=c+2*nx2-1; cnt--; pa--,pb--,pc--){
       *pa=ivs; *pb =ivp; *pc= iden;
   }

/*
   cnt=(nf2-2)*nx2;
   fvp *=fvp;
   len=nx2*(nz2+1+2);
   for(pa=a+len,pb=b+len,pc=c+len; cnt--; pa++,pb++,pc++){
          *pa=0.00000; *pb =fvp; *pc= fden;
   }
*/

/*
   iz1 = 1;
   for(i=nrec-1; i>0; i--){

       thickness= 0.0;
       for(j=nrec-1; j>=i; j--)
	   thickness += tth[j];

       iz2 = (int) (-(thickness+1.e-3)/(0.5*h)) + 2;
       if(iz2 < nd2) iz2=nd2;

       if(iz2 > iz1)break;

       ivp  = cc[i]*cc[i]; 
       ivs  = ss[i]*ss[i]; 
       iden = dd[i];

       for(j=iz1; j>=iz2; j--){
	   cnt=2*nl;
	   cnt=nx2;
           for(pa=a+j*nx2,pb=b+j*nx2,pc=c+j*nx2; cnt--;pa++,pb++,pc++){
	       *pa=ivs; *pb=ivp; *pc=iden;
           }
       }

       iz1=iz2-1;
    }
*/
   
   iz1=2;
/*
   for(i=nrec; i<nfil-1; i++){
*/
   for(i=nrec; i<nfil+400; i++){

       thickness= 0.0;
       for(j=nrec; j<=i; j++)
	   thickness += tth[j];

       iz2 = (int) ((thickness+1.e-3)/(0.5*h)) + 2;
       if(iz2 > nz2 + nf2) iz2=nz2+nf2;

       if(iz2 < iz1)break;

       ivp  = cc[i]*cc[i]; 
       ivs  = ss[i]*ss[i]; 
       iden = dd[i];

       for(j=iz1; j<=iz2; j++){
	   cnt=2*nl;
	   /* test only */
	   cnt=nx2;
           for(pa=a+j*nx2,pb=b+j*nx2,pc=c+j*nx2; cnt--;pa++,pb++,pc++){
	       *pa=ivs; *pb=ivp; *pc=iden;
           }
       }

       iz1=iz2+1;
    }


    free(nen);free(nna);free(ncoun);free(nray);

 
    if(readpars == 0){
        /* reading the basin model */
        fop  = open_file(fdmodel,"r");
        fscanf(fop,  "%d\n", &nll);

        if(nll < 0){
            firstlayer = 1;
            nll = -nll;
        }


        ratio   = (int   *) malloc ((nll+2)*sizeof(int)); 
        np   = (int   *) malloc ((nll+2)*sizeof(int)); 
        vp   = (float *) malloc ((nll+2)*sizeof(float)); 
        vs   = (float *) malloc ((nll+2)*sizeof(float)); 
        rho  = (float *) malloc ((nll+2)*sizeof(float)); 
        an   = matrix(0,nll-1,0,6000);
        bn   = matrix(0,nll-1,0,6000);
/*
        an   = matrix(0,nll,0,6000);
        bn   = matrix(0,nll,0,6000);
*/
        for(j=0; j<nll; j++){
            ratio[j] = 0;
            fscanf(fop, "%d %f %f %f\n",&np[j],&vp[j],&vs[j],&rho[j]);
            if(np[j] < 0){
                ratio[j] = 1;
                np[j]    = - np[j];
            }
            vp[j] *= vp[j]; vs[j] *= vs[j];
            if(np[j] > 6000){
                fprintf(stderr,"an bn allocation error, np[j] >6000\n");
                exit(-2);
            }
            for(i=0; i<np[j]; i++){
                fscanf(fop,"%f ",&f0);
                an[j][i] = f0;
            } 
            fscanf(fop,"\n");

            for(i=0; i<np[j]; i++){
                fscanf(fop,"%f ",&f0);
                bn[j][i] = f0;
            } 
            fscanf(fop,"\n");
        }
        fclose(fop);
    
    
        for(i=2*nl; i<nx2; i++){
            x   = (i-2*nl)*0.5*h;
            iz1 = 2; 
            if(firstlayer == 1)
                iz1 = (int) ((lip(0,np,an,bn,x)+1.e-3)/(0.5*h))+2 + 1;
            for(j=firstlayer; j<nll; j++){
                iz2 = (int) ((lip(j,np,an,bn,x)+1.e-3)/(0.5*h))+2;
                for(iz = iz1; iz<=iz2; iz++){
                    if(ratio[j] == 1){
                        b[iz*nx2+i] *= vp[j]; a[iz*nx2+i] *= vs[j];
                        c[iz*nx2+i] *=rho[j];
                    } else {
                        b[iz*nx2+i] = vp[j]; a[iz*nx2+i] = vs[j];
                        c[iz*nx2+i] =rho[j];
                    }
                }
               iz1 = iz2 + 1;
            }
        }

        /* assign parameters for the left region: set to those at x=0 */

        for(i=0; i<2*nl; i++){
    	    x = 0.0;
            iz1 = 2; 
            if(firstlayer == 1)
                iz1 = (int) ((lip(0,np,an,bn,x)+1.e-3)/(0.5*h))+2 + 1;
            for(j=firstlayer; j<nll; j++){
                iz2 = (int) ((lip(j,np,an,bn,x)+1.e-3)/(0.5*h))+2;
                for(iz = iz1; iz<=iz2; iz++){
                    if(ratio[j] == 1){
                        b[iz*nx2+i] *= vp[j]; a[iz*nx2+i] *= vs[j];
                        c[iz*nx2+i] *=rho[j];
                    } else {
                        b[iz*nx2+i] = vp[j]; a[iz*nx2+i] = vs[j];
                        c[iz*nx2+i] =rho[j];
                    }
                }
               iz1 = iz2 + 1;
            }
        }

/*
        free(np); free(vp); free(vs); free(rho);
        free(ratio);
        free_dmatrix(an,0,nll-1,0,6000);
        free_dmatrix(bn,0,nll-1,0,6000);
*/
        
   } else {

       scalex /= (0.5 *h); scalez /= (0.5 *h);
       thickness= 0.0;
       for(j=nrec; j<nfil-1; j++)
          thickness += tth[j];

       af = open_file(pvel,"rb");
       bf = open_file(svel,"rb");
       cf = open_file(den,"rb");
       iz1 = (int) ((thickness+1.e-3)/(0.5*h)) + 2;
       for(i=nfil-1; i<nfil+400; i++){

           thickness= 0.0;
           for(j=nrec; j<=i; j++){
	       thickness += tth[j];
}

           iz2 = (int) ((thickness+1.e-3)/(0.5*h)) + 2;
           if(iz2 > nz2 + nf2 -5) iz2=nz2+nf2-5;

           if(iz2 < iz1)break;

           /* constructing for a model input */
           ivp  = cc[i]; 
           ivs  = ss[i]; 
           iden = dd[i];
           x0   = 0;
           x1   = (nx2-2*nh-2) - (2*nl + 2) -5 +1;
           z0   = -1;
           z1   = iz2 - iz1 +2;

	   for(iz=iz1; iz<=iz2; iz++){
                for(ii=2*nl+2; ii<x1+2*nl+2; ii++){
                    fread(&ivp,sizeof(float),1,af);
                    fread(&ivs,sizeof(float),1,bf);
                    fread(&iden,sizeof(float),1,cf);
		    a[iz*nx2 + ii] = ivs * ivs;
		    b[iz*nx2 + ii] = ivp *ivp;
		    c[iz*nx2 + ii] = iden;
                }
           }
	   iz1 = iz2+1;
        }
        fclose(af); fclose(bf); fclose(cf);

   }
    free(cc); free(ss); free(dd); free(tth);


   /* check for stability  */
   vpmax=0.;
   for(i=nx2*nd2;i<nx2*(nz2+nf2+1);i++){
      if(sqrt(b[i])>vpmax){
	  vpmax=sqrt(b[i]);
	  j    = i;
      }

   }
   fac=1.4142*dt*vpmax/h;
   if(fac>=1.0){
      fprintf(stderr,"fac= %4.2f< instability !\n",fac);
      fprintf(stderr,"vpmax= %4.2f at iz=%d ix=%d\n",vpmax, j/nx2/2, j-2*nx2*(j/nx2/2));
      exit(-1);
   }
   else{
      fprintf(stderr,"fac= %4.2f, stable\n",fac);
   }

   /* convert to b = lambda, a = mu, and multiply by grid factor */
   fac=dt/h;
   for(i=nx2*nd2;i<nx2*(nz2+nf2+1);i++){
      b[i]=(b[i]-2.0*a[i])*c[i]*fac;
      a[i]*=c[i]*fac;
      c[i] =1.0/c[i]*fac;
   }

/* Encode media parameters in table */
   fprintf(stderr,"encode media parameters \n");
   *ntable = (float *) malloc (5*NTSIZE*sizeof(float));
   na = *ntable;
   len=0;

    nz =(nz2+nf2+1)/2-1;
    nx = nx2/2;
   /* check for homogeneous regions */


   for(iz=nd+1; iz < nz; iz++){
       for(ix=1; ix<nx-1; ix++){
	   i=nx*iz+ix;
	   pa = a+iz*2*nx2+ix*2; 
	   pb = b+iz*2*nx2+ix*2; 
	   pc = c+iz*2*nx2+ix*2; 

	   ld    = pb[1]; 
	   ld2mu = pb[1]+2.0*pa[1];
	   mu    = pa[nx2];
	   Bij   = pc[0];
	   Bi1j1 = pc[nx2+1];

           /* following for the forth order portion */

	   if(   (ix>1 && ix<nl-1 && iz>nd+1 && iz<nz-1)
	       ||(ix>nl+1 && ix<nx-2 && iz>nd+1 && iz<-1)
	       ||(ix>nl+1 && ix<nx-2 && iz>1 && iz<nz-1)){
	       
	       ld     *=(float) (9.0/8.0);
	       ld2mu  *=(float) (9.0/8.0);
	       mu     *=(float) (9.0/8.0);
	       Bij    *=(float) (9.0/8.0);
	       Bi1j1  *=(float) (9.0/8.0);
           }

	   test=0;

	   if(readpars == 0) { 
 /* skip this if inputing random model */

	   for(j2=0; j2<4; j2++){
	       if(j2==0 && ix==1)continue;
	       if(j2==1 && (ix==1 || iz==nd+1)) continue;
	       if(j2==2 && iz==nd+1)continue;
	       if(j2==3 && (ix==nx-2 || iz==nd+1))continue;
	       if (j2==0) pt[0] = ptable[i-1];
	       else       pt[0] = ptable[i-nx+j2-2];

	       if(ta(0)==mu && tb(0)==ld && tc(0)==Bij && td(0)==ld2mu && te(0)==Bi1j1){
		   test=1;
		   ptable[i]=pt[0];
	       }
            }
            for(j=len-10; j>=0 && test==0; j -=5){
		pt[0]=j;
	       if(ta(0)==mu && tb(0)==ld && tc(0)==Bij && td(0)==ld2mu && te(0)==Bi1j1){
		   test=1;
		   ptable[i]=pt[0];
	       }
            }
            }
	    if(test==0){
		pt[0]=len;
		ta(0)=mu;
		tb(0)=ld;
		tc(0)=Bij;
		td(0)=ld2mu;
		te(0)=Bi1j1;
		ptable[i]=pt[0];

		len +=5;
		if(len>=5*NTSIZE){
		    NTSIZE *= 2;
		    syze = NTSIZE * 5 * sizeof(float);
		    fprintf(stderr, 
		    "redimensioning ntable to %d entries. (%d bytes, %.2f Mbytes)\n" 
		    ,NTSIZE, syze, syze/1.e6);
		    *ntable = (float *) realloc((char *)na,syze);
		    if(*ntable==NULL){
			fprintf(stderr, "redimensioning ntable fails (try again)\n");
		        NTSIZE *= 0.55;
		        syze = NTSIZE * 5 * sizeof(float);
		        fprintf(stderr, 
		        "redimensioning ntable to %d entries. (%d bytes, %.2f Mbytes)\n" 
		        ,NTSIZE, syze, syze/1.e6);
		        *ntable = (float *) realloc((char *)na,syze);
		        if(*ntable==NULL){
			    fprintf(stderr, "Memory re-allocation for ntable fails (quit)\n");
			    exit (-3);
	                }
	            }
		    if(*ntable != na) copy(na,*ntable,len);
		    na = *ntable;
		    fflush(stderr);
	        }
            }
        }	    
    }

   fprintf(stderr,"finish the forth order portion, now start edges\n");
   /* tabulate edges and                                      */
   /* set the edges to (1-v*dt/h)/(1+v*dt/h) for the absorbing*/
   /* boundary conditions, a[] for u comp, b[] for w comp     */
   fprintf(stderr,"do edges\n");
   for(k=0;k<nx2+2*(nz-nd)-2;k++)
   {
      /* left bottom corner */
      if(k==0) {
         i =nd*nx+k;
	 j =nd2*nx2+2*k;
         side=1;
      }
      /* bottom boundary 1<=ix<=nx-2 */
      else if(k<nx-1) {
         i = nd*nx+k;
	 j = nd2*nx2+2*k;
         side=0;
      }
      /* right boundary ix=nx-1 iz=nd--nz */
      else if(k<nx+nz-nd){ 
         i = nx*(k-nx+2+nd)-1;
	 j = 4*nx*(k-nx+2+nd)-nx2-2;
         side=1;
      }
      /* top boundary iz=nz, ix=1--nx-2 */
      else if(k<nx2+nz-nd-2){ 
         i = k-nx-nz+nd+nx*nz+1;
	 j = 2*(k-nx-nz+nd)+4*nx*nz+2;
         side=0;
      }
      /* left boundary ix=0 iz=nd--nz */
      else { 
         i = (k-nx2-nz+nd+2+nd+1)*nx;
         j = 2*(k-nx2-nz+nd+2+nd+1)*nx2;
         side=1;
      }
      pa = a + j;
      pb = b + j;
      pc = c + j;
      ld = pb[1];
      ld2mu= pb[1]+2.0*pa[1];
      mu =pa[nx2];
      test=0;

      if(side==0){
         facu=sqrt(pa[0]*pc[0]);
         facw=sqrt((pb[1]+2.*pa[1])*pc[1]);
         if(pa[0]<0.001){
             facu=1.0;
         }
      }
      else{
         facu=sqrt((pb[0]+2.*pa[0])*pc[0]);
         facw=sqrt(pa[nx2]*pc[nx2]);
         if(pa[nx2]<0.001){
             facw=1.0;
         }
      }
      tdtest= (1.0-facu)/(1.0+facu);
      tetest= (1.0-facw)/(1.0+facw);
      
/*
      for(j=len-5;j>=0;j-=5){
         pt[0]=j;
         if(tdtest==tc(0)&&tetest==te(0)&& mu==ta(0)&& ld==tb(0)&& ld2mu==td(0)) 
         {
            test=1;
            ptable[i]=pt[0];
         }
      }
*/
     
      if(test==0){
         pt[0]=len;
         tc(0)=tdtest;
         te(0)=tetest;
	 ta(0)=mu;
         tb(0)=ld;
         td(0)=ld2mu;
         ptable[i]=pt[0];
         len+=5;
         if(len>5*NTSIZE){
	    NTSIZE += nx2+2*(nz-nd);
	    syze = NTSIZE * 5 * sizeof(float);
	    fprintf(stderr,
	    "redimensioning ntable to %d entries. (%d bytes, %.2f Mbytes)\n"
	    ,NTSIZE,syze,syze/1000000.);
	    *ntable = (float *) realloc((char *)na,syze);
	    if (*ntable==NULL) {
	       fprintf(stderr,"Memory re-allocation for ntable fails!\n");
	       exit(-3);
	    }
	    if (*ntable != na) copy(na,*ntable,len);
	    na = *ntable;
	    fflush(stderr);
         }
      }
   }
   fprintf(stderr,"done with table formation, %d types of media! \n",len/5);
   NTSIZE = len/5 + 1;
   syze = NTSIZE * 5 * sizeof(float);
   fprintf(stderr,
   "final ntable size: %d entries. (%d bytes, %.2f Mbytes)\n"
   ,NTSIZE,syze,syze/1000000.);


   *ntable = (float *) realloc((char *)na,syze);

   if (*ntable != na) copy(na,*ntable,len);
   fflush(stderr);

   free(ppa);
   free(ppb);
   free(ppc);

   return(syze);
}

float lip(layer,np,an,bn,x)
int layer, *np; float **an, **bn; float x;
{
    float dep; int i;
    if( x < an[layer][0] || x > an[layer][np[layer]-1] ){
        fprintf(stderr,"fdmodel out of range\n");
        exit(-3);
    }
    for(i=0; i<np[layer]-1; i++)
        if(x >=an[layer][i]  && x <an[layer][i+1]){
            dep  = (bn[layer][i+1]-bn[layer][i])/(an[layer][i+1]-an[layer][i]);
            dep  = bn[layer][i] + dep*(x-an[layer][i]); 
        }
    return dep;
}
            
void copy(a,b,n)
register float *a, *b;
register int n;

/* REQUIRES: a[0..n-1], b[0..n-1]                           */
/* MODIFIES: b[i],      i=0..n-1                            */
/* ENSURES:  b[i]=a[i]  i=0..n-1                            */

{
   while (n--) *b++ = *a++;
}



void atten_taper(nh,ass,taper)
int nh;
float ass, *taper;

/* REQUIRES: ass Positive                                      */
/* MODIFIES: taper[0..nh-1]                                    */
/* ENSURES:  taper[0]=1, taper[nh-1]=ass, taper exp function   */

{
    int n;
    double alpha, log10(), exp();

/*    alpha=0.25/(nh*nh)*log10((double) (ass))/log10(2.73); */
    alpha=0.5/nh*log10((double) (ass))/log10(2.73);
    
    for(n=0; n<2*nh; n++)
        taper[n]=exp((double) (alpha*n));
}

void atten_wave(p,nx,nz,nd,nh,atten,taper)
register struct vel_stress *p;
register int nh,nx,nz,nd,atten;
register float *taper;
{
  int ix, iz, offset;
  float *ptaper;
  
  for (ix=0; ix<nx; ix++)
      for (iz=nd; iz<=nz; iz++){ 
           if(atten>0 && iz>nz-nh){
               ptaper = taper +2*(iz-nz+nh-1);
               offset = iz*nx+ix;
               p[offset].u   *= ptaper[0];
               p[offset].t11 *= ptaper[0];
               p[offset].t22 *= ptaper[0];
               p[offset].w   *= ptaper[1];
               p[offset].t12 *= ptaper[1];
           }
           if(atten >1 && ix>=nx-nh){
               ptaper = taper +2*(ix-nx+nh);
               offset = iz*nx+ix;
               p[offset].u   *= ptaper[0];
               p[offset].t11 *= ptaper[1];
               p[offset].t22 *= ptaper[1];
               p[offset].w   *= ptaper[1];
               p[offset].t12 *= ptaper[0];
           }

           if(atten >2 && iz<nd+nh1){
               ptaper = taper +2*(nd+nh1-1-iz);
               offset = iz*nx+ix;
               p[offset].u   *= ptaper[1];
               p[offset].t11 *= ptaper[1];
               p[offset].t22 *= ptaper[1];
               p[offset].w   *= ptaper[0];
               p[offset].t12 *= ptaper[0];
           }

           if(atten >3 && ix<nh){
               ptaper = taper +2*(nh-1-ix);
               offset = iz*nx+ix;
               p[offset].u   *= ptaper[1];
               p[offset].t11 *= ptaper[0];
               p[offset].t22 *= ptaper[0];
               p[offset].w   *= ptaper[0];
               p[offset].t12 *= ptaper[1];
           }
       }
}


void absorb(p1,pt,na,offset,tang,norm,cnt,before)
register struct vel_stress *p1;
register float *na;
register int *pt;
register int tang, norm, cnt, before;
int offset;

/* REQUIRES p1[offset+tang*(cnt-1)] are at the edges               */
/* MODIFIES p1[offset+tang*(cnt-1)].u .w                           */
/* ENSURES: if before=BEFORE: operate for the previous step values
	    if before=AFTER : finish the calculation of absorbing 
			      boundary condition                   */


{
   p1 += offset;
   pt += offset;
   while(cnt--)
   {
      if(before){
          if(tc(0)>0.00001 )
	      p1[0].u  = p1[norm].u +tc(0)*p1[0].u;
          if(te(0)>0.00001)
              p1[0].w  = p1[norm].w +te(0)*p1[0].w;
      }
      else {
          if(tc(0)>0.00001 )
              p1[0].u -= tc(0)*p1[norm].u;
          if(te(0)>0.00001)
              p1[0].w -= te(0)*p1[norm].w;
      }

      if(tc(0)<0.00001)p1[0].u = p1[norm].u;
      if(te(0)<0.00001)p1[0].w = p1[norm].w;

      p1 += tang;
      pt += tang;
   }
}

#define Pi               3.1415926
void AFilter(nh,filter)
int nh;
float *filter;

{
    int n;
    double cos();

    for(n=0; n<2*nh; n++)
        filter[n]=0.5*(1.0+cos(((double) (n)/ (double) (2*nh))*Pi));

}



kirrecord(struct vel_stress *p, int *ptable, float *na, int iz, 
	  float h, float dt, int nx, int nl, int nl_skip, int nh, int kdx, 
	  FILE *kirout, int source)

/* source = 0   -> P wave only (grad div u)x is recorded */ 
/* source = 1   -> S wave only (curl curl u)x is recorded */ 
/* source = 2   -> P and S waves are  recorded */ 

{
    register int i;
    register struct vel_stress *pp;
    register float *hold;
    register int *pt;
    int knx, offset;

    knx = (int) ((nx-nl -nh -3 +0.01)/kdx) -nl_skip;
    hold= (float *) malloc(sizeof(float)*(knx));

    if(source==0 || source ==2){    		
	for(i=0; i< knx; i++)		/* (gra div u)x */
	{
	    offset   = nl+1 + (i+1+nl_skip)*kdx +iz*nx;
	    pt       = ptable + offset;
	    pp       = p      + offset;
	    
	    hold[i]  = (pp[  +1].t11 - pp[   0].t11)* (9.0/8.0)
	              -(pp[  +2].t11 - pp[  -1].t11)* (1.0/24.0); 
	    hold[i] += (pp[  +1].t22 - pp[  -0].t22)* (9.0/8.0)
	              -(pp[  +2].t22 - pp[  -1].t22)* (1.0/24.0); 
	    hold[i] *=  9.0/8.0*0.5/(ta(0)+tb(0))*dt/h/h;          
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* z-derivative */
	{
	    offset   = nl+1 + (i+1+nl_skip)*kdx +iz*nx;
	    pt       = ptable + offset;
	    pp       = p      + offset;
	    
	    hold[i]  = (pp[ nx+1].t11 - pp[ nx-1].t11)
	              -(pp[-nx+1].t11 - pp[-nx-1].t11);
	    hold[i] += (pp[ nx+1].t22 - pp[ nx-1].t22)
	              -(pp[-nx+1].t22 - pp[-nx-1].t22);
	    hold[i] *=  9.0/8.0*0.125/(ta(0)+tb(0))*dt/h/h/h;          
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* (gra div u)z */
	{
	    offset   = nl+1 + (i+1+nl_skip)*kdx +iz*nx;
	    pt       = ptable + offset;
	    pp       = p      + offset;
	    
	    hold[i]  = (pp[  nx].t11 - pp[  0].t11)* (9.0/8.0)
	              -(pp[2*nx].t11 - pp[-nx].t11)* (1.0/24.0); 
	    hold[i] += (pp[  nx].t22 - pp[  0].t22)* (9.0/8.0)
	              -(pp[2*nx].t22 - pp[-nx].t22)* (1.0/24.0); 
	    hold[i] *=  9.0/8.0*0.5/(ta(0)+tb(0))*dt/h/h;          
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* z-derivative */
	{
	    offset   = nl+1 + (i+1+nl_skip)*kdx +iz*nx;
	    pt       = ptable + offset;
	    pp       = p      + offset;
	    
	    hold[i]  = (pp[nx].t11 - 2.0*pp[0].t11 + pp[-nx].t11)
	              +(pp[nx].t22 - 2.0*pp[0].t22 + pp[-nx].t22); 
	    hold[i] *=  9.0/8.0*0.5/(ta(0)+tb(0))*dt/h/h/h;          
	}
	fwrite(hold,sizeof(float),knx,kirout);
    }
    if(source==1 || source ==2){	
	for(i=0; i< knx; i++)		/* (curl curl u)x */
	{
	    offset   = nl+1 + (i+1+nl_skip)*kdx +iz*nx;
	    pt       = ptable + offset;
	    pp       = p      + offset;
	    
            hold[i]  = ( ( pp[1+nx].w-pp[-1+nx].w
                          -pp[1-nx].w+pp[-1-nx].w)/4.0
                        -( pp[nx].u-2.*pp[0].u+pp[-nx].u
                          +pp[2*nx].u-2.*pp[nx].u+pp[0].u
                          +pp[1+nx].u-2.*pp[1].u+pp[1-nx].u
                          +pp[1+2*nx].u-2.*pp[1+nx].u+pp[1].u)/4.0)/h/h; 
	    
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* z-derivative */
	{
	    offset   = nl+1 + (i+1+nl_skip)*kdx +iz*nx;
	    pt       = ptable + offset;
	    pp       = p      + offset;

	    hold[i] =  ( (  pp[1+nx].w-2.*pp[1].w+pp[1-nx].w 
	                  -(pp[nx-1].w-2.*pp[-1].w+pp[-1-nx].w)) 
	                -(  pp[2*nx].u-2.*pp[nx].u+pp[0].u + 
	                    pp[2*nx+1].u-2.*pp[nx+1].u+pp[1].u 
	                  -(pp[nx].u-2.*pp[0].u+pp[-nx].u +
	                    pp[nx+1].u-2.*pp[1].u+pp[-nx+1].u) ))/2./h/h/h;
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* (curl curl u)z */
	{
	    offset   = nl+1 + (i+1+nl_skip)*kdx +iz*nx;
	    pt       = ptable + offset;
	    pp       = p      + offset;

	    hold[i] = ( (pp[nx+1].u-pp[nx].u-pp[1].u+pp[0].u)/2.
	               -(pp[1].w-2.*pp[0].w+pp[-1].w))/h/h;
	    
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* z-derivative */
	{
	    offset   = nl+1 + (i+1+nl_skip)*kdx +iz*nx;
	    pt       = ptable + offset;
	    pp       = p      + offset;

            hold[i] = ( ((pp[2*nx+1].u-pp[2*nx].u-pp[nx+1].u+pp[nx].u)
                        -(pp[1].u-pp[0].u-pp[-nx+1].u+pp[-nx].u))/4. 
                       +((pp[1+nx].w-2.*pp[nx].w+pp[nx-1].w)
                        -(pp[1-nx].w-2.*pp[-nx].w+pp[-nx-1].w))/2.)/h/h/h;
	
	}
	fwrite(hold,sizeof(float),knx,kirout);
    }	
}

