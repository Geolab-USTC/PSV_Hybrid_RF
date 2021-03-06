
/******************************************************************/
/*                                                                */
/* GRT for P-SV system                                            */
/*                                                                */
/* Written by: Lianxing Wen, Seismological Lab. Caltech           */
/* Last modified: Sunday, Feb. 10, 1997                           */
/*                                                                */
/******************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "getpar.h"

#include "nrutil.h"

#define RAY_LENGTH      5000 
#define RAY_NUMBER      2000

#define DIRECT          0
#define REFLECTED       1
#define TOTAL           2
#define DOWN           -1
#define UP              1


int     nx;              /* x-dimension of mesh */
int     nz;              /* z-dimension of mesh */
int     nf         =10;  /* number of fluid meshs */
int     nl;              /* left region of mesh */ 
int     nt;              /* number of time steps */

float   h;               /* spatial mesh interval */
float   dp;              /* time digization interval */
float   xmin;            /* the left edge of the FD region */

/* filenames for rayfile and Green's function */
char rayfile[80];

/* following for the GRT */
float tstart = 0.0;
char raymodel[80], greenfile[80];

/* following for the Gauss sources */
float sgm, ts, *source, *dsource;
float theta, dip, lamda, azmuth;


main(ac,av)
int ac; char **av;

{   
    FILE *fop_model, *fop_green, *open_file();
    int n34, istress, nstress, isource; 
    int ispecial;
    int nrec, nfil, zlay, nfil0;
    int ix, iz, i, i0, mad1, n, j, k;
    float  xx, thickness, dep;
    int nso, nconv, l2n, nx1, ncent; float fscl;
    float *so, *dso;
    float *green; 

/* following for the model */
    int jo, nb, jo0, nb0;
    float *cc,  *ss,  *dd, *tth;
    float *cc0, *ss0, *dd0,*tth0;
/* following for the ray descriptions */
    int lfinal, lfinal0; 
    int *nen, *ncoun, *na, *nray; 
    int *nen0,*ncoun0,*na0,*nray0; 
    int bounce = 3;

    int iflat = 1;     /* =1 flat 0: no (spherical geometry)*/

    int horizon_profile = 1;  /* horizon profile is added */
    int vertical_profile = 1;  /* vertical profile is added */

    int grt_(), gauss_(), convt0_(), convt_();

    setpar(ac,av);
    mstpar("raymodel",  "s", raymodel);
    mstpar("greenfile", "s", greenfile);
    getpar("tstart",    "f", &tstart);

    mstpar("h",         "f", &h);
    mstpar("dp",        "f", &dp);
    mstpar("nt",        "d", &nt);
    mstpar("nx",        "d", &nx);
    mstpar("nf",        "d", &nf);
    mstpar("nl",        "d", &nl);
    mstpar("xmin",      "f", &xmin);
    
    getpar("flat",      "d", &iflat);
    getpar("bounce",    "d", &bounce);

    mstpar("nl",        "d", &nl);
    mstpar("sgm",       "f", &sgm);
    mstpar("ts",        "f", &ts);
    mstpar("theta",     "f", &theta);
    mstpar("dip",       "f", &dip);
    mstpar("lamda",     "f", &lamda);
    mstpar("azmuth",    "f", &azmuth);

    getpar("horizon_profile",  "d", &horizon_profile);
    getpar("vertical_profile", "d", &vertical_profile);
    endpar();

/* calculate the radiation patterns */
    so = (float *)malloc(5*sizeof(float));
    fault_(&theta,&dip,&lamda,&azmuth,so);

/* FFT the Gauss source function */
    nso = (int) (2.*ts/dp);
    source =(float *) malloc (nso*sizeof(float));
    gauss_(&sgm,&dp,&ts,source,&nso);

    dso =(float *) malloc (40200*sizeof(float));
    convt0_(&nt,&nso,source,dso,&dp,&nconv,&l2n,&nx1,&ncent,&fscl);
    free(source);

/* Read models from the file model */
    read_model(iflat,raymodel,&jo,&nb,&lfinal,&cc,&ss,&dd,&tth,
	       &nen,&ncoun,&na,&nray);
    for(n=0; n<lfinal; n++)	       
        ncoun[n]=((na[nen[0]-2]<=nb) ? UP : DOWN);

/* The layer number where n=3,4 are */
    nrec = na[nen[0]-2];

    fop_green = open_file(greenfile,"wb");

    green =(float *) malloc ((nt+3)*sizeof(float));


/* Calculate the initial responses at m=0 */

/* Nfil is the layer number where fluid or free surface is */
    nfil = nfinal(jo,ss,nrec);
    fprintf(stderr,"nfil = %d\n",nfil);
    
/* Evaluate the thickness of the FD region; modify the thickness 
   of the last layer, or order to fit the actual region into the 
   FD grids, nz is the number of total HALF grids in the solid 
*/

    thickness=0.0;
    if(nfil>nrec){
        for(n=nrec; n<nfil-1; n++)
	    thickness += tth[n];
    } else {	    
        for(n=nfil; n<nrec-1; n++)
	    thickness += tth[n];
    }	    

    nz = (int) ((thickness+1.e-3)/(0.5*h)); 
    n =  ((nfil>nrec) ? nfil-2 : nfil);
    tth[n]= 0.5*nz*h- (thickness-tth[n])+0.0001+0.25*h;
    thickness = nz*0.5*h + 0.25 *h;

/* This set of points are derived from the initial model and rays */

    cc0   = (float *) malloc ((jo+1)*sizeof(float));
    ss0   = (float *) malloc ((jo+1)*sizeof(float));
    dd0   = (float *) malloc ((jo+1)*sizeof(float));
    tth0  = (float *) malloc ((jo+1)*sizeof(float));

    nen0   = (int *) malloc (RAY_NUMBER * sizeof(int));
    ncoun0 = (int *) malloc (RAY_NUMBER * sizeof(int));
    na0    = (int *) malloc (RAY_LENGTH * sizeof(int));
    nray0  = (int *) malloc (RAY_LENGTH * sizeof(int));

    mad1=nx+nl;
    isource=1;
    nz = (nz+1+2*nf)/2-1;
    iflat = 1;   /* take care of GRT calculation:model has been flattened*/

/*  fwrite(&nz,sizeof(int),1,fop_green);   */ 

/* Calculate the initial responses at n=3 and 4 */
    if(horizon_profile){
        ispecial = 0;
        for(n34=3; n34<=3; n34++){

            dep =(n34-3)*0.5*h+1.e-7;
            zlay = modify_model(dep,nrec,nfil,jo,cc,ss,dd,tth,
                            &jo0,cc0,ss0,dd0,tth0);
                                   
            nfil0 = nfinal(jo0,ss0,nrec);

    
            for(istress=0; istress<=0; istress++){

/*              nstress =0        w   velocity 
                nstress =12       T12 stress 
                nstress =1        u   velocity
                nstress =22       T22 stress 
*/

                if(n34==4 && istress==0)nstress=12;
                if(n34==4 && istress==1)nstress=0;
                if(n34==3 && istress==0)nstress=1;
                if(n34==3 && istress==1)nstress=22;
        
                i0=(istress==0) ? nl+1 : nl; 	    
            
	        printf("Getting the Rays for REFLECTED FOR LEFT REGION\n");
                getrays(lfinal,nen,na,nray,ncoun,zlay,nrec,nfil0,jo,nb, 
         	    &lfinal0,nen0,ncoun0,na0,nray0,jo0,&nb0,bounce,REFLECTED);
	        printf("End Getting the Rays for REFLECTED FOR LEFT REGION\n");
         	    
                for(i=0; i<i0; i++){	    
            
		    printf("REFLECTED i=%d",i);
                    xx=xmin+(i-nl)*h+0.5*istress*h;
                
                    grt_(&isource,&iflat,&ispecial,&nb0,&jo0,cc0,ss0,dd0,tth0,
		      &lfinal0,nen0,na0,nray0,ncoun0,
                      &dp,&nt,&xx,&tstart,&nstress,so,green);
                    convt_(green,&nt,green,dso,&nconv,&l2n,&nx1,
                       &ncent,&fscl);  

                    for(j=0;j<nt;j++) green[j] *= (-1);
                      
                    fwrite(green,sizeof(float),nt,fop_green);     

                }                

	        printf("Getting the Rays for DIRECT for right REGION\n");
                getrays(lfinal,nen,na,nray,ncoun,zlay,nrec,nfil0,jo,nb, 
         	    &lfinal0,nen0,ncoun0,na0,nray0,jo0,&nb0,bounce,DIRECT);
	        printf("End Getting the Rays for DIRECT for right REGION\n");

                for(i=i0; i<mad1; i++){
		    printf("DIRECT i=%d",i);
            
                    xx=xmin+(i-nl)*h+0.5*istress*h;
        
                    grt_(&isource,&iflat,&ispecial,&nb0,&jo0,cc0,ss0,dd0,tth0,
		     &lfinal0,nen0,na0,nray0,ncoun0,
                     &dp,&nt,&xx,&tstart,&nstress,so,green);
                    convt_(green,&nt,green,dso,&nconv,&l2n,&nx1,
                       &ncent,&fscl);  

                     
                    fwrite(green,sizeof(float),nt,fop_green);     
                }
            }
        }
    }


    if(vertical_profile){
    for(n34=3; n34<=3; n34++){

        xx = ( (n34==3) ? xmin : xmin+0.5*h);

	for(istress=0; istress<=0; istress++){

	        
/*          nstress =0        w   velocity 
            nstress =12       T12 stress 
            nstress =1        u   velocity
            nstress =11       T11 stress 
*/
            if(istress==1 && n34==4) nstress =0;
            if(istress==0 && n34==3) nstress =1;
            if(istress==0 && n34==4) nstress =11;
            if(istress==1 && n34==3) nstress =12;

            for(iz=0; iz<=nz; iz++){

                printf("Total iz=%d", iz);
                dep =iz*h+istress*0.5*h+1.e-7;
		ispecial = (( dep > thickness) ? 0 : 0 );
                zlay = modify_model(dep,nrec,nfil,jo,cc,ss,dd,tth,
                                   &jo0,cc0,ss0,dd0,tth0);
                                   

                nfil0 = nfinal(jo0,ss0);

                getrays(lfinal,nen,na,nray,ncoun,zlay,nrec,nfil0,jo,nb, 
	                &lfinal0,nen0,ncoun0,na0,nray0,jo0,&nb0,bounce,TOTAL);
	                
                grt_(&isource,&iflat,&ispecial,&nb0,&jo0,cc0,ss0,dd0,tth0,
		     &lfinal0,nen0,na0,nray0,ncoun0,
                     &dp,&nt,&xx,&tstart,&nstress,so,green);
                convt_(green,&nt,green,dso,&nconv,&l2n,&nx1,
                       &ncent,&fscl);  

                fwrite(green,sizeof(float),nt,fop_green);     
            }
        }
    }
    }
    
    fclose(fop_green);
}    

/* getrays determins the ray parameters for the responses 
   at zlay layer. The input is the ray group for the n=3,4, 
       which are: lfinal, nen, na, nray 
                  zlay: the layer where the receive is
                  nrec: the layer where n=3,4 lines are
                  nfil: the layer where solid-liquid interface is
                  type =0: direct wavefield
                       =1: reflected wavefield
                       =2: total wavefield
*/                  

getrays(lfinal,nen,na,nray,ncoun,zlay,nrec,nfil,jo,nb,
        lfinal0,nen0,ncoun0,na0,nray0,jo0,nb0,bounce,raytype)
int lfinal, *nen, *na, *nray, *ncoun; 
int zlay, nrec, nfil, bounce, raytype, jo, nb, jo0, *nb0;
int *lfinal0, *nen0, *na0, *nray0, *ncoun0;

{
    int ray, k, n, j, rays1, rays2, nenn, rays, layer2, incre=0;
    int P_SV, SV_P;

    int **rtype, *type;
    int nbounce, i, m, num;

    *nb0=nb;
    if(nfil<nrec+1 && jo != jo0){
        *nb0=nb+1; nrec++; incre++;
    }

    n=k=0;
    for (ray=0; ray<lfinal; ray++){
        
        P_SV = nray[nen[ray]-2];
	SV_P = ( (P_SV==3) ? 5 : 3);

	switch(raytype){

	case DIRECT:
	    rays1=0; rays2=0; break;
	case REFLECTED:     
	    rays1=1; rays2=1; break;
	case TOTAL:     
	    rays1=0; rays2=1; break;

        }

/*	if(zlay==nfil) rays2=0;  */

	rtype = imatrix(0,bounce,1,power(2,bounce+2));
	type  = ivector(0,bounce);

        for(rays=rays1; rays<=rays2; rays++){

            if(rays == 0){
                printf("DIRECT WAVE\n");
                nenn=0;

                /* keep the ray paths to n=3, 4 */
                for(j=0; j<nen[ray]-1; j++){
                    na0[k]   =na[j]+incre;    
                    nray0[k] =nray[j];    
                    k++;
                    nenn++;
                } 

                /* calculate the direct wave from the source */
		nbounce  = 1;
		type[0] = P_SV;
		nen0[n] = nenn;

		k  = MultiBounces(nrec,nfil,zlay,k,nbounce,type,
				   na0,nray0,&ncoun0[n],&nen0[n]);

                 
                n++;
            }     
            else { 
            
                /* reflected waves from the solid-liquid boundary */
                for(nbounce=2; nbounce<=bounce; nbounce++){

		    if((zlay == nfil) && ((nbounce % 2) == 0)) 
			continue;     /* skip the uparriving ray */
		    if((zlay == nrec) && ((nbounce % 2) != 0)) 
			continue;     /* skip the downarriving ray */
		     
                    num = BounceRays(nbounce,rtype,P_SV,SV_P);

		    for(j=1; j<=num; j++){

                        nenn=0;

                        /* keep the ray paths to n=3, 4 */
                       for(m=0; m<nen[ray]-1; m++){
                           na0[k]   =na[m]+incre;    
                           nray0[k] =nray[m];    
                           k++;
                           nenn++;
                        } 

			for(i=0; i<nbounce; i++){
			    type[i]  = rtype[i][j];
			    printf("%d ",type[i]);     
                         }
			 printf("\n"); 
		    
		        nen0[n] = nenn;
		        k  = MultiBounces(nrec,nfil,zlay,k,nbounce,type,
				   na0,nray0,&ncoun0[n],&nen0[n]);


                        n++;
                    }
                }
            }            
        }
	free_imatrix(rtype,0,bounce,1,power(2,bounce+2));
	free_ivector(type,0,bounce);
    }

    *lfinal0=n;
}            

/*
    int BounceRays(int nbounce, int ** rtype, P_SV SV_P

    rtype[i][j]  i -> ray segments j -> number of rays

*/

int BounceRays(nbounce,rtype,P_SV,SV_P)
int nbounce, P_SV, SV_P;
int **rtype;
{
    int i, j, m, num;

    /* first fundamental ray */ 
    for(i=0; i<nbounce; i++)
         rtype[i][1] = P_SV;

    num = 1;
    for(i=1; i< nbounce; i++){

        for(j=1; j<=num; j++){

	    for(m=0; m<nbounce; m++)
	        rtype[m][j+num] = rtype[m][j];

            rtype[nbounce - i][j+num] = SV_P;
        }
	num *= 2;

/*
        for(j=1; j<=num; j++){ 
	    printf("j=%d ", j);
	    for(m=0; m<nbounce; m++)
		printf("%d ",rtype[m][j]);
	    printf("\n");
        }
*/

    }

    return num;
}

int MultiBounces(nrec,nfil,zlay,k,bounce,type,na,nray,ncoun,nen)
int nrec, nfil, zlay; int k, bounce, *type;
int *na, *nray, *ncoun, *nen;
{
    int j, i; int start, end, temp; int going;
    int nen0, k0;


    start = ((nfil > nrec) ? nrec + 1 : nrec);
    end   = ((nfil > nrec) ? nfil - 1 : nfil +1);
    going = ((nfil > nrec) ? DOWN     : UP);


    nen0 = *nen;  k0 = k;
    i = 0;
    while ( i < bounce){

	if(i == bounce-1)
	    /* last segement of the ray */
	    end  = ((going == UP) ? zlay+1 : zlay);

        if(going == DOWN){
	    for(j=start; j<=end; j++){
	        na[k0]   = j;
	        nray[k0] = type[i];
		if((nfil > nrec) && (j >= nfil))
		    nray[k0] = 5;    /* fixed to P wave */
		if((nfil < nrec) && (j <= nfil))
		    nray[k0] = 5;    /* fixed to P wave */
	        k0++; nen0++;
            }
	    going = UP;
        } else { 
	    for(j=start; j>=end; j--){
	        na[k0]   = j;
	        nray[k0] = type[i];
		if((nfil > nrec) && (j >= nfil))
		    nray[k0] = 5;    /* fixed to P wave */
		if((nfil < nrec) && (j <= nfil))
		    nray[k0] = 5;    /* fixed to P wave */
	        k0++; nen0++;
            }
	    going = DOWN;
        }
	i++;
	temp = start; start = end; end = temp; 
    }
    na[k0]   = 1;
    nray[k0] = 5; 
    k0++; nen0++;
    *ncoun  = -going;
    *nen    =  nen0;
    return k0;
}


/* modify the model parameters for each point at m=0*/
/* returns the layer number of the receiver */ 

int modify_model(dep,nrec,nfil,jo,cc,ss,dd,tth,
                 jo0,cc0,ss0,dd0,tth0)
float dep; int nrec, nfil; 
int jo; float *cc, *ss, *dd, *tth;
int *jo0; float *cc0, *ss0, *dd0, *tth0;

{   
    int layer, zlay;
    float thickness=0.0;

    if(nrec<nfil){
        for(layer=nrec; layer<nfil-1; layer++){
	    thickness += tth[layer];
	    if(thickness >= dep) break;
        }
        zlay = ( (dep<0.00001) ? nrec : layer+1);
    } else {
        for(layer=nrec-2; layer>nfil-1; layer--){
	    thickness += tth[layer];
	    if(thickness >= dep) break;
        }
        zlay = ( (dep<0.00001) ? nrec-1 : layer+1);
    }

    if(dep > thickness)zlay = nfil;

    if(fabs(thickness-dep) <0.0001 || zlay==nrec || zlay==nfil){
        for(layer=0; layer<jo; layer++){
            cc0[layer]  =  cc[layer];
            ss0[layer]  =  ss[layer];
            dd0[layer]  =  dd[layer];
            tth0[layer] = tth[layer];
        }            
	if(zlay==nfil)tth0[nfil-1]=dep-thickness;
        *jo0 =jo; 
    }    
    else {
        for(layer=0; layer<zlay; layer++){
            cc0[layer]  =  cc[layer];
            ss0[layer]  =  ss[layer];
            dd0[layer]  =  dd[layer];
            tth0[layer] = tth[layer];
        }
        
    
         cc0[zlay]  =  cc[zlay-1]+0.00001; 
         ss0[zlay]  =  ss[zlay-1]+0.00001; 
         dd0[zlay]  =  dd[zlay-1]+0.00001; 
         tth0[zlay]   = (nrec<nfil) ? thickness-dep : tth[zlay-1]-(thickness-dep); 
         tth0[zlay-1] = (nrec>nfil) ? thickness-dep : tth[zlay-1]-(thickness-dep); 

         for(layer=zlay+1; layer<jo+1; layer++){
            cc0[layer]  =  cc[layer-1];
            ss0[layer]  =  ss[layer-1];
            dd0[layer]  =  dd[layer-1];
            tth0[layer] = tth[layer-1];
        }

        *jo0 =jo+1;
    }     

    return zlay;
}

int power(int a, int n)
{
    int j, num;

    num = a;
    for(j=2; j<=n; j++)
	num *= a;

    return num;
}



