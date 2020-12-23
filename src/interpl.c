#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void interpl(infile,rec_len,nst,rdt,mtotal,ns,sy,sy0,nx0,nz0)
FILE *infile; int nst, mtotal, *ns; float rdt, **sy, *sy0;
long int rec_len;
int nx0, nz0;

{
    float *weight, *green;
    int next, interp_const();
    int ns0;
    int i, j, size_f;
    int stress;
    float nst0;

    size_f = sizeof(float);

    green  = (float *) malloc(mtotal*sizeof(float));
    weight = (float *) malloc(4*sizeof(float));

    for(stress=0; stress < 2; stress++){
	nst0      = nst + stress * 0.5;
        next = interp_const(nst0,rdt,&ns0,ns,weight);

        if(next < 0){
    	    fread(green,size_f,mtotal,infile);
            for(i=0; i<mtotal; i++)
	        sy[i][ns[3]] = green[i];
        }

	if(stress == 0){
	    for(i=0; i<nx0; i++){
	        sy0[i] = 0.0;
	        for(j=0; j<4; j++)
	            sy0[i] += weight[j]*sy[i][ns[j]];
            }
	    for(i=3*nx0; i<4*nx0 + nz0; i++){
	        sy0[i] = 0.0;
	        for(j=0; j<4; j++)
	            sy0[i] += weight[j]*sy[i][ns[j]];
            }
	    for(i=4*nx0+3*nz0; i<4*nx0+4*nz0; i++){
	        sy0[i] = 0.0;
	        for(j=0; j<4; j++)
	            sy0[i] += weight[j]*sy[i][ns[j]];
            }
        } else {
	    for(i=nx0; i<3*nx0; i++){
	        sy0[i] = 0.0;
	        for(j=0; j<4; j++)
	            sy0[i] += weight[j]*sy[i][ns[j]];
            }
	    for(i=4*nx0+nz0; i<4*nx0+3*nz0; i++){
	        sy0[i] = 0.0;
	        for(j=0; j<4; j++)
	            sy0[i] += weight[j]*sy[i][ns[j]];
            }
       }
/*
        for(i=0; i<mtotal; i++){
	    sy0[i] = 0.0;
	    for(j=0; j<4; j++)
	        sy0[i] += weight[j]*sy[i][ns[j]];
        }
*/
    }
    free(green); free(weight);
}

int interp_const(nst,rdt,ns0,ns,weight)
float nst; int *ns0, *ns; float rdt;
float *weight;

{
    float t, dty1, dty2, dty3;
    int ns1; int next;
    int j;

    t     = nst*rdt;

    *ns0  = (int) (t+0.00001);

    next = 1;
    if(*ns0 > (int) ((nst-0.5)*rdt+0.00001)){
	ns1 = ns[0];
	for(j=0; j<3; j++)
	    ns[j] = ns[j+1];
        ns[3] = ns1;
	next = -1;
    }

    dty1      = t - (int) (t+0.00001);
    dty2      = dty1*dty1;
    dty3      = dty1*dty2;

    weight[0]  = -0.5*dty3 +     dty2 - 0.5 * dty1;
    weight[1]  =  3.0*dty3 - 2.5*dty2 + 1.0;
    weight[2]  = -3.0*dty3 + 2.0*dty2 + 0.5 * dty1;
    weight[3]  =  0.5*(dty3 - dty2);

    return next;
}


