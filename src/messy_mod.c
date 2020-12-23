#include	<stdlib.h>
#include	<stdio.h>
#include	<sys/file.h>
#include	<fcntl.h>
#include	<math.h>

/*
	messy_mod takes an input file (normally named "run".gen,
	where "run" is the code name for the FD run).  From this,
	it constructs three files, pvelocity, svelocity, and density.
	The first line of the input file is the background medium
	specified as pvelocity, svelocity, density.

	In the input file, you may specify homogeneous polygons:
		nvert pvel svel den ipat
			x0     z0
			x1     z1
			 .      .
			 .      .
			xnvert znvert
			x0     z0
	where xn,zn are the x,z coordinates of the vertices; ipat
	is a pattern designator for other programs which draw the
	model; and nvert is the number of vertices (2 < nvert < 50).

	You may also specify polygons with gradients:
		-nvert dum dum dum ipat
			x0     z0     p0     s0     d0
			x1     z1     p0     s0     d0
			 .      .      .      .      .
			 .      .      .      .      .
			xnvert znvert pnvert snvert dnvert
			x0     z0     p0     s0     d0
	xn, zn, ipat are as above, nvert must be preceded by a '-' (this
	is used as a flag in the program); dum indicates a dummy float
	which must be present; and pn, sn, dn are the pvelocity, svelocity
	and density at the vertex.  (3 <= nvert <= 4; if more vertices
	are required, break region into 3 and 4 vertex polys).

	A final option is a random medium within a polygon of an
	arbitrary number of vertices.  Such a region is specified by
	a negative ipat.  The line immediately thereafter contains
	statistical properties of the region.  The medium specified
	in the first line is the mean for the distribution.
		nvert pvel svel den -ipat
		variation scalex scalez seed
			x0     z0
			x1     z1
			 .      .
			 .      .
			xnvert znvert
			x0     z0
	where variation is the RMS variation of the distribution, scale
	is the scale length of inhomogeneities, and seed is the initial
	seed for generating random values.  The P-velocity field is
	generated using these parameters, and the S-velocity field is
	generated from the P-velocity field using a constant ratio as
	determined by the ratio implicit in the first line.  The density
	field only receives 30% of the variation.

	Required parameters:
		model="run".grad
		pvelocity="run".pvel
		svelocity="run".svel
		density="run".den
		nnx=nnx
		nnz=nnz
	Optional parameter:
		shift=1 (default is 0):
			this parameter tells grademodel to shift
			coordinates by half a grid point before
			calculating interpolated media parameters,
			allowing exact coorespondence to the offset
			grid scheme used to fill the polygons.
*/

double ap, bp, cp, dp,
       as, bs, cs, ds,
       ad, bd, cd, dd;

#define Pfunc(X,Z) ( dp * X * Z + cp * Z + bp * X + ap )
#define Sfunc(X,Z) ( ds * X * Z + cs * Z + bs * X + as )
#define Dfunc(X,Z) ( dd * X * Z + cd * Z + bd * X + ad )

#define F (float)

struct vertex
{
  int xv, zv;
  struct vertex *f;
  struct vertex *b;
};

struct gradvertex
{
  int xv, zv;
  float pv, sv, dv;
  struct gradvertex *f;
  struct gradvertex *b;
};


struct gradvertex *gradbase;
struct vertex *base;

struct complex *buf;

void cfft(), cfft2d(), getrand(), gauss_mod();

int *crosses;
int shift = 0;

int messy_mod(int nnx, int nnz, float *pvel, float *svel, float *density, char *model)
{
float *a, *b, *c;
  struct vertex *base, *head, *p;
  struct gradvertex *gradbase, *ghead, *gp;
  FILE  *af, *bf, *cf, *open_file();
  FILE *mf, *fopen();
  float aback,bback,cback,aval,bval,cval;
  int ipat, seed, xvmax, xvmin, xdim, zdim;
  float rms_var, alphax, alphaz;
  int x, z;
  float p1, s1, d1;
  float *pa, *pb, *pc;
  float *ra, *rb, *rc;
  float *fa, *fb, *fc;
  int npts, xstr, xend, i;
  int zvmin, zvmax, ncross;
  int n;
  
  gradbase = (struct gradvertex *) malloc (20000*sizeof(struct gradvertex));
  base     = (struct vertex *) malloc (20000*sizeof(struct vertex));
  crosses  = (int *) malloc (40000*sizeof(int));
  a = pvel; b = svel; c = density;
  if( (mf= fopen(model,"r")) == NULL)
    {
      fprintf(stderr,"cannot open model=%s\n",model);
      exit(-1);
    }
  
  n=fscanf(mf,"%f %f %f",&aback,&bback,&cback);
  for(i=0; i< (nnx*nnz); i++)
    {
      a[i]= aback;
      b[i]= bback;
      c[i]= cback;
    }
  
  while( (n=fscanf(mf,"%d %f %f %f %d",&npts,&aval,&bval,&cval,&ipat)) != -1 )
    {
      if (npts < 0) {	/* gradient poly */
	npts = -npts;
	fprintf(stderr,"grad poly: %d vertices\n",npts);
	if(npts < 3)
	  {
	    fprintf(stderr,"polynomial with < 3 vertices\n");
	    exit(-1);
	  }
	if(npts > 4)
	  {
	    fprintf(stderr,"gradient polynomial with > 4 vertices\n");
	    fprintf(stderr,"not supported, please break up into 3 and 4 vertex polys\n");
	    exit(-1);
	  }
	
	/* Set up a double linked ring of vertices. The pointer
	 * to the ring (head) is set to the node with the maximum
	 * x-value so that intersect() will not eliminate 'head'
	 * while casting off vertices.
	 */
	zvmin= 2000000; zvmax= -2000000;
	ghead= gradbase;
	for(i=0; i<npts; i++)
	  {
	    gp= &gradbase[i];
	    gp->f= gp+1; gp->b= gp-1;
	    fscanf(mf,"%d %d %f %f %f",&x,&z,&p1,&s1,&d1);
	    gp->xv= x; gp->zv= z;
	    gp->pv= p1; gp->sv= s1; gp->dv= d1;
	    if(gp->zv > zvmax) { zvmax= gp->zv; ghead= gp; }
	    if(gp->zv < zvmin)   zvmin= gp->zv;
	  }
	gp->f= gradbase; gradbase->b= gp; /* close the ring */
	fscanf(mf,"%d %d %f %f %f",&x,&z,&p1,&s1,&d1);
	/* check for closure */
	if(gradbase->xv != x || gradbase->zv !=  z ||
	   gradbase->pv != p1 || gradbase->sv != s1 ||
	   gradbase->dv != d1)
	  {
	    fprintf(stderr,"grad poly does not close\n");
	    fprintf(stderr,"first vertex: %d %d %f %f %f\n",
		    gradbase->xv,gradbase->zv,
		    gradbase->pv,gradbase->sv,gradbase->dv);
	    fprintf(stderr,"last  vertex: %d %d %f %f %f\n",
		    x,z,p1,s1,d1);
	    exit(-1);
	  }
	interpol(gradbase, npts);
	if(zvmax > (nnz-1)) zvmax= (nnz-1);
	if(zvmin < 0 ) zvmin= 0;
	gp= ghead;
	do gp->zv = 2* gp->zv + 1; while( (gp=gp->f) != ghead);
	for(z= zvmin; z <= zvmax; z++)
	  {
	    ncross= gintersect(2*z,crosses,ghead);
	    sort(crosses,ncross);
	    pa= a + z*nnx;
	    pb= b + z*nnx;
	    pc= c + z*nnx;
	    for(i=0; i< ncross; i += 2)
	      {
		xstr= crosses[i];
		if (xstr<0) xstr=0;
		xend= crosses[i+1];
		if (xend>=nnx) xend = nnx-1;
		for(x=xstr; x<=xend; x++)
		  {
		    pa[x] = Pfunc(x,z);
		    pb[x] = Sfunc(x,z);
		    pc[x] = Dfunc(x,z);
		  }
	      }
	  }
      }
      else if (ipat < 0) {	/* random poly */
	fprintf(stderr,"random poly: %d vertices\n",npts);
	if(npts < 3)
	  {
	    fprintf(stderr,"polynomial with < 3 vertices");
	    exit(-1);
	  }
	fscanf(mf,"%f %f %f %d",&rms_var,&alphax,&alphaz,&seed);
	
	/* Set up a double linked ring of vertices. The pointer
	 * to the ring (head) is set to the node with the maximum
	 * x-value so that intersect() will not eliminate 'head'
	 * while casting off vertices.
	 */
	zvmin= 2000000; zvmax= -2000000;
	xvmin= 2000000; xvmax= -2000000;
	head= base;
	for(i=0; i<npts; i++)
	  {
	    p= &base[i];
	    p->f= p+1; p->b= p-1;
	    fscanf(mf,"%d %d",&x,&z);
	    p->xv= x; p->zv= z;
	    if(p->zv > zvmax) { zvmax= p->zv; head= p; }
	    if(p->zv < zvmin)   zvmin= p->zv;
	    if(p->xv > xvmax)   xvmax= p->xv;
	    if(p->xv < xvmin)   xvmin= p->xv;
	  }
	p->f= base; base->b= p; /* close the ring */
	fscanf(mf,"%d %d",&x,&z);
	/* check for closure */
	if(base->xv != x || base->zv !=  z)
	  {
	    fprintf(stderr,"polygon does not close");
	    fprintf(stderr,"first vertex: %d %d\n",
		    base->xv,base->zv);
	    fprintf(stderr,"last  vertex: %d %d\n", x,z);
	    exit(-1);
	  }
	if(zvmax > (nnz-1)) zvmax= (nnz-1);
	if(zvmin < 0 ) zvmin= 0;
	xdim = xvmax - xvmin;
	xdim++;
	xdim = pow2(xdim);
	zdim = zvmax - zvmin;
	zdim++;
	zdim = pow2(zdim);
	fa = (float *) malloc(sizeof(float)*xdim*zdim);
	fb = (float *) malloc(sizeof(float)*xdim*zdim);
	fc = (float *) malloc(sizeof(float)*xdim*zdim);
        ra = fa; rb = fb; rc = fc;
	gauss_mod(xdim,zdim,seed,rms_var,alphax,alphaz,aval,bval,cval,ra,rb,rc);
	ra -= xvmin;
	rb -= xvmin;
	rc -= xvmin;
	p= head;
	do p->zv = 2* p->zv + 1; while( (p=p->f) != head);
	for(z= zvmin; z <= zvmax; z++)
	  {
	    ncross= intersect(2*z,crosses,head);
	    sort(crosses,ncross);
	    pa= a + z*nnx;
	    pb= b + z*nnx;
	    pc= c + z*nnx;
	    for(i=0; i< ncross; i += 2)
	      {
		xstr= crosses[i];
		if (xstr<0) xstr=0;
		xend= crosses[i+1];
		if (xend>=nnx) xend = nnx-1;
		for(x=xstr; x<=xend; x++)
		  {
		    pa[x] = aback + (ra[x]-aback);
		    pb[x] = bback + (rb[x]-bback);
		    pc[x] = cback + (rc[x]-cback);
		  }
	      }
	    ra += xdim;
	    rb += xdim;
	    rc += xdim;
	  }
	free(fa);
	free(fb);
	free(fc);
      }
      else {		/* homogeneous poly */
	fprintf(stderr,"uniform poly: %d vertices\n",npts);
	if(npts < 3)
	  {
	    fprintf(stderr,"polynomial with < 3 vertices");
	    exit(-1);
	  }
	
	/* Set up a double linked ring of vertices. The pointer
	 * to the ring (head) is set to the node with the maximum
	 * x-value so that intersect() will not eliminate 'head'
	 * while casting off vertices.
	 */
	zvmin= 2000000; zvmax= -2000000;
	head= base;
	for(i=0; i<npts; i++)
	  {
	    p= &base[i];
	    p->f= p+1; p->b= p-1;
	    fscanf(mf,"%d %d",&x,&z);
	    p->xv= x; p->zv= z;
	    if(p->zv > zvmax) { zvmax= p->zv; head= p; }
	    if(p->zv < zvmin)   zvmin= p->zv;
	  }
	p->f= base; base->b= p; /* close the ring */
	fscanf(mf,"%d %d",&x,&z);
	/* check for closure */
	if(base->xv != x || base->zv !=  z)
	  {
	    fprintf(stderr,"polygon does not close");
	    fprintf(stderr,"first vertex: %d %d\n",
		    base->xv,base->zv);
	    fprintf(stderr,"last  vertex: %d %d\n", x,z);
	    exit(-1);
	  }
	if(zvmax > (nnz-1)) zvmax= (nnz-1);
	if(zvmin < 0 ) zvmin= 0;
	p= head;
	do p->zv = 2* p->zv + 1; while( (p=p->f) != head);
	for(z= zvmin; z <= zvmax; z++)
	  {
	    ncross= intersect(2*z,crosses,head);
	    sort(crosses,ncross);
	    pa= a + z*nnx;
	    pb= b + z*nnx;
	    pc= c + z*nnx;
	    for(i=0; i< ncross; i += 2)
	      {
		xstr= crosses[i];
		if (xstr<0) xstr=0;
		xend= crosses[i+1];
		if (xend>=nnx) xend = nnx-1;
		for(x=xstr; x<=xend; x++)
		  {
		    pa[x] = aval;
		    pb[x] = bval;
		    pc[x] = cval;
		  }
	      }
	  }
      }
    }
  for(x=0; x<nnx; x++){
      for(z=0; z<nnz; z++)
	  fprintf(stdout,"%f ",a[nnx*z+x]);
      fprintf(stdout,"\n");
  }
  free(base); free(gradbase); free(crosses);
  fclose(mf);
  return 1;
}

intersect(z,crosses,head)
register int z; register int *crosses; struct vertex *head;
{
  register int ncross, x; register struct vertex *p;
  
  ncross= 0;
  p= head;
  do
    {
      if(p->zv > z && p->b->zv > z) continue;
      if(p->zv < z && p->b->zv < z)
	{
	  if(p->f->zv > z) continue;
	  /* eliminate vertex */
	  p->b->f= p->f;
	  p->f->b= p->b;
	  continue;
	}
      x= solve(z,p->zv,p->xv,p->b->zv,p->b->xv);
      crosses[ncross++]= x;
    }	while( (p=p->f) != head );
  return(ncross);
}

gintersect(z,crosses,head)
register int z; register int *crosses; struct gradvertex *head;
{
  register int ncross, x; register struct gradvertex *p;

  ncross= 0;
  p= head;
  do
    {
      if(p->zv > z && p->b->zv > z) continue;
      if(p->zv < z && p->b->zv < z)
	{
	  if(p->f->zv > z) continue;
	  /* eliminate vertex */
	  p->b->f= p->f;
	  p->f->b= p->b;
	  continue;
	}
      x= solve(z,p->zv,p->xv,p->b->zv,p->b->xv);
      crosses[ncross++]= x;
    }	while( (p=p->f) != head );
  return(ncross);
}

sort(vec,n)
register int *vec; register int n;
{
  /* sort the elements of vec into ascending order. */
  register int *above, *below, *last;
  register int temp;
  last = vec + (n-1);
  for(above=vec; above<last; above++)
    {
      for(below=above+1; below<=last; below++)
	{
	  if(*above > *below)
	    {
	      temp=   *above;
	      *above= *below;
	      *below= temp;
	    }
	}
    }
}

solve(pnot,p1,q1,p2,q2)
register int pnot,p1,q1,p2,q2;
{
  /* floating point version */
  double invslope;
  register int qnot;
  if(pnot==p1) return(q1);
  if(pnot==p2) return(q2);
  if(q1==q2) return(q1);
  invslope= (q1-q2)/( (double) (p1-p2));
  qnot= (pnot-p1)*invslope + (double) q1 + 0.5;
  return(qnot);
}

interpol(verts,nvert)
/*
	fits a polynomial surface to specified values of
	p-velocity, s-velocity and density (in turn) for
	three (plane) or four specified points.  Uses
	gaussian elimination - partial pivoting subroutines.
*/
struct gradvertex *verts;
int nvert;
{
  double **a, b[4];
  int *rows;
  double **init_pivot();
  int solve_pivot();
  void pivot();
  void backsub();
  void reaugment();
  
  switch (nvert) {
    
  case 3:
    /*
      Form of equation:
      | a00 a01 a02 |   | ap |   | a03 |
      | a10 a11 a12 | * | bp | = | a13 |
      | a20 a21 a22 |   | cp |   | a23 |
      Equation derived from
      ap + bp * xv + cp * zv + dp * xv * zv = pv,
      one equation for each vertex.  (For just three
      vertices, dp = 0.)  an0 = 1.  an3 are the media
      parameters specified at each grid point.
      */
    
    a = init_pivot(nvert,&rows);
    a[0][0] = 1; a[0][1] = verts[0].xv;
    a[0][2] = verts[0].zv + shift * 0.5;
    a[1][0] = 1; a[1][1] = verts[1].xv;
    a[1][2] = verts[1].zv + shift * 0.5;
    a[2][0] = 1; a[2][1] = verts[2].xv;
    a[2][2] = verts[2].zv + shift * 0.5;
    /* do p-velocity first */
    a[0][3] = verts[0].pv;
    a[1][3] = verts[1].pv;
    a[2][3] = verts[2].pv;
    if (solve_pivot(a, rows, nvert) < 0) {
      fprintf(stderr,"Singular matrix, check input vertices!\n");
      exit(-1);
    }
    backsub(a, nvert);
    ap = a[0][3];
    bp = a[1][3];
    cp = a[2][3];
    dp = 0.0;
    /* now do s-velocity using reaugmentation of a */
    b[0] = verts[0].sv;
    b[1] = verts[1].sv;
    b[2] = verts[2].sv;
    reaugment(a, b, rows, nvert);
    backsub(a, nvert);
    as = a[0][3];
    bs = a[1][3];
    cs = a[2][3];
    ds = 0.0;
    /* finally, do density */
    b[0] = verts[0].dv;
    b[1] = verts[1].dv;
    b[2] = verts[2].dv;
    reaugment(a, b, rows, nvert);
    backsub(a, nvert);
    ad = a[0][3];
    bd = a[1][3];
    cd = a[2][3];
    dd = 0.0;
    break;
	
  case 4:
    /*
      Form of equation:
      | a00 a01 a02 a03 |   | ap |   | a04 |
      | a10 a11 a12 a13 | * | bp | = | a14 |
      | a20 a21 a22 a23 |   | cp |   | a24 |
      | a30 a31 a32 a33 |   | cp |   | a24 |
      Equation derived from
      ap + bp * xv + cp * zv + dp * xv * zv = pv,
      one equation for each vertex.  (For just three
      vertices, dp = 0.)  an0 = 1.  an4 are the media
      parameters specified at each grid point.
      */
    
    a = init_pivot(nvert,&rows);
    a[0][0] = 1; a[0][1] = verts[0].xv;
    a[0][2] = verts[0].zv + shift * 0.5;
    a[0][3] = verts[0].xv * a[0][2];
    a[1][0] = 1; a[1][1] = verts[1].xv;
    a[1][2] = verts[1].zv + shift * 0.5;
    a[1][3] = verts[1].xv * a[1][2];
    a[2][0] = 1; a[2][1] = verts[2].xv;
    a[2][2] = verts[2].zv + shift * 0.5;
    a[2][3] = verts[2].xv * a[2][2];
    a[3][0] = 1; a[3][1] = verts[3].xv;
    a[3][2] = verts[3].zv + shift * 0.5;
    a[3][3] = verts[3].xv * a[3][2];
    /* do p-velocity first */
    a[0][4] = verts[0].pv;
    a[1][4] = verts[1].pv;
    a[2][4] = verts[2].pv;
    a[3][4] = verts[3].pv;
    if (solve_pivot(a, rows, nvert) < 0) {
      fprintf(stderr,"Singular matrix, check input vertices!\n");
      exit(-1);
    }
    backsub(a, nvert);
    ap = a[0][4];
    bp = a[1][4];
    cp = a[2][4];
    dp = a[3][4];
    /* now do s-velocity using reaugmentation of a */
    b[0] = verts[0].sv;
    b[1] = verts[1].sv;
    b[2] = verts[2].sv;
    b[3] = verts[3].sv;
    reaugment(a, b, rows, nvert);
    backsub(a, nvert);
    as = a[0][4];
    bs = a[1][4];
    cs = a[2][4];
    ds = a[3][4];
    /* finally, do density */
    b[0] = verts[0].dv;
    b[1] = verts[1].dv;
    b[2] = verts[2].dv;
    b[3] = verts[3].dv;
    reaugment(a, b, rows, nvert);
    backsub(a, nvert);
    ad = a[0][4];
    bd = a[1][4];
    cd = a[2][4];
    dd = a[3][4];
    break;
	
  default:
    fprintf(stderr,"Error in interpol:  %d vertices unsupported\n",
	    nvert);
    exit(-1);
  }
}

#define		PI		F 3.1415926536

struct complex { float re, im; };

void
gauss_mod(nnx,ny,seed,amp,alphax,alphaz,aval,bval,cval,ra,rb,rc)
float amp, alphax, alphaz, aval, bval, cval, *ra, *rb, *rc;
int nnx, ny, seed;
{
  float background, ratio;
  int k, l;
  float var, *y, df, c1, c2, nnx2, ny2;
  struct complex *x;
  
  background = aval;
  y = ra;
  x  = (struct complex *) malloc(sizeof(struct complex)*ny*nnx);
  /* READ 0 MEAN, 1 VARIANCE RANDOM ARRAY */
  getrand(x,2*nnx*ny,seed);
  
  /* ESTABLISH CONSTANTS */
  alphax *= PI;
  alphaz *= PI;
  alphax *= -alphax;
  alphaz *= -alphaz;
  alphax /= F 4.0;
  alphaz /= F 4.0;
  df = ( F 2.0 * PI);
  df *= df / F (nnx * ny);
  df = sqrt(df);
  nnx2 = F (nnx * nnx);
  ny2 = F (ny * ny);
  
  /* ZERO THE ZERO FREQ AND NYQUIST */
  x[0].re = F 0.0;
  x[0].im = F 0.0;
  x[nnx/2*(ny+1)].re = F 0.0;
  x[nnx/2*(ny+1)].im = F 0.0;
  for (k=0; k<ny; k++) {
    x[k*nnx+nnx/2].re = F 0.0;
    x[k*nnx+nnx/2].im = F 0.0;
  }
  for (l=0; l<nnx; l++) {
    x[ny*nnx/2+l].re = F 0.0;
    x[ny*nnx/2+l].im = F 0.0;
  }
  
  /* ADJUST TOP ROW */
  for (l=1; l<nnx/2; l++) {
    c1 = df * exp(alphax * F (l*l) / nnx2);
    x[l].re *= c1;
    x[l].im *= c1;
    x[nnx-l].re *= c1;
    x[nnx-l].im *= c1;
  }
  
  /* ADJUST LEFT COLUMN */
  for (k=1; k<ny/2; k++) {
    c1 = df * exp(alphaz * F (k*k) / ny2);
    x[k*nnx].re *= c1;
    x[k*nnx].im *= c1;
    x[(ny-k)*nnx].re *= c1;
    x[(ny-k)*nnx].im *= c1;
  }
  
  /* ADJUST THE FOUR QUANDRANTS */
  for (k=1; k<ny/2; k++) {
    c2 = F (k * k) / ny2;
    for (l=1; l<nnx/2; l++) {
      c1 = df * exp(alphax * F (l*l) / nnx2 + alphaz * c2);
      x[k*nnx+l].re *= c1;
      x[k*nnx+l].im *= c1;
      x[k*nnx+nnx-l].re *= c1;
      x[k*nnx+nnx-l].im *= c1;
      x[(ny-k)*nnx+l].re *= c1;
      x[(ny-k)*nnx+l].im *= c1;
      x[(ny-k)*nnx+nnx-l].re *= c1;
      x[(ny-k)*nnx+nnx-l].im *= c1;
    }
  }
  
  /* FFT */
  buf= (struct complex *)malloc(ny*sizeof(struct complex));
  if (buf == NULL) {
    fprintf(stderr,"cannot allocate memory for cfft2d\n");
    exit(-1);
  }
  cfft2d(x,nnx,ny);
  
  /* MAKE VARIANCE OF OUTPUT amp */
  var = F 0.0;
  for (l=0; l<nnx*ny; l++) {
    y[l] = x[l].re;
    var += y[l] * y[l];
  }
  var = sqrt(var / (nnx * ny)) / amp;
  
  /* make densities (at 30% of variance) */
  ratio = 0.3 * cval / aval;
  ratio /= var;
  for (l=0; l<nnx*ny; l++)
    rc[l] = y[l] * ratio + cval;
  
  /* make p-velocities */
  for (l=0; l<nnx*ny; l++)
    y[l] = y[l] / var + background;
  
  /* make s-velocities */
  ratio = bval / aval;
  for (l=0; l<nnx*ny; l++)
    rb[l] = ra[l] * ratio;
  free(x);
  free(buf);
}

static unsigned long state1[64] = { 0x9e1974e3, 0x1139f111, 0x9697cf7c, 0xae8f8c56,
                           0xb68981ec, 0xced50563, 0x51e0fde0, 0x7919536f,
                           0x29bf2258, 0x32e0de56, 0x73b34052, 0x0cd46c0e,
                           0x6372ce11, 0xb6acd734, 0x05be5c24, 0x0688c427,
                           0x35c51da4, 0xc1e1aec2, 0x0cd79145, 0xc7753992,
                           0xf5c44de4, 0xb489cf57, 0x396a5272, 0x74be9f2f,
                           0x30fd01fe, 0x19209141, 0x3f4257eb, 0x379c4eec,
                           0x29cb87e5, 0x7ec6962e, 0xc42c2f5e, 0x6cc9c4cf,
                           0x0b8f7ebc, 0xebfad6c7, 0x5a9db56d, 0xfdecda8e,
                           0x2ac9c5b2, 0x54c33fb5, 0x23b14a24, 0x7c38fe5c,
                           0xbb2d7111, 0x162541ec, 0x6586b906, 0xf108d75d,
                           0xfbdfce8f, 0x34894fff, 0x03b68ce9, 0x243bb116,
                           0x194bb8d3, 0xe4df3cb5, 0xb0f38722, 0x9aec9a03,
                           0xe41e6afb, 0xd0a85280, 0x6665b88c, 0x0e7a62c6,
                           0x634d78a0, 0x25d231b0, 0x358e023e, 0x3a52526f,
                           0xd9b8f072, 0xae393f93, 0xa90342f3, 0x12fa295b };

#define LEN 256

void
getrand(x,n,seed)
float *x;
int n, seed;
{
  int i, j;
  float sum, var, ave;
  
  /* SUM OVER 12 SETS OF RANDOM NUMBERS TO GET NORMAL DISTRIBUTION */
  /* initstate((unsigned)1,state1,LEN);
     setstate(state1); */
  srand(seed);
  for (i=0; i<n; i++)
    x[i] = F rand();
  for (j=0; j<11; j++)
    for (i=0; i<n; i++)
      x[i] += F rand();
  /* REMOVE MEAN */
  sum = F 0.0;
  for (i=0; i<n; i++)
    sum += x[i];
  ave = sum / F n;
  /* NORMALIZE VARIANCE TO 1.0 */
  var = F 0.0;
  for (i=0; i<n; i++) {
    x[i] -= ave;
    var += x[i] * x[i];
  }
  var = sqrt(var / F n);
  for (i=0; i<n; i++)
    x[i] /= var;
}

int
pow2(n)
register int n;
{
  register int k=1;
  
  while (k<n) k *= 2;
  return(k);
}

void
cfft2d(x,nnx,ny)
struct complex *x;
int nnx, ny;
{
  register struct complex *px;
  register int ix, iy;
  
  px= x;
  for (iy=0; iy<ny; iy++) {
    cfft(px,nnx,1);
    px += nnx;
  }
  
  for (ix=0; ix<nnx; ix++) {
    px = x + ix;
    for (iy=0; iy<ny; iy++) {
      buf[iy].re = px->re;
      buf[iy].im = px->im;
      px += nnx;
    }
    cfft(buf,ny,1);
    px = x + ix;
    for (iy=0; iy<ny; iy++) {
      px->re = buf[iy].re;
      px->im = buf[iy].im;
      px += nnx;
    }
  }
}

static float sin_table[] =
{
  1.0000000e+00,	/* sin(pi/2) */
  7.0710678e-01,	/* sin(pi/4) */
  3.8268343e-01,	/* sin(pi/8) */
  1.9509032e-01,	/* sin(pi/16) */
  9.8017140e-02,	/* sin(pi/32) */
  4.9067674e-02,	/* sin(pi/64) */
  2.4541228e-02,	/* sin(pi/128) */
  1.2271538e-02,	/* sin(pi/256) */
  6.1358846e-03,	/* sin(pi/512) */
  3.0679568e-03,	/* sin(pi/1024) */
  1.5339802e-03,	/* sin(pi/2048) */
  7.6699032e-04,	/* sin(pi/4096) */
  3.8349519e-04,	/* sin(pi/8192) */
  1.9174760e-04,	/* sin(pi/16384) */
  9.5873799e-05,        /* sin(pi/32768) */
  4.7936899e-05,        /* sin(pi/65536) */
  2.3968499e-05,        /* sin(pi/131072) */
  1.1984224e-05,        /* sin(pi/262144) */
  5.9921125e-06,        /* sin(pi/524288) */
  2.9960562e-06,        /* sin(pi/1048576) */
  1.4980281e-06         /* sin(pi/2097152) */
};

/*
 * cfft - radix 2 FFT for complex data
 *
 * n is the number of complex points.
 * x is the single precision (not double!!) complex vector.
 * isign is the sign of the transform.
 *
 * the routine does no normalization
 */
void
cfft(x,n,isign)
struct complex *x;
int n,isign;
{
  register struct complex *px, *qx, *rx;
  struct complex *limit, *qlimit, dtemp;
  float cn, sn, cd, sd, temp, real, imag;
  int m, j, istep;
  float *psintab;
  extern float sin_table[];
  
  limit= x + n;
  j= 0;
  for(px=x; px<limit; px++)
    {
      if(px < (qx= x+j))
	{	dtemp= *qx; *qx= *px; *px= dtemp;	}
      m = n>>1;
      while( m>=1 && j>=m )
	{ j-= m; m>>= 1;    }
      j+= m;
    }
  rx= x+1;
  for(px=x; px<limit; px+= 2, rx+= 2)
    {
      temp= rx->re;
      rx->re= px->re -temp;
      px->re += temp;
      temp= rx->im;
      rx->im= px->im -temp;
      px->im += temp;
    }
  j=2;
  psintab= sin_table;
  while( j < n )
    {
      istep= j<<1;
      sd= *psintab++;
      temp= *psintab;
      cd= F 2.0 * temp * temp;
      cn= F 1.0;
      sn= F 0.0;
      if( isign < 0 ) sd= -sd;
      qlimit= x+j;
      for(qx=x; qx< qlimit; qx++)
	{
	  for(px=qx; px<limit; px+= istep)
	    {
	      rx= px + j;
	      real= cn * rx->re - sn * rx->im;
	      imag= sn * rx->re + cn * rx->im;
	      rx->re = px->re - real;
	      rx->im = px->im - imag;
	      px->re += real;
	      px->im += imag;
	    }
	  temp= cd * cn + sd * sn;
	  sn += (sd * cn - cd * sn);
	  cn -= temp;
	}
      j= istep;
    }
}
