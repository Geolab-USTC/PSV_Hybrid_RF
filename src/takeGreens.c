#include <stdio.h>

main(ac,av)
int ac; char **av;
{

    int m34, n34, ix, iz;
    int nt; int nx, nz, nl, nf, nd;
    long int rec_len;
    int i;

    char greenfile[40], outfile[40];
    float green[10000];
    FILE *infile, *open_file();
    void read_record();

    setpar(ac,av);
    getpar("ix",    "d", &ix);
    getpar("nx",    "d", &nx);
    getpar("nl",    "d", &nl);
    getpar("nf",    "d", &nf);
    getpar("nd",    "d", &nd);
    getpar("iz",    "d", &iz);
    getpar("m34",    "d", &m34);
    getpar("n34",    "d", &n34);
    getpar("nt",    "d", &nt);
    getpar("greenfile", "s",greenfile);

    infile=open_file(greenfile,"r");

    nx = nx+nl;

    if(m34)
        rec_len = (long) (ix*nt*sizeof(float));

    if(n34)
        rec_len = (long) ((4*nx+iz)*nt*sizeof(float));
	 
    read_record(infile,rec_len,green,nt,sizeof(float),(long)(sizeof(float)));

    for(i=0; i<nt; i++)
        fprintf(stdout,"%e\n",green[i]);
    fclose(stdout);
}






