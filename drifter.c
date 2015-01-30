#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include <interp2d.h>
#include "/media/jason/Program/Library/interp2d-master/interp2d.h"
#include "/media/jason/Program/Library/interp2d-master/interp2d_spline.h"
#include "/media/jason/Program/Library/interp2d-master/interp2d.c"
#include "/media/jason/Program/Library/interp2d-master/interp2d_spline.c"
#include "/media/jason/Program/Library/interp2d-master/bicubic.c"

struct drifter{
    double x, y;
    double u, v;
};


void velinterp (double xdr, double ydr, int xdim, int ydim, double *u, double *v, double *ptudr, double *ptvdr);

void drift(double *fields, int xdim, int ydim, unsigned int *lat, unsigned *lon){

    int n = 2; // number of the drifters

//    clock_t start, finish;
//    double duration;

    int i, j, tdim;
    tdim = xdim * ydim;
	
    double *u, *v, *P;
	
    u = calloc(tdim, sizeof(double));
    v = calloc(tdim, sizeof(double));
    P = calloc(tdim, sizeof(double));

    for(i = 0; i < tdim; i++){
        u[lat[i] * ydim + lon[i]] = fields[i];
		v[lat[i] * ydim + lon[i]] = fields[i + tdim];
        P[lat[i] * ydim + lon[i]] = fields[i + 2 * tdim];
		
    }

	struct drifter *ptdrifter = calloc(n,sizeof(struct drifter));
	
	ptdrifter[0].x = 5.0;
	ptdrifter[0].y = 7.0;
	
	int mi = 0;
	do{
		mi++;
	}while(lat[mi] == ptdrifter[0].x && lon[mi] == ptdrifter[0].y);
	
	ptdrifter[1].x = 10.0;
	ptdrifter[1].y = 3.0;
	
	
	// need initialize the positions of the n drifters
	
	for (i = 0; i < n; i++){
		double xdr, ydr;  // temp memory for drifter position 
		double udr, vdr;
		double *ptudr = &udr;
		double *ptvdr = &vdr;
		xdr = ptdrifter[i].x;
		ydr = ptdrifter[i].y;
		velinterp (xdr, ydr, xdim, ydim, u, v, ptudr, ptvdr);
		ptdrifter[i].u = udr;
		ptdrifter[i].v = vdr;
	}

	for (i = 0; i < n; i++){
		printf("%i drifter, u = %f, v = %f\n", i, ptdrifter[i].u, ptdrifter[i].v);
		printf("%i drifter, u = %f, v = %f\n", i, u[mi], v[mi]);
	}
	
}


void velinterp (double xdr, double ydr, int xdim, int ydim, double *u, double *v, double *ptudr, double *ptvdr){
	
	int i, j;
	double *xarr, *yarr;
	xarr = calloc(xdim, sizeof(double));
	yarr = calloc(ydim, sizeof(double));
	
	for (i = 0; i < xdim; i++)
		xarr[i] = i * 1.0;
	
	for (j = 0; j < ydim; j++)
		yarr[j]= j * 1.0;

	gsl_interp_accel *xa, *ya;
	
	const interp2d_type* T = interp2d_bicubic;
	
	xa = gsl_interp_accel_alloc();
    ya = gsl_interp_accel_alloc();

    interp2d_spline* interp_u = interp2d_spline_alloc(T, xdim, ydim);
	interp2d_spline* interp_v = interp2d_spline_alloc(T, xdim, ydim);
	
	interp2d_spline_init(interp_u, xarr, yarr, u, xdim, ydim);
	interp2d_spline_init(interp_v, xarr, yarr, v, xdim, ydim);
	
	*ptudr = interp2d_spline_eval(interp_u, xdr, ydr, xa, ya);
	*ptvdr = interp2d_spline_eval(interp_v, xdr, ydr, xa, ya);
	
	free(xarr);
	free(yarr);
	gsl_interp_accel_free(xa);
    gsl_interp_accel_free(ya);
    interp2d_spline_free(interp_u);
	interp2d_spline_free(interp_v);
	
}