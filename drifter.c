#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "shallow_water.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include "interp2d.h"
#include "interp2d_spline.h"
#include "interp2d.c"
#include "interp2d_spline.c"
#include "bicubic.c"

void velinterp (struct drifter *ptdrifter, int ndr, int xdim, int ydim, double *u, double *v);

void drift(struct drifter *ptdrifter, size_t ndr, double *fields, int xdim, int ydim, double dx, double dy, unsigned int *lat, unsigned *lon, double dt){


//    clock_t start, finish;
//    double duration;

    int i, tdim;
    tdim = xdim * ydim;
	
    double *u, *v, *P;
	
    u = calloc((xdim + 1) * (ydim + 1), sizeof(double));
    v = calloc((xdim + 1) * (ydim + 1), sizeof(double));
    P = calloc((xdim + 1) * (ydim + 1), sizeof(double));

    for(i = 0; i < tdim; i++){
        u[lat[i] * xdim + lon[i] + lat[i]] = fields[i];
        v[lat[i] * xdim + lon[i] + lat[i]] = fields[i + tdim];
        //P[lat[i] * xdim + lon[i] + lat[i]] = fields[i + 2 * tdim];
    }

    for (i = 0; i < xdim; i++){
        u[xdim * ydim + xdim + i] = u[i];
        v[xdim * ydim + xdim + i] = v[i];
    }

    for (i = 0; i <= ydim; i++){
        u[i * xdim  + i + xdim] = u[i * xdim + i];
        v[i * xdim  + i + xdim] = v[i * xdim + i];
    }


    //printf("u =  %f, %f \n", u[0], u[15]);


	velinterp(ptdrifter, ndr, xdim + 1, ydim + 1, u, v);

	for (i = 0; i < ndr; i++){
		ptdrifter[i].x += ptdrifter[i].u * dt * 1.0 / dx;
		ptdrifter[i].y += ptdrifter[i].v * dt * 1.0 / dy;
		ptdrifter[i].x -= xdim * floor(ptdrifter[i].x / xdim);
        ptdrifter[i].y -= ydim * floor(ptdrifter[i].y / ydim);
    }

    free(u);
    free(v);
    free(P);

}


void velinterp (struct drifter *ptdrifter, int ndr, int xdim, int ydim, double *u, double *v){

	int i, j;
	double *xarr, *yarr;
	xarr = calloc(xdim, sizeof(double));
	yarr = calloc(ydim, sizeof(double));
	
	for (i = 0; i < xdim; i++)
		xarr[i] = i * 1.0;
	
	for (j = 0; j < ydim; j++)
		yarr[j] = j * 1.0;

	gsl_interp_accel *xa, *ya;
	
	const interp2d_type* T = interp2d_bicubic;
	
	xa = gsl_interp_accel_alloc();
    ya = gsl_interp_accel_alloc();

    interp2d_spline* interp_u = interp2d_spline_alloc(T, xdim, ydim);
	interp2d_spline* interp_v = interp2d_spline_alloc(T, xdim, ydim);
	
	interp2d_spline_init(interp_u, xarr, yarr, u, xdim, ydim);
	interp2d_spline_init(interp_v, xarr, yarr, v, xdim, ydim);

	double xdr, ydr;

	for (i = 0; i < ndr; i++) {
		xdr = ptdrifter[i].x;
		ydr = ptdrifter[i].y;
		ptdrifter[i].u = interp2d_spline_eval(interp_u, xdr, ydr, xa, ya);
		ptdrifter[i].v = interp2d_spline_eval(interp_v, xdr, ydr, xa, ya);
	}

	free(xarr);
	free(yarr);
	gsl_interp_accel_free(xa);
    gsl_interp_accel_free(ya);
    interp2d_spline_free(interp_u);
	interp2d_spline_free(interp_v);
	
}
