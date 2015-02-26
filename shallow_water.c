#include <stdlib.h>
#include <stdio.h>

#include "shallow_water.h"


int main (int argc, char ** argv)
{
    int i, j;

    /*
    * the geometry of the water layer
    * tdim is the total number of grids, i.e. tdim = xdim * ydim
    */

    int xdim, ydim, tdim;
    double dx, dy;


    /*
    * coordinate of neighbors
    */

    int **neighbors;
    /*
    * latitude and longitude values of each point
    */

    unsigned int *lat, *lon;

    /*
    * order of data (read in and print out the layer)
    */

    unsigned int *print_out_order;

    /*
    * temporal aspects of the simulation
    */

    long int ncycle;
    double dt = 5.1;
    long int n_iter = 25800;
    int write_out_interval = 5000;


    /*
    * physical aspects
    *
    * f0:           coriolis parameter (1/sec)
    * beta:         linear beta for coriolis force (1/meter sec)
    * forcingTerm:  wind stress amplitude, "tau_0" (N/meter**2)
    * dissipationTerm:  viscosity coefficient (meters**2/dec)
    * RayleighFriction: Rayleigh friction parameter (1/sec)
    */

    double f0 = 5.0e-5;
    double beta = 2e-11;
    double forcingTerm = 0.0005;
    double dissipationTerm = 0.00005;
    double RayleighFriction = 5e-8;

    /*
    * upper layer equilibrium depth (meter)
    */

    double EquilibriumDepth = 50000.0;
    double A;

    /*
    * physcal parameters and the forcing term
    */

    double *parameters;
    double *forcing;
    double *x_forcing, *y_forcing, *z_forcing;

    /*
    * fields holds the values of the u, v and P
    */
    double *fields;
    double *u, *v, *P;

    /*
    * fields_msr holds the measured data
    */

//    double *fields_msr;
//    double *u_msr, *v_msr, *P_msr;

    /*
    * read parameters from the command line
    */

    if (argc == 2){

    get_parameters(argv[1], &xdim, &ydim, &dx, &dy, &n_iter, &dt, &write_out_interval, &f0, &beta, &forcingTerm, &dissipationTerm, &RayleighFriction, &EquilibriumDepth, &A);

    }

    else{

    fprintf(stderr, "Error: You should give a file containing the parameters as an argument! \n");
    exit(1);

    }

    print_parameters(xdim, ydim, dx, dy, n_iter, dt, write_out_interval, f0, beta, forcingTerm, dissipationTerm, RayleighFriction, EquilibriumDepth, A);

    /*
    * generate some physical parameters
    */

    tdim = xdim * ydim;

    /*
    * allocate arrays containing geometrical information
    */

    lat = calloc(tdim, sizeof(unsigned int));
    lon = calloc(tdim, sizeof(unsigned int));
    print_out_order = calloc(tdim, sizeof(unsigned int));

    neighbors = calloc(tdim, sizeof(int*));
    for (i = 0; i< tdim; i++){
        neighbors[i] = calloc(4,sizeof(int));
    }

    initialize_geometry(xdim, ydim, neighbors, lat, lon, print_out_order);


    /*
    * memory for physical parameters and forcing
    */
    parameters = calloc(5, sizeof(double));

    parameters[0] = f0;
    parameters[1] = beta;
    parameters[2] = forcingTerm;
    parameters[3] = dissipationTerm;
    parameters[4] = RayleighFriction;


    forcing = calloc(3 * tdim, sizeof(double));

    x_forcing = forcing;
    y_forcing = forcing + tdim;
    z_forcing = forcing + 2 * tdim;

    /*
    * memory of fields measured
    */

    fields = calloc(3 * tdim, sizeof(double));

    u = fields;
    v = fields + tdim;
    P = fields + 2 * tdim;

    /*
    * memory of fields measured
    */

//    fields_msr = calloc(3 * tdim, sizeof(double));
//
//    u_msr = fields_msr;
//    v_msr = fields_msr + tdim;
//    P_msr = fields_msr + 2 * tdim;

    double *fields_dot;
    fields_dot = calloc(3 * tdim, sizeof(double));


    /*
    * initialize fields
    */
    initialize_fields (u, v, P, x_forcing, y_forcing, z_forcing, xdim, ydim, dx, dy, EquilibriumDepth, A, neighbors, lat, lon);

    ncycle=0;
    
    print_field(u, "u", ncycle, xdim, ydim, print_out_order);
    print_field(v, "v", ncycle, xdim, ydim, print_out_order);
    print_field(P, "P", ncycle, xdim, ydim, print_out_order);


    int ndr = 196; // number of drifters
    int sndr = 14;

    struct drifter *model_ptdrifter = calloc(ndr,sizeof(struct drifter));
    struct drifter *data_ptdrifter = calloc(ndr,sizeof(struct drifter));

    for(i =0; i< sndr; i++) {
        for(j = 0; j < sndr; j++) {
            model_ptdrifter[i * sndr + j].x = data_ptdrifter[i * sndr + j].x = 1 + j;
            model_ptdrifter[i * sndr + j].y = data_ptdrifter[i * sndr + j].y = 1 + i;
        }
    }
//    for(i = 0; i < ndr; i++){
//        printf("%i, %f, %f\n", i, model_ptdrifter[i].x, data_ptdrifter[i].y);
//    }


    /*
    * loop through time steps to do the actual simulation
    */

    if(tsyn > 0) {
        for (ncycle = 1; ncycle < n_iter; ncycle++) {
            if (ncycle >= tstart && ncycle < tstart + tsyn) {
                //jacobian(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon, print_out_order, ncycle, write_out_interval);
                jacobiandrifter(model_ptdrifter, data_ptdrifter, ndr, fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, dt, neighbors, lat, lon, print_out_order, ncycle, write_out_interval);
                printf("ncycle!! = %li \n", ncycle);
            }
            else {
                RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);
                print_field(u, "u", ncycle, xdim, ydim, print_out_order);
                print_field(v, "v", ncycle, xdim, ydim, print_out_order);
                print_field(P, "P", ncycle, xdim, ydim, print_out_order);
                printf("ncycle = %li \n", ncycle);
            }
        }
    }
    else if(tsyn == 0)
    {
        for (ncycle = 1; ncycle < n_iter; ncycle++) {
                RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);
                print_field(u, "u", ncycle, xdim, ydim, print_out_order);
                print_field(v, "v", ncycle, xdim, ydim, print_out_order);
                print_field(P, "P", ncycle, xdim, ydim, print_out_order);
                printf("ncycle = %li \n", ncycle);
            }
        }


	printf("\n");


    free(model_ptdrifter);
    free(data_ptdrifter);

    return (0);
}


