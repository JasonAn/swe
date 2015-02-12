#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "shallow_water.h"


void jacobiandrifter(struct drifter *ptdrifter, struct drifter *ptmsrdrifter, size_t ndr, double *fields_dot, double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, int **n, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order, const long int ncycle, int write_out_interval){


    /*
    *   time tool
    */

    clock_t start, finish;
    double duration;


    int i, j,  tdim;
    tdim = xdim * ydim;

    /*
     * memory of fields
     */

    double *u, *v, *P;

    u = fields;
    v = fields + tdim;
    P = fields + 2 * tdim;

    /*
     * coordinate of neighbors
     */

    int **neighbors;

    neighbors = calloc(tdim, sizeof(int*));
    for (i = 0; i< tdim; i++){
        neighbors[i] = calloc(4,sizeof(int));
    }

    initialize_geometry(xdim, ydim, neighbors, lat, lon, print_out_order);

    /*
    * memory for saving the fields
    */

    double *fields_temp;
    fields_temp = calloc(3 * tdim, sizeof(double));
    memcpy(fields_temp, fields, 3 * tdim * sizeof(double));
    const double * fields_save = fields_temp;
    free(fields_temp);

    double *fields_msr;
    fields_msr = calloc(3 * tdim, sizeof(double));

    for (i = 0; i < Dm; i++){
        read_field(fields_msr, "u", xdim, ydim, ncycle - 1, print_out_order);
        read_field(fields_msr + tdim, "v", xdim, ydim, ncycle - 1, print_out_order);
        read_field(fields_msr + 2 * tdim, "P", xdim, ydim, ncycle - 1, print_out_order);

    }


    double ***delaytensor_model, ***delaytensor_data;

    delaytensor_model = calloc(ndr, sizeof(double**));
    delaytensor_data = calloc(ndr, sizeof(double**));


    for(i = 0; i < ndr; i++ ){
        delaytensor_model[i] = calloc(Dm, sizeof(double*));
        delaytensor_data[i] = calloc(Dm, sizeof(double*));
        for(j = 0; j < Dm; j++){
            delaytensor_model[i][j] = calloc(2, sizeof(double));
            delaytensor_data[i][j] = calloc(2, sizeof(double));
        }
    }


    driftdelay(ptdrifter, delaytensor_model, ndr, fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, lat, lon, neighbors, dt, print_out_order,ncycle);
    driftdelay(ptmsrdrifter, delaytensor_data, ndr, fields_dot, fields_msr, parameters, forcing, xdim, ydim, dx, dy, lat, lon, neighbors, dt,print_out_order,ncycle);

    double *drifterdiff_0, *drifterdiff_1;
    drifterdiff_0 = calloc(Dm * ndr, sizeof(double));
    drifterdiff_1 = calloc(Dm * ndr, sizeof(double));

    for(i = 0; i < Dm; i++){
        for (j = 0; j < ndr; j++) {
            drifterdiff_0[i * ndr + j] = delaytensor_data[j][i][0] - delaytensor_model[j][i][0];
            drifterdiff_1[i * ndr + j] = delaytensor_data[j][i][1] - delaytensor_model[j][i][1];
        }
    }


    int ddim; // ddim: 0 - Dm
    int ele = 0; // ele: 0 - Dtot

    double ****jacobian;
    jacobian= calloc(3 * tdim, sizeof(double ***));

    for(ele = 0; ele < 3 * tdim; ele++){
        jacobian[ele] = calloc(ndr, sizeof(double **));
        for (i = 0; i < ndr; i++){
            jacobian[ele][i] = calloc(Dm, sizeof(double *));
            for(ddim = 0; ddim < Dm; ddim++){
                jacobian[ele][i][ddim] = calloc(2, sizeof(double));
            }
        }
    }

    double *fieldspl;
    fieldspl = calloc(3 * tdim, sizeof(double));

    double *fieldsmi;
    fieldsmi = calloc(3 * tdim, sizeof(double));

    // perturbation
    double *psi;
    psi = calloc(3 * tdim, sizeof(double));


    double ***delaytensorpl, ***delaytensormi;

    delaytensorpl = calloc(ndr, sizeof(double**));
    delaytensormi = calloc(ndr, sizeof(double**));


    for(i = 0; i < ndr; i++ ){
        delaytensorpl[i] = calloc(Dm, sizeof(double*));
        delaytensormi[i] = calloc(Dm, sizeof(double*));
        for(j = 0; j < Dm; j++){
            delaytensorpl[i][j] = calloc(2, sizeof(double));
            delaytensormi[i][j] = calloc(2, sizeof(double));
        }
    }

    start = clock();

    for (ele = 0; ele < 3 * tdim; ele++){
        printf("ele = %i, ncycle = %i\n", ele, ncycle);
        
        memcpy(fieldspl, fields_save, 3 * tdim * sizeof(double));
        memcpy(fieldsmi, fields_save, 3 * tdim * sizeof(double));

        fieldspl[ele] = fieldspl[ele] + eps * fields_save[ele];
        fieldsmi[ele] = fieldsmi[ele] - eps * fields_save[ele];

        driftdelay(ptdrifter, delaytensorpl, ndr, fields_dot, fieldspl, parameters, forcing, xdim, ydim, dx, dy, lat, lon, neighbors, dt, print_out_order,ncycle);
        driftdelay(ptdrifter, delaytensormi, ndr, fields_dot, fieldsmi, parameters, forcing, xdim, ydim, dx, dy, lat, lon, neighbors, dt, print_out_order,ncycle);

            for (i = 0; i < ndr; i++){
                for(ddim = 0; ddim < Dm; ddim++){
                    jacobian[ele][i][ddim][0] = (delaytensorpl[i][ddim][0] - delaytensormi[i][ddim][0])/ (2 * eps * fields_save[ele]);
                    jacobian[ele][i][ddim][1] = (delaytensorpl[i][ddim][1] - delaytensormi[i][ddim][1])/ (2 * eps * fields_save[ele]);
                }
            }

    }

    finish = clock();

//    print_jacobian(jac, "jac", ncycle, xdim, ydim, svd_m, svd_n);
//    fflush(stdout);

    double *jacT_0, *jacT_1;

    int svd_m = Dm * ndr, svd_n = 3 * tdim;

    jacT_0 = calloc(svd_m * svd_n, sizeof(double));
    jacT_1 = calloc(svd_m * svd_n, sizeof(double));


    for(ele = 0; ele < tdim;  ele++){
        for (ddim = 0; ddim < Dm;  ddim++){
            for (i = 0; i < ndr; i ++){
                jacT_0[ele * svd_m + ddim * ndr + i] = jacobian[ele][i][ddim][0];
                jacT_1[ele * svd_m + ddim * ndr + i] = jacobian[ele][i][ddim][1];
            }
        }
    }
	

//    print_jacobian(jacT, "jacT", ncycle, xdim, ydim, svd_n, svd_m);
// svd

    double *nudge_0 = calloc(3 * tdim, sizeof(double));
    double *nudge_1 = calloc(3 * tdim, sizeof(double));

    coupling(svd_m, svd_n, jacT_0, drifterdiff_0, nudge_0);
    coupling(svd_m, svd_n, jacT_1, drifterdiff_1, nudge_1);

//next step

    memcpy(fields, fields_save, 3 * tdim * sizeof(double));

//    for (i = 0; i < 3 * tdim; i++){
//        fields[i] = fields_save[i];
//    }

//    double *fields_next;
//
//    fields_next = calloc(tdim, sizeof(double));

//    read_field(fields_next, xdim, ydim, ncycle, print_out_order);

    RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);

    for (i = 0; i < 2 * tdim; i++) {
        fields[i] += 0.5 * (nudge_0[i] + nudge_1[i]);
    }
    for (i = 2 * tdim; i < 3 * tdim; i++){
        fields[i] += 1.5 * (nudge_0[i] + nudge_1[i]);
    }

    print_field(u, "u", ncycle, xdim, ydim, print_out_order);
    print_field(v, "v", ncycle, xdim, ydim, print_out_order);
    print_field(P, "P", ncycle, xdim, ydim, print_out_order);
    fflush(stdout);

    duration = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("time = %f \n", duration);

    for (i = 0; i< tdim; i++){
        free(neighbors[i]);
    }

    free(neighbors);

    free(fields_save);
//    free(fields_ori);
    free(fields_msr);

    free(fieldspl);
    free(fieldsmi);

    free(psi);

    free(jacT_0);
    free(jacT_1);

    free(drifterdiff_0);
    free(drifterdiff_1);


}


void driftdelay(struct drifter *ptdrifter, double ***delaytensor, size_t ndr,  double *fields_dot, double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, unsigned int *lat, unsigned *lon, int **neighbors, double dt, unsigned int *print_out_order, long int ncycle){

    int i, j, k;
    

    for(i = 0; i < ndr; i++){
        delaytensor[i][0][0] = ptdrifter[i].x;
        delaytensor[i][0][1] = ptdrifter[i].y;
    }

    for(i = 1; i < Dm; i++){
        for (j = 0; j < tau; j++){
            RK4(fields_dot, fields, parameters, forcing, xdim, ydim, dx, dy, neighbors, lat, lon, print_out_order, ncycle);
            drift(ptdrifter, ndr, fields, xdim, ydim, dx, dy, lat, lon, dt);
        }

        for(k = 0; k < ndr; k++) {
            delaytensor[k][i][0] = ptdrifter[i].x;
            delaytensor[k][i][1] = ptdrifter[i].y;
        }

    }

}


void coupling(int svd_m, int svd_n, double *jacT, double *fields_diff, double *temp2) {

    int info, lwork;
    double wkopt;
    double *work;

    double *svd_s, *svd_u, *svd_vt;

    svd_s = calloc(svd_n, sizeof(double));
    svd_u = calloc(svd_m * svd_m, sizeof(double));
    svd_vt = calloc(svd_n * svd_n, sizeof(double));


    int lda = svd_m, ldu = svd_m, ldvt = svd_n;

    lwork = -1;
    dgesvd_("All", "All", &svd_m, &svd_n, jacT, &lda, svd_s, svd_u, &ldu, svd_vt, &ldvt, &wkopt, &lwork, &info);

    lwork = (int) wkopt;

    work = (double *) malloc(lwork * sizeof(double));

/* Compute SVD */
//        printf("part two");
    dgesvd_("All", "All", &svd_m, &svd_n, jacT, &lda, svd_s, svd_u, &ldu, svd_vt, &ldvt, work, &lwork, &info);

/* Check for convergence */
    if (info > 0) {
        printf("The algorithm computing SVD failed to converge.\n");
        exit(1);
    }

//        print_svd(jac_s, "jac_s", ncycle, xdim, ydim, 1, svd_n, 1 );
////      printf("%i, %i, %i, %i \n" , ldu, ldvt, svd_m, svd_n);
//        print_svd(jac_u, "jac_u", ncycle, xdim, ydim, svd_m, svd_n, ldu);
//        print_svd(jac_vt, "jac_v", ncycle, xdim, ydim, svd_n, svd_n, ldvt);


    free((void *) work);


//inverse

    double *svd_sinv;
    svd_sinv = calloc(svd_n, sizeof(double));

    int i;
    for (i = 0; i < svd_n; i++) {
        if (svd_s[i] > 0.1)
            svd_sinv[i] = 1.0 / svd_s[i];
    }

    double *temp;
    temp = calloc(svd_n, sizeof(double));

    int svd_k = 1;
    double alpha = 1.0, beta = 0.0;

    dgemm_("T", "N", &svd_n, &svd_k, &svd_m, &alpha, svd_u, &svd_m, fields_diff, &svd_m, &beta, temp, &svd_n);
//    print_svd(temp, "temp1", ncycle, xdim, ydim, svd_n, 1, svd_n);


//temp and temp2 are used for SVD, implicit vectors.

    for (i = 0; i < svd_n; i++) {
        temp[i] = temp[i] * svd_sinv[i];
    }

//    print_svd(temp, "temp2", ncycle, xdim, ydim, svd_n, 1, svd_n);

    dgemm_("T", "N", &svd_n, &svd_k, &svd_n, &alpha, svd_vt, &svd_n, temp, &svd_n, &beta, temp2, &svd_n);

    free(svd_s);
    free(svd_u);
    free(svd_vt);

    free(svd_sinv);

    free(temp);

}
