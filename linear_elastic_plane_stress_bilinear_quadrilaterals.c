/*
 * Linear elastisc FEM (2-D plane stress bilinear quadrilaterials)
 * Disclaimer: Testing stage + re-learnig C at the moment.
 * gcc -Wall linear_elastic_plane_stress_bilinear_quadrilaterals.c -llapacke -llapack -lblas -lcblas -lm 
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cblas.h>
#include <lapacke.h>

int print_array_col_major(char *msg, double *array, int nrow, int ncol);
int create_elasticity_matrix_isotropic_plane_stress(double *C, double *pe_modulus, double *pnu);
int create_element_stiffness_matrix(double *Ke, int *num_gp, double *X);

int main (void) {
    /* Array Format Explanation (uncomment '//' to show examples)
     * ======================== */

    /* Example illustrating the array format with a 5 row 3 column matrix */
    // double mat[5*3] = {1,2,3,4,5,
    //                    1,3,5,2,4,
    //                    1,4,2,5,3};
    // int nr = 5;  /* num rows */
    // mat[1 + 2*nr] = 4.5;  /* access mat[2,3] so mat[1][2] */
    /* try and ignore the '+' and the '*nr'; so convert mat[1 + 2*nr] to mat[1][2] in your head.*/
    // print_array_col_major("mat", mat, 5, 3);
   
    // double rowvec[3] = {1, 2, 3};  /* [1x3] */
    // print_array_col_major("rowvec:", rowvec, 1, 3);

    // double colvec[3] = {1, 2, 3};  /* [3x1] */
    // print_array_col_major("colvec:", colvec, 3, 1);

    /* Actual Program
     * ============== */

    double length_x = 2;
    double length_y = 2;
    double thickness = 0.1;

    double e_modulus = 10000;
    double nu = 0.3;

    int num_ele_x = 2;
    int num_ele_y = 2;

    int num_gp = 1;


    int num_ele = num_ele_x * num_ele_y;
    int num_dof = 2*((num_ele_x + 1)*(num_ele_y + 1));

    double length_ele_x = length_x/num_ele_x;
    double length_ele_y = length_y/num_ele_y;

    printf("Number of Elements: %d\n", num_ele);
    printf("Number of DOFs: %d\n", num_dof);

    double C[3*3];
    create_elasticity_matrix_isotropic_plane_stress(C, &e_modulus, &nu);
    print_array_col_major("Elasticity Matrix, C:", C, 3, 3);

    double a[2*2] = {1,1,2,2};
    double b[2*2] = {2, 3, 4, 6};
    double c[2*2];
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    2, 2, 2, 1.0, a, 2, b, 2, 1.0, c, 2);
    print_array_col_major("res c:", c, 2, 2);

    double Ke[8*8];
    double X[4*2];
    create_element_stiffness_matrix(Ke, &num_gp, X);
    

    exit(0);
}

/* Functions
 * --------- */

int print_array_col_major(char *msg, double *array, int nrow, int ncol) {
   /* Print column-major arrays of one or two dimensions as 'normal' row-major arrays. */
   /* TODO: Make output prettier */
    printf("%s\n", msg);
    printf("[%d x %d]\n", nrow, ncol);

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            printf("%10E \t", array[i + j*nrow]);
        }
        printf("\n");
    }
    printf("\n");

    return 0;
}

int create_elasticity_matrix_isotropic_plane_stress(double *C, double *pe_modulus, double *pnu) {
    int nr = 3;
    int ne = 9;
    int mul_factor = *pe_modulus/(1.0 - pow(*pnu, 2.0));

    for (int i = 0; i < ne; i++) {
        C[i] = 0.0;
    }

    C[0 + 0*nr] = 1;  /* C[0][0] = 1 */
    C[0 + 1*nr] = *pnu;
    C[1 + 0*nr] = *pnu;
    C[1 + 1*nr] = 1;
    C[2 + 2*nr] = (1.0 - *pnu)/2.0;

    for (int i = 0; i < ne; i++) {
        C[i] *= mul_factor;
    }
    
    return 0;
}

int create_element_stiffness_matrix(double *Ke, int *num_gp, double *X) {
    /* TODO: Extend to allow for 4 Gauss Points */
    lapack_int info_dgetrf_Jinv = 0;
    lapack_int info_dgetri_Jinv = 0;
    
    double xi = 0;
    double eta = 0;
    double weight = 4;

    /* Construct P, containing partial derivatives of the shape functions */
    double P[2*4];  /* P is a [2x4] matrix */
    int nrp = 2;
    P[0 + 0*nrp] = -(1.0/4.0)*(1 - eta);  /* P[0][0] = dN1dxi */
    P[0 + 1*nrp] =  (1.0/4.0)*(1 - eta);  /* P[0][1] = dN2dxi */
    P[0 + 2*nrp] =  (1.0/4.0)*(1 + eta);
    P[0 + 3*nrp] = -(1.0/4.0)*(1 + eta);
    P[1 + 0*nrp] = -(1.0/4.0)*(1 - xi );  /* P[1][0] = dN1deta*/
    P[1 + 1*nrp] = -(1.0/4.0)*(1 + xi );  /* P[1][1] */
    P[1 + 2*nrp] =  (1.0/4.0)*(1 + xi );
    P[1 + 3*nrp] =  (1.0/4.0)*(1 - xi );
    print_array_col_major("Partial derivatives of the shape functions, P:",
    P, 2, 4);

    /* Uncomment '//' for testing */
    int nrx = 4;
    X[0 + 0*nrx] = 0;  /* X[0][0] = 0.0 */
    X[0 + 1*nrx] = 0;  /* X[0][1] = 0.0 */
    X[1 + 0*nrx] = 0.8;
    X[1 + 1*nrx] = 0;
    X[2 + 0*nrx] = 1;
    X[2 + 1*nrx] = 1;
    X[3 + 0*nrx] = 0;
    X[3 + 1*nrx] = 1;
    print_array_col_major("Nodal Coordinates, X:", X, 4, 2);

    /* Calculate the jacobian matrix, J, via J = P*X */
    double J[2*2];
    for (int i = 0; i < 4; i++) {
        J[i] = 0.0;
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    2, 2, 4, 1.0, P, 2, X, 4, 1.0, J, 2);
    print_array_col_major("Jacobian Matrix, J:", J, 2, 2);

    /* Calculate the inverse of the jacobian matrix, J^-1 */
    int ipiv[2];
    double Jinv[2*2];
    memcpy(Jinv, J, sizeof(J));
    print_array_col_major("Inverse of the Jacobian Matrix, J^-1:",
    Jinv, 2, 2);
    info_dgetrf_Jinv = LAPACKE_dgetrf(LAPACK_COL_MAJOR, 2, 2, Jinv, 2, ipiv);
    /* printf("info_dgetrf_Jinv: %d\n", info_dgetrf_Jinv); */
    info_dgetri_Jinv = LAPACKE_dgetri(LAPACK_COL_MAJOR, 2, Jinv, 2, ipiv);
    /* printf("info_dgetri_Jinv: %d\n", info_dgetri_Jinv); */ 
    print_array_col_major("Inverse of the Jacobian Matrix, J^-1:",
    Jinv, 2, 2);

    /* Calculate the strain-displacement matrix, B */
    /* First we calculate dNidx and dNidy for each shape function,
     * then assenble B from dNidx and dNidy starting with an empty B matrix. */
    


    return 0;
}
