/*
 * Linear elastisc FEM (2-D plane stress bilinear quadrilaterials)
 * Disclaimer: Testing stage + re-learnig C at the moment.
 * gcc -Wall linear_elastic_plane_stress_bilinear_quadrilaterals.c -llapacke -llapack -lblas -lcblas -lm 
 */

#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cblas.h>
#include <lapacke.h>

int print_array_col_major(char *msg, double *array, int nrow, int ncol);
int print_array_col_major_int(char *msg, int *array, int nrow, int ncol);
int create_elasticity_matrix_isotropic_plane_stress(double *C, double *pe_modulus, double *pnu);
int create_element_stiffness_matrix(double *Ke, double *C, int *num_gp, double *X, double *thickness);
int create_global_stiffness_matrix(double *K, int* nrk, int *num_ele, int *num_ele_x, double *len_ele_x,
double *len_ele_y, int *num_gp, double *thickness, double *C);
int apply_bcs(double *Kmod, int *bc, double *F, int *num_ele_x, int *num_ele_y, int *num_dof, double *len_ele_x);
int solve_for_displscements(double *Kmod, double *U, int *num_dof);

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

    double len_x = 2;
    double len_y = 2;
    double thickness = 0.1;

    double e_modulus = 10000;
    double nu = 0.3;

    int num_ele_x = 2;
    int num_ele_y = 2;

    int num_gp = 1;


    int num_ele = num_ele_x * num_ele_y;
    int num_dof = 2*((num_ele_x + 1)*(num_ele_y + 1));

    double len_ele_x = len_x/num_ele_x;
    double len_ele_y = len_y/num_ele_y;

    printf("Number of Elements: %d\n", num_ele);
    printf("Number of DOFs: %d\n", num_dof);

    double C[3*3];
    create_elasticity_matrix_isotropic_plane_stress(C, &e_modulus, &nu);
    print_array_col_major("Elasticity Matrix, C:", C, 3, 3);

    /*
    double Ke[8*8];
    double X[4*2];
    create_element_stiffness_matrix(Ke, C, &num_gp, X, &thickness);
    print_array_col_major("Element Stiffness Matrix, Ke", Ke, 8, 8);
    */

    double K[num_dof*num_dof];
    for (int i = 0; i < num_dof*num_dof; i++) {
        K[i] = 0.0;
    }
    create_global_stiffness_matrix(
        K, &num_dof, &num_ele, &num_ele_x, &len_ele_x, &len_ele_y, &num_gp, &thickness, C);
    print_array_col_major("Global Stiffness Matrix without B.C.'s, K:", K, num_dof, num_dof);

    double Kmod[num_dof*num_dof];
    memcpy(Kmod, K, sizeof(K));
    int bc[num_dof];
    for (int i = 0; i < num_dof; i++) {
        bc[i] = 0;
    }
    double F[num_dof];
    for (int i = 0; i < num_dof; i++) {
        F[i] = 0.0;
    }
    apply_bcs(Kmod, bc, F, &num_ele_x, &num_ele_y, &num_dof, &len_ele_x);
    print_array_col_major("Modified stiffness matrix, Kmod:", Kmod, num_dof, num_dof);
    print_array_col_major("Node force vector, F:", F, num_dof, 1);

    double U[num_dof];
    memcpy(U, F, sizeof(F));
    solve_for_displscements(Kmod, U, &num_dof);  /* Kmod*U = F */
    print_array_col_major("Displacement vector, U:", U, num_dof, 1);
    double solution_fem = U[2*(num_ele_x/2) + 1 + 0];  /* solution_fem = U[...][0] so in the middle at the bottom.*/ 
    double solution_ref = -0.3459;
    double epsilon = fabs(solution_fem - solution_ref)/solution_ref;
    printf("\nFEM Solution, u = %f\n", solution_fem);
    printf("Reference Solution, u = %f\n", solution_ref);
    printf("Epsilon Error = %f\n", epsilon);

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

int print_array_col_major_int(char *msg, int *array, int nrow, int ncol) {
   /* Print column-major arrays of one or two dimensions as 'normal' row-major arrays. */
    printf("%s\n", msg);
    printf("[%d x %d]\n", nrow, ncol);

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            printf("%10d \t", array[i + j*nrow]);
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

int create_element_stiffness_matrix(double *Ke, double *C, int *num_gp, double *X, double *thickness) {
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
    // int nrx = 4;
    // X[0 + 0*nrx] = 0;  /* X[0][0] = 0.0 */
    // X[0 + 1*nrx] = 0;  /* X[0][1] = 0.0 */
    // X[1 + 0*nrx] = 0.8;
    // X[1 + 1*nrx] = 0;
    // X[2 + 0*nrx] = 1;
    // X[2 + 1*nrx] = 1;
    // X[3 + 0*nrx] = 0;
    // X[3 + 1*nrx] = 1;
    // print_array_col_major("Nodal Coordinates, X:", X, 4, 2);

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
    info_dgetrf_Jinv = LAPACKE_dgetrf(LAPACK_COL_MAJOR, 2, 2, Jinv, 2, ipiv);
    /* printf("info_dgetrf_Jinv: %d\n", info_dgetrf_Jinv); */
    info_dgetri_Jinv = LAPACKE_dgetri(LAPACK_COL_MAJOR, 2, Jinv, 2, ipiv);
    /* printf("info_dgetri_Jinv: %d\n", info_dgetri_Jinv); */ 
    print_array_col_major("Inverse of the Jacobian Matrix, J^-1:",
    Jinv, 2, 2);

    /* Calculate the strain-displacement matrix, B */
    /* First we calculate dNidx and dNidy for each shape function,
     * then assenble B from dNidx and dNidy starting with an empty B matrix. */

    /* dNidx[2x4] = Jinv[2x2]*P[2*4] */
    double dNidx[2*4];
    for (int i = 0; i < 8; i++) {
        dNidx[i] = 0.0;
    }
    int nrdNidx = 2;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    2, 4, 2, 1.0, Jinv, 2, P, 2, 1.0, dNidx, 2);
    print_array_col_major("dNidx:", dNidx, 2, 4);
    
    /* Construct B */
    double B[3*8];
    for (int i = 0; i < 24; i++) {
        B[i] = 0.0;
    }
    int nrb = 3;
    for (int i = 0; i < 4; i++) {
       B[0 + (0+i*nrdNidx)*nrb] = dNidx[0 + i*nrdNidx];  /* First loop B[0][0] = dNidx[0][0] */ 
       B[2 + (0+i*nrdNidx)*nrb] = dNidx[1 + i*nrdNidx];  /* First loop B[2][0] = dNidx[1][0] */ 
       B[1 + (1+i*nrdNidx)*nrb] = dNidx[1 + i*nrdNidx];  /* First loop B[1][1] = dNidx[1][0] */ 
       B[2 + (1+i*nrdNidx)*nrb] = dNidx[0 + i*nrdNidx];  /* First loop B[2][1] = dNidx[0][0] */ 
    }
    print_array_col_major("Strain-Displacement Matrix, B:", B, 3, 8);

   /* Next we will calculate Ke via the integral of Btrans*C*B [8x8] over the volume */

    double BtransC[8*3];
    for (int i = 0; i < 24; i++) {
        BtransC[i] = 0.0;
    }
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
    8, 3, 3, 1.0, B, 8, C, 3, 1.0, BtransC, 8);
    print_array_col_major("BtransC:", BtransC, 8, 3);

    double BtransCB[8*8];
    for (int i = 0; i < 8*8; i++) {
        BtransCB[i] = 0.0;
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
    8, 8, 3, 1.0, BtransC, 8, B, 3, 1.0, BtransCB, 8);
    print_array_col_major("Matrix Btrans*C*B, BtransCB:", BtransCB, 8, 8);  /* Bug here? */

    /* Next we want to calculate Ke via 1-point numerical integration:
     * Ke = Btrans*C*B*Jdet*thickness*weight */
     
    /* Calculate the determinant of the Jacobi Matrix */
    double Jdet = 0.0;
    int nrj = 2;
    Jdet = J[0 + 0*nrj]*J[1 + 1*nrj] - J[0 + 1*nrj]*J[1 + 0*nrj];  /* det(A) = a*d - b*c */

    for (int i = 0; i < 8*8; i++) {
        Ke[i] = BtransCB[i]*Jdet*(*thickness)*weight; 
    }

    return 0;
}

int create_global_stiffness_matrix(
    double *K, int *nrk, int *num_ele, int *num_ele_x, double *len_ele_x, double *len_ele_y, int *num_gp, double *thickness, double *C) {
    /* Here comes the geometry of the model into play. */
    int eftab[8*(*num_ele)];
    int nreftab = 8;
    for (int i = 0; i < 8*(*num_ele); i++) {
        eftab[i] = 0;
    }

    double X[4*2];
    int nrx = 4;

    double Ke[8*8];
    int nrke = 8;

    for (int i = 0; i < *num_ele; i++) {
        int row = (i/(*num_ele_x)) + 1;
        int col = (i%(*num_ele_x)) + 1;

        eftab[0 + i*nreftab] = 2*((row - 1)*(*num_ele_x + 1) + col) - 1;  /* First loop eftab[0][0] = ... */
        eftab[1 + i*nreftab] = 2*((row - 1)*(*num_ele_x + 1) + col);
        eftab[2 + i*nreftab] = 2*((row - 1)*(*num_ele_x + 1) + col + 1) - 1;
        eftab[3 + i*nreftab] = 2*((row - 1)*(*num_ele_x + 1) + col + 1);
        eftab[4 + i*nreftab] = 2*(row*(*num_ele_x + 1) + col + 1) - 1;
        eftab[5 + i*nreftab] = 2*(row*(*num_ele_x + 1) + col + 1);
        eftab[6 + i*nreftab] = 2*(row*(*num_ele_x + 1) + col) - 1;
        eftab[7 + i*nreftab] = 2*(row*(*num_ele_x + 1) + col);

        X[0 + 0*nrx] = (*len_ele_x)*(col - 1);  /* x1: X[0][0] = ... */
        X[0 + 1*nrx] = (*len_ele_y)*(row - 1);  /* y1: X[0][1] = ... */
        X[1 + 0*nrx] = (*len_ele_x)*col;
        X[1 + 1*nrx] = (*len_ele_y)*(row - 1);  /* y2: X[1][1] = ... */
        X[2 + 0*nrx] = (*len_ele_x)*col; 
        X[2 + 1*nrx] = (*len_ele_y)*row; 
        X[3 + 0*nrx] = (*len_ele_x)*(col - 1);
        X[3 + 1*nrx] = (*len_ele_y)*row;
        printf("i = %d\n", i);
        print_array_col_major("Nodal Coordinates, X:", X, 4, 2);

        if (i == 0) {
            create_element_stiffness_matrix(Ke, C, num_gp, X, thickness);
            print_array_col_major("Element Stiffness Matrix, Ke", Ke, 8, 8);
        }

        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                K[eftab[j + i*nreftab] - 1 + (eftab[k + i*nreftab] - 1)*(*nrk)] += Ke[j + k*nrke];
            }
        }

    }
    /* Debugging */
    print_array_col_major_int("Last element freedom table, eftab:", eftab, 8, *num_ele);

    return 0;
}

int apply_bcs(
    double *Kmod, int *bc, double *F, int *num_ele_x, int *num_ele_y, int *num_dof, double *len_ele_x) {
    int start;
    int end;

    /* Dirichlet (Displacement) Boundary Conditions */
    for (int i = 0; i < (*num_ele_y + 1); i++) {
        bc[2*(*num_ele_x + 1)*i + 0] = 1;  /* First loop bc[0][0] = 1 */
        bc[(2*(*num_ele_x + 1)*i + 1) + 0] = 1;  /* First loop bc[1][0] = 1 */
        bc[2*(*num_ele_x + 1)*i + 2*(*num_ele_x) + 0] = 1;
        bc[(2*(*num_ele_x + 1)*i + 2*(*num_ele_x) + 1) + 0] = 1;
    }

    /* print_array_col_major_int("Boundary Conditions, bc:", bc, *num_dof, 1); */

    for (int i = 0; i < *num_dof; i++) {
        if (bc[i] == 1) {
            for (int j = 0; j < *num_dof; j++) {
                Kmod[i + j*(*num_dof)] = 0.0;  /* First loop Kmod[0][0] = 0.0 */
                Kmod[j + i*(*num_dof)] = 0.0;
            }
            Kmod[i + i*(*num_dof)] = 1.0;
        }
    }

    start = 2*(*num_ele_x + 1)*(*num_ele_y) + 4; 
    end = 2*(*num_ele_x + 1)*(*num_ele_y) + 2*(*num_ele_x) + 1;
    /* printf("start, %d, end, %d\n", start, end); */

    for (int i = start; i < end; ) {
        F[i - 1] = -1*(*len_ele_x);
        i += 2;
    }

    return 0;
}

int solve_for_displscements(double *Kmod, double *U, int *num_dof) {
    /* Kmod*U = F */
    lapack_int info_dgels_U = 0;

    info_dgels_U = LAPACKE_dgels(
        LAPACK_COL_MAJOR, 'N', *num_dof, *num_dof, 1, Kmod, *num_dof, U, *num_dof);
    printf("info_dgels_U: %d\n", info_dgels_U);

    return 0; 
}
