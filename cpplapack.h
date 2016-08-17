/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Copyright (C) 2014 Claudio Berti



-----------------------------------------------------------------------------
Claudio Berti, Ph.D.,
Department of Molecular Biophysics and Physiology,
Rush University Medical Center
and
University of Bologna, Bologna, Italy

Tel: 312.942.6756
Fax: 312.942.8711
E-Mail: brt.cld@gmail.com
-----------------------------------------------------------------------------		



This file is part of BROWNIES.

BROWNIES is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or 
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Private, research, and institutional use is free under the condition of
acknowledge the author and contributors citing the following papers:

- Berti, C.; Furini, S.; Gillespie, PACO: PArticle COunting Method To Enforce
Concentrations in Dynamic Simulations. J. Chem. Theory Comput. (in press)

- Berti, C.; Furini, S.; Gillespie, D.; Boda, D.; Eisenberg, R. S.; 
Sangiorgi, E; Fiegna, A 3-D Brownian dynamics simulator for the study of ion 
permeation through membrane pores J. Chem. Theory Comput. 2014, 10, 2911-2926

- Berti, C.; Gillespie, D.; Bardhan, J. P.; Eisenberg, R. S.; Fiegna, C.,
Comparison of three-dimensional Poisson solution methods for particle-based 
simulation and inhomogeneous dielectrics. Phys. Rev. E 2012, 86, 011912

- Berti, C.; Gillespie, D.; Eisenberg, R.; Fiegna, C. Particle-based 
simulation of charge transport in discrete-charge nano-scale systems: the 
electrostatic problem. Nanoscale Research Letters 2012, 7, 135



You may distribute modified versions of this code UNDER THE CONDITION THAT
THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE SAME FILE REMAIN UNDER 
COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE MADE 
FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE 
MODIFICATIONS.  Distribution of this code as part of a commercial system
is permissible ONLY BY DIRECT ARRANGEMENT WITH THE AUTHOR. (If you are 
not directly supplying this code to a customer, and you are instead 
telling them how they can obtain it for free, then you are not required
to make any arrangement with me.) 

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


/////////////////////////////////////////////////////////////////////
//		BLAS
/////////////////////////////////////////////////////////////////////
extern "C" void dgemv_(char *TRANS,int *M, int *N, double *ALPHA, double *A, int *LDA, double *X, int *INCX, double *BETA,double *Y, int *INCY);

extern "C" double dnrm2_(int *N,double *X, int *INCX);
/////////////////////////////////////////////////////////////////////
//		LAPACK
/////////////////////////////////////////////////////////////////////
extern "C" void dsyevd_(char *JOBZ,char *UPLO,int *N,double *A,int *LDA,double *W
		,double *WORK,int *LWORK
		,int *IWORK,int *LIWORK,int *INFO );

// DGELSY computes the minimum-norm solution to a real linear least
// squares problem: minimize || A * X - B ||
extern "C" void dgelsy_(int *M, int *N, int *NRHS, double  *A, int *LDA,
		double *B, int *LDB,int *JPVT, double *RCOND, int *RANK,
		double *WORK, int *LWORK, int *INFO );
// DGELS solves overdetermined or underdetermined real linear systems
// involving an M-by-N matrix A, or its transpose, using a QR or LQ
// factorization of A.  It is assumed that A has full rank.
extern "C" void dgels_(char *TRANS, int *M,int *N, int *NRHS,
		double *A, int *LDA, double *B,int *LDB,
		double *WORK,int *LWORK,int *INFO);
		
// DGETRF computes an LU factorization of a general M-by-N matrix A
// using partial pivoting with row interchanges.
//
// The factorization has the form
//    A = P * L * U
// where P is a permutation matrix, L is lower triangular with unit
// diagonal elements (lower trapezoidal if m > n), and U is upper
// triangular (upper trapezoidal if m < n).
//
// This is the right-looking Level 3 BLAS version of the algorithm.
extern "C" void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);



// DGETRI computes the inverse of a matrix using the LU factorization
// computed by DGETRF.
//
// This method inverts U and then computes inv(A) by solving the system
// inv(A)*L = inv(U) for inv(A).
extern "C" void dgetri_(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO );


// DGEMV  performs one of the matrix-vector operations
// y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
//  where alpha and beta are scalars, x and y are vectors and A is an
// m by n matrix.
extern "C" void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *x, int *incx, double *beta, double *y, int *incy);


// DGEMM  performs one of the matrix-matrix operations
//     C := alpha//op( A )//op( B ) + beta//C,
//  where  op( X ) is one of
//     op( X ) = X   or   op( X ) = X',
//  alpha and beta are scalars, and A, B and C are matrices, with op( A )
//  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
extern "C" void   dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);


//  DNRM2 returns the euclidean norm of a vector via the function
//  name, so that
//  DNRM2 := sqrt( x'*x )
extern "C" double dnrm2_(int *n,double *x, int* incx);


//  DGESV computes the solution to a real system of linear equations
//     A // X = B,
//  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
//
//  The LU decomposition with partial pivoting and row interchanges is
//  used to factor A as
//     A = P // L // U,
//  where P is a permutation matrix, L is unit lower triangular, and U is
//  upper triangular.  The factored form of A is then used to solve the
//  system of equations A // X = B.
extern "C" void   dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info );


// copies a vector, x, to a vector, y.
// uses unrolled loops for increments equal to one.
extern "C" void   dcopy_(int *n, double *x, int *incx, double *y, int *incy);
