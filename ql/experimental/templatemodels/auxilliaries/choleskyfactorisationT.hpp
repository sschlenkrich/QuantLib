/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016 Sebastian Schlenkrich

*/



#ifndef quantlib_templateauxilliaries_choleskyfactorisation_hpp
#define quantlib_templateauxilliaries_choleskyfactorisation_hpp

//#include <ql/types.hpp>
//#include <boost/function.hpp>

#include <stdlib.h>
#include <math.h>


namespace TemplateAuxilliaries {

    // This method is unsafe! Memory allocation need to be ensured by user
    // A(i,j) (and L) is stored row-wise as A[i*n+j]
    // L(i,j) is lower triangular matrix with L L^T = A
    template <typename Type> 
    inline void cholesky(std::vector< Type >& A, std::vector< Type >& L, size_t n) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < (i+1); ++j) {
                Type s = 0;
                for (size_t k = 0; k < j; ++k) s += L[i * n + k] * L[j * n + k];
                if (i==j) {
                    if (A[i * n + i] < s) throw std::exception();
                    L[i * n + j] = sqrt(A[i * n + i] - s);
                }
                else L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));
            }
        } 
    }

    template <typename Type>
    std::vector< std::vector< Type > > cholesky(const std::vector< std::vector< Type > >& A) {
        std::vector<Type> arrayA(A.size()*A.size());
        std::vector<Type> arrayL(A.size()*A.size());
        for (size_t i = 0; i < A.size(); ++i) {
            if (A.size()!=A[i].size()) throw std::exception();
            for (size_t j = 0; j < A[i].size(); ++j) arrayA[i*A[i].size() + j] = A[i][j];
        }
        cholesky(arrayA, arrayL, A.size());
        std::vector< std::vector< Type > > L(A.size(), std::vector< Type >(A.size(), 0.0));
        for (size_t i = 0; i < L.size(); ++i) {
            for (size_t j = 0; j <= i; ++j) L[i][j] = arrayL[i*L[i].size() + j];
        }
        return L;
    }


    // alternative implementation of Cholesky decomposition
    template<class T>
    void performCholesky(std::vector< std::vector<T> >& matrix, size_t dimIn, bool flexible) {
        size_t dim = dimIn;
        for (size_t i = 0; i < dim; i++) {
            for (size_t j = i; j < dim; j++) {
                if (abs(matrix[i][j] - matrix[j][i])>QL_EPSILON)
                    QL_FAIL(std::string("A symmetrix matrix is necessary to apply Cholesky Decomposition"));
            }
        }

        //Perform Decomposition as described in script of Trottenberg (S. 47):
        //Because script delivers Upper right matrix, we interchange indizes in order to end up with lower left matrix.

        //First step:
        matrix[0][0] = sqrt(matrix[0][0]);
        if (abs(matrix[0][0]) < QL_EPSILON && !flexible)
            QL_FAIL("No positive definite Correlation Matrix because rank is not full.");
        for (size_t i = 1; i < dim; i++) {
            matrix[i][0] = abs(matrix[0][0]) < QL_EPSILON ? 0.0 : matrix[i][0] / matrix[0][0];
        }

        //Now iterate:
        for (size_t i = 1; i < dim; i++) {
            for (size_t k = 0; k < i; k++) {
                matrix[i][i] = matrix[i][i] - matrix[i][k] * matrix[i][k];
            }
            matrix[i][i] = sqrt(matrix[i][i]);
            if (abs(matrix[i][i]) < QL_EPSILON && !flexible)
                QL_FAIL(std::string("No positive definite Correlation Matrix because rank is not full."));
            if (matrix[i][i] != matrix[i][i] && !flexible) {
                QL_FAIL(std::string("No positive definite Correlation Matrix as diagonal square entry is negative."));
            }
            matrix[i][i] = (matrix[i][i] != matrix[i][i]) ? 0.0 : matrix[i][i];
            for (size_t j = 0; j < dim; j++) {
                if (j >= i + 1) {
                    for (size_t k = 0; k < i; k++) {
                        matrix[j][i] = matrix[j][i] - matrix[i][k] * matrix[j][k];
                    }
                    matrix[j][i] = (abs(matrix[i][i]) < QL_EPSILON || matrix[i][i] != matrix[i][i]) ? 0.0 : matrix[j][i] / matrix[i][i];
                }
                else if (j<i) {
                    matrix[j][i] = 0;
                }
                else {
                    //Nothing because diagonale.
                }
            }
        }
    }


}

#endif  /* ifndef quantlib_templateauxilliaries_qrfactorisation_hpp */
