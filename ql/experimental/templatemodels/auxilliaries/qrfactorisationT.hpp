/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016 Sebastian Schlenkrich

*/



#ifndef quantlib_templateauxilliaries_qrfactorisation_hpp
#define quantlib_templateauxilliaries_qrfactorisation_hpp

//#include <ql/types.hpp>
//#include <boost/function.hpp>


namespace TemplateAuxilliaries {

    template <typename Type> 
    inline void givensrotation(std::vector<Type>&                  x,
                               const size_t                        i, 
                               const size_t                        j, // indicies leq # rows of x
                               const Type                          a,
                               const Type                          b ) {
        //
        //  x <- G * x  where
        //       [ 1         ]
        //       [   c  -s   ] <- i       c = a / (a²+b²)^{1/2}
        //  G  = [     1     ]
        //       [   s   c   ] <- j       s = b / (a²+b²)^{1/2}
        //       [         1 ]
        //
        Type c, s, den, y, w;
        if (a==0.0) {
            c = 0.0;
            s = (b>0) ? (1.0) : (-1.0);  // msign(b);
        } else {
            den = sqrt(a*a+b*b);
            if (den<1.0E-100) {
                c = 1.0;
                s = 0.0;
            } else {
                c = a / den;
                s = b / den;
            }
        }
        // multiply rows
        y    = x[i];
        w    = x[j];
        x[i] = c*y - s*w;
        x[j] = s*y + c*w;
    }


    template <typename Type> 
    inline void givensrotation(std::vector<std::vector<Type> >&    M,
                               const size_t                        i, 
                               const size_t                        j, // indicies leq # rows of M
                               const Type                          a,
                               const Type                          b ) {
        //
        //  M <- G * M  where
        //       [ 1         ]
        //       [   c  -s   ] <- i       c = a / (a²+b²)^{1/2}
        //  G  = [     1     ]
        //       [   s   c   ] <- j       s = b / (a²+b²)^{1/2}
        //       [         1 ]
        //
        Type c, s, den, y, w;
        if (a==0.0) {
            c = 0.0;
            s = (b>0) ? (1.0) : (-1.0);  // msign(b);
        } else {
            den = sqrt(a*a+b*b);
            if (den<1.0E-100) {
                c = 1.0;
                s = 0.0;
            } else {
                c = a / den;
                s = b / den;
            }
        }
        // multiply rows
        for (size_t k=0; k<M[i].size(); ++k) {
            y = M[i][k];
            w = M[j][k];
            M[i][k] = c*y - s*w;
            M[j][k] = s*y + c*w;
        }
    }


    template <typename Type>                                            // in place QR solution of ||Mx - b||->min
    inline void qrsolveles( std::vector<std::vector<Type> >&   M,     // input dim1-by-dim2 matrix, output overwritten
                            std::vector<Type>&                 b ) {  // input dim-1 vector b, output first dim2 elements x
        size_t dim1 = M.size();
        QL_REQUIRE(dim1>0,"QR factorisation wrong dimensions");
        size_t dim2 = M[0].size();
        QL_REQUIRE(dim2>0,"QR factorisation wrong dimensions");
        QL_REQUIRE(dim1>=dim2,"QR factorisation wrong dimensions");

        // factorisation
        for (size_t j=0; j<dim2; ++j) {  // iterate cols of M
            for (size_t i=dim1-1; i>j; --i) { // iterate rows of M
                if (M[i][j] != 0) {
                    givensrotation(b,j,i,M[j][j],-M[i][j]);
                    givensrotation(M,j,i,M[j][j],-M[i][j]);
                    M[i][j] = 0;
                }
            }
        }

        // back-subtitution
        for (size_t i=dim2; i>0; --i) {
            Type sum = 0.0;
            for (size_t j=i; j<dim2; ++j) sum += M[i-1][j]*b[j];
            b[i-1] = (b[i-1] - sum)/M[i-1][i-1];
        }
    }

    template <typename Type>
    inline std::vector<std::vector<Type> > qrinverse(const std::vector<std::vector<Type> >&  M) {
        std::vector<std::vector<Type> > X;
        for (size_t i = 0; i < M.size(); ++i) {
            QL_REQUIRE(M[i].size()==M.size(), "QR dimensions miss-match");
            std::vector<std::vector<Type> > A(M);  // we need to copy due to in-place factorisation
            std::vector<Type> b(M.size(), 0.0);
            b[i] = 1.0;
            qrsolveles(A, b);
            X.push_back(b);
        }
        std::vector<std::vector<Type> > Minv(M.size(), std::vector<Type>(M.size(), 0.0));
        for (size_t i = 0; i < M.size(); ++i)  // we need to transpose X coz matrix is stored row-wise
            for (size_t j = 0; j < M.size(); ++j)
                Minv[i][j] = X[j][i];
        return Minv;
    }

}

#endif  /* ifndef quantlib_templateauxilliaries_qrfactorisation_hpp */
