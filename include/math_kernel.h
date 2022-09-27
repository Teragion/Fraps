#ifndef MATH_KERNEL_H
#define MATH_KERNEL_H

#include "array2.h"
#include "vec.h"

#include <mkl_cblas.h>

#include <omp.h>

/**
 * @brief computes dot product of vectors
 * 
 * @tparam Vec 
 * @tparam Scalar=Vec::value_type 
 * @param a 
 * @param b 
 * @return double 
 */
template<typename Vector, typename Scalar=typename Vector::value_type>
double dot_product(const Vector& a, const Vector& b) {
   double res;
   for (int i = 0; i < Vector::size(); i++) {
      res += a(i) * b(i);
   }
   return res;
}

// overloading with BLAS
template<unsigned int N>
double dot_product(const Vec<N, double>& a, const Vec<N, double>& b) {
   double res;
   res = cblas_ddot(a.dim(), a.cdata(), 1, b.cdata(), 1);
   return res;
}

/**
 * @brief a = a * scale
 * 
 * @param a 
 * @param scale 
 */
template<class T, class S>
void mat_scale(Array2<T, Array1<T> > &a, const S& scale) {
#pragma omp parallel for
   for (int j = 0; j < a.nj; ++j) {
#pragma omp simd
      for (int i = 0; i < a.ni; ++i) {
         a(i, j) *= scale;
      }
   }
}

/**
 * @brief c = a + b
 * 
 * @tparam T Scalar type 
 * @param a matrix
 * @param b matrix
 * @param c matrix
 */
template<class T>
void mat_sum(const Array2<T, Array1<T> > &a, const Array2<T, Array1<T> > &b, Array2<T, Array1<T> > &c) {
#pragma omp parallel for
   for (int j = 0; j < a.nj; ++j) {
#pragma omp simd
      for (int i = 0; i < a.ni; ++i) {
         c(i, j) = a(i, j) + b(i, j);
      }
   }
}

/**
 * @brief c = a * alpha + b
 * 
 * @param a matrix
 * @param b matrix
 * @param c matrix
 * @param alpha scalar
 */
template<class T, class S>
void saxpy(const Array2<T, Array1<T> > &a, const Array2<T, Array1<T> > &b, Array2<T, Array1<T> > &c, const S& alpha) {
#pragma omp parallel for
   for (int j = 0; j < a.nj; ++j) {
#pragma omp simd
      for (int i = 0; i < a.ni; ++i) {
         c(i, j) = a(i, j) * alpha + b(i, j);
      }
   }
}

#endif // MATH_KERNEL_H
