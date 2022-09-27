#ifndef MAT_H
#define MAT_H

#include "vec.h"

/**
 * @brief minimum wrapper around vec
 * 
 * @tparam N row dim 
 * @tparam M col dim
 * @tparam T 
 */
template<unsigned int N, unsigned int M, class T>
struct Mat {
    typedef Vec<M * N, T> vector_type;

    Vec<M * N, T> v;

    Mat<N, M, T>() {}

    Mat<N, M, T>(T value_for_all) : v(value_for_all) {}

    template<class S>
    Mat<N, M, T>(const S* source) {
        for (unsigned int i = 0; i < M * N; ++i) v[i] = (T)source[i];
    }

    template<class S>
    explicit Mat<N, M, T>(const Mat<N, M, T>& source) {
        for (unsigned int i = 0; i < M * N; ++i) v[i] = (T)source[i];
    }

    T& operator()(int i) {
        return v(i);
    }

    const T& operator()(int i) const {
        return v(i);
    }

    T& operator()(int i, int j) {
        return v(i + j * N);
    }    
    
    const T& operator()(int i, int j) const {
        return v(i + j * N);
    }

    const T* cdata() const {
        return v.cdata();
    }

    T* data() {
        return v.data();
    }

    void set_zero() {
        v.set_zero();
    }

    void swap(Mat<N, M, T>& m) {
        v.swap(m.v);
    }

    static constexpr unsigned int size() { return vector_type::size(); }
};

#endif // MAT_H
