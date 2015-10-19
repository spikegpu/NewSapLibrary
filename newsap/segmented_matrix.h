#ifndef SEGMENTED_MATRIX_H
#define SEGMENTED_MATRIX_H

#include <newsap/common.h>

#include <cusp/array1d.h>

namespace newsap {

template <typename T, typename MemorySpace>
class SegmentedMatrix {
private:
    typedef typename cusp::array1d<T, MemorySpace>          Vector;
    typedef typename cusp::array1d<T, cusp::host_memory>    VectorH;
    typedef typename cusp::array1d<int, MemorySpace>        IntVector;
    typedef typename cusp::array1d<int, cusp::host_memory>  IntVectorH;

    unsigned    m_num_partitions;
    IntVectorH  m_num_rows_host;
    IntVectorH  m_num_columns_host;
    IntVectorH  m_A_offsets;
    IntVectorH  m_B_offsets;
    IntVectorH  m_C_offsets;
    Vector      m_A;
    Vector      m_B;
    Vector      m_C;

public:
    SegmentedMatrix(
        unsigned            num_partitions,
        const IntVectorH&   num_rows,
        const IntVectorH&   num_columns,
        const IntVectorH&   A_offsets,
        const IntVectorH&   B_offsets,
        const IntVectorH&   C_offsets,
        const VectorH&      A,
        const VectorH&      B,
        const VectorH&      C
    );
};

template <typename T, typename MemorySpace>
SegmentedMatrix<T, MemorySpace>::SegmentedMatrix(
        unsigned            num_partitions,
        const IntVectorH&   num_rows,
        const IntVectorH&   num_columns,
        const IntVectorH&   A_offsets,
        const IntVectorH&   B_offsets,
        const IntVectorH&   C_offsets,
        const VectorH&      A,
        const VectorH&      B,
        const VectorH&      C
    ): m_num_partitions(num_partitions),
       m_num_rows_host(num_rows),
       m_num_columns_host(num_columns),
       m_A_offsets(A_offsets),
       m_B_offsets(B_offsets),
       m_C_offsets(C_offsets),
       m_A(A),
       m_B(B),
       m_C(C)
{
}


}

#endif // SEGMENTED_MATRIX_H
