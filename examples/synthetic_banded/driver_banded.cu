#include <iostream>
#include <cmath>
#include <stdlib.h>

#include <newsap/common.h>

#ifdef USE_OLD_CUSP
#  include <cusp/blas.h>
#else
#  include <cusp/blas/blas.h>
#endif
#include <cusp/array1d.h>
#include <cusp/print.h>

#include <newsap/segmented_matrix.h>

typedef double REAL;
typedef typename cusp::array1d<REAL, cusp::device_memory>             Vector;
typedef typename cusp::array1d<REAL, cusp::host_memory>               VectorH;
typedef typename cusp::array1d<int, cusp::device_memory>              IntVector;
typedef typename cusp::array1d<int, cusp::host_memory>                IntVectorH;
typedef typename newsap::SegmentedMatrix<REAL, cusp::device_memory>   Matrix;

using std::cout;
using std::cerr;
using std::cin;
using std::endl;

// -----------------------------------------------------------------------------
// Macro to obtain a random number between two specified values
// -----------------------------------------------------------------------------
#define RAND(L, H) ((L) + ((H) - (L)) * (float)rand() / (float)RAND_MAX)


// -----------------------------------------------------------------------------
// Definitions for SimpleOpt and SimpleGlob
// -----------------------------------------------------------------------------
#include <SimpleOpt/SimpleOpt.h>

// ID values to identify command line arguments
enum {OPT_HELP, OPT_BAND};

// Table of CSimpleOpt::Soption structures. Each entry specifies:
// - the ID for the option (returned from OptionId() during processing)
// - the option as it should appear on the command line
// - type of the option
// The last entry must be SO_END_OF_OPTIONS
CSimpleOptA::SOption g_options[] = {
	{ OPT_BAND,          "-b",                   SO_MULTI   },
	{ OPT_BAND,          "--banded-synthetic",   SO_MULTI   },
	{ OPT_HELP,          "-?",                   SO_NONE    },
	{ OPT_HELP,          "-h",                   SO_NONE    },
	{ OPT_HELP,          "--help",               SO_NONE    },
	SO_END_OF_OPTIONS
};

void ShowUsage();
bool
GetProblemSpecs(int             argc, 
                char**          argv,
                int&            N,
                int&            k,
                REAL&           d);

void
GetSegmentedMatrices(
    int N,
    int k,
    REAL d,
    IntVectorH& num_rows,
    IntVectorH& num_columns,
    IntVectorH& A_offsets,
    IntVectorH& B_offsets,
    IntVectorH& C_offsets,
    VectorH&    A,
    VectorH&    B,
    VectorH&    C
);

int main(int argc, char **argv) {
    int            pN;
    int            pk;
    REAL           pd;

	if (!GetProblemSpecs(argc, argv, pN, pk, pd)) {
        ShowUsage();
        return 1;
    }

    unsigned num_partitions = pN / pk;

    IntVectorH       num_rows, num_columns;
    IntVectorH       A_offsets, B_offsets, C_offsets;
    VectorH          subA, subB, subC;

    GetSegmentedMatrices(
        pN,
        pk,
        pd,
        num_rows,
        num_columns,
        A_offsets,
        B_offsets,
        C_offsets,
        subA,
        subB,
        subC
    );
    cusp::print(A_offsets);
    cusp::print(subA);
    cusp::print(B_offsets);
    cusp::print(subB);
    cusp::print(C_offsets);
    cusp::print(subC);

    Matrix A(num_partitions, num_rows, num_columns, A_offsets, B_offsets, C_offsets, subA, subB, subC);

    return 0;
}

void ShowUsage()
{
	cout << "Usage:  driver_mm [OPTIONS]" << endl;
	cout << endl;
	cout << " -b SIZE BW DD" << endl;
	cout << " --banded-synthetic SIZE BW DD" << endl;
	cout << "        Use a synthetic banded matrix of size SIZE, half-bandwidth BW," << endl;
	cout << "        and degree of diagonal dominance DD." << endl;
}

// -----------------------------------------------------------------------------
// GetProblemSpecs()
//
// This function parses the specified program arguments and sets up the problem
// to be solved.
// -----------------------------------------------------------------------------
bool
GetProblemSpecs(int             argc, 
                char**          argv,
                int&            N,
                int&            k,
                REAL&           d)
{
    N = -1;
    k = -1;
    d = -1.0;
	// Create the option parser and pass it the program arguments and the array
	// of valid options. Then loop for as long as there are arguments to be
	// processed.
	CSimpleOptA args(argc, argv, g_options);

	while (args.Next()) {
		// Exit immediately if we encounter an invalid argument.
		if (args.LastError() != SO_SUCCESS) {
			cout << "Invalid argument: " << args.OptionText() << endl;
			ShowUsage();
			return false;
		}

		// Process the current argument.
		switch (args.OptionId()) {
			case OPT_HELP:
				return false;
			case OPT_BAND:
				{
					char **mArgs = args.MultiArg(3);
					if (!mArgs) {
						return false;
					}
					N = atoi(mArgs[0]);
					k = atoi(mArgs[1]);
					d = atof(mArgs[2]);

					break;
				}
		}
	}

    if (N < 0 || k < 0 || d < 0) {
        ShowUsage();
        return false;
    }

	return true;
}

void
GetSegmentedMatrices(
    int N,
    int k,
    REAL d,
    IntVectorH& num_rows,
    IntVectorH& num_columns,
    IntVectorH& A_offsets,
    IntVectorH& B_offsets,
    IntVectorH& C_offsets,
    VectorH&    A,
    VectorH&    B,
    VectorH&    C
) {
    int num_partitions = N / k;
    int remainder = N % num_partitions;

    num_rows.resize(num_partitions);
    num_columns.resize(num_partitions);
    A_offsets.resize(num_partitions + 1);
    B_offsets.resize(num_partitions + 1);
    C_offsets.resize(num_partitions + 1);

    int cur_a_offset = 0;
    int cur_b_offset = 0;
    int cur_c_offset = 0;

    for (int i = 0; i < num_partitions; i++) {
        A_offsets[i] = cur_a_offset;
        B_offsets[i] = cur_b_offset;
        C_offsets[i] = cur_c_offset;
        if (i < remainder) {
            num_rows[i] = num_columns[i] = k + 1;
            cur_a_offset += (k + 1) * (k + 1);
            cur_b_offset += ((i > 0) ? ((k + 1) * (k + 1)) : 0);
            cur_c_offset += (k + 1) * ((i < remainder - 1) ? (k + 1) : k);
        } else {
            num_rows[i] = num_columns[i] = k;
            cur_a_offset += k * k;
            cur_b_offset += k * ((i > remainder) ? k : (i == 0 ? 0 : (k + 1)));
            cur_c_offset += k * ((i == (num_partitions - 1)) ? 0 : k);
        }
    }

    A_offsets[num_partitions] = cur_a_offset;
    B_offsets[num_partitions] = cur_b_offset;
    C_offsets[num_partitions] = cur_c_offset;

    A.resize(cur_a_offset, REAL(0));
    B.resize(cur_b_offset, REAL(0));
    C.resize(cur_c_offset, REAL(0));

    for (int i = 0; i < num_partitions; i++) {
        int start_row, end_row;
        int a_num_columns, b_num_columns, c_num_columns;
        if (i < remainder) {
            start_row = (k + 1) * i;
            end_row = start_row + (k + 1);
            a_num_columns = k + 1;
            b_num_columns = (i > 0 ? (k + 1) : 0);
            c_num_columns = (i < remainder - 1 ? (k + 1) : k);
        } else {
            start_row = k * i + remainder;
            end_row = start_row + k;
            a_num_columns = k;
            b_num_columns = (i > remainder ? k : (i == 0 ? 0 : (k + 1)));
            c_num_columns = (i < (num_partitions - 1) ? k : 0);
        }

        REAL my_sum = REAL(0);
        for (int j = start_row; j < end_row; j++) {
            int column_min = j - k;
            int column_max = j + k;
            if (column_min < 0) {
                column_min = 0;
            }
            if (column_max >= N) {
                column_max = N - 1;
            }

            for (int l = column_min; l <= column_max; l++) {
                REAL value = RAND(-10.0, 10.0);
                my_sum += (value < 0 ? -value : value);

                if (l < start_row) {
                    B[B_offsets[i] + (j - start_row + 1) * b_num_columns + l - start_row] = value;
                } else if (l >= end_row) {
                    C[C_offsets[i] + (j - start_row) * c_num_columns + l - end_row] = value;
                } else {
                    A[A_offsets[i] + (j - start_row) * a_num_columns + l - start_row] = value;
                }
            }
            A[A_offsets[i] + (j - start_row) * a_num_columns + j - start_row] = my_sum * d;
        }
    }
}
