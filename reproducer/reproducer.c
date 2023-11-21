#include <mpi.h>
#include <stdio.h>

// Compile  : mpicc -DBB5_WITH_BUG=1 $FILE_NAME -o $OUTPUT_NAME
// Run (v1) : mpirun -n 1 $OUTPUT_NAME
// Run (v2) : $OUTPUT_NAME

#define BB5_CHECK(mpi_code) if((mpi_code) != MPI_SUCCESS) { printf("MPI failed."); }

int main (int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int mpi_code;

    MPI_Datatype inner_type;
    MPI_Datatype intermediate_type;
    MPI_Datatype outer_type;

    int inner_count = 1;
    int inner_counts[] = {10};
    MPI_Aint inner_offsets[] = {16 * sizeof(int)};

    mpi_code = MPI_Type_create_hindexed(
        inner_count,
        inner_counts,
        inner_offsets,
        MPI_INT,
        &inner_type
    );                                                                  BB5_CHECK(mpi_code)

#if BB5_WITH_BUG == 1
    mpi_code = MPI_Type_dup(inner_type, &intermediate_type);            BB5_CHECK(mpi_code)
    mpi_code = MPI_Type_free(&inner_type);                              BB5_CHECK(mpi_code)
#else
    intermediate_type = inner_type;
#endif

    int outer_count = 1;
    int outer_counts[] = {1};
    MPI_Aint outer_offset[] = {0};

    mpi_code = MPI_Type_create_struct(
        outer_count,
        outer_counts,
        outer_offset,
        &intermediate_type,
        &outer_type
    );                                                                  BB5_CHECK(mpi_code)

    mpi_code = MPI_Type_free(&intermediate_type);                       BB5_CHECK(mpi_code)
    mpi_code = MPI_Type_free(&outer_type);                              BB5_CHECK(mpi_code)

    MPI_Finalize();
}
