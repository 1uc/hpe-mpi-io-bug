# Double free in `MPI_Type_free`.

## Reproducer
The reproducer can be found in the subfolder `reproducer`.

## Acceptance Test
To confirm that the fix solves the issue we've adjust an example from the HDF5
community for reading irregular hyperslabs. The code including a CMake file can
be found in the subfolder `acceptance_test`.

