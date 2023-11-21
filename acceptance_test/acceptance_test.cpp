// This reproducer of a double free bug in HPE-MPI is based on the following
// HDF5 example:
//   https://support.hdfgroup.org/ftp/HDF5/examples/parallel/Hyperslab_by_chunk.c

#include <mpi.h>
#include <hdf5.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <iostream>

#define H5FILE_NAME     "collective.h5"
#define DATASETNAME     "dset" 
#define NX     1024
#define RANK   1

int
main (int argc, char **argv)
{
    // MPI variables
    int mpi_size, mpi_rank;
    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);  
 
    // -- Setup --------------------------------------------------------------------
    // We create a simple file with one, one-dimensional dataset `dset`.
    auto fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, comm, info);

    // Create a new file collectively and release property list identifier.
    hid_t file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    H5Pclose(fapl_id);

    // Create the dataspace for the dataset.
    hsize_t dims[RANK] = {NX};
    hid_t filespace = H5Screate_simple(RANK, dims, NULL); 
    std::vector<int> array(NX);
    for(size_t i = 0; i < array.size(); ++i) {
        array[i] = i;
    }

    // Create the dataset with default properties.
    hid_t dset_id = H5Dcreate(
        file_id, DATASETNAME, H5T_NATIVE_INT, filespace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
    );

    H5Dwrite(dset_id, H5T_NATIVE_INT, filespace, filespace,
                      H5P_DEFAULT, array.data());

    // -- Reproducer ----------------------------------------------------------------
    // We will read a irregular hyperslab from `dset` using collective MPI-IO.
    std::vector<std::array<hsize_t, 2>> ranges{{0, 16}, {100, 16}, {128, 16}};
    hsize_t memspace_size = 0;
    for(const auto& range : ranges) {
        memspace_size += range[1];
    }

    hid_t memspace = H5Screate_simple(RANK, &memspace_size, NULL);

    // Select hyperslab in the file.
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &ranges[0][0], NULL, &ranges[0][1], NULL);
    for(size_t i = 1; i < ranges.size(); ++i) {
        H5Sselect_hyperslab(filespace, H5S_SELECT_OR, &ranges[i][0], NULL, &ranges[i][1], NULL);
    }

    std::vector<int> subarray(memspace_size, -1);

    // Create property list for collective dataset write.
    hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);

    // This call fails when using HPE-MPI, other implementations work flawlessly.
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, dxpl_id, subarray.data());

    // -- Acceptance tests -------------------------------------------------------
    uint32_t local_cause;
    uint32_t global_cause;
    H5Pget_mpio_no_collective_cause(dxpl_id, &local_cause, &global_cause);

    if(local_cause != 0 || global_cause != 0) {
        std::cout << "Broken: local_cause = " << local_cause << ", global_cause = " << global_cause << std::endl;
    }

    {
        size_t i = 0;
        for(const auto& range : ranges) {
            for(size_t j = range[0]; j < range[0] + range[1]; ++j) {
                if(subarray[i] != array[j]) {
                    std::cout << "Broken: i = " << i << ", j = " << j << ", subarray[i] = " << subarray[i] << std::endl;
                }
                ++i;
            }
        }
    }

    // -- Clean up ---------------------------------------------------------------
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(dxpl_id);
    H5Fclose(file_id);
 
    MPI_Finalize();

    return 0;
}
