# Set-up Instructions

1. Clone this repo, the run `git submodule init` and `git submodule update` to fill Petsc3 submodule directory
2. Compile fftw and PETSc3 rom their root directories: fftw-2.1.5 with `./configure --with-pic && make`; Petsc3 with `./configure && make all check`
3. In all makefile.opt files, (CurveLab-2.1.3/fdct3d_mpi/src/makefile.opt, CurveLab-2.1.3/makefile.opt) set `FFTW_DIR =  <curveletTF_directory>/fftw-2.1.5` (or wherever your fftw installation is)
4. In CurveLab-2.1.3/fdct3d_mpi/src/makefile.opt, set `PETSC_DIR =  <curveletTF_directory>/petsc` (or wherever your PETSc3 installation is)
5. In CurveLab-2.1.3 build and compile 2D and 3D (in-core, out-of-core) Curvelet TF implementations with `make lib`, then test CPP installations with `make test`
6. To compile MPI 3D CTF, go into CurveLab-2.1.3/fdct3d_mpi/src/makefile.opt, then run 'make libfdct3d.a' to compile the library, then test CPP installation with
`mpirun -np <num_processors> ./test -options_sml`

More info on fftw 2.1.5 here: https://www.fftw.org/download.html
More info on PETSc3 here: https://petsc.org/release/overview/

# 3D MPI-CTF Custom Program Usage

A custom program, `custom.cpp`, is provided in `CurveLab-2.1.3/fdct3d_mpi/src`, which is built off of `test.cpp`. To use the program: 

1. Configure the option file `options_cstm`:

   a. Set `-m,-n,-p` to input volume dimensions
   
   b. Set `-b` to be desired block size (block size 1 is highly recommended for multi-scale decomposition)
   
   c. Set `-nbscales` to be desired number of scales
   
   d. Set `-nbdstz_coarse` to be number of discretizations starting from 2nd coarest scale
   
3. In `custom.cpp`, on line 23, replace "data3d.bin" with your file (either txt or bin), then comment out lines corresponding to input
4. In `custom.cpp`, set `s_select` variable to be desired scale from 0 to (nscales - 2) for multi-scale. Otherwise, set `s_select = -1`
5. Compile `custom.cpp` with `make custom`
6. Run `mpirun -np <num_processors> ./custom -options_cstm`.

In running `custom`, each process will return an output file `outdata_3d_mpi_<proc>.raw` which can be concatenated into a single raw file by running the included scripts `concat_raw.sh`. `concat_raw.sh` takes in the number of files `n_proc` and returns a single output file `outdata_3d_mpi.raw`. Run `chmod +x concat_raw.sh`, then run `.concat_raw.sh <n_proc>`. 

Likewise, excess files can be removed with `rm_raw.sh` which will remove any excess raw files produced by `custom`and can be set up and executed in in the same way as `concat_raw.sh`
