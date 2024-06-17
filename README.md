# Set-up Instructions

1. Compile fftw and PETSc3 rom their root directories: fftw-2.1.5 with `./configure --with-pic && make`; Petsc3 with `./configure && make all check`
2. In all makefile.opt files, (CurveLab-2.1.3/fdct3d_mpi/src/makefile.opt, CurveLab-2.1.3/makefile.opt) set `FFTW_DIR =  <curveletTF_directory>/fftw-2.1.5` (or wherever your fftw installation is)
3. In CurveLab-2.1.3/fdct3d_mpi/src/makefile.opt, set `PETSC_DIR =  <curveletTF_directory>/petsc` (or wherever your PETSc3 installation is)
4. In CurveLab-2.1.3 build and compile 2D and 3D (in-core, out-of-core) Curvelet TF implementations with `make lib`, then test CPP installations with `make test`
5. To compile MPI 3D CTF, go into CurveLab-2.1.3/fdct3d_mpi/src/makefile.opt, then run 'make libfdct3d.a' to compile the library, then test CPP installation with
`mpirun -np <num_processors> ./test -options_sml`

More info on fftw 2.1.5 here: https://www.fftw.org/download.html
More info on PETSc3 here: https://petsc.org/release/overview/
