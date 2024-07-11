/* Using Lexing Ying's test.cpp FDCT3D_mpi (Fast 3d Curvelet Transform) as base*/
#include <iostream>
#include <vector>
#include <fstream>  
#include <math.h> 
#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

int main(int argc, char** argv)
{
  // Load data from file (lines 12-38)

  // uncomment lines 13-19, comment 23-36 to read in data from text file
  /*std::ifstream file;
  file.open("data3d.txt");
  double d;    
  while (file >> d) {
    v.push_back(d);
    v_out.push_back(0.0f);
  }*/
  
  // comment lines 13-19, uncomment 23-36 to read in data from binary file (default)
  std::ifstream file("data3d.bin", std::ios::binary);

  if(!file)
  {
    std::cout << "There was a problem opening file data3d.txt for reading";
        exit (0);
  }

  file.seekg(0, std::ios::end);
  std::streampos fileSize = file.tellg();
  file.seekg(0, std::ios::beg);
  size_t numFloats = fileSize / sizeof(float);
  std::vector <float> v(numFloats);
  file.read(reinterpret_cast<char*>(v.data()), fileSize);

  file.close(); // don't comment

  // MPI Init => Begin parallelization
  PetscInitialize(&argc,&argv,"options",NULL);  //PetscTruth flg = PETSC_FALSE;
  srand48(0);
  
  PetscInt num_entries = v.size();
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  PetscInt num_entries_proc = num_entries/mpisize;

  std::vector <float> v_out(num_entries_proc,0);
 
  iC( PetscPrintf(MPI_COMM_WORLD, "mpisize %d\n", mpisize) );  
  iC( PetscPrintf(MPI_COMM_WORLD, "entries per proc %d\n", num_entries_proc) );
  iC( PetscPrintf(MPI_COMM_SELF, "mpirank %d\n", mpirank) );
  

  // global vars: x,y,z dimensions of input data, block (dim) size b, number of scales
  PetscBool flg = PETSC_FALSE;
  int m; iC( PetscOptionsGetInt(NULL, NULL, "-m", &m, &flg) ); iA(flg==PETSC_TRUE);
  int n;  iC( PetscOptionsGetInt(NULL, NULL, "-n", &n, &flg) ); iA(flg==PETSC_TRUE);
  int p;  iC( PetscOptionsGetInt(NULL, NULL, "-p", &p, &flg) ); iA(flg==PETSC_TRUE);
  int b; iC( PetscOptionsGetInt(NULL, NULL, "-b", &b, &flg) ); iA(flg==PETSC_TRUE);
  int nbscales;  iC( PetscOptionsGetInt(NULL,NULL, "-nbscales", &nbscales, &flg) ); iA(flg==PETSC_TRUE);
  int nbdstz_coarse;  iC( PetscOptionsGetInt(NULL,NULL, "-nbdstz_coarse", &nbdstz_coarse, &flg) ); iA(flg==PETSC_TRUE);
  
  // CpxNumTnsBlkd is global input data 
  CpxNumTnsBlkd X, zeros;
  // NumTns keeps track of processor specific info for each block of blocked CTF
  BolNumTns newtnsexists(m/b,n/b,p/b);
  IntNumTns newtnsowners(m/b,n/b,p/b);
  iC( fdct3d_partition_cpxnumtnsblkd_z(m,n,p,b, newtnsexists, newtnsowners) );
  X.setup(m,n,p,b, newtnsowners);
  zeros.setup(m,n,p,b, newtnsowners);
  
  //1. Load data into X
  iC( PetscPrintf(MPI_COMM_WORLD, "b is %i\n", b) )
  iC( PetscPrintf(MPI_COMM_WORLD, "num data_entries: %i\n", v.size()) );
  int cnt = 0;
  int e = X.e();  int f = X.f();  int g = X.g();
  //cerr << e << "," << f << "," << g << endl;
  for(int k=0; k<g; k++)	 for(int j=0; j<f; j++)		for(int i=0; i<e; i++) { // for each block in X
  //for(int i=0; i<e; i++)	 for(int j=0; j<f; j++)		for(int k=0; k<g; k++) {
	  if(X.owners()(i,j,k)==mpirank && zeros.owners()(i,j,k)==mpirank) {
		  CpxNumTns& Xblk = X.block(i,j,k);
      CpxNumTns& zerosBlk = zeros.block(i,j,k);
		  for(int koff=0; koff<b; koff++)		  for(int joff=0; joff<b; joff++)			 for(int ioff=0; ioff<b; ioff++) {
      //for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++)			 for(int koff=0; koff<b; koff++) {
		    int ind = (num_entries_proc * mpirank) + cnt;
        Xblk(ioff,joff,koff) = cpx(v[ind], 0.0d);
        zerosBlk(ioff,joff,koff) = cpx(0, 0); 
        cnt++;
		  }
	  }
  }
  
  
  double xene = X.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "energy %e\n", xene) );  //  if(mpirank==0)	 cerr<<"X energy  "<<xene<<endl;
  

  time_t tm0, tm1;  tm0 = time(NULL);
  //2. forward --> C is curvelet coeff data structure, W is meant to hold wrapped data (the data mapped to the rectangle centered around the origin)
  CpxCrvletPrtd C;  CpxNumTnsBlkd W;
  iC( fdct3d_forward(m,n,p, nbscales,nbdstz_coarse, X, C, W) );
  tm1 = time(NULL);  iC( PetscPrintf(MPI_COMM_WORLD, "FORWARD %e\n", difftime(tm1,tm0)) );  tm0 = tm1;
  double cene = C.globalenergy();
  double wene = W.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "CW energy %e %e %e\n", cene, wene, cene+wene) );

  //2.5 multi-scale decomposition
  int s_select = -1; // -1 for all scales, 0 <= s_select < nscales - 2 otherwise
  vector< vector<int> >& c = C.nx();
  for (int s=0; s<c.size(); s++){
    if (s_select != -1 && s != s_select) for(int w=0; w<c[s].size(); w++){
      if(C.owners()[s][w]==mpirank){
          CpxNumTns& Cblk = C.block(s,w);
          for(int k1=0; k1<Cblk.m(); k1++)		  for(int k2=0; k2<Cblk.n(); k2++)			 for(int k3=0; k3<Cblk.p(); k3++) Cblk(k1,k2,k3) = cpx(0.0f,0.0f);
      } 
    }
  }
  
  //3. inverses
  CpxNumTnsBlkd Y;
  iC( fdct3d_inverse(m,n,p, nbscales,nbdstz_coarse, C,W,Y) );
  tm1 = time(NULL);  iC( PetscPrintf(MPI_COMM_WORLD, "INVERSE %e\n", difftime(tm1,tm0)) );  tm0 = tm1;
  double yene = Y.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "Y energy %e\n", yene) );
  
  // 3.5 Load output data into to v_out
  cnt = 0;
  for(int k=0; k<g; k++)	 for(int j=0; j<f; j++)		for(int i=0; i<e; i++) {
  //for(int i=0; i<e; i++)	 for(int j=0; j<f; j++)		for(int k=0; k<g; k++) { 
	  if(X.owners()(i,j,k)==mpirank) {
		iA(Y.owners()(i,j,k)==mpirank); 		//CHECK THEY HAVE THE SAME DISTRIBUTION
		CpxNumTns& Xblk = X.block(i,j,k);		CpxNumTns& Yblk = Y.block(i,j,k);  
		for(int koff=0; koff<b; koff++)		  for(int joff=0; joff<b; joff++)			 for(int ioff=0; ioff<b; ioff++) {
    //for(int ioff=0; ioff<b; ioff++)		  for(int joff=0; joff<b; joff++)			 for(int koff=0; koff<b; koff++) {
      v_out[cnt] = real(Yblk(ioff,joff,koff));
      cnt++;
      Yblk(ioff,joff,koff) -= Xblk(ioff,joff,koff);
      
		}
	 }
  }

  // compute global error between input X and non-decomposed output Y
  double eene = Y.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "E energy %e\n", eene) ); 

  
  //4. Write data from v_out to output file (lines 152-167)
  // uncomment lines 152-158, comment lines 161-166 to load out to text file
  /*ofstream myfile("outdata_3d_mpi_" + std::to_string(mpirank) + ".txt");
  if (myfile.is_open()){
    for (int c = 0; c < num_entries_proc; c++){
      myfile << v_out[c] << "\n" ;
    }
    myfile.close();
  }*/

  // comment lines 152-158, uncomment lines 161-166 to load out to raw binary file (default)
  ofstream myfile("outdata_3d_mpi_" + std::to_string(mpirank) + ".raw");
  if (!myfile) {
        std::cerr << "Error opening file for writing." << std::endl;
        return 1;
  }
  myfile.write(reinterpret_cast<const char*>(v_out.data()), v_out.size() * sizeof(float));
  
  myfile.close(); // don't comment out
  
  PetscFinalize();
  return 0;
}
