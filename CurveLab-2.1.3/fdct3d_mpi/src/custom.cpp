/* Using Lexing Ying's test.cpp FDCT3D_mpi (Fast 3d Curvelet Transform) as base
*/
#include <iostream>
#include <vector>
#include <fstream>  
#include <math.h> 
#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

int main(int argc, char** argv)
{
  PetscInitialize(&argc,&argv,"options",NULL);  //PetscTruth flg = PETSC_FALSE;
  srand48(0);
  
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
 
  iC( PetscPrintf(MPI_COMM_WORLD, "mpisize %d\n", mpisize) );  
  iC( PetscPrintf(MPI_COMM_SELF, "mpirank %d\n", mpirank) );

  PetscBool flg = PETSC_FALSE;
  int m; iC( PetscOptionsGetInt(NULL, NULL, "-m", &m, &flg) ); iA(flg==PETSC_TRUE);
  int n;  iC( PetscOptionsGetInt(NULL, NULL, "-n", &n, &flg) ); iA(flg==PETSC_TRUE);
  int p;  iC( PetscOptionsGetInt(NULL, NULL, "-p", &p, &flg) ); iA(flg==PETSC_TRUE);
  int b; iC( PetscOptionsGetInt(NULL, NULL, "-b", &b, &flg) ); iA(flg==PETSC_TRUE);
  int nbscales;  iC( PetscOptionsGetInt(NULL,NULL, "-nbscales", &nbscales, &flg) ); iA(flg==PETSC_TRUE);
  int nbdstz_coarse;  iC( PetscOptionsGetInt(NULL,NULL, "-nbdstz_coarse", &nbdstz_coarse, &flg) ); iA(flg==PETSC_TRUE);
  
  CpxNumTnsBlkd X, zeros;
  BolNumTns newtnsexists(m/b,n/b,p/b);
  IntNumTns newtnsowners(m/b,n/b,p/b);
  iC( fdct3d_partition_cpxnumtnsblkd_z(m,n,p,b, newtnsexists, newtnsowners) );
  X.setup(m,n,p,b, newtnsowners);
  zeros.setup(m,n,p,b, newtnsowners);
  
  //1. generate data
  int e = X.e();  int f = X.f();  int g = X.g();
  std::ifstream file;
  file.open("data3d.txt");
  if(!file)
  {
    std::cout << "There was a problem opening file data3d.txt for reading";
        exit (0);
  }
  std::vector <double> v,v_out;
  double d;    
  while (file >> d) {
    v.push_back(d);
    v_out.push_back(d);
  }
  file.close();

  iC( PetscPrintf(MPI_COMM_WORLD, "b is %i\n", d) )
  iC( PetscPrintf(MPI_COMM_WORLD, "num data_entries: %i\n", v.size()) );
  
  int cnt = 0;
  for(int k=0; k<g; k++)	 for(int j=0; j<f; j++)		for(int i=0; i<e; i++) {
	  if(X.owners()(i,j,k)==mpirank && zeros.owners()(i,j,k)==mpirank) {
      //cerr << X.owners()(i,j,k) << "," <<  cnt << endl;
		  CpxNumTns& Xblk = X.block(i,j,k);
      CpxNumTns& zerosBlk = zeros.block(i,j,k);
		  for(int koff=0; koff<b; koff++)		  for(int joff=0; joff<b; joff++)			 for(int ioff=0; ioff<b; ioff++) {
		    int ind = (mpirank+1) * cnt;
        //cerr << mpirank << "," <<  cnt << "," << ind << endl;
        Xblk(ioff,joff,koff) = cpx(v[ind], 0); //cpx(drand48(), drand48());
        zerosBlk(ioff,joff,koff) = cpx(0, 0); //cpx(drand48(), drand48());
        cnt++;
		  }
	  }
  }
  
  
  double xene = X.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "energy %e\n", xene) );  //  if(mpirank==0)	 cerr<<"X energy  "<<xene<<endl;
  

  time_t tm0, tm1;  tm0 = time(NULL);
  //2. forward
  CpxCrvletPrtd C;  CpxNumTnsBlkd W; CpxCrvletPrtd C0;  CpxNumTnsBlkd W0;
  iC( fdct3d_forward(m,n,p, nbscales,nbdstz_coarse, X, C, W) );
  //iC( fdct3d_forward(m,n,p, nbscales,nbdstz_coarse, zeros, C0,W0) );
  tm1 = time(NULL);  iC( PetscPrintf(MPI_COMM_WORLD, "FORWARD %e\n", difftime(tm1,tm0)) );  tm0 = tm1;
  double cene = C.globalenergy();
  double wene = W.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "CW energy %e %e %e\n", cene, wene, cene+wene) );

  //C0[0] = C[0];

  
  //3. inverse
  CpxNumTnsBlkd Y;
  iC( fdct3d_inverse(m,n,p, nbscales,nbdstz_coarse, C,W,Y) );
  tm1 = time(NULL);  iC( PetscPrintf(MPI_COMM_WORLD, "INVERSE %e\n", difftime(tm1,tm0)) );  tm0 = tm1;
  double yene = Y.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "Y energy %e\n", yene) );
  
  //4. compute difference

  ofstream myfile("outdata_3d_mpi1.txt");
  cnt = 0;
  if (myfile.is_open()){
    for(int k=0; k<g; k++)	 for(int j=0; j<f; j++)		for(int i=0; i<e; i++) {
	    if(X.owners()(i,j,k)==mpirank) {
        //cerr << X.owners()(i,j,k) << "," <<  cnt << endl;
		    iA(Y.owners()(i,j,k)==mpirank); 		//CHECK THEY HAVE THE SAME DISTRIBUTION
		    CpxNumTns& Xblk = X.block(i,j,k);	CpxNumTns& Yblk = Y.block(i,j,k);
		    for(int koff=0; koff<b; koff++)		  for(int joff=0; joff<b; joff++)			 for(int ioff=0; ioff<b; ioff++) {
          int ind = (mpirank+1) * cnt;
          //cerr << mpirank << "," << real(Yblk(ioff,joff,koff))-v[ind]<< endl;
          v_out[ind] = real(Yblk(ioff,joff,koff));
          //iC( PetscPrintf(MPI_COMM_WORLD, "%e\n", real(Yblk(i,j,k))) );
          cnt++;
		    }
	    }
    }
    myfile.close();
  }
  else std::cout << "Unable to open file";
  //cerr<< cnt<< endl;
  double eene = Y.globalenergy();
  iC( PetscPrintf(MPI_COMM_WORLD, "E energy %e\n", eene) );
  
  
  /*if(mpisize==1) {
	 for(int s=0; s<nbscales-1; s++) {
		int nw = C.owners()[s].size();
		for(int w=0; w<nw/2; w++) {
		  CpxNumTns& A = C.block(s,w);		  CpxNumTns& B = C.block(s,w+nw/2);
		  double maxerr = 0;
		  for(int i=0; i<A.m(); i++)			 for(int j=0; j<A.n(); j++)				for(int k=0; k<A.p(); k++)
			 maxerr = max(maxerr, abs(A(i,j,k)-conj(B(i,j,k))));
		  cerr<<s<<" "<<w<<" "<<maxerr<<endl;
		}
	 }
  }*/
  

  PetscFinalize();
  return 0;
}
