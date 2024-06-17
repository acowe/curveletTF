/* COPIED OVER FROM test.cpp by Lexing Ying
*/

#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

#include <iostream>  
#include <vector>
#include <fstream>  

int optionsCreate(const char* optfile, map<string,string>& options)
{
  options.clear();
  ifstream fin(optfile); assert(fin.good());
  string name;  fin>>name;
  while(fin.good()) {
	 char cont[100];	 fin.getline(cont, 99);
	 options[name] = string(cont);
	 fin>>name;
  }
  fin.close();
  return 0;
}

int main(int argc, char** argv)
{
  time_t tm0, tm1;
  
  int m = 64;
  int n = 64;
  int p = 64;
  int nbscales = 4;//floor(log2(n))-2;
  int nbdstz_coarse = 8;

  std::ifstream file;
    file.open("data3d.txt");
    if(!file)
    {
        std::cout << "There was a problem opening file data3d.txt for reading";
        exit (0);
    }

    std::vector <double> v;
    double d;    
    
    while (file >> d)
        v.push_back(d);
        
    file.close();

  int cnt = 0;
  CpxNumTns x(m,n,p);
  for(int k=0; k<p; k++){
    for(int j=0; j<n; j++){
        for(int i=0; i<m; i++){
            x(i,j,k) = cpx(v[cnt], 0);
            cnt++;
        }
    }			  
  }

  tm0 = time(NULL);
  
  vector< vector<double> > fxs,fys,fzs;
  vector< vector<int> > nxs,nys,nzs;
  fdct3d_param(m, n, p, nbscales, nbdstz_coarse, fxs,fys,fzs,nxs,nys,nzs);
  tm1 = time(NULL);  cout<<"fdct3d_param "<<difftime(tm1,tm0)<<" seconds"<<endl;  tm0 = tm1;
  
  //1. fdct3d
  CpxCrvletOcr c("tmpc");
  CpxNumTns w;
  fdct3d_forward(m, n, p, nbscales, nbdstz_coarse, x, c,w);
  tm1 = time(NULL);  cout<<"fdct3d_forward "<<difftime(tm1,tm0)<<" seconds"<<endl;  tm0 = tm1;
  
  //2. ifdct3d
  //fdct3d_inverse(m, n, p, nbscales, nbdstz_coarse, c,w, x);
  //tm1 = time(NULL);  cout<<"fdct3d_inverse "<<difftime(tm1,tm0)<<" seconds"<<endl;  tm0 = tm1;
  
  CpxNumTns newx(x); clear(newx);
  fdct3d_inverse(m, n, p, nbscales, nbdstz_coarse, c,w, newx);
  tm1 = time(NULL);  cout<<"fdct3d_inverse "<<difftime(tm1,tm0)<<" seconds"<<endl;  tm0 = tm1;  //cerr<<energy(newx)<<endl;
  
  double mv = 0.0;
  for(int i=0; i<m; i++)	 for(int j=0; j<n; j++)		for(int k=0; k<p; k++)		  mv = max(mv, abs(newx(i,j,k)-x(i,j,k)));
  cerr<<"max error "<<mv<<endl;

  ofstream myfile("outdata_3d_outcore.txt");
    if (myfile.is_open()){
        for(int k=0; k<n; k++){
	        for(int j=0; j<n; j++){
		        for(int i=0; i<n; i++) {
		            myfile << real(x(i,j,k)) << "\n" ;
                    //myfile << imag(x_2(i,j,k)) << "\n" ;
		        }
            } 
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";


  return 0;
}
