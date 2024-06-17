#include <iostream>  
#include <vector>
#include <fstream>  
#include <math.h> 
#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

int main(){
    int n = 64; //64
    int nbscales = floor(log2(n))-2;
    int nbdstz_coarse = 8; 
    int allcurvelets = 0; 
    CpxNumTns x = CpxNumTns(n,n,n); //mex2cpp(prhs[6], x);
    CpxNumTns zeros = CpxNumTns(n,n,n);
    vector< vector<CpxNumTns> > c;  //vector<int> extra;
    vector< vector<CpxNumTns> > c0;  //vector<int> extra;

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
    for(int k = 0; k < n; k++){
        for(int j = 0; j < n; j++){
            for(int i = 0; i < n; i++){
                x(i,j,k) = cpx(v[cnt],0);
                zeros(i,j,k) = cpx(0,0);
                cnt++;
            }
        }
    }

    fdct3d_forward(n,n,n,nbscales, nbdstz_coarse, allcurvelets, x, c);
    fdct3d_forward(n,n,n,nbscales, nbdstz_coarse, allcurvelets, zeros, c0);
    
    std::cout << c.size() << std::endl;
    //std::cout << c[3].size() << std::endl;
    //cpx cpx1 = (1,0);
    //setvalue(c[2][1],cpx1);
    c0[0]= c[0];
    

    
    CpxNumTns x_2;
    fdct3d_inverse(n,n,n, nbscales, nbdstz_coarse, allcurvelets, c0, x_2);

    std::cout << x_2._data[0] << std::endl;

    ofstream myfile("outdata_3d_scale1_reals.txt");
    if (myfile.is_open()){
        for(int k=0; k<n; k++){
	        for(int j=0; j<n; j++){
		        for(int i=0; i<n; i++) {
		            myfile << real(x_2(i,j,k)) << "\n" ;
                    //myfile << imag(x_2(i,j,k)) << "\n" ;
		        }
            } 
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";
    
    return 0;
}
