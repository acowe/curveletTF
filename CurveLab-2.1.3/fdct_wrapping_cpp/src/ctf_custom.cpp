#include <iostream>  
#include <vector>
#include <fstream>  
#include <math.h> 
#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"


int main(){
    int n = 256; //64, 256,
    int nbscales = floor(log2(n))-3;
    int nbangles_coarse = 16; 
    int allcurvelets = 0; 
    fdct_wrapping_ns::CpxNumMat x = fdct_wrapping_ns::CpxNumMat(n,n); //mex2cpp(prhs[6], x);
    fdct_wrapping_ns::CpxNumMat zeros = fdct_wrapping_ns::CpxNumMat(n,n);
    vector< vector<fdct_wrapping_ns::CpxNumMat> > c;  //vector<int> extra;
    vector< vector<fdct_wrapping_ns::CpxNumMat> > c0;  //vector<int> extra;

    std::ifstream file;
    file.open("data2d.txt");
    if(!file)
    {
        std::cout << "There was a problem opening file data2d.txt for reading";
        exit (0);
    }

    std::vector <double> v;
    double d;    
    
    while (file >> d)
        v.push_back(d);
        
    file.close();

    int cnt = 0;
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            x(i,j) = fdct_wrapping_ns::cpx(v[cnt],0);
            zeros(i,j) = fdct_wrapping_ns::cpx(0,0);
            cnt++;
        }
    }
    

    fdct_wrapping(n, n, nbscales, nbangles_coarse, allcurvelets, x, c);
    //fdct_wrapping(n, n, nbscales, nbangles_coarse, allcurvelets, zeros, c0);
    
    //std::cout << c.size() << std::endl;
    std::cout << c[3][0].m() << std::endl;
    std::cout << c[3][0].n() << std::endl;
    
    
    /*fdct_wrapping_ns::cpx cpx1(1,0);
    std::cout << cpx1 << std::endl;
    c0[5][0](65,22) = cpx1;
    std::cout << c0[5][0].m() << std::endl;
    std::cout << c0[5][0].n() << std::endl;
    //c0[2][1535]= c[2][1535];*/
    
    fdct_wrapping_ns::CpxNumMat x_2;
    ifdct_wrapping(n, n, nbscales, nbangles_coarse, allcurvelets, c, x_2);

    std::cout << x_2.m() << std::endl;
    std::cout << x_2.n() << std::endl;
    //std::cout << x_2(511,511) << std::endl;

    std::ofstream myfile("outdata_2d_allreals.txt");
    if (myfile.is_open()){
	    for(int i=0; i<n; i++){
		    for(int j=0; j<n; j++) {
		        myfile << std::real(x_2(i,j)) << "\n" ;
                //myfile << std::imag(x_2(i,j)) << "\n" ;
		    }
        } 
        myfile.close();
    }
    else std::cout << "Unable to open file";
    
    return 0;
}