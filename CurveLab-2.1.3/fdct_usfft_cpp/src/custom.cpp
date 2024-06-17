#include <iostream>  
#include <vector>
#include <fstream>  
#include <math.h> 
#include "fdct_usfft.hpp"
#include "fdct_usfft_inline.hpp"


int main(){
    int n = 256; //64, 256,
    int nbscales = floor(log2(n))-3;
    int nbangles_coarse = 16; 
    int allcurvelets = 0; 
    fdct_usfft_ns::CpxNumMat x = fdct_usfft_ns::CpxNumMat(n,n); //mex2cpp(prhs[6], x);
    fdct_usfft_ns::CpxNumMat zeros = fdct_usfft_ns::CpxNumMat(n,n);
    vector< vector<fdct_usfft_ns::CpxNumMat> > c;  //vector<int> extra;
    vector< vector<fdct_usfft_ns::CpxNumMat> > c0;  //vector<int> extra;

    /*std::ifstream file;
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
        
    file.close();*/

    int cnt = 0;
    
    for(int j = 0; j < n; j++){
        for(int i = 0; i < n; i++){
            //x(i,j) = fdct_wrapping_ns::cpx(v[cnt],0);
            zeros(i,j) = fdct_usfft_ns::cpx(0,0);
            cnt++;
        }
    }
    

    //fdct_wrapping(n, n, nbscales, nbangles_coarse, allcurvelets, x, c);
    fdct_usfft(n, n, nbscales, nbangles_coarse, allcurvelets, zeros, c0);
    
    //std::cout << c.size() << std::endl;
    //std::cout << c[3].size() << std::endl;
    
    fdct_usfft_ns::cpx cpx1(1,0);
    std::cout << cpx1 << std::endl;
    c0[5][0](65,22) = cpx1;
    //std::cout << c0[5][0] << std::endl;
    //c0[2][1535]= c[2][1535];
    
    fdct_usfft_ns::CpxNumMat x_2, x_3;
    ifdct_usfft(n, n, nbscales, nbangles_coarse, allcurvelets, c0, x_2);
    afdct_usfft(n, n, nbscales, nbangles_coarse, allcurvelets, c0, x_3);

    std::cout << x_2(259,254) << std::endl;
    std::cout << x_3(259,254) << std::endl;
    

    std::ofstream myfile("outdata_2d_basic.txt");
    if (myfile.is_open()){
       
	    for(int j=0; j<n; j++){
		    for(int i=0; i<n; i++) {
		        myfile << std::real(x_2(i,j)) << "\n" ;
                myfile << std::imag(x_2(i,j)) << "\n" ;
		    }
        } 
        myfile.close();
    }
    else std::cout << "Unable to open file";
    
    return 0;
}