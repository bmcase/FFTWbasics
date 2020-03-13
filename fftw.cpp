#include <fftw3.h>
#include <iostream>
// see documentation tutorial at
// http://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs 
// compile with g++ -g -o2 fftw.cpp -o foo -lfftw3 -lm
int main() {
	fftw_complex *in, *out;
	fftw_plan p_forward, p_backward;
	//...
	
	int N; 
	N = 64;  //fft will do evaluation at the 64th roots of unity 
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	//in is the coeffs ordered from constant term up to highest degree
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	std::cout << "sizeof(fftw_complex) " << sizeof(fftw_complex) << "\n"; 
	std::cout << "sizeof(*in) " << sizeof(*in) << "\n";  //sizeof gives the number of bits
	std::cout << "sizeof(*out) " << sizeof(*out) << "\n"; 
	

	p_forward = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);
	p_backward = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD , FFTW_MEASURE);
	//...put something in in
	int deg = 8; 
	//double acoef [deg] = { 1,2,3,4,5,6,7,8 }; 
	double acoef [deg] = { 1,1,0,0,0,0,0,0 };  // 1 + 1 x + ... + 2 x^7
	double bcoef [deg] = { 2,2,0,0,0,0,0,0 };  
	// product is (1 + x )(2 + 2x) = 2 + 4x + 2x^2 which is what the output shows up to small truncation errors
	
	//fft a
	for(int i = 0; i < N; i++){
		if (i < deg){
			in[i][0] = acoef[i];  //set to the coefs of a
			in[i][1] = 0; //set imaginary part to zero
		}else {
			in[i][0] = 0; // set higher powers to have zero coef
			in[i][1] = 0; //set imaginary part to zero
		}
	}
	
	fftw_execute(p_forward); 
	
	double afft[N][2]; 
	
	std::cout << "sizeof(acoef)" << deg<< "\n"; 
	for(int i = 0; i < N; i++){
		afft[i][0] = out[i][0]; 
		afft[i][1] = out[i][1];  
		std::cout << out[i][0] << " + ";
		std::cout << out[i][1] << " i \n"; 
	}
	
	//fftb
	for(int i = 0; i < N; i++){
		if (i < deg){
			in[i][0] = bcoef[i];  //set to the coefs of a
			in[i][1] = 0; //set imaginary part to zero
		}else {
			in[i][0] = 0; // set higher powers to have zero coef
			in[i][1] = 0; //set imaginary part to zero
		}
	}
	fftw_execute(p_forward); 
	
	double bfft[N][2]; 
	std::cout << "bfft" << "\n"; 
	for(int i = 0; i < N; i++){
		bfft[i][0] = out[i][0];
		bfft[i][1] = out[i][1];  
		std::cout << out[i][0] << " + ";
		std::cout << out[i][1] << " i \n"; 
	}
	
	//multiply ffta and fftb
	for(int i = 0; i < N; i++){  //(a+bi)(c+di) = (ac - bd) + (bc + ad)i
		in[i][0]  = afft[i][0] * bfft[i][0] - afft[i][1] * bfft[i][1]; 
		in[i][1] =  afft[i][1] * bfft[i][0] + afft[i][0] * bfft[i][1]; 
	}
	fftw_execute(p_backward); 
	std::cout << "a * b poly \n" ; 
	for(int i = 0; i < N; i++){
		std::cout << out[i][0] / N<< " + " << out[i][1] / N<<" i\n";   // Because 
		//FFTW computes an unnormalized DFT. Thus, computing a forward followed
		// by a backward transform (or vice versa) results in the original array scaled by N
		}
	//...
	
	fftw_destroy_plan(p_forward);
	fftw_free(in); fftw_free(out);
	
}