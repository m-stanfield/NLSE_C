#ifndef FFT_COMPLEX_H_INCLUDED
#define FFT_COMPLEX_H_INCLUDED

#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <math.h>
#include <algorithm>
#include<fstream>

using std::vector;
using std::polar;
using std::cout;
using std::endl;
using std::string;

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

using std::complex;
using std::size_t;
using std::uintmax_t;


// Private function prototypes
static size_t reverseBits(size_t val, int width);

namespace Fft {

	/*
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function. The inverse transform does not perform scaling, so it is not a true inverse.
	 */
	void transform(std::vector<std::complex<double> > &vec, bool inverse);


	/*
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	 */
	void transformRadix2(std::vector<std::complex<double> > &vec, bool inverse);


	/*
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */
	void transformBluestein(std::vector<std::complex<double> > &vec, bool inverse);


	/*
	 * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
	 */
	std::vector<std::complex<double> > convolve(
		std::vector<std::complex<double> > xvec,
		std::vector<std::complex<double> > yvec);

}

#endif // FFT_COMPLEX_H_INCLUDED
