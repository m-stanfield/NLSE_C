

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include "FFT_complex.h"

using std::complex;
using std::size_t;
using std::uintmax_t;
using std::vector;


// Private function prototypes
static size_t reverseBits(size_t val, int width);


void Fft::transform(vector<complex<double> > &vec, bool inverse) {
	size_t n = vec.size();
	if (n == 0)
		return;
	else if ((n & (n - 1)) == 0)  // Is power of 2
		transformRadix2(vec, inverse);
	else  // More complicated algorithm for arbitrary sizes
		transformBluestein(vec, inverse);
}


void Fft::transformRadix2(vector<complex<double> > &vec, bool inverse) {
	// Length variables
	size_t n = vec.size();
	int levels = 0;  // Compute levels = floor(log2(n))
	for (size_t temp = n; temp > 1U; temp >>= 1)
		levels++;
	if (static_cast<size_t>(1U) << levels != n)
		throw std::domain_error("Length is not a power of 2");

	// Trigonometric table
	vector<complex<double> > expTable(n / 2);
	for (size_t i = 0; i < n / 2; i++)
		expTable[i] = std::polar(1.0, (inverse ? 2 : -2) * M_PI * i / n);

	// Bit-reversed addressing permutation
	for (size_t i = 0; i < n; i++) {
		size_t j = reverseBits(i, levels);
		if (j > i)
			std::swap(vec[i], vec[j]);
	}

	// Cooley-Tukey decimation-in-time radix-2 FFT
	for (size_t size = 2; size <= n; size *= 2) {
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (size_t i = 0; i < n; i += size) {
			for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				complex<double> temp = vec[j + halfsize] * expTable[k];
				vec[j + halfsize] = vec[j] - temp;
				vec[j] += temp;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
}


void Fft::transformBluestein(vector<complex<double> > &vec, bool inverse) {
	// Find a power-of-2 convolution length m such that m >= n * 2 + 1
	size_t n = vec.size();
	size_t m = 1;
	while (m / 2 <= n) {
		if (m > SIZE_MAX / 2)
			throw std::length_error("Vector too large");
		m *= 2;
	}

	// Trigonometric table
	vector<complex<double> > expTable(n);
	for (size_t i = 0; i < n; i++) {
		uintmax_t temp = static_cast<uintmax_t>(i) * i;
		temp %= static_cast<uintmax_t>(n) * 2;
		double angle = (inverse ? M_PI : -M_PI) * temp / n;
		expTable[i] = std::polar(1.0, angle);
	}

	// Temporary vectors and preprocessing
	vector<complex<double> > avec(m);
	for (size_t i = 0; i < n; i++)
		avec[i] = vec[i] * expTable[i];
	vector<complex<double> > bvec(m);
	bvec[0] = expTable[0];
	for (size_t i = 1; i < n; i++)
		bvec[i] = bvec[m - i] = std::conj(expTable[i]);

	// Convolution
	vector<complex<double> > cvec = convolve(std::move(avec), std::move(bvec));

	// Postprocessing
	for (size_t i = 0; i < n; i++)
		vec[i] = cvec[i] * expTable[i];
}


vector<complex<double> > Fft::convolve(
		vector<complex<double> > xvec,
		vector<complex<double> > yvec) {

	size_t n = xvec.size();
	if (n != yvec.size())
		throw std::domain_error("Mismatched lengths");
	transform(xvec, false);
	transform(yvec, false);
	for (size_t i = 0; i < n; i++)
		xvec[i] *= yvec[i];
	transform(xvec, true);
	for (size_t i = 0; i < n; i++)  // Scaling (because this FFT implementation omits it)
		xvec[i] /= static_cast<double>(n);
	return xvec;
}


static size_t reverseBits(size_t val, int width) {
	size_t result = 0;
	for (int i = 0; i < width; i++, val >>= 1)
		result = (result << 1) | (val & 1U);
	return result;
}
