#ifndef PULSE_H
#define PULSE_H

#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <math.h>
#include <algorithm>
#include<fstream>
#include "FFT_complex.h"

using std::complex;
using std::vector;
using std::polar;
using std::cout;
using std::endl;
using std::string;

class Pulse
{
    public:
        Pulse();
        Pulse(Pulse&);
        virtual ~Pulse();

        void set_freq(double freq_range, double central_freq, int bins);
        void set_gamma(double gamma);
        void set_gaussian_spectrum(double fwhm, double central_freq, double energy);
        void set_energy(double energy);
        void set_taylor_phase(vector<double> coef);
        void add_taylor_phase(vector<double> coef);
        void set_spectral_fields(vector<complex<double>> field);
        void set_temporal_fields(vector<complex<double>> field);

        void save(string fileName);

        int get_bins(void);
        vector<double> get_freq(void);
        vector<double> get_time(void);

        vector<complex<double>> get_spectral_fields(void);
        vector<double> get_spectral_intensity(void);
        vector<double> get_spectral_phase(void);

        vector<complex<double>> get_temporal_fields(void);
        vector<double> get_temporal_intensity(void);
        vector<double> get_temporal_phase(void);

        double calc_energy(void);
        double calc_peak_power(void);




        double get_energy(void);
        double get_central_freq(void);

    protected:

    private:

        vector<complex<double>> m_spectral_fields;
        vector<complex<double>> m_temporal_fields;

        vector<double> m_freq;
        vector<double> m_time;

        double m_energy;
        double m_central_freq;
        int m_bins;


        vector<double> m_abs_square(vector<complex<double>> x);
        vector<double> m_get_phase(vector<complex<double>> x);

        vector<complex<double>> m_fft(vector<complex<double>> x);
        vector<complex<double>> m_ifft(vector<complex<double>> x);
};

#endif // PULSE_H
