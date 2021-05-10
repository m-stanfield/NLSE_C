#include "Pulse.h"
#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <math.h>
#include <algorithm>
#include<fstream>
#include "FFT_complex.h"
#include "Pulse.h"

using std::complex;
using std::vector;
using std::polar;
using std::cout;
using std::endl;
using std::string;

Pulse::Pulse()
{
    //ctor
}

Pulse::Pulse(Pulse& pulse)
{
    m_spectral_fields = pulse.get_spectral_fields();
    m_temporal_fields = pulse.get_temporal_fields();

    m_freq = pulse.get_freq();
    m_time = pulse.get_time();

    m_energy = 1.0*pulse.get_energy();
    m_central_freq = 1.0*pulse.get_central_freq();
    m_bins = 1.0*pulse.get_bins();


    //ctor
}

Pulse::~Pulse()
{
    //dtor
}


void Pulse::save(string fileName) {

    vector<double> intensity = get_spectral_intensity();
    vector<double> phase = get_spectral_phase();
	std::ofstream file;
	file.open(fileName);
    file<<"Freq" << "\t" << "Int " <<"\t" << "Phase" <<endl;

	for(int i=0;i<m_freq.size();++i){
		file<< m_freq[i] << "\t" << intensity[i] << "\t" << phase[i] <<endl;
	}
	file.close();


}

void Pulse::set_spectral_fields(vector<complex<double>> field) {
	//Sets the spectral fields of the pulse to a pre-defined arbitrary vector
	//Contains both phase and amplitude infomation



    if (field.size() == m_freq.size())
    {
        if (m_spectral_fields.size() != m_freq.size())
        {

            for (int i = 0; i < m_freq.size(); i++) {
                m_spectral_fields.push_back(field[i]);
                m_temporal_fields.push_back(field[i]);
            }
        } else {
            for (int i = 0; i < m_freq.size(); i++) {
                m_spectral_fields[i] = 1.0*field[i];
                m_temporal_fields[i] = 1.0*field[i];

            }
        }


        cout << m_spectral_fields[0] << endl;
        cout << m_temporal_fields[0] << endl;
        Fft::transform(m_temporal_fields,false);



    }


}

void Pulse::set_temporal_fields(vector<complex<double>> field) {
	//Sets the temporal fields of the pulse to a pre-defined arbitrary vector
	//Contains both phase and amplitude infomation

    if (field.size() == m_freq.size())
    {
        if (m_spectral_fields.size() != m_freq.size())
        {

            for (int i = 0; i < m_freq.size(); i++) {
                m_spectral_fields.push_back(field[i]);
                m_temporal_fields.push_back(field[i]);
            }
        } else {
            for (int i = 0; i < m_freq.size(); i++) {
                m_spectral_fields[i] = 1.0*field[i];
                m_temporal_fields[i] = 1.0*field[i];

            }
        }



        Fft::transform(m_spectral_fields,true);



    }

}


void Pulse::set_freq(double freq_range, double central_freq, int bins) {
    // sets frequency array and calculates the corresponding temporal axis
    double freq_step = freq_range/bins;

    vector<double> freq;
    vector<double> time;


    for (int i = 0; i < bins; i++) {
       // cout << i << "    " << (i*freq_step - 0.5*central_freq) << endl;
        freq.push_back((i + 0.5)*freq_step - 0.5*freq_range + central_freq);
        time.push_back((i-0.5*bins)/freq_range);
    }
    m_freq = freq;
    m_central_freq = central_freq;
    m_bins = bins;

	m_time = time;
}



void Pulse::set_gaussian_spectrum(double fwhm, double central_freq, double energy) {
    //sets the intensity spectrum to be a gaussian located on the central frequency with a set fwhm
	//neesd to still set the energy of the pulse after defining the spectrum(needs implementaiton of fft first)

    vector<complex<double>> fields;
    for (auto i = m_freq.begin(); i != m_freq.end(); i++) {
        fields.push_back(sqrt(exp( -4*log(2)*pow(*i - central_freq, 2)/pow(fwhm,2))));
    }

    set_spectral_fields(fields);
    set_energy(energy);


//    m_spectral_fields.clear();
//    for (auto i = m_freq.begin(); i != m_freq.end(); i++) {
//        m_spectral_fields.push_back(sqrt(exp( -4*log(2)*pow(*i - central_freq, 2)/pow(fwhm,2))));
//    }
    //set_energy(energy);
}

void Pulse::set_energy(double energy) {
    //sets energy
	//currently not working. Requires implementaiton of FFT
    m_energy = energy;
    double sum = 0;

    vector<double> intensity = get_temporal_intensity();
    vector<double> phase = get_temporal_phase();
    vector<complex<double>> fields = get_temporal_fields();

    double time_step = (m_time[1] - m_time[0]); //converting normalized time from being ps(10e-12) to seconds
    for (int i = 0; i < m_time.size(); i++) {
        sum += intensity[i]*time_step;
    }

    for (int i = 0; i < m_time.size(); i++) {
        fields[i] = polar(sqrt(intensity[i]*(energy/sum)),phase[i]);
    }
    set_temporal_fields(fields);
}

double Pulse::calc_energy() {
    double sum = 0;

    vector<double> intensity = get_temporal_intensity();


    double time_step = (m_time[1] - m_time[0])*pow(10,-12); //converting normalized time from being ps(10e-12) to seconds
    for (int i = 0; i < m_time.size(); i++) {
        sum += intensity[i]*time_step;
    }

    return sum;
}

double Pulse::calc_peak_power() {


    vector<double> intensity = get_temporal_intensity();
    double peak_power = 0.0;

    for (int i = 0; i < m_time.size(); i++) {
        if (intensity[i] > peak_power) {
            peak_power = 1.0*intensity[i];
        }
    }

    return peak_power;
}


void Pulse::set_taylor_phase(vector<double> coef) {
	//sets the pahse to be defined by a Taylor Series. Can go to arbitrary high coef.
	//Starts on quadratic term.

    double ang_freq;
    long factorial;
    double phase;

    vector<complex<double>> temp_fields;
    for (int i = 0; i < m_freq.size(); i++) {
        ang_freq = 2*M_PI*(m_freq[i]-m_central_freq);
        phase = 0;
        factorial = 1;
        for (int j = 0; j < coef.size(); j++) {
            factorial = factorial*(j + 2);
            phase = phase + coef[j]/(factorial)*pow(ang_freq,j+2);
        }
        temp_fields.push_back(polar(abs(m_spectral_fields[i]),phase));
    }
    set_spectral_fields(temp_fields);

}


void Pulse::add_taylor_phase(vector<double> coef) {
    double ang_freq;
    long factorial;
    double phase;

    vector<double> initial_phase = get_spectral_phase();
    vector<complex<double>> temp_fields;
    for (int i = 0; i < m_freq.size(); i++) {
        ang_freq = 2*M_PI*(m_freq[i]-m_central_freq);
        phase = 0;
        factorial = 1;
        for (int j = 0; j < coef.size(); j++) {
            factorial = factorial*(j + 2);
            phase = phase + coef[j]/(factorial)*pow(ang_freq,j+2);
        }
        temp_fields.push_back(polar(abs(m_spectral_fields[i]),phase + initial_phase[i]));
    }
    set_spectral_fields(temp_fields);

}


vector<double> Pulse::get_freq(void) {
    //returns frequency array
    return m_freq;
}

vector<double> Pulse::get_time(void) {
    //returns frequency array
    return m_time;
}

vector<complex<double>> Pulse::get_spectral_fields(void) {
    //returns complex spectral fields
    return m_spectral_fields;
}

double Pulse::get_energy(void) {
    //returns energy value
    return m_energy;
}

double Pulse::get_central_freq(void){
    //returns central frequency
    return m_central_freq;
}

vector<double> Pulse::get_spectral_intensity(void) {
    return m_abs_square(m_spectral_fields);
}

int Pulse::get_bins(void) {
    return m_bins;
}

vector<double> Pulse::get_spectral_phase(void) {
    return m_get_phase(m_spectral_fields);

}


vector<complex<double>> Pulse::get_temporal_fields(void) {
//write script to take FFT and output tmeporal fields
    return m_temporal_fields;
}

vector<double> Pulse::get_temporal_intensity(void) {
    return m_abs_square(get_temporal_fields());
}

vector<double> Pulse::get_temporal_phase(void) {
    return m_get_phase(get_temporal_fields());
}

vector<double> Pulse::m_abs_square(vector<complex<double>> x) {
    vector<double> temp_vec;

    for (auto i = x.begin(); i != x.end(); ++i) {
        temp_vec.push_back(pow(abs(*i),2));
    }

    return temp_vec;
}

vector<double> Pulse::m_get_phase(vector<complex<double>> x) {
    vector<double> temp_vec;

    for (auto i = x.begin(); i != x.end(); ++i) {
        temp_vec.push_back(arg(*i));
    }
    return temp_vec;
}
