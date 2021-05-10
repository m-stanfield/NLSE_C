#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <math.h>
#include <algorithm>
#include<fstream>
#include "FFT_complex.h"
#include "Pulse.h"
#include "Media.h"

using std::complex;
using std::vector;
using std::polar;
using std::arg;

using std::cout;
using std::endl;
using std::string;


//FUNCTION DEFINITIONS

//vector<double> abs_square(vector<complex<double>> x);


// CLASS DEFINITION



//CLASS FUNCTIONS



//TESTING FUNCTIONS
//
void test_phase_set_get(Pulse pulse) {
    cout << "\n\nTESTING PHASE SET AND GET\n\n" << endl;
    vector<double> freq = pulse.get_freq();
    vector<double> phase = pulse.get_spectral_phase();
    vector<double> intensity = pulse.get_spectral_intensity();

    cout << "Index\t\tFreq\t\tInt.\t\tAngle" << endl;

    for (int i = 0; i < phase.size(); i++) {
        cout << i << "\t\t" << freq[i] << "\t\t"<< intensity[i] << "\t\t" << phase[i] << endl;
    }
}
void test_freq_set_get(Pulse pulse) {

    cout << "\n\nTESTING FREQ\n\n" << endl;
    vector<double> freq = pulse.get_freq();

    cout << "List of Frequencies" << endl;
    for ( auto i = freq.begin(); i != freq.end(); ++i) {
        cout << *i << endl;
    }
    cout << "Total Size: " << freq.size() << endl;
}


void test_gaussian_set(Pulse pulse) {

    cout << "\n\nTESTING INTENSITY SETTING AND RETRIEVING\n\n" << endl;
    vector<complex<double>> fields = pulse.get_spectral_fields();

    vector<double> freq = pulse.get_freq();
    vector<double> intensity = pulse.get_spectral_intensity();
    cout << "Element Number\t\tFrequency\t\tValue" << endl;
    for (int i = 0; i < freq.size(); i++) {
        cout << i << "\t\t" << freq[i] << "\t\t" << intensity[i] << endl;
    }
}



void test_temporal_set(Pulse pulse) {

    cout << "\n\nTESTING Temporal INTENSITY SETTING AND RETRIEVING\n\n" << endl;
    vector<complex<double>> fields = pulse.get_temporal_fields();


    vector<double> freq = pulse.get_freq();
    vector<double> intensity = pulse.get_temporal_intensity();
    cout << "Element Number\t\tFrequency\t\tValue" << endl;
    for (int i = 0; i < freq.size(); i++) {
        cout << i << "\t\t" << freq[i] << "\t\t" << intensity[i] << endl;
    }


    cout << "\n\nTESTING Temporal PHASE SET AND GET\n\n" << endl;

    vector<double> phase = pulse.get_temporal_phase();
   // vector<double> intensity = pulse.get_temporal_intensity();

    cout << "Index\t\tFreq\t\tInt.\t\tAngle" << endl;

    for (int i = 0; i < phase.size(); i++) {
        cout << i << "\t\t" << freq[i] << "\t\t"<< intensity[i] << "\t\t" << phase[i] << endl;
    }
}


void test_fft(Pulse pulse) {
    vector<complex<double>> temporal_field = pulse.get_temporal_fields();
    vector<complex<double>> spectral_field = pulse.get_spectral_fields();
    vector<double> freq = pulse.get_freq();

    for (int i = 0; i < freq.size(); i++) {
        cout << temporal_field[i] << endl;
    }

    Fft::transform(temporal_field,true);

    for (int i = 0; i < freq.size(); i++) {
        cout << spectral_field[i]/abs(spectral_field[i]) << "\t" << temporal_field[i]/abs(temporal_field[i]) << endl;
    }

}

const double C_SOL = 2.99792458e8;

int main() {


    double freq_range = 2000.0*pow(10,12);         //Frequency Range in Hz
    double central_freq = 374.7*pow(10,12);        //Central Frequency of Spectrum in Hz
    double n2 = 2.5*pow(10,-20);         //Kerr index in m^2/W
    double beam_dia = 0.012;             //e^(-2) intensity beam diameter in meters.
    int bins = pow(2,12);               //Number of bins used in system(neesd power of 2)
    double length = 0.005;              //Length in meters
    int number_steps = 100;              //number of steps to split each media into
    double beta2 = 0*36.0;                //GVD of material, in fs^2/mm


    double energy = 0.001;              //Energy of system, in Joules


    double gamma = 2*M_PI*(central_freq/C_SOL)*(2/(M_PI*pow(0.5*beam_dia,2)))*n2;

    cout << "gamma: " << gamma <<endl;
    cout << "Central Freq: " << central_freq << endl;

    double fwhm = 15.0*pow(10,12);
    vector<double> phase_coef {0*pow(10,-15*2),0*pow(10,-15*3),pow(10,-15*4)};  //Taylor phase coef starting at quadratic term, in seconds^n

    Pulse pulse;
    Media media;


    media.set_gamma(gamma);  //
    media.set_length(length);
    media.set_number_steps(number_steps);
    media.set_beta2(beta2*pow(10,-15*2)*1000);
    media.set_beta3(0.0);
    media.set_beta4(0.0);
    media.set_B(0.0);

    pulse.set_freq(freq_range,central_freq,bins);
    pulse.set_gaussian_spectrum(fwhm,central_freq, 1.0);



    pulse.set_energy(energy);

    vector<complex<double>> fields = pulse.get_spectral_fields();

    for (int i = 0; i < pulse.get_bins(); i++) {
        fields[i] = polar(abs(fields[i])*10,arg(fields[i]));
    }

    pulse.set_spectral_fields(fields);


    Pulse pulse_out = pulse;
    cout << "Energy Before: " << pulse_out.get_energy() << endl;
    //pulse_out = media.propegate(pulse_out);
    cout << "Energy After: " << pulse_out.get_energy() << endl;
    cout << "Max B: " << media.get_B() << endl;
    cout << "Mid Freq: " << pulse_out.get_freq()[1000] << endl;


    pulse.save("initial.txt");
    pulse_out.save("spm.txt");



    return 0;

}
