#define _USE_MATH_DEFINES
#include "Media.h"
#include "Pulse.h"

#include <math.h>
Media::Media()
{
    //ctor
}

Media::~Media()
{
    //dtor
}

//Functions to set values
void Media::set_length(double length) {
    m_length = length;
}
void Media::set_gamma(double gamma) {
    m_gamma = gamma;
}
void Media::set_number_steps(int number_steps) {
    m_number_steps = number_steps;
}

void Media::set_beta2(double beta2) {
    m_beta2 = beta2;
}
void Media::set_beta3(double beta3) {
    m_beta3 = beta3;
}
void Media::set_beta4(double beta4) {
    m_beta4 = beta4;
}

//Function to read values

double Media::get_length(void) {
    return m_length;
}

double Media::get_gamma(void) {
    return m_gamma;
}

int Media::get_number_steps(void) {
    return m_number_steps;
}

double Media::get_beta2(void) {
    return m_beta2;
}

double Media::get_beta3(void) {
    return m_beta3;
}

double Media::get_beta4(void) {
    return m_beta4;
}

//Functions to run simulations

double Media::get_B(void) {
    return m_Bint;
}

void Media::set_B(double B_int) {
    m_Bint = B_int;
}

Pulse Media::m_propegation_step(Pulse pulse) {

    return pulse;
};

Pulse Media::propegate(Pulse pulse_out) {
    //meat and potatoes right here
    double step_size = m_length/m_number_steps;
    vector<double> phase_coef {step_size*0.5*m_beta2,step_size*0.5*m_beta3,step_size*0.5*m_beta4};
    vector<complex<double>> t_field;
    vector<double> initial_phase;
    double phase_shift;
    double max_phase = 0.0;

    double central_freq = pulse_out.get_central_freq();

    for (int i = 0; i < m_number_steps; i++) {
        //do half steps off dispersion


        pulse_out.add_taylor_phase(phase_coef);

        t_field = pulse_out.get_temporal_fields();
        initial_phase = pulse_out.get_temporal_phase();
        for (int i = 0; i < pulse_out.get_bins(); i++) {
            phase_shift = m_gamma*pow(abs(t_field[i]),2)*step_size;

            if (phase_shift > max_phase) {
                max_phase = phase_shift;
            }

            t_field[i] = polar(abs(t_field[i]),initial_phase[i] + phase_shift);
        }
        pulse_out.set_temporal_fields(t_field);
        pulse_out.add_taylor_phase(phase_coef);
        m_Bint += max_phase;
    }

    cout << "Central Freq Prop: " << pulse_out.get_central_freq() << endl;

  //  pulse_out.set_energy(100.0);

    return pulse_out;
}



