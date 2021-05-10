#ifndef MEDIA_H
#define MEDIA_H

#include "Pulse.h"
class Media
{
    public:
        Media();
        virtual ~Media();

        void set_length(double length);
        void set_gamma(double gamma);
        void set_number_steps(int number_steps);
        void set_beta2(double beta2);
        void set_beta3(double beta3);
        void set_beta4(double beta4);
        void set_B(double B_int);

        Pulse propegate(Pulse);


        double get_length(void);
        double get_gamma(void);
        int get_number_steps(void);
        double get_beta2(void);
        double get_beta3(void);
        double get_beta4(void);
        double get_B(void);

    protected:

    private:
        double m_length;
        double m_gamma;
        double m_number_steps;

        double m_beta2;
        double m_beta3;
        double m_beta4;
        double m_Bint;

        Pulse m_propegation_step(Pulse);

};

#endif // MEDIA_H
