/***************************************************************
 * Return (0,1) normally distributed random variable
 * (C) Numerical Recipes in C
 * 
 ***************************************************************/

#include <stdlib.h>
#include <math.h>
#include "normal.h"

Normal::Normal(Param &p)
{
    x1[0] = (unsigned short) (p.random_seed >> 16);
    x1[1] = (unsigned short) (p.random_seed & 0xFFFF);
    x1[2] = (unsigned short) (0x330E);

    x2[0] = (unsigned short) (p.random_seed2 >> 16);
    x2[1] = (unsigned short) (p.random_seed2 & 0xFFFF);
    x2[2] = (unsigned short) (0x330E);
}

Normal::Normal(long seed1, long seed2)
{
    x1[0] = (unsigned short) (seed1 >> 16);
    x1[1] = (unsigned short) (seed1 & 0xFFFF);
    x1[2] = (unsigned short) (0x330E);

    x2[0] = (unsigned short) (seed2 >> 16);
    x2[1] = (unsigned short) (seed2 & 0xFFFF);
    x2[2] = (unsigned short) (0x330E);
}

Normal::~Normal() {}

double Normal::gaussdev(int reset)
{
    static int iset = 0;
    static double gset = 0;
    double fac, rsq, v1, v2;

    if (reset) {
        iset = 0; /*reinit*/
    }

    if (iset == 0) {
        do {
            v1 = 2.0 * Param::ERand(x1) - 1.0;
            v2 = 2.0 * Param::ERand(x2) - 1.0;

            rsq = v1 * v1 + v2 * v2;
        }
        while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq) / rsq);
        //  printf("VS: %lf\t%lf\t%lf\t%lf\n", v1, v2, rsq, fac);*/
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    }
    else {
        iset = 0;
        return gset;
    }
}

double Normal::draw(double mu, double sigma)
{
    return gaussdev(0) * sigma + mu;
}

