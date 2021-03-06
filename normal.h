/*Declarations for normal number generator*/

#ifndef _NORMAL_H_
#define _NORMAL_H_

#include "Param.h"

class Normal {
public:
    Normal(Param &p);
    Normal(long seed1, long seed2);
    ~Normal();

    double draw(double mu, double sigma);

private:
    double gaussdev(int);
    unsigned short x1[3], x2[3];
};

#endif

