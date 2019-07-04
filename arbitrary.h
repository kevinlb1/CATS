#ifndef _ARBITRARY_H
#define _ARBITRARY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Distribution.h"
#include "regions.h"

class Arbitrary : public Regions {
public:

    // constructor
    Arbitrary(Param &p, char * argv[], int argc);

    // destructor
    ~Arbitrary();

protected:

    // build the arbitrary adjacency map
    virtual int BuildMap();

    // decide how much weight to add if location a and b are in the same bid
    virtual double weightFromLocations(int a, int b);

};

#endif // _ARBITRARY_H

