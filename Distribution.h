// Distribution.h: implementation of the Distribution abstract class.
//
//////////////////////////////////////////////////////////////////////

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H
#endif

#include "Param.h"

class BidSet;

class Distribution {
public:
    Distribution() {
        probOfParams = 1.0;
    }

    // the main function to be overridden in subclasses
    virtual BidSet *generate(Param &p) = 0;
    virtual char *outputSettings(bool val) = 0;

    virtual ~Distribution() {
    }

    virtual void randomizeParams(Param &p) = 0;
    double probOfParams;

protected:
    // given a polynomial-defined distribution dist,
    // and a normalized estimate dist_estimate
    // this procedure fills the array params with
    // parameter settings drawn from that distribution.

    // mins and maxs are the minimum and maximum possible
    // values of each param

    void randomizeParamsWithDistribution(PolyModel *dist, int num_params, 
            double *param_mins, double *param_maxs, double *params, Param &p);
};

#endif

