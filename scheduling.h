#ifndef _SCHEDULING_H
#define _SCHEDULING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "Distribution.h"
#include "Param.h"
#include "BidSet.h"
#include "bid.h"

#define sqr(x) ((x)*(x))
#define cub(x) ((x)*(x)*(x))
#define abs(x) ((x) > 0 ? (x) : -(x))
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

class Scheduling : public Distribution {
public:
    // constructor
    Scheduling(Param &p, char *argv[], int argc);

    // destructor
    ~Scheduling();

    // output values of all variables related to the distribution
    char *outputSettings(bool tofile);

    // generate the bidset
    BidSet *generate(Param &p);
    void randomizeParams(Param &p);

    // display usage information
    static void usage();

protected:
    static const double add_low, add_high,
        dev_low, dev_high,
        prob_low, prob_high;
    static int max_len_low, max_len_high;

    // parameters
    int max_length;
    double deviation;
    double prob_additional_deadline;
    double additivity;

    // create a bid on the range start... end with value val
    Bid *makeBid(int start, int end, double val);
    char *output_buffer;
    static double sched_param_mins[];
    static double sched_param_maxs[];
};

#endif // _SCHEDULING_H

