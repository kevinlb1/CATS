#ifndef _MATCHING_H
#define _MATCHING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "bid.h"
#include "BidSet.h"
#include "Distribution.h"

#define sqr(x) ((x)*(x))
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

class Matching : public Distribution {
public:
    // constructor
    Matching(Param &p, char * argv [], int argc);

    // destructor
    ~Matching();

    // generate the bidset
    BidSet *generate(Param &p);
    void randomizeParams(Param &p);

    // output command-line usage
    static void usage();

    // output distribution parameters, either to the screen or to a file
    char *outputSettings(bool tofile);

protected:
    // constants
    static const int MAX_LONGEST_FLIGHT_LENGTH; // = 10
    static const int Z_MAX_DEV; // = 10
    static const int NUMCITIES; // = 4

    // parameters
    int m_max_airport_price;
    int m_early_takeoff_deviation;
    int m_late_takeoff_deviation;
    int m_early_land_deviation;
    int m_late_land_deviation;
    double m_deviation;
    double m_delay_coeff;
    double m_amount_late_coeff;

    // locals
    int num_times;

    typedef struct {
        double x, y;
    } point;

    point *cities;
    double **edges;
    double *cost;
    double max_l;
    char *output_buffer;

    // helper functions
    inline double longestFlightLength(double x) {
        return min(MAX_LONGEST_FLIGHT_LENGTH, x - 1);
    }

    // compute the L2 distance between two cities
    inline double RealDistance(int city1, int city2) {
        return sqrt(sqr(cities[city2].x - cities[city1].x) +
                sqr(cities[city2].y - cities[city1].y));
    }

    // initialize the cities.
    // regardless of the number of goods, there are only 4 cities, wich are the real
    // ones used in the FAA application

    void SetupCities();
};

#endif // _MATCHING_H

