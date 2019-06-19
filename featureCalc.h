#ifndef _FEATURE_CALC_H
#define _FEATURE_CALC_H

#ifdef USE_CPLEX
#include <cplex.h>
#else
#include "lp_solve_4.0/lpkit.h"
#include "lp_solve_4.0/patchlevel.h"
#endif

#include "BidSet.h"

struct FeatureCalc {

    struct LPData {
#ifdef USE_CPLEX
        CPXENVptr env;
        CPXLPptr lp;

        double *lb;
        double *ub;

        double *rhs;
        char *sense;

        int *rmatbeg;
        int *rmatind;
        double *rmatval;
        double *obj;
#endif

        int status;
        double *xvals;

        LPData();
        void deallocate();
        void allocate(BidSet *b);
    } lpd;

    static const char *featureNames[];
    static const int numFeatures = 32;
    double featureVals[numFeatures];
    int featureNum;

    int num_bids;
    int num_goods;

    static const char *realistFeatureNames[];
    static const int numRealismFeatures = 3;
    double realismFeatureVals[numFeatures];
    int realismFeatureNum;

    void calcFeatures(BidSet &b);
    float testSmallWorldBidGraph(int bids, int **bidGraph, int calc);
    void LPfeatures(BidSet& bids);
    double LPsolve(BidSet *b);

    void writeFeatures(const char *filename);
    void writeFeatNames(const char *filename);
};

#endif
