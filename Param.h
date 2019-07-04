// Param.h: interface for the Param class.
//
//////////////////////////////////////////////////////////////////////

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#ifndef _PARAM_H_
#define _PARAM_H_

#include <string.h>
#include "polyModel.h"

#ifndef _WIN32
#ifndef stricmp
#define stricmp strcasecmp
#endif
#endif

class Param {

public:
    static const char* CATS_VERSION_STRING; // = "2.2"
    static const int dist_count;
    static const char* dname[];
    static const char *defaultModelFilenamesFile, *defaultDistDist, *uniformHybridDist;

    enum dist_type {
        ARBITRARY, ARBITRARY_NPV, ARBITRARY_UPV, MATCHING, PATHS, REGIONS, REGIONS_NPV, 
        REGIONS_UPV, SCHEDULING, L1, L2, L3, L4, L5, L6, L7, L8, UNASSIGNED
    };

    dist_type distribution;
    static double ERand(unsigned short x[3]);
    static int Round(double val);

    // functions
    void Seed();
    static double DRand();
    static long LRand();
    static double DRand(double minval, double maxval);
    static long LRand(long minval, long maxval);

    Param(char** argv, int argc);
    virtual ~Param();
    char *outputSettings(bool tofile);
    int output_frequency;
    bool *argRecognized;

    // general variables
    int num_runs;
    bool cplex_output;
    int random_seed, random_seed2;
    const char *filename;
    const char *cplex_filename;
    const char *feat_filename;
    bool calculating;
    bool converting;

    // if non-NULL, param_dist
    // gives an estimate to the true dist based on parameters

    // param_feat_polys define a distribution based on
    // min of feature polynomials

    int num_feat_polys;
    const char *model_filenames_file;
    const char *dist_dist_file;
    PolyModel **param_dists;
    PolyModel *feat_polys[10];
    bool no_normalization;
    double *dist_weight;

//  const char *feature_norm_vect_file;
    double *feature_norm_means;
    double *feature_norm_devs;

    int numWeightedSamples;
    int num_bids;
    bool random_bids;
    int min_bids, max_bids;   

    int num_goods;
    bool random_goods;
    int min_goods, max_goods;
    bool remove_dominated_bids;
    bool output_parameter_settings;
    bool integer_prices;
    bool parse_error;
    bool verbatim;
    bool random_parameters;
    int bid_alpha; // for rounding to integer

private:
    char *output_settings_buffer;
    void parseModelFiles();
    void parseDistDistFile();
    void usage();
    void paramDistUsage();
};

#endif

