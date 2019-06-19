// Param.cpp: implementation of the Param class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <assert.h>
#ifndef _WIN32
#include <unistd.h>
#endif

#define max(a,b) (a>b?a:b)

// this must preceed including Param.h
// #define DNAME_DEF
#include "Param.h"
#include "Legacy.h"
#include "regions.h"
#include "arbitrary.h"
#include "scheduling.h"
#include "paths.h"
#include "matching.h"
#include "featureCalc.h"

// static members

const char* Param::CATS_VERSION_STRING = "2.1";

const int Param::dist_count = 17; // not counting "unassigned"
const char* Param::dname[] = {"arbitrary", "arbitrary-npv", "arbitrary-upv", "matching", "paths", "regions", "regions-npv", "regions-upv", "scheduling", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "unassigned"};

const char* Param::defaultModelFilenamesFile = "default_model_file";
const char* Param::defaultDistDist = "default_dist_dist";
const char* Param::uniformHybridDist = "uniform_hybrid";

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// parse the parameters: done in the constructor
// all distribution-specific parameters are stored in a table as strings, and parsed in the distribution itself

Param::Param(char** argv, int argc)
{
    // initalize variables to defaults
    distribution = UNASSIGNED;

    num_runs = 1;
    cplex_output = false;
    random_seed = 0; //- make results replicatable
    random_seed2 = 0; //- make results replicatable
    filename = NULL;
    cplex_filename = NULL;
    feat_filename = NULL;

    num_feat_polys = 0;
    model_filenames_file = NULL;

    dist_dist_file = NULL;
    dist_weight = NULL;

    argRecognized = new bool[argc];

    param_dists = new PolyModel*[dist_count];
    for (int i = 0; i < dist_count; i++)
        param_dists[i] = NULL;

    //    feature_norm_vect_file = 0;

    numWeightedSamples = 50;

    output_settings_buffer = NULL;
    num_bids = 0;
    num_goods = 0;
    remove_dominated_bids = true;
    output_parameter_settings = true;
    integer_prices = false;
    parse_error = false;
    verbatim = false;
    output_frequency = 20; // no parameter for this one
    random_parameters = false;
    random_goods = false;
    random_bids = false;
    bid_alpha = 1000;
    converting = false;
    calculating = false;

    bool default_hard = false;
    bool needsHelp = false;

    // no parameters passed
    if (argc == 1) {
        usage();
        parse_error = true;
        return;
    }

    argRecognized[0] = true;
    for (int i = 1; i < argc; i++)
        argRecognized[i] = false;

    // parse parameters
    for (int i = 1; i < argc; i++) {
        // now read in values
        if (argv[i][0] == '-' || argv[i][0] == '/') {
            char *in = &(argv[i][1]); // skip the first character

            // model_file_help
            if (!stricmp(in, "model_file_help")) {
                paramDistUsage();
                exit(0);
            }

            else if (!stricmp(in, "default_hard")) {
                argRecognized[i] = true;
                default_hard = true;
            }

                // distribution
            else if (!stricmp(in, "d")) {
                argRecognized[i] = true;
                if (i + 1 >= argc ||
                    distribution != UNASSIGNED ||
                    dist_dist_file != NULL) {
                    usage();
                    parse_error = true;
                    return;
                }

                for (int t = 0; t < dist_count; t++)
                    if (stricmp(argv[i + 1], dname[t]) == 0)
                        distribution = (dist_type) t;

                if (distribution == UNASSIGNED) {
                    printf("Unknown distribution: %s\n", argv[i + 1]);
                    parse_error = true;
                    return;
                }
                argRecognized[i + 1] = true;

                i++;
            }

                // dist_dist
            else if (!stricmp(in, "dist_dist")) {
                argRecognized[i] = true;
                if (i + 1 >= argc ||
                    distribution != UNASSIGNED ||
                    dist_dist_file != NULL) {
                    usage();
                    parse_error = true;
                    return;
                }

                dist_dist_file = argv[i + 1];
                argRecognized[i + 1] = true;
                i++;
            }

                // num_runs
            else if (!stricmp(in, "n")) {
                argRecognized[i] = true;
                if (i + 1 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                num_runs = atoi(argv[i + 1]);
                if (!num_runs) {
                    usage();
                    parse_error = true;
                    return;
                }
                i++;
                argRecognized[i] = true;
            }

                // cplex_output
            else if (!stricmp(in, "cplex")) {
                argRecognized[i] = true;
                cplex_output = true;
            }

                // bid_alpha
            else if (!stricmp(in, "bid_alpha")) {
                argRecognized[i] = true;
                if (i + 1 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                bid_alpha = atoi(argv[i + 1]);
                if (!bid_alpha) {
                    usage();
                    parse_error = true;
                    return;
                }
                i++;
                argRecognized[i] = true;
            }

                // seed1
            else if (!stricmp(in, "seed")) {
                argRecognized[i] = true;
                if (i + 1 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                random_seed = atoi(argv[i + 1]);
                if (!random_seed) {
                    usage();
                    parse_error = true;
                    return;
                }
                i++;
                argRecognized[i] = true;
            }

                // seed2
            else if (!stricmp(in, "seed2")) {
                argRecognized[i] = true;
                if (i + 1 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                random_seed2 = atoi(argv[i + 1]);
                if (!random_seed2) {
                    usage();
                    parse_error = true;
                    return;
                }
                i++;
                argRecognized[i] = true;
            }

                // bids
            else if (!stricmp(in, "bids")) {
                argRecognized[i] = true;
                if (i + 1 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                num_bids = atoi(argv[i + 1]);
                if (!num_bids) {
                    usage();
                    parse_error = true;
                    return;
                }
                i++;
                argRecognized[i] = true;
            }

                // random bids
            else if (!stricmp(in, "random_bids")) {
                argRecognized[i] = true;
                if (i + 2 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                random_bids = true;
                min_bids = atoi(argv[++i]);
                argRecognized[i] = true;
                max_bids = atoi(argv[++i]);
                argRecognized[i] = true;
                if (!min_bids || !max_bids) {
                    usage();
                    parse_error = true;
                    return;
                }
            }

                // goods
            else if (!stricmp(in, "goods")) {
                argRecognized[i] = true;
                if (i + 1 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                num_goods = atoi(argv[i + 1]);
                if (!num_goods) {
                    usage();
                    parse_error = true;
                    return;
                }
                i++;
                argRecognized[i] = true;
            }

                // random goods
            else if (!stricmp(in, "random_goods")) {
                argRecognized[i] = true;
                if (i + 2 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                random_goods = true;
                min_goods = atoi(argv[++i]);
                argRecognized[i] = true;
                max_goods = atoi(argv[++i]);
                argRecognized[i] = true;
                if (!min_goods || !max_goods) {
                    usage();
                    parse_error = true;
                    return;
                }
            }

                // no_dom_check
            else if (!stricmp(in, "no_dom_check")) {
                argRecognized[i] = true;
                remove_dominated_bids = false;
            }

                // random_parameters
            else if (!stricmp(in, "random_parameters")) {
                argRecognized[i] = true;
                random_parameters = true;
            }

                // model_filenames
            else if (!stricmp(in, "model_filenames")) {
                argRecognized[i] = true;
                if (i + 1 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                model_filenames_file = argv[i + 1];
                i++;
                argRecognized[i] = true;
            }

                // filename
            else if (!stricmp(in, "filename")) {
                argRecognized[i] = true;
                if (i + 1 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                filename = argv[i + 1];
                cplex_filename = filename; // in the future this could be a separate parameter
                i++;
                argRecognized[i] = true;
            }

                // verbatim
            else if (!stricmp(in, "verbatim")) {
                argRecognized[i] = true;
                verbatim = true;
            }

                // integer_prices
            else if (!stricmp(in, "int_prices")) {
                argRecognized[i] = true;
                integer_prices = true;
            }

                // output_parameter_settings
            else if (!stricmp(in, "no_output")) {
                argRecognized[i] = true;
                output_parameter_settings = false;
            }

                // cats2cplex
            else if (!stricmp(in, "cats2cplex")) {
                argRecognized[i] = true;
                if (i + 2 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                filename = argv[i + 1];
                cplex_filename = argv[i + 2];
                converting = true;
                verbatim = true;
                argRecognized[i + 1] = true;
                argRecognized[i + 2] = true;
                i += 2;
            }

                // calc_features
            else if (!stricmp(in, "calc_features")) {
                argRecognized[i] = true;
                if (i + 2 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                filename = argv[i + 1];
                feat_filename = argv[i + 2];
                calculating = true;
                argRecognized[i + 1] = true;
                argRecognized[i + 2] = true;
                i += 2;
            }

            else if (!stricmp(in, "num_weighted_samples")) {
                argRecognized[i] = true;
                if (i + 1 >= argc) {
                    usage();
                    parse_error = true;
                    return;
                }
                numWeightedSamples = atoi(argv[i + 1]);
                if (numWeightedSamples == 0) {
                    printf("Error: num_weighted_samples must be given a positive integer argument.\n");
                    exit(1);
                }
                i++;
                argRecognized[i] = true;
            }

            else if (!stricmp(in, "no_normalization")) {
                argRecognized[i] = true;
                no_normalization = true;
            }

            else if (!stricmp(in, "uniform_hybrid")) {
                argRecognized[i] = true;
                dist_dist_file = uniformHybridDist;
            }

                // help

            else if (!stricmp(in, "h") || !stricmp(in, "help") || !stricmp(in, "?")) {
                argRecognized[i] = true;
                needsHelp = true;
            }
        }
    }

    if (needsHelp) {
        usage();
        parse_error = true;
        return;
    }

    if (converting || calculating) return;

    if (distribution == L1 || distribution == L5) {
        printf("\nWarning: the L1 and L5 distributions may take extremely long to generate\n");
        printf("  any significant number of non-dominated bids.  Are you sure you want\n");
        printf("  to continue? (y/n):");

        char response;

        cin >> response;
        while (response != 'n' && response != 'N' &&
               response != 'y' && response != 'Y') {
            printf("Please enter \'y\' or \'n\'.\n");
            cin >> response;
        }

        if (response == 'n' || response == 'N')
            exit(0);
    }

    if (default_hard) {
        if (num_bids == 0 && num_goods == 0 && !random_goods && !random_bids) {
            num_bids = 1000;
            num_goods = 256;
        }
        else if (num_bids != 1000 || num_goods != 256 || random_goods || random_bids) {
            printf("Error: the built-in hardness models are for 1000 bids and 256 goods.\n");
            printf("Try again without specifying # of goods or bids and these values will be used.\n");
            exit(1);
        }

        model_filenames_file = defaultModelFilenamesFile;
        if (distribution == UNASSIGNED) {
            if (dist_dist_file != NULL) {
                printf("Error: a distribution over distributions cannot be used\n");
                printf("with the -default_hard flag.  Choose a specific distribution\n");
                printf("or do not specify one, and the hybrid distribution will be used.\n");
                exit(1);
            }

            dist_dist_file = defaultDistDist;
        }
        else if (distribution == PATHS || distribution == L1 ||
                 distribution == L5 || distribution == L8) {
            printf("Error: the distributions L1, L5, L8 and PATHS do not currently have\n");
            printf("built-in hardness models.  Choose a different distribution, or do not\n");
            printf("specify one, and the hybrid distribution over all others will be used.\n");
            exit(1);
        }
    }

    if (model_filenames_file != NULL)
        parseModelFiles();

    if (dist_dist_file != NULL)
        parseDistDistFile();

    // check for errors in the input
    if ((!random_bids && num_bids == 0) ||
        (!random_goods && num_goods == 0) ||
        (distribution == UNASSIGNED && dist_weight == NULL)) {
        printf("Error: number of bids and goods, and type of distribution must be entered\n");
        printf("on command line.  Type \"cats -help\" for more information.\n");
        exit(1);
    }

    if (num_runs > 1)
        verbatim = false;

    if (!filename) {
        filename = "";
        cplex_filename = "";
    }

    // other initialization

    if (!random_seed)
        random_seed = (long) time(NULL) + (long) getpid();

#ifdef DEBUG
    printf("Seed: %d\n", random_seed);
#endif

    if (!random_seed2) {
        // wait a bit of time in case we just took seed1 from the clock
        for (int t = 0; t < random_seed % 100000; t++);
        random_seed2 = time(NULL) % 99873 - (long) getpid();
    }

    Seed();
}

Param::~Param()
{
    if (output_settings_buffer) delete[] output_settings_buffer;
    if (param_dists) {
        for (int i = 0; i < dist_count; i++)
            if (param_dists[i])
                delete param_dists[i];
        delete[] param_dists;
    }

    if (dist_weight)
        delete[] dist_weight;

    for (int i = 0; i < num_feat_polys; i++)
        delete feat_polys[i];

    delete[] argRecognized;
}

void Param::parseModelFiles()
{

    ifstream modelFilenamesStream(model_filenames_file);
    if (!modelFilenamesStream) {
        printf("Error: could not open model_filenames file \"%s\"\n", model_filenames_file);
        exit(1);
    }

    char buf[10000];

    while (!modelFilenamesStream.eof()) {

        modelFilenamesStream >> buf;

        if (!stricmp(buf, "-feat_poly")) {
            if (modelFilenamesStream.eof()) {
                paramDistUsage();
                parse_error = true;
                return;
            }
            if (num_feat_polys == 10) {
                printf("Error: max 10 feat_polys exceeded\n");
                parse_error = true;
                return;
            }
            modelFilenamesStream >> buf;
            if (strlen(buf) == 0) {
                paramDistUsage();
                parse_error = true;
                return;
            }
            feat_polys[num_feat_polys++] = new PolyModel(buf);

        }
        else if (!stricmp(buf, "-param_dist")) {
            if (modelFilenamesStream.eof()) {
                paramDistUsage();
                parse_error = true;
                return;
            }
            modelFilenamesStream >> buf;
            int distIndex = -1;
            for (int i = 0; i < dist_count; i++)
                if (!stricmp(buf, dname[i])) {
                    distIndex = i;
                    break;
                }

            if (distIndex == -1 ||
                distIndex == Param::REGIONS ||
                distIndex == Param::ARBITRARY) {
                paramDistUsage();
                parse_error = true;
                return;
            }

            if (modelFilenamesStream.eof()) {
                paramDistUsage();
                parse_error = true;
                return;
            }

            modelFilenamesStream >> buf;
            if (strlen(buf) == 0) {
                paramDistUsage();
                parse_error = true;
                return;
            }
            param_dists[distIndex] = new PolyModel(buf);

        }
    }
}

void Param::parseDistDistFile()
{

    ifstream distDistFileStream(dist_dist_file);
    if (!distDistFileStream) {
        printf("Error: could not open dist_dist file \"%s\"\n", dist_dist_file);
        exit(1);
    }

    char buf[10000];

    dist_weight = new double[dist_count];

    double weightSum = 0;

    for (int i = 0; i < dist_count; i++) {
        if (i == ARBITRARY || i == REGIONS) {
            dist_weight[i] = 0;
            continue;
        }

        if (distDistFileStream.eof()) {
            printf("Error: not enough entries in dist_dist file.  Expected %d\n", dist_count - 2);
            printf("Type: \"cats -help\" for more info on the -dist_dist flag.\n");
            exit(1);
        }

        distDistFileStream >> buf;
        dist_weight[i] = atof(buf);

        if (dist_weight[i] < 0) {
            printf("Error: dist_dist weights must be >= 0.\n");
            printf("Type: \"cats -help\" for more info on the -dist_dist flag.\n");
            exit(1);
        }
        else
            weightSum += dist_weight[i];
    }

    if (!distDistFileStream.eof()) {
        printf("Error: extra entries in dist_dist file.  Expected %d\n", dist_count - 2);
        printf("Type: \"cats -help\" for more info on the -dist_dist flag.\n");
        paramDistUsage();
        exit(1);
    }

    if (dist_weight[L1] > 0 || dist_weight[L5] > 0) {
        printf("\nWarning: the L1 and L5 distributions may take extremely long\n");
        printf("  to generate any significant number of non-dominated bids, but\n");
        printf("  they have been assigned positive weight in the dist_dist.\n");
        printf("  Are you sure you want to continue? (y/n):");

        char response;

        cin >> response;
        while (response != 'n' && response != 'N' &&
               response != 'y' && response != 'Y') {
            printf("Please enter \'y\' or \'n\'.\n");
            cin >> response;
        }

        if (response == 'n' || response == 'N')
            exit(0);
    }


    if (weightSum < 0.000000001) {
        printf("Error: some dist. must have positive weight in dist_dist.\n");
        printf("Type: \"cats -help\" for more info on the -dist_dist flag.\n");
        exit(1);
    }

    //  fprintf(stderr, "Sum: %g Weights: ",weightSum);
    // normalize
    for (int i = 0; i < dist_count; i++) {
        //      fprintf(stderr,"%g\t", dist_weight[i]);
        dist_weight[i] /= weightSum;
    }
    //  fprintf(stderr, "\n");

}


// allow the use of a different random number generator for different platforms

double Param::DRand()
{
#ifndef _WIN32
    return drand48();
#else
    return (double) rand() / RAND_MAX;
#endif
}

// random number on a range

double Param::DRand(double minval, double maxval)
{
    return Param::DRand() * (maxval - minval) + minval;
}

long Param::LRand()
{
#ifndef _WIN32
    return lrand48();
#else
    return rand();
#endif
}

// random number on a range

long Param::LRand(long minval, long maxval)
{
    return Param::LRand() % (maxval - minval + 1) + minval;
}

void Param::Seed()
{
#ifndef _WIN32
    srand48(random_seed);
#else
    srand(random_seed);
#endif
}

// write out the parameter settings

char *Param::outputSettings(bool tofile)
{
    if (!output_settings_buffer)
        output_settings_buffer = new char[500];

    char dist_dist_str[500], param_dist_str[500], feat_poly_str[500];

    if (dist_weight) {
        sprintf(dist_dist_str, "Distribution over distributions taken from file \"%s\".", dist_dist_file);
        if (tofile)
            strcat(dist_dist_str, "; ");
        else
            strcat(dist_dist_str, "\n");
    }
    else
        dist_dist_str[0] = 0;

    if (param_dists[distribution]) {
        sprintf(param_dist_str, "Parameters generated from distribution in file \"%s\".", param_dists[distribution]->filename);
        if (tofile)
            strcat(param_dist_str, "; ");
        else
            strcat(param_dist_str, "\n");
    }
    else
        param_dist_str[0] = 0;

    if (num_feat_polys > 0) {
        sprintf(feat_poly_str, "Problems generated from distribution over features in files: ");
        if (!tofile)
            strcat(feat_poly_str, "\n   ");
        for (int i = 0; i < num_feat_polys; i++) {
            strcat(feat_poly_str, feat_polys[i]->filename);
            if (i < num_feat_polys - 1)
                if (tofile)
                    strcat(feat_poly_str, ", ");
                else
                    strcat(feat_poly_str, "\n   ");
        }
        if (tofile)
            strcat(feat_poly_str, "; ");
        else
            strcat(feat_poly_str, "\n");
    }
    else
        feat_poly_str[0] = 0;

    if (tofile) sprintf(output_settings_buffer, "%% Goods: %d; Bids: %d; Distribution: %s; Runs: %d; %s%s%sSeed: %d", num_goods, num_bids, dname[distribution], num_runs, dist_dist_str, param_dist_str, feat_poly_str, random_seed);
    else sprintf(output_settings_buffer, "Goods: %d\nBids: %d\nDistribution: %s\nRuns: %d\n%s%s%sSeed: %d\n", num_goods, num_bids, dname[distribution], num_runs, dist_dist_str, param_dist_str, feat_poly_str, random_seed);

    return output_settings_buffer;
}

int Param::Round(double val)
{
    return (int) (val + (double) 0.5);
}

// support for ERand under windows compilers that can't handle it

double Param::ERand(unsigned short x[3])
{
#ifndef _WIN32
    return erand48(x);
#else
    return Param::DRand();
#endif
}

void Param::usage()
{
    // 80 cols:  12345678901234567890123456789012345678901234567890123456789012345678901234567890
    int count = printf("CATS v%s (http://robotics.stanford.edu/CATS)\n", CATS_VERSION_STRING);
    int count2 = printf("Kevin Leyton-Brown, Mark Pearson, Galen Andrew, Yoav Shoham; Stanford University\n");
    int t; // -- NOTE: MSVC uses older syntax, so cannot declare it withn a forloop :(
    for (t = 0; t < max(count, count2); t++) printf("=");
    printf("\n\n");

    if (distribution == UNASSIGNED) {
        printf("What follows are the general parameters for CATS.\n");
        printf("To see parameters for a specific distribution, select a distribution with -help\n\n");
    }

    printf("Required Parameters (no default values):\n");
    printf("  -d [no default]: selects a distribution.  Valid options (without quotes) are:");
    count = 80;

    for (t = 0; t < dist_count; t++) {
        if (count >= 70) {
            printf("\n     ");
            count = 0;
        }
        count += printf("\"%s\"; ", dname[t]);
    }

    printf("\n       --OR--\n");
    printf("  -dist_dist [filename]: specifies a distribution over distributions\n");
    printf("       --OR--\n");
    printf("  -uniform_hybrid: uses the hybrid distribution selecting each built-in distribution\n");
    printf("     except for L1 and L5 with equal probability.\n\n");

    printf("  -goods [no default]: number of goods.\n");
    printf("       --OR--\n");
    printf("  -random_goods [no default] [no default]: min/max number of goods.\n\n");

    printf("  -bids [no default]: number of bids to generate.\n");
    printf("       --OR--\n");
    printf("  -random_bids [no default] [no default]: min/max number of bids.\n\n");

    printf("The argument to -dist_dist is a file containing a vector of weights given\n");
    printf("to each distribution in the following order:\n");
    for (int i = 0; i < dist_count; i++)
        if (i != REGIONS &&
            i != ARBITRARY)
            printf("[%s] ", dname[i]);
    printf("\n(Note that regions and arbtirary *must* be split into -upv and -npv.)\n");
    printf("If the weights do not sum to 1, they will be normalized.\n\n");

    printf("Optional Parameters:\n");
    printf("  -n [1]: number of instance files that you want CATS to generate\n");
    printf("  -seed [taken from clock]: random seed\n");
    printf("  -seed2 [taken from clock]: random seed 2: used only for normal distributions\n");
    printf("  -bid_alpha [1000]: multiply bids by this before rounding, when -int_prices is set\n");
    printf("  -num_weighted_samples [50]: number of samples to take in generating from");
    printf("      distribution over features\n");

    printf("\nOptional Flags (those with arguments have no default values):\n");
    printf("  -cplex: write a CPLEX-compatible .lp file as well as a CATS data file\n");
    printf("  -no_dom_check: turn off removal of dominated bids.  The default CATS behavior is to\n");
    printf("      keep generating bids until the given # of non-dominated bids have been made.\n");
    printf("  -int_prices: write integer prices in bids instead of real-valued prices\n");
    printf("  -filename:  set filename prefix for generated instances\n");
    printf("  -no_output: do not write unnecessary output to the screen\n");
    printf("  -verbatim: do not append indices to the filename.  works only with -n 1\n");
    printf("  -random_parameters: choose parameters uniformly from their ranges\n");
    printf("  -model_filenames: file containing filenames of polynomially defined distributions\n");
    printf("      over distribution parameters and/or polynomial functions of instance\n");
    printf("      features for weighted sampling\n");
    printf("  -model_file_help: print help screen for model file features)\n");
    printf("  -no_normalization: turns off normalization of all polynomially defined\n");
    printf("      If the polynomials do not represent proper probability density functions\n");
    printf("      over their ranges, behavior is undefined.\n");
    printf("  -default_hard: weight problems for hardness using the built-in hardness model\n");
    printf("  -help, -h, -?: print this usage screen\n");

    // more information about supported distributions:
    if (distribution != UNASSIGNED) {
        printf("\n");
        switch (distribution) {
        case L1:
        case L2:
        case L3:
        case L4:
        case L5:
        case L6:
        case L7:
        case L8:
            Legacy::usage();
            break;

        case ARBITRARY:
        case ARBITRARY_NPV:
        case ARBITRARY_UPV:
            Regions::usage(true);
            break;

        case MATCHING:
            Matching::usage();
            break;

        case PATHS:
            Paths::usage();
            break;

        case REGIONS:
        case REGIONS_NPV:
        case REGIONS_UPV:
            Regions::usage(false);
            break;

        case SCHEDULING:
            Scheduling::usage();
            break;

        default:
            assert(false);
        }
    }

    printf("\nCATS also has utilities for converting from CATS output to CPLEX output\n");
    printf("  and for computing features of an instance.\n");
    printf("These are run with the following flags:\n");
    printf("  -cats2cplex [input filename] [output filename]\n");
    printf("  -calc_features [input filename] [output filename]\n\n");
}

void Param::paramDistUsage()
{
    printf("The -dist_file flag specifies a file containing paths to files with polynomially defined\n");
    printf("distributions over distribution parameters or a polynomial\n");
    printf("weighting function over instance features.\n\n");

    printf("The format is a list of the following flags:\n");
    printf("  -param_dist [dist_name] [filename]: specifies polynomial distribution\n");
    printf("      over parameters for given distribution.  Note that\n");
    printf("      regions and arbitrary cannot be specified without the -npv or -upv\n");
    printf("      suffix, since the two price distributions take different parameters.\n");
    printf("  -feat_poly [filename]: specifies feature weight function-- if several\n");
    printf("      are specified, the weight will be the minimum of those.\n\n");

    printf("Arguments to -param_dist and -feat_poly must be polynomial model files\n");
    printf("with the following entries (separated by whitespace):\n");
    printf("the degree of the polynomial, the number of variables, then the coefficients\n");
    printf("ordered first according to the degree of the term (inc. constant - degree 0)\n");
    printf("and then by the constituent variables, listed from largest to smallest.\n");
    printf("For example, the second-order coefficients would be in the following order,\n");
    printf("where c_mn is the coefficient of the term with variables m and n:\n");
    printf("c_00 c_10 c_11 c_20 c_21 c_22 c_30 c_31 c_32 c_33 ... c_d0 c_d1 ... c_dd\n");

}
