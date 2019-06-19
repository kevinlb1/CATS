// remember doing NORMAL for all these sampled Parameters

#include <assert.h>

// main file for all distributions
#include "BidSet.h"
#include "Param.h"
#include "Legacy.h"
#include "regions.h"
#include "arbitrary.h"
#include "scheduling.h"
#include "matching.h"
#include "paths.h"
#include "featureCalc.h"

#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))

void chooseDistribution(Param &p, Distribution *&d, int argc, char **argv)
{
    if (p.dist_weight != NULL) {
        double choiceVal = Param::DRand(0, 1.0);
        p.distribution = (Param::dist_type)0;
        choiceVal -= p.dist_weight[(int) p.distribution];
        while (choiceVal > 0) {
            p.distribution = Param::dist_type(p.distribution + 1);
            assert(p.distribution < Param::UNASSIGNED);
            choiceVal -= p.dist_weight[(int) p.distribution];
        }
    }

    switch (p.distribution) {

    case Param::ARBITRARY:
    case Param::ARBITRARY_NPV:
    case Param::ARBITRARY_UPV:
        d = new Arbitrary(p, argv, argc);
        break;

    case Param::MATCHING:
        d = new Matching(p, argv, argc);
        break;

    case Param::PATHS:
        d = new Paths(p, argv, argc);
        break;

    case Param::REGIONS:
    case Param::REGIONS_NPV:
    case Param::REGIONS_UPV:
        d = new Regions(p, argv, argc);
        break;

    case Param::SCHEDULING:
        d = new Scheduling(p, argv, argc);
        break;

    case Param::L1:
    case Param::L2:
    case Param::L3:
    case Param::L4:
    case Param::L5:
    case Param::L6:
    case Param::L7:
    case Param::L8:
        d = new Legacy(p, argv, argc);
        break;

    default:
        assert(false); // shouldn't be possible- but lets keep the compiler happy!
    }

    for (int i = 0; i < argc; i++) {
        if (!p.argRecognized[i]) {
            printf("Error: unknown parameter %s\n", argv[i]);
            exit(1);
        }
    }

    if (p.random_parameters || p.param_dists[p.distribution])
        d->randomizeParams(p);

}

int main(int argc, char *argv[])
{

    // parse Parameters
    Param p(argv, argc);
    if (p.parse_error) return 1;

    if (p.converting) {
        BidSet temp(p.filename);
        temp.writeCPLEXFile(p);
        return (0);
    }

    if (p.calculating) {
        FeatureCalc f;
        BidSet temp(p.filename);
        f.calcFeatures(temp);
        f.writeFeatures(p.feat_filename);
        f.writeFeatNames("featurenames");
        return 0;
    }

    Distribution *d;

    FeatureCalc f;

    /*  
    fprintf(stdout, "LEMAN:\n");
    p.feat_polys[2]->print();
    fprintf(stdout, "=======================================\n");
     */

    //cout.flush();

    for (int t = 0; t < p.num_runs; t++) {
        if (p.random_bids)
            p.num_bids = Param::LRand(p.min_bids, p.max_bids);

        if (p.random_goods)
            p.num_goods = Param::LRand(p.min_goods, p.max_goods);

        BidSet *b;

        if (p.num_feat_polys) { // generate according to distribution

            // generate samples from estimate
            Distribution **dists = new Distribution*[p.numWeightedSamples];
            BidSet **samples = new BidSet*[p.numWeightedSamples];
            double *scores = new double[p.numWeightedSamples];
            Param::dist_type* distTypes = new Param::dist_type[p.numWeightedSamples];
            double sumSquaredUnweightedScore = 0;
            double totalScore = 0;

            for (int sampleNum = 0; sampleNum < p.numWeightedSamples; sampleNum++) {

                /*  
                fprintf(stdout, "LEMAN %d:\n", sampleNum);
                p.feat_polys[2]->print();
                fprintf(stdout, "=======================================\n");
                 */

                chooseDistribution(p, d, argc, argv);

                // produce bids
                BidSet *sample = d->generate(p);

                // make integer prices
                if (p.integer_prices) {
                    for (int t = 0; t < sample->numBids(); t++)
                        sample->getBid(t)->amount = Param::Round(p.bid_alpha * sample->getBid(t)->amount);
                }

                f.calcFeatures(*sample);

                /*	
                fprintf(stderr, "####### FEATS: ");
                for(int i=0; i<f.numFeatures; i++)
                  fprintf(stderr, "%g, ", f.featureVals[i]);
                fprintf(stderr, "\n");
                 */

                /* // normalize features
                   for (int i=0; i<f.numFeatures; i++) {
                   f.featureVals[i] -= p.feature_norm_means[i];
                   f.featureVals[i] /= p.feature_norm_devs[i];
                   }*/

                double unweightedScore = p.feat_polys[0]->evalAt(f.featureVals);
                for (int i = 1; i < p.num_feat_polys; i++) {
                    double nextScore = p.feat_polys[i]->evalAt(f.featureVals);
                    if (nextScore < unweightedScore)
                        unweightedScore = nextScore;
                    // ==	  fprintf(stderr, "##### NEXT %g\n", nextScore);
                }

                unweightedScore = MAX(unweightedScore, 0);

                // == fprintf(stderr, "####### Score: %g\n", unweightedScore);

                sumSquaredUnweightedScore += sqr(unweightedScore);

                double factor = (p.dist_weight ? p.dist_weight[(int) p.distribution] : 1.0);
                scores[sampleNum] = unweightedScore / (d->probOfParams * factor);
                totalScore += scores[sampleNum];

                dists[sampleNum] = d;
                samples[sampleNum] = sample;
                distTypes[sampleNum] = p.distribution;

            }

            // NOTE this might be modified to stop when sum reaches
            // a certain point instead of the same # of samples every time
            /*
              if (sumSquaredUnweightedScore < 3)
              printf("Warning: sum squared unweighted scores low in sampling from polynomial-defined distribution: %f\n", sumSquaredUnweightedScore);
             */

            double draw = Param::DRand(0.0, totalScore);
            int selectedIndex = -1;
            do {
                assert(selectedIndex + 1 < p.numWeightedSamples);
                draw -= scores[++selectedIndex];
            }
            while (draw > 0);

            // == fprintf(stderr, "####### Selected %d: %g\n", selectedIndex, scores[selectedIndex]);
            d = dists[selectedIndex];
            b = samples[selectedIndex];
            //  p.distribution = distTypes[selectedIndex];

            for (int i = 0; i < p.numWeightedSamples; i++)
                if (i != selectedIndex) {
                    delete dists[i];
                    delete samples[i];
                }

            delete[] dists;
            delete[] samples;
            delete[] scores;
            delete[] distTypes;

        }
        else {

            chooseDistribution(p, d, argc, argv);

            // produce bids
            b = d->generate(p);

            // make integer prices
            if (p.integer_prices) {
                for (int t = 0; t < b->numBids(); t++)
                    b->getBid(t)->amount = Param::Round(p.bid_alpha * b->getBid(t)->amount);
            }
        }

        // output
        if (p.output_parameter_settings)
            printf("CATS running %d of %d....\n\n%s\n%s\n", t + 1, p.num_runs, p.outputSettings(false), d->outputSettings(false));

        // write CATS output
        b->writeOutput(p, d, t);

        // write CPLEX output if desired
        if (p.cplex_output)
            b->writeCPLEXFile(p, t);

        //    printf("Num components in bid graph: %d\n", b->numComponents());

        delete d;
        delete b;
    }

    return 0;
}
