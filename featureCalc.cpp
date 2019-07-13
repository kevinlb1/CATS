#include <iostream>
#include <fstream>
#include <set>
#include <cmath>
#include <algorithm>
#include <unistd.h>
#include <assert.h>

using namespace std;

#include "BidSet.h"
#include "featureCalc.h"

const char *FeatureCalc::featureNames[] = {"bg_abs_price_dev", "bg_ppg_dev", "bg_psqrg_dev",
    "bg_edge_dense", "bg_avg_deg", "bg_deg_dev",
    "bg_max_deg", "bg_min_deg", "bp_avg_good_deg",
    "bp_good_deg_dev", "bp_max_good_deg", "bp_min_good_deg",
    "bp_avg_bid_deg", "bp_bid_deg_dev", "bp_max_bid_deg",
    "bp_min_bid_deg", "bg_deg_q1", "bg_deg_q2", "bg_deg_q3",
    "lp_l1", "lp_avg", "lp_l2", "lp_l2_avg", "lp_linf",
    "cc_bg_radius", "cc_bg_diameter",
    "cc_bg_ecc_avg", "cc_bg_ecc_dev", "cc_clustering",
    "cc_path", "cc_cp_ratio", "cc_clust_dev"};

//const char *FeatureCalc::realismFeatureNames[] = {"
void FeatureCalc::calcFeatures(BidSet &bids)
{
    bids.makeBidGraph();
    num_bids = bids.numBids();
    num_goods = bids.numGoods();

    // -- bid graph properties
    int n = bids.numBids();
    long nEdges = 0;
    double sqEdges = 0;
    double sumPrices = 0, sumPSquared = 0;
    double sumPricesPG = 0, sumPPGSquared = 0;
    double sumPricesPSQRG = 0, sumPSQRGSquared = 0;
    int maxDegree = -1, minDegree = n + 1;
    int maxBidSize = -1, minBidSize = bids.numGoods() + bids.numDummyGoods() + 1;
    double avgBidSize = 0, sumSQBidSize = 0;

    int i;
    int* degArray = new int[n];

    for (i = 0; i < n; i++) {
        Bid* b = bids.getBid(i);

        // -- bid size stats
        int deg = b->num_goods;
        avgBidSize += deg;
        sumSQBidSize += deg*deg;
        maxBidSize = (deg > maxBidSize ? deg : maxBidSize);
        minBidSize = (deg < minBidSize ? deg : minBidSize);

        // -- prices
        sumPrices += b->amount;
        sumPSquared += (b->amount) * (b->amount);
        double ppg = (b->amount) / double(b->num_goods);
        sumPricesPG += ppg;
        sumPPGSquared += ppg*ppg;
        ppg = (b->amount) / sqrt(double(b->num_goods));
        sumPricesPSQRG += ppg;
        sumPSQRGSquared += ppg*ppg;

        // -- bid graph
        //deg = b->conflicts.size();
        deg = b->conflicts;
        degArray[i] = deg;
        nEdges += deg;
        sqEdges += deg * deg;
        maxDegree = (deg > maxDegree ? deg : maxDegree);
        minDegree = (deg < minDegree ? deg : minDegree);
    }

    // -- prices
    sumPrices /= double(n);
    sumPSquared /= double(n);
    sumPricesPG /= double(n);
    sumPPGSquared /= double(n);
    sumPricesPSQRG /= double(n);
    sumPSQRGSquared /= double(n);
    double priceDev = sqrt(sumPSquared - sumPrices * sumPrices);
    double ppgDev = sqrt(sumPPGSquared - sumPricesPG * sumPricesPG);
    double psqrgDev = sqrt(sumPSQRGSquared - sumPricesPSQRG * sumPricesPSQRG);

    // -- bid graph
    double connectivity = double(nEdges) / double(n);
    double conDev = sqrt(sqEdges / double(n) - connectivity * connectivity);
    nEdges /= 2;
    double density = nEdges / double( n * (n - 1) / 2);

    // -- bid size stats
    avgBidSize /= double(n);
    double bidSizeDev = sqrt(sumSQBidSize / double(n) - avgBidSize * avgBidSize);

    // -- calculate a few percentile statistics for the bid graph
    sort(degArray, degArray + n);
    int quotOne = degArray [n / 4 - 1];
    int median = degArray [n / 2];
    int quotThree = degArray [(3 * n) / 4 - 1];
    delete[] degArray;

    // -- bipartite grap properties
    n = bids.getGoodMap()->size();
    int maxGoodDegree = -1, minGoodDegree = bids.numBids() + 1;
    double avgGoodDeg = 0, sqGoodDeg = 0;

    for (i = 0; i < n; i++) {
        int deg = (*bids.getGoodMap())[i].size();
        avgGoodDeg += deg;
        sqGoodDeg += deg*deg;
        maxGoodDegree = (deg > maxGoodDegree ? deg : maxGoodDegree);
        minGoodDegree = (deg < minGoodDegree ? deg : minGoodDegree);
    }

    avgGoodDeg /= double(n);
    double goodDegDev = sqrt(sqGoodDeg / double(n) - avgGoodDeg * avgGoodDeg);

    // -- store results
    featureNum = 0;

    // -- prices
    featureVals[featureNum++] = priceDev;
    featureVals[featureNum++] = ppgDev;
    featureVals[featureNum++] = psqrgDev;

    // -- bid graph
    featureVals[featureNum++] = density;
    featureVals[featureNum++] = connectivity;
    featureVals[featureNum++] = conDev;
    featureVals[featureNum++] = maxDegree;
    featureVals[featureNum++] = minDegree;

    // -- Good nodes in bipartite
    featureVals[featureNum++] = avgGoodDeg;
    featureVals[featureNum++] = goodDegDev;
    featureVals[featureNum++] = maxGoodDegree;
    featureVals[featureNum++] = minGoodDegree;

    // -- bid nodes in bipartite
    featureVals[featureNum++] = avgBidSize;
    featureVals[featureNum++] = bidSizeDev;
    featureVals[featureNum++] = maxBidSize;
    featureVals[featureNum++] = minBidSize;

    // -- percentiles
    featureVals[featureNum++] = quotOne;
    featureVals[featureNum++] = median;
    featureVals[featureNum++] = quotThree;

    LPfeatures(bids);
    testSmallWorldBidGraph(bids.numBids(), bids.bidGraph, 1/*calc avg min path*/);

    assert(featureNum == numFeatures);

}

/* adapted/copied from testbidgraph.c */
float FeatureCalc::testSmallWorldBidGraph(int bids, int **bidGraph, int calcavgpath)
{
    int i, j, k, nn, neighbours, n1, n2;
    float avgpath, avgclustering, almost, stdevclust;
    int neighbour[6000], deg[6000];
    float clust[6000];
    long int eccentricity[6000];
    long int diameter = 0;
    long int radius = 0;
    double eccaverage = 0, eccdev = 0;

    /* Calculate the average clustering index of the graph */
    avgclustering = 0.0;

    for (i = 0; i < bids; i++) {
        neighbours = 0;

        for (j = 0; j < bids; j++) {
            if (bidGraph[i][j]) {
                neighbour[neighbours++] = j;
            }
        }

        deg[i] = neighbours;

        if (neighbours > 1) {
            nn = 0;
            almost = (neighbours * (neighbours - 1)) / 2;

            for (n1 = 0; n1 < (neighbours - 1); n1++) {
                for (n2 = n1 + 1; n2 < neighbours; n2++) {
                    if (bidGraph[neighbour[n1]][neighbour[n2]]) nn++;
                }
            }

            clust[i] = (((float) nn) / almost);
            avgclustering += clust[i];
        }
        else {
            if (neighbours == 1) {
                clust[i] = 1.0;
                avgclustering += 1.0;
            }
            else clust[i] = 0.0;
        }
    }

    avgclustering /= ((float) bids);
    stdevclust = 0.0;

    for (i = 0; i < bids; i++) {
        stdevclust += ((clust[i] - avgclustering)*(clust[i] - avgclustering));
    }

    stdevclust = sqrt(stdevclust / ((float) bids));

    for (i = 0; i < bids; i++) {
        eccentricity[i] = 0;

        for (j = 0; j < bids; j++) {
            if (i != j) {
                if (bidGraph[i][j] == 0) bidGraph[i][j] = bids * bids;
            }
        }
    }

    if (calcavgpath) {
        //  printf( "calculating shortest paths\n" );
        avgpath = 0.0;

        //  floyd walsh all pairs shortest path algorithm to compute shortest distances
        for (k = 0; k < (bids - 1); k++) {
            for (i = 0; i < bids; i++) {
                for (j = i + 1; j < bids; j++) {
                    bidGraph[i][j] = min(bidGraph[i][j], bidGraph[i][k] + bidGraph[k][j]);
                    bidGraph[j][i] = bidGraph[i][j];

                    if ((k == (bids - 2)) && (i < j)) {
                        avgpath += bidGraph[i][j];

                        // == ^^ Should be inside the if, outside for hist. reasons
                        if (bidGraph[i][j] != bids * bids) {
                            if (eccentricity[i] < bidGraph[i][j]) {
                                eccentricity[i] = bidGraph[i][j];
                            }

                            if (eccentricity[j] < bidGraph[i][j]) {
                                eccentricity[j] = bidGraph[i][j];
                            }
                        }
                    }
                } // for j

            } // for i

        } //for k

        avgpath /= ((float) (bids * (bids - 1) / 2));
        radius = bids * bids + 173;
        diameter = 0;

        for (i = 0; i < bids; i++) {
            if (eccentricity[i] < radius) {
                radius = eccentricity[i];
            }

            if (eccentricity[i] > diameter) {
                diameter = eccentricity[i];
            }

            eccaverage += eccentricity[i];
            eccdev += eccentricity[i] * eccentricity[i];
        }

        eccaverage /= (double) bids;
        eccdev = sqrt(eccdev / (double) bids - eccaverage * eccaverage);
    }
    else {
        avgpath = 1.0;
    }

    /*
    fprintf( fp, "%ld, %ld, %lf, %lf, %f,  %f, %f, %f, %f, %f, ",
             radius, diameter, eccaverage, eccdev, avgclustering, avgpath, 
             avgclustering/avgpath, avgdeg, stdevclust, stdevdeg  );
     */

    featureVals[featureNum++] = radius;
    featureVals[featureNum++] = diameter;
    featureVals[featureNum++] = eccaverage;
    featureVals[featureNum++] = eccdev;
    featureVals[featureNum++] = avgclustering;
    featureVals[featureNum++] = avgpath;
    featureVals[featureNum++] = (avgclustering / avgpath);
    featureVals[featureNum++] = stdevclust;

    return (avgclustering / avgpath);
}

FeatureCalc::LPData::LPData()
{
#ifdef USE_CPLEX
    env = NULL;
    lp = NULL;
#endif
    status = 0;
}

void FeatureCalc::LPData::allocate(BidSet *b)
{
#ifdef USE_CPLEX
    int bids = b->numBids();
    int goods = b->numGoods() + b->numDummyGoods();
    int goodsMentioned = 0;

    for (int i = 0; i < b->numBids(); i++) {
        goodsMentioned += b->getBid(i)->num_goods;
    }

    obj = new double[bids];
    lb = new double[bids];
    ub = new double[bids];
    xvals = new double[bids];

    int rowsInLPMat = goods;
    int nonzeroInLPMat = goodsMentioned;

    rhs = new double[rowsInLPMat];
    sense = new char[rowsInLPMat];
    rmatbeg = new int[rowsInLPMat];
    rmatind = new int[nonzeroInLPMat];
    rmatval = new double[nonzeroInLPMat];
#else
    xvals = new double [b->numBids()];
#endif
}

void FeatureCalc::LPData::deallocate()
{
#ifdef USE_CPLEX
    delete[] lb;
    delete[] ub;
    delete[] rhs;
    delete[] sense;
    delete[] rmatbeg;
    delete[] rmatind;
    delete[] rmatval;
    delete[] obj;
#endif
    delete[] xvals;
}

void FeatureCalc::LPfeatures(BidSet& b)
{
    // BidSet b = *bids;
    lpd.allocate(&b);
    LPsolve(&b);

    double sum = 0;
    double l2 = 0;
    double max_dist = -1;
    double *dists = new double[b.numBids()];

    for (int i = 0; i < b.numBids(); i++) {
        double dist = min(abs(lpd.xvals[i]), abs(lpd.xvals[i] - 1));
        sum += dist;
        l2 += dist*dist;
        max_dist = max(max_dist, dist);
        dists[i] = dist;
    }

    sort(dists, dists + b.numBids());
    delete[] dists;

    featureVals[featureNum++] = sum;
    featureVals[featureNum++] = sum / double(b.numBids());
    featureVals[featureNum++] = sqrt(l2);
    featureVals[featureNum++] = sqrt(l2 / double(b.numBids()));
    featureVals[featureNum++] = max_dist;
    lpd.deallocate();

    return;
}

#ifdef USE_CPLEX

// == LP Features using cplex
double FeatureCalc::LPsolve(BidSet *b)
{
    /* Initialize the CPLEX environment */
    lpd.env = CPXopenCPLEX(&lpd.status);

    while (lpd.env == NULL) {
        fprintf(stderr, "Could not open CPLEX environment. Trying again in 3 sec.\n");
        sleep(3);
        lpd.env = CPXopenCPLEX(&lpd.status);
    }

    int bids = b->numBids();
    int goods = b->numGoods() + b->numDummyGoods();
    int goodsMentioned = 0;
    int solstat;
    double objval;

    for (int i = 0; i < b->numBids(); i++) {
        goodsMentioned += b->getBid(i)->num_goods;
    }

    int rowsInLPMat = goods;
    int colsInLPMat = bids;
    int nonzeroInLPMat = goodsMentioned;
    int currIndex = 0;

    for (int i = 0; i < bids; i++) {
        lpd.obj[i] = b->getBid(i)->amount;
        lpd.lb[i] = 0;
        lpd.ub[i] = 1;
    }

    /*lpd.status = CPXsetintparam (lpd.env, CPX_PARAM_SCRIND, CPX_ON);
    if (lpd.status) {
       fprintf (stderr,
                "Failure to turn on screen indicator, error %d.\n", lpd.status);
       exit(1);
    }*/

    /* Turn on data checking */
    /*
    lpd.status = CPXsetintparam (lpd.env, CPX_PARAM_DATACHECK, CPX_ON);
    if (lpd.status) {
       fprintf (stderr,
                "Failure to turn on data checking, error %d.\n", lpd.status);
       exit(1);
    } */

    lpd.lp = CPXcreateprob(lpd.env, &lpd.status, "lpex1");

    if (lpd.lp == NULL) {
        fprintf(stderr, "Failed to create LP.\n");
        CPXcloseCPLEX(&lpd.env);
        return -1;
    }

    CPXchgobjsen(lpd.env, lpd.lp, CPX_MAX); /* Problem is maximization */
    lpd.status = CPXnewcols(lpd.env, lpd.lp, colsInLPMat, lpd.obj, lpd.lb, lpd.ub, NULL, NULL);
    int goodCount = 0;

    for (int j = 0; j < goods; j++) {
        lpd.rmatbeg[goodCount] = currIndex;

        for (int k = 0; k < bids; k++) {
            if (b->getBid(k)->contains(j)) {
                lpd.rmatind[currIndex] = k;
                lpd.rmatval[currIndex++] = 1;
            }
        }

        lpd.sense[goodCount] = 'L';
        lpd.rhs[goodCount] = 1;
        goodCount++;
    }

    lpd.status = CPXaddrows(lpd.env, lpd.lp, 0, rowsInLPMat, nonzeroInLPMat, 
            lpd.rhs, lpd.sense, lpd.rmatbeg,
            lpd.rmatind, lpd.rmatval, NULL, NULL);

    /*
     lpd.status = CPXwriteprob (lpd.env, lpd.lp, "try1.lp", NULL);
     if (lpd.status) {
        fprintf (stderr, "Failed to write LP to disk.\n");
        exit(1);
     }*/

    lpd.status = CPXlpopt(lpd.env, lpd.lp);
    if (lpd.status) {
        fprintf(stderr, "Failed to optimize LP.\n");
    }

    lpd.status = CPXsolution(lpd.env, lpd.lp, &solstat, &objval, lpd.xvals, NULL, NULL, NULL);
    if (lpd.lp != NULL) {
        lpd.status = CPXfreeprob(lpd.env, &lpd.lp);

        if (lpd.status) {
            fprintf(stderr, "CPXfreeprob failed, error code %d.\n", lpd.status);
        }
    }

    lpd.lp = NULL;
    //  printf("lp solve = %f\n", objval);
    CPXcloseCPLEX(&lpd.env);

    return objval;
}

#else
// == LP Features using LP_SOLVE
double FeatureCalc::LPsolve(BidSet *b)
{
    int bids = b->numBids();
    int goods = b->numGoods() + b->numDummyGoods();
    lprec* lp = make_lp(0, bids);

    if (lp == NULL) {
        fprintf(stderr, "ERROR: Failed to init lp_solve\n");
        lpd.status = -1;
        return -1;
    }

    REAL obj[bids + 1];
    REAL cons[bids + 1];

    // -- set up bounds
    for (int i = 1; i <= bids; i++) {
        obj[i] = b->getBid(i - 1)->amount;
        set_upbo(lp, i, 1);
        set_lowbo(lp, i, 0);
    }

    // -- add constraints
    for (int j = 0; j < goods; j++) {
        for (int i = 0; i <= bids; i++) {
            cons[i] = 0;
        }

        REAL rhs = 1.0;

        for (int k = 0; k < bids; k++) {
            if (b->getBid(k)->contains(j)) {
                cons[k + 1] = 1;
            }
        }

        if (add_constraint(lp, cons, LE, rhs) != TRUE) {
            fprintf(stderr, "ERROR: add_constraint\n");
            delete_lp(lp);
            lpd.status = -1;
            return -1;
        }
    }

    // -- objective
    set_obj_fn(lp, obj);
    set_maxim(lp);

    // -- solve
    int res = solve(lp);

    // -- get results
    if (res != OPTIMAL) {
        fprintf(stderr, "ERROR in lp solve!\n");
        lpd.status = -1;
    }

    REAL objval;

    if (res == OPTIMAL) {
        objval = get_objective(lp);
    }
    else {
        objval = -1;
    }

    if (res == OPTIMAL) {
        if (!get_variables(lp, lpd.xvals)) {
            fprintf(stderr, "ERROR in LP comp!\n");
            lpd.status = -1;
        }
    }

    delete_lp(lp);
    return objval;
}
#endif

void FeatureCalc::writeFeatures(const char *filename)
{
    ofstream outfile(filename);

    if (!outfile) {
        printf("Error: could not open file \"%s\".\n", filename);
        exit(1);
    }

    outfile << num_bids << ',' << num_goods << ',';

    for (int i = 0; i < numFeatures; i++) {
        outfile << featureVals[i] << (i < numFeatures - 1 ? "," : "");
    }

    outfile << endl;
    outfile.close();
}

void FeatureCalc::writeFeatNames(const char *filename)
{
    ofstream outfile(filename);

    if (!outfile) {
        printf("Error: could not open file \"%s\".\n", filename);
        exit(1);
    }

    outfile << "bids" << ',' << "goods" << ',';

    for (int i = 0; i < numFeatures; i++) {
        outfile << featureNames[i] << (i < numFeatures - 1 ? "," : "");
    }

    outfile << endl;
    outfile.close();
}

