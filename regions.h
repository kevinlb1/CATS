// TODO: 
// - read in rows and cols as a parameter

#ifndef _REGIONS_H
#define _REGIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "normal.h"
#include "bid.h"
#include "BidSet.h"

class Regions : public Distribution {
public:
    // constructor: parse input
    Regions(Param &p, char * argv[], int argc);
    virtual ~Regions();

    // the main function that generates all the bids
    BidSet *generate(Param &p);
    void randomizeParams(Param &p);

    // output information on how to use the distribution
    static void usage(bool arbitrary);

    // output values of all variables related to the distribution
    char *outputSettings(bool tofile);
    static const int Z_NUMGOODS; // = 500

protected:
    int num_goods; // temp value that stores p->num_goods

    // parameters

    int max_good_value;
    int max_substitutible_bids;
    double remove_neighbor;
    double additional_neighbor;
    double additional_location;
    double jump_prob;
    double additivity;
    double deviation;
    double budget_factor;
    double resale_factor;
    double normal_mean;
    double normal_stdev;
    bool normal_prices;
    bool arbitrary;
    int effectiveDist;

    // local structs
    struct item {
        int numneigh;
        int arraySize;
        int * neighbor;
        double c; /* common value */
        double p; /* private value */
        double pn; /* intermediate calculations */
        double *d; /* array for use in "arbitrary" subclass; not used in "regions" */

        static const int arrayIncSize; // = 5

        item() {
            numneigh = 0;
            arraySize = 0;
            neighbor = NULL;
            d = NULL;
        }

        ~item() {
            if (neighbor) {
                delete[] neighbor;
            }
            if (d) {
                delete[] d;
            }
        }

        void addNeighbor(int n) {
            if (numneigh == arraySize) {
                int *newNeigh = new int[arraySize += arrayIncSize];

                if (neighbor) {
                    for (int i = 0; i < numneigh; i++) {
                        newNeigh[i] = neighbor[i];
                    }

                    delete[] neighbor;
                }

                neighbor = newNeigh;
            }

            neighbor[numneigh++] = n;
        }

    };

    struct bid_type_for_generating { //formerly "B"
        int *items;
        int num_items;
        double c; /* total common value for this bid */
        double value; /* value for this bid */

        bid_type_for_generating() {
            items = new int[Z_NUMGOODS];
        }

        ~bid_type_for_generating() {
            delete[] items;
        }
    };

    // locals
    item *location;
    int nrows, ncols;
    int numgoods;
    char *output_buffer;
    Normal *norm;

    // helper functions
    inline double S(double x) {
        return pow(x, 1 + additivity);
    }

    bool IsConnectedLocations(int loc1, int loc2);
    void ConnectLocations(int loc1, int loc2);
    void Connect(int i, int j, int x, int y);
    bool IsConnected(int i, int j, int x, int y);

    // NUMBER OF ROWS AND COLS SHOULD BE READ IN DIRECTLY RATHER THAN NUMBER OF GOODS
    // note: this function modifies the number of goods as floor(sqrt(num_goods))**2
    // build the real estate map

    virtual int BuildMap();
    void DisplayMap(char *edgefile, char *cityfile);

    // return the value of the bid. This amount is private + common value
    double Value(bid_type_for_generating &b);

    // return only the common value of the bid

    double CommonValue(bid_type_for_generating &b);

    int InBidAlready(bid_type_for_generating &b, int new_item);

    // decide how much weight to add if location a and b are in the same bid

    virtual double weightFromLocations(int a, int b);
    void Add_Good_to_Bundle(bid_type_for_generating *b);
    static int bid_compare(const void *a, const void *b);
    Bid *MakeGeneratedBidIntoRealBid(bid_type_for_generating *b);
    int Identical(bid_type_for_generating &a, bid_type_for_generating &b);
}; // Object Regions

#endif // _REGIONS_H

