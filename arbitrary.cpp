#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "Distribution.h"
#include "regions.h"
#include "arbitrary.h"

// build the arbitrary adjacency map
int Arbitrary::BuildMap()
{
    int i, j;
    //	numbids = 0;
    //	numdummy = 0;

    for (i = 0; i < num_goods; i++) {
        if (location[i].d == NULL) {
            location[i].d = new double[num_goods];
        }

        location[i].d[i] = 0;
        location[i].numneigh = 0; // NULL;  // only used in the "regions" superclass
        location[i].neighbor = NULL; // only used in the "regions" superclass
    }

    for (i = 0; i < num_goods; i++) {
        for (j = i + 1; j < num_goods; j++) {
            do {
                location[i].d[j] = Param::DRand();
                location[j].d[i] = location[i].d[j];
            }
            while (!location[i].d[j]); // make sure nothing has zero probability
        }
    }

    return 0;
}

// decide how much weight to add if location a and b are in the same bid
double Arbitrary::weightFromLocations(int a, int b)
{
    return location[b].pn * location[b].d[a];
}

// constructor
Arbitrary::Arbitrary(Param &p, char * argv[], int argc) : Regions(p, argv, argc)
{
    jump_prob = 0.0; // this is never supported under arbitrary, regardless of settings
}

// destructor
Arbitrary::~Arbitrary()
{
    // No special cleanup necessary: everything is taken care of in Regions
}

