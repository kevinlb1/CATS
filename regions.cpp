// TODO: 
// - read in rows and cols as a parameter

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "normal.h"
#include "bid.h"
#include "BidSet.h"
#include "regions.h"

const int Regions::Z_NUMGOODS = 500;
const int Regions::item::arrayIncSize = 5;

bool Regions::IsConnectedLocations(int loc1, int loc2)
{
    int i;

    for (i = 0; i < location[loc1].numneigh; i++) {
        if (location[loc1].neighbor[i] == loc2) {
            return true;
        }
    }

    return false;
}

void Regions::ConnectLocations(int loc1, int loc2)
{
    if (IsConnectedLocations(loc1, loc2)) {
        return;
    }

    location[loc1].addNeighbor(loc2);
    location[loc2].addNeighbor(loc1);
}

void Regions::Connect(int i, int j, int x, int y)
{
    if ((i >= 0) && (i < nrows) &&
        (x >= 0) && (x < nrows) &&
        (j >= 0) && (j < ncols) &&
        (y >= 0) && (y < ncols)) {
        ConnectLocations(i * ncols + j, x * ncols + y);
    }
}

bool Regions::IsConnected(int i, int j, int x, int y)
{
    if ((i >= 0) && (i < nrows) &&
        (x >= 0) && (x < nrows) &&
        (j >= 0) && (j < ncols) &&
        (y >= 0) && (y < ncols)) {
        return IsConnectedLocations(i * ncols + j, x * ncols + y);
    }
    else {
        return false;
    }
}

// NUMBER OF ROWS AND COLS SHOULD BE READ IN DIRECTLY RATHER THAN NUMBER OF GOODS
// note: this function modifies the number of goods as floor(sqrt(num_goods))**2
// build the real estate map
int Regions::BuildMap()
{
    int i, j, tmp;
    nrows = (int) floor(sqrt(num_goods));
    ncols = (int) floor(sqrt(num_goods));
    numgoods = nrows * ncols;

    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            if (Param::DRand() > remove_neighbor) {
                Connect(i, j, i + 1, j);
            }

            if (Param::DRand() > remove_neighbor) {
                Connect(i, j, i, j + 1);
            }

            if (Param::DRand() < additional_neighbor) {
                tmp = Param::LRand() % 2;

                switch (tmp) {
                case 0:
                    if (!IsConnected(i, j, i + 1, j + 1)) Connect(i, j, i + 1, j + 1);
                    break;

                case 1:
                   if (!IsConnected(i + 1, j, i, j + 1)) Connect(i + 1, j, i, j + 1);
                    break;
                }
            }
        }
    }

    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            if (location[i * ncols + j].numneigh == 0) {
                return 1;
            }
        }
    }

    return 0;
}

void Regions::DisplayMap(char * edgefile, char * cityfile)
{
    int i, j;

    FILE * fp;
    fp = fopen(edgefile, "w");

    if (fp == NULL) {
        printf("error opening display file for writing\n");
        exit(2);
    }

    for (i = 0; i < num_goods; i++) {
        for (j = 0; j < location[i].numneigh; j++) {
            fprintf(fp, "%d %d\n%d %d\n\n", i / ncols, i % ncols,
                    location[i].neighbor[j] / ncols, location[i].neighbor[j] % ncols);
        }
    }

    fclose(fp);

    fp = fopen(cityfile, "w");

    if (fp == NULL) {
        printf("error opening display file for writing\n");
        exit(2);
    }

    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            fprintf(fp, "%d %d\n", i, j);
        }
    }

    fclose(fp);
}

// return the value of the bid.  This amount is private + common value
double Regions::Value(bid_type_for_generating &b)
{
    int i;
    double total = 0;

    for (i = 0; i < b.num_items; i++) {
        total += location[b.items[i]].p + location[b.items[i]].c;
    }

    total += S(b.num_items);
    return total;
}

// return only the common value of the bid
double Regions::CommonValue(bid_type_for_generating &b)
{
    int i;
    double total = 0;

    for (i = 0; i < b.num_items; i++) {
        total += location[b.items[i]].c;
    }

    return total;
}

int Regions::InBidAlready(bid_type_for_generating &b, int new_item)
{
    int i;

    for (i = 0; i < b.num_items; i++) {
        if (b.items[i] == new_item) {
            return 1;
        }
    }

    return 0;
}

// decide how much weight to add if location a and b are in the same bid
double Regions::weightFromLocations(int a, int b)
{
    if (IsConnectedLocations(a, b)) {
        return location[b].pn;
    }
    else {
        return 0.0;
    }
}

void Regions::Add_Good_to_Bundle(bid_type_for_generating * b)
{
    int i, j, new_item;
    double prob; //,s;

    if (b->num_items >= num_goods) {
        return;
    }

    // JUMP: pick a good uniformly at random
    if (Param::DRand() < jump_prob) {
        do {
            new_item = Param::LRand() % num_goods;
        }
        while (InBidAlready(*b, new_item));

        b->items[b->num_items] = new_item;
        b->num_items++;
    }
    // add a good based on the weights
    else {
        // build the distribution
        double *p = new double[num_goods];
        double sum = 0.0;

        for (i = 0; i < num_goods; i++) {
            if (InBidAlready(*b, i)) {
                p[i] = 0.0;
            }
            else {
                prob = 0.0;

                if (b->num_items > Z_NUMGOODS) {
                    printf("invalid num_items: increase Z_NUMGOODS\n");
                    exit(1);
                }

                for (j = 0; j < b->num_items; j++) {
                    prob += weightFromLocations(b->items[j], i); // different for Regions, Arbitrary
                }

                p[i] = prob;
                sum += prob;
            }
        }

        // check that we have something to add
        if (sum <= 0) {
            return;
        }

        // normalize
        for (i = 0; i < num_goods; i++) {
            p[i] /= sum;
        }

        // draw a value
        double rand_num = Param::DRand();
        sum = 0.0;

        for (new_item = 0; rand_num > p[new_item] + sum && new_item < num_goods; new_item++) {
            sum += p[new_item];
        }

        if (new_item >= num_goods) {
            for (new_item = num_goods - 1; location[new_item].pn == 0; new_item--); // handle rounding error
        }

        b->items[b->num_items] = new_item;
        b->num_items++;
        delete[] p;
    }
}

int Regions::bid_compare(const void * a, const void * b)
{
    bid_type_for_generating * bid_a = (bid_type_for_generating *) a;
    bid_type_for_generating * bid_b = (bid_type_for_generating *) b;

    if (bid_a->value > bid_b->value) {
        return -1;
    }
    if (bid_a->value < bid_b->value) {
        return 1;
    }

    return 0;
}

Bid *Regions::MakeGeneratedBidIntoRealBid(bid_type_for_generating *b)
{
    Bid *newb = new Bid(b->value);

    for (int i = 0; i < b->num_items; i++) {
        newb->addGood(b->items[i]);
    }

    return newb;
}

int Regions::Identical(bid_type_for_generating &a, bid_type_for_generating &b)
{
    int i, j, foundit;

    if (a.num_items != b.num_items) {
        return 0;
    }

    for (i = 0; i < a.num_items; i++) {
        foundit = 0;

        for (j = 0; (j < a.num_items) && !foundit; j++) {
            if (a.items[i] == b.items[j]) {
                foundit = 1;
            }
        }

        if (!foundit) {
            return 0;
        }
    }

    return 1;
}

// the main function that generates all the bids
BidSet *Regions::generate(Param &p)
{
    int i;

    for (i = 0; i < num_goods; i++) {
        location[i].numneigh = 0;
        //  location[i].neighbor=NULL;

        //  location[i].d = NULL;
    }

    // make the map
    int z;

    do {
        z = BuildMap();
        //  make sure that every good is connected to at least one thing
        //  for (i = 0; i < num_goods; i++)
        //  {
        //	    bool con = false;
        //
        //		for (int j = 0; j < num_goods; j++)
        //		{
        //			if (i != j && IsConnectedLocations(i, j))
        //			{
        //		    	con = true;
        //				break;
        //				}
        //			}
        //
        //			if (!con)
        //			{
        //				z = 1;
        //				break;
        //			}
        //		}
    }
    while (z > 0);

    /* DisplayMap("real_e", "real_c"); */
    int g, one_valid, numsubs, seenbefore, j;
    double total, budget, min_resale_value;
    double minP = normal_mean + 1000 * normal_stdev;

    BidSet *bids = new BidSet(num_goods, p.num_bids);
    bid_type_for_generating b, *bs;
    bs = new bid_type_for_generating[num_goods];

    for (i = 0; i < num_goods; i++) {
        location[i].c = Param::DRand() * (max_good_value - 1) + 1;
    }

    // make the bids
    for (int t = 0; bids->numBids() < p.num_bids; t++) // was a while loop without the "t"
    {
restart_generation:
        total = 0;
        b.c = 0;
        b.value = 0;
        b.num_items = 0;

        for (i = 0; i < num_goods; i++) {
            // make the private values for the bidder: either uniform or normal
            if (!normal_prices) {
                location[i].p = -deviation * max_good_value +
                        (Param::DRand() * 2 * deviation * max_good_value);

                location[i].pn = (location[i].p + deviation * max_good_value) /
                        2 * deviation * max_good_value;

                total += location[i].pn;
            }
            else {
                location[i].p = norm->draw(normal_mean, normal_stdev);

                total += location[i].p;

                if (location[i].p < minP) {
                    minP = location[i].p;
                }
            }
        }

        // normalize
        if (!normal_prices) {
            for (i = 0; i < num_goods; i++) {
                location[i].pn /= total;
            }
        }
        else {
            total -= num_goods * minP;

            for (i = 0; i < num_goods; i++) {
                location[i].pn = location[i].p - minP;
                location[i].pn /= total;
            }
        }

        // now make the actual bids for the bidder
        b.num_items = 0;

        // pick the first good for this bidder
        double sum = 0.0;
        double rand_num = Param::DRand();

        for (g = 0; rand_num > sum + location[g].pn && g < num_goods; g++) 
            // this relies on sum_i location[i].pn == 1 
        {
            sum += location[g].pn;
        }

        if (g >= num_goods) {
            for (g = num_goods - 1; location[g].pn == 0; g--); // handle rounding error
        }

        b.items[b.num_items] = g;
        b.num_items++;

        // add additional goods as required
        int i = 0;

        while (Param::DRand() <= additional_location) {
            Add_Good_to_Bundle(&b);
        }

        // calculate value of the bid
        b.c = CommonValue(b);
        b.value = Value(b);

        if (b.value <= 0) {
            goto restart_generation;
        }

        budget = budget_factor * b.value;
        min_resale_value = resale_factor * b.c;
        one_valid = 0;

        // make substitutible bids
        for (i = 0; i < b.num_items; i++) {
            bs[i].num_items = 1;
            bs[i].items[0] = b.items[i];

            while (bs[i].num_items < b.num_items) {
                Add_Good_to_Bundle(&bs[i]);
            }

            bs[i].value = Value(bs[i]);
            bs[i].c = CommonValue(bs[i]);

            if ((bs[i].value >= 0) && (bs[i].value <= budget) &&
                (bs[i].c >= min_resale_value) && !Identical(b, bs[i])) {
                one_valid = 1;
            }
        }

        // add the generated bid(s) to the bidset
        if (!one_valid) {
            bids->add(MakeGeneratedBidIntoRealBid(&b), p.remove_dominated_bids);
            /* printf("Only one valid.\n"); */
        }
        // if there's more than one bid, put them in order so that the most
        // competitive ones are selected.
        // make sure that multiple identical bids are not substitutes
        else {
            /* printf("Multiple valid.\n");*/
            bids->addXor(MakeGeneratedBidIntoRealBid(&b));
            numsubs = 0;
            qsort(bs, b.num_items, sizeof (bid_type_for_generating), bid_compare);

            for (i = 0; i < b.num_items; i++) {
                if (numsubs >= max_substitutible_bids) {
                    break;
                }

                if ((bs[i].value >= 0) && (bs[i].value <= budget) &&
                    (bs[i].c >= min_resale_value) && !Identical(b, bs[i])) {
                    seenbefore = 0;

                    for (j = 0; (j < i) && !seenbefore; j++) {
                        if (Identical(bs[j], bs[i])) {
                            seenbefore = 1;
                        }
                    }

                    if (!seenbefore) {
                        bids->addXor(MakeGeneratedBidIntoRealBid(&(bs[i])));
                        numsubs++;
                        //  t++;
                    }
                }
            }

            bids->doneAddingXor(p.remove_dominated_bids);
        }

        // output progress
        if (p.output_parameter_settings && (!(t % p.output_frequency) || bids->numBids() >= p.num_bids)) {
            int count = printf("Number of %sbids: %d  Number of bidders: %d",
                               (p.remove_dominated_bids ? "non-dominated " : ""), bids->numBids(), t);

            if (bids->numBids() < p.num_bids) {
                for (int counter = 0; counter < count; counter++) {
                    printf("\b");
                }
            }
            else {
                printf("\n");
            }
        }
    }

    delete[] bs;
    return bids;
}

// output information on how to use the distribution
void Regions::usage(bool arbitrary)
{
    int count = printf("Usage for %s Distribution:\n", arbitrary ? "Arbitrary" : "Regions");

    for (int t = 0; t < count; t++) {
        printf("=");
    }

    printf(" -goodvl : sets value of max good value.\n");
    printf(" -maxbid : sets value of max substitutable bids.\n");
    printf(" -remnpr : sets value of removing hv neighbor.\n");
    printf(" -addnpr : sets value of adding diag neighbor.\n");
    printf(" -addloc : sets value of additional location.\n");

    if (!arbitrary) {
        printf(" -jump   : sets value of jump prob.\n");
    }

    printf(" -addit  : sets value of additivity.\n");
    printf(" -dev    : sets value of deviation.\n");
    printf(" -mean   : sets value of normal_mean.\n");
    printf(" -stdev  : sets value of normal_stdev.\n");
    printf(" -budget : sets value of budget_factor.\n");
    printf(" -resale : sets value of resale_factor.\n");

    printf("\n To specify uniform or normal price distributions, add \"-upv\" or \"-npv\" to the distribution name.\n");
    //  printf(" -run    : produce an output file and set the threshold.\n");
}

// constructor: parse input
Regions::Regions(Param &p, char * argv[], int argc)
{
    output_buffer = NULL;
    int i;

    // set initial values
    max_good_value = 100;
    max_substitutible_bids = 5;
    remove_neighbor = 0.1;
    additional_neighbor = 0.2;
    additional_location = 0.9;
    jump_prob = 0.05;
    additivity = 0.2;
    deviation = 0.5;
    budget_factor = 1.5;
    resale_factor = 0.5;
    normal_mean = 0.0;
    normal_stdev = 30.0;
    num_goods = p.num_goods;

    // effectiveDist is always -npv or -upv, even
    // if the real distribution does not specify
    effectiveDist = p.distribution;

    if (p.distribution == Param::REGIONS ||
        p.distribution == Param::ARBITRARY) {
        if (Param::LRand(0, 1) == 1) {
            normal_prices = true;
            effectiveDist += 1;
        }
        else {
            normal_prices = false;
            effectiveDist += 2;
        }
    }
    else if (p.distribution == Param::REGIONS_NPV ||
             p.distribution == Param::ARBITRARY_NPV) {
        normal_prices = true;
    }
    else {
        assert(p.distribution == Param::REGIONS_UPV ||
               p.distribution == Param::ARBITRARY_UPV);
        normal_prices = false;
    }

    if (p.distribution == Param::REGIONS ||
        p.distribution == Param::REGIONS_NPV ||
        p.distribution == Param::REGIONS_UPV) {
        arbitrary = false;
    }
    else {
        arbitrary = true;
    }

    // allocation
    location = new item[p.num_goods];
    norm = new Normal(p);

    // parse parameters, ignore those that don't apply
    for (i = 1; i < argc; i++) {
        // now read in values
        if (argv[i][0] == '-' || argv[i][0] == '/') {
            char *in = &(argv[i][1]); // skip the first character

            if (!stricmp("goodvl", in)) {
                p.argRecognized[i] = true;
                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                max_good_value = atoi(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("maxbid", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                max_substitutible_bids = atoi(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("remnpr", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                remove_neighbor = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("addnpr", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                additional_neighbor = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("addloc", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                additional_location = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("jump", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                jump_prob = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("addit", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                additivity = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("dev", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                deviation = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("mean", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                normal_mean = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("stdev", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                normal_stdev = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("budget", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                budget_factor = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("resale", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage(arbitrary);
                }

                resale_factor = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
        }
    }
}

Regions::~Regions()
{
    delete[] location;
    delete norm;

    if (output_buffer) {
        delete[] output_buffer;
    }
}

// output values of all variables related to the distribution
char *Regions::outputSettings(bool tofile)
{
    // allocate variables
    char comment[3];

    if (tofile) {
        comment[0] = '%';
        comment[1] = ' ';
        comment[2] = 0;
    }
    else {
        comment[0] = 0;
    }

    char buffer1[80], buffer2[80], buffer3[80], buffer4[80], buffer5[80];
    char buffer6[80], buffer7[80], buffer8[80], buffer9[80], buffer10[80];
    char buffer11[80], buffer12[80]; //,buffer13[80],buffer14[80],buffer15[80];

    // generate output
    sprintf(buffer1, "%s%s Distribution Parameters:\n", comment, arbitrary ? "Arbitrary" : "Regions");
    sprintf(buffer2, "%sMax Substitutable Bids = %d\n", comment, max_substitutible_bids);
    sprintf(buffer3, "%sRemoving HV Neighbor   = %g\n", comment, remove_neighbor);
    sprintf(buffer4, "%sAdding DIAG Neighbor   = %g\n", comment, additional_neighbor);
    sprintf(buffer5, "%sAdditional Location    = %g\n", comment, additional_location);

    if (!arbitrary) {
        sprintf(buffer6, "%sJump Probability       = %g\n", comment, jump_prob);
    }
    else {
        buffer6[0] = 0;
    }

    sprintf(buffer7, "%sAdditivity             = %g\n", comment, additivity);
    sprintf(buffer8, "%sBudget Factor          = %g\n", comment, budget_factor);
    sprintf(buffer9, "%sResale Factor          = %g\n", comment, resale_factor);

    if (!normal_prices) {
        sprintf(buffer10, "%sUniform Private Value:\n", comment);
        sprintf(buffer11, "%s  Deviation            = %g\n", comment, deviation);
        sprintf(buffer12, "%s  Max Good Value       = %d\n", comment, max_good_value);
    }
    else {
        sprintf(buffer10, "%sNormal Private Value:\n", comment);
        sprintf(buffer11, "%s  Mean                 = %g\n", comment, normal_mean);
        sprintf(buffer12, "%s  Standard Deviation   = %g\n", comment, normal_stdev);
    }

    // prepare the string and return it
    int length = strlen(buffer1) + strlen(buffer2) + strlen(buffer3) + strlen(buffer4) 
        + strlen(buffer5) + strlen(buffer6) + strlen(buffer7) + strlen(buffer8) + strlen(buffer9) 
        + strlen(buffer10) + strlen(buffer11) + strlen(buffer12)
        /*+strlen(buffer13)+strlen(buffer14)+strlen(buffer15)*/ + 20;

    if (output_buffer) {
        delete[] output_buffer;
    }

    output_buffer = new char[length];
    sprintf(output_buffer, "%s%s%s%s%s%s%s%s%s%s%s%s\n", buffer1, buffer2, buffer3, buffer4, 
            buffer5, buffer6, buffer7, buffer8, buffer9, buffer10, buffer11, 
            buffer12/*,buffer13,buffer14,buffer15*/);

    return output_buffer;
}

void Regions::randomizeParams(Param &p)
{
    if (p.param_dists[effectiveDist] != NULL) {
        double params[7];

        // set these originally for regions, uniform value dist
        // override for ARBITRARY, and/or normal value dist
        double param_mins[] = {1, 0, 0, 0.5, 0, -0.2, 0};
        double param_maxs[] = {10, 0.5, 1, 0.99, 0.2, 1, 1};

        // override for normal value dist
        if (normal_prices) {
            param_mins[6] = 5;
            param_maxs[6] = 50;
        }

        if (p.distribution == Param::ARBITRARY ||
            p.distribution == Param::ARBITRARY_NPV ||
            p.distribution == Param::ARBITRARY_UPV) {
            // no "jump_prob" in arbitrary
            param_mins[4] = param_mins[5];
            param_mins[5] = param_mins[6];
            param_maxs[4] = param_maxs[5];
            param_maxs[5] = param_maxs[6];

            randomizeParamsWithDistribution(p.param_dists[(int) effectiveDist], 6, 
                    param_mins, param_maxs, params, p);
        }
        else {
            randomizeParamsWithDistribution(p.param_dists[(int) effectiveDist], 7, 
                    param_mins, param_maxs, params, p);
        }

        max_substitutible_bids = (int) floor(params[0] + 0.5);
        remove_neighbor = params[1];
        additional_neighbor = params[2];
        additional_location = params[3];

        if (p.distribution == Param::ARBITRARY ||
            p.distribution == Param::ARBITRARY_NPV ||
            p.distribution == Param::ARBITRARY_UPV) {
            additivity = params[4];

            if (p.distribution == Param::REGIONS_UPV) {
                deviation = params[5];
            }
            else {
                normal_stdev = params[5];
            }
        }
        else {
            jump_prob = params[4];
            additivity = params[5];

            if (p.distribution == Param::REGIONS_UPV) {
                deviation = params[6];
            }
            else {
                normal_stdev = params[6];
            }
        }
    }
    else {
        max_substitutible_bids = Param::LRand(1, 10); // 5
        remove_neighbor = Param::DRand(0, 0.5); // 0.1
        additional_neighbor = Param::DRand(0, 1); // 0.2
        additional_location = Param::DRand(0.5, 0.99); // 0.9
        jump_prob = Param::DRand(0, 0.2); // 0.05
        additivity = Param::DRand(-0.2, 1); // 0.2
        //  budget_factor = Param::DRand(min,max);          // 1.5
        //	resale_factor = Param::DRand(min,max);          // 0.5
        deviation = Param::DRand(0, 1); // 0.5
        //	max_good_value = Param::DRand(min,max);         // 100
        //	normal_mean =  // leave this one alone          // 0
        normal_stdev = Param::DRand(5, 50); // 30
    }
}

