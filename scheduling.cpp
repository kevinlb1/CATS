#include <stdio.h>

#include <stdlib.h>

#include <math.h>

#include <time.h>



#include "Distribution.h"

#include "Param.h"

#include "BidSet.h"

#include "bid.h"

#include "scheduling.h"



double Scheduling::sched_param_mins[] = {3, 0, 0, -0.2};

double Scheduling::sched_param_maxs[] = {0, 1.0, 1.0, 1.0};

// note sched_param_maxs[0], max_length is set in constructor

// because it depends on number of goods



// create a bid on the range start..end with value val

Bid *Scheduling::makeBid(int start, int end, double val)
{

    Bid *b = new Bid(val);

    for (int i = 0; i <= end - start; i++)

        b->addGood(start + i);

    return b;

}



// generate the bidset

BidSet *Scheduling::generate(Param &p)
{

    int l, d_1, cur_max_deadline, new_d, start;

    double dev;



    BidSet *bids = new BidSet(p.num_goods, p.num_bids);



    // make bids

    for (int t = 0; bids->numBids() < p.num_bids; t++) {



        // initialize

        l = 1 + (Param::LRand() % max_length);

        d_1 = l - 1 + Param::LRand() % (p.num_goods - l + 1);

        dev = 1 - deviation + 2 * Param::DRand() * deviation;

        cur_max_deadline = -1;

        new_d = d_1;



        do {

            // make a deadline

            for (start = 0; start < p.num_goods; start++) {

                if ((start + l - 1 <= new_d) && (cur_max_deadline < start + l - 1))

                    bids->addXor(makeBid(start, start + l - 1, dev * pow((double) l, 1 + additivity) * d_1 / new_d));

            }

            cur_max_deadline = new_d;

            if (cur_max_deadline == p.num_goods - 1) break;

            new_d = cur_max_deadline + 1 + Param::LRand() % (p.num_goods - cur_max_deadline - 1);

        }
        while (Param::DRand() <= prob_additional_deadline);



        // finish up the bids; refine

        bids->doneAddingXor(p.remove_dominated_bids);



        // output progress

        if (p.output_parameter_settings && (!(t % p.output_frequency) || bids->numBids() >= p.num_bids)) {

            int count = printf("Number of %sbids: %d  Number of bidders: %d",

                               (p.remove_dominated_bids ? "non-dominated " : ""), bids->numBids(), t);

            if (bids->numBids() < p.num_bids)

                for (int counter = 0; counter < count; counter++)

                    printf("\b");

            else

                printf("\n");

        }



    }



    return bids;

}



// display usage information

void Scheduling::usage()
{



    int count = printf("Usage for Scheduling Distribution:\n");

    for (int t = 0; t < count; t++) printf("=");



    printf(" -maxlength         : max number of timeslots required by a bidder\n");

    printf(" -deviation         : (uniform) deviation in value\n");

    printf(" -prob_add_deadline : probability that an additional deadline will be added\n");

    printf(" -additivity        : degree of super/subadditivity in valuations\n");

}



// constructor

Scheduling::Scheduling(Param &p, char *argv[], int argc)
{

    output_buffer = NULL;



    int i;



    // defaults

    max_length = 10;

    deviation = 0.5;

    prob_additional_deadline = 0.7;

    additivity = 0.2;



    // parse parameters. Ignore those that don't apply

    for (i = 1; i < argc; i++) {

        // now read in values

        if (argv[i][0] == '-' || argv[i][0] == '/') {



            if (!stricmp("-max_length", argv[i])) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                max_length = atoi(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }

            else if (!stricmp("-deviation", argv[i])) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                deviation = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }

            else if (!stricmp("-prob_add_deadline", argv[i])) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                prob_additional_deadline = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }

            else if (!stricmp("-additivity", argv[i])) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                additivity = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }

        }

    }



    // this must be set here, because it depends on num_goods

    sched_param_maxs[0] = (int) (0.5 * p.num_goods);



}

void Scheduling::randomizeParams(Param &p)
{

    if (p.param_dists[(int) Param::SCHEDULING] != NULL) {



        double params[4];



        randomizeParamsWithDistribution(p.param_dists[(int) Param::SCHEDULING], 4, sched_param_mins, sched_param_maxs, params, p);



        max_length = (int) params[0];

        deviation = params[1];

        prob_additional_deadline = params[2];

        additivity = params[3];



    }
    else {

        max_length = Param::LRand((int) sched_param_mins[0], (int) sched_param_maxs[0]); // 10

        deviation = Param::DRand(sched_param_mins[1], sched_param_maxs[1]); // 0.5

        prob_additional_deadline = Param::DRand(sched_param_mins[2], sched_param_maxs[2]); // 0.7

        additivity = Param::DRand(sched_param_mins[3], sched_param_maxs[3]); // 0.2

    }

}

Scheduling::~Scheduling()
{

    if (output_buffer)

        delete[] output_buffer;

}



// output values of all variables related to the distribution

char *Scheduling::outputSettings(bool tofile)
{

    // allocate variables

    char comment[3];

    if (tofile) {
        comment[0] = '%';
        comment[1] = ' ';
        comment[2] = 0;
    }

    else comment[0] = 0;

    char buffer1[80], buffer2[80], buffer3[80], buffer4[80], buffer5[80];



    // generate output

    sprintf(buffer1, "%sScheduling Distribution Parameters:\n", comment);

    sprintf(buffer2, "%smax_length                = %d\n", comment, max_length);

    sprintf(buffer3, "%sdeviation                 = %f\n", comment, deviation);

    sprintf(buffer4, "%sprob_additional_deadline  = %f\n", comment, prob_additional_deadline);

    sprintf(buffer5, "%sadditivity                = %f\n", comment, additivity);



    // prepare the string and return it

    int length = strlen(buffer1) + strlen(buffer2) + strlen(buffer3) + strlen(buffer4) +

            strlen(buffer5) + 20;



    if (output_buffer)

        delete[] output_buffer;

    output_buffer = new char[length];



    sprintf(output_buffer, "%s%s%s%s%s\n", buffer1, buffer2, buffer3, buffer4, buffer5);

    return output_buffer;

}

