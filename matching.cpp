#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "bid.h"
#include "BidSet.h"
#include "Distribution.h"
#include "matching.h"

/* * Constants & Parameters **************************************** */
//#define NUMGOODS 20
//#define NUMBIDS 100
//#define m_max_airport_price 5
//#define m_deviation 0.5
//#define m_early_takeoff_deviation 1
//#define m_late_takeoff_deviation 2
//#define m_early_land_deviation 1
//#define m_late_land_deviation 2
//#define m_delay_coeff 0.9
//#define m_amount_late_coeff 0.75

const int Matching::MAX_LONGEST_FLIGHT_LENGTH = 10;
const int Matching::Z_MAX_DEV = 10;
const int Matching::NUMCITIES = 4;

// initialize the cities.
// regardless of the number of goods, there are only 4 cities, wich are the real
// ones used in the FAA application

void Matching::SetupCities()
{
    int i, j;

    cities[0].x = -87.75;
    cities[0].y = 41.98333333;

    cities[1].x = -77.03333333;
    cities[1].y = 38.85;

    cities[2].x = -73.783333;
    cities[2].y = 40.65;

    cities[3].x = -73.8666666;
    cities[3].y = 40.76666666;

    for (i = 0; i < NUMCITIES; i++) {
        for (j = 0; j < NUMCITIES; j++) {
            edges[i][j] = -1;
        }
    }

    for (i = 0; i < NUMCITIES; i++) {
        for (j = 0; j < NUMCITIES; j++) {
            if (i != j) edges[i][j] = RealDistance(i, j);
        }
    }

    edges[2][3] = -1;
    edges[3][2] = -1;

    max_l = -1;

    for (i = 0; i < NUMCITIES; i++) {
        for (j = i + 1; j < NUMCITIES; j++) {
            max_l = max(max_l, edges[i][j]);
        }
    }
}

// generate the bidset
BidSet *Matching::generate(Param &p)
{
    // variables
    int city1, city2, min_flight_length, start_time,
            takeoff, land, amount_late, delay;
    double l, dev;

    BidSet *bids = new BidSet(p.num_goods, p.num_bids);

    // initialize
    SetupCities();

    if (p.num_goods % NUMCITIES) {
        printf("The MATCHING distribution only works with num_goods %% %d == 0.\n", NUMCITIES);
        exit(1);
    }

    num_times = (int) floor(p.num_goods / NUMCITIES);

    for (int i = 0; i < NUMCITIES; i++) {
        cost[i] = Param::DRand() * m_max_airport_price;
    }

    // make bids
    for (int t = 0; bids->numBids() < p.num_bids; t++) {
        city1 = Param::LRand() % NUMCITIES;

        do {
            city2 = Param::LRand() % NUMCITIES;
        }
        while (!((city1 != city2) && (edges[city1][city2] > 0)));

        l = RealDistance(city1, city2);
        min_flight_length = Param::Round(longestFlightLength(num_times) * l / max_l);

        if (num_times == min_flight_length) {
            start_time = 1;
        }
        else {
            start_time = 1 + Param::LRand() % (num_times - min_flight_length);
        }

        dev = 1 - m_deviation + 2 * Param::DRand() * m_deviation;

        // takeoff time
        for (takeoff = max(1, start_time - m_early_takeoff_deviation);
                takeoff <= min(num_times, start_time + m_late_takeoff_deviation);
                takeoff++) {
            // land time
            for (land = takeoff + min_flight_length;
                    land <= min(start_time + min_flight_length + m_late_land_deviation, num_times);
                    land++) {
                // make one XOR bid for each combination
                amount_late = min(land - (start_time + min_flight_length), 0);
                delay = max(land - takeoff - min_flight_length, 0);

                Bid *b = new Bid(dev * (cost[city1] + cost[city2]) *
                                 pow(m_delay_coeff, delay) * pow(m_amount_late_coeff, amount_late));

                b->addGood((takeoff - 1) * NUMCITIES + city1);
                b->addGood((land - 1) * NUMCITIES + city2);

                bids->addXor(b);
            }
        }

        // that's it for this bidder
        bids->doneAddingXor(p.remove_dominated_bids);

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

    return bids;
}

// output command-line usage
void Matching::usage()
{
    int count = printf("Usage for Matching Distribution:\n");

    for (int t = 0; t < count; t++) {
        printf("=");
    }

    printf(" -price  : sets value of MAX AIRPORT PRICE.\n");
    printf(" -etkdev : sets value of EARLY TAKEOFF DEVIATION.\n");
    printf(" -ltkdev : sets value of LATE TAKEOFF DEVIATION.\n");
    printf(" -elndev : sets value of EARLY LANDING DEVIATION.\n");
    printf(" -llndev : sets value of LATE LANDING DEVIATION.\n");
    printf(" -dev    : sets value of DEVIATION.\n");
    printf(" -delay  : sets value of DELAY COEFF.\n");
    printf(" -lateco : sets value of AMOUNT LATE COEFF.\n");
}

// constructor
Matching::Matching(Param &p, char * argv [], int argc)
{
    output_buffer = NULL;
    int i;

    // create arrays
    cities = new point[NUMCITIES];
    edges = new double*[NUMCITIES];

    for (i = 0; i < NUMCITIES; i++) {
        edges[i] = new double[NUMCITIES];
    }

    cost = new double[NUMCITIES];

    // initialize
    m_max_airport_price = 5;
    m_early_takeoff_deviation = 1;
    m_late_takeoff_deviation = 2;
    m_early_land_deviation = 1;
    m_late_land_deviation = 2;
    m_deviation = 0.5;
    m_delay_coeff = 0.9;
    m_amount_late_coeff = 0.75;

    // parse parameters. Ignore those that don't apply
    for (i = 1; i < argc; i++) {
        // now read in values
        if (argv[i][0] == '-' || argv[i][0] == '/') {
            char *in = &(argv[i][1]); // skip the first character
            if (!stricmp("price", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage();
                }

                m_max_airport_price = atoi(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("etkdev", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage();
                }

                m_early_takeoff_deviation = atoi(argv[i + 1]);
                if (m_early_takeoff_deviation > Z_MAX_DEV) {
                    m_early_takeoff_deviation = Z_MAX_DEV;
                }

                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("ltkdev", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage();
                }

                m_late_takeoff_deviation = atoi(argv[i + 1]);
                if (m_late_takeoff_deviation > Z_MAX_DEV) {
                    m_late_takeoff_deviation = Z_MAX_DEV;
                }

                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("elndev", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage();
                }

                m_early_land_deviation = atoi(argv[i + 1]);
                if (m_early_land_deviation > Z_MAX_DEV) {
                    m_early_land_deviation = Z_MAX_DEV;
                }

                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("llndev", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage();
                }

                m_late_land_deviation = atoi(argv[i + 1]);
                if (m_late_land_deviation > Z_MAX_DEV) {
                    m_late_land_deviation = Z_MAX_DEV;
                }

                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("dev", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage();
                }

                m_deviation = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("delay", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage();
                }

                m_delay_coeff = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
            else if (!stricmp("lateco", in)) {
                p.argRecognized[i] = true;

                if (i + 1 >= argc) {
                    usage();
                }

                m_amount_late_coeff = atof(argv[i + 1]);
                i++;
                p.argRecognized[i] = true;
            }
        }
    }

    if (m_late_takeoff_deviation < m_early_takeoff_deviation) {
        m_late_takeoff_deviation = m_early_takeoff_deviation;
    }

    if (m_late_land_deviation < m_early_land_deviation) {
        m_late_land_deviation = m_early_land_deviation;
    }
}

// destructor
Matching::~Matching()
{
    int i;

    // delete arrays
    delete[] cities;

    for (i = 0; i < NUMCITIES; i++) {
        delete[] edges[i];
    }

    delete[] edges;
    delete[] cost;

    if (output_buffer) {
        delete[] output_buffer;
    }
}

// output distribution parameters, either to the screen or to a file
char *Matching::outputSettings(bool tofile)
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
    char buffer6[80], buffer7[80], buffer8[80], buffer9[80];

    // generate output
    sprintf(buffer1, "%sMatching Distribution Parameters:\n", comment);
    sprintf(buffer2, "%sMax airport price       = %d\n", comment, m_max_airport_price);
    sprintf(buffer3, "%sEarly takeoff deviation = %d\n", comment, m_early_takeoff_deviation);
    sprintf(buffer4, "%sLate takeoff deviation  = %d\n", comment, m_late_takeoff_deviation);
    sprintf(buffer5, "%sEarly landing deviation = %d\n", comment, m_early_land_deviation);
    sprintf(buffer6, "%sLate landing deviation  = %d\n", comment, m_late_land_deviation);
    sprintf(buffer7, "%sDeviation               = %f\n", comment, m_deviation);
    sprintf(buffer8, "%sDelay coefficient       = %f\n", comment, m_delay_coeff);
    sprintf(buffer9, "%sAmount late coefficient = %f\n", comment, m_amount_late_coeff);

    // prepare the string and return it
    int length = strlen(buffer1) + strlen(buffer2) + strlen(buffer3) + strlen(buffer4) +
            strlen(buffer5) + strlen(buffer6) + strlen(buffer7) + strlen(buffer8) + strlen(buffer9) + 20;

    if (output_buffer) {
        delete[] output_buffer;
    }

    output_buffer = new char[length];
    sprintf(output_buffer, "%s%s%s%s%s%s%s%s%s\n", buffer1, buffer2, buffer3, buffer4, buffer5, buffer6,
            buffer7, buffer8, buffer9);

    return output_buffer;
}

void Matching::randomizeParams(Param &p)
{
    do {
        if (p.param_dists[(int) Param::MATCHING] != NULL) {
            double params[7];
            double param_mins[] = {0, 0, 0, 0, 0, 0.5, 0.5};
            double param_maxs[] = {Param::Round(0.05 * p.num_goods / 4.0),
                Param::Round(0.1 * p.num_goods / 4.0),
                Param::Round(0.05 * p.num_goods / 4.0),
                Param::Round(0.1 * p.num_goods / 4.0),
                1, 0.95, 0.95};

            randomizeParamsWithDistribution(p.param_dists[(int) Param::MATCHING], 7, 
                    param_mins, param_maxs, params, p);
            m_early_takeoff_deviation = (int) floor(params[0] + 0.5);
            m_late_takeoff_deviation = (int) floor(params[1] + 0.5);
            m_early_land_deviation = (int) floor(params[2] + 0.5);
            m_late_land_deviation = (int) floor(params[3] + 0.5);
            m_deviation = params[4];
            m_delay_coeff = params[5];
            m_amount_late_coeff = params[6];
        }
        else {
            //	m_max_airport_price = Param::LRand(min,max);		// 5
            m_early_takeoff_deviation = Param::LRand(0, Param::Round(0.05 * p.num_goods / 4.0)); // 1
            m_late_takeoff_deviation = Param::LRand(0, Param::Round(0.1 * p.num_goods / 4.0)); // 2
            m_early_land_deviation = Param::LRand(0, Param::Round(0.05 * p.num_goods / 4.0)); // 1
            m_late_land_deviation = Param::LRand(0, Param::Round(0.1 * p.num_goods / 4.0)); // 2
            m_deviation = Param::DRand(0, 1); // 0.5
            m_delay_coeff = Param::DRand(0.5, 0.95); // 0.9
            m_amount_late_coeff = Param::DRand(0.5, 0.95); // 0.75
        }
    }
    while ((m_early_takeoff_deviation ? 1 : 0) + (m_late_takeoff_deviation ? 1 : 0)
            + (m_early_land_deviation ? 1 : 0) + (m_late_land_deviation ? 1 : 0) 
            < 3 && m_early_takeoff_deviation + m_late_land_deviation == 0);

    p.num_goods -= p.num_goods % 4;
}

