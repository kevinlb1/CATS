// Legacy.cpp: implementation of the Legacy class.

//

//////////////////////////////////////////////////////////////////////



#include <stdio.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <strings.h>


#include "Legacy.h"



//////////////////////////////////////////////////////////////////////

// Construction/Destruction

//////////////////////////////////////////////////////////////////////

Legacy::Legacy(Param &p, char * argv[], int argc)

{

    norm = NULL;

    output_buffer = NULL;



    // first set default values

    num_goods_type = DECAY;

    binom_cutoff = 0.2;

    q = 5;

    const_goods = 3;

    mu_goods = 4;

    sigma_goods = 1;

    pricing_type = LINEAR_RANDOM;

    high_fixed = 1000;

    low_fixed = 0;

    high_linearly = 1000;

    low_linearly = 1;

    mu_price = 16;

    sigma_price = 3;

    alpha = 0.55;



    // read through the parameters.  Ignore those that don't apply

    for (int i = 1; i < argc; i++) {

        // now read in values

        if (argv[i][0] == '-' || argv[i][0] == '/') {

            char *in = &(argv[i][1]); // skip the first character



            // num_goods_type

            if (!strcasecmp(in, "num_goods_type")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                num_goods_type = (goods_enum) atoi(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



            // binom_cutoff

            if (!strcasecmp(in, "binom_cutoff")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                binom_cutoff = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // q

            else if (!strcasecmp(in, "q")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                q = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // const_goods

            else if (!strcasecmp(in, "const_goods")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                const_goods = atoi(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // mu_goods

            else if (!strcasecmp(in, "mu_goods")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                mu_goods = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // sigma_goods

            else if (!strcasecmp(in, "sigma_goods")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                sigma_goods = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // pricing_type

            else if (!strcasecmp(in, "pricing_type")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                pricing_type = (pricing_enum) atoi(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // high_fixed

            else if (!strcasecmp(in, "high_fixed")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                high_fixed = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // low_fixed

            else if (!strcasecmp(in, "low_fixed")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                low_fixed = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // high_linearly

            else if (!strcasecmp(in, "high_linearly")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                high_linearly = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // low_linearly

            else if (!strcasecmp(in, "low_linearly")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                low_linearly = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // mu_price

            else if (!strcasecmp(in, "mu_price")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                mu_price = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // sigma_price

            else if (!strcasecmp(in, "sigma_price")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                sigma_price = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }



                // alpha

            else if (!strcasecmp(in, "alpha")) {

                p.argRecognized[i] = true;

                if (i + 1 >= argc) usage();

                alpha = atof(argv[i + 1]);

                i++;

                p.argRecognized[i] = true;

            }

        }

    }



    switch (p.distribution) {

        // L1

    case Param::L1:

        num_goods_type = RANDOM;

        pricing_type = FIXED_RANDOM;

        break;



        // L2

    case Param::L2:

        num_goods_type = RANDOM;

        pricing_type = LINEAR_RANDOM;

        break;



        // L3

    case Param::L3:

        num_goods_type = CONSTANT;

        pricing_type = FIXED_RANDOM;

        break;



        // L4

    case Param::L4:

        num_goods_type = DECAY;

        pricing_type = LINEAR_RANDOM;

        break;



        // L5

    case Param::L5:

        num_goods_type = NORMAL_GOODS;

        pricing_type = NORMAL_PRICE;

        break;



        // L6

    case Param::L6:

        num_goods_type = EXPONENTIAL;

        pricing_type = LINEAR_RANDOM;

        break;



        // L7

    case Param::L7:

        num_goods_type = BINOMIAL;

        pricing_type = LINEAR_RANDOM;

        break;



        // L8

    case Param::L8:

        num_goods_type = CONSTANT;

        pricing_type = QUADRATIC;

        break;



    default:

        printf("Error: specific legacy distribution (L1 - L8) must be chosen.\n");

        usage();

        exit(1);

    }



    if (num_goods_type == NORMAL_GOODS ||

        pricing_type == NORMAL_PRICE)

        norm = new Normal(p);



}

void Legacy::randomizeParams(Param &p)
{



    if (p.param_dists[(int) p.distribution] != NULL) {



        // arrays have 4 elts to accomodate L5 with 4 params

        double params[4];

        double param_mins[4], param_maxs[4];



        // all dists (except L5) take low_linearly or low_fixed

        // as the second param

        param_mins[1] = 0;

        param_maxs[1] = 1;



        switch (p.distribution) {

        case Param::L1:

            // add one to mins, maxs because L1 has only one param

            randomizeParamsWithDistribution(p.param_dists[(int) p.distribution], 1, param_mins + 1, param_maxs + 1, params, p);



            low_linearly = params[0];

            high_linearly = 2 - low_linearly;

            break;



        case Param::L2:

            // add one to mins, maxs because L2 has only one param

            randomizeParamsWithDistribution(p.param_dists[(int) p.distribution], 1, param_mins + 1, param_maxs + 1, params, p);



            low_linearly = params[0];

            high_linearly = 2 - low_linearly;

            break;



        case Param::L3:

        case Param::L8:

            // const_goods



            param_mins[0] = 3;

            param_maxs[0] = 10;



            randomizeParamsWithDistribution(p.param_dists[(int) p.distribution], 2, param_mins, param_maxs, params, p);



            const_goods = (int) params[0];

            low_fixed = params[1];



            break;



        case Param::L4:

            // alpha

            param_mins[0] = 0.05;

            param_maxs[0] = 0.95;



            randomizeParamsWithDistribution(p.param_dists[(int) p.distribution], 2, param_mins, param_maxs, params, p);



            alpha = params[0];

            low_linearly = params[1];

            high_linearly = 2 - low_linearly;

            break;



        case Param::L5:

            // mu_goods, mu_price, sigma_goods, sigma_price



            param_mins[0] = 2;

            param_maxs[0] = 10;

            param_mins[1] = 10;

            param_maxs[1] = 20;

            param_mins[2] = 0.0001;

            param_maxs[2] = 10;

            param_mins[3] = 0.0001;

            param_maxs[3] = 40;



            randomizeParamsWithDistribution(p.param_dists[(int) p.distribution], 4, param_mins, param_maxs, params, p);



            mu_goods = params[0];

            mu_price = params[1];

            sigma_goods = params[2];

            sigma_price = params[3];

            break;



        case Param::L6:



            param_mins[0] = 0.5;

            param_maxs[0] = 1.0;



            randomizeParamsWithDistribution(p.param_dists[(int) p.distribution], 2, param_mins, param_maxs, params, p);



            q = params[0];

            low_linearly = params[1];

            high_linearly = 2 - low_linearly;

            break;



        case Param::L7:

            param_mins[0] = 0.005;

            param_maxs[0] = 0.4;



            randomizeParamsWithDistribution(p.param_dists[(int) p.distribution], 2, param_mins, param_maxs, params, p);



            binom_cutoff = params[0];

            low_linearly = params[1];

            high_linearly = 2 - low_linearly;

            break;



        default:

            assert(false); // shouldn't be here if its not a LEG dist. (note L8 is with L3)

        }



        high_fixed = 1.0;



    }
    else {

        // L3 / L8

        const_goods = Param::LRand(3, 10);



        // L4

        alpha = Param::DRand(0.5, 0.95);



        // L5

        mu_goods = Param::LRand(2, 10);

        mu_price = Param::LRand(10, 20);

        sigma_goods = Param::DRand(0.0001, mu_goods);

        sigma_price = Param::DRand(0.0001, 2 * mu_price);



        // L6

        q = 1.0 / Param::DRand(0.5, 1);



        // L7

        binom_cutoff = Param::DRand(0.005, 0.4);



        double dev = Param::DRand(0, 1);



        // price deviations

        high_linearly = 1 + dev;

        low_linearly = 1 - dev;

        low_fixed = 1.0 - dev;

        high_fixed = 1.0;

    }

}

void Legacy::usage()

{

    // 80 cols:  12345678901234567890123456789012345678901234567890123456789012345678901234567890

    int count = printf("Usage for Legacy Distribution:\n");

    for (int t = 0; t < count; t++) printf("=");



    printf("\n\nNumber of Goods:\n");

    printf("  -num_goods_type [5]: {1=Binom.;2=Expon.;3=Random;4=Constant;5=Decay;6=Normal}\n");

    printf("  -binom_cutoff [0.2]: probability of adding each good under binomial\n");

    printf("  -q [5]: parameter for exponential\n");

    printf("  -const_goods [3]: number of constant goods\n");

    printf("  -alpha [0.55]: probability of adding another good for decay\n");

    printf("  -mu_goods [4]: mean for normally distributed goods\n");

    printf("  -sigma_goods [1]: stdev for normally distributed goods\n");



    printf("\nPricing:\n");

    printf("  -pricing_type [2]: {1=FixedRandom; 2=LinearRandom; 3=Normal; 4=Quadratic}\n");

    printf("  -low_fixed [0]: smallest fixed value\n");

    printf("  -high_fixed [1000]: largest fixed value\n");

    printf("  -low_linearly [1]: smallest linear value\n");

    printf("  -high_linearly [1000]: largest linear value\n");

    printf("  -mu_price [16]: mean for normally distributed prices\n");

    printf("  -sigma_price [3]: stdev for normally distributed prices\n");



    printf("\nPredefined Legacy Distributions:\n");

    printf("  L1:  Random; FixedRandom, low_fixed=0, high_fixed=1\n");

    printf("  L2:  Random; LinearRandom, low_linearly=0, high_linearly=1\n");

    printf("  L3:  Constant, const_goods=3; FixedRandom, low_fixed=0, high_fixed=1\n");

    printf("  L4:  Decay, alpha=0.55; LinearRandom, low_linearly=0, high_linearly=1\n");

    printf("  L5:  Normal, mu_goods=4, sigma_goods=1; Normal, mu_price=16, sigma_price=3\n");

    printf("  L6:  Exponential, q=5; LinearRandom, low_linearly=0.5, high_linearly=1.5\n");

    printf("  L7:  Binom,binom_cutoff=0.2;LinearRandom,low_linearly=0.5,high_linearly=1.5\n");

    printf("  L8:  Constant, const_goods=3; Quadratic\n");

    //

    //	printed_usage = true;

}

Legacy::~Legacy()

{

    if (norm) delete norm;

    if (output_buffer) delete[] output_buffer;

}



// generate the set of bids

BidSet *Legacy::generate(Param &p)

{

    BidSet *bidset = new BidSet(p.num_goods, p.num_bids);



    // keep adding bids until the right number has been generated

    for (int t = 1; bidset->numBids() < p.num_bids; t++) {

        int num = generateNumGoods(p);

        double price = generatePrice(num, p);

        Bid *bid = generateBid(num, p.num_goods, price); // this selects num goods uniformly without replacement

        bidset->add(bid, p.remove_dominated_bids);



        if (p.output_parameter_settings && (!(t % p.output_frequency) || bidset->numBids() == p.num_bids)) {

            int count = printf("Number of %sbids: %d  Number of tries: %d",

                               (p.remove_dominated_bids ? "non-dominated " : ""), bidset->numBids(), t);

            if (bidset->numBids() < p.num_bids)

                for (int counter = 0; counter < count; counter++)

                    printf("\b");

            else

                printf("\n");

        }



#ifdef OUTPUT_DOMINATION_DATA

        if (t >= 35000) break;

#endif

    }



    if (p.output_parameter_settings) printf("\n");



    return bidset;

}



// Decide how many goods to include in the bid

int Legacy::generateNumGoods(Param &p)

{

    int num_goods = 0;



    switch (num_goods_type) {



        // binomial

    case BINOMIAL:

    {

        do {

            num_goods = 0;

            for (int j = 0; j < p.num_goods; j++) {

                double r = Param::DRand();

                if (r < binom_cutoff)

                    num_goods++;

            }

        }
        while (num_goods == 0);

        break;

    }



        // exponential

    case EXPONENTIAL:

    {

        // a normal exponential distribution is num_goods = - 1/lambda * ln(1-rand()))

        // here we call 1/lambda q, and we make sure that num_goods < p->num_goods

        // note that "log" in c is actually "ln"

        do num_goods = 1 + (int) floor(-q * log(1.0 - (1.0 - exp(-p.num_goods / q)) * Param::DRand()));

        while (num_goods == 0);

        break;

    }



        // Sandholm Random: pick a random number of unique goods

    case RANDOM:

    {

        do num_goods = (int) (Param::DRand() * p.num_goods) + 1;

        while (num_goods > p.num_goods || num_goods == 0); // I think this will never happen anyway

        break;

    }



        // Sandholm Uniform: always the same number of goods

    case CONSTANT:

    {

        num_goods = const_goods;

        break;

    }



        // Sandholm Decay: keep adding more goods with probability p

    case DECAY:

    {

        do for (num_goods = 0; num_goods < p.num_goods && Param::DRand() <= alpha; num_goods++);

            while (num_goods == 0);

        break;

    }



        // normal

    case NORMAL_GOODS:

    {

        do num_goods = Param::Round(norm->draw(mu_goods, sigma_goods));

        while (num_goods == 0 || num_goods > p.num_goods);

        break;

    }

    }



    return num_goods;

}



// generate a price for the bid

double Legacy::generatePrice(int n, const Param &p)

{

    double price = 0.0;



    switch (pricing_type) {

    case FIXED_RANDOM:

    {

        price = Param::DRand() * (high_fixed - low_fixed) + low_fixed;

        break;

    }

    case LINEAR_RANDOM:

    {

        price = Param::DRand() * n * (high_linearly - low_linearly) + low_linearly * n;

        break;

    }

    case NORMAL_PRICE:

    {

        price = norm->draw(mu_price, sigma_price);

        break;

    }

    case QUADRATIC:

    {

        // needs a concept of bidders.  I'm not implementing this one yet

        break;

    }

    }



    return price;

}





// return a bid with n goods and price p

Bid *Legacy::generateBid(int goods_in_bid, int total_goods, double price)

{

    Bid *temp = new Bid(price);



    bool *assigned = new bool[total_goods];

    int n;



    for (n = 0; n < total_goods; n++) assigned[n] = false;



    for (n = goods_in_bid; n > 0; n--) {

        int g = Param::LRand(0, total_goods - 1);

        while (assigned[g])

            g = Param::LRand(0, total_goods - 1);



        temp->addGood(g);

        assigned[g] = true;

    }



    delete[] assigned;



    return temp;

}

char *Legacy::outputSettings(bool tofile)

{

    const char *goods[] = {"", "Binomial", "Exponential", "Random", "Constant", "Decay", "Normal"};

    const char *pricing[] = {"", "Fixed Random", "Linear Random", "Normal", "Quadratic"};

    char comment[3];

    if (tofile) {
        comment[0] = '%';
        comment[1] = ' ';
        comment[2] = 0;
    }

    else comment[0] = 0;

    char buffer1[150], buffer2[100], buffer3[100], buffer4[100];



    sprintf(buffer1, "%sLegacy distribution parameters:\n%sNumber of goods chosen according to %s distribution",

            comment, comment, goods[num_goods_type]);

    switch (num_goods_type) {

    case BINOMIAL:

        sprintf(buffer2, " (binom_cutoff = %g)", binom_cutoff);

        break;

    case EXPONENTIAL:

        sprintf(buffer2, " (q = %g)", q);

        break;

    case CONSTANT:

        sprintf(buffer2, " (const_goods = %d)", const_goods);

        break;

    case DECAY:

        sprintf(buffer2, " (alpha = %g)", alpha);

        break;

    case NORMAL_GOODS:

        sprintf(buffer2, " (mu_goods = %g, sigma_goods = %g)", mu_goods, sigma_goods);

        break;

    default:

        buffer2[0] = 0;

    }

    sprintf(buffer3, ".\n%sPrice chosen according to %s distribution", comment, pricing[pricing_type]);

    switch (pricing_type) {

    case FIXED_RANDOM:

        sprintf(buffer4, " (low_fixed = %g, high_fixed = %g)", low_fixed, high_fixed);

        break;

    case LINEAR_RANDOM:

        sprintf(buffer4, " (low_linearly = %g, high_linearly = %g)", low_linearly, high_linearly);

        break;

    case NORMAL_PRICE:

        sprintf(buffer4, " (mu_price = %g, sigma_price = %g)", mu_price, sigma_price);

        break;

    default:

        buffer4[0] = 0;

    }



    // put it all together



    if (output_buffer)

        delete[] output_buffer;



    output_buffer = new char[strlen(buffer1) + strlen(buffer2) + strlen(buffer3) + strlen(buffer4) + 20];



    sprintf(output_buffer, "%s%s%s%s%s%s.\n", buffer1, (tofile || !buffer2[0] ? "" : "\n  "), buffer2, buffer3,

            (tofile || !buffer4[0] ? "" : "\n  "), buffer4);



    return output_buffer;

}

