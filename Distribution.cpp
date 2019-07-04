#include "Distribution.h"


#include <assert.h>

#include <stdio.h>

#include <iostream>

#include <stdlib.h>


#include "polyModel.h"


#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))


void Distribution::randomizeParamsWithDistribution(PolyModel *dist, int num_params, double *param_mins, double *param_maxs, double *params, Param &p)

{


    if (dist->numVars != num_params) {

        printf("Error: parameter model has incorrect number of params for distribution.\n");

        exit(1);

    }


    PolyModel normalizedDist(*dist);


    if (!p.no_normalization) {


        // adjust estimate to have minimum zero

        if (!dist->minimumKnown) {

            if (p.output_parameter_settings)

                printf("Normalizing distribution from \"%s\"...\n", dist->filename);

            dist->findMinimum(param_mins, param_maxs);

        }


        if (!dist->isConstant)

            normalizedDist.add(-dist->minimum);


        // integrate to get normalizing constant

        PolyModel copyForNorm(normalizedDist);

        for (int i = 0; i < copyForNorm.numVars; i++) {

            assert(param_mins[i] <= param_maxs[i]);

            if (param_mins[i] == param_maxs[i])

                copyForNorm.instantiate(i, param_mins[i]);

            else

                copyForNorm.integrateOut(i, param_mins[i], param_maxs[i]);

        }


        double divideFactor = copyForNorm.constantTerm();


        normalizedDist.multiplyBy(1.0 / divideFactor);

    }


    for (int i = 0; i < num_params; i++) {


        if (param_mins[i] == param_maxs[i]) {

            params[i] = param_mins[i];

            normalizedDist.instantiate(i, param_mins[i]);

            continue;

        }


        PolyModel copy(normalizedDist);


        // integrate out non-instantiated vars other than this

        for (int j = i + 1; j < num_params; j++) {

            if (param_mins[j] == param_maxs[j])

                copy.instantiate(j, param_mins[j]);

            else

                copy.integrateOut(j, param_mins[j], param_maxs[j]);

        }


        double integralBase = copy.indefIntegralAt(i, param_mins[i]);

        double integral = copy.indefIntegralAt(i, param_maxs[i]) - integralBase;


        double targetCDFval = Param::DRand(0, 1);

        double search_min = param_mins[i];

        double search_max = param_maxs[i];

        bool found = false;

        double guess;


        // find val for this var

        while (!found) {

            guess = (search_min + search_max) / 2;


            double integralAtGuess = copy.indefIntegralAt(i, guess) - integralBase;


            double cdfAtGuess = integralAtGuess / integral;


            if (cdfAtGuess > targetCDFval + 0.000001)

                search_max = guess;

            else if (cdfAtGuess < targetCDFval - 0.000001)

                search_min = guess;

            else

                found = true;

        }


        // instantiate that var with new val

        normalizedDist.instantiate(i, guess);

        params[i] = guess;

    }


    // what's left in normalizedDist after everything is

    // instantiated is pdf at these params

    assert(normalizedDist.numFreeVars == 0);

    probOfParams = normalizedDist.constantTerm();

}

