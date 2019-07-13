#include "polyModel.h"

#include <fstream>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

// #ifdef LINUX
// #include <values.h>
// #else
// #include <limits.h>
// #endif

#include "Param.h"

const int PolyModel::EMPTY = -1;

PolyModel::PolyModel(const char *fname)
{
    filename = new char[strlen(fname) + 1];
    strcpy(filename, fname);

    minimumKnown = false;
    argMin = NULL;

    std::ifstream coeffFile(filename);

    if (!coeffFile) {
        printf("Could not open coefficient file \"%s\".\n", filename);
        exit(1);
    }

    // ==    fprintf(stderr, "OPENED: %s\n", filename);
    coeffFile >> degree;

    if (degree < 1) {
        printf("Invalid number of degrees (%d) specified in coefficient file \"%s\".\n", numVars, filename);
        exit(1);
    }

    coeffFile >> numVars;

    if (numVars < 1) {
        printf("Invalid number of parameters (%d) specified in coefficient file \"%s\".\n", 
                numVars, filename);
        exit(1);
    }

    numFreeVars = numVars;
    numCoeffs = simplicalNK(degree, numVars + 1);
    coeffs = new double[numCoeffs];
    varIsFree = new bool[numVars];

    for (int i = 0; i < numVars; i++) {
        varIsFree[i] = true;
    }

    isConstant = true;

    // read coeffs
    for (int i = 0; i < numCoeffs; i++) {
        if (coeffFile.eof()) {
            printf("Too few coefficients (%d) in file \"%s\".  Expected %d.\n", i, filename, numCoeffs);
            exit(1);
        }

        coeffFile >> coeffs[i];

        if (i > 0 && coeffs[i] != 0) {
            isConstant = false;
        }

        // eat whitespace
        char next;
        while (coeffFile.get(next)) {
            if (!isspace(next)) {
                coeffFile.putback(next);
                break;
            }
        }
    }

    if (!coeffFile.eof()) {
        printf("Too many coefficients specified in file \"%s\".  Expected %d.\n", filename, numCoeffs);
        exit(1);
    }
}

PolyModel::PolyModel(const PolyModel &orig)
{
    degree = orig.degree;
    numVars = orig.numVars;
    numFreeVars = orig.numFreeVars;

    filename = new char[strlen(orig.filename) + 1];
    strcpy(filename, orig.filename);

    numCoeffs = orig.numCoeffs;
    coeffs = new double[numCoeffs];
    varIsFree = new bool[numVars];

    for (int i = 0; i < numCoeffs; i++) {
        coeffs[i] = orig.coeffs[i];
    }

    for (int i = 0; i < numVars; i++) {
        varIsFree[i] = orig.varIsFree[i];
    }

    minimumKnown = orig.minimumKnown;
    minimum = orig.minimum;

    if (orig.argMin) {
        argMin = new double[numVars];

        for (int i = 0; i < numVars; i++) {
            argMin[i] = orig.argMin[i];
        }
    }
    else {
        argMin = NULL;
    }
}

PolyModel::~PolyModel()
{
    delete[] coeffs;
    delete[] varIsFree;
    delete[] filename;

    if (argMin) {
        delete[] argMin;
    }
}

double PolyModel::evalAt(double pos[])
{
    int vars[degree];

    for (int i = 0; i < degree; i++) {
        vars[i] = EMPTY;
    }

    double answer = 0;

    do {
        double val = coeffOf(vars);

        for (int i = 0; i < degree && vars[i] != EMPTY; i++) {
            val *= pos[vars[i]];
        }

        answer += val;
    }
    while (increment(vars));

    return answer;
}

double PolyModel::derivAt(double pos[], int var)
{
    int vars[degree];

    for (int i = 0; i < degree; i++) {
        vars[i] = EMPTY;
    }

    double answer = 0;

    do {
        int varDegree = 0;

        for (int i = 0; i < degree; i++) {
            if (vars[i] == var) {
                varDegree++;
            }
            else if (vars[i] < var) {
                break;
            }
        }

        if (varDegree == 0) {
            continue;
        }

        double derivOfTerm = coeffOf(vars) * varDegree;

        for (int i = 0; i < degree && vars[i] != EMPTY; i++) {
            // skips first instance of var, which has been "derivated out"
            if (vars[i] != var || (i != 0 && vars[i - 1] == var)) {
                derivOfTerm *= pos[vars[i]];
            }
        }

        answer += derivOfTerm;
    }
    while (increment(vars));

    return answer;
}

double PolyModel::indefIntegralAt(int var, double value)
{
    // bad things happen if this function is misused...
    assert(numFreeVars == 1);

    for (int i = 0; i < numVars; i++) {
        if (i == var) {
            assert(varIsFree[i] == true);
        }
        else {
            assert(varIsFree[i] == false);
        }
    }

    int vars[degree];

    for (int i = 0; i < degree; i++) {
        vars[i] = EMPTY;
    }

    double answer = 0;

    for (int i = 0; i <= degree; i++) {
        answer += coeffOf(vars) * pow(value, i + 1) / (i + 1);

        if (i < degree) {
            vars[i] = var;
        }
    }

    return answer;
}

// returns next term in poly
bool PolyModel::increment(int *vars)
{
    int lowestFreeVar = 0, highestFreeVar = numVars - 1;

    while (!varIsFree[lowestFreeVar]) {
        if (lowestFreeVar == numVars - 1) {
            return false;
        }
        else {
            lowestFreeVar++;
        }
    }

    while (!varIsFree[highestFreeVar]) {
        highestFreeVar--;
    }

    int pos = degree - 1;

    while (vars[pos] == EMPTY && pos > 0) {
        pos--;
    }

    if (pos == 0) {
        if (vars[pos] == EMPTY) {
            vars[pos] = lowestFreeVar;
            return true;
        }
    }
    else {
        while (pos > 0 && vars[pos - 1] == vars[pos])
            pos--;
    }

    if (pos == 0 && vars[pos] == highestFreeVar) {
        while (pos < degree && vars[pos] != EMPTY) {
            vars[pos++] = lowestFreeVar;
        }

        if (pos == degree) {
            return false;
        }
        else {
            vars[pos] = lowestFreeVar;
        }
    }
    else {
        vars[pos]++;

        if (vars[pos] > highestFreeVar) {
            return false;
        }
        else {
            while (!varIsFree[vars[pos]]) {
                vars[pos]++;
            }
        }

        while (++pos < degree && vars[pos] != EMPTY) {
            vars[pos] = lowestFreeVar;
        }
    }

    return true;
}

int PolyModel::indexOf(int *vars)
{
    int lastEmpty = degree;
    int lastVar = EMPTY;
    int index = 0;

    for (int pos = degree - 1; pos >= 0; pos--) {
        assert(vars[pos] >= lastVar);

        if (vars[pos] == EMPTY) {
            continue;
        }

        if (lastVar == EMPTY) {
            lastEmpty = pos + 1;
            index += simplicalNK(pos, numVars + 1);
        }

        index += simplicalNK(lastEmpty - pos, vars[pos]);
        lastVar = vars[pos];
    }

    return index;
}

double &PolyModel::coeffOf(int *vars)
{
    return coeffs[indexOf(vars)];
}

double PolyModel::constantTerm()
{
    return coeffs[0];
}

void PolyModel::integrateOut(int var, double min, double max)
{
    assert(var < numVars && var >= 0);
    assert(varIsFree[var] == true);

    int vars[degree];

    for (int i = 0; i < degree; i++) {
        vars[i] = EMPTY;
    }

    double oldCoeffs[numCoeffs];

    for (int i = 0; i < numCoeffs; i++) {
        oldCoeffs[i] = coeffs[i];
        coeffs[i] = 0;
    }

    do {
        int varDegree = 0;
        int integratedTerm[degree];

        for (int i = 0; i < degree; i++) {
            if (vars[i] == var) {
                varDegree++;
            }
            else {
                integratedTerm[i - varDegree] = vars[i];
            }
        }

        for (int i = degree - varDegree; i < degree; i++) {
            integratedTerm[i] = EMPTY;
        }

        coeffOf(integratedTerm) +=
                (oldCoeffs[indexOf(vars)] / (varDegree + 1)) * (pow(max, varDegree + 1) - pow(min, 
                            varDegree + 1));
    }
    while (increment(vars));

    varIsFree[var] = false;
    numFreeVars--;
}

void PolyModel::instantiate(int var, double value)
{
    assert(var < numVars && var >= 0);
    assert(varIsFree[var] == true);

    int vars[degree];

    for (int i = 0; i < degree; i++) {
        vars[i] = EMPTY;
    }

    double oldCoeffs[numCoeffs];

    for (int i = 0; i < numCoeffs; i++) {
        oldCoeffs[i] = coeffs[i];
        coeffs[i] = 0;
    }

    do {
        int varDegree = 0;
        int instantiatedTerm[degree];

        for (int i = 0; i < degree; i++) {
            if (vars[i] == var) {
                varDegree++; 
            }
            else {
                instantiatedTerm[i - varDegree] = vars[i];
            }
        }

        for (int i = degree - varDegree; i < degree; i++) {
            instantiatedTerm[i] = EMPTY;
        }

        coeffOf(instantiatedTerm) += oldCoeffs[indexOf(vars)] * pow(value, varDegree);
    }
    while (increment(vars));

    varIsFree[var] = false;
    numFreeVars--;
}

void PolyModel::multiplyBy(double constant)
{
    for (int i = 0; i < numCoeffs; i++) {
        coeffs[i] *= constant;
    }
}

void PolyModel::add(double constant)
{
    coeffs[0] += constant;
}

void PolyModel::findMinimum(double mins[], double maxs[])
{
    const double learnRate = 0.001;
    const int numRestarts = 50;

    double *width = new double[numVars];

    for (int i = 0; i < numVars; i++) {
        assert(mins[i] <= maxs[i]);
        width[i] = maxs[i] - mins[i];
    }

    double *pos = new double[numVars];
    double *step = new double[numVars];
    double valAtPos;

    argMin = new double[numVars];
    double lowestVal;
    bool lowestValValid = false;

    for (int i = 0; i < numRestarts; i++) {
        for (int j = 0; j < numVars; j++) {
            pos[j] = Param::DRand(mins[j], maxs[j]);
        }

        valAtPos = evalAt(pos);

        bool done = false;

        while (!done) {
            for (int j = 0; j < numVars; j++) {
                step[j] = -(derivAt(pos, j) * learnRate * width[j]);
            }

            for (int j = 0; j < numVars; j++) {
                pos[j] += step[j];

                if (pos[j] < mins[j]) {
                    pos[j] = mins[j];
                }
                else if (pos[j] > maxs[j]) {
                    pos[j] = maxs[j];
                }
            }

            double valAtLastPos = valAtPos;
            valAtPos = evalAt(pos);

            if (valAtPos >= valAtLastPos) {
                done = true;
            }
        }

        if (!lowestValValid || valAtPos < lowestVal) {
            for (int j = 0; j < numVars; j++) {
                argMin[j] = pos[j];
            }

            lowestVal = valAtPos;
            lowestValValid = true;
        }
    }

    minimumKnown = true;
    minimum = lowestVal;

    delete[] width;
    delete[] pos;
    delete[] step;
}

int PolyModel::nChooseK(int n, int k)
{
    int result = 1;

    for (int i = k + 1; i <= n; i++) {
        result *= i;
    }

    for (int i = 1; i <= (n - k); i++) {
        result /= i;
    }

    return result;
}

inline int PolyModel::simplicalNK(int n, int k)
{
    if (n == 0) {
        return 1;
    }
    else if (k == 0) {
        return 0;
    }
    else {
        return nChooseK(n + k - 1, k - 1);
    }
}

void PolyModel::print()
{
    int vars[degree];

    for (int i = 0; i < degree; i++) {
        vars[i] = EMPTY;
    }

    do {
        for (int i = 0; i < degree; i++) {
            if (vars[i] != EMPTY) {
                printf("%d", vars[i]);
            }
            else {
                break; 
            }
        }

        printf(":\t");
        printf("%f\n", coeffOf(vars));
    }
    while (increment(vars));
}

void PolyModel::testIncrement()
{
    int *vars = new int[degree];

    for (int i = 0; i < degree; i++) {
        vars[i] = EMPTY;
    }

    while (increment(vars)) {
        for (int i = 0; i < degree; i++) {
            if (vars[i] != EMPTY) {
                printf("%d", vars[i]); 
            }
        }

        printf(": %d\n", indexOf(vars));
    }
}

