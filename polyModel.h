#ifndef _QUAD_POLY_H
#define _QUAD_POLY_H

class PolyModel {
public:

    static const int EMPTY;

    int degree;
    int numVars;
    int numFreeVars;
    int numCoeffs;

    double *argMin;
    double minimum;
    bool minimumKnown;
    bool isConstant;

    char *filename;

private:
    double *coeffs;
    bool *varIsFree;

public:

    PolyModel(const char* filename);
    PolyModel(const PolyModel &orig);
    ~PolyModel();

    double evalAt(double pos[]);
    double derivAt(double pos[], int var);
    double indefIntegralAt(int var, double val);

    int indexOf(int *vars);
    double &coeffOf(int *vars);
    double constantTerm();

    void integrateOut(int var, double min, double max);
    void instantiate(int var, double value);

    void multiplyBy(double constant);
    void add(double constant);

    void findMinimum(double mins[], double maxs[]);

    void print();
    void testIncrement();

private:
    int nChooseK(int n, int k);
    int simplicalNK(int n, int k);
    bool increment(int *vars);

};

#endif
