#include "lpkit.h"
#include <stdio.h>

int main(int argc, char *argv[]) {
    lprec *lp;
    FILE *fpin, *fpout;
    int i;

    for (i = 1; i < argc; i++)
        if (argv[i][0] == '-') {
            printf("mps to lp file converter\n");
            printf("Usage: mps2lp [inputfile.mps [outputfile.lp]] [<inputfile.mps] [>outputfile.lp]\n");
            return (1);
        }

    if (argc >= 2) {
        fpin = fopen(argv[1], "r");
        if (fpin == NULL) {
            fprintf(stderr, "Unable to open input file %s\n", argv[1]);
            return (2);
        }
    } else
        fpin = stdin;

    if (argc >= 3) {
        fpout = fopen(argv[2], "w");
        if (fpout == NULL) {
            fprintf(stderr, "Unable to open output file %s\n", argv[2]);
            return (3);
        }
    } else
        fpout = stdout;

    fprintf(stderr, "reading mps file\n");
    lp = read_mps(fpin, FALSE);
    if (fpin != stdin)
        fclose(fpin);

    if (lp == NULL) {
        fprintf(stderr, "Unable to read mps file\n");
        return (4);
    } else {
        fprintf(stderr, "writing lp file\n");
        if (!write_LP(lp, fpout)) {
            fprintf(stderr, "Unable to write lp file\n");
            return (5);
        }
    }
    if (fpout != stdout)
        fclose(fpout);

    return (0);
}
