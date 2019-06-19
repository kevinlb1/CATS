#include <string.h>
#include <time.h>
#include "lpkit.h"
#include "lpglob.h"
#include "patchlevel.h"

#if defined FORTIFY

int EndOfPgr() {
    exit(1);
    return (0);
}
#endif

void print_help(char *argv[]) {
    printf("Usage of %s version %d.%d.%d.%d:\n", argv[0], MAJORVERSION, MINORVERSION, RELEASE, BUILD);
    printf("%s [options] [[<]input_file]\n", argv[0]);
    printf("List of options:\n");
    printf("-h\t\tprints this message\n");
    printf("-v <level>\tverbose mode, gives flow through the program.\n");
    printf("\t\t if level not provided (-v) then -v2 (NORMAL) is taken.\n");
    printf("\t -v0: CRITICALSTOP\n");
    printf("\t -v1: CRITICAL\n");
    printf("\t -v2: SEVERE\n");
    printf("\t -v3: IMPORTANT (default)\n");
    printf("\t -v4: NORMAL\n");
    printf("\t -v5: DETAILED\n");
    printf("\t -v6: FULL\n");
    printf("-d\t\tdebug mode, all intermediate results are printed,\n\t\tand the branch-and-bound decisions\n");
    printf("-p\t\tprint the values of the dual variables\n");
    printf("-b <bound>\tspecify a lower bound for the objective function\n\t\tto the program. If close enough, may speed up the\n\t\tcalculations.\n");
    printf("-i\t\tprint all intermediate valid solutions.\n\t\tCan give you useful solutions even if the total run time\n\t\tis too long\n");
    printf("-e <number>\tspecifies the epsilon which is used to determine whether a\n\t\tfloating point number is in fact an integer.\n\t\tShould be < 0.5\n");
    printf("-c\t\tduring branch-and-bound, take the ceiling branch first\n");
    printf("-ca\t\tduring branch-and-bound, the algorithm chooses branch\n");
    printf("-B <rule>\tspecify branch-and-bound rule\n");
    printf("\t -B0: Select Lowest indexed non-integer column (default)\n");
    printf("\t -B1: Select Random non-integer column\n");
    printf("\t -B2: Select Largest deviation from an integer value\n");
    printf("\t -B3: Select Best ???\n");
    printf("\t -B4: Select Median value deviation from an integer value\n");
    printf("\t -B5: Select Greedy ???\n");
    printf("-s <mode>\tuse automatic problem scaling.\n");
    printf("\t  -s:\n");
    printf("\t -s0: Numerical range-based scaling\n");
    printf("\t -s1: Geometric scaling\n");
    printf("\t -s2: Curtis-reid scaling\n");
    printf("-sp\t\talso do power scaling.\n\t\tThis option must come AFTER -s <mode> option\n");
    printf("-sl\t\talso do Lagrange scaling.\n\t\tThis option must come AFTER -s <mode> option\n");
    printf("-si\t\talso do Integer scaling.\n\t\tThis option must come AFTER -s <mode> option\n");
    printf("-I\t\tprint info after reinverting\n");
    printf("-t\t\ttrace pivot selection\n");
    printf("-lp\t\tread from LP file (default)\n");
    printf("-mps\t\tread from MPS file\n");
    printf("-degen\t\tuse perturbations to reduce degeneracy,\n\t\tcan increase numerical instability\n");
    printf("-trej <Trej>\tset minimum pivot value\n");
    printf("-parse_only\tparse input file but do not calculate (ie check)\n");
    printf("-presolve\tpresolve problem before start optimizing\n");
    printf("-improve <level>\titerative improvement level\n");
    printf("\t -improve0: none (default)\n");
    printf("\t -improve1: FTRAN only\n");
    printf("\t -improve2: BTRAN only\n");
    printf("\t -improve3: FTRAN + BTRAN\n");
    printf("-time\t\tPrint CPU time to parse input and to calculate result\n");
    printf("-min\t\tMinimize the lp problem (overrules setting in file)\n");
    printf("-max\t\tMaximize the lp problem (overrules setting in file)\n");
    printf("-S <detail>\tPrint solution. If detail ommited, then -S2 is used.\n");
    printf("\t -S0: Print nothing\n");
    printf("\t -S1: Only objective value\n");
    printf("\t -S2: Objective value + variables (default)\n");
    printf("\t -S3: Objective value + variables + constraints\n");
    printf("\t -S4: Objective value + variables + constraints + duals\n");
    printf("\t -S5: Objective value + variables + constraints + duals + lp model\n");
    printf("\t -S6: Objective value + variables + constraints + duals + lp model + lp scales\n");
}

void print_cpu_times(const char *info) {
    static clock_t last_time = 0;
    clock_t new_time;

    new_time = clock();
    fprintf(stderr, "CPU Time for %s: %gs (%gs total since program start)\n",
            info, (new_time - last_time) / (double) CLOCKS_PER_SEC,
            new_time / (double) CLOCKS_PER_SEC);
    last_time = new_time;
}

int main(int argc, char *argv[]) {
    lprec *lp;
    char *filen;
    int i;
    short verbose = IMPORTANT /* CRITICALSTOP */;
    MYBOOL debug = FALSE;
    MYBOOL print_sol = FALSE;
    MYBOOL PRINT_DUALS = FALSE;
    MYBOOL floor_first = TRUE;
    short scaling = 0;
    short print_at_invert = FALSE;
    MYBOOL tracing = FALSE;
    short mps = FALSE;
    MYBOOL anti_degen = FALSE;
    short print_timing = FALSE;
    short parse_only = FALSE;
    MYBOOL do_presolve = FALSE;
    short objective = 0;
    short PRINT_SOLUTION = 2;
    MYBOOL improve = IMPROVE_NONE;
    MYBOOL bb_rule = FIRST_SELECT;
    MYBOOL scalemode = MMSCALING;
    int result;
    REAL obj_bound = (REAL) DEF_INFINITE;
    REAL epsilon = (REAL) DEF_EPSILON;
    REAL epspivot = (REAL) DEF_EPSPIVOT;
    FILE *fpin = stdin;

    /* read command line arguments */

#if defined FORTIFY
    Fortify_EnterScope();
#endif

    for (i = 1; i < argc; i++) {
        if (strncmp(argv[i], "-v", 2) == 0) {
            if (argv[i][2])
                verbose = (short) atoi(argv[i] + 2);
            else
                verbose = NORMAL;
        } else if (strcmp(argv[i], "-d") == 0)
            debug = TRUE;
        else if (strcmp(argv[i], "-i") == 0)
            print_sol = TRUE;
        else if (strcmp(argv[i], "-c") == 0)
            floor_first = FALSE;
        else if (strcmp(argv[i], "-ca") == 0)
            floor_first = AUTOMATIC;
        else if (strncmp(argv[i], "-B", 2) == 0) {
            if (argv[i][2])
                bb_rule = (MYBOOL) atoi(argv[i] + 2);
            else
                bb_rule = FIRST_SELECT;
        } else if (strcmp(argv[i], "-b") == 0)
            obj_bound = atof(argv[++i]);
        else if (strcmp(argv[i], "-e") == 0) {
            epsilon = atof(argv[++i]);
            if ((epsilon <= 0.0) || (epsilon >= 0.5)) {
                fprintf(stderr, "Invalid epsilon %g; 0 < epsilon < 0.5\n",
                        (double) epsilon);
                exit(EXIT_FAILURE);
            }
        } else if (strcmp(argv[i], "-p") == 0)
            PRINT_DUALS = TRUE;
        else if (strcmp(argv[i], "-h") == 0) {
            print_help(argv);
            exit(EXIT_SUCCESS);
        } else if (strcmp(argv[i], "-sp") == 0)
            scalemode = (MYBOOL) (scalemode | POWERSCALE);
        else if (strcmp(argv[i], "-sl") == 0)
            scalemode = (MYBOOL) (scalemode | LAGRANGESCALE);
        else if (strcmp(argv[i], "-s2") == 0)
            scaling = 2;
        else if (strcmp(argv[i], "-si") == 0)
            scalemode = (MYBOOL) (scalemode | INTEGERSCALE);
        else if (strncmp(argv[i], "-s", 2) == 0) {
            scaling = 1;
            if (argv[i][2])
                scalemode = (MYBOOL) atoi(argv[i] + 2);
            else
                scalemode = MMSCALING;
        } else if (strcmp(argv[i], "-I") == 0)
            print_at_invert = TRUE;
        else if (strcmp(argv[i], "-t") == 0)
            tracing = TRUE;
        else if (strncmp(argv[i], "-S", 2) == 0) {
            if (argv[i][2])
                PRINT_SOLUTION = (short) atoi(argv[i] + 2);
            else
                PRINT_SOLUTION = 2;
        } else if (strncmp(argv[i], "-improve", 8) == 0) {
            if (argv[i][8])
                improve = (MYBOOL) atoi(argv[i] + 8);
            else
                improve = 0;
        } else if (strcmp(argv[i], "-mps") == 0)
            mps = TRUE;
        else if (strcmp(argv[i], "-lp") == 0)
            mps = FALSE;
        else if (strcmp(argv[i], "-degen") == 0)
            anti_degen = TRUE;
        else if (strcmp(argv[i], "-time") == 0) {
            if (clock() == -1)
                fprintf(stderr, "CPU times not available on this machine\n");
            else
                print_timing = TRUE;
        } else if (strcmp(argv[i], "-trej") == 0)
            epspivot = atof(argv[++i]);
        else if (strcmp(argv[i], "-parse_only") == 0)
            /* only useful for parser software development */
            parse_only = TRUE;
        else if (strcmp(argv[i], "-presolve") == 0)
            do_presolve = TRUE;
        else if (strcmp(argv[i], "-min") == 0)
            objective = -1;
        else if (strcmp(argv[i], "-max") == 0)
            objective = 1;
        else if (fpin == stdin) {
            filen = argv[i];
            if (*filen == '<')
                filen++;
            if ((fpin = fopen(filen, "r")) == NULL) {
                fprintf(stderr, "Error, Unable to open input file '%s'\n",
                        argv[i]);
                print_help(argv);
                exit(EXIT_FAILURE);
            }
        } else {
            filen = argv[i];
            if (*filen != '>') {
                fprintf(stderr, "Error, Unrecognized command line argument '%s'\n",
                        argv[i]);
                print_help(argv);
                exit(EXIT_FAILURE);
            }
        }
    }

    if (mps)
        lp = read_mps(fpin, verbose);
    else /* standard lp_solve syntax expected */
        lp = read_lp_file(fpin, verbose, "lp");

    if (fpin != stdin)
        fclose(fpin);

    if (print_timing)
        print_cpu_times("Parsing input");

    if (parse_only)
        exit(0);

    if (lp == NULL) {
        fprintf(stderr, "Unable to read model.\n");
        exit(EXIT_FAILURE);
    }

    if (objective != 0) {
        if (objective == 1)
            set_maxim(lp);
        else
            set_minim(lp);
    }

    if (PRINT_SOLUTION >= 5)
        print_lp(lp);

    set_print_sol(lp, print_sol);
    set_epsilon(lp, epsilon);
    set_epspivot(lp, epspivot);
    set_print_duals(lp, PRINT_DUALS);
    set_debug(lp, debug);
    set_floor_first(lp, floor_first);
    set_print_at_invert(lp, print_at_invert);
    set_trace(lp, tracing);
    if (obj_bound != DEF_INFINITE)
        set_obj_bound(lp, obj_bound);
    set_anti_degen(lp, anti_degen);
    set_do_presolve(lp, do_presolve);
    set_improve(lp, improve);
    set_scalemode(lp, scalemode);
    set_bb_rule(lp, bb_rule);

    if (scaling == 1)
        auto_scale(lp);
    if (scaling == 2)
        scaleCR(lp);

    if (PRINT_SOLUTION >= 6)
        print_scales(lp);

    result = solve(lp);

    if (print_timing)
        print_cpu_times("solving");

    if (result == OPTIMAL) {
        if (PRINT_SOLUTION >= 1)
            print_objective(lp);

        if (PRINT_SOLUTION >= 2)
            print_solution(lp);

        if (PRINT_SOLUTION >= 3)
            print_constraints(lp);

        if (PRINT_SOLUTION >= 4)
            print_duals(lp);

        if (tracing)
            fprintf(stderr,
                "Branch & Bound depth: %d\nNodes processed: %d\nSimplex pivots: %d\n",
                get_max_level(lp), get_total_nodes(lp), get_total_iter(lp));
    } else if (result == INFEASIBLE) {
        if (PRINT_SOLUTION >= 1)
            printf("This problem is infeasible\n");
    } else if (result == UNBOUNDED) {
        if (PRINT_SOLUTION >= 1)
            printf("This problem is unbounded\n");
    } else if (result == FAILURE) {
        if (PRINT_SOLUTION >= 1)
            printf("lp_solve failed\n");
    }

    delete_lp(lp);

#if defined FORTIFY
    Fortify_LeaveScope();
#endif

    return (result);
}
