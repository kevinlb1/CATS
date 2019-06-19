/*
  Main header file of the LP_SOLVE toolkit.

  Original by Jeroen Dirks, 21-2-95
  Maintained by Michel Berkelaar

  Starting at version 3.0, LP_Solve is released under the LGPL license.
  For full information, see the enclosed file LGPL.txt.

  See CHANGELOG file for a tracking of changes

 */

#ifndef __LPKIT_H__
#define __LPKIT_H__

/* let's please C++ users */
#ifdef __cplusplus
extern "C" {
#endif

#if defined _WINDOWS
#include <windows.h>
#else
#if defined WINAPI
#undef WINAPI
#endif
#define WINAPI
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hash.h"
#include "fortify.h"

#ifndef NULL
#define NULL          0
#endif

#define FALSE           0
#define TRUE            1
#define AUTOMATIC       2

#define ISREAL          0
#define ISINTEGER       1
#define ISSEMI          2
#define ISSOS           4
#define ISSOSTEMPINT    8

#define ROWNAMEMASK     "r_%d"
#define COLNAMEMASK     "v_%d"

    /* REPORT defines */
#define CRITICALSTOP    0
#define CRITICAL        1
#define SEVERE          2
#define IMPORTANT       3
#define NORMAL          4
#define DETAILED        5
#define FULL            6

    /* SOS constraint defines */
#define SOS1            1
#define SOS2            2
#define SOS3           -1
#define SOS_START_SIZE 10     /* start size of SOS_list array; realloced if needed */

    /* MESSAGE defines */
#define MSG_NONE             0
#define MSG_PRESOLVE         1
#define MSG_ITERATION        2
#define MSG_INVERT           4
#define MSG_LPFEASIBLE       8
#define MSG_LPEQUAL         16
#define MSG_LPBETTER        32
#define MSG_MILPFEASIBLE    64
#define MSG_MILPEQUAL      128
#define MSG_MILPBETTER     256
#define MSG_MILPSTRATEGY   512

#ifndef DEFNUMINV
#define DEFNUMINV      30
#endif

#ifndef HASHSIZE
#define HASHSIZE       10007  /* A prime number! */
#endif

#ifndef INITIAL_MAT_SIZE
#define INITIAL_MAT_SIZE 10000
#endif

#ifndef DELTACOLALLOC
#define DELTACOLALLOC  50
#endif

    /* Improvement defines */
#define IMPROVE_NONE    0
#define IMPROVE_FTRAN   1
#define IMPROVE_BTRAN   2
#define IMPROVE_FULL    (IMPROVE_FTRAN + IMPROVE_BTRAN);

    /* Scaling defines */
#define MMSCALING       0
#define GEOSCALING      1
#define POWERSCALE      2
#define CURTISREIDSCALE 4
#define LAGRANGESCALE   8
#define INTEGERSCALE   16

    /* B&B strategies */
#define FIRST_NI 0
#define RAND_NI  1

#define FIRST_SELECT    0
#define RAND_SELECT 1
#define WORST_SELECT    2
#define BEST_SELECT     3
#define MEDIAN_SELECT   4
#define GREEDY_SELECT   5
    /* unused ...
    #define AUTO_SELECT     6
    #define USER_SELECT     7
     */

    /* solve status values */
#define UNKNOWN        -5
#define NOTRUN         -4
#define USERABORT      -3
#define TIMEOUT        -2
#define IGNORED        -1
#define OPTIMAL      0
#define MILP_FAIL    1
#define INFEASIBLE   2
#define UNBOUNDED    3
#define FAILURE      4
#define RUNNING      5

    /* lag_solve extra status values */
#define FEAS_FOUND    6
#define NO_FEAS_FOUND  7
#define BREAK_BB        8

    /* Status values internal to the solver */
#define SWITCH_TO_PRIMAL         9
#define SINGULAR_BASIS          10
#define LOST_PRIMAL_FEASIBILITY 11
#define OUT_OF_MEMORY           12

#define LE      0
#define EQ      1
#define GE      2
#define OF      3

    /* Solver parameters and tolerances */
#ifndef RESIZEFACTOR
#define RESIZEFACTOR       1.25
#endif

#ifndef DEFNUMINV
#define DEFNUMINV            30 /* maximum number of pivots before reinversion */
#endif
#ifndef DEF_MAXRELAX
#define DEF_MAXRELAX          4 /* maximum number of non-BB relaxations in MILP */
#endif

#ifndef DEF_MAXSINGULARITIES
#define DEF_MAXSINGULARITIES  10 /* maximum number of singularities in inversion */
#endif

#ifndef LAG_SINGULARLIMIT
#define LAG_SINGULARLIMIT     5 /* Number of times the objective does not change
                                    before it is assumed that the Lagrangean constraints
                                    are non-binding, and therefore impossible to converge;
                                    upper iteration limit is divided by this threshold */
#endif

#ifndef DEF_INFINITE
#define DEF_INFINITE  1e24     /* limit for dynamic range */
#endif

#ifndef DEF_NEGRANGE
#define DEF_NEGRANGE     0     /* downward limit for expanded variable range
                                   before the variable is split into positive and
                                   negative components */
#endif

#ifndef DEF_EPSPIVOT
#define DEF_EPSPIVOT  1e-5     /* 5 pivot reject (try others first)  */
#endif

#ifndef DEF_EPSEL
#ifdef ORGPARAM
#define DEF_EPSEL     1e-8    /* for rounding other values (vectors) to 0 */
#else
#define DEF_EPSEL     1e-11    /* 10 for rounding other values (vectors) to 0 */
#endif
#endif

#ifndef DEF_EPSB
#ifdef ORGPARAM
#define DEF_EPSB   5.01e-7    /* 8 for rounding RHS values to 0;
                                   determine infeasibility basis */
#else
#define DEF_EPSB    5.01e-9    /* for rounding RHS values to 0 determine
				   infeasibility basis */
#endif
#endif

#ifndef DEF_EPSD
#ifdef ORGPARAM
#define DEF_EPSD      1e-6    /* for rounding reduced costs to zero */
#else
#define DEF_EPSD      1e-9    /* 7 for rounding reduced costs to 0 */
#endif
#endif

#ifndef RANDSCALE
#define RANDSCALE     100      /* Randomization scaling range */
#endif

#ifndef DEF_PERTURB
#define DEF_PERTURB   1e-5     /* Perturbation scalar for degenerative problems;
                                   must at least be RANDSCALE greater than DEF_EPSB */
#endif

#ifndef SOLUTIONEPS
#define SOLUTIONEPS   1e-5     /* Margin of error for solution bounds */
#endif

#ifndef DEF_EPSILON
#ifdef ORGPARAM
#define DEF_EPSILON  1e-3     /* to determine if a float value is integer */
#else
#define DEF_EPSILON  1e-6     /* to determine if a float value is integer */
#endif
#endif

#ifndef DEF_LAGACCEPT
#define DEF_LAGACCEPT 1e-3     /* Default Lagrangian convergence acceptance criterion */
#endif

#ifndef MINSCALAR
#define MINSCALAR    1e-10     /* Smallest allowed scaling adjustment */
#endif

#ifndef MAXSCALAR
#define MAXSCALAR     1e+10    /* Largest allowed scaling adjustment */
#endif

#define MINTIMEPIVOT  5e-2      /* Minimum time per pivot for pivot optimization purposes */
#define SCALINGEPS    1e-2      /* Scaling convergence criterion */

#define my_abs(x)       ((x) < 0 ? -(x) : (x))
#define my_min(x, y)    ((x) < (y) ? (x) : (y))
#define my_max(x, y)    ((x) > (y) ? (x) : (y))
#define my_if(t, x, y)  ((t) ? (x) : (y))

#define MAX_WARN_COUNT 20

#ifdef CHECK
#define my_round(val, eps) { \
    REAL absv; \
        absv = ((val) < 0 ? -(val) : (val)); \
        if(absv < (eps)) \
          val = 0; \
    if(Warn_count < MAX_WARN_COUNT) \
      { \
        if(absv > 0.5 * (eps) && absv < 2 * (eps)) \
          { \
        Warn_count++; \
        report(NULL, NORMAL, \
            "Warning Value close to epsilon V: %e E: %e\n", \
            (double)absv, (double)(eps)); \
        if(Warn_count == MAX_WARN_COUNT) \
          { \
            report(NULL, NORMAL, \
                "*** Surpressing further rounding warnings\n"); \
          } \
          } \
      } \
}

#else
#define my_round(val,eps) if (((val) < 0 ? -(val) : (val)) < (eps)) val = 0;
#endif

#ifndef REAL /* to allow -DREAL=<float type> while compiling */
#define REAL double
#endif

#ifndef LREAL
#define LREAL long double
#endif

#ifndef RREAL
#define RREAL long double
#endif

#ifndef MYBOOL
    /*#define MYBOOL  unsigned short */
#define MYBOOL     unsigned char
#endif

#ifndef STATUS
#define STATUS   short
#endif

#define NAMELEN 25
#define MAXSTRL (NAMELEN-1)
#define STD_ROW_NAME_PREFIX "r_"

#define MALLOC(ptr, nr)\
  ((((nr) == 0) || ((ptr = malloc((size_t)((nr) * sizeof(*ptr)))) == NULL)) ? \
   (void *) report(NULL, CRITICAL, "malloc of %d bytes failed on line %d of file %s",\
           (nr) * sizeof(*ptr), __LINE__, __FILE__), (ptr = (void *) 0) : \
   ptr\
  )

#define CALLOC(ptr, nr)\
  ((((nr) == 0) || ((ptr = calloc((size_t)(nr), sizeof(*ptr))) == NULL)) ? \
   (void *) report(NULL, CRITICAL, "calloc of %d bytes failed on line %d of file %s",\
           (nr) * sizeof(*ptr), __LINE__, __FILE__), (ptr = (void *) 0) : \
   ptr\
  )

#define REALLOC(ptr, nr)\
  ((((nr) == 0) || ((ptr = realloc(ptr, (size_t)((nr) * sizeof(*ptr)))) == NULL)) ? \
   (void *) report(NULL, CRITICAL, "realloc of %d bytes failed on line %d of file %s",\
           (nr) * sizeof(*ptr), __LINE__, __FILE__), (ptr = (void *) 0) : \
   ptr\
  )

#define FREE(ptr) if (ptr != NULL) {free(ptr), ptr = NULL;} else

#define MALLOCCPY(nptr, optr, nr)\
  (MALLOC(nptr, nr),(nptr != NULL) ? memcpy(nptr, optr, (size_t)((nr) * sizeof(*optr))) : 0, nptr)

#define MEMCPY(nptr, optr, nr)\
  memcpy(nptr, optr, (size_t)((nr) * sizeof(*optr)));

    typedef char nstring[NAMELEN];

    typedef struct _matrec {
        int row_nr;
        REAL value;
    } matrec;

    typedef struct _column {
        int row;
        REAL value;
        struct _column *next;
    } column;

    typedef struct _constraint_name {
        char name[NAMELEN];
        int row;
        struct _constraint_name *next;
    } constraint_name;

    typedef struct _bound {
        REAL upbo;
        REAL lowbo;
    } bound;

    typedef struct _tmp_store_struct {
        nstring name;
        int row;
        REAL value;
        REAL rhs_value;
        short relat;
    } tmp_store_struct;

    typedef struct _rside /* contains relational operator and rhs value */ {
        int row;
        REAL value;
        REAL range_value;
        struct _rside *next;
        short relat;
        short range_relat;
    } rside;

    /* SOS storage structure */

#define LINEARSEARCH 0

    typedef struct _SOSrec {
        int tagorder;
        char *name;
        short type;
        int size;
        int priority;
        int *members;
        REAL *weights;
        int *membersSorted;
        int *membersMapped;
    } SOSrec;

    /* Prototypes for call-back functions*/
    typedef int WINAPI abortfunc(void *lp, void *userhandle);
    typedef void WINAPI logfunc(void *lp, void *userhandle, char *buf);
    typedef void * WINAPI msgfunc(void *lp, void *userhandle, int message);

    /* fields indicated with ## may be modified directly */

    /* pointers will have their array size in the comments */

    typedef struct _lprec {
        char *lp_name; /* the name of the lp */
        short verbose; /* ## Verbose flag */
        MYBOOL print_duals; /* ## Set TRUE to print duals in PrintSolution */
        MYBOOL print_sol; /* ## Set TRUE to print optimal solution */
        MYBOOL debug; /* ## Set TRUE to print extra debug information */
        short print_at_invert; /* ## Print information at every reinversion */
        MYBOOL trace; /* ## Print information on simplex progression */
        MYBOOL anti_degen; /* ## Set TRUE to do perturbations to avoid cycling */
        MYBOOL do_presolve; /* ## Set TRUE to perform matrix presolving */

        int rows; /* Nr of constraint rows in the problem */
        int rows_alloc; /* The allocated memory for Rows sized data */
        int columns; /* The number of columns (= variables) */
        int columns_alloc;
        int sum; /* The size of the variables + the slacks */
        int sum_alloc;

        MYBOOL names_used; /* Flag to indicate if names for rows and columns are used */

        char **row_name; /* rows_alloc+1 */
        char **col_name; /* columns_alloc+1 */

        /* Row[0] of the sparce matrix is the objective function */

        int non_zeros; /* The number of elements in the sparse matrix*/
        int mat_alloc; /* The allocated size for matrix sized structures */
        matrec *mat; /* mat_alloc :The sparse matrix */
        int *col_end; /* columns_alloc+1 :Cend[i] is the index of the
                                   first element after column i.
                		   column[i] is stored in elements
		                   col_end[i-1] to col_end[i]-1 */
        int *col_no; /* mat_alloc :From Row 1 on, col_no contains the
                		   column nr. of the nonzero elements, row by row */
        MYBOOL row_end_valid; /* true if row_end & col_no are valid */
        int *row_end; /* rows_alloc+1 :row_end[i] is the index of the
                		   first element in Colno after row i */
        REAL *orig_rh; /* rows_alloc+1 :The RHS after scaling & sign
		                   changing, but before 'Bound transformation' */
        REAL *rh; /* rows_alloc+1 :As orig_rh, but after Bound transformation */
        RREAL *rhs; /* rows_alloc+1 :The RHS of the current simplex tableau */
        MYBOOL *must_be_int; /* sum_alloc+1 :TRUE if variable must be Integer */
        REAL *orig_upbo; /* sum_alloc+1 :Bound before transformations */
        REAL *orig_lowbo; /*  "       "                   */
        REAL *upbo; /*  " " :Upper bound after transformation & B&B work */
        REAL *lowbo; /*  "       "  :Lower bound after transformation & B&B work */

        MYBOOL basis_valid; /* TRUE is the basis is still valid */
        int *bas; /* rows_alloc+1 :The basis column list */
        MYBOOL *basis; /* sum_alloc+1 : basis[i] is TRUE if the column is in the basis */
        MYBOOL *lower; /*  "       "  :TRUE if the variable is at its
		                   lower bound (or in the basis), it is FALSE
                		   if the variable is at its upper bound */

        MYBOOL eta_valid; /* TRUE if current Eta structures are valid */
        int eta_alloc; /* The allocated memory for Eta non-zero values */
        int eta_size; /* The number of Eta columns */
        int num_inv; /* Number of pivots since last refactorization */
        int max_num_inv; /* ## Number of pivots between reinversions of the ETA matrix */
        REAL *eta_value; /* eta_alloc : Structure containing values of Eta */
        int *eta_row_nr; /*  "     "  : Structure containing row indexes of Eta */
        int *eta_col_end; /* rows_alloc + MaxNumInv : eta_col_end[i] is
                		   the start index of the next Eta column */

        MYBOOL bb_rule; /* ## Set rule for selecting B&B variables:
                                      FIRST_SELECT : Lowest indexed non-integer column
                                      RAND_SELECT  : Random non-integer column
                                      WORST_SELECT : Largest deviation from an integer value
                                      MEDIAN_SELECT: Median value deviation from an integer value */

        REAL obj_bound; /* ## Set initial "at least better than" guess for objective function
                                   (can in particular speed up B&B iterations) */
        int iter; /* Number of iterations in the simplex solver (LP) */
        int total_iter; /* Number of iterations (B&B) (Integer LP) */
        int max_level; /* The Deepest B&B level of the last solution */
        int total_nodes; /* Total number of nodes processed in B&B */
        REAL *solution; /* sum_alloc+1 : Solution array of the next to optimal LP,
                                   Index   0           : Objective function value,
                                   Indeces 1..rows     : Slack variable values,
                		   Indeced rows+1..sum : Variable values */
        REAL *best_solution; /* sum_alloc+1 : Solution array of optimal 'Integer' LP,
                                   structured as the solution array above */
        REAL *duals; /* sum_alloc+1 :The dual variables/reduced costs of the last LP */

        MYBOOL maximise; /* TRUE if the goal is to maximise the objective function */
        MYBOOL floor_first; /* FALSE: B&B does ceiling bound first
                                   TRUE:  B&B does floor bound first
                                   AUTOMATIC: B&B desides which branch to take first */
        MYBOOL *ch_sign; /* rows_alloc+1 :TRUE if the Row in the matrix has changed sign
                                   (a`x > b, x>=0) is translated to
	          		    s + -a`x = -b with x>=0, s>=0) */

        MYBOOL scaling_used; /* TRUE if scaling is used */
        MYBOOL columns_scaled; /* TRUE if the columns are scaled too, Only use
                	 	   if all variables are non-integer */
        REAL *scale; /* sum_alloc+1:0..Rows the scaling of the Rows,
                		   Rows+1..Sum the scaling of the columns */

        int nr_lagrange; /* Nr. of Langrangian relaxation constraints */
        REAL **lag_row; /* NumLagrange, columns+1:Pointer to pointer of rows */
        REAL *lag_rhs; /* NumLagrange :Pointer to pointer of Rhs */
        REAL *lambda; /* NumLagrange :Lambda Values */
        MYBOOL *lag_con_type; /* NumLagrange :TRUE if constraint type EQ */
        REAL lag_bound; /* The Lagrangian lower bound */

        MYBOOL valid; /* Has this lp pased the 'test' */
        REAL infinite; /* ## limit for dynamic range */
        REAL epsilon; /* ## to determine if a float value is integer */
        REAL epsb; /* ## for rounding RHS values to 0/infeasibility */
        REAL epsd; /* ## for rounding reduced costs to zero */
        REAL epsel; /* ## for rounding other values (vectors) to 0 */
        hashstruct *rowname_hashtab; /* hash table to store row names */
        hashstruct *colname_hashtab; /* hash table to store column names */

        REAL *dualsfrom; /* sum_alloc+1 :The sensitivity on dual variables/reduced costs of the
                		   last LP */
        REAL *dualstill; /* sum_alloc+1 :The sensitivity on dual variables/reduced costs of the
                		   last LP */
        REAL *objfrom; /* columns_alloc+1 :The sensitivity on object function of the
		                   last LP */
        REAL *objtill; /* columns_alloc+1 :The sensitivity on object function of the
                                   last LP */

        int orig_rows; /* Number of actual problem rows before deletions */
        int orig_columns; /* Number of actual problem columns before working columns */
        short spx_status; /* Simplex solver feasibility/mode code */
        short lag_status; /* Extra status variable for lag_solve */
        int solutioncount; /* number of equal-valued solutions found (up to solutionlimit) */
        int solutionlimit; /* upper number of equal-valued solutions kept track of */
        int num_refact; /* Number of times the basis was refactored */
        MYBOOL scalemode; /* ## Set 0 for max-min, 1 for geometric */
        MYBOOL improve; /* ## Set to non-zero for iterative improvement */
        MYBOOL lag_trace; /* ## Print information on Lagrange progression */
        MYBOOL piv_rule; /* ## Set rule for selecting row and column entering/leaving */
        MYBOOL break_at_first; /* ## Set TRUE to stop at first feasible solution */
        short bb_floorfirst; /* ## Set TRUE for B&B to set variables to floor bound first;
                                      conversely with FALSE, the ceiling value is set first */
        REAL break_at_value; /* ## Set value for the objective function deemed
                                      to be sufficiently good in an integer problem */
        int int_count; /* Number of integers required */
        int *var_is_free; /* columns+1: Index of twin variable if variable is free */

        REAL Extrad;
        MYBOOL doIterate; /* Perform iteration at next opportunity */
        MYBOOL doInvert; /* Force basis reinversion immediateluy */
        MYBOOL justInverted; /* The basis was recently reinverted */
        MYBOOL wasprocessed; /* The solve preprocessing was performed */

        MYBOOL Break_bb; /* Solver working variable */
        int Level; /* Solver B&B working variable (recursion depth) */

        REAL lag_accept; /* The Lagrangial convergence acceptance criterion */

        REAL negrange; /* ## limit for negative variable range */

        REAL epsperturb; /* ## perturbation scalar */

        REAL epspivot; /* ## Pivot reject tolerance (try others first) */

        /* Time/timer variables */
        long sectimeout;
        double timestart;
        double timeend;
        double time_refactstart; /* Time since start of last refactorization-pivots cyle */

        /* Message processing callbacks */
        abortfunc *abort;
        void *aborthandle;
        logfunc *writelog;
        void *loghandle;
        logfunc *debuginfo;
        msgfunc *usermessage;
        int msgmask;
        void *msghandle;
#if 0
        int *var_to_orig; /* sum_alloc+1 : Mapping of variables from solution to
				   best_solution to account for removed variables and
				   rows during presolve; a non-positive value indicates
                                   that the constraint or variable was removed */
        int *orig_to_var;
#endif
        int sc_count; /* Number of semi-continuous variables */
        int sos_alloc; /* Size allocated to specially ordered sets (SOS1, SOS2...) */
        int sos_count; /* Number of specially ordered sets (SOS1, SOS2...) */
        int sos_ints; /* Number of integers in SOS'es above */
        REAL *var_is_sc; /* sum_columns+1 : TRUE if variable is semi-continuous;
                                   value replaced by conventional lower bound during solve */
        SOSrec **sos_list; /* Array of pointers to SOS lists */
        int sos_vars; /* Number of variables in the sos_nodes list */
        int *sos_priority; /* Priority-sorted list of variables (no duplicates) */
    } lprec;



    /* function interface for the user */

    void lp_solve_version(int *majorversion,
            int *minorversion,
            int *release,
            int *build);
    void set_magic(int code, int param);

    lprec *make_lp(int rows, int columns);
    /* create and initialise a lprec structure
       defaults:
       Empty (Rows * Columns) matrix,
       Minimize the objective function
       constraints all type <=
       Upperbounds all Infinite
       no integer variables
       floor first in B&B
       no scaling
       default basis */

    lprec *read_LP(char *input, short verbose, char *lp_name);
    lprec *read_lp_file(FILE *input, short verbose, char *lp_name);
    /* create and read an .lp file from input (input must be open) */

    void delete_lp(lprec *lp);
    /* Remove problem from memory */

    lprec *copy_lp(lprec *lp);
    /* copy a lp structure */

    int set_mat(lprec *lp, int row, int column, REAL value);
    /* fill in element (Row,Column) of the matrix
       Row in [0..Rows] and Column in [1..Columns] */

    int set_obj_fn(lprec *lp, REAL *row);
    /* set the objective function (Row 0) of the matrix */
    int str_set_obj_fn(lprec *lp, char *row);
    /* The same, but with string input */

    int add_constraint(lprec *lp, REAL *row, short constr_type, REAL rh);
    /* Add a constraint to the problem,
       row is the constraint row,
       rh is the right hand side,
       constr_type is the type of constraint (LE (<=), GE(>=), EQ(=)) */
    int str_add_constraint(lprec *lp, char *row_string, short constr_type, REAL rh);
    /* The same, but with string input */

    int del_constraint(lprec *lp, int del_row);
    /* Remove constrain nr del_row from the problem */

    int add_lag_con(lprec *lp, REAL *row, short con_type, REAL rhs);
    /* add a Lagrangian constraint of form Row' x contype Rhs */
    int str_add_lag_con(lprec *lp, char *row, short con_type, REAL rhs);
    /* The same, but with string input */

    int add_column(lprec *lp, REAL *column);
    /* Add a Column to the problem */
    int str_add_column(lprec *lp, char *col_string);
    /* The same, but with string input */

    int del_column(lprec *lp, int column);
    /* Delete a column */

    int set_upbo(lprec *lp, int column, REAL value);
    /* Set the upperbound of a variable */

    int set_lowbo(lprec *lp, int column, REAL value);
    /* Set the lowerbound of a variable */

    int set_uprange(lprec *lp, int row, REAL value);
    /* Set the upper range of a constraint */

    int set_lowrange(lprec *lp, int row, REAL value);
    /* Set the lower range of a constraint */

    int set_int(lprec *lp, int column, short must_be_int);
    /* Set the type of variable, if must_be_int = TRUE then the variable must be integer */

    int is_int(lprec *lp, int column);

    int set_semicont(lprec *lp, int column, short must_be_sc);

    int is_semicont(lprec *lp, int column);

    void set_verbose(lprec *lp, short verbose);
    short get_verbose(lprec *lp);

    void set_timeout(lprec *lp, long sectimeout);
    long get_timeout(lprec *lp);

    void set_print_duals(lprec *lp, MYBOOL print_duals);
    MYBOOL is_print_duals(lprec *lp);

    void set_print_sol(lprec *lp, MYBOOL print_sol);
    MYBOOL is_print_sol(lprec *lp);

    void set_debug(lprec *lp, MYBOOL debug);
    MYBOOL is_debug(lprec *lp);

    void set_print_at_invert(lprec *lp, short print_at_invert);
    short is_print_at_invert(lprec *lp);

    void set_trace(lprec *lp, MYBOOL trace);
    MYBOOL is_trace(lprec *lp);

    void set_anti_degen(lprec *lp, MYBOOL anti_degen);
    MYBOOL is_anti_degen(lprec *lp);

    void set_do_presolve(lprec *lp, MYBOOL do_presolve);
    MYBOOL is_do_presolve(lprec *lp);

    void set_max_num_inv(lprec *lp, int max_num_inv);
    int get_max_num_inv(lprec *lp);

    void set_bb_rule(lprec *lp, MYBOOL bb_rule);
    MYBOOL get_bb_rule(lprec *lp);

    void set_obj_bound(lprec *lp, REAL obj_bound);
    REAL get_obj_bound(lprec *lp);

    void set_floor_first(lprec *lp, MYBOOL floor_first);
    MYBOOL get_floor_first(lprec *lp);

    void set_infinite(lprec *lp, REAL infinite);
    REAL get_infinite(lprec *lp);

    void set_epsilon(lprec *lp, REAL epsilon);
    REAL get_epsilon(lprec *lp);

    void set_epsb(lprec *lp, REAL epsb);
    REAL get_epsb(lprec *lp);

    void set_epsd(lprec *lp, REAL epd);
    REAL get_epsd(lprec *lp);

    void set_epsel(lprec *lp, REAL epsel);
    REAL get_epsel(lprec *lp);

    void set_scalemode(lprec *lp, MYBOOL scalemode);
    MYBOOL get_scalemode(lprec *lp);

    void set_improve(lprec *lp, MYBOOL improve);
    MYBOOL is_improve(lprec *lp);

    void set_lag_trace(lprec *lp, MYBOOL lag_trace);
    MYBOOL is_lag_trace(lprec *lp);

    void set_piv_rule(lprec *lp, MYBOOL piv_rule);
    MYBOOL get_piv_rule(lprec *lp);

    void set_break_at_first(lprec *lp, MYBOOL break_at_first);
    MYBOOL is_break_at_first(lprec *lp);

    void set_bb_floorfirst(lprec *lp, short bb_floorfirst);
    short is_bb_floorfirst(lprec *lp);

    void set_break_at_value(lprec *lp, REAL break_at_value);
    REAL get_break_at_value(lprec *lp);

    void set_negrange(lprec *lp, REAL negrange);
    REAL get_negrange(lprec *lp);

    void set_epsperturb(lprec *lp, REAL epsperturb);
    REAL get_epsperturb(lprec *lp);

    void set_epspivot(lprec *lp, REAL epspivot);
    REAL get_epspivot(lprec *lp);

    int get_max_level(lprec *lp);

    int get_total_nodes(lprec *lp);

    int get_total_iter(lprec *lp);

    REAL get_objective(lprec *lp);

    int get_variables(lprec *lp, REAL *var);

    int get_ptr_variables(lprec *lp, REAL **var);

    int get_constraints(lprec *lp, REAL *constr);

    int get_ptr_constraints(lprec *lp, REAL **constr);

    int get_sensitivity_rhs(lprec *lp, REAL *duals, REAL *dualsfrom, REAL *dualstill);

    int get_ptr_sensitivity_rhs(lprec *lp, REAL **duals, REAL **dualsfrom, REAL **dualstill);

    int get_sensitivity_obj(lprec *lp, REAL *objfrom, REAL *objtill);

    int get_ptr_sensitivity_obj(lprec *lp, REAL **objfrom, REAL **objtill);

    int get_Nrows(lprec *lp);

    int get_Ncolumns(lprec *lp);

    /* Add SOS constraints */
    int add_SOS(lprec *lp, char *name, short sostype, int priority, int count, int *sosvars, REAL *weights);
    int is_SOS_var(lprec *lp, int column);

    int set_lp_name(lprec *lp, char *name);

    int set_rh(lprec *lp, int row, REAL value);
    REAL get_rh(lprec *lp, int row);
    /* Set the right hand side of a constraint row */

    int set_rh_range(lprec *lp, int row, REAL deltavalue);
    REAL get_rh_range(lprec *lp, int row);
    /* Set the RHS range; i.e. the lower and upper bounds of a constraint row */

    void set_rh_vec(lprec *lp, REAL *rh);
    /* Set the right hand side vector */
    int str_set_rh_vec(lprec *lp, char *rh_string);
    /* The same, but with string input */

    void set_maxim(lprec *lp);
    /* maximise the objective function */

    void set_minim(lprec *lp);
    /* minimize the objective function */

    int set_constr_type(lprec *lp, int row, short con_type);
    /* Set the type of constraint in row Row (LE, GE, EQ) */

    int get_constr_type(lprec *lp, int row);
    /* get the type of constraint in row Row (LE, GE, EQ) */

    int set_row_name(lprec *lp, int row, char *new_name);
    /* Set the name of a constraint row, make sure that the name has < 25 characters */

    char *get_row_name(lprec *lp, int row);

    int set_col_name(lprec *lp, int column, char *new_name);
    /* Set the name of a varaible column, make sure that the name has < 25 characters */

    char *get_col_name(lprec *lp, int column);

    REAL scale(lprec *lp, REAL *myrowscale, REAL *mycolscale);
    /* Automatic scaling of the problem */

    REAL auto_scale(lprec *lp);
    /* Automatic scaling of the problem */

    int scaleCR(lprec *lp);
    /* Curtis-Reid scaling */

    void unscale(lprec *lp);
    /* Remove all scaling from the problem */

    /* Set/Get basis for a re-solved system */ /* Added by KE */
    void set_basis(lprec *lp, int *bascolumn);
    void get_basis(lprec *lp, int *bascolumn);

    int solve(lprec *lp);
    /* Solve the problem */

    int lag_solve(lprec *lp, REAL start_bound, int num_iter, short verbose);
    /* Do NumIter iterations with Lagrangian relaxation constraints */

    void reset_basis(lprec *lp);
    /* Reset the basis of a problem, can be usefull in case of degeneracy - JD */

    REAL mat_elm(lprec *lp, int row, int column);
    /* get a single element from the matrix */

    int get_row(lprec *lp, int row_nr, REAL *row);
    /* fill row with the row row_nr from the problem */

    int get_column(lprec *lp, int col_nr, REAL *column);
    /* fill column with the column col_nr from the problem */

    int get_reduced_costs(lprec *lp, REAL *rc);
    /* get the reduced costs vector */

    int is_feasible(lprec *lp, REAL *values);
    /* returns TRUE if the vector in values is a feasible solution to the lp */

    int column_in_lp(lprec *lp, REAL *column);
    /* returns TRUE if column is already present in lp. (Does not look at bounds
       and types, only looks at matrix values */

    lprec *read_MPS(char *input, short verbose);
    lprec *read_mps(FILE *input, short verbose);
    /* read a MPS file */

    int write_mps(lprec *lp, char *output);
    int write_MPS(lprec *lp, FILE *output);
    /* write a MPS file to output */

    int write_lp(lprec *lp, char *output);
    int write_LP(lprec *lp, FILE *output);
    /* write a LP file to output */

    void print_lp(lprec *lp);
    /* Print the current problem, only usefull in very small (test) problems.
      Shows the effect of scaling */

    void print_objective(lprec *lp);
    /* Print the value of the objective */

    void print_solution(lprec *lp);
    /* Print the values of the variables */

    void print_constraints(lprec *lp);
    /* Print the values of the constraints */

    void print_duals(lprec *lp);
    /* Print the dual variables of the solution */

    void print_scales(lprec *lp);
    /* If scaling is used, print the scaling factors */

    int print_file(char *filename);
    /* file where results are printed to. Default stderr. If NULL then back stderr */

    /* Allow the user to define an interruption callback function */
    void put_abortfunc(lprec *lp, abortfunc newabort, void *aborthandle);

    /* Allow the user to define a logging function */
    void put_logfunc(lprec *lp, logfunc newlog, void *loghandle);

    /* Allow the user to define an event-driven message/reporting */
    void put_msgfunc(lprec *lp, msgfunc newmsg, void *msghandle);


    /* functions used internaly by the lp toolkit */
    lprec *make_lpext(int rows, int columns, int non_zeros, int mat_alloc, char *lp_name);
    int report(lprec *lp, short level, char *format, ...);
    int inc_mat_space(lprec *lp, int max_extra);
    int inc_row_space(lprec *lp);
    int inc_col_space(lprec *lp);
    void unscale_columns(lprec *lp);
    void btran(lprec *lp, REAL *row, REAL roundzero);
    int presolve(lprec *lp);
    MYBOOL isvalid(lprec *lp);
    REAL get_mat(lprec *lp, int row, int column);
    int set_upbo(lprec *lp, int column, REAL value);
    REAL get_upbo(lprec *lp, int column);
    int set_lowbo(lprec *lp, int column, REAL value);
    REAL get_lowbo(lprec *lp, int column);
    int set_bounds(lprec *lp, int column, REAL lower, REAL upper);
    int preprocess(lprec *lp);
    void postprocess(lprec *lp);
    int set_matrix(lprec *lp, int Row, int Column, REAL Value, MYBOOL doscale);
    int SOS_infeasible(lprec *lp, int sosindex);
    int SOS_get_type(lprec *lp, int sosindex);
    int SOS_is_member(lprec *lp, int sosindex, int column);
    int SOS_fix_unmarked(lprec *lp, int variable, int sosindex, REAL *bound, REAL value, MYBOOL isupper, int *diffcount);
    int append_SOSrec(lprec *lp, SOSrec *SOS, int size, int *variables, REAL *weights);
    MYBOOL SOS_set_marked(lprec *lp, int sosindex, int column, MYBOOL islive);
    REAL get_mat_raw(lprec *lp, int row, int column);
    MYBOOL SOS_unmark(lprec *lp, int sosindex, int column, MYBOOL islive);
    MYBOOL SOS_is_member_of_type(lprec *lp, int column, short sostype);
    MYBOOL SOS_is_active(lprec *lp, int sosindex, int column);
    int SOS_is_satisfied(lprec *lp, int sosindex, REAL *solution);
    MYBOOL SOS_can_mark(lprec *lp, int sosindex, int column);

    /* define yyparse() to make compilers happy. There should be some system
       include file for this */
    int yyparse(void);
    void yyerror(char *);
    void check_decl(int);
    int yywrap();

#ifdef __cplusplus
};
#endif

#endif
