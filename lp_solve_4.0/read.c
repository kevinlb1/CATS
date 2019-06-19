/*
   ============================================================================
   NAME : read.c

   PURPOSE : translation of lp-problem and storage in sparse matrix

   SHORT : Subroutines for yacc program to store the input in an intermediate
   data-structure. The yacc and lex programs translate the input.  First the
   problemsize is determined and the date is read into an intermediate
   structure, then readinput fills the sparse matrix.

   USAGE : call yyparse(); to start reading the input.  call readinput(); to
   fill the sparse matrix.
   ============================================================================
   Rows : contains the amount of rows + 1. Rows-1 is the amount of constraints
   (no bounds) Rows also contains the rownr 0 which is the objective function

   Columns : contains the amount of columns (different variable names found in
   the constraints)

   Nonnuls : contains the amount of nonnuls = sum of different entries of all
   columns in the constraints and in the objectfunction

   Hash_tab : contains all columnnames on the first level of the structure the
   row information is kept under each column structure in a linked list (also
   the objective funtion is in this structure) Bound information is also
   stored under under the column name

   First_rside : points to a linked list containing all relational operators
   and the righthandside values of the constraints the linked list is in
   reversed order with respect to the rownumbers
   ============================================================================ */
#include <string.h>
#include <limits.h>
#include "lpkit.h"
#include "lpglob.h"

short *relat;
int Verbose;
constraint_name *First_constraint_name;
rside *First_rside, *rs;
tmp_store_struct tmp_store;
short Ignore_decl;
hashstruct *Hash_tab;

/*
 * error handling routine for yyparse()
 */
void yyerror(char *string) {
    report(NULL, CRITICALSTOP, string);
}

int yywrap() /* supply a default yywrap() function */ {
    return (1);
}

void check_decl(int within_int_decl) {
    if (within_int_decl) {
        Ignore_decl = FALSE;
    } else {
        if (Verbose >= NORMAL)
            report(NULL, NORMAL, "Unknown declaration specifier on line %d, ignored",
                yylineno);
        Ignore_decl = TRUE;
    }
}

void add_int_var(char *name) {
    hashelem *hp;

    if ((hp = findhash(name, Hash_tab)) == NULL) {
        if (Verbose >= NORMAL)
            report(NULL, NORMAL,
                "Unknown variable %s declared integer on line %d, ignored",
                name, yylineno);
    } else if (hp->must_be_int) {
        if (Verbose >= NORMAL)
            report(NULL, NORMAL, "Variable %s declared integer more than once on line %d",
                name, yylineno);
    } else
        hp->must_be_int = TRUE;
}

/*
 * initialisation of hashstruct and globals.
 */
int init_read(void) {
    int ok = FALSE;

    Rows = 0;
    Non_zeros = 0;
    Columns = 0;
    if (CALLOC(First_rside, 1) != NULL) {
        rs = First_rside;
        rs->value = rs->range_value = 0;
        /* first row (nr 0) is always the objective function */
        rs->relat = OF;
        rs->range_relat = -1;
        if ((Hash_tab = create_hash_table(HASHSIZE)) == NULL) {
            FREE(First_rside);
        } else
            ok = TRUE;
    }
    return (ok);
} /* init */

/*
 * searchs in column-list (p is pointer to first element of column-list)
 * for column->row = row.
 * getrow() returns a pointer to this column structure.
 * If not found a NULL-pointer is returned
 */
static column *getrow(column *p,
        int row) {
    for (; p != NULL; p = p->next)
        if (p->row == row)
            return (p);
    return (p);
} /* getrow */

/*
 * Creates a bound record.
 * Set lowbo = 0 and upbo = Infinite
 *
 */
static bound *create_bound_rec(void) {
    bound *bp;

    if (CALLOC(bp, 1) != NULL) {
        bp->upbo = DEF_INFINITE;
        bp->lowbo = 0;
    }
    return (bp);
} /* create_bound_rec */

/*
 * clears the tmp_store variable after all information has been copied
 */
void null_tmp_store(void) {
    tmp_store.value = 0;
    tmp_store.rhs_value = 0;
}

/*
 * variable : pointer to text array with name of variable
 * row      : the rownumber of the constraint
 * value    : value of matrixelement
 *            A(row, variable).
 * Sign     : (global)  determines the sign of value.
 * store()  : stores value in matrix
 *	      A(row, variable). If A(row, variable) already contains data,
 *	      value is added to the existing value.
 */
static int store(char *variable,
        int row,
        REAL value) {
    hashelem *h_tab_p;
    column *col_p;

    if (value == 0) {
        if (Verbose >= NORMAL)
            report(NULL, NORMAL,
                "(store) Warning, variable %s has an effective coefficient of 0 on line %d. Ignored.",
                variable, yylineno);
        return (TRUE);
    }

    if ((h_tab_p = findhash(variable, Hash_tab)) == NULL) {
        if (((h_tab_p = puthash(variable, Hash_tab)) == NULL) ||
                (CALLOC(h_tab_p->col, 1) == NULL)
                ) return (FALSE);
        Columns++; /* counter for calloc of final array */

        Non_zeros++; /* for calloc of final arrays */
        h_tab_p->col->row = row;
        h_tab_p->col->value = value;
    } else if ((col_p = getrow(h_tab_p->col, row)) == NULL) {
        if (CALLOC(col_p, 1) == NULL)
            return (FALSE);
        Non_zeros++; /* for calloc of final arrays */
        col_p->value = value;
        col_p->row = row;
        col_p->next = h_tab_p->col;
        h_tab_p->col = col_p;
    } else
        col_p->value += value;
    return (TRUE);
} /* store */

/*
 * store relational operator given in yylex[0] in the rightside list.
 * Also checks if it constaint was a bound and if so stores it in the
 * boundslist
 */
void store_re_op(void) {
    short tmp_relat;

    switch (yytext[0]) {

        case '=':
            tmp_relat = EQ;
            break;

        case '>':
            tmp_relat = GE;
            break;

        case '<':
            tmp_relat = LE;
            break;

        default:
            report(NULL, CRITICALSTOP, "Error: unknown relational operator %s on line %d",
                    yytext, yylineno);
            return;
            break;
    }

    if (Lin_term_count > 1) /* it is not a bound */
        rs->relat = tmp_relat;
    else if (Lin_term_count == 0) { /* it is a range */
        if (rs == NULL) { /* range before row, already reported */
            if (Verbose >= CRITICAL)
                report(NULL, CRITICAL, "Error on line %d: range for undefined row.",
                    yylineno);
        } else if (rs->range_relat != -1) {
            if (Verbose >= CRITICAL)
                report(NULL, CRITICAL, "Error on line %d: There was already a range for this row.",
                    yylineno);
        } else if (tmp_relat == rs->relat) {
            if (Verbose >= CRITICAL)
                report(NULL, CRITICAL, "Error on line %d: relational operator for range is the same as relation operator for equation.",
                    yylineno);
            rs->range_relat = -2;
        } else
            rs->range_relat = tmp_relat;
    } else /* could be a bound */
        tmp_store.relat = tmp_relat;
} /* save_re_op */

/*
 * store RHS value in the rightside structure
 * if type = true then
 */
void rhs_store(REAL value) {
    if (Lin_term_count > 1) /* not a bound */
        rs->value += value;
    else if (Lin_term_count == 0) { /* a range */
        if (rs == NULL); /* range before row, already reported */
        else if (rs->range_relat < 0) /* was a bad range; ignore */;
        else if (((rs->relat == LE) && (rs->range_relat == GE) &&
                (rs->value < value)) ||
                ((rs->relat == GE) && (rs->range_relat == LE) &&
                (rs->value > value)) ||
                ((rs->relat == EQ) || (rs->range_relat == EQ))) {
            if (Verbose >= CRITICAL)
                report(NULL, CRITICAL, "Error on line %d: range restriction is conflicting with equation on previous line",
                    yylineno);
            rs->range_relat = -2;
        } else
            rs->range_value += value;
    } else /* a bound */
        tmp_store.rhs_value += value;
} /* RHS_store */

/*
 * store all data in the right place
 * count the amount of lineair terms in a constraint
 * only store in data-structure if the constraint is not a bound
 */
int var_store(char *var, int row, REAL value) {
    if (strlen(var) > MAXSTRL) {
        if (Verbose >= CRITICAL)
            report(NULL, CRITICAL, "Variable name '%s' too long, at most %d characters allowed",
                var, MAXSTRL);
        return (FALSE);
    }
    /* also in a bound the same var name can occur more than once. Check for
       this. Don't increment Lin_term_count */

    if (Lin_term_count != 1 || strcmp(tmp_store.name, var) != 0)
        Lin_term_count++;

    /* always store objective function with rownr == 0. */
    if (row == 0)
        return (store(var, row, value));

    if (Lin_term_count == 1) { /* don't store yet. could be a bound */
        strcpy(tmp_store.name, var);
        tmp_store.row = row;
        tmp_store.value += value;
        return (TRUE);
    }

    if (Lin_term_count == 2) { /* now you can also store the first variable */
        rside *rp;

        /* make space for the rhs information */
        if (CALLOC(rp, 1) == NULL)
            return (FALSE);
        rp->next = First_rside;
        First_rside = rs = rp;
        rs->row = row;
        rs->value = tmp_store.rhs_value;
        rs->relat = tmp_store.relat;
        rs->range_relat = -1;

        if (tmp_store.value != 0) {
            if (!store(tmp_store.name, tmp_store.row, tmp_store.value))
                return (FALSE);
        } else {
            if (Verbose >= NORMAL)
                report(NULL, NORMAL,
                    "Warning, variable %s has an effective coefficient of 0 on line %d. Ignored.",
                    tmp_store.name, yylineno);
        }

        null_tmp_store();
    }

    return (store(var, row, value));
} /* var_store */

/*
 * store the information in tmp_store because it is a bound
 */
int store_bounds(void) {
    if (tmp_store.value != 0) {
        hashelem *h_tab_p;
        REAL boundvalue;

        if ((h_tab_p = findhash(tmp_store.name, Hash_tab)) == NULL) {
            /* a new columnname is found, create an entry in the hashlist */
            if ((h_tab_p = puthash(tmp_store.name, Hash_tab)) == NULL)
                return (FALSE);
            Columns++; /* counter for calloc of final array */
            /* create a place to store bounds information */
            if ((h_tab_p->bnd = create_bound_rec()) == NULL)
                return (FALSE);
        } else if (h_tab_p->bnd == NULL)
            /* create a place to store bounds information */
            if ((h_tab_p->bnd = create_bound_rec()) == NULL)
                return (FALSE);

        /* else bound_rec already exists */

        if (tmp_store.value < 0) { /* divide by negative number, */
            /* relational operator may change */
            if (tmp_store.relat == GE)
                tmp_store.relat = LE;
            else if (tmp_store.relat == LE)
                tmp_store.relat = GE;
        }
        /* Check sanity of bound; all variables should be positive */
        boundvalue = tmp_store.rhs_value / tmp_store.value;
        if (((tmp_store.relat == EQ) && (boundvalue < 0))
                || ((tmp_store.relat == LE) && (boundvalue < 0))) { /* Error */
            if (Verbose >= CRITICAL)
                report(NULL, CRITICAL, "Error on line %d: variables must always be non-negative",
                    yylineno);
            return (FALSE);
        }

        if ((tmp_store.relat == GE) && (boundvalue <= 0)) /* Warning */
            if (Verbose >= NORMAL)
                report(NULL, NORMAL,
                    "Warning on line %d: useless bound; variables are always >= 0",
                    yylineno);

        /* bound seems to be sane, add it */
        if ((tmp_store.relat == GE) || (tmp_store.relat == EQ)) {
            if (h_tab_p->bnd->lowbo < boundvalue)
                h_tab_p->bnd->lowbo = boundvalue;
            else
                if (Verbose >= NORMAL)
                report(NULL, NORMAL, "Ineffective lower bound on line %d, ignored",
                    yylineno);
        }
        if ((tmp_store.relat == LE) || (tmp_store.relat == EQ)) {
            if (h_tab_p->bnd->upbo > boundvalue)
                h_tab_p->bnd->upbo = boundvalue;
            else
                if (Verbose >= NORMAL)
                report(NULL, NORMAL, "Ineffective upper bound on line %d, ignored",
                    yylineno);
        }

        /* check for empty range */
        if (h_tab_p->bnd->upbo < h_tab_p->bnd->lowbo) {
            if (Verbose >= CRITICAL)
                report(NULL, CRITICAL, "Error: bound on line %d contradicts earlier bounds, exiting",
                    yylineno);
            return (FALSE);
        }
    } else /* tmp_store.value = 0 ! */ {
        if (Verbose >= CRITICAL)
            report(NULL, CRITICAL, "Error, variable %s has an effective coefficient of 0 in bound on line %d. Exiting.",
                tmp_store.name, yylineno);
        return (FALSE);
    }

    null_tmp_store();
    return (TRUE);
} /* store_bounds */

int add_constraint_name(char *name, int row) {
    constraint_name *cnp;

    if (!First_constraint_name) { /* first time only */
        if (CALLOC(First_constraint_name, 1) == NULL)
            return (FALSE);
        cnp = First_constraint_name;
        rs = NULL;
    } else {
        cnp = First_constraint_name;
        while (cnp) {
            if (strcmp(cnp->name, name) == 0)
                break;
            cnp = cnp->next;
        }
        if (cnp) {
            row = cnp->row;

            rs = First_rside;
            while ((rs != NULL) && (rs->row != row))
                rs = rs->next;
        } else {
            cnp = First_constraint_name;
            if (CALLOC(First_constraint_name, 1) == NULL)
                return (FALSE);
            First_constraint_name->next = cnp;
            cnp = First_constraint_name;
            rs = NULL;
        }
    }
    strcpy(cnp->name, name);
    cnp->row = row;
    return (TRUE);
}

/*
 * transport the data from the intermediate structure to the sparse matrix
 * and free the intermediate structure
 */
int readinput(lprec *lp) {
    int i, j, index, nn_ind;
    column *cp, *tcp; /* tcp (temporary cp) points to memory-space to free */
    hashelem *hp, *thp;
    bound *bp;
    rside *rp;
    constraint_name *cnp;
    REAL *row;

    if (MALLOC(row, 1 + Rows) == NULL)
        return (FALSE);

    /* fill names with the rownames */
    for (cnp = First_constraint_name; cnp; cnp = cnp->next)
        if (!set_row_name(lp, cnp->row, cnp->name)) {
            free(row);
            return (FALSE);
        }

    for (i = Rows; i >= 0; i--) {
        rp = First_rside;
        if ((rp->range_relat >= 0) && (rp->range_value == rp->value)) {
            rp->relat = EQ;
            rp->range_relat = EQ;
        }
        relat[i] = rp->relat;
        lp->orig_rh[i] = rp->value;
        if (rp->range_relat >= 0)
            lp->orig_upbo[i] = my_abs(rp->range_value - rp->value);
        First_rside = rp->next;
        free(rp); /* free memory when data has been read */
    }

    /* change upperbound to zero if the relational operator is the equal sign */
    for (i = 1; i <= Rows; i++)
        if (relat[i] == EQ)
            lp->orig_upbo[i] = 0;

    /*
      for(i = 0; i <= Rows; i++)
        if((!lp->names_used) || (lp->row_name[i] == NULL) || (strcmp(lp->row_name[i], "")==0)) {
          char holdstr[NAMELEN];
          sprintf(holdstr,ROWNAMEMASK,i);
          if (!set_row_name(lp,i,holdstr)) {
            free(row);
            return(FALSE);
          }
        }
     */

    /* start reading the Hash_list structure */
    index = 0;
    nn_ind = 0;

    for (i = 0; i < Hash_tab->size; i++) {
        hp = Hash_tab->table[i];
        while (hp != NULL) {
            /* put an index in the cend array when a new name is found */
            lp->col_end[index++] = nn_ind;

            /* check if it must be an integer variable */
            if (hp->must_be_int) {
                /* lp->must_be_int[Rows + index]=TRUE; */
                set_int(lp, index, TRUE);
            }
            /* check for bound */
            if (hp->bnd != NULL) {
                bp = hp->bnd;
                lp->orig_lowbo[Rows + index] = bp->lowbo;
                lp->orig_upbo[Rows + index] = bp->upbo;
                free(bp); /* free memory when data has been read*/
            }

            /* copy name of column variable */
            /* lp->col_name[index] = hp->name; */
            if (!set_col_name(lp, index, hp->name)) {
                free(row);
                return (FALSE);
            }

            for (j = 0; j <= Rows; j++)
                row[j] = 0.0;

            /* put matrix values in intermediate row */
            cp = hp->col;
            while (cp != NULL) {
                row[cp->row] = cp->value;
                tcp = cp;
                cp = cp->next;
                free(tcp); /* free memory when data has been read */
            }
            thp = hp;
            hp = hp->next;
            free(thp->name);
            free(thp); /* free memory when data has been read */

            /* put matrix values in sparse matrix */
            /* this makes sure that values are in order in the sparse matrix. This is a requirement */
            for (j = 0; j <= Rows; j++)
                if (row[j]) {
                    lp->mat[nn_ind].row_nr = j;
                    lp->mat[nn_ind].value = row[j];
                    nn_ind++;
                }
        }
        Hash_tab->table[i] = NULL;
    }
    lp->col_end[index] = nn_ind;

    /* the following should be replaced by a call to the MPS print routine MB */

#if 0
    if (Verbose) {
        int j;

        printf("\n");
        printf("**********Data read**********\n");
        printf("Rows    : %d\n", Rows);
        printf("Columns : %d\n", Columns);
        printf("Nonnuls : %d\n", Non_zeros);
        printf("NAME          LPPROB\n");
        printf("ROWS\n");
        for (i = 0; i <= Rows; i++) {
            if (relat[i] == LE)
                printf(" L  ");
            else if (relat[i] == EQ)
                printf(" E  ");
            else if (relat[i] == GE)
                printf(" G  ");
            else if (relat[i] == OF)
                printf(" N  ");
            printf("%s\n", get_row_name(lp, i));
        }

        printf("COLUMNS\n");
        j = 0;
        for (i = 0; i < Non_zeros; i++) {
            if (i == lp->col_end[j])
                j++;
            printf("    %-8s  %-8s  %g\n", get_col_name(lp, j),
                    get_row_name(lp, lp->mat[i].row_nr), (double) lp->mat[i].value);
        }

        printf("RHS\n");
        for (i = 0; i <= Rows; i++) {
            printf("    RHS       %-8s  %g\n", get_row_name(lp, i),
                    (double) lp->orig_rh[i]);
        }

        printf("RANGES\n");
        for (i = 1; i <= Rows; i++)
            if ((lp->orig_upbo[i] != lp->infinite) && (lp->orig_upbo[i] != 0)) {
                printf("    RGS       %-8s  %g\n", get_row_name(lp, i),
                        (double) lp->orig_upbo[i]);
            } else if ((lp->orig_lowbo[i] != 0)) {
                printf("    RGS       %-8s  %g\n", get_row_name(lp, i),
                        (double) -lp->orig_lowbo[i]);
            }

        printf("BOUNDS\n");
        for (i = Rows + 1; i <= Rows + Columns; i++) {
            if ((lp->orig_lowbo[i] != 0) && (lp->orig_upbo[i] < lp->infinite) &&
                    (lp->orig_lowbo[i] == lp->orig_upbo[i])) {
                printf(" FX BND       %-8s  %g\n", get_col_name(lp, i - Rows),
                        (double) lp->orig_upbo[i]);
            } else {
                if (lp->orig_upbo[i] < lp->infinite)
                    printf(" UP BND       %-8s  %g\n", get_col_name(lp, i - Rows),
                        (double) lp->orig_upbo[i]);
                if (lp->orig_lowbo[i] > 0)
                    printf(" LO BND       %-8s  %g\n", get_col_name(lp, i - Rows),
                        (double) lp->orig_lowbo[i]);
            }
        }

        printf("ENDATA\n");
    }
#endif

    free(row);
    return (TRUE);
} /* readinput */

lprec *read_LP(char *input, short verbose, char *lp_name) {
    FILE *fpin;
    lprec *lp = NULL;

    if ((fpin = fopen(input, "r")) != NULL) {
        lp = read_lp_file(fpin, verbose, lp_name);
        fclose(fpin);
    }
    return (lp);
}

lprec *read_lp_file(FILE *input, short verbose, char *lp_name) {
    lprec *lp = NULL;
    int i;

    Verbose = verbose;

    yyin = input;
    Maximise = TRUE;
    yyparse();

    Rows--;

    if (CALLOC(relat, Rows + 1) != NULL) {
        lp = make_lpext(Rows, Columns, Non_zeros, Non_zeros, lp_name);
        if (lp != NULL) {
            lp->verbose = verbose;

            if (!readinput(lp)) {
                delete_lp(lp);
                lp = NULL;
            } else {
                if (Maximise)
                    set_maxim(lp);

                for (i = 1; i <= Rows; i++) {
                    REAL a = lp->orig_upbo[i]; /* keep upper bound (range) */
                    set_constr_type(lp, i, relat[i]);
                    lp->orig_upbo[i] = a; /* restore upper bound (range) */
                }

                /* lets free the temporary list of constraint names */
                while (First_constraint_name) {
                    constraint_name *cp;

                    cp = First_constraint_name;
                    First_constraint_name = First_constraint_name->next;
                    free(cp);
                }

                free_hash_table(Hash_tab);
            }
        }
        free(relat);
    }
    return (lp);
}
