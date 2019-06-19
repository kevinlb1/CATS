#include <stdio.h>

#ifndef __LPGLOB_H__
#define __LPGLOB_H__

/* Globals for parser */

extern int Rows;
extern int Columns;
extern int Non_zeros;

extern FILE *yyin;
extern FILE *lpfilename;
extern short Maximise;
extern short *relat;
extern int yylineno;
extern int yyleng;
extern int Lin_term_count;
extern int Sign;
extern constraint_name *First_constraint_name;
/* I hate #ifdefs, but there seems to be no "standard" way to do this */
#if defined(__hpux) || defined(__apollo) || defined(_AIX) || defined(_OSF_SOURCE)
/* for HP/UX, Apollo, AIX, DEC OSF  */
extern unsigned char yytext[];
#else
/* For other computers */
extern char yytext[];
#endif

extern rside *First_rside;
extern short Ignore_decl;

extern tmp_store_struct tmp_store;

#endif
