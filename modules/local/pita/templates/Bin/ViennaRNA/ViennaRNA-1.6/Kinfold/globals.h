/*
  Last changed Time-stamp: <2005-02-16 13:16:43 ivo>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: globals.h,v 1.2 2005/02/16 17:00:48 ivo Exp $
*/

#ifndef GLOBDEFS_H
#define GLOBDEFS_H

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

typedef struct _GlobVars {
  int len;
  int num;
  int maxS;
  int steps;
  float cut;
  float Temp;
  float startE;
  float stopE;
  float currE;
  double grow;
  int    glen;
  double time;
  double simTime;
} GlobVars;

typedef struct _GlobArrays {
  char *ParamFile;
  char *ProgramName;
  char *BaseName;      /* output file basename */
  char *farbe;         /* sequence */
  char *farbe_full;    /* full sequence (for chain growth simulation) */
  char *startform;     /* start structure */
  char **stopform;     /* stop structure(s) */
  char *currform;      /* current structure */
  char *prevform;      /* current structure of previous time step */
  float *sE;           /* energy(s) of stop structure(s) */
  unsigned short subi[3]; /* seeds for random-number-generator */
} GlobArrays;

typedef struct _GlobToogles {
  int Par;
  int seed;
  int dangle;
  int logML;
  int noLP;
  int noShift;
  int start;
  int stop;
  int silent;
  int lmin;
  int fpt;
  int mc;
  int verbose;
} GlobToggles;

void decode_switches(int argc, char *argv[]);
void clean_up_globals(void);
void log_prog_params(FILE *FP);
void log_start_stop(FILE *FP);

GlobVars GSV;
GlobArrays GAV;
GlobToggles GTV;

#endif


/* End of file */
