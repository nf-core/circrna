/*
  Last changed Time-stamp: <2006-01-16 11:28:43 ivo>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: nachbar.c,v 1.6 2006/01/16 10:38:11 ivo Exp $
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include "globals.h"
#include "assert.h"
#include "fold_vars.h"
#include "energy_const.h"
#include "utils.h"
#include "cache_util.h"
#include "baum.h"

static char UNUSED rcsid[]="$Id: nachbar.c,v 1.6 2006/01/16 10:38:11 ivo Exp $";

/* arrays */
static short *neighbor_list=NULL;
static float *bmf; /* boltzmann weight of structure */

/* variables */
/*  static double highestE = -1000.0; */
/*  static double OhighestE = -1000.0; */
/*  static char *highestS, *OhighestS; */
static int lmin = 1;
static int top = 0;
static int is_from_cache = 0;
/*  static double meanE = 0.0; */
static double totalflux = 0.0;
static double Zeit = 0.0;
static double zeitInc = 0.0;
static double _RT = 0.6;

/* public functiones */
void ini_nbList(int chords);
void update_nbList(int i, int j, int iE);
int sel_nb(void);
void clean_up_nbList(void);
extern void update_tree(int i, int j);

/* privat functiones */
static void reset_nbList(void);
static void grow_chain(void);
static FILE *logFP=NULL;

/**/
void ini_nbList(int chords) {
  char logFN[256];

  _RT = (((temperature + K0) * GASCONST) / 1000.0);
  if (neighbor_list!=NULL) return;
  /*
    list for move coding
    make room for 2*chords neighbors (safe bet)
  */
  if (chords == 0) chords = 1;
  neighbor_list = (short *)calloc(4*chords, sizeof(short));
  assert(neighbor_list != NULL);
  /*
    list for Boltzmann-factors
  */
  bmf = (float *)calloc(2*chords, sizeof(double));
  assert(bmf != NULL);

  /* open log-file */
  logFP = fopen(strcat(strcpy(logFN, GAV.BaseName), ".log"), "a+");
  assert(logFP != NULL);

  /* log initial condition */
  log_prog_params(logFP);
  log_start_stop(logFP);
}

/**/
void update_nbList(int i, int j, int iE) {
  double E, dE, p;

  E = (double)iE/100.;
  neighbor_list[2*top] = (short )i;
  neighbor_list[2*top+1] = (short )j;

  /* compute rates and some statistics */
  /*    meanE += E; */
  dE = E-GSV.currE;
  if( GTV.mc ) {
    /* metropolis rule */
    if (dE < 0) p = 1;
    else p = exp(-(dE / _RT));
  }
  else  /* kawasaki rule */
    p = exp(-0.5 * (dE / _RT));

  totalflux += p;
  bmf[top++] = (float )p;
  if (dE < 0) lmin = 0;
  if ((dE == 0) && (lmin==1)) lmin = 2;
}

/**/
void get_from_cache(cache_entry *c) {
  top = c->top;
  totalflux = c->flux;
  GSV.currE = c->energy;
  lmin = c->lmin;
  memcpy(neighbor_list, c->neighbors, 2*top*sizeof(short));
  memcpy(bmf, c->rates, top*sizeof(float));
  is_from_cache = 1;
}

/**/
void put_in_cache(void) {
  cache_entry *c;

  if ((c = (cache_entry *) malloc(sizeof(cache_entry)))==NULL) {
    fprintf(stderr, "out of memory\n"); exit(255);
  }
/*    c->structure = strdup(GAV.currform); */
  c->structure = (char *) calloc(GSV.len+1, sizeof(char));
  strcpy(c->structure, GAV.currform);
  c->neighbors = (short *) malloc(top*2*sizeof(short));
  memcpy(c->neighbors,neighbor_list,top*2*sizeof(short));
  c->rates = (float *) malloc(top*sizeof(float));
  memcpy(c->rates, bmf, top*sizeof(float));
  c->top = top;
  c->lmin = lmin;
  c->flux = totalflux;
  c->energy = GSV.currE;
  write_cache(c);
}

/*============*/

int sel_nb(void) {

  char trans, **s;
  int next;
  double pegel = 0.0, schwelle = 0.0, zufall = 0.0;
  int found_stop=0;

  /* before we select a move, store current conformation in cache */
  /* ... unless it just came from there */
  if ( !is_from_cache ) put_in_cache();
  is_from_cache = 0;

  /* draw 2 different a random number */
  schwelle = urn();
  while ( zufall==0 ) zufall = urn();

  /* advance internal clock */
  if (totalflux>0)
    zeitInc = (log(1. / zufall) / totalflux);
  else {
    if (GSV.grow>0) zeitInc=GSV.grow;
    else zeitInc = GSV.time;
  }

  Zeit += zeitInc;
  if (GSV.grow>0 && GSV.len < strlen(GAV.farbe_full)) grow_chain();

  /* meanE /= (double)top; */

  /* normalize boltzmann weights */
  schwelle *=totalflux;

  /* and choose a neighbour structure next */
  for (next = 0; next < top; next++) {
    pegel += bmf[next];
    if (pegel > schwelle) break;
  }

  /* in case of rounding errors */
  if (next==top) next=top-1;

  /*
    process termination contitiones
  */
  /* is current structure identical to a stop structure ?*/
  for (found_stop = 0, s = GAV.stopform; *s; s++) {
    if (strcmp(*s, GAV.currform) == 0) {
      found_stop = (s - GAV.stopform) + 1;
      break;
    }
  }

  if ( ((found_stop > 0) && (GTV.fpt == 1)) || (Zeit > GSV.time) ) {
    /* met condition to stop simulation */

    /* this gose to stdout */
    if ( !GTV.silent ) {
      printf("%s %6.2f %10.3f", GAV.currform, GSV.currE, Zeit);
      if ( GTV.verbose ) printf("%4d _ %d", top, lmin);
      if (found_stop) printf(" X%d\n", found_stop);/* found a stop structure */
      else printf(" O\n"); /* time for simulation is exceeded */
      fflush(stdout);
    }

    /* this goes to log */
    /* comment log steps of simulation as well !!! %6.2f  round */
    if ( found_stop ) fprintf(logFP," X%02d %12.3f\n", found_stop, Zeit);
    else fprintf(logFP," O   %12.3f %d %s\n", Zeit, lmin, GAV.currform);
    GAV.subi[0] = xsubi[0];
    GAV.subi[1] = xsubi[1];
    GAV.subi[2] = xsubi[2];
    fprintf(logFP, "(%5hu %5hu %5hu)", GAV.subi[0], GAV.subi[1], GAV.subi[2]);
    fflush(logFP);

    Zeit = 0.0;
    /*  highestE = OhighestE = -1000.0; */
    reset_nbList();
    return(1);
  }
  else {
    /* continue simulation */
    int flag = 0;
    if( (!GTV.silent) && (GSV.currE <= GSV.stopE+GSV.cut) ) {
      if (!GTV.lmin || (lmin==1 && strcmp(GAV.prevform, GAV.currform) != 0)) {
	char format[64];
	flag = 1;
	sprintf(format, "%%-%ds %%6.2f %%10.3f", strlen(GAV.farbe_full));
	printf(format, GAV.currform, GSV.currE, Zeit);
      }
      if ( flag && GTV.verbose ) {
	int ii, jj;
	if (next<0) trans='g'; /* growth */
	else {
	  ii = neighbor_list[2*next];
	  jj = neighbor_list[2*next+1];
	  if (abs(ii) < GSV.len) {
	    if ((ii > 0) && (jj > 0)) trans = 'i';
	    else if ((ii < 0) && (jj < 0)) trans = 'd';
	    else if ((ii > 0) && (jj < 0)) trans = 's';
	    else trans = 'S';
	  } else {
	    if ((ii > 0) && (jj > 0)) trans = 'I';
	    else trans = 'D';
	  }
	}
	printf("%4d %c %d", top, trans, lmin);
      }
      if (flag) printf("\n");
    }
  }


  /* store last lmin seen, so we can avoid printing the same lmin twice */
  if (lmin==1)
    strcpy(GAV.prevform, GAV.currform);

#if 0
  if (lmin==1) {
    /* went back to previous lmin */
    if (strcmp(GAV.prevform, GAV.currform) == 0) {
      if (OhighestE < highestE) {
	highestE = OhighestE;  /* delete loop */
	strcpy(highestS, OhighestS);
      }
    } else {
      strcpy(GAV.prevform, GAV.currform);
      OhighestE = 10000.;
    }
  }

  if ( strcmp(GAV.currform, GAV.startform)==0 ) {
    OhighestE = highestE = -1000.;
    highestS[0] = 0;
  }

  /* log highes energy */
  if (GSV.currE > highestE) {
    OhighestE = highestE;
    highestE = GSV.currE;
    strcpy(OhighestS, highestS);
    strcpy(highestS, GAV.currform);
  }
#endif

  if (next>=0) update_tree(neighbor_list[2*next], neighbor_list[2*next+1]);
  else {
    clean_up_rl(); ini_or_reset_rl();
  }

  reset_nbList();
  return(0);
}

/*==========================*/
static void reset_nbList(void) {

  top = 0;
  totalflux = 0.0;
  /*    meanE = 0.0; */
  lmin = 1;
}

/*======================*/
void clean_up_nbList(void){

  free(neighbor_list);
  free(bmf);
  fprintf(logFP,"\n");
  fclose(logFP);
}

/*======================*/
static void grow_chain(void){
  int newl;
  /* note Zeit=0 corresponds to chain length GSV.glen */
  if (Zeit<(GSV.len+1-GSV.glen) * GSV.grow) return;
  newl = GSV.len+1;
  Zeit = (newl-GSV.glen) * GSV.grow;
  top=0; /* prevent structure move in sel_nb */

  if (GSV.len<newl) {
    strncpy(GAV.farbe, GAV.farbe_full, newl);
    GAV.farbe[newl] = '\0';
    strcpy(GAV.startform, GAV.currform);
    strcat(GAV.startform, ".");

    GSV.len = newl;
  }
}
