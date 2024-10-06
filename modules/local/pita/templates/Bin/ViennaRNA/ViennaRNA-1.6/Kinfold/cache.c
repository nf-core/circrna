/*
  Last changed Time-stamp: <2001-08-02 14:59:19 xtof>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: cache.c,v 1.2 2005/02/16 17:00:48 ivo Exp $
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "cache_util.h"

/*
  modify cache_f(), cache_comp() and the typedef of cache_entry
  in cache_utils.h to suit your application
*/

/* PUBLIC FUNCTIONES */
cache_entry *lookup_cache (char *x);
int write_cache (cache_entry *x);
/*  void delete_cache (cache_entry *x); */
void kill_cache();
void initialize_cache();

/* PRIVATE FUNCTIONES */
/*  static int cache_comp(cache_entry *x, cache_entry *y); */
static unsigned cache_f (char *x);

/* #define CACHESIZE 67108864 -1 */ /* 2^26 -1   must be power of 2 -1 */
/* #define CACHESIZE 33554432 -1 */ /* 2^25 -1   must be power of 2 -1 */
/* #define CACHESIZE 16777216 -1 */ /* 2^24 -1   must be power of 2 -1 */ 
/* #define CACHESIZE  4194304 -1 */ /* 2^22 -1   must be power of 2 -1 */
#define CACHESIZE  1048576 -1  /* 2^20 -1   must be power of 2 -1 */
/* #define CACHESIZE   262144 -1 */ /* 2^18 -1   must be power of 2 -1 */
/* next is default */
/* #define CACHESIZE    65536 -1 */ /* 2^16 -1   must be power of 2 -1 */
/* #define CACHESIZE    16384 -1 */ /* 2^14 -1   must be power of 2 -1 */
/* #define CACHESIZE     4096 -1 */ /* 2^12 -1   must be power of 2 -1 */

static cache_entry *cachetab[CACHESIZE+1];
static char UNUSED rcsid[] ="$Id: cache.c,v 1.2 2005/02/16 17:00:48 ivo Exp $";
unsigned long collisions=0;

/* stolen from perl source */
char coeff[] = {
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1};

/* key must not be longer than 128 */
#pragma inline (cache_f)
unsigned cache_f(char *x) { 
  register char *s;
  register int i;
  register unsigned cache;

  s = x;

  for (i=0,    cache = 0;
       /* while */ *s;
       s++,           i++ , cache *= 5 ) {
    cache += *s * coeff[i];
  }

   /* divide through CACHESIZE for normalization */
  return ((cache) & (CACHESIZE));
}


/* returns NULL unless x is in the cache */
cache_entry *lookup_cache (char *x) {
  int cacheval;
  cache_entry *c;

  cacheval=cache_f(x);
  if ((c=cachetab[cacheval]))
    if (strcmp(c->structure,x)==0) return c;
  
  return NULL;
}

/* returns 1 if x already was in the cache */
int write_cache (cache_entry *x) {
  int cacheval;
  cache_entry *c;
  
  cacheval=cache_f(x->structure);
  if ((c=cachetab[cacheval])) {
    free(c->structure);
    free(c->neighbors);
    free(c->rates);
    free(c);
  }
  cachetab[cacheval]=x;
  return 0;
}

/**/
void initialize_cache () { }

/**/
void kill_cache () {
  int i;
  
  for (i=0;i<CACHESIZE+1;i++) {
    if ( cachetab[i] ) {
      free (cachetab[i]->structure);
      free (cachetab[i]->neighbors);
      free (cachetab[i]->rates);
      free (cachetab[i]);
    }
    cachetab[i]=NULL;
  }
}

#if 0
/**/
static int cache_comp(cache_entry *x, cache_entry *y) {
  return strcmp(((cache_entry *)x)->structure, ((cache_entry *)y)->structure);
}
#endif

/* End of file */
