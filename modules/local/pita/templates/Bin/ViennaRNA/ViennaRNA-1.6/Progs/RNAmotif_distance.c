/************************************************************************/
/*   		  Distances of Secondary Structures                     */
/************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "dist_vars.h"
#include "RNAstruct.h"
#include "treedist.h"
#include "stringdist.h"
#include "utils.h"


#define MAXNUM      1000    /* max number of structs for distance matrix */
#define PRIVATE     static
#define PUBLIC



PRIVATE void usage(void);
PRIVATE int  parse_input(char *line);



PRIVATE int  types = 1;         /* number of distance types */
PRIVATE char ttype[10] = "f";   /* distance types */
PRIVATE int  n = 0;             /* number of structures read so far */


int main(int argc, char *argv[])
{
  int       tree_types = 0, ttree;
  int       string_types = 0, tstr;
  char     *line=NULL, *xstruc, *cc;
  int       i, j, tt, type;
  int       it, is; 
  float     dist;
  int       comp_to_first = 0;
  int       max = MAXNUM;

  
  edit_backtrack = 0;
  
  /* reading parameters */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (argv[i][1] == 'D') {
	strncpy(ttype, argv[i]+2, 9);
	types = (int)strlen(ttype);
      }
      else if (argv[i][1] == 'F')
	comp_to_first = 1;
      else if (argv[i][1] == 'M') {
	max = atoi(argv[i+1]) + 10;
	i++;
      }
      else
	usage();
    }
    else
      usage();
  }

  Tree     *T[5][max];
  swString *S[5][max];
  char     *P[max];   /* structures for base pair distances */
  
  /* reading sequences from stdin and calculating distances */
  line = get_line(stdin);
  type = parse_input(line);
  while (type != 999) {

    while (type == 0) {
      if (line != NULL) 
	free(line);

      line = get_line(stdin);
      type = parse_input(line);
    }

    tree_types = 0;
    string_types = 0;
    for(tt = 0; tt < types; tt++) {

      switch(ttype[tt]){
      case 'f' :
      case 'F' :
	xstruc = expand_Full(line);
	if(islower(ttype[tt]))  /* tree_edit */
	  T[tree_types++][n] = make_tree(xstruc);
	if(isupper(ttype[tt])) /* string edit */
	  S[string_types++][n] = Make_swString(xstruc);
	free(xstruc);
	break;

      case 'P':
	P[n] = strdup(line);
	break;

      case 'h' :
      case 'H' :
	xstruc = b2HIT(line);
	if(islower(ttype[tt])) 
	  T[tree_types++][n] = make_tree(xstruc);
	if(isupper(ttype[tt])) 
	  S[string_types++][n] = Make_swString(xstruc);
	free(xstruc);
	break;

      case 'c' :
      case 'C' :    
	cc = b2C(line);
	xstruc = expand_Shapiro(cc);
	free(cc);
	if(islower(ttype[tt]))
	  T[tree_types++][n] = make_tree(xstruc);
	if(isupper(ttype[tt]))
	  S[string_types++][n] = Make_swString(xstruc);
	free(xstruc);   
	break;

      case 'w' :
      case 'W' :          
	xstruc = b2Shapiro(line);
	if(islower(ttype[tt])) 
	  T[tree_types++][n] = make_tree(xstruc);
	if(isupper(ttype[tt]))
	  S[string_types++][n] = Make_swString(xstruc);
	free(xstruc); 
	break;

      default: 
	nrerror("Unknown distance type");
      }
    }
    n++;

    if (n >= 2) {
      it = 0;
      is = 0;
      if (comp_to_first) {
	for (i=0; i < types; i++) {
	  if(islower(ttype[i])) {
	    dist = tree_edit_distance(T[it][0], T[it][1]);
	    free_tree(T[it][1]);
	    it++;
	  }
	  else if (ttype[i]=='P') {
	    dist = (float)bp_distance(P[0], P[1]);
	    free(P[1]);
	  }
	  else /* isupper(ttype[i]) */ {
	    dist = string_edit_distance(S[is][0], S[is][1]);
	    free(S[is][1]);
	    is++;
	  }
	  printf("%c: %g  \n", ttype[i], dist);
	}
	n = 1;
      }
      else {
	for (i=0; i < types; i++) {
	  if(islower(ttype[i])) {
	    dist = tree_edit_distance(T[it][0], T[it][1]);
	    free_tree(T[it][0]); 
	    free_tree(T[it][1]);
	    it++;
	  }
	  else if (ttype[i]=='P') {
	    dist = (float)bp_distance(P[0], P[1]);
	    free(P[0]); 
	    free(P[1]);
	  }
	  else /* isupper(ttype[i]) */ {
	    dist = string_edit_distance(S[is][0], S[is][1]);
	    free(S[is][0]); 
	    free(S[is][1]);
	    is++;
	  }
	  printf("%c: %g  \n", ttype[i], dist);
	}
	n = 0;
      }
    }
    fflush(stdout);

    free(line);
    line = get_line(stdin);
    type = parse_input(line);
  }

  if (comp_to_first) {
    it = 0;
    is = 0;
    for (i=0; i < types; i++) {
      if(islower(ttype[i]))
	free_tree(T[it++][0]);
      else if (ttype[i]=='P')
	free(P[0]);
      else /* isupper(ttype[i]) */
	free(S[is++][1]);
    }
  }
  
  return 0;
}

/*--------------------------------------------------------------------------*/

PRIVATE int parse_input(char *line)
{
  int i,o;

  if (line == NULL) 
    return 999;

  i=o=0;
  while( line[i] ){
    switch(line[i]) {
    case '(' :
      o++;
      i++;
      break;
    case '.' :
      i++;
      break;
    case ')' : 
      i++;
      o--;
      if(o<0) return 0;
      break;
    default:
      return 0;
    }
  }
  if (o>0) 
    return 0;
  return 1;
}

/*--------------------------------------------------------------------------*/

PRIVATE void usage(void)
{
  nrerror("usage: RNAmotif_distance [-D[fhwcFHWCP]] [-F] [-M <num>]");
}
