/* Last changed Time-stamp: <2005-07-23 16:50:24 ivo> */
/*                
	     Compute duplex structure of two RNA strands

			   c Ivo L Hofacker
			  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "duplex.h"
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"

extern void  read_parameter_file(const char fname[]);
extern int subopt_sorted;
static void  print_struc(duplexT const *dup);
void parse_mapping_array (int **map_array, int *array_size, char *arguments_string);

/*@unused@*/
static char rcsid[] = "$Id: RNAduplex.c,v 1.4 2005/07/24 08:35:15 ivo Exp $";

static char  scale[] = "....,....1....,....2....,....3....,....4"
                       "....,....5....,....6....,....7....,....8";

static void usage(void);

#define PRIVATE static

int  array_size_i = 0;
int  array_size_j = 0;
int  *array_i = NULL;
int  *array_j = NULL;
int five_prime_length = 0;
int auto_five_prime_length = 0;
int E5p, E3p;
int	debug_mode = 0;
int	all_pairs = 0;

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  char *s1, *s2, *line;
  char  fname[13];
  char  tmp_string[255];
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  int   i, l, sym, r;
  double deltaf;
  double kT, sfact=1.07;
  int   pf=0, istty, delta=-1;
  int noconv=0;
  int l1, l2;
   
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-') 
      switch ( argv[i][1] )
	{
	case 'T':  if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &temperature);
	  if (!r) usage();
	  break;
	case 'p':
	  fprintf(stderr, "partition function folding not yet implemented\n");
	  usage();
	  pf=1;
	  if (argv[i][2]!='\0')
	    (void) sscanf(argv[i]+2, "%d", &do_backtrack);
	  break;
	case 'n':
	  if ( strcmp(argv[i], "-noGU")==0) noGU=1;
	  if ( strcmp(argv[i], "-noCloseGU")==0) no_closingGU=1;
	  if ( strcmp(argv[i], "-noLP")==0) noLonelyPairs=1;
	  if ( strcmp(argv[i], "-nsp") ==0) {
	    if (i==argc-1) usage();
	    ns_bases = argv[++i];
	  }
	  if ( strcmp(argv[i], "-noconv")==0) noconv=1;
	  break;
	case '4':
	  tetra_loop=0;
	  break;
	case 'e':
	  if (i>=argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &deltaf);
	  if (r!=1) usage();
	  delta = (int) (0.1+deltaf*100);
	  break;
	case 's': subopt_sorted=1;
	  break;
	case 'i':
	  if (i>=argc-1) usage();
	  r=sscanf(argv[++i], "%s", tmp_string);
	  if (r!=1) usage();
	  
	  parse_mapping_array (&array_i, &array_size_i, tmp_string);
	  break;

	case 'j':
	  if (i>=argc-1) usage();
	  r=sscanf(argv[++i], "%s", tmp_string);
	  if (r!=1) usage();
	  
	  parse_mapping_array (&array_j, &array_size_j, tmp_string);
	  break;

	case '5':
	  if (i>=argc-1) usage();
	  r=sscanf(argv[++i], "%d", &five_prime_length);
	  if (r!=1) usage();
	  if (five_prime_length == 0) auto_five_prime_length = 1;
	  break;

	case 'D':
	  if (i>=argc-1) usage();
	  r=sscanf(argv[++i], "%d", &debug_mode);
	  if (r!=1) usage();
	  break;

	case 'a':
	  if (i>=argc-1) usage();
	  r=sscanf(argv[++i], "%d", &all_pairs);
	  if (r!=1) usage();
	  break;
	
	case 'd': dangles=0;
	  if (argv[i][2]!='\0') {
	    r=sscanf(argv[i]+2, "%d", &dangles);
	    if (r!=1) usage();
	  }
	  break;
	case 'P':
	  if (i==argc-1) usage();
	  ParamFile = argv[++i];
	  break;
	default: usage();
	} 
  }

  if (ParamFile != NULL)
    read_parameter_file(ParamFile);
   
  if (ns_bases != NULL) {
    nonstandards = space(33);
    c=ns_bases;
    i=sym=0;
    if (*c=='-') {
      sym=1; c++;
    }
    while (*c!='\0') {
      if (*c!=',') {
	nonstandards[i++]=*c++;
	nonstandards[i++]=*c;
	if ((sym)&&(*c!=*(c-1))) {
	  nonstandards[i++]=*c;
	  nonstandards[i++]=*(c-1);
	}
      }
      c++;
    }
  }
  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
	
  do {				/* main loop: continue until end of file */
    duplexT mfe, *subopt;
    if (istty) {
      printf("\nInput two sequences (one line each); @ to quit\n");
      printf("%s\n", scale);
    }
    fname[0]='\0';
    
    if ((line = get_line(stdin))==NULL) break;
    /* skip empty lines, comment lines, name lines */
    while (line && ((*line=='*')||(*line=='\0')||(*line=='>'))) {
      printf("%s\n", line); free(line);
      if ((line = get_line(stdin))==NULL) break;
    } 
    if ((line ==NULL) || (strcmp(line, "@") == 0)) break;
    
    s1 = (char *) space(strlen(line)+1);
    (void) sscanf(line,"%s",s1);  free(line);
    
    if ((line = get_line(stdin))==NULL) break;
    /* skip comment lines and get filenames */
    while ((*line=='*')||(*line=='\0')||(*line=='>')) {
      printf("%s\n", line); free(line);
      if ((line = get_line(stdin))==NULL) break;
    } 
    if ((line ==NULL) || (strcmp(line, "@") == 0)) break;
    
    s2 = (char *) space(strlen(line)+1);
    (void) sscanf(line,"%s",s2); free(line);


	/* get force_binding strings */
	if (!all_pairs)
	{	
		if ((line = get_line(stdin))==NULL) break;
		(void) sscanf(line,"%s",tmp_string); free(line);
		parse_mapping_array (&array_i, &array_size_i, tmp_string);
	
		if ((line = get_line(stdin))==NULL) break;
		(void) sscanf(line,"%s",tmp_string); free(line);
		parse_mapping_array (&array_j, &array_size_j, tmp_string);
	}

    for (l = 0; l < strlen(s1); l++) {
      s1[l] = toupper(s1[l]);
      if (!noconv && s1[l] == 'T') s1[l] = 'U';
    }
    for (l = 0; l < strlen(s2); l++) {
      s2[l] = toupper(s2[l]);
      if (!noconv && s2[l] == 'T') s2[l] = 'U';
    }
    if (istty)
      printf("lengths = %d,%d\n", strlen(s1), strlen(s2));

    /* initialize_fold(length); */
    update_fold_params();
    if (delta>=0) {
      duplexT *sub;
      subopt = duplex_subopt(s1, s2, delta, 0);
      for (sub=subopt; sub->i >0; sub++) {
	print_struc(sub);
	free(sub->structure);
      }
      free(subopt);
    }
    else {
      mfe = duplexfold(s1, s2);
      
      /* print_struc(&mfe); */
      
      l1 = strchr(mfe.structure, '&') - mfe.structure;
      l2 = strlen(mfe.structure) - l1 - 1;
      
	  printf ("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%5.2f\t%5.2f\t%5.2f\n", 
  			mfe.structure,
  			mfe.i+1-l1,
	 		mfe.i,
	 		mfe.j,
	 		mfe.j+l2-1,
	 		l1,
	 		l2,
	 		mfe.energy, (float) (E5p) / 100., (float) (E3p) / 100.);
      
      
      free(mfe.structure);
    }
    (void) fflush(stdout);
       
    (void) fflush(stdout);
    free(s1); free(s2);
  } while (1);
  return 0;
}

static void print_struc(duplexT const *dup) {
  int l1;
  l1 = strchr(dup->structure, '&')-dup->structure;
  

  printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", dup->structure, dup->i+1-l1,
	 dup->i, dup->j, dup->j+strlen(dup->structure)-l1-2, dup->energy);
}
    
static void usage(void)
{
  nrerror("usage:\n"
	  "RNAduplex [-e range] [-s]\n"
	  "          [-T temp] [-4] [-d] [-noGU] [-noCloseGU]\n" 
	  "          [-noLP] [-P paramfile] [-nsp pairs] [-noconv]\n");
}

void parse_mapping_array (int **map_array, int *array_size, char *arguments_string)
{
	char *p;
	int  value;
	int  cur_array_size = 0;
	
	if (*map_array) free (*map_array);
	
	*map_array = (int *) space(sizeof(int) * (50));

	p = strtok (arguments_string,",");
	while (p != NULL)
	{
		value = atoi (p);
		/* printf ("position %i: string=%s, value=%i\n", cur_array_size, p, value); */
		(*map_array)[cur_array_size++] = value;
		p = strtok (NULL, ",");
	}
	
	*array_size = cur_array_size;
}
