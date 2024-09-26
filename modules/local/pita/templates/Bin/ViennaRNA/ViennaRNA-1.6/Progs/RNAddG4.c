/* Last changed Time-stamp: <2005-09-22 10:10:55 ivo> */
/*
 * Ineractive Access to folding Routines
 * 
 * c Ivo L Hofacker Vienna RNA package
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "PS_dot.h"
#include "utils.h"
#include "duplex.h"

static void  print_struc(duplexT const *dup);
double LogOfSumOfExps (double* exps, int num_exps);

/* @unused@ */
static char     rcsid[] = "$Id: RNAfold.c,v 1.18 2005/10/30 20:07:12 ivo Exp $";

#define PRIVATE static
#define MIN_NUMBER_TO_EXPONENTIATE -100

static char     scale1[] = "....,....1....,....2....,....3....,....4";
static char     scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE void    usage(void);


int  array_size_i = 0;
int  array_size_j = 0;
int  *array_i;
int  *array_j;
int five_prime_length = 0;
int E5p, E3p;
int	debug_mode = 0;

/*--------------------------------------------------------------------------*/

int 
main(int argc, char *argv[])
{
	char           *string, *line;
	char           *structure = NULL, *constraint = NULL, *sub_target = NULL;
	char           *ParamFile = NULL;
	char           *ns_bases = NULL, *c;
	int             i, length, l, sym, r;
	double          energy, min_en;
	double          kT, sfact = 1.07;
	int             noconv = 0;
	
	duplexT 		hybrid;
	
	double			dG0, dG1, dG2, P;
	double			ddG[50];
	double			dG1array[50];
	double			dG2array[50];
	
	int				n_exps = 0;

	
	int				restrict_from	= 0;	/* Where to start and end restriction of               */
	int				restrict_to		= 0;	/* mRNA folding so that miRNA can bind                 */
	int				target_start	= 0;	/* Base where target starts                            */
	int				upstream_rest	= 0;	/* Number of bases to be restricted upstream of target */

	string = NULL;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'T':
				if (argv[i][2] != '\0')
					usage();
				if (i == argc - 1)
					usage();
				r = sscanf(argv[++i], "%lf", &temperature);
				if (!r)
					usage();
				break;
			case 'n':
				if (strcmp(argv[i], "-noGU") == 0)
					noGU = 1;
				if (strcmp(argv[i], "-noCloseGU") == 0)
					no_closingGU = 1;
				if (strcmp(argv[i], "-noLP") == 0)
					noLonelyPairs = 1;
				if (strcmp(argv[i], "-nsp") == 0) {
					if (i == argc - 1)
						usage();
					ns_bases = argv[++i];
				}
				if (strcmp(argv[i], "-noconv") == 0)
					noconv = 1;
				break;
			case '4':
				tetra_loop = 0;
				break;
			case 'e':
				if (i == argc - 1)
					usage();
				r = sscanf(argv[++i], "%d", &energy_set);
				if (!r)
					usage();
				break;

			case 'f':
				if (i == argc - 1)
					usage();
				r = sscanf(argv[++i], "%i", &restrict_from);
				if (!r)
					usage();
				break;

			case 't':
				if (i == argc - 1)
					usage();
				r = sscanf(argv[++i], "%i", &restrict_to);
				if (!r)
					usage();
				break;

			case 's':
				if (i == argc - 1)
					usage();
				r = sscanf(argv[++i], "%i", &target_start);
				if (!r)
					usage();
				target_start--; /* zero-based */
				break;

			case 'u':
				if (i == argc - 1)
					usage();
				r = sscanf(argv[++i], "%i", &upstream_rest);
				if (!r)
					usage();
				break;
				
			case 'S':
				if (i == argc - 1)
					usage();
				r = sscanf(argv[++i], "%lf", &sfact);
				if (!r)
					usage();
				break;
			case 'd':
				dangles = 0;
				if (argv[i][2] != '\0') {
					r = sscanf(argv[i] + 2, "%d", &dangles);
					if (r != 1)
						usage();
				}
				break;
			case 'P':
				if (i == argc - 1)
					usage();
				ParamFile = argv[++i];
				break;
			default:
				usage();
			}
	}

	if (ParamFile != NULL)
		read_parameter_file(ParamFile);

	if (ns_bases != NULL) {
		nonstandards = space(33);
		c = ns_bases;
		i = sym = 0;
		if (*c == '-') {
			sym = 1;
			c++;
		}
		while (*c != '\0') {
			if (*c != ',') {
				nonstandards[i++] = *c++;
				nonstandards[i++] = *c;
				if ((sym) && (*c != *(c - 1))) {
					nonstandards[i++] = *c;
					nonstandards[i++] = *(c - 1);
				}
			}
			c++;
		}
	}

	/* Some global parameters for the folding routines */
	do_backtrack = 0;
	dangles = 2;
	fold_constrained = 0;
	
	/* main loop: continue until end of file */
	
	do {
	
		n_exps = 0;
		
		/* Read and prepare sequence ***************************/
		
		if ((line = get_line(stdin)) == NULL)
			break;

		string = (char *) space(strlen(line) + 1);
		(void) sscanf(line, "%s", string);
		free(line);
		length = (int) strlen(string);

		structure 	= (char *) space((unsigned) length + 1);
		constraint 	= (char *) space((unsigned) length + 1);
		
		sub_target 	= (char *) space(restrict_to + 1);

		
		for (l = 0; l < length; l++) {
			string[l] = toupper(string[l]);
			if (!noconv && string[l] == 'T')
				string[l] = 'U';
		}

		/* Fold with no constraint ********************************************/

		min_en = fold(string, structure);		
		printf ("%s\t", structure);

		fold_constrained = 0;
		
		kT = (temperature + 273.15) * 1.98717 / 1000.;	/* in Kcal */
		pf_scale = exp(-(sfact * min_en) / kT / length);
		init_pf_fold(length);

		dG0 = pf_fold(string, structure);
		
		free_pf_arrays();

		/* Create constraint sequence from length restrict_from to restrict_to */
		/* and calculate the fold for each constrain						   */
		
		int i,j;
		
		for (i = restrict_from; i <= restrict_to; i++)
		{
			/* Build constraint string */
			
			memset (constraint, ' ', length);
			
			for (j = (target_start + upstream_rest); j > (target_start - i); j--)
			{
				constraint[j] = 'x';
			}
			
			/* for (j = 0 - upstream_rest; j < i; j++)
			{
				constraint[j+target_start] = 'x';
			}*/
			
			if (i < restrict_to)
				constraint[target_start-i] = '|';
		
			fold_constrained = 1;
			strncpy(structure, constraint, length + 1);
	
			/* printf ("%d (from=%d; to=%d; upstream=%d): Constraint is: **%s**\n", i, restrict_from, restrict_to, upstream_rest, structure); */
			
			/* Fold with constraint *******************************************/
			
			min_en = fold(string, structure);
			printf ("%s", structure);
			if (i < restrict_to)
				printf (";");
			
			pf_scale = exp(-(sfact * min_en) / kT / length);		
			init_pf_fold(length);
	
			strncpy(structure, constraint, length + 1);
			dG1 = pf_fold(string, structure);
			
			/* Call RNAduplex to calculate the hybridization energy */
			
			/* strncpy (sub_target, string+target_start, i);
			sub_target[i]='\x0'; */
			
			
			update_fold_params();

			/* For future use: calcualte pf-like by summing many suboptimal structs  
			if (delta>=0) {
				duplexT *sub;
				subopt = duplex_subopt(s1, s2, delta, 0);
				for (sub=subopt; sub->i >0; sub++) {
					print_struc(sub);
					free(sub->structure);
				}
				free(subopt);
			}
			*/
 
			/* hybrid = duplexfold(sub_target, miRNA);
      		/* print_struc(&hybrid);
      		free(hybrid.structure); */
      		
      		printf ("\tno_hybrid_struct");
      		
      		dG2 = 0; /*hybrid.energy; */
      		
      		dG1array[n_exps] = dG1;
      		dG2array[n_exps] = dG2;
      		
      		ddG[n_exps++]    = - (dG1 + dG2 - dG0);
      		/* printf ("dG0=%g, dG1=%g, dG2=%g, ddG=%g\t", dG0, dG1, dG2, -ddG[n_exps-1]); */
		}
		
		double ddGsum = -LogOfSumOfExps (ddG, n_exps);
		
		P = 1 / ( 1 + exp(ddGsum / kT) );	/* Probability of being bound */
		
		printf ("\t%g\t", dG0);
		
		for (i=0; i<n_exps-1; i++)
			printf ("%g;", dG1array[i]);
			
		printf ("%g\t", dG1array[n_exps-1]);
		
		for (i=0; i<n_exps-1; i++)
			printf ("%g;", dG2array[i]);
		
		printf ("%g\t%g\t%g\n", dG2array[n_exps-1], ddGsum, P);

		/* Free everything */

		free_pf_arrays();
		
		if (length >= 2000) free(base_pair);
		(void) fflush(stdout);
		free(string);
		free(structure);
		free(constraint);
	} while (1);
	return 0;
}

static void print_struc(duplexT const *dup) {
  int l1;
  l1 = strchr(dup->structure, '&')-dup->structure;
  printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", dup->structure, dup->i+1-l1,
	 dup->i, dup->j, dup->j+strlen(dup->structure)-l1-2, dup->energy);
}


double LogOfSumOfExps (double* exps, int num_exps)
{

   // compute log(exp(a_1)+exp(a_2)+...exp(a_n)) using:
   // max(a_1,a_2,..a_n) + log(1+exp(a_2-max(a_1,a_2,..a_n))+...exp(a_n-max(a_1,a_2,..a_n)))

	double max_index = 0;
	double max = exps[0];
	int i;
	
	for (i = 1; i < num_exps; ++i)
	{
		if (exps[i] > max)
		{
			max = exps[i];
			max_index = i;
		} 
	}
	
	double remaining_summand = 1.0;
	for (i = 0; i < num_exps; ++i)
	{
		if (i != max_index)
		{
			double difference = exps[i] - max;
			if (difference > MIN_NUMBER_TO_EXPONENTIATE)
				remaining_summand += exp(difference);
		} 
	}
	
	double result = max + log(remaining_summand);
	
	return result;
}


PRIVATE void 
usage(void)
{
	nrerror("usage:\n"
			"RNAfold [-p[0]] [-C] [-T temp] [-4] [-d[2|3]] [-noGU] [-noCloseGU]\n"
			"        [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale]\n"
			"        [-noconv] \n\n"
			" RNAddG special switches:\n"
			"		-f <from_contraint>\n"
			"		-t <to_constraint>\n"	);
}
