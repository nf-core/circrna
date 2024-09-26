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
/* @unused@ */
static char     rcsid[] = "$Id: RNAfold.c,v 1.18 2005/10/30 20:07:12 ivo Exp $";

#define PRIVATE static

static char     scale1[] = "....,....1....,....2....,....3....,....4";
static char     scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE void    usage(void);

/*--------------------------------------------------------------------------*/

int 
main(int argc, char *argv[])
{
	char           *string, *line;
	char           *structure = NULL, *cstruc = NULL;
	char           *ParamFile = NULL;
	char           *ns_bases = NULL, *c;
	int             i, length, l, sym, r;
	double          energy, min_en;
	double          kT, sfact = 1.07;
	int             noconv = 0;
	
	float			dG2;
	double			dG0, dG1, ddG, P;

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
	
		/* Read and prepare sequence and constraint ***************************/
		
		if ((line = get_line(stdin)) == NULL)
			break;

		string = (char *) space(strlen(line) + 1);
		(void) sscanf(line, "%s", string);
		free(line);
		length = (int) strlen(string);

		structure = (char *) space((unsigned) length + 1);

		/* Always read folding constraint */
		if ((cstruc = get_line(stdin)) == NULL)
			fprintf(stderr, "RNAdGG error: constraints missing\n");

		if ((line = get_line(stdin)) == NULL)
			fprintf (stderr, "RNAdGG error: dG2 missing\n");
		else
		{
			sscanf (line, "%f", &dG2);
			free (line);
		}
		
		for (l = 0; l < length; l++) {
			string[l] = toupper(string[l]);
			if (!noconv && string[l] == 'T')
				string[l] = 'U';
		}

		/* Fold with no constraint ********************************************/

		min_en = fold(string, structure);		
		printf ("%s\t", structure);
		
		kT = (temperature + 273.15) * 1.98717 / 1000.;	/* in Kcal */
		pf_scale = exp(-(sfact * min_en) / kT / length);
		init_pf_fold(length);

		dG0 = pf_fold(string, structure);
		
		free_pf_arrays();

		/* Fold with constraint ***********************************************/

		fold_constrained = 1;
		strncpy(structure, cstruc, length + 1);

		min_en = fold(string, structure);		
		printf ("%s\t", structure);
		
		pf_scale = exp(-(sfact * min_en) / kT / length);		
		init_pf_fold(length);

		strncpy(structure, cstruc, length + 1);
		dG1 = pf_fold(string, structure);

		free_pf_arrays();
		
		/* Calculate and output ***********************************************/
		
		ddG = dG1 + dG2 - dG0;
		P = 1 / ( 1 + exp(ddG / kT) );	/* Probability of being bound */
		
		printf ("%5.2f\t%5.2f\t%5.2f\t%5.2f\t%g\n", dG0, dG1, dG2, ddG, P);
		
		/* Free everything */
		
		if (cstruc != NULL)
			free(cstruc);
		if (length >= 2000)
			free(base_pair);
		(void) fflush(stdout);
		free(string);
		free(structure);
	} while (1);
	return 0;
}

PRIVATE void 
usage(void)
{
	nrerror("usage:\n"
		"RNAfold [-p[0]] [-C] [-T temp] [-4] [-d[2|3]] [-noGU] [-noCloseGU]\n"
	"        [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale]\n"
		"        [-noconv] \n");
}
