/*
 * RNA folding statistics program Based on RNAsubopt.
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "part_func.h"
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "subopt.h"
extern void     read_parameter_file(const char fname[]);
extern int      st_back;
/* @unused@ */
static char UNUSED rcsid[] = "$Id: RNAstats.c,v 1.0 2006/01/02 Exp $";

#define ENTROPY(X)	(X==0 ? 0 : (-X*log2f(X)))

#define PRIVATE static

static char     scale[] = "....,....1....,....2....,....3....,....4"
						  "....,....5....,....6....,....7....,....8";

extern double   print_energy;
PRIVATE void    usage(void);
float *bp_entropy(int length);
extern char    *pbacktrack(char *sequence);
/*--------------------------------------------------------------------------*/

int 
main(int argc, char *argv[])
{
	char           *line;
	char           *sequence;
	char		   *ss = NULL;
	char           *ParamFile = NULL;
	char           *ns_bases = NULL, *c;
	int             i, length, l, sym, r;
	double          deltaf, deltap = 0;
	int             delta = 100;
	int             statType = 1;
	int				slideOffset = 1;
	int				slideWindow = 200;
	
	do_backtrack = 1;
	dangles = 2;
		
	/* Read command line arguments ********************************************/
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-')
			switch (argv[i][1]) {
			case 'T':
				if (argv[i][2] != '\0')
					usage();
				if (i == argc - 1)
					usage();
				r = sscanf(argv[++i], "%lf", &temperature);
				if (r != 1)
					usage();
				break;
			case 'S':
				if (argv[i][2] != '\0')
					usage();
				if (i == argc - 1)
					usage();
				(void) sscanf(argv[++i], "%d", &statType);
				break;
			case 'o':
				if (argv[i][2] != '\0')
					usage();
				if (i == argc - 1)
					usage();
				(void) sscanf(argv[++i], "%d", &slideOffset);
				break;
			case 'w':
				if (argv[i][2] != '\0')
					usage();
				if (i == argc - 1)
					usage();
				(void) sscanf(argv[++i], "%d", &slideWindow);
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
				break;
			case '4':
				tetra_loop = 0;
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
			case 's':
				subopt_sorted = 1;
				break;
			case 'l':
				if (strcmp(argv[i], "-logML") == 0) {
					logML = 1;
					break;
				} else
					usage();
				break;
			case 'e':
				if (i >= argc - 1)
					usage();
				if (strcmp(argv[i], "-ep") == 0)
					r = sscanf(argv[++i], "%lf", &deltap);
				else {
					r = sscanf(argv[++i], "%lf", &deltaf);
					delta = (int) (0.1 + deltaf * 100);
				}
				if (r != 1)
					usage();
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
		while (*c) {
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


	/* Main Loop. Process each line as a sequence *****************************/

	int		mergePos = (int)(floor(((float)(slideWindow - slideOffset)) / 2)) + 1;
	
	do {			
	
		/* Read sequence, make sure it's long enough and allocate memory for structure */
		
		if ((line = get_line(stdin)) == NULL)
			break;
	
		if (slideWindow > strlen(line))
		{
			fprintf(stderr, "Sequence length (%d) should be at least window-length (%d)", strlen(line), slideWindow);
			nrerror ("Quitting.");
		}
		
		/* printf ("Sliding window %d\n", slideWindow); */
		
		sequence = (char *) space(strlen(line) + 1);
		(void) sscanf(line, "%s", sequence);
		free(line);
		
		length = (int) strlen(sequence);

		for (l = 0; l < length; l++)
			sequence[l] = toupper(sequence[l]);
		
		
		/* Some deltap calculations */
		
		if (logML != 0 || dangles == 1 || dangles == 3)
			if (deltap <= 0)
				deltap = delta / 100. + 0.001;
		if (deltap > 0)
			print_energy = deltap;
		
		int             i, j;
		double          mfe, kT;

		int				startPos;
		char			rememberBp;
		
		st_back = 1;
			
		/* Construct the sub-sequence to be processed, based on the window size
		   and offset. Then calculate the required statistics for each such sub-sequence */

		ss = (char *) space((unsigned)(slideWindow + 1));
		   
		for (startPos = 0; startPos <= (length - slideWindow); startPos += slideOffset)
		{	
			/* Crop the sequence at window size */
			
			rememberBp = sequence[startPos + slideWindow];
			sequence[startPos + slideWindow] = 0;
			
			/* Initial folding of the sequence */
				
			mfe = fold((char*)(sequence + startPos), ss);
			
			kT = (temperature + 273.15) * 1.98717 / 1000.;	/* in Kcal */
			pf_scale = exp(-(1.03 * mfe) / kT / length);

			(void) pf_fold((char*)(sequence + startPos), NULL);
			
			float 	*bpProfile;
			float	*bpEntropy;
			float	basePairedP;
			char	na;
			
			/* Use the Make_bp_profile function to calculate the probability of
			   each bp to be paired. Then either print out this probability 
			   or the entropy of this probability */

			bpProfile = (float *)Make_bp_profile (slideWindow);
			bpEntropy = (float *)bp_entropy (slideWindow);
	
			/* P[i*3+0] unpaired, P[i*3+1] upstream, P[i*3+2] downstream p */
			
			/* First time only -- pad with data before mergePos */
			if (!startPos)
			{
				for (j = 1; j < mergePos; j++)
				{
					basePairedP = (1-bpProfile[j*3]);
					na = sequence[startPos + j - 1];
					
					printf("%d\t%d\tPairness\t%.2f\t%c\n", j, j, basePairedP, na);
					printf("%d\t%d\tPairness_Entropy\t%.2f\t%c\n", j, j, ENTROPY(basePairedP) + ENTROPY((1-basePairedP)), na);
					printf("%d\t%d\tEntropy\t%.2f\t%c\n", j, j,  bpEntropy[j], na);
					/* Not printing dG */
				}
			}
				
			/* Always -- print center data */
			for (j = mergePos; j < mergePos + slideOffset; j++)
			{
				basePairedP = (1-bpProfile[j*3]);
				na = sequence[startPos + j - 1];
				
				printf("%d\t%d\tPairness\t%.2f\t%c\n", j+startPos, j+startPos, basePairedP, na);
				printf("%d\t%d\tPairness_Entropy\t%.2f\t$c\n", j+startPos, j+startPos, ENTROPY(basePairedP) + ENTROPY((1-basePairedP)), na);
				printf("%d\t%d\tEntropy\t%.2f\t%c\n", j+startPos, j+startPos,  bpEntropy[j], na);
			}
			
			/* Always -- print dG data (but only once) */
			printf("%d\t%d\tdG\t%.2f\n", startPos+mergePos, (startPos+mergePos+slideOffset-1), -mfe);

			/* Last time only -- pad with data after mergePos */
			if (startPos + slideOffset > length - slideWindow)
			{
				for (j = mergePos + slideOffset; j < slideWindow+1; j++)
				{
					basePairedP = (1-bpProfile[j*3]);
					na = sequence[startPos + j - 1];
					
					printf("%d\t%d\tPairness\t%.2f\t%c\n", j+startPos, j+startPos, basePairedP, na);
					printf("%d\t%d\tPairness_Entropy\t%.2f\t%c\n", j+startPos, j+startPos, ENTROPY(basePairedP) + ENTROPY((1-basePairedP)), na);
					printf("%d\t%d\tEntropy\t%.2f\t%c\n", j+startPos, j+startPos,  bpEntropy[j], na);
				}
			}
				
			free (bpProfile);
			free (bpEntropy);
	
			/* Restore the bp at the position where we placed the "0" beofre */
			sequence[startPos + slideWindow] = rememberBp;
		}
			
		(void) fflush(stdout);

		free_pf_arrays();
		free(ss);
		free(sequence);
		
	} while (1);
	return 0;
}

/******************************************************************************/

float *bp_entropy(int length)
{
	int i,j;
	int L=3;
	float *P;
	float entropy;
	
	P =  (float *) space((length)*sizeof(float));

	for (i=1; i<length+1; i++)
	{
		entropy = 0;
		for (j=1; j<i; j++)
		{
			entropy += ENTROPY(pr[iindx[j]-i]);
		}
			
    	for (j=i+1; j<=length; j++)
    	{
			entropy += ENTROPY(pr[iindx[i]-j]);
		}
		
		P[i-1] = entropy;
	}
	
	return (float *) P;
}
 
 
PRIVATE void 
usage(void)
{
	nrerror("usage: \n\n"
		
		 "RNAstats [-e range] [-ep prange] [-s] [-logML]\n"
         "[-C] [-T temp] [-4] [-d[2]] [-noGU] [-noCloseGU]\n"
         "[-noLP] [-P paramfile] [-nsp pairs]\n"
         "[-S statType] [-w windowSize] [-o offset] \n\n"

        "-S statistics type:\n"
        "    1 = Pairing probability\n"
        "    2 = Pairness entropy\n"
        "    3 = Base entropy\n"
        "    4 = delta G\n\n"
            
        " -w windowSize\n"
        " 	Size of sliding window used when folding the sequence\n\n"
         	
        " -o offset\n"
        " 	Offset (in bases) in which the sliding window should be moved in\n"
        " 	every step.\n\n");

}
