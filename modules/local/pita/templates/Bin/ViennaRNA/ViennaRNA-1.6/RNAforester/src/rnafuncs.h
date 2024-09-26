#ifndef _RNA_FUNCS_H
#define _RNA_FUNCS_H

#include <list>
#include <string>

#include "types.h"

using namespace std;
	
class RNAFuncs
{
 public:
  struct SquigglePlotOptions
  {
	  bool hideBaseNumbers;
	  Uint baseNumInterval;
	  bool greyColors;
	  bool generatePNG;
	  bool generateJPG;
          bool generateFIG;
	  double scale;
  };

	static bool isRNAString(const string &str);
	static bool isViennaString(const string &str, Ulong &basePairCount, Ulong &maxDepth);
	static void drawRNAStructure(const string &seq, const string &structure, const string &filename_prefix, const string &structname, const list<pair<Uint,Uint> > &regions, const SquigglePlotOptions &options);
    static void drawRNAAlignment(const string &structure, const string &altStructure,  const string &seq1, const string &seq2, const string &strname1, const string &strname2, const string &filename_prefix, const bool atX, const SquigglePlotOptions &options);
	static void generateRNAAlignmentXML(const string &structure, const string &altStructure, const string &seq1, const string &seq2, const string &strname1, const string &strname2, ostream &s);
	static void printAli(const string &name1, const string &name2, const string &seq1, const string &seq2, const string &str1, const string &str2);
	static Uint treeSize(const string &viennaStr);

};

#endif

