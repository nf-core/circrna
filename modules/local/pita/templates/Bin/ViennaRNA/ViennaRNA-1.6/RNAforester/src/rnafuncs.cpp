#include <algorithm>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>

#include "debug.h"
#include "misc.h"
#include "rnafuncs.h"
#include "utils.h"

#ifndef WIN32
#include "config.h"
#endif

#ifndef HAVE_LIBRNA
#undef HAVE_LIBG2
#endif

#include <assert.h>

#ifdef HAVE_LIBG2
#include <g2.h>
#include <g2_PS.h>
#include <g2_FIG.h>

#ifdef HAVE_LIBGD
#include <g2_gd.h>
#endif
#endif

#define COLOR_RED       0
#define COLOR_GREEN     1
#define COLOR_BLUE      2
#define COLOR_YELLOW    3
#define COLOR_MAGENTA   4
#define COLOR_TURKIS    5
#define NUM_COLORS      6

#define COLOR_DEF_RED       1,0,0
#define COLOR_DEF_GREEN     0,1,0
#define COLOR_DEF_BLUE      0,0,1
#define COLOR_DEF_YELLOW    1,1,0
#define COLOR_DEF_MAGENTA   0,1,1
#define COLOR_DEF_TURKIS    1,0,1

bool RNAFuncs::isRNAString(const string &str)
{	
  string::const_iterator it;
  for(it=str.begin(); it!=str.end(); it++)
    {
      switch(tolower(*it))
	{
	case 'a':
	case 'c':
	case 'g':
	case 'u':
	case 't':
	case '.':
	case '-':
	  break;
	default:
	  return false;
	}
    }

  return true;
}


bool RNAFuncs::isViennaString(const string &str, Ulong &basePairCount, Ulong &maxDepth)
{
  long brackets=0;
  maxDepth=0;
  basePairCount=0;       

  string::const_iterator it;
  for(it=str.begin(); it!=str.end(); it++)
    {
      switch(*it)
	{
	case '.':
	  break;
	case '(':
	  brackets++;
	  basePairCount++;
	  if(brackets>0)	
	    maxDepth=max(static_cast<Ulong>(maxDepth),static_cast<Ulong>(brackets));
	  break;
	case ')':
	  if(brackets)
	    {
	      brackets--;
	      break;
	    }
	default:
	  return false;
	}
    }

  // equal nr of opening and closing brackets ?
  if(brackets)
    return false;

  return true;
}

#ifdef HAVE_LIBG2  // This features require the g2 library
void RNAFuncs::drawRNAStructure(const string &seq, const string &structure, const string &filename_prefix, const string &structname, const list<pair<Uint,Uint> > &regions, const SquigglePlotOptions &options)
{
  const double base_fontsize=12;

  float *X, *Y,min_X=0,max_X=0,min_Y=0,max_Y=0;
  Uint i;
  short *pair_table;
  int id_PS,id_FIG=0,id;
  int color_black, *colors;
  double xpos,ypos;
  char buf[10];
  string filename;
  double *points;
  int numPoints;

#ifdef HAVE_LIBGD
  int id_PNG=0,id_JPG=0;
#endif

  X = new float[structure.size()];
  Y = new float[structure.size()];
  points=new double[2*structure.size()];
  colors=new int[NUM_COLORS];

  assert(seq.size() == structure.size());

  pair_table = make_pair_table(structure.c_str());
  i = naview_xy_coordinates(pair_table, X, Y);
  if(i!=structure.size())
    cerr << "strange things happening in squigglePlot ..." << endl;
    
  // scale image
  for(i=0;i<structure.size();i++)
    {
      X[i]*=static_cast<float>(options.scale);
      Y[i]*=static_cast<float>(options.scale);
    }  

  // calculate image dimensions
  for(i=0;i<structure.size();i++)
    {
      min_X=min(min_X,X[i]);
      max_X=max(max_X,X[i]);
      min_Y=min(min_Y,Y[i]);
      max_Y=max(max_Y,Y[i]);
    }

  // add a border to image size
  min_X-=10;
  max_X+=10;
  min_Y-=10;
  max_Y+=10;

  //id_PS  = g2_open_PS("ali.ps", g2_A4, g2_PS_port);
  filename=filename_prefix + ".ps";
  id_PS  = g2_open_EPSF((char*)filename.c_str());
  g2_set_coordinate_system(id_PS,-min_X,-min_Y,1,1);

  if(options.generateFIG)
    {
      filename=filename_prefix + ".fig";
      id_FIG=g2_open_FIG((char*)filename.c_str());
      g2_set_coordinate_system(id_FIG,-min_X,-min_Y,1,1);
    }
#ifdef HAVE_LIBGD
  if(options.generatePNG)
    {
      filename=filename_prefix + ".png";
      id_PNG=g2_open_gd((char*)filename.c_str(),(int)(max_X-min_X),(int)(max_Y-min_Y),g2_gd_png);
      g2_set_coordinate_system(id_PNG,-min_X,-min_Y,1,1);
    }
   if(options.generateJPG)
     {
       filename=filename_prefix + ".jpg"; 
       id_JPG=g2_open_gd((char*)filename.c_str(),(int)(max_X-min_X),(int)(max_Y-min_Y),g2_gd_jpeg);
       g2_set_coordinate_system(id_PS,-min_X,-min_Y,1,1);
     }
#endif

  id     = g2_open_vd();
  g2_attach(id,id_PS);

  if(options.generateFIG)
    {
      g2_attach(id,id_FIG);
    }

#ifdef HAVE_LIBGD
  if(options.generatePNG)
    g2_attach(id,id_PNG);
  if(options.generateJPG)
    g2_attach(id,id_JPG);
#endif

  // cout << "min_X: " << min_X <<",max_X: " << max_X << ",min_Y: " << min_Y << "max_Y: " << max_Y << endl; 
  //  g2_set_coordinate_system(id_PS,595/2.0,842/2.0,0.5,0.5);

  // define colors
  if(options.greyColors)
    {
      color_black=g2_ink(id_PS,0,0,0);
      for(i=0;i<NUM_COLORS;i++)
	colors[i]=g2_ink(id_PS,0,0,0);

      if(options.generateFIG)
	{
	  color_black=g2_ink(id_FIG,0,0,0);
	  for(i=0;i<NUM_COLORS;i++)
	    colors[i]=g2_ink(id_FIG,0,0,0);
	}
      
#ifdef HAVE_LIBGD
      if(options.generatePNG)
	{
	  color_black=g2_ink(id_PNG,0,0,0);
	  for(i=0;i<NUM_COLORS;i++)
	    colors[i]=g2_ink(id_PNG,0,0,0);	  
	}

      if(options.generateJPG)
	{
	  color_black=g2_ink(id_JPG,0,0,0);
	  for(i=0;i<NUM_COLORS;i++)
	    colors[i]=g2_ink(id_JPG,0,0,0);	  
	}     
#endif
    }
  else
    {
      color_black=g2_ink(id,0,0,0);
      colors[COLOR_RED]=g2_ink(id_PS,COLOR_DEF_RED);
      colors[COLOR_GREEN]=g2_ink(id_PS,COLOR_DEF_GREEN);
      colors[COLOR_BLUE]=g2_ink(id_PS,COLOR_DEF_BLUE);
      colors[COLOR_YELLOW]=g2_ink(id_PS,COLOR_DEF_YELLOW);
      colors[COLOR_MAGENTA]=g2_ink(id_PS,COLOR_DEF_MAGENTA);
      colors[COLOR_TURKIS]=g2_ink(id_PS,COLOR_DEF_TURKIS);

      if(options.generateFIG)
	{
	  color_black=g2_ink(id_FIG,0,0,0);
	  colors[COLOR_RED]=g2_ink(id_FIG,COLOR_DEF_RED);
	  colors[COLOR_GREEN]=g2_ink(id_FIG,COLOR_DEF_GREEN);
	  colors[COLOR_BLUE]=g2_ink(id_FIG,COLOR_DEF_BLUE);
	  colors[COLOR_YELLOW]=g2_ink(id_FIG,COLOR_DEF_YELLOW);
	  colors[COLOR_MAGENTA]=g2_ink(id_FIG,COLOR_DEF_MAGENTA);
	  colors[COLOR_TURKIS]=g2_ink(id_FIG,COLOR_DEF_TURKIS);
	}
      
#ifdef HAVE_LIBGD
      if(options.generatePNG)
	{
	  color_black=g2_ink(id_PNG,0,0,0);
	  colors[COLOR_RED]=g2_ink(id_PNG,COLOR_DEF_RED);
	  colors[COLOR_GREEN]=g2_ink(id_PNG,COLOR_DEF_GREEN);
	  colors[COLOR_BLUE]=g2_ink(id_PNG,COLOR_DEF_BLUE);
	  colors[COLOR_YELLOW]=g2_ink(id_PNG,COLOR_DEF_YELLOW);
	  colors[COLOR_MAGENTA]=g2_ink(id_PNG,COLOR_DEF_MAGENTA);
	  colors[COLOR_TURKIS]=g2_ink(id_PNG,COLOR_DEF_TURKIS);
	}

      if(options.generateJPG)
	{
	  color_black=g2_ink(id_JPG,0,0,0);
	  colors[COLOR_RED]=g2_ink(id_JPG,COLOR_DEF_RED);
	  colors[COLOR_GREEN]=g2_ink(id_JPG,COLOR_DEF_GREEN);
	  colors[COLOR_BLUE]=g2_ink(id_JPG,COLOR_DEF_BLUE);
	  colors[COLOR_YELLOW]=g2_ink(id_JPG,COLOR_DEF_YELLOW);
	  colors[COLOR_MAGENTA]=g2_ink(id_JPG,COLOR_DEF_MAGENTA);
	  colors[COLOR_TURKIS]=g2_ink(id_JPG,COLOR_DEF_TURKIS);
	}
#endif
    }
  
  for(i=0;i<structure.size();i++)
    {
      // connection to next base
      if(i<structure.size()-1)
	{
	 // if the connected bases are only in one of the structures
	  // use the appropriate color
	  g2_line(id,X[i],Y[i],X[i+1],Y[i+1]);
	  
	  // draw circles at line endpoints
	  g2_filled_circle(id,X[i],Y[i],0.7*options.scale);       // circles are drawn twice, but thats ok ...
	}
    }

  // draw pairings
  // !!! pair_table indexing begins at 1 !!!
  for(i=0;i<structure.size();i++)
    {
	if((unsigned short)pair_table[i+1]>i+1)
	  {	    	    
	    // pairs in both structures
	    g2_line(id,X[i],Y[i],X[pair_table[i+1]-1],Y[pair_table[i+1]-1]);
	  }
    }

  // draw regions
  g2_set_line_width(id,0.4);
  list<pair<Uint,Uint> >::const_iterator it;  
  Uint regionNr=0;

  for(it=regions.begin();it!=regions.end();it++)
    {
      double center_x=0,center_y=0;

      // fill coordinate list
	for(i=0;i < it->second;i++)
	  {
	    points[2*i]=X[it->first+i];
	    points[2*i+1]=Y[it->first+i];

	    center_x+=X[it->first+i];
	    center_y+=Y[it->first+i];
	  }
	numPoints=it->second-1;
	center_x/=numPoints;   // center of gravity
	center_y/=numPoints;
	
	g2_pen(id,colors[regionNr % NUM_COLORS]);
	g2_poly_line(id,numPoints,points);
	sprintf(buf,"%d",regionNr+1);
	g2_string(id,center_x,center_y,buf);   // draw region number

	regionNr++;
    }
  
  // mark 5' end

  g2_pen(id,color_black);
  g2_string(id,X[0]-20,Y[0],"5'");

  g2_set_font_size(id,base_fontsize*options.scale);
  g2_set_line_width(id,0.2);
  
  // draw sequence
  for(i=0;i<structure.size();i++)
    {
     g2_pen(id,color_black); // match 
     xpos=X[i]-(base_fontsize*options.scale)/2;
     ypos=Y[i]-4;
     sprintf(buf,"%c",seq[i]);
     g2_string(id,xpos,ypos,buf);

     /*
     if(!options.hideBaseNumbers)
       {
	 // draw base number
	 if(basenr % options.baseNumInterval == 0)
	   {	 
	     sprintf(buf,"%d",basenr);
	     g2_string(id,xpos-20,ypos,buf);
	   }
	   }*/     
    }

    // draw structure name
    g2_string(id,min_X,max_Y-10,(char*)structname.c_str());
        
    g2_flush(id);
    g2_close(id);

    free(pair_table);
    DELETE(X);
    DELETE(Y);
    DELETE(points);
    DELETE(colors);
}
#endif

#ifdef HAVE_LIBG2  // This features require the g2 library
void RNAFuncs::drawRNAAlignment(const string &structure, const string &altStructure, const string &seq1, const string &seq2, const string &strname1, const string &strname2, const string &filename_prefix, const bool atX, const SquigglePlotOptions &options)
{
  const double base_fontsize=12;

  string base,structname;
  float *X, *Y,min_X=0,max_X=0,min_Y=0,max_Y=0;
  Uint i;
  short *pair_table,*alt_pair_table;
  int id_PS,id;
  int ps_color_black,ps_color_red,ps_color_blue,ps_color_green;
  Uint basenr_x=0, basenr_y=0;
  double xpos,ypos;
  char buf[10];
  bool isDel,isIns;
  string filename;
#ifdef HAVE_LIBGD
  int id_PNG=0,id_JPG=0;
  int png_color_black,png_color_red,png_color_blue,png_color_green;
  int jpg_color_black,jpg_color_red,jpg_color_blue,jpg_color_green;
#endif

  X = new float[structure.size()];
  Y = new float[structure.size()];

  assert(seq1.size() == seq2.size());
  assert(seq1.size() == structure.size());
  assert(seq1.size() == altStructure.size());

  pair_table = make_pair_table(structure.c_str());
  alt_pair_table = make_pair_table(altStructure.c_str());
  i = naview_xy_coordinates(pair_table, X, Y);
  if(i!=structure.size())
    cerr << "strange things happening in squigglePlot ..." << endl;
    
  // scale image
  for(i=0;i<structure.size();i++)
    {
      X[i]*=static_cast<float>(options.scale);
      Y[i]*=static_cast<float>(options.scale);
    }  

  // calculate image dimensions
  for(i=0;i<structure.size();i++)
    {
      min_X=min(min_X,X[i]);
      max_X=max(max_X,X[i]);
      min_Y=min(min_Y,Y[i]);
      max_Y=max(max_Y,Y[i]);
    }

  // add a border to image size
  min_X-=10;
  max_X+=10;
  min_Y-=10;
  max_Y+=10;

  //id_PS  = g2_open_PS("ali.ps", g2_A4, g2_PS_port);
  filename=filename_prefix + ".ps";
  id_PS  = g2_open_EPSF((char*)filename.c_str());
  g2_set_coordinate_system(id_PS,-min_X,-min_Y,1,1);

#ifdef HAVE_LIBGD
  if(options.generatePNG)
    {
      filename=filename_prefix + ".png";
      id_PNG=g2_open_gd((char*)filename.c_str(),(int)(max_X-min_X),(int)(max_Y-min_Y),g2_gd_png);
      g2_set_coordinate_system(id_PNG,-min_X,-min_Y,1,1);
    }
   if(options.generateJPG)
     {
       filename=filename_prefix + ".jpg"; 
       id_JPG=g2_open_gd((char*)filename.c_str(),(int)(max_X-min_X),(int)(max_Y-min_Y),g2_gd_jpeg);
       g2_set_coordinate_system(id_PS,-min_X,-min_Y,1,1);
     }
#endif

  id     = g2_open_vd();
  g2_attach(id,id_PS);
#ifdef HAVE_LIBGD
  if(options.generatePNG)
    g2_attach(id,id_PNG);
  if(options.generateJPG)
    g2_attach(id,id_JPG);
#endif

  // cout << "min_X: " << min_X <<",max_X: " << max_X << ",min_Y: " << min_Y << "max_Y: " << max_Y << endl; 
  //  g2_set_coordinate_system(id_PS,595/2.0,842/2.0,0.5,0.5);
  g2_set_line_width(id,0.2);

  // mark 5' end
  g2_string(id,X[0]-20,Y[0],"5'");

  // define colors
  if(options.greyColors)
    {
      ps_color_black=g2_ink(id_PS,0,0,0);
      ps_color_red=g2_ink(id_PS,0,0,0);
      ps_color_blue=g2_ink(id_PS,0.7,0.7,0.7);
      ps_color_green=g2_ink(id_PS,0.4,0.4,0.4);

#ifdef HAVE_LIBGD
      if(options.generatePNG)
	{
	  png_color_black=g2_ink(id_PNG,0,0,0);
	  png_color_red=g2_ink(id_PNG,0,0,0);
	  png_color_blue=g2_ink(id_PNG,0.7,0.7,0.7);
	  png_color_green=g2_ink(id_PNG,0.4,0.4,0.4);
	}

       if(options.generateJPG)
	{
	  jpg_color_black=g2_ink(id_JPG,0,0,0);
	  jpg_color_red=g2_ink(id_JPG,0,0,0);
	  jpg_color_blue=g2_ink(id_JPG,0.7,0.7,0.7);
	  jpg_color_green=g2_ink(id_JPG,0.4,0.4,0.4);  
	}     
#endif
    }
  else
    {
      ps_color_black=g2_ink(id_PS,0,0,0);
      ps_color_red=g2_ink(id_PS,1,0,0);
      ps_color_blue=g2_ink(id_PS,0,0,0.5);
      ps_color_green=g2_ink(id_PS,0,0.5,0);

#ifdef HAVE_LIBGD
      if(options.generatePNG)
	{
	  png_color_black=g2_ink(id_PNG,0,0,0);
	  png_color_red=g2_ink(id_PNG,1,0,0);
	  png_color_blue=g2_ink(id_PNG,0,0,0.5);
	  png_color_green=g2_ink(id_PNG,0,0.5,0);
	}

      if(options.generateJPG)
	{
	  jpg_color_black=g2_ink(id_JPG,0,0,0);
	  jpg_color_red=g2_ink(id_JPG,1,0,0);
	  jpg_color_blue=g2_ink(id_JPG,0,0,0.5);
	  jpg_color_green=g2_ink(id_JPG,0,0.5,0);
	}
#endif
    }  

  // draw sequence
  g2_set_font_size(id,base_fontsize*options.scale);
    for(i=0;i<structure.size();i++)
   {
     isDel = false;
     isIns = false;

     //base
     g2_pen(id,ps_color_black); // match 
     base = "";
     if(seq1[i]!='-')
       {
	 base += seq1[i];
	 basenr_x++;
       }
     else
       {
	 g2_pen(id,ps_color_green); // insertion
	 isIns=true;
       }
     
     if(seq2[i]!='-')
       {
	 base += seq2[i];
	 basenr_y++;
       }
     else
       {
	 g2_pen(id,ps_color_blue); // deletion 
	 isDel=true;
       }
     
     if(base.size()==2 && base[0]!=base[1])
       g2_pen(id,ps_color_red); // mismatch  
     else
       base = base[0];          // show duplicate base only once
     
     xpos=X[i]-base.length()*(base_fontsize*options.scale)/2;
     ypos=Y[i]-4;
     g2_string(id,xpos,ypos,(char*)base.c_str());

     if(!options.hideBaseNumbers)
       {
	 // draw base number
	 if(!isIns && basenr_x % options.baseNumInterval == 0)
	   {	 
	     g2_pen(id,ps_color_blue);
	     sprintf(buf,"%d",basenr_x);
	     g2_string(id,xpos-20,ypos,buf);
	   }
	 if(!isDel && basenr_y % options.baseNumInterval == 0)
	   {	 
	     g2_pen(id,ps_color_green);
	     sprintf(buf,"%d",basenr_y);
	     g2_string(id,xpos+20,ypos,buf);
	   }
       }
     
     // connection to next base
     if(i<structure.size()-1)
       {
	 // if the connected bases are only in one of the structures
	 // use the appropriate color
	 if(seq1[i]=='-' && seq1[i+1]=='-')
	   g2_pen(id,ps_color_green);
	 else
	   if(seq2[i]=='-' && seq2[i+1]=='-')
	     g2_pen(id,ps_color_blue);
	   else	      
	     g2_pen(id,ps_color_black);

	 g2_line(id,X[i],Y[i],X[i+1],Y[i+1]);

	 // draw circles at line endpoints
	 if(seq1[i]=='-')
	   g2_pen(id,ps_color_green);
	 else
	   if(seq2[i]=='-')
	     g2_pen(id,ps_color_blue);
	   else	      
	     g2_pen(id,ps_color_black); 
	 g2_filled_circle(id,X[i],Y[i],0.7*options.scale);       // circles are drawn twice, but thats ok ...

	 if(seq1[i+1]=='-')
	   g2_pen(id,ps_color_green);
	 else
	   if(seq2[i+1]=='-')
	     g2_pen(id,ps_color_blue);
	   else	      
	     g2_pen(id,ps_color_black);
	 g2_filled_circle(id,X[i+1],Y[i+1],0.7*options.scale);
       }
   }
   
    // draw pairings
    // !!! pair_table indexing begins at 1 !!!
    for(i=0;i<structure.size();i++)
      {
	if((unsigned short)pair_table[i+1]>i+1 && (unsigned short)alt_pair_table[i+1]>i+1)
	  {	    	    
	    // pairs in both structures
	    g2_pen(id,ps_color_black);	   	    
	    g2_line(id,X[i],Y[i],X[pair_table[i+1]-1],Y[pair_table[i+1]-1]);
	  }
	else
	  {
	    if((unsigned short)pair_table[i+1]>i+1 && (unsigned short)alt_pair_table[i+1]<=i+1)
	      {	    	    
		// pairs only in first structure
		if(atX)
		  g2_pen(id,ps_color_blue);
		else
		  g2_pen(id,ps_color_green);

		g2_line(id,X[i],Y[i],X[pair_table[i+1]-1],Y[pair_table[i+1]-1]);
	      }
	    else
	      {
		if((unsigned short)pair_table[i+1]<=i+1 && (unsigned short)alt_pair_table[i+1]>i+1)
		  {	    	    
		    // pairs only in second structure
		    if(atX)
		      g2_pen(id,ps_color_green);	   	    
		    else
		      g2_pen(id,ps_color_blue);

		    double dashes=2.0;
		    g2_set_dash(id,1,&dashes);
		    g2_line(id,X[i],Y[i],X[alt_pair_table[i+1]-1],Y[alt_pair_table[i+1]-1]);
		    dashes=0;
		    g2_set_dash(id,1,&dashes);
		  }
	      }
	  }
      }


    // draw structure names
    // x-at-y or y-at-x
    if(atX)
      {
	g2_pen(id,ps_color_green);
	g2_string(id,min_X,max_Y-10,(char*)strname2.c_str());
	g2_pen(id,ps_color_black);
	g2_string(id,min_X,max_Y-20,"at");
	g2_pen(id,ps_color_blue);
	g2_string(id,min_X,max_Y-30,(char*)strname1.c_str());
      }
    else
      {
	g2_pen(id,ps_color_blue);
	g2_string(id,min_X,max_Y-10,(char*)strname1.c_str());
	g2_pen(id,ps_color_black);
	g2_string(id,min_X,max_Y-20,"at");
	g2_pen(id,ps_color_green);
	g2_string(id,min_X,max_Y-30,(char*)strname2.c_str());
      }	   

    g2_flush(id);
    g2_close(id);

    free(pair_table);
    free(alt_pair_table);
    DELETE(X);
    DELETE(Y);
}
#endif

#ifdef HAVE_LIBRNA
void RNAFuncs::generateRNAAlignmentXML(const string &structure, const string &altStructure, const string &seq1, const string &seq2, const string &strname1, const string &strname2, ostream &s)
{
  string base;
  Uint i;
  float *X, *Y;
  short *pair_table,*alt_pair_table;
  Uint basenr_x=1, basenr_y=1;
  string filename;

  X = new float[structure.size()];
  Y = new float[structure.size()];

  assert(seq1.size() == seq2.size());
  assert(seq1.size() == structure.size());
  assert(seq1.size() == altStructure.size());

  // calculate coordinates
  pair_table = make_pair_table(structure.c_str());
  alt_pair_table = make_pair_table(altStructure.c_str());
  i = naview_xy_coordinates(pair_table, X, Y);
  if(i!=structure.size())
    cerr << "strange things happening in squigglePlot ..." << endl;
    
  // generate XML
  s << "<alignment xname=\"" << strname1 << "\" yname=\"" << strname2 << "\">" << endl;

  // seqalignment
  s << "  <seqalignment>" << endl;
  for(i=0;i<seq1.length();i++)
    {
      s << "    <base id=\"" << i+1 << "\" xbasenr=\"" << basenr_x << "\" ybasenr=\"" << basenr_y << "\" xbase=\"" << seq1[i] <<"\" ybase=\"" << seq2[i] << "\" />" << endl;
      if(seq1[i]!='-')
	basenr_x++;

      if(seq2[i]!='-')
	basenr_y++;      
    }
  s << "  </seqalignment>" << endl;

  // pairs
  s << "  <pairs>" << endl;
  for(i=0;i<structure.size();i++)
    {
      // both
      if((unsigned short)pair_table[i+1]>i+1 && (unsigned short)alt_pair_table[i+1]>i+1) 
	{
	  s << "    <pair drawbaseid1=\"" << i+1 << "\" drawbaseid2=\"" << pair_table[i+1] << "\"/>" << endl;
	}
      else
	{
	  if((unsigned short)pair_table[i+1]>i+1)
	    {
	      s << "    <pair drawbaseid1=\"" << i+1 << "\" drawbaseid2=\"" << pair_table[i+1] << "\" structid=\"1\"/>" << endl;
	    }
	  else
	    {
	    if((unsigned short)alt_pair_table[i+1]>i+1)
	      s << "    <pair drawbaseid1=\"" << i+1 << "\" drawbaseid2=\"" << alt_pair_table[i+1] << "\" structid=\"2\"/>" << endl;
	    }
	}
    }
  s << "  </pairs>" << endl;

  // coordinates
  // x structure
  s << "  <structure id=\"1\">" << endl;
  s << "    <drawbases>" << endl;
  for(i=0;i<seq1.length();i++)
    {
      s << "      <drawbase id=\"" << i+1 << "\" xcoor=\"" << X[i] << "\" ycoor=\"" << -Y[i] << "\"/>" << endl;
    }
  s << "    </drawbases>" << endl;
  s << "  </structure>" << endl;
  
  DELETE(X);
  DELETE(Y);
  X = new float[altStructure.size()];
  Y = new float[altStructure.size()];

  // y structure
  i = naview_xy_coordinates(alt_pair_table, X, Y);
  if(i!=structure.size())
    cerr << "strange things happening in squigglePlot ..." << endl;

  s << "  <structure id=\"2\">" << endl;
  s << "    <drawbases>" << endl;
  for(i=0;i<seq1.length();i++)
    {
      s << "      <drawbase id=\"" << i+1 << "\" xcoor=\"" << X[i] << "\" ycoor=\"" << -Y[i] << "\"/>" << endl;
    }
  s << "    </drawbases>" << endl;
  s << "  </structure>" << endl;
  s << "</alignment>" << endl;

  free(pair_table);
  free(alt_pair_table);
  DELETE(X);
  DELETE(Y);
}
#endif

void RNAFuncs::printAli(const string &name1, const string &name2, const string &seq1, const string &seq2, const string &str1, const string &str2)
{
  Uint i;
  string info;

  for(i=0;i<seq1.length();i++)		 
    info += seq1[i]==seq2[i] ? "*" : " ";		  		

  for(i=0;i<seq1.length();i+=55)
    {
      cout << setw(20) << setfill(' ') << left << name1.substr(0,20) << setw(5) << " " << seq1.substr(i,55) << endl;
      cout << setw(20) << setfill(' ') << left << name2.substr(0,20) << setw(5) << " " << seq2.substr(i,55) << endl;
      cout << setw(25) << " " << info.substr(i,55) << endl;
    }

  cout << endl;

  info="";
  for(i=0;i<str1.length();i++)		 
    info += str1[i]==str2[i] ? "*" : " ";

  for(i=0;i<str1.length();i+=55)
    {
      cout << setw(20) << setfill(' ') << left << name1.substr(0,20) << setw(5) << " " << str1.substr(i,55) << endl;
      cout << setw(20) << setfill(' ') << left << name2.substr(0,20) << setw(5) << " " << str2.substr(i,55) << endl;
      cout << setw(25) << " " << info.substr(i,55) << endl;
    }

  cout << endl;

}

Uint RNAFuncs::treeSize(const string &viennaStr)
{
   Ulong basePairCount=0,maxDepth=0;

   RNAFuncs::isViennaString(viennaStr,basePairCount,maxDepth);
   // the size of the forests is determined by the number of bases and basepairs	
   return viennaStr.length() + basePairCount;  
};
