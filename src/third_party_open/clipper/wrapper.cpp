#include "stdafx.h"
#include "clipper.hpp" 
#include <iostream>
#include <fstream>
using namespace ClipperLib;
using namespace std;
//from clipper.hpp ... 
//typedef long long long64; 
//struct IntPoint {long64 X; long64 Y;}; 
//typedef std::vector<IntPoint> Polygon; 
//typedef std::vector<Polygon> Polygons; 
extern "C" int wrapclip() {
	//set up the subject and clip polygons ... 
	long64 fact;
	Polygons subj(1), clip(1), solution;
	fact = 1;
	subj[0].push_back(IntPoint(180*fact,200*fact));
	subj[0].push_back(IntPoint(260*fact,200*fact));
	subj[0].push_back(IntPoint(260*fact,150*fact));
	subj[0].push_back(IntPoint(180*fact,150*fact));
	 for ( size_t counter = 0; counter <= solution[0].size()-1; counter++)  // size_t intstead of int otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
		 outputFILE <<  solution[0][counter].X/fact << ' ' <<   solution[0][counter].Y/fact <<"\n" ;

	//subj[1].push_back(IntPoint(215*fact,160*fact));
	//subj[1].push_back(IntPoint(230*fact,190*fact));
	//subj[1].push_back(IntPoint(200*fact,190*fact));

	clip[0].push_back(IntPoint(190*fact,210*fact));
	clip[0].push_back(IntPoint(240*fact,210*fact));
	clip[0].push_back(IntPoint(240*fact,130*fact));
	clip[0].push_back(IntPoint(190*fact,130*fact));

	//get the intersection of the subject and clip polygons ... 
	Clipper clpr; 
	clpr.AddPolygons(subj, ptSubject); 
	clpr.AddPolygons(clip, ptClip); 
	  clpr.Execute(ctIntersection, solution, pftEvenOdd, pftEvenOdd); 
	//clpr.Execute(ctUnion, solution, pftEvenOdd, pftEvenOdd);  
	//printf("%d %d\n" ,subj[0].push_back(IntPoint));
	// printf("%lld %llu\d" ,subj[0].IntPoint[]);
	// cout <<  solution[0].size()  <<"\n"  ;
	// cout <<  solution[0][1] ;

	// cout <<"\n"  ;
	 ofstream outputFILE; 
	 ofstream fileSUBJ; 
	 ofstream fileCLIP; 
     outputFILE.open ("solution.txt");
     fileSUBJ.open ("SUBJ.txt");
     fileCLIP.open ("CLIP.txt");
 
	 for ( size_t counter = 0; counter <= solution[0].size()-1; counter++)  // size_t intstead of int otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
		 outputFILE <<  solution[0][counter].X/fact << ' ' <<   solution[0][counter].Y/fact <<"\n" ;
	       
	 for ( size_t counter = 0; counter <= subj[0].size()-1; counter++)  // size_t intstead of int otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
		 fileSUBJ <<  subj[0][counter].X/fact << ' ' <<   subj[0][counter].Y/fact <<"\n" ;
	 for ( size_t counter = 0; counter <= clip[0].size()-1; counter++)  // size_t intstead of int otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
		 fileCLIP <<  clip[0][counter].X/fact << ' ' <<   clip[0][counter].Y/fact <<"\n" ;

	//finally draw the intersection polygons ... 
	//DrawPolygons(img->Bitmap, sol, 0x40808080); 
	return 0;
}