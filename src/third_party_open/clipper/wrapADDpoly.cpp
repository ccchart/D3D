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
extern "C" void wrapAddPolySubj( int Nnod[1], double xFOR[],double yFOR[], double absMAXx, double absMAXy) {
   // declarations
   int NnodCpp ;
   NnodCpp = Nnod[0];
   Polygons subj(1); //I add 1 polygon for each call of this subroutine
   double factx = fact/absMAXx;
   double facty = fact/absMAXy;

    for ( int counter = 0; counter <= NnodCpp-1; counter++) { 
     subj[0].push_back(IntPoint(long64(xFOR[counter]*factx+0.5),long64(yFOR[counter]*facty+0.5))); }//add +0.5 couse conversion from double to long long implies truncation, i.e. 29.7 becomes 29
   //add the subject polygon to Clipper ... 
   Clipper clpr; 
   clpr.AddPolygons(subj, ptSubject); 
}

 
extern "C" void wrapAddPolyClip( int Nnod[1], double xFOR[],double yFOR[], double absMAXx, double absMAXy) {
   // declarations
   int NnodCpp ;
   NnodCpp = Nnod[1];
   Polygons clip(NnodCpp);
   double factx = fact/absMAXx;
   double facty = fact/absMAXy;

    for ( size_t counter = 0; counter <= clip.size()-1; counter++) { 
     clip[0].push_back(IntPoint(long64(xFOR[counter]*factx+0.5),long64(yFOR[counter]*facty+0.5))); }//add +0.5 couse conversion from double to long long implies truncation, i.e. 29.7 becomes 29
   //add the Clip polygon to Clipper ... 
   Clipper clpr; 
   clpr.AddPolygons(clip, ptClip);  
}

 
extern "C" void wrapClipperRes(int oper[1],double POLYresX[13][100],double POLYresY[13][100],int VERTres[13], double absMAXx, double absMAXy) {
   //  oper = 1: intersection 2: union 3:
   //
   //
   // declarations
   int typeOPER ;
   typeOPER = oper[0];
   Polygons solution;
   Clipper clpr; 
   double factx = fact/absMAXx;
   double facty = fact/absMAXy;
//
    switch (typeOPER)
    {   case 1:  //get the intersection of the subject and clip polygons ... 
           clpr.Execute(ctIntersection, solution, pftNonZero , pftNonZero); 
           break;
        case 2:  //get the union of the subject and clip polygons ... 
           clpr.Execute(ctUnion, solution, pftNonZero, pftNonZero); 
           break;
        case 3:  //get the difference of the subject and clip polygons ... 
           clpr.Execute(ctDifference, solution, pftNonZero, pftNonZero); 
           break;
        case 4:
           clpr.Execute(ctXor , solution, pftNonZero, pftNonZero); 
           break;
        default:
           cout   << "Clipping operation not valid";
          //return();
         // exit(0);
    }
// COPY solution to arrays POLunion, VERTinters
    for ( size_t Np = 0; Np <= solution.size()-1; Np++) {
       VERTres[Np] = solution[Np].size();
       for ( size_t Nv = 0; Nv <= solution[Np].size()-1; Nv++) { // size_t intstead of itn otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
          POLYresX[Np][Nv] = solution[Np][Nv].X/factx;
          POLYresY[Np][Nv] = solution[Np][Nv].Y/facty;
       } 
    }

}
