#include "stdafx.h"
#include "clipper.hpp" 
#include <iostream>
#include <fstream>
using namespace ClipperLib;
using namespace std;

extern "C" void wrapUNI_intCEL(int Ndry[1],double xclip[5],double yclip[5], double POLYStoBEjoinedX[1][20000000],double POLYStoBEjoinedY[1][20000000],int NPOLYtoBEjoined[1], int VERTtoBEjoined[1],  double POLYSunionX[1][20000000],double POLYSunionY[1][20000000], int NPOLYSunion[1], int VERTunion[1], double POLYintersX[1][20000000],double POLYintersY[1][20000000], int NPOLYSinters[1], int VERTinters[1], double absMAXx[1], double absMAXy[1]) { 
   //  oper = 1: intersection 2: union 3:
   //
   //
   // declarations
   int NdryC = Ndry[0];
   int NPOLYtoBEjoinedC = NPOLYtoBEjoined[0];
   //int NPOLYSintersC = NPOLYSinters[0];
   Polygons subj(NPOLYtoBEjoinedC),solution,Punion,clip(1);
   Clipper clpr; 
   double factx = fact/absMAXx[0];
   double facty = fact/absMAXy[0];
  // factOK = 
   // COPY POLYStoBEjoinedX,POLYStoBEjoinedY to arrays POLYStoBEjoinedC, VERTinters
    for ( int Np = 0; Np < NPOLYtoBEjoinedC; Np++) {
       for ( int Nv = 0; Nv <VERTtoBEjoined[Np]; Nv++) { // size_t intstead of itn otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
           subj[Np].push_back(IntPoint(long64(POLYStoBEjoinedX[Np][Nv]*factx+0.5),long64(POLYStoBEjoinedY[Np][Nv]*facty+0.5)));  //+0.5 because when converting to integer the number is truncated, adding 0.5 make 1000.7 become 1001 instead of 1000
       } 
    }
   // COPY fortran clip to C++ polygon clip 
   for ( int Nv = 0; Nv < NdryC; Nv++) {
         clip[0].push_back(IntPoint(long64(xclip[Nv]*factx+0.5),long64(yclip[Nv]*facty+0.5)));  //[Nv].X = xclip[Nv];
         //clip[0][Nv].Y = yclip[Nv];
   }
   //clpr.Clear(); //clear all the polygons in clipper
   clpr.AddPolygons(subj, ptSubject);   //add polygons to be joined to clipper
   clpr.Execute(ctUnion, Punion, pftNonZero, pftNonZero);  //get the union of the subjects (already loaded in clipper) and clip polygons ... 
   clpr.Clear(); //clear all the polygons in clipper
   clpr.AddPolygons(Punion, ptSubject); //add the union to the subject
   clpr.AddPolygons(clip, ptClip);    //add the cell to the clipper
   clpr.Execute(ctIntersection, solution, pftNonZero , pftNonZero);   //get the intersection of the subjects and clip polygons ..

   //ofstream outputFILE; 
   //ofstream fileSUBJ; 
   //ofstream fileCLIP; 
   //ofstream filePunion; 
   //outputFILE.open ("solution.txt");
   //fileSUBJ.open ("SUBJ.txt");
   //fileCLIP.open ("CLIP.txt");
   //filePunion.open ("union.txt");
   //fileSUBJ.precision(20);
   //fileCLIP.precision(20);
   //outputFILE.precision(20);

   //// print  subject polygons on text file
   //fileSUBJ << subj.size()<<"\n"; //numb of polygons
   //for ( size_t Np = 0; Np < subj.size(); Np++) {
   //   fileSUBJ << subj[Np].size()<<"\n"; //numb of verices
   //   for ( size_t counter = 0; counter < subj[Np].size(); counter++)  // size_t intstead of itn otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
   //      fileSUBJ <<  double(subj[Np][counter].X)/factx << ' ' <<   double(subj[Np][counter].Y)/facty  <<"\n" ; 
   //}
   //fileSUBJ << flush;

   //// print  clip polygons on text file
   //fileCLIP << clip.size()<<"\n"; //numb of polygons
   //for ( size_t Np = 0; Np < clip.size(); Np++) {
   //   fileCLIP << clip[Np].size()<<"\n"; //numb of verices
   //   for ( size_t counter = 0; counter <clip[Np].size(); counter++)  // size_t intstead of itn otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
   //     fileCLIP <<  double(clip[Np][counter].X)/factx << ' ' <<   double(clip[Np][counter].Y)/facty  <<"\n" ; 
   //}
   //fileCLIP << flush;

   // // print  union polygons on text file
   //filePunion << Punion.size()<<"\n"; //numb of polygons
   //for ( size_t Np = 0; Np < Punion.size(); Np++) {
   //   filePunion << Punion[Np].size()<<"\n"; //numb of verices
   //   for ( size_t counter = 0; counter < Punion[Np].size(); counter++)  // size_t intstead of itn otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
   //     filePunion <<  double(Punion[Np][counter].X)/factx << ' ' <<   double(Punion[Np][counter].Y)/facty  <<"\n" ; 
   //}
   //filePunion << flush;

   //// print  solution polygons on text file
   //if (solution.size()>0) {
   //   outputFILE << solution.size()<<"\n"; //numb of polygons
   //   for ( size_t Np = 0; Np < solution.size(); Np++) {
   //      outputFILE << solution[Np].size()<<"\n"; //numb of verices
   //      for ( size_t counter = 0; counter <solution[Np].size(); counter++)  // size_t intstead of itn otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
   //      outputFILE <<  double(solution[Np][counter].X)/factx << ' ' <<   double(solution[Np][counter].Y)/facty  <<"\n" ; 
   //   }
   //   outputFILE << flush;
   //}
// COPY union to arrays POLYunion, VERTunion
   NPOLYSunion[0] = Punion.size();
   if (NPOLYSunion[0]>0) {      
      for ( size_t Np = 0; Np < Punion.size(); Np++) {
         VERTunion[Np] = Punion[Np].size();
         for ( size_t Nv = 0; Nv < Punion[Np].size(); Nv++) { // size_t intstead of itn otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
            POLYSunionX[Np][Nv] = double(Punion[Np][Nv].X)/factx;
            POLYSunionY[Np][Nv] = double(Punion[Np][Nv].Y)/facty;
         }
      }    
   } 

// COPY solution to arrays POLYinters, VERTinters
   NPOLYSinters[0] = solution.size();
   if (NPOLYSinters[0]>0) {      
      for ( size_t Np = 0; Np < solution.size(); Np++) {
         VERTinters[Np] = solution[Np].size();
         for ( size_t Nv = 0; Nv < solution[Np].size(); Nv++) { // size_t intstead of itn otherwise you have warning C4018: '<=' : signed/unsigned mismatch. (since it doesn't make sense for containers to contain a negative number of elements).
            POLYintersX[Np][Nv] = double(solution[Np][Nv].X)/factx;
            POLYintersY[Np][Nv] = double(solution[Np][Nv].Y)/facty;
         }
      }    
   } 
}




