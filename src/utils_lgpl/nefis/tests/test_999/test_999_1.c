//---- LGPL --------------------------------------------------------------------
//                                                                              
// Copyright (C)  Stichting Deltares, 2011.                                     
//                                                                              
// This library is free software; you can redistribute it and/or                
// modify it under the terms of the GNU Lesser General Public                   
// License as published by the Free Software Foundation version 2.1.                 
//                                                                              
// This library is distributed in the hope that it will be useful,              
// but WITHOUT ANY WARRANTY; without even the implied warranty of               
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
// Lesser General Public License for more details.                              
//                                                                              
// You should have received a copy of the GNU Lesser General Public             
// License along with this library; if not, see <http://www.gnu.org/licenses/>. 
//                                                                              
// contact: delft3d.support@deltares.nl                                         
// Stichting Deltares                                                           
// P.O. Box 177                                                                 
// 2600 MH Delft, The Netherlands                                               
//                                                                              
// All indications and logos of, and references to, "Delft3D" and "Deltares"    
// are registered trademarks of Stichting Deltares, and remain the property of  
// Stichting Deltares. All rights reserved.                                     
//                                                                              
//------------------------------------------------------------------------------
#include <stdio.h>
#include <string.h>
#include "nefis.h"
int main()
{
  long fds;
  char coding = ' ';
  long error, i, j, k, l, m, n;
  char *  cel_name;
  char *  grp_name;
  long grpdms[5], grpord[5], usrord[5], uindex[5][3];
  float buffer[28660];
  char errstr[134];
  float elap_w, cpu2 = 0., cpu1 = 0.;
  char names[100][16+1];

  fds = -1;

  error = Opndef( &fds, "nefis_ex.def", coding);
  if (error!=0) {error = Neferr( 1, errstr);}

  error = Opndat( &fds, "nefis_ex.dat", coding);
  if (error!=0) {error = Neferr( 1, errstr);}

  strcpy(names[0], "ELEM_R_4_DIM_1");
  strcpy(names[1], "ELEM_R_4");
  cel_name = strdup("CEL_TEST_3");
  error= Defcel( &fds, cel_name, 2, names);
  if (error!=0) {error = Neferr( 1, errstr);}

  grpdms[0]= 4;
  grpdms[1]= 5;
  grpdms[2]= 6;
  grpdms[3]= 7;
  grpdms[4]= 8;

  grpord[0]= 1;
  grpord[1]= 2;
  grpord[2]= 3;
  grpord[3]= 4;
  grpord[4]= 5;

  grp_name = strdup("GRP_TEST_3A");
  error= Defgrp( &fds, grp_name, cel_name, 5, grpdms, grpord);
  if (error!=0) {error = Neferr( 1, errstr);}

  grpord[0]= 5;
  grpord[1]= 4;
  grpord[2]= 3;
  grpord[3]= 2;
  grpord[4]= 1;
  grp_name = strdup("GRP_TEST_3B");
  error= Defgrp( &fds, grp_name, cel_name, 5, grpdms, grpord);
  if (error!=0) {error = Neferr( 1, errstr);}

  grpord[0]= 5;
  grpord[1]= 3;
  grpord[2]= 1;
  grpord[3]= 2;
  grpord[4]= 4;
  grp_name = strdup("GRP_TEST_3C");
  error= Defgrp( &fds, grp_name, cel_name, 5, grpdms, grpord);
  if (error!=0) {error = Neferr( 1, errstr);}

  grpdms[0]= 0;

  grpord[0]= 1;
  grpord[1]= 2;
  grpord[2]= 3;
  grpord[3]= 4;
  grpord[4]= 5;

  grp_name = strdup("GRP_TEST_3D");
  error= Defgrp( &fds, grp_name, cel_name, 5, grpdms, grpord);
  if (error!=0) {error = Neferr( 1, errstr);}
/*--------------------------------------------------------------------- */
  error= Credat( &fds, "DATAGRP_TEST_3A", "GRP_TEST_3A");
  if (error!=0) {error = Neferr( 1, errstr);}

  error= Credat( &fds, "DATAGRP_TEST_3B", "GRP_TEST_3B");
  if (error!=0) {error = Neferr( 1, errstr);}

  error= Credat( &fds, "DATAGRP_TEST_3C", "GRP_TEST_3C");
  if (error!=0) {error = Neferr( 1, errstr);}

  error= Credat( &fds, "DATAGRP_TEST_3D", "GRP_TEST_3D");
  if (error!=0) {error = Neferr( 1, errstr);}
/*--------------------------------------------------------------------- */
      fprintf(stdout,"Initialisation NEFIS files [sec] %g\n",cpu2-cpu1);

      uindex[0][2] = 1;
      uindex[1][2] = 1;
      uindex[2][2] = 1;
      uindex[3][2] = 1;
      uindex[4][2] = 1;
      usrord[0] = 1;
      usrord[1] = 2;
      usrord[2] = 3;
      usrord[3] = 4;
      usrord[4] = 5;

      fprintf(stdout,"\n");
      fprintf(stdout,"6720 schrijfopdrachten (6720 cellen) van 16 bytes elk\n");
      for ( i=1; i<=8; i++) {
          uindex[4][0] = i;
          uindex[4][1] = i;
          for ( j=1; j<=7; j++) {
             uindex[3][0] = j;
             uindex[3][1] = j;
             for ( k=1; k<=6; k++) {
                uindex[2][0] = k;
                uindex[2][1] = k;
                for ( l=1; l<=5; l++) {
                    uindex[1][0] = l;
                    uindex[1][1] = l;
                    for ( m=1; m<=4; m++) {
                        uindex[0][0] = m;
                        uindex[0][1] = m;
                        for ( n=0; n<4; n++) {
                            buffer[n]= (float) (i*j*k*l*m*(n+1));
                        }

                        error= Putelt( &fds, "DATAGRP_TEST_3A", "*", uindex, usrord, buffer);
                        if (error!=0) {error = Neferr( 1, errstr);}
      }}}}}
      elap_w = cpu2-cpu1;
      fprintf(stdout,"DATAGRP_TEST_3A written in [sec] %#13.5g\n", cpu2-cpu1);

	  error = Clsnef( &fds );
      error = Crenef( &fds, "nefis_ex.dat", "nefis_ex.def", ' ', 'u');

      fprintf(stdout,"\n");
      fprintf(stdout,"6720 schrijfopdrachten (6720 cellen) van 16 bytes elk\n");
      for ( i=1; i<=8; i++) {
          uindex[4][0] = i;
          uindex[4][1] = i;
          for ( j=1; j<=7; j++) {
             uindex[3][0] = j;
             uindex[3][1] = j;
             for ( k=1; k<=6; k++) {
                uindex[2][0] = k;
                uindex[2][1] = k;
                for ( l=1; l<=5; l++) {
                    uindex[1][0] = l;
                    uindex[1][1] = l;
                    for ( m=1; m<=4; m++) {
                        uindex[0][0] = m;
                        uindex[0][1] = m;
                        for ( n=0; n<4; n++) {
                            buffer[n]= (float) (i*j*k*l*m*(n+1)*2);
                        }

                        error= Putelt( &fds, "DATAGRP_TEST_3D", "*", uindex, usrord, buffer);
                        if (error!=0) {error = Neferr( 1, errstr);}
      }}}}}
      elap_w = elap_w + cpu2-cpu1;
      fprintf(stdout,"DATAGRP_TEST_3D written in [sec] %#13.5g\n", cpu2-cpu1);

	  fprintf(stdout,"Elapsed time %#13.5g\n", elap_w);

	  error = Clsnef( &fds );


  return(0);
}
