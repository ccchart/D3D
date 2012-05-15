//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2012.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
// $Id$
// $HeadURL$
#if defined(WIN32) || defined (WIN64)
#  include <io.h>
#  include <wtypes.h>
#elif defined (salford32)
#  include <io.h>
#  include <windows.h>
#endif
#include <string.h>

#if HAVE_CONFIG_H
#   include "config.h"
#endif
#include "version_number.h"
#define CAT(a, b) a ## b
#define FUNC_CAT(a, b) CAT(a, b)

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

extern  char * FUNC_CAT( version_getFileVersionString_, MOD_NAME)(void);
extern  char * FUNC_CAT( version_getFullVersionString_, MOD_NAME)(void);
extern  char * FUNC_CAT( version_getCompanyString_, MOD_NAME)(void);
extern  char * FUNC_CAT( version_getVersionNumberString_, MOD_NAME)(void);
extern  char * FUNC_CAT( version_getProgramNameString_, MOD_NAME)(void);
extern  char * FUNC_CAT( version_getShortProgramNameString_, MOD_NAME)(void);
extern  char * FUNC_CAT( version_getSvnRevisionString_, MOD_NAME)(void);
extern  char * FUNC_CAT( version_getFeatureNumberString_, MOD_NAME)(void);


/*
 * FTN_CAPITAL   : dvf6, salford
 * FTN_UNDERSCORE: sgi, sun, cygwin, linux
 * FTN_SMALL     : hp
 */

#if HAVE_CONFIG_H
#   define FTN_CALL  /* nothing */
#   define F90_GETFILEVERSIONSTRING      FC_FUNC( FUNC_CAT(getfileversionstring,F90_MOD_NAME), FUNC_CAT(GETFILEVERSIONSTRING,F90_MOD_NAME))
#   define F90_GETFULLVERSIONSTRING      FC_FUNC( FUNC_CAT(getfullversionstring,F90_MOD_NAME), FUNC_CAT(GETFULLVERSIONSTRING,F90_MOD_NAME))
#   define F90_GETCOMPANYSTRING          FC_FUNC( FUNC_CAT(getcompanystring,F90_MOD_NAME), FUNC_CAT(GETCOMPANYSTRING,F90_MOD_NAME))
#   define F90_GETVERSIONNUMBERSTRING    FC_FUNC( FUNC_CAT(getversionnumberstring,F90_MOD_NAME), FUNC_CAT(GETVERSIONNUMBERSTRING,F90_MOD_NAME))
#   define F90_GETPROGRAMNAMESTRING      FC_FUNC( FUNC_CAT(getprogramnamestring,F90_MOD_NAME), FUNC_CAT(GETPROGRAMNAMESTRING,F90_MOD_NAME))
#   define F90_GETSHORTPROGRAMNAMESTRING FC_FUNC( FUNC_CAT(getshortprogramnamestring,F90_MOD_NAME), FUNC_CAT(GETSHORTPROGRAMNAMESTRING,F90_MOD_NAME))
#   define F90_GETSVNREVISIONSTRING      FC_FUNC( FUNC_CAT(getsvnrevisionstring,F90_MOD_NAME), FUNC_CAT(GETSVNREVISIONSTRING,F90_MOD_NAME))
#   define F90_GETFEATURENUMBERSTRING    FC_FUNC( FUNC_CAT(getfeaturenumberstring,F90_MOD_NAME), FUNC_CAT(GETFEATURENUMBERSTRING,F90_MOD_NAME))
#else
/* WIN32 or WIN64 */
#   define FTN_CALL
#   define F90_GETFILEVERSIONSTRING      FUNC_CAT( GETFILEVERSIONSTRING_, MOD_NAME)
#   define F90_GETFULLVERSIONSTRING      FUNC_CAT( GETFULLVERSIONSTRING_, MOD_NAME)
#   define F90_GETCOMPANYSTRING          FUNC_CAT( GETCOMPANYSTRING_, MOD_NAME)
#   define F90_GETVERSIONNUMBERSTRING    FUNC_CAT( GETVERSIONNUMBERSTRING_, MOD_NAME)
#   define F90_GETPROGRAMNAMESTRING      FUNC_CAT( GETPROGRAMNAMESTRING_, MOD_NAME)
#   define F90_GETSHORTPROGRAMNAMESTRING FUNC_CAT( GETSHORTPROGRAMNAMESTRING_, MOD_NAME)
#   define F90_GETSVNREVISIONSTRING      FUNC_CAT( GETSVNREVISIONSTRING_, MOD_NAME)
#   define F90_GETFEATURENUMBERSTRING    FUNC_CAT( GETFEATURENUMBERSTRING_, MOD_NAME)
#endif

/*==========================================================================*/
 void FTN_CALL F90_GETFILEVERSIONSTRING( char * str, int length_str )
{
  char * string;
  int i;
  for (i=0; i<length_str; i++) {str[i] = ' ';}
  string = FUNC_CAT( version_getFileVersionString_, MOD_NAME)();
  i  = min((int) length_str, (int) strlen(string));
  strncpy(str, string, i);
}

/*==========================================================================*/
 void FTN_CALL F90_GETFULLVERSIONSTRING( char * str, int length_str )
{
   char * string;
  int i;
  for (i=0; i<length_str; i++) {str[i] = ' ';}
  string = FUNC_CAT( version_getFullVersionString_, MOD_NAME)();
  i  = min((int) length_str, (int) strlen(string));
  strncpy(str, string, i);
}

/*==========================================================================*/
 void FTN_CALL F90_GETCOMPANYSTRING( char * str, int length_str )
{
  char * string;
  int i;
  for (i=0; i<length_str; i++) {str[i] = ' ';}
  string = FUNC_CAT( version_getCompanyString_, MOD_NAME)();
  i  = min(length_str, (int) strlen(string));
  strncpy(str, string, i);
}

/*==========================================================================*/
 void FTN_CALL F90_GETVERSIONNUMBERSTRING( char * str, int length_str )
{
  char * string;
  int i;
  for (i=0; i<length_str; i++) {str[i] = ' ';}
  string = FUNC_CAT( version_getVersionNumberString_, MOD_NAME)();
  i  = min(length_str, (int) strlen(string));
  strncpy(str, string, i);
}

/*==========================================================================*/
 void FTN_CALL F90_GETPROGRAMNAMESTRING( char * str, int length_str )
{
  char * string;
  int i;
  for (i=0; i<length_str; i++) {str[i] = ' ';}
  string = FUNC_CAT( version_getProgramNameString_, MOD_NAME)();
  i  = min(length_str, (int) strlen(string));
  strncpy(str, string, i);
}

/*==========================================================================*/
 void FTN_CALL F90_GETSHORTPROGRAMNAMESTRING( char * str, int length_str )
{
  char * string;
  int i;
  for (i=0; i<length_str; i++) {str[i] = ' ';}
  string = FUNC_CAT( version_getShortProgramNameString_, MOD_NAME)();
  i  = min(length_str, (int) strlen(string));
  strncpy(str, string, i);
}

/*==========================================================================*/
 void FTN_CALL F90_GETSVNREVISIONSTRING( char * str, int length_str )
{
  char * string;
  int i;
  for (i=0; i<length_str; i++) {str[i] = ' ';}
  string = FUNC_CAT( version_getSvnRevisionString_, MOD_NAME)();
  i  = min(length_str, (int) strlen(string));
  strncpy(str, string, i);
}
/*==========================================================================*/
 void FTN_CALL F90_GETFEATURENUMBERSTRING( char * str, int length_str )
{
  char * string;
  int i;
  for (i=0; i<length_str; i++) {str[i] = ' ';}
  string = FUNC_CAT( version_getFeatureNumberString_, MOD_NAME)();
  i  = min(length_str, (int) strlen(string));
  strncpy(str, string, i);
}
