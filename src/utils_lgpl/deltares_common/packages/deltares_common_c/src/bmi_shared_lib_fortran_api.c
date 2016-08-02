//---- LGPL --------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2016.
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
// $Id$
// $HeadURL$
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "so_fortran_api.h"

#ifndef min
#  define min(a,b) (a)<(b) ? (a) : (b)
#  define max(a,b) (a)>(b) ? (a) : (b)
#endif

#if defined(WIN32)
#  include <windows.h>
#elif defined(salford32)
#  include <windows.h>
#elif defined(HAVE_CONFIG_H)
#  include <dlfcn.h>
#endif

#if defined(WIN32)
#  define BMI_INITIALIZE        BMI_INITIALIZE
#  define BMI_UPDATE            BMI_UPDATE
#  define BMI_FINALIZE          BMI_FINALIZE
#  define BMI_GET_START_TIME    BMI_GET_START_TIME
#  define BMI_GET_END_TIME      BMI_GET_END_TIME
#  define BMI_GET_CURRENT_TIME  BMI_GET_CURRENT_TIME
#  define BMI_GET_VAR_TYPE      BMI_GET_VAR_TYPE
#  define BMI_GET_VAR_RANK      BMI_GET_VAR_RANK
#  define BMI_GET_VAR_SHAPE     BMI_GET_VAR_SHAPE
#  define BMI_SET_VAR           BMI_SET_VAR
#  define BMI_GET_VAR           BMI_GET_VAR
#  define STDCALL
#elif defined(HAVE_CONFIG_H)
#   include "config.h"
#  define BMI_INITIALIZE        FC_FUNC(initialize,BMI_INITIALIZE)
#  define BMI_UPDATE            FC_FUNC(update,BMI_UPDATE)
#  define BMI_FINALIZE          FC_FUNC(finalize,BMI_FINALIZE)
#  define BMI_GET_START_TIME    FC_FUNC(get_start_time,BMI_GET_START_TIME)
#  define BMI_GET_END_TIME      FC_FUNC(get_end_time,BMI_GET_END_TIME)
#  define BMI_GET_CURRENT_TIME  FC_FUNC(get_current_time,BMI_GET_CURRENT_TIME)
#  define BMI_GET_VAR_TYPE      FC_FUNC(get_var_type,BMI_GET_VAR_TYPE)
#  define BMI_GET_VAR_RANK      FC_FUNC(get_var_rank,BMI_GET_VAR_RANK)
#  define BMI_GET_VAR_SHAPE     FC_FUNC(get_var_shape,BMI_GET_VAR_SHAPE)
#  define BMI_SET_VAR           FC_FUNC(set_var,BMI_SET_VAR)
#  define BMI_GET_VAR           FC_FUNC(get_var,BMI_GET_VAR)
#  define STDCALL
#endif

/*
*
* Connection routine between F90 (main) -> C (interface) -> F90 (DLL).
* Special attention to the WINAPI define, which is needed if the DLL is written in F90
* Support methods are implemented in shared_library_fortran_api:
* .  strFcpy
* .  RemoveTrailingBlanks_dll
*
*/

#if defined(WIN32)
typedef HMODULE DllHandle;
#elif defined(HAVE_CONFIG_H)
typedef void * DllHandle;
#endif

typedef struct {
	DllHandle   dllHandle;
} SharedDLL;


long STDCALL BMI_INITIALIZE(long * sharedDLLHandle,
	char   * config_file,
	int      config_file_len)
{

	long error = -1;

	char * c_config_file = strFcpy(config_file, config_file_len);
	RemoveTrailingBlanks_dll(c_config_file);

	char * fun_name = "initialize";
	SharedDLL * sharedDLL = (SharedDLL *)(*sharedDLLHandle);

	typedef void * (STDCALL * MyProc)(char *);
	MyProc proc;

#if defined(WIN32)
	proc = (MyProc)GetProcAddress(sharedDLL->dllHandle, fun_name);
#elif defined(HAVE_CONFIG_H)
	proc = (MyProc)dlsym(sharedDLL->dllHandle, fun_name);
#endif

	if (proc != NULL)
	{
		error = 0;
		(void *)(*proc)(c_config_file);
	}
	free(c_config_file); c_config_file = NULL;

	return error;
}