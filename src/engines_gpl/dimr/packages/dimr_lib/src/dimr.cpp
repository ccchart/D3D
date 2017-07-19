//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2017.
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
// $Id: dimr.cpp 962 2011-10-31 21:52:47Z elshoff $
// $HeadURL: $
//------------------------------------------------------------------------------
//  dimr Main Program
//
//  Irv.Elshoff@Deltares.NL
//  6 mar 13
//------------------------------------------------------------------------------

#define DIMR_LIB


#include "dimr.h"
#include "dimr_lib_version.h"

#if defined(HAVE_CONFIG_H)
#include "config.h"
#include <dlfcn.h>
#include <libgen.h>
#endif
#include <expat.h>
#include <limits.h>
#include <time.h>

#if defined (MEMCHECK)
#include <mcheck.h>
#endif


#include <typeinfo>
#include "dimr_bmi_utils.h"
using namespace std;

#include <string>
#include <sstream>
#include <math.h>
#include <mpi.h>

#if defined (WIN32)
//#  include "getopt.h"
#  include <Strsafe.h>
#  include <windows.h>
#  include <direct.h>
#include <errno.h>
#  include <io.h>
#  include <sys/stat.h>
#  define strdup _strdup
#  define chdir _chdir
#  define getcwd _getcwd
#  define dup2 _dup2
#else
#  include <unistd.h>
#endif

// static added to prevent name conflicts on Linux.
static Dimr * thisDimr = NULL;     // global pointer to single object instance

Dimr::Dimr(void) {
    FILE * logFile               = stdout;
    ready                        = false;
    exePath                      = NULL;
    exeName                      = NULL;
    clock                        = new Clock ();
    logMask                      = Log::MAJOR;
    feedbackMask                 = Log::MAJOR;
    log                          = new Log (logFile, clock, logMask);
    config                       = NULL;
    mainArgs                     = NULL;
    slaveArg                     = NULL;
    control                      = NULL;
    componentsList.numComponents = 0;
    couplersList.numCouplers     = 0;
    use_mpi                      = false;
    my_rank                      = 0;
    numranks                     = 1;
    configfile                   = NULL;
    done                         = false;
    // Initialize redirectFile. Default: switched On (!=NULL)
    const char * filename        = "dimr_redirected.log";
    int len                      = strlen(filename);
    redirectFile                 = (char *) malloc((len+1)*sizeof(char));
    strncpy(redirectFile, (const char*)filename, len);
    redirectFile[len]            = '\0';
}


extern "C" {
//------------------------------------------------------------------------------
DllExport void set_logger_callback(WriteCallback writeCallBack) {
    if (thisDimr == NULL) {
        thisDimr = new Dimr();
    }
    thisDimr->log->SetWriteCallBack(writeCallBack);
}
    
//------------------------------------------------------------------------------
DllExport void set_dimr_logger(Log * loggerFromDimrExe) {
    if (thisDimr == NULL) {
        thisDimr = new Dimr();
    }
    thisDimr->log = loggerFromDimrExe;
}

BMI_API void set_logger(Logger logger){
	if (thisDimr == NULL) {
		thisDimr = new Dimr();
	}
	
	thisDimr->log->SetExternalLogger(logger);
	//Update or add the extertnal logger function in the kernels
	for (int i = 0; i < thisDimr->componentsList.numComponents; i++) {
		if (thisDimr->componentsList.components[i].type == COMP_TYPE_FLOW1D)
		{
			if (thisDimr->componentsList.components[i].setLogger != NULL) {
				thisDimr->componentsList.components[i].setLogger(logger);
			}
			double level = (double)DimrBmiUtils::convertDimrLogLevelToLogLevel(thisDimr->logMask);
			thisDimr->componentsList.components[i].dllSetVar("debugLevel", (const void *)&level);
		}
	}
}
//------------------------------------------------------------------------------
BMI_API int initialize(const char * configfile) {
    int nSettingsSet, nParamsSet;
    if (thisDimr == NULL) {
        thisDimr = new Dimr();
    }
    if (thisDimr->redirectFile != NULL) {
        // RedirectFile must be including the full path:
        // - Get the basename (platform dependent implementation)
        // - if (redirectfile == basename) then
        //       Make copy of redirectfile
        //       Put CWD in redirectfile
        //       redirectfile = redirectfile + / + copy
        char *fileBasename = new char[MAXSTRING];
#if defined(HAVE_CONFIG_H)
        fileBasename = strdup (basename (thisDimr->redirectFile));
        const char *dirSeparator = "/";
#else
        char * ext = new char[5];
        _splitpath (thisDimr->redirectFile, NULL, NULL, fileBasename, ext);
        StringCbCatA (fileBasename, MAXSTRING, ext);
        delete [] ext;
        const char *dirSeparator = "\\";
#endif
        if (strcmp(thisDimr->redirectFile,fileBasename) == 0) {
            char *filenameCopy = new char[MAXSTRING];
            strcpy(filenameCopy, thisDimr->redirectFile);

            delete [] thisDimr->redirectFile;
            thisDimr->redirectFile                 = (char *) malloc((MAXSTRING)*sizeof(char));
            if (!getcwd(thisDimr->redirectFile, MAXSTRING)) {
                throw new Exception (true, "ERROR obtaining the current working directory (init)");
            }
            strcat(thisDimr->redirectFile,dirSeparator);
            strcat(thisDimr->redirectFile,filenameCopy);
            delete [] filenameCopy;
        }
        // Redirection to file is currently handled in the logger by writing directly to the specified file
        thisDimr->log->redirectFile = thisDimr->redirectFile;
        printf("DIMR messages are redirected to file \"%s\"\n", thisDimr->redirectFile);
        fflush(stdout);
        // Create an empty file
        FILE * fp = fopen(thisDimr->redirectFile,"w+");
        fclose(fp);
        delete [] fileBasename;
    }

    thisDimr->log->Write(Log::MAJOR, thisDimr->my_rank, getfullversionstring_dimr_lib());
    thisDimr->log->Write(Log::MAJOR, thisDimr->my_rank, "dimr_dll:initialize(%s)", configfile);
    //
    //
    // Read XML configuration file into tree structure
    thisDimr->configfile = configfile;
    FILE * conf;
    if (strcmp (thisDimr->configfile, "-") == 0)
        conf = stdin;
    else {
        conf = fopen (thisDimr->configfile, "r");
        if (conf == NULL){
            thisDimr->log->Write(Log::MAJOR, thisDimr->my_rank, "Cannot open configuration file \"%s\"", thisDimr->configfile);
            throw new Exception(true, "Cannot open configuration file \"%s\"", thisDimr->configfile);
        }
    }

    thisDimr->config = new XmlTree (conf);
    thisDimr->config->ExpandEnvironmentVariables();
    fclose (conf);
    //
    // Build controlBlock administration by scanning the XmlTree
    thisDimr->log->Write(Log::MAJOR, thisDimr->my_rank, "Build controlBlock administration by scanning the XmlTree");
    thisDimr->scanConfigFile();
    //
    // ToDo: check whether a core dump is requested on abort; if so set global variable for Dimr_CoreDump
    //
    // This is a good time to attach to the processes in case you want to debug
    thisDimr->processWaitFile();
    //
    // Build connection with dlls
    thisDimr->connectLibs();
    // Init the timers before calling the dllInitialize routines!
    thisDimr->timersInit();
    //
    // Initialize the components in the first controlBlock only
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        thisDimr->runParallelInit(&(thisDimr->control->subBlocks[0]));
    } else {
        // Start block

        // Hack for WAVE:
        if (thisDimr->control->subBlocks[0].unit.component->type == COMP_TYPE_WAVE) {
            int *waveModePtr        = NULL;
            const char *key = "mode";
            (thisDimr->control->subBlocks[0].unit.component->dllGetVar) (key, (void**)(&waveModePtr));
            *waveModePtr = 0;
        }

        chdir(thisDimr->control->subBlocks[0].unit.component->workingDir);
        thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "%s.Initialize(%s)", thisDimr->control->subBlocks[0].unit.component->name, thisDimr->control->subBlocks[0].unit.component->inputFile);
        nSettingsSet = thisDimr->control->subBlocks[0].unit.component->dllSetKeyVals(thisDimr->control->subBlocks[0].unit.component->settings);
        thisDimr->timerStart(thisDimr->control->subBlocks[0].unit.component);
        thisDimr->control->subBlocks[0].unit.component->result = (thisDimr->control->subBlocks[0].unit.component->dllInitialize) (thisDimr->control->subBlocks[0].unit.component->inputFile);
        thisDimr->timerEnd(thisDimr->control->subBlocks[0].unit.component);
        nParamsSet = thisDimr->control->subBlocks[0].unit.component->dllSetKeyVals(thisDimr->control->subBlocks[0].unit.component->parameters );
        (thisDimr->control->subBlocks[0].unit.component->dllGetStartTime) (&thisDimr->control->subBlocks[0].tStart);
        (thisDimr->control->subBlocks[0].unit.component->dllGetEndTime) (&thisDimr->control->subBlocks[0].tEnd);
        (thisDimr->control->subBlocks[0].unit.component->dllGetTimeStep) (&thisDimr->control->subBlocks[0].tStep);
        (thisDimr->control->subBlocks[0].unit.component->dllGetCurrentTime) (&thisDimr->control->subBlocks[0].tCur);
    }
    // all ok (no exceptions)
    return 0;
}

//------------------------------------------------------------------------------
BMI_API int update(double tStep) {
    thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr_lib:update");
    // Execute update on the first controlBlock only
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        thisDimr->runControlBlock(&(thisDimr->control->subBlocks[0]), tStep, GLOBAL_PHASE_UPDATE);
    } else {
        // Start block
        chdir(thisDimr->control->subBlocks[0].unit.component->workingDir);
        thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "%s.Update(%6.1f)", thisDimr->control->subBlocks[0].unit.component->name, tStep);
        thisDimr->timerStart(thisDimr->control->subBlocks[0].unit.component);
        (thisDimr->control->subBlocks[0].unit.component->dllUpdate) (tStep);
        thisDimr->timerEnd(thisDimr->control->subBlocks[0].unit.component);
        (thisDimr->control->subBlocks[0].unit.component->dllGetCurrentTime) (&thisDimr->control->subBlocks[0].tCur);
    }
	return 0;
}

//------------------------------------------------------------------------------
BMI_API int finalize(void) {
    thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr_lib:finalize");
    // Execute finalize on the first controlBlock and
    // initialize, step, finalize on all other controlBlocks
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        thisDimr->runParallelFinish(&(thisDimr->control->subBlocks[0]));
    } else {
        // Start block
        chdir(thisDimr->control->subBlocks[0].unit.component->workingDir);
        thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "%s.Finalize()", thisDimr->control->subBlocks[0].unit.component->name);
        thisDimr->timerStart(thisDimr->control->subBlocks[0].unit.component);
        (thisDimr->control->subBlocks[0].unit.component->dllFinalize) ();
        thisDimr->timerEnd(thisDimr->control->subBlocks[0].unit.component);
        fflush(stdout);
    }
    thisDimr->timersFinish();
    for (int i = 1 ; i < thisDimr->control->numSubBlocks ; i++) {
        thisDimr->runControlBlock(&(thisDimr->control->subBlocks[i]),999999999.0, GLOBAL_PHASE_FINISH);
    }

    if (thisDimr->redirectFile != NULL) {
        thisDimr->log->redirectFile = NULL;
        printf("Finished: redirecting DIMR messages to file \"%s\"", thisDimr->redirectFile);
        fflush(stdout);
    }
	return 0;
}

//------------------------------------------------------------------------------
BMI_API void get_start_time(double * tStart) {
    thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr_lib:get_start_time");
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        *tStart = thisDimr->control->subBlocks[0].subBlocks[thisDimr->control->subBlocks[0].masterSubBlockId].tStart;
    } else {
        // Start block
        *tStart = thisDimr->control->subBlocks[0].tStart;
    }
}

//------------------------------------------------------------------------------
BMI_API void get_end_time(double * tEnd) {
    thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr_lib:get_end_time");
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        *tEnd = thisDimr->control->subBlocks[0].subBlocks[thisDimr->control->subBlocks[0].masterSubBlockId].tEnd;
    } else {
        // Start block
        *tEnd = thisDimr->control->subBlocks[0].tEnd;
    }
}

//------------------------------------------------------------------------------
BMI_API void get_time_step(double * tStep) {
    thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr_lib:get_time_step");
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        *tStep = thisDimr->control->subBlocks[0].subBlocks[thisDimr->control->subBlocks[0].masterSubBlockId].tStep;
    } else {
        // Start block
        *tStep = thisDimr->control->subBlocks[0].tStep;
    }
}

//------------------------------------------------------------------------------
BMI_API void get_current_time(double * tCur) {
    thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr_lib:get_current_time");
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        *tCur = thisDimr->control->subBlocks[0].subBlocks[thisDimr->control->subBlocks[0].masterSubBlockId].tCur;
    } else {
        // Start block
        *tCur = thisDimr->control->subBlocks[0].tCur;
    }
}

//------------------------------------------------------------------------------
BMI_API void get_var(const char * key, void ** ref) {
    char           * componentName = new char[thisDimr->MAXSTRING];
    const char     * slash         = strstr(key, "/");
    const char     * sourceName;
    dimr_component * compPtr       = NULL;
    double         * sourceVarPtr  = NULL; // In a coupler, when the source component is get/setting by ref, the pointer to
                                           // the resulting parameter is stored in sourceVarPtr
                                           // In this get_var, sourceVarPtr is always undefined.
                                           // NULL flags that it has to be retrieved.
    int              sourceProcess  = 0;   // With multiple possible source processes, sourceProcess flags what process is actualy delivering the value
                                           // Only relevant for parallel calculations. Possibly not working yet?

    thisDimr->log->Write (Log::MINOR, thisDimr->my_rank, "dimr_lib:get_var");
    // Assumption: "key" has the structure "componentName/group/id/parameter"
    if (slash == NULL) {
        // No component name specified in "key"
        *ref == NULL;
        delete[] componentName;
        return;
    }
    // componentName is everything before the first / in key
    strncpy(componentName, key, slash-key);
    componentName[slash-key] = '\0';
    // sourceName is everything behind the first / in key
    sourceName = slash + 1;
    if (strlen(sourceName) < 1) {
        throw new Exception(true, "dimr::get_var: No parameter specified. Expecting \"componentName/parameterName\"\n");
    }
    // Search componentName in the list of components of thisDimr
    for (int i = 0 ; i < thisDimr->componentsList.numComponents ; i++) {
        if (strcmp(thisDimr->componentsList.components[i].name, componentName) == 0) {
            compPtr = &thisDimr->componentsList.components[i];
            break;
        }
    }
    if (compPtr == NULL) {
        throw new Exception(true, "dimr::get_var: Unrecognized component \"%s\". Expecting \"componentName/parameterName\"\n", componentName);
    }
    // Get the pointer to the variable being asked for and put it in argument "ref"
    double * transfer = new double [compPtr->numProcesses];

    thisDimr->getAddress(sourceName, compPtr->type, compPtr->dllGetVar, &sourceVarPtr, compPtr->processes, compPtr->numProcesses, transfer);
    *ref = thisDimr->send(sourceName, compPtr->type, sourceVarPtr, compPtr->processes, compPtr->numProcesses, transfer);
    
    delete[] transfer;
    delete[] componentName;
}

//------------------------------------------------------------------------------
BMI_API void set_var(const char * key, const void * value) {
    char           * componentName = new char[thisDimr->MAXSTRING];
    const char     * slash = strstr(key, "/");
    const char     * targetName;
    dimr_component * compPtr = NULL;
    double         * targetVarPtr = NULL; // In a coupler, when the target component is get/setting by ref, the pointer to
                                          // the resulting parameter is stored in sourceVarPtr
                                          // In this set_var, targetVarPtr is always undefined.
                                          // NULL flags that it has to be retrieved via dllSetVar.
    int              sourceProcess = 0;   // With multiple possible source processes, sourceProcess flags what process is actualy delivering the value
                                          // Only relevant for parallel calculations. Possibly not working yet?

    // thisDimr->log is not initialized when set_var is called before initialize
    if (thisDimr == NULL) {
        thisDimr = new Dimr();
    }
    thisDimr->log->Write (Log::MINOR, thisDimr->my_rank, "dimr_lib:set_var");
    // Catch special keywords for Dimr_dll itself
    if (strcmp(key, "useMPI") == 0) {
        thisDimr->use_mpi = *(bool *)value;
    } else if (strcmp(key, "numRanks") == 0) {
        thisDimr->numranks = *(int *)value;
    } else if (strcmp(key, "myRank") == 0) {
        thisDimr->my_rank = *(int *)value;
    } else if (strcmp(key, "debugLevel") == 0) {
        thisDimr->logMask = *(Log::Mask *)value;
        thisDimr->log->SetMask(thisDimr->logMask);
    } else if (strcmp(key, "feedbackLevel") == 0) {
        thisDimr->feedbackMask = *(Log::Mask *)value;
        thisDimr->log->SetFeedbackLevel(thisDimr->feedbackMask);
    } else if (strcmp(key, "redirectFile") == 0) {
        // value is a char*
        // Special value: "stdout/stderr" => switch off redirection to file by setting to NULL
        // Else: value is the name of the file to redirect to
        if (strcmp((const char*)value,"stdout/stderr") == 0 && thisDimr->redirectFile != NULL) {
            free(thisDimr->redirectFile);
            thisDimr->redirectFile = NULL;
        } else {
            if (thisDimr->redirectFile != NULL) {
                free(thisDimr->redirectFile);
            }
            int len = strlen((const char*)value);
            thisDimr->redirectFile = (char *) malloc((len+1)*sizeof(char));
            strncpy(thisDimr->redirectFile, (const char*)value, len);
            thisDimr->redirectFile[len] = '\0';
        }
    } else {
        // Assumption: "key" has the structure "componentName/group/id/parameter"
        if (slash == NULL) {
            // No component name specified in "key"
            throw new Exception(true, "dimr::set_var: Unrecognized keyword \"%s\"\n", key);
        }
        // componentName is everything before the first / in key
        strncpy(componentName, key, slash-key);
        componentName[slash-key] = '\0';
        // targetName is everything behind the first / in key
        targetName = slash + 1;
        if (strlen(targetName) < 1) {
            throw new Exception(true, "dimr::set_var: No parameter specified. Expecting \"componentName/parameterName\"\n");
        }
        // Search componentName in the list of components of thisDimr
        for (int i = 0 ; i < thisDimr->componentsList.numComponents ; i++) {
            if (strcmp(thisDimr->componentsList.components[i].name, componentName) == 0) {
                compPtr = &thisDimr->componentsList.components[i];
                break;
            }
        }
        if (compPtr == NULL) {
            throw new Exception(true, "dimr::set_var: Unrecognized component \"%s\". Expecting \"componentName/parameterName\"\n", componentName);
        }
        // Send value to the receiving component
        thisDimr->receive (targetName, 
                               compPtr->type,
                               compPtr->dllSetVar,
                               compPtr->dllGetVar,
                               targetVarPtr,
                               compPtr->processes,
                               compPtr->numProcesses,
                               -1,
                               value);
    }
    delete[] componentName;
}


} // extern "C"


//------------------------------------------------------------------------------
int dimr_component::dllSetKeyVals (keyValueLL * kv) {
        // Pass parameters for the first controll block's component: parameters
        int count = 0;
        while(kv){
            if (dllSetVar!=NULL){
               (dllSetVar) (kv->key,(void*)kv->val);
            } else {
               if (dllGetVar!=NULL){
                  (dllGetVar) (kv->key,(void**)kv->val);
               }
            }
           kv = kv->nextkv;
           count++;
        }
        return count;
}


//
//------------------------------------------------------------------------------
//  Destructor
Dimr::~Dimr (void) {

    if (done)
        return;

    // to do:  (void) FreeLibrary(handle);
    freeLibs();

    log->Write (Log::MAJOR, thisDimr->my_rank, "dimr shutting down normally");

#if defined(HAVE_CONFIG_H)
    free (exeName);
#else
    delete [] exeName;
#endif

    delete clock;
    delete config;
    free (exePath);
    delete log;
    delete [] mainArgs;
    // componentsList
    delete [] componentsList.components;
    // couplersList
    for (int i = 0 ; i < couplersList.numCouplers ; i++) {
        delete [] couplersList.couplers[i].items;
    }
    delete [] couplersList.couplers;
    // control
    deleteControlBlock(*(control));
    free(control);
    done = true;
}



//------------------------------------------------------------------------------
void Dimr::deleteControlBlock (dimr_control_block cb) {
    if (cb.numSubBlocks > 0) {
        for (int i = 0 ; i < cb.numSubBlocks ; i++) {
            // Recursively delete all subBlocks
            deleteControlBlock(cb.subBlocks[i]);
        }
        delete [] cb.subBlocks;
    }
}




// WARNING: dimr is not BMI compliant yet!
// tStep is not used in runControlBlock
//------------------------------------------------------------------------------
void Dimr::runControlBlock (dimr_control_block * cb, double tStep, int phase) {
    if (cb->type == CT_PARALLEL) {
        log->Write (Log::MAJOR, my_rank, "PARALLEL:");
        //
        //
        // Initialize loop
        if (phase == GLOBAL_PHASE_FINISH) {
            runParallelInit(cb);
        }
        //
        //
        // Update loop
        runParallelUpdate(cb, tStep);
        //
        //
        // Finalize loop
        if (phase == GLOBAL_PHASE_FINISH) {
            runParallelFinish(cb);
        }
    } else if (cb->type == CT_START) {
        runStartBlock(cb, tStep, phase);
    }
    fflush(stdout);
}



//------------------------------------------------------------------------------
void Dimr::runStartBlock (dimr_control_block * cb, double tStep, int phase) {
    log->Write (Log::MAJOR, my_rank, "START:");

    chdir(cb->unit.component->workingDir);
    if (phase == GLOBAL_PHASE_FINISH) {
        log->Write (Log::MAJOR, my_rank, "%s.Initialize(%s)", cb->unit.component->name, cb->unit.component->inputFile);
        timerStart(cb->unit.component);
        cb->unit.component->result = (cb->unit.component->dllInitialize) (cb->unit.component->inputFile);
        timerEnd(cb->unit.component);
        (cb->unit.component->dllGetStartTime) (&cb->tStart);
        (cb->unit.component->dllGetEndTime) (&cb->tEnd);
    }
    cb->tStep = tStep;
    log->Write (Log::MAJOR, my_rank, "%s.Update(%6.1f)", cb->unit.component->name, cb->tStep);
    timerStart(cb->unit.component);
    (cb->unit.component->dllUpdate) (cb->tStep);
    timerEnd(cb->unit.component);
    if (phase == GLOBAL_PHASE_FINISH) {
        log->Write (Log::MAJOR, my_rank, "%s.Finalize()", cb->unit.component->name);
        timerStart(cb->unit.component);
        (cb->unit.component->dllFinalize) ();
        timerEnd(cb->unit.component);
    }
    fflush(stdout);
}



//------------------------------------------------------------------------------
void Dimr::runParallelInit (dimr_control_block * cb) {
    int ierr;
    MPI_Group mpiGroupWorld;
    MPI_Group mpiGroupComp;

    if (use_mpi) {
        ierr = MPI_Comm_group(MPI_COMM_WORLD, &mpiGroupWorld);
        if (ierr != MPI_SUCCESS) {
            throw new Exception(true, "runParallelInit: cannot obtain MPI world group. Code: %d.", ierr);
        }
    }

    // set masterSubBlockId
    for (int i = 0 ; i < cb->numSubBlocks ; i++) {
        if (cb->subBlocks[i].type == CT_START) {
            if (cb->masterSubBlockId == -1) {
                cb->masterSubBlockId = i;
                log->Write (Log::MINOR, my_rank, "Master: %s", cb->subBlocks[cb->masterSubBlockId].unit.component->name);
            } else {
                throw new Exception (true, "runParallelInit: a parallel block cannot have more than one start element.");
            }
        }
    }
    if (cb->masterSubBlockId == -1) {
        throw new Exception (true, "runParallelInit: a parallel block must have at least one start element.");
    }

    // Hack:
    // The masterComponent must be initialized first
    // Wave can only be initialized after the flow component
    dimr_component * masterComponent = cb->subBlocks[cb->masterSubBlockId].unit.component;

    // Create an MPI subgroup and subcommunicator and pass it on to the component.
    if (use_mpi && masterComponent->mpiCommVar != NULL  && masterComponent->numProcesses > 1) { // TODO: consider removing the numproc>1 check.
        ierr = MPI_Group_incl(mpiGroupWorld, masterComponent->numProcesses, masterComponent->processes, &mpiGroupComp);
        if (ierr != MPI_SUCCESS) {
            throw new Exception(true, "runParallelInit: cannot create a subgroup of %d processes for component \"%s\". Code: %d.", masterComponent->name, ierr);
        }
        // Needs to be called by *all* ranks:
        ierr = MPI_Comm_create(MPI_COMM_WORLD, mpiGroupComp, &masterComponent->mpiComm);
        if (ierr != MPI_SUCCESS) {
            throw new Exception(true, "runParallelInit: cannot create a subcommunicator of %d processes for component \"%s\". Code: %d.", masterComponent->name, ierr);
        }
        if (masterComponent->onThisRank) {
            MPI_Fint *fComm;
            masterComponent->dllGetVar(masterComponent->mpiCommVar, (void**)(&fComm));
            if (fComm == NULL) {
                throw new Exception(true, "runParallelInit: cannot obtain reference to communicator handle \"%s\" from component \"%s\".", masterComponent->mpiCommVar, masterComponent->name);
            }
            *fComm = MPI_Comm_c2f(masterComponent->mpiComm);
        }
    }
    if (masterComponent->onThisRank) {
        chdir(masterComponent->workingDir);
        log->Write (Log::MAJOR, my_rank, "%s.Initialize(%s)", masterComponent->name, masterComponent->inputFile);
        timerStart(masterComponent);
        masterComponent->result = (masterComponent->dllInitialize) (masterComponent->inputFile);
        timerEnd(masterComponent);
        (masterComponent->dllGetStartTime) (&cb->subBlocks[cb->masterSubBlockId].tStart);
        (masterComponent->dllGetEndTime) (&cb->subBlocks[cb->masterSubBlockId].tEnd);
        (masterComponent->dllGetTimeStep) (&cb->subBlocks[cb->masterSubBlockId].tStep);
        (masterComponent->dllGetCurrentTime) (&cb->subBlocks[cb->masterSubBlockId].tCur);
    }



    //
    // Then initialize the other components and couplers
    for (int i = 0 ; i < cb->numSubBlocks ; i++) {
        if (i != cb->masterSubBlockId) {
            // CT_STARTGROUP
            // First all components
            for (int j = 0 ; j < cb->subBlocks[i].numSubBlocks ; j++) {
                if (cb->subBlocks[i].subBlocks[j].type == CT_START) {
                    dimr_component * thisComponent = cb->subBlocks[i].subBlocks[j].unit.component;

                    if (thisComponent->onThisRank) { // TODO: AvD/AM: if FM is not start, but startblock, we need all the MPI stuff here as well: make a generic initializeComponent helper routine.

                        // Hack for WAVE:
                        if (thisComponent->type == COMP_TYPE_WAVE) {
                            int *waveModePtr        = NULL;
                            const char *key = "mode";
                            (thisComponent->dllGetVar) (key, (void**)(&waveModePtr));
                            *waveModePtr = 1;
                        }

                        chdir(thisComponent->workingDir);
                        log->Write (Log::MAJOR, my_rank, "%s.Initialize(%s)", thisComponent->name, thisComponent->inputFile);
                        timerStart(thisComponent);
                        thisComponent->result = (thisComponent->dllInitialize) (thisComponent->inputFile);
                        timerEnd(thisComponent);
                    }
                }
            }
            // Then all couplers
            for (int j = 0 ; j < cb->subBlocks[i].numSubBlocks ; j++) {
                if (cb->subBlocks[i].subBlocks[j].type != CT_START) {
                    dimr_coupler   * thisCoupler = cb->subBlocks[i].subBlocks[j].unit.coupler;
                    for (int k = 0 ; k < thisCoupler->numItems ; k++) {
                        if (thisCoupler->sourceComponent->type == COMP_TYPE_RTC   ||
                            thisCoupler->sourceComponent->type == COMP_TYPE_WANDA ||
                            thisCoupler->sourceComponent->type == COMP_TYPE_FLOW1D2D) {
                            // RTCTools/Wanda: impossible to autodetect which partition will deliver this source var
                            // Assumption: there is only one RTC-partition
                            thisCoupler->items[k].sourceProcess = thisCoupler->sourceComponent->processes[0];
                        } else {
                            // For each item: get the pointers to the variables inside the dlls to be exchanged
                            // Currently this does not work for RTC-Tools
                            //
                            // Source variable
                            // autodetect which (single!) partition will deliver this source var
                            int *  sources = (int *)malloc(thisCoupler->sourceComponent->numProcesses * sizeof(int));
                            int * gsources = (int *)malloc(thisCoupler->sourceComponent->numProcesses * sizeof(int));
                            for (int m = 0 ; m < thisCoupler->sourceComponent->numProcesses ; m++) {
                                sources[m] = 0;
                                if (my_rank == thisCoupler->sourceComponent->processes[m]) {
                                    // Also for RTCTools_BMI: this is a dummy getvar call, just to check whether it works for this partition
                                    log->Write (Log::MINOR, my_rank, "%s.getVar(%s)", thisCoupler->sourceComponentName, thisCoupler->items[k].sourceName);
                                    (thisCoupler->sourceComponent->dllGetVar) (thisCoupler->items[k].sourceName, (void**)(&thisCoupler->items[k].sourceVarPtr));
                                    if (thisCoupler->items[k].sourceVarPtr != NULL) {
                                        // Yes, this partition can deliver the source var
                                        sources[m] = 1;
                                    }
                                }
                            }
                            // Do not call MPI_Allreduce when the number of partitions is 1. It will cause a crash on free(gsources)
                            if (numranks > 1) {
                                int ierr = MPI_Allreduce(sources, gsources, thisCoupler->sourceComponent->numProcesses, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                            } else {
                                for (int m = 0 ; m < thisCoupler->sourceComponent->numProcesses ; m++) {
                                    gsources[m] = sources[m];
                                }
                            }
                            thisCoupler->items[k].sourceProcess = -1;
                            for (int m = 0 ; m < thisCoupler->sourceComponent->numProcesses ; m++) {
                                if (gsources[m] == 1) {
                                    if (thisCoupler->items[k].sourceProcess == -1) {
                                        // First partition that can deliver the source var
                                        thisCoupler->items[k].sourceProcess = m;
                                    } else {
                                        // Second/Third/... partition that can deliver the source var
                                        // Produce a warning
                                        // The "if (my_rank == m)" avoids multiple identical messages
                                        if (my_rank == m) {
                                            log->Write(Log::WARN, my_rank, "WARNING: coupler %s: item %d: \"%s\" will be delivered by partition %d. Ignoring deliverance by partition %d",
                                                thisCoupler->name, k, thisCoupler->items[k].sourceName, thisCoupler->items[k].sourceProcess, m);
                                        }
                                    }
                                }
                            }
                            free(sources);
                            free(gsources);
                        }

                        // Target variable

                        if (thisCoupler->targetComponent->type == COMP_TYPE_RTC     ||
                            thisCoupler->targetComponent->type == COMP_TYPE_WANDA   ||
                            thisCoupler->targetComponent->type == COMP_TYPE_FLOW1D2D) {
                            // nothing
                        } else {
                            // Target variable
                            // autodetect which (possibly multiple!) partition(s) will accept this target var
                            int *  targets = (int *)malloc(thisCoupler->targetComponent->numProcesses * sizeof(int));
                            int * gtargets = (int *)malloc(thisCoupler->targetComponent->numProcesses * sizeof(int));
                            for (int m = 0 ; m < thisCoupler->targetComponent->numProcesses; m++) {
                                targets[m] = 0;
                                if (my_rank == thisCoupler->targetComponent->processes[m]) {
                                        log->Write (Log::MINOR, my_rank, "%s.getVar(%s)", thisCoupler->targetComponentName, thisCoupler->items[k].targetName);
                                        (thisCoupler->targetComponent->dllGetVar) (thisCoupler->items[k].targetName, (void**)(&thisCoupler->items[k].targetVarPtr));
                                        if (thisCoupler->items[k].targetVarPtr != NULL) {
                                            // Yes, this partition can accept the target var
                                            targets[m] = 1;
                                        }
                                }
                            }
                            // Do not call MPI_Allreduce when the number of partitions is 1. It will cause a crash on free(gtargets)
                            if (numranks > 1) {
                                int ierr = MPI_Allreduce(targets, gtargets, thisCoupler->targetComponent->numProcesses, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                            } else {
                                for (int m = 0 ; m < thisCoupler->targetComponent->numProcesses ; m++) {
                                    gtargets[m] = targets[m];
                                }
                            }
                            thisCoupler->items[k].targetProcess = -1;
                            for (int m = 0 ; m < thisCoupler->targetComponent->numProcesses ; m++) {
                                if (gtargets[m] == 1) {
                                    if (thisCoupler->items[k].targetProcess == -1) {
                                        // First partition that can accept the target var
                                        thisCoupler->items[k].targetProcess = m;
                                    } else {
                                        // Second/Third/... partition that can accept the target var
                                        // Produce a warning
                                        // The "if (my_rank == m)" avoids multiple identical messages
                                        if (my_rank == m) {
                                            log->Write(Log::WARN, my_rank, "WARNING: coupler %s: item %d: \"%s\" will be delivered to multiple partitions: %d and %d.",
                                                thisCoupler->name, k, thisCoupler->items[k].targetName, thisCoupler->items[k].targetProcess, m);
                                        }
                                    }
                                }
                            }
                            if (thisCoupler->items[k].targetProcess == -1) {
                                throw new Exception (true, "Coupler %s: item %d: \"%s\" is not accepted by any of the partitions.",
                                    thisCoupler->name, k, thisCoupler->items[k].targetName);
                            }
                            free(targets);
                            free(gtargets);
                        }
                    }
                }
            }
        }
    }
}



//------------------------------------------------------------------------------
void Dimr::runParallelUpdate (dimr_control_block * cb, double tStep) {
    dimr_control_block * masterComponent = &cb->subBlocks[cb->masterSubBlockId];
    if (!masterComponent->unit.component->onThisRank) {
        throw new Exception (true, "runParallelUpdate: not supported yet: master component \"%s\" should run on all processes.", masterComponent->unit.component->name);
        // TODO: AvD/AM: is this allowed: master component not on *all* dimr processes? NOT YET, but yes we want it, e.g. 3xFM, 7xWAVE. Rethink.
    }
    // Initialize time parameters
    double * currentTime = &masterComponent->tCur;
    // masterComponent->tStart : simulation start time as obtained from the masterComponent
    // masterComponent->tEnd   : simulation end   time as obtained from the masterComponent
    (masterComponent->unit.component->dllGetStartTime) (&masterComponent->tStart);
    (masterComponent->unit.component->dllGetEndTime) (&masterComponent->tEnd);
    masterComponent->tNext = min(*currentTime + tStep, masterComponent->tEnd);
    if (*currentTime == masterComponent->tStart) {
        // Set the currentTime and nextTime in all other components, relative to currentTime
        for (int i = 0 ; i < cb->numSubBlocks ; i++) {
            if (i != cb->masterSubBlockId) {
                // subBlock.tCur:
                //     componentType = WAVE:
                //         Assumption: Reference date is equal to masterComponent's reference date
                //         => : tCur = 0.0
                //     componentType = RTC:
                //         tCur is not used (update is always called with argument -1)
                cb->subBlocks[i].tCur  = 0.0;
                //
                // subBlock.tNext:
                //     tNext must be tStart, also when tStart=0,
                //     to distinguish "start at begin of simulatione" and "after one time step"
                cb->subBlocks[i].tNext = *currentTime + cb->subBlocks[i].tStart;
                cb->subBlocks[i].tEnd  = *currentTime + cb->subBlocks[i].tEnd;
            }
        }
    }
    //
    // TIME LOOP
    //
    while (*currentTime < masterComponent->tNext) {
        //
        // define tStep
        // Start with maximum value defined by masterComponent
        double tStep = masterComponent->tNext - *currentTime;
        // tStep is the minimum allowed time step over all followers
        for (int i = 0 ; i < cb->numSubBlocks ; i++) {
            if (i != cb->masterSubBlockId) {
                // follower.tNext is the next point in time that this follower must be executed
                // For this follower, the allowed timestep is "follower.tNext - currentTime"
                double tStepFollower = cb->subBlocks[i].tNext - *currentTime;
                if (tStepFollower == 0.0) {
                    // This follower is already active
                    if (i < cb->masterSubBlockId) {
                        // Before MasterComponent:
                        // Check the step size instead of tNext
                        tStep = min(tStep,cb->subBlocks[i].tStep);
                    } else {
                        // After MasterComponent:
                        // When the MasterComponent is going to run, it will increase currentTime
                        // This does not match with executing this follower at the current time
                        throw new Exception (true, "runParallelUpdate: Zero timestep, needed for block %d, is not possible.", i);
                    }
                } else {
                    // This follower is not active yet
                    tStep = min(tStep,tStepFollower);
                }
            }
        }

        //
        // Update all components in the specified order
        for (int i = 0 ; i < cb->numSubBlocks ; i++) {
            // Sync all partitions to execute the same component
            // This ensures that all flow calculations are finished before a wave calculation is started
            if (use_mpi) {
                int ierr = MPI_Barrier(MPI_COMM_WORLD);
            }

            if (i == cb->masterSubBlockId) {
                // masterComponent
                chdir(masterComponent->unit.component->workingDir);
                log->Write (Log::MAJOR, my_rank, "%10.1f:    %s.Update(%10.1f)", *currentTime, masterComponent->unit.component->name, tStep);
                timerStart(masterComponent->unit.component);
                (masterComponent->unit.component->dllUpdate) (tStep);
                timerEnd(masterComponent->unit.component);
                *currentTime = *currentTime + tStep;
            } else {
                // CT_STARTGROUP
                if (*currentTime >= cb->subBlocks[i].tNext) {
                    // Yes, it's time to execute all childBlocks of this subBlock
                    for (int j = 0 ; j < cb->subBlocks[i].numSubBlocks ; j++) {
                        if (cb->subBlocks[i].subBlocks[j].type == CT_START) {
                            // Component
                            dimr_component * thisComponent = cb->subBlocks[i].subBlocks[j].unit.component;
                            if (!thisComponent->onThisRank) {
                                continue;
                            }
                            //
                            // tUpdate is the timeInterval since the last time this component was executed
                            // This is not the overall tStep!
                            double tUpdate;
                            tUpdate = *currentTime - cb->subBlocks[i].tCur;
                            // Hack: Always call RTCTools.Update with argument -1.0
                            if (cb->subBlocks[i].subBlocks[j].unit.component->type == COMP_TYPE_RTC) {
                                tUpdate = -1.0;
                            }
                            // Hack: Wanda: The first time that WANDA.update(tUpdate) is called, tUpdate must be 0.0
                            //              Otherwise, WANDA will hang in case tUpdate is very big (FLOW start time big, related to ref time)
                            //              It doesn't matter for the WANDA results; with this coupling, WANDA doesn't have a time knowledge at all
                            if (cb->subBlocks[i].subBlocks[j].unit.component->type == COMP_TYPE_WANDA) {
                                if (cb->subBlocks[i].tCur == 0.0) {
                                    tUpdate = 0.0;
                                }
                            }
                            // Update
                            chdir(thisComponent->workingDir);
                            log->Write (Log::MAJOR, my_rank, "%10.1f:    %s.Update(%10.1f)", *currentTime, thisComponent->name, tUpdate);
                            timerStart(thisComponent);
                            (thisComponent->dllUpdate) (tUpdate);
                            timerEnd(thisComponent);
                        } else {
                            // Coupler
                            dimr_coupler * thisCoupler = cb->subBlocks[i].subBlocks[j].unit.coupler;
                            double * transferValuePtr;
                            log->Write (Log::MINOR, my_rank, "%10.1f:    %s.communicate", *currentTime, thisCoupler->name);
                            for (int k = 0 ; k < thisCoupler->numItems ; k++) {
                                log->Write (Log::DETAIL, my_rank, "    %s -> %s", thisCoupler->items[k].sourceName, thisCoupler->items[k].targetName);
                                // Getting and Setting of data is split to enable
                                // transferring data, possibly inbetween different partitions
                                // TODO: This does not work for arrays (yet), only scalar double
                                //
                                // Getting data:
                                double * transfer = new double [thisCoupler->sourceComponent->numProcesses];

                                // addresses eventually updated
                                getAddress(thisCoupler->items[k].sourceName, thisCoupler->sourceComponent->type, thisCoupler->sourceComponent->dllGetVar, &(thisCoupler->items[k].sourceVarPtr), 
                                   thisCoupler->sourceComponent->processes, thisCoupler->sourceComponent->numProcesses, transfer);
                                // the value of sourceVarPtr is sent here to all MPI processes
                                transferValuePtr = send(thisCoupler->items[k].sourceName, thisCoupler->sourceComponent->type, thisCoupler->items[k].sourceVarPtr, 
                                   thisCoupler->sourceComponent->processes, thisCoupler->sourceComponent->numProcesses, transfer);

                                delete[] transfer;

                                // Optional TODO: assuming one source(partition):
                                //    if there is only one target (partition), and all partitions can act as source partition:
                                //        choose the target partition to act as source partition
                                //        no MPI_Bcast needed
                                //
                                receive (thisCoupler->items[k].targetName,
                                         thisCoupler->targetComponent->type,
                                         thisCoupler->targetComponent->dllSetVar,
                                         thisCoupler->targetComponent->dllGetVar,
                                         thisCoupler->items[k].targetVarPtr,
                                         thisCoupler->targetComponent->processes,
                                         thisCoupler->targetComponent->numProcesses,
                                         thisCoupler->items[k].targetProcess,
                                         transferValuePtr);
                            }
                        }
                    }
                    // All child blocks are executed. Update tCur and define the next time step that this block has to be executed.
                    cb->subBlocks[i].tCur  = *currentTime;
                    cb->subBlocks[i].tNext = cb->subBlocks[i].tNext + cb->subBlocks[i].tStep;
                    if (cb->subBlocks[i].tNext > cb->subBlocks[i].tEnd)
                        // This subBlock does not have to be executed anymore
                        // Force this by giving it a nextTime > simulationEndTime
                        cb->subBlocks[i].tNext = 2.0 * masterComponent->tEnd;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
// set "value" in the component target location
void Dimr::receive(const char * name,
   int          compType,
   BMI_SETVAR   dllSetVar,
   BMI_GETVAR   dllGetVar,
   double     * targetVarPtr,
   int        * processes,
   int          nProc,
   int          targetProcess,
   const void * transferValuePtr) {
   // Assumption: transferValuePtr points to a (scalar) double
   // TODO: allow more types of pointers
   //
   // NaN check:
   if (*(double*)(transferValuePtr) == *(double*)(transferValuePtr)) {
      for (int m = 0; m < nProc; m++) {
         if (my_rank == processes[m]) {
            log->Write(Log::DETAIL, my_rank, "Dimr::receive(%s) %20.10f", name, *(double*)(transferValuePtr));
            if (compType == COMP_TYPE_RTC
               || compType == COMP_TYPE_RR
               || compType == COMP_TYPE_FLOW1D
               || compType == COMP_TYPE_FLOW1D2D
               || compType == COMP_TYPE_WANDA) {
               if (dllSetVar == NULL) {
                  throw new Exception(true, "ABORT: Dimr::receive: set_var function not defined while processing %s", name);
               }
               (dllSetVar)(name, (const void *)transferValuePtr);
               if (compType == COMP_TYPE_RTC) {
                  // target = rtc
                  // SetVar(name, value) sets variable named "name" to "value" at the current time (t = n)
                  // But, in case of IMPLICIT method, this should be the next time (t = n + 1)
                  // Clean solution: do not use IMPLICIT method in RTCTools (To Do for Stef Hummel)
                  // Shortcut hack:
                  // 1. Set value at current time (t=n)
                  // AND 2. Set also value at next time (t=n+1),
                  // by adding a * at the end of the targetName. This is a signal to RTCTools that this
                  // value belongs to t=n+1
                  // This "shortcut hack" is identical to what currently is used in Delta Shell via the
                  // OpenMI interface with RTCTools
                  char *targetName = new char[strlen(name) + 2];
                  sprintf(targetName, "%s*\0", name);
                  (dllSetVar)(targetName, (const void *)transferValuePtr);
                  delete[] targetName;
               }
            }
            else {
               // target is a component that uses direct pointer access to the actual variable
               // When doing a dllSetVar, targetVarPtr is not defined. First do a "getAddress" to get it defined
               // NOTE: transferValuePtr contains the value calculated by RTCTools: it has been already communicated 
               // to all processes by a previous send call. 

               if (targetVarPtr == NULL)
               {
                  double * transfer = new double[nProc];
                  //here we get the address (e.g. weir levels)
                  getAddress(name, compType, dllGetVar, &targetVarPtr, processes, nProc, transfer);
                  delete[] transfer;
               }
 
               //  }
               // Know we know the location, we need to set it in targetVarPtr
               //here we do the check
               if (targetVarPtr == NULL)
               {
                  if (targetProcess == -1 || targetProcess == my_rank) 
                  {
                     // targetProcess=-1: no process can accept this item
                     // targetProcess=my_rank: this process is registered to be able to accept this item but something goes wrong
                     throw new Exception(true, "ABORT: Dimr::receive: get_var function not defined while processing %s", name);
                  }
               }
               else
               {
                  // Here write the RTC value in FM (direct access to memory using a pointer). 
                  // Here we already know targetVarPtr is not null. 
                     *(targetVarPtr) = *(double *)transferValuePtr;
               }
            }
         }
      }
   }
}

//------------------------------------------------------------------------------
// This function gets the address from the components. 
// Does not do any transfer of the values stored in sourceVarPtr. 

void Dimr::getAddress(
   const char * name,
   int          compType,
   BMI_GETVAR   dllGetVar,
   double    ** sourceVarPtr,
   int        * processes,
   int          nProc,
   double     * transfer)
{
   for (int m = 0; m < nProc; m++) {
      if (my_rank == processes[m]) {
         log->Write(Log::DETAIL, my_rank, "Dimr::getAddress (%s)", name);
         if (compType == COMP_TYPE_RTC
            || compType == COMP_TYPE_RR
            || compType == COMP_TYPE_FLOW1D
            || compType == COMP_TYPE_FLOW1D2D
            || *sourceVarPtr == NULL) {
            // These components only returns a new pointer to a copy of the double value, so call it each time.
            // sourceVarPtr=NULL: getVar not yet called for this parameter, probably because "send" is being called
            //                    via the toplevel "get_var"
            if (dllGetVar == NULL) {
               throw new Exception(true, "ABORT: get_var function not defined while processing %s", name);
            }
            (dllGetVar)(name, (void**)(sourceVarPtr));
         }
         else if (compType == COMP_TYPE_WANDA) {
            if (dllGetVar == NULL) {
               throw new Exception(true, "ABORT: get_var function not defined while processing %s", name);
            }
            // Wanda does not use pointers to internal structures:
            // - Use the DIMR-transfer array
            // - Note the missing & in the dllGetVar call, when comparing with the dllGetVar call above
            *sourceVarPtr = transfer;
            (dllGetVar)(name, (void**)(*sourceVarPtr));
         }
         else 
         {
            // Other components already have direct pointer access to the actual variable.
         }
      }
   }
}

//------------------------------------------------------------------------------
// This function sends the sourceVarPtr value to all processes 
double* Dimr::send(
   const char * name,
   int          compType,
   double     * sourceVarPtr,
   int        * processes,
   int          nProc,
   double     * transfer)
{
   double * reducedTransfer = new double[nProc];

   for (int m = 0; m < nProc; m++) 
   {
      transfer[m] = -999000.0;
      if (my_rank == processes[m]) 
      {
         log->Write(Log::DETAIL, my_rank, "Dimr::send (%s)", name);
         if (sourceVarPtr != NULL && compType != COMP_TYPE_WANDA) 
         {
            transfer[m] = *sourceVarPtr;
         }
      }
   }

   //now you have all source values in the transfer array, we can reduce them
   double maxValue = -999000.0;
   // Do not call MPI_Allreduce when the number of partitions is 1. 
   if (numranks > 1)
   {
      // NOTE: here the transfer array has a defined value only at position==my_rank and if *sourceVarPtr != NULL.
      // We perform a reduction operation to collect the maximum values at each position and afterwards take the maximum of all values.
      // We could also use only one double instead of a transfer array and perform the reduction on a single scalar (so avoiding the loop below) 
      int ierr = MPI_Allreduce(transfer, reducedTransfer, nProc, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      for (int m = 0; m < nProc; m++)
      {
         if (reducedTransfer[m] > maxValue) 
            maxValue = reducedTransfer[m];
      }
   }
   else
   {
      maxValue = transfer[0];
   }

   delete [] reducedTransfer;

   //the reduced value now is set and return to all ranks
   transferValue = maxValue;
   return &(transferValue);
}


//------------------------------------------------------------------------------
void Dimr::runParallelFinish (dimr_control_block * cb) {
    if (use_mpi) {
        int ierr = MPI_Barrier(MPI_COMM_WORLD);
    }

    for (int i = 0 ; i < cb->numSubBlocks ; i++) {
        if (cb->subBlocks[i].type == CT_START) {
            if (!cb->subBlocks[i].unit.component->onThisRank) {
                continue;
            }
            chdir(cb->subBlocks[i].unit.component->workingDir);
            log->Write (Log::MAJOR, my_rank, "    %s.Finalize()", cb->subBlocks[i].unit.component->name);
            timerStart(cb->subBlocks[i].unit.component);
            (cb->subBlocks[i].unit.component->dllFinalize) ();
            timerEnd(cb->subBlocks[i].unit.component);

            if (use_mpi && cb->subBlocks[i].unit.component->mpiComm != NULL) {
                MPI_Group mpiGroupComp;
                int ierr = MPI_Comm_group(cb->subBlocks[i].unit.component->mpiComm, &mpiGroupComp);
                ierr = MPI_Group_free(&mpiGroupComp);
                ierr = MPI_Comm_free(&cb->subBlocks[i].unit.component->mpiComm);
            }
        } else {
            // CT_STARTGROUP
            for (int j = 0 ; j < cb->subBlocks[i].numSubBlocks ; j++) {
                if (cb->subBlocks[i].subBlocks[j].type == CT_START) {
                    if (!cb->subBlocks[i].subBlocks[j].unit.component->onThisRank) {
                        continue;
                    }
                    chdir(cb->subBlocks[i].subBlocks[j].unit.component->workingDir);
                    log->Write (Log::MAJOR, my_rank, "    %s.Finalize()", cb->subBlocks[i].subBlocks[j].unit.component->name);
                    timerStart(cb->subBlocks[i].subBlocks[j].unit.component);
                    (cb->subBlocks[i].subBlocks[j].unit.component->dllFinalize) ();
                    timerEnd(cb->subBlocks[i].subBlocks[j].unit.component);
                    }
            }
        }
    }
}



//------------------------------------------------------------------------------
void Dimr::scanConfigFile (void) {

    XmlTree * rootXml     = config->Lookup ("/dimrConfig");
    if (rootXml == NULL)
        throw new Exception (true, "Configuration file \"%s\" does not have a <dimrConfig> root element", configfile);
    XmlTree * fileversion = rootXml->Lookup ("documentation/fileVersion");
    if (fileversion == NULL)
        throw new Exception (true, "Configuration file \"%s\" does not have a deltaresHydro documentation->fileVersion element", configfile);

    // Check version number
    const char * versionnr = fileversion->charData;
    float versionnumber;
    int intRead = sscanf(versionnr, "%f", &versionnumber);
    if (intRead != 1)
        throw new Exception (true, "Configuration file \"%s\" does not have a version number", configfile);
    if ((int)floor(versionnumber) != 1)
        throw new Exception (true, "Configuration file \"%s\": Version number (%3.2f) must have main version 2", configfile,versionnumber);

    XmlTree * controlXml  = rootXml->Lookup ("control");
    if (controlXml == NULL)
        throw new Exception (true, "Configuration file \"%s\" does not have a deltaresHydro control element", configfile);
    // Allocate the control structure and check its size
    control = (dimr_control_block *) malloc(sizeof(dimr_control_block));
    control->numSubBlocks = 0;
    control->subBlocks    = NULL;
    control->masterSubBlockId = -1;
    // First scan the config file for all components and couplers (= units)
    scanUnits(rootXml);
    // Then scan the control part
    // References are added to the list of components/couplers
    scanControl(controlXml, control);
}



//------------------------------------------------------------------------------
void Dimr::scanUnits(XmlTree * rootXml) {
    // Init
    componentsList.numComponents = 0;
    componentsList.components    = NULL;
    couplersList.numCouplers     = 0;
    couplersList.couplers        = NULL;
    // Scan
    for (int i = 0 ; i < rootXml->numChildren ; i++) {
        if (strcmp(rootXml->children[i]->name, "component") == 0) {
            componentsList.numComponents++;
            if (componentsList.components == NULL) {
                componentsList.components = (dimr_component*)malloc(componentsList.numComponents * sizeof(dimr_component));
            } else {
                componentsList.components = (dimr_component*)realloc(componentsList.components, 
                                                                          componentsList.numComponents * sizeof(dimr_component));
                if (componentsList.components == NULL) {
                    throw new Exception (true, "Allocation error in scanUnits (component)");
                }
            }
            scanComponent(rootXml->children[i], &(componentsList.components[componentsList.numComponents - 1]));
        }
        if (strcmp(rootXml->children[i]->name, "coupler") == 0) {
            couplersList.numCouplers++;
            if (couplersList.couplers == NULL) {
                couplersList.couplers = (dimr_coupler*)malloc(couplersList.numCouplers * sizeof(dimr_coupler));
            } else {
                couplersList.couplers = (dimr_coupler*)realloc(couplersList.couplers, couplersList.numCouplers * sizeof(dimr_coupler));
                if (couplersList.couplers == NULL) {
                    throw new Exception (true, "Allocation error in scanUnits (coupler)");
                }
            }
            scanCoupler(rootXml->children[i], &(couplersList.couplers[couplersList.numCouplers - 1]));
        }
    }
}



//------------------------------------------------------------------------------
void Dimr::scanComponent(XmlTree * xmlComponent, dimr_component * newComp) {
    // Needed for path handling
#if defined(HAVE_CONFIG_H)
    const char *dirSeparator = "/";
#else
    const char *dirSeparator = "\\";
#endif
    char *curPath = new char[MAXSTRING];
    if (!getcwd(curPath, MAXSTRING))
        throw new Exception (true, "ERROR obtaining the current working directory (scan)");
    //
    //
    newComp->name = xmlComponent->GetAttrib("name");
    // Element library
    XmlTree * libraryElement = xmlComponent->Lookup ("library");
    if (libraryElement == NULL)
        throw new Exception (true, "Component \"%s\" does not contain a library element", newComp->name);
    newComp->library = libraryElement->charData;
    int libLen = strlen(newComp->library);
    char *libNameLowercase= new char[libLen+1];
    strncpy(libNameLowercase, newComp->library, libLen);
    libNameLowercase[libLen] = '\0';
    for (int i=0; i < libLen; i++) {
        libNameLowercase[i] = (tolower(libNameLowercase[i]));
    }
    if (strstr(libNameLowercase, "rtc") != NULL) {
        newComp->type = COMP_TYPE_RTC;
    } else if (strstr(libNameLowercase, "wave") != NULL) {
        newComp->type = COMP_TYPE_WAVE;
    } else if (strstr(libNameLowercase, "fm") != NULL){
        newComp->type = COMP_TYPE_FM;
    } else if (strstr(libNameLowercase, "cf_dll") != NULL){
        newComp->type = COMP_TYPE_FLOW1D;
    } else if (strstr(libNameLowercase, "rr_dll") != NULL){
        newComp->type = COMP_TYPE_RR;
    } else if (strstr(libNameLowercase, "wandaengine_native") != NULL){
        newComp->type = COMP_TYPE_WANDA;
    } else if (strstr(libNameLowercase, "flow2d3d") != NULL){
        newComp->type = COMP_TYPE_FLOW2D3D;
    } else if (strstr(libNameLowercase, "flow1d2d") != NULL){
       newComp->type = COMP_TYPE_FLOW1D2D;
    } else if (strstr(libNameLowercase, "delwaq") != NULL){
       newComp->type = COMP_TYPE_DELWAQ;
    } else if (strstr(libNameLowercase, "dimr_testcomponent") != NULL){
       newComp->type = COMP_TYPE_TEST;
    }
    else {
        throw new Exception (true, "Name of library, \"%s\", is not recognized", newComp->library);
    }
    delete [] libNameLowercase;

    // Element process (optional)
    XmlTree * processElement = xmlComponent->Lookup ("process");
    if (processElement != NULL) {
        // Store process rank numbers in component's processes array.
        char_to_ints(processElement->charData, &(newComp->processes), &(newComp->numProcesses));

        // Check whether this process' rank is also in components configured processes array.
        newComp->onThisRank = false;         // Not found (yet): only active on other ranks.
        for (int i=0; i < newComp->numProcesses; i++) {
            if (newComp->processes[i] >= numranks) {
                throw new Exception(true, "Component \"%s\" configured for process #%d, but max running MPI rank is only %d.",
                    newComp->name, newComp->processes[i], numranks-1);
            } else  if (newComp->processes[i] == my_rank) {
                newComp->onThisRank = true;  // Found:     active.
            }
        }

    } else {
        // No <process> specified, default: only run on rank #0.
        log->Write (Log::DETAIL, my_rank, "INFO: \"<process>\" not specified for component \"%s\". Assuming it only runs on rank #0.", newComp->name);
        newComp->numProcesses = 1;
        char *defaultProc = "0";
        char_to_ints(defaultProc, &(newComp->processes), &(newComp->numProcesses));

        newComp->onThisRank = (my_rank==0);
    }

    // parameter (optional)
    int nparameter = xmlComponent->Lookup ("parameter",0,newComp->parameters);

    // setting (optional)
    int nsetting = xmlComponent->Lookup ("setting",0,newComp->settings);

    // Element mpiCommunicator (optional)
    XmlTree * commElement = xmlComponent->Lookup ("mpiCommunicator");
    if (commElement != NULL) {
        // Store communicator var name in component.
        newComp->mpiCommVar = commElement->charData;
    }

    // Element inputFile (optional?)
    XmlTree * inputFileElement = xmlComponent->Lookup ("inputFile");
    if (inputFileElement == NULL) {
        log->Write (Log::MAJOR, my_rank, "WARNING: No inputFile specified for component %s.", newComp->name);
        }
    else {
        newComp->inputFile = inputFileElement->charData;
        }
    // Element workingDir
    XmlTree * workingDirElement = xmlComponent->Lookup ("workingDir");
    if (workingDirElement == NULL) {
        newComp->workingDir = (char *) &curPath;
        log->Write (Log::MAJOR, my_rank, "WARNING: No workingDir specified for component %s.", newComp->name);
        log->Write (Log::MAJOR, my_rank, "         workingDir is set to %s", newComp->workingDir);
        }
    else {
        newComp->workingDir = workingDirElement->charData;
        }
    // Is workingDir a valid relative path?
    char *combinedPath= new char[MAXSTRING];
    sprintf(combinedPath, "%s%s%s", curPath, dirSeparator, newComp->workingDir);
    if (chdir(combinedPath)) {
        // CombinedPath is not correct. May be just workingDir?
        delete [] combinedPath;
        // Is workingDir a valid absolute path?
        if (chdir(newComp->workingDir)) {
            throw new Exception (true, "Component \"%s\" has an invalid workingDir \"%s\"", newComp->name, newComp->workingDir);
        }
    }
    else {
        // newComp->workingDir was a pointer to workingDirElement->charData; now it will point to the new combinedPath
        newComp->workingDir = combinedPath;

    }
    chdir(curPath);

    // Hack for RTC-Tools:
    // inputFile is an input directory
    // "." means workingDir
    if (strcmp(newComp->inputFile, ".") == 0) {
        newComp->inputFile = newComp->workingDir;
    }
    delete [] curPath;
}



//------------------------------------------------------------------------------
void Dimr::scanCoupler(XmlTree * xmlCoupler, dimr_coupler * newCoup) {
    newCoup->name = xmlCoupler->GetAttrib("name");
    // Element sourceComponent
    XmlTree * sourceComponent = xmlCoupler->Lookup ("sourceComponent");
    if (sourceComponent == NULL)
        throw new Exception (true, "The coupler \"%s\" does not contain a sourceComponent element", newCoup->name);
    newCoup->sourceComponentName = sourceComponent->charData;
    // Add reference to the actual component acting as source
    newCoup->sourceComponent = getComponent(newCoup->sourceComponentName);
    // Element targetComponent
    XmlTree * targetComponent = xmlCoupler->Lookup ("targetComponent");
    if (targetComponent == NULL)
        throw new Exception (true, "The coupler \"%s\" does not contain a targetComponent element", newCoup->name);
    newCoup->targetComponentName = targetComponent->charData;
    // Add reference to the actual component acting as target
    newCoup->targetComponent = getComponent(newCoup->targetComponentName);
    // Items
    newCoup->numItems = 0;
    newCoup->items    = NULL;
    for (int j = 0 ; j < xmlCoupler->numChildren ; j++) {
        if (strcmp(xmlCoupler->children[j]->name, "item") == 0) {
            // Create the item
            newCoup->numItems++;
            if (newCoup->items == NULL) {
                newCoup->items = (dimr_couple_item*)malloc(newCoup->numItems * sizeof(dimr_couple_item));
            } else {
                newCoup->items = (dimr_couple_item*)realloc(newCoup->items, newCoup->numItems * sizeof(dimr_couple_item));
                if (newCoup->items == NULL) {
                    throw new Exception (true, "Allocation error in scanUnits (couple unit)");
                }
            }
            dimr_couple_item *newItem = &(newCoup->items[newCoup->numItems - 1]);

            // Read sourceName
            XmlTree * xmlSource = xmlCoupler->children[j]->Lookup ("sourceName");
            if (xmlSource == NULL)
                throw new Exception (true, "The coupler \"%s\", item %d, does not contain a sourceName element", newCoup->name, newCoup->numItems);
            newItem->sourceName = xmlSource->charData;
            if (newItem->sourceName == NULL)
                throw new Exception (true, "Item %d of coupler \"%s\" does not contain a source::name element", newCoup->numItems, newCoup->name);

            // Read targetName
            XmlTree * xmlTarget = xmlCoupler->children[j]->Lookup ("targetName");
            if (xmlTarget == NULL)
                throw new Exception (true, "The coupler \"%s\", item %d, does not contain a targetName element", newCoup->name, newCoup->numItems);
            newItem->targetName = xmlTarget->charData;
            if (newItem->targetName == NULL)
                throw new Exception (true, "Item %d of coupler \"%s\" does not contain a target::name element", newCoup->numItems, newCoup->name);

            // source/targetVarPtr will be set in runParallelInit
            newItem->sourceVarPtr = NULL;
            newItem->targetVarPtr = NULL;
        }
    }
}


//------------------------------------------------------------------------------
void Dimr::scanControl(XmlTree * controlBlockXml, dimr_control_block * controlBlock) {
    if (strcmp(controlBlockXml->name, "control") == 0) {
        controlBlock->type = CT_SEQUENTIAL;
    } else if (strcmp(controlBlockXml->name, "parallel") == 0) {
        controlBlock->type = CT_PARALLEL;
    } else if (strcmp(controlBlockXml->name, "start") == 0) {
        controlBlock->type = CT_START;
        controlBlock->unit.component = getComponent(controlBlockXml->GetAttrib("name"));
        controlBlock->unit.coupler = NULL;
    } else if (strcmp(controlBlockXml->name, "startGroup") == 0) {
        controlBlock->type = CT_STARTGROUP;
        XmlTree * timeElt = controlBlockXml->Lookup ("time");
        if (timeElt == NULL)
            throw new Exception (true, "The startGroup component \"%s\" does not contain a time element", controlBlockXml->name);
        int intRead = sscanf(timeElt->charData, "%lf %lf %lf", &(controlBlock->tStart), &(controlBlock->tStep), &(controlBlock->tEnd));
        if (intRead != 3)
            throw new Exception (true, "Cannot find tStart, tStep, tEnd");
    } else if (strcmp(controlBlockXml->name, "coupler") == 0) {
        controlBlock->type = CT_COUPLER;
        controlBlock->unit.component = NULL;
        controlBlock->unit.coupler = getCoupler(controlBlockXml->GetAttrib("name"));
    }
    controlBlock->numSubBlocks = 0;
    controlBlock->subBlocks    = NULL;
    controlBlock->masterSubBlockId = -1;
    for (int i = 0 ; i < controlBlockXml->numChildren ; i++) {
        if (   strcmp(controlBlockXml->children[i]->name, "parallel"  ) == 0
            || strcmp(controlBlockXml->children[i]->name, "start"     ) == 0
            || strcmp(controlBlockXml->children[i]->name, "startGroup") == 0
            || strcmp(controlBlockXml->children[i]->name, "coupler"   ) == 0) {
            controlBlock->numSubBlocks++;
            if (controlBlock->subBlocks == NULL) {
                controlBlock->subBlocks = (dimr_control_block*)malloc(controlBlock->numSubBlocks * sizeof(dimr_control_block));
            } else {
                controlBlock->subBlocks = (dimr_control_block*)realloc(controlBlock->subBlocks, controlBlock->numSubBlocks * sizeof(dimr_control_block));
                if (controlBlock->subBlocks == NULL) {
                    throw new Exception (true, "Allocation error in scanControl");
                }
            }
            scanControl(controlBlockXml->children[i], &(controlBlock->subBlocks[controlBlock->numSubBlocks - 1]));
        }
        fflush(stdout);
    }
}



//------------------------------------------------------------------------------
// Search for a named component in the list of components
dimr_component * Dimr::getComponent(const char * compName) {
    for (int i = 0 ; i < componentsList.numComponents ; i++) {
        if (strcmp(componentsList.components[i].name, compName) == 0) {
            return &(componentsList.components[i]);
        }
    }
}



//------------------------------------------------------------------------------
// Search for a named coupler in the list of couplers
dimr_coupler * Dimr::getCoupler(const char * coupName) {
    for (int i = 0 ; i < couplersList.numCouplers ; i++) {
        if (strcmp(couplersList.couplers[i].name, coupName) == 0) {
            return &(couplersList.couplers[i]);
        }
    }
}



//------------------------------------------------------------------------------
void Dimr::connectLibs (void) {

    //          linux windows   mac
    // lib        so    dll     dylib
    // module     so    dll     so

#if defined (OSX)
    // Macintosh:VERY SIMILAR TO LINUX
    throw new Exception (true, "ABORT: %s has not be ported to Apple Mac OS/X yet", exeName);
#endif
#if defined (HAVE_CONFIG_H)
    char *err;
#endif

    // do for all libraries
    for (int i = 0 ; i < componentsList.numComponents ; i++) {
        if (!componentsList.components[i].onThisRank) {
            log->Write (Log::DETAIL, my_rank, "Not necessary to load component library \"%s\" on this rank.", componentsList.components[i].library);
            continue;
        }

#if defined (HAVE_CONFIG_H)
        char * lib = new char[strlen (componentsList.components[i].library) + 3+3+1];
        sprintf (lib, "lib%s%s", componentsList.components[i].library, D3D_PLUGIN_EXT);
        if (   strchr (componentsList.components[i].library, '/' ) != NULL 
            || strchr (componentsList.components[i].library, '\\') != NULL 
            || strchr (componentsList.components[i].library, '.' ) != NULL) {
            throw new Exception (true, "Invalid component library name \"%s\"\n", lib, -1);
        }
#else
        char * lib = new char[strlen (componentsList.components[i].library) + 4+1];
        sprintf (lib, "%s.dll", componentsList.components[i].library);
        if (   strchr (componentsList.components[i].library, '/' ) != NULL 
            || strchr (componentsList.components[i].library, '\\') != NULL 
            || strchr (componentsList.components[i].library, '.' ) != NULL) {
            throw new Exception (true, "Invalid component library name \"%s\"\n", lib, -1);
        }
#endif

        log->Write (Log::DETAIL, my_rank, "Loading library \"%s\"", lib);

#if defined (HAVE_CONFIG_H)
        dlerror(); /* clear error code */
        void * dllhandle = dlopen (lib, RTLD_LAZY);
        componentsList.components[i].libHandle = dllhandle;
        #define GETPROCADDRESS dlsym
        #define GetLastError dlerror
        #define Sleep sleep
#else
        SetLastError(0); /* clear error code */
        HINSTANCE dllhandle = LoadLibrary (LPCSTR(lib));
        componentsList.components[i].libHandle = dllhandle;
        #define GETPROCADDRESS GetProcAddress
#endif

        if (dllhandle == NULL) {

#if defined (HAVE_CONFIG_H)
            if ((err = dlerror()) != NULL)
                throw new Exception (true, "Cannot load component library \"%s\". Error: %s\n", lib, err);
#else
            if (GetLastError() == 193)
                throw new Exception (true, "Cannot load component library \"%s\". Return code: %d\n    Most probably a 32bit - 64bit conflict.", lib, GetLastError());
            else
                throw new Exception (true, "Cannot load component library \"%s\". Return code: %d", lib, GetLastError());
#endif
        }

        // Collect BMI entry points
        componentsList.components[i].dllInitialize = (BMI_INITIALIZE) GETPROCADDRESS (dllhandle, BmiInitializeEntryPoint);
        if (componentsList.components[i].dllInitialize == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiInitializeEntryPoint, lib, GetLastError());
        }

        componentsList.components[i].dllUpdate = (BMI_UPDATE) GETPROCADDRESS (dllhandle, BmiUpdateEntryPoint);
        if (componentsList.components[i].dllUpdate == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiUpdateEntryPoint, lib, GetLastError());
        }

        componentsList.components[i].dllFinalize = (BMI_FINALIZE) GETPROCADDRESS (dllhandle, BmiFinalizeEntryPoint);
        if (componentsList.components[i].dllFinalize == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiFinalizeEntryPoint, lib, GetLastError());
        }

        componentsList.components[i].dllGetStartTime = (BMI_GETSTARTTIME) GETPROCADDRESS (dllhandle, BmiGetStartTimeEntryPoint);
        if (componentsList.components[i].dllGetStartTime == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetStartTimeEntryPoint, lib, GetLastError());
        }

        componentsList.components[i].dllGetEndTime = (BMI_GETENDTIME) GETPROCADDRESS (dllhandle, BmiGetEndTimeEntryPoint);
        if (componentsList.components[i].dllGetEndTime == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetEndTimeEntryPoint, lib, GetLastError());
        }

        componentsList.components[i].dllGetTimeStep = (BMI_GETTIMESTEP) GETPROCADDRESS (dllhandle, BmiGetTimeStepEntryPoint);
        if (componentsList.components[i].dllGetStartTime == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetStartTimeEntryPoint, lib, GetLastError());
        }

        componentsList.components[i].dllGetCurrentTime = (BMI_GETCURRENTTIME) GETPROCADDRESS (dllhandle, BmiGetCurrentTimeEntryPoint);
        if (componentsList.components[i].dllGetCurrentTime == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetCurrentTimeEntryPoint, lib, GetLastError());
        }

        componentsList.components[i].dllGetAttribute = (BMI_GETATTRIBUTE) GETPROCADDRESS (dllhandle, BmiGetAttributeEntryPoint);
        if (componentsList.components[i].dllGetAttribute == NULL) {
            log->Write (Log::DETAIL, my_rank, "No GetAttribute entry point in %s !", componentsList.components[i].library);
        }    
//      If GetAttribute is optional in a lib, no need to throw an exception
//      if (componentsList.components[i].dllGetStartTime == NULL) {
//          throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetAttributeEntryPoint, lib, GetLastError());
//        }

		if (   componentsList.components[i].type == COMP_TYPE_FM
            || componentsList.components[i].type == COMP_TYPE_RTC 
            || componentsList.components[i].type == COMP_TYPE_RR 
            || componentsList.components[i].type == COMP_TYPE_FLOW1D 
            || componentsList.components[i].type == COMP_TYPE_FLOW1D2D
            || componentsList.components[i].type == COMP_TYPE_DELWAQ
            || componentsList.components[i].type == COMP_TYPE_TEST
            || componentsList.components[i].type == COMP_TYPE_WANDA) {
            // RTC-Tools: setVar is used
            componentsList.components[i].dllSetVar = (BMI_SETVAR) GETPROCADDRESS (dllhandle, BmiSetVarEntryPoint);
            if (componentsList.components[i].dllSetVar == NULL) {
                throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiSetVarEntryPoint, lib, GetLastError());
            }
        } else {
            componentsList.components[i].dllSetVar = NULL;
        }

		if (componentsList.components[i].type == COMP_TYPE_FLOW1D)
		{
			componentsList.components[i].setLogger = (BMI_SET_LOGGER)GETPROCADDRESS(dllhandle, BmiSetLogger);
			if (componentsList.components[i].setLogger == NULL) {
				throw new Exception(true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiSetLogger, lib, GetLastError());
			}
			componentsList.components[i].setLogger(&_log);
			double level = (double)DimrBmiUtils::convertDimrLogLevelToLogLevel(this->logMask);
			componentsList.components[i].dllSetVar("debugLevel", (const void *)&level);

		}

        // Not implemented yet in Delwaq:
	        if (componentsList.components[i].type != COMP_TYPE_DELWAQ) { 
            componentsList.components[i].dllGetVar = (BMI_GETVAR) GETPROCADDRESS (dllhandle, BmiGetVarEntryPoint);
            if (componentsList.components[i].dllGetVar == NULL) {
                throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetVarEntryPoint, lib, GetLastError());
            }
        }

        delete [] lib;
    }
    if (my_rank == 0) {
        printComponentVersionStrings(Log::MAJOR);            // List component version to log 
    }
}

//void Dimr::printComponentVersionStrings (Log::Mask my_mask) {
void Dimr::printComponentVersionStrings (unsigned int my_mask) {
    const char * version = "version";
    char * versionstr = new char[thisDimr->MAXSTRING];
    log->Write (my_mask, my_rank, "");
    log->Write (my_mask, my_rank, "Version Information of Components");
    log->Write (my_mask, my_rank, "=================================");
    for (int i=0;i<componentsList.numComponents;i++){
       strcpy(versionstr,"");
       if (componentsList.components[i].dllGetAttribute!=NULL){
          componentsList.components[i].dllGetAttribute(version, versionstr);
       } 
       if (strlen(versionstr)==0){
          strcpy(versionstr,"Unknown");
       }
       log->Write (my_mask, my_rank, "%-35s: %s", componentsList.components[i].library, versionstr);
    }
    log->Write (my_mask, my_rank, "---------------------------------");
    log->Write (my_mask, my_rank, "");
    delete[] versionstr;
}


//void _log(int level, char * msg) {
void _log(Level level, const char * msg) {
	if (thisDimr == NULL) {
		thisDimr = new Dimr();
	}	
	thisDimr->log->Write(DimrBmiUtils::convertBMILogLevelToDimrLogLevel(level), thisDimr->my_rank, msg);
}

//------------------------------------------------------------------------------
void Dimr::freeLibs (void) {

    //          linux windows   mac
    // lib        so    dll     dylib
    // module     so    dll     so

#if defined (OSX)
    // Macintosh:VERY SIMILAR TO LINUX
    throw new Exception (true, "ABORT: %s has not be ported to Apple Mac OS/X yet", exeName);
#endif
#if defined (HAVE_CONFIG_H)
    char *err;
#endif

    // do for all libraries
    for (int i = 0 ; i < componentsList.numComponents ; i++) {
        if (!componentsList.components[i].onThisRank) {
            continue;
        }

        log->Write (Log::DETAIL, my_rank, "Freeing library \"%s\"", componentsList.components[i].library);
#if defined (HAVE_CONFIG_H)
        dlerror(); /* clear error code */
        int ierr = dlclose(componentsList.components[i].libHandle);
        if ((err = dlerror()) != NULL) {
            throw new Exception (true, "Cannot free component library \"%s\". Error: %s\n",  componentsList.components[i].library, err);
        }
#else
        DWORD ierr;
        SetLastError(0); /* clear error code */
        bool success = FreeLibrary(componentsList.components[i].libHandle);
        if ((ierr = GetLastError()) != 0) {
            throw new Exception (true, "Cannot free component library \"%s\". Return code: %d.", componentsList.components[i].library, ierr);
        }
#endif

    }
}


//------------------------------------------------------------------------------
void Dimr::processWaitFile (void) {
    // The following waitFile code is introduced for
    // debugging parallel runs. It should NOT be used for any other purpose!

    XmlTree * rootXml     = config->Lookup ("/dimrConfig");
    const char * waitFile = rootXml->GetElement ("waitFile");
    if (waitFile != NULL) {
        log->Write (Log::MAJOR, my_rank, "Waiting for file \"%s\" to appear...", waitFile);
        fflush (stdout);
        FILE * f;
        do {
            f = fopen (waitFile, "r");
            Sleep (1000);
        } while (f == NULL);
        log->Write (Log::MAJOR, my_rank, "Found waitfile \"%s\".", waitFile);
        fclose (f);
    }
    ready = 1;
 }



//------------------------------------------------------------------------------
void Dimr::timersInit (void) {
    for (int i = 0 ; i < thisDimr->componentsList.numComponents ; i++) {
        componentsList.components[i].timerSum   =  0;
        componentsList.components[i].timerStart =  0;
    }
}



//------------------------------------------------------------------------------
void Dimr::timerStart (dimr_component * thisComponent) {
    thisComponent->timerStart = clock->Epoch();
}



//------------------------------------------------------------------------------
void Dimr::timerEnd (dimr_component * thisComponent) {
    Clock::Timestamp curtime = clock->Epoch();
    thisComponent->timerSum   =  curtime - thisComponent->timerStart + thisComponent->timerSum;
    thisComponent->timerStart =  0;
}



//------------------------------------------------------------------------------
void Dimr::timersFinish (void) {
    log->Write (Log::MAJOR, my_rank, "TIMER INFO:\n");
    for (int i = 0 ; i < thisDimr->componentsList.numComponents ; i++) {
        componentsList.components[i].timerStart = 0;
        log->Write (Log::MAJOR, my_rank, "%s\t: %d.%d sec", componentsList.components[i].name, 
                          componentsList.components[i].timerSum/1000000,
                          componentsList.components[i].timerSum%1000000);
        componentsList.components[i].timerSum   =  0.0;
    }
}




//------------------------------------------------------------------------------
void Dimr::char_to_ints(char * line, int ** iarr, int * count) {
    std::stringstream stream(line);
    // TODO: support also:
    // The processes may be specified as a space separated list with series compressed using colons
    //       e.g. 16:31

    int np = 0;
    int p;
    while (1) { // NOTE: there's no checking on a valid first number yet.
        stream >> p;
        if (!stream) {
            break;
        }
        np++;
    }
    *count = np;
    *iarr = (int *)malloc(*count * sizeof(int));

    stream.clear();
    stream.seekg(0);
    stream.str(line);
    np = 0;
    while (1) { // NOTE: there's no checking on a valid first number yet.
        stream >> (*iarr)[np];
        if (!stream) {
            break;
        }
        np++;
    }
}

