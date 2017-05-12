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
#  define strdup _strdup
#  define chdir _chdir
#  define getcwd _getcwd
#else
#  include <unistd.h>
#endif

// static added to prevent name conflicts on Linux.
static Dimr * thisDimr = NULL;     // global pointer to single object instance

Dimr::Dimr(void) {
    FILE *      logFile = stdout;
    this->ready = false;
    this->exePath = NULL;
    this->exeName = NULL;
    this->clock = new Clock ();
    this->logMask = Log::MAJOR;
    this->log = new Log (logFile, this->clock, this->logMask);
    this->config = NULL;
    this->mainArgs = NULL;
    this->slaveArg = NULL;
    this->control = NULL;
    this->componentsList.numComponents = 0;
    this->couplersList.numCouplers = 0;
    this->use_mpi = false;
    this->my_rank = 0;
    this->numranks = 1;
    this->configfile = NULL;
    this->done = false;
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
DllExport void set_logger(Log * loggerFromDimrExe) {
	if (thisDimr == NULL) {
		thisDimr = new Dimr();
	}
	thisDimr->log = loggerFromDimrExe;
}


//------------------------------------------------------------------------------
DllExport int initialize(const char * configfile) {
	int nSettingsSet, nParamsSet;
    if (thisDimr == NULL) {
        thisDimr = new Dimr();
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
DllExport void update (double tStep) {
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
}

//------------------------------------------------------------------------------
DllExport void finalize (void) {
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
}

//------------------------------------------------------------------------------
DllExport void get_start_time (double * tStart) {
    thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr_lib:get_start_time");
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        *tStart = thisDimr->control->subBlocks[0].subBlocks[thisDimr->control->subBlocks[0].masterSubBlockId].tStart;
    } else {
        // Start block
        *tStart = thisDimr->control->subBlocks[0].tStart;
    }
}

//------------------------------------------------------------------------------
DllExport void get_end_time (double * tEnd) {
    thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr_lib:get_end_time");
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        *tEnd = thisDimr->control->subBlocks[0].subBlocks[thisDimr->control->subBlocks[0].masterSubBlockId].tEnd;
    } else {
        // Start block
        *tEnd = thisDimr->control->subBlocks[0].tEnd;
    }
}

//------------------------------------------------------------------------------
DllExport void get_time_step (double * tStep) {
    thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr_lib:get_time_step");
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        *tStep = thisDimr->control->subBlocks[0].subBlocks[thisDimr->control->subBlocks[0].masterSubBlockId].tStep;
    } else {
        // Start block
        *tStep = thisDimr->control->subBlocks[0].tStep;
    }
}

//------------------------------------------------------------------------------
DllExport void get_current_time (double * tCur) {
    thisDimr->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr_lib:get_current_time");
    if (thisDimr->control->subBlocks[0].type == CT_PARALLEL) {
        *tCur = thisDimr->control->subBlocks[0].subBlocks[thisDimr->control->subBlocks[0].masterSubBlockId].tCur;
    } else {
        // Start block
        *tCur = thisDimr->control->subBlocks[0].tCur;
    }
}

//------------------------------------------------------------------------------
DllExport void get_var (const char * key, void ** ref) {
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
    *ref = thisDimr->send (sourceName, 
                           compPtr->type,
                           compPtr->dllGetVar, 
                           &sourceVarPtr,
                           compPtr->processes,
                           compPtr->numProcesses,
                           sourceProcess);
    delete[] componentName;
}

//------------------------------------------------------------------------------
DllExport void set_var (const char * key, const void * value) {
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
			if (this->dllSetVar!=NULL){
		       (this->dllSetVar) (kv->key,(void*)kv->val);
			} else {
			   if (this->dllGetVar!=NULL){
		          (this->dllGetVar) (kv->key,(void**)kv->val);
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

    if (this->done)
        return;

    // to do:  (void) FreeLibrary(handle);
    freeLibs();

    this->log->Write (Log::MAJOR, thisDimr->my_rank, "dimr shutting down normally");

#if defined(HAVE_CONFIG_H)
    free (this->exeName);
#else
    delete [] this->exeName;
#endif

    delete this->clock;
    delete this->config;
    free (this->exePath);
    delete this->log;
    delete [] this->mainArgs;
    // componentsList
    delete [] this->componentsList.components;
    // couplersList
    for (int i = 0 ; i < this->couplersList.numCouplers ; i++) {
        delete [] this->couplersList.couplers[i].items;
    }
    delete [] this->couplersList.couplers;
    // control
    deleteControlBlock(*(this->control));
    free(this->control);
    this->done = true;
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
        this->log->Write (Log::MAJOR, my_rank, "PARALLEL:");
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
    this->log->Write (Log::MAJOR, my_rank, "START:");

    chdir(cb->unit.component->workingDir);
    if (phase == GLOBAL_PHASE_FINISH) {
        this->log->Write (Log::MAJOR, my_rank, "%s.Initialize(%s)", cb->unit.component->name, cb->unit.component->inputFile);
        this->timerStart(cb->unit.component);
        cb->unit.component->result = (cb->unit.component->dllInitialize) (cb->unit.component->inputFile);
        this->timerEnd(cb->unit.component);
        (cb->unit.component->dllGetStartTime) (&cb->tStart);
        (cb->unit.component->dllGetEndTime) (&cb->tEnd);
    }
    cb->tStep = tStep;
    this->log->Write (Log::MAJOR, my_rank, "%s.Update(%6.1f)", cb->unit.component->name, cb->tStep);
    this->timerStart(cb->unit.component);
    (cb->unit.component->dllUpdate) (cb->tStep);
    this->timerEnd(cb->unit.component);
    if (phase == GLOBAL_PHASE_FINISH) {
        this->log->Write (Log::MAJOR, my_rank, "%s.Finalize()", cb->unit.component->name);
        this->timerStart(cb->unit.component);
        (cb->unit.component->dllFinalize) ();
        this->timerEnd(cb->unit.component);
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
                this->log->Write (Log::MINOR, my_rank, "Master: %s", cb->subBlocks[cb->masterSubBlockId].unit.component->name);
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
        this->log->Write (Log::MAJOR, my_rank, "%s.Initialize(%s)", masterComponent->name, masterComponent->inputFile);
        this->timerStart(masterComponent);
        masterComponent->result = (masterComponent->dllInitialize) (masterComponent->inputFile);
        this->timerEnd(masterComponent);
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
                        this->log->Write (Log::MAJOR, my_rank, "%s.Initialize(%s)", thisComponent->name, thisComponent->inputFile);
                        this->timerStart(thisComponent);
                        thisComponent->result = (thisComponent->dllInitialize) (thisComponent->inputFile);
                        this->timerEnd(thisComponent);
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
                                    this->log->Write (Log::MINOR, my_rank, "%s.getVar(%s)", thisCoupler->sourceComponentName, thisCoupler->items[k].sourceName);
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
                                            this->log->Write(Log::WARN, my_rank, "WARNING: coupler %s: item %d: \"%s\" will be delivered by partition %d. Ignoring deliverance by partition %d",
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
                                        this->log->Write (Log::MINOR, my_rank, "%s.getVar(%s)", thisCoupler->targetComponentName, thisCoupler->items[k].targetName);
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
                                            this->log->Write(Log::WARN, my_rank, "WARNING: coupler %s: item %d: \"%s\" will be delivered to multiple partitions: %d and %d.",
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
                this->log->Write (Log::MAJOR, my_rank, "%10.1f:    %s.Update(%10.1f)", *currentTime, masterComponent->unit.component->name, tStep);
                this->timerStart(masterComponent->unit.component);
                (masterComponent->unit.component->dllUpdate) (tStep);
                this->timerEnd(masterComponent->unit.component);
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
                            this->log->Write (Log::MAJOR, my_rank, "%10.1f:    %s.Update(%10.1f)", *currentTime, thisComponent->name, tUpdate);
                            this->timerStart(thisComponent);
                            (thisComponent->dllUpdate) (tUpdate);
                            this->timerEnd(thisComponent);
                        } else {
                            // Coupler
                            dimr_coupler * thisCoupler = cb->subBlocks[i].subBlocks[j].unit.coupler;
                            double * transferValuePtr;
                            this->log->Write (Log::MINOR, my_rank, "%10.1f:    %s.communicate", *currentTime, thisCoupler->name);
                            for (int k = 0 ; k < thisCoupler->numItems ; k++) {
                                this->log->Write (Log::DETAIL, my_rank, "    %s -> %s", thisCoupler->items[k].sourceName, thisCoupler->items[k].targetName);
                                // Getting and Setting of data is split to enable
                                // transferring data, possibly inbetween different partitions
                                // TODO: This does not work for arrays (yet), only scalar double
                                //
                                // Getting data:
                                transferValuePtr = this->send (thisCoupler->items[k].sourceName, 
                                                         thisCoupler->sourceComponent->type,
                                                         thisCoupler->sourceComponent->dllGetVar, 
                                                         &(thisCoupler->items[k].sourceVarPtr),
                                                         thisCoupler->sourceComponent->processes,
                                                         thisCoupler->sourceComponent->numProcesses,
                                                         thisCoupler->items[k].sourceProcess);
                                // Optional TODO: assuming one source(partition):
                                //    if there is only one target (partition), and all partitions can act as source partition:
                                //        choose the target partition to act as source partition
                                //        no MPI_Bcast needed
                                //
                                this->receive (thisCoupler->items[k].targetName,
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
double * Dimr::send (const char * name, 
                     int          compType, 
                     BMI_GETVAR   dllGetVar, 
                     double    ** sourceVarPtr,    // Output parameter! Needed when doing a dllSetVar via "get_var, var=value". See "send" call in Dimr::receive
                     int        * processes,
                     int          nProc,
                     int          sourceProcess) {
    // Array "transfer" is used to collect a scalar double from all partitions
    // For now: The resulting scalar double is stored in the class parameter "transferValeu"
    // TODO: Point to variable in related kernel? How when running in parallel?
    double * transfer       = (double *)malloc(nProc * sizeof(double));   
    double * reducedTransfer = (double *)malloc(nProc * sizeof(double));
    // put data to be transferred from this partitions source location into the transfer array
    for (int m = 0 ; m < nProc; m++) 
	{
		// source provider component not available on this rank.
		// Define a missing value (large negative number for later reduction operations)
		transfer[m] = -999000.0;

        if (my_rank == processes[m]) 
		{
            this->log->Write (Log::DETAIL, my_rank, "Send(%s)", name);
			if (    compType == COMP_TYPE_RTC
				|| compType == COMP_TYPE_RR
				|| compType == COMP_TYPE_FLOW1D
				|| compType == COMP_TYPE_FLOW1D2D
                || *sourceVarPtr == NULL          ) {
                // These components only returns a new pointer to a copy of the double value, so call it each time.
                // sourceVarPtr=NULL: getVar not yet called for this parameter, probably because "send" is being called
                //                    via the toplevel "get_var"
                if (dllGetVar == NULL) {
                    throw new Exception (true, "ABORT: get_var function not defined while processing %s", name);
                }
                (dllGetVar) (name, (void**)(sourceVarPtr));
            } else if (compType == COMP_TYPE_WANDA) {
                if (dllGetVar == NULL) {
                    throw new Exception (true, "ABORT: get_var function not defined while processing %s", name);
                }
                // Wanda does not use pointers to internal structures:
                // - Use the DIMR-transfer array
                // - Note the missing & in the dllGetVar call, when comparing with the dllGetVar call above
                *sourceVarPtr = transfer;
                (dllGetVar) (name, (void**)(*sourceVarPtr));
            } else {
                // Other components already have direct pointer access to the actual variable.
            }
            // The sourceVarPtr may be NULL for my_rank
            // Wanda: the value is already in the transfer array
            if (*sourceVarPtr != NULL && compType != COMP_TYPE_WANDA)
			{
                transfer[m] = *(*sourceVarPtr);
            }
        } 
    }

	double maxValue = -999000.0;
	// Do not call MPI_Allreduce when the number of partitions is 1. Big chance that it will cause a crash
	if (numranks > 1) 
	{
		// NOTE: here the transfer array has a defined value only at position==my_rank and if *sourceVarPtr != NULL.
		// We perform a reduction operation to collect the maximum values at each position and afterwards take the maximum of all values.
		// We could also use only one double.instead of a transfer array and perform the reduction on a single scalar (so avoiding the loop below) 
		int ierr = MPI_Allreduce(transfer, reducedTransfer, nProc, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		
		for (int m = 0; m < nProc; m++)
		{
			if (reducedTransfer[m]> maxValue) maxValue = reducedTransfer[m];
		}
    }
	else
	{
		maxValue = transfer[0];
	}
	
	transferValue = maxValue;
    free(transfer);
	free(reducedTransfer);
	return &(transferValue);
}



//------------------------------------------------------------------------------
// set "value" in the comopnents target location
void Dimr::receive (const char * name,
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
        for (int m = 0 ; m < nProc; m++) {
            if (my_rank == processes[m]) {
                this->log->Write (Log::DETAIL, my_rank, "Dimr::receive(%s) %20.10f", name, *(double*)(transferValuePtr));
				if (   compType == COMP_TYPE_RTC 
                    || compType == COMP_TYPE_RR 
                    || compType == COMP_TYPE_FLOW1D 
					|| compType == COMP_TYPE_FLOW1D2D
                    || compType == COMP_TYPE_WANDA   ) {
                    if (dllSetVar == NULL) {
                        throw new Exception (true, "ABORT: Dimr::receive: set_var function not defined while processing %s", name);
                    }
                    (dllSetVar) (name,(const void *)transferValuePtr);
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
						(dllSetVar) (targetName, (const void *)transferValuePtr);
						delete[] targetName;
					}
                } else {
                    // target is a component that uses direct pointer access to the actual variable
                    // When doing a dllSetVar, targetVarPtr is not defined. First do a "send" to get it defined.
                    // WARNING: the test "if (targetVarPtr == NULL)" is not allowed:
                    //          When target is defined for 1 partition and undefined for the other partitions,
                    //          send will be called for all but 1 partitions. This screws the MPI administration,
                    //          resulting in a crash.
                    //          If the test "if (targetVarPtr == NULL)" is prefered, a more detailed administration is needed.
                    //   if (targetVarPtr == NULL) {

					// transferValuePtr contains the value calculated by RTCTools. 
					// Here we already know transferValuePtr is not null.
					double * valueToTransfer = ( double *)transferValuePtr;
					thisDimr->send(name,
                                    compType,
					                dllGetVar, 
									&valueToTransfer,
					                processes,
					                nProc,
					                0);
                    //  }
                    if (targetVarPtr == NULL) 
					   {
                        if (targetProcess == -1 || targetProcess == my_rank) {
                            // targetProcess=-1: no process can accept this item
                            // targetProcess=my_rank: this process is registered to be able to accept this item but something goes wrong
                            throw new Exception (true, "ABORT: Dimr::receive: get_var function not defined while processing %s", name);
                        }
                    } 
					else 
					{
						// Here writes the RTC value in FM (direct access to memory using a pointer). 
						// Here we already know targetVarPtr is not null. 
						*(targetVarPtr) = transferValue;
                    }
                }
            }
        }
    }
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
            this->log->Write (Log::MAJOR, my_rank, "    %s.Finalize()", cb->subBlocks[i].unit.component->name);
            this->timerStart(cb->subBlocks[i].unit.component);
            (cb->subBlocks[i].unit.component->dllFinalize) ();
            this->timerEnd(cb->subBlocks[i].unit.component);

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
                    this->log->Write (Log::MAJOR, my_rank, "    %s.Finalize()", cb->subBlocks[i].subBlocks[j].unit.component->name);
                    this->timerStart(cb->subBlocks[i].subBlocks[j].unit.component);
                    (cb->subBlocks[i].subBlocks[j].unit.component->dllFinalize) ();
                    this->timerEnd(cb->subBlocks[i].subBlocks[j].unit.component);
                    }
            }
        }
    }
}



//------------------------------------------------------------------------------
void Dimr::scanConfigFile (void) {

	XmlTree * rootXml     = this->config->Lookup ("/dimrConfig");
    if (rootXml == NULL)
        throw new Exception (true, "Configuration file \"%s\" does not have a <dimrConfig> root element", this->configfile);
    XmlTree * fileversion = rootXml->Lookup ("documentation/fileVersion");
    if (fileversion == NULL)
        throw new Exception (true, "Configuration file \"%s\" does not have a deltaresHydro documentation->fileVersion element", this->configfile);

	// Check version number
    const char * versionnr = fileversion->charData;
    float versionnumber;
    int intRead = sscanf(versionnr, "%f", &versionnumber);
    if (intRead != 1)
        throw new Exception (true, "Configuration file \"%s\" does not have a version number", this->configfile);
    if ((int)floor(versionnumber) != 1)
        throw new Exception (true, "Configuration file \"%s\": Version number (%3.2f) must have main version 2", this->configfile,versionnumber);

	XmlTree * controlXml  = rootXml->Lookup ("control");
    if (controlXml == NULL)
        throw new Exception (true, "Configuration file \"%s\" does not have a deltaresHydro control element", this->configfile);
    // Allocate the control structure and check its size
    this->control = (dimr_control_block *) malloc(sizeof(dimr_control_block));
    this->control->numSubBlocks = 0;
    this->control->subBlocks    = NULL;
    this->control->masterSubBlockId = -1;
    // First scan the config file for all components and couplers (= units)
    scanUnits(rootXml);
    // Then scan the control part
    // References are added to the list of components/couplers
    scanControl(controlXml, this->control);
}



//------------------------------------------------------------------------------
void Dimr::scanUnits(XmlTree * rootXml) {
    // Init
    this->componentsList.numComponents = 0;
    this->componentsList.components    = NULL;
    this->couplersList.numCouplers     = 0;
    this->couplersList.couplers        = NULL;
    // Scan
    for (int i = 0 ; i < rootXml->numChildren ; i++) {
        if (strcmp(rootXml->children[i]->name, "component") == 0) {
            this->componentsList.numComponents++;
            if (this->componentsList.components == NULL) {
                this->componentsList.components = (dimr_component*)malloc(this->componentsList.numComponents * sizeof(dimr_component));
            } else {
                this->componentsList.components = (dimr_component*)realloc(this->componentsList.components, 
                                                                          this->componentsList.numComponents * sizeof(dimr_component));
                if (this->componentsList.components == NULL) {
                    throw new Exception (true, "Allocation error in scanUnits (component)");
                }
            }
            scanComponent(rootXml->children[i], &(this->componentsList.components[this->componentsList.numComponents - 1]));
        }
        if (strcmp(rootXml->children[i]->name, "coupler") == 0) {
            this->couplersList.numCouplers++;
            if (this->couplersList.couplers == NULL) {
                this->couplersList.couplers = (dimr_coupler*)malloc(this->couplersList.numCouplers * sizeof(dimr_coupler));
            } else {
                this->couplersList.couplers = (dimr_coupler*)realloc(this->couplersList.couplers, this->couplersList.numCouplers * sizeof(dimr_coupler));
                if (this->couplersList.couplers == NULL) {
                    throw new Exception (true, "Allocation error in scanUnits (coupler)");
                }
            }
            scanCoupler(rootXml->children[i], &(this->couplersList.couplers[this->couplersList.numCouplers - 1]));
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
    char curPath[FILENAME_MAX];
    if (!getcwd(curPath, sizeof(curPath)))
        throw new Exception (true, "ERROR obtaining the current working directory");
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
        this->log->Write (Log::MAJOR, my_rank, "WARNING: \"<process>\" not specified for component \"%s\". Assuming it only runs on rank #0.", newComp->name);
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
        this->log->Write (Log::MAJOR, my_rank, "WARNING: No inputFile specified for component %s.", newComp->name);
        }
    else {
        newComp->inputFile = inputFileElement->charData;
        }
    // Element workingDir
    XmlTree * workingDirElement = xmlComponent->Lookup ("workingDir");
    if (workingDirElement == NULL) {
        newComp->workingDir = (char *) &curPath;
        this->log->Write (Log::MAJOR, my_rank, "WARNING: No workingDir specified for component %s.", newComp->name);
        this->log->Write (Log::MAJOR, my_rank, "         workingDir is set to %s", newComp->workingDir);
        }
    else {
        newComp->workingDir = workingDirElement->charData;
        }
    // Is workingDir a valid relative path?
    char *combinedPath= new char[FILENAME_MAX];
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
    newCoup->sourceComponent = this->getComponent(newCoup->sourceComponentName);
    // Element targetComponent
    XmlTree * targetComponent = xmlCoupler->Lookup ("targetComponent");
    if (targetComponent == NULL)
        throw new Exception (true, "The coupler \"%s\" does not contain a targetComponent element", newCoup->name);
    newCoup->targetComponentName = targetComponent->charData;
    // Add reference to the actual component acting as target
    newCoup->targetComponent = this->getComponent(newCoup->targetComponentName);
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
        controlBlock->unit.component = this->getComponent(controlBlockXml->GetAttrib("name"));
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
        controlBlock->unit.coupler = this->getCoupler(controlBlockXml->GetAttrib("name"));
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
            this->scanControl(controlBlockXml->children[i], &(controlBlock->subBlocks[controlBlock->numSubBlocks - 1]));
        }
        fflush(stdout);
    }
}



//------------------------------------------------------------------------------
// Search for a named component in the list of components
dimr_component * Dimr::getComponent(const char * compName) {
    for (int i = 0 ; i < this->componentsList.numComponents ; i++) {
        if (strcmp(this->componentsList.components[i].name, compName) == 0) {
            return &(this->componentsList.components[i]);
        }
    }
}



//------------------------------------------------------------------------------
// Search for a named coupler in the list of couplers
dimr_coupler * Dimr::getCoupler(const char * coupName) {
    for (int i = 0 ; i < this->couplersList.numCouplers ; i++) {
        if (strcmp(this->couplersList.couplers[i].name, coupName) == 0) {
            return &(this->couplersList.couplers[i]);
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
    throw new Exception (true, "ABORT: %s has not be ported to Apple Mac OS/X yet", this->exeName);
#endif
#if defined (HAVE_CONFIG_H)
    char *err;
#endif

    // do for all libraries
    for (int i = 0 ; i < this->componentsList.numComponents ; i++) {
        if (!this->componentsList.components[i].onThisRank) {
            this->log->Write (Log::DETAIL, my_rank, "Not necessary to load component library \"%s\" on this rank.", this->componentsList.components[i].library);
            continue;
        }

#if defined (HAVE_CONFIG_H)
        char * lib = new char[strlen (this->componentsList.components[i].library) + 3+3+1];
        sprintf (lib, "lib%s%s", this->componentsList.components[i].library, D3D_PLUGIN_EXT);
        if (   strchr (this->componentsList.components[i].library, '/' ) != NULL 
            || strchr (this->componentsList.components[i].library, '\\') != NULL 
            || strchr (this->componentsList.components[i].library, '.' ) != NULL) {
            throw new Exception (true, "Invalid component library name \"%s\"\n", lib, -1);
        }
#else
        char * lib = new char[strlen (this->componentsList.components[i].library) + 4+1];
        sprintf (lib, "%s.dll", this->componentsList.components[i].library);
        if (   strchr (this->componentsList.components[i].library, '/' ) != NULL 
            || strchr (this->componentsList.components[i].library, '\\') != NULL 
            || strchr (this->componentsList.components[i].library, '.' ) != NULL) {
            throw new Exception (true, "Invalid component library name \"%s\"\n", lib, -1);
        }
#endif

        this->log->Write (Log::DETAIL, my_rank, "Loading library \"%s\"", lib);

#if defined (HAVE_CONFIG_H)
        dlerror(); /* clear error code */
        void * dllhandle = dlopen (lib, RTLD_LAZY);
        this->componentsList.components[i].libHandle = dllhandle;
        #define GETPROCADDRESS dlsym
        #define GetLastError dlerror
        #define Sleep sleep
#else
        SetLastError(0); /* clear error code */
        HINSTANCE dllhandle = LoadLibrary (LPCSTR(lib));
        this->componentsList.components[i].libHandle = dllhandle;
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
        this->componentsList.components[i].dllInitialize = (BMI_INITIALIZE) GETPROCADDRESS (dllhandle, BmiInitializeEntryPoint);
        if (this->componentsList.components[i].dllInitialize == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiInitializeEntryPoint, lib, GetLastError());
        }

        this->componentsList.components[i].dllUpdate = (BMI_UPDATE) GETPROCADDRESS (dllhandle, BmiUpdateEntryPoint);
        if (this->componentsList.components[i].dllUpdate == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiUpdateEntryPoint, lib, GetLastError());
        }

        this->componentsList.components[i].dllFinalize = (BMI_FINALIZE) GETPROCADDRESS (dllhandle, BmiFinalizeEntryPoint);
        if (this->componentsList.components[i].dllFinalize == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiFinalizeEntryPoint, lib, GetLastError());
        }

        this->componentsList.components[i].dllGetStartTime = (BMI_GETSTARTTIME) GETPROCADDRESS (dllhandle, BmiGetStartTimeEntryPoint);
        if (this->componentsList.components[i].dllGetStartTime == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetStartTimeEntryPoint, lib, GetLastError());
        }

        this->componentsList.components[i].dllGetEndTime = (BMI_GETENDTIME) GETPROCADDRESS (dllhandle, BmiGetEndTimeEntryPoint);
        if (this->componentsList.components[i].dllGetEndTime == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetEndTimeEntryPoint, lib, GetLastError());
        }

        this->componentsList.components[i].dllGetTimeStep = (BMI_GETTIMESTEP) GETPROCADDRESS (dllhandle, BmiGetTimeStepEntryPoint);
        if (this->componentsList.components[i].dllGetStartTime == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetStartTimeEntryPoint, lib, GetLastError());
        }

        this->componentsList.components[i].dllGetCurrentTime = (BMI_GETCURRENTTIME) GETPROCADDRESS (dllhandle, BmiGetCurrentTimeEntryPoint);
        if (this->componentsList.components[i].dllGetCurrentTime == NULL) {
            throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetCurrentTimeEntryPoint, lib, GetLastError());
        }

        this->componentsList.components[i].dllGetAttribute = (BMI_GETATTRIBUTE) GETPROCADDRESS (dllhandle, BmiGetAttributeEntryPoint);
        if (this->componentsList.components[i].dllGetAttribute == NULL) {
            this->log->Write (Log::DETAIL, my_rank, "No GetAttribute entry point in %s !", this->componentsList.components[i].library);
		}	
//      If GetAttribute is optional in a lib, no need to throw an exception
//      if (this->componentsList.components[i].dllGetStartTime == NULL) {
//          throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetAttributeEntryPoint, lib, GetLastError());
//        }

        if (   this->componentsList.components[i].type == COMP_TYPE_FM
            || this->componentsList.components[i].type == COMP_TYPE_RTC 
            || this->componentsList.components[i].type == COMP_TYPE_RR 
            || this->componentsList.components[i].type == COMP_TYPE_FLOW1D 
            || this->componentsList.components[i].type == COMP_TYPE_FLOW1D2D
            || this->componentsList.components[i].type == COMP_TYPE_DELWAQ
            || this->componentsList.components[i].type == COMP_TYPE_TEST
            || this->componentsList.components[i].type == COMP_TYPE_WANDA) {
            // RTC-Tools: setVar is used
            this->componentsList.components[i].dllSetVar = (BMI_SETVAR) GETPROCADDRESS (dllhandle, BmiSetVarEntryPoint);
            if (this->componentsList.components[i].dllSetVar == NULL) {
                throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiSetVarEntryPoint, lib, GetLastError());
            }
        } else {
            this->componentsList.components[i].dllSetVar = NULL;
		}

        // Not implemented yet in Delwaq:
        if (this->componentsList.components[i].type != COMP_TYPE_DELWAQ) { 
            this->componentsList.components[i].dllGetVar = (BMI_GETVAR) GETPROCADDRESS (dllhandle, BmiGetVarEntryPoint);
            if (this->componentsList.components[i].dllGetVar == NULL) {
                throw new Exception (true, "Cannot find function \"%s\" in library \"%s\". Return code: %d", BmiGetVarEntryPoint, lib, GetLastError());
            }
        }

        delete [] lib;
    }
    if (my_rank == 0) {
        this->printComponentVersionStrings(Log::MAJOR);			// List component version to log 
    }
}

//void Dimr::printComponentVersionStrings (Log::Mask my_mask) {
void Dimr::printComponentVersionStrings (unsigned int my_mask) {
    const char * version = "version";
	char * versionstr = new char[thisDimr->MAXSTRING];
    this->log->Write (my_mask, my_rank, "");
    this->log->Write (my_mask, my_rank, "Version Information of Components");
    this->log->Write (my_mask, my_rank, "=================================");
	for (int i=0;i<this->componentsList.numComponents;i++){
	   strcpy(versionstr,"");
	   if (this->componentsList.components[i].dllGetAttribute!=NULL){
          this->componentsList.components[i].dllGetAttribute(version, versionstr);
	   } 
	   if (strlen(versionstr)==0){
	      strcpy(versionstr,"Unknown");
	   }
       this->log->Write (my_mask, my_rank, "%-35s: %s", this->componentsList.components[i].library, versionstr);
	}
    this->log->Write (my_mask, my_rank, "---------------------------------");
    this->log->Write (my_mask, my_rank, "");
    delete[] versionstr;
}

//------------------------------------------------------------------------------
void Dimr::freeLibs (void) {

    //          linux windows   mac
    // lib        so    dll     dylib
    // module     so    dll     so

#if defined (OSX)
    // Macintosh:VERY SIMILAR TO LINUX
    throw new Exception (true, "ABORT: %s has not be ported to Apple Mac OS/X yet", this->exeName);
#endif
#if defined (HAVE_CONFIG_H)
    char *err;
#endif

    // do for all libraries
    for (int i = 0 ; i < this->componentsList.numComponents ; i++) {
        if (!this->componentsList.components[i].onThisRank) {
            continue;
        }

        this->log->Write (Log::DETAIL, my_rank, "Freeing library \"%s\"", this->componentsList.components[i].library);
#if defined (HAVE_CONFIG_H)
        dlerror(); /* clear error code */
        int ierr = dlclose(this->componentsList.components[i].libHandle);
        if ((err = dlerror()) != NULL) {
            throw new Exception (true, "Cannot free component library \"%s\". Error: %s\n",  this->componentsList.components[i].library, err);
        }
#else
        DWORD ierr;
        SetLastError(0); /* clear error code */
        bool success = FreeLibrary(this->componentsList.components[i].libHandle);
        if ((ierr = GetLastError()) != 0) {
            throw new Exception (true, "Cannot free component library \"%s\". Return code: %d.", this->componentsList.components[i].library, ierr);
        }
#endif

    }
}


//------------------------------------------------------------------------------
void Dimr::processWaitFile (void) {
    // The following waitFile code is introduced for
    // debugging parallel runs. It should NOT be used for any other purpose!

    XmlTree * rootXml     = this->config->Lookup ("/dimrConfig");
    const char * waitFile = rootXml->GetElement ("waitFile");
    if (waitFile != NULL) {
        this->log->Write (Log::MAJOR, my_rank, "Waiting for file \"%s\" to appear...", waitFile);
        fflush (stdout);
        FILE * f;
        do {
            f = fopen (waitFile, "r");
            Sleep (1000);
        } while (f == NULL);
        this->log->Write (Log::MAJOR, my_rank, "Found waitfile \"%s\".", waitFile);
        fclose (f);
    }
    this->ready = 1;
 }



//------------------------------------------------------------------------------
void Dimr::timersInit (void) {
    for (int i = 0 ; i < thisDimr->componentsList.numComponents ; i++) {
        this->componentsList.components[i].timerSum   =  0;
        this->componentsList.components[i].timerStart =  0;
    }
}



//------------------------------------------------------------------------------
void Dimr::timerStart (dimr_component * thisComponent) {
    thisComponent->timerStart = this->clock->Epoch();
}



//------------------------------------------------------------------------------
void Dimr::timerEnd (dimr_component * thisComponent) {
    Clock::Timestamp curtime = this->clock->Epoch();
    thisComponent->timerSum   =  curtime - thisComponent->timerStart + thisComponent->timerSum;
    thisComponent->timerStart =  0;
}



//------------------------------------------------------------------------------
void Dimr::timersFinish (void) {
    this->log->Write (Log::MAJOR, my_rank, "TIMER INFO:\n");
    for (int i = 0 ; i < thisDimr->componentsList.numComponents ; i++) {
        this->componentsList.components[i].timerStart = 0;
        this->log->Write (Log::MAJOR, my_rank, "%s\t: %d.%d sec", this->componentsList.components[i].name, 
                          this->componentsList.components[i].timerSum/1000000,
                          this->componentsList.components[i].timerSum%1000000);
        this->componentsList.components[i].timerSum   =  0.0;
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

