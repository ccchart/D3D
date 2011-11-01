//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011.
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
// $Id: flowol.cpp 962 2011-10-31 21:52:47Z elshoff $
// $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/branches/research/Deltares/20110420_OnlineVisualisation/src/engines_gpl/flow2d3d/packages/flow2d3d/src/flowol.cpp $
//-------------------------------------------------------------------------------
//  d_hydro Flow2FLOW Component
//  Interface to DelftOnline - IMPLEMENTATION
//
//  Irv.Elshoff@deltares.NL
//  6 jun 11
//-------------------------------------------------------------------------------


// #define WITH_FLOWOL



#ifdef WIN32
#   define strcasecmp _stricmp
#endif


#include "flow2d3d.h"
#include "flowol.h"


//-------------------------------------------------------------------------------


FlowOL::FlowOL (
    DeltaresHydro * dh,
    XmlTree * config
    ) {

#if defined (WITH_FLOWOL)
    this->dh = dh;
    this->numSubdomains = 0;

    // Get DOL options from the configuration tree

    char * jrePath = config->GetAttrib ("JREPath");
    if (jrePath == NULL)
        throw new Exception (true, "JREPath not specified in DelftOnline section of the configuration file");

    char * jarPath = config->GetAttrib ("JARPath");
    if (jarPath == NULL)
        throw new Exception (true, "JARPath not specified in DelftOnline section of the configuration file");

    char * urlFile = config->GetAttrib ("URLFile");
    if (urlFile == NULL)
        throw new Exception (true, "URLFile not specified in DelftOnline section of the configuration file");

    char * jvmOpts = config->GetAttrib ("JVMOpts");
    bool allowStart = ! config->GetBoolAttrib ("wait");
    bool allowControl = config->GetBoolAttrib ("control");
    allowControl |= ! allowStart;

    DOL::Verbosity verbosity = DOL::SILENT;
    char * verb = config->GetAttrib ("verbosity");
    if (verb != NULL) {
        if (strcasecmp (verb, "silent") == 0)
            verbosity = DOL::SILENT;
        else if (strcasecmp (verb, "error") == 0)
            verbosity = DOL::ERROR;
        else if (strcasecmp (verb, "info") == 0)
            verbosity = DOL::INFO;
        else if (strcasecmp (verb, "trace") == 0)
            verbosity = DOL::TRACE;
        else
            throw new Exception (true, "Invalid verbosity value \"%s\" in DelftOnline section of the configuration file", verb);
        }

    // Initialize the Java Virtual Machine

    try {
        this->dh->log->Write (Log::MINOR, "Starting Java virtual machine");
        this->java = new JavaLaunch (jrePath, jarPath, jvmOpts);
        }
    catch (JavaLaunch::Exception * ex) {
        throw new Exception (true, "Java initialization fails: %s", ex->message);
        }

    // Initialize DelftOnline and write handle to URL file

    try {
        this->dh->log->Write (Log::MINOR, "Creating DOL server");
        this->dol = new DOL::Server (allowStart, true, verbosity, NULL, this->java);

        char * url = this->dol->Handle ();
        this->dh->log->Write (Log::MAJOR, "DOL handle is \"%s\"", url);

        FILE * f = fopen (urlFile, "w");
        if (f == NULL)
            throw new Exception ("Cannot create URL file \"%s\": %s", urlFile, strerror (errno));

        fprintf (f, "'%s'\n", url);
        fclose (f);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "DOL initialization fails: %s", ex->message);
        }
#endif
    }


FlowOL::~FlowOL (
    void
    ) {

#if defined (WITH_FLOWOL)
    this->dh->log->Write (Log::MINOR, "Shutting down DOL");

    delete this->dol;
    delete this->java;
#endif
    }


//-------------------------------------------------------------------------------


void
FlowOL::RegisterSubdomain (
    const char * name
    ) {

#if defined (WITH_FLOWOL)
    if (this->dol == NULL) return;

    try {
        this->dol->RegisterThread  (this->numSubdomains, name);
#ifdef WIN32
        this->dol->CreateDirectoryA (name);
#else
        this->dol->CreateDirectory (name);
#endif
        this->dol->ChangeDirectory (name);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "DOL subdomain registration fails: %s", ex->message);
        }
#endif
    }


void
FlowOL::UnregisterSubdomain (
    void
    ) {

    // No-op
    }


//-------------------------------------------------------------------------------
//  Fortran-to-C++ string conversion function


static char *
createCString (
    const char *    fortString,
    int             fortLen
    ) {

    // Look for the last non-space character in the Fortran string (padded with spaces)

    int end;
    for (end = fortLen-1 ; end >= 0 ; end--)
        if (fortString[end] != ' ')
            break;

    if (end < 0) end = 0;       // all spaces becomes empty string

    // Allocate a C buffer for the truncated string

    end++;
    char * cString = new char [end + 1 /* for EOL */];

    // Copy Fortran to C string

    for (int i = 0 ; i < end ; i++)
        cString[i] = fortString[i];

    cString[end] = '\0';

    return cString;
    }


//-------------------------------------------------------------------------------
//  Publication and timestep functions (called from Fortran)


extern "C" {    // closing brace at end of file

void STDCALL
FLOWOL_ArrayShape (
#if defined(MIXED_STR_LEN_ARG)
    const char *    ftn_name,
    int             len_name,
    int *           dimensionality,
    int             dimensions []
#else
    const char *    ftn_name,
    int *           dimensionality,
    int             dimensions [],
    int             len_name
#endif
    ) {

#if defined (WITH_FLOWOL)
    FlowOL * flowol = FLOW2D3D->flowol;
    if (flowol == NULL)
        return;

    char * name = createCString (ftn_name, len_name);
    char * dir  = "";

    try {
        flowol->dol->ArrayShape (dir, name, *dimensionality, dimensions);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "DOL ArrayShape fails: %s", ex->message);
        }

    delete [] name;
#endif
    }


void STDCALL
FLOWOL_ChangeDirectory (
    const char *    ftn_dirname,
    int             len_dirname
    ) {

#if defined (WITH_FLOWOL)
    FlowOL * flowol = FLOW2D3D->flowol;
    if (flowol == NULL)
        return;

    char * dirname = createCString (ftn_dirname, len_dirname);

    try {
        flowol->dol->ChangeDirectory (dirname);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "Cannot set DOL description: %s", ex->message);
        }

    delete [] dirname;
#endif
    }


void STDCALL
FLOWOL_Publish_c (
#if defined(MIXED_STR_LEN_ARG)
    const char *    ftn_name,
    int             len_name,
    const char *    ftn_description,
    int             len_description,
    const char *    ftn_units,
    int             len_units,
    const char *    ftn_definedon,
    int             len_definedon,
    const char *    ftn_arrayshape,
    int             len_arrayshape,
    int *           basetype,
    void *          address,
    int *           inout
#else
    const char *    ftn_name,
    const char *    ftn_description,
    const char *    ftn_units,
    const char *    ftn_definedon,
    const char *    ftn_arrayshape,
    int *           basetype,
    void *          address,
    int *           inout,
    int             len_name,
    int             len_description,
    int             len_units,
    int             len_definedon,
    int             len_arrayshape
#endif
    ) {

#if defined (WITH_FLOWOL)
    FlowOL * flowol = FLOW2D3D->flowol;
    if (flowol == NULL)
        return;

    char * name        = createCString (ftn_name, len_name);
    char * description = createCString (ftn_description, len_description);
    char * units       = createCString (ftn_units, len_units);
    char * definedon   = createCString (ftn_definedon, len_definedon);
    char * arrayshape  = createCString (ftn_arrayshape, len_arrayshape);
    char * dir = "";

    try {
        flowol->dol->Publish (dir, name, description, units, definedon, arrayshape, *basetype, address, *inout);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "DOL Publish fails: %s", ex->message);
        }

    delete [] name;
    delete [] description;
    delete [] units;
    delete [] definedon;
    delete [] arrayshape;
#endif
    }


void STDCALL
FLOWOL_Publish_string_c (
#if defined(MIXED_STR_LEN_ARG)
    const char *    ftn_name,
    int             len_name,
    const char *    ftn_description,
    int             len_description,
    const char *    ftn_units,
    int             len_units,
    const char *    ftn_definedon,
    int             len_definedon,
    const char *    ftn_arrayshape,
    int             len_arrayshape,
    int *           basetype,
    void *          address,
    int             len_addres,
    int *           inout
#else
    const char *    ftn_name,
    const char *    ftn_description,
    const char *    ftn_units,
    const char *    ftn_definedon,
    const char *    ftn_arrayshape,
    int *           basetype,
    void *          address,
    int *           inout,
    int             len_name,
    int             len_description,
    int             len_units,
    int             len_definedon,
    int             len_arrayshape,
    int             len_address
#endif
    ) {

#if defined (WITH_FLOWOL)
    FlowOL * flowol = FLOW2D3D->flowol;
    if (flowol == NULL)
        return;

#if defined(MIXED_STR_LEN_ARG)
    FLOWOL_Publish_c (ftn_name, len_name, ftn_description, len_description, ftn_units, len_units, ftn_definedon, len_definedon, ftn_arrayshape, len_arrayshape, basetype, address, inout);
#else
    FLOWOL_Publish_c (ftn_name, ftn_description, ftn_units, ftn_definedon, ftn_arrayshape, basetype, address, inout, len_name, len_description, len_units, len_definedon, len_arrayshape);
#endif
#endif
    }


void STDCALL
FLOWOL_PublishFunction (
#if defined(MIXED_STR_LEN_ARG)
    const char *    ftn_name,
    int             len_name,
    const char *    ftn_description,
    int             len_description,
    void *          address,
    void *          dataptr
#else
    const char *    ftn_name,
    const char *    ftn_description,
    void *          address,
    void *          dataptr,
    int             len_name,
    int             len_description
#endif
    ) {

#if defined (WITH_FLOWOL)
    FlowOL * flowol = FLOW2D3D->flowol;
    if (flowol == NULL)
        return;

    char * name = createCString (ftn_name, len_name);
    char * description = createCString (ftn_description, len_description);
    char * dir = "";

    try {
        flowol->dol->PublishFunction (dir, name, description, DOL::FORTRAN, (int (STDCALL *)(void *, const int *)) address, dataptr);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "DOL PublishFunction fails: %s", ex->message);
        }

    delete [] name;
    delete [] description;
#endif
    }


void STDCALL
FLOWOL_SetDescription (
#if defined(MIXED_STR_LEN_ARG)
    const char *    ftn_prog,
    int             len_prog,
    const char *    ftn_model,
    int             len_model
#else
    const char *    ftn_prog,
    const char *    ftn_model,
    int             len_prog,
    int             len_model
#endif
    ) {

#if defined (WITH_FLOWOL)
    FlowOL * flowol = FLOW2D3D->flowol;
    if (flowol == NULL)
        return;

    char * prog = createCString (ftn_prog, len_prog);
    char * model = createCString (ftn_model, len_model);

    char * desc = new char [strlen (prog) + strlen (model) + 10];
    sprintf (desc, "%s :: %s", prog, model);

    delete [] prog;
    delete [] model;

    try {
        flowol->dol->SetDescription (desc);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "Cannot set DOL description: %s", ex->message);
        }

    delete [] desc;
#endif
    }


void STDCALL
FLOWOL_Timestep (
    int *   timestep
    ) {

#if defined (WITH_FLOWOL)
    FlowOL * flowol = FLOW2D3D->flowol;
    if (flowol == NULL || flowol->dol == NULL)
        return;

    try {
        flowol->dol->PassMilestone ((DOL::Milestone) *timestep);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "DOL fails: %s", ex->message);
        }
#endif
    }


//-------------------------------------------------------------------------------
//  Retraction functions (called from Fortran)


void STDCALL
FLOWOL_Retract (
    const char *    ftn_name,
    int             len_name
    ) {

#if defined (WITH_FLOWOL)
    FlowOL * flowol = FLOW2D3D->flowol;
    if (flowol == NULL)
        return;

    char * name = createCString (ftn_name, len_name);

    try {
        flowol->dol->Retract (NULL, name);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "Cannot Retract DOL DataElement: %s", ex->message);
        }

    delete [] name;
#endif
    }


void STDCALL
FLOWOL_RetractArrayShape (
    const char *    ftn_name,
    int             len_name
    ) {

#if defined (WITH_FLOWOL)
    FlowOL * flowol = FLOW2D3D->flowol;
    if (flowol == NULL)
        return;

    char * name = createCString (ftn_name, len_name);

    try {
        flowol->dol->RetractArrayShape (NULL, name);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "Cannot Retract DOL ArrayShape: %s", ex->message);
        }

    delete [] name;
#endif
    }


void STDCALL
FLOWOL_RetractFunction (
    const char *    ftn_name,
    int             len_name
    ) {

#if defined (WITH_FLOWOL)
    FlowOL * flowol = FLOW2D3D->flowol;
    if (flowol == NULL)
        return;

    char * name = createCString (ftn_name, len_name);

    try {
        flowol->dol->RetractFunction (NULL, name);
        }

    catch (DOL::Exception * ex) {
        throw new Exception (true, "Cannot Retract DOL Function: %s", ex->message);
        }

    delete [] name;
#endif
    }

}




