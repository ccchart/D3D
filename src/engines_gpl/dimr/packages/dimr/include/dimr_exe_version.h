#ifndef DIMR_EXE_VERSION
#define DIMR_EXE_VERSION

#define CAT(a, b) a ## b
#define FUNC_CAT(a, b) CAT(a, b)

#define MOD_NAME         DIMR_EXE
#define modname_program  "DIMR_EXE"
#if HAVE_CONFIG_H
#   define F90_MOD_NAME   FC_FUNC(dimr, DIMR)
#else
#   define F90_MOD_NAME   MOD_NAME
#endif

#define modname_major "2"
#define modname_minor "00"
#define modname_revision "00"
#define modname_build "000000"

#define modname_company "Deltares"
#define modname_company_url "http://www.deltares.nl"

#define modname_sourcecode_url "@(#) $HeadURL$"

/*=================================================== DO NOT MAKE CHANGES BELOW THIS LINE ===================================================================== */

static char modname_version       [] = {modname_major "." modname_minor "." modname_revision "." modname_build};
static char modname_version_short [] = {modname_major "." modname_minor};
static char modname_version_full  [] = {modname_company ", " modname_program " Version " modname_major "." modname_minor "." modname_revision "." modname_build ", " __DATE__ ", " __TIME__ ""};
static char modname_url           [] = {modname_sourcecode_url};

char * getversionstring_dimr_exe(void);
char * getfullversionstring_dimr_exe(void);
char * getshortversionstring_dimr_exe(void);
char * geturlstring_dimr_exe(void);
char * getversionidstring_dimr_exe(void);

#endif /* DIMR_EXE_VERSION */

