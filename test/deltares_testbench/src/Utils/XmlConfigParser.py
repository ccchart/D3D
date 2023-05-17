'''
Description: Xml Configuration parser
-----------------------------------------------------
Copyright (C)  Stichting Deltares, 2013
'''

import copy, re

import logging
from typing import Any, Iterable, List

from lxml import etree
from src.Config.Credentials import Credentials
from src.Config.LocalPaths import LocalPaths
from src.Config.Location import Location
from src.Config.ProgramConfig import ProgramConfig
from src.Config.Type import PathType, FileType, PresenceType
from src.Config.TestCaseConfig import TestCaseConfig
from src.Config.FileCheck import FileCheck
from src.Config.FileCheck import Parameter
from src.Config.FileCheck import SkipLine


def loop(dict: dict[str, Any], key: str) -> List:
    if (key in dict):
        if type(dict[key]) == list:
            return dict[key]
        elif type(dict[key]) == dict:
            return list(dict[key].values())
        else:
            return [dict[key]]
    else:
        return []


def branch(xml_tree: etree._ElementTree, prefix: str) -> dict[str, Any]:
    new_tree = {"txt": xml_tree.text}

    for ch in xml_tree.getchildren():
        if (type(ch.tag) == str):
            branch_name = ch.tag.replace(prefix, "")

            if branch_name not in new_tree:
                new_tree[branch_name] = []

            new_tree[branch_name].append(branch(ch, prefix))

    for key, val in xml_tree.attrib.items():
        new_tree[key.replace(prefix, "")] = [val]

    return new_tree


# Parse the xml configuration file
class XmlConfigParser(object):

    def load(self, path: str, rstr: str, cred: Credentials) -> tuple[LocalPaths, List[ProgramConfig], list]:
        """load the config file

        Args:
            path (str): (relative or absolute), overwritten roots (if any)
            rstr (str): ??
            cred (Credentials): credential from the command line

        Returns:
            tuple[LocalPaths, List[ProgramConfig], list]: local paths, program configs and test case configs
        """
        self.__path: str = path
        self.__rstr: str = rstr
        self.__validate__()
        self.__initialize__()
        self.__credentials.append(cred)

        return self.__parse__()

    def __maketree__(self, path):
        global schema
        parser = etree.XMLParser(remove_blank_text=True, attribute_defaults=False)
        parsed_tree: etree._ElementTree = etree.parse(path, parser)

        parsed_tree.xinclude()
        root_node: etree._Element = parsed_tree.getroot()

        schema = root_node.nsmap[None]
        prefix = "{%s}" % schema
        root_name = root_node.tag.replace(prefix, "")

        parsed_tree = branch(root_node, prefix)
        return (parsed_tree, schema, root_name)

    def __validate__(self):
        """validate Xml file format"""
        # parser = make_parser()
        # parser.setContentHandler(ContentHandler())
        # parser.parse(self.__path)
        pass

    def __initialize__(self):
        """initialize defaults"""
        self.__credentials: List[Credentials] = []
        self.__locations: List[Location] = []
        self.__program_configs: List[ProgramConfig] = []
        self.__default_cases: List[TestCaseConfig] = []

    # parse the xml file
    # output: loglevel, local paths, program configs and test case configs
    def __parse__(self):
        xml_doc, schema, root_name = self.__maketree__(self.__path)

        config_tags = xml_doc["config"]

        for config_tag in config_tags:
            #
            # The following for-loop should be deleted (when preparations are finished)
            self.__credentials =  self.__credentials + list(self.__parse_credentials(config_tag))
            local_paths = self.__parse_local_paths(config_tag)
            self.__locations = list(self.__parse_locations(config_tag))

        for programs in loop(xml_doc, "programs"):
            if root_name == "deltaresTestbench_v3":
                for program in programs["program"]:
                    programInstance = self.__fillProgram__(program)
                    if programInstance is not None:
                        self.__program_configs.append(programInstance)

        for defaultCases in loop(xml_doc, "defaultTestCases"):
            for case in defaultCases["testCase"]:
                self.__default_cases.append(self.__fillCase__(case))

        result = []
        caseNr = -1
        for cases in loop(xml_doc, "testCases"):
            for case in loop(cases, "testCase"):
                caseNr = caseNr + 1
                try:
                    self.__default_cases.append(self.__fillCase__(case))
                    result.append(self.__fillCase__(case))
                except:
                    logging.warning("Something is wrong with test case: " + str(cases['testCase'][caseNr]['path'][0]['txt']) + ", test case will be ignored")
        return local_paths, self.__program_configs, result

    def __parse_locations(self, config_tag: dict[str, Any]) -> Iterable[Location]:
        for locations_tags in loop(config_tag, "locations"):
            for location_tag in locations_tags["location"]:
                yield self.__fillLocation__(location_tag)

    def __parse_local_paths(self, config_tag) -> LocalPaths:
        local_paths = LocalPaths()

        get_text = lambda d: str(d[0]["txt"])

        for lcl in loop(config_tag, "localPaths"):
            local_paths.setCasesPath(get_text(lcl["testCasesDir"]))
            local_paths.setEnginesPath(get_text(lcl["enginesDir"]))
            local_paths.setReferencePath(get_text(lcl["referenceDir"]))

        return local_paths

    def __parse_credentials(self, config_tag) -> Iterable[Credentials]:

        for credentials_tag in loop(config_tag, "credentials"):
            for credential_tag in credentials_tag["credential"]:
                new_credentials = Credentials()
                new_credentials.name = str(credential_tag["name"][0])
                new_credentials.username = str(credential_tag["username"][0]["txt"])
                new_credentials.password = str(credential_tag["password"][0]["txt"])
                yield new_credentials

    # fill network path (including default)
    # input: XmlElement
    # output: Location
    def __fillLocation__(self, element: dict[str, Any]) -> Location:
        if not "ref" in element and not "name" in element:
            return None
        if not "ref" in element:
            new_location = Location()
            new_location.name = str(element["name"][0])
            if "credential" in element:
                c = self.__getCredentials__(str(element["credential"][0]["ref"][0]))
                if not c:
                    raise XmlError("invalid credential reference value in " + new_location.name)
                new_location.credentials = c
            # overwrite roots if specified
            newroot = self.__getOverwritePaths__(self.__rstr, new_location.name, "root")
            if newroot:
                new_location.root = newroot
            else:
                new_location.root = str(element["root"][0]["txt"].strip())
        else:
            new_location = copy.deepcopy(self.__getLocations__(element["ref"][0]))
            if not new_location:
                raise XmlError("invalid network path reference value in " + element["ref"][0])
        if "type" in element:
            if str(element["type"][0]).lower() == "input":
                new_location.type = PathType.INPUT
            if str(element["type"][0]).lower() == "check":
                new_location.type = PathType.CHECK
            if str(element["type"][0]).lower() == "reference":
                new_location.type = PathType.REFERENCE
        if "path" in element:
            #  overwrite paths if specified
            newpath = self.__getOverwritePaths__(self.__rstr, new_location.name, "path")
            if newpath:
                new_location.setPath(newpath)
            else:
                new_location.setPath(str(element["path"][0]["txt"]))
        if "version" in element:
            new_location.setVersion = str(element["version"][0]["txt"])
        if "from" in element:
            # overwrite from if specified
            newfrom = self.__getOverwritePaths__(self.__rstr, new_location.name, "from")
            if newfrom:
                new_location.from_path = newfrom
            else:
                # Remove leading/trailing slashes; they mess up the building of the path
                new_location.from_path = str(element["from"][0]["txt"]).strip('/\\')
        if "to" in element:
            # overwrite to if specified
            newto = self.__getOverwritePaths__(self.__rstr, new_location.name, "to")
            if newto:
                new_location.to_path = newto
            else:
                new_location.to_path = str(element["to"][0]["txt"])
        return new_location

    # fill program (including defaults)
    # input: XmlElement
    # output: Program
    def __fillProgram__(self, element):
        p = ProgramConfig()
        if "ignore" in element:  # ignore program for this case [RL666]
            if (element["ignore"][0].lower() == "true"):
                return None
        if "name" in element:
            p.setName(str(element["name"][0]))
        if "programStringRemoveQuotes" in element and str(element["programStringRemoveQuotes"][0]).lower() == "true":
            p.setProgramRemoveQuotes(True)
        if "shellStringRemoveQuotes" in element and str(element["shellStringRemoveQuotes"][0]).lower() == "true":
            p.setShellRemoveQuotes(True)
        if "ignoreStandardError" in element and str(element["ignoreStandardError"][0]).lower() == "true":
            p.setIgnoreStandardError(True)
        if "ignoreReturnValue" in element and str(element["ignoreReturnValue"][0]).lower() == "true":
            p.setIgnoreReturnValue(True)
        if "logOutputToFile" in element and str(element["logOutputToFile"][0]).lower() == "true":
            p.setLogOutputToFile(True)
        if "storeOutput" in element and str(element["storeOutput"][0]).lower() == "true":
            p.setStoreOutput(True)
        if "addSearchPaths" in element and str(element["addSearchPaths"][0]).lower() == "true":
            p.setAddSearchPaths(True)
        if "excludeSearchPathsContaining" in element:
            p.setExcludeSearchPathsContaining(str(element["excludeSearchPathsContaining"][0]))
        if "ref" in element:
            p.setName(str(element["ref"][0]))
        if "seq" in element:
            p.setSequence(int(element["seq"][0]))
        if "delay" in element:
            p.setDelay(float(element["delay"][0]))
        if "maxRunTime" in element:
            p.setMaxRunTime(float(element["maxRunTime"][0]["txt"]))
        if "path" in element:
            # overwrite path if specified
            newpath = self.__getOverwritePaths__(self.__rstr, p.getName(), "path")
            if newpath:
                p.setPath(newpath)
            else:
                p.setPath(str(element["path"][0]["txt"]))
        if "workingDirectory" in element:
            p.setWorkingDirectory(str(element["workingDirectory"][0]["txt"]))
            # logging.debug (p.getWorkingDirectory())
        for e in loop(element, "location"):
            nwp = self.__fillLocation__(e)
            if nwp:
                nwpExists = False
                for enp in p.getLocations():
                    if enp.name == nwp.name and enp.type() == nwp.type():
                        enp = nwp
                        nwpExists = True
                if not nwpExists:
                    p.getLocations().append(nwp)
        for el in loop(element, "shell"):
            shellProgram = self.__getPrograms__(str(el["ref"][0]))
            if shellProgram == None:
                raise XmlError("Can not find shell program '" + str(el["ref"][0]) + \
                               "'. Is this program defined in the config.xml file before being used as shell?")
            else:
                p.setShell(shellProgram)
        for el in loop(element, "arguments"):
            for package in el["argument"]:
                p.getArguments().append(str(package["txt"]))
        for el in loop(element, "modules"):
            for module in el["module"]:
                p.getModules().append(str(module["txt"]))
        for el in loop(element, "environments"):
            for env in el["environment"]:
                # append search paths if necessary
                if str(env["name"][0]).lower() == "%path%":
                    p.getSearchPaths().append(str(env["txt"]))
                else:
                    p.getEnvironmentVariables()[str(env["name"][0])] = [str(env["type"][0]), str(env["txt"])]
        return p

    # fill file checks
    # input: XmlElement
    # output: list of FileCheck
    def __fillFileCheck__(self, element):
        defined_file_types = {"ascii": FileType.ASCII,
                              "nefis": FileType.NEFIS,
                              "his": FileType.HIS,
                              "netcdf": FileType.NETCDF,
                              "numbertext": FileType.NUMBERTEXT,
                              "dseriesregression": FileType.DSERIESREGRESSION,
                              "dseriesverification": FileType.DSERIESVERIFICATION,
                              "timeseries_pi": FileType.PITIMESERIES,
                              "timeseries_csv": FileType.CSVTIMESERIES}
        defined_presence_types = {"present": PresenceType.PRESENT,
                              "absent": PresenceType.ABSENT}

        fc = FileCheck()
        skiplines = {}
        skipline = []
        i = -1

        fc.setName(str(element["name"][0]))
        if "ignore" in element and str(element["ignore"][0]).lower() == "true":
            fc.setIgnore(True)
        else:
            if "type" in element:
                typename = str(element["type"][0]).lower()
                if typename in defined_file_types:
                    fc.setType(defined_file_types[typename])
                else:
                    fc.setType(FileType.NONE)

            if "presence" in element:
                presencename = str(element["presence"][0]).lower()
                if presencename in defined_presence_types:
                    fc.setPresence(defined_presence_types[presencename])
                else:
                    fc.setPresence(PresenceType.NONE)

        parameters = {}
        for el in loop(element, "parameters"):
            params = []
            name = ""
            if "name" in el:
                name = str(el["name"][0])
            for param in loop(el,"parameter"):
                # parameters MUST have a name when this is a nefis file, otherwise it is optional.
                # But name is not allowed to be empty!
                if name == "":
                    if fc.getType() == FileType.NEFIS:
                        raise OSError(
                            "In config file, checkfile " + fc.getName() + " has type nefis but field <parameters> has no name attribute")
                    name = "parameters"
                p = Parameter()
                p.setName(str(param["name"][0]))
                if "location" in param:
                    p.setLocation(str(param["location"][0]))
                if "tolerance" in param:
                    p.setTolerance(float(param["tolerance"][0]))
                if "toleranceAbsolute" in param:
                    p.setToleranceAbsolute(float(param["toleranceAbsolute"][0]))
                if "toleranceRelative" in param:
                    p.setToleranceRelative(float(param["toleranceRelative"][0]))
                if "ignore" in param and str(param["ignore"][0]).lower() == "true":
                    p.setIgnore(True)
                params.append(p)
            parameters[name] = params
        for el in loop(element, "skipline"):
            i = i + 1
            p = SkipLine()
            p.setName(str(element["skipline"][i]["txt"]))
            skipline.append(p)
            skiplines['skipline'] = skipline

        fc.setParameters(parameters)
        fc.setSkipLines(skiplines)
        return fc

    # #######################################################################################################################################

    # fill cases (including default)
    # input: XmlElement
    # output: TestCaseConfig
    def __fillCase__(self, element):
        if not "ref" in element:
            c = TestCaseConfig()
            if "maxRunTime" not in element:
                raise XmlError("no maximum run time specified for " + c.getName())
        else:
            c = copy.deepcopy(self.__getCase__(str(element["ref"][0])))
            if "programs" in element:
                c.setProgramConfigs([])
        c.setName(str(element["name"][0]))
        if "ignore" in element:
            if str(element["ignore"][0]).lower() == "true":
                c.setIgnore(True)
        for e in loop(element, "location"):
            nwp = self.__fillLocation__(e)
            if nwp:
                nwpExists = False
                for enp in c.getLocations():
                    if enp.type == nwp.type:
                        enp = nwp
                        nwpExists = True
                if not nwpExists:
                    c.getLocations().append(nwp)
        # add case path
        if "path" in element:
            # overwrite path if specified
            newpath = self.__getOverwritePaths__(self.__rstr, c.getName(), "path")
            if newpath:
                c.setPath(newpath)
            else:
                c.setPath(str(element["path"][0]["txt"]))
        if "maxRunTime" in element:
            c.setMaxRunTime(float(element["maxRunTime"][0]["txt"]))
            for el in element["maxRunTime"]:
                if "OverruleRefMaxRunTime" in el and str(el["OverruleRefMaxRunTime"][0]).lower() == "true":
                    c.setOverruleRefMaxRunTime(True)
        for el in loop(element, "programs"):
            for program in loop(el, "program"):
                programInstance = self.__fillProgram__(program)
                if programInstance is not None:
                    c.getProgramConfigs().append(programInstance)
        for el in loop(element, "errors"):
            for error in el["error"]:
                c.getErrors().append(str(error["txt"]))
        for el in loop(element, "checks"):
            for check in el["file"]:
                c.getChecks().append(self.__fillFileCheck__(check))
        for el in loop(element, "shellarguments"):
            for shellargument in el["shellargument"]:
                c.getShellArguments().append(str(shellargument["txt"]))
        if "shell" in element:
            localShellName = str(element["shell"][0]["txt"])
            c.setShell(self.__getPrograms__(localShellName))
        return c

        # parse roots

    def __getOverwritePaths__(self, rstr, who, what):
        if rstr == None or rstr == "":
            return None
        pathparts = rstr.split(",")
        for part in pathparts:
            actual = PathParts()
            call, path = part.split("=")
            # only check calls for my name (in config)
            name = re.findall(r'(?<=\[)(.*?)(?=\])', call)[0]
            if not who == name:
                continue
            # only check calls with correct name (root, form, to or path)
            if not str(call).startswith(what):
                continue
            # found it
            return path
        # found nothing
        return None

    # get reference value if exists
    def __getCredentials__(self, name):
        for credential in self.__credentials:
            if credential.name == name:
                return credential
        return None

    # get reference value if exists
    def __getLocations__(self, name):
        for location in self.__locations:
            if location.name == name:
                return location
        return None

    # get reference value if exists
    def __getPrograms__(self, name):
        for program in self.__program_configs:
            if program.getName() == name:
                return program
        return None

    # get reference value if exists
    def __getCase__(self, name):
        for case in self.__default_cases:
            if case.getName() == name:
                return case
        return None


# Parsable for paths overrides
class PathParts:
    def __init__(self):
        self.__name = None
        self.__root = None
        self.__from = None
        self.__to = None
        self.__path = None

    def getName(self):
        return self.__name

    def setName(self, name):
        self.__name = name

    def getRoot(self):
        return self.__root

    def setRoot(self, root):
        self.__root = root

    def getFrom(self):
        return self.__from

    def setFrom(self, _from):
        self.__from = _from

    def getTo(self):
        return self.__to

    def setTo(self, to):
        self.__to = to

    def getPath(self):
        return self.__path

    def setPath(self, path):
        self.__path = path


# custom error for Xml handler
class XmlError(Exception):
    def __init__(self, value):
        self.__value = value

    def __str__(self):
        return repr(self.__value)
