'''
Description: Local Paths Data Class
-----------------------------------------------------
Copyright (C)  Stichting Deltares, 2013
'''

# contains locations to given root directories
class LocalPaths(object):

   def __init__(self):
      self.__casespath: str = "cases"
      self.__enginespath: str = "engines"
      self.__referencepath: str = "references"

   def getCasesPath(self) -> str:
      return self.__casespath

   def setCasesPath(self, value: str):
      self.__casespath = value

   def getEnginesPath(self) -> str:
      return self.__enginespath

   def setEnginesPath(self, value: str):
      self.__enginespath = value

   def setReferencePath(self, value: str):
      self.__referencepath = value

   def getReferencePath(self) -> str:
      return self.__referencepath
