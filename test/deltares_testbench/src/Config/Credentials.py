'''
Description: Credentials Data Class
-----------------------------------------------------
Copyright (C)  Stichting Deltares, 2013
'''

# (network) credentials
class Credentials(object):

   def __init__(self):
      self.__name: str = ""
      self.__username: str = ""
      self.__password: str = ""

   @property
   def name(self) -> str:
      return self.__name

   @name.setter
   def name(self, value: str):
      self.__name = value

   @property
   def username(self) -> str:
      return self.__username

   @username.setter
   def username(self, value: str):
      self.__username = value

   @property
   def password(self) -> str:
      return self.__password

   @password.setter
   def password(self, value):
      self.__password = value