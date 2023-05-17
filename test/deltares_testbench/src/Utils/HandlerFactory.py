'''
Description: Handler Factory
-----------------------------------------------------
Copyright (C)  Stichting Deltares, 2013
'''

import os, logging
from distutils import dir_util
from src.Config.Credentials import Credentials
from src.Utils.IHandler import IHandler
from src.Utils.Paths import Paths
from src.Config.HandlerType import HandlerType
from src.Utils.ResolveHandler import ResolveHandler
from src.Utils.FtpHandler import FTPHandler
from src.Utils.SvnHandler import SvnHandler
from src.Utils.HttpHandler import HTTPHandler
from src.Utils.LocalNetHandler import LocalNetHandler
from src.Utils.Unzipper import Unzipper


# Chooses which type of handler is used for upload and download actions
class HandlerFactory(object):

    @classmethod
    def get_handler(self, to_path: str, credentials: Credentials=None, autocommit: bool=False)-> IHandler:
        """Creates handler based on destination path

        Args:
            to_path (str): destination path
            credentials (Credentials, optional): credentials needed for connection.
                                                 Defaults to None.
            autocommit (bool): use auto commit (for svn)

        Raises:
            AttributeError: if handler could not be detected

        Returns:
            IHandler: Specific handler
        """
        handler_type = ResolveHandler().detect(to_path, credentials)
        handler: IHandler = None

        if handler_type == HandlerType.WEB:
            logging.debug("using HTTP handler for %s", to_path)
            handler = HTTPHandler()
        if handler_type == HandlerType.SVN:
            logging.debug("using SVN handler for %s", to_path)
            handler = SvnHandler(autocommit)
        if handler_type == HandlerType.FTP:
            logging.debug("using FTP handler for %s", to_path)
            handler = FTPHandler()
        if handler_type == HandlerType.NET or handler_type == HandlerType.PATH:
            logging.debug("using LocalNet handler for %s", to_path)
            handler = LocalNetHandler()
        if handler_type == HandlerType.NONE:
            raise AttributeError("upload :: no type specified")

        return handler

    def prepare_upload(self, from_path: str, to_path: str, credentials: Credentials=None):

        rfp = Paths().rebuildToLocalPath(from_path)

        handler = self.get_handler(to_path, credentials)
        handler.prepare_upload(rfp, to_path, credentials)

    # Upload data to location
    # input: from, to (assumes this is network) and optional credentials
    def upload(self, from_path, to_path, credentials=None, autocommit=False):
        rfp = Paths().rebuildToLocalPath(from_path)
        handler = self.get_handler(to_path, credentials, autocommit)
        handler.upload(rfp, to_path, credentials)

    # Download data from location
    # input: from, to and optional credentials
    def download(self, from_path, to_path, credentials=None, version=None, unzip=False, autocommit=False):

        rtp = Paths().rebuildToLocalPath(to_path)
        if os.path.exists(rtp):
            dir_util.remove_tree(rtp)
        os.makedirs(rtp)

        handler = self.get_handler(from_path, credentials, autocommit)
        try:
            handler.download(from_path, rtp, credentials, version)
        except Exception as e:
            raise e

        if unzip:
            Unzipper().recursive(rtp)
