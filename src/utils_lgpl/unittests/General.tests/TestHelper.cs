﻿using System;
using System.IO;

namespace General.tests
{
    public static class TestHelper
    {
        public static void DeleteIfExists(string filename)
        {
            if (File.Exists(filename))
            {
                File.Delete(filename);
            }
        }

        public static string CreateLocalCopy(string filename)
        {
            string targetFilePath = TestDirectoryPath() + @"\" + filename;
            string sourceFilePath = TestFilesDirectoryPath() + @"\" + filename;
            File.Copy(sourceFilePath, targetFilePath, true);
            return targetFilePath;
        }

        public static string TestDirectoryPath()
        {
            FileInfo fileInfo = new FileInfo(AppDomain.CurrentDomain.BaseDirectory);
            return fileInfo.FullName;
        }

        public static string TestFilesDirectoryPath()
        {
            FileInfo fileInfo = new FileInfo(AppDomain.CurrentDomain.BaseDirectory);
            string path = fileInfo.Directory.Parent.Parent.Parent.Parent.FullName + @"\test_data";
            return path;
        }

        public static string GetLibraryPath(string libName)
        {
            FileInfo fileInfo = new FileInfo(AppDomain.CurrentDomain.BaseDirectory);
            string path = fileInfo.Directory.Parent.Parent.Parent.Parent.Parent.Parent.FullName;
            bool is64bit = Environment.Is64BitProcess;
            string prefix = @"\";
            // If 64-bit process, load 64-bit DLL otherwise load the 32 bit dll 
            if (is64bit)
            {
                prefix = @"\x64\";
            }
            path = path + @"\" + libName + @"\packages\" + libName + @"\dll" + prefix + NativeLibrary.mode + @"\" + libName + @".dll";
            return path;
        }

        public static IntPtr LoadLibrary(string dllToLoad)
        {
            //uint ierr = GetLastError(); //use it to check the errors
            IntPtr _ptr = NativeLibrary.LoadLibrary(dllToLoad);
            return _ptr;
        }

        //this function sets the path in case the dll depends on other dlls. It is costumized for io_netcdf (e.g. @"\ifort12")
        public static void SetSharedPath(string libName)
        {
            FileInfo fileInfo = new FileInfo(AppDomain.CurrentDomain.BaseDirectory);
            string path = fileInfo.Directory.Parent.Parent.Parent.Parent.Parent.Parent.Parent.FullName;
            bool is64bit = Environment.Is64BitProcess;
            string prefix = @"\win32";
            // If 64-bit process, load 64-bit DLL otherwise load the 32 bit dll 
            if (is64bit)
            {
                prefix = @"\x64";
            }
            path = path + @"\third_party_open\" + libName + @"\src" + prefix + @"\2005\libsrc\" + NativeLibrary.mode;
            var envpath = Environment.GetEnvironmentVariable("PATH");
            if (envpath != null && envpath.Contains(path)) return;

            envpath = envpath + ";" + path;
            Environment.SetEnvironmentVariable("PATH", envpath, EnvironmentVariableTarget.Process);
        }
    }
}
