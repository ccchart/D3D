﻿using System;
using System.IO;
using NUnit.Framework;
using System.Runtime.InteropServices;
using General.tests;

// The build of this test is disabled by default because it requires NUnit.
// If you decide to build the test make sure to install Nunit in your solution.


namespace ECModuleTests
{
    public class ECModuleTests
    {
       //The Constructor loads the library
        static ECModuleTests()
        {
            string filename = TestHelper.GetLibraryPath(Ec_ModuleLibWrapper.LibDetails.LIB_NAME);
            __libecmodule = TestHelper.LoadLibrary(filename);
            Assert.That(__libecmodule, Is.Not.Null);
        }

        //pointer to the loaded dll
        public static IntPtr __libecmodule;

        // Create the netcdf files
        [Test]
        [NUnit.Framework.Category("TestTriangleInterpolation")]
        public void TestTriangleInterpolation()
        {
            string tempPath = TestHelper.TestDirectoryPath();
            string pathIn = TestHelper.TestFilesDirectoryPath() + @"\inTestTriangleInterpolation.txt";
            Assert.IsTrue(File.Exists(pathIn));
            string pathOut = TestHelper.TestFilesDirectoryPath() + @"\outTestTriangleInterpolation.txt";
            Assert.IsTrue(File.Exists(pathOut));

            string[] pathInLines = File.ReadAllLines(pathIn);
            string[] firstLineIn = pathInLines[0].Split(' '.ToString().ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
            int nsample = Convert.ToInt16(firstLineIn[0].Trim(' '));
            double[] sx = new double[nsample];
            double[] sy = new double[nsample];
            double[] sv = new double[nsample];


            /// Store the values
            for (int i = 0; i < nsample; i++)
            {
                string[] val = pathInLines[i+1].Split(' '.ToString().ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                sx[i] = Convert.ToDouble(val[0]);
                sy[i] = Convert.ToDouble(val[1]);
                sv[i] = Convert.ToDouble(val[2]);
            }
            
            string[] pathOutLines = System.IO.File.ReadAllLines(pathOut);
            string[] firstLineOut = pathOutLines[0].Split( ' '.ToString().ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
            int nGridPoints = Convert.ToInt16(firstLineOut[0]);
            double[] dx = new double[nGridPoints];
            double[] dy = new double[nGridPoints];
            double[] res = new double[nGridPoints];


            /// Store the values
            for (int i = 0; i < nGridPoints; i++)
            { 
                string[] val = pathOutLines[i+1].Split(' '.ToString().ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                dx[i] = Convert.ToDouble(val[0]);
                dy[i] = Convert.ToDouble(val[1]);
                res[i] = Convert.ToDouble(val[2]);
            }

            IntPtr cptr_sx = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nsample);
            IntPtr cptr_sy = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nsample);
            IntPtr cptr_sv = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nsample);


            IntPtr cptr_dx = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nGridPoints);
            IntPtr cptr_dy = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nGridPoints);
            IntPtr cptr_res = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nGridPoints);

            Marshal.Copy(sx, 0,cptr_sx, nsample);
            Marshal.Copy(sy, 0, cptr_sy, nsample);
            Marshal.Copy(sv, 0, cptr_sv, nsample);

            Marshal.Copy(dx, 0, cptr_dx, nGridPoints);
            Marshal.Copy(dy, 0, cptr_dy, nGridPoints);


            var wrapper = new Ec_ModuleLibWrapper();
            int ierr = wrapper.triang(ref cptr_sx, ref cptr_sy, ref cptr_sv, ref nsample,
                    ref cptr_dx, ref cptr_dy, ref nGridPoints, ref cptr_res);
            Assert.That(ierr, Is.EqualTo(0));

            //check the interpolation results
            double[] resInterp = new double[nGridPoints];
            Marshal.Copy(cptr_res, resInterp, 0, nGridPoints);

            for (int i = 0; i < nGridPoints; i++)
            {
                Assert.That(resInterp[i], Is.EqualTo(res[i]).Within(1e-6));
            }

        }

        // Create the netcdf files
        [Test]
        [NUnit.Framework.Category("TestAverageInterpolation")]
        public void TestAverageInterpolation()
        {
            string tempPath = TestHelper.TestDirectoryPath();
            string pathIn = TestHelper.TestFilesDirectoryPath() + @"\inTestAverageInterpolation.txt";
            Assert.IsTrue(File.Exists(pathIn));
            string pathOut = TestHelper.TestFilesDirectoryPath() + @"\outTestAverageInterpolation.txt";
            Assert.IsTrue(File.Exists(pathOut));

            string[] pathInLines = File.ReadAllLines(pathIn);
            string[] firstLineIn = pathInLines[0].Split(' '.ToString().ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
            int nsample = Convert.ToInt16(firstLineIn[0].Trim(' '));
            double[] sx = new double[nsample];
            double[] sy = new double[nsample];
            double[] sv = new double[nsample];


            /// Store the values
            for (int i = 0; i < nsample; i++)
            {
                string[] val = pathInLines[i + 1].Split(' '.ToString().ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                sx[i] = Convert.ToDouble(val[0]);
                sy[i] = Convert.ToDouble(val[1]);
                sv[i] = Convert.ToDouble(val[2]);
            }

            string[] pathOutLines = System.IO.File.ReadAllLines(pathOut);
            string[] firstLineOut = pathOutLines[0].Split(' '.ToString().ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
            int nGridPoints = Convert.ToInt16(firstLineOut[0]);
            double[] dx = new double[nGridPoints];
            double[] dy = new double[nGridPoints];
            double[] res = new double[nGridPoints];


            /// Store the values
            for (int i = 0; i < nGridPoints; i++)
            {
                string[] val = pathOutLines[i + 1].Split(' '.ToString().ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                dx[i] = Convert.ToDouble(val[0]);
                dy[i] = Convert.ToDouble(val[1]);
                res[i] = Convert.ToDouble(val[2]);
            }

            IntPtr cptr_sx = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nsample);
            IntPtr cptr_sy = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nsample);
            IntPtr cptr_sv = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nsample);


            IntPtr cptr_dx = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nGridPoints);
            IntPtr cptr_dy = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nGridPoints);
            IntPtr cptr_res = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nGridPoints);

            Marshal.Copy(sx, 0, cptr_sx, nsample);
            Marshal.Copy(sy, 0, cptr_sy, nsample);
            Marshal.Copy(sv, 0, cptr_sv, nsample);

            Marshal.Copy(dx, 0, cptr_dx, nGridPoints);
            Marshal.Copy(dy, 0, cptr_dy, nGridPoints);

            // parameters to define!
            int jsferic = 0;
            int jasfer3D = 0;
            int c_numc = 10;
            int c_n6 = 12;
            int cptr_meth = 0;
            int cptr_nmin = 0;
            double cptr_csize = 0;

            IntPtr cptr_cxx = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nGridPoints);
            IntPtr cptr_cyy = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nGridPoints);
            IntPtr cptr_cnp = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nGridPoints);

            var wrapper = new Ec_ModuleLibWrapper();
            int ierr = wrapper.averaging(
            ref cptr_sx,
            ref cptr_sy,
            ref cptr_sv,
            ref nsample,
            ref cptr_dx,
            ref cptr_dy,
            ref cptr_cxx,
            ref cptr_cyy,
            ref cptr_cnp,
            ref c_numc,
            ref c_n6,
            ref cptr_res,
            ref cptr_meth,
            ref cptr_nmin,
            ref cptr_csize,
            ref jsferic,
            ref jasfer3D);

            Assert.That(ierr, Is.EqualTo(0));

            //check the interpolation results
            double[] resInterp = new double[nGridPoints];
            Marshal.Copy(cptr_res, resInterp, 0, nGridPoints);

            for (int i = 0; i < nGridPoints; i++)
            {
                Assert.That(resInterp[i], Is.EqualTo(res[i]).Within(1e-6));
            }

        }

    }
}
