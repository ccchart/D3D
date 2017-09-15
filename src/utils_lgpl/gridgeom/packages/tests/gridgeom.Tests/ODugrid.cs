﻿using System;
using System.Runtime.InteropServices;
using NUnit.Framework;
using NUnit.Framework.Internal;
using General.tests;
using UGrid.tests;
using System.IO;

//To compile and run this test you add NUnit to your visual studio solution and 
// enable the build in configuration manager
// todo: enable pre-build events to generate the version number

namespace gridgeom.Tests
{
    public class ODugrid
    {
        //Constructor loads the library
        static ODugrid()
        {
            string filename = TestHelper.GetLibraryPath(LibDetails.LIB_NAME);
            _gridgeom_libptr = TestHelper.LoadLibrary(filename);
            //we should chek the pointer is not null
            Assert.That(_gridgeom_libptr, Is.Not.Null);

            //load netcdf for reading in meshgeom
            TestHelper.SetSharedPath(LibDetails.NETCDF_DEP);
            filename = TestHelper.GetLibraryPath(LibDetails.NETCDF_LIB_NAME);
            _netcdf_libptr = TestHelper.LoadLibrary(filename);
            Assert.That(_netcdf_libptr, Is.Not.Null);
        }

        //pointer to the loaded dll
        public static IntPtr _gridgeom_libptr;
        //pointer to the loaded dll
        public static IntPtr _netcdf_libptr;

        [Test]
        [NUnit.Framework.Category("ougridTests")]
        public void convert1DUgridToXY()
        {
            //dimension info
            int nmeshpoints = 10;
            int nbranches = 3;
            int ngeopoints = 9;

            //branches info
            double[] branchlengths = { 4.0, 3.0, 3.0 };
            int[] nbranchgeometrynodes = { 3, 3, 3 };

            //geometry info
            double[] geopointsX = { 1.0, 3.0, 5.0, 5.0, 5.0, 5.0, 5.0, 7.0, 8.0 };
            double[] geopointsY = { 4.0, 4.0, 4.0, 1.0, 2.0, 4.0, 4.0, 4.0, 4.0 };

            //mesh geometry
            int[] branchids = { 1, 1, 1, 1, 2, 2, 2, 3, 3, 3 };
            double[] branchoffsets = { 0.0, 2.0, 3.0, 4.0, 0.0, 1.5, 3.0, 0.0, 1.5, 3.0 };

            // Create the netcdf files
            double[] meshXCoords = { 1.0, 3.0, 4.0, 5.0, 5.0, 5.0, 5.0, 5.0, 6.5, 8.0 };
            double[] meshYCoords = { 4.0, 4.0, 4.0, 4.0, 1.0, 2.5, 4.0, 4.0, 4.0, 4.0 };

            IntPtr c_branchids = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nmeshpoints);
            IntPtr c_branchoffsets = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nmeshpoints);
            IntPtr c_geopointsX = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * ngeopoints);
            IntPtr c_geopointsY = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * ngeopoints);
            IntPtr c_nbranchgeometrynodes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nbranches);
            IntPtr c_branchlengths = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nbranches);
            IntPtr c_meshXCoords = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nmeshpoints);
            IntPtr c_meshYCoords = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nmeshpoints);
            try
            {
                var wrapper = new GridGeomLibWrapper();
                Marshal.Copy(branchids, 0, c_branchids, nmeshpoints);
                Marshal.Copy(branchoffsets, 0, c_branchoffsets, nmeshpoints);
                Marshal.Copy(geopointsX, 0, c_geopointsX, ngeopoints);
                Marshal.Copy(geopointsY, 0, c_geopointsY, ngeopoints);
                Marshal.Copy(nbranchgeometrynodes, 0, c_nbranchgeometrynodes, nbranches);
                Marshal.Copy(branchlengths, 0, c_branchlengths, nbranches);

                //call the function and assert for validity 
                int ierr = wrapper.ggeo_get_xy_coordinates(ref c_branchids, ref c_branchoffsets, ref c_geopointsX,
                    ref c_geopointsY, ref c_nbranchgeometrynodes, ref c_branchlengths, ref c_meshXCoords,
                    ref c_meshYCoords, ref nbranches, ref ngeopoints, ref nmeshpoints);
                Assert.That(ierr, Is.EqualTo(0));

                double[] rc_meshXCoords = new double[nmeshpoints];
                double[] rc_meshYCoords = new double[nmeshpoints];
                Marshal.Copy(c_meshXCoords, rc_meshXCoords, 0, nmeshpoints);
                Marshal.Copy(c_meshYCoords, rc_meshYCoords, 0, nmeshpoints);

                // check conversion is correct
                for (int i = 0; i < nmeshpoints; i++)
                {
                    Assert.That(rc_meshXCoords[i], Is.EqualTo(meshXCoords[i]));
                    Assert.That(rc_meshYCoords[i], Is.EqualTo(meshYCoords[i]));
                }
            }
            finally
            {
                Marshal.FreeCoTaskMem(c_branchids);
                Marshal.FreeCoTaskMem(c_branchoffsets);
                Marshal.FreeCoTaskMem(c_geopointsX);
                Marshal.FreeCoTaskMem(c_geopointsY);
                Marshal.FreeCoTaskMem(c_nbranchgeometrynodes);
                Marshal.FreeCoTaskMem(c_branchlengths);
                Marshal.FreeCoTaskMem(c_meshXCoords);
                Marshal.FreeCoTaskMem(c_meshYCoords);
            }
        }

        [Test]
        [NUnit.Framework.Category("ougridTests")]
        public void largeNumberOfPoints()
        {
            //dimension info
            int s_nmeshpoints = 6;
            int s_nbranches = 2;
            int s_ngeopoints = 4;
            int repetition = 100000;

            //branches info
            double[] s_branchlengths = { 2.0, 2.0 };
            int[] s_nbranchgeometrynodes = { 2, 2 };

            //geometry info
            double[] s_geopointsX = { 0.0, 2.0, 2.0, 2.0 };
            double[] s_geopointsY = { 0.0, 0.0, 0.0, -2.0 };

            //mesh geometry
            int[] s_branchids = { 1, 1, 1, 2, 2, 2 };
            double[] s_branchoffsets = { 0.0, 1.0, 2.0, 0.0, 1.0, 2.0 };

            //mesh coordinates
            double[] s_meshXCoords = { 0.0, 1.0, 2.0, 2.0, 2.0, 2.0 };
            double[] s_meshYCoords = { 0.0, 0.0, 0.0, 0.0, -1.0, -2.0 };

            //repeat small structure
            int nmeshpoints = s_nmeshpoints * repetition;
            int nbranches = s_nbranches * repetition;
            int ngeopoints = s_ngeopoints * repetition;

            double[] branchlengths = new double[nbranches];
            int[] nbranchgeometrynodes = new int[nbranches];
            double[] geopointsX = new double[ngeopoints];
            double[] geopointsY = new double[ngeopoints];
            int[] branchids = new int[nmeshpoints];
            double[] branchoffsets = new double[nmeshpoints];
            double[] meshXCoords = new double[nmeshpoints];
            double[] meshYCoords = new double[nmeshpoints];

            //generate a large number of point, by repeating the structure above
            int brid = 0;
            int geoid = 0;
            int mid = 0;
            for (int i = 0; i < repetition; i++)
            {

                for (int j = 0; j < s_nbranches; j++)
                {
                    branchlengths[brid] = s_branchlengths[j];
                    nbranchgeometrynodes[brid] = s_nbranchgeometrynodes[j];
                    brid = brid + 1;
                }
                for (int j = 0; j < s_ngeopoints; j++)
                {
                    geopointsX[geoid] = s_geopointsX[j] + i; //add positive offset
                    geopointsY[geoid] = s_geopointsY[j] - i; //add negative offset
                    geoid = geoid + 1;
                }
                for (int j = 0; j < s_nmeshpoints; j++)
                {
                    meshXCoords[mid] = s_meshXCoords[j] + i; //add positive offset
                    meshYCoords[mid] = s_meshYCoords[j] - i; //add negative offset
                    branchids[mid] = s_branchids[j] + i * 2; //add branch ids
                    branchoffsets[mid] = s_branchoffsets[j];
                    mid = mid + 1;
                }
            }

            IntPtr c_branchids = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nmeshpoints);
            IntPtr c_branchoffsets = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nmeshpoints);
            IntPtr c_geopointsX = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * ngeopoints);
            IntPtr c_geopointsY = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * ngeopoints);
            IntPtr c_nbranchgeometrynodes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nbranches);
            IntPtr c_branchlengths = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nbranches);
            IntPtr c_meshXCoords = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nmeshpoints);
            IntPtr c_meshYCoords = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nmeshpoints);
            try
            {
                var wrapper = new GridGeomLibWrapper();
                Marshal.Copy(branchids, 0, c_branchids, nmeshpoints);
                Marshal.Copy(branchoffsets, 0, c_branchoffsets, nmeshpoints);
                Marshal.Copy(geopointsX, 0, c_geopointsX, ngeopoints);
                Marshal.Copy(geopointsY, 0, c_geopointsY, ngeopoints);
                Marshal.Copy(nbranchgeometrynodes, 0, c_nbranchgeometrynodes, nbranches);
                Marshal.Copy(branchlengths, 0, c_branchlengths, nbranches);

                //call the function 
                int ierr = wrapper.ggeo_get_xy_coordinates(ref c_branchids, ref c_branchoffsets, ref c_geopointsX,
                    ref c_geopointsY, ref c_nbranchgeometrynodes, ref c_branchlengths, ref c_meshXCoords, ref c_meshYCoords, ref nbranches, ref ngeopoints, ref  nmeshpoints);
                Assert.That(ierr, Is.EqualTo(0));

                double[] rc_meshXCoords = new double[nmeshpoints];
                double[] rc_meshYCoords = new double[nmeshpoints];
                Marshal.Copy(c_meshXCoords, rc_meshXCoords, 0, nmeshpoints);
                Marshal.Copy(c_meshYCoords, rc_meshYCoords, 0, nmeshpoints);

                // check conversion is correct
                for (int i = 0; i < nmeshpoints; i++)
                {
                    Assert.That(rc_meshXCoords[i], Is.EqualTo(meshXCoords[i]));
                    Assert.That(rc_meshYCoords[i], Is.EqualTo(meshYCoords[i]));
                }
            }
            finally
            {
                Marshal.FreeCoTaskMem(c_branchids);
                Marshal.FreeCoTaskMem(c_branchoffsets);
                Marshal.FreeCoTaskMem(c_geopointsX);
                Marshal.FreeCoTaskMem(c_geopointsY);
                Marshal.FreeCoTaskMem(c_nbranchgeometrynodes);
                Marshal.FreeCoTaskMem(c_branchlengths);
                Marshal.FreeCoTaskMem(c_meshXCoords);
                Marshal.FreeCoTaskMem(c_meshYCoords);
            }
        }

        /// <summary>
        /// Read a netcdf file and populate meshgeom datastructure
        /// </summary>
        [Test]
        [NUnit.Framework.Category("readFile")]
        public void createLinksFrom1d2dFile()
        {
            //mesh2d
            int twoddim = 2;
            int twodnumnode = 121;
            int twodnumedge = 220;
            int twodnumface = 100;
            int twodmaxnumfacenodes = 4;
            int twodnumlayer = 0;
            int twodlayertype = 0;

            //mesh1d
            int oneddim = 1;
            int onednumnode = 21;
            int onednumedge = 20;
            int onednumface = 0;
            int onedmaxnumfacenodes = 0;
            int onednumlayer = 0;
            int onedlayertype = 0;
            int nt_nbranches = 1;
            int nt_ngeometry = 4;

            int numnodes = 10;
            int numedge = 1;
            //string c_path = TestHelper.TestDirectoryPath() + @"\write1d.nc";
            string c_path = TestHelper.TestFilesDirectoryPath() + @"\1d2d_ugrid_net.nc";
            Assert.IsTrue(File.Exists(c_path));
            int ioncid = 0; //file variable 
            int mode = 0;   //create in read mode
            var wrapperNetcdf = new IoNetcdfLibWrapper();
            int iconvtype = 2;
            double convversion = 0.0;
            var ierr = wrapperNetcdf.ionc_open(c_path, ref mode, ref ioncid, ref iconvtype, ref convversion);
            Assert.That(ierr, Is.EqualTo(0));

            int meshid = 1;
            var meshtwoddim = new meshgeomdim();
            ierr = wrapperNetcdf.ionc_get_meshgeom_dim(ref ioncid, ref meshid, ref meshtwoddim);
            Assert.That(ierr, Is.EqualTo(0));

            Assert.That(meshtwoddim.dim, Is.EqualTo(twoddim));
            Assert.That(meshtwoddim.numnode, Is.EqualTo(twodnumnode));
            Assert.That(meshtwoddim.numedge, Is.EqualTo(twodnumedge));
            Assert.That(meshtwoddim.numface, Is.EqualTo(twodnumface));
            Assert.That(meshtwoddim.maxnumfacenodes, Is.EqualTo(twodmaxnumfacenodes));
            Assert.That(meshtwoddim.numlayer, Is.EqualTo(twodnumlayer));
            Assert.That(meshtwoddim.layertype, Is.EqualTo(twodlayertype));

            //You need to know in advance the number of mesh points
            var meshtwod = new meshgeom();
            meshtwod.nodex = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * twodnumnode);
            meshtwod.nodey = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * twodnumnode);
            meshtwod.nodez = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * twodnumnode);
            meshtwod.edge_nodes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * meshtwoddim.numedge * 2);

            //var gridGeomWrapper = new GridGeomLibWrapper();
            bool includeArrays = true;
            ierr = wrapperNetcdf.ionc_get_meshgeom(ref ioncid, ref  meshid, ref meshtwod, ref  includeArrays);
            double[] rc_twodnodex = new double[twodnumnode];
            double[] rc_twodnodey = new double[twodnumnode];
            double[] rc_twodnodez = new double[twodnumnode];
            Marshal.Copy(meshtwod.nodex, rc_twodnodex, 0, twodnumnode);
            Marshal.Copy(meshtwod.nodey, rc_twodnodey, 0, twodnumnode);
            Marshal.Copy(meshtwod.nodez, rc_twodnodez, 0, twodnumnode);

            // mesh1d 
            int meshonedid = 2;
            var meshoneddim = new meshgeomdim();
            ierr = wrapperNetcdf.ionc_get_meshgeom_dim(ref ioncid, ref meshonedid, ref meshoneddim);
            Assert.That(ierr, Is.EqualTo(0));

            Assert.That(meshoneddim.dim, Is.EqualTo(oneddim));
            Assert.That(meshoneddim.numnode, Is.EqualTo(onednumnode));
            Assert.That(meshoneddim.numedge, Is.EqualTo(onednumedge));
            Assert.That(meshoneddim.numlayer, Is.EqualTo(onednumlayer));
            Assert.That(meshoneddim.layertype, Is.EqualTo(onedlayertype));

            var meshoned = new meshgeom();
            meshoned.nodex = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * onednumnode);
            meshoned.nodey = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * onednumnode);
            meshoned.nodez = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * onednumnode);
            meshoned.edge_nodes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * onednumedge * 2);

            meshoned.branchids = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * onednumnode);
            meshoned.nbranchgeometrynodes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nt_nbranches);
            meshoned.branchoffsets = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * onednumnode);
            meshoned.geopointsX = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nt_ngeometry);
            meshoned.geopointsY = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nt_ngeometry);
            meshoned.branchlengths = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nt_nbranches);

            ierr = wrapperNetcdf.ionc_get_meshgeom(ref ioncid, ref meshonedid, ref meshoned, ref includeArrays);
            Assert.That(ierr, Is.EqualTo(0));

            double[] rc_onednodex = new double[onednumnode];
            double[] rc_onednodey = new double[onednumnode];
            double[] rc_onednodez = new double[onednumnode];
            Marshal.Copy(meshoned.nodex, rc_onednodex, 0, onednumnode);
            Marshal.Copy(meshoned.nodey, rc_onednodey, 0, onednumnode);
            Marshal.Copy(meshoned.nodez, rc_onednodez, 0, onednumnode);

            //ggeo_convert to fill in kn matrix, so we can call make1D2Dinternalnetlinks_dll
            var wrapperGridgeom = new GridGeomLibWrapper();
            ierr = wrapperGridgeom.ggeo_convert(ref meshoned, ref meshoneddim);
            Assert.That(ierr, Is.EqualTo(0));
            ierr = wrapperGridgeom.ggeo_convert(ref meshtwod, ref meshtwoddim);
            Assert.That(ierr, Is.EqualTo(0));
            //convert ugrid to xy coordinates
            ierr = wrapperGridgeom.ggeo_get_xy_coordinates(ref meshoned.branchids, ref meshoned.branchoffsets, ref meshoned.geopointsX,
                ref meshoned.geopointsY, ref meshoned.nbranchgeometrynodes, ref meshoned.branchlengths, ref meshoned.nodex, ref meshoned.nodey, ref meshoneddim.nt_nbranches, ref meshoneddim.nt_ngeometry, ref  meshoneddim.numnode);
            Assert.That(ierr, Is.EqualTo(0));
            
            //call make1d2dlinks, no argument needed (all in memory)
            ierr = wrapperGridgeom.ggeo_make1D2Dinternalnetlinks();

            //free arrays
            Marshal.FreeCoTaskMem(meshtwod.nodex);
            Marshal.FreeCoTaskMem(meshtwod.nodey);
            Marshal.FreeCoTaskMem(meshtwod.nodez);
            Marshal.FreeCoTaskMem(meshoned.edge_nodes);
            
            Marshal.FreeCoTaskMem(meshoned.nodex);
            Marshal.FreeCoTaskMem(meshoned.nodey);
            Marshal.FreeCoTaskMem(meshoned.nodez);
            Marshal.FreeCoTaskMem(meshtwod.edge_nodes);

        }

        /// <summary>
        /// In this test we read a 2d grid from file and we provide the 1d discretization points.
        /// We emulate the case where  Delta Shell has a 2d file opened and the user creates a 1d mesh and wants to generate the links.
        /// 1D discretization and links are saved in memory and written afterwards. 
        /// </summary>
        [Test]
        [NUnit.Framework.Category("readFile")]
        public void createLinksFrom2dFile()
        {
            //mesh2d
            int twoddim = 2;
            int twodnumnode = 16;
            int twodnumedge = 24;
            int twodnumface = 9;
            int twodmaxnumfacenodes = 4;
            int twodnumlayer = 0;
            int twodlayertype = 0;

            //mesh1d
            //discretization points information
            int nmeshpoints = 4;
            int nbranches = 1;
            int[] branchids = { 1, 1, 1, 1 };
            double[] meshXCoords = { -6, 5, 23, 34 };
            double[] meshYCoords = { 22, 16, 16, 7 };

            //links
            int[] arrayfrom = { 2, 8 };
            int[] arrayto = { 2, 3 };

            //1. open the file with the 2d mesh
            string c_path = TestHelper.TestFilesDirectoryPath() + @"\2d_ugrid_net.nc";
            Assert.IsTrue(File.Exists(c_path));
            int ioncid = 0; //file variable 
            int mode = 0;   //create in read mode
            var wrapperNetcdf = new IoNetcdfLibWrapper();
            int iconvtype = 2;
            double convversion = 0.0;
            var ierr = wrapperNetcdf.ionc_open(c_path, ref mode, ref ioncid, ref iconvtype, ref convversion);
            Assert.That(ierr, Is.EqualTo(0));

            //2. get the 2d mesh id
            int meshid = 1;
            ierr = wrapperNetcdf.ionc_get_2d_mesh_id(ref ioncid, ref meshid);
            Assert.That(ierr, Is.EqualTo(0));

            //3. get the dimensions of the 2d mesh
            var meshtwoddim = new meshgeomdim();
            ierr = wrapperNetcdf.ionc_get_meshgeom_dim(ref ioncid, ref meshid, ref meshtwoddim);
            Assert.That(ierr, Is.EqualTo(0));

            Assert.That(meshtwoddim.dim, Is.EqualTo(twoddim));
            Assert.That(meshtwoddim.numnode, Is.EqualTo(twodnumnode));
            Assert.That(meshtwoddim.numedge, Is.EqualTo(twodnumedge));
            Assert.That(meshtwoddim.numface, Is.EqualTo(twodnumface));
            Assert.That(meshtwoddim.maxnumfacenodes, Is.EqualTo(twodmaxnumfacenodes));
            Assert.That(meshtwoddim.numlayer, Is.EqualTo(twodnumlayer));
            Assert.That(meshtwoddim.layertype, Is.EqualTo(twodlayertype));

            //4. allocate the arrays in meshgeom for storing the 2d mesh coordinates, edge_nodes
            var meshtwod = new meshgeom();
            meshtwod.nodex = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * twodnumnode);
            meshtwod.nodey = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * twodnumnode);
            meshtwod.nodez = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * twodnumnode);
            meshtwod.edge_nodes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * meshtwoddim.numedge * 2);

            //5. get the meshgeom arrays
            bool includeArrays = true;
            ierr = wrapperNetcdf.ionc_get_meshgeom(ref ioncid, ref meshid, ref meshtwod, ref includeArrays);
            Assert.That(ierr, Is.EqualTo(0));
            double[] rc_twodnodex = new double[twodnumnode];
            double[] rc_twodnodey = new double[twodnumnode];
            double[] rc_twodnodez = new double[twodnumnode];
            Marshal.Copy(meshtwod.nodex, rc_twodnodex, 0, twodnumnode);
            Marshal.Copy(meshtwod.nodey, rc_twodnodey, 0, twodnumnode);
            Marshal.Copy(meshtwod.nodez, rc_twodnodez, 0, twodnumnode);

            //6. allocate the 1d arrays for storing the 1d coordinates and edge_nodes
            IntPtr c_meshXCoords = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nmeshpoints);
            IntPtr c_meshYCoords = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nmeshpoints);
            IntPtr c_branchids = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nmeshpoints);

            Marshal.Copy(branchids, 0, c_branchids, nmeshpoints);
            Marshal.Copy(meshXCoords, 0, c_meshXCoords, nmeshpoints);
            Marshal.Copy(meshYCoords, 0, c_meshYCoords, nmeshpoints);

            //7. fill kn (Herman datastructure) for creating the links
            var wrapperGridgeom = new GridGeomLibWrapper();
            ierr = wrapperGridgeom.ggeo_convert_1d_arrays(ref c_meshXCoords, ref c_meshYCoords, ref c_branchids, ref nbranches, ref nmeshpoints);
            Assert.That(ierr, Is.EqualTo(0));
            ierr = wrapperGridgeom.ggeo_convert(ref meshtwod, ref meshtwoddim);
            Assert.That(ierr, Is.EqualTo(0));

            //9. make the links
            ierr = wrapperGridgeom.ggeo_make1D2Dinternalnetlinks();
            Assert.That(ierr, Is.EqualTo(0));

            //10. get the number of links
            int n1d2dlinks = 0;
            ierr = wrapperGridgeom.ggeo_get_links_count(ref n1d2dlinks);
            Assert.That(ierr, Is.EqualTo(0));

            //11. get the links: arrayfrom = 2d cell index, arrayto = 1d node index 
            IntPtr c_arrayfrom = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * n1d2dlinks); //2d cell number
            IntPtr c_arrayto = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * n1d2dlinks); //1d node
            ierr = wrapperGridgeom.ggeo_get_links(ref c_arrayfrom, ref c_arrayto, ref n1d2dlinks);
            Assert.That(ierr, Is.EqualTo(0));


            int[] rc_arrayfrom = new int[n1d2dlinks];
            int[] rc_arrayto = new int[n1d2dlinks];
            Marshal.Copy(c_arrayfrom, rc_arrayfrom, 0, n1d2dlinks);
            Marshal.Copy(c_arrayto, rc_arrayto, 0, n1d2dlinks);
            for (int i = 0; i < n1d2dlinks; i++)
            {
                Assert.That(rc_arrayfrom[i], Is.EqualTo(arrayfrom[i]));
                Assert.That(rc_arrayto[i], Is.EqualTo(arrayto[i]));
            }
            //for writing the links look io_netcdf ionc_def_mesh_contact, ionc_put_mesh_contact 

            //Free 2d arrays
            Marshal.FreeCoTaskMem(meshtwod.nodex);
            Marshal.FreeCoTaskMem(meshtwod.nodey);
            Marshal.FreeCoTaskMem(meshtwod.nodez);
            Marshal.FreeCoTaskMem(meshtwod.edge_nodes);

            //Free 1d arrays
            Marshal.FreeCoTaskMem(c_meshXCoords);
            Marshal.FreeCoTaskMem(c_meshYCoords);
            Marshal.FreeCoTaskMem(c_branchids);

            //Free from and to arrays describing the links 
            Marshal.FreeCoTaskMem(c_arrayfrom);
            Marshal.FreeCoTaskMem(c_arrayto);
        }

    }
}
