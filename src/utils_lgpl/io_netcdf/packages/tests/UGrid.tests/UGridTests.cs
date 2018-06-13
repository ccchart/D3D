﻿using System;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;
using NUnit.Framework;
using General.tests;
using System.Threading;

// The build of this test is disabled by default because it requires NUnit.
// If you decide to build the test make sure to install Nunit in your solution.

namespace UGrid.tests
{
    public class UGridTests
    {
        //Constructor loads the library
        static UGridTests()
        {
            TestHelper.SetSharedPath(IoNetcdfLibWrapper.LibDetails.LIB_DEP);
            string filename = TestHelper.GetLibraryPath(IoNetcdfLibWrapper.LibDetails.LIB_NAME);
            m_libnetcdf = TestHelper.LoadLibrary(filename);
            //we should chek the pointer is not null
            Assert.That(m_libnetcdf, Is.Not.Null);
            filename = TestHelper.GetLibraryPath(GridGeomLibWrapper.LibDetails.LIB_NAME);
            m_libgridgeom = TestHelper.LoadLibrary(filename);
            Assert.That(m_libgridgeom, Is.Not.Null);
        }

        //pointer to the loaded dll
        public static IntPtr m_libnetcdf;
        public static IntPtr m_libgridgeom;

        //network name
        private string networkName = "network";

        //dimension info
        private int nNodes = 4;
        private int nBranches = 3;
        private int nGeometry = 9;
        private int startIndex = 1;

        //node info
        private double[] nodesX = { 1.0, 5.0, 5.0, 8.0 };
        private double[] nodesY = { 4.0, 4.0, 1.0, 4.0 };
        private string[] nodesids = { "node1", "node2", "node3", "node4" };
        private string[] nodeslongNames = { "nodelong1", "nodelong2", "nodelong3", "nodelong4" };
        private int[] sourcenodeid = { 1, 2, 2 };
        private int[] targetnodeid = { 2, 3, 4 };

        //branches info
        private double[] branchlengths = { 4.0, 3.0, 3.0 };
        private int[] nbranchgeometrypoints = { 3, 3, 3 };
        private string[] branchids = { "branch1", "branch2", "branch3" };
        private string[] branchlongNames = { "branchlong1", "branchlong2", "branchlong3" };
        private int[] branch_order = { -1, -1, -1 };

        //geometry info
        double[] geopointsX = { 1.0, 3.0, 5.0, 5.0, 5.0, 5.0, 5.0, 7.0, 8.0 };
        double[] geopointsY = { 4.0, 4.0, 4.0, 4.0, 2.0, 1.0, 4.0, 4.0, 4.0 };

        //mesh name
        private string meshname = "1dmesh";

        //mesh dimension
        private int nmeshpoints = 8;
        private int nedgenodes  = 7; // nmeshpoints - 1 

        //mesh geometry
        private int nmesh1dPoints = 8;
        private int[] branchidx = { 1, 1, 1, 1, 2, 2, 3, 3 };
        private double[] offset = { 0.0, 2.0, 3.0, 4.0, 1.5, 3.0, 1.5, 3.0 };

        private string[] meshnodeids = { "node1_branch1", "node2_branch1", "node3_branch1", "node4_branch1", "node1_branch2",
                                         "node2_branch2", "node3_branch2", "node1_branch3", "node2_branch3", "node3_branch3" };

        private string[] meshnodelongnames = { "node1_branch1_long_name", "node2_branch1_long_name", "node3_branch1_long_name", "node4_branch1_long_name", "node1_branch2_long_name",
                                               "node2_branch2_long_name", "node3_branch2_long_name", "node1_branch3_long_name", "node2_branch3_long_name", "node3_branch3_long_name" };

        //netcdf file specifications 
        private int iconvtype = 2;

        private double convversion = 0.0;

        //mesh links
        private string linkmeshname = "links";

        private int nlinks = 3;
        private int linkmesh1 = 1;
        private int linkmesh2 = 2;
        private int locationType1 = 1;
        private int locationType2 = 1;
        private int[] mesh1indexes = { 1, 2, 3 };
        private int[] mesh2indexes = { 1, 2, 3 };
        private string[] linksids = { "link1", "link2", "link3" };
        private string[] linkslongnames = { "linklong1", "linklong2", "linklong3" };
        private int[] contacttype = { 3, 3, 3 };

        // mesh2d
        private int numberOf2DNodes = 5;
        private int numberOfFaces = 2;
        private int numberOfMaxFaceNodes = 4;
        private double[] mesh2d_nodesX = { 0, 10, 15, 10, 5 };
        private double[] mesh2d_nodesY = { 0, 0, 5, 10, 5 };
        private double[,] mesh2d_face_nodes = { { 1, 2, 5, -999 }, { 2, 3, 4, 5 } };

        //function to check mesh1d data
        private void read1dnetwork
            (
            int ioncid, 
            int networkid, 
            ref IoNetcdfLibWrapper wrapper,
            ref int l_nnodes,
            ref int l_nbranches,
            ref int l_nGeometry,
            int l_startIndex,
            ref StringBuilder l_networkName,
            ref IoNetcdfLibWrapper.interop_charinfo[] l_nodesinfo,
            ref IoNetcdfLibWrapper.interop_charinfo[] l_branchinfo,
            ref double[] l_nodesX,
            ref double[] l_nodesY,
            ref int[] l_sourcenodeid,
            ref int[] l_targetnodeid,
            ref double[] l_branchlengths,
            ref int[] l_nbranchgeometrypoints,
            ref double[] l_geopointsX,
            ref double[] l_geopointsY,
            ref int[] l_branch_order
            )
        {
            //1. Get the mesh name
            int ierr = -1;
            ierr = wrapper.ionc_get_1d_network_name(ref ioncid, ref networkid, l_networkName);
            Assert.That(ierr, Is.EqualTo(0));

            //2. Get the node count
            ierr = wrapper.ionc_get_1d_network_nodes_count(ref ioncid, ref networkid, ref l_nnodes);
            Assert.That(ierr, Is.EqualTo(0));
            IntPtr c_nodesX = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nnodes);
            IntPtr c_nodesY = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nnodes);

            //3. Get the number of branches
            ierr = wrapper.ionc_get_1d_network_branches_count(ref ioncid, ref networkid, ref l_nbranches);
            Assert.That(ierr, Is.EqualTo(0));
            IntPtr c_sourcenodeid = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nbranches);
            IntPtr c_targetnodeid = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nbranches);
            IntPtr c_branchlengths = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nbranches);
            IntPtr c_nbranchgeometrypoints = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nbranches);
            IntPtr c_branch_order = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nbranches);

            //4. Get the number of geometry points
            ierr = wrapper.ionc_get_1d_network_branches_geometry_coordinate_count(ref ioncid, ref networkid, ref l_nGeometry);
            Assert.That(ierr, Is.EqualTo(0));
            IntPtr c_geopointsX = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nGeometry);
            IntPtr c_geopointsY = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nGeometry);

            //5. Get nodes info and coordinates
            ierr = wrapper.ionc_read_1d_network_nodes(ref ioncid, ref networkid, ref c_nodesX, ref c_nodesY,
                l_nodesinfo, ref l_nnodes);
            Assert.That(ierr, Is.EqualTo(0));
            Marshal.Copy(c_nodesX, l_nodesX, 0, l_nodesX.Length);
            Marshal.Copy(c_nodesY, l_nodesY, 0, l_nodesY.Length);

            //6. Get the branch info and coordinates
            ierr = wrapper.ionc_get_1d_network_branches(ref ioncid, ref networkid, ref c_sourcenodeid, ref c_targetnodeid,
                    ref c_branchlengths, l_branchinfo, ref c_nbranchgeometrypoints, ref l_nbranches, ref l_startIndex);
            Assert.That(ierr, Is.EqualTo(0));
            Marshal.Copy(c_targetnodeid, l_targetnodeid, 0, l_targetnodeid.Length);
            Marshal.Copy(c_sourcenodeid, l_sourcenodeid, 0, l_sourcenodeid.Length);
            Marshal.Copy(c_branchlengths, l_branchlengths, 0, l_branchlengths.Length);
            Marshal.Copy(c_nbranchgeometrypoints, l_nbranchgeometrypoints, 0, l_nbranchgeometrypoints.Length);

            //7. Get the 1d branch geometry
            ierr = wrapper.ionc_read_1d_network_branches_geometry(ref ioncid, ref networkid, ref c_geopointsX,
                    ref c_geopointsY, ref l_nGeometry);
            Assert.That(ierr, Is.EqualTo(0));
            Marshal.Copy(c_geopointsX, l_geopointsX, 0, l_geopointsX.Length);
            Marshal.Copy(c_geopointsY, l_geopointsY, 0, l_geopointsY.Length);

            //8. Get the branch order 
            ierr = wrapper.ionc_get_1d_network_branchorder(ref ioncid, ref networkid, ref c_branch_order, ref nBranches);
            Assert.That(ierr, Is.EqualTo(0));

            Marshal.Copy(c_branch_order, l_branch_order, 0, l_branch_order.Length);

            //free all pointers
            Marshal.FreeCoTaskMem(c_nodesX);
            Marshal.FreeCoTaskMem(c_nodesY);
            Marshal.FreeCoTaskMem(c_sourcenodeid);
            Marshal.FreeCoTaskMem(c_targetnodeid);
            Marshal.FreeCoTaskMem(c_branchlengths);
            Marshal.FreeCoTaskMem(c_nbranchgeometrypoints);
            Marshal.FreeCoTaskMem(c_geopointsX);
            Marshal.FreeCoTaskMem(c_geopointsY);
            Marshal.FreeCoTaskMem(c_branch_order);
        }

        //function to check mesh1d data
        private void check1dmesh(int ioncid, int meshid, ref IoNetcdfLibWrapper wrapper)
        {
            //mesh variables
            IntPtr c_branchidx = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nmeshpoints);
            IntPtr c_offset = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nmeshpoints);
            //links variables
            IntPtr c_mesh1indexes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nlinks);
            IntPtr c_mesh2indexes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nlinks);
            IntPtr c_contacttype = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nlinks);
            try
            {
                //1. Get the mesh name
                int ierr = -1;
                var rmeshName = new StringBuilder(IoNetcdfLibWrapper.LibDetails.MAXSTRLEN);
                ierr = wrapper.ionc_get_mesh_name(ref ioncid, ref meshid, rmeshName);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.That(rmeshName.ToString().Trim(), Is.EqualTo(meshname));

                //2. Get the number of mesh points
                int rnmeshpoints = -1;
                ierr =
                    wrapper.ionc_get_1d_mesh_discretisation_points_count(ref ioncid, ref meshid,
                        ref rnmeshpoints);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.That(rnmeshpoints, Is.EqualTo(nmeshpoints));

                //3. Get the coordinates of the mesh points
                IoNetcdfLibWrapper.interop_charinfo[] nodeinfo = new IoNetcdfLibWrapper.interop_charinfo[nmesh1dPoints];

                ierr = wrapper.ionc_get_1d_mesh_discretisation_points(ref ioncid, ref meshid, ref c_branchidx,
                    ref c_offset, nodeinfo, ref rnmeshpoints, ref startIndex);
                Assert.That(ierr, Is.EqualTo(0));

                for (int i = 0; i < nmesh1dPoints; i++)
                {
                    string tmpstring = new string(nodeinfo[i].ids);
                    Assert.That(tmpstring.Trim(), Is.EqualTo(meshnodeids[i]));
                    tmpstring = new string(nodeinfo[i].longnames);
                    Assert.That(tmpstring.Trim(), Is.EqualTo(meshnodelongnames[i]));
                }

                int[] rc_branchidx = new int[rnmeshpoints];
                double[] rc_offset = new double[rnmeshpoints];
                Marshal.Copy(c_branchidx, rc_branchidx, 0, rnmeshpoints);
                Marshal.Copy(c_offset, rc_offset, 0, rnmeshpoints);
                for (int i = 0; i < rnmeshpoints; i++)
                {
                    Assert.That(rc_branchidx[i], Is.EqualTo(branchidx[i]));
                    Assert.That(rc_offset[i], Is.EqualTo(offset[i]));
                }

                //4. Get the number of links. TODO: add a mesh 2d for this to work!
                //int linkmesh = 1;
                //int r_nlinks = -1;
                //ierr = wrapper.ionc_get_contacts_count(ref ioncid, ref linkmesh, ref r_nlinks);
                //Assert.That(ierr, Is.EqualTo(0));
                //Assert.That(r_nlinks, Is.EqualTo(nlinks));
                //IoNetcdfLibWrapper.interop_charinfo[] linksinfo = new IoNetcdfLibWrapper.interop_charinfo[nlinks];

                ////5. Get the links values
                //ierr = wrapper.ionc_get_mesh_contact(ref ioncid, ref linkmesh, ref c_mesh1indexes, ref c_mesh2indexes, ref c_contacttype,
                //    linksinfo, ref nlinks, ref startIndex);
                //Assert.That(ierr, Is.EqualTo(0));
                //int[] rc_contacttype = new int[nlinks];
                //int[] rc_mesh1indexes = new int[nlinks];
                //int[] rc_mesh2indexes = new int[nlinks];
                //Marshal.Copy(c_mesh1indexes, rc_mesh1indexes, 0, nlinks);
                //Marshal.Copy(c_mesh2indexes, rc_mesh2indexes, 0, nlinks);
                //Marshal.Copy(c_contacttype, rc_contacttype, 0, nlinks);
                //for (int i = 0; i < nlinks; i++)
                //{
                //    string tmpstring = new string(linksinfo[i].ids);
                //    Assert.That(tmpstring.Trim(), Is.EqualTo(linksids[i]));
                //    tmpstring = new string(linksinfo[i].longnames);
                //    Assert.That(tmpstring.Trim(), Is.EqualTo(linkslongnames[i]));
                //    Assert.That(rc_mesh1indexes[i], Is.EqualTo(mesh1indexes[i]));
                //    Assert.That(rc_mesh2indexes[i], Is.EqualTo(mesh2indexes[i]));
                //    Assert.That(rc_contacttype[i], Is.EqualTo(contacttype[i]));
                //}

                //6. Get the written nodes ids
                //StringBuilder varname = new StringBuilder("node_ids");
                //IoNetcdfLibWrapper.interop_charinfo[] nodeidsvalues = new IoNetcdfLibWrapper.interop_charinfo[nmesh1dPoints];

                //ierr = wrapper.ionc_get_var_chars(ref ioncid, ref meshid, varname, nodeidsvalues, ref nmesh1dPoints);
                //Assert.That(ierr, Is.EqualTo(0));
                //for (int i = 0; i < nmesh1dPoints; i++)
                //{
                //    string tmpstring = new string(nodeidsvalues[i].ids);
                //    Assert.That(tmpstring.Trim(), Is.EqualTo(meshnodeids[i]));
                //}
            }
            finally
            {
                Marshal.FreeCoTaskMem(c_branchidx);
                Marshal.FreeCoTaskMem(c_offset);
                Marshal.FreeCoTaskMem(c_mesh1indexes);
                Marshal.FreeCoTaskMem(c_mesh2indexes);
            }
        }

        private void check2dmesh(int ioncid, int meshid, ref IoNetcdfLibWrapper wrapper)
        {
            IntPtr c_nodesX = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * numberOf2DNodes);
            IntPtr c_nodesY = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * numberOf2DNodes);
            IntPtr c_face_nodes =Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * numberOfFaces * numberOfMaxFaceNodes);
            try
            {

                int nnodes = -1;
                int ierr = wrapper.ionc_get_node_count(ref ioncid, ref meshid, ref nnodes);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.That(nnodes, Is.EqualTo(5));

                int nedge = -1;
                ierr = wrapper.ionc_get_edge_count(ref ioncid, ref meshid, ref nedge);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.That(nedge, Is.EqualTo(6));

                int nface = -1;
                ierr = wrapper.ionc_get_face_count(ref ioncid, ref meshid, ref nface);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.That(nface, Is.EqualTo(numberOfFaces));

                int maxfacenodes = -1;
                ierr = wrapper.ionc_get_max_face_nodes(ref ioncid, ref meshid, ref maxfacenodes);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.That(maxfacenodes, Is.EqualTo(numberOfMaxFaceNodes));

                //get all node coordinates
                ierr = wrapper.ionc_get_node_coordinates(ref ioncid, ref meshid, ref c_nodesX, ref c_nodesY, ref nnodes);
                Assert.That(ierr, Is.EqualTo(0));

                double[] rc_nodeX = new double[numberOf2DNodes];
                double[] rc_nodeY = new double[numberOf2DNodes];
                Marshal.Copy(c_nodesX, rc_nodeX, 0, nnodes);
                Marshal.Copy(c_nodesY, rc_nodeY, 0, nnodes);
                for (int i = 0; i < nnodes; i++)
                {
                    Assert.That(rc_nodeX[i], Is.EqualTo(mesh2d_nodesX[i]));
                    Assert.That(rc_nodeY[i], Is.EqualTo(mesh2d_nodesY[i]));
                }

                //Check face nodes
                int fillvalue = -1;
                ierr = wrapper.ionc_get_face_nodes(ref ioncid, ref meshid, ref c_face_nodes, ref nface, ref maxfacenodes, ref fillvalue, ref startIndex);
                Assert.That(ierr, Is.EqualTo(0));
                int[] rc_face_nodes = new int[nface * maxfacenodes];
                Marshal.Copy(c_face_nodes, rc_face_nodes, 0, nface * maxfacenodes);
                int ind = 0;
                for (int i = 0; i < nface; i++)
                {
                    for (int j = 0; j < maxfacenodes; j++)
                    {
                        Assert.That(rc_face_nodes[ind], Is.EqualTo(mesh2d_face_nodes[i, j]));
                        ind += 1;
                    }
                }

            }
            finally
            {
                Marshal.FreeCoTaskMem(c_nodesX);
                Marshal.FreeCoTaskMem(c_nodesY);
                Marshal.FreeCoTaskMem(c_face_nodes);
            }
        }

        private void write1dmesh(
            int ioncid, 
            ref StringBuilder l_networkName,
            ref StringBuilder l_meshname,
            ref IoNetcdfLibWrapper wrapper,
            ref IoNetcdfLibWrapper.interop_charinfo[] meshnodeidsinfo,
            int l_nmeshpoints,
            int l_nedgenodes,
            int l_nBranches,
            int l_nlinks,
            ref double[] l_branchoffset, 
            ref double[] l_branchlength, 
            ref int[] l_branchidx, 
            ref int[] l_sourcenodeid, 
            ref int[] l_targetnodeid, 
            int l_startIndex 
            )
        {
            // Mesh variables
            IntPtr c_branchidx    = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int))    * l_nmeshpoints);
            IntPtr c_edgenodes    = Marshal.AllocCoTaskMem(2* Marshal.SizeOf(typeof(int)) * l_nedgenodes);
            IntPtr c_sourcenodeid = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int))    * l_nBranches);
            IntPtr c_targetnodeid = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int))    * l_nBranches);
            IntPtr c_branchoffset = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nmeshpoints);
            IntPtr c_branchlength = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nBranches);
            // Links variables
            IntPtr c_mesh1indexes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nlinks);
            IntPtr c_mesh2indexes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nlinks);
            IntPtr c_contacttype  = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nlinks);
            try
            {
                string tmpstring;
                int meshid = -1;

                //1. Create: the assumption here is that l_nedgenodes is known (we could move this calculation inside ionc_create_1d_mesh)
                int ierr = wrapper.ionc_create_1d_mesh(ref ioncid, l_networkName.ToString(), ref meshid, l_meshname.ToString(), ref l_nmeshpoints, ref l_nedgenodes);
                Assert.That(ierr, Is.EqualTo(0));

                //2. Create the edge nodes (the algorithm is in gridgeom.dll, not in ionetcdf.dll)
                Marshal.Copy(l_branchidx,    0, c_branchidx, l_nmeshpoints);
                Marshal.Copy(l_sourcenodeid, 0, c_sourcenodeid, l_nBranches);
                Marshal.Copy(l_targetnodeid, 0, c_targetnodeid, l_nBranches);
                Marshal.Copy(l_branchoffset, 0, c_branchoffset, l_nmeshpoints);
                Marshal.Copy(l_branchlength, 0, c_branchlength, l_nBranches);
                
                var gridwrapper = new GridGeomLibWrapper();
                ierr = gridwrapper.ggeo_create_edge_nodes(ref c_branchoffset, ref c_branchlength, ref c_branchidx, ref c_sourcenodeid, ref c_targetnodeid, ref c_edgenodes,ref l_nBranches, ref l_nmeshpoints, ref l_nedgenodes, ref l_startIndex);
                Assert.That(ierr, Is.EqualTo(0));
                
                //4. Write the discretization points
                ierr = wrapper.ionc_put_1d_mesh_discretisation_points(ref ioncid, ref meshid, ref c_branchidx, ref c_branchoffset, ref c_edgenodes, meshnodeidsinfo, ref l_nmeshpoints, ref l_nedgenodes, ref l_startIndex);
                Assert.That(ierr, Is.EqualTo(0));
            }
            finally
            {
                Marshal.FreeCoTaskMem(c_branchidx);
                Marshal.FreeCoTaskMem(c_branchoffset);
                Marshal.FreeCoTaskMem(c_mesh1indexes);
                Marshal.FreeCoTaskMem(c_mesh2indexes);
                Marshal.FreeCoTaskMem(c_edgenodes);
                Marshal.FreeCoTaskMem(c_sourcenodeid);
                Marshal.FreeCoTaskMem(c_targetnodeid);
                Marshal.FreeCoTaskMem(c_contacttype);
            }
        }

        public void write1d2dlinks(
            int ioncid,
            int l_linkmesh1,
            int l_linkmesh2,
            int l_locationType1,
            int l_locationType2,
            ref StringBuilder l_linkmeshname,
            ref IoNetcdfLibWrapper wrapper,
            int l_nlinks,
            ref int[] l_mesh1indexes,
            ref int[] l_mesh2indexes,
            ref int[] l_contacttype,
            int l_startIndex)
        {
            IntPtr c_mesh1indexes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nlinks);
            IntPtr c_mesh2indexes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nlinks);
            IntPtr c_contacttype = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nlinks);
            try
            {
                // 1. Define the mesh with the contacts
                int linkmesh = -1;
                int ierr = wrapper.ionc_def_mesh_contact(ref ioncid, ref linkmesh, l_linkmeshname.ToString(), ref l_nlinks, ref l_linkmesh1,
                    ref l_linkmesh2, ref l_locationType1, ref l_locationType2);
                Assert.That(ierr, Is.EqualTo(0));

                // 2. Put the contacts in thge defined mesh
                Marshal.Copy(l_mesh1indexes, 0, c_mesh1indexes, l_nlinks);
                Marshal.Copy(l_mesh2indexes, 0, c_mesh2indexes, l_nlinks);
                Marshal.Copy(l_contacttype, 0, c_contacttype, l_nlinks);
                IoNetcdfLibWrapper.interop_charinfo[] linksinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nlinks];

                for (int i = 0; i < l_nlinks; i++)
                {
                    string tmpstring = "linkid";
                    tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                    linksinfo[i].ids = tmpstring.ToCharArray();
                    tmpstring = "linklongname";
                    tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                    linksinfo[i].longnames = tmpstring.ToCharArray();
                }

                ierr = wrapper.ionc_put_mesh_contact(ref ioncid, ref linkmesh, ref c_mesh1indexes, ref c_mesh2indexes, ref c_contacttype, linksinfo, ref l_nlinks, ref l_startIndex);
                Assert.That(ierr, Is.EqualTo(0));
            }
            finally
            {
                Marshal.FreeCoTaskMem(c_mesh1indexes);
                Marshal.FreeCoTaskMem(c_mesh2indexes);
                Marshal.FreeCoTaskMem(c_contacttype);
            }
        }


        // writes a network from the arrays 
        private void write1dnetwork(
            int ioncid, 
            int networkid, 
            ref IoNetcdfLibWrapper wrapper, 
            int l_nnodes, 
            int l_nbranches, 
            int l_nGeometry,
            ref IoNetcdfLibWrapper.interop_charinfo[] l_nodesinfo,
            ref IoNetcdfLibWrapper.interop_charinfo[] l_branchinfo,
            ref double[]  l_nodesX,
            ref double[]  l_nodesY,
            ref int[]     l_sourcenodeid,
            ref int[]     l_targetnodeid,
            ref double[]  l_branchlengths,
            ref int[]     l_nbranchgeometrypoints,
            ref double[]  l_geopointsX,
            ref double[]  l_geopointsY,
            ref int[]     l_branch_order)
        {
            IntPtr c_nodesX = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nnodes);
            IntPtr c_nodesY = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nnodes);
            IntPtr c_sourcenodeid = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nbranches);
            IntPtr c_targetnodeid = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nbranches);
            IntPtr c_branchlengths = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nbranches);
            IntPtr c_nbranchgeometrypoints = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nbranches);
            IntPtr c_geopointsX = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nGeometry);
            IntPtr c_geopointsY = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nGeometry);
            IntPtr c_branch_order = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nbranches);
            try
            {
                int ierr = -1;
                string tmpstring;
                //1. Write 1d network network nodes
                Marshal.Copy(l_nodesX, 0, c_nodesX, l_nnodes);
                Marshal.Copy(l_nodesY, 0, c_nodesY, l_nnodes);
            
                ierr = wrapper.ionc_write_1d_network_nodes(ref ioncid, ref networkid, ref c_nodesX, ref c_nodesY,
                    l_nodesinfo, ref l_nnodes);
                Assert.That(ierr, Is.EqualTo(0));

                //2. Write 1d network branches
                Marshal.Copy(l_sourcenodeid, 0, c_sourcenodeid, l_nbranches);
                Marshal.Copy(l_targetnodeid, 0, c_targetnodeid, l_nbranches);
                Marshal.Copy(l_branchlengths, 0, c_branchlengths, l_nbranches);
                Marshal.Copy(l_nbranchgeometrypoints, 0, c_nbranchgeometrypoints, l_nbranches);
                ierr = wrapper.ionc_put_1d_network_branches(ref ioncid, ref networkid, ref c_sourcenodeid,
                    ref c_targetnodeid, l_branchinfo, ref c_branchlengths, ref c_nbranchgeometrypoints, ref l_nbranches, ref startIndex);
                Assert.That(ierr, Is.EqualTo(0));

                //3. Write 1d network geometry
                Marshal.Copy(l_geopointsX, 0, c_geopointsX, l_nGeometry);
                Marshal.Copy(l_geopointsY, 0, c_geopointsY, l_nGeometry);
                ierr = wrapper.ionc_write_1d_network_branches_geometry(ref ioncid, ref networkid, ref c_geopointsX,
                    ref c_geopointsY, ref l_nGeometry);
                Assert.That(ierr, Is.EqualTo(0));

                //4. Define the branch order 
                Marshal.Copy(l_branch_order, 0, c_branch_order, l_nbranches);
                ierr = wrapper.ionc_put_1d_network_branchorder(ref ioncid, ref networkid, ref c_branch_order, ref l_nbranches);
                Assert.That(ierr, Is.EqualTo(0));
            }
            finally
            {
                Marshal.FreeCoTaskMem(c_nodesX);
                Marshal.FreeCoTaskMem(c_nodesY);
                Marshal.FreeCoTaskMem(c_sourcenodeid);
                Marshal.FreeCoTaskMem(c_targetnodeid);
                Marshal.FreeCoTaskMem(c_branchlengths);
                Marshal.FreeCoTaskMem(c_nbranchgeometrypoints);
                Marshal.FreeCoTaskMem(c_geopointsX);
                Marshal.FreeCoTaskMem(c_geopointsY);
                Marshal.FreeCoTaskMem(c_branch_order);
            }
        }

        private void addglobalattributes(int ioncid, ref IoNetcdfLibWrapper wrapper)
        {
            string tmpstring;
            IoNetcdfLibWrapper.interop_metadata metadata;
            tmpstring = "Deltares";
            tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.metadatasize, ' ');
            metadata.institution = tmpstring.ToCharArray();
            tmpstring = "Unknown";
            tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.metadatasize, ' ');
            metadata.source = tmpstring.ToCharArray();
            tmpstring = "Unknown";
            tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.metadatasize, ' ');
            metadata.references = tmpstring.ToCharArray();
            tmpstring = "Unknown";
            tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.metadatasize, ' ');
            metadata.version = tmpstring.ToCharArray();
            tmpstring = "Unknown";
            tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.metadatasize, ' ');
            metadata.modelname = tmpstring.ToCharArray();
            int ierr = wrapper.ionc_add_global_attributes(ref ioncid, ref metadata);
            Assert.That(ierr, Is.EqualTo(0));
        }

        private void getNetworkid(int ioncid, ref int networkid, ref IoNetcdfLibWrapper wrapper)
        {
            //get the number of networks
            int nnumNetworks = -1;
            int ierr = wrapper.ionc_get_number_of_networks(ref ioncid, ref nnumNetworks);
            Assert.That(ierr, Is.EqualTo(0));
            Assert.That(nnumNetworks, Is.EqualTo(1));
            // get the networks ids 
            IntPtr c_networksids = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nnumNetworks);
            ierr = wrapper.ionc_get_network_ids(ref ioncid, ref c_networksids, ref nnumNetworks);
            Assert.That(ierr, Is.EqualTo(0));
            int[] rc_networksids = new int[nnumNetworks];
            Marshal.Copy(c_networksids, rc_networksids, 0, nnumNetworks);
            networkid = rc_networksids[0];
        }

        private void getMeshid(int ioncid, ref int meshid, int meshType, ref IoNetcdfLibWrapper wrapper)
        {
            // get the number meshes 
            int numMeshes = -1;
            int ierr = wrapper.ionc_get_number_of_meshes(ref ioncid, ref meshType, ref numMeshes);
            Assert.That(ierr, Is.EqualTo(0));
            Assert.That(numMeshes, Is.EqualTo(1));
            // get the mesh id
            IntPtr c_meshids = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * numMeshes);
            ierr = wrapper.ionc_ug_get_mesh_ids(ref ioncid, ref meshType, ref c_meshids, ref numMeshes);
            Assert.That(ierr, Is.EqualTo(0));
            int[] rc_meshids = new int[numMeshes];
            Marshal.Copy(c_meshids, rc_meshids, 0, numMeshes);
            meshid = rc_meshids[0];
        }

        // Create the netcdf files
        [Test]
        [NUnit.Framework.Category("UGRIDTests")]
        public void create1dUGridNetworkAndMeshNetcdf()
        {
            //1. Create a netcdf file 
            int ioncid = 0; //file variable 
            int mode = 1; //create in write mode
            var ierr = -1;
            string tmpstring; //temporary string for several operations
            string c_path = TestHelper.TestDirectoryPath() + @"\write1d.nc";
            TestHelper.DeleteIfExists(c_path);
            Assert.IsFalse(File.Exists(c_path));
            var wrapper = new IoNetcdfLibWrapper();

            //2. make a local copy of the variables 
            int l_nnodes = nNodes;
            int l_nbranches = nBranches;
            int l_nGeometry = nGeometry;
            double[] l_nodesX = nodesX;
            double[] l_nodesY = nodesY;
            int[] l_sourcenodeid = sourcenodeid;
            int[] l_targetnodeid = targetnodeid;
            double[] l_branchlengths = branchlengths;
            int[] l_nbranchgeometrypoints = nbranchgeometrypoints;
            double[] l_geopointsX = geopointsX;
            double[] l_geopointsY = geopointsY;
            int[] l_branch_order = branch_order;

            IoNetcdfLibWrapper.interop_charinfo[] l_nodesinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nnodes];
            IoNetcdfLibWrapper.interop_charinfo[] l_branchinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nbranches];

            for (int i = 0; i < l_nnodes; i++)
            {
                tmpstring = nodesids[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                l_nodesinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = nodeslongNames[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                l_nodesinfo[i].longnames = tmpstring.ToCharArray();
            }

            for (int i = 0; i < l_nbranches; i++)
            {
                tmpstring = branchids[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                l_branchinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = branchlongNames[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                l_branchinfo[i].longnames = tmpstring.ToCharArray();
            }

            //4. Create the file, will not add any dataset 
            ierr = wrapper.ionc_create(c_path, ref mode, ref ioncid);
            Assert.That(ierr, Is.EqualTo(0));
            Assert.IsTrue(File.Exists(c_path));

            //5. For reading the grid later on we need to add metadata to the netcdf file. 
            //   The function ionc_add_global_attributes adds to the netCDF file the UGRID convention
            addglobalattributes(ioncid, ref wrapper);

            //6. Create a 1d network
            int networkid = -1;
            StringBuilder l_networkName = new StringBuilder(networkName);
            ierr = wrapper.ionc_create_1d_network(ref ioncid, ref networkid, l_networkName.ToString(), ref nNodes, ref nBranches,
                ref nGeometry);
            Assert.That(ierr, Is.EqualTo(0));

            StringBuilder l_meshname = new StringBuilder(meshname);
            int l_nmeshpoints = nmeshpoints;
            int l_nedgenodes = nedgenodes;
            int l_nBranches = nBranches;
            int l_nlinks = nlinks;
            double[] l_branchoffset = offset;
            double[] l_branchlength = branchlengths;
            int[] l_branchidx = branchidx;
            int l_startIndex = startIndex;

            //3. Create the node branchidx, offsets, meshnodeidsinfo
            IoNetcdfLibWrapper.interop_charinfo[] meshnodeidsinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nmeshpoints];
            for (int i = 0; i < l_nmeshpoints; i++)
            {
                tmpstring = meshnodeids[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                meshnodeidsinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = meshnodelongnames[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                meshnodeidsinfo[i].longnames = tmpstring.ToCharArray();
            }

            //7. Write 1d network and mesh
            write1dnetwork(ioncid,
                networkid,
                ref wrapper,
                l_nnodes,
                l_nbranches,
                l_nGeometry,
                ref l_nodesinfo,
                ref l_branchinfo,
                ref l_nodesX,
                ref l_nodesY,
                ref l_sourcenodeid,
                ref l_targetnodeid,
                ref l_branchlengths,
                ref l_nbranchgeometrypoints,
                ref l_geopointsX,
                ref l_geopointsY,
                ref l_branch_order
                );

            write1dmesh(
                ioncid,
                ref l_networkName,
                ref l_meshname,
                ref wrapper,
                ref meshnodeidsinfo,
                l_nmeshpoints,
                l_nedgenodes,
                l_nBranches,
                l_nlinks,
                ref l_branchoffset,
                ref l_branchlength,
                ref l_branchidx,
                ref l_sourcenodeid,
                ref l_targetnodeid,
                l_startIndex);

            //8. Close the file
            ierr = wrapper.ionc_close(ref ioncid);
            Assert.That(ierr, Is.EqualTo(0));
        }

        ////// read the netcdf file created in the test above
        [Test]
        [NUnit.Framework.Category("UGRIDTests")]
        public void read1dUGRIDNetcdf()
        {
            IntPtr c_meshidsfromnetworkid = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)));
            try
            {
                //1. Open a netcdf file 
                string c_path = TestHelper.TestDirectoryPath() + @"\write1d.nc";
                Assert.IsTrue(File.Exists(c_path));
                int ioncid = 0; //file variable 
                int mode = 0; //create in read mode
                var wrapper = new IoNetcdfLibWrapper();
                var ierr = wrapper.ionc_open(c_path, ref mode, ref ioncid, ref iconvtype, ref convversion);
                Assert.That(ierr, Is.EqualTo(0));

                //2. Get the 1D network and mesh ids
                // network
                int networkid = -1;
                getNetworkid(ioncid, ref networkid, ref wrapper);
                Assert.That(networkid, Is.EqualTo(1));

                // 1d mesh mesh
                int meshType = 1;
                int meshid = -1;
                getMeshid(ioncid, ref meshid, meshType, ref wrapper);
                Assert.That(meshid, Is.EqualTo(1));

                ierr = wrapper.ionc_get_1d_mesh_id(ref ioncid, ref meshid);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.That(meshid, Is.EqualTo(1));

                //3. create local variables
                int l_nnodes = -1;
                int l_nbranches = -1;
                int l_nGeometry = -1;
                int l_startIndex = startIndex;
                StringBuilder l_networkName = new StringBuilder(IoNetcdfLibWrapper.LibDetails.MAXSTRLEN);
                IoNetcdfLibWrapper.interop_charinfo[] l_nodesinfo = new IoNetcdfLibWrapper.interop_charinfo[nNodes];
                IoNetcdfLibWrapper.interop_charinfo[] l_branchinfo = new IoNetcdfLibWrapper.interop_charinfo[nNodes];
                double[] l_nodesX = new double[nNodes];
                double[] l_nodesY = new double[nNodes];
                int[] l_sourcenodeid = new int[nBranches];
                int[] l_targetnodeid = new int[nBranches];
                double[] l_branchlengths = new double[nBranches];
                int[] l_nbranchgeometrypoints = new int[nBranches];
                double[] l_geopointsX = new double[nGeometry];
                double[] l_geopointsY = new double[nGeometry];
                int[] l_branch_order = new int[nBranches];

                read1dnetwork(
                    ioncid,
                    networkid,
                    ref wrapper,
                    ref l_nnodes,
                    ref l_nbranches,
                    ref l_nGeometry,
                    l_startIndex,
                    ref l_networkName,
                    ref l_nodesinfo,
                    ref l_branchinfo,
                    ref l_nodesX,
                    ref l_nodesY,
                    ref l_sourcenodeid,
                    ref l_targetnodeid,
                    ref l_branchlengths,
                    ref l_nbranchgeometrypoints,
                    ref l_geopointsX,
                    ref l_geopointsY,
                    ref l_branch_order);

                //4. check the read values
                Assert.That(l_networkName.ToString().Trim(), Is.EqualTo(networkName));
                Assert.That(l_nnodes, Is.EqualTo(l_nnodes));
                Assert.That(l_nbranches, Is.EqualTo(nBranches));
                Assert.That(l_nGeometry, Is.EqualTo(nGeometry));

                for (int i = 0; i < l_nnodes; i++)
                {
                    string tmpstring = new string(l_nodesinfo[i].ids);
                    Assert.That(tmpstring.Trim(), Is.EqualTo(nodesids[i]));
                    tmpstring = new string(l_nodesinfo[i].longnames);
                    Assert.That(tmpstring.Trim(), Is.EqualTo(nodeslongNames[i]));
                    Assert.That(l_nodesX[i], Is.EqualTo(nodesX[i]));
                    Assert.That(l_nodesY[i], Is.EqualTo(nodesY[i]));
                }

                for (int i = 0; i < l_nbranches; i++)
                {
                    string tmpstring = new string(l_branchinfo[i].ids);
                    Assert.That(tmpstring.Trim(), Is.EqualTo(branchids[i]));
                    tmpstring = new string(l_branchinfo[i].longnames);
                    Assert.That(tmpstring.Trim(), Is.EqualTo(branchlongNames[i]));
                    Assert.That(l_targetnodeid[i], Is.EqualTo(targetnodeid[i]));
                    Assert.That(l_sourcenodeid[i], Is.EqualTo(sourcenodeid[i]));
                    Assert.That(l_branchlengths[i], Is.EqualTo(branchlengths[i]));
                    Assert.That(l_nbranchgeometrypoints[i], Is.EqualTo(nbranchgeometrypoints[i]));
                    Assert.That(l_branch_order[i], Is.EqualTo(branch_order[i]));
                }

                for (int i = 0; i < l_nGeometry; i++)
                {
                    Assert.That(l_geopointsX[i], Is.EqualTo(geopointsX[i]));
                    Assert.That(l_geopointsY[i], Is.EqualTo(geopointsY[i]));
                }

                // LC: REFACTOR NEEDED!
                check1dmesh(ioncid, meshid, ref wrapper);

                //5. count the meshes associated with this network
                int nmeshids = -1;
                ierr = wrapper.ionc_count_mesh_ids_from_network_id(ref ioncid, ref networkid, ref nmeshids);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.That(nmeshids, Is.EqualTo(1));

                int[] meshidsfromnetworkid = new int[nmeshids];
                ierr = wrapper.ionc_get_mesh_ids_from_network_id(ref ioncid, ref networkid, ref nmeshids, ref c_meshidsfromnetworkid);
                Assert.That(ierr, Is.EqualTo(0));
                Marshal.Copy(c_meshidsfromnetworkid, meshidsfromnetworkid, 0, nmeshids);
                Assert.That(meshidsfromnetworkid[0], Is.EqualTo(1));

                //6. get the network id from the mesh id
                networkid = -1;
                ierr = wrapper.ionc_get_network_id_from_mesh_id(ref ioncid, ref meshid, ref networkid);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.That(networkid, Is.EqualTo(1));

                //7. Close the file
                ierr = wrapper.ionc_close(ref ioncid);
                Assert.That(ierr, Is.EqualTo(0));
            }
            finally
            {
                Marshal.FreeCoTaskMem(c_meshidsfromnetworkid);
            }
        }

        // Deltashell creates a new file to write the 1d geometry and mesh as in the first test create1dUGRIDNetcdf
        // and clones the 2d mesh data read from a file produced by RGFgrid. 
        [Test]
        [NUnit.Framework.Category("UGRIDTests")]
        public void Clones2dMesh()
        {
            var wrapper = new IoNetcdfLibWrapper();

            //1. RGF grid creates a 2d mesh. The info is in memory, here simulated by opening a file containing a mesh2d
            // and by reading all data in
            string sourcetwod_path = TestHelper.CreateLocalCopy("Custom_Ugrid.nc");
            Assert.IsTrue(File.Exists(sourcetwod_path));
            int sourcetwodioncid = -1; //file id 
            int sourcetwomode = 0; //read mode
            int ierr = wrapper.ionc_open(sourcetwod_path, ref sourcetwomode, ref sourcetwodioncid, ref iconvtype,
                ref convversion);
            Assert.That(ierr, Is.EqualTo(0));

            // make a local copy of the variables 
            int l_nnodes = nNodes;
            int l_nbranches = nBranches;
            int l_nGeometry = nGeometry;
            double[] l_nodesX = nodesX;
            double[] l_nodesY = nodesY;
            int[] l_sourcenodeid = sourcenodeid;
            int[] l_targetnodeid = targetnodeid;
            double[] l_branchlengths = branchlengths;
            int[] l_nbranchgeometrypoints = nbranchgeometrypoints;
            double[] l_geopointsX = geopointsX;
            double[] l_geopointsY = geopointsY;
            int[] l_branch_order = branch_order;

            IoNetcdfLibWrapper.interop_charinfo[] l_nodesinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nnodes];
            IoNetcdfLibWrapper.interop_charinfo[] l_branchinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nbranches];

            for (int i = 0; i < l_nnodes; i++)
            {
                string tmpstring = "";
                tmpstring = nodesids[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                l_nodesinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = nodeslongNames[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                l_nodesinfo[i].longnames = tmpstring.ToCharArray();
            }

            for (int i = 0; i < l_nbranches; i++)
            {
                string tmpstring = "";
                tmpstring = branchids[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                l_branchinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = branchlongNames[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                l_branchinfo[i].longnames = tmpstring.ToCharArray();
            }

            //2. Now we create a new empty file where to save 1d and 2d meshes
            int targetioncid = -1; //file id  
            int targetmode = 1; //create in write mode
            string target_path = TestHelper.TestDirectoryPath() + "/target.nc";
            TestHelper.DeleteIfExists(target_path);
            Assert.IsFalse(File.Exists(target_path));

            //3. Create the file
            ierr = wrapper.ionc_create(target_path, ref targetmode, ref targetioncid);
            Assert.That(ierr, Is.EqualTo(0));
            Assert.IsTrue(File.Exists(target_path));

            //4. Add global attributes in the file
            addglobalattributes(targetioncid, ref wrapper);

            //5. Get the id of the 2d mesh in the RGF grid file (Custom_Ugrid.nc)
            int meshType = 2;
            int sourcemesh2d = -1;
            getMeshid(sourcetwodioncid, ref sourcemesh2d, meshType, ref wrapper);
            Assert.That(sourcemesh2d, Is.EqualTo(1));

            //6. Create 1d geometry and mesh in the new file (target.nc)
            int networkid = -1;
            ierr = wrapper.ionc_create_1d_network(ref targetioncid, ref networkid, networkName, ref nNodes,
                ref nBranches, ref nGeometry);
            Assert.That(ierr, Is.EqualTo(0));

            StringBuilder l_networkName = new StringBuilder(networkName);
            StringBuilder l_meshname = new StringBuilder(meshname);
            StringBuilder l_linkmeshname = new StringBuilder(linkmeshname);
            int l_nmeshpoints = nmeshpoints;
            int l_nedgenodes = nedgenodes;
            int l_nBranches = nBranches;
            int l_nlinks = nlinks;
            double[] l_branchoffset = offset;
            double[] l_branchlength = branchlengths;
            int[] l_branchidx = branchidx;
            int l_startIndex = startIndex;

            //3. Create the node branchidx, offsets, meshnodeidsinfo
            IoNetcdfLibWrapper.interop_charinfo[] meshnodeidsinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nmeshpoints];
            for (int i = 0; i < l_nmeshpoints; i++)
            {
                string tmpstring = "";
                tmpstring = meshnodeids[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                meshnodeidsinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = meshnodelongnames[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                meshnodeidsinfo[i].longnames = tmpstring.ToCharArray();
            }


            //6. Write the 1d data in the new file (1d geometry, mesh and links)
            write1dnetwork(targetioncid,
                networkid,
                ref wrapper,
                l_nnodes,
                l_nbranches,
                l_nGeometry,
                ref l_nodesinfo,
                ref l_branchinfo,
                ref l_nodesX,
                ref l_nodesY,
                ref l_sourcenodeid,
                ref l_targetnodeid,
                ref l_branchlengths,
                ref l_nbranchgeometrypoints,
                ref l_geopointsX,
                ref l_geopointsY,
                ref l_branch_order);

            write1dmesh(
                targetioncid,
                ref l_networkName,
                ref l_meshname,
                ref wrapper,
                ref meshnodeidsinfo,
                l_nmeshpoints,
                l_nedgenodes,
                l_nBranches,
                l_nlinks,
                ref l_branchoffset,
                ref l_branchlength,
                ref l_branchidx,
                ref l_sourcenodeid,
                ref l_targetnodeid,
                l_startIndex);

            //7. Clone the 2d mesh definitions in the new file
            int target2dmesh = -1;
            ierr = wrapper.ionc_clone_mesh_definition(ref sourcetwodioncid, ref targetioncid, ref sourcemesh2d,
                ref target2dmesh);
            Assert.That(ierr, Is.EqualTo(0));

            //8. Clone the 2d mesh data
            ierr = wrapper.ionc_clone_mesh_data(ref sourcetwodioncid, ref targetioncid, ref sourcemesh2d,
                ref target2dmesh);
            Assert.That(ierr, Is.EqualTo(0));

            //9. Close all files 
            ierr = wrapper.ionc_close(ref sourcetwodioncid);
            Assert.That(ierr, Is.EqualTo(0));
            ierr = wrapper.ionc_close(ref targetioncid);
            Assert.That(ierr, Is.EqualTo(0));

            //10. Now open the file with cloned meshes and check if the data written there are correct
            targetioncid = -1; //file id  
            targetmode = 0; //open in write mode
            ierr = wrapper.ionc_open(target_path, ref targetmode, ref targetioncid, ref iconvtype, ref convversion);
            Assert.That(ierr, Is.EqualTo(0));

            //11. Check 2 meshes are present
            int nmesh = -1;
            ierr = wrapper.ionc_get_mesh_count(ref targetioncid, ref nmesh);
            Assert.That(ierr, Is.EqualTo(0));
            Assert.That(nmesh, Is.EqualTo(2));

            //12. Get the mesh  and network ids
            // network
            int source1dnetwork = -1;
            getNetworkid(targetioncid, ref source1dnetwork, ref wrapper);
            Assert.That(networkid, Is.EqualTo(1));

            // 1d mesh mesh
            meshType = 1;
            int source1dmesh = -1;
            getMeshid(targetioncid, ref source1dmesh, meshType, ref wrapper);
            Assert.That(source1dmesh, Is.EqualTo(1));

            // 2d mesh mesh
            meshType = 2;
            int source2dmesh = -1;
            getMeshid(targetioncid, ref source2dmesh, meshType, ref wrapper);
            Assert.That(source2dmesh, Is.EqualTo(2));

            StringBuilder ll_networkName = new StringBuilder(IoNetcdfLibWrapper.LibDetails.MAXSTRLEN);

            //13. Check all 1d and 2d data
            read1dnetwork(
                targetioncid,
                networkid,
                ref wrapper,
                ref l_nnodes,
                ref l_nbranches,
                ref l_nGeometry,
                l_startIndex,
                ref ll_networkName,
                ref l_nodesinfo,
                ref l_branchinfo,
                ref l_nodesX,
                ref l_nodesY,
                ref l_sourcenodeid,
                ref l_targetnodeid,
                ref l_branchlengths,
                ref l_nbranchgeometrypoints,
                ref l_geopointsX,
                ref l_geopointsY,
                ref l_branch_order);

            //4. check the read values
            Assert.That(l_networkName.ToString().Trim(), Is.EqualTo(networkName));
            Assert.That(l_nnodes, Is.EqualTo(l_nnodes));
            Assert.That(l_nbranches, Is.EqualTo(nBranches));
            Assert.That(l_nGeometry, Is.EqualTo(nGeometry));

            for (int i = 0; i < l_nnodes; i++)
            {
                string tmpstring = new string(l_nodesinfo[i].ids);
                Assert.That(tmpstring.Trim(), Is.EqualTo(nodesids[i]));
                tmpstring = new string(l_nodesinfo[i].longnames);
                Assert.That(tmpstring.Trim(), Is.EqualTo(nodeslongNames[i]));
                Assert.That(l_nodesX[i], Is.EqualTo(l_nodesX[i]));
                Assert.That(l_nodesY[i], Is.EqualTo(nodesY[i]));
            }

            for (int i = 0; i < l_nbranches; i++)
            {
                string tmpstring = new string(l_branchinfo[i].ids);
                Assert.That(tmpstring.Trim(), Is.EqualTo(branchids[i]));
                tmpstring = new string(l_branchinfo[i].longnames);
                Assert.That(tmpstring.Trim(), Is.EqualTo(branchlongNames[i]));
                Assert.That(l_targetnodeid[i], Is.EqualTo(targetnodeid[i]));
                Assert.That(l_sourcenodeid[i], Is.EqualTo(sourcenodeid[i]));
                Assert.That(l_branchlengths[i], Is.EqualTo(branchlengths[i]));
                Assert.That(l_nbranchgeometrypoints[i], Is.EqualTo(nbranchgeometrypoints[i]));
                Assert.That(l_branch_order[i], Is.EqualTo(branch_order[i]));
            }

            for (int i = 0; i < l_nGeometry; i++)
            {
                Assert.That(l_geopointsX[i], Is.EqualTo(geopointsX[i]));
                Assert.That(l_geopointsY[i], Is.EqualTo(geopointsY[i]));
            }

            //LC refactor needed also for this part!
            check1dmesh(targetioncid, source1dmesh, ref wrapper);
            check2dmesh(targetioncid, source2dmesh, ref wrapper);

            //14. Close the file
            ierr = wrapper.ionc_close(ref targetioncid);
            Assert.That(ierr, Is.EqualTo(0));
        }

        // Create a standalone network
        [Test]
        [NUnit.Framework.Category("UGRIDTests")]
        public void create1dNetwork()
        {
            //1. Create a netcdf file 
            int ioncid = 0; //file variable 
            int mode = 1; //create in write mode
            var ierr = -1;
            string c_path = TestHelper.TestDirectoryPath() + @"\write1dNetwork.nc";
            TestHelper.DeleteIfExists(c_path);
            Assert.IsFalse(File.Exists(c_path));
            var wrapper = new IoNetcdfLibWrapper();

            // make a local copy of the variables 
            int l_nnodes = nNodes;
            int l_nbranches = nBranches;
            int l_nGeometry = nGeometry;
            double[] l_nodesX = nodesX;
            double[] l_nodesY = nodesY;
            int[] l_sourcenodeid = sourcenodeid;
            int[] l_targetnodeid = targetnodeid;
            double[] l_branchlengths = branchlengths;
            int[] l_nbranchgeometrypoints = nbranchgeometrypoints;
            double[] l_geopointsX = geopointsX;
            double[] l_geopointsY = geopointsY;
            int[] l_branch_order = branch_order;


            IoNetcdfLibWrapper.interop_charinfo[] l_nodesinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nnodes];
            IoNetcdfLibWrapper.interop_charinfo[] l_branchinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nbranches];

            for (int i = 0; i < l_nnodes; i++)
            {
                string tmpstring = nodesids[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                l_nodesinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = nodeslongNames[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                l_nodesinfo[i].longnames = tmpstring.ToCharArray();
            }

            for (int i = 0; i < l_nbranches; i++)
            {
                string tmpstring = branchids[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                l_branchinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = branchlongNames[i];
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                l_branchinfo[i].longnames = tmpstring.ToCharArray();
            }

            //2. Create the file, will not add any dataset 
            ierr = wrapper.ionc_create(c_path, ref mode, ref ioncid);
            Assert.That(ierr, Is.EqualTo(0));
            Assert.IsTrue(File.Exists(c_path));

            //3. For reading the grid later on we need to add metadata to the netcdf file. 
            //   The function ionc_add_global_attributes adds to the netCDF file the UGRID convention
            addglobalattributes(ioncid, ref wrapper);

            //4. Create a 1d network
            int networkid = -1;
            ierr = wrapper.ionc_create_1d_network(ref ioncid, ref networkid, networkName, ref nNodes,
                ref nBranches, ref nGeometry);
            Assert.That(ierr, Is.EqualTo(0));

            //5. Write 1d network and mesh
            write1dnetwork(ioncid,
                networkid,
                ref wrapper,
                l_nnodes,
                l_nbranches,
                l_nGeometry,
                ref l_nodesinfo,
                ref l_branchinfo,
                ref l_nodesX,
                ref l_nodesY,
                ref l_sourcenodeid,
                ref l_targetnodeid,
                ref l_branchlengths,
                ref l_nbranchgeometrypoints,
                ref l_geopointsX,
                ref l_geopointsY,
                ref l_branch_order);
            //6. Close the file
            ierr = wrapper.ionc_close(ref ioncid);
            Assert.That(ierr, Is.EqualTo(0));
        }

        // read the standalone network
        [Test]
        [NUnit.Framework.Category("UGRIDTests")]
        public void read1dNetwork()
        {
            //1. Open a netcdf file 
            string c_path = TestHelper.TestDirectoryPath() + @"\write1dNetwork.nc";
            Assert.IsTrue(File.Exists(c_path));
            int ioncid = 0; //file variable 
            int mode = 0; //create in read mode
            var wrapper = new IoNetcdfLibWrapper();
            var ierr = wrapper.ionc_open(c_path, ref mode, ref ioncid, ref iconvtype, ref convversion);
            Assert.That(ierr, Is.EqualTo(0));

            //2. Get the 1D network and mesh ids
            int networkid = -1;
            ierr = wrapper.ionc_get_1d_network_id(ref ioncid, ref networkid);
            Assert.That(ierr, Is.EqualTo(0));
            Assert.That(networkid, Is.EqualTo(1));
            int meshid = -1;
            ierr = wrapper.ionc_get_1d_mesh_id(ref ioncid, ref meshid);

            // Mesh should not be found!
            Assert.That(ierr, Is.EqualTo(-1));
            Assert.That(meshid, Is.EqualTo(-1));

            //3. Check if all 1d data written in the file are correct
            int l_nnodes = -1;
            int l_nbranches = -1;
            int l_nGeometry = -1;
            int l_startIndex = startIndex;
            StringBuilder l_networkName = new StringBuilder(IoNetcdfLibWrapper.LibDetails.MAXSTRLEN);
            IoNetcdfLibWrapper.interop_charinfo[] l_nodesinfo = new IoNetcdfLibWrapper.interop_charinfo[nNodes];
            IoNetcdfLibWrapper.interop_charinfo[] l_branchinfo = new IoNetcdfLibWrapper.interop_charinfo[nNodes];
            double[] l_nodesX = new double[nNodes];
            double[] l_nodesY = new double[nNodes];
            int[] l_sourcenodeid = new int[nBranches];
            int[] l_targetnodeid = new int[nBranches];
            double[] l_branchlengths = new double[nBranches];
            int[] l_nbranchgeometrypoints = new int[nBranches];
            double[] l_geopointsX = new double[nGeometry];
            double[] l_geopointsY = new double[nGeometry];
            int[] l_branch_order = new int[nBranches];

            read1dnetwork(
                ioncid,
                networkid,
                ref wrapper,
                ref l_nnodes,
                ref l_nbranches,
                ref l_nGeometry,
                l_startIndex,
                ref l_networkName,
                ref l_nodesinfo,
                ref l_branchinfo,
                ref l_nodesX,
                ref l_nodesY,
                ref l_sourcenodeid,
                ref l_targetnodeid,
                ref l_branchlengths,
                ref l_nbranchgeometrypoints,
                ref l_geopointsX,
                ref l_geopointsY,
                ref l_branch_order);

            //4. check the read values
            Assert.That(l_networkName.ToString().Trim(), Is.EqualTo(networkName));
            Assert.That(l_nnodes, Is.EqualTo(l_nnodes));
            Assert.That(l_nbranches, Is.EqualTo(nBranches));
            Assert.That(l_nGeometry, Is.EqualTo(nGeometry));

            for (int i = 0; i < l_nnodes; i++)
            {
                string tmpstring = new string(l_nodesinfo[i].ids);
                Assert.That(tmpstring.Trim(), Is.EqualTo(nodesids[i]));
                tmpstring = new string(l_nodesinfo[i].longnames);
                Assert.That(tmpstring.Trim(), Is.EqualTo(nodeslongNames[i]));
                Assert.That(l_nodesX[i], Is.EqualTo(l_nodesX[i]));
                Assert.That(l_nodesY[i], Is.EqualTo(nodesY[i]));
            }

            for (int i = 0; i < l_nbranches; i++)
            {
                string tmpstring = new string(l_branchinfo[i].ids);
                Assert.That(tmpstring.Trim(), Is.EqualTo(branchids[i]));
                tmpstring = new string(l_branchinfo[i].longnames);
                Assert.That(tmpstring.Trim(), Is.EqualTo(branchlongNames[i]));
                Assert.That(l_targetnodeid[i], Is.EqualTo(targetnodeid[i]));
                Assert.That(l_sourcenodeid[i], Is.EqualTo(sourcenodeid[i]));
                Assert.That(l_branchlengths[i], Is.EqualTo(branchlengths[i]));
                Assert.That(l_nbranchgeometrypoints[i], Is.EqualTo(nbranchgeometrypoints[i]));
                Assert.That(l_branch_order[i], Is.EqualTo(branch_order[i]));
            }

            for (int i = 0; i < l_nGeometry; i++)
            {
                Assert.That(l_geopointsX[i], Is.EqualTo(geopointsX[i]));
                Assert.That(l_geopointsY[i], Is.EqualTo(geopointsY[i]));
            }
            //4. Close the file
            ierr = wrapper.ionc_close(ref ioncid);
        }

        /* 
        1)	Allocates the arrays defining the network (nodes, branches, geometry points, and all ids/longnames). Here I used 1000001 nodes, 1000000 branches, 3000000 geometry points 
        2)	Opens the first netcdf file “sewer_system.nc”
        3)	Creates and writes the network
        4)	Closes the file
        5)	Opens the file
        6)	Reads the arrays back in
        7)	Closes the file
        8)	Allocates the arrays for the augmented network(1000001 nodes, 1000002 branches, 3000003 geometry points)
        9)	Copies the values read from the “sewer_system.nc” in the new arrays
        10)	Adds the “stranger” node, branch and geometry points
        11)	Opens a second netcdf file “LargeSewerSystemSecondTest.nc”
        12)	Creates and writes the network “sewer_system_with_the_stranger”
        13)	Closes the file
        */
        [Test]
        [NUnit.Framework.Category("UGRIDTests")]
        public void LargeSewerSystem()
        {
            int stackSize = 1024 * 1024 * 16; //LC: in C# the default stack size is 1 MB, increase it to something larger for this test!
            Thread th = new Thread(() =>
            {
                //1. Allocates the arrays defining the network 
                int firstCaseNumberOfNodes = 100000; //5000 limit win64, without stack increase
                string tmpstring;

                int l_nnodes = firstCaseNumberOfNodes + 1;
                int l_nbranches = firstCaseNumberOfNodes;
                int l_nGeometry = firstCaseNumberOfNodes * 3;
                double[] l_nodesX = new double[l_nnodes];
                double[] l_nodesY = new double[l_nnodes];
                double[] l_branchlengths = new double[l_nbranches];
                int[] l_nbranchgeometrypoints = new int[l_nbranches];

                int[] l_sourcenodeid = new int[l_nbranches];
                int[] l_targetnodeid = new int[l_nbranches];
                int[] l_branch_order = new int[l_nbranches];

                double[] l_geopointsX = new double[l_nGeometry];
                double[] l_geopointsY = new double[l_nGeometry];

                IoNetcdfLibWrapper.interop_charinfo[] l_nodesinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nnodes];
                IoNetcdfLibWrapper.interop_charinfo[] l_branchinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nbranches];

                for (int i = 0; i < l_nnodes; i++)
                {
                    tmpstring = "node_id_" + i.ToString();
                    tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                    l_nodesinfo[i].ids = tmpstring.ToCharArray();
                    tmpstring = "node_longname_" + i.ToString();
                    tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                    l_nodesinfo[i].longnames = tmpstring.ToCharArray();
                    l_nodesX[i] = i;
                    l_nodesY[i] = i;
                }

                for (int i = 0; i < l_nbranches; i++)
                {
                    tmpstring = "branch_id_" + i.ToString();
                    tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                    l_branchinfo[i].ids = tmpstring.ToCharArray();
                    tmpstring = "branch_longname_" + i.ToString();
                    tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                    l_branchinfo[i].longnames = tmpstring.ToCharArray();
                    l_branchlengths[i] = 1.414;
                    l_nbranchgeometrypoints[i] = 1;
                    l_sourcenodeid[i] = i;
                    l_targetnodeid[i] = i + 1;
                    l_branch_order[i] = i;
                }

                for (int i = 0; i < l_nGeometry; i++)
                {
                    l_geopointsX[i] = i + 1.0 / 2.0;
                    l_geopointsX[i] = i + 1.0 / 2.0;
                }

                //2
                int ioncid = 0;     // file variable 
                int mode = 1;       // create in write mode
                var ierr = -1;
                string c_path = TestHelper.TestDirectoryPath() + @"\LargeSewerSystem.nc";
                TestHelper.DeleteIfExists(c_path);
                Assert.IsFalse(File.Exists(c_path));
                var wrapper = new IoNetcdfLibWrapper();

                ierr = wrapper.ionc_create(c_path, ref mode, ref ioncid);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.IsTrue(File.Exists(c_path));

                //The function ionc_add_global_attributes adds to the netCDF file the UGRID convention
                addglobalattributes(ioncid, ref wrapper);

                string l_network_name = "sewer_system";
                int networkid = -1;
                ierr = wrapper.ionc_create_1d_network(ref ioncid, ref networkid, l_network_name, ref l_nnodes,
                    ref l_nbranches, ref l_nGeometry);
                Assert.That(ierr, Is.EqualTo(0));

                //3. Write 1d network and mesh
                write1dnetwork(ioncid,
                    networkid,
                    ref wrapper,
                    l_nnodes,
                    l_nbranches,
                    l_nGeometry,
                    ref l_nodesinfo,
                    ref l_branchinfo,
                    ref l_nodesX,
                    ref l_nodesY,
                    ref l_sourcenodeid,
                    ref l_targetnodeid,
                    ref l_branchlengths,
                    ref l_nbranchgeometrypoints,
                    ref l_geopointsX,
                    ref l_geopointsY,
                    ref l_branch_order);

                //4. Close the file
                ierr = wrapper.ionc_close(ref ioncid);
                Assert.That(ierr, Is.EqualTo(0));

                //5. Open the file in readmode
                mode = 0; //read
                ierr = wrapper.ionc_open(c_path, ref mode, ref ioncid, ref iconvtype, ref convversion);
                Assert.That(ierr, Is.EqualTo(0));

                //Get the 1D network id
                getNetworkid(ioncid, ref networkid, ref wrapper);
                Assert.That(networkid, Is.EqualTo(1));

                //Create local variables
                int l_startIndex = startIndex;
                StringBuilder l_networkName = new StringBuilder(IoNetcdfLibWrapper.LibDetails.MAXSTRLEN);
                //6. Read the arrays back in
                read1dnetwork(
                    ioncid,
                    networkid,
                    ref wrapper,
                    ref l_nnodes,
                    ref l_nbranches,
                    ref l_nGeometry,
                    l_startIndex,
                    ref l_networkName,
                    ref l_nodesinfo,
                    ref l_branchinfo,
                    ref l_nodesX,
                    ref l_nodesY,
                    ref l_sourcenodeid,
                    ref l_targetnodeid,
                    ref l_branchlengths,
                    ref l_nbranchgeometrypoints,
                    ref l_geopointsX,
                    ref l_geopointsY,
                    ref l_branch_order);

                //7. Close the file
                ierr = wrapper.ionc_close(ref ioncid);
                Assert.That(ierr, Is.EqualTo(0));

                //8. Allocate arrays for augmented network
                int secondCaseNumberOfNodes = firstCaseNumberOfNodes + 1;

                int sl_nnodes = secondCaseNumberOfNodes + 1;
                int sl_nbranches = secondCaseNumberOfNodes;
                int sl_nGeometry = secondCaseNumberOfNodes * 3;
                double[] sl_nodesX = new double[sl_nnodes];
                double[] sl_nodesY = new double[sl_nnodes];


                double[] sl_branchlengths = new double[sl_nbranches];
                int[] sl_nbranchgeometrypoints = new int[sl_nbranches];
                int[] sl_sourcenodeid = new int[sl_nbranches];
                int[] sl_targetnodeid = new int[sl_nbranches];
                int[] sl_branch_order = new int[sl_nbranches];

                double[] sl_geopointsX = new double[sl_nGeometry];
                double[] sl_geopointsY = new double[sl_nGeometry];

                IoNetcdfLibWrapper.interop_charinfo[] sl_nodesinfo = new IoNetcdfLibWrapper.interop_charinfo[sl_nnodes];
                IoNetcdfLibWrapper.interop_charinfo[] sl_branchinfo = new IoNetcdfLibWrapper.interop_charinfo[sl_nbranches];

                //9. Copy the values of sewer_system in the arrays of the augmented network
                for (int i = 0; i < l_nnodes; i++)
                {
                    sl_nodesinfo[i].ids = l_nodesinfo[i].ids;
                    sl_nodesinfo[i].longnames = l_nodesinfo[i].longnames;
                }

                for (int i = 0; i < l_nbranches; i++)
                {
                    sl_branchinfo[i].ids = l_branchinfo[i].ids;
                    sl_branchinfo[i].longnames = l_branchinfo[i].longnames;
                }

                Array.Copy(l_nodesX, sl_nodesX, l_nodesX.Length);
                Array.Copy(l_nodesY, sl_nodesY, l_nodesY.Length);
                Array.Copy(l_branchlengths, sl_branchlengths, l_branchlengths.Length);
                Array.Copy(l_nbranchgeometrypoints, sl_nbranchgeometrypoints, l_nbranchgeometrypoints.Length);
                Array.Copy(l_sourcenodeid, sl_sourcenodeid, l_sourcenodeid.Length);
                Array.Copy(l_targetnodeid, sl_targetnodeid, l_targetnodeid.Length);
                Array.Copy(l_branch_order, sl_branch_order, l_branch_order.Length);
                Array.Copy(l_geopointsX, sl_geopointsX, l_geopointsX.Length);
                Array.Copy(l_geopointsY, sl_geopointsY, l_geopointsY.Length);

                //10. Add the stranger..
                tmpstring = "i_am_the_stranger_nodeid";
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                sl_nodesinfo[l_nnodes].ids  = tmpstring.ToCharArray();
                tmpstring = "i_am_the_stranger_nodelongname";
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                sl_nodesinfo[l_nnodes].longnames = tmpstring.ToCharArray();
                sl_nodesX[l_nnodes] = -1.0;
                sl_nodesY[l_nnodes] = -1.0;


                tmpstring = "i_am_the_stranger_branchid";
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                sl_branchinfo[l_nbranches].ids = tmpstring.ToCharArray();
                tmpstring = "i_am_the_stranger_branchlongname";
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                sl_branchinfo[l_nbranches].longnames = tmpstring.ToCharArray();
                sl_branchlengths[l_nbranches] =-1.0;
                sl_nbranchgeometrypoints[l_nbranches] = 1;
                sl_sourcenodeid[l_nbranches] =  -1;
                sl_targetnodeid[l_nbranches] =  -1;
                sl_branch_order[l_nbranches] =  -1;

                sl_geopointsX[l_nGeometry] = -1.0;
                sl_geopointsY[l_nGeometry] = -1.0;

                //11. Creates the second file, will not add any dataset 
                c_path = TestHelper.TestDirectoryPath() + @"\LargeSewerSystemSecondTest.nc";
                TestHelper.DeleteIfExists(c_path);
                Assert.IsFalse(File.Exists(c_path));
                mode = 1;       // create in write mode
                ierr = wrapper.ionc_create(c_path, ref mode, ref ioncid);
                Assert.That(ierr, Is.EqualTo(0));
                Assert.IsTrue(File.Exists(c_path));
                addglobalattributes(ioncid, ref wrapper);

                //12. Creates and write 1d network
                string sl_network_name = "sewer_system_with_the_stranger";
                ierr = wrapper.ionc_create_1d_network(ref ioncid, ref networkid, sl_network_name, ref sl_nnodes,
                    ref sl_nbranches, ref sl_nGeometry);
                Assert.That(ierr, Is.EqualTo(0));

                write1dnetwork(ioncid,
                    networkid,
                    ref wrapper,
                    sl_nnodes,
                    sl_nbranches,
                    sl_nGeometry,
                    ref sl_nodesinfo,
                    ref sl_branchinfo,
                    ref sl_nodesX,
                    ref sl_nodesY,
                    ref sl_sourcenodeid,
                    ref sl_targetnodeid,
                    ref sl_branchlengths,
                    ref sl_nbranchgeometrypoints,
                    ref sl_geopointsX,
                    ref sl_geopointsY,
                    ref sl_branch_order);

                //13. Close the second file
                ierr = wrapper.ionc_close(ref ioncid);
                Assert.That(ierr, Is.EqualTo(0));
            },
        stackSize);
            th.Start();
            th.Join();

        }

        // Load 2D, create 1D, create links 1D-2D, save them to file for later processing
        [Test]
        [NUnit.Framework.Category("CreateNetInputFile")]
        public void CreateNetInputFile()
        {

            //1. Load 2d file 
            var wrapperNetcdf = new IoNetcdfLibWrapper();
            string sourcetwod_path = TestHelper.CreateLocalCopy("2d_net_river.nc");
            Assert.IsTrue(File.Exists(sourcetwod_path));
            int sourcetwodioncid = -1; //file id 
            int sourcetwomode = 0; //read mode
            int ierr = wrapperNetcdf.ionc_open(sourcetwod_path, ref sourcetwomode, ref sourcetwodioncid, ref iconvtype,
                ref convversion);
            Assert.That(ierr, Is.EqualTo(0));

            //2. Define the network
            int l_nnodes = 2;
            int l_nbranches = 1;
            int l_nGeometry = 25;
            double[] l_nodesX = { 293.78, 538.89 };
            double[] l_nodesY = { 27.48, 956.75 };
            int[] l_sourcenodeid = { 1 };
            int[] l_targetnodeid = { 2 };
            double[] l_branchlengths = { 1165.29 };
            int[] l_nbranchgeometrypoints = { 25 };
            double[] l_geopointsX =
            {
                293.78,
                278.97,
                265.31,
                254.17,
                247.44,
                248.30,
                259.58,
                282.24,
                314.61,
                354.44,
                398.94,
                445.00,
                490.60,
                532.84,
                566.64,
                589.08,
                600.72,
                603.53,
                599.27,
                590.05,
                577.56,
                562.97,
                547.12,
                530.67,
                538.89
            };
            double[] l_geopointsY =
            {
                27.48,
                74.87,
                122.59,
                170.96,
                220.12,
                269.67,
                317.89,
                361.93,
                399.39,
                428.84,
                450.76,
                469.28,
                488.89,
                514.78,
                550.83,
                594.93,
                643.09,
                692.60,
                742.02,
                790.79,
                838.83,
                886.28,
                933.33,
                980.17,
                956.75
            };
            int[] l_branch_order = { -1 };
            int l_startIndex = startIndex;

            IoNetcdfLibWrapper.interop_charinfo[] l_nodesinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nnodes];
            IoNetcdfLibWrapper.interop_charinfo[] l_branchinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nbranches];

            for (int i = 0; i < l_nnodes; i++)
            {
                string tmpstring = "";
                tmpstring = "nodesids";
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                l_nodesinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = "nodeslongNames";
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                l_nodesinfo[i].longnames = tmpstring.ToCharArray();
            }

            for (int i = 0; i < l_nbranches; i++)
            {
                string tmpstring = "";
                tmpstring = "branchids";
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                l_branchinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = "branchlongNames";
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                l_branchinfo[i].longnames = tmpstring.ToCharArray();
            }

            //3. Now we create a new empty file where to save 1d and 2d meshes
            int targetioncid = -1; //file id  
            int targetmode   =  1;    //create in write mode
            string target_path = TestHelper.TestDirectoryPath() + "/river1_full_net.nc";
            TestHelper.DeleteIfExists(target_path);
            Assert.IsFalse(File.Exists(target_path));

            ierr = wrapperNetcdf.ionc_create(target_path, ref targetmode, ref targetioncid);
            Assert.That(ierr, Is.EqualTo(0));
            Assert.IsTrue(File.Exists(target_path));

            addglobalattributes(targetioncid, ref wrapperNetcdf);

            //4. Get the 2d mesh id from the existing file 
            int meshType = 2;
            int sourcemesh2d = -1;
            getMeshid(sourcetwodioncid, ref sourcemesh2d, meshType, ref wrapperNetcdf);
            Assert.That(sourcemesh2d, Is.EqualTo(1));

            //5. Create 1d geometry and mesh in the new file (target.nc)
            int networkid = -1;
            ierr = wrapperNetcdf.ionc_create_1d_network(ref targetioncid, ref networkid, networkName, ref l_nnodes, ref l_nbranches, ref l_nGeometry);
            Assert.That(ierr, Is.EqualTo(0));

            //6. define the network variables
            StringBuilder l_networkName = new StringBuilder(networkName);
            StringBuilder l_meshname = new StringBuilder(meshname);
            int l_nmeshpoints = 25;
            int l_nedgenodes = 24;
            int l_nlinks = 10;
            double[] l_branchoffset =
            {
                0.00,
                49.65,
                99.29,
                148.92,
                198.54,
                248.09,
                297.62,
                347.15,
                396.66,
                446.19,
                495.80,
                545.44,
                595.08,
                644.63,
                694.04,
                743.52,
                793.07,
                842.65,
                892.26,
                941.89,
                991.53,
                1041.17,
                1090.82,
                1140.46,
                1165.29
            };

            int[] l_branchidx =
            {
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1
            };

            IoNetcdfLibWrapper.interop_charinfo[] meshnodeidsinfo = new IoNetcdfLibWrapper.interop_charinfo[l_nmeshpoints];
            for (int i = 0; i < l_nmeshpoints; i++)
            {
                string tmpstring = "";
                tmpstring = "meshnodeids";
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.idssize, ' ');
                meshnodeidsinfo[i].ids = tmpstring.ToCharArray();
                tmpstring = "meshnodelongnames";
                tmpstring = tmpstring.PadRight(IoNetcdfLibWrapper.longnamessize, ' ');
                meshnodeidsinfo[i].longnames = tmpstring.ToCharArray();
            }


            //7. Write the 1d data in the new file (1d geometry, mesh)
            write1dnetwork(targetioncid,
                networkid,
                ref wrapperNetcdf,
                l_nnodes,
                l_nbranches,
                l_nGeometry,
                ref l_nodesinfo,
                ref l_branchinfo,
                ref l_nodesX,
                ref l_nodesY,
                ref l_sourcenodeid,
                ref l_targetnodeid,
                ref l_branchlengths,
                ref l_nbranchgeometrypoints,
                ref l_geopointsX,
                ref l_geopointsY,
                ref l_branch_order);

            write1dmesh(
                targetioncid,
                ref l_networkName,
                ref l_meshname,
                ref wrapperNetcdf,
                ref meshnodeidsinfo,
                l_nmeshpoints,
                l_nedgenodes,
                l_nbranches,
                l_nlinks,
                ref l_branchoffset,
                ref l_branchlengths,
                ref l_branchidx, 
                ref l_sourcenodeid,
                ref l_targetnodeid,
                l_startIndex);

            //8. Clone the 2d mesh definitions in the new file
            int target2dmesh = -1;
            ierr = wrapperNetcdf.ionc_clone_mesh_definition(ref sourcetwodioncid, ref targetioncid, ref sourcemesh2d,
                ref target2dmesh);
            Assert.That(ierr, Is.EqualTo(0));

            //9. Clone the 2d mesh data
            ierr = wrapperNetcdf.ionc_clone_mesh_data(ref sourcetwodioncid, ref targetioncid, ref sourcemesh2d,
                ref target2dmesh);
            Assert.That(ierr, Is.EqualTo(0));

            //10. Close all target and source files
            ierr = wrapperNetcdf.ionc_close(ref sourcetwodioncid);
            Assert.That(ierr, Is.EqualTo(0));
            ierr = wrapperNetcdf.ionc_close(ref targetioncid);
            Assert.That(ierr, Is.EqualTo(0));

            //----------------------------------------------------------------------//
            // Create the links from the file containing 2d and 1d ugrid meshes
            //----------------------------------------------------------------------//
            
            //11. Open file again
            int mode = 1; // weite mode (need to write the links)
            ierr = wrapperNetcdf.ionc_open(target_path, ref mode, ref targetioncid, ref iconvtype, ref convversion);

            //12. Get 2d mesh id
            meshType = 2;
            int mesh2d = -1;
            getMeshid(targetioncid, ref mesh2d, meshType, ref wrapperNetcdf);

            //13. Get 2d mesh dimensions
            var meshtwoddim = new meshgeomdim();
            ierr = wrapperNetcdf.ionc_get_meshgeom_dim(ref targetioncid, ref mesh2d, ref meshtwoddim);
            Assert.That(ierr, Is.EqualTo(0));

            //14. Get the 2d mesh arrays
            var meshtwod = new meshgeom();
            meshtwod.nodex = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * meshtwoddim.numnode);
            meshtwod.nodey = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * meshtwoddim.numnode);
            meshtwod.nodez = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * meshtwoddim.numnode);
            meshtwod.edge_nodes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * meshtwoddim.numedge * 2);
            bool includeArrays = true;
            int start_index = 1;
            ierr = wrapperNetcdf.ionc_get_meshgeom(ref targetioncid, ref mesh2d, ref meshtwod, ref start_index, ref includeArrays);
            Assert.That(ierr, Is.EqualTo(0));

            //15. Using the existing arrays in memory (delta shell scenario), convert 1d into herman datastructure
            IntPtr c_meshXCoords = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nmeshpoints);
            IntPtr c_meshYCoords = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nmeshpoints);
            IntPtr c_branchoffset = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nmeshpoints);
            IntPtr c_branchlength = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * l_nbranches);
            IntPtr c_branchids = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nmeshpoints);
            IntPtr c_sourcenodeid = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nbranches);
            IntPtr c_targetnodeid = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * l_nbranches);


            Marshal.Copy(l_geopointsX, 0, c_meshXCoords, l_nmeshpoints); // in this case geopoints == meshpoints
            Marshal.Copy(l_geopointsY, 0, c_meshYCoords, l_nmeshpoints);
            Marshal.Copy(l_branchoffset, 0, c_branchoffset, l_nmeshpoints);
            Marshal.Copy(l_branchlengths, 0, c_branchlength, l_nbranches);
            Marshal.Copy(l_branchidx, 0, c_branchids, l_nmeshpoints);
            Marshal.Copy(sourcenodeid, 0, c_sourcenodeid, l_nbranches);
            Marshal.Copy(targetnodeid, 0, c_targetnodeid, l_nbranches);

            //16. fill kn (Herman datastructure) for creating the links
            var wrapperGridgeom = new GridGeomLibWrapper();
            ierr = wrapperGridgeom.ggeo_convert_1d_arrays
                (
                ref c_meshXCoords,
                ref c_meshYCoords,
                ref c_branchoffset,
                ref c_branchlength,
                ref c_branchids,
                ref c_sourcenodeid,
                ref c_targetnodeid,
                ref l_nbranches,
                ref l_nmeshpoints,
                ref l_startIndex
                );
            Assert.That(ierr, Is.EqualTo(0));
            ierr = wrapperGridgeom.ggeo_convert(ref meshtwod, ref meshtwoddim);
            Assert.That(ierr, Is.EqualTo(0));

            //17. make the links
            ierr = wrapperGridgeom.ggeo_make1D2Dinternalnetlinks();
            Assert.That(ierr, Is.EqualTo(0));

            //18. get the number of links
            int n2dl1dinks = 0;
            ierr = wrapperGridgeom.ggeo_get_links_count(ref n2dl1dinks);
            Assert.That(ierr, Is.EqualTo(0));

            //19. get the links: arrayfrom = 2d cell index, arrayto = 1d node index 
            IntPtr c_arrayfrom = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * n2dl1dinks); //2d cell number
            IntPtr c_arrayto = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * n2dl1dinks); //1d node
            ierr = wrapperGridgeom.ggeo_get_links(ref c_arrayfrom, ref c_arrayto, ref n2dl1dinks);
            Assert.That(ierr, Is.EqualTo(0));

            int[] rc_arrayfrom = new int[n2dl1dinks];
            int[] rc_arrayto = new int[n2dl1dinks];
            Marshal.Copy(c_arrayfrom, rc_arrayfrom, 0, n2dl1dinks);
            Marshal.Copy(c_arrayto, rc_arrayto, 0, n2dl1dinks);

            int l_linkmesh1   =  mesh2d;
            int l_linkmesh2 =  -1;
            getMeshid(targetioncid, ref l_linkmesh2, 1, ref wrapperNetcdf);
            int l_locationType1 = 1;
            int l_locationType2 = 1;
            StringBuilder l_linkmeshname = new StringBuilder("2d1dlinks");
            int[] contacttype = new int[n2dl1dinks];
            for (int i = 0; i < n2dl1dinks; i++)
            {
                contacttype[i] = 3;
            }

            //20. write 2d-1d links to file.
            write1d2dlinks(
                targetioncid,
                l_linkmesh1,
                l_linkmesh2,
                l_locationType1,
                l_locationType2,
                ref l_linkmeshname,
                ref wrapperNetcdf,
                n2dl1dinks,
                ref rc_arrayfrom,
                ref rc_arrayto,
                ref contacttype,
                l_startIndex);

            //21. Close all files 
            ierr = wrapperNetcdf.ionc_close(ref targetioncid);
            Assert.That(ierr, Is.EqualTo(0));

            //22. Free memory  
               //Free 2d arrays
            Marshal.FreeCoTaskMem(meshtwod.nodex);
            Marshal.FreeCoTaskMem(meshtwod.nodey);
            Marshal.FreeCoTaskMem(meshtwod.nodez);
            Marshal.FreeCoTaskMem(meshtwod.edge_nodes);
              //Free 1d arrays
            Marshal.FreeCoTaskMem(c_meshXCoords);
            Marshal.FreeCoTaskMem(c_meshYCoords);
            Marshal.FreeCoTaskMem(c_branchoffset);
            Marshal.FreeCoTaskMem(c_branchlength);
            Marshal.FreeCoTaskMem(c_branchids);
            Marshal.FreeCoTaskMem(c_sourcenodeid);
            Marshal.FreeCoTaskMem(c_targetnodeid);
              //Free links arrays
            Marshal.FreeCoTaskMem(c_arrayfrom);
            Marshal.FreeCoTaskMem(c_arrayto);
        }

        // Test: get only 1d network using get mesh geom
        [Test]
        [NUnit.Framework.Category("Read1dNetworkUsingGetMeshGeom")]
        public void Load1dNetworkUsingGetMeshGeom()
        {
            //1.Open a netcdf file
            string c_path = TestHelper.TestDirectoryPath() + @"\write1dNetwork.nc";
            Assert.IsTrue(File.Exists(c_path));
            int ioncid = 0; //file variable 
            int mode = 0; //create in read mode
            var wrapper = new IoNetcdfLibWrapper();
            var ierr = wrapper.ionc_open(c_path, ref mode, ref ioncid, ref iconvtype, ref convversion);
            Assert.That(ierr, Is.EqualTo(0));

            //2. Get 1d mesh dimensions
            var networkdim = new meshgeomdim();
            int l_meshid = 0; //set an invalid index
            int nnumNetworks = -1;
            ierr = wrapper.ionc_get_number_of_networks(ref ioncid, ref nnumNetworks);
            Assert.That(ierr, Is.EqualTo(0));
            Assert.That(nnumNetworks, Is.EqualTo(1));
            IntPtr c_networkids = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * nnumNetworks);
            
            //3. Get a valid networkid
            ierr = wrapper.ionc_get_network_ids(ref ioncid, ref c_networkids, ref nnumNetworks);
            Assert.That(ierr, Is.EqualTo(0));
            int[] l_networkid = new int[nnumNetworks];
            Marshal.Copy(c_networkids, l_networkid, 0, nnumNetworks);

            //4. Get network dimensions using meshgeom
            ierr = wrapper.ionc_get_meshgeom_dim(ref ioncid, ref l_meshid, ref l_networkid[0], ref networkdim);
            Assert.That(ierr, Is.EqualTo(0));

            //5. Allocate memory
            var network = new meshgeom();
            network.nnodex = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * networkdim.nnodes);
            network.nnodey = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * networkdim.nnodes);
            network.nedge_nodes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * networkdim.nbranches);
            network.nbranchlengths = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * networkdim.nbranches);
            network.nbranchgeometrynodes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * networkdim.nbranches);
            network.ngeopointx = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * networkdim.ngeometry);
            network.ngeopointy = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(double)) * networkdim.ngeometry);
            network.nbranchorder = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * networkdim.nbranches);

            var includeArrays = true;
            int start_index = 1;
            ierr = wrapper.ionc_get_meshgeom(ref ioncid, ref l_meshid, ref l_networkid[0], ref network, ref start_index, ref includeArrays);
            Assert.That(ierr, Is.EqualTo(0));

            //2. Define the network
            int l_nnodes = 2;
            int l_nbranches = 1;
            int l_nGeometry = 25;
            double[] l_nodesX = new double[networkdim.nnodes];
            double[] l_nodesY = new double[networkdim.nnodes];
            int[] l_nedge_nodes = new int[networkdim.nbranches];
            double[] l_nbranchlengths = new double[networkdim.nbranches];
            double[] l_ngeopointx = new double[networkdim.ngeometry];
            double[] l_ngeopointy = new double[networkdim.ngeometry];

            Marshal.Copy(network.nnodex, l_nodesX, 0, networkdim.nnodes);
            Marshal.Copy(network.nnodey, l_nodesY, 0, networkdim.nnodes);
            Marshal.Copy(network.nedge_nodes, l_nedge_nodes, 0, networkdim.nbranches);
            Marshal.Copy(network.nbranchlengths, l_nbranchlengths, 0, networkdim.nbranches);
            Marshal.Copy(network.ngeopointx, l_ngeopointx, 0, networkdim.ngeometry);
            Marshal.Copy(network.ngeopointy, l_ngeopointy, 0, networkdim.ngeometry);

        }

        //open a file test
        [Test]
        [NUnit.Framework.Category("UGRIDTests")]
        public void openAfile()
        {
            var wrapper = new IoNetcdfLibWrapper();

            // Open a file 
            string c_path = @"D:\carniato\LUCA\ENGINES\FM\FMSource\vanRob\Custom_Ugrid.nc";
            Assert.IsTrue(File.Exists(c_path));
            int ioncid = -1; //file id 
            int sourcetwomode = 0; //read mode
            int ierr = wrapper.ionc_open(c_path, ref sourcetwomode, ref ioncid, ref iconvtype,
                ref convversion);
            Assert.That(ierr, Is.EqualTo(0));
            
            // test functions on costum file 
            int meshid = 0;
            ierr = wrapper.ionc_get_2d_mesh_id(ref ioncid, ref meshid);
            Assert.That(ierr, Is.EqualTo(0));
            int nmaxfacenodes = 0;

            int nnodes = -1;
            ierr = wrapper.ionc_get_node_count(ref ioncid, ref meshid, ref nnodes);
            Assert.That(ierr, Is.EqualTo(0));

            int nedge = -1;
            ierr = wrapper.ionc_get_edge_count(ref ioncid, ref meshid, ref nedge);
            Assert.That(ierr, Is.EqualTo(0));

            int nface = -1;
            ierr = wrapper.ionc_get_face_count(ref ioncid, ref meshid, ref nface);
            Assert.That(ierr, Is.EqualTo(0));

            int maxfacenodes = -1;
            ierr = wrapper.ionc_get_max_face_nodes(ref ioncid, ref meshid, ref maxfacenodes);
            Assert.That(ierr, Is.EqualTo(0));

            IntPtr c_face_nodes = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * nface * maxfacenodes);

            int fillvalue = -1;
            int startIndex = 0;
            ierr = wrapper.ionc_get_face_nodes(ref ioncid, ref meshid, ref c_face_nodes, ref nface, ref maxfacenodes, ref fillvalue, ref startIndex);

            int[] rc_face_nodes = new int[nface * maxfacenodes];
            Marshal.Copy(c_face_nodes, rc_face_nodes, 0, nface * maxfacenodes);
        }
    }
}

