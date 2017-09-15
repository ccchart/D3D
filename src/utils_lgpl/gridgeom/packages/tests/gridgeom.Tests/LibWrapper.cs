﻿using System;
using System.Runtime.InteropServices;
using General.tests;

namespace gridgeom.Tests
{
    //this class wraps the single library functions
    public class GridGeomLibWrapper
    {
        #region ggeo_functions
        /// <summary>
        /// Converts 1d ugrid coordinates to x-y coordinates
        /// </summary>
        /// <param name="c_branchids">The branch ids</param>
        /// <param name="c_branchoffsets">The branch offsets</param>
        /// <param name="c_geopointsX">The x coordiantes of the geometrical points</param>
        /// <param name="c_geopointsY">The y coordiantes of the geometrical points</param>
        /// <param name="c_nbranchgeometrynodes">The number of geometrical points per branch</param>
        /// <param name="c_branchlengths">The branch lengths</param>
        /// <param name="c_meshXCoords">The x coordinates of the mesh points</param>
        /// <param name="c_meshYCoords">The y coordinates of the mesh points</param>
        /// <param name="nbranches">The number of branches</param>
        /// <param name="ngeopoints">The number of geometrical points</param>
        /// <param name="nmeshnodes">The number of mesh nodes</param>
        /// <returns></returns>
        [DllImport(LibDetails.LIB_DLL_NAME, EntryPoint = "ggeo_get_xy_coordinates", CallingConvention = CallingConvention.Cdecl)]
        public static extern int ggeo_get_xy_coordinates_dll(
            [In] ref IntPtr c_branchids,
            [In] ref IntPtr c_branchoffsets,
            [In] ref IntPtr c_geopointsX,
            [In] ref IntPtr c_geopointsY,
            [In] ref IntPtr c_nbranchgeometrynodes,
            [In] ref IntPtr c_branchlengths,
            [In, Out] ref IntPtr c_meshXCoords,
            [In, Out] ref IntPtr c_meshYCoords,
            [In] ref int nbranches,
            [In] ref int ngeopoints,
            [In] ref int nmeshnodes);

        /// <summary>
        /// Use meshgeom to fill kn matrix
        /// </summary>
        /// <param name="meshgeom"> The data structure containing the mesh information</param>
        /// <param name="meshgeomdim">The data structure containing the mesh dimensions</param>
        /// <returns></returns>
        [DllImport(LibDetails.LIB_DLL_NAME, EntryPoint = "ggeo_convert", CallingConvention = CallingConvention.Cdecl)]
        public static extern int ggeo_convert_dll([In, Out] ref meshgeom meshgeom, [In] ref meshgeomdim meshgeomdim);

        /// <summary>
        /// Makes the 1d/2d links (results are stored in memory)
        /// </summary>
        /// <returns></returns>
        [DllImport(LibDetails.LIB_DLL_NAME, EntryPoint = "ggeo_make1D2Dinternalnetlinks", CallingConvention = CallingConvention.Cdecl)]
        public static extern int ggeo_make1D2Dinternalnetlinks_dll();

        /// <summary>
        /// Use 1d array to fill kn matrix
        /// </summary>
        /// <param name="c_meshXCoords">The x coordinates of the mesh points</param>
        /// <param name="c_meshYCoords">The y coordinates of the mesh points</param>
        /// <param name="c_branchids">The branch ids</param>
        /// <param name="nbranches">The number of branches</param>
        /// <param name="nmeshnodes">The number of mesh nodes</param>
        /// <returns></returns>
        [DllImport(LibDetails.LIB_DLL_NAME, EntryPoint = "ggeo_convert_1d_arrays", CallingConvention = CallingConvention.Cdecl)]
        public static extern int ggeo_convert_1d_arrays_dll([In] ref IntPtr c_meshXCoords, [In] ref IntPtr c_meshYCoords, [In] ref IntPtr c_branchids, [In] ref int nbranches, [In] ref int nmeshnodes);


        /// <summary>
        /// Gets the number of 1d-2d links produced by ggeo_make1D2Dinternalnetlinks_dll
        /// </summary>
        /// <param name="nlinks">The number of links</param>
        /// <returns></returns>
        [DllImport(LibDetails.LIB_DLL_NAME, EntryPoint = "ggeo_get_links_count", CallingConvention = CallingConvention.Cdecl)]
        public static extern int ggeo_get_links_count_dll([In, Out] ref int nlinks);

        /// <summary>
        /// Gets the number the 1d-2d links produced by ggeo_make1D2Dinternalnetlinks_dll
        /// </summary>
        /// <param name="arrayfrom">The cell indexes where the links start</param>
        /// <param name="arrayto">The node indexes where the links end</param>
        /// <param name="nlinks">The number of links</param>
        /// <returns></returns>
        [DllImport(LibDetails.LIB_DLL_NAME, EntryPoint = "ggeo_get_links", CallingConvention = CallingConvention.Cdecl)]
        public static extern int ggeo_get_links_dll([In, Out] ref IntPtr arrayfrom, [In, Out] ref IntPtr arrayto, [In] ref int nlinks);

        #endregion ggeo_functions


        public int ggeo_get_xy_coordinates(
            ref IntPtr c_branchids,
            ref IntPtr c_branchoffsets,
            ref IntPtr c_geopointsX,
            ref IntPtr c_geopointsY,
            ref IntPtr c_nbranchgeometrynodes,
            ref IntPtr c_branchlengths,
            ref IntPtr c_meshXCoords,
            ref IntPtr c_meshYCoords,
            ref int nbranches,
            ref int ngeopoints,
            ref int nmeshnodes
        )
        {
          int ierr = ggeo_get_xy_coordinates_dll(
                ref  c_branchids,
                ref  c_branchoffsets,
                ref  c_geopointsX,
                ref  c_geopointsY,
                ref  c_nbranchgeometrynodes,
                ref  c_branchlengths,
                ref  c_meshXCoords,
                ref  c_meshYCoords,
                ref  nbranches,
                ref  ngeopoints,
                ref  nmeshnodes);
            return ierr;
        }


        public int ggeo_convert(ref meshgeom c_meshgeom, ref meshgeomdim c_meshgeomdim )
        {
            int ierr = ggeo_convert_dll(ref  c_meshgeom, ref  c_meshgeomdim );
            return ierr;
        }

        public int ggeo_make1D2Dinternalnetlinks()
        {
            int ierr = ggeo_make1D2Dinternalnetlinks_dll();
            return ierr;
        }

        public int ggeo_convert_1d_arrays(ref IntPtr c_meshXCoords, ref IntPtr c_meshYCoords, ref IntPtr c_branchids, ref int  nbranches,ref int nmeshnodes)
        {
            int ierr = ggeo_convert_1d_arrays_dll(ref c_meshXCoords, ref c_meshYCoords, ref c_branchids, ref nbranches, ref nmeshnodes);
            return ierr;
        }

        public int ggeo_get_links_count(ref int nbranches)
        {
            int ierr = ggeo_get_links_count_dll(ref nbranches);
            return ierr;
        }

        public int ggeo_get_links(ref IntPtr arrayfrom, ref IntPtr arrayto, ref int nlinks)
        {
            int ierr = ggeo_get_links_dll(ref arrayfrom, ref arrayto, ref nlinks);
            return ierr;
        }

    }
}
