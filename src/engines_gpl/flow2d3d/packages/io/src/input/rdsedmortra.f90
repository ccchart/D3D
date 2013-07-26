subroutine rdsedmortra(lundia    ,error     ,lsal      ,ltem      ,lsed      , &
                     & lsedtot   ,lstsci    ,ltur      ,facdss    ,namcon    , &
                     & iopsus    ,filnam    ,mmax      ,nmax      ,nmaxus    , &
                     & nmmax     ,nto       ,nambnd    ,lsec      ,gdp       )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2013.                                     
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  $Id$
!  $HeadURL$
!!--description-----------------------------------------------------------------
!
! Read sediment, morphology and transport parameters from the input files
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use properties
    use m_rdmor
    use m_rdsed
    use m_rdtrafrm
    !
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    integer                   , parameter    :: NPARDEF = 20
!
! Global variables
!
    integer                                  , intent(in)  :: lsal    !  Description and declaration in dimens.igs
    integer                                  , intent(in)  :: lsed    !  Description and declaration in iidim.f90
    integer                                  , intent(in)  :: lsedtot !  Description and declaration in iidim.f90
    integer                                  , intent(in)  :: lstsci  !  Description and declaration in iidim.f90
    integer                                  , intent(in)  :: ltem    !  Description and declaration in dimens.igs
    integer                                  , intent(in)  :: ltur    !  Description and declaration in iidim.f90
    integer                                                :: lundia  !  Description and declaration in inout.igs
    logical                                  , intent(out) :: error   !!  Flag=TRUE if an error is encountered
    real(fp)      , dimension(lsed)                        :: facdss  !  Description and declaration in rjdim.f90
    character(20) , dimension(lstsci + ltur)               :: namcon  !  Description and declaration in ckdim.f90
    integer                                                :: iopsus
    character(*)                                           :: filnam
    integer                                  , intent(in)  :: mmax
    integer                                  , intent(in)  :: nmax
    integer                                  , intent(in)  :: nmaxus
    integer                                  , intent(in)  :: nmmax
    integer                                  , intent(in)  :: nto
    integer                                  , intent(in)  :: lsec
    character(20) , dimension(nto)                         :: nambnd  !  Description and declaration in ckdim.f90
!
! Local variables
!
    integer                                  :: i
    real(fp)                                 :: fwfacmor
    character(256)                           :: filmor
    character(256)                           :: filsed
    character(256)                           :: filtrn
    character(256)                           :: string
    type(tree_data)               , pointer  :: mor_ptr
    type(tree_data)               , pointer  :: sed_ptr
    integer, dimension(2,NPARDEF)            :: ipardef
    real(fp), dimension(NPARDEF)             :: rpardef
!
!! executable statements -------------------------------------------------------
!
    !
    ! Read name of default transport formula
    !
    filtrn = ' '
    call prop_get_string(gdp%mdfile_ptr,'*','TraFrm',filtrn)
    !
    call initrafrm(lundia    ,error     ,lsedtot   ,gdp%gdtrapar)
    if (error) goto 999
    !
    ! Read name of sediment input file
    !
    filsed = ' '
    call prop_get_string(gdp%mdfile_ptr, '*', 'Filsed', filsed)
    !
    ! Sediment input has been placed in input_tree in subroutine dimsedconst
    ! get pointer
    !
    call tree_get_node_by_name( gdp%input_tree, 'Sediment Input', sed_ptr )
    !
    ! Read data from that file
    !
    call rdsed(lundia    ,error     ,lsal      ,ltem      ,lsed      , &
             & lsedtot   ,lstsci    ,ltur      ,facdss    ,namcon    , &
             & iopsus    ,gdp%d%nmlb,gdp%d%nmub,filsed    , &
             & sed_ptr   ,gdp%gdsedpar,gdp%gdtrapar)
    if (error) goto 999
    !
    ! Read name of morphology input file
    !
    filmor = ' '
    call prop_get_string(gdp%mdfile_ptr, '*', 'Filmor', filmor)
    !
    ! Create Morphology branch in input tree
    !
    call tree_create_node(gdp%input_tree, 'Morphology Input', mor_ptr )
    call tree_put_data(mor_ptr, transfer(trim(filmor),node_value), 'STRING' )
    !
    ! Read data from that file
    !
    call rdmor(lundia     ,error     ,filmor    ,lsec      ,lsedtot    , &
             & nmaxus     ,nto       , &
             & nambnd     ,gdp%gdinttim%julday  ,mor_ptr   ,gdp%gdsedpar, &
             &gdp%gdmorpar,fwfacmor  ,gdp%gdmorlyr, gdp%griddim)
    if (error) goto 999
    !
    ! FWFac is independent of sediment transport. Use the value historically
    ! specified in mor file only if it hasn't been specified in the mdf file.
    !
    string = ' '
    call prop_get(gdp%mdfile_ptr, '*', 'FWFac' , string)
    if (string == ' ') then
       gdp%gdnumeco%fwfac = fwfacmor
    endif
    !
    ! Some other parameters are transport formula specific. Use the value
    ! historically specified in mor file as default.
    !
    ipardef = 0
    rpardef = 0.0_fp
    call setpardef(ipardef, rpardef, NPARDEF, -1, 1, gdp%gdmorpar%iopsus)
    call setpardef(ipardef, rpardef, NPARDEF, -1, 2, gdp%gdmorpar%aksfac)
    call setpardef(ipardef, rpardef, NPARDEF, -1, 4, gdp%gdmorpar%rwave)
    call setpardef(ipardef, rpardef, NPARDEF, -1, 4, gdp%gdmorpar%rdc)
    call setpardef(ipardef, rpardef, NPARDEF, -1, 5, gdp%gdmorpar%rdw)
    call setpardef(ipardef, rpardef, NPARDEF, -1, 6, gdp%gdmorpar%iopkcw)
    call setpardef(ipardef, rpardef, NPARDEF, -1, 7, gdp%gdmorpar%epspar)
    call setpardef(ipardef, rpardef, NPARDEF, -2, 1, gdp%gdmorpar%iopsus)
    call setpardef(ipardef, rpardef, NPARDEF, -2, 2, gdp%gdmorpar%pangle)
    call setpardef(ipardef, rpardef, NPARDEF, -2, 3, gdp%gdmorpar%fpco)
    call setpardef(ipardef, rpardef, NPARDEF, -2, 4, gdp%gdmorpar%subiw)
    call setpardef(ipardef, rpardef, NPARDEF, -2, 5, gdp%gdmorpar%epspar)
    !
    call rdtrafrm(lundia    ,error     ,filtrn    ,lsedtot   , &
                & ipardef   ,rpardef   ,NPARDEF   ,gdp%gdtrapar, &
                & gdp%gdsedpar%sedtyp  ,gdp%gdsedpar%sedblock  , &
                & gdp%griddim)
    if (error) goto 999
    !
    !--------------------------------------------------------------------------
    !
    ! Echo sediment and transport parameters
    !
    call echosed(lundia    ,error     ,lsed      ,lsedtot   ,facdss    , &
               & iopsus    ,gdp%gdsedpar, gdp%gdtrapar)
    if (error) goto 999
    !
    ! Echo morphology parameters
    !
    call echomor(lundia    ,error     ,lsec      ,lsedtot   ,nto        , &
               & nambnd    ,gdp%gdsedpar, gdp%gdmorpar)
    if (error) goto 999
    !
    ! Read scour and echo parameters
    !
    call rdscour(lundia    ,error     ,nmmax     ,gdp       )
    !
    ! If Van Rijn 2004 transport formula is used (iform = -2), switch on the
    ! bed roughness height predictor. By default this predictor is set to the
    ! Van Rijn 2004 formulations; give a warning if this has been set to a
    ! different predictor by the user.
    !
    do i = 1, lsedtot
       if (gdp%gdtrapar%iform(i) == -2) then
          if (gdp%gdbedformpar%bdfrpt /= 0) then
             call prterr(lundia, 'U190', 'Van Rijn 2004 transport formula combined with different bedform roughness predictor')
          endif
          gdp%gdbedformpar%lfbedfrmrou = .true.
          exit
       endif
    enddo
999 continue
    if (error) call d3stop(1, gdp)
end subroutine rdsedmortra
