subroutine initdredge(gdp       )
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011.                                     
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
!!--description-----------------------------------------------------------------
! NONE
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use precision
    use globaldata
    !
    implicit none
    !
    type(globdat),target :: gdp
    !
    ! The following list of pointer parameters is used to point inside the gdp structure
    !
    real(fp)      , dimension(:,:) , pointer :: link_percentage
    real(fp)      , dimension(:)   , pointer :: link_distance
    real(fp)      , dimension(:,:) , pointer :: link_sum
    real(fp)      , dimension(:)   , pointer :: dzdred
    real(fp)      , dimension(:)   , pointer :: refplane
    real(fp)      , dimension(:,:) , pointer :: voldred
    real(fp)      , dimension(:)   , pointer :: voldune
    real(fp)      , dimension(:)   , pointer :: totvoldred
    real(fp)      , dimension(:)   , pointer :: globalareadred
    real(fp)      , dimension(:,:) , pointer :: voldump
    real(fp)      , dimension(:,:) , pointer :: percsupl
    real(fp)      , dimension(:)   , pointer :: totvoldump
    real(fp)      , dimension(:)   , pointer :: localareadump
    real(fp)      , dimension(:)   , pointer :: globalareadump
    real(fp)      , dimension(:)   , pointer :: globaldumpcap
    real(fp)                       , pointer :: alpha_dh
    integer                        , pointer :: dredge_domainnr
    integer                        , pointer :: dredge_ndomains
    integer                        , pointer :: nadred
    integer                        , pointer :: nadump
    integer                        , pointer :: nasupl
    integer                        , pointer :: nalink
    integer       , dimension(:,:) , pointer :: link_def
    logical                        , pointer :: tsmortime
    logical                        , pointer :: firstdredge
    logical                        , pointer :: use_dunes
    character(256)                 , pointer :: dredgefile
    character( 80), dimension(:)   , pointer :: dredge_areas
    character( 80), dimension(:)   , pointer :: dump_areas
    type (dredtype), dimension(:)  , pointer :: dredge_prop
    type (dumptype), dimension(:)  , pointer :: dump_prop
    !
    integer :: istat
!
!! executable statements -------------------------------------------------------
!
    link_percentage   => gdp%gddredge%link_percentage
    link_distance     => gdp%gddredge%link_distance
    link_sum          => gdp%gddredge%link_sum
    dzdred            => gdp%gddredge%dzdred
    refplane          => gdp%gddredge%refplane
    voldred           => gdp%gddredge%voldred
    voldune           => gdp%gddredge%voldune
    totvoldred        => gdp%gddredge%totvoldred
    globalareadred    => gdp%gddredge%globalareadred
    voldump           => gdp%gddredge%voldump
    percsupl          => gdp%gddredge%percsupl
    totvoldump        => gdp%gddredge%totvoldump
    localareadump     => gdp%gddredge%localareadump
    globalareadump    => gdp%gddredge%globalareadump
    globaldumpcap     => gdp%gddredge%globaldumpcap
    alpha_dh          => gdp%gddredge%alpha_dh
    dredge_domainnr   => gdp%gddredge%dredge_domainnr
    dredge_ndomains   => gdp%gddredge%dredge_ndomains
    nadred            => gdp%gddredge%nadred
    nadump            => gdp%gddredge%nadump
    nasupl            => gdp%gddredge%nasupl
    nalink            => gdp%gddredge%nalink
    link_def          => gdp%gddredge%link_def
    tsmortime         => gdp%gddredge%tsmortime
    firstdredge       => gdp%gddredge%firstdredge
    use_dunes         => gdp%gddredge%use_dunes
    dredgefile        => gdp%gddredge%dredgefile
    dredge_areas      => gdp%gddredge%dredge_areas
    dump_areas        => gdp%gddredge%dump_areas
    dredge_prop       => gdp%gddredge%dredge_prop
    dump_prop         => gdp%gddredge%dump_prop
    !
    nullify(gdp%gddredge%link_percentage)
    nullify(gdp%gddredge%link_distance)
    nullify(gdp%gddredge%link_sum)
    nullify(gdp%gddredge%dzdred)
    nullify(gdp%gddredge%refplane)
    nullify(gdp%gddredge%voldred)
    nullify(gdp%gddredge%totvoldred)
    nullify(gdp%gddredge%globalareadred)
    nullify(gdp%gddredge%voldune)
    nullify(gdp%gddredge%percsupl)
    nullify(gdp%gddredge%totvoldump)
    nullify(gdp%gddredge%localareadump)
    nullify(gdp%gddredge%globalareadump)
    nullify(gdp%gddredge%globaldumpcap)
    nullify(gdp%gddredge%voldump)
    !
    dredge_domainnr = 0
    dredge_ndomains = 0
    nadred  = 0
    nadump  = 0
    nasupl  = 0
    nalink  = 0
    alpha_dh = 0.5_fp
    !
    nullify(gdp%gddredge%link_def)
    !
    tsmortime = .false.
    firstdredge = .true.
    use_dunes = .false.
    !
    nullify(gdp%gddredge%dredge_areas)
    nullify(gdp%gddredge%dump_areas)
    dredgefile   = ' '
    !
    nullify(gdp%gddredge%dredge_prop)
    nullify(gdp%gddredge%dump_prop)
end subroutine initdredge
