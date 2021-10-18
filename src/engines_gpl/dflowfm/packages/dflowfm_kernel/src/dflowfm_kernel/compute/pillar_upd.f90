!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2021.                                
!                                                                               
!  This file is part of Delft3D (D-Flow Flexible Mesh component).               
!                                                                               
!  Delft3D is free software: you can redistribute it and/or modify              
!  it under the terms of the GNU Affero General Public License as               
!  published by the Free Software Foundation version 3.                         
!                                                                               
!  Delft3D  is distributed in the hope that it will be useful,                  
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU Affero General Public License for more details.                          
!                                                                               
!  You should have received a copy of the GNU Affero General Public License     
!  along with Delft3D.  If not, see <http://www.gnu.org/licenses/>.             
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D",                  
!  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting 
!  Deltares, and remain the property of Stichting Deltares. All rights reserved.
!                                                                               
!-------------------------------------------------------------------------------

! $Id$
! $HeadURL$

  ! =================================================================================================
  ! =================================================================================================
  subroutine pillar_upd()
    use m_flowexternalforcings, only: Cpil
    use m_flowgeom            , only: lnx, ln, dx
    use m_flow                , only: u1, v, advi
    use m_flowparameters      , only: japillar
    implicit none
    integer          :: L, k1, k2
    double precision :: CpilL, uv

    if (japillar == 1) then
      do L = 1,lnx
        k1 = ln(1,L)
        k2 = ln(2,L)
        CpilL = ( Cpil(k1) + Cpil(k2) ) * 0.5d0
        uv = sqrt( u1(L) * u1(L) + v(L) * v(L) )
        advi(L) = advi(L) + CpilL * uv / dx(L)
      enddo
    else if (japillar == 3) then
      do L = 1,lnx
        if (Cpil(L) == 0d0) cycle
        CpilL = Cpil(L)
        uv = sqrt( u1(L) * u1(L) + v(L) * v(L) )
        advi(L) = advi(L) + CpilL * uv / dx(L)
      enddo
    endif

  end subroutine pillar_upd
  
  ! =================================================================================================
  ! =================================================================================================
  subroutine sealock_upd()
    use m_flowgeom            , only: lnx, ln, dx
    use m_flow                , only: sa1, hs, vol1, s1
    use m_flowparameters      , only: jasealock
    use m_flowexternalforcings, only: sealock, nsealocksg, qstss
    use gridoperations
    use m_transport           , only: numconst
    use m_zsf
    implicit none
    integer                 :: L, k, k1, k2, kb, kt, m, isorsin
    double precision        :: salsum, volsum
    type(zsf_param_t)       :: p
    type(zsf_results_t)     :: results
    type(zsf_aux_results_t) :: aux_results
    integer                 :: err_code
    double precision        :: discharge_from_lake
    double precision        :: discharge_to_lake
    double precision        :: discharge_from_sea
    double precision        :: discharge_to_sea
 
    do m = 1,nsealocksg
       call zsf_param_default(p)
       
       k1 = sealock(m)%ksea_probe
       k2 = sealock(m)%klake_probe
       
       call getkbotktop(k1, kb, kt)
       volsum = 0d0
       salsum = 0d0
       do k = kb, kt
          salsum = salsum + sa1(k) * vol1(k)
          volsum = volsum + vol1(k) 
       enddo
       p%salinity_sea = p%salinity_sea / volsum
       p%head_sea = s1(k1)
    
       call getkbotktop(k2, kb, kt)
       volsum = 0d0
       salsum = 0d0
       do k = kb, kt
          salsum = salsum + sa1(k) * vol1(k)
          volsum = volsum + vol1(k) 
       enddo
       p%salinity_lake = salsum / volsum
       p%head_lake = s1(k2)
    
       err_code = zsf_calc_steady (p, results, aux_results)
       
       discharge_to_lake   = results%discharge_to_lake    !TODO: Needs to be checked is it equalent to DFM? Check t_cycle
       discharge_from_lake = results%discharge_from_lake
       discharge_to_sea    = results%discharge_to_sea
       discharge_from_sea  = results%discharge_from_sea
    !   qstss(??) = volume_to_lake / dts
       isorsin = sealock(m)%sorsin_index(1)
       qstss((1+numconst)*(isorsin-1) + 1) = results%discharge_from_sea

       isorsin = sealock(m)%sorsin_index(2)
       qstss((1+numconst)*(isorsin-1) + 1) = results%discharge_from_lake

    enddo
    

  end subroutine sealock_upd
