!!  Copyright (C)  Stichting Deltares, 2012-2014.
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License version 3,
!!  as published by the Free Software Foundation.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program. If not, see <http://www.gnu.org/licenses/>.
!!
!!  contact: delft3d.support@deltares.nl
!!  Stichting Deltares
!!  P.O. Box 177
!!  2600 MH Delft, The Netherlands
!!
!!  All indications and logos of, and references to registered trademarks
!!  of Stichting Deltares remain the property of Stichting Deltares. All
!!  rights reserved.

      subroutine phcarb  ( pmsa   , fl     , ipoint , increm , noseg  ,
     &                    noflux , iexpnt , iknmrk , noq1   , noq2   ,
     &                    noq3   , noq4   )
!>\file
!>       Integrated calculation of pH and CO2 system
 
C***********************************************************************
C
C     Author  : Lucas Guerrero - Anna de Kluijver
C     Editing : Michel Jeuken
C     Date    : 030714             Version : 0.01 for pH-carb in delwaq
C
C     History :
C
C     Date    Author          Description
C     ------  --------------  -----------------------------------
C     030714  LG - AdK - MJ   first version
C
C***********************************************************************
C
C     Description of the module :
C
C        General water quality module for DELWAQ:
C        Integrated calculation of pH and CO2 system
C        Restrictions: ...
C
C Name    T   L I/O   Description                                    Uni
C ----    --- -  -    -------------------                             --
C AHPLUS  R*4 1 L  activity of H+                                 [mol/kg]
C ALK     R*4 1 L  Alkalinity as HCO3-                            [mmol/kg]
C ALKA    R*4 1 I  Alkalinity as HCO3-                            [g/m3]
C TICM    R*4 1 L  total carbonate concentration                  [mmol/kg]
C TIC     R*4 1 I  total carbonate concentration                  [gC/m3]
C BT      R*4 1 L  total borate                                   [mmol/kg]
C Ca      R*4 1 L  calcium concentration                          [mmol/kg]
C MMTOM   R*4 1 L  conversion from mmol to mol                    [mmol/mol]
C M3TOL   R*4 1 L  conversion from 1/m3 to 1/l                    [m3/l]
C MC      R*4 1 L  from gC to mol C (molar weight)                [g/mol]
C MCO2    R*4 1 L  from gCO2 to mol CO2 (molar weight)            [g/mol]
C MHCO3   R*4 1 L  fron gHCO3 to mol HCO3 (molar weight)          [g/mol]
C MB      R*4 1 L  from gB to mol B (molar weight)                [g/mol]
C K1      R*4 1 L  dissociation constant CO2-HCO3                 [mol/kg]
C K2      R*4 1 L  dissociation constant HCO3-CO3                 [mol/kg]
C KB      R*4 1 L  dissociation constant BOH3-BOH4                [mol/kg]
C KW      R*4 1 L  dissociation constant H2O                      [mol2/kg2]
C Kcal    R*4 1 L  solubility constant calcite                    [mol2/kg2]
C Karg    R*4 1 L  solubility constant aragonite                  [mol2/kg2]
C K0      R*4 1 L  solubility of CO2 in seawater                  [mol/l*atm]
C PH      R*4 1 O  pH                                             [-]
C CO2     R*4 1 O  concentration CO2 in water                     [gCO2/m3]
C pCO2    R*4 1 O  pressure of CO2	                          [�atm]
C HCO3    R*4 1 O  concentration of HCO3                          [gC/m3]
C CO3	  R*4 1 O  concentration of CO3                           [gC/m3]
C Satcal  R*4 1 O  saturation state of calcite                    [-]
C Satarg  R*4 1 O  saturation state of aragonite                  [-]
C BOH4    R*4 1 O  concentration of borate		          [gB/m3]
C SAL     R*4 1 I  salinity                                       [g/kg]
C TEMP    R*4 1 I  temperature                                    [oC]
C TEMPK   R*4 1 L  temperature in Kelvin                          [K]
C RHOH2O  R*4 1 L  water density                                  [kg/l]
C-----------------------------------------------------------------------

C     Logical Units : -

C     Modules called : -
C
      USE MOD_CHEMCONST
      USE MOD_ACBW_PHSOLVERS

      IMPLICIT NONE
      
      REAL     PMSA  ( * ) , FL    (*)
      DOUBLE PRECISION AHPLUSD, P_VAL
      
      INTEGER  ILUMON
      INTEGER  IPOINT( * ) , INCREM(*) , NOSEG , ISEG, NOFLUX,
     +         IEXPNT(4,*) , IKNMRK(*) , NOQ1, NOQ2, NOQ3, NOQ4
      INTEGER  IP1 , IP2 , IP3 , IP4 , IP5 , IP6 , IP7 , IP8, IP9, IP10,
     +         IP11, IP12, IP13, IP14
      
      LOGICAL,SAVE  :: FIRST = .TRUE.
C
C     Local declarations, constants in source
C
      REAL, PARAMETER :: MC     =    12.0
      REAL, PARAMETER :: MCO2   =    44.0
      REAL, PARAMETER :: MHCO3  =    61.0
      REAL, PARAMETER :: MB     =    10.8
      REAL, PARAMETER :: MMTOM  =     1.0E-3
      REAL, PARAMETER :: M3TOL  =     1.0E-3
      REAL, PARAMETER :: KELVIN =   273.15
      
      REAL            :: SAL, TEMP, TIC, ALKA, PH_MIN, PH_MAX
      REAL            :: PH_OLD, TEMPK, LNKW, KW
      REAL            :: LNK0, K0, LNK1, K1, LNK2, K2, LNKB, KB
      REAL            :: LOGCAL, KCAL, LOGARG, KARG
      REAL            :: RHOH2O, TICM, ALK, BT, CA
      REAL            :: AHPLUS, PH
      REAL            :: CO2, pCO2water, HCO3, CO3, BOH4, SATCAL, SATARG
      
      integer, save   :: nr_mes = 0     ! message count negative total carbonate
      integer, save   :: nrmes2 = 0     ! message count negative salinity
      integer, save   :: nrmes3 = 0     ! message count high salinity
      integer, save   :: nrmes4 = 0     ! message count negative alkalinity
      integer, save   :: nrmes5 = 0     ! message count negative H+
C
      IP1   = IPOINT( 1)
      IP2   = IPOINT( 2)
      IP3   = IPOINT( 3)
      IP4   = IPOINT( 4)
      IP5   = IPOINT( 5)
      IP6   = IPOINT( 6)
      IP7   = IPOINT( 7)
      IP8   = IPOINT( 8)
      IP9   = IPOINT( 9)
      IP10  = IPOINT( 10)
      IP11  = IPOINT( 11)
      IP12  = IPOINT( 12)
      IP13  = IPOINT( 13)
      IP14  = IPOINT( 14)
C
C     Loop over de segmenten
C
      DO 9000 ISEG = 1 , NOSEG
C
C     Eerste kenmerk actief of inactief segment
C
!!    CALL DHKMRK(1,IKNMRK(ISEG),IKMRK1)
C
C     Alleen actieve en bodem segmenten behandelen
C
!!    IF (IKMRK1.EQ.1) THEN
      IF (BTEST(IKNMRK(ISEG),0)) THEN
C
C     Map PMSA on local variables
C
      SAL     = PMSA(IP1 )
      TIC     = PMSA(IP2 )
      ALKA    = PMSA(IP3 )
      TEMP    = PMSA(IP4 )
      PH_MIN  = PMSA(IP5 )
      PH_MAX  = PMSA(IP6 )
      
C     Try to get the old value, this is a good initial value for our solvers.
      PH_OLD  = PMSA(IP7 )
C     Because the value might not exist, if it is outside the range start neutral with pH = 7.0      
      IF (PH_OLD.LT.PH_MIN .OR. PH_OLD.GT.PH_MAX) THEN
          PH_OLD = 7.0e0
      ENDIF
C
C     Error messages

      IF ( TIC .LT. 1E-30 ) THEN
         CALL GETMLU(ILUMON)
         IF ( NR_MES .LT. 10 ) THEN
            NR_MES = NR_MES + 1
            WRITE ( ILUMON , * ) 'WARNING :total carbonate <= 0',
     +                           ' segment=',ISEG,' conc=',TIC
         ENDIF
         IF ( NR_MES .EQ. 10 ) THEN
            NR_MES = NR_MES + 1
            WRITE(ILUMON,*) ' 10 WARNINGS on total carbonate'
            WRITE(ILUMON,*) ' No further messages on total carbonate'
         ENDIF
         TIC = 1E-30
      ENDIF
      IF ( SAL .LT. 1E-30 ) THEN
         CALL GETMLU(ILUMON)
         IF ( NRMES2 .LT. 10 ) THEN
            NRMES2 = NRMES2 + 1
            WRITE ( ILUMON , * ) 'WARNING :salinity <= 0',
     +                           ' segment=',ISEG,' conc=',SAL
         ENDIF
         IF ( NRMES2 .EQ. 10 ) THEN
            NRMES2 = NRMES2 + 1
            WRITE(ILUMON,*) ' 10 WARNINGS on salinity'
            WRITE(ILUMON,*) ' No further messages on salinity'
         ENDIF
         SAL = 1E-30
      ENDIF
      IF ( SAL .GT. 50. ) THEN
         CALL GETMLU(ILUMON)
         IF ( NRMES4 .LT. 10 ) THEN
            NRMES4 = NRMES4 + 1
            WRITE ( ILUMON , * ) 'WARNING :salinity => 50.',
     +                           ' segment=',ISEG,' conc=',SAL
         ENDIF
         IF ( NRMES4 .EQ. 10 ) THEN
            NRMES4 = NRMES4 + 1
            WRITE(ILUMON,*) ' 10 WARNINGS on salinity'
            WRITE(ILUMON,*) ' No further messages on salinity'
         ENDIF
         SAL = 50.
      ENDIF
      IF ( ALKA .LT. 1E-30 ) THEN
         CALL GETMLU(ILUMON)
         IF ( NRMES3 .LT. 10 ) THEN
            NRMES3 = NRMES3 + 1
            WRITE ( ILUMON , * ) 'WARNING: alkalinity <= 0',
     +                           ' segment=',ISEG,' conc=',ALKA
         ENDIF
         IF ( NRMES3 .EQ. 10 ) THEN
            NRMES3 = NRMES3 + 1
            WRITE(ILUMON,*) ' 10 WARNINGS on alkalinity'
            WRITE(ILUMON,*) ' No further messages on alkalinity'
         ENDIF
         ALKA = 1E-30
      ENDIF
      IF (TEMP .LE. -KELVIN) THEN
        WRITE (ILUMON,*) ' WARNING: Temperature drops below 0 Kelvin',
     +   ' segment=',ISEG,' Temp set to 15 oC (288.15 K)'
        TEMP = 15.0E0
      ENDIF
C
C---- Process formulation ---------------------------------------
C ********************************
C Dissociation constants depending on temperature and salinity
      TEMPK   = TEMP + KELVIN

C Dissociation constant of water. DOE (1994), Zeebe and Wolf-Gladrow (2001). Total pH scale. [mol^2/kg^2 solution]
      LNKW = 148.96502 - 13847.26/TEMPK - 23.6521 * log(TEMPK) + 
     +       (118.67/TEMPK - 5.977 + 1.0495 * log(TEMPK)) *
     +       SAL**0.5 - 0.01615 * SAL
      
      KW  = exp (LNKW)
     
      IF (SAL .LT. 5.0e0) THEN
C Dissociation constant of carbonic acid. [mol/kg H2O] (Molality)
C Roy et al (1993), Millero (1995). Salinities below 5 (freshwater). Total pH scale.
         LNK1 = 290.9097 - 14554.21/TEMPK - 45.0575*log(TEMPK) +
     +          (-228.39774 + 9714.36839/TEMPK + 34.485796*log(TEMPK)) *
     +          SAL**0.5 + (54.20871 - 2310.48919/TEMPK -
     +          8.19515*log(TEMPK)) * SAL + (-3.969101 + 170.22169/TEMPK +
     +          0.603627*log(TEMPK))* SAL**1.5 - 0.00258768* SAL**2.0

C --- Unit of K1 and K2 [mol/kg H2O]. To convert to [mol/kg solution] a fraction is added (Millero 1995)
         LNK1 = LNK1 + log(1.0E0-SAL*0.001005)

         K1  = exp (LNK1)

C Roy et al (1993), Millero (1995). Salinities below 5 (freshwater). Total pH scale.
         LNK2 = 207.6548 - 11843.79/TEMPK - 33.6485*log(TEMPK) +
     +          (-167.69908 + 6551.35253/TEMPK + 25.928788*log(TEMPK)) *
     +          SAL**0.5 + (39.75854 - 1566.13883/TEMPK -
     +          6.171951*log(TEMPK)) * SAL + (-2.892532 + 116.270079/TEMPK
     +       + 0.45788501*log(TEMPK))*SAL**1.5 - 0.00613142* SAL**2.0

         LNK2 = LNK2 + log(1.0E0-SAL*0.001005)

         K2 = exp (LNK2)
      
      ELSE

C Roy et al (1993), Millero (1995). Salinities 5-45 (seawater). Total pH scale.
         LNK1 = 2.83655 - 2307.1266/TEMPK - 1.5529413 * log(TEMPK) +
     +         (-0.20760841 - 4.0484/TEMPK) * SAL**0.5 + 0.08468345 * SAL - 
     +          0.00654208 * SAL**1.5
      
         LNK1 = LNK1 + log(1.0E0-SAL*0.001005)
      
         K1  = exp (LNK1)

C Roy et al (1993), Millero (1995). Salinities 5-45 (seawater). Total pH scale.
         LNK2 = -9.226508 - 3351.6106/TEMPK - 0.2005743 * log(TEMPK) +
     +         (-0.106901773 - 23.9722/TEMPK) * SAL**0.5 + 0.1130822 * SAL -
     +           0.00846934 * SAL**1.5

         LNK2 = LNK2 + log(1.0E0-SAL*0.001005)

         K2 = exp (LNK2)

      END IF

C Dissociation constant of boric acid. [mol/kg solution]
C Dickson (1990b). Salinities 5-45 (seawater). Total pH scale.
      LNKB = (-8966.90 - 2890.53 * SAL**0.5 - 77.942 * SAL +
     +         1.728 * SAL**1.5 - 0.0996 * SAL**2)/TEMPK + (148.0248 + 
     +         137.1942 * SAL**0.5 + 1.62142 * SAL) + (-24.4344 - 
     +         25.085 * SAL**0.5 - 0.2474 * SAL) * log(TEMPK) + 
     +         0.053105 * SAL**0.5 * TEMPK
     
      KB = exp (LNKB)

C Solubility constants of calcium carbonate.  Mucci (1983), Millero (1995), Zeebe and Wolf-Gladrow (2001). [mol^2/kg^2 solution]
C Solubility constant of calcite
      LOGCAL = -171.9065 - 0.077993 * TEMPK + 2839.319/TEMPK + 
     +          71.595 * log10(TEMPK) + (-0.77712 + 0.0028426 * TEMPK + 
     +          178.34/TEMPK) * SAL**0.5 - 0.07711 * SAL + 0.0041249 * 
     +          SAL**1.5
     
      Kcal = 10**LOGCAL
     
C Solubility constant of aragonite
      LOGARG = -171.945 - 0.077993 * TEMPK + 2903.293/TEMPK + 
     +          71.595 * log10(TEMPK) + (-0.068393 + 0.0017276 * 
     +          TEMPK + 88.135/TEMPK) * SAL**0.5 -0.10018 * SAL + 
     +          0.0059415 * SAL**1.5

      Karg = 10**LOGARG

C Solubility of CO2 in seawater (Weiss, 1974)
      LNK0 = -58.0931 + 90.5069 / (TEMPK/100.0) + 22.2940 *
     +        log(TEMPK/100.0) + SAL * (0.027766 - 0.025888 * (TEMPK/100.0) +
     +        0.0050578 * (TEMPK/100.0)**2)
      
      K0 = exp (LNK0)

C ******************************
C Conversion of [g/m3] to [mol/kg solvent]
C Density seawater [kg/l]
      RHOH2O = (1000. + 0.7 * SAL / (1.0-SAL/1000.)
     +       - 0.0061 * (TEMP-4.0) * (TEMP-4.0))/1000.

C Total inorganic carbon [mmol/kg]
      TICM = 1000. * (TIC  * M3TOL) / MC / RHOH2O

C Alkalinity [mmol/kg]
      ALK  = 1000. * (ALKA * M3TOL) / MHCO3 / RHOH2O

C Total borate. Millero (1982) [mmol/kg]
      BT = 1000. * 0.000416 * (SAL/35.0)
      
C Calcium concentration. Millero (1982) [mmol/kg]
      Ca = 1000. * 0.01028 * (SAL/35.0)

C
C Solve using solvesaphe from Munhoven (2013)
C

C Set temperature and salinity
      CALL SETUP_API4PHTOT(DBLE(TEMPK), DBLE(SAL), 1.0D0)
C First try the fast poly solver
      AHPLUS = SNGL(SOLVE_ACBW_POLYFAST(DBLE(ALK), DBLE(TICM), DBLE(BT)))
      IF (AHPLUS .LT. 0.0d0) THEN
C If not succesfull try the normal poly solver
          AHPLUS = SNGL(SOLVE_ACBW_POLY(DBLE(ALK), DBLE(TICM), DBLE(BT)))
          IF (AHPLUS .LT. 0.0d0) THEN
C If still not succesfull use the robust solver
              AHPLUS = SNGL(SOLVE_ACBW_GENERAL(DBLE(ALK), DBLE(TICM), DBLE(BT)))
          ENDIF
      ENDIF    
      PH = -LOG10(AHPLUS)
      PH = MAX(PH_MIN,PH)
      PH = MIN(PH_MAX,PH)
      AHPLUS = 10.**(-PH)

      
C ---- SPECCARB ----
C ---- Calculation of CO2 ----
      CO2  = ((TICM * MMTOM * AHPLUS**2)/(AHPLUS**2 + K1 * AHPLUS + K1 * K2))
     +      * MCO2 * RHOH2O / M3TOL
 
C ---- Calculation of pCO2 ----
      pCO2water = (CO2/K0 * 1000000) * M3TOL / MCO2

C ---- Calculation of HCO3 ----
      HCO3 = ((TICM * MMTOM * K1 * AHPLUS)/(AHPLUS**2 + K1 * AHPLUS + K1 * K2))
     +      * MC * RHOH2O / M3TOL
      
C ---- Calculation of CO3 ---- 
      CO3  = ((TICM * MMTOM * K1 * K2)/(AHPLUS**2 + K1 * AHPLUS + K1 * K2))
     +      * MC * RHOH2O / M3TOL
      
C ---- SPECBOR ----
C ---- Calculation of BOH4 ----      
      BOH4 = (BT * MMTOM * KB/(AHPLUS + KB))
     +      * MB * RHOH2O / M3TOL

C ---- Saturation state of calcium carbonate ----
C ---- Calculation of saturation state of calcite ----
      Satcal = (Ca * MMTOM * CO3 /Kcal)* M3TOL / MC / RHOH2O
      
C ---- Calculation of saturation state of aragonite ----
      Satarg = (Ca * MMTOM * CO3 /Karg)* M3TOL / MC / RHOH2O

C
C---- Output --------------------
C
      PMSA(IP7)  = PH
      PMSA(IP8)  = CO2
      PMSA(IP9)  = pCO2water
      PMSA(IP10) = HCO3
      PMSA(IP11) = CO3
      PMSA(IP12) = Satcal
      PMSA(IP13) = Satarg
      PMSA(IP14) = BOH4
C
      ENDIF
C
C---- Pointers ophogen ( altijd buiten de if's op kenmerk) in de segment loop
C
10    IP1    = IP1    + INCREM (  1 )
      IP2    = IP2    + INCREM (  2 )
      IP3    = IP3    + INCREM (  3 )
      IP4    = IP4    + INCREM (  4 )
      IP5    = IP5    + INCREM (  5 )
      IP6    = IP6    + INCREM (  6 )
      IP7    = IP7    + INCREM (  7 )
      IP8    = IP8    + INCREM (  8 )
      IP9    = IP9    + INCREM (  9 )
      IP10   = IP10   + INCREM (  10 )
      IP11   = IP11   + INCREM (  11 )
      IP12   = IP12   + INCREM (  12 )
      IP13   = IP13   + INCREM (  13 )
      IP14   = IP14   + INCREM (  14 )
c
 9000 CONTINUE
c
C
      RETURN
      END
