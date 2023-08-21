import os
from pathlib import Path

descrp = "This script creates a global.f90 file with an asked iniital NPLS number"


s = int(input("Number of initial random EV (10 - 3000): "))

global_f90 = f"""
! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! This program is based on the work of those that gave us the INPE-CPTEC-PVM2 model

module types
   implicit none

   ! FOR THE GNU FORTRAN COMPILER
   integer,parameter,public :: l_1 = 2  ! standard Logical type
   integer,parameter,public :: i_2 = 2  ! 16 bits integer
   integer,parameter,public :: i_4 = 4  ! 32 bits integer
   integer,parameter,public :: r_4 = 4  ! 32 bits float
   integer,parameter,public :: r_8 = 8  ! 64 bits float

end module types


module global_par
   
   use types
   
   implicit none
   real(r_4),parameter,public :: q10 = 1.4
   real(r_4),parameter,public :: h = 1.0                         ! soil layer thickness (meters)
   real(r_4),parameter,public :: diffu = 1.036800e14             ! soil thermal diffusivity (m2/mes)
   real(r_4),parameter,public :: tau = (h**2)/(2.0*diffu)        ! e-folding times (months)
   real(r_4),parameter,public :: rcmax = 5000.0                  ! ResistÊncia estomática máxima s/m
   real(r_4),parameter,public :: rcmin = 100                     ! ResistÊncia estomática mínima s/m
   real(r_8),parameter,public :: cmin = 1.0D-6                   ! Minimum to survive kg m-2
   ! real(r_4),parameter,public :: wmax = 500.0                    ! Maximum water soil capacity (Kg m-2)
   real(r_8),parameter,public :: vcmax = 50.0D-5 
   real(r_8),parameter,public :: csru = 0.5D0                    ! Root attribute
   real(r_8),parameter,public :: alfm = 1.391D0                  ! Root attribute
   real(r_8),parameter,public :: gm = 3.26D0 * 86400D0           ! (*86400 transform s/mm to dia/mm)
   real(r_8),parameter,public :: sapwood = 0.05D0                ! Fraction of wood tissues that are sapwood
   real(r_4),parameter,public :: ks = 0.25                       ! P Sorption
   integer(i_4),parameter,public :: npls = {s}                  ! Number of Plant Life Strategies-PLSs simulated (Defined at compile time)
   integer(i_4),parameter,public :: ntraits = 17                 ! Number of traits for each PLS
   real(r_8),parameter, public :: ncl = (1.0/60.0)          !(gN/gC) used in maintenance respiration (from Scheiter & Higgins 2008)
   real(r_8),parameter, public :: ncf = (1.0/60.0)          !(gN/gC)
   real(r_8),parameter, public :: ncs = (1.0/330.0)         !(gN/gC) from Smith et al 2001

   !allometric parameters
   real(r_8), parameter, public :: dwood = 0.74*1.D6 !provisory
   real(r_8), parameter, public :: pi = 3.1415926536
   real(r_8), parameter, public :: k_allom2 = 29.
   real(r_8), parameter, public :: k_allom3 = 0.4
   real(r_8), parameter, public :: klatosa =  8000.0
   real(r_8), parameter, public :: sla_allom = 0.023!15.36 !0.023!provisory
   real(r_8), parameter, public :: ltor = 0.77302587552347657

   !allocation parameters
   real(r_8), parameter, public :: tol = 0.0000001
   real(r_8), parameter, public :: reinickerp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)
   real(r_8), parameter, public :: nseg = 20
   real(r_8), parameter, public :: xacc =  0.1     !x-axis precision threshold for the allocation solution
   real(r_8), parameter, public :: yacc =  1.e-10  !y-axis precision threshold for the allocation solution

   !turnover parameters - see parameters in Bucley & Roberts  2006 - Tree physiology
   !provisory
   real(r_8), parameter, public :: l_turnover = 1./4. !Sitch et al 2003
   real(r_8), parameter, public :: r_turnover = 1./2. !Sitch et al 2003
   real(r_8), parameter, public :: s_turnover = 1./20. !0.05 !Sitch et al 2003
   !ATTENTION TO HEARt TURNOVER
   real(r_8), parameter, public :: h_turnover = 1./60.!1/50. !Sitch et al 2003


end module global_par


module photo_par
   use types, only : r_8
   implicit none


   real(r_8),public, parameter ::       &
        a   = 0.8300D0       ,&          !Photosynthesis co-limitation coefficient
        a2  = 0.930D0        ,&          !Photosynthesis co-limitation coefficient
        p3  = 21200.0D0      ,&          !Atmospheric oxygen concentration (Pa)
        p4  = 0.080D0        ,&          !Quantum efficiency (mol electrons/Ein)
        p5  = 0.150D0        ,&          !Light scattering rate
        p6  = 2.0D0          ,&          !Parameter for jl
        p7  = 0.50D0         ,&          !Ratio of light limited photosynthesis to Rubisco carboxylation
        p8  = 5200.0D0       ,&          !Photo-respiration compensation point
        p9  = 0.570D0        ,&          !Photosynthesis co-limitation coefficient
        p10 = 0.100D0        ,&          !Q10 function
        p11 = 25.0D0         ,&          !Q10 function reference temperature (oC)
        p12 = 30.0D0         ,&          !Michaelis-Menten constant for CO2 (Pa)
        p13 = 2.100D0        ,&          !Michaelis-Menten constant for CO2
        p14 = 30000.0D0      ,&          !Michaelis-Menten constant for O2 (Pa)
        p15 = 1.20D0         ,&          !Michaelis-Menten constant for O2
        p19 = 0.90D0         ,&          !Maximum ratio of internal to external CO2
        p20 = 0.10D0         ,&          !Critical humidity deficit (kg/kg)
        p25 = 9.1D-5         ,&          !Maximum gross photosynthesis rate (molCO2/m2/s)
        p26 = 0.50D0         ,&          !light extinction coefficient for IPAR/sun (0.5/sen90)
        p27 = 1.50D0         ,&          !light extinction coefficient for IPAR/shade (0.5/sen20)
        alphap = 0.0913D0    ,&          ! 0.0913 parameter for v4m. Hard to explain. See Chen et al. 1994
        vpm25 =  85.0D0      ,&          ! µmol m-2 s-1 PEPcarboxylase CO2 saturated rate of carboxilation at 25°C
        h_vpm = 185075.0D0   ,&          ! Arrhenius eq. constant
        s_vpm = 591.0D0      ,&          ! Arrhenius eq. constant
        r_vpm = 8.314D0      ,&          ! Arrhenius eq. constant
        e_vpm = 60592.0D0    ,&          ! Arrhenius eq. constant
        kp25 = 82.0D0                    ! µmol mol-1 (ppm)  MM constant PEPcase at
end module photo_par"""



root = Path(os.getcwd()).resolve()
with open(f"{root}/global.f90", 'w') as fh:
    fh.write(global_f90)
    