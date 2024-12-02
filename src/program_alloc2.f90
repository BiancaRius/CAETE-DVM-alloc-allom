!Provisory program for testingn allocation2 module

program program_alloc
    
    use types 
    use alloc2

    !declaring variables:

    real(r_4) :: npp = 1.
    !C that enters in the logic (kgC/m2) for each PLS (attention it is not individual yet)
    real(r_8) :: scl1 = 1.5 !kgC/m2 initial C leaf input 
    real(r_8) :: sca1 = 0.0 !kgC/m2 initial C awood input
    real(r_8) :: scf1 = 1.  !kgC/m2 initial C root input    
    real(r_8) :: scs1 = 1.  !kgC/m2 initial C sapwood input  
    real(r_8) :: sch1 = 1.  !kgC/m2 initial C heart input  
    real(r_8) :: sct1 = 1.  !kgC/m2 initial C storage input  
    real(r_8) :: dens1 = 1.  !kgC/m2 initial C storage input  




    !C that exit in the logic (kgC/m2) for each PLS (attention it is not individual yet)
    real(r_8) :: scl2 
    real(r_8) :: sca2
    real(r_8) :: scf2
    real(r_8) :: scs2
    real(r_8) :: sch2
    real(r_8) :: sct2


    ! sca1 - sum of heartwood and sapwood
    sca1 = scs1 + sch1

    
    call allocation2(npp, scl1, sca1, scf1, scs1, sch1, sct1, dens1, &
    scl2, sca2, scf2, scs2, sch2, sct2)


end program program_alloc