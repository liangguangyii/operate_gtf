module math

    implicit none
    contains
    
    ! ONLY FOR OPENSHELL SYSTEMS
    !************************************************************
    !calculate the density matrix(ONLY SUM UP OCCUPIED orbits)
    ! that means, two density
    ! Total SCF density is the sum of the two densities
    ! Spin SCF density is the difference of the two densities
    subroutine genDensity(basis_num, alpha_elec_num, beta_elec_num, aMOcoeff, bMOcoeff,&
        Ptot, Pspin, Palpha, Pbeta)
        integer,  intent(in):: basis_num, alpha_elec_num, beta_elec_num
        real*8, allocatable, intent(in):: aMOcoeff(:,:), bMOcoeff(:,:)
        
        real*8, allocatable, intent(out):: Ptot(:,:), Pspin(:,:), Palpha(:,:), Pbeta(:,:)
        
        integer:: iBasis, jBasis, iMO
        real*8:: rtemp1, rtemp2
        
        allocate(Palpha(basis_num, basis_num))
        allocate(Pbeta(basis_num, basis_num))
        allocate(Ptot(basis_num, basis_num))
        allocate(Pspin(basis_num, basis_num))  
        
        do iBasis = 1, basis_num
            do jBasis = iBasis, basis_num
                rtemp1 = 0D0
                do iMO = 1, alpha_elec_num
                    rtemp1 = rtemp1 + aMOcoeff(iMO, iBasis) * aMOcoeff(iMO, jBasis)
                end do
                Palpha(iBasis, jBasis) = rtemp1
            end do
        end do
    
        do iBasis = 1, basis_num
            do jBasis = iBasis, basis_num
                rtemp2 = 0D0
                do iMO = 1, beta_elec_num
                    rtemp2 = rtemp2 + bMOcoeff(iMO, iBasis) * bMOcoeff(iMO, jBasis)
                end do
                Pbeta(iBasis, jBasis) = rtemp2
            end do
        end do    
    
        Ptot = Palpha + Pbeta
        Pspin = Palpha - Pbeta        
        
    end subroutine genDensity
    
end module math
