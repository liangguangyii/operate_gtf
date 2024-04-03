program main
    use type_defination
    use variable_defination
    use fileIO
    use math
    implicit none
    
    integer:: fileNum, ifile, iMO, iBasis, jBasis
    
    character(len=200):: chartemp
    character(len=1), allocatable:: operatorList(:)
    character(len=200), allocatable:: filenameList(:)
    
    type(atom_type), allocatable::atoms(:)
    type(gtf_type), allocatable::gtfs(:)
        
    integer:: elec_num, alpha_elec_num, beta_elec_num, atoms_num, basis_num, MO_num, shell_num, primitive_shell_num
    integer, allocatable:: shell_type(:), shell_contraction(:), shell2atom(:)
    real*8, allocatable:: primitive_exponents(:), contraction_coeff(:), SPcontraction_coeff(:)
    real*8, allocatable:: aMOcoeff(:,:), bMOcoeff(:,:)
    
    
    real*8, allocatable:: Palpha(:,:), Pbeta(:,:), Ptot(:,:), Pspin(:,:)
    
    !the input file format
    ! file num
    ! filename1
    ! operator1, filename2
    ! operator2, filename3
    ! operator3, filename4
    ! ...
    
    open(unit=10, file='input.txt', status='old', action='read')
        read(10,*) fileNum
        allocate(filenameList(fileNum+1), operatorList(fileNum))
        
        read(10,"(a)") filenameList(1)
        
        do ifile=1, fileNum
            read(10,"(a)") chartemp 
            read(chartemp,*) operatorList(ifile), filenameList(ifile+1)
        end do
        
    close(10)
    
    call readfch(filenameList(1), elec_num, alpha_elec_num, beta_elec_num, atoms_num, basis_num, MO_num, shell_num, primitive_shell_num,&
        atoms, gtfs, shell_type, shell_contraction, shell2atom, primitive_exponents, contraction_coeff, SPcontraction_coeff,&
        aMOcoeff, bMOcoeff)
    
    
    call genDensity(basis_num, alpha_elec_num, beta_elec_num, aMOcoeff, bMOcoeff,&
        Ptot, Pspin, Palpha, Pbeta)
    
    !do iMO = 1, MO_num
    !    rtemp1 = 0D0
    !    rtemp2 = 0D0
    !    do iBasis = 1, basis_num
    !        rtemp1 = rtemp1 + aMOcoeff(iMO, iBasis) ** 2
    !        rtemp2 = rtemp2 + bMOcoeff(iMO, iBasis) ** 2
    !    end do
    !    write(*,*) iMO, rtemp1, rtemp2
    !end do
    


    
end program main