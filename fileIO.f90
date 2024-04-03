module fileIO
    use type_defination
    use variable_defination
    
    implicit none
    contains
    
     subroutine locate_label(fileID, label, find_Boolean, maxline, rewind_Boolean)
        integer,intent(in):: fileID
        integer::error_Boolean, i
        integer,optional:: find_Boolean, rewind_Boolean, maxline
        character(len=200)::cha_temp_200
        character(len=*)::label
        
        
        if (.not.present(rewind_Boolean)) then
            rewind(fileID)
        elseif (rewind_Boolean==1) then
            rewind(fileID)
        end if
    
        if (.not.present(maxline)) then
            do while (.true.)
                read(fileID, "(a)", iostat=error_Boolean) cha_temp_200
                if (error_Boolean/=0) exit
                if (index(cha_temp_200, label)/=0) then
                    backspace(fileID)
                    if (present(find_Boolean)) find_Boolean=1
                    return
                end if 
            end do
        else
            do i=1, maxline
                read(fileID, "(a)", iostat=error_Boolean) cha_temp_200
                if (error_Boolean/=0) exit
                if (index(cha_temp_200, label)/=0) then
                    backspace(fileID)
                    if (present(find_Boolean)) find_Boolean=1
                    return
                end if
            end do
        end if
        if (present(find_Boolean)) find_Boolean=0
    end subroutine locate_label
     
    subroutine readfch(filename, elec_num, alpha_elec_num, beta_elec_num, atoms_num, basis_num, MO_num, shell_num, primitive_shell_num,&
        atoms, gtfs, shell_type, shell_contraction, shell2atom, primitive_exponents, contraction_coeff, SPcontraction_coeff,&
        aMOcoeff, bMOcoeff)
        
        character(len=200), intent(in):: filename
        
        type(atom_type), allocatable, intent(out)::atoms(:)
        type(gtf_type), allocatable, intent(out)::gtfs(:)
        
        integer, intent(out):: elec_num, alpha_elec_num, beta_elec_num, atoms_num, basis_num, MO_num, shell_num, primitive_shell_num
        integer, allocatable, intent(out):: shell_type(:), shell_contraction(:), shell2atom(:)
        real*8, allocatable, intent(out):: primitive_exponents(:), contraction_coeff(:), SPcontraction_coeff(:)
        real*8, allocatable, intent(out):: aMOcoeff(:,:), bMOcoeff(:,:)
        
        
        integer:: i, iMO, iBasis, beta_Boolean=1, find_Boolean
        
        !initialize the parameters
        call init_fch()
        
        open(10, file=filename, status='old', action='read')
        
            !electron number
            call locate_label(10,'Number of electrons')
            read(10,"(49x,i12)") elec_num
            read(10,"(49x,i12)") alpha_elec_num
            read(10,"(49x,i12)") beta_elec_num
    
            !Atomic index in peridic table
            call locate_label(10,'Atomic numbers')
            read(10,"(49x,i12)") atoms_num
            allocate(atoms(atoms_num))
            read(10,"(6i12)") (atoms(i)%atom_index, i=1,atoms_num)
    
            !Atomic name
            do i=1,atoms_num
            atoms(i)%atom_name=periodic_table(atoms(i)%atom_index)
            end do
    
            !Atomic cartesian coordinates(standard orientation)
            call locate_label(10,'Current cartesian coordinates')
            read(10,*)
            read(10,"(5(1PE16.8))") (atoms(i)%x, atoms(i)%y, atoms(i)%z, i=1,atoms_num)
    
            !Atomic charge
            call locate_label(10,'Nuclear charges')
            read(10,*)
            read(10,"(5(1PE16.8))") (atoms(i)%charge, i=1,atoms_num)
    
            !Number of basis functions
            call locate_label(10,'Number of basis functions')
            read(10,"(49x,i12)") basis_num
    
            !Number of independent functions
            call locate_label(10,'Number of independent functions')
            read(10,"(49x,i12)") MO_num
    
            !Number of Shell(Contracted)
            call locate_label(10,'Shell types')
            read(10,"(49x,i12)") shell_num
            allocate(shell_type(shell_num))
            read(10,"(6i12)") (shell_type(i), i=1, shell_num)
    
            !Number of primitives per shell
            call locate_label(10,'Number of primitives per shell')
            read(10,*)
            allocate(shell_contraction(shell_num))
            read(10,"(6i12)") (shell_contraction(i), i=1, shell_num)
    
            !Mapping Shell to atoms(the index of atoms(:), not the index in perdic table)
            call locate_label(10,'Shell to atom map')
            read(10,*)
            allocate(shell2atom(shell_num))
            read(10,"(6i12)") (shell2atom(i), i=1, shell_num)
    
            !Parameters of primitive shells
            call locate_label(10,'Primitive exponents')
            read(10,"(49x,i12)") primitive_shell_num
            allocate(primitive_exponents(primitive_shell_num))
            read(10,"(5(1PE16.8))") (primitive_exponents(i), i=1, primitive_shell_num)
    
            call locate_label(10,'Contraction coefficients')
            read(10,*)
            allocate(contraction_coeff(primitive_shell_num))
            read(10,"(5(1PE16.8))") (contraction_coeff(i), i=1, primitive_shell_num)
    
            call locate_label(10,'P(S=P) Contraction coefficients', find_Boolean)
            if (find_Boolean==1) then
                read(10,*)
                allocate(SPcontraction_coeff(primitive_shell_num))
                read(10,"(5(1PE16.8))") (SPcontraction_coeff(i), i=1, primitive_shell_num)
            end if
        
        
            allocate(aMOcoeff(MO_num, basis_num))
            allocate(bMOcoeff(MO_num, basis_num))
    
            call locate_label(10,'lpha MO coefficients')
            read(10,*)
            read(10,"(5(1PE16.8))") ((aMOcoeff(iMO, iBasis), iBasis=1, basis_num), iMO=1, MO_num)
    
            !For restrict orbital, beta=alpha
            call locate_label(10,'eta MO coefficients', find_Boolean)
            if (find_Boolean==1) then
                read(10,*)
                read(10,"(5(1PE16.8))") ((bMOcoeff(iMO, iBasis), iBasis=1, basis_num), iMO=1, MO_num)
            else
                beta_Boolean=0
                bMOcoeff=aMOcoeff
            end if
        
            !!!!!!!!!!!!!!!!!!!!!!!!!TODO: add transformation from spherical to cartisian
        
        close(10)
     
    end subroutine readfch
    
end module fileIO
