module pdb_objects

    type atom
        character(:), allocatable, private :: atomName
        real, dimension(3), private :: coordinates
    contains
        procedure :: display_atom
        procedure :: init_atom
        generic :: write(formatted) => display_atom
    end type atom

    type molecule
        type(atom), dimension(:), allocatable, private :: atoms
        integer, private :: nbAtoms
    contains
        procedure :: display_molecule
        procedure :: init_molecule
        procedure :: remove_molecule
        procedure :: get_atom_size
        procedure :: get_atom
        generic :: write(formatted) => display_molecule
    end type molecule

contains
    ! Choisir 'a' au lieu de 'at' comme nom ici cause un SEGFAULT... ?
    subroutine init_atom(at, atomNameIn, coordinates)
        class(atom), intent(inout) :: at
        character(len = *), intent(in) :: atomNameIn
        real, dimension(3), intent(in) :: coordinates

        at%atomName = atomNameIn
        at%coordinates = coordinates
    end subroutine init_atom

    ! Choisir 'a' au lieu de 'at' comme nom ici cause un SEGFAULT... ?
    subroutine display_atom(at, unit, iotype, v_list, iostat, iomsg)
        class(atom), intent(in) :: at
        integer, intent(in) :: unit             ! Internal unit to write to.
        character(*), intent(in) :: iotype      ! LISTDIRECTED or DTxxx
        integer, intent(in) :: v_list(:)        ! parameters from fmt spec.
        integer, intent(out) :: iostat          ! non zero on error, etc.
        character(*), intent(inout) :: iomsg    ! define if iostat non zero.
        write (unit, 999, IOSTAT = iostat, IOMSG = iomsg) at%atomName, at%coordinates(1), at%coordinates(2), at%coordinates(3)
        999 format(a4, 1x, 3(f8.3, 1x), /)
    end subroutine display_atom

    subroutine translate_atom(at, translationVector)
        class(atom), intent(inout) :: at
        real, dimension(:), intent(in) :: translationVector

        at%coordinates = at%coordinates + translationVector
    end subroutine translate_atom

    subroutine init_molecule(m, size)
        class(molecule), intent(inout) :: m
        integer, intent(in) :: size
        integer :: ok

        allocate(m%atoms(size), stat = ok)

        if(ok /= 0) then
            print '(a)', 'Error during the array allocation ==> aborting'
            stop 666
        end if

        m%nbAtoms = 0
    end subroutine init_molecule

    subroutine add_atom(m, at)
        class(molecule), intent(inout) :: m
        type(atom), intent(in) :: at

        m%nbAtoms = m%nbAtoms + 1
        m%atoms(m%nbAtoms) = at
    end subroutine add_atom

    subroutine display_molecule(m, unit, iotype, v_list, iostat, iomsg)
        class(molecule), intent(in) :: m
        integer, intent(in) :: unit             ! Internal unit to write to.
        character(*), intent(in) :: iotype      ! LISTDIRECTED or DTxxx
        integer, intent(in) :: v_list(:)        ! parameters from fmt spec.
        integer, intent(out) :: iostat          ! non zero on error, etc.
        character(*), intent(inout) :: iomsg    ! define if iostat non zero.
        write (unit, *, IOSTAT = iostat, IOMSG = iomsg) m%atoms
    end subroutine display_molecule

    subroutine remove_molecule(m)
        class(molecule), intent(inout) :: m
        integer :: ok

        deallocate(m%atoms, stat = ok)

        if(ok /= 0) then
            print '(a)', 'Error during the atoms array deallocation'
            stop 666
        end if
    end subroutine remove_molecule

    integer function get_atom_size(m) result(atomSize)
        class(molecule), intent(in) :: m

        atomSize = m%nbAtoms
    end function get_atom_size

    type(atom) function get_atom(m, atomIndex) result(at)
        class(molecule), intent(in) :: m
        integer, intent(in) :: atomIndex

        at = m%atoms(atomIndex)
    end function get_atom

    character function get_atom_name(at) result(atomName)
        type(atom), intent(in) :: at

        atomName = at%atomName
    end function get_atom_name

    subroutine translate_molecule(originalMolecule, translatedMolecule, translationVector)
        class(molecule), intent(in) :: originalMolecule
        class(molecule), intent(inout) :: translatedMolecule
        real, dimension(:), intent(in) :: translationVector
        integer :: atomIndex

        print *, 'Translation vector:', translationVector

        call init_molecule(translatedMolecule, originalMolecule%nbAtoms)

        do atomIndex = 1, originalMolecule%nbAtoms
            call add_atom(translatedMolecule, originalMolecule%atoms(atomIndex))
            call translate_atom(translatedMolecule%atoms(atomIndex), translationVector)
        end do
    end subroutine translate_molecule

end module pdb_objects