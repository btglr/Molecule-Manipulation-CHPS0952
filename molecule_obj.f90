module molecule_obj
    use atom_obj

    type molecule
        type(atom), dimension(:), allocatable, private :: atoms
        integer, private :: numberOfAtoms
    contains
        procedure :: displayMolecule
        procedure :: initMolecule
        procedure :: removeMolecule
        procedure :: getNumberOfAtoms
        procedure :: getAtom
        procedure :: setAtom
        procedure :: computeRMSD
        generic :: write(formatted) => displayMolecule
    end type molecule
contains
    subroutine initMolecule(m, size)
        class(molecule), intent(inout) :: m
        integer, intent(in) :: size
        integer :: ok

        allocate(m%atoms(size), stat = ok)

        if(ok /= 0) then
            print '(a)', "Error during the array allocation ==> aborting"
            stop 666
        end if

        m%numberOfAtoms = 0
    end subroutine initMolecule

    subroutine addAtom(m, at)
        class(molecule), intent(inout) :: m
        type(atom), intent(in) :: at

        m%numberOfAtoms = m%numberOfAtoms + 1
        m%atoms(m%numberOfAtoms) = at
    end subroutine addAtom

    subroutine displayMolecule(m, unit, iotype, v_list, iostat, iomsg)
        class(molecule), intent(in) :: m
        integer, intent(in) :: unit             ! Internal unit to write to.
        character(*), intent(in) :: iotype      ! LISTDIRECTED or DTxxx
        integer, intent(in) :: v_list(:)        ! parameters from fmt spec.
        integer, intent(out) :: iostat          ! non zero on error, etc.
        character(*), intent(inout) :: iomsg    ! define if iostat non zero.
        write (unit, *, IOSTAT = iostat, IOMSG = iomsg) m%atoms
    end subroutine displayMolecule

    subroutine removeMolecule(m)
        class(molecule), intent(inout) :: m
        integer :: ok

        deallocate(m%atoms, stat = ok)

        if(ok /= 0) then
            print '(a)', "Error during the atoms array deallocation"
            stop 666
        end if
    end subroutine removeMolecule

    subroutine translateMolecule(originalMolecule, translatedMolecule, translationVector)
        class(molecule), intent(in) :: originalMolecule
        class(molecule), intent(inout) :: translatedMolecule
        real, dimension(:), intent(in) :: translationVector
        integer :: atomIndex

        print *, "Translation vector:", translationVector

        call initMolecule(translatedMolecule, originalMolecule%numberOfAtoms)

        do atomIndex = 1, originalMolecule%numberOfAtoms
            call addAtom(translatedMolecule, originalMolecule%atoms(atomIndex))
            call translateAtom(translatedMolecule%atoms(atomIndex), translationVector)
        end do
    end subroutine translateMolecule

    integer function getNumberOfAtoms(m) result(numberOfAtoms)
        class(molecule), intent(in) :: m

        numberOfAtoms = m%numberOfAtoms
    end function getNumberOfAtoms

    type(atom) function getAtom(m, atomIndex) result(at)
        class(molecule), intent(in) :: m
        integer, intent(in) :: atomIndex

        at = m%atoms(atomIndex)
    end function getAtom

    subroutine setAtom(m, atomIndex, at)
        class(molecule), intent(inout) :: m
        integer, intent(in) :: atomIndex
        type(atom), intent(in) :: at

        m%atoms(atomIndex) = at
    end subroutine setAtom

    real function computeRMSD(m1, m2, type) result(rmsd)
        class(molecule), intent(in) :: m1
        class(molecule), intent(in) :: m2
        character(*), intent(in) :: type

        real, dimension(:, :), allocatable :: difference
        integer :: numberOfAtoms, numberOfHeavyAtoms, atomIndex
        logical :: sameMolecule

        if (getNumberOfAtoms(m1) /= getNumberOfAtoms(m2)) then
            print '(a)', "Error: the molecules do not have the same number of atoms ==> aborting"
            stop 20
        end if

        numberOfAtoms = getNumberOfAtoms(m1)

        sameMolecule = .TRUE.

        do atomIndex = 1, numberOfAtoms
            if (getAtomName(getAtom(m1, atomIndex)) /= getAtomName(getAtom(m2, atomIndex))) then
                print *, getAtomName(getAtom(m1, atomIndex)), getAtomName(getAtom(m2, atomIndex))
                sameMolecule = .FALSE.
                exit
            end if
        end do

        if (sameMolecule.eqv..FALSE.) then
            print '(a)', "Error: the molecules are not identical ==> aborting"
            stop 20
        end if

        if (type == "standard") then
            print '(a25,i8)', "Number of atoms:", numberOfAtoms
            allocate(difference(numberOfAtoms, 3))

            do atomIndex = 1, numberOfAtoms
                difference(atomIndex, :) = getCoordinates(getAtom(m1, atomIndex)) - getCoordinates(getAtom(m2, atomIndex))
            end do
        else if (type == "heavy") then
            numberOfHeavyAtoms = 0

            do atomIndex = 1, numberOfAtoms
                if (getAtomName(getAtom(m1, atomIndex)) /= "H") then
                    numberOfHeavyAtoms = numberOfHeavyAtoms + 1
                end if
            end do

            allocate(difference(numberOfHeavyAtoms, 3))

            do atomIndex = 1, numberOfHeavyAtoms
                difference(atomIndex, :) = getCoordinates(getAtom(m1, atomIndex)) - getCoordinates(getAtom(m2, atomIndex))
            end do

            numberOfAtoms = numberOfHeavyAtoms
            print '(a25,i8)', "Number of heavy atoms:", numberOfAtoms
        else
            stop 20
        end if

        print '(a25,f8.0)', "Sum:", sum(difference ** 2)

        rmsd = sqrt(sum(difference ** 2) / numberOfAtoms)
    end function computeRMSD
end module molecule_obj