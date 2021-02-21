module atom_obj
    type atom
        character(:), allocatable, private :: atomName
        real, dimension(3), private :: coordinates
    contains
        procedure :: displayAtom
        procedure :: initAtom
        procedure :: atomMinusAtom
        generic :: write(formatted) => displayAtom
        generic :: operator(-) => atomMinusAtom
    end type atom
contains
    subroutine initAtom(at, atomNameIn, coordinates)
        class(atom), intent(inout) :: at
        character(len = *), intent(in) :: atomNameIn
        real, dimension(3), intent(in) :: coordinates

        at%atomName = atomNameIn
        at%coordinates = coordinates
    end subroutine initAtom

    subroutine displayAtom(at, unit, iotype, v_list, iostat, iomsg)
        class(atom), intent(in) :: at
        integer, intent(in) :: unit             ! Internal unit to write to.
        character(*), intent(in) :: iotype      ! LISTDIRECTED or DTxxx
        integer, intent(in) :: v_list(:)        ! parameters from fmt spec.
        integer, intent(out) :: iostat          ! non zero on error, etc.
        character(*), intent(inout) :: iomsg    ! define if iostat non zero.
        write (unit, 999, IOSTAT = iostat, IOMSG = iomsg) at%atomName, at%coordinates
        999 format(a2, 1x, 3(f8.3, 1x), /)
    end subroutine displayAtom

    subroutine translateAtom(at, translationVector)
        class(atom), intent(inout) :: at
        real, dimension(:), intent(in) :: translationVector

        call setCoordinates(at, at%coordinates + translationVector)
    end subroutine translateAtom

    character function getAtomName(at) result(atomName)
        type(atom), intent(in) :: at

        atomName = at%atomName
    end function getAtomName

    function getCoordinates(at) result(coordinates)
        class(atom), intent(in) :: at
        real, dimension(3) :: coordinates

        coordinates = at%coordinates
    end function getCoordinates

    subroutine setCoordinates(at, newCoordinates)
        class(atom), intent(inout) :: at
        real, dimension(3), intent(in) :: newCoordinates

        at%coordinates = newCoordinates
    end subroutine setCoordinates

    function atomMinusAtom(at1, at2) result(coordinates)
        class(atom), intent(in) :: at1, at2
        real, dimension(3) :: coordinates

        coordinates = getCoordinates(at1) - getCoordinates(at2)
    end function atomMinusAtom

    subroutine rotateAtom(at, rotationMatrix, tail)
        class(atom), intent(inout) :: at
        real, dimension(3,3), intent(in) :: rotationMatrix
        real, dimension(3), intent(in) :: tail

        real, dimension(3) :: newCoordinates, currentCoordinates, translatedPoint

        currentCoordinates = getCoordinates(at)
        translatedPoint = currentCoordinates - tail
        translatedPoint = matMul(rotationMatrix, translatedPoint)
        newCoordinates = translatedPoint + tail

        call setCoordinates(at, newCoordinates)
    end subroutine rotateAtom
end module atom_obj