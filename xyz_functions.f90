module xyz_functions
    use atom_obj
    use molecule_obj

contains

    subroutine writeXYZ(m, filename)
        class(molecule), intent(in) :: m
        character(*), intent(in) :: fileName

        type(atom) :: at
        type(transformation), dimension(:), allocatable :: transformations
        character(len = 8) :: atomName, nbAtomsC
        character(len = 75) :: stringTransformation
        character(len = 2048) :: comment, tmp
        integer :: atomIndex, nbAtoms, ok, transformationIndex, unit
        logical :: exist
        real, dimension(3) :: atomCoordinates

        unit = 11
        inquire(file = fileName, exist = exist)

        if (exist) then
            open(unit, file = fileName, status = 'old', action = 'write', iostat = ok)
        else
            open(unit, file = fileName, status = 'new', action = 'write', iostat = ok)
        end if

        if(ok /= 0) then
            print '(a, 4x, a)', 'Error occurred while opening file', fileName
            stop 20
        end if

        nbAtoms = getNumberOfAtoms(m)
        transformations = getTransformations(m)

        comment = ''

        do transformationIndex = size(transformations), 1, -1
            write(stringTransformation, *) getTransformation(m, transformationIndex)
            write(tmp, '(a75,a)') stringTransformation, trim(comment)
            comment = tmp
        end do

        write(tmp, '(i4,a)') size(transformations), trim(comment)
        comment = tmp

        write(comment, '(a2048)') comment
        write(nbAtomsC, '(i8)') nbAtoms

        write(unit, '(a)') adjustl(nbAtomsC)
        write(unit, '(a)') comment

        do atomIndex = 1, nbAtoms
            at = getAtom(m, atomIndex)
            atomName = getAtomName(at)
            atomCoordinates = getCoordinates(at)

            write(unit, '(a2, 2x, 3(f8.3, 2x))') adjustl(atomName), atomCoordinates(1), atomCoordinates(2), atomCoordinates(3)
        end do

        close(unit)

        deallocate(transformations)
    end subroutine writeXYZ

    subroutine readXYZ(m, filename, unit)
        class(molecule), intent(inout) :: m
        character(*), intent(in) :: filename
        integer, intent(in) :: unit

        type(atom) :: currentAtom
        type(transformation) :: t
        character(len = 2) :: atomName
        character(len = 20) :: transformationType
        character(len = 75) :: stringTransformation
        character(len = 128) :: line
        character(len = 2048) :: comment, stringTransformations
        integer :: end, ok, numberOfAtoms, numberOfTransformations, transformationIndex
        integer, dimension(2) :: atomsIndices
        logical :: exist
        real :: x, y, z, angle
        real, dimension(3) :: translationVector

        inquire(file = filename, exist = exist)

        if (exist) then
            open(unit, file = filename, status = 'old', iostat = ok)
        else
            print '(a, 1x, a, a)', 'File', adjustl(trim(filename)), ' doesn''t exist'
            stop 20
        end if

        if(ok /= 0) then
            print '(a, 1x, a)', 'Error occurred while opening file', adjustl(trim(filename))
            stop 20
        end if

        read(unit, '(i8)', iostat = end) numberOfAtoms

        call initMolecule(m, numberOfAtoms)

        ! La ligne du commentaire
        read(unit, '(a2048)', iostat = end) comment

        ! Lecture des transformations du commentaire
        read(comment, '(i4, a)') numberOfTransformations
        read(comment, '(4x, a)') stringTransformations

        do transformationIndex = 1, numberOfTransformations
            read(stringTransformations((((transformationIndex - 1) * 75) + 1):len(stringTransformations)), '(a75)') &
            stringTransformation
            read(stringTransformation, '(1x, a20)') transformationType

            if (adjustl(trim(transformationType)) == 'Translation') then
                read(stringTransformation, '(22x, 3(f8.3, 1x))') translationVector
                call initTranslation(t, translationVector)
            else if (adjustl(trim(transformationType)) == 'Global rotation') then
                read(stringTransformation, '(32x, 2(i5, 1x), 10x, f8.3)') atomsIndices, angle
                call initGlobalRotation(t, atomsIndices, angle)
            else if (adjustl(trim(transformationType)) == 'Internal rotation') then
                read(stringTransformation, '(32x, 2(i5, 1x), 10x, f8.3)') atomsIndices, angle
                call initInternalRotation(t, atomsIndices, angle)
            else
                stop 20
            end if

            call addTransformation(m, t)
        end do

        do
            read(unit, '(a)', iostat = end) line
            if(end /= 0) then
                exit
            else
                read(line, '(a2, 2x, 3(f8.3, 2x))') atomName, x, y, z
                call initAtom(currentAtom, atomName, (/x, y, z/))
                call addAtom(m, currentAtom)
            end if
        end do

        close(unit)
    end subroutine readXYZ

end module xyz_functions