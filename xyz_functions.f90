module xyz_functions
    use atom_obj
    use molecule_obj

contains

    subroutine writeXYZ(m, filename)
        class(molecule), intent(in) :: m
        character(*), intent(in) :: fileName

        type(atom) :: at
        character(len = 8) :: atomName, nbAtomsC
        character(len = 512) :: comment
        integer :: atomIndex, nbAtoms, ok, unit
        logical :: exist
        real :: globalRotationAngle, internalRotationAngle
        real, dimension(3) :: atomCoordinates, translationVector

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
        translationVector = getTranslationVector(m)
        globalRotationAngle = getGlobalRotationAngle(m)
        internalRotationAngle = getInternalRotationAngle(m)

        write(comment, '(a20, 3(f8.3), 1x, a30, f8.3, 1x, a30, f5.3)') &
                'Translation: ', translationVector, &
                'Global rotation angle: ', globalRotationAngle, &
                'Internal rotation angle: ', internalRotationAngle
        write (nbAtomsC, '(i8)') nbAtoms

        write(unit, '(a)') adjustl(nbAtomsC)
        write(unit, '(a)') comment

        do atomIndex = 1, nbAtoms
            at = getAtom(m, atomIndex)
            atomName = getAtomName(at)
            atomCoordinates = getCoordinates(at)

            write(unit, '(a2, 2x, 3(f8.3, 2x))') adjustl(atomName), atomCoordinates(1), atomCoordinates(2), atomCoordinates(3)
        end do

        close(unit)

    end subroutine writeXYZ

    subroutine readXYZ(m, filename, unit)
        class(molecule), intent(inout) :: m
        character(*), intent(in) :: filename
        integer, intent(in) :: unit

        type(atom) :: currentAtom
        character(len = 2) :: atomName
        character(len = 128) :: line
        integer :: ok, numberOfAtoms, end
        logical :: exist
        real :: globalRotationAngle, internalRotationAngle, x, y, z
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
        read(unit, '(20x,3(f8.3),31x,f8.3,31x,f5.3)', iostat = end) translationVector, globalRotationAngle, internalRotationAngle

        call setTranslationVector(m, translationVector)
        call setGlobalRotationAngle(m, globalRotationAngle)
        call setInternalRotationAngle(m, internalRotationAngle)

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