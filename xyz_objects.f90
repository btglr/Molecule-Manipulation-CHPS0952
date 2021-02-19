module xyz_handler
    use atom_obj
    use molecule_obj

contains

    subroutine writeXYZ(m, filename)
        class(molecule), intent(in) :: m
        character(*), intent(in) :: fileName

        character(len = 8) :: nbAtomsC, atomName
        character(len = 512) :: comment
        integer :: ok, atomIndex, unit, nbAtoms
        logical :: exist
        type(atom) :: at
        real, dimension(3) :: atomCoordinates, translationVector, rotationVector

        unit = 11
        inquire(file = fileName, exist = exist)

        if (exist) then
            open(unit, file = fileName, status = 'old', action = 'write', iostat = ok)
        else
            open(unit, file = fileName, status = 'new', action = 'write', iostat = ok)
        end if

        if(ok /= 0) then
            print '(a, 4x, a)', "Error occurred while opening file", fileName
            stop 20
        end if

        nbAtoms = getNumberOfAtoms(m)
        translationVector = getTranslationVector(m)
        rotationVector = getRotationVector(m)

        write(comment, '(a15,x,3(f8.3),x,a15,x,3(f8.3))') "Translation: ", translationVector, ", Rotation: ", rotationVector
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
        character(*), intent(in) :: fileName
        integer, intent(in) :: unit

        character(len = 128) :: line
        character(len = 2) :: atomName
        integer :: ok, numberOfAtoms, end
        logical :: exist
        real :: x, y, z
        type(atom) :: currentAtom
        real, dimension(3) :: translationVector, rotationVector

        inquire(file = fileName, exist = exist)

        if (exist) then
            open(unit, file = fileName, status = 'old', iostat = ok)
        else
            print '(a, 4x, a)', "File ", fileName, " doesn't exist"
            stop 20
        end if

        if(ok /= 0) then
            print '(a, 4x, a)', "Error occurred while opening file", fileName
            stop 20
        end if

        read(unit, '(i8)', iostat = end) numberOfAtoms

        call initMolecule(m, numberOfAtoms)

        ! La ligne du commentaire
        read(unit, '(16x,3(f8.3),17x,3(f8.3))', iostat = end) translationVector, rotationVector

        call setTranslationVector(m, translationVector)
        call setRotationVector(m, rotationVector)

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

end module xyz_handler