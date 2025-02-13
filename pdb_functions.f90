module pdb_functions
    use atom_obj
    use molecule_obj

contains

    subroutine readPDB(m, pdbFile)
        class(molecule), intent(inout) :: m
        character(len = 128), intent(in) :: pdbFile

        type(atom) :: currentAtom
        character(len = 2) :: elementSymbol
        character(len = 128) :: line
        integer :: arraySize, end, idx, ok, unit
        real :: x, y, z

        unit = 10
        open(unit, file = pdbFile, iostat = ok, status = 'old')

        if(ok /= 0) then
            print '(a, 4x, a)', 'Error during opening', pdbFile
            stop 20
        end if

        ! Lecture du nombre de lignes commençant par 'ATOM'
        arraySize = 0
        do
            read(unit, *, iostat = end) line
            if(end /= 0)then
                exit
            else
                if(line(1:6) == 'ATOM') then
                    arraySize = arraySize + 1
                end if
            end if
        end do

        print '(a40, i8)', 'Number of atoms: ', arraySize

        rewind(unit)
        call initMolecule(m, arraySize)

        do
            read(unit, '(a)', iostat = end) line
            if(end /= 0) then
                exit
            else
                if(line(1:6) == 'ATOM') then
                    read(line, '(6x, i5, 19x, 3(f8.3), 23x, a2)') idx, x, y, z, elementSymbol
                    call initAtom(currentAtom, elementSymbol, ([x, y, z]))
                    call addAtom(m, currentAtom)
                end if
            end if
        end do

        close(unit)
    end subroutine readPDB

end module pdb_functions