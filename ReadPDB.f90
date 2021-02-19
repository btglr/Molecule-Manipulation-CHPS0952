program ReadPDBObject
    use atom_obj
    use molecule_obj
    use xyz_handler
    implicit none

    integer :: i, end, ok, idx, arraySize, unit
    character(len = 128) :: inputFile, outputFile, line, newFileName
    character(len = 4) :: atomName, elementSymbol
    real :: x, y, z
    type(atom) :: currentAtom
    type(atom), dimension(:), allocatable :: atoms
    type(molecule) :: currentMolecule
    real, dimension(3) :: translationVector

    if(iargc() /= 2) then
        print '(a)', "Please provide a PDB file to open and an XYZ file to write to"
        stop 10
    end if

    call getarg(1, inputFile)
    call getarg(2, outputFile)

    inputFile = trim(inputFile)
    outputFile = trim(outputFile)

    print '(a, a)', "File to read = ", inputFile

    unit = 10
    open(unit, file = inputFile, iostat = ok, status = 'old')

    if(ok /= 0) then
        print '(a, 4x, a)', "Error during opening", inputFile
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

    print *, "Number of atoms:", arraySize

    rewind(unit)
    call initMolecule(currentMolecule, arraySize)

    do
        read(unit, '(a)', iostat = end) line
        if(end /= 0) then
            exit
        else
            if(line(1:6) == 'ATOM') then
                read(line, '(6x, i5, 2x, a4, 14x, 3(f8.3), 22x, a2)') idx, atomName, x, y, z, elementSymbol
                call initAtom(currentAtom, elementSymbol, (/x, y, z/))
                call addAtom(currentMolecule, currentAtom)
            end if
        end if
    end do

!    print *, currentMolecule

    call random_number(translationVector)
    call translateMolecule(currentMolecule, translationVector)

    call writeXYZ(currentMolecule, outputFile)

    call removeMolecule(currentMolecule)

    close(unit)
end program ReadPDBObject