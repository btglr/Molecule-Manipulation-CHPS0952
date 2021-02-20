program ReadPDBObject
    use atom_obj
    use molecule_obj
    use xyz_functions
    use pdb_functions
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

    call readPDB(currentMolecule, inputFile)
!    print *, currentMolecule

    call random_number(translationVector)
    call translateMolecule(currentMolecule, translationVector)
    call rotateMolecule(currentMolecule, 180.0)

    call writeXYZ(currentMolecule, outputFile)

    call removeMolecule(currentMolecule)

    close(unit)
end program ReadPDBObject