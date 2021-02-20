program ReadPDBObject
    use atom_obj
    use molecule_obj
    use xyz_functions
    use pdb_functions
    use utilities
    implicit none

    integer :: i, end, ok, idx, arraySize, unit, numberOfFiles, currentOutputFile
    character(len = 128) :: inputFile, outputFile, basename, line, newFileName, outputDirectory
    character(len = 4) :: atomName, elementSymbol, numberOfFilesChar, currentFileChar
    real :: x, y, z, angleInDegrees
    type(atom) :: currentAtom
    type(atom), dimension(:), allocatable :: atoms
    type(molecule) :: currentMolecule, originalMolecule
    real, dimension(3) :: translationVector
    integer :: stat

    if (iargc() < 1) then
        print '(a)', "Please provide a PDB file to open"
        stop 10
    end if

    call getarg(1, inputFile)

    if (iargc() >= 2) then
        call getarg(2, numberOfFilesChar)
        read(numberOfFilesChar, *, iostat = stat) numberOfFiles

        if (stat /= 0) then
            print '(a)', "Second argument was not of type integer"
            stop 10
        end if
    else
        numberOfFiles = 1
    end if

    outputDirectory = 'XYZ'
    call execute_command_line('mkdir -p ' // adjustl(trim(outputDirectory)))

    ! Write the number of files back in the string for later use instead of using numberOfFilesChar directly from
    ! getarg in case the provided number isn't correctly formatted (0050, ...)
    write (numberOfFilesChar, '(i4)') numberOfFiles

    inputFile = trim(inputFile)
    basename = trim(getNameFromPath(inputFile))

    print '(a, a)', "File to read: ", inputFile
    print '(a, a)', "Number of output XYZ files to generate: ", adjustl(numberOfFilesChar)

    call readPDB(originalMolecule, inputFile)

    do currentOutputFile = 1, numberOfFiles
        write (currentFileChar, '(i4.4)') currentOutputFile
        outputFile = adjustl(trim(outputDirectory)) // '/' // trim(basename) // "_" // trim(currentFileChar) // ".xyz"

        currentMolecule = originalMolecule

        call random_number(translationVector)
        call random_number(angleInDegrees)

        call translateMolecule(currentMolecule, translationVector)
        call rotateMolecule(currentMolecule, angleInDegrees * 360.0)
        call writeXYZ(currentMolecule, outputFile)
    end do

    call removeMolecule(originalMolecule)
end program ReadPDBObject