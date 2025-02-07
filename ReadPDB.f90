program ReadPDBObject
    use atom_obj
    use molecule_obj
    use xyz_functions
    use pdb_functions
    use utilities
    implicit none

    type(molecule) :: currentMolecule, originalMolecule
    type(VdWManager) :: manager
    character(len = 4) :: currentFileChar, numberOfFilesChar
    character(len = 128) :: basename, inputFile, outputDirectory, outputFile
    integer :: currentOutputFile, numberOfFiles, stat

    if (iargc() < 1) then
        print '(a)', 'Please provide a PDB file to open'
        stop 10
    end if

    call getarg(1, inputFile)

    if (iargc() >= 2) then
        call getarg(2, numberOfFilesChar)
        read(numberOfFilesChar, *, iostat = stat) numberOfFiles

        if (stat /= 0) then
            print '(a)', 'Second argument was not of type integer'
            stop 10
        end if

        if (numberOfFiles <= 0) then
            print '(a)', 'Second argument should be a strictly positive integer'
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

    print '(a40, a)', 'Number of output XYZ files to generate: ', adjustl(numberOfFilesChar)

    print '(/, a, 1x, a)', 'Reading PDB file', adjustl(trim(inputFile))

    call readPDB(originalMolecule, inputFile)

    print '(/, a)', 'Reading Van der Waals radii'
    call readVdW(manager, 'VdW_radii.txt')

    do currentOutputFile = 1, numberOfFiles
        print '(/, a, i4.4, a, i4.4)', 'Creating output file ', currentOutputFile, ' of ', numberOfFiles
        write (currentFileChar, '(i4.4)') currentOutputFile

        currentMolecule = originalMolecule

        call printMenu(currentMolecule, manager)
        call checkTopology(currentMolecule, manager)

        if (isValidTopology(currentMolecule)) then
            outputFile = adjustl(trim(outputDirectory)) // '/' // trim(basename) // '_' // trim(currentFileChar) // '.xyz'
        else
            outputFile = adjustl(trim(outputDirectory)) // '/NON_VALIDE_' // trim(basename) // '_' // &
                    trim(currentFileChar) // '.xyz'
        end if

        call writeXYZ(currentMolecule, outputFile)

        call removeMolecule(currentMolecule)
    end do

    call removeVdWManager(manager)
    call removeMolecule(originalMolecule)
end program ReadPDBObject