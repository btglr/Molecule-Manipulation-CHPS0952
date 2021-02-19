program MoleculeManipulation
    use atom_obj
    use molecule_obj
    use xyz_handler
    use molecule_func

    implicit none

    character(len = 128) :: firstXYZ, secondXYZ
    type(molecule) :: firstMolecule, secondMolecule

    if(iargc() /= 2) then
        print '(a)', "Please provide two XYZ files to read"
        stop 10
    end if

    call getarg(1, firstXYZ)
    call getarg(2, secondXYZ)

    firstXYZ = trim(firstXYZ)
    secondXYZ = trim(secondXYZ)

    call readXYZ(firstMolecule, firstXYZ, 12)
    call readXYZ(secondMolecule, secondXYZ, 13)

    call rotateMolecule(firstMolecule, 180.0)

    call writeXYZ(firstMolecule, "test.xyz", "Test")
end program MoleculeManipulation