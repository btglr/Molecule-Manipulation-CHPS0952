program RMSD
    use atom_obj
    use molecule_obj
    use xyz_handler

    implicit none

    character(len = 128) :: firstXYZ, secondXYZ
    type(molecule) :: firstMolecule, secondMolecule
    real :: rmsdValue

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

    rmsdValue = computeRMSD(firstMolecule, secondMolecule, type = "standard")

    print '(a25,f8.3)', "RMSD Standard:", rmsdValue

    rmsdValue = computeRMSD(firstMolecule, secondMolecule, type = "heavy")

    print '(a25,f8.3)', "RMSD Heavy:", rmsdValue
end program RMSD