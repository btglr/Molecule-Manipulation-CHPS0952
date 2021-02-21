module molecule_obj
    use atom_obj
    use vdw_obj
    use vdw_manager

    type molecule
        type(atom), dimension(:), allocatable, private :: atoms
        integer, private :: numberOfAtoms
        real, dimension(3), private :: translationVector
        real, private :: rotationAngle
        logical, private :: validTopology
    contains
        procedure :: displayMolecule
        procedure :: initMolecule
        procedure :: removeMolecule
        procedure :: getNumberOfAtoms
        procedure :: getAtom
        procedure :: setAtom
        procedure :: computeRMSD
        generic :: write(formatted) => displayMolecule
    end type molecule
contains
    subroutine initMolecule(m, size)
        class(molecule), intent(inout) :: m
        integer, intent(in) :: size
        integer :: ok

        allocate(m%atoms(size), stat = ok)

        if(ok /= 0) then
            print '(a)', "Error during the array allocation ==> aborting"
            stop 666
        end if

        m%numberOfAtoms = 0
        m%translationVector = 0.0
        m%rotationAngle = 0.0
        m%validTopology = .TRUE.
    end subroutine initMolecule

    subroutine addAtom(m, at)
        class(molecule), intent(inout) :: m
        type(atom), intent(in) :: at

        m%numberOfAtoms = m%numberOfAtoms + 1
        m%atoms(m%numberOfAtoms) = at
    end subroutine addAtom

    subroutine displayMolecule(m, unit, iotype, v_list, iostat, iomsg)
        class(molecule), intent(in) :: m
        integer, intent(in) :: unit             ! Internal unit to write to.
        character(*), intent(in) :: iotype      ! LISTDIRECTED or DTxxx
        integer, intent(in) :: v_list(:)        ! parameters from fmt spec.
        integer, intent(out) :: iostat          ! non zero on error, etc.
        character(*), intent(inout) :: iomsg    ! define if iostat non zero.
        write (unit, *, IOSTAT = iostat, IOMSG = iomsg) m%atoms
    end subroutine displayMolecule

    subroutine removeMolecule(m)
        class(molecule), intent(inout) :: m
        integer :: ok

        deallocate(m%atoms, stat = ok)

        if(ok /= 0) then
            print '(a)', "Error during the atoms array deallocation"
            stop 666
        end if
    end subroutine removeMolecule

    subroutine translateMolecule(m, translationVector)
        class(molecule), intent(inout) :: m
        real, dimension(:), intent(in) :: translationVector
        integer :: atomIndex

        print *, "Translation vector:", translationVector

        m%translationVector = m%translationVector + translationVector

        do atomIndex = 1, m%numberOfAtoms
            call translateAtom(m%atoms(atomIndex), translationVector)
        end do
    end subroutine translateMolecule

    subroutine furthestAtoms(m, firstAtom, secondAtom)
        class(molecule), intent(in) :: m
        type(atom), intent(inout) :: firstAtom, secondAtom
        real, dimension(3) :: currentCoords, otherCoords
        integer :: numberOfAtoms, atomIndex, otherAtomIndex, firstAtomIndex, secondAtomIndex
        real :: maxDistance, distance

        numberOfAtoms = getNumberOfAtoms(m)
        maxDistance = -1

        do atomIndex = 1, numberOfAtoms
            currentCoords = getCoordinates(getAtom(m, atomIndex))

            do otherAtomIndex = atomIndex + 1, numberOfAtoms
                otherCoords = getCoordinates(getAtom(m, otherAtomIndex))
                distance = norm2(currentCoords - otherCoords)

                if (distance > maxDistance) then
                    maxDistance = distance

                    firstAtom = getAtom(m, atomIndex)
                    secondAtom = getAtom(m, otherAtomIndex)

                    firstAtomIndex = atomIndex
                    secondAtomIndex = otherAtomIndex
                end if
            end do
        end do

        print '(a35,x,i8,x,i8)', "Two furthest atoms: ", firstAtomIndex, secondAtomIndex
        print '(a35,x,f8.3,a)', "Distance between the two atoms: ", maxDistance, " Å"
    end subroutine furthestAtoms

    subroutine rotateMoleculeGlobally(m, angleInDegrees)
        class(molecule), intent(inout) :: m
        real, intent(in) :: angleInDegrees

        type(atom) :: firstAtom, secondAtom, newAtom
        real, dimension(3) :: u, unorm
        real :: angleInRadians, PI
        integer :: atomIndex, i, j
        integer, dimension(3, 3) :: identityMatrix
        real, dimension(3, 3) :: wRodrigues, rotationMatrix

        PI = 4.D0 * DATAN(1.D0)
        theta = angleInDegrees * PI / 180

        call furthestAtoms(m, firstAtom, secondAtom)

        u = firstAtom - secondAtom
        unorm = u / norm2(u)

        forall(i = 1:3, j = 1:3) identityMatrix(i, j) = (i / j) * (j / i)

        ! https://mathworld.wolfram.com/RodriguesRotationFormula.html
        wRodrigues = reshape([0.0, unorm(3), -unorm(2), -unorm(3), 0.0, unorm(1), unorm(2), -unorm(1), 0.0], &
                shape(wRodrigues))

        rotationMatrix = identityMatrix + sin(theta) * wRodrigues + (1.0 - cos(theta)) * &
                matMul(wRodrigues, wRodrigues)

        m%rotationAngle = m%rotationAngle + angleInDegrees

        do atomIndex = 1, getNumberOfAtoms(m)
            call rotateAtom(m%atoms(atomIndex), rotationMatrix, getCoordinates(secondAtom))
        end do

        print '(a35,x,3(f8.3))', "Axis vector: ", u
        print '(a35,x,3(f8.3))', "Normalized Axis vector: ", unorm
    end subroutine rotateMoleculeGlobally

    subroutine rotateMoleculeInternally(m, angleInDegrees)
        class(molecule), intent(inout) :: m
        real, intent(in) :: angleInDegrees

        integer :: atomIndex, carbonIndex, numberOfCarbonBonds, selectedBond
        integer, dimension(3, 3) :: identityMatrix
        integer, dimension(:, :), allocatable :: carbonBonds
        real :: theta, PI, random
        real, dimension(3) :: u, unorm
        real, dimension(3, 3) :: wRodrigues, rotationMatrix
        type(atom) :: firstAtom, secondAtom

        if (angleInDegrees > 10.0) then
            print '(a)', 'Internal rotation should not be greater than 10 degrees'
            stop 10
        end if

        PI = 4.D0 * DATAN(1.D0)
        theta = angleInDegrees * PI / 180

        call getCarbonBonds(m, carbonBonds, numberOfCarbonBonds)
        call random_number(random)

        ! Select a random bond between the first and last - 1
        selectedBond = 1 + floor((numberOfCarbonBonds - 1) * random)

        print '(a35, x, i8)', 'Selected bond: ', selectedBond
        print '(a35, x, 2(i8))', 'Atoms of selected bond: ', carbonBonds(selectedBond, :)

        firstAtom = getAtom(m, carbonBonds(selectedBond, 1))
        secondAtom = getAtom(m, carbonBonds(selectedBond, 2))

        u = firstAtom - secondAtom
        unorm = u / norm2(u)

        forall(i = 1:3, j = 1:3) identityMatrix(i, j) = (i / j) * (j / i)

        ! https://mathworld.wolfram.com/RodriguesRotationFormula.html
        wRodrigues = reshape([0.0, unorm(3), -unorm(2), -unorm(3), 0.0, unorm(1), unorm(2), -unorm(1), 0.0], &
                shape(wRodrigues))

        rotationMatrix = identityMatrix + sin(theta) * wRodrigues + (1.0 - cos(theta)) * &
                matMul(wRodrigues, wRodrigues)

        do atomIndex = carbonBonds(selectedBond, 2), m%numberOfAtoms
            call rotateAtom(m%atoms(atomIndex), rotationMatrix, getCoordinates(secondAtom))
        end do

        print '(a35,x,3(f8.3))', 'Axis vector: ', u
        print '(a35,x,3(f8.3))', 'Normalized Axis vector: ', unorm
    end subroutine rotateMoleculeInternally

    real function computeRMSD(m1, m2, type) result(rmsd)
        class(molecule), intent(in) :: m1
        class(molecule), intent(in) :: m2
        character(*), intent(in) :: type

        real, dimension(:, :), allocatable :: difference
        integer :: numberOfAtoms, numberOfHeavyAtoms, atomIndex
        logical :: sameMolecule

        if (getNumberOfAtoms(m1) /= getNumberOfAtoms(m2)) then
            print '(a)', "Error: the molecules do not have the same number of atoms ==> aborting"
            stop 20
        end if

        numberOfAtoms = getNumberOfAtoms(m1)

        sameMolecule = .TRUE.

        do atomIndex = 1, numberOfAtoms
            if (getAtomName(getAtom(m1, atomIndex)) /= getAtomName(getAtom(m2, atomIndex))) then
                print *, getAtomName(getAtom(m1, atomIndex)), getAtomName(getAtom(m2, atomIndex))
                sameMolecule = .FALSE.
                exit
            end if
        end do

        if (sameMolecule.eqv..FALSE.) then
            print '(a)', "Error: the molecules are not identical ==> aborting"
            stop 20
        end if

        if (type == "standard") then
            print '(a25,i8)', "Number of atoms:", numberOfAtoms
            allocate(difference(numberOfAtoms, 3))

            do atomIndex = 1, numberOfAtoms
                difference(atomIndex, :) = getCoordinates(getAtom(m1, atomIndex)) - getCoordinates(getAtom(m2, atomIndex))
            end do
        else if (type == "heavy") then
            numberOfHeavyAtoms = 0

            do atomIndex = 1, numberOfAtoms
                if (getAtomName(getAtom(m1, atomIndex)) /= "H") then
                    numberOfHeavyAtoms = numberOfHeavyAtoms + 1
                end if
            end do

            allocate(difference(numberOfHeavyAtoms, 3))

            do atomIndex = 1, numberOfHeavyAtoms
                difference(atomIndex, :) = getCoordinates(getAtom(m1, atomIndex)) - getCoordinates(getAtom(m2, atomIndex))
            end do

            numberOfAtoms = numberOfHeavyAtoms
            print '(a25,i8)', "Number of heavy atoms:", numberOfAtoms
        else
            stop 20
        end if

        print '(a25,f8.0)', "Sum:", sum(difference ** 2)

        rmsd = sqrt(sum(difference ** 2) / numberOfAtoms)
    end function computeRMSD

    integer function getNumberOfAtoms(m) result(numberOfAtoms)
        class(molecule), intent(in) :: m

        numberOfAtoms = m%numberOfAtoms
    end function getNumberOfAtoms

    type(atom) function getAtom(m, atomIndex) result(at)
        class(molecule), intent(in) :: m
        integer, intent(in) :: atomIndex

        at = m%atoms(atomIndex)
    end function getAtom

    subroutine setAtom(m, atomIndex, at)
        class(molecule), intent(inout) :: m
        integer, intent(in) :: atomIndex
        type(atom), intent(in) :: at

        m%atoms(atomIndex) = at
    end subroutine setAtom

    function getTranslationVector(m) result(translationVector)
        class(molecule), intent(in) :: m

        real, dimension(3) :: translationVector

        translationVector = m%translationVector
    end function getTranslationVector

    subroutine setTranslationVector(m, translationVector)
        class(molecule), intent(inout) :: m
        real, dimension(3), intent(in) :: translationVector

        m%translationVector = translationVector
    end subroutine setTranslationVector

    type(real) function getRotationAngle(m) result(rotationVector)
        class(molecule), intent(in) :: m

        rotationAngle = m%rotationAngle
    end function getRotationAngle

    subroutine setRotationAngle(m, rotationAngle)
        class(molecule), intent(inout) :: m
        real, intent(in) :: rotationAngle

        m%rotationAngle = rotationAngle
    end subroutine setRotationAngle

    function filterByElement(m, element) result(atoms)
        class(molecule), intent(in) :: m
        character(*), intent(in) :: element

        integer :: atomsOfElement, atomIndex, ok, atomOfElementIndex
        type(atom), dimension(:), allocatable :: atoms

        atomsOfElement = 0

        do atomIndex = 1, m%numberOfAtoms
            if (adjustl(getAtomName(getAtom(m, atomIndex))) == adjustl(element)) then
                atomsOfElement = atomsOfElement + 1
            end if
        end do

        allocate(atoms(atomsOfElement), stat = ok)

        if(ok /= 0) then
            print '(a)', "Error during the array allocation ==> aborting"
            stop 666
        end if

        atomOfElementIndex = 1
        do atomIndex = 1, m%numberOfAtoms
            if (adjustl(getAtomName(getAtom(m, atomIndex))) == adjustl(element)) then
                atoms(atomOfElementIndex) = getAtom(m, atomIndex)
                atomOfElementIndex = atomOfElementIndex + 1
            end if
        end do
    end function filterByElement

    subroutine getCarbonBonds(m, carbonBonds, numberOfCarbonBonds)
        class(molecule), intent(in) :: m
        integer, dimension(:, :), allocatable, intent(inout) :: carbonBonds
        integer, intent(inout) :: numberOfCarbonBonds

        integer :: atomIndex, carbonIndex

        numberOfCarbonBonds = 0

        do atomIndex = 1, m%numberOfAtoms - 1
            if (getAtomName(m%atoms(atomIndex)) == getAtomName(m%atoms(atomIndex + 1))) then
                numberOfCarbonBonds = numberOfCarbonBonds + 1
            end if
        end do

        print '(a35, x, i8)', 'Number of carbon bonds: ', numberOfCarbonBonds

        allocate(carbonBonds(numberOfCarbonBonds, 2))

        carbonIndex = 1
        do atomIndex = 1, m%numberOfAtoms - 1
            if (getAtomName(m%atoms(atomIndex)) == getAtomName(m%atoms(atomIndex + 1))) then
                carbonBonds(carbonIndex, 1) = atomIndex
                carbonBonds(carbonIndex, 2) = atomIndex + 1

                carbonIndex = carbonIndex + 1
            end if
        end do
    end subroutine getCarbonBonds

    subroutine checkTopology(m, manager)
        class(molecule), intent(inout) :: m
        type(VdWManager), intent(in) :: manager

        integer :: atomIndex, otherAtomIndex
        real :: sumRadii, distance
        type(VdWRadius) :: radius

        do atomIndex = 1, m%numberOfAtoms
            do otherAtomIndex = 1, m%numberOfAtoms
                distance = norm2(getCoordinates(m%atoms(atomIndex)) - getCoordinates(m%atoms(otherAtomIndex)))
                sumRadii = getVdWRadius(getVdWObject(manager, getAtomName(m%atoms(atomIndex)))) + &
                        getVdWRadius(getVdWObject(manager, getAtomName(m%atoms(otherAtomIndex))))

                if (distance > 0 .AND. distance / sumRadii < 0.35) then
                    m%validTopology = .FALSE.

                    print '(a35,x,i8)', 'First atom: ', atomIndex
                    print *, m%atoms(atomIndex)
                    print '(a35,x,i8)', 'Second atom: ', otherAtomIndex
                    print *, m%atoms(otherAtomIndex)
                    print '(a35,x,f8.3)', 'Distance: ', distance
                    print '(a35,x,f8.3)', 'Sum of radii: ', sumRadii
                    print '(a35,x,f8.3)', 'Ratio: ', distance / sumRadii
                    print '(a35)', 'Invalid topology'
                    exit
                end if
            end do

            if (m%validTopology .eqv. .FALSE.) then
                exit
            end if
        end do
    end subroutine checkTopology

end module molecule_obj