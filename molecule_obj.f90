module molecule_obj
    use atom_obj
    use vdw_obj
    use vdw_manager
    use transformation_obj

    type molecule
        type(atom), dimension(:), allocatable, private :: atoms
        type(transformation), dimension(:), allocatable, private :: transformations, tmpTransformations
        integer, private :: numberOfAtoms, numberOfTransformations
        logical, private :: validTopology
        real, private :: globalRotationAngle, internalRotationAngle
        real, dimension(3), private :: translationVector
    contains
        procedure :: initMolecule
        procedure :: addAtom
        procedure :: displayMolecule
        procedure :: removeMolecule
        procedure :: translateMolecule
        procedure :: furthestAtoms
        procedure :: rotateMoleculeGlobally
        procedure :: rotateMoleculeInternally
        procedure :: computeRMSD
        procedure :: getNumberOfAtoms
        procedure :: getAtom
        procedure :: setAtom
        procedure :: filterByElement
        procedure :: getCarbonBonds
        procedure :: checkTopology
        procedure :: isValidTopology
        procedure :: getTransformations
        generic :: write(formatted) => displayMolecule
    end type molecule
contains
    subroutine initMolecule(m, size)
        class(molecule), intent(inout) :: m
        integer, intent(in) :: size

        integer :: ok

        allocate(m%atoms(size), stat = ok)

        if(ok /= 0) then
            print '(a)', 'Error during the array allocation ==> aborting'
            stop 666
        end if

        m%numberOfAtoms = 0
        m%numberOfTransformations = 0
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
            print '(a)', 'Error during the atoms array deallocation'
            stop 666
        end if

        if (m%numberOfTransformations > 0) then
            deallocate(m%transformations, stat = ok)

            if(ok /= 0) then
                print '(a)', 'Error during the molecule transformations deallocation'
                stop 666
            end if
        end if
    end subroutine removeMolecule

    subroutine translateMolecule(m, translationVector)
        class(molecule), intent(inout) :: m
        real, dimension(:), intent(in) :: translationVector

        type(transformation) :: t
        integer :: atomIndex

        print '(/, a40)', '=== Translating molecule ==='
        print '(/, a40, 3(f8.3))', 'Translation vector: ', translationVector

        do atomIndex = 1, m%numberOfAtoms
            call translateAtom(m%atoms(atomIndex), translationVector)
        end do

        call initTranslation(t, translationVector)
        call addTransformation(m, t)
    end subroutine translateMolecule

    subroutine furthestAtoms(m, firstAtomIndex, secondAtomIndex)
        class(molecule), intent(in) :: m
        integer, intent(inout) :: firstAtomIndex, secondAtomIndex

        integer :: atomIndex, numberOfAtoms, otherAtomIndex
        real :: distance, maxDistance
        real, dimension(3) :: currentCoords, otherCoords

        numberOfAtoms = getNumberOfAtoms(m)
        maxDistance = -1

        do atomIndex = 1, numberOfAtoms
            currentCoords = getCoordinates(getAtom(m, atomIndex))

            do otherAtomIndex = atomIndex + 1, numberOfAtoms
                otherCoords = getCoordinates(getAtom(m, otherAtomIndex))
                distance = norm2(currentCoords - otherCoords)

                if (distance > maxDistance) then
                    maxDistance = distance

                    firstAtomIndex = atomIndex
                    secondAtomIndex = otherAtomIndex
                end if
            end do
        end do

        print '(a40, i8, 1x, i8)', 'Two furthest atoms: ', firstAtomIndex, secondAtomIndex
        print '(a40, f8.3, a)', 'Distance between the two atoms: ', maxDistance, ' Å'
    end subroutine furthestAtoms

    subroutine rotateMoleculeGlobally(m, angleInDegrees)
        class(molecule), intent(inout) :: m
        real, intent(in) :: angleInDegrees

        type(atom) :: firstAtom, secondAtom
        type(transformation) :: t
        integer :: atomIndex, firstAtomIndex, i, j, secondAtomIndex
        integer, dimension(3, 3) :: identityMatrix
        real :: PI, theta
        real, dimension(3) :: u, unorm
        real, dimension(3, 3) :: rotationMatrix, wRodrigues

        print '(/, a40)', '=== Rotating molecule globally ==='
        print '(/, a40, f7.3)', 'Angle of rotation: ', angleInDegrees

        PI = 4.D0 * datan(1.D0)
        theta = angleInDegrees * PI / 180

        call furthestAtoms(m, firstAtomIndex, secondAtomIndex)

        firstAtom = getAtom(m, firstAtomIndex)
        secondAtom = getAtom(m, secondAtomIndex)

        u = firstAtom - secondAtom
        unorm = u / norm2(u)

        forall(i = 1:3, j = 1:3) identityMatrix(i, j) = (i / j) * (j / i)

        ! https://mathworld.wolfram.com/RodriguesRotationFormula.html
        wRodrigues = reshape([0.0, unorm(3), -unorm(2), -unorm(3), 0.0, unorm(1), unorm(2), -unorm(1), 0.0], &
                shape(wRodrigues))

        rotationMatrix = identityMatrix + sin(theta) * wRodrigues + (1.0 - cos(theta)) * &
                matMul(wRodrigues, wRodrigues)

        do atomIndex = 1, getNumberOfAtoms(m)
            call rotateAtom(m%atoms(atomIndex), rotationMatrix, getCoordinates(secondAtom))
        end do

        print '(a40, 1x, 3(f8.3))', 'Axis vector: ', u
        print '(a40, 1x, 3(f8.3))', 'Normalized Axis vector: ', unorm

        call initGlobalRotation(t, [firstAtomIndex, secondAtomIndex], angleInDegrees)
        call addTransformation(m, t)
    end subroutine rotateMoleculeGlobally

    subroutine rotateMoleculeInternally(m, angleInDegrees)
        class(molecule), intent(inout) :: m
        real, intent(in) :: angleInDegrees

        type(atom) :: firstAtom, secondAtom
        type(transformation) :: t
        integer :: atomIndex, i, j, numberOfCarbonBonds, selectedBond
        integer, dimension(3, 3) :: identityMatrix
        integer, dimension(:, :), allocatable :: carbonBonds
        real :: PI, random, theta
        real, dimension(3) :: u, unorm
        real, dimension(3, 3) :: rotationMatrix, wRodrigues

        if (angleInDegrees > 10.0) then
            print '(a)', 'Internal rotation should not be greater than 10 degrees'
            stop 10
        end if

        print '(/, a40)', '=== Rotating molecule internally ==='
        print '(/, a40, f7.3)', 'Angle of rotation: ', angleInDegrees

        PI = 4.D0 * datan(1.D0)
        theta = angleInDegrees * PI / 180

        call getCarbonBonds(m, carbonBonds, numberOfCarbonBonds)
        call random_number(random)

        ! Select a random bond between the first and last - 1
        selectedBond = 1 + floor((numberOfCarbonBonds - 1) * random)

        print '(a40, 1x, i8)', 'Selected bond: ', selectedBond
        print '(a40, 1x, 2(i8))', 'Atoms of selected bond: ', carbonBonds(selectedBond, :)

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

        print '(a40, 1x, 3(f8.3))', 'Axis vector: ', u
        print '(a40, 1x, 3(f8.3))', 'Normalized Axis vector: ', unorm

        call initInternalRotation(t, [carbonBonds(selectedBond, 1), carbonBonds(selectedBond, 2)], angleInDegrees)
        call addTransformation(m, t)

        deallocate(carbonBonds)
    end subroutine rotateMoleculeInternally

    real function computeRMSD(m1, m2, type) result(rmsd)
        class(molecule), intent(in) :: m1
        class(molecule), intent(in) :: m2
        character(*), intent(in) :: type

        integer :: atomIndex, numberOfAtoms, numberOfHeavyAtoms
        logical :: sameMolecule
        real, dimension(:, :), allocatable :: difference

        if (getNumberOfAtoms(m1) /= getNumberOfAtoms(m2)) then
            print '(a)', 'Error: the molecules do not have the same number of atoms ==> aborting'
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
            print '(a)', 'Error: the molecules are not identical ==> aborting'
            stop 20
        end if

        if (type == 'standard') then
            print '(a40, i12)', 'Number of atoms:', numberOfAtoms
            allocate(difference(numberOfAtoms, 3))

            do atomIndex = 1, numberOfAtoms
                difference(atomIndex, :) = getCoordinates(getAtom(m1, atomIndex)) - getCoordinates(getAtom(m2, atomIndex))
            end do
        else if (type == 'heavy') then
            numberOfHeavyAtoms = 0

            do atomIndex = 1, numberOfAtoms
                if (getAtomName(getAtom(m1, atomIndex)) /= 'H') then
                    numberOfHeavyAtoms = numberOfHeavyAtoms + 1
                end if
            end do

            allocate(difference(numberOfHeavyAtoms, 3))

            do atomIndex = 1, numberOfHeavyAtoms
                difference(atomIndex, :) = getCoordinates(getAtom(m1, atomIndex)) - getCoordinates(getAtom(m2, atomIndex))
            end do

            numberOfAtoms = numberOfHeavyAtoms
            print '(a40, i12)', 'Number of heavy atoms:', numberOfAtoms
        else
            stop 20
        end if

        print '(a40, f12.0)', 'Sum:', sum(difference ** 2)

        rmsd = sqrt(sum(difference ** 2) / numberOfAtoms)

        deallocate(difference)
    end function computeRMSD

    pure elemental integer function getNumberOfAtoms(m) result(numberOfAtoms)
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

    function filterByElement(m, element) result(atoms)
        class(molecule), intent(in) :: m
        character(*), intent(in) :: element

        type(atom), dimension(:), allocatable :: atoms
        integer :: atomIndex, atomsOfElement, atomOfElementIndex, ok

        atomsOfElement = 0

        do atomIndex = 1, m%numberOfAtoms
            if (adjustl(getAtomName(getAtom(m, atomIndex))) == adjustl(element)) then
                atomsOfElement = atomsOfElement + 1
            end if
        end do

        allocate(atoms(atomsOfElement), stat = ok)

        if(ok /= 0) then
            print '(a)', 'Error during the array allocation ==> aborting'
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

        print '(a40, 1x, i8)', 'Number of carbon bonds: ', numberOfCarbonBonds

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

        integer :: firstAtomOfPair, secondAtomOfPair
        real :: distance, sumRadii

        print '(/, a40, /)', '=== Checking topology of molecule ==='

        do firstAtomOfPair = 1, m%numberOfAtoms
            do secondAtomOfPair = firstAtomOfPair + 1, m%numberOfAtoms
                distance = norm2(getCoordinates(m%atoms(firstAtomOfPair)) - getCoordinates(m%atoms(secondAtomOfPair)))
                sumRadii = getVdWRadius(getVdWObject(manager, getAtomName(m%atoms(firstAtomOfPair)))) + &
                        getVdWRadius(getVdWObject(manager, getAtomName(m%atoms(secondAtomOfPair))))

                if (distance > 0 .AND. distance / sumRadii < 0.35) then
                    m%validTopology = .FALSE.

                    print '(a50)', '/!\ INVALID TOPOLOGY /!\'
                    print '(a40, 1x, i8)', 'First atom: ', firstAtomOfPair
                    print *, m%atoms(firstAtomOfPair)
                    print '(a40, 1x, i8)', 'Second atom: ', secondAtomOfPair
                    print *, m%atoms(secondAtomOfPair)
                    print '(a40, 1x, f8.3, a)', 'Distance: ', distance, ' Å'
                    print '(a40, 1x, f8.3, a)', 'Sum of radii: ', sumRadii, ' Å'
                    print '(a40, 1x, f8.3, a)', 'Ratio: ', distance / sumRadii, ' < 0.35'
                    exit
                end if
            end do

            if (m%validTopology .eqv. .FALSE.) then
                exit
            end if
        end do
    end subroutine checkTopology

    pure elemental logical function isValidTopology(m) result(validTopology)
        class(molecule), intent(in) :: m

        validTopology = m%validTopology
    end function isValidTopology

    pure function getTransformations(m) result(ts)
        class(molecule), intent(in) :: m

        type(transformation), dimension(size(m%transformations)) :: ts

        ts = m%transformations
    end function getTransformations

    pure function getTransformation(m, index) result(t)
        class(molecule), intent(in) :: m
        integer, intent(in) :: index

        type(transformation) :: t

        t = m%transformations(index)
    end function getTransformation

    subroutine addTransformation(m, t)
        class(molecule), intent(inout) :: m
        type(transformation), intent(in) :: t

        type(transformation), dimension(:), allocatable :: tmpTransformations

        if (m%numberOfTransformations > 0) then
            allocate(tmpTransformations(m%numberOfTransformations + 1))
            tmpTransformations(1:m%numberOfTransformations) = m%transformations(1:m%numberOfTransformations)
            call move_alloc(tmpTransformations, m%transformations)
        else
            allocate(m%transformations(m%numberOfTransformations + 1))
        end if

        m%numberOfTransformations = m%numberOfTransformations + 1
        m%transformations(m%numberOfTransformations) = t
    end subroutine addTransformation

end module molecule_obj