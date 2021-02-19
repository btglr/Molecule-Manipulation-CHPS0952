module molecule_func
    use molecule_obj

contains

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

    subroutine rotateMolecule(m, angleInDegrees)
        class(molecule), intent(inout) :: m
        real, intent(in) :: angleInDegrees

        type(atom) :: firstAtom, secondAtom, newAtom
        real, dimension(3) :: u, unorm
        real :: angleInRadians, PI
        integer :: atomIndex

        PI = 4.D0 * DATAN(1.D0)
        angleInRadians = angleInDegrees * PI / 180

        call furthestAtoms(m, firstAtom, secondAtom)

        u = firstAtom - secondAtom
        unorm = u / norm2(u)

        m%rotationVector = m%rotationVector + u

        do atomIndex = 1, getNumberOfAtoms(m)
            newAtom = getAtom(m, atomIndex)
            call rotateAtom(newAtom, u, unorm, angleInRadians, getCoordinates(secondAtom))
            call setAtom(m, atomIndex, newAtom)
        end do

        print '(a35,x,3(f8.3))', "Axis vector: ", u
        print '(a35,x,3(f8.3))', "Normalized Axis vector: ", unorm
        print *, firstAtom
        print *, secondAtom
    end subroutine rotateMolecule

end module molecule_func