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

end module molecule_func