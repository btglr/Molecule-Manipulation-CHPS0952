module utilities
    use molecule_obj
    use vdw_manager
contains
    function getNameFromPath(path) result(name)
        character(len = 128), intent(in) :: path

        character(len = 1) :: pathSepUnix = '/'
        character(len = 1) :: pathSepWindows = '\'
        character(len = 128) name
        integer :: sepPos, dotPos

        ! With back = .true., search from right to left
        sepPos = scan(path, pathSepUnix, back = .true.)

        if (sepPos == 0) then
            sepPos = scan(path, pathSepWindows, back = .true.)
        end if

        if (sepPos > 0) then
            name = path(sepPos + 1:len(path))
            dotPos = scan(name, '.', back = .false.)

            if (dotPos > 0) then
                name = name(1:dotPos - 1)
            end if
        end if
    end function getNameFromPath

    subroutine printMenu(mo, ma)
        type(molecule), intent(inout) :: mo
        type(VdwManager), intent(in) :: ma

        integer :: secondChoice, transformationChoice
        real :: angle, lowerBound, upperBound
        real, dimension(3) :: translationVector

        transformationChoice = -1

        do while (transformationChoice /= 0)
            print '(/, a40, /)', '=== Choose an action ==='
            print '(a22)', '1. Translate molecule'
            print '(a28)', '2. Rotate molecule globally'
            print '(a30)', '3. Rotate molecule internally'
            print '(a18)', '4. Check topology'
            print '(a8, /)', '0. Stop'

            read '(i1)', transformationChoice

            secondChoice = 0

            if (transformationChoice == 1) then
                print '(/, a40, /)', '=== Choose an action ==='
                print '(a32)', '1. Translate by a random vector'
                print '(a31, /)', '2. Translate by a given vector'

                read '(i1)', secondChoice
            else if (transformationChoice == 2 .OR. transformationChoice == 3) then
                print '(/, a40, /)', '=== Choose an action ==='
                print '(a28)', '1. Rotate by a random angle'
                print '(a27, /)', '2. Rotate by a given angle'

                read '(i1)', secondChoice
            else if (transformationChoice == 4) then
                call checkTopology(mo, ma)
            end if

            if (secondChoice == 2) then
                if (transformationChoice == 1) then
                    print '(/, a20)', 'Translation vector: '
                    read '(3(f8.3))', translationVector
                else if (transformationChoice == 2 .OR. transformationChoice == 3) then
                    lowerBound = 0.0
                    if (transformationChoice == 2) then
                        upperBound = 360.0
                    else
                        upperBound = 10.0
                    end if

                    angle = 0.0

                    do while (.NOT. isValidAngle(lowerBound, upperBound, angle))
                        print '(/, a14, f8.3, a, f8.3, a)', 'Angle > ', lowerBound, ' and < ', upperBound, ': '
                        read '(f8.3)', angle
                    end do
                end if
            else if (secondChoice == 1) then
                if (transformationChoice == 1) then
                    call random_number(translationVector)
                else if (transformationChoice == 2 .OR. transformationChoice == 3) then
                    call random_number(angle)

                    if (transformationChoice == 2) then
                        angle = angle * 360.0
                    end if

                    if (transformationChoice == 3) then
                        angle = angle * 10.0
                    end if
                end if
            end if

            if (transformationChoice == 1) then
                call translateMolecule(mo, translationVector)
            else if (transformationChoice == 2) then
                call rotateMoleculeGlobally(mo, angle)
            else if (transformationChoice == 3) then
                call rotateMoleculeInternally(mo, angle)
            end if
        end do
    end subroutine printMenu

    pure logical elemental function isValidAngle(lowerBound, upperBound, value) result(isValid)
        real, intent(in) :: lowerBound, upperBound, value

        isValid = .FALSE.

        if (value > lowerBound .AND. value < upperBound) then
            isValid = .TRUE.
        end if
    end function isValidAngle
end module utilities