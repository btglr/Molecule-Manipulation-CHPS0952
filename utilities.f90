module utilities
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
end module utilities