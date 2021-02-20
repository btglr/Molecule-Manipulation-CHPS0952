module vdw_manager
    use vdw_obj

    type VdWManager
        type(VdWRadius), dimension(:), allocatable, private :: radii
        integer :: numberOfRadii
    contains
        procedure :: printRadii
        procedure :: initVdWManager
        procedure :: addVdWRadius
        procedure :: readVdW
        generic :: write(formatted) => printRadii
    end type VdWManager
contains
    subroutine printRadii(manager, unit, iotype, v_list, iostat, iomsg)
        class(VdWManager), intent(in) :: manager
        integer, intent(in) :: unit             ! Internal unit to write to.
        character(*), intent(in) :: iotype      ! LISTDIRECTED or DTxxx
        integer, intent(in) :: v_list(:)        ! parameters from fmt spec.
        integer, intent(out) :: iostat          ! non zero on error, etc.
        character(*), intent(inout) :: iomsg    ! define if iostat non zero.
        write (unit, '(/)', IOSTAT = iostat, IOMSG = iomsg)
        write (unit, *, IOSTAT = iostat, IOMSG = iomsg) manager%radii
    end subroutine printRadii

    subroutine initVdWManager(manager, size)
        class(VdWManager), intent(inout) :: manager
        integer, intent(in) :: size
        integer :: ok

        allocate(manager%radii(size), stat = ok)

        if(ok /= 0) then
            print '(a)', "Error during the array allocation ==> aborting"
            stop 666
        end if

        manager%numberOfRadii = 0
    end subroutine initVdWManager

    subroutine addVdWRadius(manager, radius)
        class(VdWManager), intent(inout) :: manager
        type(VdWRadius), intent(in) :: radius

        manager%numberOfRadii = manager%numberOfRadii + 1
        manager%radii(manager%numberOfRadii) = radius
    end subroutine addVdWRadius

    subroutine readVdW(manager, vdwFile)
        class(VdWManager), intent(inout) :: manager
        character(*), intent(in) :: vdwFile

        type(VdWRadius) :: currentRadius
        character(len = 2) :: element
        character(len = 10) :: line
        integer :: i, end, ok, arraySize, unit

        unit = 20
        open(unit, file = vdwFile, iostat = ok, status = 'old')

        if(ok /= 0) then
            print '(a, 4x, a)', "Error during opening", vdwFile
            stop 20
        end if

        ! Read the number of non-comment lines
        arraySize = 0
        do
            read(unit, *, iostat = end) line
            if (end /= 0) then
                exit
            else
                if (line(1:1) /= '#') then
                    arraySize = arraySize + 1
                end if
            end if
        end do

        rewind(unit)
        call initVdWManager(manager, arraySize)

        do
            read(unit, '(a)', iostat = end) line
            if(end /= 0) then
                exit
            else
                if(line(1:1) /= '#') then
                    read(line, '(a2, f6.4)') element, radius
                    call initVdWRadius(currentRadius, element, radius)
                    call addVdWRadius(manager, currentRadius)
                end if
            end if
        end do

        close(unit)
    end subroutine readVdW

    function getVdWRadius(manager, element) result(radius)
        class(VdWManager), intent(in) :: manager
        character(*), intent(in) :: element

        type(VdWRadius) :: radius
        integer :: i

        do i = 1, manager%numberOfRadii
            if (adjustl(getVdWRadiusElement(manager%radii(i))) == element) then
                radius = manager%radii(i)
                exit
            end if
        end do
    end function getVdWRadius
end module vdw_manager