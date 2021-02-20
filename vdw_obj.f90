module vdw_obj
    type VdWRadius
        character(:), allocatable, private :: element
        real, private :: radius
    contains
        procedure :: printRadius
        procedure :: initVdWRadius
        generic :: write(formatted) => printRadius
    end type VdWRadius
contains
    subroutine initVdWRadius(vr, element, radius)
        class(VdWRadius), intent(inout) :: vr
        character(len = *), intent(in) :: element
        real, intent(in) :: radius

        vr%element = element
        vr%radius = radius
    end subroutine initVdWRadius

    subroutine printRadius(vr, unit, iotype, v_list, iostat, iomsg)
        class(VdWRadius), intent(in) :: vr
        integer, intent(in) :: unit             ! Internal unit to write to.
        character(*), intent(in) :: iotype      ! LISTDIRECTED or DTxxx
        integer, intent(in) :: v_list(:)        ! parameters from fmt spec.
        integer, intent(out) :: iostat          ! non zero on error, etc.
        character(*), intent(inout) :: iomsg    ! define if iostat non zero.
        write (unit, 999, IOSTAT = iostat, IOMSG = iomsg) vr%element, vr%radius
        999 format(a2, x, f6.4, /)
    end subroutine printRadius
end module vdw_obj