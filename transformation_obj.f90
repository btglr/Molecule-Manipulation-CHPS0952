module transformation_obj
    type transformation
        integer, private :: transformationType
        integer, dimension(2), private :: atomsIndices
        real, private :: angle
        real, dimension(3), private :: translationVector
    contains
        procedure :: printTransformation
        procedure :: initGlobalRotation
        procedure :: initInternalRotation
        procedure :: initTranslation
        generic :: write(formatted) => printTransformation
    end type transformation
contains
    subroutine printTransformation(t, unit, iotype, v_list, iostat, iomsg)
        class(transformation), intent(in) :: t
        integer, intent(in) :: unit             ! Internal unit to write to.
        character(*), intent(in) :: iotype      ! LISTDIRECTED or DTxxx
        integer, intent(in) :: v_list(:)        ! parameters from fmt spec.
        integer, intent(out) :: iostat          ! non zero on error, etc.
        character(*), intent(inout) :: iomsg    ! define if iostat non zero.

        if (t%transformationType == 0) then
            write (unit, 1000, IOSTAT = iostat, IOMSG = iomsg) "Translation", t%translationVector
            1000 format(a20, 1x, 3(f8.3, 1x))
        else
            if (t%transformationType == 1) then
                write (unit, 1001, IOSTAT = iostat, IOMSG = iomsg) "Global rotation", "Atoms: ", t%atomsIndices, &
                        "Angle: ", t%angle
            else
                write (unit, 1001, IOSTAT = iostat, IOMSG = iomsg) "Internal rotation", "Atoms: ", t%atomsIndices, &
                        "Angle: ", t%angle
            end if

            1001 format(a20, 1x, a10, 2(i5, 1x), a10, f8.3)
        end if
    end subroutine printTransformation

    subroutine initGlobalRotation(t, atomsIndices, angle)
        class(transformation), intent(inout) :: t
        integer, dimension(2), intent(in) :: atomsIndices
        real, intent(in) :: angle

        t%transformationType = 1
        t%atomsIndices = atomsIndices
        t%angle = angle
    end subroutine initGlobalRotation

    subroutine initInternalRotation(t, atomsIndices, angle)
        class(transformation), intent(inout) :: t
        integer, dimension(2), intent(in) :: atomsIndices
        real, intent(in) :: angle

        t%transformationType = 2
        t%atomsIndices = atomsIndices
        t%angle = angle
    end subroutine initInternalRotation

    subroutine initTranslation(t, translationVector)
        class(transformation), intent(inout) :: t
        real, dimension(3), intent(in) :: translationVector

        t%transformationType = 0
        t%translationVector = translationVector
    end subroutine initTranslation
end module transformation_obj