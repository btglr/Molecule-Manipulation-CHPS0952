program VdWTest
    use vdw_obj
    use vdw_manager

    implicit none

    type(VdWManager) :: manager

    call readVdW(manager, 'VdW_radii.txt')
    print *, manager
end program VdWTest