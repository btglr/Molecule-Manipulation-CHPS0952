program VdWTest
    use vdw_obj
    use vdw_manager

    implicit none

    type(VdWManager) :: manager

    call readVdW(manager, 'VdW_radii.txt')

    print '(a)', 'Printing all radii'
    print *, manager

    print '(a, 1x, a2)', 'Getting one radius from element', 'Ba'
    print *, getVdWObject(manager, 'Ba')
end program VdWTest