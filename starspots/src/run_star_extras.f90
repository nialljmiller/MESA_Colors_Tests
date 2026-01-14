! x_ctrl(1) = fspot (spot filling factor)
! x_ctrl(2) = xspot (temperature contrast T_spot/T_photosphere)

module run_star_extras

    use star_lib
    use star_def
    use const_def
    use math_lib

    implicit none

contains

    subroutine extras_controls(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        ! Set the starspot gradr factor routine
        s% other_gradr_factor => starspot_gradr_factor

    end subroutine extras_controls


    subroutine starspot_gradr_factor(id, ierr)
        ! Modifies radiative gradient to account for starspot coverage
        ! Following Somers & Pinsonneault 2015, Somers et al. 2020
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type (star_info), pointer :: s
        real(dp) :: fspot, xspot, xspot_local, factor
        integer :: k

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        ! Get spot parameters from x_ctrl
        fspot = s% x_ctrl(1)  ! Filling factor
        xspot = s% x_ctrl(2)  ! Temperature contrast at surface

        ! Apply gradr modification throughout the star
        do k = 1, s% nz
            ! Depth-dependent xspot following MESA VI eq. 44
            ! xspot(r) = 1 - (1 - xspot) * T_ambient(r) / T(r)
            if (s% T(k) > 0d0) then
                xspot_local = 1d0 - (1d0 - xspot) * s% T(s% nz) / s% T(k)
            else
                xspot_local = xspot
            end if

            ! gradr factor = 1 / [fspot * xspot^4 + (1 - fspot)]
            factor = fspot * pow4(xspot_local) + (1d0 - fspot)
            if (factor > 0d0) then
                s% gradr_factor(k) = 1d0 / factor
            else
                s% gradr_factor(k) = 1d0
            end if
        end do

    end subroutine starspot_gradr_factor


    subroutine extras_startup(id, restart, ierr)
        integer, intent(in) :: id
        logical, intent(in) :: restart
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

    end subroutine extras_startup


    integer function extras_start_step(id)
        integer, intent(in) :: id
        integer :: ierr
        type (star_info), pointer :: s

        extras_start_step = 0
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

    end function extras_start_step


    integer function extras_check_model(id)
        integer, intent(in) :: id
        integer :: ierr
        type (star_info), pointer :: s

        extras_check_model = keep_going
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

    end function extras_check_model


    integer function how_many_extra_history_columns(id)
        integer, intent(in) :: id

        how_many_extra_history_columns = 2

    end function how_many_extra_history_columns


    subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
        integer, intent(in) :: id, n
        character(len=maxlen_history_column_name) :: names(n)
        real(dp) :: vals(n)
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        ! Record spot parameters in history
        names(1) = 'fspot'
        vals(1) = s% x_ctrl(1)

        names(2) = 'xspot'
        vals(2) = s% x_ctrl(2)

    end subroutine data_for_extra_history_columns


    integer function how_many_extra_profile_columns(id)
        integer, intent(in) :: id

        how_many_extra_profile_columns = 0

    end function how_many_extra_profile_columns


    subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
        integer, intent(in) :: id, n, nz
        character(len=maxlen_profile_column_name) :: names(n)
        real(dp) :: vals(nz, n)
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

    end subroutine data_for_extra_profile_columns


    integer function extras_finish_step(id)
        integer, intent(in) :: id
        integer :: ierr
        type (star_info), pointer :: s

        extras_finish_step = keep_going
        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

    end function extras_finish_step


    subroutine extras_after_evolve(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

    end subroutine extras_after_evolve

end module run_star_extras
