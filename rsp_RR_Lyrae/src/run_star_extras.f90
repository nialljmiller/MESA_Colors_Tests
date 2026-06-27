! ***********************************************************************
!
!   Copyright (C) 2018-2019 The MESA Team
!
! ***********************************************************************

module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff

      implicit none

      include "test_suite_extras_def.inc"

      logical :: need_to_write_LINA_data

      integer, parameter :: mode_settle = 1
      integer, parameter :: mode_colors = 2

      integer :: last_checked_period = -1
      integer :: stable_period_count = 0
      integer :: colors_start_period = -1

      contains

      include "test_suite_extras.inc"


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_start_step => extras_start_step
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns
      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call test_suite_startup(s, restart, ierr)
         if (ierr /= 0) return

         last_checked_period = -1
         stable_period_count = 0
         colors_start_period = -1

         if (.not. restart) then
            need_to_write_LINA_data = len_trim(s% x_character_ctrl(10)) > 0
         else
            need_to_write_LINA_data = .false.
         end if
      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr, io, i
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = keep_going

         if (need_to_write_LINA_data) then
            io = 61
            open(io,file=trim(s% x_character_ctrl(10)),status='unknown')
            write(io, '(99d16.5)') s% RSP_mass, s% RSP_L, s% RSP_Teff, &
               (s% rsp_LINA_periods(i), s% rsp_LINA_growth_rates(i), i=1, s% RSP_nmodes)
            close(io)
            write(*,*) 'write ' // trim(s% x_character_ctrl(10))
            need_to_write_LINA_data = .false.
         end if
      end function extras_start_step


      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr, mode, min_periods, required_stable_periods
         integer :: max_settle_periods, target_color_periods
         integer :: nperiod, rel_period
         real(dp) :: grekm_tol, grekm_abs, rel_run_E_err
         type (star_info), pointer :: s

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         mode = s% x_integer_ctrl(1)
         nperiod = s% rsp_num_periods

         select case (mode)

         case (mode_settle)
            min_periods = max(0, s% x_integer_ctrl(2))
            required_stable_periods = max(1, s% x_integer_ctrl(3))
            max_settle_periods = max(min_periods, s% x_integer_ctrl(4))
            grekm_tol = s% x_ctrl(2)
            if (grekm_tol <= 0d0) grekm_tol = 1d-3

            ! Only check once per completed pulsation period.
            if (nperiod <= last_checked_period) return
            last_checked_period = nperiod

            if (nperiod < min_periods) then
               if (mod(nperiod, 25) == 0) then
                  write(*,'(A,I0,A,I0)') 'settling RSP: period ', nperiod, &
                     ' of minimum ', min_periods
               end if
               return
            end if

            grekm_abs = abs(s% rsp_GREKM)
            if (grekm_abs <= grekm_tol) then
               stable_period_count = stable_period_count + 1
            else
               stable_period_count = 0
            end if

            write(*,'(A,I0,A,1PE12.4,A,I0,A,I0)') 'settle check: period ', nperiod, &
               ' abs(rsp_GREKM)=', grekm_abs, ' stable_count=', &
               stable_period_count, '/', required_stable_periods

            if (stable_period_count >= required_stable_periods) then
               call report_energy_error(s)
               write(*,'(A)') 'RSP settling criterion satisfied; saving settled.mod and terminating stage 1.'
               extras_finish_step = terminate
               return
            end if

            if (max_settle_periods > 0 .and. nperiod >= max_settle_periods) then
               call report_energy_error(s)
               write(*,'(A,I0,A)') 'WARNING: reached maximum settling period ', &
                  max_settle_periods, ' before satisfying stability criterion.'
               write(*,'(A)') 'Saving current model as settled.mod anyway so stage 2 can proceed.'
               extras_finish_step = terminate
               return
            end if

         case (mode_colors)
            target_color_periods = max(1, s% x_integer_ctrl(5))

            if (colors_start_period < 0) then
               colors_start_period = nperiod
               write(*,'(A,I0,A,I0,A)') 'Colors stage starts at completed period ', &
                  colors_start_period, '; target is +', target_color_periods, ' periods.'
               return
            end if

            rel_period = nperiod - colors_start_period
            if (nperiod > last_checked_period) then
               last_checked_period = nperiod
               if (mod(rel_period, 10) == 0) then
                  write(*,'(A,I0,A,I0)') 'Colors production: completed additional periods ', &
                     rel_period, '/', target_color_periods
               end if
            end if

            if (rel_period >= target_color_periods) then
               call report_energy_error(s)
               write(*,'(A,I0,A)') 'Completed Colors production run of ', &
                  target_color_periods, ' additional periods; saving colors_final.mod.'
               extras_finish_step = terminate
               return
            end if

         case default
            ! No custom termination mode requested.
            return

         end select
      end function extras_finish_step


      subroutine report_energy_error(s)
         type (star_info), pointer :: s
         real(dp) :: rel_run_E_err
         if (s% total_energy /= 0d0) then
            rel_run_E_err = s% cumulative_energy_error/s% total_energy
            write(*,*) 'rel_run_E_err', rel_run_E_err
            if (abs(rel_run_E_err) > 1d-5) then
               write(*,*) '*** WARNING: large rel_run_E_error ***', rel_run_E_err
            end if
         end if
      end subroutine report_energy_error


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns

end module run_star_extras
