program tsunami

  ! Tsunami simulator.
  !
  ! Solves the non-linear 2-d shallow water equation system:
  !
  !     du/DT + u du/DX + v du/DY + G dh/DX = 0
  !     dv/DT + u dv/DX + v dv/DY + G dh/DY = 0
  !     dh/DT + d(hu)/DX + d(hv)/DY = 0
  !
  ! This version is parallelized and uses derived types.

  !! Used module
  use iso_fortran_env, only: int32, real32


  !! Define variables
  implicit none

  character(len=*), parameter :: RUNFILE = 'run.def'

  integer(int32) :: im, jm
  integer(int32) :: num_time_steps
  real(real32) :: decay

  !! Define MPI communicator
  
  ! Call the subroutine to read parameters from file
  ! im = getin_p(RUNFILE, 'X_LENGTH')
  ! jm = getin_p(RUNFILE, 'Y_LENGTH')
  ! num_time_steps =  getin_p(RUNFILE, 'TIME_STEPS')
  ! decay = getin_p(RUNFILE, 'DECAY')
  
  ! Display the read parameters
  ! print *, "Grid size (X_LENGTH, Y_LENGTH): ", im, jm
  ! print *, "Time steps (TIME_STEPS): ", num_time_steps
  ! print *, "Model parameter (DECAY): ", decay

  !! Read configs


  !! Get the grids distributed on all processors


  !! Prepare for the simulation

  ! type(Field) :: h, hm, u, v

  ! real(real32) :: hmin, hmax, hmean

  ! u = Field('u', [im, jm])
  ! v = Field('v', [im, jm])
  ! h = Field('h', [im, jm])
  ! hm = Field('h_mean', [im, jm])


  ! ! initialize a gaussian blob in the center
  ! call h % init_gaussian(decay, ic, jc)

  ! hm = 10.

  ! call h % write(0)

  ! !! Time loop

  ! time_loop: do n = 1, num_time_steps

  !   ! compute u at next time step
  !   u = u - (u * diffx(u) / DX + v * diffy(u) / DY &
  !     + G * diffx(h) / DX) * DT

  !   ! compute v at next time step
  !   v = v - (u * diffx(v) / DX + v * diffy(v) / DY &
  !     + G * diffy(h) / DY) * DT

  !   ! compute h at next time step
  !   h = h - (diffx(u * (hm + h)) / DX &
  !          + diffy(v * (hm + h)) / DY) * DT

  !   hmin = minval(h % data)
  !   call co_min(hmin, 1)

  !   hmax = maxval(h % data)
  !   call co_max(hmax, 1)

  !   hmean = sum(h % data(h % lb(1):h % ub(1),h % lb(2):h % ub(2))) &
  !         / size(h % data(h % lb(1):h % ub(1),h % lb(2):h % ub(2)))
  !   call co_sum(hmean, 1)
  !   hmean = hmean / num_images()

  !   if (this_image() == 1) &
  !     print '(a, i5, 3(f10.6))', &
  !       'step, min(h), max(h), mean(h):', &
  !        n, hmin, hmax, hmean

  !   call h % write(n)

  ! end do time_loop

  !! Close everything

end program tsunami
