module mod_io

  ! Provides input/output procedures.

  use iso_fortran_env, only: int32, real32

  implicit none

  private
  public :: write_field

contains

  subroutine write_field(field, fieldname, time)
    ! Writes a field into a binary file.
    real(real32), intent(in) :: field(:,:)
    character(*), intent(in) :: fieldname
    integer(int32), intent(in) :: time
    integer(int32) :: fileunit, record_length
    character(100) :: filename, timestr
    write(timestr, '(i4.4)') time
    filename = 'tsunami_' // fieldname // '_' // trim(timestr) // '.dat'
    record_length = storage_size(field) / 8 * size(field)
    open(newunit=fileunit, file=filename, access='direct', recl=record_length)
    write(unit=fileunit, rec=1) field
    close(fileunit)
  end subroutine write_field



  subroutine getin_p(filename, param_name, param_value)
      character(len=*), intent(in) :: filename, param_name
      integer, intent(inout) :: param_value
      integer :: ierr
      character(100) :: line
      logical :: found
      
      ! Open the file
      open(unit=10, file=filename, status='old', action='read', iostat=ierr)
      if (ierr /= 0) then
          print *, "Error opening file ", filename
          stop
      end if
      
      ! Read parameters from file
      found = .false.
      do
          read(10, '(A)', iostat=ierr) line
          if (ierr /= 0) exit
          
          if (index(line, trim(param_name)) /= 0) then
              read(line, *) param_value
              found = .true.
              exit
          end if
      end do
      
      ! Close the file
      close(10)
      
      if (.not. found) then
          print *, "Parameter ", trim(param_name), " not found in file ", filename
          stop
      end if
      
  end subroutine getin_p

end module mod_io
