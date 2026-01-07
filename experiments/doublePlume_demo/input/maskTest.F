      program test_mask
  
      implicit none
      integer, parameter :: nx=50, ny=11
      REAL*8:: mask(nx,ny)
      integer :: iounit, ios

      open(unit=10, file="plumeMask.bin", form="unformatted", 
     &  access="direct", recl=nx*ny*8, status="old",
     & convert = "big_endian", iostat=ios)

      if (ios /= 0) then
        print *, "Open failed, ios=", ios
        stop
      end if

      read(10, rec=1, iostat=ios) mask
      if (ios /= 0) then
        print *, "Read failed, ios=", ios
        stop
      end if

      print *, "First 10 values:"
      print *, mask
      end program test_mask
