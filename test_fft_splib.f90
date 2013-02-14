      program test
!     Purpose:
!     Test routines in splib that use FFT, since it now uses FFTW3.
!     For theory related to these functions FCCHGL, FVCHGL, FCCHGA
!     consult splib's manual written by D. Funaro available online.
!
      implicit none
      integer, parameter :: n=7
      real(kind=8), dimension(0:n) :: et,qn,qn1,co,co1
      real(kind=8), dimension(n) :: cs,dz,qz
      real(kind=8), dimension(0:n-1) :: co3,co4
      ! Locals
      integer :: i
      
      call zechgl(n,et)
      qn=exp(et)
      call cochgl(n,qn,co)
      call fcchgl(n,qn,co1)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Chebyshev-Fourier coefficients at Gauss-Lobatto points:'
      write ( *, '(a)' ) ' '

      write (*,'(a)') '          Matrix:       FFT:'      
      do i = 0, n
        write(*,'(2x,i4,2x,3g14.6)') i, co(i), co1(i)
      end do
      write ( *, '(a)' ) ' '
      
      call fvchgl(n,co,qn1)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Chebyshev interpolant values at Gauss-Lobatto points:'
      write ( *, '(a)' ) ' '

      write (*,'(a)') '          Original:     FFT:'       
      do i = 0, n
        write(*,'(2x,i4,2x,3g14.6)') i, qn(i), qn1(i)
      end do
      write ( *, '(a)' ) ' '

      call zechga(n,cs,dz)
      qz=exp(cs)
      call cochga(n,qz,co3)
      call fcchga(n,qz,co4)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Chebyshev-Fourier coefficients at Gauss points:'
      write ( *, '(a)' ) ' '

      write (*,'(a)') '          Matrix:       FFT:'    
      do i = 0, n-1
        write(*,'(2x,i4,2x,3g14.6)') i, co3(i), co4(i)
      end do
      write ( *, '(a)' ) ' '
  
      end program
  
  
