      program test
!     Purpose:
!     Show that we can reproduce the results of the examples from NAG documents:
!     'NAG Library Routine Document: C06EBF',
!     where C06EBF and C06EAF are used to perform the FFT.
!     Our approach:
!     1) Using FFTW3 library to perform forward and backward FFT.
!     2) Output results were scaled with 1/sqrt(n) for forward and
!        with 1/n for backward transform to be able to match results from 
!        NAG routines.
!     3) Simple routines to transform complex sequences from Hermitian from
!        to complex conjugate form and vice versa.
!        Hermite form is used to reduce storage.
!
      implicit none
      include "fftw3.f"
      integer, parameter :: n=7
      double precision, dimension(0:n-1) :: x,xrec
      double complex, dimension(0:n-1) :: y
!     fftw related
      double complex, dimension(0:n-1) :: out
      integer*8 :: plan_forward
      integer*8 :: plan_backward
!     locals
      integer :: i

      ! Test when n is odd, n=7
      x = reshape([0.34907, 0.54890, 0.74776, 0.94459, 1.13850, 1.32850, 1.51370],shape(x))

      ! Test when n is even, n=6
!      x = reshape([0.34907, 0.54890, 0.74776, 0.94459, 1.13850, 1.32850],shape(x))

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Given Hermitian sequence:'
      write ( *, '(a)' ) ' '
      do i = 0, n-1
        write(*,'(2x,i4,2x,1g14.6)') i, x(i)
      end do

      call hermite_to_complex(x,n,y)
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Complex conjugate of the given Hermitian sequence:'
      write ( *, '(a)' ) ' '
      do i = 0, n-1
        write(*,'(2x,i4,2x,2g14.6)') i, y(i)
      end do

!
!  make a plan for the fft, and forward transform the data.
!
      call dfftw_plan_dft_1d_ ( plan_forward, n, y, out, &
        fftw_forward, fftw_estimate )

      call dfftw_execute_ ( plan_forward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Output fft coefficients: '
      write ( *, '(a)' ) '  NOTE: imaginary part is zero, and is discarted in the output!'
      write ( *, '(a)' ) ' '

!  note we have to divide by sqrt(n) to obtain the same output as in c06ebf of nag library
      do i = 0, n-1
        write ( *, '(2x,i4,2x,1g14.6)' ) i, dble(out(i))/sqrt(dble(n))
      end do
!
!  Now take this real sequence and do IFFT to it to recover original sequence
!  We have two ways of doing this:
!  Since the sequence 'out' can be seen as a real sequence, because it's imag
!  part is zero we can perform real FFT, using r2c transform in fftw3.
!  In that case we use dble(out) as the input, and whole out as the output
      call dfftw_plan_dft_r2c_1d_ ( plan_backward, n, dble(out), out, &
       fftw_backward, fftw_estimate )
!  ...or we can use the usual complex FFT, in that case erase comments '!#$'.
!#$      call dfftw_plan_dft_1d_ ( plan_backward, n, out, out, &
!#$        fftw_backward, fftw_estimate )

      call dfftw_execute_ ( plan_backward )

!  Convert back to Herite form
      call complex_to_hermite_r2c(out,n,xrec)
!#$      call complex_to_hermite(out,n,xrec)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recovered input data:'
      write ( *, '(a)' ) ' '

      write (*,'(a)') '          Restored:    Original:'
      do i = 0, n-1
        write ( *, '(2x,i4,2x,3g14.6)' ) i, xrec(i)/dble(n), x(i)
      end do

!
!  discard the information associated with the plans.
!
      call dfftw_destroy_plan_ ( plan_forward )
      call dfftw_destroy_plan_ ( plan_backward )

      end program test

      subroutine hermite_to_complex(x,n,q)
!     Purpose:
!     Form a complex conjugate from a Hermitian sequence of 'n' data values.
!     Funaro used FFT routines from NAG library C06EBF and C06EAF which use 
!     data in so called Hermitian sequence.
!
!     Author: Nikola Mirkov email:nmirkov@vinca.rs
!
      implicit none
      integer, intent(in) :: n
      double precision, dimension(0:n-1), intent(in) :: x
      double complex, dimension(0:n-1), intent(out) :: q 
!     locals
      integer :: j,nj,n2
      double precision, dimension(0:n-1) :: a,b
      
      a(0) = x(0)
      b(0) = 0.0d0
      n2 = int(ceiling((n-1)/2.))
      do j = 1, n2
        nj = n - j
        a(j) = x(j)
        a(nj) = x(j)
        b(j) = x(nj)
        b(nj) = -x(nj)
      enddo
      if (mod(n,2).eq.0) then
        a(n2+1) = x(n2+1)
        b(n2+1) = 0.0e0
      end if
      
      q=cmplx(a,b)

      return
      end 

      subroutine complex_to_hermite(q,n,x)
!     Purpose:
!     Form a Hermitian sequence of 'n' data values from a complex conjugate.
!     Funaro used FFT routines from NAG library C06EBF and C06EAF which use 
!     data in so called Hermitian sequence
!
!     Author: Nikola Mirkov email:nmirkov@vinca.rs
!
      implicit none
      integer, intent(in) :: n
      double complex, dimension(0:n-1), intent(in) :: q 
      double precision, dimension(0:n-1), intent(out) :: x
!     locals
      integer :: j,nj,n2

      n2 = int(ceiling((n-1)/2.))
      do j = 0, n2-1
        nj = n-j-1
        x(j) = dble(q(j))
        x(nj) = -dimag(q(nj))
      enddo
      x(n2) = dble(q(n2))  
      return
      end

      subroutine complex_to_hermite_r2c(q,n,x)
!     Purpose:
!     Form a Hermitian sequence of n data values from a complex conjugate.
!     Funaro used FFT routines from NAG library C06EBF and C06EAF which use 
!     data in so called Hermitian sequence.
!     NOTE: We use this if we performed real to compex ifft in fftw3, therefore r2c.
!
!     Author: Nikola Mirkov email:nmirkov@vinca.rs
!
      implicit none
      integer, intent(in) :: n
      double complex, dimension(0:n-1), intent(in) :: q 
      double precision, dimension(0:n-1), intent(out) :: x
!     locals
      integer :: j,nj,n2

      n2 = int(ceiling((n-1)/2.))
      do j = 0, n2-1
        nj = n-j-1
        x(j) = dble(q(j))
        x(nj) = -dimag(q(j+1))
      enddo
      x(n2) = dble(q(n2)) 
      return
      end
