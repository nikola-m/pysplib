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
      n2 = (n-1)/2
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

      n2 = (n-1)/2
      do j = 0, n2-1
        nj = n-j-1
        x(j) = dble(q(j))
        x(nj) = -dimag(q(nj))
      enddo
      if (mod(n,2).ne.0) x(n2) = dble(q(n2)) 
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

      n2 = (n-1)/2
      do j = 0, n2-1
        nj = n-j-1
        x(j) = dble(q(j))
        x(nj) = -dimag(q(j+1))
      enddo
      if (mod(n,2).ne.0) x(n2) = dble(q(n2)) 
      return
      end
