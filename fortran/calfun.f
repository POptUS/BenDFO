     subroutine calfun(n,x,f)
      integer n
      double precision f, x(n)
c     *********
c
c     Subroutine calfun
c
c     This subroutine returns a function value as used in:
c
c     Benchmarking Derivative-Free Optimization Algorithms
c     Jorge J. More' and Stefan M. Wild
c     SIAM J. Optimization, Vol. 20 (1), pp.172-191, 2009.
c
c     The latest version of this subroutine is always available at
c     http://www.mcs.anl.gov/~more/dfo/
c     The authors would appreciate feedback and experiences from numerical
c     studies conducted using this subroutine.
c
c     The subroutine returns the function value f(x)
c
c     The subroutine statement is
c
c       subroutine calfun(n,x,f)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c       f is an output that contains the function value at x.
c
c
c     Additional problem descriptors are passed through the common block
c     calfun_int containing:
c       m a positive integer (length of output from dfovec).
c          m must not exceed n.
c       nprob is a positive integer that defines the number of the problem.
c          nprob must not exceed 22.
c       job is a positive integer specifying the type of problem desired:
c           1 corresponds to smooth problems
c           2 corresponds to piecewise-smooth problems
c           3 corresponds to deterministically noisy problems
c           4 corresponds to stochastically noisy problems
c       nseed is a random number generator seed needed when job=4
c
c     To store the evaluation history, additional variables are passed 
c     through the common blocks calfun_int and calfun_fevals. These
c     may be commented out if a user desires. They are:
c       nfev is a non-negative integer containing the number of function 
c          evaluations done so far (nfev=0 is a good default).
c          nfev should be less than 1500 unless fevals is modified.
c          after calling calfun, nfev will be incremented by one.
c       np is a counter for the test problem number. np=1 is a good
c          default if only a single problem/run will be done.
c          np should be no bigger than 100 unless fevals is modified.
c       fevals(1500,100) is an array containing the history of function
c          values, the entry fevals(nfev+1,np) being updated here.
c
c     Argonne National Laboratory
c     Jorge More' and Stefan Wild. January 2008.
c     **********

      double precision dasum, dnrm2

c     Function value array
      double precision phi, epsf, xmax, xnorm1, xnorm2
      double precision y(100)
      double precision fvec(1500)

      integer q, nprob, nfev, np, nseed, job
      common /calfun_int/ q, nprob, nfev, np, nseed, job
      double precision fevals(1500,100)
      common /calfun_fevals/ fevals

      if (job .eq. 1) then
         call dfovec(q,n,x,fvec,nprob)
         f = dnrm2(q,fvec,1)**2
      else if (job .eq. 2) then
         call dfovec(q,n,x,fvec,nprob)
         f = dasum(q,fvec,1)
         if (nprob .eq. 8  .or. nprob .eq.  9 .or. nprob .eq. 13 .or.
     +       nprob .eq. 16 .or. nprob .eq. 17 .or. nprob .eq. 18) then
            do i = 1, n
               y(i) = max(x(i),0.0d0)
            end do
            call dfovec(q,n,y,fvec,nprob)
            f = dasum(q,fvec,1)
         end if
      else if (job .eq. 3) then
         call dfovec(q,n,x,fvec,nprob)
         epsf = 1.0d-3
         do i = 1, q
            fvec(i) = fvec(i)*(1 + epsf*(2*surn01(nseed)-1))
         end do
         f = dnrm2(q,fvec,1)**2
      else if (job .eq. 4) then
         call dfovec(q,n,x,fvec,nprob)
         epsf = 1.0d-3
         xmax = 0.0
         do i = 1, n
            xmax = max(xmax,abs(x(i)))
         end do
         xnorm1 = dasum(n,x,1)
         xnorm2 = dnrm2(n,x,1)
         phi = 0.9*sin(100*xnorm1)*cos(100*xmax) + 0.1*cos(xnorm2)
         phi = phi*(4*phi**2 -3)
         f = (1 + epsf*phi)*dnrm2(q,fvec,1)**2
      else
         write (*,*) "Parameter JOB does not satisfy 1 <= JOB <= 4"
      end if

      call wallclock(wctime2)

      nfev = nfev + 1
      fevals(nfev,np) = f

      end