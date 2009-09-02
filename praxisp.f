c This Fortran version of PRAXIS 
c (supplied by Linda Kaufman, lck@research.att.com)
c calls Port utilities to get scratch storage.
c
c R. Brent (rpb@phys4.anu.oz.au), 31 July 1990.
c
      subroutine fminn(t0,h0,n,x,f,fmin)
      real t0,h0,x(n),f,fmin
      external f
c
c     fminn returns the minimum of the function f(x,n) of n variables
c     using the principal axis method.  The gradient of the function is
c     not required.
c
c     For a description of the algorithm, see chapter 7 of
c     "Algorithms for finding zeros and extrema of functions without
c     calculating derivatives" by Richard P Brent (Prentice-Hall, 1973).
c
c     The original Algol W version of praxis was translated into Fortran
c     by Linda Kaufman.
c
c     The input parameters are
c
c     t0       is a tolerance.  fminn attempts to return fmin=f(x)
c              such that if x0 is the true local minimum near x, then
c              norm(x-x0) ) t0 + squareroot(machep)*norm(x).
c              where machep is the machine precision
c     h0       is the maximum step size.  h0 should be set to about the
c              maximum distance from the initial guess to the minimum.
c              (if h0 is set too large or too small, the initial rate of
c              convergence may be slow.)
c     n        (at least two) is the number of variables upon which
c              the function depends.
c     x        is an array containing on entry a guess of the point of
c              minimum, on return the estimated point of minimum.
c     f(x,n)   is the function to be minimized.  f should be a real
c              function declared external in the calling program.
c
c output parameters are
c
c     fmin     the value of the function at the minimum
c     x        the solution of the problem
c     the approximating quadratic form is
c              q(x*) = f(x,n) + (1/2) * (x*-x)-transpose * a * (x*-x)
c     where x is the best estimate of the minimum and a is
c              inverse(v-transpose) * d * inverse(v)
c     (v(*,*) is the matrix of search directions\ d(*) is the array
c     of second differences).  if f has continuous second derivatives
c     near x0, a will tend to the hessian of f at x0 as x approaches x0.
c
c     It is assumed that on floating-point underflow the result is set
c     to zero.
c     The user should observe the comment on heuristic numbers 
c
c     Error condition*
c       n.lt.1   fatal
c
c Extra storage required: 7n+n*n real storage locations
c
       logical illc
       integer ktm
       real scbd
       common /cstak/ d
       double precision d(500)
       real r(1000)
       equivalence (d(1),r(1))
       call enter(1)
       if (n.lt.1) call seterr(12hfminn-n.lt.1,12,1,2)
       id=istkgt(n,3)
       iy=istkgt(n,3)
       iz= istkgt(n,3)
       iq1=istkgt(n,3)
       iq0=istkgt(n,3)
       ie=istkgt(n,3)
       iv=istkgt(n*n,3)
c      Heuristic numbers
c      If the axes are badly scaled,set scbd to 10, otherwise
c      set it to 1. If the problem is known to be ill conditioned,
c      set illc to true, otherwise false.
c      ktm is the number of iterations without improvement
c      before the algorithm terminates. ktm is 4 is
c      very cautious/usually ktm=1 is sufficient.
       scbd=1.0
       illc=.false.
       ktm=1
       eps= r1mach(4)
       call praxis(t0,eps,h0,n,x,f,fmin,scbd,ktm,illc,
     1 r(id),r(iy),r(iz),r(iq0),r(iq1),r(iv),r(ie))
       call leave
       return
       end
c
       subroutine praxis(t0,machep,h0,n,x,f,fmin,scbd,ktm,
     1 illc,d,y,z,q0,q1,v,e)
      real t0,machep,h0,x(n),f,fmin
      external f
      logical illc
      integer nl,nf,kl,kt,ktm
      real s,sl,dn,dmin,fx,f1,lds,ldt,t,h,sf,df,qf1,qd0,qd1,qa,qb,qc
      real m2,m4,small,vsmall,large,vlarge,scbd,ldfac,t2,dni,value
c
      real d(n),y(n),z(n),q0(n),q1(n),v(n,n)
      real e(n)
       logical kenplt
       common /kenpt/ kenplt
      common /global/ fx,ldt,dmin,nf,nl
     .       /q/ qa,qb,qc,qd0,qd1,qf1
      idim=n
      kenplt=.false.
c
c.....initialization.....
c     machine dependent numbers[
c
      small=machep*machep
      vsmall=small*small
      large=1.0/small
      vlarge=1.0/vsmall
      m2=sqrt(machep)
      m4=sqrt(m2)
c
c
      ldfac=0.01
      if (illc) ldfac=0.1
      kt=0
      nl=0
      nf=1
      fx=f(x,n)
      qf1=fx
      t=small+abs(t0)
      t2=t
      dmin=small
      h=h0
      if (h.lt.100.0*t) h=100.0*t
      ldt=h
c.....the first set of search directions v is the identity matrix.....
      do 20 i=1,n
           do 10 j=1,n
10              v(i,j)=0.0
20         v(i,i)=1.0
      d(1)=0.0
      qd0=0.0
      do 30 i=1,n
           q0(i)=x(i)
30         q1(i)=x(i)
c
c.....the main loop starts here.....
40    sf=d(1)
      d(1)=0.0
      s=0.0
c
c.....minimize along the first direction v(*,1).
c     fx must be passed to min by value.
      value=fx
      call min(n,1,2,d(1),s,value,.false.,f,x,t,machep,h,
     1v,q0,q1)
      if (s.gt.0.0) go to 50
           do 45 i=1,n
45              v(i,1)=-v(i,1)
50    if (sf.gt.0.9*d(1).and.0.9*sf.lt.d(1)) go to 70
           do 60 i=2,n
60              d(i)=0.0
c
c.....the inner loop starts here.....
70    do 170 k=2,n
           do 75 i=1,n
75              y(i)=x(i)
           sf=fx
           if (kt.gt.0) illc=.true.
80         kl=k
           df=0.0
c
c.....a random step follows (to avoid resolution valleys).
c     praxis assumes that random returns a random number uniformly
c     distributed in (0,1).
c
           if(.not.illc) go to 95
                do 90 i=1,n
                    s=(0.1*ldt+t2*(10.0**kt))*(uni(0)-0.5)
                     z(i)=s
                     do 85 j=1,n
85                        x(j)=x(j)+s*v(j,i)
90              continue
                fx=f(x,n)
                nf=nf+1
c
c.....minimize along the ?non-conjugate? directions v(*,k),...,v(*,n)
c
95         do 105 k2=k,n
                sl=fx
                s = 0.0
                value=fx
                call min(n,k2,2,d(k2),s,value,.false.,f,x,t,machep,h,
     1           v,q0,q1)
          write(6,900)(x(i),i=1,n)
           write(6,901)fx,nf
900       format(" current x",5e15.7)
901      format(" functions value is",e15.7," after evaluation",i5)
                if (illc) go to 97
                     s=sl-fx
                     go to 99
97              s=d(k2)*((s+z(k2))**2)
99              if (df.gt.s) go to 105
                     df=s
                     kl=k2
105        continue
           if (illc.or.(df.ge.abs((100.*machep)*fx))) go to 110
c
c.....if there was not much improvement on the first try, set
c     illc=true and start the inner loop again.....
c
           illc=.true.
           go to 80
110        continue
c
c.....minimize along the ?conjugate? directions v(*,1),...,v(*,k-1)
c
           km1=k-1
           do 120 k2=1,km1
           s=0
           value=fx
           call min(n,k2,2,d(k2),s,value,.false.,f,x,t,machep,h,
     1      v,q0,q1)
         write(6,900)(x(i),i=1,n)
         write(6,901)fx,nf
120        continue
           f1=fx
           fx=sf
           lds=0
           do 130 i=1,n
                sl=x(i)
                x(i)=y(i)
                sl=sl-y(i)
                y(i)=sl
130             lds=lds+sl*sl
           lds=sqrt(lds)
           if (lds.le.small) go to 160
c
c.....discard direction v(*,kl).
c     if no random step was taken, v(*,kl) is the ?non-conjugate?
c     direction along which the greatest improvement was made.....
c
           klmk=kl-k
           if (klmk.lt.1) go to 141
           do 140 ii=1,klmk
                i=kl-ii
                do 135 j=1,n
135                  v(j,i+1)=v(j,i)
140             d(i+1)=d(i)
141        d(k)=0
           do 145 i=1,n
145             v(i,k)=y(i)/lds
c
c.....minimize along the new ?conjugate? direction v(*,k), which is
c     the normalized vector[  (new x) - (0ld x).....
c
           value=f1
           call min(n,k,4,d(k),lds,value,.true.,f,x,t,machep,h,
     1     v,q0,q1)
          write(6,902)(x(i),i=1,n)
          write(6,901)fx,nf
902      format(" x at end of inner loop",4e15.7)
         kenplt=.true.
        fff=f(x,n)
         kenplt=.false.
          if (lds.gt.0.0) go to 160
                lds=-lds
                do 150 i=1,n
150                  v(i,k)=-v(i,k)
160        ldt=ldfac*ldt
           if (ldt.lt.lds) ldt=lds
           t2 = 0.0
           do 165 i=1,n
165             t2=t2+x(i)**2
           t2=m2*sqrt(t2)+t
c
c.....see whether the length of the step taken since starting the
c     inner loop exceeds half the tolerance.....
c
           if (ldt.gt.(0.5*t2)) kt=-1
           kt=kt+1
           if (kt.gt.ktm) go to 400
170   continue
c.....the inner loop ends here.
c
c     try quadratic extrapolation in case we are in a curved valley.
c
171   call quad(n,f,x,t,machep,h,v,q0,q1)
      dn=0.0
      do 175 i=1,n
           d(i)=1.0/sqrt(d(i))
           if (dn.lt.d(i)) dn=d(i)
175   continue
      do 180 j=1,n
           s=d(j)/dn
           do 180 i=1,n
180             v(i,j)=s*v(i,j)
c
c.....scale the axes to try to reduce the condition number.....
c
      if (scbd.le.1.0) go to 200
           s=vlarge
           do 185 i=1,n
                sl=0.0
                do 182 j=1,n
182                  sl=sl+v(i,j)*v(i,j)
                z(i)=sqrt(sl)
                if (z(i).lt.m4) z(i)=m4
                if (s.gt.z(i)) s=z(i)
185        continue
           do 195 i=1,n
                sl=s/z(i)
                z(i) = 1.0/sl
                if (z(i).le.scbd) go to 189
                     sl = 1.0/scbd
                     z(i)=scbd
189             do 190 j=1,n
190                  v(i,j)=sl*v(i,j)
195        continue
c
c.....calculate a new set of orthogonal directions before repeating
c     the main loop.
c     first transpose v for minfit[
c
200   do 220 i=2,n
           im1=i-1
           do 210 j=1,im1
                s=v(i,j)
                v(i,j)=v(j,i)
210             v(j,i)=s
220   continue
c
c.....call minfit to find the singular value decomposition of v.
c     this gives the principal values and principal directions of the
c     approximating quadratic form without squaring the condition
c     number.....
c
      call minfit(idim,n,machep,vsmall,v,d,e)
c
c.....unscale the axes.....
c
      if (scbd.le.1.0) go to 250
           do 230 i=1,n
                s=z(i)
                do 225 j=1,n
225                  v(i,j)=s*v(i,j)
230        continue
           do 245 i=1,n
                s=0.0
                do 235 j=1,n
235                  s=s+v(j,i)**2
                s=sqrt(s)
                d(i)=s*d(i)
                s=1.0/s
                do 240 j=1,n
240                  v(j,i)=s*v(j,i)
245        continue
c
c
250   do 270 i=1,n
           dni=dn*d(i)
           if (dni.gt.large) go to 265
                if (dni.lt.small) go to 260
                     d(i)=1.0/(dni*dni)
                     go to 270
260             d(i)=vlarge
                go to 270
265        d(i)=vsmall
270   continue
c
c.....sort the eigenvalues and eigenvectors.....
c
      call sort(idim,n,d,v)
      dmin=d(n)
      if (dmin.lt.small) dmin=small
      illc=.false.
      if (m2*d(1).gt.dmin) illc=.true.
c.....the main loop ends here.....
c
      go to 40
c
c.....return.....
c
400       continue
      fmin=fx
      return
      end
      subroutine minfit(m,n,machep,tol,ab,q,e)
      real machep
      dimension ab(m,n),q(n)
c...an improved version of minfit (see golub and reinsch, 1969)
c   restricted to m=n,p=0.
c   the singular values of the array ab are returned in q and ab is
c   overwritten with the orthogonal matrix v such that u.diag(q) = ab.v,
c   where u is another orthogonal matrix.
         real e(n)
c...householder*s reduction to bidiagonal form...
      if (n.eq.1) go to 200
      eps = machep
      g = 0.0
      x = 0.0
      do 11 i=1,n
         e(i) = g
         s = 0.0
         l = i + 1
         do 1 j=i,n
1           s = s + ab(j,i)**2
         g = 0.0
         if (s.lt.tol) go to 4
            f = ab(i,i)
           g = sqrt(s)
            if (f.ge.0.0) g = -g
            h = f*g - s
            ab(i,i)=f-g
            if (l.gt.n) go to 4
            do 3 j=l,n
               f = 0.0
               do 2 k=i,n
2                 f = f + ab(k,i)*ab(k,j)
               f = f/h
               do 3 k=i,n
3                 ab(k,j) = ab(k,j) + f*ab(k,i)
4        q(i) = g
         s = 0.0
         if (i.eq.n) go to 6
         do 5 j=l,n
5           s = s + ab(i,j)*ab(i,j)
6        g = 0.0
         if (s.lt.tol) go to 10
            if (i.eq.n) go to 16
            f = ab(i,i+1)
16          g = sqrt(s)
            if (f.ge.0.0) g = -g
            h = f*g - s
            if (i.eq.n) go to 10
            ab(i,i+1) = f - g
            do 7 j=l,n
7              e(j) = ab(i,j)/h
            do 9 j=l,n
               s = 0.0
               do 8 k=l,n
8                 s = s + ab(j,k)*ab(i,k)
               do 9 k=l,n
9                 ab(j,k) = ab(j,k) + s*e(k)
10       y = abs(q(i)) + abs(e(i))
11       if (y.gt.x) x = y
c...accumulation of right-hand transformations...
      ab(n,n) = 1.0
      g = e(n)
      l = n
      do 25 ii=2,n
         i = n - ii + 1
         if (g.eq.0.0) go to 23
         h = ab(i,i+1)*g
         do 20 j=l,n
20          ab(j,i) = ab(i,j)/h
         do 22 j=l,n
            s = 0.0
            do 21 k=l,n
21             s = s + ab(i,k)*ab(k,j)
            do 22 k=l,n
22             ab(k,j) = ab(k,j) + s*ab(k,i)
23       do 24 j=l,n
            ab(i,j) = 0.0
24          ab(j,i) = 0.0
         ab(i,i) = 1.0
         g = e(i)
25       l = i
c...diagonalization of the bidiagonal form...
100   eps = eps*x
      do 150 kk=1,n
         k = n - kk + 1
         kt = 0
101      kt = kt + 1
         if (kt.le.30) go to 102
            e(k) = 0.0
             call seterr(22hfminn-problem with svd,22,2,2)
             return
102      do 103 ll2=1,k
            l2 = k - ll2 + 1
            l = l2
            if (abs(e(l)).le.eps) go to 120
            if (l.eq.1) go to 103
            if (abs(q(l-1)).le.eps) go to 110
103         continue
c...cancellation of e(l) if l'1...
110      c = 0.0
         s = 1.0
         do 116 i=l,k
            f = s*e(i)
            e(i) = c*e(i)
            if (abs(f).le.eps) go to 120
            g = q(i)
c...q(i) = h = dsqrt(g*g + f*f)...
            if (abs(f).lt.abs(g)) go to 113
            if (f) 112,111,112
111         h = 0.0
            go to 114
112         h = abs(f)*sqrt(1.0 + (g/f)**2)
            go to 114
113         h = abs(g)*sqrt(1.0+ (f/g)**2)
114         q(i) = h
            if (h.ne.0.0) go to 115
               g = 1.0
               h = 1.0
115         c = g/h
116         s = -f/h
c...test for convergence...
120      z = q(k)
         if (l.eq.k) go to 140
c...shift from bottom 2*2 minor...
         x = q(l)
         y = q(k-1)
         g = e(k-1)
         h = e(k)
         f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0*h*y)
         g = sqrt(f*f + 1.0)
         temp = f - g
         if (f.ge.0.0) temp = f + g
         f = ((x - z)*(x + z) + h*(y/temp - h))/x
c...next qr transformation...
         c = 1.0
         s = 1.0
         lp1 = l + 1
         if (lp1.gt.k) go to 133
         do 132 i=lp1,k
            g = e(i)
            y = q(i)
            h = s*g
            g = g*c
            if (abs(f).lt.abs(h)) go to 123
            if (f) 122,121,122
121         z = 0.0
            go to 124
122         z = abs(f)*sqrt(1.0+ (h/f)**2)
            go to 124
123         z = abs(h)*sqrt(1.0+ (f/h)**2)
124         e(i-1) = z
            if (z.ne.0.0) go to 125
               f = 1.0
               z = 1.0
125         c = f/z
            s = h/z
            f = x*c + g*s
            g = -x*s + g*c
            h = y*s
            y = y*c
            do 126 j=1,n
               x = ab(j,i-1)
               z = ab(j,i)
               ab(j,i-1) = x*c + z*s
126            ab(j,i) = -x*s + z*c
            if (abs(f).lt.abs(h)) go to 129
            if (f) 128,127,128
127         z = 0.0
            go to 130
128         z = abs(f)*sqrt(1.0+ (h/f)**2)
            go to 130
129         z = abs(h)*sqrt(1.0+ (f/h)**2)
130         q(i-1) = z
            if (z.ne.0.0) go to 131
               f = 1.0
               z = 1.0
131         c = f/z
            s = h/z
            f = c*g + s*y
132         x = -s*g + c*y
133      e(l) = 0.0
         e(k) = f
         q(k) = x
         go to 101
c...convergence[  q(k) is made non-negative...
140      if (z.ge.0.0) go to 150
         q(k) = -z
         do 141 j=1,n
141         ab(j,k) = -ab(j,k)
150      continue
      return
200   q(1) = ab(1,1)
      ab(1,1) = 1.0
      return
      end
      subroutine min(n,j,nits,d2,x1,f1,fk,f,x,t,machep,h,
     1v,q0,q1)
      external f
      logical fk
      real machep,x(n),ldt
      dimension v(n,n),q0(n),q1(n)
      common /global/ fx,ldt,dmin,nf,nl
     .       /q/ qa,qb,qc,qd0,qd1,qf1
c...the subroutine min minimizes f from x in the direction v(*,j) unless
c   j is less than 1, when a quadratic search is made in the plane
c   defined by q0,q1,x.
c   d2 is either zero or an approximation to half f?.
c   on entry, x1 is an estimate of the distance from x to the minimum
c   along v(*,j) (or, if j=0, a curve).  on return, x1 is the distance
c   found.
c   if fk=.true., then f1 is flin(x1).  otherwise x1 and f1 are ignored
c   on entry unless final fx is greater than f1.
c   nits controls the number of times an attempt will be made to halve
c   the interval.
      logical dz
      real m2,m4
      small = machep**2
      m2=sqrt(machep)
      m4 = sqrt(m2)
      sf1 = f1
      sx1 = x1
      k = 0
      xm = 0.0
      fm = fx
      f0 = fx
      dz = d2.lt.machep
c...find the step size...
      s = 0.0
      do 1 i=1,n
1        s = s + x(i)**2
      s = sqrt(s)
      temp = d2
      if (dz) temp = dmin
      t2 = m4*sqrt(abs(fx)/temp + s*ldt) + m2*ldt
      s = m4*s + t
      if (dz.and.t2.gt.s) t2 = s
      t2=amax1(t2,small)
      t2=amin1(t2,.01*h)
      if (.not.fk.or.f1.gt.fm) go to 2
      xm = x1
      fm = f1
2     if (fk.and.abs(x1).ge.t2) go to 3
      temp = 1.0
      if (x1.lt.0.0) temp = -1.0
      x1=temp*t2
      f1 = flin(n,j,x1,f,x,nf,v,q0,q1)
3     if (f1.gt.fm) go to 4
      xm = x1
      fm = f1
4     if (.not.dz) go to 6
c...evaluate flin at another point and estimate the second derivative...
      x2 = -x1
      if (f0.ge.f1) x2 = 2.0*x1
      f2 = flin(n,j,x2,f,x,nf,v,q0,q1)
      if (f2.gt.fm) go to 5
         xm = x2
         fm = f2
5     d2 = (x2*(f1 - f0)-x1*(f2 - f0))/((x1*x2)*(x1 - x2))
c...estimate the first derivative at 0...
6     d1 = (f1 - f0)/x1 - x1*d2
      dz = .true.
c...predict the minimum...
      if (d2.gt.small) go to 7
         x2 = h
         if (d1.ge.0.0) x2 = -x2
         go to 8
7        x2 = (-0.5*d1)/d2
8     if (abs(x2).le.h) go to 11
         if (x2) 9,9,10
9        x2 = -h
         go to 11
10       x2 = h
c...evaluate f at the predicted minimum...
11    f2 = flin(n,j,x2,f,x,nf,v,q0,q1)
      if (k.ge.nits.or.f2.le.f0) go to 12
c...no success, so try again...
         k = k + 1
         if (f0.lt.f1.and.(x1*x2).gt.0.0) go to 4
         x2 = 0.5*x2
         go to 11
c...increment the one-dimensional search counter...
12    nl = nl + 1
      if (f2.le.fm) go to 13
      x2 = xm
      go to 14
13    fm = f2
c...get a new estimate of the second derivative...
14    if (abs(x2*(x2-x1)).le.small) go to 15
         d2 = (x2*(f1-f0) - x1*(fm-f0))/((x1*x2)*(x1 - x2))
         go to 16
15       if (k.gt.0) d2 = 0.0
16    if (d2.le.small) d2 = small
      x1 = x2
      fx = fm
      if (sf1.ge.fx) go to 17
         fx = sf1
         x1 = sx1
c...update x for linear but not parabolic search...
17    if (j.eq.0) return
      do 18 i=1,n
18       x(i) = x(i) + x1*v(i,j)
      return
      end
      real function flin(n,j,l,f,x,nf,v,q0,q1)
      real l,x(n)
      dimension v(n,n),q0(n),q1(n)
c...flin is the function of one real variable l that is minimized
c   by the subroutine min...
      common /q/ qa,qb,qc,qd0,qd1,qf1
      common /cstack/ d(500)
      double precision d(500)
      real  t(1000)
      equivalence(d(1),t(1))
       call enter(1)
      it =istkgt(n,3)
      if (j .eq. 0) go to 2
c...the search is linear...
      it1=it-1
      do 1 i=1,n
          itt=it1+i
1        t(itt) = x(i) + l*v(i,j)
      go to 4
c...the search is along a parabolic space curve...
2     qa = (l*(l - qd1))/(qd0*(qd0 + qd1))
      qb = ((l + qd0)*(qd1 - l))/(qd0*qd1)
      qc = (l*(l + qd0))/(qd1*(qd0 + qd1))
      it1=it-1
      do 3 i=1,n
           itt=it1+i
3        t(itt) = (qa*q0(i) + qb*x(i)) + qc*q1(i)
c...the function evaluation counter nf is incremented...
4     nf = nf + 1
      flin = f(t(it),n)
      call leave
      return
      end
      subroutine sort(m,n,d,v)
      real d(n),v(m,n)
c...sorts the elements of d(n) into descending order and moves the
c   corresponding columns of v(n,n).
c   m is the row dimension of v as declared in the calling program.
      real s
      if (n.eq.1) return
      nm1 = n - 1
      do 3 i = 1,nm1
         k=i
         s = d(i)
         ip1 = i + 1
         do 1 j = ip1,n
            if (d(j) .le. s) go to 1
            k = j
            s = d(j)
1           continue
         if (k .le. i) go to 3
         d(k) = d(i)
         d(i) = s
         do 2 j = 1,n
            s = v(j,i)
            v(j,i) = v(j,k)
2           v(j,k) = s
3        continue
      return
      end
      subroutine quad(n,f,x,t,machep,h,v,q0,q1)
      external f
c...quad looks for the minimum of f along a curve defined by q0,q1,x...
      real x(n),machep,ldt,l
      dimension v(n,n),q0(n),q1(n)
      common /global/ fx,ldt,dmin,nf,nl
     .       /q/ qa,qb,qc,qd0,qd1,qf1
      s = fx
      fx = qf1
      qf1 = s
      qd1 = 0.0
      do 1 i=1,n
         s = x(i)
         l = q1(i)
         x(i) = l
         q1(i) = s
1        qd1 = qd1+ (s-l)**2
      qd1 = sqrt(qd1)
      l = qd1
      s = 0.0
      if (qd0.le.0.0.or.qd1.le.0.0.or.nl.lt.3*n*n) go to 2
      value=qf1
      call min(n,0,2,s,l,value,.true.,f,x,t,machep,h,v,q0,q1)
      qa = (l*(l-qd1))/(qd0*(qd0+qd1))
      qb = ((l+qd0)*(qd1-l))/(qd0*qd1)
      qc = (l*(l+qd0))/(qd1*(qd0+qd1))
      go to 3
2     fx = qf1
      qa = 0.0
      qb = qa
      qc = 1.0
3     qd0 = qd1
      do 4 i=1,n
         s = q0(i)
         q0(i) = x(i)
4        x(i) = (qa*s + qb*x(i)) + qc*q1(i)
      return
      end

