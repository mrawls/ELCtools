          program checkfit
c
c   This will use the command line arguments.
c
c
c
          implicit double precision(a-h,o-z)
          dimension xmodel(9990000),ymodel(9990000)
          dimension xkeep(9990000)
          dimension ykeep(9990000)
          dimension cxmodel(9990000),cymodel(9990000)
          character*40 filein
          character*80 line
          real*8 zzy,zzx
c
          M=iargc()    ! find out how many arguments are on the command line
c
          call getvariables(M,filein)
c
          open(unit=20,file=filein)
c
          ymin=900999999990.d0
          ymax=-1000.d0
          icount=0
          Ndot=0
          Nmodel=0
          do 10 i=1,9990000                            ! read in the model here
            read(20,16,end=15)line
            call countdot(line,Ndot)
            if(Ndot.eq.2)then
              backspace(20)
              read(20,*,end=15)zzx,zzy
              if(zzy.lt.0.0d0)go to 10
              icount=icount+1
              if(zzy.lt.ymin)ymin=zzy
            endif
 10       continue
 15       Nmodel=icount
16        format(a80)
          close(20)
c
          write(*,*)'N=',Nmodel,ymin
          open(unit=20,file=filein)
c
c          ymin=900999999990.d0
c          ymax=-1000.d0
          icount=0
          do 120 i=1,9990000                       ! read in the model here
            read(20,16,end=150)line
            call countdot(line,Ndot)
            if(Ndot.eq.2)then
              backspace(20)
              read(20,*,end=150)zzx,zzy
              if(zzy.lt.0.0d0)go to 120
              if(zzy.lt.ymin+50.0d0)then
                icount=icount+1
                xmodel(icount)=(zzx)
                ymodel(icount)=(zzy)
              endif
            endif
 120      continue
 150      Nmodel=icount

          close(20)
c

          Kmodel=Nmodel
          call sort2(Nmodel,ymodel,xmodel)
c
          chimin=ymodel(1)
          icount=1
c
          open(unit=40,file='newlowerenv.points',status='unknown')

          do 500 i=1,Nmodel
            cxmodel(i)=xmodel(i)
            cymodel(i)=ymodel(i)
            if(ymodel(i).lt.chimin+50.0d0)write(40,*)xmodel(i),ymodel(i)-chimin
500       continue
          close(40)
c
          icount=0
          dist=0.0d0
          do 501 i=2,Nmodel
            if(ymodel(i).gt.chimin+50.0d0)go to 502
c            dist=dsqrt((cxmodel(i)-cxmodel(i-1))**2+
c     @                (cymodel(i)-cymodel(i-1))**2)
c            if(dist.gt.1.0d-9)then
              icount=icount+1
              xmodel(icount)=cxmodel(i)
              ymodel(icount)=cymodel(i)
c            endif
501        continue
cc
502       Nmodel=icount
          icount=1



          xsmall=xmodel(1)
          write(*,888)xmodel(1),ymodel(1),Kmodel,filein
888       format(f16.9,2x,f10.3,2x,'N=',i7,1x,a40)
          write(33,*)xmodel(1),ymodel(1)
          pmax=-9.9d40
          pmin=9.0d40
          xkeep(1)=xmodel(1)
          ykeep(1)=ymodel(1)
          do 30 i=2,Nmodel
c              do 29 j=1,i
c                if(ymodel(i).gt.chimin+50.0d0)go to 30
                if(xmodel(i).gt.pmax)then
                  pmax=xmodel(i)
                  ymax=ymodel(i)
                  icount=icount+1
                  xkeep(icount)=pmax
                  ykeep(icount)=ymax
                endif
                if(xmodel(i).lt.pmin)then
                  pmin=xmodel(i)
                  ymin=ymodel(i)
                 icount=icount+1
                 xkeep(icount)=pmin
                 ykeep(icount)=ymin
                endif
c29            continue
c              write(33,*)pmin,ymin
c              write(33,*)pmax,ymax
 30       continue
c

          if(icount.gt.1)call sort2(icount,xkeep,ykeep)
          open(unit=33,file='newlowerenv.out',status='unknown')
          open(unit=40,file='newlowerenv.off',status='unknown')
c
          write(33,*)xkeep(1),ykeep(1)
          siglow=xsmall
          sighigh=xsmall
          chilow=99.99d0
          chihigh=99.99d0
          xmin=1.0d33
          xmax=-1.0d33
          ilowflag=0
          ihighflag=0
          do 50 i=2,icount
            if(xkeep(i).gt.xmax)xmax=xkeep(i)
            if(xkeep(i).lt.xmin)xmin=xkeep(i)
            dist=dsqrt((xkeep(i)-xkeep(i-1))**2+(ykeep(i)-ykeep(i-1))**2)
            if(dist.gt.1.0d-6)write(33,*)xkeep(i),ykeep(i)
            fred=ykeep(i)-chimin
            if(dist.gt.1.0d-6)write(40,*)xkeep(i),fred
            if(xkeep(i).lt.xsmall)then
              if((fred.gt.1.0d0).and.(fred.lt.chilow))then
                chilow=fred
                siglow=xkeep(i)
                ilowflag=99
              endif
            endif
            if(xkeep(i).gt.xsmall)then
              if((fred.gt.1.0d0).and.(fred.lt.chihigh))then
                chihigh=fred
                sighigh=xkeep(i)
                ihighflag=99
              endif
            endif

50        continue
c
          if(ilowflag.eq.0)siglow=xmin
          if(ihighflag.eq.0)sighigh=xmax
          write(*,889)siglow,chilow,sighigh,chihigh
          write(*,890)dabs(xsmall-siglow),dabs(xsmall-sighigh)
          if(dabs(xsmall-siglow).gt.dabs(xsmall-sighigh))then
            write(*,891)xsmall,dabs(xsmall-siglow),filein
          else
            write(*,891)xsmall,dabs(xsmall-sighigh),filein
          endif
889       format(4x,2(f16.9,2x,f7.3,2x))
890       format(4x,2(f16.9,2x))
891       format(f16.9,2x,f16.9,2x,a40/)
          close(33)
          close(40)
c
c    Now we have an array xinter and yinter containing the model resampled
c    at the same x-values as the data.  Perform the chi^2 fit.
c
          end
c
c    ===================================================
c
          subroutine getvariables(M,filein)
c
          implicit double precision(a-h,o-z)
          character*40 filein
          integer M
c
          if(M.eq.0)then   !case of no arguments
c
            write(*,*)'enter the file with parameter and chi^2 '
            read(*,100)filein
c
          endif
c
          if(M.eq.1)then        !case of only one argument
c
            call getarg(1,filein)     ! assign filein the value of the entered
                                      ! variable 
          endif
c
c
 100      format(a40)
          return
          end

c
c
c &&&&&&&&&&&&&&&&&&&&
c
      SUBROUTINE SORT2(N,RA,RB)
          implicit double precision(a-h,o-z)
      DIMENSION RA(N),RB(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
      END
c
c
c
c
      SUBROUTINE SORT3(N,RA,RB,rc)
          implicit double precision(a-h,o-z)
      DIMENSION RA(N),RB(N),rc(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
          rrc=rc(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          rrc=rc(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          rc(IR)=RC(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            rc(1)=rrc
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            rc(i)=rc(j)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
        rc(i)=rrc
      GO TO 10
      END
c
c  -----------------------------------------------------
c
          subroutine getchi(N,y,err,ydum,chisq)
c
          implicit double precision(a-h,o-z)
          dimension y(N),ydum(N),err(N)
          chisq=0.0d0
c
          do 10 i=1,N
            chisq=chisq+(y(i)-ydum(i))*(y(i)-ydum(i))/(err(i)*err(i))
 10       continue
          return
          end
c
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine offset(N,x,y,yoff,off)
c
c   Will offset the y values by offset and return new array yoff.
c
          implicit double precision(a-h,o-z)
          dimension x(N),y(N),yoff(N)
c
          do 10 i=1,N
            yoff(i)=y(i)+off
 10       continue
c
          return
          end
c
c  *********************************************************
c
          subroutine getmean(N,x,y,average)
c
          implicit double precision(a-h,o-z)
          dimension x(N),y(N)
c
          summ=0.0d0
          average=0.0d0
c
          do 10 i=1,N
            summ=summ+y(i)
c            write(*,*)x(i),y(i),summ
 10       continue
c 
          average=summ/dble(N)
c
c          write(*,*)'***'
c          write(*,*)average,N
 
          return
          end
c
c  &&&&&&&&&&&&&&&&&
c
         subroutine spline(x,y,n,yp1,ypn,y2)
c
c   November 12, 1999
c
c   This is a spline interpolation routine taken from NUMERICAL RECIPES.
c
         integer n,NMAX
         REAL*8 yp1,ypn,x(n),y(n),y2(n)
         parameter(NMAX=50000)
         integer i,k
         REAL*8 p,qn,sig,un,u(NMAX)

         if(yp1.gt..99d30)then
           y2(1)=0.0d0
           u(1)=0.0d0
         else
           y2(1)=-0.5d0
           u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
         endif
         do 11 i=2,n-1
           sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
           p=sig*y2(i-1)+2.0d0
           y2(i)=(sig-1.0d0)/p
           u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     #       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
 11      continue
         if(ypn.gt..99e30)then
           qn=0.0d0
           un=0.0d0
         else
           qn=0.5d0
           un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
         endif
         y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
         do 12 k=n-1,1,-1
           y2(k)=y2(k)*y2(k+1)+u(k)
 12      continue
         return
         end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
         subroutine splint(xa,ya,y2a,n,x,y)
c
c   November 12, 1999
c
c   This is a spline interpolation routine taken from Numerical Recipes.
c
         integer n
         real*8 x,y,xa(n),y2a(n),ya(n)
         integer k,khi,klo
         real*8 a,b,h
         klo=1
         khi=n
 1       if(khi-klo.gt.1)then
           k=(khi+klo)/2
           if(xa(k).gt.x)then
             khi=k
           else
             klo=k
           endif
           go to 1
         endif
         h=xa(khi)-xa(klo)
         if(h.eq.0.0d0)stop
         a=(xa(khi)-x)/h
         b=(x-xa(klo))/h
         y=a*ya(klo)+b*ya(khi)+
     $     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
         return
         end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
         subroutine countdot(line,Ndot)
c
         character*80 line
         integer Ndot
c
         Ndot=0
         k=lnblnk(line)
c
         do 10 i=1,k
           if(line(i:i).eq.'.')Ndot=Ndot+1
10       continue
c
         return
         end

       
