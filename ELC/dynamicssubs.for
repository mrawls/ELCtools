c          program test
cc
c          implicit double precision(a-h,o-z)
cc
c          PARAMETER (NDIM=30)
c          DIMENSION PP(NDIM),QQ(NDIM),QQBC(30),PPBC(30)
c          dimension rmass(10),posarray(400000,30),velarray(400000,30)
c          dimension zzq(400000,30,6),QQm2(3),PPm2(3),Esec(10000)
c          dimension odetime(400000),timeinterp(6),Eprim(10000)
c          dimension rIBCinp(10,6)
cc
c          pie=4.0d0*datan(1.0d0)
c          Nbody=3
c          h=4.0d0/400.0d0
c          tstart=-35.0d0
c          tend=3000.0d0    !1200
c
c          rmass(1)=0.7d0
c          rmass(2)=0.2d0
c          rmass(3)=0.000005d0 
c
c          rIBCinp(1,6)=-34.3437d0     !  Tconj
c          rIBCinp(1,1)=4.0d0    !  period
c          rIBCinp(1,2)=0.0d0
c          rIBCinp(1,3)=0.0d0
c          rIBCinp(1,4)=1.5769d0  !  finc
c          rIBCinp(1,5)=0.000d0     ! Omega
c
c          rIBCinp(2,1)=1600.0d0   ! period 3
c          rIBCinp(2,2)=-0.05d0     ! ecos for body 3
c          rIBCinp(2,3)=0.1d0    !    esin for body 3
c          rIBCinp(2,6)=-25.0d0    ! Tconj body 3
c          rIBCinp(2,4)=1.5769d0  !  finc for body 3  
c          rIBCinp(2,5)=0.0d0   !  Omega for boduy 3
c
c          call getIBC(Nbody,rmass,rIBCinp,QQ,PP,tstart)
c
c          write(*,101)QQ(1)
c          write(*,101)QQ(2)
c          write(*,101)QQ(3)
c          write(*,101)QQ(4)
c          write(*,101)QQ(5)
c          write(*,101)QQ(6)
c          write(*,101)QQ(7)
c          write(*,101)QQ(8)
c          write(*,101)QQ(9)
c 
c          write(*,*)' '
c          write(*,101)PP(1)
c          write(*,101)PP(2)
c          write(*,101)PP(3)
c          write(*,101)PP(4)
c          write(*,101)PP(5)
c          write(*,101)PP(6)
c          write(*,101)PP(7)
c          write(*,101)PP(8)
c          write(*,101)PP(9)
cc 
c 101      format(1(f19.15,2x))c
c
c          call solveorbit(Nbody,h,tstart,tend,Nstep,QQ,PP,rmass,
c     @         posarray,velarray,odetime,zzq,timeinterp)
c
c          call findprimaryeclipse(Nbody,posarray,velarray,odetime,
c     @       zzq,timeinterp,nstep,Nprim,Eprim,Nsec,Esec)
c
c          eb=1.0d0
c          write(*,*)Nprim,Nsec
c          do jj=1,Nprim
c            write(40,*)jj,Eprim(jj),eb
c          enddo
c          do jj=1,Nsec
c            write(41,*)jj,Esec(jj),eb
c          enddo
c
c          jlo=100
c          do i=1,20000
c            timein=dble(i)/1000.0d0
c            timein=timein+1.64d0
c            call LTTsky(odetime,zzq,Nstep,timeinterp,1,timein,
c     @        xout,yout,velout,ii)
c            call LTTsky(odetime,zzq,Nstep,timeinterp,2,timein,
c     @        xouts,youts,velouts,ii)
c            dist=dsqrt((xout-xouts)**2+(yout-youts)**2)
c            write(37,444)timein,dist/0.0046491d0
c          enddo

c          do i=1,20000
c            timein=dble(i)/1000.0d0
c            timein=timein+761.62d0
c            call LTTsky(odetime,zzq,Nstep,timeinterp,1,timein,
c     @        xout,yout,velout,ii)
c            call LTTsky(odetime,zzq,Nstep,timeinterp,2,timein,
c     @        xouts,youts,velouts,ii)
cc
c            dist=dsqrt((xout-xouts)**2+(yout-youts)**2)
c            write(38,444)timein,dist/0.0046491d0
c          enddo
cc
c
c          end
c
c  #!!!%(@^%#%*(&(#@#*$^@^$##*#@&(%)
c
          subroutine cartKep(Nbody,ibody,rmass,QQ,PP,Tref,timein,eccout,
     @      argout,rinclout,Omegaout,astarout,rmeanout,trueanomout,
     @      periodout)
c
c    converts Cartesian x,y,z to Keplerian
c
          implicit real*8 (a-h,o-z)

          dimension QQ(30),PP(30),rmass(10)
          dimension QQBC(30),PPBC(30)
c
c   Here is the Gravitational Constant in AU**2/solar_mass/day**2 units
c
          pie=4.0d0*datan(1.0d0)
          G=(0.01720209895d0)**2
c
          if(ibody.eq.2)G=G*(rmass(1)+rmass(2))
          if(ibody.eq.3)G=G*(rmass(1)+rmass(2)+rmass(3))
c
c   invert the move-to-barycenter operation
c
          call movetoastrocentric(Nbody,rmass,QQ,PP,QQBC,PPBC)
c
c    find the position and velocity of the body
c
          jndex=(ibody-1)*3
          x=QQBC(jndex+1)
          y=QQBC(jndex+2)
          z=QQBC(jndex+3)
          vx=PPBC(jndex+1)
          vy=PPBC(jndex+2)
          vz=PPBC(jndex+3)
c
c   find the direction of the angular momentum
c   vector
c
          rxv_x=y*vz-z*vy
          rxv_y=z*vx-x*vz
          rxv_z=x*vy-y*vx
          h=sqrt(rxv_x**2+rxv_y**2+rxv_z**2)
c
c   find the radius and velocity 
c
          r=sqrt(x**2+y**2+z**2)
          vs=vx**2+vy**2+vz**2
c
          rdotv=x*vx+y*vy+z*vz
          rdot=rdotv/r
c
          rhatx=x/r
          rhaty=y/r
          rhatz=z/r
c
          vxh_x=vy*rxv_z-vz*rxv_y
          vxh_y=vz*rxv_x-vx*rxv_z
          vxh_z=vx*rxv_y-vy*rxv_x
c
c          e_x=vxh_x/G-rhat_x
c          e_y=vxh_y/G-rhat_y
c          e_z=vxh_z/G-rhat_z
cc
c          write(*,*)sqrt(e_x**2+e_y**2+e_z**2)
c
          hs=h*h
          parm=hs/G
c
          rinclout=acos(rxv_z/h)
c
          if((rxv_x.ne.0.0d0).or.(rxv_y.ne.0.0d0))then
            Omegaout=atan2(rxv_x,-rxv_y)
          else
            Omegaout=0.0d0
          endif
c
          
          ecostrueanom=parm/r-1.0d0
          esintrueanom=rdot*h/G
          eccout=sqrt(ecostrueanom * ecostrueanom + 
     @         esintrueanom * esintrueanom)
c
          if((esintrueanom.ne.0.0d0).or.
     @         (ecostrueanom.ne.0.0d0))then
            trueanomout=atan2(esintrueanom,ecostrueanom)
          else
            trueanomout=0.0d0
          endif
c
          cosnode=cos(Omegaout)
          sinnode=sin(Omegaout)
c
          rcosu=x*cosnode+y*sinnode
          rsinu=(y*cosnode-x*sinnode)/cos(rinclout)
c
c   u is the argument of latitude
c          
          if((rsinu.ne.0.0d0).or.(rcosu.ne.0.0d0))then
            u=atan2(rsinu,rcosu)
          else
            u=0.0d0
          endif
c
          argout=u-trueanomout
          astarout=1.0d0/(2.0d0/r-vs/G)          
          if(argout.lt.0.0d0)argout=argout+2.0d0*pie
c
          eccanom=2.0d0*atan(sqrt((1.0d0-eccout)/
     @        (1.0d0+eccout))*tan(trueanomout/2.0d0))
          rmeanout=eccanom-eccout*sin(ecanom)
c
          if(astarout.gt.0.0d0)then
            rn=sqrt(G/astarout**3)
            tp=-rmeanout/rn
          endif
          periodout=abs(2.0d0*pie*tp/rmeanout)

          return
          end
c
c   !!%(*&^%^$**!)!)))!@@@@&@&%*$$
c
          subroutine getIBC(Nbody,rmass,rIBCinp,QQ,PP,tstart,itconj,
     @          isw28)

          implicit double precision(a-h,o-z)
c
          PARAMETER (NDIM=30)
          DIMENSION PP(NDIM),QQ(NDIM),QQBC(30),PPBC(30)
          dimension rmass(10)
          dimension QQm2(3),PPm2(3)
          dimension rIBCinp(10,6),tempmass(10)
          
          pie=4.0d0*datan(1.0d0)
          rm1=rmass(1)
          rm2=rmass(2)
          period=rIBCinp(1,1)
          ecos=rIBCinp(1,2)
          esin=rIBCinp(1,3)
          argper=datan2(esin,ecos)
          ecc=dsqrt(ecos*ecos+esin*esin)
          Tconj=rIBCinp(1,6)
          rmean=TrueAtoMeanA(pie/2.0d0-argper,ecc)
          if(isw28.eq.2)then
            rmean=TrueAtoMeanA(3.0d0*pie/2.0d0-argper,ecc)
          endif
          T0B=Tconj-period*rmean/(2.0d0*pie)
          rmeanbin=2.0d0*pie*(tstart-T0B)/period
          finc=rIBCinp(1,4)
          Omega=rIBCinp(1,5)
c
          call M2toM1centric(rM1,rM2,period,ecc,rmeanbin,
     @       argper,finc,Omega,QQm2,PPm2)
c
          QQ(1)=0.0d0
          QQ(2)=0.0d0
          QQ(3)=0.0d0
          QQ(4)=QQm2(1)
          QQ(5)=QQm2(2)
          QQ(6)=QQm2(3)
          PP(1)=0.0d0
          PP(2)=0.0d0
          PP(3)=0.0d0
          PP(4)=PPm2(1)
          PP(5)=PPm2(2)
          PP(6)=PPm2(3)
c
          QQBC(1)=0.0d0
          QQBC(2)=0.0d0
          QQBC(3)=0.0d0
          QQBC(4)=QQm2(1)
          QQBC(5)=QQm2(2)
          QQBC(6)=QQm2(3)
          PPBC(1)=0.0d0
          PPBC(2)=0.0d0
          PPBC(3)=0.0d0
          PPBC(4)=PPm2(1)
          PPBC(5)=PPm2(2)
          PPBC(6)=PPm2(3)
c
          tempmass(1)=rmass(1)
          tempmass(2)=rmass(2)
          tempmass(3)=0.0d0

          total=rM1+rM2
          do kk=2,Nbody-1
            period=rIBCinp(kk,1)
            ecos=rIBCinp(kk,2)
            esin=rIBCinp(kk,3)
            argper=datan2(esin,ecos)
            ecc=dsqrt(ecos*ecos+esin*esin)
            Tconj=rIBCinp(kk,6)
            rmean=TrueAtoMeanA(pie/2.0d0-argper,ecc)
            if(itconj.eq.2)then
              rmean=TrueAtoMeanA(3.0d0*pie/2.0d0-argper,ecc)
            endif
            T0B=Tconj-period*rmean/(2.0d0*pie)
            rmeanbin=2.0d0*pie*(tstart-T0B)/period
            finc=rIBCinp(kk,4)
            Omega=rIBCinp(kk,5)
            planetm=rmass(kk+1)
c
            call M2toM1centric(total,planetm,period,ecc,rmeanbin,
     @         argper,finc,Omega,QQm2,PPm2)
            total=total+planetm
c
            QQBC(3*kk+1)=QQm2(1)
            QQBC(3*kk+2)=QQm2(2)
            QQBC(3*kk+3)=QQm2(3)
            PPBC(3*kk+1)=PPm2(1)
            PPBC(3*kk+2)=PPm2(2)
            PPBC(3*kk+3)=PPm2(3)
          enddo
c
 100      format(3(1pe23.15,2x))
          call newmovetobarycenter(Nbody,rmass,QQBC,PPBC,QQ,PP)
c
          return
          end
c
c   !@#$%^&*&^%$#@#$%^%$#@!
c
          subroutine findprimaryeclipse(Nbody,posarray,velarray,odetime,
     @       zzq,timeinterp,nstep,Nprim,Eprim,Nsec,Esec,Ndyn,Distprim,
     @       Distsec,reff1,reff2,separ,ibody1,ibody2,Nmaxeclipse,
     @       durprim1,dursec1,durprim2,dursec2)
c
          implicit double precision (a-h,o-z)

          dimension posarray(Ndyn,30),velarray(Ndyn,30)
          dimension zzq(6,60,Ndyn)
          dimension odetime(Ndyn),timeinterp(6)
          dimension xx(15),dd(15),Distprim(Nmaxeclipse)
          dimension Distsec(Nmaxeclipse)
          dimension Eprim(Nmaxeclipse),Esec(Nmaxeclipse)
          dimension durprim1(Nmaxeclipse),dursec1(Nmaxeclipse)
          dimension durprim2(Nmaxeclipse),dursec2(Nmaxeclipse)
c
          rsunperAU=0.00465116d0     !solar radius in AU
          AC=1.49597870691D11/86400.0d0  ! AU/day --> m/sec
          speed=299792458.0d0/AC    !speed of light in AU/day
c
          jlo=100
          Neclipse1=0
          Neclipse2=0
          ttiny=1.0d-14
          sumdist=(reff1+reff2)*separ
          if(sumdist.le.0.0d0)sumdist=ttiny
          do i=1,Nstep-1
            if(posarray(i,ibody1*3).lt.posarray(i,ibody2*3))then !body1 in back
              timein=odetime(i)
              call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,timein,
     @        xp,yp,vp,jlo,Ndyn)
              call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,timein,
     @        xs,ys,vs,jlo,Ndyn)
c
              dx=xp-xs  !dxp-dxs   !posarray(i,1)-posarray(i,4)
              dy=yp-ys   !dyp-dys   !posarray(i,2)-posarray(i,5)
              vx=velarray(i,3*ibody1-2)-velarray(i,3*ibody2-2)
              vy=velarray(i,3*ibody1-1)-velarray(i,3*ibody2-1)
              deriv=dx*vx+dy*vy
              xx(1)=odetime(i)
              dd(1)=deriv
c
              timein=odetime(i+1)
              call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,timein,
     @        xp,yp,vp,jlo,Ndyn)
              call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,timein,
     @        xs,ys,vs,jlo,Ndyn)
              dx1=xp-xs   !dxp-dxs  !posarray(i+1,1)-posarray(i+1,4)
              dy1=yp-ys   !dyp-dys  !posarray(i+1,2)-posarray(i+1,5)
              dist1=dsqrt(dx*dx+dy*dy)
              vx1=velarray(i+1,3*ibody1-2)-velarray(i+1,3*ibody2-2)
              vy1=velarray(i+1,3*ibody1-1)-velarray(i+1,3*ibody2-1)
              deriv1=dx1*vx1+dy1*vy1
              xx(2)=odetime(i+1)
              dd(2)=deriv1
c
              if((deriv.lt.0.0d0).and.(deriv1.ge.0.0d0))then  !sign change    
                do jj=3,15
                  diff=dd(jj-1)-dd(jj-2)
                  if(diff.eq.0.0d0)then
                    xx(jj)=xx(jj-1)
                  else
                    xx(jj)=xx(jj-1)-dd(jj-1)*((xx(jj-1)-xx(jj-2))
     @                   /(diff))
                  endif
                  if(dabs(xx(jj)-xx(jj-1)).lt.ttiny)go to 98
                  timemid=xx(jj)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,
     @              timemid,xp,yp,vp,jlo,Ndyn)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,
     @              timemid,xs,ys,vs,jlo,Ndyn)
                  call divdiff(odetime,zzq,Nstep,timeinterp,3*ibody1-2,
     @              timemid,xout1,veloutx1,jlo,Ndyn)
                  call divdiff(odetime,zzq,Nstep,timeinterp,3*ibody1-1,
     @              timemid,yout1,velouty1,jlo,Ndyn)
c                  call divdiff(odetime,zzq,Nstep,timeinterp,ibody1,
c     @              timemid,zout1,veloutz1,jlo,Ndyn)
                  call divdiff(odetime,zzq,Nstep,timeinterp,3*ibody2-2,
     @              timemid,xout2,veloutx2,jlo,Ndyn)
                  call divdiff(odetime,zzq,Nstep,timeinterp,3*ibody2-1,
     @              timemid,yout2,velouty2,jlo,Ndyn)
c                  call divdiff(odetime,zzq,Nstep,timeinterp,ibody2,
c     @              timemid,zout2,veloutz2,jlo,Ndyn)
                  dx=xp-xs  !dxp-dxs   !xout1-xout2
                  dy=yp-ys  !dyp-dys   !yout1-yout2
                  vx=veloutx1-veloutx2
                  vy=velouty1-velouty2
                  deriv1=dx*vx+dy*vy
                  dd(jj)=deriv1
                enddo
 98             Neclipse1=Neclipse1+1
                Eprim(Neclipse1)=timemid
                Distprim(Neclipse1)=dsqrt(dx*dx+dy*dy)/sumdist
                Nprim=Neclipse1
                if(Neclipse1.gt.Nmaxeclipse)then
                  write(*,*)'Nmaxeclipse is too small'
                  stop
                endif
c
c   find ingress and egress points
c
                tmidsave=timemid
                aaa=timemid
                bbb=timemid-10.0d0
                do kback=-1,-150,-1
                  timein=odetime(i+kback)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,
     @                timein,xp,yp,vp,jlo,Ndyn)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,
     @                timein,xs,ys,vs,jlo,Ndyn)
                  ddist=sqrt((xp-xs)**2+(yp-ys)**2)
                  ffunc=ddist-sumdist
                  if(ffunc.lt.0.0d0)then
                    aaa=timein
                  else
                    bbb=timein
                    go to 66
                  endif
                enddo
c
 66             do kbis=1,20
                  tmid=0.5d0*(bbb+aaa)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,
     @                 tmid,xp,yp,vp,jlo,Ndyn)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,
     @                 tmid,xs,ys,vs,jlo,Ndyn)
                  ddist=sqrt((xp-xs)**2+(yp-ys)**2)
                  ffunc=ddist-sumdist
                  if(ffunc.lt.0.0d0)then
                    aaa=tmid
                  else
                    bbb=tmid
                  endif
                enddo
c
                durprim1(Neclipse1)=tmid
c
                aaa=tmidsave
                bbb=tmidsave+10.0d0
                do kback=1,150
                  timein=odetime(i+kback)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,
     @                timein,xp,yp,vp,jlo,Ndyn)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,
     @                timein,xs,ys,vs,jlo,Ndyn)
                  ddist=sqrt((xp-xs)**2+(yp-ys)**2)
                  ffunc=ddist-sumdist
                  if(ffunc.lt.0.0d0)then
                    aaa=timein
                  else
                    bbb=timein
                    go to 67
                  endif
                enddo
c
 67             do kbis=1,20
                  tmid=0.5d0*(bbb+aaa)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,
     @                 tmid,xp,yp,vp,jlo,Ndyn)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,
     @                 tmid,xs,ys,vs,jlo,Ndyn)
                  ddist=sqrt((xp-xs)**2+(yp-ys)**2)
                  ffunc=ddist-sumdist
                  if(ffunc.lt.0.0d0)then
                    aaa=tmid
                  else
                    bbb=tmid
                  endif
                enddo
c
                durprim2(Neclipse1)=tmid
c
              endif
            endif
c
            if(posarray(i,ibody2*3).lt.posarray(i,ibody1*3))then !body2 in back
              timein=odetime(i)
              call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,timein,
     @           xp,yp,vp,jlo,Ndyn)
              call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,timein,
     @           xs,ys,vs,jlo,Ndyn)
              dx=xp-xs  !posarray(i,1)-posarray(i,4)
              dy=yp-ys  !posarray(i,2)-posarray(i,5)
              dist=dsqrt(dx*dx+dy*dy)
              vx=velarray(i,3*ibody1-2)-velarray(i,3*ibody2-2)
              vy=velarray(i,3*ibody1-1)-velarray(i,3*ibody2-1)
              deriv=dx*vx+dy*vy
              xx(1)=odetime(i)
              dd(1)=deriv
c
              timein=odetime(i+1)
              call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,timein,
     @        xp,yp,vp,jlo,Ndyn)
              call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,timein,
     @        xs,ys,vs,jlo,Ndyn)
              dx1=xp-xs   !posarray(i+1,1)-posarray(i+1,4)
              dy1=yp-ys   !posarray(i+1,2)-posarray(i+1,5)
              dist1=dsqrt(dx*dx+dy*dy)
              vx1=velarray(i+1,3*ibody1-2)-velarray(i+1,3*ibody2-2)
              vy1=velarray(i+1,3*ibody1-1)-velarray(i+1,3*ibody2-1)
              deriv1=dx1*vx1+dy1*vy1
              xx(2)=odetime(i+1)
              dd(2)=deriv1
c
              if((deriv.lt.0.0d0).and.(deriv1.ge.0.0d0))then  !sign change    
                do jj=3,15
                  diff=dd(jj-1)-dd(jj-2)
                  if(diff.eq.0.0d0)then
                    xx(jj)=xx(jj-1)
                  else
                    xx(jj)=xx(jj-1)-dd(jj-1)*((xx(jj-1)-xx(jj-2))
     @                   /(diff))
                  endif
                  if(dabs(xx(jj)-xx(jj-1)).lt.ttiny)go to 99
                  timemid=xx(jj)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,
     @                timemid,xp,yp,vp,jlo,Ndyn)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,
     @                timemid,xs,ys,vs,jlo,Ndyn)
                  call divdiff(odetime,zzq,Nstep,timeinterp,3*ibody1-2,
     @                timemid,xout1,veloutx1,jlo,Ndyn)
                  call divdiff(odetime,zzq,Nstep,timeinterp,3*ibody1-1,
     @                timemid,yout1,velouty1,jlo,Ndyn)

c                  call divdiff(odetime,zzq,Nstep,timeinterp,3,timemid,
c     @              zout1,veloutz1,jlo,Ndyn)

                  call divdiff(odetime,zzq,Nstep,timeinterp,3*ibody2-2,
     @              timemid,xout2,veloutx2,jlo,Ndyn)
                  call divdiff(odetime,zzq,Nstep,timeinterp,3*ibody2-1,
     @              timemid,yout2,velouty2,jlo,Ndyn)

c                  call divdiff(odetime,zzq,Nstep,timeinterp,6,timemid,
c     @              zout2,veloutz2,jlo,Ndyn)

                  dx=xp-xs  !xout1-xout2
                  dy=yp-ys  !yout1-yout2
                  dist=dx*dx+dy*dy
                  vx=veloutx1-veloutx2
                  vy=velouty1-velouty2
                  deriv1=dx*vx+dy*vy
                  dd(jj)=deriv1
                  if(deriv1.lt.0.0d0)then
                    aa=timemid
                  else
                    bb=timemid
                  endif
                enddo
 99             Neclipse2=Neclipse2+1
                Esec(Neclipse2)=timemid
                Distsec(Neclipse2)=dsqrt(dx*dx+dy*dy)/sumdist
                Nsec=Neclipse2
                if(Neclipse2.gt.Nmaxeclipse)then
                  write(*,*)'Nmaxeclipse is too small'
                  stop
                endif
c
c   find ingress and egress points
c
                tmidsave=timemid
                aaa=timemid
                bbb=timemid-10.0d0
                do kback=-1,-150,-1
                  jdx=i+kback
                  if(jdx.lt.1)jdx=1
                  timein=odetime(jdx)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,
     @                timein,xp,yp,vp,jlo,Ndyn)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,
     @                timein,xs,ys,vs,jlo,Ndyn)
                  ddist=sqrt((xp-xs)**2+(yp-ys)**2)
                  ffunc=ddist-sumdist
                  if(ffunc.lt.0.0d0)then
                    aaa=timein
                  else
                    bbb=timein
                    go to 366
                  endif
                enddo
c
366             do kbis=1,20
                  tmid=0.5d0*(bbb+aaa)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,
     @                 tmid,xp,yp,vp,jlo,Ndyn)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,
     @                 tmid,xs,ys,vs,jlo,Ndyn)
                  ddist=sqrt((xp-xs)**2+(yp-ys)**2)
                  ffunc=ddist-sumdist
                  if(ffunc.lt.0.0d0)then
                    aaa=tmid
                  else
                    bbb=tmid
                  endif
                enddo
c
                dursec1(Neclipse2)=tmid
c
                aaa=tmidsave
                bbb=tmidsave+10.0d0
                do kback=1,150
                  timein=odetime(i+kback)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,
     @                timein,xp,yp,vp,jlo,Ndyn)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,
     @                timein,xs,ys,vs,jlo,Ndyn)
                  ddist=sqrt((xp-xs)**2+(yp-ys)**2)
                  ffunc=ddist-sumdist
                  if(ffunc.lt.0.0d0)then
                    aaa=timein
                  else
                    bbb=timein
                    go to 367
                  endif
                enddo
c
367             do kbis=1,20
                  tmid=0.5d0*(bbb+aaa)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody1,
     @                 tmid,xp,yp,vp,jlo,Ndyn)
                  call LTTsky(odetime,zzq,Nstep,timeinterp,ibody2,
     @                 tmid,xs,ys,vs,jlo,Ndyn)
                  ddist=sqrt((xp-xs)**2+(yp-ys)**2)
                  ffunc=ddist-sumdist
                  if(ffunc.lt.0.0d0)then
                    aaa=tmid
                  else
                    bbb=tmid
                  endif
                enddo
c
                dursec2(Neclipse2)=tmid
c
              endif
            endif
          enddo
c
c 100      format(f13.6,2x,f14.11)

          return
          end
c
c   #!!!%(@^%#%*(&(@#*$^@^$##*#@&(%)
c
          subroutine solveorbit(Nbody,h,tstart,tend,Nstep,QQ,PP,rmass,
     @       posarray,velarray,odetime,zzq,timeinterp,Ndyn,iGR,tideparm)
c
c   This will integrate the equations of motion using a 12th order   
c   Gaussian Runga Kutta scheme devised by E. Hairer.
c
c   Nbody = number of bodies, limit 10.  Body 1 is the primary of the binary
c   and body 2 is the secondary.
c
c   h = the stepsize in days.  It should be about P_binary/400 or smaller
c
c   tstart,tend =  The start and end times of the integration, in days
c
c   Nstep = the number of integration steps, which is computed below
c 
c   QQ = the array with the positions of each body in AU:
c        QQ(1) = x-coordinate of primary
c        QQ(2) = y-coordinate of primary
c        QQ(3) = z-coordinate of primary
c        QQ(4) = x-coordinate of secondary
c
c   PP = the array with the velocities, in AU/day
c
c   rmass = the array with the masses in solar masses
c
c   posarray, velarray = arrays of dimension (Nstep,30) with the positions
c     and velocities of all of the bodies
c
c   odetime = array of size Nstep with the times at each point
c
c   zzq = array of dimension (6,30,nstep) interior node points from
c    the ode solver.  These will be used to interpolate intermediate
c    values
c
c   timeinterp = array of dimension 6 that has the times of the interior
c     node points
c
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NDIM=30)
      dimension oldPP(NDIM),oldQQ(NDIM)
      DIMENSION PP(NDIM),QQ(NDIM),IPAR(20),RPAR(20),DPAR(10)
      dimension rmass(10),posarray(Ndyn,30),velarray(Ndyn,30)
      dimension zzq(6,60,Ndyn),DY(60),oldzzq(6,60,Ndyn)
      dimension smallzzq(6,60,2),smallposarray(2,30)
      dimension smallodetime(2),smallvelarray(2,30)
      dimension odetime(Ndyn),timeinterp(6),METH(10)
      dimension tideparm(10)
      EXTERNAL EQUA,SOLFIX,DEQUA
c
C --- CHOOSE THE PROBLEM
C  IPROB = 1 : KEPLER PROBLEM, ECCENTRICITY IN RPAR(1)
C  IPROB = 2 : HARMONIC OSCILLATOR
C  IPROB = 3 : PENDULUM
C  IPROB = 4 : OUTER SOLAR SYSTEM
      IPROB=4
      IPAR(11)=IPROB
c      IF (IPROB.EQ.1) RPAR(11)=0.6D0
c      CALL PDATA(X,XEND,NDIM,N,Q,P,QEX,PEX,RPAR,IPAR)
C --- CHOOSE THE THE METHOD
C --- GAUSS METHOD OF ORDER 2*METH
      MMETH=6
c      WRITE (6,*) '     METHOD, PROBLEM  ',METH,IPROB
c      H=4.0d0

      x=tstart
      xend=tend
      NSTEP=int((XEND-X)/H)
      if(Nstep.gt.Ndyn)then
        write(*,*)'Error:  Number of integration steps > Ndyn'
        stop
      endif

      DO I=1,10
        RPAR(I)=0.0D0
        IPAR(I)=0
      END DO
      IPAR(12)=0
      IOUT=0

      if(iGR.le.0)then
        CALL GNI_IRK2(Nbody*3,EQUA,NSTEP,X,PP,QQ,XEND,
     &     MMETH,SOLFIX,IOUT,RPAR,IPAR,rmass,nbody,posarray,velarray,
     @     odetime,zzq,timeinterp,Ndyn,h,iGR)
      endif
c
      Ntimes=Nstep
      if(iGR.ge.1)then
c
        x=tstart
        xend=tend
        Ntimes=Nstep
        Ndynsmall=2
        do kkk=1,Ntimes
          xend=x+h
          Nstep=1
          xsave=x
          xendsave=xend
          do jj=1,NDIM
            oldPP(jj)=PP(jj)
            oldQQ(jj)=QQ(jj)
          enddo
          CALL GNI_IRK2(Nbody*3,EQUA,NSTEP,X,oldPP,oldQQ,XEND,
     &       MMETH,SOLFIX,IOUT,RPAR,IPAR,rmass,nbody,
     @       smallposarray,smallvelarray,
     @       smallodetime,smallzzq,timeinterp,Ndynsmall,h,iGR)
c
           do i=1,6*Nbody
             do j=1,6
               oldzzq(j,i,1)=smallzzq(j,i,1)
             enddo
           enddo
c
           x=xsave
           xend=xendsave
           NEQU=Nbody*6
           NN=60
           NPER=100
           METH(1)=6
           METH(2)=0
           METH(3)=0
           METH(4)=0
           METH(5)=0
           METH(6)=0
           METH(7)=0
           METH(8)=0
           METH(9)=0
           METH(10)=0
           DX=X
           DXEND=XEND
           DPAR(1)=2.0D-20
           IOUT=0
           do kk=1,3*Nbody
             DY(kk)=QQ(kk)
             DY(kk+3*Nbody)=PP(kk)
           enddo
c
           CALL GRKAAD(NEQU,DEQUA,Nstep,DX,DY,DXEND,METH,IOUT,DPAR,
     @       rmass,nbody,smallposarray,smallvelarray,
     @       smallodetime,smallzzq,timeinterp,
     @       Ndynsmall,h,iGR,oldzzq,tideparm)

           do kk=1,3*Nbody
             QQ(kk)=DY(kk)
             PP(kk)=DY(kk+3*Nbody)
           enddo
c
           do i=1,3*Nbody
             do j=1,6
               zzq(j,i,kkk)=smallzzq(j,i,1)
             enddo
           enddo
c
           x=x+h
c
           odetime(kkk)=smallodetime(1)
c
           do i=1,3*Nbody
             posarray(kkk,i)=smallposarray(1,i)           
             velarray(kkk,i)=smallvelarray(1,i)           
           enddo
        enddo   !kkk over time steps
      endif
c
      Nstep=Ntimes
      call newtondiff(odetime,zzq,Nstep,timeinterp,Nbody,timein,
     @        xout,velout,ii,Ndyn)
c
      return
      END
C
c  @!#$%!^!&!*!(!)!_!)!(!*!*&!&!^!^%!!!#!#!!
c
          subroutine goback(Nbody,h,tstart,tend,Nstep,QQ,PP,rmass,
     @       posarray,velarray,odetime,zzq,timeinterp,Ndyn,iGR,tideparm)
c
c   This will integrate the equations of motion using a 12th order   
c   Gaussian Runga Kutta scheme devised by E. Hairer.
cc
c   It will be used in cases where tstart is no equal to tref
c
c   Nbody = number of bodies, limit 10.  Body 1 is the primary of the binary
c   and body 2 is the secondary.
c
c   h = the stepsize in days.  It should be about P_binary/400 or smaller
c
c   tstart,tend =  The start and end times of the integration, in days
c
c   Nstep = the number of integration steps, which is computed below
c 
c   QQ = the array with the positions of each body in AU:
c        QQ(1) = x-coordinate of primary
c        QQ(2) = y-coordinate of primary
c        QQ(3) = z-coordinate of primary
c        QQ(4) = x-coordinate of secondary
c
c   PP = the array with the velocities, in AU/day
c
c   rmass = the array with the masses in solar masses
c
c   posarray, velarray = arrays of dimension (Nstep,30) with the positions
c     and velocities of all of the bodies
c
c   odetime = array of size Nstep with the times at each point
c
c   zzq = array of dimension (6,30,Nstep) interior node points from
c    the ode solver.  These will be used to interpolate intermediate
c    values
c
c   timeinterp = array of dimension 6 that has the times of the interior
c     node points
c
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NDIM=30)
      dimension oldPP(NDIM),oldQQ(Ndim)
      DIMENSION PP(NDIM),QQ(NDIM),IPAR(20),RPAR(20),DPAR(10)
      dimension rmass(10),posarray(Ndyn,30),velarray(Ndyn,30)
      dimension zzq(6,60,Ndyn),meth(10),DY(60),oldzzq(6,60,Ndyn)
      dimension odetime(Ndyn),timeinterp(6)
      dimension smallzzq(6,60,2),smallposarray(2,30)
      dimension smallodetime(2),smallvelarray(2,30)
      dimension tideparm(10)
      EXTERNAL EQUA,SOLFIX,DEQUA
c
C --- CHOOSE THE PROBLEM
C  IPROB = 1 : KEPLER PROBLEM, ECCENTRICITY IN RPAR(1)
C  IPROB = 2 : HARMONIC OSCILLATOR
C  IPROB = 3 : PENDULUM
C  IPROB = 4 : OUTER SOLAR SYSTEM
      IPROB=4
      IPAR(11)=IPROB
c      IF (IPROB.EQ.1) RPAR(11)=0.6D0
c      CALL PDATA(X,XEND,NDIM,N,Q,P,QEX,PEX,RPAR,IPAR)
C --- CHOOSE THE THE METHOD
C --- GAUSS METHOD OF ORDER 2*METH
      MMETH=6
c      WRITE (6,*) '     METHOD, PROBLEM  ',METH,IPROB
c      H=4.0d0

      x=tstart
      xend=tend
      NSTEP=int((XEND-X)/H)
      Nstep=abs(nstep)
      if(Nstep.gt.Ndyn)then
        write(*,*)'Error:  Number of integration steps > Ndyn'
        stop
      endif

      DO I=1,10
        RPAR(I)=0.0D0
        IPAR(I)=0
      END DO
      IPAR(12)=0
      IOUT=0

      if(iGR.le.0)then
        CALL GNI_IRK2(Nbody*3,EQUA,NSTEP,X,PP,QQ,XEND,
     &     MMETH,SOLFIX,IOUT,RPAR,IPAR,rmass,nbody,posarray,velarray,
     @     odetime,zzq,timeinterp,Ndyn,h,iGR)
      endif
c
c      if(iGR.eq.1)then
c
c        do kk=1,Ndim
c          PP(kk)=-PP(kk)
c        enddo
c
c        x=tstart
c        xend=x+abs(tstart-tend)
c        Ntimes=Nstep
c        Ndynsmall=2
c        h=abs(h)
c        do kkk=1,Ntimes
c          xend=x+h
c          Nstep=1
c          xsave=x
c          xendsave=xend
c          do jj=1,NDIM
c            oldPP(jj)=PP(jj)
c            oldQQ(jj)=QQ(jj)
c          enddo
cc
c          CALL GNI_IRK2(Nbody*3,EQUA,NSTEP,X,oldPP,oldQQ,XEND,
c     &       MMETH,SOLFIX,IOUT,RPAR,IPAR,rmass,nbody,
c     @       smallposarray,smallvelarray,
c     @       smallodetime,smallzzq,timeinterp,Ndynsmall,h,iGR)
cc
c           do i=1,3*Nbody
c             do j=1,6
c               oldzzq(j,i,1)=smallzzq(j,i,1)
c             enddo
c           enddo
cc
c           x=xsave
c           xend=xendsave
c           NEQU=Nbody*6
c           NN=60
c           NPER=100
c           METH(1)=6
c           METH(2)=0
c           METH(3)=0
c           METH(4)=0
c           METH(5)=0
c           METH(6)=0
c           METH(7)=0
c           METH(8)=0
c           METH(9)=0
c           METH(10)=0
c           DX=X
c           DXEND=XEND
c           DPAR(1)=2.0D-20
c           IOUT=0
c           do kk=1,3*Nbody
c             DY(kk)=QQ(kk)
c             DY(kk+3*Nbody)=PP(kk)
c           enddo
cc
c           CALL GRKAAD(NEQU,DEQUA,Nstep,DX,DY,DXEND,METH,IOUT,DPAR,
c     @       rmass,nbody,smallposarray,smallvelarray,
c     @       smallodetime,smallzzq,timeinterp,
c     @       Ndynsmall,h,iGR,oldzzq)
c
c           do kk=1,3*Nbody
c             QQ(kk)=DY(kk)
c             PP(kk)=DY(kk+3*Nbody)
c           enddo
cc
c           do i=1,3*Nbody
c             do j=1,6
c               zzq(j,i,kkk)=smallzzq(j,i,1)
c             enddo
c           enddo
cc
c           x=x+h
cc 
c           odetime(kkk)=smallodetime(1)
cc
c           do i=1,3*Nbody
c             posarray(kkk,i)=smallposarray(1,i)           
c             velarray(kkk,i)=smallvelarray(1,i)           
c           enddo
c        enddo   !kkk over time steps
c        do kk=1,3*Nbody
c          PP(kk)=-PP(kk)
c        enddo
c      endif
c
      if(iGR.ge.1)then
         NEQU=Nbody*6
         NN=60
         NPER=100
         METH(1)=6
         METH(2)=0
         METH(3)=0
         METH(4)=0
         METH(5)=0
         METH(6)=0
         METH(7)=0
         METH(8)=0
         METH(9)=0
         METH(10)=0
         DX=X
         DXEND=XEND
         DPAR(1)= 2.0d-16 !2.0D-20
         IOUT=0
         do kk=1,3*Nbody
           DY(kk)=QQ(kk)
           DY(kk+3*Nbody)=PP(kk)
         enddo
         CALL oldGRKAAD(NEQU,DEQUA,Nstep,DX,DY,DXEND,METH,IOUT,DPAR,
     @     rmass,nbody,posarray,velarray,odetime,zzq,timeinterp,Ndyn,h,
     @     iGR,tideparm)
         do kk=1,3*Nbody
           QQ(kk)=DY(kk)
           PP(kk)=DY(kk+3*Nbody)
         enddo
      endif

      return
      END
c
c
c      return
c      END
C
c   #!!!%(@^%#%*(&(#@#
c
      SUBROUTINE SOLFIX (NR,XOLD,X,P,Q,N,RPAR,IPAR,rmass,nbody)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(N),P(N)
      DIMENSION IPAR(20),RPAR(20)
      dimension rmass(10)
      CALL HAMILTON (N,Q,P,HAMIL,IPAR,rmass,nbody)
      IF (NR.EQ.0) THEN
        RPAR(12)=HAMIL
        RPAR(13)=0.0D0
      ELSE
        RPAR(13)=MAX(RPAR(13),ABS(RPAR(12)-HAMIL))
      END IF
      RETURN
      END
C
C  IPROB = 1 : KEPLER PROBLEM, ECCENTRICITY IN RPAR(1)
C  IPROB = 2 : HARMONIC OSCILLATOR
C  IPROB = 3 : PENDULUM
C  IPROB = 4 : OUTER SOLAR SYSTEM
C
      SUBROUTINE PDATA(X,XEND,NDIM,N,Q,P,QEX,PEX,RPAR,IPAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q(NDIM),P(NDIM),QEX(NDIM),PEX(NDIM)
      DIMENSION IPAR(*),RPAR(*)
      IPROB=IPAR(11)
      IF (IPROB.EQ.1) THEN
        N=2
        PI=4*ATAN(1.D0)
        ECCENT=RPAR(11)
        X=0.D0
        XEND=2*PI
        Q(1)=1-ECCENT
        Q(2)=0.d0
        P(1)=0.d0
        P(2)=SQRT((1+ECCENT)/(1-ECCENT))
        QEX(1)=1-ECCENT
        QEX(2)=0.d0
        PEX(1)=0.d0
        PEX(2)=SQRT((1+ECCENT)/(1-ECCENT))
      END IF
      IF (IPROB.EQ.2) THEN
        N=1
        PI=4*ATAN(1.D0)
        X=0.D0
        XEND=2*PI
        Q(1)=0.0D0
        P(1)=1.0D0
        QEX(1)=0.0D0
        PEX(1)=1.0D0
      END IF
      IF (IPROB.EQ.3) THEN
        N=1
        PI=4*ATAN(1.D0)
        X=0.D0
        XEND=2*PI
        Q(1)=0.0D0
        P(1)=1.0D0
        QEX(1)=-0.443944662290259D0
        PEX(1)= 0.897846803312944D0
      END IF
      IF (IPROB.EQ.4) THEN
        N=18
        X=0.D0
        XEND=500000.D0
        Q(1)=-3.5023653D0
        Q(2)=-3.8169847D0
        Q(3)=-1.5507963D0
        Q(4)=9.0755314D0
        Q(5)=-3.0458353D0
        Q(6)=-1.6483708D0
        Q(7)=8.3101420D0
        Q(8)=-16.2901086D0
        Q(9)=-7.2521278D0
        Q(10)=11.4707666D0
        Q(11)=-25.7294829D0
        Q(12)=-10.8169456D0
        Q(13)=-15.5387357D0
        Q(14)=-25.2225594D0
        Q(15)=-3.1902382D0
        Q(16)=0.0D0
        Q(17)=0.0D0
        Q(18)=0.0D0
        P(1)=0.00565429D0
        P(2)=-0.00412490D0
        P(3)=-0.00190589D0
        P(4)=0.00168318D0
        P(5)=0.00483525D0
        P(6)=0.00192462D0
        P(7)=0.00354178D0
        P(8)=0.00137102D0
        P(9)=0.00055029D0
        P(10)=0.00288930D0
        P(11)=0.00114527D0
        P(12)=0.00039677D0
        P(13)=0.00276725D0
        P(14)=-0.00170702D0
        P(15)=-0.00136504D0
        P(16)=0.0D0
        P(17)=0.0D0
        P(18)=0.0D0
        QEX( 1)=   0.7766584086800482D+01
        QEX( 2)=   0.2531065754551048D+00
        QEX( 3)=  -0.9410571402013185D-01
        QEX( 4)=  -0.5564967162844037D+01
        QEX( 5)=   0.1674849740822012D+01
        QEX( 6)=   0.9767232069533176D+00
        QEX( 7)=   0.1963899572895227D+02
        QEX( 8)=   0.8958504552286460D+01
        QEX( 9)=   0.3611839157057347D+01
        QEX(10)=   0.2493570870305177D+02
        QEX(11)=   0.1769518676153705D+02
        QEX(12)=   0.6583785164549242D+01
        QEX(13)=   0.3178592511375764D+02
        QEX(14)=   0.3863618958160644D+02
        QEX(15)=   0.3192794169732889D+01
        QEX(16)=   0.3084118473380683D+01
        QEX(17)=  -0.1227726356581642D+01
        QEX(18)=  -0.6162537634647217D+00
        PEX( 1)=  -0.2495503201917009D-02
        PEX( 2)=   0.6896467194473328D-02
        PEX( 3)=   0.3007950247474123D-02
        PEX( 4)=  -0.2255335935351989D-02
        PEX( 5)=  -0.4905913854771086D-02
        PEX( 6)=  -0.1938473641716708D-02
        PEX( 7)=  -0.2186170231167942D-02
        PEX( 8)=   0.2817177012110666D-02
        PEX( 9)=   0.1262882639181183D-02
        PEX(10)=  -0.2148728705895163D-02
        PEX(11)=   0.2128650077635786D-02
        PEX(12)=   0.9248501411662923D-03
        PEX(13)=  -0.1675173186229401D-02
        PEX(14)=   0.1011833320388655D-02
        PEX(15)=   0.8231800038576520D-03
        PEX(16)=   0.9417379703028725D-05
        PEX(17)=  -0.7855256238249194D-05
        PEX(18)=  -0.3646926313230521D-05
      END IF
      RETURN
      END
c
      SUBROUTINE EQUA(N,X,Q,F,RPAR,IPAR,rmass,nbody)
      IMPLICIT REAL*8 (A-H,O-Z) 
      DOUBLE PRECISION  D(10,10),K,rmass(10)
      DIMENSION Q(N),F(N)
      DIMENSION IPAR(*),RPAR(*)
c      IPROB=IPAR(11)
      IPAR(12)=IPAR(12)+1
      iprob=4
      IF (IPROB.EQ.4) THEN
        K=(0.01720209895d0)**2
        DO I=1,Nbody-1   !5
          I1=3*(I-1)+1
          DO J=I+1,Nbody   !6
            J1=3*(J-1)+1
            D(I,J)=(SQRT((Q(I1)-Q(J1))**2+(Q(I1+1)-Q(J1+1))**2+
     @        (Q(I1+2)-Q(J1+2))**2))**3
            D(J,I)=D(I,J)
          END DO
        END DO

        DO I=1,Nbody   !    6
          I1=3*(I-1)+1
          F(I1)=0.0D0
          F(I1+1)=0.0D0
          F(I1+2)=0.0D0
          DO J=1,Nbody   !6
            IF (J.NE.I) THEN
              J1=3*(J-1)+1
              F(I1)=F(I1)+rmass(J)*(Q(J1)-Q(I1))/D(I,J)
              F(I1+1)=F(I1+1)+rmass(J)*(Q(J1+1)-Q(I1+1))/D(I,J)
              F(I1+2)=F(I1+2)+rmass(J)*(Q(J1+2)-Q(I1+2))/D(I,J)
            END IF
          END DO
          F(I1)=K*F(I1)
          F(I1+1)=K*F(I1+1)
          F(I1+2)=K*F(I1+2)
        END DO
      END IF
      RETURN
      END 
c
      SUBROUTINE HAMILTON (N,Q,P,HAMIL,IPAR,rmass,Nbody)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DOUBLE PRECISION  D(10,10),K,rmass(10),Q(N),P(N),Y(100)
      DIMENSION IPAR(*)
      
      IPROB=IPAR(11)
c      IF (IPROB.EQ.1) THEN
c        HAMIL=0.5d0*(P(1)**2+P(2)**2)-1/SQRT(Q(1)**2+Q(2)**2)
c      END IF
c      IF (IPROB.EQ.2) THEN
c        HAMIL=0.5d0*(P(1)**2+Q(1)**2)
c      END IF
c      IF (IPROB.EQ.3) THEN
c        HAMIL=0.5d0*P(1)**2-COS(Q(1))
c      END IF
      iprob=4
      IF (IPROB.EQ.4) THEN
        DO I=1,N
           Y(I)=Q(I)
           Y(I+N)=P(I)
        END DO
c        K=2.95912208286D-4
c        M(1)=0.000954786104043D0
c        M(2)=0.000285583733151D0
c        M(3)=0.0000437273164546D0
c        M(4)=0.0000517759138449D0
c        M(5)=1.0D0/1.3D8
c        M(6)=1.00000597682D0
c
c   change the Gaussian constant to match Don
c
        K=(0.01720209895d0)**2

        DO I=1,Nbody-1   ! 5
          I1=3*(I-1)+1
          DO J=I+1,Nbody   ! 6
            J1=3*(J-1)+1
            D(I,J)=SQRT((Y(I1)-Y(J1))**2+(Y(I1+1)-Y(J1+1))**2+
     @         (Y(I1+2)-Y(J1+2))**2)
            D(J,I)=D(I,J)
          END DO
        END DO

        HAMIL=0.0D0
        DO I=1,Nbody  !  6
           I1=N+3*(I-1)+1
           HAMIL=HAMIL+rmass(I)*(Y(I1)**2+Y(I1+1)**2+Y(I1+2)**2)
        END DO
        HAMIL=HAMIL/2.0D0
        POT=0.0D0
        DO I=2,Nbody   !6
           DO J=1,I-1
              POT=POT+rmass(I)*rmass(J)/D(I,J)
           END DO
        END DO
        HAMIL=HAMIL-K*POT
      END IF
      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE GNI_IRK2(N,FCN,NSTEP,X,P,Q,XEND,
     &    METH,SOLFIX,IOUT,RPAR,IPAR,rmass,nbody,posarray,velarray,
     @    odetime,zzq,timeinterp,Ndyn,h,iGR)
C-----------------------------------------------------------------------
C                 VERSION OF SEPTEMBER 4,2002  
C  E-MAIL CONTACT ADDRESS : Ernst.Hairer@math.unige.ch
C-----------------------------------------------------------------------
C  SOLVES SECOND ORDER ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM
C                       Q'' = F(X,Q)
C  BASED ON THE SYMPLECTIC AND SYMMETRIC GAUSS (IRK) METHODS
C  DESCRIBED IN SECTIONS II.1, VIII.6 OF THE BOOK:
C
C      E. HAIRER, C. LUBICH, G. WANNER, GEOMETRIC NUMERICAL INTEGRATION,
C         STRUCTURE-PRESERVING ALGORITHMS FOR ODES.
C         SPRINGER SERIES IN COMPUT. MATH. 31, SPRINGER 2002.
C
C  AND IN THE PUBLICATION
C
C      E. HAIRER, M. HAIRER, GNI-CODES - MATLAB PROGRAMS FOR
C         GEOMETRIC NUMERICAL INTEGRATION.
C
C  INPUT..
C     N           DIMENSION OF Q AND F(X,Q) 
C
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING F(X,Q):
C                    SUBROUTINE FCN(N,X,Q,F,RPAR,IPAR)
C                    REAL*8 Q(N),F(N)
C                    F(1)=...   ETC.
C
C     NSTEP       NUMBER OF INTEGRATION STEPS
C                    CONSTANT STEP SIZE, H=(XEND-X)/NSTEP
C
C     X           INITIAL X-VALUE
C     P(N)        INITIAL VELOCITY VECTOR
C     Q(N)        INITIAL POSITION VECTOR
C     XEND        FINAL X-VALUE
C
C     METH        NUMBER OF STAGES OF THE GAUSS METHOD
C                    FOR THE MOMENT ONLY POSSIBLE VALUES: 2,4,6.
C
C     SOLFIX      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION. 
C                 IF IOUT=1, IT IS CALLED AFTER EVERY STEP.
C                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. 
C                    SUBROUTINE SOLFIX (NR,XOLD,X,P,Q,N,IRTRN,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(N),CONT(LRC)
C                      ....  
C                 SOLFIX FURNISHES THE SOLUTION "Q,P" AT THE NR-TH
C                    GRID-POINT "X" (INITIAL VALUE FOR NR=0).
C                 "XOLD" IS THE PRECEEDING GRID-POINT.
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, RETURN TO THE CALLING PROGRAM.
C     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLFIX:
C                    IOUT=0: SUBROUTINE IS NEVER CALLED
C                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
C
C     RPAR(LR)    REAL PARAMETER ARRAY; LR MUST BE AT LEAST LR=10
C                    RPAR(1),...,RPAR(10) SERVE AS PARAMETERS FOR
C                    THE CODE. FURTHER VALUES CAN BE USED FOR DEFINING
C                    PARAMETERS IN THE PROBLEM
C     IPAR(LI)    INTEGER PARAMETER ARRAY; LI MUST BE AT LEAST LI=10
C                    IPAR(1),...,IPAR(10) SERVE AS PARAMETERS FOR
C                    THE CODE. FURTHER VALUES CAN BE USED FOR DEFINING
C                    PARAMETERS IN THE PROBLEM
C
C  OUTPUT..
C     P(N)        SOLUTION (VELOCITY) AT XEND
C     Q(N)        SOLUTION (POSITION) AT XEND
C-----------------------------------------------------------------------
C     SOPHISTICATED SETTING OF PARAMETERS 
C-----------------------------------------------------------------------
C    RPAR(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
C    IPAR(1)   NITMAX, MAXIMAL NUMER OF FIXED POINT ITERAT., DEFAULT 50
C-----------------------------------------------------------------------
      PARAMETER (NDGL=500,NSD=6,NMD=3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(NDGL*NSD),EQ(NDGL),EP(NDGL),YH(NDGL),QQ(NDGL)
      DIMENSION C(NSD),AA(NSD,NSD),E(NSD,NSD+NMD),B(NSD),BC(NSD)
      DIMENSION SM(NMD),AM(NSD+NMD),AAT(NSD,NSD)
      DIMENSION Q(N),P(N),FS(NDGL),PS(NDGL),ZQ(NDGL,NSD)
      DIMENSION IPAR(*),RPAR(*),odetime(Ndyn)
      dimension rmass(10),posarray(Ndyn,30),velarray(Ndyn,30)
      dimension zzq(6,60,Ndyn),timeinterp(6)
      dimension FM(30,6),result(30,6),fakeB(6),fakeC(6)
      EXTERNAL FCN
      EXTERNAL SOLFIX
C -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0  
      IF (RPAR(1).EQ.0.0D0) THEN
         UROUND=1.0D-16
      ELSE
         UROUND=RPAR(1)
      END IF
C -------- NITMAX, MAXIMAL NUMER OF FIXED POINT ITERATIONS
      IF (IPAR(1).EQ.0) THEN
         NITMAX=50
      ELSE
         NITMAX=IPAR(1)
      END IF
C --------
      NS=METH
      uround=1.d-19
c
c      H=(XEND-X)/NSTEP
c
      CALL COEFG(NS,C,B,BC,NSD,AA,E,NM,SM,NMD,AM,H)
c
      do jj=1,6
        timeinterp(jj)=c(jj)
      enddo
      IF (IOUT.NE.0) CALL SOLFIX (0,X,X,P,Q,N,RPAR,IPAR,rmass,
     @       nbody)
      CALL FCN(N,X,Q,FS,RPAR,IPAR,rmass,nbody)
      DO IS=1,NS
        FAC=C(IS)**2/2
        DO I=1,N
          ZQ(I,IS)=C(IS)*P(I)+FAC*FS(I)
        END DO
      END DO
      DO I=1,N
        PS(I)=P(I)
        EQ(I)=0.0D0
        EP(I)=0.0D0
      END DO
C --- LOOP FOR THE ITERATIONS

      DO ISTEP=1,NSTEP
        IF (ISTEP.GT.1) CALL STARTB (FCN,N,X,P,Q,NS,NDGL,FS,PS,
     &              ZQ,NSD,E,YH,NM,SM,AM,F,C,RPAR,IPAR,rmass,nbody,h)
C --- FIXED POINT ITERATION

          do i=1,N
            do is=1,6
c              zzq(istep,i,is)=Q(I)+zq(i,is)   
              zzq(is,i,istep)=Q(I)+zq(i,is)   
            enddo
          enddo        
        odetime(istep)=x
        do ii=1,N
          posarray(istep,ii)=Q(ii)
          velarray(istep,ii)=P(ii)
        enddo


        NITER=0
        DYNOLD=0.0D0
  40    CONTINUE
        CALL RKNITE(FCN,N,NS,X,Q,P,NSD,AA,C,NDGL,QQ,ZQ,F,DYNO,RPAR,IPAR,
     &        rmass,nbody)
        NITER=NITER+1
        if(iGR.gt.0.and.niter.ge.2)go to 50
        IF (DYNOLD.LT.DYNO.AND.DYNO.LT.10*UROUND) GOTO 50
        IF (NITER.GE.NITMAX) THEN
          WRITE (6,*) 'NO CONVERGENCE OF ITERATION',DYNO
          IF (DYNO.GT.UROUND*1.D6) RETURN
        END IF
        IF (DYNO.GT.0.1D0*UROUND) GOTO 40
 50     CONTINUE
c
        call gaussQ(6,fakeC,fakeB,6,AAT)
        do JS=1,6
          do j=1,N
            FM(j,JS)=F(j+(JS-1)*N) 
          enddo
        enddo
 
        do i=1,N
          do j=1,6
            dotprodij=0.0d0
            do k=1,6
              dotprodij=dotprodij+FM(i,k)*AAT(j,k)
            enddo
            result(i,j)=dotprodij
          enddo
        enddo
c
        do ii=1,N
          do is=1,6
            zzq(is,ii+N,istep)=P(ii)  +result(ii,is)*h
          enddo
        enddo

c
        do i=1,N
          do is=1,6
c            zzq(istep,i,is)=Q(I)+zq(i,is)   
            zzq(is,i,istep)=Q(I)+zq(i,is)   
          enddo
        enddo        

c
C --- UPDATE OF THE SOLUTION

        X=X+H
        DO I=1,N
          EQI=EQ(I)
          QI=Q(I)
          AY=QI
          EQI=EQI+H*P(I)
          QI=AY+EQI
          EQI=EQI+(AY-QI)
          DO IS=1,NS
c            zzq(istep,i,is)=zq(i,is)+QQ(IS)
            AY=QI
            EQI=EQI+F(I+(IS-1)*N)*BC(IS)
            QI=AY+EQI
            EQI=EQI+(AY-QI)
c            zzq(istep,i,is)=QQ(I)+zq(i,is)   
          END DO
          AY=Q(I)
          Q(I)=QI
          EQ(I)=EQI
          EPI=EP(I)
          PI=P(I)
          DO IS=1,NS/2
            AY=PI
            EPI=EPI+(F(I+(IS-1)*N)+F(I+(NS-IS)*N))*B(IS)
            PI=AY+EPI
            EPI=EPI+(AY-PI)
          END DO
          P(I)=PI
          EP(I)=EPI
        END DO

        IF (IOUT.NE.0) CALL SOLFIX (ISTEP,X-H,X,P,Q,N,RPAR,IPAR,
     $        rmass,nbody)
c
c  the times were here
c
      END DO

      RETURN
      END
C
      SUBROUTINE STARTB (FCN,N,X,P,Q,NS,NDGL,FS,PS,ZQ,NSD,E,YH,
     &                   NM,SM,AM,F,C,RPAR,IPAR,rmass,nbody,h)
C ----------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ZQ(NDGL,NS),E(NSD,NS+NM),PS(N),F(N*NS),C(NS)
      REAL*8 AM(NS+NM),SM(NM),P(N),Q(N),YH(NDGL),FS(NDGL)
      DIMENSION IPAR(*),RPAR(*)
      EXTERNAL FCN
      dimension rmass(10)
        NS1=NS+1
        NS2=NS+2
        NSM=NS+NM
        DO I=1,N
          SAV=0.0D0
          DO JS=1,NS
            SAV=SAV+AM(JS)*ZQ(I,JS)
          END DO
          YH(I)=SAV+AM(NS1)*PS(I)+AM(NS2)*P(I)+Q(I)
          DO IS=1,NS
            SAV=0.0D0
            DO JS=1,NS
              SAV=SAV+E(IS,JS)*F(I+(JS-1)*N)
            END DO
            ZQ(I,IS)=SAV+E(IS,NS1)*FS(I)
          END DO
        END DO
C
        CALL FCN(N,X+H,Q,FS,RPAR,IPAR,rmass,nbody)
        CALL FCN(N,X+H*SM(NM),YH,F,RPAR,IPAR,rmass,nbody)
        DO I=1,N
          PS(I)=P(I)
          DO IS=1,NS
            ZQIIS=ZQ(I,IS)+E(IS,NS2)*FS(I)+E(IS,NSM)*F(I)
            ZQ(I,IS)=ZQIIS+C(IS)*P(I)
          END DO
        END DO
      RETURN
      END
C
      SUBROUTINE RKNITE(FCN,N,NS,X,Q,P,NSD,AA,C,NDGL,QQ,ZQ,
     &                  F,DYNO,RPAR,IPAR,rmass,nbody)
C ----------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ZQ(NDGL,NS),C(NS),Q(N),P(N)
      REAL*8 AA(NSD,NS),F(NDGL*NS),QQ(NDGL)
      DIMENSION IPAR(*),RPAR(*)
      dimension rmass(10)
      EXTERNAL FCN
C ---
      DO JS=1,NS
        DO J=1,N
          QQ(J)=Q(J)+ZQ(J,JS)
        END DO
        CALL FCN(N,X+C(JS),QQ,F(1+(JS-1)*N),RPAR,IPAR,rmass,nbody)
      END DO
C ---
      DYNO=0.D0
      DO I=1,N
        DNOM=MAX(1.0D-1,ABS(Q(I)))
        DO IS=1,NS
          SSUM=C(IS)*P(I)
          DO JS=1,NS
            SSUM=SSUM+AA(IS,JS)*F(I+(JS-1)*N)
          END DO
          DYNO=DYNO+((SSUM-ZQ(I,IS))/DNOM)**2
          ZQ(I,IS)=SSUM
        END DO
      END DO
      DYNO=DSQRT(DYNO/(NS*N))
C ---
      RETURN
      END
C
      SUBROUTINE COEFG(NS,C,B,BC,NSD,AA,E,NM,SM,NMD,AM,HSTEP)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 C(NS),AA(NSD,NS),E(NSD,NS+NMD),B(NS)
      REAL*8 AM(NSD+NMD),SM(NMD),BC(NS)
 
      NM=3
      IF (NS.EQ.2) THEN
         C(1)= 0.21132486540518711775D+00
         C(2)= 0.78867513459481288225D+00
         B(1)= 0.50000000000000000000D+00
         B(2)= 0.50000000000000000000D+00
         BC(1)= 0.39433756729740644113D+00
         BC(2)= 0.10566243270259355887D+00
         AA(1,1)= 0.41666666666666666667D-01
         AA(1,2)=-0.19337567297406441127D-01
         AA(2,1)= 0.26933756729740644113D+00
         AA(2,2)= 0.41666666666666666667D-01
         E(1,1)=-0.28457905077110526160D-02
         E(1,2)=-0.63850024471784160410D-01
         E(1,3)= 0.48526095198694517563D-02
         E(1,4)= 0.11305688530429939012D+00
         E(1,5)=-0.28884580475413403312D-01
         E(2,1)= 0.41122751744511433137D-01
         E(2,2)=-0.18654814888622834132D+00
         E(2,3)=-0.18110185277445209332D-01
         E(2,4)= 0.36674109449368040786D+00
         E(2,5)= 0.10779872188955481745D+00
         SM(1)= 0.00000000000000000000D+00
         SM(2)= 0.10000000000000000000D+01
         SM(3)= 0.16000000000000000000D+01
         AM(1)= 0.25279583039343438291D+02
         AM(2)=-0.86907830393434382912D+01
         AM(3)=-0.80640000000000000000D+00
         AM(4)= 0.29184000000000000000D+01
         AM(5)= 0.00000000000000000000D+00
      END IF
      IF (NS.EQ.4) THEN
         C(1)= 0.69431844202973712388D-01
         C(2)= 0.33000947820757186760D+00
         C(3)= 0.66999052179242813240D+00
         C(4)= 0.93056815579702628761D+00
         B(1)= 0.17392742256872692869D+00
         B(2)= 0.32607257743127307131D+00
         B(3)= 0.32607257743127307131D+00
         B(4)= 0.17392742256872692869D+00
         BC(1)= 0.16185132086231030665D+00
         BC(2)= 0.21846553629538057030D+00
         BC(3)= 0.10760704113589250101D+00
         BC(4)= 0.12076101706416622036D-01
         AA(1,1)= 0.40381914508467311298D-02
         AA(1,2)=-0.32958609449446961650D-02
         AA(1,3)= 0.26447829520668538006D-02
         AA(1,4)=-0.97672296325588161023D-03
         AA(2,1)= 0.43563580902396261254D-01
         AA(2,2)= 0.13818951406296126013D-01
         AA(2,3)=-0.43401341944349953440D-02
         AA(2,4)= 0.14107297391595337720D-02
         AA(3,1)= 0.10586435263357640763D+00
         AA(3,2)= 0.10651836096505307395D+00
         AA(3,3)= 0.13818951406296126013D-01
         AA(3,4)=-0.17580153590805494993D-02
         AA(4,1)= 0.14879849619263780300D+00
         AA(4,2)= 0.19847049885237718995D+00
         AA(4,3)= 0.81671359795877570687D-01
         AA(4,4)= 0.40381914508467311298D-02
         E(1,1)=-0.21272768296134340207D-01
         E(1,2)= 0.11059138674756969912D-01
         E(1,3)= 0.38999255049973564023D-02
         E(1,4)=-0.43986226789008967612D-01
         E(1,5)= 0.13581590305438849621D-01
         E(1,6)= 0.39922421675314269059D-01
         E(1,7)=-0.79369058065113002021D-03
         E(2,1)=-0.75671119283734809953D-02
         E(2,2)= 0.10209394000843457002D-01
         E(2,3)=-0.12880197817980892596D-01
         E(2,4)=-0.56381316813776501277D-01
         E(2,5)= 0.37440782682669799960D-02
         E(2,6)= 0.11522469441011273193D+00
         E(2,7)= 0.21035877343246316334D-02
         E(3,1)=-0.39890571772473709759D+00
         E(3,2)= 0.26819725655216894347D+00
         E(3,3)=-0.82551711648854471247D-01
         E(3,4)=-0.85516559106259630212D+00
         E(3,5)= 0.24433810515772642570D+00
         E(3,6)= 0.10234155624049009806D+01
         E(3,7)= 0.25115745967236579242D-01
         E(4,1)=-0.40964796048052939224D+00
         E(4,2)= 0.29949323098224574487D+00
         E(4,3)=-0.13867460566101912494D+00
         E(4,4)=-0.98859300714628940382D+00
         E(4,5)= 0.24671351779481625627D+00
         E(4,6)= 0.12912760231350872304D+01
         E(4,7)= 0.13241134766742798418D+00
         SM(1)= 0.00000000000000000000D+00
         SM(2)= 0.10000000000000000000D+01
         SM(3)= 0.16500000000000000000D+01
         AM(1)= 0.10806374869244001787D+04
         AM(2)=-0.66008818661284690206D+03
         AM(3)= 0.61810154357557529566D+03
         AM(4)=-0.31341427826212857229D+03
         AM(5)=-0.10187174765625000000D+02
         AM(6)= 0.31173050390625000000D+02
         AM(7)= 0.00000000000000000000D+00
      END IF
      IF (NS.EQ.6) THEN
         C(1)= 0.33765242898423986094D-01
         C(2)= 0.16939530676686774317D+00
         C(3)= 0.38069040695840154568D+00
         C(4)= 0.61930959304159845432D+00
         C(5)= 0.83060469323313225683D+00
         C(6)= 0.96623475710157601391D+00
         B(1)= 0.85662246189585172520D-01
         B(2)= 0.18038078652406930378D+00
         B(3)= 0.23395696728634552369D+00
         B(4)= 0.23395696728634552369D+00
         B(5)= 0.18038078652406930378D+00
         B(6)= 0.85662246189585172520D-01
         BC(1)= 0.82769839639769234611D-01
         BC(2)= 0.14982512785597570103D+00
         BC(3)= 0.14489179419935320895D+00
         BC(4)= 0.89065173086992314743D-01
         BC(5)= 0.30555658668093602753D-01
         BC(6)= 0.28924065498159379092D-02
         AA(1,1)= 0.90625420195651151857D-03
         AA(1,2)=-0.72859711612531400024D-03
         AA(1,3)= 0.79102695861167691135D-03
         AA(1,4)=-0.70675390218535384182D-03
         AA(1,5)= 0.45647714224056921122D-03
         AA(1,6)=-0.14836147050330408643D-03
         AA(2,1)= 0.11272367531794365387D-01
         AA(2,2)= 0.39083482447840698486D-02
         AA(2,3)=-0.14724868010943911900D-02
         AA(2,4)= 0.10992669056588431310D-02
         AA(2,5)=-0.67689040729401428165D-03
         AA(2,6)= 0.21677950347174141516D-03
         AA(3,1)= 0.30008019623627547434D-01
         AA(3,2)= 0.36978289259468146662D-01
         AA(3,3)= 0.65490339168957822692D-02
         AA(3,4)=-0.16615098173008262274D-02
         AA(3,5)= 0.84753461862041607649D-03
         AA(3,6)=-0.25877462623437421721D-03
         AA(4,1)= 0.49900269650650898941D-01
         AA(4,2)= 0.82003427445271620462D-01
         AA(4,3)= 0.54165111295060067982D-01
         AA(4,4)= 0.65490339168957822692D-02
         AA(4,5)=-0.11352871017627472322D-02
         AA(4,6)= 0.28963081055952389031D-03
         AA(5,1)= 0.68475836671617248304D-01
         AA(5,2)= 0.11859257878058808400D+00
         AA(5,3)= 0.10635984886129551097D+00
         AA(5,4)= 0.47961474042181382443D-01
         AA(5,5)= 0.39083482447840698486D-02
         AA(5,6)=-0.34600839001342442657D-03
         AA(6,1)= 0.79729071619449992615D-01
         AA(6,2)= 0.14419100392702230613D+00
         AA(6,3)= 0.13628542646896576408D+00
         AA(6,4)= 0.81956586217401900627D-01
         AA(6,5)= 0.23736460480774324642D-01
         AA(6,6)= 0.90625420195651151857D-03
         E(1,1)=-0.16761132335280609813D-01
         E(1,2)= 0.10201050166615899799D-01
         E(1,3)=-0.58593121685075943100D-02
         E(1,4)=-0.11907383391366998251D-03
         E(1,5)= 0.10615611118132982241D-01
         E(1,6)=-0.30692054230989138447D-01
         E(1,7)= 0.10615182045216224925D-01
         E(1,8)= 0.22586707045496892369D-01
         E(1,9)=-0.16931992776201068110D-04
         E(2,1)= 0.10671755276327262128D-01
         E(2,2)=-0.51098203653251450913D-02
         E(2,3)= 0.16062647299186369205D-03
         E(2,4)= 0.64818802653621866868D-02
         E(2,5)=-0.12132386914873895089D-01
         E(2,6)=-0.99709737725909584834D-02
         E(2,7)=-0.70287093442894942752D-02
         E(2,8)= 0.31243249755879001843D-01
         E(2,9)= 0.31763603839792897936D-04
         E(3,1)=-0.40875203230945019464D+00
         E(3,2)= 0.28214948905763253599D+00
         E(3,3)=-0.22612660499718519054D+00
         E(3,4)= 0.13640993962985420478D+00
         E(3,5)= 0.15888529591697266925D+00
         E(3,6)=-0.11667863471317749710D+01
         E(3,7)= 0.25224964119340060668D+00
         E(3,8)= 0.10440940643938620983D+01
         E(3,9)= 0.33914722176493324285D-03
         E(4,1)=-0.29437531285359759661D+01
         E(4,2)= 0.20017220470127690267D+01
         E(4,3)=-0.15383035791443948798D+01
         E(4,4)= 0.78114323215109899716D+00
         E(4,5)= 0.13930345104184182146D+01
         E(4,6)=-0.75958281612589849630D+01
         E(4,7)= 0.18220129530415584951D+01
         E(4,8)= 0.62663163493155487560D+01
         E(4,9)= 0.54279630166374655267D-02
         E(5,1)=-0.79572842006457093076D+01
         E(5,2)= 0.53527892762707449170D+01
         E(5,3)=-0.40049139768467199697D+01
         E(5,4)= 0.18326058141135591515D+01
         E(5,5)= 0.39753886181058367500D+01
         E(5,6)=-0.19423696478604790213D+02
         E(5,7)= 0.49362128400107292627D+01
         E(5,8)= 0.15601708062381928560D+02
         E(5,9)= 0.32142123424873719685D-01
         E(6,1)=-0.78463118056075171475D+01
         E(6,2)= 0.53580869574441241664D+01
         E(6,3)=-0.41476905275607763365D+01
         E(6,4)= 0.21275912797813913113D+01
         E(6,5)= 0.37642416878253538582D+01
         E(6,6)=-0.20329681631523484613D+02
         E(6,7)= 0.48515418060343387549D+01
         E(6,8)= 0.16604467346259915039D+02
         E(6,9)= 0.84559690262225766975D-01
         SM(1)= 0.00000000000000000000D+00
         SM(2)= 0.10000000000000000000D+01
         SM(3)= 0.17500000000000000000D+01
         AM(1)= 0.58080578375796358720D+05
         AM(2)=-0.33214989339522861968D+05
         AM(3)= 0.28376088288312020853D+05
         AM(4)=-0.27923430684614999462D+05
         AM(5)= 0.29743005589491042677D+05
         AM(6)=-0.15525927919158826444D+05
         AM(7)=-0.27700591278076171875D+03
         AM(8)= 0.73086943817138671875D+03
         AM(9)= 0.00000000000000000000D+00
      END IF
C
            HSTEP2=HSTEP*HSTEP
            DO IS=1,NS
              B(IS)=HSTEP*B(IS)
              BC(IS)=HSTEP2*BC(IS)
              C(IS)=HSTEP*C(IS)
              DO JS=1,NS
                AA(IS,JS)=HSTEP2*AA(IS,JS)
                E(IS,JS)=HSTEP2*E(IS,JS)
              END DO
            END DO
            DO IM=1,NM
              DO IS=1,NS
                E(IS,NS+IM)=HSTEP2*E(IS,NS+IM)
              END DO
              AM(NS+IM)=HSTEP*AM(NS+IM)
            END DO
C
      RETURN
      END
c
c  #$%^#^#%#*#&#^#%#$#%#^#&#&#^#%#$#
c
          subroutine LTTsky(odetime,zzq,Nstep,timeinterp,ibody,timein,
     @        xout,yout,velout,ii,Ndyn)
c
          implicit double precision(a-h,o-z)
          dimension odetime(Ndyn),zzq(6,60,Ndyn),timeinterp(6)
          dimension xx(10),ff(10)
c
          rsunperAU=0.00465116d0     !solar radius in AU
          AC=1.49597870691D11/86400.0d0  ! AU/day --> m/sec
          speed=299792458.0d0/AC    !speed of light in AU/day

          jndex=3*(ibody-1)
          call divdiffzonly(odetime,zzq,Nstep,timeinterp,jndex+3,timein,
     @        zout1,velz1,ii,Ndyn)
          zsave=z1out

          Tr=-zout1/speed
          xx(1)=-Tr
          call divdiffzonly(odetime,zzq,Nstep,timeinterp,jndex+3,timein+xx(1),
     @        zout1,velz1,ii,Ndyn)
          ff(1)=zout1-speed*(xx(1))

          xx(2)=Tr          
          call divdiffzonly(odetime,zzq,Nstep,timeinterp,jndex+3,timein+xx(2),
     @        zout1,velz1,ii,Ndyn)
          ff(2)=zout1-speed*(xx(2))
c
          deltaT=0.0d0
          do kk=3,10
            tol=abs(ff(kk-1)-ff(kk-2))
            if(tol.lt.1.0d-13)go to 68
            xx(kk)=(xx(kk-2)*ff(kk-1)-xx(kk-1)*ff(kk-2))/(ff(kk-1)-
     @               ff(kk-2))
            call divdiff(odetime,zzq,Nstep,timeinterp,jndex+3,timein+
     @          xx(kk),zout1,velz1,ii,Ndyn)
            ff(kk)=zout1-speed*(xx(kk))
          enddo
c 68       write(*,*)xx(1)*86400.0d0,ff(1),ff(2)
c          write(*,*)xx(kk-1)*86400.0d0,ff(kk-1),kk
c
 68       deltaT=xx(kk-1)

c          deltaT=-Tr     !this gives nearly the same as the iteration above
 
c          deltaT=0.0d0

 69       call divdiffzonly(odetime,zzq,Nstep,timeinterp,jndex+1,timein+deltaT,
     @        xout1,velx1,ii,Ndyn)
          call divdiffzonly(odetime,zzq,Nstep,timeinterp,jndex+2,timein+deltaT,
     @        yout1,vely1,ii,Ndyn)

c          call divdiff(odetime,zzq,Nstep,timeinterp,jndex+1,timein,
c     @        xold,velx1,ii,Ndyn)
c          call divdiff(odetime,zzq,Nstep,timeinterp,jndex+2,timein,
c     @        yold,vely1,ii,Ndyn)
c          write(*,6969)ibody,timein,86400.0d0*deltaT,velx1,vely1
c          write(*,6970)(xout1-(xold+velx1*deltaT))/rsunperAU,(yout1-yold)/rsunperAU
 6969     format(i1,2x,4(f15.8,1x))
 6970     format(10x,2(f15.10,2x))
c
          velout=velz1
c          zout1=dsqrt(xout1**2+yout1**2+zout1**2)
          xout=xout1 ! +velx1*zsave/speed  ! adding this hardly matters
          yout=yout1 ! +vely1*zsave/speed
          xout=xout/rsunperAU
          yout=yout/rsunperAU
          velout=velout*AC/1000.0d0
c
          return
          end
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine newtondiff(odetime,zzq,Nstep,timeinterp,Nbody,timein,
     @        xout,velout,ii,Ndyn)
c
          implicit double precision(a-h,o-z)
          dimension odetime(Ndyn),zzq(6,60,Ndyn),timeinterp(6)
          dimension xx(6),yy(6),d(6,6),A(6)
c
c
c          call hunt(odetime,Nstep,timein,ii)
c
          do ii=1,Nstep
            do kkk=1,Nbody*3
              do jj=1,6
                xx(jj)=odetime(ii)+timeinterp(jj)
                d(jj,1)=zzq(jj,kkk,ii)
              enddo
c
c   Construct the divided difference table          
c
              do j=2,6 
                do i=j,6
                  d(i,j)=(d(i,j-1)-d(i-1,j-1))/(xx(i)-xx(i-j+1))
                enddo
              enddo
c
              do i=1,6
                zzq(i,kkk,ii)=d(i,i)
              enddo
c
            enddo ! 1,Nbody*3
          enddo   !1, Ndyn
          return
          end
c
c     #!!!%(@^%#%(*&(#@#*$@^$##&#@
c
          subroutine divdiff(odetime,zzq,Nstep,timeinterp,ibody,timein,
     @        xout,velout,ii,Ndyn)
c
          implicit double precision(a-h,o-z)
          dimension odetime(Ndyn),zzq(6,60,Ndyn),timeinterp(6)
          dimension xx(6),yy(6),d(6,6),A(6)
c
c
          call hunt(odetime,Nstep,timein,ii)

          if(ii.le.0)ii=1

          do jj=1,6
            xx(jj)=odetime(ii)+timeinterp(jj)
          enddo

cc            d(jj,1)=zzq(ii,ibody,jj)
c            d(jj,1)=zzq(jj,ibody,ii)
c          enddo
cc
cc   Construct the divided difference table          
cc
c          do j=2,6 
c            do i=j,6
c              d(i,j)=(d(i,j-1)-d(i-1,j-1))/(xx(i)-xx(i-j+1))
c            enddo
c          enddo
cc
          do i=1,6
            A(i)=zzq(i,ibody,ii)
          enddo
c
          xout=A(6)          
          velout=xout
c
          do i=5,2,-1
            xout=A(i)+(timein-xx(i))*xout
            velout=xout+(timein-xx(i-1))*velout
          enddo
c
          xout=A(1)+(timein-xx(1))*xout
c
          return
          end
c
c     #!!!%(@^%#%(*&(#@#*$@^$##&#@
c
          subroutine divdiffzonly(odetime,zzq,Nstep,timeinterp,ibody,
     @      timein,xout,velout,ii,Ndyn)
c
          implicit double precision(a-h,o-z)
          dimension odetime(Ndyn),zzq(6,60,Ndyn),timeinterp(6)
          dimension xx(6),yy(6),d(6,6),A(6)
c
c
          call hunt(odetime,Nstep,timein,ii)

          if(ii.le.0)ii=1

          do jj=1,6
            xx(jj)=odetime(ii)+timeinterp(jj)
          enddo
c
cc            d(jj,1)=zzq(ii,ibody,jj)
c            d(jj,1)=zzq(jj,ibody,ii)
c          enddo
cc
cc   Construct the divided difference table          
cc
c          do j=2,6 
c            do i=j,6
c              d(i,j)=(d(i,j-1)-d(i-1,j-1))/(xx(i)-xx(i-j+1))
c            enddo
c          enddo
cc
          do i=1,6
            A(i)=zzq(i,ibody,ii)
          enddo
c
          xout=A(6)          
c
c    skip the velocity
c
c          velout=xout
c
          do i=5,2,-1
            xout=A(i)+(timein-xx(i))*xout
c            velout=xout+(timein-xx(i-1))*velout
          enddo
c
          xout=A(1)+(timein-xx(1))*xout
c
          return
          end
c
c     #!!!%(@^%#%(*&(#@#*$@^$##&#@
c
          subroutine divdiffonestep(x,zzq,Nstep,timeinterp,ibody,timein,
     @        xout,velout,ii,Ndyn)
c
          implicit double precision(a-h,o-z)
          dimension zzq(6,60,Ndyn),timeinterp(6)
          dimension xx(6),yy(6),d(6,6),A(6)
c
c
c          call hunt(odetime,Nstep,timein,ii)
c
          if(ii.le.0)ii=1
c
          ii=1
          do jj=1,6
            xx(jj)=timeinterp(jj) +x
            d(jj,1)=zzq(jj,ibody,ii)
          enddo
cc
cc   Construct the divided difference table          
cc
          do j=2,6 
            do i=j,6
              d(i,j)=(d(i,j-1)-d(i-1,j-1))/(xx(i)-xx(i-j+1))
            enddo
          enddo
c
          do i=1,6
            A(i)=d(i,i)
          enddo
c
          xout=A(6)          
          velout=xout
c
          do i=5,2,-1
            xout=A(i)+(timein-xx(i))*xout
            velout=xout+(timein-xx(i-1))*velout
          enddo
c
          xout=A(1)+(timein-xx(1))*xout
c
          return
          end
c
c     #!!!%(@^%#%(*&(#@#*$@^$##&#@
c
          subroutine M2toM1centric(rM1,rM2,period,ecc,rMeanAnom,
     @       argper,finc,Omega,QQm2,PPm2)
c
c   this routine will compute the xyz M1-centric coordinates and velocity
c   of M2 given these orbital parameters.
c
c   rM1 and rM2 are in solar masses, period is in days.  The angles
c   argper, finc, and Omega are in degrees
c
c   The mean anomaly is a phase, and will be in radians
c
          implicit double precision(A-h,o-z)
c
          dimension QQm2(3),PPm2(3),PeriX(3)
c
c   Here is the Gravitational Constant in AU**2/solar_mass/day**2 units
c
          pie=4.0d0*datan(1.0d0)
          G=(0.01720209895d0)**2
c
c   Find the separation
c          
          a=(G*(rM1+rM2)*Period*Period/(4.0d0*pie*pie))**(1.0d0/3.0d0)
          p=a*(1.0d0-ecc**2)
          h=dsqrt(G*(rM1+rM2)*p)
c
c   Find the true anomaly f
c
c          write(*,444)(rM1+rM2)*G,rMeanAnom*180.0d0/pie
c          write(*,555)(rM1+rM2)*G,a,period,ecc
c          write(*,556)finc*180.0d0/pie,argper*180.0d0/pie,Omega*180.0d0/pie,rMeanAnom*180.0d0/pie
c          write(*,*)' '
c 555      format(4(1pe16.9,3x))
c 556      format(4(1pe16.9,3x))
c
c 444      format(2(f17.12,3x),'   total mass, rmean')
c          write(*,445)ecc,finc*180.0d0/pie,argper*180.0d0/pie,Omega*180.0d0/pie
c 445      format(f10.8,2x,3(f12.7,2x),'  ecc, finc, argper, Omega')

          cc=rMeanAnom
          cc=dmod(cc,2.0d0*pie)
          cc=cc-2.0d0*pie
          bigE=cc
          if(bigE.lt.2.0d0*pie)bigE=bigE+2.0d0*pie
c
          do i=1,10
            bigE=bigE-(bigE-ecc*dsin(bigE)-cc)/(1.0d0-ecc*dcos(bigE))
          enddo
c
          f=2.0d0*datan(dsqrt((1.0d0+ecc)/(1.0d0-ecc))*dtan(bigE/2.0d0))
          r=a*(1.0d0-ecc*dcos(bigE))
c
c     position in orbital plane
c 
          QQm2(1)=r*dcos(f)
          QQm2(2)=r*dsin(f)
          QQm2(3)=0.0d0
c
          PeriX(1)=a*(1.0d0-ecc)
          PeriX(2)=0.0d0
          PeriX(3)=0.0d0
c
c     velocity in orbital plane
c
          PPm2(1)=h*(-dsin(f)/r+ecc/p*dsin(f)*dcos(f))
          PPm2(2)=h*(dcos(f)/r+ecc/p*dsin(f)*dsin(f))
          PPm2(3)=0.0d0
c
c         write(*,123)QQm2(1),QQm2(2),QQ
c         write(*,123)PPm2(1),PPm2(2)
c123      format(2(f16.9,2x),'    old')

c
c   rotate through the longitude angle, then the inclination angle, and then
c   through argper angle
c          
          fincr=finc  !*pie/180.0d0
          argperrad=argper  !*pie/180.0d0
          Omegarad=Omega   !*pie/180.0d0
c
c          wfcos=dcos(argperrad+f)
c          wfsin=dsin(argperrad+f)
          Omcos=dcos(Omegarad)
          Omsin=dsin(Omegarad)
          finccos=dcos(fincr)
          fincsin=dsin(fincr)
c          
          argcos=dcos(argperrad)
          argsin=dsin(argperrad)
c
c   rotate velocity vector by -argper
c
          t1=QQm2(1)*argcos-QQm2(2)*argsin
          t2=QQm2(1)*argsin+QQm2(2)*argcos
          t3=QQm2(3)
          QQm2(1)=t1
          QQm2(2)=t2
          QQm2(3)=t3
c
c   rotate the result by finc
c
          t1=QQm2(1)
          t2=QQm2(2)*finccos-QQm2(3)*fincsin
          t3=QQm2(2)*fincsin+QQm2(3)*finccos
          QQm2(1)=t1
          QQm2(2)=t2
          QQm2(3)=t3
c
c   rotate by -Omega
c 
          t1=QQm2(1)*Omcos-QQm2(2)*Omsin
          t2=QQm2(1)*Omsin+QQm2(2)*Omcos
          t3=QQm2(3)
          QQm2(1)=t1
          QQm2(2)=t2
          QQm2(3)=t3
c
c   rotate velocity vector by -argper
c
          t1=PPm2(1)*argcos-PPm2(2)*argsin
          t2=PPm2(1)*argsin+PPm2(2)*argcos
          t3=PPm2(3)
          PPm2(1)=t1
          PPm2(2)=t2
          PPm2(3)=t3
c
c   rotate the result by finc
c
          t1=PPm2(1)
          t2=PPm2(2)*finccos-PPm2(3)*fincsin
          t3=PPm2(2)*fincsin+PPm2(3)*finccos
          PPm2(1)=t1
          PPm2(2)=t2
          PPm2(3)=t3
c
c   rotate by -Omega
c 
          t1=PPm2(1)*Omcos-PPm2(2)*Omsin
          t2=PPm2(1)*Omsin+PPm2(2)*Omcos
          t3=PPm2(3)
          PPm2(1)=t1
          PPm2(2)=t2
          PPm2(3)=t3
c

c          write(*,333)QQm2(1),QQm2(2),QQm2(3)
c          write(*,333)PPm2(1),PPm2(2),PPm2(3)

c 333      format(3(1pe16.9,2x))

          return
          end
c
c   #!!!%(@^%#%*(&(#@#*$^@^$##*#@&(%)#!!!%(@^%#%(*&(#@#*$@^$##&#@
c
          subroutine newM2toM1centric(rM1,rM2,period,ecc,rMeanAnom,
     @       argper,finc,Omega,QQm2,PPm2)
c
c   this routine will compute the xyz M1-centric coordinates and velocity
c   of M2 given these orbital parameters.
c
c   rM1 and rM2 are in solar masses, period is in days.  The angles
c   argper, finc, and Omega are in degrees
c
c   The mean anomaly is a phase, and will be in radians
c
          implicit double precision(A-h,o-z)
c
          dimension QQm2(3),PPm2(3),PeriX(3)
c
c   Here is the Gravitational Constant in AU**2/solar_mass/day**2 units
c
          pie=4.0d0*datan(1.0d0)
          G=(0.01720209895d0)**2
c
c   Find the separation
c          
          a=(G*(rM1+rM2)*Period*Period/(4.0d0*pie*pie))**(1.0d0/3.0d0)
          p=a*(1.0d0-ecc**2)
          h=dsqrt(G*(rM1+rM2)*p)
c
c   Find the true anomaly f
c
          cc=rMeanAnom
          cc=dmod(cc,2.0d0*pie)
          cc=cc-2.0d0*pie
          bigE=cc
          if(bigE.lt.2.0d0*pie)bigE=bigE+2.0d0*pie
c
          do i=1,10
            bigE=bigE-(bigE-ecc*dsin(bigE)-cc)/(1.0d0-ecc*dcos(bigE))
          enddo
c
          f=2.0d0*datan(dsqrt((1.0d0+ecc)/(1.0d0-ecc))*dtan(bigE/2.0d0))
          r=a*(1.0d0-ecc*dcos(bigE))
c
c     position in orbital plane
c 
          QQm2(1)=a*(dcos(bigE)-ecc)    !r*dcos(f)
          QQm2(2)=a*dsqrt(1.0d0-ecc*ecc)*dsin(bigE)  !r*dsin(f)
          QQm2(3)=0.0d0
c
          PeriX(1)=a*(1.0d0-ecc)
          PeriX(2)=0.0d0
          PeriX(3)=0.0d0
c
c     velocity in orbital plane
c
          PPm2(1)=h*(-dsin(f)/r+ecc/p*dsin(f)*dcos(f))
          PPm2(2)=h*(dcos(f)/r+ecc/p*dsin(f)*dsin(f))
          PPm2(3)=0.0d0
c
          rmeanmotion=h/(a*a*dsqrt(1.0d0-ecc*ecc))  !dsqrt(G/(a*a*a))
          PPm2(1)=-a*rmeanmotion*dsin(bigE)/(1.0d0-ecc*dcos(bigE))
          PPm2(2)=dsqrt(1.0d0-ecc*ecc)*a*rmeanmotion*dcos(bigE)/(1.0d0-
     @          ecc*dcos(bigE))
c
c          write(*,123)QQm2(1),QQm2(2)
c          write(*,123)PPm2(1),PPm2(2)
c 123      format(2(f16.9,2x),'    new')
c
c   rotate through the longitude angle, then the inclination angle, and then
c   through argper angle
c          
          fincr=finc  !*pie/180.0d0
          argperrad=argper  !*pie/180.0d0
          Omegarad=Omega   !*pie/180.0d0
c
c          wfcos=dcos(argperrad+f)
c          wfsin=dsin(argperrad+f)
          Omcos=dcos(Omegarad)
          Omsin=dsin(Omegarad)
          finccos=dcos(fincr)
          fincsin=dsin(fincr)
c          
          argcos=dcos(argperrad)
          argsin=dsin(argperrad)
c
c   rotate velocity vector by -argper
c
          t1=QQm2(1)*argcos-QQm2(2)*argsin
          t2=QQm2(1)*argsin+QQm2(2)*argcos
          t3=QQm2(3)
          QQm2(1)=t1
          QQm2(2)=t2
          QQm2(3)=t3
c
c   rotate the result by finc
c
          t1=QQm2(1)
          t2=QQm2(2)*finccos-QQm2(3)*fincsin
          t3=QQm2(2)*fincsin+QQm2(3)*finccos
          QQm2(1)=t1
          QQm2(2)=t2
          QQm2(3)=t3
c
c   rotate by -Omega
c 
          t1=QQm2(1)*Omcos-QQm2(2)*Omsin
          t2=QQm2(1)*Omsin+QQm2(2)*Omcos
          t3=QQm2(3)
          QQm2(1)=t1
          QQm2(2)=t2
          QQm2(3)=t3
c
c   rotate velocity vector by -argper
c
          t1=PPm2(1)*argcos-PPm2(2)*argsin
          t2=PPm2(1)*argsin+PPm2(2)*argcos
          t3=PPm2(3)
          PPm2(1)=t1
          PPm2(2)=t2
          PPm2(3)=t3
c
c   rotate the result by finc
c
          t1=PPm2(1)
          t2=PPm2(2)*finccos-PPm2(3)*fincsin
          t3=PPm2(2)*fincsin+PPm2(3)*finccos
          PPm2(1)=t1
          PPm2(2)=t2
          PPm2(3)=t3
c
c   rotate by -Omega
c 
          t1=PPm2(1)*Omcos-PPm2(2)*Omsin
          t2=PPm2(1)*Omsin+PPm2(2)*Omcos
          t3=PPm2(3)
          PPm2(1)=t1
          PPm2(2)=t2
          PPm2(3)=t3
c
          return
          end
c
c   #!!!%(@^%#%*(&(#@#*$^@^$##*#@&(%)
c
          subroutine newmovetobarycenter(Nbody,rmass,QQ,PP,QQBC,PPBC)
c
          implicit double precision (a-h,o-z)
c
          dimension rmass(10),QQ(30),PP(30),QQBC(30),PPBC(30)
          dimension eta(0:9),tempmass(0:9),tempx(0:9),tempy(0:9)
          dimension tempz(0:9),tempvx(0:9),tempvy(0:9),tempvz(0:9)
          dimension tempmx(0:9),tempmy(0:9),tempmz(0:9)
c
          totalmass=0.0d0
          do i=1,Nbody
            totalmass=totalmass+rmass(i)
          enddo
c
          p0x=0.0d0
          p0y=0.0d0
          p0z=0.0d0
          do i=0,Nbody-1
            tempmass(i)=rmass(i+1)
            tempx(i)=QQ(3*i+1)
            tempy(i)=QQ(3*i+2)
            tempz(i)=QQ(3*i+3)
            tempvx(i)=PP(3*i+1)
            tempvy(i)=PP(3*i+2)
            tempvz(i)=PP(3*i+3)
          enddo
c
          do i=0,Nbody-1
            summ=0.0d0
            do j=0,i
              summ=summ+tempmass(j)
            enddo
            eta(i)=summ
          enddo
c
          tempmx(0)=totalmass*tempvx(0)
          tempmy(0)=totalmass*tempvy(0)
          tempmz(0)=totalmass*tempvz(0)
          p0x=tempmx(0)
          p0y=tempmy(0)
          p0z=tempmz(0)
c
c
          do i=1,Nbody-1
            tempmx(i)=eta(i-1)/eta(i)*tempmass(i)*tempvx(i)
            tempmy(i)=eta(i-1)/eta(i)*tempmass(i)*tempvy(i)
            tempmz(i)=eta(i-1)/eta(i)*tempmass(i)*tempvz(i)
          enddo

          cmx=0.0d0
          cmy=0.0d0
          cmz=0.0d0
          vx=0.0d0
          vy=0.0d0
          vz=0.0d0
          do i=0,0
            summx=0.0d0
            summy=0.0d0
            summz=0.0d0
            vx=0.0d0
            vy=0.0d0
            vz=0.0d0
            do j=1,Nbody-1
              summx=summx+tempmass(j)/eta(j)*tempx(j)
              summy=summy+tempmass(j)/eta(j)*tempy(j)
              summz=summz+tempmass(j)/eta(j)*tempz(j)
              vx=vx+tempmass(i)/eta(j-1)*tempmx(j)
              vy=vy+tempmass(i)/eta(j-1)*tempmy(j)
              vz=vz+tempmass(i)/eta(j-1)*tempmz(j)
            enddo
            QQBC(3*i+1)=-summx
            QQBC(3*i+2)=-summy
            QQBC(3*i+3)=-summz
            PPBC(3*i+1)=(p0x*tempmass(i)/eta(Nbody-1)-vx)/tempmass(i)
            PPBC(3*i+2)=(p0y*tempmass(i)/eta(Nbody-1)-vy)/tempmass(i)
            PPBC(3*i+3)=(p0z*tempmass(i)/eta(Nbody-1)-vz)/tempmass(i)
          enddo
          do i=1,Nbody-2
            summx=0.0d0
            summy=0.0d0
            summz=0.0d0
            vx=0.0d0
            vy=0.0d0
            vz=0.0d0
            do j=i+1,Nbody-1
              summx=summx+tempmass(j)/eta(j)*tempx(j)
              summy=summy+tempmass(j)/eta(j)*tempy(j)
              summz=summz+tempmass(j)/eta(j)*tempz(j)
              vx=vx+tempmass(i)/eta(j-1)*tempmx(j)
              vy=vy+tempmass(i)/eta(j-1)*tempmy(j)
              vz=vz+tempmass(i)/eta(j-1)*tempmz(j)
            enddo
            QQBC(3*i+1)=QQ(3*i+1)*eta(i-1)/eta(i)-summx
            QQBC(3*i+2)=QQ(3*i+2)*eta(i-1)/eta(i)-summy
            QQBC(3*i+3)=QQ(3*i+3)*eta(i-1)/eta(i)-summz
            PPBC(3*i+1)=(p0x*tempmass(i)/eta(Nbody-1)+tempmx(i)-vx)/tempmass(i)
            PPBC(3*i+2)=(p0y*tempmass(i)/eta(Nbody-1)+tempmy(i)-vy)/tempmass(i)
            PPBC(3*i+3)=(p0z*tempmass(i)/eta(Nbody-1)+tempmz(i)-vz)/tempmass(i)
          enddo
          do i=Nbody-1,Nbody-1
            QQBC(3*i+1)=eta(Nbody-2)/eta(Nbody-1)*tempx(i)
            QQBC(3*i+2)=eta(Nbody-2)/eta(Nbody-1)*tempy(i)
            QQBC(3*i+3)=eta(Nbody-2)/eta(Nbody-1)*tempz(i)
            PPBC(3*i+1)=(p0x*tempmass(i)/eta(Nbody-1)+tempmx(i))/tempmass(i)
            PPBC(3*i+2)=(p0y*tempmass(i)/eta(Nbody-1)+tempmy(i))/tempmass(i)
            PPBC(3*i+3)=(p0z*tempmass(i)/eta(Nbody-1)+tempmz(i))/tempmass(i)
          enddo
          return
          end
c
c  %#$#^%&%*#&$($^#!^#%!#%^@*@&%!&%!^%!^%!^%!^%!
c
          subroutine movetoastrocentric(Nbody,rmass,QQ,PP,QQBC,PPBC)
c
c   inverts newmovetobarycenter
c
          implicit double precision (a-h,o-z)
c
          dimension rmass(10),QQ(30),PP(30),QQBC(30),PPBC(30)
          dimension eta(0:9),tempmass(0:9),tempx(0:9),tempy(0:9)
          dimension tempz(0:9),tempvx(0:9),tempvy(0:9),tempvz(0:9)
          dimension tempmx(0:9),tempmy(0:9),tempmz(0:9)
c
          totalmass=0.0d0
          do i=1,Nbody
            totalmass=totalmass+rmass(i)
          enddo
c
          p0x=0.0d0
          p0y=0.0d0
          p0z=0.0d0
          do i=0,Nbody-1
            tempmass(i)=rmass(i+1)
c            tempx(i)=QQ(3*i+1)
c            tempy(i)=QQ(3*i+2)
c            tempz(i)=QQ(3*i+3)
c            tempvx(i)=PP(3*i+1)
c            tempvy(i)=PP(3*i+2)
c            tempvz(i)=PP(3*i+3)
          enddo
c
          do i=0,Nbody-1
            summ=0.0d0
            do j=0,i
              summ=summ+tempmass(j)
            enddo
            eta(i)=summ
          enddo
c
          do i=Nbody-1,Nbody-1
            QQBC(3*i+1)=eta(Nbody-1)/eta(Nbody-2)*QQ(3*i+1)
            QQBC(3*i+2)=eta(Nbody-1)/eta(Nbody-2)*QQ(3*i+2)
            QQBC(3*i+3)=eta(Nbody-1)/eta(Nbody-2)*QQ(3*i+3)
            PPBC(3*i+1)=PP(3*i+1)*tempmass(i)
            PPBC(3*i+2)=PP(3*i+2)*tempmass(i)
            PPBC(3*i+3)=PP(3*i+3)*tempmass(i)
          enddo
c
          shiftx=QQ(1)
          shiftvx=PP(1)
          QQBC(1)=QQ(1)-shiftx
          QQBC(4)=QQ(4)-shiftx
          PPBC(1)=PP(1)-shiftvx
          PPBC(4)=PP(4)-shiftvx

          shifty=QQ(2)
          shiftvy=PP(2)
          QQBC(2)=QQ(2)-shifty
          QQBC(5)=QQ(5)-shifty
          PPBC(2)=PP(2)-shiftvy
          PPBC(5)=PP(5)-shiftvy

          shiftz=QQ(3)
          shiftvz=PP(3)
          QQBC(3)=QQ(3)-shiftz
          QQBC(6)=QQ(6)-shiftz
          PPBC(3)=PP(3)-shiftvz
          PPBC(6)=PP(6)-shiftvz

c          tempmx(0)=totalmass*tempvx(0)
c          tempmy(0)=totalmass*tempvy(0)
c          tempmz(0)=totalmass*tempvz(0)
c          p0x=tempmx(0)
c          p0y=tempmy(0)
c          p0z=tempmz(0)
c          do i=1,Nbody-1
c            tempmx(i)=eta(i-1)/eta(i)*tempmass(i)*tempvx(i)
c            tempmy(i)=eta(i-1)/eta(i)*tempmass(i)*tempvy(i)
c            tempmz(i)=eta(i-1)/eta(i)*tempmass(i)*tempvz(i)
c          enddo
c
c          cmx=0.0d0
c          cmy=0.0d0
c          cmz=0.0d0
c          vx=0.0d0
c          vy=0.0d0
c          vz=0.0d0
c          do i=0,0
c            summx=0.0d0
c            summy=0.0d0
c            summz=0.0d0
c            vx=0.0d0
c            vy=0.0d0
c            vz=0.0d0
c            do j=1,Nbody-1
c              summx=summx+tempmass(j)/eta(j)*tempx(j)
c              summy=summy+tempmass(j)/eta(j)*tempy(j)
c              summz=summz+tempmass(j)/eta(j)*tempz(j)
c              vx=vx+tempmass(i)/eta(j-1)*tempmx(j)
c              vy=vy+tempmass(i)/eta(j-1)*tempmy(j)
c              vz=vz+tempmass(i)/eta(j-1)*tempmz(j)
c            enddo
c            QQBC(3*i+1)=-summx
c            QQBC(3*i+2)=-summy
c            QQBC(3*i+3)=-summz
c            PPBC(3*i+1)=(p0x*tempmass(i)/eta(Nbody-1)-vx)/tempmass(i)
c            PPBC(3*i+2)=(p0y*tempmass(i)/eta(Nbody-1)-vy)/tempmass(i)
c            PPBC(3*i+3)=(p0z*tempmass(i)/eta(Nbody-1)-vz)/tempmass(i)
c          enddo
c          do i=1,Nbody-2
c            summx=0.0d0
c            summy=0.0d0
c            summz=0.0d0
c            vx=0.0d0
c            vy=0.0d0
c            vz=0.0d0
c            do j=i+1,Nbody-1
c              summx=summx+tempmass(j)/eta(j)*tempx(j)
c              summy=summy+tempmass(j)/eta(j)*tempy(j)
c              summz=summz+tempmass(j)/eta(j)*tempz(j)
c              vx=vx+tempmass(i)/eta(j-1)*tempmx(j)
c              vy=vy+tempmass(i)/eta(j-1)*tempmy(j)
c              vz=vz+tempmass(i)/eta(j-1)*tempmz(j)
c            enddo
c            QQBC(3*i+1)=QQ(3*i+1)*eta(i-1)/eta(i)-summx
c            QQBC(3*i+2)=QQ(3*i+2)*eta(i-1)/eta(i)-summy
c            QQBC(3*i+3)=QQ(3*i+3)*eta(i-1)/eta(i)-summz
c            PPBC(3*i+1)=(p0x*tempmass(i)/eta(Nbody-1)+tempmx(i)-vx)/tempmass(i)
c            PPBC(3*i+2)=(p0y*tempmass(i)/eta(Nbody-1)+tempmy(i)-vy)/tempmass(i)
c            PPBC(3*i+3)=(p0z*tempmass(i)/eta(Nbody-1)+tempmz(i)-vz)/tempmass(i)
c          enddo

          return
          end


c
c    %%%%%%%%%%%%%%%%%%%%%%%%
c
          integer function ifloor(x)
c
c    returns the nearest integer less than or equal to x
c
          implicit double precision (a-h,o-z)
c
          i=int(x)
c
          if(dble(i).gt.x)i=i-1
c
          ifloor=i

          return
          end
c                                                     
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&     
c                                                 
          double precision function TrueAtoMeanA(Theta,ecc)
c
          implicit double precision (a-h,o-z)
c
          pie=4.0d0*datan(1.0d0)
c
          n=ifloor(Theta/(2.0d0*pie))
          Psi=Theta-2.0d0*pie*dble(n+1)
          Psippi=Psi+pie
          bigE=2.0d0*datan(dsqrt((1.0d0-ecc)/
     @         (1.0d0+ecc))*dtan(Psi/2.0d0))
          sss=0.0d0
          if(Psippi.gt.0.0d0)sss=1.0d0
          if(Psippi.lt.0.0d0)sss=-1.0d0
          TrueAtoMeanA=bigE-ecc*dsin(bigE)+
     @      2.0d0*pie*(dble(n)+0.5d0+sss/2.0d0)
          return
          end
c
c    @#!!@%$@#%$!^&&%&%&$!^#%#!%#DHFHJR!&^&!#%!@$%@!%#!
c   
c   #!!!%(@^%#%*(&(#@#*$^@^$##*#@&(%)
c
          subroutine movetobarycenter(Nbody,rmass,QQ,PP,QQBC,PPBC)
c
          implicit double precision (a-h,o-z)
c
          dimension rmass(10),QQ(30),PP(30),QQBC(30),PPBC(30)
c
          totalmass=0.0d0
          do i=1,Nbody
            totalmass=totalmass+rmass(i)
          enddo
c
          cmx=0.0d0
          cmy=0.0d0
          cmz=0.0d0
          vx=0.0d0
          vy=0.0d0
          vz=0.0d0
          do i=0,Nbody-1
            cmx=cmx+QQ(3*i+1)*rmass(i+1)
            cmy=cmy+QQ(3*i+2)*rmass(i+1)
            cmz=cmz+QQ(3*i+3)*rmass(i+1)
            vx=vx+PP(3*i+1)*rmass(i+1)
            vy=vy+PP(3*i+2)*rmass(i+1)
            vz=vz+PP(3*i+3)*rmass(i+1)
          enddo
          cmx=cmx/totalmass
          cmy=cmy/totalmass
          cmz=cmz/totalmass
          vx=vx/totalmass
          vy=vy/totalmass
          vz=vz/totalmass
c
          do i=0,Nbody-1
            QQBC(3*i+1)=QQ(3*i+1)-cmx
            QQBC(3*i+2)=QQ(3*i+2)-cmy
            QQBC(3*i+3)=QQ(3*i+3)-cmz
            PPBC(3*i+1)=PP(3*i+1)-vx
            PPBC(3*i+2)=PP(3*i+2)-vy
            PPBC(3*i+3)=PP(3*i+3)-vz
          enddo
c
          if(Nbody.le.2)return
c          do i=Nbody-1,Nbody-1
c            QQBC(3*i+1)=QQ(3*i+1)-cmx
c            QQBC(3*i+2)=QQ(3*i+2)-cmy
c            QQBC(3*i+3)=QQ(3*i+3)-cmz
c            PPBC(3*i+1)=PP(3*i+1)!-vx
c            PPBC(3*i+2)=PP(3*i+2)!-vy
c            PPBC(3*i+3)=PP(3*i+3)!-vz
c          enddo
c
c         do j=1,Nbody*3
c           write(*,100)QQ(j),QQBC(j)
c         enddo
c100      format(2(f20.12,3x))
c         write(*,*)'  '
c         do j=1,Nbody*3
c           write(*,100)PP(j),PPBC(j)
c         enddo
c         write(*,*)'----------------------------'

          return
          end
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c   first-order routines
C
C
      SUBROUTINE SOLFID (NR,XOLD,X,Y,N,IRTRN,iGR,rmass,Nbody,tideparm)
      PARAMETER (NSD=11,NMAX=30)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(N),DANG(3),rmass(10),tideparm(10)
      REAL*4 XDES,YDES,AN,HAMEX,YMEAN,YDEVI
      REAL*4 ERRA1,ERRA2,ERRA3
c      REAL*16 QHAMEX
      REAL*8 QHAMEX
      COMMON /INTERP/DHAMEX,DANGEX(3),XDES(5000),YDES(5000),AN,IST
      COMMON /DESTAT/YMEAN(5000),YDEVI(5000)
        IF (MOD(NR,51).EQ.0) THEN
           CALL DHAMIL(N,Y,DHAM,iGR,rmass,Nbody,tideparm)
           CALL DMOMEN(N,Y,DANG,iGR,rmass,Nbody,tideparm)
           ERRH=(DHAM-DHAMEX)/DHAMEX
           ERRA1=(DANG(1)-DANGEX(1))/DANGEX(1)
           ERRA2=(DANG(2)-DANGEX(2))/DANGEX(2)
           ERRA3=(DANG(3)-DANGEX(3))/DANGEX(3)
c           write (6,*) '  error ',nr,ERRH,ERRA1,ERRA2,ERRA3
           ist=ist+1
           XDES(ist)=X
           YDES(ist)=ERRH
           YMEAN(ist)=YMEAN(ist)+ERRH
           YDEVI(ist)=YDEVI(ist)+ERRH**2
        END IF
      RETURN
      END
c
        SUBROUTINE DEQUA(N,X,Y,F,iGR,rmass,Nbody,tideparm)
        IMPLICIT REAL*8 (A-H,O-Z) 
        REAL*8  D(10,10),M(10),rmass(10)
        DIMENSION Y(N),F(N),tideparm(10)
c
c   force equations for first-order system, from E. Hairer
c
c   if iGR = 1, then add GR corrections
c   if iGR = 2, then add tidal corrections
c
        AK=2.95912208286D-4
        AK=(0.01720209895d0)**2
        AC=1.49597870691d11/86400.0d0   !!   converts AU/day into m/sec
        SpLt=299792458.0d0/AC    !!   Speed of light in AU/day
c
        do jj=1,Nbody
          M(jj)=rmass(jj)
        enddo
        DO I=1,Nbody-1     
          I1=3*(I-1)+1
          DO J=I+1,Nbody     
	    J1=3*(J-1)+1
	    D(I,J)=(SQRT((Y(I1)-Y(J1))**2+(Y(I1+1)-Y(J1+1))**2+
     %	      (Y(I1+2)-Y(J1+2))**2))**3
	     D(J,I)=D(I,J)
	  END DO
	END DO
        DO I=1,Nbody   
	  I1=3*(I-1)+1
	  F(Nbody*3+I1)=0.0D0       
	  F(Nbody*3+I1+1)=0.0D0     
          F(Nbody*3+I1+2)=0.0D0     
	  DO J=1,Nbody    
            IF(J.NE.I)THEN
	      J1=3*(J-1)+1
	      F(Nbody*3+I1)=F(Nbody*3+I1)+M(J)*(Y(J1)-Y(I1))/D(I,J)
	      F(Nbody*3+I1+1)=F(Nbody*3+I1+1)+
     @                   M(J)*(Y(J1+1)-Y(I1+1))/D(I,J)
	      F(Nbody*3+I1+2)=F(Nbody*3+I1+2)+
     @                   M(J)*(Y(J1+2)-Y(I1+2))/D(I,J)
            END IF
	  END DO
	  F(Nbody*3+I1+0)=AK*F(Nbody*3+I1+0)
	  F(Nbody*3+I1+1)=AK*F(Nbody*3+I1+1)
	  F(Nbody*3+I1+2)=AK*F(Nbody*3+I1+2)
	  F(I1+0)=Y(Nbody*3+I1+0)
	  F(I1+1)=Y(Nbody*3+I1+1)
	  F(I1+2)=Y(Nbody*3+I1+2)
        END DO
c
        if(iGR.eq.1)then !! add GR correction to the force equations
          rm21=rmass(1)+rmass(2)
          rnu=rmass(1)*rmass(2)/rm21**2
          x21_x=y(4)-y(1)
          x21_y=y(5)-y(2)
          x21_z=y(6)-y(3)
          jndex2=4*3  
          jndex1=3*3  
          v21_x=y(jndex2+1)-y(jndex1+1)
          v21_y=y(jndex2+2)-y(jndex1+2)
          v21_z=y(jndex2+3)-y(jndex1+3)
          r21norm=sqrt(x21_x*x21_x+x21_y*x21_y+x21_z*x21_z)
          r21dot=(x21_x*v21_x+x21_y*v21_y+x21_z*v21_z)/r21norm
          v21dot=(v21_x*v21_x+v21_y*v21_y+v21_z*v21_z)
          xv21dot=(x21_x*v21_x+x21_y*v21_y+x21_z*v21_z)
c 
          agr=(1.0d0+3.0d0*rnu)*v21dot-2.0d0*(2.0d0+
     @         rnu)*AK*rm21/r21norm-1.5d0*rnu*(r21dot)**2
          bgr=-2.0d0*(2.0d0-rnu)*r21dot*r21norm
          cgrx=-(AK/(r21norm**3*SpLt**2))*(agr*x21_x+bgr*v21_x)
          cgry=-(AK/(r21norm**3*SpLt**2))*(agr*x21_y+bgr*v21_y)
          cgrz=-(AK/(r21norm**3*SpLt**2))*(agr*x21_z+bgr*v21_z)
c
          F(jndex1+1)=F(jndex1+1)-rmass(2)*cgrx
          F(jndex1+2)=F(jndex1+2)-rmass(2)*cgry
          F(jndex1+3)=F(jndex1+3)-rmass(2)*cgrz
          F(jndex2+1)=F(jndex2+1)+rmass(1)*cgrx
          F(jndex2+2)=F(jndex2+2)+rmass(1)*cgry
          F(jndex2+3)=F(jndex2+3)+rmass(1)*cgrz
        endif
c
        if(iGR.eq.2)then !! add tidal correction to the force equations
          rm21=rmass(1)+rmass(2)
          rnu=rmass(1)*rmass(2)/rm21**2
          x21_x=y(4)-y(1)
          x21_y=y(5)-y(2)
          x21_z=y(6)-y(3)
          jndex2=4*3  
          jndex1=3*3  
          v21_x=y(jndex2+1)-y(jndex1+1)
          v21_y=y(jndex2+2)-y(jndex1+2)
          v21_z=y(jndex2+3)-y(jndex1+3)
          r21norm=sqrt(x21_x*x21_x+x21_y*x21_y+x21_z*x21_z)
          r21dot=(x21_x*v21_x+x21_y*v21_y+x21_z*v21_z)/r21norm
          v21dot=(v21_x*v21_x+v21_y*v21_y+v21_z*v21_z)
          xv21dot=(x21_x*v21_x+x21_y*v21_y+x21_z*v21_z)
c
c  omega below is hut*(2*pi/period) where hut is given by
c
c     hut=(1.0d0+7.5d0*ecc*ecc+5.625d0*ecc**4+
c       0.3125d0*ecc**6)/((1.0d0+3.0d0*ecc*ecc+
c       3.0d0/8.0d0*ecc**4)*dsqrt((1.0d0-ecc*ecc)**3))
c 
          omega1x=tideparm(5)   ! x-component usually zero
          omega1y=tideparm(6)   ! y-component = omega*sin(i)
          omega1z=tideparm(7)   ! z-component = -omega*cos(i)
          omega2x=tideparm(8)   ! x-component usually zero
          omega2y=tideparm(9)   ! y-component = omega*sin(i)
          omega2z=tideparm(10)  ! z-component = -omega*cos(i)
          omega1norm=sqrt(omega1x**2+omega1y**2+omega1z**2)
          omega2norm=sqrt(omega2x**2+omega2y**2+omega2z**2)
          tidek1=tideparm(1) !  coefficient on the order of 0.005
          tidek2=tideparm(2) !  coefficient on the order of 0.005
          frac1=tideparm(3)  !  this is the radius in AU
          frac2=tideparm(4)  !  this is the radius in AU
c
c   term for star 1.  Note the factor is 6 and not 12
c
          aqp=5.0d0*((omega1x*x21_x+omega1y*x21_y+
     @        omega1z*x21_z)/r21norm)**2-omega1norm**2-
     @        6.0d0*AK*rmass(2)/r21norm**3
          bqp=-2.0d0*(omega1x*x21_x+omega1y*x21_y+omega1z*x21_z)
          qpt=tidek1*((frac1/r21norm)**5) 
          qptx=qpt*(aqp*x21_x+bqp*omega1x)
          qpty=qpt*(aqp*x21_y+bqp*omega1y)
          qptz=qpt*(aqp*x21_z+bqp*omega1z)
c
          F(jndex1+1)=F(jndex1+1)-qptx*(1.0d0+rmass(2)/rmass(1))
          F(jndex1+2)=F(jndex1+2)-qpty*(1.0d0+rmass(2)/rmass(1))
          F(jndex1+3)=F(jndex1+3)-qptz*(1.0d0+rmass(2)/rmass(1))
c
c   don't modify star 2 forces
c
c          F(jndex2+1)=F(jndex2+1)+rmass(1)*qptx/rmass(2)
c          F(jndex2+2)=F(jndex2+2)+rmass(1)*qpty/rmass(2)
c          F(jndex2+3)=F(jndex2+3)+rmass(1)*qptz/rmass(2)
c
c  do the same for star 2.  
c  Note the coefficient is 6 and not 12
c
          aqp=5.0d0*((omega2x*x21_x+omega2y*x21_y+
     @        omega2z*x21_z)/r21norm)**2-omega2norm**2-
     @        6.0d0*AK*rmass(1)/r21norm**3
          bqp=-2.0d0*(omega2x*x21_x+omega2y*x21_y+omega2z*x21_z)
          qpt=tidek2*((frac2/r21norm)**5)
          qptx=qpt*(aqp*x21_x+bqp*omega2x)
          qpty=qpt*(aqp*x21_y+bqp*omega2y)
          qptz=qpt*(aqp*x21_z+bqp*omega2z)
c
c   don't modify the forces for star 1
c
c          F(jndex1+1)=F(jndex1+1)-rmass(1)*qptx
c          F(jndex1+2)=F(jndex1+2)-rmass(1)*qpty
c          F(jndex1+3)=F(jndex1+3)-rmass(1)*qptz
          F(jndex2+1)=F(jndex2+1)+qptx*(1.0d0+rmass(1)/rmass(2))
          F(jndex2+2)=F(jndex2+2)+qpty*(1.0d0+rmass(1)/rmass(2))
          F(jndex2+3)=F(jndex2+3)+qptz*(1.0d0+rmass(1)/rmass(2))

        endif
c
        if(iGR.eq.3)then !! add tidal correction and GR to the force equations
          rm21=rmass(1)+rmass(2)
          rnu=rmass(1)*rmass(2)/rm21**2
          x21_x=y(4)-y(1)
          x21_y=y(5)-y(2)
          x21_z=y(6)-y(3)
          jndex2=4*3  
          jndex1=3*3  
          v21_x=y(jndex2+1)-y(jndex1+1)
          v21_y=y(jndex2+2)-y(jndex1+2)
          v21_z=y(jndex2+3)-y(jndex1+3)
          r21norm=sqrt(x21_x*x21_x+x21_y*x21_y+x21_z*x21_z)
          r21dot=(x21_x*v21_x+x21_y*v21_y+x21_z*v21_z)/r21norm
          v21dot=(v21_x*v21_x+v21_y*v21_y+v21_z*v21_z)
          xv21dot=(x21_x*v21_x+x21_y*v21_y+x21_z*v21_z)
c 
          agr=(1.0d0+3.0d0*rnu)*v21dot-2.0d0*(2.0d0+
     @         rnu)*AK*rm21/r21norm-1.5d0*rnu*(r21dot)**2
          bgr=-2.0d0*(2.0d0-rnu)*r21dot*r21norm
          cgrx=-(AK/(r21norm**3*SpLt**2))*(agr*x21_x+bgr*v21_x)
          cgry=-(AK/(r21norm**3*SpLt**2))*(agr*x21_y+bgr*v21_y)
          cgrz=-(AK/(r21norm**3*SpLt**2))*(agr*x21_z+bgr*v21_z)
c
          F(jndex1+1)=F(jndex1+1)-rmass(2)*cgrx
          F(jndex1+2)=F(jndex1+2)-rmass(2)*cgry
          F(jndex1+3)=F(jndex1+3)-rmass(2)*cgrz
          F(jndex2+1)=F(jndex2+1)+rmass(1)*cgrx
          F(jndex2+2)=F(jndex2+2)+rmass(1)*cgry
          F(jndex2+3)=F(jndex2+3)+rmass(1)*cgrz
c
c  omega below is hut*(2*pi/period) where hut is given by
c
c     hut=(1.0d0+7.5d0*ecc*ecc+5.625d0*ecc**4+
c       0.3125d0*ecc**6)/((1.0d0+3.0d0*ecc*ecc+
c       3.0d0/8.0d0*ecc**4)*dsqrt((1.0d0-ecc*ecc)**3))
c 
          omega1x=tideparm(5)   ! x-component usually zero
          omega1y=tideparm(6)   ! y-component = omega*sin(i)
          omega1z=tideparm(7)   ! z-component = -omega*cos(i)
          omega2x=tideparm(8)   ! x-component usually zero
          omega2y=tideparm(9)   ! y-component = omega*sin(i)
          omega2z=tideparm(10)  ! z-component = -omega*cos(i)
          omega1norm=sqrt(omega1x**2+omega1y**2+omega1z**2)
          omega2norm=sqrt(omega2x**2+omega2y**2+omega2z**2)
          tidek1=tideparm(1) !  coefficient on the order of 0.005
          tidek2=tideparm(2) !  coefficient on the order of 0.005
          frac1=tideparm(3)  !  this is the radius in AU
          frac2=tideparm(4)  !  this is the radius in AU
c
c   term for star 1.  Note the factor is 6 and not 12
c
          aqp=5.0d0*((omega1x*x21_x+omega1y*x21_y+
     @        omega1z*x21_z)/r21norm)**2-omega1norm**2-
     @        6.0d0*AK*rmass(2)/r21norm**3
          bqp=-2.0d0*(omega1x*x21_x+omega1y*x21_y+omega1z*x21_z)
          qpt=tidek1*((frac1/r21norm)**5) 
          qptx=qpt*(aqp*x21_x+bqp*omega1x)
          qpty=qpt*(aqp*x21_y+bqp*omega1y)
          qptz=qpt*(aqp*x21_z+bqp*omega1z)
c
          F(jndex1+1)=F(jndex1+1)-qptx*(1.0d0+rmass(2)/rmass(1))
          F(jndex1+2)=F(jndex1+2)-qpty*(1.0d0+rmass(2)/rmass(1))
          F(jndex1+3)=F(jndex1+3)-qptz*(1.0d0+rmass(2)/rmass(1))
c
c   don't modify star 2 forces
c
c          F(jndex2+1)=F(jndex2+1)+rmass(1)*qptx/rmass(2)
c          F(jndex2+2)=F(jndex2+2)+rmass(1)*qpty/rmass(2)
c          F(jndex2+3)=F(jndex2+3)+rmass(1)*qptz/rmass(2)
c
c  do the same for star 2.  
c  Note the coefficient is 6 and not 12
c
          aqp=5.0d0*((omega2x*x21_x+omega2y*x21_y+
     @        omega2z*x21_z)/r21norm)**2-omega2norm**2-
     @        6.0d0*AK*rmass(1)/r21norm**3
          bqp=-2.0d0*(omega2x*x21_x+omega2y*x21_y+omega2z*x21_z)
          qpt=tidek2*((frac2/r21norm)**5)
          qptx=qpt*(aqp*x21_x+bqp*omega2x)
          qpty=qpt*(aqp*x21_y+bqp*omega2y)
          qptz=qpt*(aqp*x21_z+bqp*omega2z)
c
c   don't modify the forces for star 1
c
c          F(jndex1+1)=F(jndex1+1)-rmass(1)*qptx
c          F(jndex1+2)=F(jndex1+2)-rmass(1)*qpty
c          F(jndex1+3)=F(jndex1+3)-rmass(1)*qptz
          F(jndex2+1)=F(jndex2+1)+qptx*(1.0d0+rmass(1)/rmass(2))
          F(jndex2+2)=F(jndex2+2)+qpty*(1.0d0+rmass(1)/rmass(2))
          F(jndex2+3)=F(jndex2+3)+qptz*(1.0d0+rmass(1)/rmass(2))
c
        endif
c
        if(igr.ge.4)return
c
        RETURN
        END 
c
        SUBROUTINE DHAMIL(N,Y,HAM,iGR,rmass,Nbody,tideparm)
        IMPLICIT REAL*8 (A-H,O-Z) 
        DOUBLE PRECISION  D(6,6),M(6)
        DIMENSION Y(N),rmass(10),tideparm(10)
        common /COUNT/NFCN,NPROB
          nfcn=nfcn+1
        Nprob=4
        write(*,*)'DHAMIL'
        IF (NPROB.EQ.4) THEN
	       AK=2.95912208286D-4
               AK=(0.01720209895d0)**2
           do jj=1,Nbody
             M(jj)=rmass(jj)
           enddo
           DO I=1,Nbody-1    !5
             I1=3*(I-1)+1
             DO J=I+1,Nbody    !6
               J1=3*(J-1)+1
               D(I,J)=SQRT((Y(I1)-Y(J1))**2+(Y(I1+1)-Y(J1+1))**2+
     *	           (Y(I1+2)-Y(J1+2))**2)
	           D(J,I)=D(I,J)
            END DO
          END DO
          HAM=0.0D0
          DO I=1,Nbody    !6
	        I1=Nbody*3+3*(I-1)+1
            HAM=HAM+M(I)*(Y(I1)**2+Y(I1+1)**2+Y(I1+2)**2)
          END DO
          HAM=HAM/2
          POT=0.0D0
          DO I=2,Nbody    !6
            DO J=1,I-1
              POT=POT+M(I)*M(J)/D(I,J)
            END DO
          END DO
          HAM=HAM-AK*POT
        END IF
        RETURN
        END 
c
        SUBROUTINE DMOMEN(N,Y,DANG,iGR,rmass,Nbody,tideparm)
        IMPLICIT REAL*8 (A-H,O-Z) 
        REAL*8 Y(N),DANG(3),M(6),rmass(10),tideparm(10)
c
        do jj=1,Nbody
          M(jj)=rmass(jj)
        enddo
        ND2=N/2
        DO ID=1,3
          DANG(ID)=0.0D0
          DO I=1,6
            ID1=MOD(ID,3)+1
            ID2=MOD(ID+1,3)+1
            I1=(I-1)*3+ID1
            I2=(I-1)*3+ID2
            DANG(ID)=DANG(ID)+M(I)*(Y(I1)*Y(ND2+I2)-Y(I2)*Y(ND2+I1))
          END DO
        END DO
        RETURN
        END 
c
      SUBROUTINE GRKAAD(N,FCN,NSTEP,X,Y,XEND,METH,IOUT,RPAR,rmass,nbody,
     @     posarray,velarray,odetime,zzq,timeinterp,Ndyn,h,iGR,oldzzq,
     @     tideparm)
C ----------------------------------------------------------
      PARAMETER (NDGL=50,NSD=15)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),YH(NDGL),YCS(NDGL),F(NDGL*NSD),Z(NDGL*NSD)
      DIMENSION RPAR(10),METH(10)
      DIMENSION C(NSD),A(NSD,NSD),AP(NSD,NSD),B(NSD),BP(NSD)
      dimension rmass(10),posarray(Ndyn,30),velarray(Ndyn,30)
      dimension zzq(6,60,Ndyn),oldzzq(6,60,Ndyn)
      dimension odetime(Ndyn),timeinterp(6),tideparm(10)
      EXTERNAL FCN

      NS=METH(1)
      NSS=NS*N
      ICOS=METH(2)
      ITSW=METH(3)
      ISYM=METH(4)
      IF (ICOS.NE.0) WRITE (6,*) ' METH(2) = ',ICOS, ' NOT ALLOWED'
      IF (ISYM.NE.0) WRITE (6,*) ' METH(4) = ',ISYM, ' NOT ALLOWED'
 
      do i=1,NSD
        do j=1,NSD
          A(j,i)=0.0d0
          AP(j,i)=0.0d0
          B(i)=0.0d0
          BP(i)=0.0d0
          C(i)=0.0d0
        enddo
      enddo

      CALL GAUSPD(NS,C,B,BP,NSD,A,AP)

      xstart=x

      do jj=1,6
        timeinterp(jj)=c(jj)*h
      enddo

      IF (IOUT.NE.0) CALL SOLFID (0,X,X,Y,N,IRTRN,iGR,rmass,Nbody,
     @      tideparm)
      DO I=1,N
        YCS(I)=0.0D0
      END DO
C ---

      timein=0.0d0
      DO ISTEP=1,NSTEP
        CALL FCN(N,X,Y,F,iGR,rmass,Nbody,tideparm)
        DO I=1,N/2
          FFS=H*F(I)
          YYI=Y(I)
          DO IS=1,6     !NS
              Z(I+(IS-1)*N)=oldzzq(is,i,istep)
c              Z(I+(IS-1)*N)=YYI+C(IS)*FFS
          END DO
        END DO
c
        DO I=N/2+1,N
          FFS=H*F(I)
          YYI=Y(I)
          DO IS=1,6     !NS
            xoff=x  
            timein=H*C(IS)    +x
c            call divdiffonestep(xoff,oldzzq,Nstep,
c     @           timeinterp,I-N/2,
c     $           timein,xout,velout,ii,Ndyn)
            Z(I+(IS-1)*N)=oldzzq(is,i,istep)   !velout
          END DO
        END DO

        odetime(istep)=x  !xstart+dble(istep-1)*h
        do ii=1,N/2
          posarray(istep,ii)=y(ii)
          velarray(istep,ii)=y(ii+N/2)
        enddo

        DYMIN=100.0d0
        DYNO=10.0d0
        ITER=0
        EPS=RPAR(1)
        ITSW=0
        IF (ITSW.EQ.0) THEN
          icount=0
          DO WHILE ((DYNO.LT.DYMIN.AND.DYNO.NE.0.0d0))
            DYMIN=MIN(DYMIN,DYNO)
            CALL ITAAD(FCN,N,NS,NSS,X,Y,YH,NSD,A,AP,C,NDGL,F,Z,H,DYNO,
     @          iGR,rmass,nbody,tideparm)
            icount=icount+1
          END DO
        ELSE
          DO WHILE (DYNO.GT.EPS.AND.ITER.LE.50)
            ITER=ITER+1
            CALL ITAAD(FCN,N,NS,NSS,X,Y,YH,NSD,A,AP,C,NDGL,F,Z,H,DYNO,
     @         iGR,rmass,nbody,tideparm)
          END DO
          if (iter.ge.49) write (6,*) ' no convergence',iter, dyno
          DO IS=1,NS
            ISN=1+(IS-1)*N
            CALL FCN(N,X+H*C(IS),Z(ISN),F(ISN),iGR,rmass,nbody,tideparm)
          END DO
        END IF
c
c        do ii=1,N/2
c          posarray(istep,ii)=y(ii)
c          velarray(istep,ii)=y(ii+N/2)
c        enddo

C
  33    CONTINUE
        X=X+H
c        X=dble(ISTEP)*H
        DO I=1,N
          SUMA=0.0D0
          SUME=0.0D0
          DO IS=1,NS/2
            FF=F(I+(IS-1)*N)+F(I+(NS-IS)*N)
            SUMA=SUMA+B(IS)*FF
            SUME=SUME+BP(IS)*FF
          END DO
          IF (2*(NS/2).NE.NS) THEN
            IS=1+NS/2
            FF=F(I+(IS-1)*N)
            SUMA=SUMA+B(IS)*FF
            SUME=SUME+BP(IS)*FF
          END IF
          TEMP=Y(I)
          YCSI=YCS(I)+H*SUMA
          YI=TEMP+YCSI
          YCSI=YCSI+(TEMP-YI)
          YCSI=YCSI+H*SUME
          Y(I)=YI+YCSI
          YCS(I)=YCSI+(YI-Y(I))
c
          do IS=1,6
c            if(i.le.N/2)zzq(istep,i,IS)=Z(I+(IS-1)*N)  
            if(i.le.N/2)zzq(IS,i,istep)=Z(I+(IS-1)*N)  
          enddo
c          do ii=1,N/2
c            posarray(istep+1,ii)=y(ii)
c            velarray(istep+1,ii)=y(ii+N/2)
c          enddo
        END DO
        IF (IOUT.NE.0) CALL SOLFID (ISTEP,X-H,X,Y,N,IRTRN,iGR,rmass,
     @       Nbody,tideparm)
      END DO

      RETURN
      END
C
      SUBROUTINE ITAAD(FCN,N,NS,NSS,X,Y,YH,NSD,A,AP,C,NDGL,F,Z,H,DYNO,
     @         iGR,rmass,Nbody,tideparm)
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),YH(NDGL),F(NDGL*NS),Z(NDGL*NS)
      DIMENSION C(NS),A(NSD,NS),AP(NSD,NS)
      dimension rmass(10),tideparm(10)
      EXTERNAL FCN
c
      DO IS=1,NS
        ISN=1+(IS-1)*N
        CALL FCN(N,X+H*C(IS),Z(ISN),F(ISN),iGR,rmass,nbody,tideparm)
      END DO
C ---
      DYNO=0.0D0
      DO IS=1,NS
        DO I=1,N
          SSUM=0.0D0
          SUMP=0.0D0
          DO JS=NS,1,-1
            FF=F(I+(JS-1)*N)
            SSUM=SSUM+A(IS,JS)*FF
            SUMP=SUMP+AP(IS,JS)*FF
          END DO
          ISN=I+(IS-1)*N
          ZNEW=Y(I)+H*SSUM+H*SUMP
          DYNO=MAX(DYNO,ABS(ZNEW-Z(ISN)))
          Z(ISN)=ZNEW
        END DO
      END DO
      RETURN
      END
C
      SUBROUTINE GAUSPD(NS,C,B,BP,NSD,A,AP)
      PARAMETER (NSDIM=10)
      IMPLICIT DOUBLE PRECISION (A-H,O-P,R-Z)
c      IMPLICIT REAL*16 (Q)
      IMPLICIT REAL*8 (Q)
      DIMENSION C(NS),B(NS),BP(NS),A(NSD,NS),AP(NSD,NS)
      DIMENSION QC(NSDIM),QB(NSDIM),QA(NSDIM,NSDIM)
      IF (NS.EQ.2) THEN
         B(1)=0.5D0
         A(1,1)= 0.25D0
         A(1,2)=-5.0D0/128.0D0
         A(2,1)= 69.0D0/128.0D0
         A(2,2)= 0.25D0
      END IF
      IF (NS.EQ.3) THEN
         B(1) = 36.0D0/128.0D0
         B(2)=  1.0D0-2*B(1)
         A(1,1) = 18.0D0/128.0D0
         A(1,2) = -5.0D0/128.0D0
         A(1,3) = 1.0D0/128.0D0
         A(2,1) = 38.0D0/128.0D0
         A(2,2) = 28.0D0/128.0D0
         A(2,3) = -3.0D0/128.0D0
         A(3,1) = 34.0D0/128.0D0
         A(3,2) = 61.0D0/128.0D0
         A(3,3) = 18.0D0/128.0D0
      END IF
      IF (NS.EQ.4) THEN
         B(1) = 178.0D0/1024.0D0
         B(2) = 334.0D0/1024.0D0
         A(1,1) = 89.0D0/1024.0D0
         A(1,2) = -27.0D0/1024.0D0
         A(1,3) = 13.0D0/1024.0D0
         A(1,4) = -4.0D0/1024.0D0
         A(2,1) = 193.0D0/1024.0D0
         A(2,2) = 167.0D0/1024.0D0
         A(2,3) = -29.0D0/1024.0D0
         A(2,4) = 7.0D0/1024.0D0
         A(3,1) = 171.0D0/1024.0D0
         A(3,2) = 362.0D0/1024.0D0
         A(3,3) = 167.0D0/1024.0D0
         A(3,4) = -15.0D0/1024.0D0
         A(4,1) = 182.0D0/1024.0D0
         A(4,2) = 321.0D0/1024.0D0
         A(4,3) = 361.0D0/1024.0D0
         A(4,4) = 89.0D0/1024.0D0
      END IF
      IF (NS.EQ.6) THEN
         B(1) = 88.0D0/1024.0D0
         B(2) = 185.0D0/1024.0D0
         B(3) = 239.0D0/1024.0D0
         A(1,1) = 44.0D0/1024.0D0
         A(1,2) = -15.0D0/1024.0D0
         A(1,3) = 10.0D0/1024.0D0
         A(1,4) = -6.0D0/1024.0D0
         A(1,5) = 3.0D0/1024.0D0
         A(1,6) = -1.0D0/1024.0D0
         A(2,1) = 95.0D0/1024.0D0
         A(2,2) = 92.0D0/1024.0D0
         A(2,3) = -21.0D0/1024.0D0
         A(2,4) = 11.0D0/1024.0D0
         A(2,5) = -5.0D0/1024.0D0
         A(2,6) = 1.0D0/1024.0D0
         A(3,1) = 84.0D0/1024.0D0
         A(3,2) = 201.0D0/1024.0D0
         A(3,3) = 120.0D0/1024.0D0
         A(3,4) = -21.0D0/1024.0D0
         A(3,5) = 8.0D0/1024.0D0
         A(3,6) = -2.0D0/1024.0D0
         A(4,1) = 90.0D0/1024.0D0
         A(4,2) = 177.0D0/1024.0D0
         A(4,3) = 261.0D0/1024.0D0
         A(4,4) = 120.0D0/1024.0D0
         A(4,5) = -16.0D0/1024.0D0
         A(4,6) = 3.0D0/1024.0D0
         A(5,1) = 86.0D0/1024.0D0
         A(5,2) = 190.0D0/1024.0D0
         A(5,3) = 229.0D0/1024.0D0
         A(5,4) = 260.0D0/1024.0D0
         A(5,5) = 92.0D0/1024.0D0
         A(5,6) = -7.0D0/1024.0D0
         A(6,1) = 89.0D0/1024.0D0
         A(6,2) = 182.0D0/1024.0D0
         A(6,3) = 245.0D0/1024.0D0
         A(6,4) = 230.0D0/1024.0D0
         A(6,5) = 200.0D0/1024.0D0
         A(6,6) = 44.0D0/1024.0D0
      END IF
      CALL GAUSSQ(NS,QC,QB,NSDIM,QA)
      DO I=1,NS
         C(I)=QC(I)
         BP(I)=QB(I)-B(I)
         DO J=1,NS
            AP(I,J)=QA(I,J)-A(I,J)
         END DO
      END DO
      RETURN
      END
C
      SUBROUTINE GAUSSQ(NS,C,B,NSD,A)
c      IMPLICIT REAL*16 (A-H,O-Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(NS),B(NS),A(NSD,NS)
      IF (NS.EQ.2) THEN
         C(1)=  .211324865405187117745425609749D+00
         C(2)=  .788675134594812882254574390251D+00
         B(1) =  0.5D+00
         A(1,1)=  .250000000000000000000000000000D+00
         A(1,2)= -.386751345948128822545743902510D-01
         A(2,1)=  .538675134594812882254574390251D+00
         A(2,2)=  .250000000000000000000000000000D+00
      END IF
      IF (NS.EQ.3) THEN
         C(1)=  .112701665379258311482073460022D+00
         C(2)=  .500000000000000000000000000000D+00
         C(3)=  .887298334620741688517926539978D+00
         B(1)=  .277777777777777777777777777778D+00
         B(2)=  .444444444444444444444444444444D+00
         A(1,1)=  .138888888888888888888888888889D+00
         A(1,2)= -.359766675249389034563954710966D-01
         A(1,3)=  .978944401530832604958004222948D-02
         A(2,1)=  .300263194980864592438024947213D+00
         A(2,2)=  .222222222222222222222222222222D+00
         A(2,3)= -.224854172030868146602471694354D-01
         A(3,1)=  .267988333762469451728197735548D+00
         A(3,2)=  .480421111969383347900839915541D+00
         A(3,3)=  .138888888888888888888888888889D+00
      END IF
      IF (NS.EQ.4) THEN
         C(1)=  .694318442029737123880267555536D-01
         C(2)=  .330009478207571867598667120448D+00
         C(3)=  .669990521792428132401332879552D+00
         C(4)=  .930568155797026287611973244446D+00
         B(1) = 0.173927422568726928686531974611D+00
         B(2) = 0.326072577431273071313468025389D+00
         A(1,1)=  .869637112843634643432659873055D-01
         A(1,2)= -.266041800849987933133851304770D-01
         A(1,3)=  .126274626894047245150568805746D-01
         A(1,4)= -.355514968579568315691098184957D-02
         A(2,1)=  .188118117499868071650685545087D+00
         A(2,2)=  .163036288715636535656734012694D+00
         A(2,3)= -.278804286024708952241511064190D-01
         A(2,4)=  .673550059453815551539866908570D-02
         A(3,1)=  .167191921974188773171133305525D+00
         A(3,2)=  .353953006033743966537619131808D+00
         A(3,3)=  .163036288715636535656734012694D+00
         A(3,4)= -.141906949311411429641535704762D-01
         A(4,1)=  .177482572254522611843442956461D+00
         A(4,2)=  .313445114741868346798411144814D+00
         A(4,3)=  .352676757516271864626853155866D+00
         A(4,4)=  .869637112843634643432659873055D-01
      END IF
      IF (NS.EQ.6) THEN
         C(1)=  .337652428984239860938492227530D-01
         C(2)=  .169395306766867743169300202490D+00
         C(3)=  .380690406958401545684749139160D+00
         C(4)=  .619309593041598454315250860840D+00
         C(5)=  .830604693233132256830699797510D+00
         C(6)=  .966234757101576013906150777247D+00
         B(1) = 0.856622461895851725201480710863D-01
         B(2) = 0.180380786524069303784916756919D+00
         B(3) = 0.233956967286345523694935171995D+00
         A(1,1)=  .428311230947925862600740355432D-01
         A(1,2)= -.147637259971974124753725910605D-01
         A(1,3)=  .932505070647775119143888450801D-02
         A(1,4)= -.566885804948351190092125641622D-02
         A(1,5)=  .285443331509933513092928583012D-02
         A(1,6)= -.812780171264762112299135651564D-03
         A(2,1)=  .926734914303788631865122917633D-01
         A(2,2)=  .901903932620346518924583784594D-01
         A(2,3)= -.203001022932395859524940805243D-01
         A(2,4)=  .103631562402464237307199458066D-01
         A(2,5)= -.488719292803767146341420376580D-02
         A(2,6)=  .135556105548506177551787075080D-02
         A(3,1)=  .822479226128438738077716511411D-01
         A(3,2)=  .196032162333245006055759781564D+00
         A(3,3)=  .116978483643172761847467585997D+00
         A(3,4)= -.204825277456560976298590118654D-01
         A(3,5)=  .798999189966233579720442148033D-02
         A(3,6)= -.207562578486633419359528915759D-02
         A(4,1)=  .877378719744515067137433602439D-01
         A(4,2)=  .172390794624406967987712335439D+00
         A(4,3)=  .254439495032001621324794183860D+00
         A(4,4)=  .116978483643172761847467585997D+00
         A(4,5)= -.156513758091757022708430246450D-01
         A(4,6)=  .341432357674129871237641994524D-02
         A(5,1)=  .843066851341001107446302003355D-01
         A(5,2)=  .185267979452106975248330960685D+00
         A(5,3)=  .223593811046099099964215226188D+00
         A(5,4)=  .254257069579585109647429252519D+00
         A(5,5)=  .901903932620346518924583784593D-01
         A(5,6)= -.701124524079369066636422067690D-02
         A(6,1)=  .864750263608499346324472067378D-01
         A(6,2)=  .177526353208969968653987471089D+00
         A(6,3)=  .239625825335829035595856428410D+00
         A(6,4)=  .224631916579867772503496287487D+00
         A(6,5)=  .195144512521266716260289347979D+00
         A(6,6)=  .428311230947925862600740355433D-01
      END IF
      RETURN
      END
C
      SUBROUTINE oldGRKAAD(N,FCN,NSTEP,X,Y,XEND,METH,IOUT,RPAR,rmass,
     @    nbody,posarray,velarray,odetime,zzq,timeinterp,Ndyn,h,iGR,
     @    tideparm)
C ----------------------------------------------------------
c  
      PARAMETER (NDGL=50,NSD=15)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),YH(NDGL),YCS(NDGL),F(NDGL*NSD),Z(NDGL*NSD)
      DIMENSION RPAR(10),METH(10)
      DIMENSION C(NSD),A(NSD,NSD),AP(NSD,NSD),B(NSD),BP(NSD)
      dimension rmass(10),posarray(Ndyn,30),velarray(Ndyn,30)
      dimension zzq(6,60,Ndyn),tideparm(10)
      dimension odetime(Ndyn),timeinterp(6)
      EXTERNAL FCN

      NS=METH(1)
      NSS=NS*N
      ICOS=METH(2)
      ITSW=METH(3)
      ISYM=METH(4)
      IF (ICOS.NE.0) WRITE (6,*) ' METH(2) = ',ICOS, ' NOT ALLOWED'
      IF (ISYM.NE.0) WRITE (6,*) ' METH(4) = ',ISYM, ' NOT ALLOWED'
c      H=(XEND-X)/NSTEP
      CALL GAUSPD(NS,C,B,BP,NSD,A,AP)

      xstart=x
      do jj=1,6
        timeinterp(jj)=c(jj)*h
      enddo

      IF (IOUT.NE.0) CALL SOLFID (0,X,X,Y,N,IRTRN,iGR,rmass,Nbody,
     @          tideparm)
      DO I=1,N
        YCS(I)=0.0D0
      END DO
C ---

      DO ISTEP=1,NSTEP
        CALL FCN(N,X,Y,F,iGR,rmass,Nbody,tideparm)
        DO I=1,N
          FFS=H*F(I)
          YYI=Y(I)
          DO IS=1,6     !NS
            Z(I+(IS-1)*N)=YYI+C(IS)*FFS
c            if(i.le.N/2)zzq(istep,i,IS)=Z(I+(IS-1)*N)
          END DO
        END DO
        odetime(istep)=x  !xstart+dble(istep-1)*h
        do ii=1,N/2
          posarray(istep,ii)=y(ii)
          velarray(istep,ii)=y(ii+N/2)
        enddo

        DYMIN=100.D0
        DYNO=10.D0
        ITER=0
        EPS=RPAR(1)
        IF (ITSW.EQ.0) THEN
          DO WHILE (DYNO.LT.DYMIN.AND.DYNO.NE.0.0D0)
            DYMIN=MIN(DYMIN,DYNO)
            CALL ITAAD(FCN,N,NS,NSS,X,Y,YH,NSD,A,AP,C,NDGL,F,Z,H,DYNO,
     @          iGR,rmass,nbody,tideparm)
          END DO
        ELSE
          DO WHILE (DYNO.GT.EPS.AND.ITER.LE.50)
            ITER=ITER+1
            CALL ITAAD(FCN,N,NS,NSS,X,Y,YH,NSD,A,AP,C,NDGL,F,Z,H,DYNO,
     @         iGR,rmass,nbody,tideparm)
          END DO
          if (iter.ge.49) write (6,*) ' no convergence',iter, dyno
          DO IS=1,NS
            ISN=1+(IS-1)*N
            CALL FCN(N,X+H*C(IS),Z(ISN),F(ISN),iGR,rmass,nbody,tideparm)
          END DO
        END IF
c
c        do ii=1,N/2
c          posarray(istep,ii)=y(ii)
c          velarray(istep,ii)=y(ii+N/2)
c        enddo

C
 33         CONTINUE
        X=X+H
c        X=dble(ISTEP)*H
        DO I=1,N
          SUMA=0.0D0
          SUME=0.0D0
          DO IS=1,NS/2
            FF=F(I+(IS-1)*N)+F(I+(NS-IS)*N)
            SUMA=SUMA+B(IS)*FF
            SUME=SUME+BP(IS)*FF
          END DO
          IF (2*(NS/2).NE.NS) THEN
            IS=1+NS/2
            FF=F(I+(IS-1)*N)
            SUMA=SUMA+B(IS)*FF
            SUME=SUME+BP(IS)*FF
          END IF
          TEMP=Y(I)
          YCSI=YCS(I)+H*SUMA
          YI=TEMP+YCSI
          YCSI=YCSI+(TEMP-YI)
          YCSI=YCSI+H*SUME
          Y(I)=YI+YCSI
          YCS(I)=YCSI+(YI-Y(I))
c
          do IS=1,6
c            if(i.le.N/2)zzq(istep,i,IS)=Z(I+(IS-1)*N)  
            if(i.le.N/2)zzq(IS,i,istep)=Z(I+(IS-1)*N)  
          enddo
c          do ii=1,N/2
c            posarray(istep+1,ii)=y(ii)
c            velarray(istep+1,ii)=y(ii+N/2)
c          enddo
        END DO
        IF (IOUT.NE.0) CALL SOLFID (ISTEP,X-H,X,Y,N,IRTRN,iGR,rmass,
     @      Nbody,tideparm)
      END DO

      RETURN
      END
C
