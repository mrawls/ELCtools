          program mkenv
c
c   December 27, 2012
c
c   This code will read the gridloop.opt file
c   and attempt to make an env.com file
c
          implicit double precision (a-h,o-z)
          parameter(Nvmax=99)
          dimension obsparm(18),obv(18),eobv(18)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     %      savestep(Nvmax),Nstepnew(Nvmax)
          character*40 Udatafile,svar(nvmax),Hdatafile,Kdatafile,
     &             RV1file,RV2file
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          character*40 sobv(18),command
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,3),
     @      dwavey(8,3)
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          dimension powercoeff(8,9)
c
          do i=1,Nvmax
            svar(i)='none'
          enddo
c
          call getinput(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,Nref,
     %       rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,sw3,T0,
     $       idraw,iecheck,iidint,iatm,ism1,icnU,icnB,icnV,icnR,icnI,
     @       icnJ,icnH,icnK,iRVfilt,isw1,isw2,isw3,isw4,ilaw,wave,dbolx,
     @       dboly,dwavex,dwavey,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     $       ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,spot2parm,
     %       spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,
     %       temprat,idark1,idark2,isw12,isw13,isw21,isw22,isw23,isw24,
     @       bigI,bigbeta,sw23,sw24,powercoeff, sw25,sw26,sw27,sw28,
     @       sw29,sw30,contam,Tconj,beam1,beam2,isw25,isw26,isw27,isw28,
     &       isw29,isw30,isw31,isw32,isw33,isw34,ocose,osine,omegadot,
     @       contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49)

          call getloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,
     @      Nvmax,Nvar,svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)
c                                                                               
          icnRV1=icnvrt(RV1file(1:2))
          icnRV2=icnvrt(RV2file(1:2))
          isvel1=0
          isvel2=0
          if(icnRV1.ne.430)isvel1=299
          if(icnRV2.ne.430)isvel2=299
c
          Nterms=0
          do 699 i7=1,Nvar
            kkk=icnvrt(svar(i7)(1:2))
            if(kkk.eq.430)then
              go to 749
            else
              Nterms=Nterms+1
            endif
 699      continue
c                                                                               
 749      continue
c
          open(unit=33,file='mkenv.out',status='unknown')
          call writevar(Nvarmax,Nterms,svar,var,isvel1,isvel2,isw30)
          close(33)
c
          end
c
c   **********************************
c
          subroutine getloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,
     @      Nvmax,Nvar,svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)
c
c    November 12, 1999
c
c    This routine will read the input file to set up the optimizer (either
c    the 'grid search' or the Levenburg-Marquardt routine).
c
          implicit double precision (a-h,o-z)
c
c   UPDATE September 11, 2001
c
c   Change the dimensions of obs,eobv,sobv to 9
c
c   UPDATE September 21, 2008
c
c   make the dimension of obv, eobv, and sobv 11
c
c   UPDATE November 17, 2008
c
c   make the dimension of obv, epob, and sobv 17

          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     %        obv(18),eobv(18)
          character*40 Udatafile,svar(Nvmax),Hdatafile,Kdatafile,
     %        RV1file,RV2file,sobv(18)
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          character*1 bell
c
c   UPDATE September 11, 2001
c
c   Change the limit to 9
c
c   UPDATE September 21, 2008
c
c   make the limit 11
c
          do 1 i=1,18
            sobv(i)=' '
            obv(i)=-99.
            eobv(i)=-99.
 1        continue
          Nobv=0
c
          do 2 i=1,Nvmax
            svar(i)='zz'
 2        continue
c
          bell=char(7)
          ios=0
          open(unit=11,file='gridloop.opt',status='old',err=1200,
     @            iostat=ios)
c
          read(11,100)Udatafile
          read(11,100)Bdatafile
          read(11,100)Vdatafile
          read(11,100)Rdatafile
          read(11,100)Idatafile
          read(11,100)Jdatafile
          read(11,100)Hdatafile
          read(11,100)Kdatafile
          read(11,100)RV1file
          read(11,100)RV2file
c
          read(11,*)Nvar
          if(Nvar.eq.0)go to 999
          if(Nvar.gt.Nvmax)then
            write(*,200)bell
            stop
          endif
c
          do 10 i=1,Nvar
            read(11,100)svar(i)
 10       continue
c
          do 20 i=1,Nvar
            read(11,*)vstart(i),vstep(i),Nstep(i)
            var(i)=vstart(i)
 20       continue
c
c  Look here for possible observed parameters
c
 999      read(11,*,end=30,err=30)Nobv
c
          do 25 i=1,Nobv
            read(11,100)sobv(i)
 25       continue
          do 26 i=1,Nobv
            read(11,*)obv(i),eobv(i)
 26       continue

 30       close(11)
c
 1200     if(ios.ne.0)then    ! error in opening the file
            write(*,201)bell
            stop
          endif
c
 100      format(a40)
 200      format(a1,'Error:  too many variables in ''gridloop.opt''')
 201      format(a1,'Error:  unable to open file gridloop.opt')
c
          return
          end
c
c
c   **************************************
c
      INTEGER  FUNCTION  ICNVRT (ASTRING)
c
c    November 12, 1999
c
c    This function was stolen directly out of Peter Stetson's DAOPHOT IIe
c    source code.  The following are his comments:
c
C
C=======================================================================
C
C This little function is supposed to take two ASCII characters and
C express them as an integer in the range 0-(32**2-1) without 
C distinguishing upper and lower case:
C
C AA = Aa = aA = aa = 0, AB = Ab = aB = ab = 1, BA = Ba = bA = ba = 32,
C etc.
C
C Argument
C
C ASTRING is a character string containing two ASCII characters.
C
C=======================================================================
C
      IMPLICIT NONE
      CHARACTER*2 ASTRING
C
C-----------------------------------------------------------------------
C
      ICNVRT=32*MOD(ICHAR(ASTRING(1:1))-1,32)+
     .     MOD(ICHAR(ASTRING(2:2))-1,32)
c
      if(ASTRING(2:2).eq.'0')ICNVRT=ICNVRT+1000
      if(ASTRING(2:2).eq.'1')ICNVRT=ICNVRT+1000
      if(ASTRING(2:2).eq.'2')ICNVRT=ICNVRT+1000
      if(ASTRING(2:2).eq.'3')ICNVRT=ICNVRT+1000
      if(ASTRING(2:2).eq.'4')ICNVRT=ICNVRT+1000
      if(ASTRING(2:2).eq.'5')ICNVRT=ICNVRT+1000
      if(ASTRING(2:2).eq.'6')ICNVRT=ICNVRT+1000
      if(ASTRING(2:2).eq.'7')ICNVRT=ICNVRT+1000
      if(ASTRING(2:2).eq.'8')ICNVRT=ICNVRT+1000
      if(ASTRING(2:2).eq.'9')ICNVRT=ICNVRT+1000
c
      RETURN
C
      END!
C
c    End of stolen code.
c
c
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine recordloopopt(Udatafile,Bdatafile,Vdatafile,
     @      Rdatafile,Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,
     @      RV2file,Nvmax,Nvar,svar,var,vstart,vstep,Nstep,Nobv,sobv,
     @      obv,eobv,vstep1)
c
c    November 15, 1999
c
c    This routine will write an input file to set up the optimizer (either
c    the 'grid search' or the Levenburg-Marquardt routine), based on
c    the final parameters derived by the optimizer.
c
c     UPDATE September 21, 2008
c
c     make the dimension of sobv, eobv, and obv 11
c
c
          implicit double precision (a-h,o-z)
c
          dimension obv(18),eobv(18),vstep1(Nvmax)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax)
          character*40 Udatafile,svar(Nvmax),Hdatafile,Kdatafile,
     %        RV1file,RV2file,sobv(18)
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
c
          ios=0
          open(unit=11,file='gridELC.opt',status='unknown',iostat=ios)
c
          write(11,100)Udatafile
          write(11,100)Bdatafile
          write(11,100)Vdatafile
          write(11,100)Rdatafile
          write(11,100)Idatafile
          write(11,100)Jdatafile
          write(11,100)Hdatafile
          write(11,100)Kdatafile
          write(11,100)RV1file
          write(11,100)RV2file
c
          write(11,500)Nvar
c
          do 10 i=1,Nvar
            write(11,100)svar(i)
            if((svar(i)(1:2).eq.'f1').or.(svar(i)(1:2).eq.'F1'))then
              if(var(i).gt.1.0d0)var(i)=1.0d0
            endif
            if((svar(i)(1:2).eq.'f2').or.(svar(i)(1:2).eq.'F2'))then
              if(var(i).gt.1.0d0)var(i)=1.0d0
            endif
 10       continue
c
          do 20 i=1,Nvar
            if(vstep1(i)-vstart(i).ne.0.0d0)then
              r1=(var(i)-vstart(i))/(vstep1(i)-vstart(i))
              r2=(vstep1(i)-var(i))/(vstep1(i)-vstart(i))
            else
              r1=0.0d0
              r2=0.0d0
            endif
            write(11,501)var(i),vstep(i),Nstep(i),r1,r2,svar(i)(1:2),i
 20       continue
c
          if(Nobv.lt.1)go to 30
          write(11,500)Nobv
          do 25 i=1,Nobv
            write(11,100)sobv(i)
 25       continue
c
          do 26 i=1,Nobv
            write(11,502)obv(i),eobv(i)
 26       continue
c
 30       close(11)
c
 100      format(a40)
 500      format(i2)
 501      format(f18.10,1x,f18.10,1x,i7,2x,2(f10.7,1x),1x,a2,1x,i2)
 502      format(f16.9,3x,f15.10)
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine writevar(Nvarmax,Nterms,svar,var,isvel1,isvel2,
     @        isw30)
c
          implicit double precision (a-h,o-z)
c
          character*80 command,line
          character*40 svar(Nvarmax)
          character*19 string19
          character*30 outstring
          character*40 pwd
c
          dimension var(Nvarmax)
c
          command='pwd > mkenv.pwd'
          call system(command)
c
          open(unit=20,file='mkenv.pwd',status='old')
          read(20,55)line
          close(20)
c
 55       format(a80)
c
          lll=lnblnk(line)
c
          islash=1
          do i=1,lll
           if(line(i:i).eq.'/')islash=i
          enddo
c 
c          write(*,*)line(islash+1:lll)
c
c          write(*,*)'enter the working directory'
c          read(*,50)pwd
c
          pwd=line(islash+1:lll)
          lll=lnblnk(pwd)
          write(*,56)pwd
 56       format(a)

          igg=1
          write(*,*)'enter 1 for >, and 2 for >>'
          read(*,*)igg
c
 50       format(a40)
c

          iline=0
          iwrite=0
          ilength=2
          iargper=0              !flag for ecosw being assigned.
          do 10 i=1,Nterms  !Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
            if((kk.eq.1016))then          !k1, tag a1    
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'apsidal_k1'
              else
                write(33,501)i+2,pwd(1:lll),'apsidal_k1'
                isw29=99
              endif
            endif
c                                                
            if((kk.eq.1017))then          !k2, tag a2    
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'apsidal_k2'
              else
                write(33,501)i+2,pwd(1:lll),'apsidal_k2'
                isw29=99
              endif
            endif
c
            if((kk.eq.450))then          !ecos(omega), tag oc
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'ecos_omega'
              else
                write(33,501)i+2,pwd(1:lll),'ecos_omega'
                isw29=99
              endif
            endif
c
            if((kk.eq.466))then  !esin(omega), tag os
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'esin_omega'
              else
              write(33,501)i+2,pwd(1:lll),'esin_omega'
              isw29=99
              endif
            endif
c
            if(kk.eq.1208)then  !Tgrav1, tag g1
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'Tgrav1'
              else
              write(33,501)i+2,pwd(1:lll),'Tgrav1'
              endif
            endif

            if(kk.eq.110)then  !omega_dot, tag do
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'omega_dot'
              else
              write(33,501)i+2,pwd(1:lll),'omega_dot'
              endif
            endif
c
            if(kk.eq.1209)then  !Tgrav2, tag g2
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'Tgrav2'
              else
              write(33,501)i+2,pwd(1:lll),'Tgrav2'
              endif
            endif
c
            if(kk.eq.627)then          !tertperiod, tag tt
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'tertP'
              else
              write(33,501)i+2,pwd(1:lll),'tertP'
              endif
            endif
c
            if(kk.eq.628)then          !tertT0, tag tu
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'tertT0'
              else
              write(33,501)i+2,pwd(1:lll),'tertT0'
              endif
            endif
c
            if(kk.eq.617)then          !tertT0, tag tj
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'tertconj'
              else
              write(33,501)i+2,pwd(1:lll),'tertconj'
              endif
            endif
c
            if(kk.eq.629)then          !tertecos, tag tv
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'terte_cos'
              else
              write(33,501)i+2,pwd(1:lll),'terte_cos'
              endif
            endif
c
            if(kk.eq.630)then          !tertecos, tag tw
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'terte_sin'
              else
              write(33,501)i+2,pwd(1:lll),'terte_sin'
              endif
            endif
c
            if(kk.eq.631)then          !tertincl, tag tx
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'tertincl'
              else
              write(33,501)i+2,pwd(1:lll),'tertincl'
              endif
            endif
c

            if(kk.eq.632)then          !tertOmega, tag ty
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'tertOmega'
              else
              write(33,501)i+2,pwd(1:lll),'tertOmega'
              endif
            endif
c
            if(kk.eq.633)then          !tertQ, tag tz
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'tertQ'
              else
              write(33,501)i+2,pwd(1:lll),'tertQ'
              endif
            endif
c
            if(kk.eq.78)then          !contam, tag co
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'contam'
              else
              write(33,501)i+2,pwd(1:lll),'contam'
              endif
            endif
c
            if(kk.eq.1591)then          !seasonal contam S0, tag s0
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'s0'
              else
              write(33,501)i+2,pwd(1:lll),'s0'
              endif
            endif
c
            if(kk.eq.1592)then          !seasonal contam S1, tag s1
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'s1'
              else
              write(33,501)i+2,pwd(1:lll),'s1'
              endif
            endif
c
            if(kk.eq.1593)then          !seasonal contam S2, tag s2
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'s2'
              else
              write(33,501)i+2,pwd(1:lll),'s2'
              endif
            endif
c
            if(kk.eq.1594)then          !seasonal contam S3, tag s3
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'s3'
              else
              write(33,501)i+2,pwd(1:lll),'s3'
              endif
            endif
c
            if(kk.eq.1144)then          !beam1, tag e1
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'beam1'
              else
              write(33,501)i+2,pwd(1:lll),'beam1'
              endif
            endif
c
            if(kk.eq.1145)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'beam2'
              else
              write(33,501)i+2,pwd(1:lll),'beam2'
              endif
            endif
c
            if(kk.eq.492)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'primmass'
              else
              write(33,501)i+2,pwd(1:lll),'primmass'
              endif
            endif
c
            if(kk.eq.497)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'primrad'
              else
              write(33,501)i+2,pwd(1:lll),'primrad'
              endif
            endif
c
            if(kk.eq.612)then    ! temprat
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'temprat'
              else
              write(33,501)i+2,pwd(1:lll),'temprat'
              endif
            endif
c
            if(kk.eq.490)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'primK'
              else
              write(33,501)i+2,pwd(1:lll),'primK'
              endif
            endif
c
            if(kk.eq.544)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'ratrad'
              else
              write(33,501)i+2,pwd(1:lll),'ratrad'
              endif
            endif
c
            if(kk.eq.1528)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'frac1'
              else
              write(33,501)i+2,pwd(1:lll),'frac1'
              endif
            endif
c
            if(kk.eq.100)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'density'
              else
              write(33,501)i+2,pwd(1:lll),'density'
              endif
            endif
c
            if(kk.eq.1529)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'frac2'
              else
              write(33,501)i+2,pwd(1:lll),'frac2'
              endif
            endif
c
            if(kk.eq.1368)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'al1'
              else
              write(33,501)i+2,pwd(1:lll),'al1'
              endif
            endif
c
            if(kk.eq.1369)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'al2'
              else
              write(33,501)i+2,pwd(1:lll),'al2'
              endif
            endif
c
            if(kk.eq.111)then   !ecosw, use string dphi
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'dphi'
              else
              write(33,501)i+2,pwd(1:lll),'dphi'
              iargper=100
              endif
            endif
c
            if(kk.eq.610)then           !Tconj string
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'Tconj'
              else
              write(33,501)i+2,pwd(1:lll),'Tconj'
              endif
            endif
c
            if(kk.eq.1623)then           !T0 string
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'T0'
              else
              write(33,501)i+2,pwd(1:lll),'T0'
              endif
            endif
c
            if(kk.eq.484)then      !period string
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'Period'
              else
              write(33,501)i+2,pwd(1:lll),'Period'
              endif
            endif
c
            if(kk.eq.269)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'incl'
              else
              write(33,501)i+2,pwd(1:lll),'incl'
              endif
            endif
c
            if(kk.eq.8)then          !bigI
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'bigI'
              else
              write(33,501)i+2,pwd(1:lll),'bigI'
              endif
            endif
c
            if(kk.eq.1)then          !bigbeta
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'bigbeta'
              else
              write(33,501)i+2,pwd(1:lll),'bigbeta'
              endif
            endif
c
            if(kk.eq.384)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'Q'
              else
              write(33,501)i+2,pwd(1:lll),'Q'
              endif
            endif
c
            if(kk.eq.1176)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'fill1'
              else
              write(33,501)i+2,pwd(1:lll),'fill1'
              endif
            endif
c
            if(kk.eq.1177)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'fill2'
              else
              write(33,501)i+2,pwd(1:lll),'fill2'
              endif
            endif
c
            if(kk.eq.552)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'rin'
              else
              write(33,501)i+2,pwd(1:lll),'rin'
              endif
            endif
c
            if(kk.eq.1464)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'omega1'
              else
              write(33,501)i+2,pwd(1:lll),'omega1'
              endif
            endif
c
            if(kk.eq.1465)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'omega2'
              else
              write(33,501)i+2,pwd(1:lll),'omega2'
              endif
            endif
c
            if(kk.eq.558)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'rout'
              else
              write(33,501)i+2,pwd(1:lll),'rout'
              endif
            endif
c
            if(kk.eq.1626)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'Teff3'
              else
              write(33,501)i+2,pwd(1:lll),'Teff3'
              endif
            endif
c
            if(kk.eq.1210)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'g3'
              else
              write(33,501)i+2,pwd(1:lll),'g3'
              endif
            endif
c
            if(kk.eq.36)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'beta'
              else
              write(33,501)i+2,pwd(1:lll),'beta'
              endif
            endif
c
            if(kk.eq.1624)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'Teff1'
              else
              write(33,501)i+2,pwd(1:lll),'Teff1'
              endif
            endif
c
            if(kk.eq.1625)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'Teff2'
              else
              write(33,501)i+2,pwd(1:lll),'Teff2'
              endif
            endif
c
            if(kk.eq.744)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'xi'
              else
              write(33,501)i+2,pwd(1:lll),'xi'
              endif
            endif
c
            if(kk.eq.611)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'Tdisk'
              else
              write(33,501)i+2,pwd(1:lll),'Tdisk'
              endif
            endif
c
            if(kk.eq.375)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'rLx'
              else
              write(33,501)i+2,pwd(1:lll),'rLx'
              endif
            endif
c
            if(kk.eq.580)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'separ'
              else
              write(33,501)i+2,pwd(1:lll),'separ'
              endif
            endif
c
            if(kk.eq.576)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'SA3'
              else
              write(33,501)i+2,pwd(1:lll),'SA3'
              endif
            endif
c
            if(kk.eq.498)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'pshift'
              else
              write(33,501)i+2,pwd(1:lll),'pshift'
              endif
            endif
c
            if(kk.eq.130)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'ecc'
              else
              write(33,501)i+2,pwd(1:lll),'ecc'
              endif
            endif
c
            if(kk.eq.17)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'argper'
              else
              write(33,501)i+2,pwd(1:lll),'argper'
              endif
            endif
c
            if(kk.eq.1048)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'TF_spot1_star1'
              else
              write(33,501)i+2,pwd(1:lll),'TF_spot1_star1'
              endif
            endif
c
            if(kk.eq.1049)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'lat_spot1_star1'
              else
              write(33,501)i+2,pwd(1:lll),'lat_spot1_star1'
              endif
            endif
c
            if(kk.eq.1050)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'lon_spot1_star1'
              else
              write(33,501)i+2,pwd(1:lll),'lon_spot1_star1'
              endif
            endif
c
            if(kk.eq.1051)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'rad_spot1_star1'
              else
              write(33,501)i+2,pwd(1:lll),'rad_spot1_star1'
              endif
            endif
c
            if(kk.eq.1052)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'TF_spot2_star1'
              else
              write(33,501)i+2,pwd(1:lll),'TF_spot2_star1'
              endif
            endif
c
            if(kk.eq.1053)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'lat_spot2_star1'
              else
              write(33,501)i+2,pwd(1:lll),'lat_spot2_star1'
              endif
            endif
c
            if(kk.eq.1054)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'lon_spot2_star1'
              else
              write(33,501)i+2,pwd(1:lll),'lon_spot2_star1'
              endif
            endif
c
            if(kk.eq.1055)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'rad_spot2_star1'
              else
              write(33,501)i+2,pwd(1:lll),'rad_spot2_star1'
              endif
            endif
c
            if(kk.eq.1080)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'TF_spot1_star2'
              else
              write(33,501)i+2,pwd(1:lll),'TF_spot1_star2'
              endif
            endif
c
            if(kk.eq.1081)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'lat_spot1_star2'
              else
              write(33,501)i+2,pwd(1:lll),'lat_spot1_star2'
              endif
            endif
c
            if(kk.eq.1082)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'lon_spot1_star2'
              else
              write(33,501)i+2,pwd(1:lll),'lon_spot1_star2'
              endif
            endif
c
            if(kk.eq.1083)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'rad_spot1_star2'
              else
              write(33,501)i+2,pwd(1:lll),'rad_spot1_star2'
              endif
            endif
c
            if(kk.eq.1084)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'TF_spot2_star2'
              else
              write(33,501)i+2,pwd(1:lll),'TF_spot2_star2'
              endif
            endif
c
            if(kk.eq.1085)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'lat_spot2_star2'
              else
              write(33,501)i+2,pwd(1:lll),'lat_spot2_star2'
              endif
            endif
c
            if(kk.eq.1086)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'lon_spot2_star2'
              else
              write(33,501)i+2,pwd(1:lll),'lon_spot2_star2'
              endif
            endif
c
            if(kk.eq.1087)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'rad_spot2_star2'
              else
              write(33,501)i+2,pwd(1:lll),'rad_spot2_star2'
              endif
            endif
c
            if(kk.eq.1112)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'TF_spot1_disk'
              else
              write(33,501)i+2,pwd(1:lll),'TF_spot1_disk'
              endif
            endif
c
            if(kk.eq.1113)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'azi_spot1_disk'
              else
              write(33,501)i+2,pwd(1:lll),'azi_spot1_disk'
              endif
            endif
c
            if(kk.eq.1114)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'cut_spot1_disk'
              else
              write(33,501)i+2,pwd(1:lll),'cut_spot1_disk'
              endif
            endif
c
            if(kk.eq.1115)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'wid_spot1_disk'
              else
              write(33,501)i+2,pwd(1:lll),'wid_spot1_disk'
              endif
            endif
c
            if(kk.eq.1116)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'TF_spot2_disk'
              else
              write(33,501)i+2,pwd(1:lll),'TF_spot2_disk'
              endif
            endif
c
            if(kk.eq.1117)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'azi_spot2_disk'
              else
              write(33,501)i+2,pwd(1:lll),'azi_spot2_disk'
              endif
            endif
c
            if(kk.eq.1118)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'cut_spot2_disk'
              else
              write(33,501)i+2,pwd(1:lll),'cut_spot2_disk'
              endif
            endif
c
            if(kk.eq.1119)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'wid_spot2_disk'
              else
              write(33,501)i+2,pwd(1:lll),'wid_spot2_disk'
              endif
            endif
c
            if(kk.eq.1752)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x1_U'
              else
              write(33,501)i+2,pwd(1:lll),'x1_U'
              endif
            endif
c
            if(kk.eq.1753)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x1_B'
              else
              write(33,501)i+2,pwd(1:lll),'x1_B'
              endif
            endif
c
            if(kk.eq.1754)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x1_V'
              else
              write(33,501)i+2,pwd(1:lll),'x1_V'
              endif
            endif
c
            if(kk.eq.1755)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x1_R'
              else
              write(33,501)i+2,pwd(1:lll),'x1_R'
              endif
            endif
c
            if(kk.eq.1756)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x1_I'
              else
              write(33,501)i+2,pwd(1:lll),'x1_I'
              endif
            endif
c
            if(kk.eq.1757)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x1_J'
              else
              write(33,501)i+2,pwd(1:lll),'x1_J'
              endif
            endif
c
            if(kk.eq.1758)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x1_H'
              else
              write(33,501)i+2,pwd(1:lll),'x1_H'
              endif
            endif
c
            if(kk.eq.1759)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x1_K'
              else
              write(33,501)i+2,pwd(1:lll),'x1_K'
              endif
            endif
c
            if(kk.eq.1816)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x2_U'
              else
              write(33,501)i+2,pwd(1:lll),'x2_U'
              endif
            endif
c
            if(kk.eq.1817)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x2_B'
              else
              write(33,501)i+2,pwd(1:lll),'x2_B'
              endif
            endif
c
            if(kk.eq.1818)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x2_V'
              else
              write(33,501)i+2,pwd(1:lll),'x2_V'
              endif
            endif
c
            if(kk.eq.1819)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x2_R'
              else
              write(33,501)i+2,pwd(1:lll),'x2_R'
              endif
            endif
c
            if(kk.eq.1820)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x2_I'
              else
              write(33,501)i+2,pwd(1:lll),'x2_I'
              endif
            endif
c
            if(kk.eq.1821)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x2_J'
              else
              write(33,501)i+2,pwd(1:lll),'x2_J'
              endif
            endif
c
            if(kk.eq.1822)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x2_H'
              else
              write(33,501)i+2,pwd(1:lll),'x2_H'
              endif
            endif
c
            if(kk.eq.1823)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'x2_K'
              else
              write(33,501)i+2,pwd(1:lll),'x2_K'
              endif
            endif
c
            if(kk.eq.1784)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y1_U'
              else
              write(33,501)i+2,pwd(1:lll),'y1_U'
              endif
            endif
c
            if(kk.eq.1785)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y1_B'
              else
              write(33,501)i+2,pwd(1:lll),'y1_B'
              endif
            endif
c
            if(kk.eq.1786)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y1_V'
              else
              write(33,501)i+2,pwd(1:lll),'y1_V'
              endif
            endif
c
            if(kk.eq.1787)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y1_R'
              else
              write(33,501)i+2,pwd(1:lll),'y1_R'
              endif
            endif
c
            if(kk.eq.1788)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y1_I'
              else
              write(33,501)i+2,pwd(1:lll),'y1_I'
              endif
            endif
c
            if(kk.eq.1789)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y1_J'
              else
              write(33,501)i+2,pwd(1:lll),'y1_J'
              endif
            endif
c
            if(kk.eq.1790)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y1_H'
              else
              write(33,501)i+2,pwd(1:lll),'y1_H'
              endif
            endif
c
            if(kk.eq.1791)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y1_K'
              else
              write(33,501)i+2,pwd(1:lll),'y1_K'
              endif
            endif
c
            if(kk.eq.1720)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y2_U'
              else
              write(33,501)i+2,pwd(1:lll),'y2_U'
              endif
            endif
c
            if(kk.eq.1721)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y2_B'
              else
              write(33,501)i+2,pwd(1:lll),'y2_B'
              endif
            endif
c
            if(kk.eq.1722)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y2_V'
              else
              write(33,501)i+2,pwd(1:lll),'y2_V'
              endif
            endif
c
            if(kk.eq.1723)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y2_R'
              else
              write(33,501)i+2,pwd(1:lll),'y2_R'
              endif
            endif
c
            if(kk.eq.1724)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y2_I'
              else
              write(33,501)i+2,pwd(1:lll),'y2_I'
              endif
            endif
c
            if(kk.eq.1725)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y2_J'
              else
              write(33,501)i+2,pwd(1:lll),'y2_J'
              endif
            endif
c
            if(kk.eq.1726)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y2_H'
              else
              write(33,501)i+2,pwd(1:lll),'y2_H'
              endif
            endif
c
            if(kk.eq.1727)then
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'y2_K'
              else
              write(33,501)i+2,pwd(1:lll),'y2_K'
              endif
            endif
c
            if(kk.eq.649)then      !tag uj, P2tconj
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P2tconj'
              else
              write(33,501)i+2,pwd(1:lll),'P2tconj'
              endif
            endif
c
            if(kk.eq.659)then      !tag ut, P2period
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P2period'
              else
              write(33,501)i+2,pwd(1:lll),'P2period'
              endif
            endif
c
            if(kk.eq.660)then      !tag uu, P2T0
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P2T0'
              else
              write(33,501)i+2,pwd(1:lll),'P2T0'
              endif
            endif
c
            if(kk.eq.661)then      !tag uv, P2ecos
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P2ecos'
              else
              write(33,501)i+2,pwd(1:lll),'P2ecos'
              endif
            endif
c
            if(kk.eq.662)then      !tag uw, P2esin
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P2esin'
              else
              write(33,501)i+2,pwd(1:lll),'P2esin'
              endif
            endif
c
            if(kk.eq.663)then      !tag ux, P2incl
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P2incl'
              else
              write(33,501)i+2,pwd(1:lll),'P2incl'
              endif
            endif
c
            if(kk.eq.664)then      !tag uy, P2Omega
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P2Omega'
              else
              write(33,501)i+2,pwd(1:lll),'P2Omega'
              endif
            endif
c
            if(kk.eq.665)then      !tag uz, P2Q
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P2Q'
              else
              write(33,501)i+2,pwd(1:lll),'P2Q'
              endif
            endif
c
            if(kk.eq.641)then      !tag ub, P2ratrad
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P2ratrad'
              else
              write(33,501)i+2,pwd(1:lll),'P2ratrad'
              endif
            endif
c
c  
c
            if(kk.eq.681)then      !tag vj, P3tconj
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P3tconj'
              else
              write(33,501)i+2,pwd(1:lll),'P3tconj'
              endif
            endif
c
            if(kk.eq.691)then      !tag vt, P3period
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P3period'
              else
              write(33,501)i+2,pwd(1:lll),'P3period'
              endif
            endif
c
            if(kk.eq.692)then      !tag vu, P3T0
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P3T0'
              else
              write(33,501)i+2,pwd(1:lll),'P3T0'
              endif
            endif
c
            if(kk.eq.693)then      !tag vv, P3ecos
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P3ecos'
              else
              write(33,501)i+2,pwd(1:lll),'P3ecos'
              endif
            endif
c
            if(kk.eq.694)then      !tag vw, P3esin
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P3esin'
              else
              write(33,501)i+2,pwd(1:lll),'P3esin'
              endif
            endif
c
            if(kk.eq.695)then      !tag vx, P3incl
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P3incl'
              else
              write(33,501)i+2,pwd(1:lll),'P3incl'
              endif
            endif
c
            if(kk.eq.696)then      !tag vy, P3Omega
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P3Omega'
              else
              write(33,501)i+2,pwd(1:lll),'P3Omega'
              endif
            endif
c
            if(kk.eq.697)then      !tag vz, P3Q
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P3Q'
              else
              write(33,501)i+2,pwd(1:lll),'P3Q'
              endif
            endif
c
            if(kk.eq.673)then      !tag vb, P3ratrad
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P3ratrad'
              else
              write(33,501)i+2,pwd(1:lll),'P3ratrad'
              endif
            endif
c
c QQQ
c
            if(kk.eq.713)then      !tag wj, P4tconj
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P4tconj'
              else
              write(33,501)i+2,pwd(1:lll),'P4tconj'
              endif
            endif
c
            if(kk.eq.723)then      !tag wt, P4period
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P4period'
              else
              write(33,501)i+2,pwd(1:lll),'P4period'
              endif
            endif
c
            if(kk.eq.724)then      !tag wu, P4T0
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P4T0'
              else
              write(33,501)i+2,pwd(1:lll),'P4T0'
              endif
            endif
c
            if(kk.eq.725)then      !tag wv, P4ecos
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P4ecos'
              else
              write(33,501)i+2,pwd(1:lll),'P4ecos'
              endif
            endif
c
            if(kk.eq.726)then      !tag ww, P4esin
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P4esin'
              else
              write(33,501)i+2,pwd(1:lll),'P4esin'
              endif
            endif
c
            if(kk.eq.727)then      !tag wx, P4incl
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P4incl'
              else
              write(33,501)i+2,pwd(1:lll),'P4incl'
              endif
            endif
c
            if(kk.eq.728)then      !tag wy, P4Omega
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P4Omega'
              else
              write(33,501)i+2,pwd(1:lll),'P4Omega'
              endif
            endif
c
            if(kk.eq.729)then      !tag wz, P4Q
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P4Q'
              else
              write(33,501)i+2,pwd(1:lll),'P4Q'
              endif
            endif
c
            if(kk.eq.705)then      !tag wb, P4ratrad
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'P4ratrad'
              else
              write(33,501)i+2,pwd(1:lll),'P4ratrad'
              endif
            endif

 10       continue
c
c
c   Update July 29, 2005
c
c   If the variable ecosw was assigned, then also print out the argument
c   of periastron.
c
            if(iargper.gt.0)then
              i=i+1
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'argper'
              else
              write(33,501)i+2,pwd(1:lll),'argper'
              endif
            endif
c
            if(isw29.gt.0)then
              i=i+1
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'ecc'
              else
                write(33,501)i+2,pwd(1:lll),'ecc'
              endif
              i=i+1
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'argper'
              else
                write(33,501)i+2,pwd(1:lll),'argper'
              endif
            endif
c
c UPDATE June 10, 2003
c
c Add two more digits to gamma, use string 15
c
          if(isvel1.gt.0)then
              i=i+1
              if(igg.eq.1)then
                write(33,500)i+2,pwd(1:lll),'gam1'
              else
                write(33,501)i+2,pwd(1:lll),'gam1'
              endif
          endif
c
          if(isvel2.gt.0)then
              i=i+1
              write(33,500)i+2,pwd(1:lll),'gam2'
          endif
c
          if(igg.eq.1)then
            write(33,502)3,pwd(1:lll),'m1'
          else
            write(33,503)3,pwd(1:lll),'m1'
          endif
          if(igg.eq.1)then
            write(33,502)4,pwd(1:lll),'m2'
          else
            write(33,503)4,pwd(1:lll),'m2'
          endif
          if(igg.eq.1)then
            write(33,502)5,pwd(1:lll),'r1'
          else
            write(33,503)5,pwd(1:lll),'r1'
          endif
          if(igg.eq.1)then
            write(33,502)6,pwd(1:lll),'r2'
          else
            write(33,503)6,pwd(1:lll),'r2'
          endif
          if(igg.eq.1)then
            write(33,502)7,pwd(1:lll),'r3'
          else
            write(33,503)7,pwd(1:lll),'r3'
          endif
          if(igg.eq.1)then
            write(33,502)8,pwd(1:lll),'g1'
          else
            write(33,503)8,pwd(1:lll),'g1'
          endif
          if(igg.eq.1)then
            write(33,502)9,pwd(1:lll),'g2'
          else
            write(33,503)9,pwd(1:lll),'g2'
          endif
          if(igg.eq.1)then
            write(33,502)12,pwd(1:lll),'separ'
          else
            write(33,503)12,pwd(1:lll),'separ'
          endif
          if(igg.eq.1)then
            write(33,502)13,pwd(1:lll),'k1'
          else
            write(33,503)13,pwd(1:lll),'k1'
          endif
          if(igg.eq.1)then
            write(33,502)14,pwd(1:lll),'k2'
          else
            write(33,503)14,pwd(1:lll),'k2'
          endif
          if(igg.eq.1)then
            write(33,502)17,pwd(1:lll),'vrot1'
          else
            write(33,503)17,pwd(1:lll),'vrot1'
          endif
          if(igg.eq.1)then
            write(33,502)18,pwd(1:lll),'vrot2'
          else
            write(33,503)18,pwd(1:lll),'vrot2'
          endif
c
          if(isw30.ge.1)then
            if(igg.eq.1)then
              write(33,602)9,pwd(1:lll),'period2day'
            else
              write(33,603)9,pwd(1:lll),'period2day'
            endif
            if(igg.eq.1)then
              write(33,602)10,pwd(1:lll),'m2earth'
            else
              write(33,603)10,pwd(1:lll),'m2earth'
            endif
            if(igg.eq.1)then
              write(33,602)11,pwd(1:lll),'r2earth'
            else
              write(33,603)11,pwd(1:lll),'r2earth'
            endif
            if(igg.eq.1)then
              write(33,602)12,pwd(1:lll),'den2'
            else
              write(33,603)12,pwd(1:lll),'den2'
            endif
            if(igg.eq.1)then
              write(33,602)13,pwd(1:lll),'a2AU'
            else
              write(33,603)13,pwd(1:lll),'a2AU'
            endif
            if(igg.eq.1)then
              write(33,602)14,pwd(1:lll),'ecc2'
            else
              write(33,603)14,pwd(1:lll),'ecc2'
            endif
            if(igg.eq.1)then
              write(33,602)15,pwd(1:lll),'arg2deg'
            else
              write(33,603)15,pwd(1:lll),'arg2deg'
            endif
            if(igg.eq.1)then
              write(33,602)16,pwd(1:lll),'incl2deg'
            else
              write(33,603)16,pwd(1:lll),'incl2deg'
            endif
            if(igg.eq.1)then
              write(33,602)17,pwd(1:lll),'Omega2deg'
            else
              write(33,603)17,pwd(1:lll),'Omega2deg'
            endif
            if(igg.eq.1)then
              write(33,602)18,pwd(1:lll),'Tconj2'
            else
              write(33,603)18,pwd(1:lll),'Tconj2'
            endif
            if(igg.eq.1)then
              write(33,602)19,pwd(1:lll),'Tperi2'
            else
              write(33,603)19,pwd(1:lll),'Tperi2'
            endif
            if(igg.eq.1)then
              write(33,602)20,pwd(1:lll),'mutual2deg'
            else
              write(33,603)20,pwd(1:lll),'mutual2deg'
            endif
            if(igg.eq.1)then
              write(33,602)21,pwd(1:lll),'impact2'
            else
              write(33,603)21,pwd(1:lll),'impact2'
            endif
            if(igg.eq.1)then
              write(33,602)22,pwd(1:lll),'dur2hr'
            else
              write(33,603)22,pwd(1:lll),'dur2hr'
            endif
c
c
            if(igg.eq.1)then
              write(33,602)23,pwd(1:lll),'period3day'
            else
              write(33,603)23,pwd(1:lll),'period3day'
            endif
            if(igg.eq.1)then
              write(33,602)24,pwd(1:lll),'m3earth'
            else
              write(33,603)24,pwd(1:lll),'m3earth'
            endif
            if(igg.eq.1)then
              write(33,602)25,pwd(1:lll),'r3earth'
            else
              write(33,603)25,pwd(1:lll),'r3earth'
            endif
            if(igg.eq.1)then
              write(33,602)26,pwd(1:lll),'den3'
            else
              write(33,603)26,pwd(1:lll),'den3'
            endif
            if(igg.eq.1)then
              write(33,602)27,pwd(1:lll),'a3AU'
            else
              write(33,603)27,pwd(1:lll),'a3AU'
            endif
            if(igg.eq.1)then
              write(33,602)28,pwd(1:lll),'ecc3'
            else
              write(33,603)28,pwd(1:lll),'ecc3'
            endif
            if(igg.eq.1)then
              write(33,602)29,pwd(1:lll),'arg3deg'
            else
              write(33,603)29,pwd(1:lll),'arg3deg'
            endif
            if(igg.eq.1)then
              write(33,602)30,pwd(1:lll),'incl3deg'
            else
              write(33,603)30,pwd(1:lll),'incl3deg'
            endif
            if(igg.eq.1)then
              write(33,602)31,pwd(1:lll),'Omega3deg'
            else
              write(33,603)31,pwd(1:lll),'Omega3deg'
            endif
            if(igg.eq.1)then
              write(33,602)32,pwd(1:lll),'Tconj3'
            else
              write(33,603)32,pwd(1:lll),'Tconj3'
            endif
            if(igg.eq.1)then
              write(33,602)33,pwd(1:lll),'Tperi3'
            else
              write(33,603)33,pwd(1:lll),'Tperi3'
            endif
            if(igg.eq.1)then
              write(33,602)34,pwd(1:lll),'mutual3deg'
            else
              write(33,603)34,pwd(1:lll),'mutual3deg'
            endif
            if(igg.eq.1)then
              write(33,602)35,pwd(1:lll),'impact3'
            else
              write(33,603)35,pwd(1:lll),'impact3'
            endif
            if(igg.eq.1)then
              write(33,602)36,pwd(1:lll),'dur3hr'
            else
              write(33,603)36,pwd(1:lll),'dur3hr'
            endif
c 
            if(igg.eq.1)then
              write(33,602)37,pwd(1:lll),'period4day'
            else
              write(33,603)37,pwd(1:lll),'period4day'
            endif
            if(igg.eq.1)then
              write(33,602)38,pwd(1:lll),'m4earth'
            else
              write(33,603)38,pwd(1:lll),'m4earth'
            endif
            if(igg.eq.1)then
              write(33,602)39,pwd(1:lll),'r4earth'
            else
              write(33,603)39,pwd(1:lll),'r4earth'
            endif
            if(igg.eq.1)then
              write(33,602)40,pwd(1:lll),'den4'
            else
              write(33,603)40,pwd(1:lll),'den4'
            endif
            if(igg.eq.1)then
              write(33,602)41,pwd(1:lll),'a4AU'
            else
              write(33,603)41,pwd(1:lll),'a4AU'
            endif
            if(igg.eq.1)then
              write(33,602)42,pwd(1:lll),'ecc4'
            else
              write(33,603)42,pwd(1:lll),'ecc4'
            endif
            if(igg.eq.1)then
              write(33,602)43,pwd(1:lll),'arg4deg'
            else
              write(33,603)43,pwd(1:lll),'arg4deg'
            endif
            if(igg.eq.1)then
              write(33,602)44,pwd(1:lll),'incl4deg'
            else
              write(33,603)44,pwd(1:lll),'incl4deg'
            endif
            if(igg.eq.1)then
              write(33,602)45,pwd(1:lll),'Omega4deg'
            else
              write(33,603)45,pwd(1:lll),'Omega4deg'
            endif
            if(igg.eq.1)then
              write(33,602)46,pwd(1:lll),'Tconj4'
            else
              write(33,603)46,pwd(1:lll),'Tconj4'
            endif
            if(igg.eq.1)then
              write(33,602)47,pwd(1:lll),'Tperi4'
            else
              write(33,603)47,pwd(1:lll),'Tperi4'
            endif
            if(igg.eq.1)then
              write(33,602)48,pwd(1:lll),'mutual4deg'
            else
              write(33,603)48,pwd(1:lll),'mutual4deg'
            endif
            if(igg.eq.1)then
              write(33,602)49,pwd(1:lll),'impact4'
            else
              write(33,603)49,pwd(1:lll),'impact4'
            endif
            if(igg.eq.1)then
              write(33,602)50,pwd(1:lll),'dur4hr'
            else
              write(33,603)50,pwd(1:lll),'dur4hr'
            endif
c 
            if(igg.eq.1)then
              write(33,602)51,pwd(1:lll),'period5day'
            else
              write(33,603)51,pwd(1:lll),'period5day'
            endif
            if(igg.eq.1)then
              write(33,602)52,pwd(1:lll),'m5earth'
            else
              write(33,603)52,pwd(1:lll),'m5earth'
            endif
            if(igg.eq.1)then
              write(33,602)53,pwd(1:lll),'r5earth'
            else
              write(33,603)53,pwd(1:lll),'r5earth'
            endif
            if(igg.eq.1)then
              write(33,602)54,pwd(1:lll),'den5'
            else
              write(33,603)54,pwd(1:lll),'den5'
            endif
            if(igg.eq.1)then
              write(33,602)55,pwd(1:lll),'a5AU'
            else
              write(33,603)55,pwd(1:lll),'a5AU'
            endif
            if(igg.eq.1)then
              write(33,602)56,pwd(1:lll),'ecc5'
            else
              write(33,603)56,pwd(1:lll),'ecc5'
            endif
            if(igg.eq.1)then
              write(33,602)57,pwd(1:lll),'arg5deg'
            else
              write(33,603)57,pwd(1:lll),'arg5deg'
            endif
            if(igg.eq.1)then
              write(33,602)58,pwd(1:lll),'incl5deg'
            else
              write(33,603)58,pwd(1:lll),'incl5deg'
            endif
            if(igg.eq.1)then
              write(33,602)59,pwd(1:lll),'Omega5deg'
            else
              write(33,603)59,pwd(1:lll),'Omega5deg'
            endif
            if(igg.eq.1)then
              write(33,602)60,pwd(1:lll),'Tconj5'
            else
              write(33,603)60,pwd(1:lll),'Tconj5'
            endif
            if(igg.eq.1)then
              write(33,602)61,pwd(1:lll),'Tperi5'
            else
              write(33,603)61,pwd(1:lll),'Tperi5'
            endif
            if(igg.eq.1)then
              write(33,602)62,pwd(1:lll),'mutual5deg'
            else
              write(33,603)62,pwd(1:lll),'mutual5deg'
            endif
            if(igg.eq.1)then
              write(33,602)63,pwd(1:lll),'impact5'
            else
              write(33,603)63,pwd(1:lll),'impact5'
            endif
            if(igg.eq.1)then
              write(33,602)64,pwd(1:lll),'dur5hr'
            else
              write(33,603)64,pwd(1:lll),'dur5hr'
            endif
c
          endif

 500      format('awk ''{print $',i2,', $2}'' ./',a,
     @            '/generation.* > ttemp.',a)
 501      format('awk ''{print $',i2,', $2}'' ./',a,
     @            '/generation.* >> ttemp.',a)
 502      format('awk ''{print $',i2,', $2}'' ./',a,
     @            '/ELCparm.* > ttemp.',a)
 503      format('awk ''{print $',i2,', $2}'' ./',a,
     @            '/ELCparm.* >> ttemp.',a)
 602      format('awk ''{print $',i2,', $2}'' ./',a,
     @            '/ELCbody3parm.* > ttemp.',a)
 603      format('awk ''{print $',i2,', $2}'' ./',a,
     @            '/ELCbody3parm.* >> ttemp.',a)
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&
c
c   $%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%$%
c
          subroutine chistring(instring,chisq,outstring,ilength)
c
          implicit double precision (a-h,o-z)
c 
          character*(*) instring,outstring
c
          if(chisq.lt.1.0d0)then
            write(outstring,100)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 100      format(a,'=',f7.5)
c
          if(chisq.lt.10.0d0)then
            write(outstring,100)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
c
          if(chisq.lt.100.0d0)then
            write(outstring,101)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 101      format(a,'=',f8.5)
c
          if(chisq.lt.1000.0d0)then
            write(outstring,102)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 102      format(a,'=',f8.4)
c
          if(chisq.lt.10000.0d0)then
            write(outstring,103)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 103      format(a,'=',f8.3)
c
          if(chisq.lt.100000.0d0)then
            write(outstring,104)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 104      format(a,'=',f9.3)
c
          if(chisq.lt.1000000.0d0)then
            write(outstring,105)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 105      format(a,'=',f10.3)
c
          if(chisq.lt.10000000.0d0)then
            write(outstring,106)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 106      format(a,'=',f10.2)
c
          if(chisq.lt.100000000.0d0)then
            write(outstring,107)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 107      format(a,'=',f11.2)
c
          if(chisq.lt.1000000000.0d0)then
            write(outstring,108)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 108      format(a,'=',f12.2)
c
          if(chisq.lt.10000000000.0d0)then
            write(outstring,109)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 109      format(a,'=',f13.2)
c
          if(chisq.lt.100000000000.0d0)then
            write(outstring,110)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
 110      format(a,'=',f14.2)
c
          tempchi=999999999999.99d0
          write(outstring,111)instring,tempchi
          ilength=lnblnk(outstring)
 111      format(a,'=',f15.2)
c          
          return
          end
c
c    $@$@$%%&%&%@%%%%&&&&&&&&&%@@@@@@@@@$%$@%$@$@$%@$%@$%
c
          subroutine pstring(instring,iplace,chisq,outstring,ilength)
c
          implicit double precision (a-h,o-z)
c 
          character*(*) instring,outstring
c
c   iplace=1
c
          if(iplace.eq.1)then
            if(chisq.lt.-10000000.0d0)then
              tempchi=-99999999.9d0
              write(outstring,100)instring,tempchi
              ilength=lnblnk(outstring)
              return
            endif
 100        format(a,'=',f11.1)
c
            if(chisq.lt.-1000000.0d0)then
              write(outstring,101)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 101        format(a,'=',f10.1)
c
            if(chisq.lt.-100000.0d0)then
              write(outstring,102)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 102        format(a,'=',f9.1)
            if(chisq.lt.-9999.9d0)then
              write(outstring,103)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 103        format(a,'=',f8.1)
c
            if(chisq.lt.-999.9d0)then
              write(outstring,104)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 104        format(a,'=',f7.1)
c
            if(chisq.lt.-99.9d0)then
              write(outstring,105)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 105        format(a,'=',f6.1)
c
            if(chisq.lt.-9.9d0)then
              write(outstring,106)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 106        format(a,'=',f5.1)
c
            if(chisq.lt.-0.9d0)then
              write(outstring,107)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 107        format(a,'=',f4.1)
c
            if(chisq.lt.0.0d0)then
              write(outstring,108)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 108        format(a,'=',f4.1)
c
            if(chisq.lt.10.0d0)then
              write(outstring,109)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 109        format(a,'=',f3.1)
c
            if(chisq.lt.100.0d0)then
              write(outstring,110)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 110        format(a,'=',f4.1)
c
            if(chisq.lt.1000.0d0)then
              write(outstring,111)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 111        format(a,'=',f5.1)
c
            if(chisq.lt.10000.0d0)then
              write(outstring,112)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 112        format(a,'=',f6.1)
c
            if(chisq.lt.100000.0d0)then
              write(outstring,113)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 113        format(a,'=',f7.1)
c
            if(chisq.lt.1000000.0d0)then
              write(outstring,114)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 114        format(a,'=',f8.1)
c
            if(chisq.lt.10000000.0d0)then
              write(outstring,115)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 115        format(a,'=',f9.1)
c
            if(chisq.lt.100000000.0d0)then
              write(outstring,116)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 116        format(a,'=',f10.1)
c
            if(chisq.lt.1000000000.0d0)then
              write(outstring,117)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 117        format(a,'=',f11.1)
c
            tempchi=9999999999.9d0
            write(outstring,118)instring,tempchi
            ilength=lnblnk(outstring)
 118        format(a,'=',f12.1)
            return
c
          endif     ! iplace=1
c
c   iplace=2
c
          if(iplace.eq.2)then
            if(chisq.lt.-10000000.0d0)then
              tempchi=-99999999.9d0
              write(outstring,200)instring,tempchi
              ilength=lnblnk(outstring)
              return
            endif
 200        format(a,'=',f12.2)
c
            if(chisq.lt.-1000000.0d0)then
              write(outstring,201)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 201        format(a,'=',f11.2)
c
            if(chisq.lt.-100000.0d0)then
              write(outstring,202)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 202        format(a,'=',f10.2)
c
            if(chisq.lt.-9999.9d0)then
              write(outstring,203)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 203        format(a,'=',f9.2)
c
            if(chisq.lt.-999.9d0)then
              write(outstring,204)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 204        format(a,'=',f8.2)
c
            if(chisq.lt.-99.9d0)then
              write(outstring,205)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 205        format(a,'=',f7.2)
c
            if(chisq.lt.-9.9d0)then
              write(outstring,206)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 206        format(a,'=',f6.2)
c
            if(chisq.lt.-0.9d0)then
              write(outstring,207)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 207        format(a,'=',f5.2)
c
            if(chisq.lt.0.0d0)then
              write(outstring,208)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 208        format(a,'=',f5.2)
c
            if(chisq.lt.10.0d0)then
              write(outstring,209)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 209        format(a,'=',f4.2)
c
            if(chisq.lt.100.0d0)then
              write(outstring,210)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 210        format(a,'=',f5.2)
c
            if(chisq.lt.1000.0d0)then
              write(outstring,211)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 211        format(a,'=',f6.2)
c
            if(chisq.lt.10000.0d0)then
              write(outstring,212)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 212        format(a,'=',f7.2)
c
            if(chisq.lt.100000.0d0)then
              write(outstring,213)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 213        format(a,'=',f8.2)
c
            if(chisq.lt.1000000.0d0)then
              write(outstring,214)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 214        format(a,'=',f9.2)
c
            if(chisq.lt.10000000.0d0)then
              write(outstring,215)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 215        format(a,'=',f10.2)
c
            if(chisq.lt.100000000.0d0)then
              write(outstring,216)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 216        format(a,'=',f11.2)
c
            if(chisq.lt.1000000000.0d0)then
              write(outstring,217)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 217        format(a,'=',f12.2)
c
            tempchi=9999999999.99d0
            write(outstring,218)instring,tempchi
            ilength=lnblnk(outstring)
 218        format(a,'=',f13.2)
            return
c
          endif     ! iplace=2
c
c   iplace=3
c
          if(iplace.eq.3)then
            if(chisq.lt.-10000000.0d0)then
              tempchi=-99999999.9d0
              write(outstring,300)instring,tempchi
              ilength=lnblnk(outstring)
              return
            endif
 300        format(a,'=',f13.3)
c
            if(chisq.lt.-1000000.0d0)then
              write(outstring,301)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 301        format(a,'=',f12.3)
c
            if(chisq.lt.-100000.0d0)then
              write(outstring,302)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 302        format(a,'=',f11.3)
c
            if(chisq.lt.-9999.9d0)then
              write(outstring,303)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 303        format(a,'=',f10.3)
c
            if(chisq.lt.-999.9d0)then
              write(outstring,304)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 304        format(a,'=',f9.3)
c
            if(chisq.lt.-99.9d0)then
              write(outstring,305)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 305        format(a,'=',f8.3)
c
            if(chisq.lt.-9.9d0)then
              write(outstring,306)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 306        format(a,'=',f7.3)
c
            if(chisq.lt.-0.9d0)then
              write(outstring,307)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 307        format(a,'=',f6.3)
c
            if(chisq.lt.0.0d0)then
              write(outstring,308)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 308        format(a,'=',f6.3)
c
            if(chisq.lt.10.0d0)then
              write(outstring,309)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 309        format(a,'=',f5.3)
c
            if(chisq.lt.100.0d0)then
              write(outstring,310)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 310        format(a,'=',f6.3)
c
            if(chisq.lt.1000.0d0)then
              write(outstring,311)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 311        format(a,'=',f7.3)
c
            if(chisq.lt.10000.0d0)then
              write(outstring,312)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 312        format(a,'=',f8.3)
c
            if(chisq.lt.100000.0d0)then
              write(outstring,313)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 313        format(a,'=',f9.3)
c
            if(chisq.lt.1000000.0d0)then
              write(outstring,314)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 314        format(a,'=',f10.3)
c
            if(chisq.lt.10000000.0d0)then
              write(outstring,315)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 315        format(a,'=',f11.3)
c
            if(chisq.lt.100000000.0d0)then
              write(outstring,316)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 316        format(a,'=',f12.3)
c
            if(chisq.lt.1000000000.0d0)then
              write(outstring,317)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 317        format(a,'=',f13.3)
c
            tempchi=9999999999.999d0
            write(outstring,318)instring,tempchi
            ilength=lnblnk(outstring)
 318        format(a,'=',f14.3)
            return
c
          endif     ! iplace=3
c
c   iplace=4
c
          if(iplace.eq.4)then
            if(chisq.lt.-10000000.0d0)then
              tempchi=-99999999.9d0
              write(outstring,400)instring,tempchi
              ilength=lnblnk(outstring)
              return
            endif
 400        format(a,'=',f14.4)
c
            if(chisq.lt.-1000000.0d0)then
              write(outstring,401)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 401        format(a,'=',f13.4)
c
            if(chisq.lt.-100000.0d0)then
              write(outstring,402)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 402        format(a,'=',f12.4)
c
            if(chisq.lt.-9999.9d0)then
              write(outstring,403)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 403        format(a,'=',f11.4)
c
            if(chisq.lt.-999.9d0)then
              write(outstring,404)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 404        format(a,'=',f10.4)
c
            if(chisq.lt.-99.9d0)then
              write(outstring,405)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 405        format(a,'=',f9.4)
c
            if(chisq.lt.-9.9d0)then
              write(outstring,406)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 406        format(a,'=',f8.4)
c
            if(chisq.lt.-0.9d0)then
              write(outstring,407)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 407        format(a,'=',f7.4)
c
            if(chisq.lt.0.0d0)then
              write(outstring,408)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 408        format(a,'=',f7.4)
c
            if(chisq.lt.10.0d0)then
              write(outstring,409)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 409        format(a,'=',f6.4)
c
            if(chisq.lt.100.0d0)then
              write(outstring,410)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 410        format(a,'=',f7.4)
c
            if(chisq.lt.1000.0d0)then
              write(outstring,411)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 411        format(a,'=',f8.4)
c
            if(chisq.lt.10000.0d0)then
              write(outstring,412)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 412        format(a,'=',f9.4)
c
            if(chisq.lt.100000.0d0)then
              write(outstring,413)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 413        format(a,'=',f10.4)
c
            if(chisq.lt.1000000.0d0)then
              write(outstring,414)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 414        format(a,'=',f11.4)
c
            if(chisq.lt.10000000.0d0)then
              write(outstring,415)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 415        format(a,'=',f12.4)
c
            if(chisq.lt.100000000.0d0)then
              write(outstring,416)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 416        format(a,'=',f13.4)
c
            if(chisq.lt.1000000000.0d0)then
              write(outstring,417)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 417        format(a,'=',f14.4)
c
            tempchi=9999999999.9999d0
            write(outstring,418)instring,tempchi
            ilength=lnblnk(outstring)
 418        format(a,'=',f15.4)
            return
c
          endif     ! iplace=4
c
c   iplace=5
c
          if(iplace.eq.5)then
            if(chisq.lt.-10000000.0d0)then
              tempchi=-99999999.9d0
              write(outstring,500)instring,tempchi
              ilength=lnblnk(outstring)
              return
            endif
 500        format(a,'=',f15.5)
c
            if(chisq.lt.-1000000.0d0)then
              write(outstring,501)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 501        format(a,'=',f14.5)
c
            if(chisq.lt.-100000.0d0)then
              write(outstring,502)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 502        format(a,'=',f13.5)
c
            if(chisq.lt.-9999.9d0)then
              write(outstring,503)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 503        format(a,'=',f12.5)
c
            if(chisq.lt.-999.9d0)then
              write(outstring,504)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 504        format(a,'=',f11.5)
c
            if(chisq.lt.-99.9d0)then
              write(outstring,505)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 505        format(a,'=',f10.5)
c
            if(chisq.lt.-9.9d0)then
              write(outstring,506)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 506        format(a,'=',f9.5)
c
            if(chisq.lt.-0.9d0)then
              write(outstring,507)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 507        format(a,'=',f8.5)
c
            if(chisq.lt.0.0d0)then
              write(outstring,508)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 508        format(a,'=',f8.5)
c
            if(chisq.lt.10.0d0)then
              write(outstring,509)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 509        format(a,'=',f7.5)
c
            if(chisq.lt.100.0d0)then
              write(outstring,510)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 510        format(a,'=',f8.5)
c
            if(chisq.lt.1000.0d0)then
              write(outstring,511)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 511        format(a,'=',f9.5)
c
            if(chisq.lt.10000.0d0)then
              write(outstring,512)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 512        format(a,'=',f10.5)
c
            if(chisq.lt.100000.0d0)then
              write(outstring,513)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 513        format(a,'=',f11.5)
c
            if(chisq.lt.1000000.0d0)then
              write(outstring,514)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 514        format(a,'=',f12.5)
c
            if(chisq.lt.10000000.0d0)then
              write(outstring,515)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 515        format(a,'=',f13.5)
c
            if(chisq.lt.100000000.0d0)then
              write(outstring,516)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 516        format(a,'=',f14.5)
c
            if(chisq.lt.1000000000.0d0)then
              write(outstring,517)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 517        format(a,'=',f15.5)
c
            tempchi=9999999999.99999d0
            write(outstring,518)instring,tempchi
            ilength=lnblnk(outstring)
 518        format(a,'=',f16.5)
            return
c
          endif     ! iplace=5
c
c   iplace=6
c
          if(iplace.eq.6)then
            if(chisq.lt.-10000000.0d0)then
              tempchi=-99999999.9d0
              write(outstring,600)instring,tempchi
              ilength=lnblnk(outstring)
              return
            endif
 600        format(a,'=',f16.6)
c
            if(chisq.lt.-1000000.0d0)then
              write(outstring,601)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 601        format(a,'=',f15.6)
c
            if(chisq.lt.-100000.0d0)then
              write(outstring,602)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 602        format(a,'=',f14.6)
c
            if(chisq.lt.-9999.9d0)then
              write(outstring,603)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 603        format(a,'=',f13.6)
c
            if(chisq.lt.-999.9d0)then
              write(outstring,604)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 604        format(a,'=',f12.6)
c
            if(chisq.lt.-99.9d0)then
              write(outstring,605)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 605        format(a,'=',f11.6)
c
            if(chisq.lt.-9.9d0)then
              write(outstring,606)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 606        format(a,'=',f10.6)
c
            if(chisq.lt.-0.9d0)then
              write(outstring,607)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 607        format(a,'=',f9.6)
c
            if(chisq.lt.0.0d0)then
              write(outstring,608)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 608        format(a,'=',f9.6)
c
            if(chisq.lt.10.0d0)then
              write(outstring,609)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 609        format(a,'=',f8.6)
c
            if(chisq.lt.100.0d0)then
              write(outstring,610)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 610        format(a,'=',f9.6)
c
            if(chisq.lt.1000.0d0)then
              write(outstring,611)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 611        format(a,'=',f10.6)
c
            if(chisq.lt.10000.0d0)then
              write(outstring,612)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 612        format(a,'=',f11.6)
c
            if(chisq.lt.100000.0d0)then
              write(outstring,613)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 613        format(a,'=',f12.6)
c
            if(chisq.lt.1000000.0d0)then
              write(outstring,614)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 614        format(a,'=',f13.6)
c
            if(chisq.lt.10000000.0d0)then
              write(outstring,615)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 615        format(a,'=',f14.6)
c
            if(chisq.lt.100000000.0d0)then
              write(outstring,616)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 616        format(a,'=',f15.6)
c
            if(chisq.lt.1000000000.0d0)then
              write(outstring,617)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 617        format(a,'=',f16.6)
c
            tempchi=9999999999.999999d0
            write(outstring,618)instring,tempchi
            ilength=lnblnk(outstring)
 618        format(a,'=',f17.6)
            return
c
          endif     ! iplace=6
c
c   iplace=7
c
          if(iplace.eq.7)then
            if(chisq.lt.-10000000.0d0)then
              tempchi=-99999999.9d0
              write(outstring,700)instring,tempchi
              ilength=lnblnk(outstring)
              return
            endif
 700        format(a,'=',f17.7)
c
            if(chisq.lt.-1000000.0d0)then
              write(outstring,701)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 701        format(a,'=',f16.7)
c
            if(chisq.lt.-100000.0d0)then
              write(outstring,702)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 702        format(a,'=',f15.7)
c
            if(chisq.lt.-9999.9d0)then
              write(outstring,703)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 703        format(a,'=',f14.7)
c
            if(chisq.lt.-999.9d0)then
              write(outstring,704)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 704        format(a,'=',f13.7)
c
            if(chisq.lt.-99.9d0)then
              write(outstring,705)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 705        format(a,'=',f12.7)
c
            if(chisq.lt.-9.9d0)then
              write(outstring,706)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 706        format(a,'=',f11.7)
c
            if(chisq.lt.-0.9d0)then
              write(outstring,707)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 707        format(a,'=',f10.7)
c
            if(chisq.lt.0.0d0)then
              write(outstring,708)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 708        format(a,'=',f10.7)
c
            if(chisq.lt.10.0d0)then
              write(outstring,709)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 709        format(a,'=',f9.7)
c
            if(chisq.lt.100.0d0)then
              write(outstring,710)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 710        format(a,'=',f10.7)
c
            if(chisq.lt.1000.0d0)then
              write(outstring,711)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 711        format(a,'=',f11.7)
c
            if(chisq.lt.10000.0d0)then
              write(outstring,712)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 712        format(a,'=',f12.7)
c
            if(chisq.lt.100000.0d0)then
              write(outstring,713)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 713        format(a,'=',f13.7)
c
            if(chisq.lt.1000000.0d0)then
              write(outstring,714)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 714        format(a,'=',f14.7)
c
            if(chisq.lt.10000000.0d0)then
              write(outstring,715)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 715        format(a,'=',f15.7)
c
            if(chisq.lt.100000000.0d0)then
              write(outstring,716)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 716        format(a,'=',f16.7)
c
            if(chisq.lt.1000000000.0d0)then
              write(outstring,717)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 717        format(a,'=',f17.7)
c
            tempchi=9999999999.9999990d0
            write(outstring,718)instring,tempchi
            ilength=lnblnk(outstring)
 718        format(a,'=',f18.7)
            return
c
          endif     ! iplace=7
c
c   iplace=8
c
          if(iplace.eq.8)then
            if(chisq.lt.-10000000.0d0)then
              tempchi=-99999999.9d0
              write(outstring,800)instring,tempchi
              ilength=lnblnk(outstring)
              return
            endif
 800        format(a,'=',f18.8)
c
            if(chisq.lt.-1000000.0d0)then
              write(outstring,801)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 801        format(a,'=',f17.8)
c
            if(chisq.lt.-100000.0d0)then
              write(outstring,802)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 802        format(a,'=',f16.8)
c
            if(chisq.lt.-9999.9d0)then
              write(outstring,803)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 803        format(a,'=',f15.8)
c
            if(chisq.lt.-999.9d0)then
              write(outstring,804)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 804        format(a,'=',f14.8)
c
            if(chisq.lt.-99.9d0)then
              write(outstring,805)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 805        format(a,'=',f13.8)
c
            if(chisq.lt.-9.9d0)then
              write(outstring,806)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 806        format(a,'=',f12.8)
c
            if(chisq.lt.-0.9d0)then
              write(outstring,807)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 807        format(a,'=',f11.8)
c
            if(chisq.lt.0.0d0)then
              write(outstring,808)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 808        format(a,'=',f11.8)
c
            if(chisq.lt.10.0d0)then
              write(outstring,809)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 809        format(a,'=',f10.8)
c
            if(chisq.lt.100.0d0)then
              write(outstring,810)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 810        format(a,'=',f11.8)
c
            if(chisq.lt.1000.0d0)then
              write(outstring,811)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 811        format(a,'=',f12.8)
c
            if(chisq.lt.10000.0d0)then
              write(outstring,812)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 812        format(a,'=',f13.8)
c
            if(chisq.lt.100000.0d0)then
              write(outstring,813)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 813        format(a,'=',f14.8)
c
            if(chisq.lt.1000000.0d0)then
              write(outstring,814)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 814        format(a,'=',f15.8)
c
            if(chisq.lt.10000000.0d0)then
              write(outstring,815)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 815        format(a,'=',f16.8)
c
            if(chisq.lt.100000000.0d0)then
              write(outstring,816)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 816        format(a,'=',f17.8)
c
            if(chisq.lt.1000000000.0d0)then
              write(outstring,817)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 817        format(a,'=',f18.8)
c
            tempchi=9999999999.999999d0
            write(outstring,818)instring,tempchi
            ilength=lnblnk(outstring)
 818        format(a,'=',f19.8)
            return
c
          endif     ! iplace=8
c
c   iplace=9
c
          if(iplace.eq.9)then
            if(chisq.lt.-10000000.0d0)then
              tempchi=-99999999.9d0
              write(outstring,900)instring,tempchi
              ilength=lnblnk(outstring)
              return
            endif
 900        format(a,'=',f19.9)
c
            if(chisq.lt.-1000000.0d0)then
              write(outstring,901)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 901        format(a,'=',f18.9)
c
            if(chisq.lt.-100000.0d0)then
              write(outstring,902)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 902        format(a,'=',f17.9)
c
            if(chisq.lt.-9999.9d0)then
              write(outstring,903)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 903        format(a,'=',f16.9)
c
            if(chisq.lt.-999.9d0)then
              write(outstring,904)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 904        format(a,'=',f15.9)
c
            if(chisq.lt.-99.9d0)then
              write(outstring,905)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 905        format(a,'=',f14.9)
c
            if(chisq.lt.-9.9d0)then
              write(outstring,906)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 906        format(a,'=',f13.9)
c
            if(chisq.lt.-0.9d0)then
              write(outstring,907)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 907        format(a,'=',f12.9)
c
            if(chisq.lt.0.0d0)then
              write(outstring,908)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 908        format(a,'=',f12.9)
c
            if(chisq.lt.10.0d0)then
              write(outstring,909)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 909        format(a,'=',f11.9)
c
            if(chisq.lt.100.0d0)then
              write(outstring,910)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 910        format(a,'=',f12.9)
c
            if(chisq.lt.1000.0d0)then
              write(outstring,911)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 911        format(a,'=',f13.9)
c
            if(chisq.lt.10000.0d0)then
              write(outstring,912)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 912        format(a,'=',f14.9)
c
            if(chisq.lt.100000.0d0)then
              write(outstring,913)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 913        format(a,'=',f15.9)
c
            if(chisq.lt.1000000.0d0)then
              write(outstring,914)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 914        format(a,'=',f16.9)
c
            if(chisq.lt.10000000.0d0)then
              write(outstring,915)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 915        format(a,'=',f17.9)
c
            if(chisq.lt.100000000.0d0)then
              write(outstring,916)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 916        format(a,'=',f18.9)
c
            if(chisq.lt.1000000000.0d0)then
              write(outstring,917)instring,chisq
              ilength=lnblnk(outstring)
              return
            endif
 917        format(a,'=',f19.9)
c
            tempchi=9999999999.999999d0
            write(outstring,918)instring,tempchi
            ilength=lnblnk(outstring)
 918        format(a,'=',f20.9)
            return
c
          endif     ! iplace=9
c
c
c
          return
          end
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printS(chi1)
c
          implicit double precision(a-h,o-z)
c
          if(chi1.gt.9.999999999999999D20)then
            write(*,100)
            return
          endif
          if(chi1.gt.9.999999999999999D19)then
            write(*,101)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D18)then
            write(*,102)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D17)then
            write(*,103)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D16)then
            write(*,104)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D15)then
            write(*,105)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D14)then
            write(*,106)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D13)then
            write(*,107)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D12)then
            write(*,108)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D11)then
            write(*,109)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D10)then
            write(*,110)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D9)then
            write(*,111)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D8)then
            write(*,112)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D7)then
            write(*,113)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D6)then
            write(*,114)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D5)then
            write(*,115)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D4)then
            write(*,116)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D3)then
            write(*,117)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D2)then
            write(*,118)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D1)then
            write(*,119)chi1
            return
          endif
          if(chi1.gt.9.999999999999999D0)then
            write(*,120)chi1
            return
          endif
          write(*,121)chi1

100       format('S = 999999999999999999999.999999')
101       format('S = ',f28.6)
102       format('S = ',f27.6)
103       format('S = ',f26.6)
104       format('S = ',f25.6)
105       format('S = ',f24.6)
106       format('S = ',f23.6)
107       format('S = ',f22.6)
108       format('S = ',f21.6)
109       format('S = ',f20.6)
110       format('S = ',f19.6)
111       format('S = ',f18.6)
112       format('S = ',f17.6)
113       format('S = ',f16.6)
114       format('S = ',f15.6)
115       format('S = ',f14.6)
116       format('S = ',f13.6)
117       format('S = ',f12.6)
118       format('S = ',f11.6)
119       format('S = ',f10.6)
120       format('S = ',f9.6)
121       format('S = ',f8.6)

c
          return
          end
c
c
          subroutine writejump(mmm,nnn,jumpindex)
c          
          if(jumpindex.ge.10)then
            if((mmm.eq.1).and.(nnn.eq.1))then
              write(*,94)mmm,nnn,jumpindex
            endif
            if((mmm.eq.1).and.((nnn.gt.1).and.(nnn.lt.10)))then
              write(*,394)mmm,nnn,jumpindex
            endif
            if((mmm.eq.1).and.((nnn.ge.10).and.(nnn.lt.100)))then
              write(*,494)mmm,nnn,jumpindex
            endif
            if((mmm.eq.1).and.((nnn.ge.100).and.(nnn.lt.1000)))then
              write(*,594)mmm,nnn,jumpindex
            endif
            if((mmm.eq.1).and.((nnn.ge.1000).and.(nnn.lt.10000)))then
              write(*,694)mmm,nnn,jumpindex
            endif
c
            if(((mmm.gt.1).and.(mmm.lt.10)).and.((nnn.gt.1).and.
     @            (nnn.lt.10)))then
              write(*,1394)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.1).and.(mmm.lt.10)).and.((nnn.ge.10).and.
     @           (nnn.lt.100)))then
              write(*,1494)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.1).and.(mmm.lt.10)).and.((nnn.ge.100).and.
     @           (nnn.lt.1000)))then
              write(*,1594)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.1).and.(mmm.lt.10)).and.((nnn.ge.1000).and.
     @           (nnn.lt.10000)))then
              write(*,1694)mmm,nnn,jumpindex
            endif
c
            if(((mmm.gt.10).and.(mmm.lt.100)).and.((nnn.gt.1).and.
     @            (nnn.lt.10)))then
              write(*,2394)mmm,nnn,jumpindex
             endif
            if(((mmm.gt.10).and.(mmm.lt.100)).and.((nnn.ge.10).and.
     @            (nnn.lt.100)))then
              write(*,2494)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.10).and.(mmm.lt.100)).and.((nnn.ge.100).and.
     @             (nnn.lt.1000)))then
              write(*,2594)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.10).and.(mmm.lt.100)).and.((nnn.ge.1000).and.
     @              (nnn.lt.10000)))then
              write(*,2694)mmm,nnn,jumpindex
            endif
c
            if(((mmm.gt.100).and.(mmm.lt.1000)).and.((nnn.gt.1).and.
     @           (nnn.lt.10)))then
              write(*,3394)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.100).and.(mmm.lt.1000)).and.((nnn.ge.10).and.
     @           (nnn.lt.100)))then
              write(*,3494)mmm,nnn,jumpindex
            endif
           if(((mmm.gt.100).and.(mmm.lt.1000)).and.((nnn.ge.100).and.
     @             (nnn.lt.1000)))then
              write(*,3594)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.100).and.(mmm.lt.1000)).and.
     @            ((nnn.ge.1000).and.(nnn.lt.10000)))then
              write(*,3694)mmm,nnn,jumpindex
            endif
           if(((mmm.gt.1000).and.(mmm.lt.10000)).and.
     @            ((nnn.ge.1000).and.(nnn.lt.10000)))then
              write(*,4694)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.1000).and.(mmm.lt.10000)).and.
     @          ((nnn.ge.10000).and.(nnn.lt.100000)))then
              write(*,4695)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.10000).and.(mmm.lt.100000)).and.
     @           ((nnn.ge.10000).and.(nnn.lt.100000)))then
              write(*,4696)mmm,nnn,jumpindex
            endif
          endif
c
94        format(3x,i1,' jump taken out of ',i1,
     @      ' try for parameter #',i2)
394       format(3x,i1,' jump taken out of ',i1,
     @      ' tries for parameter #',i2)
494       format(3x,i1,' jump taken out of ',i2,
     @      ' tries for parameter #',i2)
594       format(3x,i1,' jump taken out of ',i3,
     @      ' tries for parameter #',i2)
694       format(3x,i1,' jump taken out of ',i4,
     @      ' tries for parameter #',i2)
1394      format(3x,i1,' jumps taken out of ',i1,
     @      ' tries for parameter #',i2)
1494      format(3x,i1,' jumps taken out of ',i2,
     @      ' tries for parameter #',i2)
1594      format(3x,i1,' jumps taken out of ',i3,
     @      ' tries for parameter #',i2)
1694      format(3x,i1,' jumps taken out of ',i4,
     @      ' tries for parameter #',i2)
2394      format(3x,i2,' jumps taken out of ',i1,
     @      ' tries for parameter #',i2)
2494      format(3x,i2,' jumps taken out of ',i2,
     @      ' tries for parameter #',i2)
2594      format(3x,i2,' jumps taken out of ',i3,
     @      ' tries for parameter #',i2)
2694      format(3x,i2,' jumps taken out of ',i4,
     @      ' tries for parameter #',i2)
3394      format(3x,i3,' jumps taken out of ',i1,
     @      ' tries for parameter #',i2)
3494      format(3x,i3,' jumps taken out of ',i2,
     @      ' tries for parameter #',i2)
3594      format(3x,i3,' jumps taken out of ',i3,
     @      ' tries for parameter #',i2)
3694      format(3x,i3,' jumps taken out of ',i4,
     @      ' tries for parameter #',i2)
4694      format(3x,i4,' jumps taken out of ',i4,
     @      ' tries for parameter #',i2)
4695      format(3x,i4,' jumps taken out of ',i5,
     @      ' tries for parameter #',i2)
4696      format(3x,i5,' jumps taken out of ',i5,
     @          ' tries for parameter #',i2)
c
c
c
          if(jumpindex.lt.10)then
            if((mmm.eq.1).and.(nnn.eq.1))then
              write(*,34)mmm,nnn,jumpindex
            endif
            if((mmm.eq.1).and.((nnn.gt.1).and.(nnn.lt.10)))then
              write(*,334)mmm,nnn,jumpindex
            endif
            if((mmm.eq.1).and.((nnn.ge.10).and.(nnn.lt.100)))then
              write(*,434)mmm,nnn,jumpindex
            endif
            if((mmm.eq.1).and.((nnn.ge.100).and.(nnn.lt.1000)))then
              write(*,534)mmm,nnn,jumpindex
            endif
            if((mmm.eq.1).and.((nnn.ge.1000).and.(nnn.lt.10000)))then
              write(*,634)mmm,nnn,jumpindex
            endif
c
            if(((mmm.gt.1).and.(mmm.lt.10)).and.((nnn.gt.1).and.
     @            (nnn.lt.10)))then
              write(*,1334)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.1).and.(mmm.lt.10)).and.((nnn.ge.10).and.
     @           (nnn.lt.100)))then
              write(*,1434)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.1).and.(mmm.lt.10)).and.((nnn.ge.100).and.
     @           (nnn.lt.1000)))then
              write(*,1534)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.1).and.(mmm.lt.10)).and.((nnn.ge.1000).and.
     @           (nnn.lt.10000)))then
              write(*,1634)mmm,nnn,jumpindex
            endif
c
            if(((mmm.gt.10).and.(mmm.lt.100)).and.((nnn.gt.1).and.
     @            (nnn.lt.10)))then
              write(*,2334)mmm,nnn,jumpindex
             endif
            if(((mmm.gt.10).and.(mmm.lt.100)).and.((nnn.ge.10).and.
     @            (nnn.lt.100)))then
              write(*,2434)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.10).and.(mmm.lt.100)).and.((nnn.ge.100).and.
     @             (nnn.lt.1000)))then
              write(*,2534)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.10).and.(mmm.lt.100)).and.((nnn.ge.1000).and.
     @              (nnn.lt.10000)))then
              write(*,2634)mmm,nnn,jumpindex
            endif
c
            if(((mmm.gt.100).and.(mmm.lt.1000)).and.((nnn.gt.1).and.
     @           (nnn.lt.10)))then
              write(*,3334)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.100).and.(mmm.lt.1000)).and.((nnn.ge.10).and.
     @           (nnn.lt.100)))then
              write(*,3434)mmm,nnn,jumpindex
            endif
           if(((mmm.gt.100).and.(mmm.lt.1000)).and.((nnn.ge.100).and.
     @             (nnn.lt.1000)))then
              write(*,3534)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.100).and.(mmm.lt.1000)).and.
     @            ((nnn.ge.1000).and.(nnn.lt.10000)))then
              write(*,3634)mmm,nnn,jumpindex
            endif
           if(((mmm.gt.1000).and.(mmm.lt.10000)).and.
     @            ((nnn.ge.1000).and.(nnn.lt.10000)))then
              write(*,4634)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.1000).and.(mmm.lt.10000)).and.
     @          ((nnn.ge.10000).and.(nnn.lt.100000)))then
              write(*,4635)mmm,nnn,jumpindex
            endif
            if(((mmm.gt.10000).and.(mmm.lt.100000)).and.
     @           ((nnn.ge.10000).and.(nnn.lt.100000)))then
              write(*,4636)mmm,nnn,jumpindex
            endif
          endif
c
34        format(3x,i1,' jump taken out of ',i1,
     @      ' try for parameter #',i1)
334       format(3x,i1,' jump taken out of ',i1,
     @      ' tries for parameter #',i1)
434       format(3x,i1,' jump taken out of ',i2,
     @      ' tries for parameter #',i1)
534       format(3x,i1,' jump taken out of ',i3,
     @      ' tries for parameter #',i1)
634       format(3x,i1,' jump taken out of ',i4,
     @      ' tries for parameter #',i1)
1334      format(3x,i1,' jumps taken out of ',i1,
     @      ' tries for parameter #',i1)
1434      format(3x,i1,' jumps taken out of ',i2,
     @      ' tries for parameter #',i1)
1534      format(3x,i1,' jumps taken out of ',i3,
     @      ' tries for parameter #',i1)
1634      format(3x,i1,' jumps taken out of ',i4,
     @      ' tries for parameter #',i1)
2334      format(3x,i2,' jumps taken out of ',i1,
     @      ' tries for parameter #',i1)
2434      format(3x,i2,' jumps taken out of ',i2,
     @      ' tries for parameter #',i1)
2534      format(3x,i2,' jumps taken out of ',i3,
     @      ' tries for parameter #',i1)
2634      format(3x,i2,' jumps taken out of ',i4,
     @      ' tries for parameter #',i1)
3334      format(3x,i3,' jumps taken out of ',i1,
     @      ' tries for parameter #',i1)
3434      format(3x,i3,' jumps taken out of ',i2,
     @      ' tries for parameter #',i1)
3534      format(3x,i3,' jumps taken out of ',i3,
     @      ' tries for parameter #',i1)
3634      format(3x,i3,' jumps taken out of ',i4,
     @      ' tries for parameter #',i1)
4634      format(3x,i4,' jumps taken out of ',i4,
     @      ' tries for parameter #',i1)
4635      format(3x,i4,' jumps taken out of ',i5,
     @      ' tries for parameter #',i1)
4636      format(3x,i5,' jumps taken out of ',i5,
     @          ' tries for parameter #',i1)

          return
          end
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printiter(imod,instring1,iter,instring2)
c
          character*(*)instring1,instring2
          character*(80) outstring1,outstring2
c
          if(imod.lt.0)return
          if(imod.ge.1000000)return
c
          l1=lnblnk(instring1)
          l1=l1+1
          if(imod.lt.10)write(outstring1,100)instring1(1:l1),imod
          if((imod.ge.10).and.(imod.lt.100))then
            write(outstring1,101)instring1(1:l1),imod
          endif
          if((imod.ge.100).and.(imod.lt.1000))then
            write(outstring1,102)instring1(1:l1),imod
          endif
          if((imod.ge.1000).and.(imod.lt.10000))then
            write(outstring1,103)instring1(1:l1),imod
          endif
          if((imod.ge.10000).and.(imod.lt.100000))then
            write(outstring1,104)instring1(1:l1),imod
          endif
          if((imod.ge.100000).and.(imod.lt.1000000))then
            write(outstring1,105)instring1(1:l1),imod
          endif
c
          if(iter.lt.0)return
          if(iter.ge.1000000)return
c
          l2=lnblnk(instring2)
          l2=l2+1
          if(iter.lt.10)write(outstring2,100)instring2(1:l2),iter
          if((iter.ge.10).and.(iter.lt.100))then
            write(outstring2,101)instring2(1:l2),iter
          endif
          if((iter.ge.100).and.(iter.lt.1000))then
            write(outstring2,102)instring2(1:l2),iter
          endif
          if((iter.ge.1000).and.(iter.lt.10000))then
            write(outstring2,103)instring2(1:l2),iter
          endif
          if((iter.ge.10000).and.(iter.lt.100000))then
            write(outstring2,104)instring2(1:l2),iter
          endif
          if((iter.ge.100000).and.(iter.lt.1000000))then
            write(outstring2,105)instring2(1:l2),iter
          endif
c
          ll1=lnblnk(outstring1)
          ll2=lnblnk(outstring2)
c
          write(*,200)outstring1(1:ll1),outstring2(1:ll2)
c
 100      format(a,i1)
 101      format(a,i2)
 102      format(a,i3)
 103      format(a,i4)
 104      format(a,i5)
 105      format(a,i6)
 200      format(a,', ',a)
          return
          end
c
c
c @#$&^*&*$*&@&%#&%#*&^%#^$@%^(&^%@$&@@@@#$&^*&*$*&@&%#&%#*&^%#^$@%^(&^%
c
          subroutine printamoebamed(nn,chi1)
c
          implicit double precision(a-h,o-z)
c
          if(nn.lt.10)then
            if(chi1.gt.9.999999999999999D20)then
              write(*,100)nn
              return
            endif
            if(chi1.gt.9.999999999999999D19)then
              write(*,101)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D18)then
              write(*,102)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D17)then
              write(*,103)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D16)then
              write(*,104)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D15)then
              write(*,105)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D14)then
              write(*,106)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D13)then
              write(*,107)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D12)then
              write(*,108)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D11)then
              write(*,109)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D10)then
              write(*,110)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D9)then
              write(*,111)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D8)then
              write(*,112)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D7)then
              write(*,113)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D6)then
              write(*,114)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D5)then
              write(*,115)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D4)then
              write(*,116)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D3)then
              write(*,117)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D2)then
              write(*,118)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D1)then
              write(*,119)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D0)then
              write(*,120)nn,chi1
              return
            endif
            write(*,121)nn,chi1
          endif !nn < 10
c
          if((nn.ge.10).and.(nn.lt.100))then
            if(chi1.gt.9.999999999999999D20)then
              write(*,200)nn
              return
            endif
            if(chi1.gt.9.999999999999999D19)then
              write(*,201)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D18)then
              write(*,202)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D17)then
              write(*,203)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D16)then
              write(*,204)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D15)then
              write(*,205)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D14)then
              write(*,206)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D13)then
              write(*,207)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D12)then
              write(*,208)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D11)then
              write(*,209)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D10)then
              write(*,210)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D9)then
              write(*,211)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D8)then
              write(*,212)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D7)then
              write(*,213)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D6)then
              write(*,214)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D5)then
              write(*,215)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D4)then
              write(*,216)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D3)then
              write(*,217)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D2)then
              write(*,218)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D1)then
              write(*,219)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D0)then
              write(*,220)nn,chi1
              return
            endif
            write(*,221)nn,chi1
          endif !nn > 10 and < 100
c
          if((nn.ge.100).and.(nn.lt.1000))then
            if(chi1.gt.9.999999999999999D20)then
              write(*,300)nn
              return
            endif
            if(chi1.gt.9.999999999999999D19)then
              write(*,301)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D18)then
              write(*,302)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D17)then
              write(*,303)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D16)then
              write(*,304)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D15)then
              write(*,305)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D14)then
              write(*,306)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D13)then
              write(*,307)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D12)then
              write(*,308)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D11)then
              write(*,309)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D10)then
              write(*,310)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D9)then
              write(*,311)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D8)then
              write(*,312)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D7)then
              write(*,313)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D6)then
              write(*,314)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D5)then
              write(*,315)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D4)then
              write(*,316)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D3)then
              write(*,317)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D2)then
              write(*,318)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D1)then
              write(*,319)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D0)then
              write(*,320)nn,chi1
              return
            endif
            write(*,321)nn,chi1
          endif !nn > 100  < 1000
c
  
100       format('med',i1,' = 999999999999999999999.999999')
101       format('med',i1,' = ',f28.6)
102       format('med',i1,' = ',f27.6)
103       format('med',i1,' = ',f26.6)
104       format('med',i1,' = ',f25.6)
105       format('med',i1,' = ',f24.6)
106       format('med',i1,' = ',f23.6)
107       format('med',i1,' = ',f22.6)
108       format('med',i1,' = ',f21.6)
109       format('med',i1,' = ',f20.6)
110       format('med',i1,' = ',f19.6)
111       format('med',i1,' = ',f18.6)
112       format('med',i1,' = ',f17.6)
113       format('med',i1,' = ',f16.6)
114       format('med',i1,' = ',f15.6)
115       format('med',i1,' = ',f14.6)
116       format('med',i1,' = ',f13.6)
117       format('med',i1,' = ',f12.6)
118       format('med',i1,' = ',f11.6)
119       format('med',i1,' = ',f10.6)
120       format('med',i1,' = ',f9.6)
121       format('med',i1,' = ',f8.6)
c
200       format('med',i2,' = 999999999999999999999.999999')
201       format('med',i2,' = ',f28.6)
202       format('med',i2,' = ',f27.6)
203       format('med',i2,' = ',f26.6)
204       format('med',i2,' = ',f25.6)
205       format('med',i2,' = ',f24.6)
206       format('med',i2,' = ',f23.6)
207       format('med',i2,' = ',f22.6)
208       format('med',i2,' = ',f21.6)
209       format('med',i2,' = ',f20.6)
210       format('med',i2,' = ',f19.6)
211       format('med',i2,' = ',f18.6)
212       format('med',i2,' = ',f17.6)
213       format('med',i2,' = ',f16.6)
214       format('med',i2,' = ',f15.6)
215       format('med',i2,' = ',f14.6)
216       format('med',i2,' = ',f13.6)
217       format('med',i2,' = ',f12.6)
218       format('med',i2,' = ',f11.6)
219       format('med',i2,' = ',f10.6)
220       format('med',i2,' = ',f9.6)
221       format('med',i2,' = ',f8.6)
c
c
300       format('med',i3,' = 999999999999999999999.999999')
301       format('med',i3,' = ',f28.6)
302       format('med',i3,' = ',f27.6)
303       format('med',i3,' = ',f26.6)
304       format('med',i3,' = ',f25.6)
305       format('med',i3,' = ',f24.6)
306       format('med',i3,' = ',f23.6)
307       format('med',i3,' = ',f22.6)
308       format('med',i3,' = ',f21.6)
309       format('med',i3,' = ',f20.6)
310       format('med',i3,' = ',f19.6)
311       format('med',i3,' = ',f18.6)
312       format('med',i3,' = ',f17.6)
313       format('med',i3,' = ',f16.6)
314       format('med',i3,' = ',f15.6)
315       format('med',i3,' = ',f14.6)
316       format('med',i3,' = ',f13.6)
317       format('med',i3,' = ',f12.6)
318       format('med',i3,' = ',f11.6)
319       format('med',i3,' = ',f10.6)
320       format('med',i3,' = ',f9.6)
321       format('med',i3,' = ',f8.6)

c
          return
          end
c
c  ##########
c

          subroutine printamoebachi(nn,chi1)
c
          implicit double precision(a-h,o-z)
c
          if(nn.lt.10)then
            if(chi1.gt.9.999999999999999D20)then
              write(*,100)nn
              return
            endif
            if(chi1.gt.9.999999999999999D19)then
              write(*,101)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D18)then
              write(*,102)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D17)then
              write(*,103)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D16)then
              write(*,104)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D15)then
              write(*,105)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D14)then
              write(*,106)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D13)then
              write(*,107)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D12)then
              write(*,108)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D11)then
              write(*,109)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D10)then
              write(*,110)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D9)then
              write(*,111)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D8)then
              write(*,112)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D7)then
              write(*,113)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D6)then
              write(*,114)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D5)then
              write(*,115)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D4)then
              write(*,116)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D3)then
              write(*,117)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D2)then
              write(*,118)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D1)then
              write(*,119)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D0)then
              write(*,120)nn,chi1
              return
            endif
            write(*,121)nn,chi1
          endif !nn < 10
c
          if(nn.ge.10)then
            if(chi1.gt.9.999999999999999D20)then
              write(*,200)nn
              return
            endif
            if(chi1.gt.9.999999999999999D19)then
              write(*,201)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D18)then
              write(*,202)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D17)then
              write(*,203)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D16)then
              write(*,204)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D15)then
              write(*,205)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D14)then
              write(*,206)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D13)then
              write(*,207)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D12)then
              write(*,208)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D11)then
              write(*,209)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D10)then
              write(*,210)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D9)then
              write(*,211)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D8)then
              write(*,212)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D7)then
              write(*,213)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D6)then
              write(*,214)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D5)then
              write(*,215)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D4)then
              write(*,216)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D3)then
              write(*,217)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D2)then
              write(*,218)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D1)then
              write(*,219)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D0)then
              write(*,220)nn,chi1
              return
            endif
            write(*,221)nn,chi1
          endif !nn < 10
  
100       format('chi',i1,' = 999999999999999999999.999999')
101       format('chi',i1,' = ',f28.6)
102       format('chi',i1,' = ',f27.6)
103       format('chi',i1,' = ',f26.6)
104       format('chi',i1,' = ',f25.6)
105       format('chi',i1,' = ',f24.6)
106       format('chi',i1,' = ',f23.6)
107       format('chi',i1,' = ',f22.6)
108       format('chi',i1,' = ',f21.6)
109       format('chi',i1,' = ',f20.6)
110       format('chi',i1,' = ',f19.6)
111       format('chi',i1,' = ',f18.6)
112       format('chi',i1,' = ',f17.6)
113       format('chi',i1,' = ',f16.6)
114       format('chi',i1,' = ',f15.6)
115       format('chi',i1,' = ',f14.6)
116       format('chi',i1,' = ',f13.6)
117       format('chi',i1,' = ',f12.6)
118       format('chi',i1,' = ',f11.6)
119       format('chi',i1,' = ',f10.6)
120       format('chi',i1,' = ',f9.6)
121       format('chi',i1,' = ',f8.6)
c
200       format('chi',i2,' = 999999999999999999999.999999')
201       format('chi',i2,' = ',f28.6)
202       format('chi',i2,' = ',f27.6)
203       format('chi',i2,' = ',f26.6)
204       format('chi',i2,' = ',f25.6)
205       format('chi',i2,' = ',f24.6)
206       format('chi',i2,' = ',f23.6)
207       format('chi',i2,' = ',f22.6)
208       format('chi',i2,' = ',f21.6)
209       format('chi',i2,' = ',f20.6)
210       format('chi',i2,' = ',f19.6)
211       format('chi',i2,' = ',f18.6)
212       format('chi',i2,' = ',f17.6)
213       format('chi',i2,' = ',f16.6)
214       format('chi',i2,' = ',f15.6)
215       format('chi',i2,' = ',f14.6)
216       format('chi',i2,' = ',f13.6)
217       format('chi',i2,' = ',f12.6)
218       format('chi',i2,' = ',f11.6)
219       format('chi',i2,' = ',f10.6)
220       format('chi',i2,' = ',f9.6)
221       format('chi',i2,' = ',f8.6)

c
          return
          end
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printiter1(imod,instring1)
c
          character*(*)instring1
          character*(80) outstring1
c
          if(imod.lt.0)return
          if(imod.ge.1000000)return
c
          l1=lnblnk(instring1)
          l1=l1+1
          if(imod.lt.10)write(outstring1,100)instring1(1:l1),imod
          if((imod.ge.10).and.(imod.lt.100))then
            write(outstring1,101)instring1(1:l1),imod
          endif
          if((imod.ge.100).and.(imod.lt.1000))then
            write(outstring1,102)instring1(1:l1),imod
          endif
          if((imod.ge.1000).and.(imod.lt.10000))then
            write(outstring1,103)instring1(1:l1),imod
          endif
          if((imod.ge.10000).and.(imod.lt.100000))then
            write(outstring1,104)instring1(1:l1),imod
          endif
          if((imod.ge.100000).and.(imod.lt.1000000))then
            write(outstring1,105)instring1(1:l1),imod
          endif
c
          ll1=lnblnk(outstring1)
c
          write(*,200)outstring1(1:ll1)
c
 100      format(a,i1)
 101      format(a,i2)
 102      format(a,i3)
 103      format(a,i4)
 104      format(a,i5)
 105      format(a,i6)
 200      format(a)
          return
          end
c
c
c @#$&^*&*$*&@&%#&%#*&^%#^$@%^(&^%@$&@@@@#$&^*&*$*&@&%#&%#*&^%#^$@%^(&^%
c

          subroutine getinput(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,
     @       omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,
     $       betarim,rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,
     &       Nref,rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,
     @       sw3,T0,idraw,iecheck,iidint,iatm,ism1,icnU,icnB,icnV,icnR,
     @       icnI,icnJ,icnH,icnK,iRVfilt,isw1,isw2,isw3,isw4,ilaw,wave,
     @       dbolx,dboly,dwavex,dwavey,ecc,argper,pshift,sw5,sw6,sw7,
     @       sw8,sw9,ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,
     @       spot2parm,spotdparm,primmass,primK,primrad,ratrad,frac1,
     @       frac2,ecosw,temprat,idark1,idark2,isw12,isw13,isw21,isw22,
     @       isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,
     @       sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,isw25,isw26,
     @       isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,ocose,
     @       osine,omegadot,contamS0,contamS1,contamS2,contamS3,
     @       sw47,sw48,sw49)
c
          implicit double precision(a-h,o-z)
c
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,3),
     @      dwavey(8,3)
c
c   RVG BUG ALERT  May 9, 2001
c
c   Dimension the spot arrays.
c
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          dimension powercoeff(8,9)
c
c   UPDATE June 22, 2002
c
c   Declare the variable bell to be character*1
c
          character*1 bell
c
          ios=0
          open(unit=1,file='ELC.inp',status='old',err=100,iostat=ios)
c
          read(1,*)Nalph1
          read(1,*)Nbet1
          read(1,*)Nalph2
          read(1,*)Nbet2
          read(1,*)fill1
          read(1,*)fill2
          read(1,*)omega1
          read(1,*)omega2
          read(1,*)dphase
          read(1,*)Q
          read(1,*)finc
          read(1,*)Teff1
          read(1,*)Teff2
          read(1,*)Tgrav1
          read(1,*)Tgrav2
          read(1,*)betarim
          read(1,*)rinner
          read(1,*)router
          read(1,*)tdisk
          read(1,*)xi
          read(1,*)Ntheta
          read(1,*)Nradius
          read(1,*)alb1
          read(1,*)alb2
          read(1,*)Nref
          read(1,*)rLx
          read(1,*)Period
          read(1,*)fm
          read(1,*)separ
          read(1,*)gamma
          read(1,*)t3
          read(1,*)g3
          read(1,*)SA3
          read(1,*)density
          read(1,*)sw1
          read(1,*)sw2
          read(1,*)sw3
          read(1,*)T0
          read(1,*)idraw
          read(1,*)iecheck
          read(1,*)iidint
          read(1,*)iatm
          read(1,*)ism1
          read(1,*)icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK
          read(1,*)iRVfilt
          read(1,*)isw1
          read(1,*)isw2
          read(1,*)isw3
          read(1,*)isw4
          read(1,*)ilaw
c
c  Load the limb darkening parameters. 
c         
          do 10 i=1,8
            read(1,*)wave(i),dbolx(i,1),dboly(i,1),dbolx(i,2),dboly(i,2),
     %         dwavex(i,1),dwavey(i,1),dwavex(i,2),dwavey(i,2)
            dwavex(i,3)=dwavex(i,1)
            dwavey(i,3)=dwavey(i,1)
 10                continue
c
c  Look for additional parameters related to the eccentricity.  Use
c  the end=...  option to allow for older input files.
c
          read(1,*,end=99)ecc
          read(1,*)argper
          read(1,*)pshift
          read(1,*)sw5
          read(1,*)sw6
          read(1,*)sw7
          read(1,*)sw8
c
c   UPDATE April 8, 2002
c
c   Correct bug:  change sw5 to sw9 below
c
          read(1,*)sw9
          read(1,*)ikeep
          read(1,*)isynch
          read(1,*)isw5
          read(1,*)isw6
          read(1,*)isw7
          read(1,*)isw8
          read(1,*)isw9
          ios=0
          read(1,*,end=101,err=101)spot1parm(1,1)
          read(1,*,end=101,err=101)spot1parm(1,2)
          read(1,*,end=101,err=101)spot1parm(1,3)
          read(1,*,end=101,err=101)spot1parm(1,4)
c              
          read(1,*,end=101,err=101)spot1parm(2,1)
          read(1,*,end=101,err=101)spot1parm(2,2)
          read(1,*,end=101,err=101)spot1parm(2,3)
          read(1,*,end=101,err=101)spot1parm(2,4)
c              
c              
          read(1,*,end=101,err=101)spot2parm(1,1)
          read(1,*,end=101,err=101)spot2parm(1,2)
          read(1,*,end=101,err=101)spot2parm(1,3)
          read(1,*,end=101,err=101)spot2parm(1,4)
               
          read(1,*,end=101,err=101)spot2parm(2,1)
          read(1,*,end=101,err=101)spot2parm(2,2)
          read(1,*,end=101,err=101)spot2parm(2,3)
          read(1,*,end=101,err=101)spot2parm(2,4)
c              
c              
          read(1,*,end=101,err=101)spotdparm(1,1)
          read(1,*,end=101,err=101)spotdparm(1,2)
          read(1,*,end=101,err=101)spotdparm(1,3)
          read(1,*,end=101,err=101)spotdparm(1,4)
               
          read(1,*,end=101,err=101)spotdparm(2,1)
          read(1,*,end=101,err=101)spotdparm(2,2)
          read(1,*,end=101,err=101)spotdparm(2,3)
          read(1,*,end=101,err=101)spotdparm(2,4)
c
c   UPDATE August 10, 2004
c
c   Add the new variables here.  Abort if there is an error
c
          read(1,*,end=101,err=101)primmass
          read(1,*,end=101,err=101)primK
          read(1,*,end=101,err=101)primrad
          read(1,*,end=101,err=101)ratrad
          read(1,*,end=101,err=101)frac1
          read(1,*,end=101,err=101)frac2
          read(1,*,end=101,err=101)ecosw
          read(1,*,end=101,err=101)temprat
          read(1,*,end=101,err=101)idark1
          read(1,*,end=101,err=101)idark2
          read(1,*,end=101,err=101)isw12
          read(1,*,end=101,err=101)isw13
c
c   UPDATE May 8, 2006
c
c   Add new variables here.  Abort if there is an error.
c
          read(1,*,end=101,err=101)isw21
          read(1,*,end=101,err=101)isw22
          read(1,*,end=101,err=101)isw23
          read(1,*,end=101,err=101)isw24
c
c   Power-series limb darkening coefficients for star 1
c
          read(1,*,end=101,err=101)(powercoeff(1,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(2,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(3,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(4,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(5,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(6,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(7,k),k=1,9)
          read(1,*,end=101,err=101)(powercoeff(8,k),k=1,9)
c
          read(1,*,end=101,err=101)bigI
          read(1,*,end=101,err=101)bigbeta
          read(1,*,end=101,err=101)sw23
          read(1,*,end=101,err=101)sw24
          read(1,*,end=101,err=101)sw25
          read(1,*,end=101,err=101)sw26
          read(1,*,end=101,err=101)sw27
          read(1,*,end=101,err=101)sw28
          read(1,*,end=101,err=101)sw29
          read(1,*,end=101,err=101)sw30
          read(1,*,end=101,err=101)contam
          read(1,*,end=101,err=101)Tconj
          read(1,*,end=101,err=101)beam1
          read(1,*,end=101,err=101)beam2

          read(1,*,end=101,err=101)isw25
          read(1,*,end=101,err=101)isw26
          read(1,*,end=101,err=101)isw27
          read(1,*,end=101,err=101)isw28
          read(1,*,end=101,err=101)isw29
          read(1,*,end=101,err=101)isw30
          read(1,*,end=101,err=101)isw31
          read(1,*,end=101,err=101)isw32
          read(1,*,end=101,err=101)isw33
          read(1,*,end=101,err=101)isw34
          read(1,*,end=101,err=101)ocose
          read(1,*,end=101,err=101)osine
          read(1,*,end=101,err=101)omegadot
          read(1,*,end=101,err=101)contamS0
          read(1,*,end=101,err=101)contamS1
          read(1,*,end=101,err=101)contamS2
          read(1,*,end=101,err=101)contamS3
          read(1,*,end=101,err=101)sw47
          read(1,*,end=101,err=101)sw48
          read(1,*,end=101,err=101)sw49

 99       close(1)
 100      if(ios.gt.0)then
            write(*,*)'error with ELC.inp'
            stop
          endif
c
c   Put this if-then block for successful completion.
c
          if(ios.eq.0)then 
            close(1)
            return
          endif

 101      bell=char(7)   !file ended too soon
          write(*,1002)bell
          close(1)
c
 1002     format(a1,'Error:  Bad entry in ELC.inp')
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&
c
