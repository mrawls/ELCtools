c
c    November 12, 1999
c
c    These subroutines will perform the necessary tasks for the optimizer
c    routines, such as variable assignment, chi^2 checking, etc.
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
            svar(i)='none'
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
c
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
            call makeloopopt(Nvmax)
          endif
c
 100      format(a40)
 200      format(a1,'Error:  too many variables in ''gridloop.opt''')
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine makeloopopt(Nvmax)
c
c   This routine will write a correctly formatted 'gridloop.opt' file
c   in the event one does not exist.
c     
          implicit double precision (a-h,o-z)

          character*1 bell
c
          bell=char(7)
          write(*,100)bell
c 
          open(unit=11,file='gridloop.opt',status='new')
c
          write(11,1000)
          write(11,1001)
          write(11,1002)
          write(11,1003)
          write(11,1004)
          write(11,1005)
          write(11,1006)
          write(11,1007)
          write(11,1008)
          write(11,1009)
c
          Nvar=Nvmax
          write(11,1010)Nvar
c
          write(11,1011)
          write(11,1012)
          write(11,1013)
          write(11,1014)
          write(11,1015)
          write(11,1016)
          write(11,1017)
          write(11,1018)
          write(11,1019)
          write(11,1020)
          write(11,1021)
          write(11,1022)
          write(11,1023)
          write(11,1024)
          write(11,1025)
          write(11,1026)
          write(11,1030)
          write(11,1031)
          write(11,1032)
          write(11,1033)
          write(11,1034)
          write(11,1035)
          write(11,1036)
          write(11,1037)
          write(11,1038)
          write(11,1039)
          write(11,1040)
          write(11,1041)
          write(11,1042)
          write(11,1043)
          write(11,1044)
          write(11,1045)
c
          write(11,2000)
          write(11,2001)
          write(11,2002)
          write(11,2003)
          write(11,3004)
          write(11,2004)
          write(11,2005)
          write(11,2006)
          write(11,2007)
          write(11,2008)
          write(11,2011)
          write(11,2012)
          write(11,2014)
          write(11,2015)
c
          idummy=10
          if(idummy.eq.10)then
            write(*,2020)
            stop
          endif
c
 100      format(a1,'Error:  file ''gridloop.opt'' not found!  ',
     %      'I''m making one up!')
c
 1000     format('put_your_file_for_U_data_here')
 1001     format('put_your_file_for_B_data_here')
 1002     format('put_your_file_for_V_data_here')
 1003     format('put_your_file_for_R_data_here')
 1004     format('put_your_file_for_I_data_here')
 1005     format('put_your_file_for_J_data_here')
 1006     format('put_your_file_for_H_data_here')
 1007     format('put_your_file_for_K_data_here')
 1008     format('put_your_file_for_RV1_data_here')
 1009     format('put_your_file_for_RV2_data_here')
 1010     format(i2,18x,'Number of variables to adjust')
 1011     format('inclination')
 1012     format('mass ratio')
 1013     format('f1 (fill 1)')
 1014     format('f2 (fill 2)')
 1015     format('o1 (omega1)')
 1016     format('o2 (omega2)')
 1017     format('rinner [inner disk radius (same units as fill2)]')
 1018     format('router [outer disk radius','
     %              (in units of Rl volume of star 2)]')
 1019     format('Tdisk (inner disk temperature)')
 1020     format('betarim (half-opening angle of disk rim in degrees)')
 1021     format('T1 (T_eff of star 1)')
 1022     format('T2 (T_eff of star 2)')
 1023     format('xi')
 1024     format('Lx (L_x/L_opt)')
 1025     format('separation (solar radii)')
 1026     format('gamma (gamma velocity in km/sec)')
 1030     format('70.0     1.0     3')
 1031     format('2.0      0.25    3')
 1032     format('0.70     0.1     3')
 1033     format('0.70     0.1     3')
 1034     format('1.0      0.1     3')
 1035     format('1.0      0.1     3')
 1036     format('0.01     0.001   3')
 1037     format('0.75     0.01    3')
 1038     format('20000.0  100.0   3')
 1039     format('4.0      0.5     3')
 1040     format('6500.0   100.0   3')
 1041     format('7500.0   100.0   3')
 1042     format('-0.75    0.05    3')
 1043     format('10.0     1.0     3')
 1044     format('5.0      0.5     3')
 1045     format('-69.0    10.0    3')
c
 2000   format('##############################################')
 2001   format('#    The first 8 entry MUST be a string indicating ',
     $      'the name of the FOLDED')
 2002   format('#    data files.  Put ''none'' if there is no file.')
 2003   format('#')
 3004   format('#    The next entry MUST be a positive integer (Nvar)')
 2004   format('#    The next Nvar entries MUST be strings.  Put the name ',
     $      'of the variable you want ')
 2005   format('#    in the innermost loop first, the name of the ', 
     $      'variable in the second') 
 2006   format('#    loop second, etc.',  
     $      ' If you do not want to loop over all variables, ')
 2007   format('#    put the string ''NOTHING'' in the unwanted positions.')
 2008   format('#    The names of the allowed variables that can be ',
     $      'put into loops are listed')
 2011   format('#    The next Nvar lines contain the starting values',
     $      ' of the variables') 
 2012   format('#    the step size, and the number of loops.')
 2014   format('#')
 2015   format('##############################################')
c
 2020   format('Edit the file ''gridloop.opt'' and restart the program.')
c
          return
          end
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

      RETURN
C
      END!
C
c    End of stolen code.
c
c
c    *************************************************************
c
          subroutine loaddata(Ndatamax,filein,Ndata,xdata,ydata,err)
c
c    November 12, 1999
c
c    This general routine will load a data file.  It assumes that there
c    are three columns:   phase    flux or RV     error
c
          implicit double precision (a-h,o-z)
c
          dimension xdata(Ndatamax),ydata(Ndatamax),err(Ndatamax)
c
          character*40 filein
          character*1 bell
c
          bell=char(7)
          ios=0
c
          Ndata=0
          open(unit=20,file=filein,status='old',err=999,iostat=ios)
c
          do 10 i=1,Ndatamax
            read(20,*,end=15)xdata(i),ydata(i),err(i)
 10       continue
 15       close(20)
c
          Ndata=i-1
c
 999      if(ios.ne.0)then
            write(*,10000)bell,filein
            stop
          endif
c
10000     format(a1,'Error:  data file not found (',a40,')')
c  
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  RVG BUG ALERT   May 8, 2001
c
c  This subroutine has been updated to include the spot parameters
c
c  NEW BUG August 2, 2001
c
c  This routine has been updated to include the period and phase zeropoint.
c
c
c  UPDATE January 15, 2002
c
c  Update the routine to include the albedos of the two stars (alb1, alb2)
c
c  UPDATE November 6, 2002
c
c  Add limb darkening coefficients dwavex and dwavey to the
c  argument list of assignvar and varassign.
c
c  UPDATE August 13, 2004
c
c  Add primmass,primK,primrad,ratrad,frac1,frac2
c
c  UPDATE JULY 29, 2005
c
c  Add ecosw and temprat to the list.
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c  UPDATE October 10, 2008
c
c  Add the density.
c
          subroutine assignvar(Nvarmax,svar,var,fill1,fill2,omega1,
     @       omega2,Q,finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,
     @       rLx,separ,gamma,t3,g3,sa3,ecc,argper,pshift,spot1parm,
     @       spot2parm,spotdparm,period,T0,alb1,alb2,dwavex,dwavey,
     @       primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @       bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2,contam,
     @       ocose,osine,isw29,tertperiod,tertt0,tertecos,tertesin,
     @       tertincl,tertOmega,tertQ,Tgrav1,Tgrav2,tertconj,omegadot,
     @       contamS0,contamS1,contamS2,contamS3,
     @       P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @       P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @       P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @       P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @       P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @       P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @       P6ratrad,
     @       P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
     @       P7ratrad,
     @       P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @       P8ratrad,sw72,sw73)
c
c   November 12, 1999
c
c   This routine will determine which variables need to be changed
c   based on the string codes in svar(1:Nvarmax)
c
c          implicit double precision (a-h,o-z)
c
          implicit none

          integer nvarmax,isw29,i,kk,icnvrt

          real*8 var
          real*8 fill1,fill2,omega1,
     @       omega2,Q,finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,
     @       rLx,separ,gamma,t3,g3,sa3,ecc,argper,pshift,spot1parm,
     @       spot2parm,spotdparm,period,T0,alb1,alb2,dwavex,dwavey,
     @       primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @       bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2,contam,
     @       ocose,osine,sw72,sw73
          real*8 tertperiod,tertt0,tertecos,tertesin,
     @       tertincl,tertOmega,tertQ,Tgrav1,Tgrav2,tertconj,omegadot,
     @       contamS0,contamS1,contamS2,contamS3,
     @       P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @       P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @       P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @       P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @       P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @       P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @       P6ratrad,
     @       P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
     @       P7ratrad,
     @       P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @       P8ratrad

          character*40 svar(Nvarmax)
          dimension var(Nvarmax)
c
c   RVG BUG ALERT   May 8, 2001
c
c   Dimension the spot parameter arrays.
c       
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
c
c
c  UPDATE November 6, 2002
c
c  dimension the limb darkening variables
c
          dimension dwavex(8,3),dwavey(8,3)
          dimension powercoeff(8,9)

          write(*,*)' '
        
          do 10 i=1,Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
c   tidal apsidal coefficients
c
            if(kk.eq.1016)then   !a1, sw72
              sw72=var(i)
              if(sw72.lt.0.0d0)sw72=0.0d0
            endif

            if(kk.eq.1017)then   !a2, sw73
              sw73=var(i)
              if(sw73.lt.0.0d0)sw73=0.0d0
            endif
c
c  planet 2 parameters
c
            if(kk.eq.649)then  !P2tconj, tag uj
              P2tconj=var(i)
            endif
c
            if(kk.eq.659)then  !P2period, tag ut
              P2period=var(i)
            endif
c
            if(kk.eq.660)then  !P2period, tag uu
              P2T0=var(i)
            endif
c
            if(kk.eq.661)then  !P2ecos, tag uv
              P2ecos=var(i)
            endif
c
            if(kk.eq.662)then  !P2esin, tag uw
              P2esin=var(i)
            endif
c
            if(kk.eq.663)then  !P2incl, tag ux
              P2incl=var(i)
            endif
c
            if(kk.eq.664)then  !P2Omega, tag uy
              P2Omega=var(i)
            endif
c
            if(kk.eq.665)then  !P2Q, tag uz
              P2Q=var(i)
            endif
c
            if(kk.eq.641)then  !P2ratrad, tag ub
              P2ratrad=var(i)
            endif
c
c planet 3 parameters
c
            if(kk.eq.681)then  !P3tconj, tag vj
              P3tconj=var(i)
            endif
c
            if(kk.eq.691)then  !P3period, tag vt
              P3period=var(i)
            endif
c
            if(kk.eq.692)then  !P3period, tag vu
              P3T0=var(i)
            endif
c
            if(kk.eq.693)then  !P3ecos, tag vv
              P3ecos=var(i)
            endif
c
            if(kk.eq.694)then  !P3esin, tag vw
              P3esin=var(i)
            endif
c
            if(kk.eq.695)then  !P3incl, tag vx
              P3incl=var(i)
            endif
c
            if(kk.eq.696)then  !P3Omega, tag vy
              P3Omega=var(i)
            endif
c
            if(kk.eq.697)then  !P3Q, tag vz
              P3Q=var(i)
            endif
c
            if(kk.eq.673)then  !P3ratrad, tag vb
              P3ratrad=var(i)
            endif
c
c  planet 4 parameters
c
            if(kk.eq.713)then  !P4tconj, tag wj
              P4tconj=var(i)
            endif
c
            if(kk.eq.723)then  !P4period, tag wt
              P4period=var(i)
            endif
c
            if(kk.eq.724)then  !P4period, tag wu
              P4T0=var(i)
            endif
c
            if(kk.eq.725)then  !P4ecos, tag wv
              P4ecos=var(i)
            endif
c
            if(kk.eq.726)then  !P4esin, tag ww
              P4esin=var(i)
            endif
c
            if(kk.eq.727)then  !P4incl, tag wx
              P4incl=var(i)
            endif
c
            if(kk.eq.728)then  !P4Omega, tag wy
              P4Omega=var(i)
            endif
c
            if(kk.eq.729)then  !P4Q, tag wz
              P4Q=var(i)
            endif
c
            if(kk.eq.705)then  !P4ratrad, tag wb
              P4ratrad=var(i)
            endif
c
c  planet 5 parameters
c
            if(kk.eq.745)then  !P5tconj, tag xj
              P5tconj=var(i)
            endif
c
            if(kk.eq.755)then  !P5period, tag xt
              P5period=var(i)
            endif
c
            if(kk.eq.756)then  !P5period, tag xu
              P5T0=var(i)
            endif
c
            if(kk.eq.757)then  !P5ecos, tag xv
              P5ecos=var(i)
            endif
c
            if(kk.eq.758)then  !P5esin, tag xw
              P5esin=var(i)
            endif
c
            if(kk.eq.759)then  !P5incl, tag xx
              P5incl=var(i)
            endif
c
            if(kk.eq.760)then  !P5Omega, tag xy
              P5Omega=var(i)
            endif
c
            if(kk.eq.761)then  !P5Q, tag xz
              P5Q=var(i)
            endif
c
            if(kk.eq.737)then  !P5ratrad, tag xb
              P5ratrad=var(i)
            endif
c
c  planet 6 parameters
c
            if(kk.eq.585)then  !P6tconj, tag sj
              P6tconj=var(i)
            endif
c
            if(kk.eq.595)then  !P6period, tag st
              P6period=var(i)
            endif
c
            if(kk.eq.596)then  !P6period, tag su
              P6T0=var(i)
            endif
c
            if(kk.eq.597)then  !P6ecos, tag sv
              P6ecos=var(i)
            endif
c
            if(kk.eq.598)then  !P6esin, tag sw
              P6esin=var(i)
            endif
c
            if(kk.eq.599)then  !P6incl, tag sx
              P6incl=var(i)
            endif
c
            if(kk.eq.600)then  !P6Omega, tag sy
              P6Omega=var(i)
            endif
c
            if(kk.eq.601)then  !P6Q, tag sz
              P6Q=var(i)
            endif
c
            if(kk.eq.577)then  !P6ratrad, tag sb
              P6ratrad=var(i)
            endif
c
c  planet 7 parameters
c
            if(kk.eq.233)then  !P7tconj, tag hj
              P7tconj=var(i)
            endif
c
            if(kk.eq.243)then  !P7period, tag ht
              P7period=var(i)
            endif
c
            if(kk.eq.244)then  !P7period, tag hu
              P7T0=var(i)
            endif
c
            if(kk.eq.245)then  !P7ecos, tag hv
              P7ecos=var(i)
            endif
c
            if(kk.eq.246)then  !P7esin, tag hw
              P7esin=var(i)
            endif
c
            if(kk.eq.247)then  !P7incl, tag hx
              P7incl=var(i)
            endif
c
            if(kk.eq.248)then  !P7Omega, tag hy
              P7Omega=var(i)
            endif
c
            if(kk.eq.249)then  !P7Q, tag hz
              P7Q=var(i)
            endif
c
            if(kk.eq.225)then  !P7ratrad, tag hb
              P7ratrad=var(i)
            endif
c
c  planet 8 parameters
c
            if(kk.eq.329)then  !P8tconj, tag kj
              P8tconj=var(i)
            endif
c
            if(kk.eq.339)then  !P8period, tag kt
              P8period=var(i)
            endif
c
            if(kk.eq.340)then  !P8period, tag ku
              P8T0=var(i)
            endif
c
            if(kk.eq.341)then  !P8ecos, tag kv
              P8ecos=var(i)
            endif
c
            if(kk.eq.342)then  !P8esin, tag kw
              P8esin=var(i)
            endif
c
            if(kk.eq.343)then  !P8incl, tag kx
              P8incl=var(i)
            endif
c
            if(kk.eq.344)then  !P8Omega, tag ky
              P8Omega=var(i)
            endif
c
            if(kk.eq.345)then  !P8Q, tag kz
              P8Q=var(i)
            endif
c
            if(kk.eq.42)then  !P8ratrad, tag kb
              P8ratrad=var(i)
            endif
c
c   seasonal contamination
c
            if(kk.eq.1591)then    !contamS0, tag s0
              contamS0=var(i)
              if(contamS0.lt.0.0d0)contamS0=0.0d0
            endif

            if(kk.eq.1592)then    !contamS1, tag s1
              contamS1=var(i)
              if(contamS1.lt.0.0d0)contamS1=0.0d0
            endif

            if(kk.eq.1593)then    !contamS2, tag s2
              contamS2=var(i)
              if(contamS2.lt.0.0d0)contamS2=0.0d0
            endif

            if(kk.eq.1594)then    !contamS3, tag s3
              contamS3=var(i)
              if(contamS3.lt.0.0d0)contamS3=0.0d0
            endif

            if(kk.eq.110)then  ! omegadot, tag do
              omegadot=var(i)
            endif
c
            if(kk.eq.1208)then  ! Tgrav1, tag g1
              Tgrav1=var(i)
            endif
c
            if(kk.eq.1209)then  ! Tgrav2, tag g2
              Tgrav2=var(i)
            endif
c
c   November 18, 2012
c
c   Add third body stuff here.
c
            if(kk.eq.617)then   !tertperiod tag tj
              tertconj=var(i)
            endif
c
            if(kk.eq.627)then   !tertperiod tag tt
              tertperiod=var(i)
            endif
c
            if(kk.eq.628)then   !tertT0 tag tu
              tertT0=var(i)
            endif
c
            if(kk.eq.629)then   !tertecos tag tv
              tertecos=var(i)
            endif
c
            if(kk.eq.630)then   !tertesin tag tw
              tertesin=var(i)
            endif
c
            if(kk.eq.631)then   !tertincl tag tx
              tertincl=var(i)
            endif
c
            if(kk.eq.632)then   !tertOmega tag ty
              tertOmega=var(i)
            endif
c
            if(kk.eq.633)then   !tertQ tag tz
              tertQ=var(i)
            endif
c


            if((isw29.gt.0).and.(kk.eq.450))then ! ecos(omega) oc
              ocose=var(i)
            endif

            if((isw29.gt.0).and.(kk.eq.466))then ! ecos(omega) os
              osine=var(i)
            endif
c
            if(kk.eq.78)then   !contam, use string co
              contam=var(i)
            endif
c
            if(kk.eq.1144)then   !beam1, use string e1
              beam1=var(i)
            endif
c
            if(kk.eq.1145)then   !beam2, use string e2
              beam2=var(i)
            endif
c
            if(kk.eq.610)then   !Tconj, use string tc
              Tconj=var(i)
            endif
c
            if(kk.eq.100)then   !density, use string de
              density=var(i)
            endif
c
            if(kk.eq.8)then  !bigI, use string ai
              bigI=var(i)
            endif
c
            if(kk.eq.1)then  !bigbeta, use string ab
              bigbeta=var(i)
            endif
c
            if(kk.eq.111)then  !ecosw, use string dphi
              ecosw=var(i)
            endif
c
            if(kk.eq.612)then  !temprat
              temprat=var(i)
            endif
c
            if(kk.eq.492)then  !primmass, use string pm
              primmass=var(i)
              if(primmass.lt.0.0d0)primmass=0.0d0
            endif
c
            if(kk.eq.497)then  !primrad, use string pr
              primrad=var(i)
              if(primrad.lt.0.0d0)primrad=0.0d0
            endif
c
            if(kk.eq.490)then  !primK, use string pk
              primK=var(i)
              if(primK.lt.0.0d0)primK=0.0d0
            endif
c
            if(kk.eq.544)then  !ratrad, use string ra
              ratrad=var(i)
              if(ratrad.lt.0.0d0)ratrad=0.0d0
            endif
c
            if(kk.eq.1528)then  !frac1, use string q1
              frac1=var(i)
              if(frac1.lt.0.0d0)frac1=0.0d0
            endif
c
            if(kk.eq.1529)then  !frac2, use string q2
              frac2=var(i)
c              if(frac2.lt.0.0d0)frac2=0.0d0
            endif
c
c
c   UPDATE January 15, 2002
c
c   Here are the assignments for the albedos
c
            if(kk.eq.1368)then  !alb1, use string l1
              alb1=var(i)
              if(alb1.lt.0.0d0)alb1=0.0d0
            endif
c
            if(kk.eq.1369)then  !alb2, use string l2
              alb2=var(i)
              if(alb2.lt.0.0d0)alb2=0.0d0
            endif
c
c   NEW BUG August 2, 2001
c
c   Here are the assignments for the period and T0
c
            if(kk.eq.484)then  !period
              period=var(i)
            endif
c
            if(kk.eq.1623)then ! T0
              T0=var(i)
            endif
c
c   RVG BUG ALERT   May 8, 2001
c
c   Here is the new block for the assignments of spot parameters.
c
            if(kk.eq.1112)then   ! temperature factor spot 1 on disk
              spotdparm(1,1)=var(i)
              if(spotdparm(1,1).le.-10.0d0)spotdparm(1,1)=-1.0d0
            endif
c
            if(kk.eq.1116)then   ! temperature factor spot 2 on disk
              spotdparm(2,1)=var(i)
              if(spotdparm(2,1).le.-10.0d0)spotdparm(2,1)=-1.0d0
            endif
c
            if(kk.eq.1113)then  ! azimuth spot 1 on disk
              spotdparm(1,2)=dmod(var(i),360.0d0)
            endif   
c
            if(kk.eq.1117)then  ! azimuth spot 2 on disk
              spotdparm(2,2)=dmod(var(i),360.0d0)
            endif   
c
            if(kk.eq.1114)then    ! cutoff radius for spot 1 on disk
              spotdparm(1,3)=var(i)
              if(spotdparm(1,3).gt.1.0d0)spotdparm(1,3)=1.0d0
              if(spotdparm(1,3).lt.0.0d0)spotdparm(1,3)=0.0d0
            endif
c
            if(kk.eq.1118)then    ! cutoff radius for spot 2 on disk
              spotdparm(2,3)=var(i)
              if(spotdparm(2,3).gt.1.0d0)spotdparm(2,3)=1.0d0
              if(spotdparm(2,3).lt.0.0d0)spotdparm(2,3)=0.0d0
            endif
c
            if(kk.eq.1115)then    ! angular width of spot 1 on disk
              spotdparm(1,4)=var(i)
              if(spotdparm(1,4).lt.0.0d0)spotdparm(1,4)=0.0d0
            endif
c
            if(kk.eq.1119)then    ! angular width of spot 2 on disk
              spotdparm(2,4)=var(i)
              if(spotdparm(2,4).lt.0.0d0)spotdparm(2,4)=0.0d0
            endif
c
c
            if(kk.eq.1080)then         ! temperature factor spot 1, star 2
              spot2parm(1,1)=var(i)
            endif
c
            if(kk.eq.1084)then         ! temperature factor spot 2, star 2
              spot2parm(2,1)=var(i)
            endif
c
            if(kk.eq.1048)then         ! temperature factor spot 1, star 1
              spot1parm(1,1)=var(i)
            endif
c
            if(kk.eq.1052)then         ! temperature factor spot 2, star 1
              spot1parm(2,1)=var(i)
            endif
c
            if(kk.eq.1081)then         ! latitude spot 1, star 2
              spot2parm(1,2)=dmod(var(i),180.0d0)
            endif
c
            if(kk.eq.1085)then         ! latitude spot 2, star 2
              spot2parm(2,2)=dmod(var(i),180.0d0)
            endif
c
            if(kk.eq.1049)then         ! latitude spot 1, star 1
              spot1parm(1,2)=dmod(var(i),180.0d0)
            endif
c
            if(kk.eq.1053)then         ! latitude spot 2, star 1
              spot1parm(2,2)=dmod(var(i),180.0d0)
            endif
c
            if(kk.eq.1082)then         ! longitude spot 1, star 2
              spot2parm(1,3)=dmod(var(i),360.0d0)
            endif
c
            if(kk.eq.1086)then         ! longitude spot 2, star 2
              spot2parm(2,3)=dmod(var(i),360.0d0)
            endif
c
            if(kk.eq.1050)then         ! longitude spot 1, star 1
              spot1parm(1,3)=dmod(var(i),360.0d0)
            endif
c
            if(kk.eq.1054)then         ! longitude spot 2, star 1
              spot1parm(2,3)=dmod(var(i),360.0d0)
            endif
c
            if(kk.eq.1083)then         ! radius spot 1, star 2
              spot2parm(1,4)=var(i)
            endif
c
            if(kk.eq.1087)then         ! radius spot 2, star 2
              spot2parm(2,4)=var(i)
            endif
c
            if(kk.eq.1051)then         ! radius spot 1, star 1
              spot1parm(1,4)=var(i)
            endif
c
            if(kk.eq.1055)then         ! radius spot 2, star 1
              spot1parm(2,4)=var(i)
            endif

            if(kk.eq.498)then
              pshift=var(i)
              if(pshift.gt.1.0d0)pshift=1.0d0
              if(pshift.lt.-1.0d0)pshift=-1.0d0
            endif
c
            if(kk.eq.269)then
              finc=var(i)
c              if(finc.gt.90.0)finc=90.0
              if(finc.lt.0.0)finc=0.0
            endif
c
            if(kk.eq.384)then
              Q=var(i)
              if(Q.le.0.0)Q=0.001
            endif
c
            if(kk.eq.130)then
              ecc=var(i)
              if(ecc.le.0.0d0)ecc=0.000d0
              if(ecc.ge.1.0d0)ecc=0.9999d0
            endif
c
            if(kk.eq.17)then
              argper=var(i)
c              if(argper.le.0.0d0)argper=0.000d0
c              if(argper.gt.360.0d0)argper=360.0d0
            endif
c
            if(kk.eq.1176)then
              fill1=var(i)
              if(fill1.gt.1.0)fill1=1.0
            endif
c
            if(kk.eq.552)then
              rinner=var(i)
              if((rinner.lt.fill2).and.(Teff2.gt.0.0))rinner=fill2
            endif
c
            if(kk.eq.1177)then
              fill2=var(i)
              if(fill2.gt.1.0)fill2=1.0
              if(teff2.gt.0.0)rinner=fill2
              if(fill2.lt.0.0)then
                fill2=0.000001
                rinner=fill2
              endif
            endif
c
            if(kk.eq.1464)then
              omega1=var(i)
            endif
c
            if(kk.eq.1465)then
              omega2=var(i)
            endif
c
            if(kk.eq.558)then
              router=var(i)
              if(router.gt.1.0)router=1.0
            endif
c
            if(kk.eq.611)then
              Tdisk=var(i)
              if(Tdisk.lt.100.0)Tdisk=100.
            endif
c
            if(kk.eq.36)then
              betarim=var(i)
              if(betarim.lt.0.0)betarim=0.0
            endif
c
            if(kk.eq.1624)then
              Teff1=var(i)
              if(Teff1.lt.100.0)Teff1=100.
            endif
c
            if(kk.eq.1625)then
              Teff2=var(i)
              if(Teff2.lt.100.0)Teff2=100.
            endif
c
            if(kk.eq.744)then
              xi=var(i)
            endif
c
            if(kk.eq.375)then
              rLx=var(i)
              if(rLx.lt.0.0)rLx=0.0
            endif
c
            if(kk.eq.580)then
              separ=var(i)
              if(separ.lt.0.01)separ=0.01
            endif
c
            if(kk.eq.192)then
              gamma=var(i)
            endif
c
            if(kk.eq.1626)then
              t3=var(i)
              if(t3.lt.0.01)t3=0.01
            endif
c
            if(kk.eq.1210)then
              g3=var(i)
              if(g3.lt.0.01)t3=0.01
            endif
c
            if(kk.eq.576)then
              SA3=var(i)
c
c   UPDATE January 16, 2001
c
c   comment out this if-then statement
c
c              if(SA3.lt.0.01)t3=0.01
            endif
c
c   UPDATE November 6, 2002
c
c   Here are the assignments for the limb darkening parameters.
c   use strings x1, x2, to x8 for the x-coefficient for star 1
c   and y1, y2, to y8 for the y-coefficient for star 1.
c
c   Use z1, z2, to z8 for the x-coefficient for star 2 and
c   use w1, w2, to w8 for the y-coefficient for star 2.
c
            if(kk.eq.1752)dwavex(1,1)=var(i)
            if(kk.eq.1753)dwavex(2,1)=var(i)
            if(kk.eq.1754)dwavex(3,1)=var(i)
            if(kk.eq.1755)dwavex(4,1)=var(i)
            if(kk.eq.1756)dwavex(5,1)=var(i)
            if(kk.eq.1757)dwavex(6,1)=var(i)
            if(kk.eq.1758)dwavex(7,1)=var(i)
            if(kk.eq.1759)dwavex(8,1)=var(i)
c
            if(kk.eq.1784)dwavey(1,1)=var(i)
            if(kk.eq.1785)dwavey(2,1)=var(i)
            if(kk.eq.1786)dwavey(3,1)=var(i)
            if(kk.eq.1787)dwavey(4,1)=var(i)
            if(kk.eq.1788)dwavey(5,1)=var(i)
            if(kk.eq.1789)dwavey(6,1)=var(i)
            if(kk.eq.1790)dwavey(7,1)=var(i)
            if(kk.eq.1791)dwavey(8,1)=var(i)
c
            if(kk.eq.1816)dwavex(1,2)=var(i)
            if(kk.eq.1817)dwavex(2,2)=var(i)
            if(kk.eq.1818)dwavex(3,2)=var(i)
            if(kk.eq.1819)dwavex(4,2)=var(i)
            if(kk.eq.1820)dwavex(5,2)=var(i)
            if(kk.eq.1821)dwavex(6,2)=var(i)
            if(kk.eq.1822)dwavex(7,2)=var(i)
            if(kk.eq.1823)dwavex(8,2)=var(i)
c
            if(kk.eq.1720)dwavey(1,2)=var(i)
            if(kk.eq.1721)dwavey(2,2)=var(i)
            if(kk.eq.1722)dwavey(3,2)=var(i)
            if(kk.eq.1723)dwavey(4,2)=var(i)
            if(kk.eq.1724)dwavey(5,2)=var(i)
            if(kk.eq.1725)dwavey(6,2)=var(i)
            if(kk.eq.1726)dwavey(7,2)=var(i)
            if(kk.eq.1727)dwavey(8,2)=var(i)
c
c   body 3 limb darkening
c
            if(kk.eq.1400)dwavex(1,3)=var(i)  !tag m1
            if(kk.eq.1401)dwavex(2,3)=var(i)
            if(kk.eq.1402)dwavex(3,3)=var(i)
            if(kk.eq.1403)dwavex(4,3)=var(i)
            if(kk.eq.1404)dwavex(5,3)=var(i)
            if(kk.eq.1405)dwavex(6,3)=var(i)
            if(kk.eq.1406)dwavex(7,3)=var(i)
            if(kk.eq.1407)dwavex(8,3)=var(i)   ! tag m8
c
            if(kk.eq.1432)dwavey(1,3)=var(i)   !tag n1
            if(kk.eq.1433)dwavey(2,3)=var(i)
            if(kk.eq.1434)dwavey(3,3)=var(i)
            if(kk.eq.1435)dwavey(4,3)=var(i)
            if(kk.eq.1436)dwavey(5,3)=var(i)
            if(kk.eq.1437)dwavey(6,3)=var(i)
            if(kk.eq.1438)dwavey(7,3)=var(i)
            if(kk.eq.1439)dwavey(8,3)=var(i)   !tag n8
c

 10       continue
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&
c
          subroutine checklcfit(islc,Nmodel,xmodel,yinput,Ndata,xdata,
     @        ydata,errbar,chisq,zero,ifixgamma,itime,res)
c
c   November 12, 1999
c
c   This routine will convert the linear model to magnitudes and adjust
c   the zeropoint to give the best chi^2.  The best chi^2 value is
c   returned, as well as the zeropoint.  Set islc=1 for light curves
c   (linear models are converted to magnitudes) and islc=0 for RV curves.
c
c
c   January 24, 2001
c
c   Added the parameter ifixgamma.  If this is 1 or more and we are
c   fitting a velocity curve then the zero point is forced to be the input
c   gamma.
c
c

          implicit double precision (a-h,o-z)
c
          dimension xmodel(Nmodel),yinput(Nmodel),xdata(Ndata),
     &      ydata(Ndata),
     %      errbar(Ndata),ymodel(2000000),res(Ndata)
          dimension yinter(2000000),y2(2000000),
     &      ydummy(2000000)
          dimension xpad(2000000),ypad(2000000)
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
c   First, convert to magnitudes.
c

           if(itime.le.1)then
             call addpad(Nmodel,xmodel,yinput,xpad,ypad)
             Mmodel=Nmodel*3
           else
             do jj=1,Nmodel
               xpad(jj)=xmodel(jj)
               ypad(jj)=yinput(jj)
             enddo
             Mmodel=Nmodel
           endif
c 
            do 10 i=1,Mmodel
              if(islc.gt.0)then
c
c  Check to see that the argument of log10 is positive.
c
                if(ypad(i).gt.0.0d0)then
                  ymodel(i)=-2.5d0*dlog10(ypad(i))
                else
                  ymodel(i)=-99.9d0
                endif
              else
                ymodel(i)=ypad(i)
              endif
 10         continue
c
c   Find the maximum and minimum y-values of the data
c
          ymin=1000.
          ymax=-1000.
          do 20 i=1,Ndata
            if(ydata(i).lt.ymin)ymin=ydata(i)
            if(ydata(i).gt.ymax)ymax=ydata(i)
 20       continue
c

c
c   Interpolate the model so that we have y-values at all observed phases.
c
          call spline(xpad,ymodel,Mmodel,0.0d0,0.0d0,y2)
c
          do 30 i=1,Ndata
            call splint(xpad,ymodel,y2,Mmodel,xdata(i),qqq)
            yinter(i)=qqq
 30       continue
c
c   Now find the optimal zero point that will give the lowest chi^2.
c
          call getmean(Ndata,ydata,dataave)
          call getmean(Ndata,yinter,rmodelave)
          savezero=zero
          zero=(dataave-rmodelave)
          step=abs(rmodelave-(dataave))/1000.0d0
c
          chi1=0.0d0
          chi2=0.0d0
          chi3=0.0d0
          small=1.0d35
          zerosmall=1.0d30
c
          do 750 i=1,40
            zero1=zero
            call offset(Ndata,yinter,ydummy,zero1)
            call getchi(Ndata,ydata,errbar,ydummy,chi1)
            if(chi1.lt.small)then
              small=chi1
              zerosmall=zero1
            endif
            fn=0.0d0
            zero2=zero+step
            call offset(Ndata,yinter,ydummy,zero2)
            call getchi(Ndata,ydata,errbar,ydummy,chi2)
            if(chi2.lt.small)then
              small=chi2
              zerosmall=zero2
            endif
            diff=chi1-chi2
            if(diff.gt.0.0d0)go to 5061
            step=-step
c  
            csave=chi1
            chi1=chi2
            chi2=csave
            zsave=zero1
            zero1=zero2
            zero2=zsave
c
 5061       fn=fn+1.0d0
c
            zero3=zero2+step
            call offset(Ndata,yinter,ydummy,zero3)
            call getchi(Ndata,ydata,errbar,ydummy,chi3)
            if(chi3.lt.small)then
              small=chi3
              zerosmall=zero3
            endif
            diff23=chi3-chi2
            if(diff23.lt.0.0d0)then
              chi1=chi2
              chi2=chi3
              zero1=zero2
              zero2=zero3
              zero=zero3
              go to 5061
            endif          
c
c         find the minimum of parabola defined by the last three points
c
            if((chi3-chi2).eq.0.0d0)go to 999
            step=step*(1.0d0/(1.0d0+(chi1-chi2)/(chi3-chi2))+0.5d0)
            zero=zero2-step
            step=step*fn/3.0d0
            call offset(Ndata,yinter,ydummy,zero)
            call getchi(Ndata,ydata,errbar,ydummy,chi4)
            if(chi4.lt.small)then
              small=chi4
              zerosmall=zero
            endif
c
 750      continue  ! loop over grid searches
c
          continue                  ! come here when delta chi is small
c
 999      zero=zerosmall
c
          if((islc.le.0).and.(ifixgamma.eq.0))zero=savezero
          call offset(Ndata,yinter,ydummy,zero)
          call getchi(Ndata,ydata,errbar,ydummy,chi3)
c
c          write(*,*)zero
          chisq=chi3
c
c   UPDATE October 31, 2011
c          
c   compute the residuals
c
          do ii=1,Ndata
            res(ii)=ydata(ii)-ydummy(ii)
          enddo
c
          return
c
          end
c
c  -----------------------------------------------------
c
          subroutine getchi(N,y,err,ydum,chisq)
c
          implicit double precision (a-h,o-z)
c
          dimension y(N),ydum(N),err(N)
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 

          chisq=0.0d0
c
c   If the variable rmed > 1, then define chisq to be the
c   absolute deviation
c
          do 10 i=1,N
            if(rmed.ge.1.0d0)then
              chisq=chisq+dabs((y(i)-ydum(i))/err(i))
            else
              chisq=chisq+(y(i)-ydum(i))*(y(i)-ydum(i))/(err(i)*err(i))
            endif
 10       continue
          return
          end
c
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine offset(N,y,yoff,off)
c
c   Will offset the y values by offset and return new array yoff.
c
          implicit double precision (a-h,o-z)
c
          dimension y(N),yoff(N)
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
          subroutine getmean(N,y,average)
c
          implicit double precision (a-h,o-z)
c
          dimension y(N)
c
          summ=0.0
          average=0.0
c
          do 10 i=1,N
            summ=summ+y(i)
 10       continue
c 
          average=summ/float(N)
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine initchi(chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,
     &            chisqH,chisqK,chisqRV1,chisqRV2,chisqRV3)
c
          implicit double precision (a-h,o-z)
c
          chisqU=0.d0
          chisqB=0.d0
          chisqV=0.d0
          chisqR=0.d0
          chisqI=0.d0
          chisqJ=0.d0
          chisqH=0.d0
          chisqK=0.d0
          chisqRV1=0.d0
          chisqRV2=0.d0
          chisqRV3=0.d0
c
          return
          end
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %       icnRV1,icnRV2,chiU,chiB,chiV,chiR,chiI,chiJ,chiH,chiK,
     &       chiRV1,chiRV2,icnRV3,chiRV3)
c
c   November 15, 1999
c
c   This routine will write chi^2 values to the screen in a compact way.
c
          implicit double precision (a-h,o-z)
c
          character*80 line
          character*30 outstring
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
c   UPDATE May 27, 2002.
c
c   If the variable rmed > 1, then we are doing median fitting.
c   In that case, use different format statements.
c
          iline=0
          iwrite=0
          
c
c          do 10 i=1,10
c            section(i)='1234567890123456'
c            section(i)='                '
c 10       continue
          line='  -->'
          ilength=5
c
          if(icnU.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medU',chiU,outstring,lll)
            else
              call chistring(' chiU',chiU,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
c
          if(icnB.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medB',chiB,outstring,lll)
            else
              call chistring(' chiB',chiB,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
c
          if(icnV.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medV',chiV,outstring,lll)
            else
              call chistring(' chiV',chiV,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
C
          if(icnR.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medR',chiR,outstring,lll)
            else
              call chistring(' chiR',chiR,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
C
          if(icnI.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medI',chiI,outstring,lll)
            else
              call chistring(' chiI',chiI,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
C
          if(icnJ.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medJ',chiJ,outstring,lll)
            else
              call chistring(' chiJ',chiJ,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
C
          if(icnH.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medH',chiH,outstring,lll)
            else
              call chistring(' chiH',chiH,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
C
          if(icnK.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medK',chiK,outstring,lll)
            else
              call chistring(' chiK',chiK,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
C
          if(icnRV1.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medRV1',chiRV1,outstring,lll)
            else
              call chistring(' chiRV1',chiRV1,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
C
          if(icnRV2.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medRV2',chiRV2,outstring,lll)
            else
              call chistring(' chiRV2',chiRV2,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
c
          if(icnRV3.ne.430)then
            iline=iline+1
            iwrite=0
            if(rmed.ge.1.0d0)then
              call chistring(' medRV3',chiRV3,outstring,lll)
            else
              call chistring(' chiRV3',chiRV3,outstring,lll)
            endif
            K=lnblnk(line)
            line=line(1:K)//outstring
            ilength=ilength+lll
            if(ilength.gt.51)then
              iwrite=100
              write(*,500)line
              ilength=5
              line='  -->'
            endif
          endif
c
          if(ilength.gt.5.and.iwrite.eq.0)write(*,500)line
c
 500      format(a80)

          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  RVG BUG ALERT  May 9, 2001
c
c  Add the spot parameters to the argument list.
c
c  NEW BUG ALERT August 2, 2001
c
c  Replace 'sw4' with T0
c
c  UPDATE August 10, 2004
c
c  Add 8 real and 4 integer variables to the argument list.
c
c  May 8, 2006
c
c  Add isw21-isw24, sw21-sw24, powercoeff to list.
c
c  UPDATE November 6, 2008
c
c  Add sw25-sw34 and isw25-isw34 below
c
          subroutine writegridout(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,
     &       omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,
     @       betarim,rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,
     @       Nref,rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,
     @       sw3,T0,idraw,iecheck,idint,iatm,ism1,icnU,icnB,icnV,icnR,
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
c
c   November 15, 1999
c
c   This routine will write a file called 'gridELC.inp' which will be
c   of the same format as ELC.inp.  The parameters written will be those
c   found by gridELC.
c
          implicit double precision (a-h,o-z)
c
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,3),
     @        dwavey(8,3)
          dimension powercoeff(8,9)
c
c   RVG BUG ALERT  May 9, 2001
c
c   Dimension the spot arrays.
c
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
c 
          open(unit=1,file='gridELC.inp',status='unknown')
c
c   Set the icn? control numbers back to zeros and 1s.
c
          if(icnU.eq.430)then
            icnU=0
          else
            icnU=1
          endif
          if(icnB.eq.430)then
            icnB=0
          else
            icnB=1
          endif
          if(icnV.eq.430)then
            icnV=0
          else
            icnV=1
          endif
          if(icnR.eq.430)then
            icnR=0
          else
            icnR=1
          endif
          if(icnI.eq.430)then
            icnI=0
          else
            icnI=1
          endif
          if(icnJ.eq.430)then
            icnJ=0
          else
            icnJ=1
          endif
          if(icnH.eq.430)then
            icnH=0
          else
            icnH=1
          endif
          if(icnK.eq.430)then
            icnK=0
          else
            icnK=1
          endif
c
          write(1,1000)Nalph1
          write(1,1001)Nbet1
          write(1,1002)Nalph2
          write(1,1003)Nbet2
          write(1,1004)fill1
          write(1,1005)fill2
          write(1,1006)omega1
          write(1,1007)omega2
          write(1,1008)dphase
          write(1,1009)Q
          write(1,1010)finc
          write(1,1011)Teff1
          write(1,1012)Teff2
          write(1,1013)Tgrav1
          write(1,1014)Tgrav2
          write(1,1015)betarim
          write(1,1016)rinner
          write(1,1017)router
          write(1,1018)tdisk
          write(1,2018)xi
          write(1,1019)Ntheta
          write(1,1020)Nradius
          write(1,1042)alb1
          write(1,1043)alb2
          write(1,2043)Nref
          write(1,1021)rLx
          write(1,1023)Period
          if(fm.gt.1.0d-4)then
            write(1,1024)fm
          else
            write(1,9024)fm
          endif
 9024     format(1pe16.9,9x,'fm')
          write(1,1025)separ
          write(1,4025)gamma
          write(1,5000)t3
          write(1,5001)g3
          write(1,5002)SA3
          write(1,5003)density
          write(1,5004)sw1
          write(1,5005)sw2
          write(1,5006)sw3
          write(1,5007)T0
          write(1,1026)idraw
          write(1,1027)iecheck
          write(1,1028)idint
          write(1,4000)iatm
          write(1,4001)ism1
          write(1,4002)icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK
          write(1,5008)iRVfilt
          write(1,5009)isw1
          write(1,5010)isw2
          write(1,5011)isw3
          write(1,5012)isw4
          write(1,3028)ilaw
c
          do 10 i=1,8
            write(1,2000)wave(i),dbolx(i,1),dboly(i,1),dbolx(i,2),
     @         dboly(i,2),dwavex(i,1),dwavey(i,1),dwavex(i,2),
     @         dwavey(i,2)
 10       continue
c
          write(1,6000)ecc
          write(1,6001)argper
          write(1,6002)pshift
          write(1,6003)sw5
          write(1,6004)sw6
          write(1,6005)sw7
          write(1,6006)sw8
          write(1,6007)sw9
          write(1,6008)ikeep
          write(1,6009)isynch
          write(1,6010)isw5
          write(1,6011)isw6
          write(1,6012)isw7
          write(1,6013)isw8
          write(1,6014)isw9
c
c   RVG BUG ALERT  May 9, 2001
c
c   Add the spot parameters here.
c
          write(1,2001)spot1parm(1,1)
          write(1,2002)spot1parm(1,2)
          write(1,2003)spot1parm(1,3)
          write(1,2004)spot1parm(1,4)
c                     
          write(1,2005)spot1parm(2,1)
          write(1,2006)spot1parm(2,2)
          write(1,2007)spot1parm(2,3)
          write(1,2008)spot1parm(2,4)
c                     
c                     
          write(1,2009)spot2parm(1,1)
          write(1,2010)spot2parm(1,2)
          write(1,2011)spot2parm(1,3)
          write(1,2012)spot2parm(1,4)
                      
          write(1,2013)spot2parm(2,1)
          write(1,2014)spot2parm(2,2)
          write(1,2015)spot2parm(2,3)
          write(1,2016)spot2parm(2,4)
c                     
c                     
          write(1,2017)spotdparm(1,1)
          write(1,22018)spotdparm(1,2)
          write(1,2019)spotdparm(1,3)
          write(1,2020)spotdparm(1,4)
                      
          write(1,2021)spotdparm(2,1)
          write(1,2022)spotdparm(2,2)
          write(1,2023)spotdparm(2,3)
          write(1,2024)spotdparm(2,4)
c
c   UPDATE August 10, 2004
c
c   Write out the 8 real and 4 new integer variables here.
c
          write(1,2025)primmass
          write(1,2026)primK
          write(1,2027)primrad
          write(1,2028)ratrad
c
          write(1,2030)frac1
          write(1,2031)frac2
          write(1,2032)ecosw
          write(1,2033)temprat
c
          write(1,2040)idark1
          write(1,2041)idark2
          write(1,2042)isw12
          write(1,3043)isw13
c
          write(1,8001)isw21
          write(1,8002)isw22
          write(1,8003)isw23
          write(1,8004)isw24
c
8001      format(i1,24x,'ialign (0 for rotation aligned with orbit)')
8002      format(i1,24x,'ifastgen (1 for fast genetic mode)')
8003      format(i1,24x,'iwriteeclipse (1 to fit eclipse times)')         
8004      format(i1,24x,
     @          'frac switch (>1 to enable ELCratio.???? files)')
c
          do 85000 kk=1,8
              write(1,85001)powercoeff(kk,1),powercoeff(kk,2),
     $       powercoeff(kk,3),powercoeff(kk,4),powercoeff(kk,5),
     $       powercoeff(kk,6),powercoeff(kk,7),powercoeff(kk,8),
     #       powercoeff(kk,9)
85000     continue
c
85001      format(9(f7.4,1x))
c
          write(1,8011)bigI
          write(1,8012)bigbeta
          write(1,8013)sw23
          write(1,8014)sw24
c
8011      format(f14.10,11x,'axis_I (inclination of rotation axis if',
     @       ' ialign=1)')
8012      format(f14.10,11x,'axis_beta (angle of rotation axis wrt to ',
     @       ' orbit if ialign=1)')          
8013      format(f18.7,7x,'t_start')
8014      format(f18.7,7x,'t_end')         
c
c  UPDATE November 6, 2008
c
c  write the new variables sw25-sw34 and isw25-isw34 here
c
          write(1,8025)sw25
          write(1,8026)sw26
          write(1,8027)sw27
          write(1,8028)sw28
          write(1,8029)sw29
          write(1,8030)sw30
          write(1,8031)contam
          write(1,8032)Tconj
          write(1,8033)beam1
          write(1,8034)beam2

          write(1,9025)isw25
          write(1,9026)isw26
          write(1,9027)isw27
          write(1,9028)isw28
          write(1,9029)isw29
          write(1,9030)isw30
          write(1,9031)isw31
          write(1,9032)isw32
          write(1,9033)isw33
          write(1,9034)isw34
c
          write(1,8040)ocose
          write(1,8041)osine
          write(1,8042)omegadot
          write(1,8043)contamS0
          write(1,8044)contamS1
          write(1,8045)contamS2
          write(1,8046)contamS3
          write(1,8047)sw47
          write(1,8048)sw48
          write(1,8049)sw49

8025      format(f13.10,12x,'asini error')
8026      format(f13.11,12x,'reference phase for disk fraction')
8027      format(f13.11,12x,
     @       'radfill1 (set to use fill1 in terms of R_eff')
8028      format(f13.11,12x,
     @       'radfill2 (set to use fill2 in terms of R_eff')
8029      format(f13.9,12x,'bin size for light curves (minutes)')
8030      format(f13.9,12x,'bin size for RV curves (minutes)')
8031      format(f16.14,9x,'Kepler contamination')
8032      format(f21.14,4x,'Tconj')
8033      format(f13.11,12x,'beam1 (Doppler boosting factor, star 1)')
8034      format(f13.11,12x,'beam2 (Doppler boosting factor, star 2)')


 8040     format(f20.16,5x,'e*cos(omega)')
 8041     format(f20.16,5x,'e*sin(omega)')
 8042     format(f17.13,8x,'omega_dot (degrees per year)')
 8043     format(f17.13,8x,'contamS0 (season 0 contamination, tag s0)')
 8044     format(f17.13,8x,'contamS1 (season 1 contamination, tag s1)')
 8045     format(f17.13,8x,'contamS2 (season 2 contamination, tag s2)')
 8046     format(f17.13,8x,'contamS3 (season 3 contamination, tag s3)')
 8047     format(f19.8,6x,'Tref for dynamical integrator')
 8048     format(f12.8,13x,'threshold to write chi^2')
 8049     format(f12.8,13x,'sw49 (currently inactive)')

9025      format(i1,24x,'X-ray foreshortening switch',
     @       ' (1 for point source)')
9026      format(i1,24x,'iGR (1 for GR, 2 for tidal, 3 for both)')
9027      format(i6,19x,'Nterms for fast analytic')
9028      format(i1,24x,'set to 1 to fit for Tconj')
9029      format(i1,24x,
     $      'set to 1 to fit for e*cos(omega), e*sin(omega)')
9030      format(i1,24x,'body 3 switch')
9031      format(i6,19x,'Ngap')
9032      format(i16,9x,'jdum (seed for markovELC, geneticELC, ',
     @        'randomELC)')
9033      format(i1,24x,'mandel (0 for Gimenez, 1 for Mandel & Agol)')
9034      format(i1,24x,'Iseason (1 for seasonal Kepler contamination)')

c
 2025     format(f20.16,5x,'primmass (star 1 mass in solar masses)')
 2026     format(f19.14,6x,'primK (K-velocity of star 1 in km/sec)')
 2027     format(f19.14,6x,'primrad (star 1 radius in solar radii)')
 2028     format(f21.14,4x,
     &          'ratrad (ratio of star 1 radius and star 2 radius)')
c
 2030     format(f18.16,7x,'frac1 (fractional radius star 1: R_1/a)')
 2031     format(f18.16,7x,'frac2 (fractional radius star 2: R_2/a)')
 2032     format(f17.14,8x,'ecosw (phase difference between eclipses)')
 2033     format(f17.14,8x,'temprat (T_2/T_1)')
c
 2040     format(i1,24x,'idark1')
 2041     format(i1,24x,'idark2')
 2042     format(i6,19x,'Npoly (0 for numerical integration)')
 3043     format(i1,24x,'ifasttrans (>0 for fast transit mode)')

          close(1)
c
 1000     format(i6,19x,'Nalph1')
 1001     format(i5,20x,'Nbet1')
 1002     format(i6,19x,'Nalph2')
 1003     format(i5,20x,'Nbet2')
 1004     format(f16.14,9x,'fill1')
 1005     format(f16.14,9x,'fill2')
 1006     format(f16.12,9x,'omega1')
 1007     format(f16.12,9x,'omega2')
 1008     format(f13.8,12x,'dphase')
 1009     format(f21.16,4x,'Q')
 1010     format(f21.16,4x,'finc')
 1011     format(f12.5,13x,'Teff1')
 1012     format(f12.5,13x,'Teff2')
 1013     format(f9.7,16x,'Tgrav1')
 1014     format(f9.7,16x,'Tgrav2')
 1015     format(f11.8,14x,'betarim')
 1016     format(f11.9,14x,'rinner')
 1017     format(f11.9,14x,'router')
 1018     format(f12.6,13x,'tdisk')
 2018     format(f12.9,13x,'xi')
 1019     format(i6,19x,'Ntheta')
 1020     format(i6,19x,'Nradius')
 1021     format(f13.8,12x,'Lx/Lopt')
 1023     format(f21.15,4x,'Period')
 1024     format(f8.5,17x,'fm')
 1025     format(f20.11,5x,'separ')
 1026     format(i1,24x,'idraw')
 1027     format(i1,24x,'iecheck')
 1028     format(i1,24x,'idint')
 1042     format(f12.10,13x,'alb1')
 1043     format(f12.10,13x,'alb2')
c
c   UPDATE April 15, 2002
c
c   Change the format statement below to i2,18x
c
 2043     format(i2,23x,'Nref')
 2000     format(f7.1,2x,4(f8.5,1x),4(f11.8,1x))
 3028     format(i2,23x,'ilaw  (1=lin. law, 2=log. law,',
     %           ' 3=sqrt law, 4=quad. law)')
 4000     format(i1,24x,'iatm')
 4001     format(i1,24x,'ism1')
 4002     format(8(i1,1x),9x,'icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK')
 4025     format(f18.12,7x,'gamma velocity')
 5000     format(f18.10,7x,'t3')
 5001     format(f15.11,10x,'g3')
 5002     format(f22.16,3x,'SA3')
 5003     format(f12.6,13x,'density in g/cc')
 5004     format(f12.6,13x,'onephase')
 5005     format(f12.6,13x,'usepot1')
 5006     format(f12.6,13x,'usepot2')
 5007     format(f21.14,4x,'T0')
 5008     format(i1,24x,'iRVfilt')
 5009     format(i1,24x,'ionephase')
 5010     format(i1,24x,'isquare')
 5011     format(i1,24x,'iusepot')
 5012     format(i1,24x,'ifixgamma')
 6000     format(f18.12,7x,'eccentricity')
 6001     format(f19.14,6x,'argument of peristron in degrees')
 6002     format(f14.11,11x,'pshift')
c
c   UPDATE April 8, 2002
c
c   Change the format statement below:
c
 6003     format(f17.11,8x,'asini (projected semimajor axis in seconds)')
c
c   UPDATE May 27, 2002
c
c   Change this format statement
c
 6004     format(f12.6,13x,'median fit (geneticELC only)')
 6005     format(f12.6,13x,'sw7 (currently inactive)')
 6006     format(f12.6,13x,'sw8 (currently inactive)')
 6007     format(f15.9,10x,'time step for itime=2')
 6008     format(i1,24x,'ikeep (1 to put eclipse at phase 0.0)')
 6009     format(i1,24x,'isynch (1 to keep rotation synchronous',
     $         ' at periastron)')
 6010     format(i1,24x,'ispotprof')
 6011     format(i1,24x,'igrav')
 6012     format(i1,24x,'itime')
c 
c  UPDATE JULY 4, 2004
c
c  The variable isw8 will be assigned to MonteCarlo, which
c  will be used in the Monte Carlo routine to determine
c  fractionally eclipsed pixels.
c
 6013     format(i6,19x,'MonteCarlo (0 for interpolation, >10 for ',
     $     'Monte Carlo integration)')
 6014     format(i6,19x,'ielite')
c
c
c  RVG BUG ALERT  May 9, 2001
c
c  Add these format statements.
c
 2001     format(f13.10, 12x,'Temperature factor spot 1, star 1')
 2002     format(f13.9, 12x,'Latitude of spot 1, star 1 (degrees)') 
 2003     format(f13.9, 12x,'Longitude of spot 1, star 1 (degrees)') 
 2004     format(f13.9, 12x,
     @        'Angular radius of spot 1, star 1 (degrees)') 

 2005     format(f13.10, 12x,'Temperature factor spot 2, star 1')
 2006     format(f13.9, 12x,'Latitude of spot 2, star 1 (degrees)') 
 2007     format(f13.9, 12x,'Longitude of spot 2, star 1 (degrees)') 
 2008     format(f13.9, 12x,
     @        'Angular radius of spot 2, star 1 (degrees)') 

 2009     format(f13.10, 12x,'Temperature factor spot 1, star 2')
 2010     format(f13.9, 12x,'Latitude of spot 1, star 2 (degrees)') 
 2011     format(f13.9, 12x,'Longitude of spot 1, star 2 (degrees)') 
 2012     format(f13.9, 12x,
     @        'Angular radius of spot 1, star 2 (degrees)') 

 2013     format(f13.10, 12x,'Temperature factor spot 2, star 2')
 2014     format(f13.9, 12x,'Latitude of spot 2, star 2 (degrees)') 
 2015     format(f13.9, 12x,'Longitude of spot 2, star 2 (degrees)') 
 2016     format(f13.9, 12x,
     @        'Angular radius of spot 2, star 2 (degrees)') 

 2017     format(f13.10, 12x,'Temperature factor spot 1, disk')
22018     format(f13.9, 12x,'Azimuth of spot 1, disk (degrees)') 
 2019     format(f13.9, 12x,'Radial cutoff of spot 1, disk (0 <= ',
     @            ' r_cut <=1)') 
 2020     format(f13.9, 12x,'Angular size of spot 1, disk (degrees)') 

 2021     format(f13.10, 12x,'Temperature factor spot 2, disk')
 2022     format(f13.9, 12x,'Azimuth of spot 2, disk (degrees)') 
 2023     format(f13.9, 12x,'Radial cutoff of spot 2, disk (0 <= r_cut',
     @        ' <=1)') 
 2024     format(f13.9, 12x,'Angular size of spot 2, disk (degrees)') 
c
c  Reset the values of icnU, etc. in case they were 430
c   Set the icn? control numbers back to zeros and 1s.
c
          if(icnU.eq.0)then
            icnU=430
          endif
          if(icnB.eq.0)then
            icnB=430
          endif
          if(icnV.eq.0)then
            icnV=430
          endif
          if(icnR.eq.0)then
            icnR=430
          endif
          if(icnI.eq.0)then
            icnI=430
          endif
          if(icnJ.eq.0)then
            icnJ=430
          endif
          if(icnH.eq.0)then
            icnH=430
          endif
          if(icnK.eq.0)then
            icnK=430
          endif

c
          return
          end
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
 501      format(f25.12,1x,f18.10,1x,i7,2x,2(f10.7,1x),1x,a2,1x,i2)
 502      format(f23.9,3x,f15.10)
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  NEW BUG August 2, 2001
c
c  Add the period and T0 to the argument list.
c
c  UPDATE January 15, 2002
c
c  Add alb1 and alb2 to the argument list
c
c  UPDATE November 6, 2002
c
c  Add the limb darkening coefficients dwavex and dwavey to the
c  argument list of writevar.  Dimension them as dwavex(8,2)
c  and dwavey(8,2).
c
          subroutine writevar(Nvarmax,svar,var,fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,
     @       rLx,separ,isvel1,gamma1,isvel2,gamma2,t3,g3,SA3,ecc,argper,
     @       pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,alb2,
     @       dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,frac2,
     @       ecosw,temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,
     @       beam2,contam,ocose,osine,isw29,tertperiod,tertt0,tertecos,
     @       tertesin,tertincl,tertOmega,tertQ,Tgrav1,Tgrav2,tertconj,
     @       omegadot,contamS0,contamS1,contamS2,contamS3,
     @       P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @       P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @       P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @       P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @       P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @       P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @       P6ratrad,
     @       P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
     @       P7ratrad,
     @       P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @       P8ratrad,gamma3,isvel3,sw72,sw73)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c  UPDATE October 10, 2008
c
c  Add the density to the list above
c
c   November 15, 1999
c
c   This routine will determine which variables need to be printed on the
c   screen based on the string codes in svar(1:Nvarmax).  The output
c   printed on the screen will be compact.
c
          implicit double precision (a-h,o-z)
c
          character*80 line
          character*40 svar(Nvarmax)
c          character*16 string16
c          character*10 string10
c          character*11 string11
c          character*12 string12
c          character*13 string13
c          character*14 string14
c          character*15 string15
c
c   UPDATE August 2, 2004
c
c   Add string 17 for period and T0.
c
c          character*17 string17
c          character*18 string18
          character*19 string19
c          character*20 string20
c          character*21 string21
c          character*22 string22
c          character*23 string23
c          character*24 string24
          character*30 outstring

          dimension powercoeff(8,9)
          dimension var(Nvarmax)
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          dimension dwavex(8,3),dwavey(8,3)
c
          line='**'
          iline=0
          iwrite=0
          ilength=2
          iargper=0              !flag for ecosw being assigned.
          do 10 i=1,Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
            if((kk.eq.450).and.(isw29.gt.0))then          !ecos(omega), tag oc
              iline=iline+1
              iwrite=0
              call pstring(' ecos(omega)',8,ocose,outstring,lll)
              K=lnblnk(line)
              line=line(1:K)//outstring
              ilength=ilength+lll
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if((kk.eq.466).and.(isw29.gt.0))then  !esin(omega), tag os
              iline=iline+1
              iwrite=0
              call pstring(' esin(omega)',8,osine,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c  tidal apsidal constants
c
            if(kk.eq.1016)then  !rk1, tag a1
              iline=iline+1
              iwrite=0
              call pstring(' rk1',6,sw72,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1017)then  !rk2, tag a2
              iline=iline+1
              iwrite=0
              call pstring(' rk2',6,sw73,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   planet 2 parameters
c
            if(kk.eq.649)then  !P2tconj, tag uj
              iline=iline+1
              iwrite=0
              call pstring(' P2tconj',6,P2tconj,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.659)then  !P2period, tag ut
              iline=iline+1
              iwrite=0
              call pstring(' P2period',9,P2period,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.660)then  !P2T0, tag uu
              iline=iline+1
              iwrite=0
              call pstring(' P2T0',6,P2T0,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.661)then  !P2ecos, tag uv
              iline=iline+1
              iwrite=0
              call pstring(' P2ecos',7,P2ecos,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.662)then  !P2esin, tag uw
              iline=iline+1
              iwrite=0
              call pstring(' P2esin',7,P2esin,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.663)then  !P2incl, tag ux
              iline=iline+1
              iwrite=0
              call pstring(' P2incl',7,P2incl,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.664)then  !P2Omega, tag uy
              iline=iline+1
              iwrite=0
              call pstring(' P2Omega',6,P2Omega,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.665)then  !P2Q, tag uz
              iline=iline+1
              iwrite=0
              call pstring(' P2Q',6,P2Q,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.641)then  !P2artrad, tag ub
              iline=iline+1
              iwrite=0
              call pstring(' P2ratrad',7,P2ratrad,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   planet 3 parameters
c
            if(kk.eq.681)then  !P3tconj, tag vj
              iline=iline+1
              iwrite=0
              call pstring(' P3tconj',6,P3tconj,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.691)then  !P3period, tag vt
              iline=iline+1
              iwrite=0
              call pstring(' P3period',9,P3period,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.692)then  !P3T0, tag vu
              iline=iline+1
              iwrite=0
              call pstring(' P3T0',6,P3T0,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.693)then  !P3ecos, tag vv
              iline=iline+1
              iwrite=0
              call pstring(' P3ecos',7,P3ecos,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.694)then  !P3esin, tag vw
              iline=iline+1
              iwrite=0
              call pstring(' P3esin',7,P3esin,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.695)then  !P3incl, tag vx
              iline=iline+1
              iwrite=0
              call pstring(' P3incl',7,P3incl,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.696)then  !P3Omega, tag vy
              iline=iline+1
              iwrite=0
              call pstring(' P3Omega',6,P3Omega,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.697)then  !P3Q, tag vz
              iline=iline+1
              iwrite=0
              call pstring(' P3Q',6,P3Q,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.673)then  !P3artrad, tag vb
              iline=iline+1
              iwrite=0
              call pstring(' P3ratrad',7,P3ratrad,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   planet 4 parameters
c
            if(kk.eq.713)then  !P4tconj, tag wj
              iline=iline+1
              iwrite=0
              call pstring(' P4tconj',6,P4tconj,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.723)then  !P4period, tag wt
              iline=iline+1
              iwrite=0
              call pstring(' P4period',9,P4period,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.724)then  !P4T0, tag wu
              iline=iline+1
              iwrite=0
              call pstring(' P4T0',6,P4T0,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.725)then  !P4ecos, tag wv
              iline=iline+1
              iwrite=0
              call pstring(' P4ecos',7,P4ecos,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.726)then  !P4esin, tag ww
              iline=iline+1
              iwrite=0
              call pstring(' P4esin',7,P4esin,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.727)then  !P4incl, tag wx
              iline=iline+1
              iwrite=0
              call pstring(' P4incl',7,P4incl,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.728)then  !P4Omega, tag wy
              iline=iline+1
              iwrite=0
              call pstring(' P4Omega',6,P4Omega,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.729)then  !P4Q, tag wz
              iline=iline+1
              iwrite=0
              call pstring(' P4Q',6,P4Q,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.705)then  !P4ratrad, tag wb
              iline=iline+1
              iwrite=0
              call pstring(' P4ratrad',7,P4ratrad,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   planet 5 parameters
c
            if(kk.eq.745)then  !P5tconj, tag xj
              iline=iline+1
              iwrite=0
              call pstring(' P5tconj',6,P5tconj,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.755)then  !P5period, tag xt
              iline=iline+1
              iwrite=0
              call pstring(' P5period',9,P5period,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.756)then  !P5T0, tag xu
              iline=iline+1
              iwrite=0
              call pstring(' P5T0',6,P5T0,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.757)then  !P5ecos, tag xv
              iline=iline+1
              iwrite=0
              call pstring(' P5ecos',7,P5ecos,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.758)then  !P5esin, tag xw
              iline=iline+1
              iwrite=0
              call pstring(' P5esin',7,P5esin,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.759)then  !P5incl, tag xx
              iline=iline+1
              iwrite=0
              call pstring(' P5incl',7,P5incl,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.760)then  !P5Omega, tag xy
              iline=iline+1
              iwrite=0
              call pstring(' P5Omega',6,P5Omega,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.761)then  !P5Q, tag xz
              iline=iline+1
              iwrite=0
              call pstring(' P5Q',6,P5Q,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.737)then  !P5ratrad, tag xb
              iline=iline+1
              iwrite=0
              call pstring(' P5ratrad',7,P5ratrad,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   planet 6 parameters
c
            if(kk.eq.585)then  !P6tconj, tag sj
              iline=iline+1
              iwrite=0
              call pstring(' P6tconj',6,P6tconj,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.595)then  !P6period, tag st
              iline=iline+1
              iwrite=0
              call pstring(' P6period',9,P6period,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.596)then  !P6T0, tag su
              iline=iline+1
              iwrite=0
              call pstring(' P6T0',6,P6T0,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.597)then  !P6ecos, tag sv
              iline=iline+1
              iwrite=0
              call pstring(' P6ecos',7,P6ecos,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.598)then  !P6esin, tag sw
              iline=iline+1
              iwrite=0
              call pstring(' P6esin',7,P6esin,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.599)then  !P6incl, tag sx
              iline=iline+1
              iwrite=0
              call pstring(' P6incl',7,P6incl,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.600)then  !P6Omega, tag sy
              iline=iline+1
              iwrite=0
              call pstring(' P6Omega',6,P6Omega,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.601)then  !P6Q, tag sz
              iline=iline+1
              iwrite=0
              call pstring(' P6Q',6,P6Q,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.577)then  !P6ratrad, tag sb
              iline=iline+1
              iwrite=0
              call pstring(' P6ratrad',7,P6ratrad,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   planet 7 parameters
c
            if(kk.eq.233)then  !P7tconj, tag hj
              iline=iline+1
              iwrite=0
              call pstring(' P7tconj',6,P7tconj,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.243)then  !P7period, tag ht
              iline=iline+1
              iwrite=0
              call pstring(' P7period',9,P7period,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.244)then  !P7T0, tag hu
              iline=iline+1
              iwrite=0
              call pstring(' P7T0',6,P7T0,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.245)then  !P7ecos, tag hv
              iline=iline+1
              iwrite=0
              call pstring(' P7ecos',7,P7ecos,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.246)then  !P7esin, tag hw
              iline=iline+1
              iwrite=0
              call pstring(' P7esin',7,P7esin,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.247)then  !P7incl, tag hx
              iline=iline+1
              iwrite=0
              call pstring(' P7incl',7,P7incl,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.248)then  !P7Omega, tag hy
              iline=iline+1
              iwrite=0
              call pstring(' P7Omega',6,P7Omega,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.249)then  !P7Q, tag hz
              iline=iline+1
              iwrite=0
              call pstring(' P7Q',6,P7Q,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.225)then  !P7ratrad, tag hb
              iline=iline+1
              iwrite=0
              call pstring(' P7ratrad',7,P7ratrad,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   planet 8 parameters
c
            if(kk.eq.329)then  !P8tconj, tag kj
              iline=iline+1
              iwrite=0
              call pstring(' P8tconj',6,P8tconj,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.339)then  !P8period, tag kt
              iline=iline+1
              iwrite=0
              call pstring(' P8period',9,P8period,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.340)then  !P8T0, tag ku
              iline=iline+1
              iwrite=0
              call pstring(' P8T0',6,P8T0,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.341)then  !P8ecos, tag kv
              iline=iline+1
              iwrite=0
              call pstring(' P8ecos',7,P8ecos,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.342)then  !P8esin, tag kw
              iline=iline+1
              iwrite=0
              call pstring(' P8esin',7,P8esin,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.343)then  !P8incl, tag kx
              iline=iline+1
              iwrite=0
              call pstring(' P8incl',7,P8incl,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.344)then  !P8Omega, tag ky
              iline=iline+1
              iwrite=0
              call pstring(' P8Omega',6,P8Omega,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.345)then  !P8Q, tag kz
              iline=iline+1
              iwrite=0
              call pstring(' P8Q',6,P8Q,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.321)then  !P8ratrad, tag kb
              iline=iline+1
              iwrite=0
              call pstring(' P8ratrad',7,P8ratrad,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif

c
c

            if(kk.eq.1208)then  !Tgrav1, tag g1
              iline=iline+1
              iwrite=0
              call pstring(' Tgrav1',5,Tgrav1,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1591)then  !contamS0, tag s0
              iline=iline+1
              iwrite=0
              call pstring(' contamS0',5,contamS0,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
c              line=line(1:K)//string24
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1592)then  !contamS1, tag s1
              iline=iline+1
              iwrite=0
              call pstring(' contamS1',5,contamS1,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
c              line=line(1:K)//string24
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1593)then  !contamS2, tag s2
              iline=iline+1
              iwrite=0
              call pstring(' contamS2',5,contamS2,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
c              line=line(1:K)//string24
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1594)then  !contamS3, tag s3
              iline=iline+1
              iwrite=0
              call pstring(' contamS3',5,contamS3,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
c              line=line(1:K)//string24
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c



            if(kk.eq.110)then  !omegadot, tag do
              iline=iline+1
              iwrite=0
c              ilength=ilength+24
c              write(string24,7832)osine
              call pstring(' omega_dot',8,omegadot,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
c              line=line(1:K)//string24
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1209)then  !Tgrav2, tag g2
              iline=iline+1
              iwrite=0
c              ilength=ilength+24
c              write(string24,7832)osine
              call pstring(' Tgrav2',5,Tgrav2,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
c              line=line(1:K)//string24
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   November 18, 2012
c
c   Third body variables here
c
            if(kk.eq.627)then          !tertperiod, tag tt
              iline=iline+1
              iwrite=0
c              ilength=ilength+20
c              write(string20,2783)tertPeriod
              call pstring(' tertP',9,tertperiod,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
c              line=line(1:K)//string20
              line=line(1:K)//outstring
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.628)then          !tertT0, tag tu
              iline=iline+1
              iwrite=0
c              ilength=ilength+21
c              write(string21,2784)tertT0
              call pstring(' tertT0',6,tertT0,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string21
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.617)then          !tertT0, tag tj
              iline=iline+1
              iwrite=0
c              ilength=ilength+21
c              write(string21,2784)tertT0
              call pstring(' tertconj',6,tertconj,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string21
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.629)then          !tertecos, tag tv
              iline=iline+1
              iwrite=0
c              ilength=ilength+21
c              write(string21,2785)tertecos
              call pstring(' terte_cos',7,tertecos,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string21
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.630)then          !tertecos, tag tw
              iline=iline+1
              iwrite=0
c              ilength=ilength+21
c              write(string21,2786)tertesin
              call pstring(' terte_sin',7,tertesin,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string21
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.631)then          !tertincl, tag tx
              iline=iline+1
              iwrite=0
c              ilength=ilength+20
c              write(string20,2787)tertincl
              call pstring(' tertincl',7,tertincl,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string20
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c

            if(kk.eq.632)then          !tertOmega, tag ty
              iline=iline+1
              iwrite=0
c              ilength=ilength+21
c              write(string21,2788)tertOmega
              call pstring(' tertOmega',6,tertOmega,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string21
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.633)then          !tertQ, tag tz
              iline=iline+1
              iwrite=0
c              ilength=ilength+18
c              write(string18,2789)tertQ
              call pstring(' tertQ',6,tertQ,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string18
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.78)then          !contam, tag co
              iline=iline+1
              iwrite=0
c              ilength=ilength+15
c              write(string15,783)contam
              call pstring(' contam',5,contam,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string15
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1144)then          !beam1, tag e1
              iline=iline+1
              iwrite=0
c              ilength=ilength+14
c              write(string14,781)beam1
              call pstring(' beam1',3,beam1,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string14
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1145)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+14
c              write(string14,782)beam2
              call pstring(' beam2',3,beam2,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string14
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.492)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+19
c              write(string19,777)primmass
              call pstring(' primmass',5,primmass,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string19
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.497)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+19
c              write(string19,778)primrad
              call pstring(' primrad',5,primrad,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string19
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.612)then    ! temprat
              iline=iline+1
              iwrite=0
c              ilength=ilength+19
c              write(string19,1778)temprat
              call pstring(' temprat',7,temprat,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string19
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.490)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+19
c              write(string19,779)primK
              call pstring(' primK',6,primK,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string19
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.544)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+19
c              write(string19,776)ratrad
              call pstring(' ratrad',6,ratrad,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string19
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1528)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+16
c              write(string16,775)frac1
              call pstring(' frac1',6,frac1,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string16
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.100)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+19
c              write(string19,780)density
              call pstring(' density',7,density,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string19
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1529)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+16
c              write(string16,774)frac2
              call pstring(' frac2',6,frac2,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string16
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c 
c   UPDATE January 15, 2002
c
c   Here are the blocks for alb1 and alb2
c
            if(kk.eq.1368)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+12
c              write(string12,2202)alb1
              call pstring(' al1',5,alb1,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string12
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1369)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+12
c              write(string12,2203)alb2
              call pstring(' al2',5,alb2,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string12
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.111)then   !ecosw, use string dphi
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3202)ecosw
              call pstring(' dphi',5,ecosw,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
              iargper=100
            endif
c

            if(kk.eq.610)then           !Tconj string
              iline=iline+1
              iwrite=0
c              ilength=ilength+19
c              write(string19,7701)Tconj
              call pstring(' Tconj',6,Tconj,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string19
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1623)then           !T0 string
              iline=iline+1
              iwrite=0
c              ilength=ilength+16
c              write(string16,701)T0
              call pstring(' T0',6,T0,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string16
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.484)then      !period string
              iline=iline+1
              iwrite=0
c              ilength=ilength+17
c              write(string17,700)period
              call pstring(' P',9,period,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string17
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.269)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+10
c              write(string10,200)finc
              call pstring(' i',4,finc,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string10
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   UPDATE May 8, 2006
c
c   Add bigI and bigbeta
c
            if(kk.eq.8)then          !bigI
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,88200)bigI
              call pstring(' bigI',4,bigI,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1)then          !bigbeta
              iline=iline+1
              iwrite=0
c              ilength=ilength+17
c              write(string17,88201)bigbeta
              call pstring(' bigbeta',4,bigbeta,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string17
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif

c
c
c   UPDATE JULY @7, 2004
c
c   Change Q to string length 14 (f11.8).
c
            if(kk.eq.384)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+14
c              write(string14,201)Q
              call pstring(' Q',8,Q,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string14
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   UPDATE AUGUST 2, 2004
c
c   Make fill1 and fill2 string length 14
c
            if(kk.eq.1176)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+14
c              write(string14,202)fill1
              call pstring(' fi1',7,fill1,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string14
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1177)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+14
c              write(string14,203)fill2
              call pstring(' fi2',7,fill2,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string14
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.552)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+12
c              write(string12,204)rinner
              call pstring(' rin',5,rinner,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string12
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1464)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+12
c              write(string12,205)omega1
              call pstring(' o1',4,omega1,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string12
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1465)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+12
c              write(string12,206)omega2
              call pstring(' o2',4,omega2,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string12
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.558)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,207)router
              call pstring(' rout',5,router,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1626)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,217)t3  
              call pstring(' Teff3',2,t3,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1210)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,218)g3  
              call pstring(' g3',3,g3,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.36)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+12
c              write(string12,211)betarim
              call pstring(' beta',3,betarim,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string12
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1624)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+16
c              write(string16,209)Teff1
              call pstring(' Teff1',2,Teff1,outstring,lll)
              K=lnblnk(line)
              line=line(1:K)//outstring
c              line=line(1:K)//string16
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1625)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+16
c              write(string16,210)Teff2
              call pstring(' Teff2',2,Teff2,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string16
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.744)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+11
c              write(string11,212)xi
              call pstring(' xi',3,xi,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string11
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.611)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,208)tdisk 
              call pstring(' Td',2,tdisk,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.375)then
              iline=iline+1
              iwrite=0
              ilength=ilength+19
              write(string19,213)rLx
              call pstring(' Td',2,tdisk,outstring,lll)
              K=lnblnk(line)
              line=line(1:K)//string19
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.580)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+14
c              write(string14,214)separ
              call pstring(' separ',6,separ,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string14
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.576)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+14
c              write(string14,219)SA3
              call pstring(' SA3',6,SA3,outstring,lll)
              K=lnblnk(line)
              line=line(1:K)//outstring
c              line=line(1:K)//string14
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.498)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+15
c              write(string15,220)pshift
              call pstring(' pshift',5,pshift,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string15
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.130)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+12
c              write(string12,221)ecc
              call pstring(' e',6,ecc,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string12
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.17)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+15
c              write(string15,222)argper
              call pstring(' argper',3,argper,outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string15
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1048)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,401)spot1parm(1,1)
              call pstring(' TF_spot1_star1',4,spot1parm(1,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1049)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,402)spot1parm(1,2)
              call pstring(' lat_spot1_star1',2,spot1parm(1,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1050)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,403)spot1parm(1,3)
              call pstring(' lon_spot1_star1',2,spot1parm(1,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1051)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,404)spot1parm(1,4)
              call pstring(' rad_spot1_star1',3,spot1parm(1,4),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1052)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,405)spot1parm(2,1)
              call pstring(' TF_spot2_star1',4,spot1parm(2,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1053)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,406)spot1parm(2,2)
              call pstring(' lat_spot2_star1',2,spot1parm(2,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1054)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,407)spot1parm(2,3)
              call pstring(' lon_spot2_star1',2,spot1parm(2,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1055)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,408)spot1parm(2,4)
              call pstring(' rad_spot2_star1',3,spot1parm(2,4),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1080)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,501)spot2parm(1,1)
              call pstring(' TF_spot1_star2',4,spot2parm(1,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1081)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,502)spot2parm(1,2)
              call pstring(' lat_spot1_star2',2,spot2parm(1,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1082)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,503)spot2parm(1,3)
              call pstring(' lon_spot1_star2',2,spot2parm(1,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1083)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,504)spot2parm(1,4)
              call pstring(' rad_spot1_star2',3,spot2parm(1,4),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1084)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,505)spot2parm(2,1)
              call pstring(' TF_spot2_star2',4,spot2parm(2,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1085)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,506)spot2parm(2,2)
              call pstring(' lat_spot2_star2',2,spot2parm(2,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1086)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,507)spot2parm(2,3)
              call pstring(' lon_spot2_star2',2,spot2parm(2,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1087)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+23
c              write(string23,508)spot2parm(2,4)
              call pstring(' rad_spot2_star2',3,spot2parm(2,4),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string23
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1112)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+22
c              write(string22,601)spotdparm(1,1)
              call pstring(' TF_spot1_disk',4,spotdparm(1,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string22
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1113)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+22
c              write(string22,602)spotdparm(1,2)
              call pstring(' azi_spot1_disk',2,spotdparm(1,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string22
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1114)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+22
c              write(string22,603)spotdparm(1,3)
              call pstring(' cut_spot1_disk',4,spotdparm(1,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string22
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1115)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+22
c              write(string22,604)spotdparm(1,4)
              call pstring(' wid_spot1_disk',2,spotdparm(1,4),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string22
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1116)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+22
c              write(string22,605)spotdparm(2,1)
              call pstring(' TF_spot2_disk',4,spotdparm(2,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string22
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1117)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+22
c              write(string22,606)spotdparm(2,2)
              call pstring(' azi_spot2_disk',2,spotdparm(2,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string22
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1118)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+22
c              write(string22,607)spotdparm(2,3)
              call pstring(' cut_spot2_disk',4,spotdparm(2,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string22
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1119)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+22
c              write(string22,608)spotdparm(2,4)
              call pstring(' wid_spot2_disk',3,spotdparm(2,4),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string22
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   UPDATE November 6, 2002
c
c   Add the limb darkening coefficients here.
c
c   x-coefficient, star 1
c
            if(kk.eq.1752)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,801)dwavex(1,1)
              call pstring(' x1(U)',3,dwavex(1,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1753)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,802)dwavex(2,1)
              call pstring(' x1(B)',3,dwavex(2,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1754)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,803)dwavex(3,1)
              call pstring(' x1(V)',3,dwavex(3,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1755)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,804)dwavex(4,1)
              call pstring(' x1(R)',3,dwavex(4,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1756)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,805)dwavex(5,1)
              call pstring(' x1(I)',3,dwavex(5,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1757)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,806)dwavex(6,1)
              call pstring(' x1(J)',3,dwavex(6,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1758)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,807)dwavex(7,1)
              call pstring(' x1(H)',3,dwavex(7,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1759)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,808)dwavex(8,1)
              call pstring(' x1(K)',3,dwavex(8,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   x-coefficient, star 2
c
            if(kk.eq.1816)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,2801)dwavex(1,2)
              call pstring(' x2(U)',3,dwavex(1,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1817)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,2802)dwavex(2,2)
              call pstring(' x2(B)',3,dwavex(2,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1818)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,2803)dwavex(3,2)
              call pstring(' x2(V)',3,dwavex(3,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1819)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,2804)dwavex(4,2)
              call pstring(' x2(R)',3,dwavex(4,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1820)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,2805)dwavex(5,2)
              call pstring(' x2(I)',3,dwavex(5,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1821)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,2806)dwavex(6,2)
              call pstring(' x2(J)',3,dwavex(6,2),
     @            outstring,lll)
              K=lnblnk(line)
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1822)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,2807)dwavex(7,2)
              call pstring(' x2(H)',3,dwavex(7,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1823)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,2808)dwavex(8,2)
              call pstring(' x2(K)',3,dwavex(8,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   y-coefficient, star 1
c
            if(kk.eq.1784)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,1801)dwavey(1,1)
              call pstring(' y1(U)',3,dwavey(1,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1785)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,1802)dwavey(2,1)
              call pstring(' y1(B)',3,dwavey(2,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1786)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,1803)dwavey(3,1)
              call pstring(' y1(V)',3,dwavey(3,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1787)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,1804)dwavey(4,1)
              call pstring(' y1(R)',3,dwavey(4,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1788)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,1805)dwavey(5,1)
              call pstring(' y1(I)',3,dwavey(5,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1789)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,1806)dwavey(6,1)
              call pstring(' y1(J)',3,dwavey(6,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1790)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,1807)dwavey(7,1)
              call pstring(' y1(H)',3,dwavey(7,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1791)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,1808)dwavey(8,1)
              call pstring(' y1(K)',3,dwavey(8,1),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   y-coefficient, star 2
c
            if(kk.eq.1720)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3801)dwavey(1,2)
              call pstring(' y2(U)',3,dwavey(1,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1721)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3802)dwavey(2,2)
              call pstring(' y2(B)',3,dwavey(2,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1722)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3803)dwavey(3,2)
              call pstring(' y2(V)',3,dwavey(3,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1723)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3804)dwavey(4,2)
              call pstring(' y2(R)',3,dwavey(4,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1724)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3805)dwavey(5,2)
              call pstring(' y2(I)',3,dwavey(5,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1725)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3806)dwavey(6,2)
              call pstring(' y2(J)',3,dwavey(6,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1726)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3807)dwavey(7,2)
              call pstring(' y2(H)',3,dwavey(7,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1727)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3808)dwavey(8,2)
              call pstring(' y2(K)',3,dwavey(8,2),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c   x-coefficient, star 3
c
            if(kk.eq.1400)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,801)dwavex(1,1)
              call pstring(' x3(U)',3,dwavex(1,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1401)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,802)dwavex(2,1)
              call pstring(' x3(B)',3,dwavex(2,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1402)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,803)dwavex(3,1)
              call pstring(' x3(V)',3,dwavex(3,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1403)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,804)dwavex(4,1)
              call pstring(' x3(R)',3,dwavex(4,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1404)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,805)dwavex(5,1)
              call pstring(' x3(I)',3,dwavex(5,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1405)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,806)dwavex(6,1)
              call pstring(' x3(J)',3,dwavex(6,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1406)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,807)dwavex(7,1)
              call pstring(' x3(H)',3,dwavex(7,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1407)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,808)dwavex(8,1)
              call pstring(' x3(K)',3,dwavex(8,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c
c   y-coefficient, star 3
c
            if(kk.eq.1432)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3801)dwavey(1,2)
              call pstring(' y3(U)',3,dwavey(1,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1433)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3802)dwavey(2,2)
              call pstring(' y3(B)',3,dwavey(2,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1434)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3803)dwavey(3,2)
              call pstring(' y3(V)',3,dwavey(3,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1435)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3804)dwavey(4,2)
              call pstring(' y3(R)',3,dwavey(4,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1436)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3805)dwavey(5,2)
              call pstring(' y3(I)',3,dwavey(5,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1437)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3806)dwavey(6,2)
              call pstring(' y3(J)',3,dwavey(6,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1438)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3807)dwavey(7,2)
              call pstring(' y3(H)',3,dwavey(7,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(kk.eq.1439)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+13
c              write(string13,3808)dwavey(8,2)
              call pstring(' y3(K)',3,dwavey(8,3),
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string13
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
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
              call getom(ecc,ecosw,argper)
              iline=iline+1
              iwrite=0
c              ilength=ilength+15
c              write(string15,222)argper
              call pstring(' argper',3,argper,
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string15
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
            if(isw29.gt.0)then
              iline=iline+1
              iwrite=0
c              ilength=ilength+12
c              write(string12,221)ecc
              call pstring(' e',7,ecc,
     @            outstring,lll)
              K=lnblnk(line)
              ilength=ilength+lll
              line=line(1:K)//outstring
c              line=line(1:K)//string12
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
              iline=iline+1
              iwrite=0
c              ilength=ilength+15
c              write(string15,222)argper
              call pstring(' argper',3,argper,outstring,lll)
              K=lnblnk(line)
              line=line(1:K)//outstring
c              line=line(1:K)//string15
              if(ilength.gt.45)then
                iwrite=100
                write(*,500)line
                ilength=2
                line=' *'
              endif
            endif
c
c UPDATE June 10, 2003
c
c Add two more digits to gamma, use string 15
c
          if(isvel1.gt.0)then
            iline=iline+1
            iwrite=0
c            ilength=ilength+15
c            write(string15,215)gamma1
              call pstring(' gam1',4,gamma1,outstring,lll)
            K=lnblnk(line)
            ilength=ilength+lll
            line=line(1:K)//outstring
c            line=line(1:K)//string15
            if(ilength.gt.48)then
              iwrite=100
              write(*,500)line
              ilength=2
                line=' *'
            endif
          endif
c
          if(isvel2.gt.0)then
            iline=iline+1
            iwrite=0
c            ilength=ilength+15
c            write(string15,216)gamma2
              call pstring(' gam2',4,gamma2,outstring,lll)
            K=lnblnk(line)
            ilength=ilength+lll
            line=line(1:K)//outstring
c            line=line(1:K)//string15
            if(ilength.gt.48)then
              iwrite=100
              write(*,500)line
              ilength=2
                line=' *'
            endif
          endif
c
          if(isvel3.gt.0)then
            iline=iline+1
            iwrite=0
c            ilength=ilength+15
c            write(string15,216)gamma3
              call pstring(' gam3',4,gamma3,outstring,lll)
            K=lnblnk(line)
            ilength=ilength+lll
            line=line(1:K)//outstring
c            line=line(1:K)//string15
            if(ilength.gt.48)then
              iwrite=100
              write(*,500)line
              ilength=2
                line=' *'
            endif
          endif
c
          if(ilength.gt.8.and.iwrite.eq.0)write(*,500)line
c
 213      format(' Lx/Lopt=',1pe10.4)
c
 500      format(a80)
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&
c
          subroutine evalfit(islc,Nmodel,xmodel,yinput,
     $      Ndata,xdata,index,yout,bestoff)
c
c   This routine will return the value of the model at a specific
c   datapoint
c
          implicit double precision (a-h,o-z)
c
          dimension xmodel(Nmodel),yinput(Nmodel),xdata(Ndata)
          dimension yinter(500000)  !should match Nmaxphase
          dimension ydummy(500000),y2(500000),ymodel(500000) 
c
          do 10 i=1,Nmodel
            if(islc.gt.0)then
c
c  Check to see that the argument of log10 is positive.
c
              if(yinput(i).gt.0.0d0)then
                ymodel(i)=-2.5d0*dlog10(yinput(i))
              else
                ymodel(i)=-99.9d0
              endif
            else
              ymodel(i)=yinput(i)
            endif
 10       continue
c
          call spline(xmodel,ymodel,Nmodel,0.0d0,0.0d0,y2)
c
          do 30 i=1,Ndata
            call splint(xmodel,ymodel,y2,Nmodel,xdata(i),qqq)
            yinter(i)=qqq
 30       continue
c
          call offset(Ndata,yinter,ydummy,bestoff)
c
          yout=ydummy(index)
c
          return
          end
c
c    &&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c   December 22, 2012
c
c   comment out
c
c       SUBROUTINE MATINV(ARRAY,NORDER,DET)
c       implicit double precision(a-h)
c       implicit double precision(o-z)
c       double precision array(20,20)
c       DIMENSION IK(20),JK(20)
c10     DET=1.d0
c11     DO 100 K=1,NORDER
cC        FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX
c       AMAX=0.d0
c21     DO 30 I=K,NORDER
c       DO 30 J=K,NORDER
c23     IF(DABS(AMAX)-DABS(ARRAY(I,J))) 24,24,30
c24     AMAX=ARRAY(I,J)
c       IK(K)=I
c       JK(K)=J
c30     CONTINUE
cC        INTERCHANGE ROWS AND COLUMNS TO PUT AMAX IN ARRAY(K,K)
c31     IF(AMAX) 41,32,41
c32     DET=0.d0
c       GOTO 140
c41     I=IK(K)
c       IF(I-K) 21,51,43
c43     DO 50 J=1,NORDER
c       SAVE=ARRAY(K,J)
c       ARRAY(K,J)=ARRAY(I,J)
c50     ARRAY(I,J)=-SAVE
c51     J=JK(K)
c       IF(J-K) 21,61,53
c53     DO 60 I=1,NORDER
c       SAVE=ARRAY(I,K)
c       ARRAY(I,K)=ARRAY(I,J)
c60     ARRAY(I,J)=-SAVE
cC        ACCUMULATE ELEMENTS OF INVERSE MATRIX
c61     DO 70 I=1,NORDER
c       IF(I-K) 63,70,63
c63     ARRAY(I,K)=-ARRAY(I,K)/AMAX
c70     CONTINUE
c71     DO 80 I=1,NORDER
c       DO 80 J=1,NORDER
c       IF(I-K) 74,80,74
c74     IF(J-K) 75,80,75
c75     ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J)
c80     CONTINUE
c81     DO 90 J=1,NORDER
c       IF(J-K) 83,90,83
c83     ARRAY(K,J)=ARRAY(K,J)/AMAX
c90     CONTINUE
c       ARRAY(K,K)=1.d0/AMAX
c100    DET=DET*AMAX
cC        RESTORE ORDERING OF MATRIX
c101    DO 130 L=1,NORDER
c       K=NORDER-L+1
c       J=IK(K)
c       IF(J-K) 111,111,105
c105    DO 110 I=1,NORDER
c       SAVE=ARRAY(I,K)
c       ARRAY(I,K)=-ARRAY(I,J)
c110    ARRAY(I,J)=SAVE
c111    I=JK(K)
c       IF(I-K) 130,130,113
c113    DO 120 J=1,NORDER
c       SAVE=ARRAY(K,J)
c       ARRAY(K,J)=-ARRAY(I,J)
c120    ARRAY(I,J)=SAVE
c130    CONTINUE
c140    RETURN
c       END
cc
cc  #########$%$%$#####$%$%$#$$$&&
c
      SUBROUTINE GAUSSJ(A,N,NP,B,M)
      implicit double precision (a-h,o-z)
      integer m,n,np,nmax
      real*8 a(np,np),b(np,np)
      PARAMETER (NMAX=50)
c      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      integer i,icol,irow,j,k,l,ll,indxc(nmax),indxr(nmax),ipiv(nmax)
      real*8 big,dum,pivinv
      irow=1
      icol=1
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                write(*,*)'Singular matrix in gaussj'
                stop
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.)then
           write(*,*) 'Singular matrix in gaussj.'
        endif
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE 
        ENDIF
24    CONTINUE
      RETURN
      END
c
c  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
          subroutine obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c  July 24, 2000
c
c  This routine will take the observed variables specified in the
c  gridloop.opt file and compute the additional chi^2
c
c  UPDATE September 11, 2001
c
c  Change the dimensions of obsparm, obv, sobv, and eobv to 9.  The duration
c  of the X-ray eclipse in degrees is in obsparm(9).
c
c
c   UPDATE September 21, 2008
c
c   Change the dimension of obsparm, obv, sobv, and eobv  to 11
c
c   UPDATE October 10, 2008
c
c   Make the dimensions of obv, obsparm, sobv, eobv to 17
c
c   UPDATE March 15, 2011
c
c   Make the dimensions 18 (add sum of radii)
   
          implicit double precision (a-h,o-z)
c
          dimension obv(19),eobv(19),obsparm(19)
          character*40 sobv(19)
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
          index=1
          ochi=0.0d0
          do 1001 i=1,Nobv
c
c  UPDATE NOVEMBER 23, 2001
c
c  make the argument sobv(i)(1:2)
c
            index=-99
            icn=icnvrt(sobv(i)(1:2))
            if(icn.eq.1400)index=1
            if(icn.eq.1401)index=5
            if(icn.eq.1560)index=2
            if(icn.eq.1561)index=6
            if(icn.eq.1208)index=3
c
c   UPDATE September 21, 2008
c
c   Here are the strings for K1 and K2
c
            if(icn.eq.1336)index=10
            if(icn.eq.1337)index=11
c
c   UPDATE October 10, 2008
c
c   Here are the strings for incl, mass, ecc, arg, t1, t2, t3
c
            if(icn.eq.269)index=12
            if(icn.eq.384)index=13
            if(icn.eq.130)index=14
            if(icn.eq.17)index=15
            if(icn.eq.1624)index=16
            if(icn.eq.1625)index=17
            if(icn.eq.1626)index=19
c
c   UPDATE March 15, 2011
c
c   Here is the string for the fractional sum (string su)
c
            if(icn.eq.596)index=18

c
c   UPDATE February 21, 2002
c
c   Change icn.eq.208  to  icn.eq.209
c
            if(icn.eq.1209)index=7
            if(icn.eq.1688)index=4
            if(icn.eq.1689)index=8
c
c  UPDATE September 11, 2001
c
c  Add this if block.
c
            if(icn.eq.740)index=9
c
c   NEW BUG ALERT  July 24, 2001
c
c   Add this if-then block
c
            if(icn.eq.104)index=-1   !di
            if(icn.eq.1112)index=-1
            if(icn.eq.1113)index=-1
            if(icn.eq.1114)index=-1
            if(icn.eq.1115)index=-1
            if(icn.eq.1116)index=-1
            if(icn.eq.1117)index=-1
            if(icn.eq.1118)index=-1
            if(icn.eq.1119)index=-1
c
c   UPDATE FEBRUARY 4, 2005
c
c
c   Add this if-then block
c
            if(icn.eq.369)index=-1   !lr
c
c   UPDATE May 27, 2002
c
c   If the variable rmed > 1, then define chisq as the absolute deviation,
c   rather than the normal chi^2.
c
c   UPDATE June 3, 2002
c
c   Add this if-then statement.   See below for details.
c
            if((index.eq.9).and.(obv(i).lt.0.0d0))go to 1234
            if(index.ge.1)then
              if(rmed.ge.1.0d0)then
                chisq=dabs((obv(i)-obsparm(index))/(eobv(i)))
              else
                chisq=(obv(i)-obsparm(index))*(obv(i)-obsparm(index))
     $            /(eobv(i)*eobv(i))
              endif
              ochi=ochi+chisq
            endif
c
c   UPDATE June 3, 2002
c
c   If the observed duration of the X-ray eclipse is negative, then
c   the number is meant as an upper limit.  That is, if obv = -10.0,
c   then the eclipse duration is less than 10 degrees.
c
 1234       if((index.eq.9).and.(obv(i).lt.0.0d0))then
c
c   Here is the case where the upper limit is less than the computed
c   eclipse duration.
c
              if(dabs(obv(i)).gt.obsparm(index))go to 1001
c
              if(rmed.ge.1.0d0)THEN
                chisq=dabs((dabs(obv(i))-obsparm(index))/(eobv(i)))
              else
                chisq=(dabs(obv(i))-obsparm(index))*
     %            (dabs(obv(i))-obsparm(index))/(eobv(i)*eobv(i))
              endif
              ochi=ochi+chisq
            endif
c
c 
c
 1001     continue
c
          return
          end
c
c  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
          subroutine wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
c  July 24, 2000
c
c  This routine will take the observed variables specified in the
c  gridloop.opt file and compute the additional chi^2
c
c
c  UPDATE September 11, 2001
c
c  Change the dimensions of obsparm, obv, sobv, and eobv to 9.  The duration
c  of the X-ray eclipse in degrees is in obsparm(9).
c
c
c   UPDATE September 21, 2008
c
c   Change the dimension of obsparm, obv, sobv, and eobv  to 11
c
c   UPDATE October 10, 2008
c
c   Change the dimension of obsparm, obv, sobv, and eobv to 17
c
c   UPDATE March 15, 2011
c
c   Make the dimension 18 (add sum of fractional radii)
c
          implicit double precision (a-h,o-z)
c
          dimension obv(19),eobv(19),obsparm(19)
          character*40 sobv(19)
c
          character*30 section(19),outstring
          character*80 line
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
          iline=0
          iwrite=0
          ilength=0
c
          do 10 i=1,10
            section(i)='12345678901234567890123456789'
            section(i)='                             '
 10       continue
          line='  -->'
c
          ochi=0.0d0
          do 1001 i=1,Nobv
c
c  UPDATE NOVEMBER 23, 2001
c
c  make the argument sobv(i)(1:2)
c
            icn=icnvrt(sobv(i)(1:2))
c
c   NEW BUG  July 24, 2001
c
c   Escape if icn=104 (disk fraction requested).
c
            if(icn.eq.104)go to 1001        !di
c
c   UPDATE FEBRUARY 4, 2005
c
c   Escape if icn=369 (luminosity ratio requested).
c
            if(icn.eq.369)go to 1001       !lr
c
c   END BUG
c
c   NEW BUG August 2, 2001
c 
c   Put an escape if one of the variables is not valid.
c
            index=-1
            if(icn.eq.1400)index=1   !m1
            if(icn.eq.1401)index=5
            if(icn.eq.1560)index=2    !r1
            if(icn.eq.1561)index=6
            if(icn.eq.1208)index=3    !g1
c
c   UPDATE September 21, 2008
c
c   Add K1 and K2
c
            if(icn.eq.1336)index=10
            if(icn.eq.1337)index=11
c
c
c   UPDATE October 10, 2008
c
c   Here are the strings for incl, mass, ecc, arg, t1, t2
c
            if(icn.eq.269)index=12
            if(icn.eq.384)index=13
            if(icn.eq.130)index=14
            if(icn.eq.17)index=15
            if(icn.eq.1624)index=16
            if(icn.eq.1625)index=17
            if(icn.eq.1626)index=19
c
c   UPDATE March 15, 2011
c
c   Add sum of fractional radii (string su)
c
            if(icn.eq.596)index=18
c
c   UPDATE February 21, 2002
c
c   Change icn.eq.208  to  icn.eq.209
c
            if(icn.eq.1209)index=7
            if(icn.eq.1688)index=4   !v1
            if(icn.eq.1689)index=8    !v2
c
c   UPDATE September 11, 2001
c
c   Add this if block
c
            if(icn.eq.740)index=9
c
c   Here is the escape.
c
            if(index.lt.0)go to 1001
c
c   UPDATE May 27, 2002
c
c   If the variable rmed > 1, then define chisq as the absolute
c   deviation, rather than the normal chi^2.
c
c
c   UPDATE June 3, 2002
c
c   If the observed duration of the X-ray eclipse is negative, then
c   the number is meant as an upper limit.  That is, if obv = -10.0,
c   then the eclipse duration is less than 10 degrees.
c
            if(index.eq.9)then
              if(obv(i).lt.0.0d0)then
c
c   Here is the case where the upper limit is less than the computed
c   eclipse duration.
c
                if(dabs(obv(i)).gt.obsparm(index))then
                  chisq=0.0d0
                  iline=iline+1
                  go to 1234
                else
                  if(rmed.ge.1.0d0)THEN
                    chisq=dabs((dabs(obv(i))-obsparm(index))/(eobv(i)))
                  else
                    chisq=(dabs(obv(i))-obsparm(index))*
     @                (dabs(obv(i))-obsparm(index))/(eobv(i)*eobv(i))
                  endif
                  iline=iline+1
                  go to 1234
                endif  
              else
                if(rmed.ge.1.0d0)THEN
                  chisq=dabs((dabs(obv(i))-obsparm(index))/(eobv(i)))
                else
                  chisq=(dabs(obv(i))-obsparm(index))*
     @                (dabs(obv(i))-obsparm(index))/(eobv(i)*eobv(i))
                endif
                iline=iline+1
                go to 1234
              endif
            endif
c
            if(rmed.ge.1.0d0)then
              chisq=dabs((obv(i)-obsparm(index))/(eobv(i)))
            else
              chisq=(obv(i)-obsparm(index))*(obv(i)-obsparm(index))
     $          /(eobv(i)*eobv(i))
            endif
            iline=iline+1
c
c   UPDATE May 27, 2002
c
c   If the variable rmed > 1, then define chisq as the absolute
c   deviation, rather than the normal chi^2.  In that case, write
c   using different format statements.
c
 1234       if(rmed.ge.1.0d0)then
              if(index.eq.1)then
                 call chistring(' med(m1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),91)chisq
                iwrite=0
              endif
              if(index.eq.2)then
                 call chistring(' med(r1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),92)chisq
                iwrite=0
              endif
              if(index.eq.3)then
                 call chistring(' med(g1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),93)chisq
                iwrite=0
              endif
              if(index.eq.4)then
                 call chistring(' med(v1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),94)chisq
                iwrite=0
              endif
              if(index.eq.5)then
                 call chistring(' med(m2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),95)chisq
                iwrite=0
              endif
              if(index.eq.6)then
                 call chistring(' med(r2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),96)chisq
                iwrite=0
              endif
              if(index.eq.7)then
                 call chistring(' med(g2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),97)chisq
                iwrite=0
              endif
              if(index.eq.8)then
                 call chistring(' med(v2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),98)chisq
                iwrite=0
              endif
              if(index.eq.9)then
                 call chistring(' med(duras)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),99)chisq
                iwrite=0
              endif
              if(index.eq.10)then
                 call chistring(' med(k1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),302)chisq
                iwrite=0
              endif
              if(index.eq.11)then
                 call chistring(' med(k2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),303)chisq
                iwrite=0
              endif
              if(index.eq.12)then
                 call chistring(' med(incl)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),992)chisq
                iwrite=0
              endif
              if(index.eq.13)then
                 call chistring(' med(Q)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),993)chisq
                iwrite=0
              endif
              if(index.eq.14)then
                 call chistring(' med(ecc)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),994)chisq
                iwrite=0
              endif
              if(index.eq.15)then
                 call chistring(' med(arg)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),995)chisq
                iwrite=0
              endif
              if(index.eq.16)then
                 call chistring(' med(t1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),996)chisq
                iwrite=0
              endif
              if(index.eq.17)then
                 call chistring(' med(t2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),997)chisq
                iwrite=0
              endif
              if(index.eq.18)then
                 call chistring(' med(sum)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),998)chisq
                iwrite=0
              endif
              if(index.eq.19)then
                 call chistring(' med(t3)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),998)chisq
                iwrite=0
              endif
            else
              if(index.eq.1)then
                 call chistring(' chi(m1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),1)chisq
                iwrite=0
              endif
              if(index.eq.2)then
                 call chistring(' chi(r1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),2)chisq
                iwrite=0
              endif
              if(index.eq.3)then
                 call chistring(' chi(g1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),3)chisq
                iwrite=0
              endif
              if(index.eq.4)then
                 call chistring(' chi(v1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),4)chisq
                iwrite=0
              endif
              if(index.eq.5)then
                 call chistring(' chi(m2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),5)chisq
                iwrite=0
              endif
              if(index.eq.6)then
                 call chistring(' chi(r2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),6)chisq
                iwrite=0
              endif
              if(index.eq.7)then
                 call chistring(' chi(g2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),7)chisq
                iwrite=0
              endif
              if(index.eq.8)then
                 call chistring(' chi(v2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),8)chisq
                iwrite=0
              endif
              if(index.eq.9)then
                 call chistring(' chi(duras)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),9)chisq
                iwrite=0
              endif
              if(index.eq.10)then
                 call chistring(' chi(k1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),300)chisq
                iwrite=0
              endif
              if(index.eq.11)then
                 call chistring(' chi(k2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),301)chisq
                iwrite=0
              endif
              if(index.eq.12)then
                 call chistring(' chi(inc)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),882)chisq
                iwrite=0
              endif
              if(index.eq.13)then
                 call chistring(' chi(Q)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),883)chisq
                iwrite=0
              endif
              if(index.eq.14)then
                 call chistring(' chi(ecc)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),884)chisq
                iwrite=0
              endif
              if(index.eq.15)then
                 call chistring(' chi(arg)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),885)chisq
                iwrite=0
              endif
              if(index.eq.16)then
                 call chistring(' chi(t1)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),886)chisq
                iwrite=0
              endif
              if(index.eq.17)then
                 call chistring(' chi(t2)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),887)chisq
                iwrite=0
              endif
              if(index.eq.18)then
                 call chistring(' chi(sum)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),888)chisq
                iwrite=0
              endif
              if(index.eq.19)then
                 call chistring(' chi(t3)',chisq,outstring,lll)
                 ilength=ilength+lll
                 section(iline)=outstring
c                write(section(iline),888)chisq
                iwrite=0
              endif
            endif
c
            K=lnblnk(line)
            line=line(1:K)//section(iline)
            ochi=ochi+chisq
c            if(iline.eq.3)then
            if(ilength.ge.55)then
              write(*,200)line
              iwrite=100
              iline=0
              ilength=0
              line='  -->'
            endif
 1001     continue
c
          if(iline.gt.0.and.iwrite.eq.0)write(*,200)line
c
 200      format(a80)
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&
c
          subroutine checkRVfit(Nmodel1,xmodel1,ymodel1,Nmodel2,
     @      xmodel2,ymodel2,Ndata1,xdata1,ydata1,err1,Ndata2,xdata2,
     @      ydata2,err2,chisq1,chisq2,zero,resRV1,resRV2)
c
c   April 11, 2001
c
c   This routine will fit two velocity curves simultaneously, using
c   a common gamma as a free parameter.  
c
          implicit double precision (a-h,o-z)
c
           dimension xmodel1(Nmodel1),ymodel1(Nmodel1),
     @      xdata1(Ndata1),ydata1(Ndata1),
     %      err1(Ndata1),ymodel2(Nmodel2),xmodel2(Nmodel2),
     @      xdata2(Ndata2),ydata2(Ndata2),err2(Ndata2)
          dimension resRV1(Ndata1),resRV2(Ndata2)  !match Nmaxphase below
          dimension yinter1(500000),y21(500000),ydummy1(500000)
          dimension yinter2(500000),y22(500000),ydummy2(500000)
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
c   Find the maximum and minimum y-values of the data
c
          ymin1=1.0d20
          ymax1=-1.0d20
          do 20 i=1,Ndata1
            if(ydata1(i).lt.ymin1)ymin1=ydata1(i)
            if(ydata1(i).gt.ymax1)ymax1=ydata1(i)
 20       continue
c
          ymin2=10000.0d0
          ymax2=-10000.0d0
          do 21 i=1,Ndata2
            if(ydata2(i).lt.ymin2)ymin2=ydata2(i)
            if(ydata2(i).gt.ymax2)ymax2=ydata2(i)
 21       continue
c
c   Interpolate the model so that we have y-values at all observed phases.
c
          call spline(xmodel1,ymodel1,Nmodel1,0.0d0,0.0d0,y21)
          call spline(xmodel2,ymodel2,Nmodel2,0.0d0,0.0d0,y22)
c
          do 30 i=1,Ndata1
            call splint(xmodel1,ymodel1,y21,Nmodel1,xdata1(i),qqq)
            yinter1(i)=qqq
 30       continue
c
          do 31 i=1,Ndata2
            call splint(xmodel2,ymodel2,y22,Nmodel2,xdata2(i),qqq)
            yinter2(i)=qqq
 31       continue
c
c   Now find the optimal zero point that will give the lowest chi^2.
c
          call getmean(Ndata1,ydata1,dataave1)
          call getmean(Ndata1,yinter1,rmodelave1)
          call getmean(Ndata2,ydata2,dataave2)
          call getmean(Ndata2,yinter2,rmodelave2)

c          savezero=zero
          zero=(dataave1-rmodelave1)
          step=abs(rmodelave1-dataave1)/1000.0d0
          if(step.lt.0.10d0)step=0.05d0

          chi1=0.0d0
          chi2=0.0d0
          chi3=0.0d0
          chisq1=0.0d0
          chisq2=0.0d0
          small=1.0d35
          zerosmall=1.0d20
c
          do 750 i=1,20
            zero1=zero
            call offset(Ndata1,yinter1,ydummy1,zero1)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi11)
            call offset(Ndata2,yinter2,ydummy2,zero1)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi12)
            chi1=chi11+chi12
            if(chi1.lt.small)then
              small=chi1
              zerosmall=zero1
            endif
            fn=0.0d0
            zero2=zero+step
            call offset(Ndata1,yinter1,ydummy1,zero2)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi21)
            call offset(Ndata2,yinter2,ydummy2,zero2)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi22)
            chi2=chi21+chi22
            if(chi2.lt.small)then
              small=chi2
              zerosmall=zero2
            endif
            diff=chi1-chi2
            if(diff.gt.0.0d0)go to 5061
            step=-step
c  
            csave=chi1
            chi1=chi2
            chi2=csave
            zsave=zero1
            zero1=zero2
            zero2=zsave
c
 5061       fn=fn+1.0d0
c
            zero3=zero2+step
            call offset(Ndata1,yinter1,ydummy1,zero3)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi31)
            call offset(Ndata2,yinter2,ydummy2,zero3)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi32)
            chi3=chi31+chi32
            if(chi3.lt.small)then
              small=chi3
              zerosmall=zero3
            endif
            diff23=chi3-chi2
            if(diff23.lt.0.0d0)then
              chi1=chi2
              chi2=chi3
              zero1=zero2
              zero2=zero3
              zero=zero3
              go to 5061
            endif          
c
c         find the minimum of parabola defined by the last three points
c
            if((chi3-chi2).eq.0.0d0)go to 999
            step=step*(1.0d0/(1.0d0+(chi1-chi2)/(chi3-chi2))+0.5d0)
            zero=zero2-step
            step=step*fn/3.0d0
            call offset(Ndata1,yinter1,ydummy1,zero)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi41)
            call offset(Ndata2,yinter2,ydummy2,zero)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi42)
            chi4=chi41+chi42
            if(chi4.lt.small)then
              small=chi4
              zerosmall=zero
            endif
c
 750      continue  ! loop over grid searches
c
          continue                  ! come here when delta chi is small
c
 999      zero=zerosmall
c
          call offset(Ndata1,yinter1,ydummy1,zero)
          call getchi(Ndata1,ydata1,err1,ydummy1,chisq1)
          call offset(Ndata2,yinter2,ydummy2,zero)
          call getchi(Ndata2,ydata2,err2,ydummy2,chisq2)
c
          do ii=1,Ndata1
            resRV1(ii)=ydata1(ii)-ydummy1(ii)
          enddo
          do ii=1,Ndata2
            resRV2(ii)=ydata2(ii)-ydummy2(ii)
          enddo

c 
          return
c
          end
c
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  NEW BUG August 2, 2001
c
c  Add the period and T0 to the argument list
c
c  UPDATE August 13, 2001
c
c  Make the variable line 200 characters.
c
c
c  UPDATE January 15, 2002
c
c  Add alb1 and alb2 to the list of variables.
c
c  UPDATE November 6, 2002
c
c  Add dwavex and dwavey (limb darkening coefficients) to the 
c  argument list of newwritevar.
c  Dimension them as dwavex(8,2), dwavey(8,2).
c
c  July 29, 2005
c
c  Add ecosw and temprat to the list.
c
          subroutine newwritevar(Nvarmax,svar,var,fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,
     @       rLx,separ,isvel1,gamma1,isvel2,gamma2,t3,g3,SA3,ecc,argper,
     @       pshift,spot1parm,spot2parm,spotdparm,period,T0,chisq,nnn,
     @       line,alb1,alb2,dwavex,dwavey,primmass,primK,primrad,ratrad,
     @       frac1,frac2,ecosw,temprat,bigI,bigbeta,powercoeff,density,
     @       sw5,Tconj,beam1,beam2,contam,ocose,osine,isw29,tertperiod,
     @       tertt0,tertecos,tertesin,tertincl,tertOmega,tertQ,Tgrav1,
     @       Tgrav2,tertconj,omegadot,contamS0,contamS1,contamS2,
     @       contamS3,
     @       P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @       P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @       P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @       P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @       P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @       P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @       P6ratrad,
     @       P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
     @       P7ratrad,
     @       P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @       P8ratrad,gamma3,isvel3,sw72,sw73)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list above.
c
c
c   November 15, 1999
c
c   This routine will determine which variables need to be printed on the
c   screen based on the string codes in svar(1:Nvarmax).  The output
c   printed on the screen will be compact.
c
          implicit double precision (a-h,o-z)
c
c   UPDATE October 28, 2002
c
c   Make the length of line 300 (was 200)
c
c
          character*1700 line       
          character*40 svar(Nvarmax)
          character*10 string10
          character*11 string11
          character*12 string12
          character*13 string13
          character*14 string14
          character*15 string15
          character*16 string16
          character*17 string17
          character*18 string18
c          character*19 string19
c          character*22 string22
c          character*23 string23
          character*7 string7
          character*8 string8
c          character*5 string5
c          character*6 string6
          character*9 string9

          dimension var(Nvarmax)
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          dimension dwavex(8,3),dwavey(8,3)
          dimension powercoeff(8,9)
c
          line=' '
          iline=0
c          iwrite=0
          ilength=1
c
          ilength=ilength+7
          write(string7,2000)nnn
          K=lnblnk(line)
          line=line(1:K)//string7
c
          ilength=ilength+7
c
c   NEW_BUG July 5, 2001
c
c   Check that the value of chi^ does not exceed 99,999,999
c
          tempchi=chisq
          if(chisq.gt.999999999999.9d0)tempchi=999999999999.9d0
          write(string18,2001)tempchi
c
c          write(string12,2001)chisq
c
          K=lnblnk(line)
          line=line(1:K)//string18
          ilength=ilength+18
c
 2000     format(' ',i6)
 2001     format(' ',f17.4)

          iargper=0

          do 10 i=1,Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
            if((kk.eq.450).and.(isw29.gt.0))then   !ecos(omega) oc
              iline=iline+1
c              iwrite=0
              ilength=ilength+15
              write(string15,7831)ocose
              K=lnblnk(line)
              line=line(1:K)//string15
            endif

            if((kk.eq.466).and.(isw29.gt.0))then   !esin(omega) os
              iline=iline+1
c              iwrite=0
              ilength=ilength+15
              write(string15,7831)osine
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
c   tidal apsidal constants
c
            if(kk.eq.1016)then   !rk1 tag a1
              iline=iline+1
              ilength=ilength+9
              write(string9,4444)sw72
              K=lnblnk(line)
              line=line(1:K)//string9
            endif
c
            if(kk.eq.1017)then   !rk2 tag a2
              iline=iline+1
              ilength=ilength+9
              write(string9,4444)sw73
              K=lnblnk(line)
              line=line(1:K)//string9
            endif

c
c   Third body stuff here
c
            if(kk.eq.627)then   !tertperiod tag tt
              iline=iline+1
              ilength=ilength+18
              write(string18,700)tertPeriod
              K=lnblnk(line)
              line=line(1:K)//string18
            endif
c
            if(kk.eq.659)then   !P2period tag ut
              iline=iline+1
              ilength=ilength+18
              write(string18,700)P2period
              K=lnblnk(line)
              line=line(1:K)//string18
            endif
c
            if(kk.eq.691)then   !P3period tag vt
              iline=iline+1
              ilength=ilength+18
              write(string18,700)P3period
              K=lnblnk(line)
              line=line(1:K)//string18
            endif
c
            if(kk.eq.723)then   !P4period tag wt
              iline=iline+1
              ilength=ilength+18
              write(string18,700)P4period
              K=lnblnk(line)
              line=line(1:K)//string18
            endif
c
            if(kk.eq.755)then   !P5period tag xt
              iline=iline+1
              ilength=ilength+18
              write(string18,700)P5period
              K=lnblnk(line)
              line=line(1:K)//string18
            endif
c
            if(kk.eq.595)then   !P6period tag st
              iline=iline+1
              ilength=ilength+18
              write(string18,700)P6period
              K=lnblnk(line)
              line=line(1:K)//string18
            endif
c
            if(kk.eq.243)then   !P7period tag ht
              iline=iline+1
              ilength=ilength+18
              write(string18,700)P7period
              K=lnblnk(line)
              line=line(1:K)//string18
            endif
c
            if(kk.eq.339)then   !P7period tag kt
              iline=iline+1
              ilength=ilength+18
              write(string18,700)P8period
              K=lnblnk(line)
              line=line(1:K)//string18
            endif
c
            if(kk.eq.628)then   !tertT0 tag tu
              iline=iline+1
              ilength=ilength+14
              write(string14,701)tertT0
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.660)then   !P2T0 tag uu
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P2T0
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.692)then   !P3T0 tag vu
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P3T0
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.724)then   !P4T0 tag wu
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P4T0
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.756)then   !P5T0 tag xu
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P5T0
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.596)then   !P6T0 tag xu
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P6T0
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.244)then   !P7T0 tag hu
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P7T0
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.340)then   !P8T0 tag ku
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P8T0
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.617)then   !tertconj tag tj
              iline=iline+1
              ilength=ilength+14
              write(string14,701)tertconj
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.649)then   !P2tconj tag uj
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P2tconj
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.681)then   !P3tconj tag vj
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P3tconj
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.713)then   !P4tconj tag wj
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P4tconj
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.745)then   !P5tconj tag xj
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P5tconj
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.585)then   !P6tconj tag sj
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P6tconj
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.233)then   !P7tconj tag hj
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P7tconj
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.329)then   !P8tconj tag kj
              iline=iline+1
              ilength=ilength+14
              write(string14,701)P8tconj
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.629)then   !tertecos tag tv
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)tertecos
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.661)then   !P2ecos tag uv
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P2ecos
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.693)then   !P3ecos tag vv
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P3ecos
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.725)then   !P4ecos tag wv
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P4ecos
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.757)then   !P5ecos tag xv
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P5ecos
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.597)then   !P6ecos tag sv
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P6ecos
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.245)then   !P7ecos tag hv
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P7ecos
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.341)then   !P8ecos tag kv
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P8ecos
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.630)then   !tertesin tag tw
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)tertesin
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.662)then   !P2esin tag uw
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P2esin
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.694)then   !P3esin tag vw
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P3esin
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.726)then   !P4esin tag ww
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P4esin
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.758)then   !P5esin tag xw
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P5esin
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.598)then   !P6esin tag sw
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P6esin
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.246)then   !P7esin tag hw
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P7esin
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.342)then   !P8esin tag kw
              iline=iline+1
              ilength=ilength+15
              write(string15,7831)P8esin
              K=lnblnk(line)
              line=line(1:K)//string15
            endif
c
            if(kk.eq.631)then   !tertincl tag tx
              iline=iline+1
              ilength=ilength+17
              write(string17,20000)tertincl
              K=lnblnk(line)
              line=line(1:K)//string17
            endif
c
            if(kk.eq.663)then   !P2incl tag ux
              iline=iline+1
              ilength=ilength+17
              write(string17,20000)P2incl
              K=lnblnk(line)
              line=line(1:K)//string17
            endif
c
            if(kk.eq.695)then   !P3incl tag vx
              iline=iline+1
              ilength=ilength+17
              write(string17,20000)P3incl
              K=lnblnk(line)
              line=line(1:K)//string17
            endif
c
            if(kk.eq.727)then   !P4incl tag wx
              iline=iline+1
              ilength=ilength+17
              write(string17,20000)P4incl
              K=lnblnk(line)
              line=line(1:K)//string17
            endif
c
            if(kk.eq.759)then   !P5incl tag xx
              iline=iline+1
              ilength=ilength+17
              write(string17,20000)P5incl
              K=lnblnk(line)
              line=line(1:K)//string17
            endif
c
            if(kk.eq.599)then   !P6incl tag sx
              iline=iline+1
              ilength=ilength+17
              write(string17,20000)P6incl
              K=lnblnk(line)
              line=line(1:K)//string17
            endif
c
            if(kk.eq.247)then   !P7incl tag hx
              iline=iline+1
              ilength=ilength+17
              write(string17,20000)P7incl
              K=lnblnk(line)
              line=line(1:K)//string17
            endif
c
            if(kk.eq.343)then   !P8incl tag kx
              iline=iline+1
              ilength=ilength+17
              write(string17,20000)P8incl
              K=lnblnk(line)
              line=line(1:K)//string17
            endif
c
            if(kk.eq.632)then   !tertOmega tag ty
              iline=iline+1
              ilength=ilength+16
              write(string16,222)tertOmega
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.664)then   !P2Omega tag uy
              iline=iline+1
              ilength=ilength+16
              write(string16,222)P2Omega
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.696)then   !P3Omega tag vy
              iline=iline+1
              ilength=ilength+16
              write(string16,222)P3Omega
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.728)then   !P4Omega tag wy
              iline=iline+1
              ilength=ilength+16
              write(string16,222)P4Omega
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.760)then   !P5Omega tag xy
              iline=iline+1
              ilength=ilength+16
              write(string16,222)P5Omega
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.600)then   !P6Omega tag sy
              iline=iline+1
              ilength=ilength+16
              write(string16,222)P6Omega
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.248)then   !P7Omega tag hy
              iline=iline+1
              ilength=ilength+16
              write(string16,222)P7Omega
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.344)then   !P8Omega tag ky
              iline=iline+1
              ilength=ilength+16
              write(string16,222)P8Omega
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.633)then   !tertQ tag tz
              iline=iline+1
              ilength=ilength+16
              write(string16,2201)tertQ
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.665)then   !P2Q tag uz
              iline=iline+1
              ilength=ilength+16
              write(string16,2201)P2Q
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.697)then   !P3Q tag vz
              iline=iline+1
              ilength=ilength+16
              write(string16,2201)P3Q
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.729)then   !P4Q tag wz
              iline=iline+1
              ilength=ilength+16
              write(string16,2201)P4Q
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.761)then   !P5Q tag xz
              iline=iline+1
              ilength=ilength+16
              write(string16,2201)P5Q
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.601)then   !P6Q tag sz
              iline=iline+1
              ilength=ilength+16
              write(string16,2201)P6Q
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.249)then   !P7Q tag hz
              iline=iline+1
              ilength=ilength+16
              write(string16,2201)P7Q
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.345)then   !P8Q tag kz
              iline=iline+1
              ilength=ilength+16
              write(string16,2201)P8Q
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.641)then   !P2ratrad tag ub
              iline=iline+1
              ilength=ilength+12
              write(string12,776)P2ratrad   
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.673)then   !P3ratrad tag vb
              iline=iline+1
              ilength=ilength+12
              write(string12,776)P3ratrad   
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.705)then   !P4ratrad tag wb
              iline=iline+1
              ilength=ilength+12
              write(string12,776)P4ratrad   
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.737)then   !P5ratrad tag xb
              iline=iline+1
              ilength=ilength+12
              write(string12,776)P5ratrad   
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.577)then   !P6ratrad tag sb
              iline=iline+1
              ilength=ilength+12
              write(string12,776)P6ratrad   
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.225)then   !P7ratrad tag hb
              iline=iline+1
              ilength=ilength+12
              write(string12,776)P7ratrad   
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.321)then   !P8ratrad tag kb
              iline=iline+1
              ilength=ilength+12
              write(string12,776)P8ratrad   
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c

            if(kk.eq.78)then   !use string co for contam
              iline=iline+1
c              iwrite=0
              ilength=ilength+9
              write(string9,783)contam
              K=lnblnk(line)
              line=line(1:K)//string9
            endif

            if(kk.eq.1591)then   !use string s0
              iline=iline+1
c              iwrite=0
              ilength=ilength+9
              write(string9,783)contamS0
              K=lnblnk(line)
              line=line(1:K)//string9
            endif

            if(kk.eq.1592)then   !use string s1
              iline=iline+1
c              iwrite=0
              ilength=ilength+9
              write(string9,783)contamS1
              K=lnblnk(line)
              line=line(1:K)//string9
            endif

            if(kk.eq.1593)then   !use string s2
              iline=iline+1
c              iwrite=0
              ilength=ilength+9
              write(string9,783)contamS2
              K=lnblnk(line)
              line=line(1:K)//string9
            endif

            if(kk.eq.1594)then   !use string s3
              iline=iline+1
c              iwrite=0
              ilength=ilength+9
              write(string9,783)contamS3
              K=lnblnk(line)
              line=line(1:K)//string9
            endif

            if(kk.eq.1144)then   !use string e1
              iline=iline+1
c              iwrite=0
              ilength=ilength+9
              write(string9,781)beam1
              K=lnblnk(line)
              line=line(1:K)//string9
            endif
c
            if(kk.eq.1145)then   !use string e2
              iline=iline+1
c              iwrite=0
              ilength=ilength+9
              write(string9,781)beam2
              K=lnblnk(line)
              line=line(1:K)//string9
            endif
c
            if(kk.eq.492)then   !use string pm
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,777)primmass
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.497)then   !use string pr
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,778)primrad
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.490)then   !use string pk
              iline=iline+1
c              iwrite=0
              ilength=ilength+13
              write(string13,779)primk
              K=lnblnk(line)
              line=line(1:K)//string13
            endif
c
            if(kk.eq.544)then   !use string ra
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,776)ratrad
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1528)then   !use string q1
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,203)frac1
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.1529)then   !use string q2
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,203)frac2
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.612)then   !temprat
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,778)temprat
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.100)then   !density
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,780)density
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c

 776        format(' ',f11.6)
 777        format(' ',f9.5)
 778        format(' ',f10.5)
 779        format(' ',f12.8)
780         format(' ',f10.7)
781         format(' ',f8.4)
783         format(' ',f8.6)
 7831       format(' ',f14.11)
c
c   UPDATE January 15, 2002
c
c   Here are the blocks for alb1 and alb2
c
c   UPDATE MARCH 20, 2004
c
c   USE string 10 in the two blocks below
c
            if(kk.eq.1368)then   !use string l1
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,202)alb1
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.1369)then   !use string l2
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,202)alb2
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.1208)then   !use string g1
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,202)Tgrav1
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.1209)then   !use string g2
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,202)Tgrav2
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.110)then   !use string do
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,699)omegadot
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c

            if(kk.eq.484)then
              iline=iline+1
c              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 15 (add two decimals in format 700 below).
c
              ilength=ilength+18
              write(string18,700)period
              K=lnblnk(line)
              line=line(1:K)//string18
            endif
c
            if(kk.eq.1623)then
              iline=iline+1
c              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 14 (add two decimals in format 701 below).
c
              ilength=ilength+14
c
c  UPDATE May 21, 2002
c
c  Change (string11,200) to (string12,701)
c
              write(string14,701)T0
              K=lnblnk(line)
c
c   UPDATE May 27, 2002
c
c   Change string11 to string12.
c
              line=line(1:K)//string14
            endif

            if(kk.eq.610)then          !Tconj
              iline=iline+1
c              iwrite=0
              ilength=ilength+14
              write(string14,701)Tconj
              K=lnblnk(line)
              line=line(1:K)//string14
            endif
c
            if(kk.eq.269)then
              iline=iline+1
c              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 9 (add one decimal in format 200 below).
c
c   UPDATE September 20, 2007
c
c   Use length 12 (add three decimals in format 20000 below).
c
              ilength=ilength+17
              write(string17,20000)finc
              K=lnblnk(line)
              line=line(1:K)//string17
            endif
c
c   UPDATE May 8, 2006
c
c   Add bigI and bigbeta
c
            if(kk.eq.8)then     !bigI, string aI
              iline=iline+1
c              iwrite=0
              ilength=ilength+9
              write(string9,200)bigI
              K=lnblnk(line)
              line=line(1:K)//string9
            endif
c
            if(kk.eq.1)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,88201)bigbeta
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.384)then
              iline=iline+1
c              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 12 (add four decimals in format 201 below).
c
              ilength=ilength+12
              write(string12,201)Q
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1176)then
              iline=iline+1
c              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 10 (add twp decimals in format 202 below).
c
              ilength=ilength+10
              write(string10,202)fill1
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.1177)then
              iline=iline+1
c              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 10 (add twp decimals in format 203 below).
c
              ilength=ilength+10
              write(string10,203)fill2
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.552)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+8
              write(string8,204)rinner
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.1464)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,205)omega1
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.1465)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,206)omega2
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.558)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+8
              write(string8,207)router
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.1626)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,217)t3  
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.1210)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,218)g3  
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.36)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+7
              write(string7,211)betarim
              K=lnblnk(line)
              line=line(1:K)//string7
            endif
c
            if(kk.eq.1624)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,209)Teff1
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.1625)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,210)Teff2
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.744)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+8
              write(string8,212)xi
              K=lnblnk(line)
              line=line(1:K)//string8
            endif
c
            if(kk.eq.611)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,208)tdisk 
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.375)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,213)rLx
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.580)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,214)separ
              K=lnblnk(line)
              line=line(1:K)//string10
            endif
c
            if(kk.eq.576)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+13
              write(string13,219)SA3
              K=lnblnk(line)
              line=line(1:K)//string13
            endif
c
c  UPDATE May 24, 2002
c
c  Make the length 9, as in ilength=ilength+9, and
c  write(string9,223), and //string9
c
            if(kk.eq.498)then
              iline=iline+1
c              iwrite=0
c
c   UPDATE October 20, 2002
c
c   Use length 13 (add four decimals in format 223 below).
c
              ilength=ilength+13
              write(string13,223)pshift
              K=lnblnk(line)
              line=line(1:K)//string13
            endif
c
            if(kk.eq.130)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,221)ecc
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.17)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+16
              write(string16,222)argper
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
            if(kk.eq.1048)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,401)spot1parm(1,1)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1049)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,402)spot1parm(1,2)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1050)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,403)spot1parm(1,3)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1051)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,404)spot1parm(1,4)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1052)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,405)spot1parm(2,1)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1053)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,406)spot1parm(2,2)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1054)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,407)spot1parm(2,3)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1055)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,408)spot1parm(2,4)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1080)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,501)spot2parm(1,1)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1081)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+1
              write(string11,502)spot2parm(1,2)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1082)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+1
              write(string11,503)spot2parm(1,3)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1083)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,504)spot2parm(1,4)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1084)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,505)spot2parm(2,1)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1085)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,506)spot2parm(2,2)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1086)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,507)spot2parm(2,3)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1087)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,508)spot2parm(2,4)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1112)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,601)spotdparm(1,1)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1113)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,602)spotdparm(1,2)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1114)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,603)spotdparm(1,3)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1115)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,604)spotdparm(1,4)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1116)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,605)spotdparm(2,1)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1117)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,606)spotdparm(2,2)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1118)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,607)spotdparm(2,3)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
            if(kk.eq.1119)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+11
              write(string11,608)spotdparm(2,4)
              K=lnblnk(line)
              line=line(1:K)//string11
            endif
c
c
c   UPDATE November 6, 2002
c
c   Add the limb darkening coefficients here.
c
c   x-coefficient, star 1
c
            if(kk.eq.1752)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,801)dwavex(1,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1753)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,802)dwavex(2,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1754)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,803)dwavex(3,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1755)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,804)dwavex(4,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1756)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,805)dwavex(5,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1757)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,806)dwavex(6,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1758)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,807)dwavex(7,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1759)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,808)dwavex(8,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
c   x-coefficient, star 2
c
            if(kk.eq.1816)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,2801)dwavex(1,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1817)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,2802)dwavex(2,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1818)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,2803)dwavex(3,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1819)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,2804)dwavex(4,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1820)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,2805)dwavex(5,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1821)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,2806)dwavex(6,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1822)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,2807)dwavex(7,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1823)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,2808)dwavex(8,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
c   y-coefficient, star 1
c
            if(kk.eq.1784)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,1801)dwavey(1,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1785)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,1802)dwavey(2,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1786)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,1803)dwavey(3,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1787)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,1804)dwavey(4,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1788)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,1805)dwavey(5,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1789)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,1806)dwavey(6,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1790)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,1807)dwavey(7,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1791)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,1808)dwavey(8,1)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
c   y-coefficient, star 2
c
            if(kk.eq.1720)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3801)dwavey(1,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1721)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3802)dwavey(2,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1722)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3803)dwavey(3,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1723)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3804)dwavey(4,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1724)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3805)dwavey(5,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1725)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3806)dwavey(6,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1726)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3807)dwavey(7,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1727)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3808)dwavey(8,2)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.111)then   ! ecosw, string dphi
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,202)ecosw
              K=lnblnk(line)
              line=line(1:K)//string10
              iargper=100
            endif
c
c   x-coefficient, star 3
c
            if(kk.eq.1400)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,801)dwavex(1,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1401)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,802)dwavex(2,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1402)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,803)dwavex(3,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1403)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,804)dwavex(4,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1404)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,805)dwavex(5,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1405)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,806)dwavex(6,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1406)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,807)dwavex(7,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1407)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,808)dwavex(8,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
c   y-coefficient, star 3
c
            if(kk.eq.1432)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3801)dwavey(1,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1433)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3802)dwavey(2,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1434)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3803)dwavey(3,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1435)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3804)dwavey(4,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1436)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3805)dwavey(5,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1437)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3806)dwavey(6,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1438)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3807)dwavey(7,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.1439)then
              iline=iline+1
c              iwrite=0
              ilength=ilength+12
              write(string12,3808)dwavey(8,3)
              K=lnblnk(line)
              line=line(1:K)//string12
            endif
c
            if(kk.eq.111)then   ! ecosw, string dphi
              iline=iline+1
c              iwrite=0
              ilength=ilength+10
              write(string10,202)ecosw
              K=lnblnk(line)
              line=line(1:K)//string10
              iargper=100
            endif

 10       continue
c
c  July 29, 2005
c
c  If the variable ecosw is assigned, then also print out the argument
c  of periastron
c 
c
            if(iargper.gt.0)then
              call getom(ecc,ecosw,argper)
              iline=iline+1
c              iwrite=0
              ilength=ilength+16
              write(string16,222)argper
              K=lnblnk(line)
              line=line(1:K)//string16
            endif
c
c   UPDATE AUGUST 10, 2005
c
c   Make the gammas string10
c
          if(isvel1.gt.0)then
            iline=iline+1
c            iwrite=0
            ilength=ilength+10
            write(string10,215)gamma1
            K=lnblnk(line)
            line=line(1:K)//string10
          endif
c
          if(isvel2.gt.0)then
            iline=iline+1
c            iwrite=0
            ilength=ilength+10
            write(string10,215)gamma2
            K=lnblnk(line)
            line=line(1:K)//string10
          endif
c
          if(isvel3.gt.0)then
            iline=iline+1
c            iwrite=0
            ilength=ilength+10
            write(string10,215)gamma3
            K=lnblnk(line)
            line=line(1:K)//string10
          endif
c
c   UPDATE November 6, 2008
c
c   if asini (sw5) > 0, record it here
c
          if(sw5.gt.0.0d0)then
            iline=iline+1
c            iwrite=0
            ilength=ilength+17
            write(string17,20000)sw5
            K=lnblnk(line)
            line=line(1:K)//string17
          endif
c
c    UPDATE November 5, 2008
c   
c    If ecc > 0, record the times of inferior and superior
c    conjunction.
c
          if(ecc.gt.0.0d0)then
            call getcontimes(finc,period,ecc,argper,T0,tconj1,tconj2)
            iline=iline+1
c            iwrite=0
            ilength=ilength+14
            write(string14,701)tconj1
            K=lnblnk(line)
            line=line(1:K)//string14
            iline=iline+1
c            iwrite=0
            ilength=ilength+14
            write(string14,701)tconj2
            K=lnblnk(line)
            line=line(1:K)//string14
          endif
c
c  
          if(isw29.gt.0)then
            iline=iline+1
c            iwrite=0
            ilength=ilength+11
            write(string11,221)ecc
            K=lnblnk(line)
            line=line(1:K)//string11
            iline=iline+1
c            iwrite=0
            ilength=ilength+16
            write(string16,222)argper
            K=lnblnk(line)
            line=line(1:K)//string16
          endif


c
c
c   UPDATE October 20, 2002
c
c   Use length 9 (add one decimal in format 200 below).
c
 200      format(' ',f8.5)
20000     format(' ',f16.11)
c
c   UPDATE October 20, 2002
c
c   Use length 12 (add four decimals in format 201 below).
c
 201      format(' ',f11.8)
 2201     format(' ',f15.7)
c
c   UPDATE October 20, 2002
c
c   Use length 10 (add twp decimals in format 202 below).
c
 202      format(' ',f9.7)
c
c   UPDATE October 20, 2002
c
c   Use length 10 (add twp decimals in format 203 below).
c
 203      format(' ',f9.7)
 204      format(' ',f7.5)
 205      format(' ',f9.6)
 206      format(' ',f9.6)
 207      format(' ',f7.5)
 208      format(' ',f9.2)
 209      format(' ',f9.2)
 210      format(' ',f9.2)
 211      format(' ',f6.3)
 212      format(' ',f7.4)
 213      format(' ',1pe10.4)
 214      format(' ',f9.4)
 215      format(' ',f9.4)
 217      format(' ',f9.2)
 218      format(' ',f9.2)
 219      format(' ',f12.7)
 221      format(' ',f10.8)
 222      format(' ',f15.11)
c
c  UPDATE May 24, 2002
c
c  Add this format statement
c
c
c   UPDATE October 20, 2002
c
c   Use length 13 (add four decimals in format 223 below).
c
 223      format(' ',f12.9)
c
c   RVG BUG ALERT   May 8, 2001
c
c   Add these statements for spots
c
 401      format(' ',f10.7)
 402      format(' ',f10.5)
 403      format(' ',f10.5)
 404      format(' ',f10.7)
 405      format(' ',f10.7)
 406      format(' ',f10.5)
 407      format(' ',f10.5)
 408      format(' ',f10.7)
 501      format(' ',f10.7)
 502      format(' ',f10.5)
 503      format(' ',f10.5)
 504      format(' ',f10.7)
 505      format(' ',f10.7)
 506      format(' ',f10.5)
 507      format(' ',f10.5)
 508      format(' ',f10.7)
 601      format(' ',f10.7)
 602      format(' ',f10.5)
 603      format(' ',f10.5)
 604      format(' ',f10.7)
 605      format(' ',f10.7)
 606      format(' ',f10.5)
 607      format(' ',f10.5)
 608      format(' ',f10.7)
c
 699      format(' ',f11.8)
 700      format(' ',f17.12)
c
c   UPDATE May 27, 2002
c
c   Change the format from f10.4 to f11.5.
c
c   UPDATE October 20, 2002
c
c   Use length 14 (add two decimals in format 701 below).
c
 701      format(' ',f13.7)
c
c
c  UPDATE November 6, 2002
c
c  Here are the format statements for the limb darkening coefficients.
c
c
 801      format(' ',f11.8)
 802      format(' ',f11.8)
 803      format(' ',f11.8)
 804      format(' ',f11.8)
 805      format(' ',f11.8)
 806      format(' ',f11.8)
 807      format(' ',f11.8)
 808      format(' ',f11.8)
c
 1801     format(' ',f11.8)
 1802     format(' ',f11.8)
 1803     format(' ',f11.8)
 1804     format(' ',f11.8)
 1805     format(' ',f11.8)
 1806     format(' ',f11.8)
 1807     format(' ',f11.8)
 1808     format(' ',f11.8)
c
 2801     format(' ',f11.8)
 2802     format(' ',f11.8)
 2803     format(' ',f11.8)
 2804     format(' ',f11.8)
 2805     format(' ',f11.8)
 2806     format(' ',f11.8)
 2807     format(' ',f11.8)
 2808     format(' ',f11.8)
c
 3801     format(' ',f11.8)
 3802     format(' ',f11.8)
 3803     format(' ',f11.8)
 3804     format(' ',f11.8)
 3805     format(' ',f11.8)
 3806     format(' ',f11.8)
 3807     format(' ',f11.8)
 3808     format(' ',f11.8)
c
88201     format(' ',f9.5)
c
 4444     format(' ',f8.6)
          return
          end
c
c   @##!!!%(@^%#%*(&(#@#*%^@^$##
c
c July 24, 2001
c
c  Here is a new subroutine that will compute the chi^2 of the disk fraction.
c
         subroutine diskchi(Nobv,sobv,obv,eobv,ochidisk,ochi,compfracs)
c
          implicit double precision (a-h,o-z)
c
          dimension obv(18),eobv(18)
c
c   UPDATE October 27, 2008
c
c   Make ochidisk have dimension 8.  This will allow the user to fit
c   for the disk fraction in more than 1 band at a time.  Also, put
c   the array compfracs into the argument list above
c
          dimension ochidisk(8),compfracs(8,3) 
          character*40 sobv(18)
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
         do 1 i=1,8
            ochidisk(i)=0.0d0
1        continue

          do 1001 i=1,Nobv
c
c  UPDATE NOVEMBER 23, 2001
c
c  make the argument sobv(i)(1:2)
c
            icn=icnvrt(sobv(i)(1:2))
c
c  UPDATE May 27, 2002
c
c  If the variable rmed > 1, then define ochidisk to be the absolute
c  deviation, rather than the normal chi^2.
c
c  UPDATE October 27, 2008
c
c  use the tags d1 to d8 to fit for the disk fraction in filter 1, filter 2,
c  etc.  String d1=1112, d2=1113, ... d9=1119
c
            if(icn.eq.1112)then
              frac=compfracs(1,2)
              if(obv(i).lt.0.0d0)then
                if(frac.gt.dabs(obv(i)))then
                  ochidisk(1)=0.0d0
                else
                  ochidisk(1)=1.d10
                endif
              else
                if(rmed.ge.1.0d0)then
                  ochidisk(1)=dabs((obv(i)-frac)/(eobv(i)))
                else
                  ochidisk(1)=(obv(i)-frac)*(obv(i)-frac)
     $              /(eobv(i)*eobv(i))
                endif
              endif
              ochi=ochi+ochidisk(1)
            endif
c
            if(icn.eq.1113)then
              frac=compfracs(2,2)
              if(obv(i).lt.0.0d0)then
                if(frac.gt.dabs(obv(i)))then
                  ochidisk(2)=0.0d0
                else
                  ochidisk(2)=1.d10
                endif
              else
                if(rmed.ge.1.0d0)then
                  ochidisk(2)=dabs((obv(i)-frac)/(eobv(i)))
                else
                  ochidisk(2)=(obv(i)-frac)*(obv(i)-frac)
     $              /(eobv(i)*eobv(i))
                endif
              endif
              ochi=ochi+ochidisk(2)
            endif
c
            if(icn.eq.1114)then
              frac=compfracs(3,2)
              if(obv(i).lt.0.0d0)then
                if(frac.gt.dabs(obv(i)))then
                  ochidisk(3)=0.0d0
                else
                  ochidisk(3)=1.d10
                endif
              else
                if(rmed.ge.1.0d0)then
                  ochidisk(3)=dabs((obv(i)-frac)/(eobv(i)))
                else
                  ochidisk(3)=(obv(i)-frac)*(obv(i)-frac)
     $              /(eobv(i)*eobv(i))
                endif
              endif
              ochi=ochi+ochidisk(3)
            endif
c
            if(icn.eq.1115)then
              frac=compfracs(4,2)
              if(obv(i).lt.0.0d0)then
                if(frac.gt.dabs(obv(i)))then
                  ochidisk(4)=0.0d0
                else
                  ochidisk(4)=1.d10
                endif
              else
                if(rmed.ge.1.0d0)then
                  ochidisk(4)=dabs((obv(i)-frac)/(eobv(i)))
                else
                  ochidisk(4)=(obv(i)-frac)*(obv(i)-frac)
     $              /(eobv(i)*eobv(i))
                endif
              endif
              ochi=ochi+ochidisk(4)
            endif
c
            if(icn.eq.1116)then
              frac=compfracs(5,2)
              if(obv(i).lt.0.0d0)then
                if(frac.gt.dabs(obv(i)))then
                  ochidisk(5)=0.0d0
                else
                  ochidisk(5)=1.d10
                endif
              else
                if(rmed.ge.1.0d0)then
                  ochidisk(5)=dabs((obv(i)-frac)/(eobv(i)))
                else
                  ochidisk(5)=(obv(i)-frac)*(obv(i)-frac)
     $              /(eobv(i)*eobv(i))
                endif
              endif
              ochi=ochi+ochidisk(5)
            endif
c
            if(icn.eq.1117)then
              frac=compfracs(6,2)
              if(obv(i).lt.0.0d0)then
                if(frac.gt.dabs(obv(i)))then
                  ochidisk(6)=0.0d0
                else
                  ochidisk(6)=1.d10
                endif
              else
                if(rmed.ge.1.0d0)then
                  ochidisk(6)=dabs((obv(i)-frac)/(eobv(i)))
                else
                  ochidisk(6)=(obv(i)-frac)*(obv(i)-frac)
     $              /(eobv(i)*eobv(i))
                endif
              endif
              ochi=ochi+ochidisk(6)
            endif
c
            if(icn.eq.1118)then
              frac=compfracs(7,2)
              if(obv(i).lt.0.0d0)then
                if(frac.gt.dabs(obv(i)))then
                  ochidisk(7)=0.0d0
                else
                  ochidisk(7)=1.d10
                endif
              else
                if(rmed.ge.1.0d0)then
                  ochidisk(7)=dabs((obv(i)-frac)/(eobv(i)))
                else
                  ochidisk(7)=(obv(i)-frac)*(obv(i)-frac)
     $              /(eobv(i)*eobv(i))
                endif
              endif
              ochi=ochi+ochidisk(7)
            endif
c
            if(icn.eq.1119)then
              frac=compfracs(8,2)
              if(obv(i).lt.0.0d0)then
                if(frac.gt.dabs(obv(i)))then
                  ochidisk(8)=0.0d0
                else
                  ochidisk(8)=1.d10
                endif
              else
                if(rmed.ge.1.0d0)then
                  ochidisk(8)=dabs((obv(i)-frac)/(eobv(i)))
                else
                  ochidisk(8)=(obv(i)-frac)*(obv(i)-frac)
     $              /(eobv(i)*eobv(i))
                endif
              endif
              ochi=ochi+ochidisk(8)
            endif
c
 1001     continue
c
          return
          end
c
c  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c
c  NEW BUG  July 24, 2001
c
c  Here is a new subroutine that will compute the chi^2 of the disk fraction.
c
          subroutine wdiskobschi(ochidisk,Nobv,sobv)
c
c  July 24, 2000
c
c  This routine will take the observed variables specified in the
c  gridloop.opt file and compute the additional chi^2
c
          implicit double precision (a-h,o-z)
c
c  UPDATE September 21, 2008
c
c  Make the dimension of obv, eobv, sobv to 11
c
          dimension ochidisk(8)
          character*40 sobv(18)
c
          character*33 section
          character*80 line
          character*33 outstring   
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
          iwrite=0
          ilength=0
c
          section='1234567890123456'
          section='                '
          line='  -->'
c
c          ochi=0.0d0 
c
c   UPDATE May 27, 2002
c
c   If rmed > 1, then we are doing median fitting.  In that case, write
c   to a different format statement.
c
c
c   UPDATE October 27, 2008
c
c   Update this routine so that one can specify which filter the disk
c   fraction is, and also have the ability to specify two or more disk
c   fractions.
c
           do 99 ii=1,Nobv
             icn=icnvrt(sobv(ii)(1:2))

             if(icn.eq.1112)then   ! string d1
               if(rmed.ge.1.0d0)then
                 call chistring(' med_Ufrac',ochidisk(1),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,21)ochidisk(1)
               else
                 call chistring(' chi_Ufrac',ochidisk(1),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,11)ochidisk(1)
               endif
               K=lnblnk(line)
               line=line(1:K)//section
               if(ilength.ge.55)then
                 write(*,200)line
                 iwrite=100
                 ilength=0
                 line='  -->'
               endif
c
             endif
c
             if(icn.eq.1113)then   ! string d2
               if(rmed.ge.1.0d0)then
                 call chistring(' med_Bfrac',ochidisk(2),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,22)ochidisk(2)
               else
                 call chistring(' chi_Bfrac',ochidisk(2),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,12)ochidisk(2)
               endif
               K=lnblnk(line)
               line=line(1:K)//section
               if(ilength.ge.55)then
                 write(*,200)line
                 iwrite=100
                 ilength=0
                 line='  -->'
               endif
c
             endif
c
             if(icn.eq.1114)then   ! string d3
               if(rmed.ge.1.0d0)then
                 call chistring(' med_Vfrac',ochidisk(3),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,23)ochidisk(3)
               else
                 call chistring(' chi_Vfrac',ochidisk(3),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,13)ochidisk(3)
               endif
               K=lnblnk(line)
               line=line(1:K)//section
               if(ilength.ge.55)then
                 write(*,200)line
                 iwrite=100
                 ilength=0
                 line='  -->'
               endif
c
             endif
c
             if(icn.eq.1115)then   ! string d4
               if(rmed.ge.1.0d0)then
                 call chistring(' med_Rfrac',ochidisk(4),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,24)ochidisk(4)
               else
                 call chistring(' chi_Rfrac',ochidisk(4),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,14)ochidisk(4)
               endif
               K=lnblnk(line)
               line=line(1:K)//section
               if(ilength.ge.55)then
                 write(*,200)line
                 iwrite=100
                 ilength=0
                 line='  -->'
               endif
c
             endif
c
             if(icn.eq.1116)then   ! string d5
               if(rmed.ge.1.0d0)then
                 call chistring(' med_Ifrac',ochidisk(5),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,25)ochidisk(5)
               else
                 call chistring(' chi_Ifrac',ochidisk(5),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,15)ochidisk(5)
               endif
               K=lnblnk(line)
               line=line(1:K)//section
               if(ilength.ge.55)then
                 write(*,200)line
                 iwrite=100
                 ilength=0
                 line='  -->'
               endif
c
             endif
c
             if(icn.eq.1117)then   ! string d6
               if(rmed.ge.1.0d0)then
                 call chistring(' med_Jfrac',ochidisk(6),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,26)ochidisk(6)
               else
                 call chistring(' chi_Jfrac',ochidisk(6),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,16)ochidisk(6)
               endif
               K=lnblnk(line)
               line=line(1:K)//section
               if(ilength.ge.55)then
                 write(*,200)line
                 iwrite=100
                 ilength=0
                 line='  -->'
               endif
c
             endif
c
             if(icn.eq.1118)then   ! string d7
               if(rmed.ge.1.0d0)then
                 call chistring(' med_Hfrac',ochidisk(7),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,27)ochidisk(7)
               else
                 call chistring(' chi_Hfrac',ochidisk(7),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,17)ochidisk(7)
               endif
               K=lnblnk(line)
               line=line(1:K)//section
               if(ilength.ge.55)then
                 write(*,200)line
                 iwrite=100
                 ilength=0
                 line='  -->'
               endif
c
             endif
c
             if(icn.eq.1119)then   ! string d8
               if(rmed.ge.1.0d0)then
                 call chistring(' med_Kfrac',ochidisk(8),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,28)ochidisk(8)
               else
                 call chistring(' chi_Kfrac',ochidisk(8),outstring,lll)
                 ilength=ilength+lll
                 section=outstring
                 iwrite=0
c                 write(section,18)ochidisk(8)
               endif
c
               K=lnblnk(line)
               line=line(1:K)//section
               if(ilength.ge.55)then
                 write(*,200)line
                 iwrite=100
                 ilength=0
                 line='  -->'
               endif
c
             endif
c

99        continue
c
          if(iwrite.eq.0)write(*,200)line
c
c   UPDATE May 27, 2002
c
c   Add this new format statement for median fitting.
c
 200      format(a80)

          return
          end
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c  UPDATE JUNE 19, 2008
c
c  Add ikeep, ecc, and argper as arguments.
c
c   NEW BUG August 2, 2001
c
c   These new routines will allow for the fitting of period and T0
c
          subroutine phasefold(Ndatamax,isavN,savx,savy,saverr,
     @      Nbin,Ndata,xdata,ydata,errdata,period,T0,ikeep,ecc,argper,
     @      itime)
c
c   This routine will take data arrays in savex,savey,saverr 
c   (time, mag, magerr), and fold them into the arrays
c   xdata,ydata,errdata.  The latter arrays are the ones normally
c   used in the optimization programs.
c 
c          implicit double precision(a-h,o-z)
c
          implicit none
c
          integer Ndatamax,isavN,Nbin,Ndata,ikeep,itime,i
c
          real*8 savx,savy,saverr,xdata,ydata,errdata,period,T0,ecc,argper
          real*8 pie,eshift,pconj,pconj2,argrad,trc,htrc,ecan,xmc,phper
          real*8 phper2,pshift,tshift,t1

          dimension savx(Ndatamax),savy(Ndatamax),saverr(Ndatamax)
          dimension xdata(Ndatamax),ydata(Ndatamax),errdata(Ndatamax)
c
          parameter(pie=3.14159265358979323d0)
c
c   UPDATE JUNE 19, 2008
c
c   Figure out the conjunction phases if the orbit is eccentric.
c
          if(itime.eq.2)then
            if(Ndata.gt.1)call sort3(Ndata,xdata,ydata,errdata)
            return
          endif

          eshift=0.0d0
          pconj=pie
          pconj2=0.0d0

          if((ecc.gt.0.0d0).and.(ikeep.gt.0))then
            argrad=argper*pie/180.0d0
c
c  Here is some code adapted from Wilson-Devinney to keep track of
c  phases needed for eccentric orbits.
c
            trc=0.5d0*pie-argrad
 1139       if(trc.lt.0.d0) trc=trc+2.0d0*pie
            if(trc.lt.0.d0) goto 1139
 1140       if(trc.ge.2.0d0*pie) trc=trc-2.0d0*pie
            if(trc.ge.2.0d0*pie) goto 1140
            htrc=0.5d0*trc
            if(dabs(0.5*pie-htrc).lt.7.d-6) goto 11101
            if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 11101
            ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))
            goto 11103
11101       ecan=pie
11103       xmc=ecan-ecc*dsin(ecan)
            if(xmc.lt.0.d0) xmc=xmc+2.0d0*pie
            phper=1.d0-xmc/(2.0d0*pie)
            pconj=(xmc+argrad)/(2.0d0*pie)-0.25d0
c
c   Make sure the conjunction phase is between 0 and 1
c
            if(pconj.gt.1.0d0)pconj=pconj-1.0d0

            trc=0.5d0*pie-argrad+pie
 3139       if(trc.lt.0.d0) trc=trc+2.0d0*pie
            if(trc.lt.0.d0) goto 3139
 3140       if(trc.ge.2.0d0*pie) trc=trc-2.0d0*pie
            if(trc.ge.2.0d0*pie) goto 3140
            htrc=0.5d0*trc
            if(dabs(0.5*pie-htrc).lt.7.d-6) goto 31101
            if(dabs(4.712388980384690d0-htrc).lt.7.d-6) goto 31101
            ecan=2.d0*datan(dsqrt((1.d0-ecc)/(1.d0+ecc))*dtan(htrc))
            goto 31103
31101       ecan=pie
31103       xmc=ecan-ecc*dsin(ecan)
            if(xmc.lt.0.d0) xmc=xmc+2.0d0*pie
            phper2=1.d0-xmc/(2.0d0*pie)
            pconj2=(xmc+argrad)/(2.0d0*pie)-0.25d0
c
c   UPDATE March 14, 2008
c
c   Make sure the conjunction phase is between 0 and 1
c
            if(pconj2.gt.1.0d0)pconj2=pconj2-1.0d0
c
            if(ikeep.eq.1)eshift=phper+pconj-0.5
            if(ikeep.eq.2)eshift=phper2+pconj2
          endif

c   Here is the case where no binning is requested.
c
          pshift=0.0d0
          tshift=0.0d0
          if(ikeep.eq.1)tshift=pshift+eshift-pconj
          if(ikeep.eq.2)tshift=pshift+eshift-pconj2

          if(nbin.eq.0)then
            do 10 i=1,isavN
              t1=(savx(i)-T0)/period+tshift
              xdata(i)=dmod(t1,1.0d0)
              ydata(i)=savy(i)
              errdata(i)=saverr(i)
 10         continue
            Ndata=isavN
c
           do 15 i=1,Ndata
             if(xdata(i).lt.0.0d0)xdata(i)=xdata(i)+1.0d0
             if(xdata(i).gt.1.0d0)xdata(i)=xdata(i)-1.0d0
             if(xdata(i).lt.0.0d0)xdata(i)=xdata(i)+1.0d0
             if(xdata(i).gt.1.0d0)xdata(i)=xdata(i)-1.0d0
             if(xdata(i).lt.0.0d0)xdata(i)=xdata(i)+1.0d0
             if(xdata(i).gt.1.0d0)xdata(i)=xdata(i)-1.0d0
 15        continue
c
c   Sort the data by phase.
c
c   UPDATE April 15, 2002
c
c   Put an if-then clause to cover the case when Ndata=1
c
c
            if(Ndata.gt.1)call sort3(Ndata,xdata,ydata,errdata)
c
c   We are done, so leave.
c
            return
          endif
c

          return
          end
c
c
c
          subroutine kopydata(Ndatamax,Ndata,xdata,ydata,errdata,isavN,
     $       savx,savy,saverr)
c
c   August 2, 2001
c
c   When the period and/or T0 is being adjusted, the data is entered
c   in time units.  We have to save the original arrays because the
c   data can be multiply folded during a single run of a program.
c
          implicit double precision(a-h,o-z)
c
          dimension savx(Ndatamax),savy(Ndatamax),saverr(Ndatamax)
          dimension xdata(Ndatamax),ydata(Ndatamax),errdata(Ndatamax)
c
          isavN=Ndata
c
          do 10 i=1,Ndata
            savx(i)=xdata(i)
            savy(i)=ydata(i)
            saverr(i)=errdata(i)
 10       continue
c
          return
c
          end
c
c  ##########################################################         
c
          subroutine wELCdata(Ndatamax,Ndata,xdata,ydata,err,fileout)
c
c  NEW BUG August 2, 2001
c
c  This routine is new.  When the period and/or T0 is being fitted for,
c  the folded data files will be written.
c
          implicit double precision(a-h,o-z)
c
          dimension xdata(Ndatamax),ydata(Ndatamax),err(Ndatamax)
c
c   UPDATE November 18, 2002
c
c   Change to character*(*) below.
c
          character*(*) fileout

          open(unit=51,file=fileout,status='unknown')
c
          do 10 i=1,Ndata
            write(51,100)xdata(i),ydata(i),err(i)
 10       continue
c
          close(51)
c
 100      format(f21.11,3x,f18.13,3x,f17.13)
c
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c   UPDATE August 13, 2001
c
c   This is a new subroutine, which is the inverse of assignvar.  That it,
c   the routine has values for the parameters (fill1, fill2, omega1, ...)
c   and it will assign them to the var string based on the variables
c   selected in gridloop.opt file.
c
c
c  UPDATE January 15, 2002
c
c  Update the routine to include the albedos of the two stars (alb1, alb2)
c
c
c  UPDATE November 6, 2002
c
c  Add limb darkening coefficients dwavex and dwavey to the
c  argument list of assignvar and varassign.
c
c  UPDATE August 13, 2004
c  
c  Add primmass,primK,primrad,ratrad,frac1,frac2
c
          subroutine varassign(Nvarmax,svar,var,fill1,fill2,omega1,
     $       omega2,Q,finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,
     @       rLx,separ,gamma,t3,g3,sa3,ecc,argper,pshift,spot1parm,
     @       spot2parm,spotdparm,period,T0,alb1,alb2,dwavex,dwavey,
     @       primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @       bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2,contam,
     @       ocose,osine,isw29,tertperiod,tertt0,tertecos,tertesin,
     @       tertincl,tertOmega,tertQ,Tgrav1,Tgrav2,tertconj,omegadot,
     @       contamS0,contamS1,contamS2,contamS3,
     @       P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @       P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @       P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @       P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @       P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @       P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @       P6ratrad,
     @       P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
     @       P7ratrad,
     @       P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @       P8ratrad,sw72,sw73)
c
c  UPDATE May 8, 2006
c
c  Add bigI, bigbeta, powercoeff to list
c
c  UPDATE October 10, 2008
c
c  Add the density to the list
c

          implicit double precision (a-h,o-z)
c
          character*40 svar(Nvarmax)
          dimension var(Nvarmax)
c       
          dimension spot1parm(2,4),spot2parm(2,4),spotdparm(2,4)
          dimension powercoeff(8,9)
c
c
c  UPDATE November 6, 2002
c
c  Dimension the limb darkening variables.
c
          dimension dwavex(8,3),dwavey(8,3)
c
          write(*,*)' '
          do 10 i=1,Nvarmax
            kk=icnvrt(svar(i)(1:2))
c
c  tidal apsidal k-coefficients
c
            if(kk.eq.1016)then  !sw72, tag a1
              var(i)=sw72
            endif

            if(kk.eq.1017)then  !sw73, tag a2
              var(i)=sw73
            endif

c
c  planet 2 parameters
c
            if(kk.eq.649)then  !P2tconj, tag uj
              var(i)=P2tconj
            endif
c
            if(kk.eq.659)then  !P2period, tag ut
              var(i)=P2period
            endif
c
            if(kk.eq.660)then  !P2period, tag uu
              var(i)=P2T0
            endif
c
            if(kk.eq.661)then  !P2ecos, tag uv
              var(i)=P2ecos
            endif
c
            if(kk.eq.662)then  !P2esin, tag uw
              var(i)=P2esin
            endif
c
            if(kk.eq.663)then  !P2incl, tag ux
              var(i)=P2incl
            endif
c
            if(kk.eq.664)then  !P2Omega, tag uy
              var(i)=P2Omega
            endif
c
            if(kk.eq.665)then  !P2Q, tag uz
              var(i)=P2Q
            endif
c
            if(kk.eq.641)then  !P2ratrad, tag ub
              var(i)=P2ratrad
            endif
c
c planet 3 parameters
c
            if(kk.eq.681)then  !P3tconj, tag vj
              var(i)=P3tconj
            endif
c
            if(kk.eq.691)then  !P3period, tag vt
              var(i)=P3period
            endif
c
            if(kk.eq.692)then  !P3period, tag vu
              var(i)=P3T0
            endif
c
            if(kk.eq.693)then  !P3ecos, tag vv
              var(i)=P3ecos
            endif
c
            if(kk.eq.694)then  !P3esin, tag vw
              var(i)=P3esin
            endif
c
            if(kk.eq.695)then  !P3incl, tag vx
              var(i)=P3incl
            endif
c
            if(kk.eq.696)then  !P3Omega, tag vy
              var(i)=P3Omega
            endif
c
            if(kk.eq.697)then  !P3Q, tag vz
              var(i)=P3Q
            endif
c
            if(kk.eq.673)then  !P3ratrad, tag vb
              var(i)=P3ratrad
            endif
c
c  planet 4 parameters
c
            if(kk.eq.713)then  !P4tconj, tag wj
              var(i)=P4tconj
            endif
c
            if(kk.eq.723)then  !P4period, tag wt
              var(i)=P4period
            endif
c
            if(kk.eq.724)then  !P4period, tag wu
              var(i)=P4T0
            endif
c
            if(kk.eq.725)then  !P4ecos, tag wv
              var(i)=P4ecos
            endif
c
            if(kk.eq.726)then  !P4esin, tag ww
              var(i)=P4esin
            endif
c
            if(kk.eq.727)then  !P4incl, tag wx
              var(i)=P4incl
            endif
c
            if(kk.eq.728)then  !P4Omega, tag wy
              var(i)=P4Omega
            endif
c
            if(kk.eq.729)then  !P4Q, tag wz
              var(i)=P4Q
            endif
c
            if(kk.eq.705)then  !P4ratrad, tag wb
              var(i)=P4ratrad
            endif
c
c  planet 5 parameters
c
            if(kk.eq.745)then  !P5tconj, tag xj
              var(i)=P5tconj
            endif
c
            if(kk.eq.755)then  !P5period, tag xt
              var(i)=P5period
            endif
c
            if(kk.eq.756)then  !P5period, tag xu
              var(i)=P5T0
            endif
c
            if(kk.eq.757)then  !P5ecos, tag xv
              var(i)=P5ecos
            endif
c
            if(kk.eq.758)then  !P5esin, tag xw
              var(i)=P5esin
            endif
c
            if(kk.eq.759)then  !P5incl, tag xx
              var(i)=P5incl
            endif
c
            if(kk.eq.760)then  !P5Omega, tag xy
              var(i)=P5Omega
            endif
c
            if(kk.eq.761)then  !P5Q, tag xz
              var(i)=P5Q
            endif
c
            if(kk.eq.737)then  !P5ratrad, tag xb
              var(i)=P5ratrad
            endif
c
c  planet 6 parameters
c
            if(kk.eq.585)then  !P6tconj, tag sj
              var(i)=P6tconj
            endif
c
            if(kk.eq.595)then  !P6period, tag st
              var(i)=P6period
            endif
c
            if(kk.eq.596)then  !P6period, tag su
              var(i)=P6T0
            endif
c
            if(kk.eq.597)then  !P6ecos, tag sv
              var(i)=P6ecos
            endif
c
            if(kk.eq.598)then  !P6esin, tag sw
              var(i)=P6esin
            endif
c
            if(kk.eq.599)then  !P6incl, tag sx
              var(i)=P6incl
            endif
c
            if(kk.eq.600)then  !P6Omega, tag sy
              var(i)=P6Omega
            endif
c
            if(kk.eq.601)then  !P6Q, tag sz
              var(i)=P6Q
            endif
c
            if(kk.eq.577)then  !P6ratrad, tag sb
              var(i)=P6ratrad
            endif
c
c  planet 7 parameters
c
            if(kk.eq.233)then  !P7tconj, tag hj
              var(i)=P7tconj
            endif
c
            if(kk.eq.243)then  !P7period, tag ht
              var(i)=P7period
            endif
c
            if(kk.eq.244)then  !P7period, tag hu
              var(i)=P7T0
            endif
c
            if(kk.eq.245)then  !P7ecos, tag hv
              var(i)=P7ecos
            endif
c
            if(kk.eq.246)then  !P7esin, tag hw
              var(i)=P7esin
            endif
c
            if(kk.eq.247)then  !P7incl, tag hx
              var(i)=P7incl
            endif
c
            if(kk.eq.248)then  !P7Omega, tag hy
              var(i)=P7Omega
            endif
c
            if(kk.eq.249)then  !P7Q, tag hz
              var(i)=P7Q
            endif
c
            if(kk.eq.225)then  !P7ratrad, tag hb
              var(i)=P7ratrad
            endif
c
c  planet 8 parameters
c
            if(kk.eq.329)then  !P8tconj, tag kj
              var(i)=P8tconj
            endif
c
            if(kk.eq.339)then  !P8period, tag kt
              var(i)=P8period
            endif
c
            if(kk.eq.340)then  !P8period, tag ku
              var(i)=P8T0
            endif
c
            if(kk.eq.341)then  !P8ecos, tag kv
              var(i)=P8ecos
            endif
c
            if(kk.eq.342)then  !P8esin, tag kw
              var(i)=P8esin
            endif
c
            if(kk.eq.343)then  !P8incl, tag kx
              var(i)=P8incl
            endif
c
            if(kk.eq.344)then  !P8Omega, tag ky
              var(i)=P8Omega
            endif
c
            if(kk.eq.345)then  !P8Q, tag kz
              var(i)=P8Q
            endif
c
            if(kk.eq.42)then  !P8ratrad, tag kb
              var(i)=P8ratrad
            endif
c
c   seasonal contamination
c
            if(kk.eq.1591)then   !contamS0, tag s0
              var(i)=contamS0
            endif

            if(kk.eq.1592)then   !contamS1, tag s1
              var(i)=contamS1
            endif

            if(kk.eq.1593)then   !contamS2, tag s2
              var(i)=contamS2
            endif

            if(kk.eq.1594)then   !contamS3, tag s3
              var(i)=contamS3
            endif

            if(kk.eq.110)then  !omegadot, tag do
              var(i)=omegadot
            endif
c
            if(kk.eq.1208)then  !Tgrav1, tag g1
              var(i)=Tgrav1
            endif
c
            if(kk.eq.1209)then  !Tgrav2, tag g2
              var(i)=Tgrav2
            endif
c

            if((isw29.gt.0).and.(kk.eq.450))then !ecos(omega) oc
              var(i)=ocose
            endif
c
            if((isw29.gt.0).and.(kk.eq.466))then !esin(omega) os
              var(i)=osine
            endif
c
c   November 18, 2012
c
c   Add third body stuff here
c
            if(kk.eq.617)then  !tertperiod  tag tj
              var(i)=tertconj
            endif
c
            if(kk.eq.627)then  !tertperiod  tag tt
              var(i)=tertperiod
            endif
c
            if(kk.eq.628)then  !tertT0  tag tu
              var(i)=tertT0
            endif
c
            if(kk.eq.629)then  !tertecos  tag tv
              var(i)=tertecos
            endif
c
            if(kk.eq.630)then  !tertesin  tag tw
              var(i)=tertesin
            endif
c
            if(kk.eq.631)then  !terteincl  tag tx
              var(i)=tertincl
            endif
c
            if(kk.eq.632)then  !tertOmega  tag ty
              var(i)=tertOmega
            endif
c
            if(kk.eq.633)then  !tertQ  tag tz
              var(i)=tertQ
            endif
c
            if(kk.eq.78)then  !contam, use string co
              var(i)=contam
            endif
c
            if(kk.eq.1144)then  !beam1, use string e1
              var(i)=beam1
            endif
c
            if(kk.eq.1145)then  !Tconj, use string e2
              var(i)=beam2
            endif
c
            if(kk.eq.610)then  !Tconj, use string tc
              var(i)=Tconj
            endif
c
            if(kk.eq.100)then  !density, use string de
              var(i)=density
            endif
c
            if(kk.eq.8)then  !bigI, use string ai
              var(i)=bigI
            endif
c
            if(kk.eq.1)then  !bigbeta, use string ab
              var(i)=bigbeta
            endif
c
            if(kk.eq.111)then  !ecosw, use string dphi
              var(i)=ecosw
            endif
c
            if(kk.eq.612)then  !temprat
              var(i)=temprat
            endif
c
            if(kk.eq.492)then  !primmass, use string pm
              var(i)=primmass
            endif
c
            if(kk.eq.497)then  !primrad, use string pr
              var(i)=primrad
            endif
c
            if(kk.eq.490)then  !primK, use string pk
              var(i)=primK
            endif
c
            if(kk.eq.544)then  !ratrad, use string ra
              var(i)=ratrad
            endif
c
            if(kk.eq.1528)then  !frac1, use string q1
              var(i)=frac1
            endif
c
            if(kk.eq.1529)then  !frac2, use string q2
              var(i)=frac2
            endif
c
c
c  UPDATE January 15, 2002
c
c  Here are the assignments for alb1 and alb2
c
            if(kk.eq.1368)then  !alb1
              var(i)=alb1
            endif
c
            if(kk.eq.1369)then  !alb2
              var(i)=alb2
            endif
c
            if(kk.eq.484)then  !period
              var(i)=period
            endif
c
            if(kk.eq.1623)then ! T0
              var(i)=T0
            endif
c
            if(kk.eq.1112)then   ! temperature factor spot 1 on disk
              var(i)=spotdparm(1,1)
            endif
c
            if(kk.eq.1116)then   ! temperature factor spot 2 on disk
              var(i)=spotdparm(2,1)
            endif
c
            if(kk.eq.1113)then  ! azimuth spot 1 on disk
              var(i)=spotdparm(1,2)
            endif   
c
            if(kk.eq.1117)then  ! azimuth spot 2 on disk
              var(i)=spotdparm(2,2)
            endif   
c
            if(kk.eq.1114)then    ! cutoff radius for spot 1 on disk
              var(i)=spotdparm(1,3)
            endif
c
            if(kk.eq.1118)then    ! cutoff radius for spot 2 on disk
              var(i)=spotdparm(2,3)
            endif
c
            if(kk.eq.1115)then    ! angular width of spot 1 on disk
              var(i)=spotdparm(1,4)
            endif
c
            if(kk.eq.1119)then    ! angular width of spot 2 on disk
              var(i)=spotdparm(2,4)
            endif
c
c
            if(kk.eq.1080)then         ! temperature factor spot 1, star 2
              var(i)=spot2parm(1,1)
            endif
c
            if(kk.eq.1084)then         ! temperature factor spot 2, star 2
              var(i)=spot2parm(2,1)
            endif
c
            if(kk.eq.1048)then         ! temperature factor spot 1, star 1
              var(i)=spot1parm(1,1)
            endif
c
            if(kk.eq.1052)then         ! temperature factor spot 2, star 1
              var(i)=spot1parm(2,1)
            endif
c
            if(kk.eq.1081)then         ! latitude spot 1, star 2
              var(i)=spot2parm(1,2)
            endif
c
            if(kk.eq.1085)then         ! latitude spot 2, star 2
              var(i)=spot2parm(2,2)
            endif
c
            if(kk.eq.1049)then         ! latitude spot 1, star 1
              var(i)=spot1parm(1,2)
            endif
c
            if(kk.eq.1053)then         ! latitude spot 2, star 1
              var(i)=spot1parm(2,2)
            endif
c
            if(kk.eq.1082)then         ! longitude spot 1, star 2
              var(i)=spot2parm(1,3)
            endif
c
            if(kk.eq.1086)then         ! longitude spot 2, star 2
              var(i)=spot2parm(2,3)
            endif
c
            if(kk.eq.1050)then         ! longitude spot 1, star 1
              var(i)=spot1parm(1,3)
            endif
c
            if(kk.eq.1054)then         ! longitude spot 2, star 1
              var(i)=spot1parm(2,3)
            endif
c
            if(kk.eq.1083)then         ! radius spot 1, star 2
              var(i)=spot2parm(1,4)
            endif
c
            if(kk.eq.1087)then         ! radius spot 2, star 2
              var(i)=spot2parm(2,4)
            endif
c
            if(kk.eq.1051)then         ! radius spot 1, star 1
              var(i)=spot1parm(1,4)
            endif
c
            if(kk.eq.1055)then         ! radius spot 2, star 1
              var(i)=spot1parm(2,4)
            endif

            if(kk.eq.498)then
              var(i)=pshift
            endif
c
            if(kk.eq.269)then
              var(i)=finc
            endif
c
            if(kk.eq.384)then
              var(i)=Q
            endif
c
            if(kk.eq.130)then
              var(i)=ecc
            endif
c
            if(kk.eq.17)then
              var(i)=argper
            endif
c
            if(kk.eq.1176)then
              var(i)=fill1
            endif
c
            if(kk.eq.552)then
              var(i)=rinner
            endif
c
            if(kk.eq.1177)then
              var(i)=fill2
            endif
c
            if(kk.eq.1464)then
              var(i)=omega1
            endif
c
            if(kk.eq.1465)then
              var(i)=omega2
            endif
c
            if(kk.eq.558)then
              var(i)=router
            endif
c
            if(kk.eq.611)then
              var(i)=Tdisk
            endif
c
            if(kk.eq.36)then
              var(i)=betarim
            endif
c
            if(kk.eq.1624)then
              var(i)=Teff1
            endif
c
            if(kk.eq.1625)then
              var(i)=Teff2
            endif
c
            if(kk.eq.744)then
              var(i)=xi
            endif
c
            if(kk.eq.375)then
              var(i)=rLx
            endif
c
            if(kk.eq.580)then
              var(i)=separ
            endif
c
            if(kk.eq.192)then
              var(i)=gamma
            endif
c
            if(kk.eq.1626)then
              var(i)=t3
            endif
c
            if(kk.eq.1210)then
              var(i)=g3
            endif
c
            if(kk.eq.576)then
              var(i)=SA3
            endif
c
            if(kk.eq.1752)var(i)=dwavex(1,1)
            if(kk.eq.1753)var(i)=dwavex(2,1)
            if(kk.eq.1754)var(i)=dwavex(3,1)
            if(kk.eq.1755)var(i)=dwavex(4,1)
            if(kk.eq.1756)var(i)=dwavex(5,1)
            if(kk.eq.1757)var(i)=dwavex(6,1)
            if(kk.eq.1758)var(i)=dwavex(7,1)
            if(kk.eq.1759)var(i)=dwavex(8,1)
c
            if(kk.eq.1784)var(i)=dwavey(1,1)
            if(kk.eq.1785)var(i)=dwavey(2,1)
            if(kk.eq.1786)var(i)=dwavey(3,1)
            if(kk.eq.1787)var(i)=dwavey(4,1)
            if(kk.eq.1788)var(i)=dwavey(5,1)
            if(kk.eq.1789)var(i)=dwavey(6,1)
            if(kk.eq.1790)var(i)=dwavey(7,1)
            if(kk.eq.1791)var(i)=dwavey(8,1)
c
            if(kk.eq.1816)var(i)=dwavex(1,2)
            if(kk.eq.1817)var(i)=dwavex(2,2)
            if(kk.eq.1818)var(i)=dwavex(3,2)
            if(kk.eq.1819)var(i)=dwavex(4,2)
            if(kk.eq.1820)var(i)=dwavex(5,2)
            if(kk.eq.1821)var(i)=dwavex(6,2)
            if(kk.eq.1822)var(i)=dwavex(7,2)
            if(kk.eq.1823)var(i)=dwavex(8,2)
c
            if(kk.eq.1720)var(i)=dwavey(1,2)
            if(kk.eq.1721)var(i)=dwavey(2,2)
            if(kk.eq.1722)var(i)=dwavey(3,2)
            if(kk.eq.1723)var(i)=dwavey(4,2)
            if(kk.eq.1724)var(i)=dwavey(5,2)
            if(kk.eq.1725)var(i)=dwavey(6,2)
            if(kk.eq.1726)var(i)=dwavey(7,2)
            if(kk.eq.1727)var(i)=dwavey(8,2)
c
c   body 3 parameters
c
            if(kk.eq.1400)var(i)=dwavex(1,3)   !tag m1
            if(kk.eq.1401)var(i)=dwavex(2,3)
            if(kk.eq.1402)var(i)=dwavex(3,3)
            if(kk.eq.1403)var(i)=dwavex(4,3)
            if(kk.eq.1404)var(i)=dwavex(5,3)
            if(kk.eq.1405)var(i)=dwavex(6,3)
            if(kk.eq.1406)var(i)=dwavex(7,3)
            if(kk.eq.1407)var(i)=dwavex(8,3)    ! tag m8
c
            if(kk.eq.1432)var(i)=dwavey(1,3)   ! tag n1
            if(kk.eq.1433)var(i)=dwavey(2,3)
            if(kk.eq.1434)var(i)=dwavey(3,3)
            if(kk.eq.1435)var(i)=dwavey(4,3)
            if(kk.eq.1436)var(i)=dwavey(5,3)
            if(kk.eq.1437)var(i)=dwavey(6,3)
            if(kk.eq.1438)var(i)=dwavey(7,3)
            if(kk.eq.1439)var(i)=dwavey(8,3)     !tag n8

 10       continue
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&
c
c  UPDATE October 12, 2001
c
c  Here is a new subroutine that will write the chi^2 values in a
c  single line
c
c  UPDATE FEBRUARY 4, 2005
c
c  Add ochilr to the end of the list.
c

          subroutine chiline(iter,chiall,ochidisk,
     @       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %       icnRV1,icnRV2,chiU,chiB,chiV,chiR,chiI,chiJ,chiH,chiK,
     &       chiRV1,chiRV2,Nobv,sobv,obv,eobv,obsparm,ochilr)
c
c   November 15, 1999
c
c   This routine will write chi^2 values to the screen in a compact way.
c
          implicit double precision (a-h,o-z)
c
c   UPDATE September 21, 2008
c
c   make the dimension of obsparm, obv, eobv, sobv 11
c
c   UPDATE October 10, 2008
c
c   Make the dimension of obsparm, obv, eobv, sobv 17
c
c   UPDATE October 27, 2008
c
c   ochidisk is now an array with dimension 8.
c
c   UPDATE March 15, 2011
c
c   Make the dimension of obsparm et al. 18
c
          dimension ochidisk(8)
          dimension obv(19),eobv(19),obsparm(19)
          character*40 sobv(19)
          character*15 section(20)
          character*80 line
c
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
          iline=1
c          iwrite=0
c
          do 10 i=1,10
            section(i)='123456789012345'
            section(i)='               '
 10       continue
          line='  -->'
c
c    UPDATE Oct 28, 2002
c
c    Add d0 to the end of the 999999 strings below.
c
c
          cccc=chiall
          if(cccc.gt.999999999.0d0)cccc=999999999.0d0
          write(section(1),600)cccc
          if(icnU.ne.430)then
            iline=iline+1
            cccc=chiU
            if(cccc.gt.999999999.0d0)cccc=999999999.0d0
            write(section(iline),600)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnB.ne.430)then
            iline=iline+1
            cccc=chiB
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),101)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnV.ne.430)then
            iline=iline+1
            cccc=chiV
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),102)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnR.ne.430)then
            iline=iline+1
            cccc=chiR
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),103)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnI.ne.430)then
            iline=iline+1
            cccc=chiI
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),104)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnJ.ne.430)then
            iline=iline+1
c            iwrite=0
            cccc=chiJ
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),105)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnH.ne.430)then
            iline=iline+1
c            iwrite=0
            cccc=chiH
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),106)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnK.ne.430)then
            iline=iline+1
c            iwrite=0
            cccc=chiK
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),107)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnRV1.ne.430)then
            iline=iline+1
c            iwrite=0
            cccc=chiRV1
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),108)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
          if(icnRV2.ne.430)then
            iline=iline+1
c            iwrite=0
            cccc=chiRV2
            if(cccc.gt.9999999.0d0)cccc=9999999.0d0
            write(section(iline),109)cccc
            K=lnblnk(line)
            line=line(1:K)//section(iline)
          endif
c
          do 1001 i=1,Nobv
c
c  UPDATE NOVEMBER 23, 2001
c
c  make the argument sobv(i)(1:2)
c
            icn=icnvrt(sobv(i)(1:2))
c
c   NEW BUG  July 24, 2001
c
c   Escape if icn=104 (disk fraction requested).
c
c            if(icn.eq.104)then
c              iline=iline+1
c              cccc=ochidisk(1)
c              if(cccc.gt.9999999.d0)cccc=9999999.0d0
c              write(section(iline),100)cccc
c              go to 1001
c            endif
c
            if(icn.eq.1112)then
              iline=iline+1
              cccc=ochidisk(1)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.1113)then
              iline=iline+1
              cccc=ochidisk(2)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.1114)then
              iline=iline+1
              cccc=ochidisk(3)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.1115)then
              iline=iline+1
              cccc=ochidisk(4)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.1116)then
              iline=iline+1
              cccc=ochidisk(5)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.1117)then
              iline=iline+1
              cccc=ochidisk(6)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.1118)then
              iline=iline+1
              cccc=ochidisk(7)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
            if(icn.eq.1119)then
              iline=iline+1
              cccc=ochidisk(8)
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif

c
c   Escape if icn=369 (luminosity ratio requested).
c
            if(icn.eq.369)then
              iline=iline+1
              cccc=ochilr
              if(cccc.gt.9999999.d0)cccc=9999999.0d0
              write(section(iline),100)cccc
              go to 1001
            endif
c
c   END BUG
c
c   NEW BUG August 2, 2001
c 
c   Put an escape if one of the variables is not valid.
c
            index=-1
            if(icn.eq.1400)index=1   !m1
            if(icn.eq.1401)index=5    !m2
            if(icn.eq.1560)index=2     !r1
            if(icn.eq.1561)index=6    !r2
            if(icn.eq.1208)index=3    !g1
c  
            if(icn.eq.1336)index=10   !k1
            if(icn.eq.1337)index=11    !k2
c
c
c   UPDATE October 10, 2008
c
c   Here are the strings for incl, mass, ecc, arg, t1, t2
c
            if(icn.eq.269)index=12
            if(icn.eq.384)index=13
            if(icn.eq.130)index=14
            if(icn.eq.17)index=15
            if(icn.eq.1624)index=16
            if(icn.eq.1625)index=17
c
c   T3
c
            if(icn.eq.1626)index=19
c
c   UPDATE March 15, 2011
c
c   Here is the string for sum of fractional radii (sum)
c
            if(icn.eq.596)index=18
c
c   UPDATE May 22, 2002
c
c   Change the 208 to 209 in the argument of if.
c
            if(icn.eq.1209)index=7
            if(icn.eq.1688)index=4    !v1
            if(icn.eq.1689)index=8    !v2
c
c   UPDATE September 11, 2001
c
c   Add this if block
c
            if(icn.eq.740)index=9     !xe
c
c   Here is the escape.
c
            if(index.lt.0)go to 1001
c
c   UPDATE June 3, 2002
c
c   There are two changes.  First, if the variable rmed > 1,
c   then we are doing median fitting.  chisq becomes the absolute
c   deviation. 
c   Second, if the observed duration of the X-ray eclipse is negative, then
c   the number is meant as an upper limit.  That is, if obv = -10.0,
c   then the eclipse duration is less than 10 degrees.
c
            if(index.eq.9)then
              if(obv(i).lt.0.0d0)then
c
c   Here is the case where the upper limit is less than the computed
c   eclipse duration.
c
                if(dabs(obv(i)).gt.obsparm(index))then
                  chisq=0.0d0
                  iline=iline+1
                  go to 1234
                else
                  if(rmed.ge.1.0d0)THEN
                    chisq=dabs((dabs(obv(i))-obsparm(index))/(eobv(i)))
                  else
                    chisq=(dabs(obv(i))-obsparm(index))*
     @                (dabs(obv(i))-obsparm(index))/(eobv(i)*eobv(i))
                  endif
                  iline=iline+1
                  go to 1234
                endif  
              else
                if(rmed.ge.1.0d0)THEN
                  chisq=dabs((dabs(obv(i))-obsparm(index))/(eobv(i)))
                else
                  chisq=(dabs(obv(i))-obsparm(index))*
     @                (dabs(obv(i))-obsparm(index))/(eobv(i)*eobv(i))
                endif
                iline=iline+1
                go to 1234
              endif
            endif
c
            if(rmed.ge.1.0d0)then
              chisq=dabs((obv(i)-obsparm(index))/(eobv(i)))
            else
              chisq=(obv(i)-obsparm(index))*(obv(i)-obsparm(index))
     $          /(eobv(i)*eobv(i))
            endif
            iline=iline+1
c
c   UPDATE October 28, 2002
c
c   Add 0.0d0 to the end of the 9999999 numbers
c
 1234       if(chisq.gt.9999999.0d0)chisq=9999999.0d0
            if(index.eq.1)write(section(iline),100)chisq
            if(index.eq.2)write(section(iline),100)chisq
            if(index.eq.3)write(section(iline),100)chisq
            if(index.eq.4)write(section(iline),100)chisq
            if(index.eq.5)write(section(iline),100)chisq
            if(index.eq.6)write(section(iline),100)chisq
            if(index.eq.7)write(section(iline),100)chisq
            if(index.eq.8)write(section(iline),100)chisq
c
            if(index.eq.10)write(section(iline),100)chisq
            if(index.eq.11)write(section(iline),100)chisq
            if(index.eq.12)write(section(iline),100)chisq
            if(index.eq.13)write(section(iline),100)chisq
            if(index.eq.14)write(section(iline),100)chisq
            if(index.eq.15)write(section(iline),100)chisq
            if(index.eq.16)write(section(iline),100)chisq
            if(index.eq.17)write(section(iline),100)chisq
            if(index.eq.18)write(section(iline),100)chisq
            if(index.eq.19)write(section(iline),100)chisq
c
c   UPDATE September 11, 2001
c
c   Add this if block
c
            if(index.eq.9)write(section(iline),100)chisq

            K=lnblnk(line)
            line=line(1:K)//section(iline)
 1001     continue
c 
          write(55,300)iter,(section(j),j=1,iline)
c
c          if(iline.gt.0.and.iwrite.eq.0)write(*,200)line
c
 100      format(' ',f12.4)
 600      format(' ',f14.4)
 101      format(' ',f12.4)
 102      format(' ',f12.4)
 103      format(' ',f12.4)
 104      format(' ',f12.4)
 105      format(' ',f12.4)
 106      format(' ',f12.4)
 107      format(' ',f12.4)
 108      format(' ',f12.4)
 109      format(' ',f12.4)
 300      format(i6,1x,20(a13))
c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&
c
c
c  UPDATE FEBRUARY 4, 2005
c
c  Here is a new subroutine that will compute the chi^2 of the luminosity
c  fraction.
c
         subroutine lrchi(Nphase,Nmaxphase,ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochilr,ochi)
c
c
          implicit double precision (a-h,o-z)
c
c   UPDATE September 21, 2008
c
c   Change the dimensions of obv, eobv, and sobv to 11
c
          dimension ymods1(Nmaxphase),ymods2(Nmaxphase),ymodd(Nmaxphase),
     &      ymods3(Nmaxphase)
          dimension obv(19),eobv(19)
          dimension xr(500000),rat(500000)   !should match Nmaxphase
          character*40 sobv(19)
c
         common /medblock/ rmed 
c
          ochilr=0.0d0
          icount=0
          do 10 i=1,Nphase           
            if(ymods1(i).gt.0.0d0)then
              icount=icount+1
              rat(icount)=ymods2(i)/ymods1(i)
              xr(icount)=ymods3(i)
            endif
 10       continue
c
          if(icount.le.2)then
            ochilr=9999999999.0d0
            ochi=ochi+ochilr
            return
          endif

c          write(*,*)'icount,rat,xr   ',icount
          call sort2(icount,rat,xr)
c
          q1=dble(icount/2)
          q2=dble(icount)/2.0
          if(q1.eq.q2)then
            rmedian=(rat(icount/2)+rat(icount/2+1))/2.0
          else
            rmedian=rat(icount/2+1)
          endif
          frac=rmedian
c  
          do 1001 i=1,Nobv
c
c  UPDATE NOVEMBER 23, 2001
c
c  make the argument sobv(i)(1:2)
c
            icn=icnvrt(sobv(i)(1:2))
c
c  UPDATE May 27, 2002
c
c  If the variable rmed > 1, then define ochidisk to be the absolute
c  deviation, rather than the normal chi^2.
c
            if(icn.eq.369)then        ! string lr = 369
              if(rmed.ge.1.0d0)then
                ochilr=dabs((obv(i)-frac)/(eobv(i)))
              else
                ochilr=(obv(i)-frac)*(obv(i)-frac)
     $            /(eobv(i)*eobv(i))
              endif
            endif
 1001     continue
c
          ochi=ochi+ochilr

          return
          end
c
c  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
c  UPDATE FEBRUARY 4, 2005
c
c  Here is a new subroutine that will compute the chi^2 of the luminosity
c  fraction.
c
          subroutine wlrobschi(ochilr)
c
c  July 24, 2000
c
c  This routine will take the observed variables specified in the
c  gridloop.opt file and compute the additional chi^2
c
          implicit double precision (a-h,o-z)
c
c  UPDATE September 21, 2008
c
c  change the dimension of obv, eobv, and sobv to 11
c
c
          character*30 section,outstring    !was 26 QQQQQQQ
          character*80 line
c
         common /medblock/ rmed 
c
          ilength=0
c
          section='1234567890123456'
          section='                '
          line='  -->'
c
c   If rmed > 1, then we are doing median fitting.  In that case, write
c   to a different format statement.
c
           if(rmed.ge.1.0d0)then
             call chistring(' med(lumrat)',ochilr,outstring,lll)
             ilength=ilength+lll
             section=outstring
           else
             call chistring(' chi(lumrat)',ochilr,outstring,lll)
             ilength=ilength+lll
             section=outstring
           endif
           K=lnblnk(line)
           line=line(1:K)//section
           write(*,200)line
c
 200      format(a80)

          return
          end
c
c &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine initNdata(NdataU,NdataB,NdataV,NdataR,
     &        NdataI,NdataJ,NdataH,NdataK,
     %        NRV1,NRV2,Nobv,NRV3)
c
          NdataU=0
          NdataB=0
          NdataV=0
          NdataR=0
          NdataI=0
          NdataJ=0
          NdataH=0
          NdataK=0
          NRV1=0
          NRV2=0
          NRV3=0
          Nobv=0
c
          return
          end
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printchi(chi1)
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

100       format('chi1 = 999999999999999999999.999999')
101       format('chi1 = ',f28.6)
102       format('chi1 = ',f27.6)
103       format('chi1 = ',f26.6)
104       format('chi1 = ',f25.6)
105       format('chi1 = ',f24.6)
106       format('chi1 = ',f23.6)
107       format('chi1 = ',f22.6)
108       format('chi1 = ',f21.6)
109       format('chi1 = ',f20.6)
110       format('chi1 = ',f19.6)
111       format('chi1 = ',f18.6)
112       format('chi1 = ',f17.6)
113       format('chi1 = ',f16.6)
114       format('chi1 = ',f15.6)
115       format('chi1 = ',f14.6)
116       format('chi1 = ',f13.6)
117       format('chi1 = ',f12.6)
118       format('chi1 = ',f11.6)
119       format('chi1 = ',f10.6)
120       format('chi1 = ',f9.6)
121       format('chi1 = ',f8.6)

c
          return
          end
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printmed(chi1)
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

100       format('med1 = 999999999999999999999.999999')
101       format('med1 = ',f28.6)
102       format('med1 = ',f27.6)
103       format('med1 = ',f26.6)
104       format('med1 = ',f25.6)
105       format('med1 = ',f24.6)
106       format('med1 = ',f23.6)
107       format('med1 = ',f22.6)
108       format('med1 = ',f21.6)
109       format('med1 = ',f20.6)
110       format('med1 = ',f19.6)
111       format('med1 = ',f18.6)
112       format('med1 = ',f17.6)
113       format('med1 = ',f16.6)
114       format('med1 = ',f15.6)
115       format('med1 = ',f14.6)
116       format('med1 = ',f13.6)
117       format('med1 = ',f12.6)
118       format('med1 = ',f11.6)
119       format('med1 = ',f10.6)
120       format('med1 = ',f9.6)
121       format('med1 = ',f8.6)

c
          return
          end
c
c  ##########
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printsmall(chi1)
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

100       format('chi_small = 99999999999999999999.999999')
101       format('chi_small = ',f28.6)
102       format('chi_small = ',f27.6)
103       format('chi_small = ',f26.6)
104       format('chi_small = ',f25.6)
105       format('chi_small = ',f24.6)
106       format('chi_small = ',f23.6)
107       format('chi_small = ',f22.6)
108       format('chi_small = ',f21.6)
109       format('chi_small = ',f20.6)
110       format('chi_small = ',f19.6)
111       format('chi_small = ',f18.6)
112       format('chi_small = ',f17.6)
113       format('chi_small = ',f16.6)
114       format('chi_small = ',f15.6)
115       format('chi_small = ',f14.6)
116       format('chi_small = ',f13.6)
117       format('chi_small = ',f12.6)
118       format('chi_small = ',f11.6)
119       format('chi_small = ',f10.6)
120       format('chi_small = ',f9.6)
121       format('chi_small = ',f8.6)

c
          return
          end
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printsmallmed(chi1)
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

100       format('med_small = 99999999999999999999.999999')
101       format('med_small = ',f28.6)
102       format('med_small = ',f27.6)
103       format('med_small = ',f26.6)
104       format('med_small = ',f25.6)
105       format('med_small = ',f24.6)
106       format('med_small = ',f23.6)
107       format('med_small = ',f22.6)
108       format('med_small = ',f21.6)
109       format('med_small = ',f20.6)
110       format('med_small = ',f19.6)
111       format('med_small = ',f18.6)
112       format('med_small = ',f17.6)
113       format('med_small = ',f16.6)
114       format('med_small = ',f15.6)
115       format('med_small = ',f14.6)
116       format('med_small = ',f13.6)
117       format('med_small = ',f12.6)
118       format('med_small = ',f11.6)
119       format('med_small = ',f10.6)
120       format('med_small = ',f9.6)
121       format('med_small = ',f8.6)

c
          return
          end
c
c  &&&&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine writebody3grid(Nalph3,Nbet3,tertperiod,
     @       tertt0,tertecos,tertesin,tertincl,tertOmega,tertQ,dwavex,
     @       dwavey,itconj,it1,it2,it3,it4,tertconj,tertratrad,hh,sw72,
     @       sw73,
     @         P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,P2Q,
     @         P2ratrad,
     @         P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,P3Omega,P3Q,
     @         P3ratrad,
     @         P4tconj,P4period,P4T0,P4ecos,P4esin,P4incl,P4Omega,P4Q,
     @         P4ratrad,
     @         P5tconj,P5period,P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,
     @         P5ratrad,   
     @         P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @         P6ratrad,
     @         P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
     @         P7ratrad,
     @         P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @         P8ratrad)
c
c    will write the correctly formatted file ELCbody3.inp and return
c    default parameters
c
          implicit double precision(a-h,o-z)
c
          dimension dwavex(8,3),dwavey(8,3)
c
          do 5 i=1,8
            dwavex(i,3)=0.635d0
            dwavey(i,3)=0.130d0
 5        continue
c
          open(unit=1,file='gridELCbody3.inp',status='unknown')
c
          write(1,1000)Nalph3
          write(1,1001)Nbet3
          write(1,3001)itconj
          write(1,3002)it1
          write(1,3003)it2
          write(1,3004)it3
          write(1,3005)it4
          write(1,3006)tertconj
          write(1,1002)tertperiod
          write(1,1003)tertt0
          write(1,1004)tertecos
          write(1,1005)tertesin
          write(1,1006)tertincl
          write(1,1007)tertOmega
          write(1,1008)tertQ
          write(1,2001)dwavex(1,3),dwavey(1,3)
          write(1,2002)dwavex(2,3),dwavey(2,3)
          write(1,2003)dwavex(3,3),dwavey(3,3)
          write(1,2004)dwavex(4,3),dwavey(4,3)
          write(1,2005)dwavex(5,3),dwavey(5,3)
          write(1,2006)dwavex(6,3),dwavey(6,3)
          write(1,2007)dwavex(7,3),dwavey(7,3)
          write(1,2008)dwavex(8,3),dwavey(8,3)
          write(1,3007)tertratrad
          write(1,3008)hh
          write(1,3009)sw72
          write(1,3010)sw73
c
          write(1,2020)P2tconj
          write(1,2021)P2period
          write(1,2022)P2T0
          write(1,2023)P2ecos
          write(1,2024)P2esin
          write(1,2025)P2incl
          write(1,2026)P2Omega
          write(1,2027)P2Q
          write(1,2028)P2ratrad

          write(1,2030)P3tconj
          write(1,2031)P3period
          write(1,2032)P3T0
          write(1,2033)P3ecos
          write(1,2034)P3esin
          write(1,2035)P3incl
          write(1,2036)P3Omega
          write(1,2037)P3Q
          write(1,2038)P3ratrad
c
          write(1,2040)P4tconj
          write(1,2041)P4period
          write(1,2042)P4T0
          write(1,2043)P4ecos
          write(1,2044)P4esin
          write(1,2045)P4incl
          write(1,2046)P4Omega
          write(1,2047)P4Q
          write(1,2048)P4ratrad
c
          write(1,2050)P5tconj
          write(1,2051)P5period
          write(1,2052)P5T0
          write(1,2053)P5ecos
          write(1,2054)P5esin
          write(1,2055)P5incl
          write(1,2056)P5Omega
          write(1,2057)P5Q
          write(1,2058)P5ratrad
c
          write(1,2060)P6tconj
          write(1,2061)P6period
          write(1,2062)P6T0
          write(1,2063)P6ecos
          write(1,2064)P6esin
          write(1,2065)P6incl
          write(1,2066)P6Omega
          write(1,2067)P6Q
          write(1,2068)P6ratrad
c
          write(1,2070)P7tconj
          write(1,2071)P7period
          write(1,2072)P7T0
          write(1,2073)P7ecos
          write(1,2074)P7esin
          write(1,2075)P7incl
          write(1,2076)P7Omega
          write(1,2077)P7Q
          write(1,2078)P7ratrad
c
          write(1,2080)P8tconj
          write(1,2081)P8period
          write(1,2082)P8T0
          write(1,2083)P8ecos
          write(1,2084)P8esin
          write(1,2085)P8incl
          write(1,2086)P8Omega
          write(1,2087)P8Q
          write(1,2088)P8ratrad

          close(1)
c
c 100      format(a1,'I can''t find the file ''ELCbody3.inp''!',
c     @     '   I''m making one up and setting default values')
c
 3001     format(i1,24x,'itconj (0=T_peri, 1=T_tran, 2=T_occul)')
 3002     format(i1,24x,'set to 1 for logarithmic mass ratios')
 3003     format(i1,24x,'set to 1 for informational output')
 3004     format(i1,24x,'it3 (currently inactive)')
 3005     format(i1,24x,'it4 (currently inactive)')
 3006     format(f21.13,4x,'tertconj             tag tj')
 3007     format(f15.13,10x,'tertratrad (P1 radius to star 1 radius,' 
     @            'tag tb)')
 3008     format(f16.8,9x,'hh (step size for dynamical integration)')
 3009     format(f10.8,15x,'rk1  (apsidal constant, star 1   tag a1)')
 3010     format(f10.8,15x,'rk2  (apsidal constant, star 2   tag a2)')

 1000     format(i4,21x,'Nalph3')
 1001     format(i4,21x,'Nbet3')
 1002     format(f22.14,3x,  'tertperiod (days)    tag tt')
 1003     format(f21.13,4x,  'tertT0               tag tu')
 1004     format(f20.17,5x, 'terte*cos(omega)     tag tv')
 1005     format(f20.17,5x, 'terte*sin(omega)     tag tw')
 1006     format(f20.15,5x,  'tertincl (degrees)   tag tx')
 1007     format(f20.15,5x,  'tertOmega (degrees)  tag ty')
 1008     format(f23.11,2x,  'tertQ (EB/body3)     tag tz')
 2001     format(2(f11.7,1x),1x,'limb darkening coefficients filter 1')
 2002     format(2(f11.7,1x),1x,'limb darkening coefficients filter 2')
 2003     format(2(f11.7,1x),1x,'limb darkening coefficients filter 3')
 2004     format(2(f11.7,1x),1x,'limb darkening coefficients filter 4')
 2005     format(2(f11.7,1x),1x,'limb darkening coefficients filter 5')
 2006     format(2(f11.7,1x),1x,'limb darkening coefficients filter 6')
 2007     format(2(f11.7,1x),1x,'limb darkening coefficients filter 7')
 2008     format(2(f11.7,1x),1x,'limb darkening coefficients filter 8')
c
 2020     format(f18.11,7x,'P2tconj              tag uj')
 2021     format(f22.14,3x,'P2period (days)      tag ut')
 2022     format(f21.13,4x,'P2T0                 tag uu')
 2023     format(f20.17,5x,'P2ecos               tag uv')
 2024     format(f20.17,5x,'P2esin               tag uw')
 2025     format(f20.15,5x,'P2incl (degrees)     tag ux')
 2026     format(f20.15,5x,'P2Omega (degrees)    tag uy')
 2027     format(f23.11,2x,'P2Q (EB/body3)       tag uz')
 2028     format(f18.11,7x,'P2ratrad             tag ub')
c
 2030     format(f18.11,7x,'P3tconj              tag vj')
 2031     format(f22.14,3x,'P3period (days)      tag vt')
 2032     format(f21.13,4x,'P3T0                 tag vu')
 2033     format(f20.17,5x,'P3ecos               tag vv')
 2034     format(f20.17,5x,'P3esin               tag vw')
 2035     format(f20.15,5x,'P3incl (degrees)     tag vx')
 2036     format(f20.15,5x,'P3Omega (degrees)    tag vy')
 2037     format(f23.11,2x,'P3Q (EB/body3)       tag vz')
 2038     format(f18.11,7x,'P3ratrad             tag vb')
c
2040      format(f18.11,7x,'P4tconj              tag wj')
2041      format(f22.14,3x,'P4period (days)      tag wt')
2042      format(f21.13,4x,'P4T0                 tag wu')
2043      format(f20.17,5x,'P4ecos               tag wv')
2044      format(f20.17,5x,'P4esin               tag ww')
2045      format(f20.15,5x,'P4incl (degrees)     tag wx')
2046      format(f20.15,5x,'P4Omega (degrees)    tag wy')
2047      format(f23.11,2x,'P4Q (EB/body3)       tag wz')
2048      format(f18.11,7x,'P4ratrad             tag wb')
c
2050      format(f18.11,7x,'P5tconj              tag xj')
2051      format(f22.14,3x,'P5period (days)      tag xt')
2052      format(f21.13,4x,'P5T0                 tag xu')
2053      format(f20.17,5x,'P5ecos               tag xv')
2054      format(f20.17,5x,'P5esin               tag xw')
2055      format(f20.15,5x,'P5incl (degrees)     tag xx')
2056      format(f20.15,5x,'P5Omega (degrees)    tag xy')
2057      format(f23.11,2x,'P5Q (EB/body3)       tag xz')
2058      format(f18.11,7x,'P5ratrad             tag xb')
c
2060      format(f18.11,7x,'P6tconj              tag sj')
2061      format(f22.14,3x,'P6period (days)      tag st')
2062      format(f21.13,4x,'P6T0                 tag su')
2063      format(f20.17,5x,'P6ecos               tag sv')
2064      format(f20.17,5x,'P6esin               tag sw')
2065      format(f20.15,5x,'P6incl (degrees)     tag sx')
2066      format(f20.15,5x,'P6Omega (degrees)    tag sy')
2067      format(f23.11,2x,'P6Q (EB/body3)       tag sz')
2068      format(f18.11,7x,'P6ratrad             tag sb')
c
2070      format(f18.11,7x,'P7tconj              tag hj')
2071      format(f22.14,3x,'P7period (days)      tag ht')
2072      format(f21.13,4x,'P7T0                 tag hu')
2073      format(f20.17,5x,'P7ecos               tag hv')
2074      format(f20.17,5x,'P7esin               tag hw')
2075      format(f20.15,5x,'P7incl (degrees)     tag hx')
2076      format(f20.15,5x,'P7Omega (degrees)    tag hy')
2077      format(f23.11,2x,'P7Q (EB/body3)       tag hz')
2078      format(f18.11,7x,'P7ratrad             tag hb')
c
2080      format(f18.11,7x,'P8tconj              tag kj')
2081      format(f22.14,3x,'P8period (days)      tag kt')
2082      format(f21.13,4x,'P8T0                 tag ku')
2083      format(f20.17,5x,'P8ecos               tag kv')
2084      format(f20.17,5x,'P8esin               tag kw')
2085      format(f20.15,5x,'P8incl (degrees)     tag kx')
2086      format(f20.15,5x,'P8Omega (degrees)    tag ky')
2087      format(f23.11,2x,'P8Q (EB/body3)       tag kz')
2088      format(f18.11,7x,'P8ratrad             tag kb')
c


          return
          end
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
 100      format(a,'=',f8.6)
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
c 101      format(a,'=',f8.5)
 101      format(a,'=',f10.7)
c
          if(chisq.lt.1000.0d0)then
            write(outstring,102)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
c 102      format(a,'=',f8.4)
 102      format(a,'=',f11.7)
c
          if(chisq.lt.10000.0d0)then
            write(outstring,103)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
c 103      format(a,'=',f8.3)
 103      format(a,'=',f12.7)
c
          if(chisq.lt.100000.0d0)then
            write(outstring,104)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
c 104      format(a,'=',f9.3)
 104      format(a,'=',f14.8)
c
          if(chisq.lt.1000000.0d0)then
            write(outstring,105)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
c 105      format(a,'=',f10.3)
 105      format(a,'=',f15.8)
c
          if(chisq.lt.10000000.0d0)then
            write(outstring,106)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
c 106      format(a,'=',f10.2)
 106      format(a,'=',f12.4)
c
          if(chisq.lt.100000000.0d0)then
            write(outstring,107)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
c 107      format(a,'=',f11.2)
 107      format(a,'=',f13.4)
c
          if(chisq.lt.1000000000.0d0)then
            write(outstring,108)instring,chisq
            ilength=lnblnk(outstring)
            return
          endif
c 108      format(a,'=',f12.2)
 108      format(a,'=',f14.4)
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
          subroutine printiter(kmodel,instring1,iter,instring2)
c
          character*(*)instring1,instring2
          character*(80) outstring1,outstring2
c
          if(kmodel.lt.0)return
          if(kmodel.ge.1000000)return
c
          l1=lnblnk(instring1)
          l1=l1+1
          if(kmodel.lt.10)write(outstring1,100)instring1(1:l1),kmodel
          if((kmodel.ge.10).and.(kmodel.lt.100))then
            write(outstring1,101)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.100).and.(kmodel.lt.1000))then
            write(outstring1,102)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.1000).and.(kmodel.lt.10000))then
            write(outstring1,103)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.10000).and.(kmodel.lt.100000))then
            write(outstring1,104)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.100000).and.(kmodel.lt.1000000))then
            write(outstring1,105)instring1(1:l1),kmodel
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
          endif !nn gt.10 and .lt 99
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
          endif !nn between 100 and 999
c
          if((nn.ge.1000).and.(nn.lt.10000))then
            if(chi1.gt.9.999999999999999D20)then
              write(*,400)nn
              return
            endif
            if(chi1.gt.9.999999999999999D19)then
              write(*,401)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D18)then
              write(*,402)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D17)then
              write(*,403)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D16)then
              write(*,404)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D15)then
              write(*,405)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D14)then
              write(*,406)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D13)then
              write(*,407)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D12)then
              write(*,408)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D11)then
              write(*,409)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D10)then
              write(*,410)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D9)then
              write(*,411)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D8)then
              write(*,412)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D7)then
              write(*,413)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D6)then
              write(*,414)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D5)then
              write(*,415)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D4)then
              write(*,416)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D3)then
              write(*,417)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D2)then
              write(*,418)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D1)then
              write(*,419)nn,chi1
              return
            endif
            if(chi1.gt.9.999999999999999D0)then
              write(*,420)nn,chi1
              return
            endif
            write(*,421)nn,chi1
          endif !nn between 1000 and 9999

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
300       format('chi',i3,' = 999999999999999999999.999999')
301       format('chi',i3,' = ',f28.6)
302       format('chi',i3,' = ',f27.6)
303       format('chi',i3,' = ',f26.6)
304       format('chi',i3,' = ',f25.6)
305       format('chi',i3,' = ',f24.6)
306       format('chi',i3,' = ',f23.6)
307       format('chi',i3,' = ',f22.6)
308       format('chi',i3,' = ',f21.6)
309       format('chi',i3,' = ',f20.6)
310       format('chi',i3,' = ',f19.6)
311       format('chi',i3,' = ',f18.6)
312       format('chi',i3,' = ',f17.6)
313       format('chi',i3,' = ',f16.6)
314       format('chi',i3,' = ',f15.6)
315       format('chi',i3,' = ',f14.6)
316       format('chi',i3,' = ',f13.6)
317       format('chi',i3,' = ',f12.6)
318       format('chi',i3,' = ',f11.6)
319       format('chi',i3,' = ',f10.6)
320       format('chi',i3,' = ',f9.6)
321       format('chi',i3,' = ',f8.6)
c
400       format('chi',i4,' = 999999999999999999999.999999')
401       format('chi',i4,' = ',f28.6)
402       format('chi',i4,' = ',f27.6)
403       format('chi',i4,' = ',f26.6)
404       format('chi',i4,' = ',f25.6)
405       format('chi',i4,' = ',f24.6)
406       format('chi',i4,' = ',f23.6)
407       format('chi',i4,' = ',f22.6)
408       format('chi',i4,' = ',f21.6)
409       format('chi',i4,' = ',f20.6)
410       format('chi',i4,' = ',f19.6)
411       format('chi',i4,' = ',f18.6)
412       format('chi',i4,' = ',f17.6)
413       format('chi',i4,' = ',f16.6)
414       format('chi',i4,' = ',f15.6)
415       format('chi',i4,' = ',f14.6)
416       format('chi',i4,' = ',f13.6)
417       format('chi',i4,' = ',f12.6)
418       format('chi',i4,' = ',f11.6)
419       format('chi',i4,' = ',f10.6)
420       format('chi',i4,' = ',f9.6)
421       format('chi',i4,' = ',f8.6)

c
          return
          end
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printiter1(kmodel,instring1)
c
          character*(*)instring1
          character*(80) outstring1
c
          if(kmodel.lt.0)return
          if(kmodel.ge.1000000)return
c
          l1=lnblnk(instring1)
          l1=l1+1
          if(kmodel.lt.10)write(outstring1,100)instring1(1:l1),kmodel
          if((kmodel.ge.10).and.(kmodel.lt.100))then
            write(outstring1,101)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.100).and.(kmodel.lt.1000))then
            write(outstring1,102)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.1000).and.(kmodel.lt.10000))then
            write(outstring1,103)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.10000).and.(kmodel.lt.100000))then
            write(outstring1,104)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.100000).and.(kmodel.lt.1000000))then
            write(outstring1,105)instring1(1:l1),kmodel
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
          subroutine getchilimb(Nvmax,Nterms,svar,dwavex,dwavey,chilimb)
c
          implicit double precision (a-h,o-z)
c
          dimension dwavex(8,3),dwavey(8,3)
c
          character*40 svar(Nvmax)
c
          chilimb=0.0d0
          do j=1,Nterms
            kk=icnvrt(svar(j)(1:2))
            if(kk.ge.1752.and.kk.le.1759)then
              do jjj=1,8
c                if(j.eq.jjj)then                                  
                  rlimbsum=dwavex(jjj,1)+dwavey(jjj,1)
                  if(rlimbsum.gt.1.0d0)chilimb=9.99d18
                  rlimbsum=dwavex(jjj,2)+dwavey(jjj,2)
                  if(rlimbsum.gt.1.0d0)chilimb=9.99d18
c                endif                                        
              enddo
            endif
            if(kk.ge.1784.and.kk.le.1791)then
              do jjj=1,8
c                if(j.eq.jjj)then                                    
                  rlimbsum=dwavex(jjj,1)+dwavey(jjj,1)
                  if(rlimbsum.gt.1.0d0)chilimb=9.99d18
                  rlimbsum=dwavex(jjj,2)+dwavey(jjj,2)
                  if(rlimbsum.gt.1.0d0)chilimb=9.99d18
c                endif                                               
              enddo
            endif
            if(kk.ge.1816.and.kk.le.1823)then
              do jjj=1,8
c                if(j.eq.jjj)then                                
                  rlimbsum=dwavex(jjj,1)+dwavey(jjj,1)
                  if(rlimbsum.gt.1.0d0)chilimb=9.99d18
                  rlimbsum=dwavex(jjj,2)+dwavey(jjj,2)
                  if(rlimbsum.gt.1.0d0)chilimb=9.99d18
c                endif                                               
              enddo
            endif
            if(kk.ge.1720.and.kk.le.1727)then
              do jjj=1,8
c                if(j.eq.jjj)then  
                  rlimbsum=dwavex(jjj,1)+dwavey(jjj,1)
                  if(rlimbsum.gt.1.0d0)chilimb=9.99d18
                  rlimbsum=dwavex(jjj,2)+dwavey(jjj,2)
                  if(rlimbsum.gt.1.0d0)chilimb=9.99d18
c                endif 
              enddo
            endif
          enddo
c
          return
          end
c
c &^&%$*$*$*$^#@^#^$^&&%$%$#@$&^&%$*$*$*$^#@^#^$^&
c
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     $     NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=3.e-16,RNMX=1.-EPS)
      INTEGER j,k,iivv(NTAB),iiyy
      SAVE iivv,iiyy
      DATA iivv /NTAB*0/, iiyy /0/
      if (idum.le.0.or.iiyy.eq.0) then
         idum=max(-idum,1)
         do j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iivv(j)=idum
         enddo
         iiyy=iivv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iiyy/NDIV
      iiyy=iivv(j)
      iivv(j)=idum
      ran1=min(AM*iiyy,RNMX)
      return
      END
c
c #!!!%(@^%#%*(&(#@#*$^@^$##*#@&(%)
c
          function gasdev(idum)
c
          implicit double precision(a-h,o-z)

          save iset,gset
          data iset/0/
c
          if(iset.eq.0)then
 1          v1=2.0d0*ran1(idum)-1.0d0
            v2=2.0d0*ran1(idum)-1.0d0
            rsq=v1**2+v2**2
            if(rsq.ge.1.0d0.or.rsq.eq.0.0d0)go to 1
            fac=dsqrt(-2.0d0*dlog(rsq)/rsq)
            gset=v1*fac
            gasdev=v2*fac
            iset=1
          else
            gasdev=gset
            iset=0
          endif
          return
          end
c
c  #!!!%(@^%#%*(&(#@#*$^@^$##*#@&(%)@**$!^(#((#&%!)
c
          subroutine monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     &      ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,
     @      RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,fracs1,
     @      fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
     @      chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,chisqRV1,
     &      chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,zeroU,resU,
     @      NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,
     @      errV,zeroV,resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
     @      xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,errJ,
     %      zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,NdataK,
     $      xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,errRV1,NRV2,
     @      xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,sobv,obv,eobv,ochi,
     &      ochidisk,ochilr,Nvmax,svar,var,saveasini,
     @      savxdataU,savydataU,saverrU,savxdataB,savydataB,
     @      saverrB,savxdataV,savydataV,saverrV,
     $      savxdataR,savydataR,saverrR,savxdataI,savydataI,
     &      saverrI,savxdataJ,savydataJ,saverrJ,
     &      savxdataH,savydataH,saverrH,savxdataK,savydataK,
     @      saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,savyRV2,saverrRV2,
     @      ifrac,ilum,i16,isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,
     @      isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,ibest,
     @      ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,small,
     @      Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,icnarray,
     @      RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
     @   parmstring,planetparm,dynparm,line,fill1,fill2,omega1,
     @    omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     @    rinner,router,tdisk,xi,alb1,alb2,rLx,Period,fm,separ,
     @    gamma,wave,dbolx,dboly,dwavex,dwavey,t3,g3,SA3,density,
     @    sw1,sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     @    primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @    bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,
     @    sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,omegadot,
     @    contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49,gaplow,
     @    gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @    P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @    P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @    P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @    P5esin,P5incl,P5Omega,P5Q,P5ratrad,P6tconj,P6period,P6T0,
     @    P6ecos,P6esin,P6incl,P6Omega,P6Q,P6ratrad,P7tconj,P7period,
     @    P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,
     @    P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,
     @    xSC,ySC,spot1parm,spot2parm,spotdparm,Nalph1,Nbet1,Nalph2,
     @    Nbet2,Ntheta,Nradius,Nref,idraw,iecheck,iidint,iatm,ism1,
     @    ilaw,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
     @    iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,isw5,isw6,isw7,
     @    isw8,isw9,idark1,idark2,isw12,isw13,isw21,isw22,isw23,isw24,
     @    isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,
     @    NSC,compfracs,tertperiod,tertt0,tertecos,tertesin,
     @    tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72,sw73,
     @    Nmaxeclipse,Tdur1,Tdur2)
c
c   January 1, 2013
c
c   This subroutine will collect all of the steps necessary to compute
c   models, fold the data, compute chi^2, and print the variables.
c
c          implicit double precision (a-h,o-z)
c
          implicit none

          integer Nphase,Nmaxphase,maxlines,maxmu,nmu,
     @      ifastflag,NdataU,NdataB,NdataR,NdataI,
     @      NdataJ,NdataH,NdataK,Nobv,Nvmax,NRV1,NRV2,NRV3,
     @      ifrac,ilum,i16,isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,
     @      isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,ibest,
     @      ifixgamma,ichilabel,Ncycle,Nobscycle,icnarray,icnRV3
          integer nalph1,nbet1,nalph2,nbet2,ntheta,nradius,Nref,idraw
          integer iecheck,iidint,iatm,ism1,ilaw,icnU,icnB,icnV,icnR
          integer icnI,icnJ,icnH,icnK,icnRV1,icnRV2,iRVfilt,isw1,isw2
          integer isw3,isw4,ikeep,isynch,isw5,isw6,isw7,isw8,isw9
          integer idark1,idark2,isw12,isw13,isw21,isw22,isw23,isw24
          integer isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32
          integer isw33,isw34,nSC,idum,Nlines,Nalph3,Nbet3,itconj,it2
          integer it3,it4,it1,NRVphase,NdataV,nbin,isvel3,lll,lnblnk
          integer ljk,lkk,ljj,ijk,Nmaxeclipse

          real*8 gamma3,chitimetot,cccc,Tdur1,Tdur2
          real*8 savP5Q,savP6Q,savP7Q,savP8Q,chisqRV3,gamma1,gamma2
          real*8 remf1,remf2,gasdev,savtertQ,savP2Q,savP3Q,savP4Q
          real*8 tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72
          real*8 rmed,tertperiod,tertT0,tertecos,tertesin,tmin,tmax
          real*8 xmod,ymodU,ymodB,ymodV,wave,dbolx,dboly,dwavex,gmin,
     &      ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,
     @      RV1,RV2,drv1,drv2,obsparm,xRVmod,fracs1,dwavey,gmax,sw73,
     @      fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
     @      chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,chisqRV1,
     &      chisqRV2,chilimb,chi1,xdataU,ydataU,errU,zeroU,resU,
     @      xdataB,ydataB,errB,zeroB,resB,xdataV,ydataV,spotdparm,
     @      errV,zeroV,resV,xdataR,ydataR,errR,zeroR,resR,spot1parm,
     @      ydataI,errI,zeroI,resI,xdataJ,ydataJ,errJ,spot2parm,
     %      zeroJ,resJ,xdataH,ydataH,errH,zeroH,resH,xsob,powercoeff,
     $      xdataK,ydataK,errK,zeroK,resK,xRV1,yRV1,errRV1,atmint1,
     @      xRV2,yRV2,errRV2,ggamma1,ggamma2,obv,eobv,ochi,atmint2,
     &      ochidisk,ochilr,var,saveasini,atmT,atmg,atmmu,atmint3,
     @      savxdataU,savydataU,saverrU,savxdataB,savydataB,atmint4,
     @      saverrB,savxdataV,savydataV,saverrV,atmint5,atmint6,atmint7,
     $      savxdataR,savydataR,saverrR,savxdataI,savydataI,atmint8,
     &      saverrI,savxdataJ,savydataJ,saverrJ,compfracs,xdataI,
     &      savxdataH,savydataH,saverrH,savxdataK,savydataK,
     @      saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,savyRV2,saverrRV2,
     @      savesep,resRV1,resRV2,thresh,small,chitimearray,gaplow,
     @      Ttimes,Tseps,obsTtimes,obsTerr,gaphigh,xSC,ySC,fill1,fill2,
     @      RV3,xRV3,yRV3,errRV3,resRV3,ggamma3,omega1,omega2,dphase

          real*8 p2ratrad,p3ratrad,p4ratrad,p5ratrad,p6ratrad,p7ratrad,p8ratrad
          real*8 Q,finc,teff1,teff2,tgrav1,tgrav2,betarim,rinner,router
          real*8 tdisk,xi,alb1,alb2,rlx,period,fm,separ,gamma,t3,g3,SA3
          real*8 density,sw1,sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8
          real*8 sw9,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,sw30
          real*8 temprat,bigI,bigbeta,sw23,sw24,sw25,sw26,sw27,sw28,sw29
          real*8 contam,tconj,beam1,beam2,ocose,osine,omegadot,contamS0
          real*8 contamS1,contamS2,contamS3,sw47,sw48,sw49,pie
          real*8 p2tconj,p2period,p2t0,p2ecos,p2esin,p2incl,p2Omega,p2Q
          real*8 p3tconj,p3period,p3t0,p3ecos,p3esin,p3incl,p3Omega,p3Q
          real*8 p4tconj,p4period,p4t0,p4ecos,p4esin,p4incl,p4Omega,p4Q
          real*8 p5tconj,p5period,p5t0,p5ecos,p5esin,p5incl,p5Omega,p5Q
          real*8 p6tconj,p6period,p6t0,p6ecos,p6esin,p6incl,p6Omega,p6Q
          real*8 p7tconj,p7period,p7t0,p7ecos,p7esin,p7incl,p7Omega,p7Q
          real*8 p8tconj,p8period,p8t0,p8ecos,p8esin,p8incl,p8Omega,p8Q

          dimension obsparm(19),obv(19),eobv(19)
          dimension xmod(Nmaxphase),ymodU(Nmaxphase),ymodB(Nmaxphase),
     $      ymodV(Nmaxphase),ymodR(Nmaxphase),ymodI(Nmaxphase),
     $      ymodJ(Nmaxphase),ymodH(Nmaxphase),ymodK(Nmaxphase),
     &      ymods1(Nmaxphase),ymods2(Nmaxphase),ymodd(Nmaxphase),
     &      RV1(Nmaxphase),RV2(Nmaxphase),ymods3(Nmaxphase),
     @      RV3(Nmaxphase)
          dimension xdataU(Ndatamax),ydataU(Ndatamax),errU(Ndatamax),
     &      ydataB(Ndatamax),errB(Ndatamax),ydataV(Ndatamax),
     @      errV(Ndatamax),ydataR(Ndatamax),errR(Ndatamax),
     @      ydataI(Ndatamax),errI(Ndatamax),ydataJ(Ndatamax),
     @      errJ(Ndatamax),ydataH(Ndatamax),errH(Ndatamax),
     $      ydataK(Ndatamax),errK(Ndatamax),yRV1(Ndatamax),
     @      errRV1(Ndatamax),yRV2(Ndatamax),errRV2(Ndatamax),
     @      xdataB(Ndatamax),xdataV(Ndatamax),xdataR(Ndatamax),
     @      xdataI(Ndatamax),xdataJ(Ndatamax),xdataH(Ndatamax),
     @      xdataK(Ndatamax),xRV1(Ndatamax),xRV2(Ndatamax),
     $      drv1(Nmaxphase),drv2(Nmaxphase),xRV3(Ndatamax),
     @      yRV3(Ndatamax),errRV3(Ndatamax)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,3),
     @       dwavey(8,3)
          dimension var(Nvmax)
          character*40 svar(Nvmax)
          character*40 sobv(Nvmax)
          dimension xRVmod(Nmaxphase)
          dimension resU(Ndatamax),resB(Ndatamax),resV(Ndatamax)
          dimension resR(Ndatamax),resI(Ndatamax),resJ(Ndatamax)
          dimension resH(Ndatamax),resK(Ndatamax),resRV1(Ndatamax)
          dimension resRV2(Ndatamax),resRV3(Ndatamax)
c
          character*1700 line     
c
          dimension spotdparm(2,4),spot1parm(2,4),spot2parm(2,4)
          dimension xsob(2)
          dimension powercoeff(8,9)
          character*1000 parmstring
          character*2000 planetparm
          character*1700 dynparm
c
          parameter (maxlines=1300,maxmu=115)  
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)

c
c   New variables are needed to allow for the fitting of period and T0.
c
          dimension savxdataU(Ndatamax),savydataU(Ndatamax),
     @      saverrU(Ndatamax),savydataB(Ndatamax),saverrB(Ndatamax),
     @      savydataV(Ndatamax),saverrV(Ndatamax),savydataR(Ndatamax),
     @      saverrR(Ndatamax),savydataI(Ndatamax),saverrI(Ndatamax),
     &      savydataJ(Ndatamax),saverrJ(Ndatamax),savydataH(Ndatamax),
     @      saverrH(Ndatamax),savydataK(Ndatamax),saverrK(Ndatamax),
     @      savyRV1(Ndatamax),saverrRV1(Ndatamax),savyRV2(Ndatamax),
     @      saverrRV2(Ndatamax),savxdataB(Ndatamax),savxdataV(Ndatamax),
     &      savxdataR(Ndatamax),savxdataI(Ndatamax),savxdataJ(Ndatamax),
     &      savxdataH(Ndatamax),savxdataK(Ndatamax),savxRV1(Ndatamax),
     @      savxRV2(Ndatamax)
c
          dimension compfracs(8,3),ochidisk(8)
c
          dimension fracs1(Nmaxphase,4),fracs2(Nmaxphase,4)
          dimension fracs3(Nmaxphase,4),fracs4(Nmaxphase,4)
          dimension fracs5(Nmaxphase,4),fracs6(Nmaxphase,4)
          dimension fracs7(Nmaxphase,4),fracs8(Nmaxphase,4)
          dimension Ncycle(34),Ttimes(34,Nmaxeclipse)
          dimension Tseps(34,Nmaxeclipse),Tdur1(34,Nmaxphase)
          dimension Nobscycle(34),obsTtimes(34,Nmaxeclipse)
          dimension obsTerr(34,Nmaxeclipse),Tdur2(34,Nmaxeclipse)
          dimension icnarray(34),chitimearray(34)
c
          dimension gaplow(9999),gaphigh(9999)
          dimension xSC(9999),ySC(9999)
c
c          common /stringblock/ parmstring,planetparm,dynparm,line
cc
c          common /realblock/ fill1,fill2,omega1,omega2,dphase,Q,finc,
c     @      Teff1,Teff2,Tgrav1,Tgrav2,betarim,rinner,router,tdisk,xi,
c     @      alb1,alb2,rLx,Period,fm,separ,gamma,wave,dbolx,dboly,dwavex,
c     @      dwavey,t3,g3,SA3,density,sw1,sw2,sw3,T0,ecc,argper,pshift,
c     @      sw5,sw6,sw7,sw8,sw9,primmass,primK,primrad,ratrad,frac1,
c     @      frac2,ecosw,temprat,bigI,bigbeta,sw23,sw24,powercoeff,sw25,
c     @      sw26,sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,ocose,
c     @      osine,omegadot,contamS0,contamS1,contamS2,contamS3,sw47,
c     @      sw48,sw49,gaplow,gaphigh,
c     @         P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,P2Q,
c     @         P2ratrad,
c     @         P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,P3Omega,P3Q,
c     @         P3ratrad,
c     @         P4tconj,P4period,P4T0,P4ecos,P4esin,P4incl,P4Omega,P4Q,
c     @         P4ratrad,
c     @         P5tconj,P5period,P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,
c     @         P5ratrad,   
c     @         P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
c     @         P6ratrad,
c     @         P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
c     @         P7ratrad,
c     @         P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
c     @         P8ratrad,xSC,ySC
cc
c          common /spotblock/ spot1parm,spot2parm,spotdparm
cc
c          common /intblock/ Nalph1,Nbet1,Nalph2,Nbet2,Ntheta,Nradius,
c     @      Nref,idraw,iecheck,iidint,iatm,ism1,ilaw,icnU,icnB,icnV,icnR,
c     @      icnI,icnJ,icnH,icnK,icnRV1,icnRV2,iRVfilt,isw1,isw2,isw3,
c     @      isw4,ikeep,isynch,isw5,isw6,isw7,isw8,isw9,idark1,idark2,
c     @      isw12,isw13,isw21,isw22,isw23,isw24,isw25,isw26,isw27,isw28,
c     @      isw29,isw30,isw31,isw32,isw33,isw34,NSC
cc
c          common /fracblock/ compfracs          
cc
c          common /thirdblock/ tertperiod,tertt0,tertecos,tertesin,
c     @        tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72,sw73
cc
          common /realatm/ atmT,atmg,atmmu,atmint1,atmint2,atmint3,
     &      atmint4,atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin

          common /intatm/  Nlines,Nmu,Nalph3,Nbet3,itconj,it1,it2,it3,
     @      it4
c
          common /medblock/ rmed 
c
          common /ranblock/ idum
c
c
          parameter(pie=3.14159265358979323d0)

          if(savesep.lt.0.0d0)separ=savesep 
c
          call assignvar(Nvmax,svar,var,fill1,fill2,omega1,
     @       omega2,Q,finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,
     @       rLx,separ,gamma,t3,g3,sa3,ecc,argper,pshift,spot1parm,
     @       spot2parm,spotdparm,period,T0,alb1,alb2,dwavex,dwavey,
     @       primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @       bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2,contam,
     @       ocose,osine,isw29,tertperiod,tertt0,tertecos,tertesin,
     @       tertincl,tertOmega,tertQ,Tgrav1,Tgrav2,tertconj,omegadot,
     @       contamS0,contamS1,contamS2,contamS3,
     @       P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @       P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @       P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @       P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @       P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @       P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @       P6ratrad,
     @       P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
     @       P7ratrad,
     @       P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @       P8ratrad,sw72,sw73)
c
c   if isw28 is greater than 0, then set the value of T0 based
c   on the time of transit Tconj.
c
          thresh=sw48
          Nbin=0
          if(isw28.gt.0)then
            if(isw29.gt.0)then
              ecc=dsqrt(ocose*ocose+osine*osine)
              argper=datan2(osine,ocose)*180.0d0/pie
              if(argper.le.0.0d0)argper=argper+360.0d0
              if(argper.gt.360.0d0)argper=argper-360.0d0
            endif
            call getT0(finc,period,ecc,argper,T0,Tconj)
          endif
c
c  Save the filling factors in case there is an eccentric orbit
c
          remf1=fill1
          remf2=fill2
c
c  Get the models.
c
c  If sw25>0, use that as an error for asini (sw5)
c
          if(sw25.gt.0.0d0)sw5=saveasini+gasdev(idum)*sw25
c
c  if we are doing the body3, and it1=1, then assume mass
c  ratios are logs
c
          if((isw30.ge.1).and.(it1.eq.1))then
            savtertQ=tertQ
            savP2Q=P2Q
            savP3Q=P3Q
            savP4Q=P4Q
            savP5Q=P5Q
            savP6Q=P6Q
            savP7Q=P7Q
            savP8Q=P8Q
            tertQ=10.0d0**(tertQ)
            P2Q=10.0d0**(P2Q)
            P3Q=10.0d0**(P3Q)
            P4Q=10.0d0**(P4Q)
            P5Q=10.0d0**(P5Q)
            P6Q=10.0d0**(P6Q)
            P7Q=10.0d0**(P7Q)
            P8Q=10.0d0**(P8Q)
          endif
c
          call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     @             ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %             ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @             xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     @             fracs7,fracs8,Ncycle,Ttimes,Tseps,RV3,
     @      parmstring,planetparm,dynparm,line,fill1,fill2,omega1,
     @    omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     @    rinner,router,tdisk,xi,alb1,alb2,rLx,Period,fm,separ,
     @    gamma,wave,dbolx,dboly,dwavex,dwavey,t3,g3,SA3,density,
     @    sw1,sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     @    primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @    bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,
     @    sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,omegadot,
     @    contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49,gaplow,
     @    gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @    P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @    P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @    P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @    P5esin,P5incl,P5Omega,P5Q,P5ratrad,P6tconj,P6period,P6T0,
     @    P6ecos,P6esin,P6incl,P6Omega,P6Q,P6ratrad,P7tconj,P7period,
     @    P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,
     @    P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,
     @    xSC,ySC,spot1parm,spot2parm,spotdparm,Nalph1,Nbet1,Nalph2,
     @    Nbet2,Ntheta,Nradius,Nref,idraw,iecheck,iidint,iatm,ism1,
     @    ilaw,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
     @    iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,isw5,isw6,isw7,
     @    isw8,isw9,idark1,idark2,isw12,isw13,isw21,isw22,isw23,isw24,
     @    isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,
     @    NSC,compfracs,tertperiod,tertt0,tertecos,tertesin,
     @    tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72,sw73,
     @    Nmaxeclipse,Tdur1,Tdur2)
c
          if((isw30.ge.1).and.(it1.eq.1))then
            tertQ=savtertQ
            P2Q=savP2Q
            P3Q=savP3Q
            P4Q=savP4Q
            P5Q=savP5Q
            P6Q=savP6Q
            P7Q=savP7Q
            P8Q=savP8Q
          endif
c
          fill1=remf1
          fill2=remf2
c
c
c   August 2, 2001
c
c   If we are fitting for the period and/or T0 (itime>0), then we
c   have to fold the saved arrays.
c
          if(isw7.gt.0)then
            if(icnU.ne.430)call phasefold(Ndatamax,isavNU,savxdataU,
     @        savydataU,saverrU,Nbin,NdataU,xdataU,ydataU,errU,
     @        period,T0,ikeep,ecc,argper,isw7)
            if(icnB.ne.430)call phasefold(Ndatamax,isavNB,savxdataB,
     @        savydataB,saverrB,Nbin,NdataB,xdataB,ydataB,errB,
     @        period,T0,ikeep,ecc,argper,isw7)
            if(icnV.ne.430)call phasefold(Ndatamax,isavNV,savxdataV,
     @        savydataV,saverrV,Nbin,NdataV,xdataV,ydataV,errV,
     @        period,T0,ikeep,ecc,argper,isw7)
            if(icnR.ne.430)call phasefold(Ndatamax,isavNR,savxdataR,
     @        savydataR,saverrR,Nbin,NdataR,xdataR,ydataR,errR,
     @        period,T0,ikeep,ecc,argper,isw7)
            if(icnI.ne.430)call phasefold(Ndatamax,isavNI,savxdataI,
     @        savydataI,saverrI,Nbin,NdataI,xdataI,ydataI,errI,
     @        period,T0,ikeep,ecc,argper,isw7)
            if(icnJ.ne.430)call phasefold(Ndatamax,isavNJ,savxdataJ,
     @        savydataJ,saverrJ,Nbin,NdataJ,xdataJ,ydataJ,errJ,
     @        period,T0,ikeep,ecc,argper,isw7)
            if(icnH.ne.430)call phasefold(Ndatamax,isavNH,savxdataH,
     @        savydataH,saverrH,Nbin,NdataH,xdataH,ydataH,errH,
     @        period,T0,ikeep,ecc,argper,isw7)
            if(icnK.ne.430)call phasefold(Ndatamax,isavNK,savxdataK,
     @        savydataK,saverrK,Nbin,NdataK,xdataK,ydataK,errK,
     @        period,T0,ikeep,ecc,argper,isw7)
            if(icnRV1.ne.430)call phasefold(Ndatamax,isavRV1,
     @        savxRV1,savyRV1,saverrRV1,Nbin,NRV1,xRV1,yRV1,errRV1,
     @        period,T0,ikeep,ecc,argper,isw7)
            if(icnRV2.ne.430)call phasefold(Ndatamax,isavRV2,
     $        savxRV2,savyRV2,saverrRV2,Nbin,NRV2,xRV2,yRV2,errRV2,
     @        period,T0,ikeep,ecc,argper,isw7)
          endif
c
c   Get the chi^2 values
c
          call initchi(chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,
     &        chisqH,chisqK,chisqRV1,chisqRV2,chisqRV3)
c
          if(icnU.ne.430)call checklcfit(9,Nphase,xmod,ymodU,NdataU,
     @       xdataU,ydataU,errU,chisqU,zeroU,ifixgamma,isw7,resU)
          if(icnB.ne.430)call checklcfit(9,Nphase,xmod,ymodB,NdataB,
     @       xdataB,ydataB,errB,chisqB,zeroB,ifixgamma,isw7,resB)
          if(icnV.ne.430)call checklcfit(9,Nphase,xmod,ymodV,NdataV,
     @       xdataV,ydataV,errV,chisqV,zeroV,ifixgamma,isw7,resV)
          if(icnR.ne.430)call checklcfit(9,Nphase,xmod,ymodR,NdataR,
     @       xdataR,ydataR,errR,chisqR,zeroR,ifixgamma,isw7,resR)
          if(icnI.ne.430)call checklcfit(9,Nphase,xmod,ymodI,NdataI,
     @       xdataI,ydataI,errI,chisqI,zeroI,ifixgamma,isw7,resI)
          if(icnJ.ne.430)call checklcfit(9,Nphase,xmod,ymodJ,NdataJ,
     @       xdataJ,ydataJ,errJ,chisqJ,zeroJ,ifixgamma,isw7,resJ)
          if(icnH.ne.430)call checklcfit(9,Nphase,xmod,ymodH,NdataH,
     @       xdataH,ydataH,errH,chisqH,zeroH,ifixgamma,isw7,resH)
          if(icnK.ne.430)call checklcfit(9,Nphase,xmod,ymodK,NdataK,
     @       xdataK,ydataK,errK,chisqK,zeroK,ifixgamma,isw7,resK)
c
          if((icnRV1.ne.430).and.(icnRV2.ne.430).and.
     @         (icnRV3.eq.430).and.(ifixgamma.ge.2))then
            call checkRVfit(NRVphase,xRVmod,RV1,NRVphase,xRVmod,
     @            RV2,NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     &            chisqRV1,chisqRV2,ggamma1,resRV1,resRV2)
            gamma1=ggamma1+gamma
            gamma2=ggamma1+gamma
            ggamma2=ggamma1
          endif
c
          if((icnRV1.ne.430).and.(icnRV2.ne.430).and.
     @         (icnRV3.ne.430).and.(ifixgamma.ge.2))then
            call checkRVfit3(NRVphase,xRVmod,RV1,
     @            RV2,RV3,NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     @            NRV3,xRV3,yRV3,errRV3,
     &            chisqRV1,chisqRV2,chisqRV3,ggamma1,resRV1,resRV2,
     @            resRV3)
            gamma1=ggamma1+gamma
            gamma2=ggamma1+gamma
            gamma3=ggamma1+gamma
          endif
c
          if((icnRV1.ne.430).and.(icnRV3.ne.430).and.
     @          (icnRV2.eq.430).and.(ifixgamma.ge.2))then
            call checkRVfit(NRVphase,xRVmod,RV1,NRVphase,xRVmod,
     @            RV3,NRV1,xRV1,yRV1,errRV1,NRV3,xRV3,yRV3,errRV3,
     &            chisqRV1,chisqRV3,ggamma1,resRV1,resRV3)
            gamma1=ggamma1+gamma
            gamma3=ggamma1+gamma
            ggamma3=ggamma1
          endif
c
          if(ifixgamma.lt.2)then
            if(icnRV1.ne.430)call checklcfit(0,NRVphase,xRVmod,RV1,
     &        NRV1,xRV1,
     %        yRV1,errRV1,chisqRV1,ggamma1,ifixgamma,isw7,resRV1)
             gamma1=ggamma1+gamma
             if(icnRV2.ne.430)call checklcfit(0,NRVphase,xRVmod,
     %         RV2,NRV2,xRV2,
     %         yRV2,errRV2,chisqRV2,ggamma2,ifixgamma,isw7,resRV2)
             gamma2=ggamma2+gamma
             if(icnRV3.ne.430)call checklcfit(0,NRVphase,xRVmod,
     %         RV3,NRV3,xRV3,
     %         yRV3,errRV3,chisqRV3,ggamma3,ifixgamma,isw7,resRV3)
             gamma3=ggamma3+gamma
          endif
c
          ochi=0.0d0
          ochilr=0.0d0
          if(Nobv.gt.0)call obschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
          if(ifrac.gt.10)call diskchi(Nobv,sobv,obv,eobv,ochidisk,
     @       ochi,compfracs)
c
          if(ilum.gt.10)call lrchi(Nphase,Nmaxphase,
     %          ymods1,ymods2,ymods3,
     %          ymodd,Nobv,sobv,obv,eobv,ochilr,ochi)
c
c   add in the eclipse times chi^2
c
          chitimetot=0.0d0
          if((isw30.ge.3).and.(isw23.ge.1))then
            call timechi(icnarray,Ncycle,Ttimes,Tseps,
     @       Nobscycle,obsTtimes,obsTerr,chitimearray,chitimetot,
     @       Nmaxeclipse)
          endif
          chi1=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+
     @            chisqK+chisqRV1+chisqRV2+ochi+chilimb+chitimetot+
     @            chisqRV3)
c
          isvel3=0
          if(icnRV3.ne.430)isvel3=1
          call writevar(Nvmax,svar,var,fill1,fill2,omega1,omega2,Q,
     @          finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,rLx,
     @          separ,isvel1,gamma1,isvel2,gamma2,t3,g3,SA3,ecc,argper,
     @          pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,
     @          alb2,dwavex,dwavey,primmass,primK,primrad,ratrad,frac1,
     @          frac2,ecosw,temprat,bigI,bigbeta,powercoeff,density,
     @          Tconj,beam1,beam2,contam,ocose,osine,isw29,tertperiod,
     @          tertt0,tertecos,tertesin,tertincl,tertOmega,tertQ,
     @          Tgrav1,Tgrav2,tertconj,omegadot,contamS0,contamS1,
     @          contamS2,contamS3,
     @          P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @       P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @       P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @       P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @       P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @       P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @       P6ratrad,
     @       P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
     @       P7ratrad,
     @       P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @       P8ratrad,gamma3,isvel3,sw72,sw73)
c
          if(sw6.ge.1.0d0)then
            call printamoebamed(ichilabel,chi1)
          else
            call printamoebachi(ichilabel,chi1)
          endif
c
          call writechi(icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %           icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,
     &           chisqRV1,chisqRV2,icnRV3,chisqRV3)              
          if(Nobv.gt.0)call wobschi(Nobv,sobv,obv,eobv,obsparm,ochi)
c
          if(ifrac.gt.10)call wdiskobschi(ochidisk,Nobv,sobv)
          if(ilum.gt.10)call wlrobschi(ochilr)
          if((isw30.ge.3).and.(isw23.ge.1))then
            call writetimechi(icnarray,chitimearray)
          endif
c
          if(ibest.eq.0)then
            if((thresh.le.0.0d0).or.((thresh.gt.0.d0).and.
     @          (chi1.lt.small+thresh)))then
              call chiline(i16,chi1,ochidisk,
     @          icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     %          icnRV1,icnRV2,chisqU,chisqB,chisqV,chisqR,chisqI,
     &          chisqJ,chisqH,chisqK,
     &          chisqRV1,chisqRV2,Nobv,sobv,obv,eobv,obsparm,ochilr)
             endif
          endif 
c
          call newwritevar(Nvmax,svar,var,fill1,fill2,omega1,omega2,
     @          Q,finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,rLx,
     @          separ,isvel1,gamma1,isvel2,gamma2,t3,g3,SA3,ecc,argper,
     @          pshift,spot1parm,spot2parm,spotdparm,period,T0,chi1,i16,
     @          line,alb1,alb2,dwavex,dwavey,primmass,primK,primrad,
     @          ratrad,frac1,frac2,ecosw,temprat,bigI,bigbeta,
     @          powercoeff,density,sw5,Tconj,beam1,beam2,contam,ocose,
     @          osine,isw29,tertperiod,tertt0,tertecos,tertesin,
     @          tertincl,tertOmega,tertQ,Tgrav1,Tgrav2,tertconj,
     @          omegadot,contamS0,contamS1,contamS2,contamS3,
     @          P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @       P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @       P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @       P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @       P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @       P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @       P6ratrad,
     @       P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
     @       P7ratrad,
     @       P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @       P8ratrad,gamma3,isvel3,sw72,sw73)
c
!$omp critical
          if(ibest.eq.0)then
            cccc=chi1
            if(cccc.gt.999999999999.0d0)cccc=999999999999.9d0
            lll=lnblnk(line)
            ljk=lnblnk(parmstring)
            lkk=lnblnk(planetparm)
            ljj=lnblnk(dynparm)
            if((thresh.le.0.0d0).or.((thresh.gt.0.d0).and.
     @            (chi1.lt.small+thresh)))then
              write(46,4646)i16,cccc,parmstring(1:ljk)
              if(isw30.gt.0)write(48,5746)i16,cccc,planetparm(1:lkk)
              if(isw30.ge.3)write(49,5746)i16,cccc,dynparm(1:ljj)
              write(45,5555)line(1:lll)
              if(isw24.ge.1)write(47,5566)i16,cccc,(compfracs(ijk,1),
     @           ijk=1,8),(compfracs(ijk,2),ijk=1,8),
     @           (compfracs(ijk,3),ijk=1,8)
            endif

 4646       format(i7,1x,f18.5,1x,a)
 5746       format(i7,1x,f18.5,1x,a)
c
c
 5555       format(a)
c
5566        format(i5,1x,f18.4,1x,24(1pe13.6,1x))
          endif
!$omp end critical

c
          return
          end
c
c  #!!!%(@^%#%*(&(#@#*$^@^$##*#@&(%)
c
          subroutine loadtimes(icnarray,Nobscycle,obsTtimes,obsTerr,
     @       NRV3,xRV3,yRV3,errRV3,Ndatamax,icnRV3,Nmaxeclipse)
c
c   This routine will load the observed eclipse times into the
c   arrays.
c          
          implicit double precision(a-h,o-z)
c
          dimension Nobscycle(34),obsTtimes(34,Nmaxeclipse)
          dimension obsTerr(34,Nmaxeclipse)
          dimension icnarray(34)
          dimension xRV3(Ndatamax),yRV3(Ndatamax),errRV3(Ndatamax)
c
          character*40 filein
          character*1 bell
c
          bell=char(7)
          ios=0
          ios1=0
          ios2=0
c
          open(unit=20,file='ELCeclipsetimes.opt',status='old',
     @       err=999,iostat=ios)
          do i=1,34
            read(20,100,err=444,iostat=ios2)filein
            icnarray(i)=icnvrt(filein(1:2))
            if(icnarray(i).ne.430)then
              open(unit=21,file=filein,status='old',err=555,iostat=ios1)
              do j=1,10000
                read(21,*,end=15)crap,obsTtimes(i,j),obsTerr(i,j)
              enddo
 15           Nobscycle(i)=j-1
              close(21)
            endif
          enddo
          ios3=0
          read(20,100,err=444,iostat=ios3)filein
          icnRV3=icnvrt(filein(1:2))
          if(icnRV3.ne.430)then
            open(unit=21,file=filein,status='old',err=555,iostat=ios1)
            do j=1,10000
              read(21,*,end=69)xRV3(j),yRV3(j),errRV3(j)
            enddo
 69         close(21)
            NRV3=j-1
          endif
c
 100      format(a40)
            

 999      if(ios.ne.0)then
            write(*,998)bell
            stop
          endif
 998      format(a1,'Error:  ELCeclipsetimes.opt not found')
c
 555      if(ios1.ne.0)then
            write(*,554)bell
            stop
          endif
 554      format(a1,'Error:  data file not found (',a40,')')
c
 444      if(ios2.ne.0)then
            write(*,443)bell
            stop
          endif
 443      format(a1,'Error:  error reading ELCeclipsetimes.opt')
          close(20)
          return
          end
c
c  #$!%@^&*!&!%!!#!)!***&%^%!@!~!_!+!_{}{!(!*(!<>!?U@@
c
          subroutine timechi(icnarray,Ncycle,Ttimes,Tseps,
     @       Nobscycle,obsTtimes,obsTerr,chitimearray,chitimetot,
     @       Nmaxeclipse)
c
c  will compute the chi^2 contributions from the eclipse
c  times
c
          implicit double precision (a-h,o-z)
c
          dimension icnarray(34),Ncycle(34),Ttimes(34,Nmaxeclipse),
     @       Tseps(34,Nmaxeclipse),Nobscycle(34),
     @       obsTtimes(34,Nmaxeclipse),
     @       obsTerr(34,Nmaxeclipse),chitimearray(34)

c          
          do i=1,34
            chitimearray(i)=0.0d0
          enddo
          chitimetot=0.0d0
c
c   go through each data set and match measured eclipses with
c   model ones
c
          do 10 i=1,34
            if(icnarray(i).ne.430)then
              Ngood=0
              do j=1,Ncycle(i)
                if(Tseps(i,j).le.1.0d0)Ngood=Ngood+1
              enddo
              if((Ngood.eq.0).and.(Ncycle(i).gt.0))then
                chitimearray(i)=1.0d20
                go to 10
              endif
              if(Ncycle(i).le.0)then
                chitimearray(i)=1.0d20
                go to 10
              endif
              ssum=0.0d0
              do j=1,Nobscycle(i)
                chismall=1.0d99
                chisave=chismall
                do k=1,Ncycle(i)
                  if(Tseps(i,k).le.1.0d0)then
                    chi=(Ttimes(i,k)-obsTtimes(i,j))/obsTerr(i,j)
                    chi=chi*chi
                    if(chi.lt.chismall)then
                      chismall=chi
                      chisave=chi
                    endif
                  endif
                enddo
                ssum=ssum+chisave
              enddo
              chitimearray(i)=ssum
            endif
 10       continue
c
          chitimetot=0.0d0
          do i=1,34
            chitimetot=chitimetot+chitimearray(i)
          enddo
c
          return
          end
c
c   @#&!^@*((#*#&^<>JSGS{WIW][uouw$@$r^!*()@&^%@)&@%!~``@@@@@
c
          subroutine writetimechi(icnarray,chitimearray)
c
          implicit double precision(a-h,o-z)
c
          dimension icnarray(34),chitimearray(34)
c
          character*80 line
          character*30 outstring
          character*10 instringmed(34),instringchi(34)
c
          common /medblock/ rmed 
c
          instringmed(1)=' med_ecl01'
          instringmed(2)=' med_ecl02'
          instringmed(3)=' med_ecl03'
          instringmed(4)=' med_ecl04'
          instringmed(5)=' med_ecl05'
          instringmed(6)=' med_ecl06'
          instringmed(7)=' med_ecl07'
          instringmed(8)=' med_ecl08'
          instringmed(9)=' med_ecl09'
          instringmed(10)=' med_ecl10'   
          instringmed(11)=' med_ecl11'   
          instringmed(12)=' med_ecl12'   
          instringmed(13)=' med_ecl13'   
          instringmed(14)=' med_ecl14'   
          instringmed(15)=' med_ecl15'   
          instringmed(16)=' med_ecl16'   
          instringmed(17)=' med_ecl17'   
          instringmed(18)=' med_ecl18'   
          instringmed(19)=' med_ecl19'   
          instringmed(20)=' med_ecl20'   
          instringmed(21)=' med_ecl21'   
          instringmed(22)=' med_ecl22'   
          instringmed(23)=' med_ecl23'   
          instringmed(24)=' med_ecl24'   
          instringmed(25)=' med_ecl25'   
          instringmed(26)=' med_ecl26'   
          instringmed(27)=' med_ecl27'   
          instringmed(28)=' med_ecl28'   
          instringmed(29)=' med_ecl29'   
          instringmed(30)=' med_ecl30'   
          instringmed(31)=' med_ecl31'   
          instringmed(32)=' med_ecl32'   
          instringmed(33)=' med_ecl33'   
          instringmed(34)=' med_ecl34'   
c
          instringchi(1)=' chi_ecl01'
          instringchi(2)=' chi_ecl02'
          instringchi(3)=' chi_ecl03'
          instringchi(4)=' chi_ecl04'
          instringchi(5)=' chi_ecl05'
          instringchi(6)=' chi_ecl06'
          instringchi(7)=' chi_ecl07'
          instringchi(8)=' chi_ecl08'
          instringchi(9)=' chi_ecl09'
          instringchi(10)=' chi_ecl10'   
          instringchi(11)=' chi_ecl11'   
          instringchi(12)=' chi_ecl12'   
          instringchi(13)=' chi_ecl13'   
          instringchi(14)=' chi_ecl14'   
          instringchi(15)=' chi_ecl15'   
          instringchi(16)=' chi_ecl16'   
          instringchi(17)=' chi_ecl17'   
          instringchi(18)=' chi_ecl18'   
          instringchi(19)=' chi_ecl19'   
          instringchi(20)=' chi_ecl20'   
          instringchi(21)=' chi_ecl21'   
          instringchi(22)=' chi_ecl22'   
          instringchi(23)=' chi_ecl23'   
          instringchi(24)=' chi_ecl24'   
          instringchi(25)=' chi_ecl25'   
          instringchi(26)=' chi_ecl26'   
          instringchi(27)=' chi_ecl27'   
          instringchi(28)=' chi_ecl28'   
          instringchi(29)=' chi_ecl29'   
          instringchi(30)=' chi_ecl30'   
          instringchi(31)=' chi_ecl31'   
          instringchi(32)=' chi_ecl32'   
          instringchi(33)=' chi_ecl33'   
          instringchi(34)=' chi_ecl34'   
c

          iline=0
          iwrite=0
          
          line='  -->'
          ilength=5
c
          do ii=1,34
            if(icnarray(ii).ne.430)then
              iline=iline+1
              iwrite=0
              if(rmed.ge.1.0d0)then
                call chistring(instringmed(ii),chitimearray(ii),
     @             outstring,lll)
              else
                call chistring(instringchi(ii),chitimearray(ii),
     @             outstring,lll)
              endif
              K=lnblnk(line)
              line=line(1:K)//outstring
              ilength=ilength+lll
              if(ilength.gt.48)then
                iwrite=100
                write(*,500)line
                ilength=5
                line='  -->'
              endif
            endif
          enddo
c
          if(ilength.gt.5.and.iwrite.eq.0)write(*,500)line
c
 500      format(a80)
c
          return
          end
c
c
c    $@$@$%%&%&%@%%%%&&&&&&&&&%@@@@@@@@@$%$@%$@$@$%@$%@$%
c
          subroutine istring(instring,ichi,outstring,ilength)
c
          implicit double precision (a-h,o-z)
c 
          character*(*) instring,outstring
c
c   1 digit positive
c
          if((ichi.ge.0).and.(ichi.lt.10))then
              write(outstring,100)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 100        format(a,'=',i1)
c
c   2 digits positive
c
          if((ichi.ge.10).and.(ichi.lt.100))then
              write(outstring,101)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 101        format(a,'=',i2)
c
c   3 digits positive
c
          if((ichi.ge.100).and.(ichi.lt.1000))then
              write(outstring,102)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 102        format(a,'=',i3)
c
c   4 digits positive
c
          if((ichi.ge.1000).and.(ichi.lt.10000))then
              write(outstring,103)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 103        format(a,'=',i4)
c
c   5 digits positive
c
          if((ichi.ge.10000).and.(ichi.lt.100000))then
              write(outstring,104)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 104        format(a,'=',i5)
c
c   6 digits positive
c
          if((ichi.ge.100000).and.(ichi.lt.1000000))then
              write(outstring,105)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 105        format(a,'=',i6)
c
c   7 digits positive
c
          if((ichi.ge.1000000).and.(ichi.lt.10000000))then
              write(outstring,106)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 106      format(a,'=',i7)
c
c   8 digits positive
c
          if((ichi.ge.10000000).and.(ichi.lt.100000000))then
              write(outstring,107)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 107      format(a,'=',i8)
c
c   9 digits positive
c
          if((ichi.ge.1000000000))then
              itemp=999999999
              write(outstring,108)instring,itemp
              ilength=lnblnk(outstring)
              return
            endif
 108      format(a,'=',i9)
c
c   1 digit negative
c
          if((ichi.gt.-10).and.(ichi.le.0))then
              write(outstring,200)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 200        format(a,'=',i2)
c
c   2 digits negative
c
          if((ichi.gt.-100).and.(ichi.le.-10))then
              write(outstring,201)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 201        format(a,'=',i3)
c
c   3 digits negative
c
          if((ichi.gt.-1000).and.(ichi.le.-100))then
              write(outstring,202)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 202        format(a,'=',i4)
c
c   4 digits negative
c
          if((ichi.gt.-10000).and.(ichi.le.-1000))then
              write(outstring,203)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 203        format(a,'=',i5)
c
c   5 digits negative
c
          if((ichi.gt.-100000).and.(ichi.le.-10000))then
              write(outstring,204)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 204        format(a,'=',i6)
c
c   6 digits negative
c
          if((ichi.gt.-1000000).and.(ichi.le.-100000))then
              write(outstring,205)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 205        format(a,'=',i7)
c
c   7 digits negative
c
          if((ichi.gt.-10000000).and.(ichi.le.-1000000))then
              write(outstring,206)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 206        format(a,'=',i8)
c
c   8 digits negative
c
          if((ichi.gt.-100000000).and.(ichi.le.-10000000))then
              write(outstring,207)instring,ichi
              ilength=lnblnk(outstring)
              return
            endif
 207        format(a,'=',i9)
c
c   9 digits negative
c
          if((ichi.le.-100000000))then
              itemp=-999999999
              write(outstring,209)instring,itemp
              ilength=lnblnk(outstring)
              return
            endif
 209        format(a,'=',i10)

c
c
c
c
          return
          end
c
c  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine printgriditer(kmodel,instring1,iter,instring2,vs,
     @           isw29)
c
          implicit real*8 (a-h,o-z)

          character*(*)instring1,instring2
          character*(80) outstring1,outstring2,outstring3
          character*(*) vs
c
          if(kmodel.lt.0)return
          if(kmodel.ge.1000000)return
c
          l1=lnblnk(instring1)
          l1=l1+1
          if(kmodel.lt.10)write(outstring1,100)instring1(1:l1),kmodel
          if((kmodel.ge.10).and.(kmodel.lt.100))then
            write(outstring1,101)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.100).and.(kmodel.lt.1000))then
            write(outstring1,102)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.1000).and.(kmodel.lt.10000))then
            write(outstring1,103)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.10000).and.(kmodel.lt.100000))then
            write(outstring1,104)instring1(1:l1),kmodel
          endif
          if((kmodel.ge.100000).and.(kmodel.lt.1000000))then
            write(outstring1,105)instring1(1:l1),kmodel
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
c  HERE
c
c
          kk=icnvrt(vs(1:2))
c
          if((kk.eq.450).and.(isw29.gt.0))then          !ecos(omega), tag oc
            outstring3=' ecos(omega)'
          endif
c
          if((kk.eq.466).and.(isw29.gt.0))then
              outstring3=' esin(omega)'
          endif
c
c  tidal apsidal constants
c
          if(kk.eq.1016)then  !rk1, tag a1
            outstring3=' rk1'
          endif

          if(kk.eq.1016)then  !rk2, tag a2
            outstring3=' rk2'
          endif
c
c   planet 2 parameters
c
          if(kk.eq.649)then  !P2tconj, tag uj
            outstring3=' P2tconj'
          endif
c
          if(kk.eq.659)then  !P2period, tag ut
            outstring3=' P2period'
          endif
c
          if(kk.eq.660)then  !P2T0, tag uu
            outstring3=' P2T0'
          endif
c
          if(kk.eq.661)then  !P2ecos, tag uv
            outstring3=' P2ecos'
          endif
c
          if(kk.eq.662)then  !P2esin, tag uw
            outstring3=' P2esin'
          endif
c
          if(kk.eq.663)then  !P2incl, tag ux
            outstring3=' P2incl'
          endif
c
          if(kk.eq.664)then  !P2Omega, tag uy
            outstring3=' P2Omega'
          endif
c
          if(kk.eq.665)then  !P2Q, tag uz
            outstring3=' P2Q'
          endif
c
          if(kk.eq.641)then  !P2artrad, tag ub
            outstring3=' P2ratrad'
          endif
c
c   planet 3 parameters
c
          if(kk.eq.681)then  !P3tconj, tag vj
            outstring3=' P3tconj'
          endif
c
          if(kk.eq.691)then  !P3period, tag vt
            outstring3=' P3period'
          endif
c
          if(kk.eq.692)then  !P3T0, tag vu
            outstring3=' P3T0'
          endif
c
          if(kk.eq.693)then  !P3ecos, tag vv
            outstring3=' P3ecos'
          endif
c
          if(kk.eq.694)then  !P3esin, tag vw
            outstring3=' P3esin'
          endif
c
          if(kk.eq.695)then  !P3incl, tag vx
            outstring3=' P3incl'
          endif
c
          if(kk.eq.696)then  !P3Omega, tag vy
            outstring3=' P3Omega'
          endif
c
          if(kk.eq.697)then  !P3Q, tag vz
            outstring3=' P3Q'
          endif
c
          if(kk.eq.673)then  !P3artrad, tag vb
            outstring3=' P3ratrad'
          endif
c
c   planet 4 parameters
c
          if(kk.eq.713)then  !P4tconj, tag wj
            outstring3=' P4tconj'
          endif
c
          if(kk.eq.723)then  !P4period, tag wt
            outstring3=' P4period'
          endif
c
          if(kk.eq.724)then  !P4T0, tag wu
            outstring3=' P4T0'
          endif
c
          if(kk.eq.725)then  !P4ecos, tag wv
            outstring3=' P4ecos'
          endif
c
          if(kk.eq.726)then  !P4esin, tag ww
            outstring3=' P4esin'
          endif
c
          if(kk.eq.727)then  !P4incl, tag wx
            outstring3=' P4incl'
          endif
c
          if(kk.eq.728)then  !P4Omega, tag wy
            outstring3=' P4Omega'
          endif
c
          if(kk.eq.729)then  !P4Q, tag wz
            outstring3=' P4Q'
          endif
c
          if(kk.eq.705)then  !P4ratrad, tag wb
            outstring3=' P4ratrad'
          endif
c
c   planet 5 parameters
c
          if(kk.eq.745)then  !P5tconj, tag xj
            outstring3=' P5tconj'
          endif
c
          if(kk.eq.755)then  !P5period, tag xt
            outstring3=' P5period'
          endif
c
          if(kk.eq.756)then  !P5T0, tag xu
            outstring3=' P5T0'
          endif
c
          if(kk.eq.757)then  !P5ecos, tag xv
            outstring3=' P5ecos'
          endif
c
          if(kk.eq.758)then  !P5esin, tag xw
            outstring3=' P5esin'
          endif
c
          if(kk.eq.759)then  !P5incl, tag xx
            outstring3=' P5incl'
          endif
c
          if(kk.eq.760)then  !P5Omega, tag xy
            outstring3=' P5Omega'
          endif
c
          if(kk.eq.761)then  !P5Q, tag xz
            outstring3=' P5Q'
          endif
c
          if(kk.eq.737)then  !P5ratrad, tag xb
            outstring3=' P5ratrad'
          endif
c
c   planet 6 parameters
c
          if(kk.eq.585)then  !P6tconj, tag sj
            outstring3=' P6tconj'
          endif
c
          if(kk.eq.595)then  !P6period, tag st
            outstring3=' P6period'
          endif
c
          if(kk.eq.596)then  !P6T0, tag su
            outstring3=' P6T0'
          endif
c
          if(kk.eq.597)then  !P6ecos, tag sv
            outstring3=' P6ecos'
          endif
c
          if(kk.eq.598)then  !P6esin, tag sw
            outstring3=' P6esin'
          endif
c
          if(kk.eq.599)then  !P6incl, tag sx
            outstring3=' P6incl'
          endif
c
          if(kk.eq.600)then  !P6Omega, tag sy
            outstring3=' P6Omega'
          endif
c
          if(kk.eq.601)then  !P6Q, tag sz
            outstring3=' P6Q'
          endif
c
          if(kk.eq.577)then  !P6ratrad, tag sb
            outstring3=' P6ratrad'
          endif
c
c   planet 7 parameters
c
          if(kk.eq.233)then  !P7tconj, tag hj
            outstring3=' P7tconj'
          endif
c
          if(kk.eq.243)then  !P7period, tag ht
            outstring3=' P7period'
          endif
c
          if(kk.eq.244)then  !P7T0, tag hu
            outstring3=' P7T0'
          endif
c
          if(kk.eq.245)then  !P7ecos, tag hv
            outstring3=' P7ecos'
          endif
c
          if(kk.eq.246)then  !P7esin, tag hw
            outstring3=' P7esin'
          endif
c
          if(kk.eq.247)then  !P7incl, tag hx
            outstring3=' P7incl'
          endif
c
          if(kk.eq.248)then  !P7Omega, tag hy
            outstring3=' P7Omega'
          endif
c
          if(kk.eq.249)then  !P7Q, tag hz
            outstring3=' P7Q'
          endif
c
          if(kk.eq.225)then  !P7ratrad, tag hb
            outstring3=' P7ratrad'
          endif
c
c   planet 8 parameters
c
          if(kk.eq.329)then  !P8tconj, tag kj
            outstring3=' P8tconj'
          endif
c
          if(kk.eq.339)then  !P8period, tag kt
            outstring3=' P8period'
          endif
c
          if(kk.eq.340)then  !P8T0, tag ku
            outstring3=' P8T0'
          endif
c
          if(kk.eq.341)then  !P8ecos, tag kv
            outstring3=' P8ecos'
          endif
c
          if(kk.eq.342)then  !P8esin, tag kw
            outstring3=' P8esin'
          endif
c
          if(kk.eq.343)then  !P8incl, tag kx
            outstring3=' P8incl'
          endif
c
          if(kk.eq.344)then  !P8Omega, tag ky
            outstring3=' P8Omega'
          endif
c
          if(kk.eq.345)then  !P8Q, tag kz
            outstring3=' P8Q'
          endif
c
          if(kk.eq.321)then  !P8ratrad, tag kb
            outstring3=' P8ratrad'
          endif

c
c

          if(kk.eq.1208)then  !Tgrav1, tag g1
            outstring3=' Tgrav1'
          endif
c
          if(kk.eq.1591)then  !contamS0, tag s0
            outstring3=' contamS0'
          endif
c
          if(kk.eq.1592)then  !contamS1, tag s1
            outstring3=' contamS1'
          endif
c
          if(kk.eq.1593)then  !contamS2, tag s2
            outstring3=' contamS2'
          endif
c
          if(kk.eq.1594)then  !contamS3, tag s3
            outstring3=' contamS3'
          endif
c
          if(kk.eq.110)then  !omegadot, tag do
            outstring3=' omega_dot'
          endif
c
          if(kk.eq.1209)then  !Tgrav2, tag g2
            outstring3=' Tgrav2'
          endif
c
c   November 18, 2012
c
c   Third body variables here
c
          if(kk.eq.627)then          !tertperiod, tag tt
            outstring3=' tertP'
          endif
c
          if(kk.eq.628)then          !tertT0, tag tu
            outstring3=' tertT0'
          endif
c
          if(kk.eq.617)then          !tertT0, tag tj
            outstring3=' tertconj'
          endif
c
          if(kk.eq.629)then          !tertecos, tag tv
            outstring3=' terte_cos'
          endif
c
          if(kk.eq.630)then          !tertecos, tag tw
            outstring3=' terte_sin'
          endif
c
          if(kk.eq.631)then          !tertincl, tag tx
            outstring3=' tertincl'
          endif
c

          if(kk.eq.632)then          !tertOmega, tag ty
            outstring3=' tertOmega'
          endif
c
          if(kk.eq.633)then          !tertQ, tag tz
            outstring3=' tertQ'
          endif
c
          if(kk.eq.78)then          !contam, tag co
            outstring3=' contam'
          endif
c
          if(kk.eq.1144)then          !beam1, tag e1
            outstring3=' beam1'
          endif
c
          if(kk.eq.1145)then
            outstring3=' beam2'
          endif
c
          if(kk.eq.492)then
            outstring3=' primmass'
          endif
c
          if(kk.eq.497)then
            outstring3=' primrad'
          endif
c
          if(kk.eq.612)then    ! temprat
            outstring3=' temprat'
          endif
c
          if(kk.eq.490)then
            outstring3=' primK'
          endif
c
          if(kk.eq.544)then
            outstring3=' ratrad'
          endif
c
          if(kk.eq.1528)then
            outstring3=' frac1'
          endif
c
          if(kk.eq.100)then
            outstring3=' density'
          endif
c
          if(kk.eq.1529)then
            outstring3=' frac2'
          endif
c 
          if(kk.eq.1368)then
            outstring3=' al1'
          endif
c
          if(kk.eq.1369)then
            outstring3=' al2'
          endif
c
          if(kk.eq.111)then   !ecosw, use string dphi
            outstring3=' dphi'
          endif
c

          if(kk.eq.610)then           !Tconj string
            outstring3=' Tconj'
          endif
c
          if(kk.eq.1623)then           !T0 string
            outstring3=' T0'
          endif
c
          if(kk.eq.484)then      !period string
            outstring3=' P'
          endif
c
          if(kk.eq.269)then
            outstring3=' i'
          endif
c
c   UPDATE May 8, 2006
c
c   Add bigI and bigbeta
c
          if(kk.eq.8)then          !bigI
            outstring3=' bigI'
          endif
c
          if(kk.eq.1)then          !bigbeta
            outstring3=' bigbeta'
          endif

c
c
c   UPDATE JULY @7, 2004
c
c   Change Q to string length 14 (f11.8).
c
          if(kk.eq.384)then
            outstring3=' Q'
          endif
c
c   UPDATE AUGUST 2, 2004
c
c   Make fill1 and fill2 string length 14
c
          if(kk.eq.1176)then
            outstring3=' fi1'
          endif
c
          if(kk.eq.1177)then
            outstring3=' fi2'
          endif
c
          if(kk.eq.552)then
            outstring3=' rin'
          endif
c
          if(kk.eq.1464)then
            outstring3=' o1'
          endif
c
          if(kk.eq.1465)then
            outstring3=' o2'
          endif
c
          if(kk.eq.558)then
            outstring3=' rout'
          endif
c
          if(kk.eq.1626)then
            outstring3=' Teff3'
          endif
c
          if(kk.eq.1210)then
            outstring3=' g3'
          endif
c
          if(kk.eq.36)then
            outstring3=' beta'
          endif
c
          if(kk.eq.1624)then
            outstring3=' Teff1'
          endif
c
          if(kk.eq.1625)then
            outstring3=' Teff2'
          endif
c
          if(kk.eq.744)then
            outstring3=' xi'
          endif
c
          if(kk.eq.611)then
            outstring3=' Td'
          endif
c
          if(kk.eq.375)then
            outstring3=' Td'
          endif
c
          if(kk.eq.580)then
            outstring3=' separ'
          endif
c
          if(kk.eq.576)then
            outstring3=' SA3'
          endif
c
          if(kk.eq.498)then
            outstring3=' pshift'
          endif
c
          if(kk.eq.130)then
            outstring3=' e'
          endif
c
          if(kk.eq.17)then
            outstring3=' argper'
          endif
c
          if(kk.eq.1048)then
            outstring3=' TF_spot1_star1'
          endif
c
          if(kk.eq.1049)then
            outstring3=' lat_spot1_star1'
          endif
c
          if(kk.eq.1050)then
            outstring3=' lon_spot1_star1'
          endif
c
          if(kk.eq.1051)then
            outstring3=' rad_spot1_star1'
          endif
c
          if(kk.eq.1052)then
            outstring3=' TF_spot2_star1'
          endif
c
          if(kk.eq.1053)then
            outstring3=' lat_spot2_star1'
          endif
c
          if(kk.eq.1054)then
            outstring3=' lon_spot2_star1'
          endif
c
          if(kk.eq.1055)then
            outstring3=' rad_spot2_star1'
          endif
c
          if(kk.eq.1080)then
            outstring3=' TF_spot1_star2'
          endif
c
          if(kk.eq.1081)then
            outstring3=' lat_spot1_star2'
          endif
c
          if(kk.eq.1082)then
            outstring3=' lon_spot1_star2'
          endif
c
          if(kk.eq.1083)then
            outstring3=' rad_spot1_star2'
          endif
c
          if(kk.eq.1084)then
            outstring3=' TF_spot2_star2'
          endif
c
          if(kk.eq.1085)then
            outstring3=' lat_spot2_star2'
          endif
c
          if(kk.eq.1086)then
            outstring3=' lon_spot2_star2'
          endif
c
          if(kk.eq.1087)then
            outstring3=' rad_spot2_star2'
          endif
c
          if(kk.eq.1112)then
            outstring3=' TF_spot1_disk'
          endif
c
          if(kk.eq.1113)then
            outstring3=' azi_spot1_disk'
          endif
c
          if(kk.eq.1114)then
            outstring3=' cut_spot1_disk'
          endif
c
          if(kk.eq.1115)then
            outstring3=' wid_spot1_disk'
          endif
c
          if(kk.eq.1116)then
            outstring3=' TF_spot2_disk'
          endif
c
          if(kk.eq.1117)then
            outstring3=' azi_spot2_disk'
          endif
c
          if(kk.eq.1118)then
            outstring3=' cut_spot2_disk'
          endif
c
          if(kk.eq.1119)then
            outstring3=' wid_spot2_disk'
          endif
c
c   UPDATE November 6, 2002
c
c   Add the limb darkening coefficients here.
c
c   x-coefficient, star 1
c
          if(kk.eq.1752)then
            outstring3=' x1(U)'
          endif
c
          if(kk.eq.1753)then
            outstring3=' x1(B)'
          endif
c
          if(kk.eq.1754)then
            outstring3=' x1(V)'
          endif
c
          if(kk.eq.1755)then
            outstring3=' x1(R)'
          endif
c
          if(kk.eq.1756)then
            outstring3=' x1(I)'
          endif
c
          if(kk.eq.1757)then
            outstring3=' x1(J)'
          endif
c
          if(kk.eq.1758)then
            outstring3=' x1(H)'
          endif
c
          if(kk.eq.1759)then
            outstring3=' x1(K)'
          endif
c
c   x-coefficient, star 2
c
          if(kk.eq.1816)then
            outstring3=' x2(U)'
          endif
c
          if(kk.eq.1817)then
            outstring3=' x2(B)'
          endif
c
          if(kk.eq.1818)then
            outstring3=' x2(V)'
          endif
c
          if(kk.eq.1819)then
            outstring3=' x2(R)'
          endif
c
          if(kk.eq.1820)then
            outstring3=' x2(I)'
          endif
c
          if(kk.eq.1821)then
            outstring3=' x2(J)'
          endif
c
          if(kk.eq.1822)then
            outstring3=' x2(H)'
          endif
c
          if(kk.eq.1823)then
            outstring3=' x2(K)'
          endif
c
c   y-coefficient, star 1
c
          if(kk.eq.1784)then
            outstring3=' y1(U)'
          endif
c
          if(kk.eq.1785)then
            outstring3=' y1(B)'
          endif
c
          if(kk.eq.1786)then
            outstring3=' y1(V)'
          endif
c
          if(kk.eq.1787)then
            outstring3=' y1(R)'
          endif
c
          if(kk.eq.1788)then
            outstring3=' y1(I)'
          endif
c
          if(kk.eq.1789)then
            outstring3=' y1(J)'
          endif
c
          if(kk.eq.1790)then
            outstring3=' y1(H)'
          endif
c
          if(kk.eq.1791)then
            outstring3=' y1(K)'
          endif
c
c   y-coefficient, star 2
c
          if(kk.eq.1720)then
            outstring3=' y2(U)'
          endif
c
          if(kk.eq.1721)then
            outstring3=' y2(B)'
          endif
c
          if(kk.eq.1722)then
            outstring3=' y2(V)'
          endif
c
          if(kk.eq.1723)then
            outstring3=' y2(R)'
          endif
c
          if(kk.eq.1724)then
            outstring3=' y2(I)'
          endif
c
          if(kk.eq.1725)then
            outstring3=' y2(J)'
          endif
c
          if(kk.eq.1726)then
            outstring3=' y2(H)'
          endif
c
          if(kk.eq.1727)then
            outstring3=' y2(K)'
          endif
c
c   x-coefficient, star 3
c
          if(kk.eq.1400)then
            outstring3=' x3(U)'
          endif
c
          if(kk.eq.1401)then
            outstring3=' x3(B)'
          endif
c
          if(kk.eq.1402)then
            outstring3=' x3(V)'
          endif
c
          if(kk.eq.1403)then
            outstring3=' x3(R)'
          endif
c
          if(kk.eq.1404)then
            outstring3=' x3(I)'
          endif
c
          if(kk.eq.1405)then
            outstring3=' x3(J)'
          endif
c
          if(kk.eq.1406)then
            outstring3=' x3(H)'
          endif
c
          if(kk.eq.1407)then
            outstring3=' x3(K)'
          endif
c
c
c   y-coefficient, star 3
c
          if(kk.eq.1432)then
            outstring3=' y3(U)'
          endif
c
          if(kk.eq.1433)then
            outstring3=' y3(B)'
          endif
c
          if(kk.eq.1434)then
            outstring3=' y3(V)'
          endif
c
          if(kk.eq.1435)then
            outstring3=' y3(R)'
          endif
c
          if(kk.eq.1436)then
            outstring3=' y3(I)'
          endif
c
          if(kk.eq.1437)then
            outstring3=' y3(J)'
          endif
c
          if(kk.eq.1438)then
            outstring3=' y3(H)'
          endif
c
          if(kk.eq.1439)then
            outstring3=' y3(K)'
          endif

          write(*,200)outstring1(1:ll1),outstring2(1:ll2),outstring3
c
c
c    
c
 100      format(a,i1)
 101      format(a,i2)
 102      format(a,i3)
 103      format(a,i4)
 104      format(a,i5)
 105      format(a,i6)
 200      format(a,', ',a,', variable',a)
          return
          end
c
c
c @#$&^*&*$*&@&%#&%#*&^%#^$@%^(&^%@$&@@@@#$&^*&*$*&@&%#&%#*&^%#^$@%^(&^%
c
          subroutine checkRVfit3(Nmodel1,xmodel1,ymodel1,ymodel2,
     @       ymodel3,Ndata1,xdata1,ydata1,err1,Ndata2,xdata2,ydata2,
     @       err2,Ndata3,xdata3,ydata3,err3,chisq1,chisq2,chisq3,zero,
     @       resRV1,resRV2,resRV3)
c
          implicit double precision (a-h,o-z)
c
           dimension xmodel1(Nmodel1),ymodel1(Nmodel1),
     @      xdata1(Ndata1),ydata1(Ndata1),
     %      err1(Ndata1),ymodel2(Nmodel1),ymodel3(Nmodel1),
     @      xdata2(Ndata2),ydata2(Ndata2),err2(Ndata2)
          dimension xdata3(Ndata3),ydata3(Ndata3),err3(Ndata3)
          dimension resRV1(Ndata1),resRV2(Ndata2),resRV3(Ndata3)  
          dimension yinter1(500000),y21(500000),ydummy1(500000)
          dimension yinter2(500000),y22(500000),ydummy2(500000)
          dimension yinter3(500000),y23(500000),ydummy3(500000)
c
c   UPDATE May 27, 2002
c
c   Add this common block:
c
         common /medblock/ rmed 
c
c   Find the maximum and minimum y-values of the data
c
          ymin1=1.0d20
          ymax1=-1.0d20
          do 20 i=1,Ndata1
            if(ydata1(i).lt.ymin1)ymin1=ydata1(i)
            if(ydata1(i).gt.ymax1)ymax1=ydata1(i)
 20       continue
c
          ymin2=1000099.0d0
          ymax2=-1000099.0d0
          do 21 i=1,Ndata2
            if(ydata2(i).lt.ymin2)ymin2=ydata2(i)
            if(ydata2(i).gt.ymax2)ymax2=ydata2(i)
 21       continue
c
          ymin3=1000099.0d0
          ymax3=-1000099.0d0
          do 22 i=1,Ndata3
            if(ydata3(i).lt.ymin3)ymin3=ydata3(i)
            if(ydata3(i).gt.ymax3)ymax3=ydata3(i)
 22      continue

c
c   Interpolate the model so that we have y-values at all observed phases.
c
          call spline(xmodel1,ymodel1,Nmodel1,0.0d0,0.0d0,y21)
          call spline(xmodel1,ymodel2,Nmodel1,0.0d0,0.0d0,y22)
          call spline(xmodel1,ymodel3,Nmodel1,0.0d0,0.0d0,y23)
c
          do 30 i=1,Ndata1
            call splint(xmodel1,ymodel1,y21,Nmodel1,xdata1(i),qqq)
            yinter1(i)=qqq
 30       continue
c
          do 31 i=1,Ndata2
            call splint(xmodel1,ymodel2,y22,Nmodel1,xdata2(i),qqq)
            yinter2(i)=qqq
 31       continue
          do 32 i=1,Ndata3
            call splint(xmodel1,ymodel3,y23,Nmodel1,xdata3(i),qqq)
            yinter3(i)=qqq
 32      continue
c
c   Now find the optimal zero point that will give the lowest chi^2.
c
          call getmean(Ndata1,ydata1,dataave1)
          call getmean(Ndata1,yinter1,rmodelave1)
          call getmean(Ndata2,ydata2,dataave2)
          call getmean(Ndata2,yinter2,rmodelave2)
          call getmean(Ndata3,ydata3,dataave3)
          call getmean(Ndata3,yinter3,rmodelave3)

c          savezero=zero
          zero=(dataave1-rmodelave1)
          step=abs(rmodelave1-dataave1)/1000.0d0
          if(step.lt.0.10d0)step=0.05d0

          chi1=0.0d0
          chi2=0.0d0
          chi3=0.0d0
          chisq1=0.0d0
          chisq2=0.0d0
          chisq3=0.0d0
          small=1.0d35
          zerosmall=1.0d20
c
          do 750 i=1,20
            zero1=zero
            call offset(Ndata1,yinter1,ydummy1,zero1)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi11)
            call offset(Ndata2,yinter2,ydummy2,zero1)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi12)
            call offset(Ndata3,yinter3,ydummy3,zero1)
            call getchi(Ndata3,ydata3,err3,ydummy3,chi13)
            chi1=chi11+chi12+chi13
            if(chi1.lt.small)then
              small=chi1
              zerosmall=zero1
            endif
            fn=0.0d0
            zero2=zero+step
            call offset(Ndata1,yinter1,ydummy1,zero2)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi21)
            call offset(Ndata2,yinter2,ydummy2,zero2)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi22)
            call offset(Ndata3,yinter3,ydummy3,zero2)
            call getchi(Ndata3,ydata3,err3,ydummy3,chi23)
            chi2=chi21+chi22+chi23
            if(chi2.lt.small)then
              small=chi2
              zerosmall=zero2
            endif
            diff=chi1-chi2
            if(diff.gt.0.0d0)go to 5061
            step=-step
c  
            csave=chi1
            chi1=chi2
            chi2=csave
            zsave=zero1
            zero1=zero2
            zero2=zsave
c
 5061       fn=fn+1.0d0
c
            zero3=zero2+step
            call offset(Ndata1,yinter1,ydummy1,zero3)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi31)
            call offset(Ndata2,yinter2,ydummy2,zero3)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi32)
            call offset(Ndata3,yinter3,ydummy3,zero3)
            call getchi(Ndata3,ydata3,err3,ydummy3,chi33)
            chi3=chi31+chi32+chi33
            if(chi3.lt.small)then
              small=chi3
              zerosmall=zero3
            endif
            diff23=chi3-chi2
            if(diff23.lt.0.0d0)then
              chi1=chi2
              chi2=chi3
              zero1=zero2
              zero2=zero3
              zero=zero3
              go to 5061
            endif          
c
c         find the minimum of parabola defined by the last three points
c
            if((chi3-chi2).eq.0.0d0)go to 999
            step=step*(1.0d0/(1.0d0+(chi1-chi2)/(chi3-chi2))+0.5d0)
            zero=zero2-step
            step=step*fn/3.0d0
            call offset(Ndata1,yinter1,ydummy1,zero)
            call getchi(Ndata1,ydata1,err1,ydummy1,chi41)
            call offset(Ndata2,yinter2,ydummy2,zero)
            call getchi(Ndata2,ydata2,err2,ydummy2,chi42)
            call offset(Ndata3,yinter3,ydummy3,zero)
            call getchi(Ndata3,ydata3,err3,ydummy3,chi43)
            chi4=chi41+chi42+chi43
            if(chi4.lt.small)then
              small=chi4
              zerosmall=zero
            endif
c
 750      continue  ! loop over grid searches
c
          continue                  ! come here when delta chi is small
c
 999      zero=zerosmall
c
          call offset(Ndata1,yinter1,ydummy1,zero)
          call getchi(Ndata1,ydata1,err1,ydummy1,chisq1)
          call offset(Ndata2,yinter2,ydummy2,zero)
          call getchi(Ndata2,ydata2,err2,ydummy2,chisq2)
          call offset(Ndata3,yinter3,ydummy3,zero)
          call getchi(Ndata3,ydata3,err3,ydummy3,chisq3)
c
          do ii=1,Ndata1
            resRV1(ii)=ydata1(ii)-ydummy1(ii)
          enddo
          do ii=1,Ndata2
            resRV2(ii)=ydata2(ii)-ydummy2(ii)
          enddo
          do ii=1,Ndata3
            resRV3(ii)=ydata3(ii)-ydummy3(ii)
          enddo
c 
          return
c
          end
c

