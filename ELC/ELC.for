         program ELC
c
c   October 6, 1999
c
c   This is a total rewrite of the Avni code.  
c
          implicit double precision(a-h,o-z)
c
          parameter(Nmaxphase=1500000)    ! was 1500000
          parameter(Nmaxeclipse=1000)
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm to 9.
c
c   UPDATE September 21, 2008
c
c   Change the dimension of obsparm to 11
c
c
c   UPDATE October 10, 2008
c
c   change the dimensions of obsparm, obv, sobv, and eobv to 11
c
          dimension obsparm(19)
          dimension xmod(Nmaxphase),ymodU(Nmaxphase),ymodB(Nmaxphase),
     $      ymodV(Nmaxphase),ymodR(Nmaxphase),ymodI(Nmaxphase),
     $      ymodJ(Nmaxphase),ymodH(Nmaxphase),ymodK(Nmaxphase),
     &      ymods1(Nmaxphase),ymods2(Nmaxphase),ymodd(Nmaxphase),
     &      RV1(Nmaxphase),RV2(Nmaxphase),ymods3(Nmaxphase),
     $      drv1(Nmaxphase),drv2(Nmaxphase),RV3(Nmaxphase)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,3),
     @       dwavey(8,3)
          dimension xRVmod(Nmaxphase)
c
c   RVG BUG ALERT  June 12, 2001
c
c   Dimension the variables needed for the model atmospheres here.
c
          parameter (maxlines=1300,maxmu=115)  ! was 960
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
c
c
c   RVG BUG ALERT   May 8, 2001
c
c   Define the variables needed for spots
c
          dimension spotdparm(2,4),spot1parm(2,4),spot2parm(2,4)
c
c
c   UPDATE JULY 7, 2004
c
c   This array is needed for the sub-random Sobel sequence, used in
c   place of ran9
c
          dimension xsob(2)
c
c   UPDATE May 8, 2006
c
c   Add this array for power-law limb darkening coefficients.  
c   8=filter index, 9=coefficient c1,c2,c3 ...
c
          dimension powercoeff(8,9)
c
c
c   UPDATE January 12, 2009
c
c   make fracs fracs1, fracs2, .... fracs8 and put them in the
c   argument of subroutine lightcurve.  In this was one can use
c   Nmaxphase in the dimension statement
c
          dimension fracs1(Nmaxphase,4),fracs2(Nmaxphase,4)
          dimension fracs3(Nmaxphase,4),fracs4(Nmaxphase,4)
          dimension fracs5(Nmaxphase,4),fracs6(Nmaxphase,4)
          dimension fracs7(Nmaxphase,4),fracs8(Nmaxphase,4)
          dimension Ncycle(34),Ttimes(34,Nmaxeclipse),
     @       Tseps(34,Nmaxeclipse),Tdur1(34,Nmaxeclipse)
          dimension compfracs(8,3),Tdur2(34,Nmaxeclipse)
c
          dimension gaplow(9999),gaphigh(9999)
          dimension xSC(9999),ySC(9999)
c
c   NEW BUG ALERT  July 13, 2001
c
c   Add a new character string and common block for a 'parameter string'
c   This string of parameters will be fed to the genetic code to make it
c   easier to compute uncertainties on the physical quantities like mass
c   and radius.
c
c   UPDATE June 14, 2002
c
c   Make the length of parmstring character*237.
c
c   UPDATE October 31, 2002
c
c   Make the length of parmstring 249
c
c   UPDATE October 22, 2008
c
c   Make the length of parmstring 259
c
          character*1000 parmstring
          character*2000 planetparm
          character*1700 dynparm,line
c
c
c   NEW BUG August 2, 2001
c
c   Change 'sw4' to 'T0'
c
c   UPDATE August 10, 2004
c
c   Add the 8 variables below.
c
c   UPDATE May 8, 2006
c
c   Add sw21-sw24, powercoeff below.
c  
c   UPDATE November 6, 2008
c
c   Add sw25-sw34
c
c          common /stringblock/ parmstring,planetparm,dynparm,line
cc
c          common /realblock/ fill1,fill2,omega1,omega2,dphase,Q,finc,
c     @      Teff1,Teff2,Tgrav1,Tgrav2,betarim,rinner,router,tdisk,xi,
c     &      alb1,alb2,rLx,Period,fm,separ,gamma,wave,dbolx,dboly,dwavex,
c     @      dwavey,t3,g3,SA3,density,sw1,sw2,sw3,T0,ecc,argper,pshift,
c     @      sw5,sw6,sw7,sw8,sw9,primmass,primK,primrad,ratrad,frac1,
c     @      frac2,ecosw,temprat,bigI,bigbeta,sw23,sw24,powercoeff,sw25,
c     @      sw26,sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,ocose,
c     @      osine,omegadot,contamS0,contamS1,contamS2,contamS3,sw47,
c     @      sw48,sw49,gaplow,gaphigh,
c     @      P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,P2Q,
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
c          common /intblock/ Nalph1,Nbet1,Nalph2,Nbet2,
c     &       Ntheta,Nradius,Nref,
c     $       idraw,iecheck,iidint,iatm,ism1,
c     $       ilaw,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
c     &       iRVfilt,isw1,isw2,isw3,isw4,
c     %       ikeep,isynch,isw5,isw6,isw7,isw8,isw9,idark1,idark2,isw12,isw13,
c     %       isw21,isw22,isw23,isw24,
c     %       isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,
c     $       NSC
cc
c          common /fracblock/ compfracs
cc
c          common /thirdblock/ tertperiod,tertt0,tertecos,tertesin,
c     @        tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72,sw73
c
          common /realatm/ atmT,atmg,atmmu,atmint1,atmint2,atmint3,
     @     atmint4,atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin

          common /intatm/  Nlines,Nmu,Nalph3,Nbet3,itconj,it1,it2,it3,
     @     it4
c
c
c   Open the parameter file and read all of the parameters. 
c   Pass the parameters to the light curve routines
c   via a common block.  
c
        ifastflag=0   !disable fast genetic mode
c
c   RVG BUG ALERT  May 9, 2001
c
c   Add the spot parameters to the arguments of getinput and recordparm
c
c   UPDATE August 10, 2004
c
c   Add the 8 real variables and 4 integer variables to the list.
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
c
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff to list above.
c
c
c   Load the atmosphere table here, rather than in the subroutine
c
c  If the flag iatm>0, then load the model atmosphere table.
c
          if(iatm.ge.1)then
            call loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &         atmint1,atmint2,atmint3,atmint4,atmint5,
     @         atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin)
          endif
c
c   November 17, 2012
c
c   If isw30 >= 1, then load the third body parameters
c
          if((isw30.ge.1).and.(isw7.ge.2))then
            call getbody3(Nalph3,Nbet3,tertperiod,tertt0,tertecos,
     @         tertesin,tertincl,tertOmega,tertQ,dwavex,dwavey,itconj,
     @         it1,it2,it3,it4,tertconj,tertratrad,hh,sw72,sw73,
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
c   if it1=1, assume the mass ratios are logs
c
            if(it1.eq.1)then    
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
          endif
c
c   UPDATE December 15, 2012
c
c   if isw31 >= 1, load the gap parameters
c
          if(isw31.ge.1)call getgap(isw31,gaplow,gaphigh) 
c
          NSC=0
          call getSC(NSC,xSC,ySC)
c
c   UPDATE JULY 7, 2004
c
c   Initialize the Sobel sequence here.
c
          nnn=-11
          call sobseq(nnn,xsob)
c
c   UPDATE January 30, 2010
c
c   if isw28 is greater than 0, then set the value of T0 based
c   on the time of transit Tconj.
c
          if(isw28.gt.0)then
            call getT0(finc,period,ecc,argper,T0,Tconj)
          endif
c
          call lightcurve(Nphase,Nmaxphase,xmod,ymodU,ymodB,
     $      ymodV,ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     %      ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     &         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,
     @         Ncycle,Ttimes,Tseps,RV3,parmstring,planetparm,dynparm,
     @    line,fill1,fill2,omega1,
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
c   Write out eclipse times if isw23 > 0 and Nbody >= 3
c
          if((isw30.ge.3).and.(isw23.ge.1))then
            call writeeclipse(Ncycle,Ttimes,Tseps,isw30,Nmaxeclipse,
     @        Tdur1,Tdur2)
          endif
c         
          call wlinmod(Nphase,xmod,ymodU,'modelU.linear',isw7)
          call wlinmod(Nphase,xmod,ymodB,'modelB.linear',isw7)
          call wlinmod(Nphase,xmod,ymodV,'modelV.linear',isw7)
          call wlinmod(Nphase,xmod,ymodR,'modelR.linear',isw7)
          call wlinmod(Nphase,xmod,ymodI,'modelI.linear',isw7)
          call wlinmod(Nphase,xmod,ymodJ,'modelJ.linear',isw7)
          call wlinmod(Nphase,xmod,ymodH,'modelH.linear',isw7)
          call wlinmod(Nphase,xmod,ymodK,'modelK.linear',isw7)
c
          zp=80.0d0
          call wmagmod(Nphase,xmod,ymodU,'modelU.mag',isw7,zp)
          call wmagmod(Nphase,xmod,ymodB,'modelB.mag',isw7,zp)
          call wmagmod(Nphase,xmod,ymodV,'modelV.mag',isw7,zp)
          call wmagmod(Nphase,xmod,ymodR,'modelR.mag',isw7,zp)
          call wmagmod(Nphase,xmod,ymodI,'modelI.mag',isw7,zp)
          call wmagmod(Nphase,xmod,ymodJ,'modelJ.mag',isw7,zp)
          call wmagmod(Nphase,xmod,ymodH,'modelH.mag',isw7,zp)
          call wmagmod(Nphase,xmod,ymodK,'modelK.mag',isw7,zp)

          call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
          call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
          call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
          if(iidint.ge.1)call wlinmod(Nphase,xmod,ymodd,'lcdisk.linear',isw7)
c
          call wlinmod(NRVphase,xRVmod,RV1,'star1.RV',isw7)
          call wlinmod(NRVphase,xRVmod,RV2,'star2.RV',isw7)
          if(isw30.ge.3)call wlinmod(NRVphase,xRVmod,RV3,'star3.RV',
     @          isw7)
          call wlinmod(NRVphase,xRVmod,dRV1,'star1.delRV',isw7)
          call wlinmod(NRVphase,xRVmod,dRV2,'star2.delRV',isw7)
c

c          
c   RVG BUG ALERT   May 16, 2001
c
c   Close the input file within the lightcurve subroutine, not here.
c
c          close(2)   ! close the output file
c
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          include 'lcsubs.for'
          include 'dynamicssubs.for'



