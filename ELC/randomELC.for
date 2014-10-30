          program loopELC
c
c   May 10, 2001
c
c   This program will read in the light curves and parameters specified
c   in the 'gridloop.opt' file and generate random sets of parameters and
c   compute the chi^2 for each point.
c
c 
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
          implicit double precision (a-h,o-z)
c
          parameter(Nmaxphase=750000,Ndatamax=200000)
          parameter(Nvmax=60)
          parameter(Nmaxeclipse=1000)
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm, obv, sobv and eobv to 9.
c
c
c   UPDATE September 21, 2008
c
c   Change the dimension of obsparm, obv, sobv, and eobv  to 11
c
c
c   UPDATE October 10, 2008
c
c   change the dimensions of obsparm, obv, sobv, and eobv to 17
c
c
c   UPDATE March 15, 2011
c
c   Make the dimensions of obsparm, obs, sobv, eobv 18
c
          dimension obsparm(19),obv(19),eobv(19)
          character*40 sobv(19)
c
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
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     %      stepsave(Nvmax)
          dimension xRVmod(Nmaxphase)
          character*40 Udatafile,svar(Nvmax),Hdatafile,Kdatafile,
     @        RV1file,RV2file
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          character*40 fakevar(Nvmax),command
          dimension parmarray(400000,Nvmax),chiarr(400000),
     @       indxchi(400000)
c
c   RVG BUG ALERT  June 12, 2001
c
c   Dimension the variables needed for the model atmospheres here.
c
          parameter (maxlines=1300,maxmu=115)
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
c   Define variables for the spots
c
          dimension spotdparm(2,4),spot1parm(2,4),spot2parm(2,4)
c
c   NEW BUG August 2, 2001
c
c   New variables are needed to allow for the fitting of period and T0.
c
          dimension savxdataU(Ndatamax),savydataU(Ndatamax),saverrU(Ndatamax),
     &      savydataB(Ndatamax),saverrB(Ndatamax),savydataV(Ndatamax),
     #      saverrV(Ndatamax),
     &      savydataR(Ndatamax),saverrR(Ndatamax),savydataI(Ndatamax),
     $      saverrI(Ndatamax),
     &      savydataJ(Ndatamax),saverrJ(Ndatamax),savydataH(Ndatamax),
     #      saverrH(Ndatamax),
     $      savydataK(Ndatamax),saverrK(Ndatamax),savyRV1(Ndatamax),
     $      saverrRV1(Ndatamax),
     %      savyRV2(Ndatamax),saverrRV2(Ndatamax),savxdataB(Ndatamax),
     $      savxdataV(Ndatamax),
     &      savxdataR(Ndatamax),savxdataI(Ndatamax),savxdataJ(Ndatamax),
     &      savxdataH(Ndatamax),savxdataK(Ndatamax),savxRV1(Ndatamax),
     #      savxRV2(Ndatamax)
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
c
c   UPDATE January 12, 2009
c
c   make fracs fracs1, fracs2, .... fracs8 and put them in the
c   argument of subroutine lightcurve.  In this was one can use
c   Nmaxphase in the dimension statement
c
          dimension compfracs(8,3),ochidisk(8)
          dimension fracs1(Nmaxphase,4),fracs2(Nmaxphase,4)
          dimension fracs3(Nmaxphase,4),fracs4(Nmaxphase,4)
          dimension fracs5(Nmaxphase,4),fracs6(Nmaxphase,4)
          dimension fracs7(Nmaxphase,4),fracs8(Nmaxphase,4)
          dimension resU(Ndatamax),resB(Ndatamax),resV(Ndatamax)
          dimension resR(Ndatamax),resI(Ndatamax),resJ(Ndatamax)
          dimension resH(Ndatamax),resK(Ndatamax),resRV1(Ndatamax)
          dimension resRV2(Ndatamax),gaplow(9999),gaphigh(9999),
     @       resRV3(Ndatamax)
c
          dimension sss(Nvmax)

          dimension Ncycle(34),Ttimes(34,Nmaxeclipse)
          dimension Tseps(34,Nmaxeclipse),Tdur1(34,Nmaxeclipse)
          dimension Nobscycle(34),obsTtimes(34,Nmaxeclipse)
          dimension obsTerr(34,Nmaxeclipse),Tdur2(34,Nmaxeclipse)
          dimension icnarray(34)
          dimension xSC(9999),ySC(9999)

c                                                 
c          character*600 line   
c                               
          character*1000 parmstring
          character*2000 planetparm
          character*1700 dynparm,line
c          character*600 pstringparm(100000)
          character*30 outstring,outstring1
c          dimension pstringchi(100000)
       
c          character*300 parmstring
c          character*30 outstring1,outstring2,outstringp
c
c   Add the period and T0 to the assignvar argument list, and change 'sw4'
c   to 'T0'
c
c   UPDATE August 10, 2004
c
c   Add the extra variables below (8 total)
c
c   UPDATE May 8, 2006
c
c   Add sw21-sw24, powercoeff below.
c
c   UPDATE November 6, 2008
c
c   Add sw25-sw34 below
c
c          common /stringblock/ parmstring,planetparm,dynparm
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
c     @      Nref,idraw,iecheck,iidint,iatm,ism1,ilaw,icnU,icnB,icnV,
c     @      icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,iRVfilt,isw1,isw2,
c     @      isw3,isw4,ikeep,isynch,isw5,isw6,isw7,isw8,isw9,idark1,
c     @      idark2,isw12,isw13,isw21,isw22,isw23,isw24,isw25,isw26,
c     @      isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,NSC
cc
c          common /fracblock/ compfracs
cc
c          common /thirdblock/ tertperiod,tertt0,tertecos,tertesin,
c     @        tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72,sw73
cc
          common /realatm/ atmT,atmg,atmmu,atmint1,atmint2,atmint3,
     &      atmint4,atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin
c
          common /intatm/  Nlines,Nmu,Nalph3,Nbet3,itconj,it1,it2,it3,
     @      it4
c
          common /medblock/ rmed 
c
          common /ranblock/ idum
c
         ifastflag=0   !disable fast genetic mode
c
c   Open the parameter file and read all of the parameters. 
c   Pass the parameters to the light curve routines
c   via a common block.  
c
c   RVG BUG ALERT  May 9, 2001
c
c   Add the spot parameters to getinput and recordparm.
c
c   UPDATE August 10, 2004
c
c   Add the 8 real and 4 integer variables below.
c
          call getinput(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $      omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     @      rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,Nref,
     %      rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,sw3,T0,
     $      idraw,iecheck,iidint,iatm,ism1,icnU,icnB,icnV,icnR,icnI,
     @      icnJ,icnH,icnK,iRVfilt,isw1,isw2,isw3,isw4,ilaw,wave,dbolx,
     @      dboly,dwavex,dwavey,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     $      ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,spot2parm,
     %      spotdparm,primmass,primK,primrad,ratrad,frac1,frac2,ecosw,
     @      temprat,idark1,idark2,isw12,isw13,isw21,isw22,isw23,isw24,
     @      bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,sw29,
     @      sw30,contam,Tconj,beam1,beam2,isw25,isw26,isw27,isw28,isw29,
     @      isw30,isw31,isw32,isw33,isw34,ocose,osine,omegadot,contamS0,
     @      contamS1,contamS2,contamS3,sw47,sw48,sw49)
c
c   Set the threshold to record models in files
c
          thresh=sw48
c
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff to list above.
c
c
c          call recordparm(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
c     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
c     $       rinner,router,tdisk,xi,
c     &       Ntheta,Nradius,alb1,alb2,Nref,
c     %       rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,sw3,T0,
c     $       idraw,iecheck,iidint,iatm,ism1,
c     %       icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,
c     &       isw2,isw3,isw4,
c     &       ilaw,wave,dbolx,dboly,dwavex,dwavey,
c     $       ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
c     $       ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,spot2parm,
c     %       spotdparm)
c
c   RVG BUG ALERT  June 12, 2001
c
c   Load the atmosphere table here, rather than in the subroutine
c
c  If the flag iatm>0, then load the model atmosphere table.
c
          if(iatm.ge.1)then
            call loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &        atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,
     &        atmint8,Tmax,Tmin,gmax,gmin)
          endif
c
c   Set the threshold to record models in files
c
          thresh=sw48
c
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
          endif
c
c   disable the dynamics output file
c
          it2=0
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
c   NEW BUG August 2, 2001
c
c   Use these flags when fitting for the period and/or T0
c
          itime=isw7
c
c   UPDATE May 27, 2002
c
c   Define this new variable.  If rmed > 1, then do median fitting.
c
          rmed=sw6
c
c   Disable the draw option
c
          idraw=0
c
c   Do we force the gamma to be the input value?
c
          ifixgamma=isw4
c
c   UPDATE April 15, 2002
c
c   Define the values of gamma1 and gamma2
c
          gamma1=gamma
          gamma2=gamma
c
c
c   UPDATE JULY 30, 2004
c
c   set Nbin = 0
c
          Nbin=0
c
c   UPDATE JULY 7, 2004
c
c   Initialize the Sobel sequence here.
c
          nnn=-11
          call sobseq(nnn,xsob)
c
c
c   UPDATE AUGUST 4, 2004
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              savesep=separ
c
c   Fix the inner disk radius at fill2.
c
          if(teff2.gt.0.0)rinner=fill2
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
          call initNdata(NdataU,NdataB,NdataV,NdataR,
     &        NdataI,NdataJ,NdataH,NdataK,
     %        NRV1,NRV2,Nobv,NRV3)
c
c   Initialize the steps.
c
          do 1 i=1,Nvmax
            Nstep(i)=1
            var(i)=0.0d0
            svar(i)='none'
            sss(i)=0.0d0
 1        continue
c
c   Get the parameters for the grid search.
c
          call getloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,
     @      Nvmax,Nvar,svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)
c
          ifrac=-99
          ilum=-99
          do 3111 i=1,Nobv
            icn=icnvrt(sobv(i)(1:2))
            if((icn.eq.104).and.(iidint.ge.1))ifrac=99 ! disk fraction = yes   
            if((icn.eq.1112).and.(iidint.ge.1))ifrac=99 ! disk fraction = yes
            if((icn.eq.1113).and.(iidint.ge.1))ifrac=99 ! disk fraction = yes 
            if((icn.eq.1114).and.(iidint.ge.1))ifrac=99 ! disk fraction = yes 
            if((icn.eq.1115).and.(iidint.ge.1))ifrac=99 ! disk fraction = yes 
            if((icn.eq.1116).and.(iidint.ge.1))ifrac=99 ! disk fraction = yes 
            if((icn.eq.1117).and.(iidint.ge.1))ifrac=99 ! disk fraction = yes 
            if((icn.eq.1118).and.(iidint.ge.1))ifrac=99 ! disk fraction = yes 
            if((icn.eq.1119).and.(iidint.ge.1))ifrac=99 ! disk fraction = yes 
            if((icn.eq.369))ilum=99 ! luminosity ratio = yes                 
 3111     continue
c
c   Figure out which data files were specified.
c
          icnU=icnvrt(Udatafile(1:2))
          icnB=icnvrt(Bdatafile(1:2))
          icnV=icnvrt(Vdatafile(1:2))
          icnR=icnvrt(Rdatafile(1:2))
          icnI=icnvrt(Idatafile(1:2))
          icnJ=icnvrt(Jdatafile(1:2))
          icnH=icnvrt(Hdatafile(1:2))
          icnK=icnvrt(Kdatafile(1:2))
          icnRV1=icnvrt(RV1file(1:2))
          icnRV2=icnvrt(RV2file(1:2))
c
c   icn? = 430 indicates that no file was specified.  Load the files listed.
c

          if(icnU.ne.430)call loaddata(Ndatamax,Udatafile,
     &        NdataU,xdataU,ydataU,errU)
          if(icnB.ne.430)call loaddata(Ndatamax,Bdatafile,
     &        NdataB,xdataB,ydataB,errB)
          if(icnV.ne.430)call loaddata(Ndatamax,Vdatafile,
     &        NdataV,xdataV,ydataV,errV)
          if(icnR.ne.430)call loaddata(Ndatamax,Rdatafile,
     &        NdataR,xdataR,ydataR,errR)
          if(icnI.ne.430)call loaddata(Ndatamax,Idatafile,
     &        NdataI,xdataI,ydataI,errI)
          if(icnJ.ne.430)call loaddata(Ndatamax,Jdatafile,
     &        NdataJ,xdataJ,ydataJ,errJ)
          if(icnH.ne.430)call loaddata(Ndatamax,Hdatafile,
     &        NdataH,xdataH,ydataH,errH)
          if(icnK.ne.430)call loaddata(Ndatamax,Kdatafile,
     &        NdataK,xdataK,ydataK,errK)
          if(icnRV1.ne.430)call loaddata(Ndatamax,RV1file,
     &        NRV1,xRV1,yRV1,errRV1)
          if(icnRV2.ne.430)call loaddata(Ndatamax,RV2file,
     &        NRV2,xRV2,yRV2,errRV2)
c
c   If we are fitting for the period and/or T0, we have to make copies
c   of the arrays.
c
          if(icnU.ne.430)call kopydata(Ndatamax,NdataU,
     $       xdataU,ydataU,errU,isavNU,
     $       savxdataU,savydataU,saverrU)
          if(icnB.ne.430)call kopydata(Ndatamax,NdataB,
     $       xdataB,ydataB,errB,isavNB,
     $       savxdataB,savydataB,saverrB)
          if(icnV.ne.430)call kopydata(Ndatamax,NdataV,
     $       xdataV,ydataV,errV,isavNV,
     $       savxdataV,savydataV,saverrV)
          if(icnI.ne.430)call kopydata(Ndatamax,NdataI,
     $       xdataI,ydataI,errI,isavNI,
     $       savxdataI,savydataI,saverrI)
          if(icnJ.ne.430)call kopydata(Ndatamax,NdataJ,
     $       xdataJ,ydataJ,errJ,isavNJ,
     $       savxdataJ,savydataJ,saverrJ)
          if(icnH.ne.430)call kopydata(Ndatamax,NdataH,
     $       xdataH,ydataH,errH,isavNH,
     $       savxdataH,savydataH,saverrH)
          if(icnK.ne.430)call kopydata(Ndatamax,NdataK,
     $       xdataK,ydataK,errK,isavNK,
     $       savxdataK,savydataK,saverrK)
          if(icnR.ne.430)call kopydata(Ndatamax,NdataR,
     $       xdataR,ydataR,errR,isavNR,
     $       savxdataR,savydataR,saverrR)
          if(icnRV1.ne.430)call kopydata(Ndatamax,NRV1,
     $       xRV1,yRV1,errRV1,isavRV1,
     $       savxRV1,savyRV1,saverrRV1)
          if(icnRV2.ne.430)call kopydata(Ndatamax,NRV2,
     $       xRV2,yRV2,errRV2,isavRV2,
     $       savxRV2,savyRV2,saverrRV2)
c
c   Sort the data files by phase just to be on the safe side.
c
c   UPDATE April 15, 2002
c
c   Put an if-then clause to cover the case when Ndata=1
c
          if((icnU.ne.430).and.(NdataU.gt.1))
     @         call sort3(NdataU,xdataU,ydataU,errU)
          if((icnB.ne.430).and.(NdataB.gt.1))
     @         call sort3(NdataB,xdataB,ydataB,errB)
          if((icnV.ne.430).and.(NdataV.gt.1))
     @         call sort3(NdataV,xdataV,ydataV,errV)
          if((icnR.ne.430).and.(NdataR.gt.1))
     @         call sort3(NdataR,xdataR,ydataR,errR)
          if((icnI.ne.430).and.(NdataI.gt.1))
     @         call sort3(NdataI,xdataI,ydataI,errI)
          if((icnJ.ne.430).and.(NdataJ.gt.1))
     @         call sort3(NdataJ,xdataJ,ydataJ,errJ)
          if((icnH.ne.430).and.(NdataH.gt.1))
     @         call sort3(NdataH,xdataH,ydataH,errH)
          if((icnK.ne.430).and.(NdataK.gt.1))
     @         call sort3(NdataK,xdataK,ydataK,errK)
          if((icnRV1.ne.430).and.(NRV1.gt.1))
     @         call sort3(NRV1,xRV1,yRV1,errRV1)
          if((icnRV2.ne.430).and.(NRV2.gt.1))
     @         call sort3(NRV2,xRV2,yRV2,errRV2)
c
c
c  If we are fitting eclipse times ((isw30.ge.3).and.(isw23.ge.1))  
c  then load the data                                      
c                    
          icnRV3=430
          if((isw30.ge.3).and.(isw23.ge.1))then
            call loadtimes(icnarray,Nobscycle,obsTtimes,obsTerr,
     @          NRV3,xRV3,yRV3,errRV3,Ndatamax,icnRV3,Nmaxeclipse)
          endif
c
          isvel1=0
          isvel2=0
          if(icnRV1.ne.430)isvel1=299
          if(icnRV2.ne.430)isvel2=299
c
          open(unit=55,file='chi.0000',status='unknown')
          open(unit=45,file='generation.0000',status='unknown')
          open(unit=46,file='ELCparm.0000',status='unknown')
          if(isw24.ge.1)open(unit=47,file='ELCratio.0000',
     $      status='unknown')
          if(isw30.ge.1)open(unit=48,file='ELCbody3parm.0000',
     @      status='unknown')
          if(isw30.ge.3)open(unit=49,file='ELCdynparm.0000',
     @      status='unknown')
c
c   Define the random parameter sets.  Nstep(1) is the number of random
c   sets to define.
c
          if(isw32.lt.0)then
            idum=isw32
          else
            idum=-1234567
          endif
          do 99 j=1,Nvar
            vx1=vstart(j)
            vx2=vstep(j)
            vxmult=(vx2-vx1)
            vxadd=vx1
            do 98 i=1,Nstep(1)
              varx=ran1(idum)*vxmult+vxadd
              parmarray(i,j)=varx
 98         continue
 99       continue
c
c   Start the looping here.  
c
          chi1=0.0d0
          small=1.0d44
          Ndattot=NdataU+NdataB+NdataV+NdataR+NdataI+NdataJ+NdataH+NdataK+
     %        NRV1+NRV2+NRV3
c
          Nterms=0
          do 699 i7=1,Nvar
            var(i7)=vstart(i7)     !initialize the variables
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
          write(*,*)'Nterms, Nmodels, Nsave',Nterms,Nstep(1),Nstep(2)
          do 750 i16=1,min(400000,Nstep(1))
c
            do 760 j=1,Nterms
              var(j)=parmarray(i16,j)
 760        continue
c
            chilimb=0.0d0
            ibest=0
            ichilabel=i16
            imod=ichilabel
c
            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,ymodR,
     @         ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,RV1,
     @         RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,fracs1,
     @         fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
     @         chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,
     @         chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,
     @         zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,
     &         xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,ydataR,errR,
     %         zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,resI,NdataJ,
     $         xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,
     @         zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,
     @         xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,
     @         Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,
     @         saveasini,savxdataU,savydataU,saverrU,savxdataB,
     @         savydataB,saverrB,savxdataV,savydataV,saverrV,savxdataR,
     @         savydataR,saverrR,savxdataI,savydataI,saverrI,savxdataJ,
     @         savydataJ,saverrJ,savxdataH,savydataH,saverrH,savxdataK,
     @         savydataK,saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,
     @         savyRV2,saverrRV2,ifrac,ilum,i16,isavNU,isavNB,isavNV,
     @         isavNR,isavNI,isavNJ,isavNH,isavNK,isavRV1,isavRV2,
     @         isvel1,isvel2,Ndatamax,ibest,ifixgamma,savesep,ichilabel,
     @         resRV1,resRV2,thresh,
     @         small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     @         icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
     @    parmstring,planetparm,dynparm,line,fill1,fill2,omega1,
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
c            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
c     @           ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
c     @           ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
c     @           xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
c     &           fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
c     &           chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi1,
c     @           NdataU,xdataU,ydataU,errU,zeroU,resU,NdataB,xdataB,
c     @           ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,errV,zeroV,
c     @           resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
c     @           xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,
c     @           errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,
c     @           NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,
c     @           errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,sobv,
c     @           obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,saveasini,
c     @           savxdataU,savydataU,saverrU,savxdataB,savydataB,
c     @           saverrB,savxdataV,savydataV,saverrV,savxdataR,
c     @           savydataR,saverrR,savxdataI,savydataI,saverrI,
c     @           savxdataJ,savydataJ,saverrJ,savxdataH,savydataH,
c     @           saverrH,savxdataK,savydataK,saverrK,savxRV1,savyRV1,
c     @           saverrRV1,savxRV2,savyRV2,saverrRV2,ifrac,ilum,imod,
c     @           isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,isavNH,
c     @           isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,ibest,
c     @           ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
c     @           small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
c     @           icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,
c     @           NRV3,Nmaxeclipse,Tdur1,Tdur2)
c                                                          
              call printiter1(i16,'model number = ')
              chiarr(i16)=chi1
c
              if(chi1.lt.small)then
                small=chi1
                do mmm=1,Nvmax
                  sss(mmm)=var(mmm)
                enddo
              endif
 750        continue
c 
c   Sort the chiarr and print out parameters for the best Nstep(2) sets.
c
          call indexx(Nstep(1),chiarr,indxchi)
c
          do 666 kk=1,Nstep(2)
c
            do 665 j=1,Nterms
              var(j)=parmarray(indxchi(kk),j)
c              Nstep(j)=2
              stepsave(j)=0.1d0*dabs(vstep(j)-vstart(j))
 665        continue
c
            chilimb=0.0d0
            ibest=0
            ichilabel=kk+i16-1
            imod=ichilabel
c
            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,ymodR,
     @         ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,RV1,
     @         RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,fracs1,
     @         fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
     @         chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,
     @         chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,
     @         zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,
     &         xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,ydataR,errR,
     %         zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,resI,NdataJ,
     $         xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,
     @         zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,
     @         xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,
     @         Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,
     @         saveasini,savxdataU,savydataU,saverrU,savxdataB,
     @         savydataB,saverrB,savxdataV,savydataV,saverrV,savxdataR,
     @         savydataR,saverrR,savxdataI,savydataI,saverrI,savxdataJ,
     @         savydataJ,saverrJ,savxdataH,savydataH,saverrH,savxdataK,
     @         savydataK,saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,
     @         savyRV2,saverrRV2,ifrac,ilum,i16,isavNU,isavNB,isavNV,
     @         isavNR,isavNI,isavNJ,isavNH,isavNK,isavRV1,isavRV2,
     @         isvel1,isvel2,Ndatamax,ibest,ifixgamma,savesep,ichilabel,
     @         resRV1,resRV2,thresh,
     @         small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     @         icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
     @     parmstring,planetparm,dynparm,line,fill1,fill2,omega1,
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
c
c            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
c     @           ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
c     @           ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
c     @           xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
c     &           fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
c     &           chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi1,
c     @           NdataU,xdataU,ydataU,errU,zeroU,resU,NdataB,xdataB,
c     @           ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,errV,zeroV,
c     @           resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
c     @           xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,
c     @           errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,
c     @           NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,
c     @           errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,sobv,
c     @           obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,saveasini,
c     @           savxdataU,savydataU,saverrU,savxdataB,savydataB,
c     @           saverrB,savxdataV,savydataV,saverrV,savxdataR,
c     @           savydataR,saverrR,savxdataI,savydataI,saverrI,
c     @           savxdataJ,savydataJ,saverrJ,savxdataH,savydataH,
c     @           saverrH,savxdataK,savydataK,saverrK,savxRV1,savyRV1,
c     @           saverrRV1,savxRV2,savyRV2,saverrRV2,ifrac,ilum,imod,
c     @           isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,isavNH,
c     @           isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,ibest,
c     @           ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
c     @           small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
c     @           icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,
c     @           NRV3,Nmaxeclipse,Tdur1,Tdur2)
c
             call printiter1(i16+kk-1,'model number = ')
             call printS(chi1)
c
          call writegridout(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,
     &       omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,
     &       betarim,rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,
     &       Nref,rLx,Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,
     @       sw3,T0,idraw,iecheck,iidint,iatm,ism1,icnU,icnB,icnV,icnR,
     $       icnI,icnJ,icnH,icnK,iRVfilt,isw1,isw2,isw3,isw4,ilaw,wave,
     &       dbolx,dboly,dwavex,dwavey,ecc,argper,pshift,sw5,sw6,sw7,
     @       sw8,sw9,ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,
     @       spot2parm,spotdparm,primmass,primK,primrad,ratrad,frac1,
     &       frac2,ecosw,temprat,idark1,idark2,isw12,isw13,isw21,isw22,
     @       isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,
     @       sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,isw25,isw26,
     @       isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,ocose,
     @       osine,omegadot,contamS0,contamS1,contamS2,contamS3,sw47,
     @       sw48,sw49)
c
          if(isw30.gt.0)call writebody3grid(Nalph3,Nbet3,tertperiod,
     @       tertt0,tertecos,tertesin,tertincl,tertOmega,tertQ,dwavex,
     @       dwavey,itconj,it1,it2,it3,it4,tertconj,tertratrad,hh,sw72,
     @       sw73,
     @       P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,P2Q,
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
          call recordloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,
     @      Nvmax,Nvar,svar,var,vstart,stepsave,Nstep,Nobv,sobv,
     @      obv,eobv,vstep)
c
            write(command,101)kk+1000
            call system(command)
            write(command,102)kk+1000
            call system(command)
            write(command,103)kk+1000
            call system(command)

 666      continue
c
c   reset the variables at their best values and print the chi^2
c
          do mmm=1,Nvmax
            var(mmm)=sss(mmm)
          enddo
c
            chilimb=0.0d0
            ibest=99
            ichilabel=i16+kk-1  !+Nstep(2)
            imod=ichilabel
c
            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     @           ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @           ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @           xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     &           fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi1,
     @           NdataU,xdataU,ydataU,errU,zeroU,resU,NdataB,xdataB,
     @           ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,errV,zeroV,
     @           resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
     @           xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,
     @           errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,
     @           NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,
     @           errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,sobv,
     @           obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,saveasini,
     @           savxdataU,savydataU,saverrU,savxdataB,savydataB,
     @           saverrB,savxdataV,savydataV,saverrV,savxdataR,
     @           savydataR,saverrR,savxdataI,savydataI,saverrI,
     @           savxdataJ,savydataJ,saverrJ,savxdataH,savydataH,
     @           saverrH,savxdataK,savydataK,saverrK,savxRV1,savyRV1,
     @           saverrRV1,savxRV2,savyRV2,saverrRV2,ifrac,ilum,imod,
     @           isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,isavNH,
     @           isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,ibest,
     @           ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
     @           small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     @           icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,
     @           NRV3,parmstring,planetparm,dynparm,line,fill1,
     @    fill2,omega1,
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
             call printiter1(imod,'model number = ')
             call printS(chi1)

          write(*,*)' '
c          write(*,1002)chiall
          call printS(chi1)
c
          if((isw30.ge.3).and.(isw23.ge.1))then
            call writeeclipse(Ncycle,Ttimes,Tseps,isw30,Nmaxeclipse,
     @           Tdur1,Tdur2)
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
          call wmagmod(Nphase,xmod,ymodU,'modelU.mag',isw7,zeroU)
          call wmagmod(Nphase,xmod,ymodB,'modelB.mag',isw7,zeroB)
          call wmagmod(Nphase,xmod,ymodV,'modelV.mag',isw7,zeroV)
          call wmagmod(Nphase,xmod,ymodR,'modelR.mag',isw7,zeroR)
          call wmagmod(Nphase,xmod,ymodI,'modelI.mag',isw7,zeroI)
          call wmagmod(Nphase,xmod,ymodJ,'modelJ.mag',isw7,zeroJ)
          call wmagmod(Nphase,xmod,ymodH,'modelH.mag',isw7,zeroH)
          call wmagmod(Nphase,xmod,ymodK,'modelK.mag',isw7,zeroK)
c
          call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
          call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
          call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
          if(iidint.ge.1)call wlinmod(Nphase,xmod,ymodd,'lcdisk.linear',
     @        isw7)
c
c          do ii=1,NRVphase
c            RV1(ii)=RV1(ii)+gam
c            RV2(ii)=RV2(ii)+gam
c          enddo
          call wRVmod(NRVphase,xRVmod,RV1,'star1.RV',isw7,ggamma1)
          call wRVmod(NRVphase,xRVmod,RV2,'star2.RV',isw7,ggamma2)
          if(isw30.ge.3)call wRVmod(NRVphase,xRVmod,RV3,'star3.RV',isw7,
     @       ggamma1)
          call wlinmod(NRVphase,xRVmod,dRV1,'star1.delRV',isw7)
          call wlinmod(NRVphase,xRVmod,dRV2,'star2.delRV',isw7)
c
c   NEW BUG August 2, 2001
c
c   If we are fitting for the period and/or the T0, write the current
c   folded light curves
c
          if(itime.gt.0)then
            if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,ydataU,
     %              errU,'ELCdataU.fold')
            if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,ydataB,
     %              errB,'ELCdataB.fold')
            if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,ydataV,
     %              errV,'ELCdataV.fold')
            if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,ydataR,
     %              errR,'ELCdataR.fold')
            if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,ydataI,
     %              errI,'ELCdataI.fold')
            if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,ydataJ,
     %              errJ,'ELCdataJ.fold')
            if(icnH.ne.430) call wELCdata(Ndatamax,NdataH,xdataH,ydataH,
     %              errH,'ELCdataH.fold')
            if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,ydataK,
     %              errK,'ELCdataK.fold')
            if(icnRV1.ne.430) call wELCdata(Ndatamax,NRV1,xRV1,yRV1,
     %              errRV1,'ELCdataRV1.fold')
            if(icnRV2.ne.430) call wELCdata(Ndatamax,NRV2,xRV2,yRV2,
     %              errRV2,'ELCdataRV2.fold')
            if(icnRV3.ne.430) call wELCdata(Ndatamax,NRV3,xRV3,yRV3,
     %              errRV3,'ELCdataRV3.fold')
          endif
c
            if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,resU,
     %              errU,'ELCresidualsU.fold')
            if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,resB,
     %              errB,'ELCresidualsB.fold')
            if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,resV,
     %              errV,'ELCresidualsV.fold')
            if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,resR,
     %              errR,'ELCresidualsR.fold')
            if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,resI,
     %              errI,'ELCresidualsI.fold')
            if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,resJ,
     %              errJ,'ELCresidualsJ.fold')
            if(icnH.ne.430)call wELCdata(Ndatamax,NdataH,xdataH,resH,
     %              errH,'ELCresidualsH.fold')
            if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,resK,
     %              errK,'ELCresidualsK.fold')
            if(icnRV1.ne.430) call wELCdata(Ndatamax,NRV1,xRV1,resRV1,
     %              errRV1,'ELCresidualsRV1.fold')
            if(icnRV2.ne.430) call wELCdata(Ndatamax,NRV2,xRV2,resRV2,
     %              errRV2,'ELCresidualsRV2.fold')
            if(icnRV3.ne.430) call wELCdata(Ndatamax,NRV3,xRV3,resRV3,
     %              errRV3,'ELCresidualsRV3.fold')
c
c   The output file is closed within the light curve subroutine now.
c
c          close(2)

          close(45)
          close(46)
          close(55)
          if(isw24.ge.1)close(47)
          if(isw30.gt.0)close(48)

 101      format('cp gridELC.opt gridELC.',i4)
 102      format('cp gridELC.inp ELC.',i4)
 103      format('cp gridELCbody3.inp ELCbody3.',i4)

          end
c
c &^&%$*$*$*$^@@^@^$^&&%&&^&%$*$*$*$^@@^@^$^&#$@%%&^&%$*$*$*$^@@^@^$^&
c
c      FUNCTION ran1(idum)
c      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
c      REAL*8 ran1,AM,EPS,RNMX
c      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
c     $     NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=3.e-16,RNMX=1.-EPS)
c      INTEGER j,k,iv(NTAB),iy
c      SAVE iv,iy
c      DATA iv /NTAB*0/, iy /0/
c      if (idum.le.0.or.iy.eq.0) then
c         idum=max(-idum,1)
c         do j=NTAB+8,1,-1
c            k=idum/IQ
c            idum=IA*(idum-k*IQ)-IR*k
c            if (idum.lt.0) idum=idum+IM
c            if (j.le.NTAB) iv(j)=idum
c         enddo
c         iy=iv(1)
c      endif
c      k=idum/IQ
c      idum=IA*(idum-k*IQ)-IR*k
c      if (idum.lt.0) idum=idum+IM
c      j=1+iy/NDIV
c      iy=iv(j)
c      iv(j)=idum
c      ran1=min(AM*iy,RNMX)
c      return
c      END
cc
c     &&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine indexx(n,arr,indx)
c
          integer n,indx(n),m,nstack
          real*8 arr(n)
          parameter(m=7,nstack=50)
          integer i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
          real*8 a
c
          do 11 j=1,n
            indx(j)=j
 11       continue
c
          jstack=0
          l=1
          ir=n
 1        if(ir-l.lt.m)then
            do 13 j=l+1,ir
              indxt=indx(j)
              a=arr(indxt)
              do 12 i=j-1,1,-1
                if(arr(indx(i)).le.a)go to 2
                indx(i+1)=indx(i)
 12           continue
              i=0
 2            indx(i+1)=indxt
 13         continue
            if(jstack.eq.0)return
            ir=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
          else
            k=(l+ir)/2
            itemp=indx(k)
            indx(k)=indx(l+1)
            indx(l+1)=itemp
            if(arr(indx(l+1)).gt.arr(indx(ir)))then
              itemp=indx(l+1)
              indx(l+1)=indx(ir)
              indx(ir)=itemp
            endif
            if(arr(indx(l)).gt.arr(indx(ir)))then
              itemp=indx(l)
              indx(l)=indx(ir)
              indx(ir)=itemp
            endif
            if(arr(indx(l+1)).gt.arr(indx(l)))then
              itemp=indx(l+1)
              indx(l+1)=indx(l)
              indx(l)=itemp
            endif
            i=l+1
            j=ir
            indxt=indx(l)
            a=arr(indxt)
 3          continue
            i=i+1
            if(arr(indx(i)).lt.a)go to 3
 4          continue
            j=j-1
            if(arr(indx(j)).gt.a)go to 4
            if(j.lt.i) go to 5
            itemp=indx(i)
            indx(i)=indx(j)
            indx(j)=itemp
            go to 3
 5          indx(l)=indx(j)
            indx(j)=indxt
            jstack=jstack+2
            if(jstack.gt.nstack)then 
              write(*,*) 'nstack too small in indexx'
              stop
            endif
            if(ir-i+1.ge.j-l)then
              istack(jstack)=ir
              istack(jstack-1)=i
              ir=j-1
            else
              istack(jstack)=j-1
              istack(jstack-1)=l
              l=i
            endif
          endif
          go to 1
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
          include 'dynamicssubs.for'









