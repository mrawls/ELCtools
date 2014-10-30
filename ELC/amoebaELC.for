      program amobaeELC
c
c   April 26, 2001
c
c   This is an optimizer code based on the 'amoeba' program
c   discussed in Numerical Recipes.  This program will call the
c   subroutines contained within 'lcsubs.for' and 'optimizesubs.for'.
c
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
          implicit double precision (a-h,o-z)
c
          parameter(ftol=1.0d-4,maxiter=150)
          parameter(Nmaxphase=750000,Ndatamax=200000)
          parameter(Nvmax=60)
          parameter(Nmaxeclipse=1000)
c
c   UPDATE March 15, 2011
c
c   Make obsparm, obv, eobv, sobv dimension 18
c
          character*40 Udatafile,svar(Nvmax),
     &             Hdatafile,Kdatafile,RV1file,RV2file
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm, obv, sobv, and eobv to 9.
c
c
c   UPDATE September 21, 2008
c
c   Change the dimension of obsparm, obv, sobv, and eobv  to 11
c
c   UPDATE October 10, 2008
c
c   change the dimensions of obsparm, obv, sobv, and eobv to 11
c
          dimension obsparm(19),obv(19),eobv(19)
          character*40 sobv(19)
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
     @      xdataB(Ndatamax),xdataV(Ndatamax),
     &      xdataR(Ndatamax),xdataI(Ndatamax),xdataJ(Ndatamax),
     &      xdataH(Ndatamax),xdataK(Ndatamax),xRV1(Ndatamax),
     @      xRV2(Ndatamax),xRV3(Ndatamax),yRV3(Ndatamax),
     @      errRV3(Ndatamax)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,3),
     @        dwavey(8,3)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     $      drv1(Nmaxphase),drv2(Nmaxphase),vstep1(Nvmax)
          dimension xRVmod(Nmaxphase)
          dimension resU(Ndatamax),resB(Ndatamax),resV(Ndatamax)
          dimension resR(Ndatamax),resI(Ndatamax),resJ(Ndatamax)
          dimension resH(Ndatamax),resK(Ndatamax),resRV1(Ndatamax)
          dimension resRV2(Ndatamax),resRV3(Ndatamax)
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
c   Here is the matrix and vector needed for the amoeba. 
c
          dimension p(Nvmax+1,Nvmax),yval(Nvmax),ptry(Nvmax),psum(Nvmax)
c
c
          dimension spotdparm(2,4),spot1parm(2,4),spot2parm(2,4)
c
c   New variables are needed to allow for the fitting of period and T0.
c
          dimension savxdataU(Ndatamax),savydataU(Ndatamax),
     $      saverrU(Ndatamax),savydataB(Ndatamax),saverrB(Ndatamax),
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
c   UPDATE OCTOBER 10, 2007
c
c   Add this array to keep track of the light curves of the individual
c   components in all band passes.
c
          dimension compfracs(8,3),ochidisk(8)
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
          dimension gaplow(9999),gaphigh(9999)
c
          dimension sss(Nvmax)
c
          character*1000 parmstring
          character*2000 planetparm
          character*1700 dynparm,line
          character*30 outstring1,outstring2,outstringp
c
          dimension Ncycle(34),Ttimes(34,Nmaxeclipse)
          dimension Tseps(34,Nmaxeclipse),Tdur1(34,Nmaxphase)
          dimension Nobscycle(34),obsTtimes(34,Nmaxeclipse)
          dimension obsTerr(34,Nmaxeclipse)
          dimension icnarray(34),Tdur2(34,Nmaxeclipse)
          dimension xSC(9999),ySC(9999)

c          common /stringblock/ parmstring,planetparm,dynparm,line
cc
c          common /realblock/ fill1,fill2,omega1,omega2,dphase,Q,finc,
c     @      Teff1,Teff2,Tgrav1,Tgrav2,betarim,rinner,router,tdisk,xi,
c     @      alb1,alb2,rLx,Period,fm,separ,gamma,wave,dbolx,dboly,dwavex,
c     @      dwavey,t3,g3,SA3,density,sw1,sw2,sw3,T0,ecc,argper,pshift,
c     @      sw5,sw6,sw7,sw8,sw9,primmass,primK,primrad,ratrad,frac1,
c     @      frac2,ecosw,temprat,bigI,bigbeta,sw23,sw24,powercoeff,sw25,
c     @      sw26,sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,ocose,
c     @      osine,sw42,contamS0,contamS1,contamS2,contamS3,sw47,
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
c     %      icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,iRVfilt,isw1,isw2,
c     @      isw3,isw4,ikeep,isynch,isw5,isw6,isw7,isw8,isw9,idark1,
c     &      idark2,isw12,isw13,isw21,isw22,isw23,isw24,isw25,isw26,
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
c   Disable median fitting (if selected by sw6 > 0)
c
          rmed=0.0d0
c
          ifastflag=0   !disable fast genetic mode
c
c   Open the parameter file and read all of the parameters. 
c   Pass the parameters to the light curve routines
c   via a common block.  
c
c   May 9, 2001
c
c   Add the spot parameters to getinput and recordparm.
c
c   UPDATE August 10, 2004
c
c   Add the 8 real variables and 4 integers to the argument list.
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
     &       temprat,idark1,idark2,isw12,isw13,isw21,isw22,isw23,isw24,
     @       bigI,bigbeta,sw23,sw24,powercoeff, sw25,sw26,sw27,sw28,
     @       sw29,sw30,contam,Tconj,beam1,beam2,isw25,isw26,isw27,isw28,
     @       isw29,isw30,isw31,isw32,isw33,isw34,ocose,osine,omegadot,
     @       contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49)
c
c   Set the threshold to record models in files
c
          thresh=sw48
c
c   Load the atmosphere table here, rather than in the subroutine
c
c  If the flag iatm>0, then load the model atmosphere table.
c
          if(iatm.ge.1)then
            call loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &         atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,
     &         atmint7,atmint8,Tmax,Tmin,gmax,gmin)
          endif
c
c   Set the threshold to record models in files
c
          thresh=sw48
c
c   November 17, 2012
c
c   If isw30 >= 1, then load the third body parameters
c
          if((isw30.ge.1).and.(isw7.ge.2))then
            call getbody3(Nalph3,Nbet3,tertperiod,tertt0,tertecos,
     @         tertesin,tertincl,tertOmega,tertQ,dwavex,dwavey,itconj,
     @         it1,it2,it3,it4,tertconj,tertratrad,hh,sw72,sw73,P2tconj,
     @         P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,P2Q,
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
c   Use these flags when fitting for the period and/or T0
c
          itime=isw7
          ifixgamma=isw4

c
          gamma1=gamma
          gamma2=gamma
c
c   Disable the draw option
c
          idraw=0
c
c   Fix the inner disk radius at fill2.
c
          if(teff2.gt.0.0d0)rinner=fill2
c
c   UPDATE AUGUST 4, 2004
c
c   Define a new variable called sepsave, which is equal to the
c   value of separ given in ELC.inp
c
          savesep=separ
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
          call initNdata(NdataU,NdataB,NdataV,NdataR,
     &        NdataI,NdataJ,NdataH,NdataK,
     %        NRV1,NRV2,Nobv,NRV3)
c 
c   Initialize the steps. 
c  
          do  i=1,nvmax
            Nstep(i)=1
            var(i)=0.0d0
            svar(i)='none'
            sss(i)=0.0d0
          enddo
c
c   Get the parameters for the amoeba.
c
          call getloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,
     &      Nvmax,Nvar,svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)
c
c   UPDATE FEBRUARY 5, 2005
c
c   Add also the ability to include the luminosity ratio in the total chi^2
c
          ifrac=-99
          ilum=-99
          do 3111 i=1,Nobv
            vstep1(i)=vstep(i)
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
c  If we are fitting eclipse times ((isw30.ge.3).and.(isw23.ge.1))
c  then load the data
c
          icnRV3=430
          if((isw30.ge.3).and.(isw23.ge.1))then
            call loadtimes(icnarray,Nobscycle,obsTtimes,obsTerr,
     @        NRV3,xRV3,yRV3,errRV3,Ndatamax,icnRV3,Nmaxeclipse)
          endif
          isvel1=0
          isvel2=0
          if(icnRV1.ne.430)isvel1=299
          if(icnRV2.ne.430)isvel2=299
c
          gamma1=gamma
          gamma2=gamma
c                      
          if(isw24.ge.1)open(unit=47,file='ELCratio.0000',
     @      status='unknown')
          if(isw30.ge.1)open(unit=48,file='ELCbody3parm.0000',
     @      status='unknown')
          if(isw30.ge.3)open(unit=49,file='ELCdynparm.0000',
     @      status='unknown')

c
c   We have to load the p and yvar matrix.  The values of p(1,1:Nvar)
c   are simply the variables in the gridloop.opt file.  The value of
c   yval(1) is the chi^2 there.
c
          chi1=0.0d0
          chilimb=0.0d0
          small=1.0d30
          Ndattot=NdataU+NdataB+NdataV+NdataR+NdataI+NdataJ+NdataH+
     @      NdataK+NRV1+NRV2+Nobv+NRV3
c
          Nterms=0
          do 698 i7=1,Nvar
            kkk=icnvrt(svar(i7)(1:2))
            if(kkk.ne.430)then
              Nterms=Nterms+1
            endif
 698      continue

 748      do 699 i7=1,Nvar      !Nterms
            p(1,i7)=vstart(i7)     
            var(i7)=vstart(i7)
 699      continue
c
          write(*,*)'Nterms, Nvar, Ndata ',Nterms,Nvar,Ndattot

          icount=0
 749      if(Nterms.eq.0)go to 69   !no valid variables specified
c
          open(unit=55,file='chi.0000',status='unknown')
          open(unit=45,file='generation.0000',status='unknown')
          open(unit=46,file='ELCparm.0000',status='unknown')
          i16=0
c
          ibest=0
          ichilabel=0
          call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
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
     @      isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,
     @      ibest,ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
     @      small,
     @      Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,icnarray,
     @      RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
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
          yval(1)=chi1
          if(chi1.lt.small)then
            small=chi1
            do mmm=1,Nvmax
              sss(mmm)=var(mmm)
            enddo
          endif
c
          i16=i16+1
          chikkk=chi1
c
c   Now we have to fill out the other Nterms dimension
c 
          do 800 nn=1,Nvar     !Nterms
            kkk=icnvrt(svar(nn)(1:2))
            if(kkk.eq.430)go to 800
c
c   Reset the var array.
c
            do 899 i7=1,Nvar     !Nterms
              var(i7)=vstart(i7)
 899        continue
c
c   Alter the nn'th term
c
            var(nn)=var(nn)+vstep(nn)
c
c   Fill up the row in the p matrix
c
            do 898 i7=1,Nvar   !Nterms
              p(nn+1,i7)=var(i7)
 898        continue
c
            ichilabel=nn
            ibest=0
            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,ymodR,
     @         ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,
     @         RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     @         fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,
     @         chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,
     &         chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,
     $         zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,
     @         xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,ydataR,errR,
     &         zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,resI,NdataJ,
     @         xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,
     @         zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,
     @         xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,
     @         Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,
     @         saveasini,savxdataU,savydataU,saverrU,savxdataB,
     &         savydataB,saverrB,savxdataV,savydataV,saverrV,savxdataR,
     &         savydataR,saverrR,savxdataI,savydataI,saverrI,savxdataJ,
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
            yval(nn+1)=chi1
c
            if(chi1.lt.small)then
              small=chi1
              do mmm=1,Nvmax
                sss(mmm)=var(mmm)
              enddo
            endif
c
            i16=i16+1
c
 800      continue
c
c   Go ahead with the amoeba.
c
 1        iter=0 
          do 12 n=1,Nvar   !Nterms
            sum=0.0d0
            do 11 m=1,Nvar+1   !Nterms+1
              sum=sum+p(m,n)
 11         continue
            psum(n)=sum
 12       continue
c
 2        ilo=1
          if(yval(1).gt.yval(2))then
            ihi=1
            inhi=2
          else
            ihi=2
            inhi=1
          endif
c
          do 13 i=1,Nvar+1      !Nterms+1
            if(yval(i).le.yval(ilo))ilo=i
            if(yval(i).gt.yval(ihi))then
              inhi=ihi
              ihi=i
            else if(yval(i).gt.yval(inhi))then
              if(i.ne.ihi)ihi=i
            endif
 13       continue
c
          rtol=2.0d0*dabs(yval(ihi)-yval(ilo))/(dabs(yval(ihi))+
     @          dabs(yval(ilo)))
c        
          call chistring('chi^2(hi)',yval(ihi),outstring1,lll1)
          call chistring('chi^2(lo)',yval(ilo),outstring2,lll2)
          call pstring('rtol',6,rtol,outstringp,lllp)
          write(*,*)' '
          write(*,457)outstringp(1:lllp),outstring1(1:lll1),
     @          outstring2(1:lll2),ihi,ilo
c
 457      format(a,',',2x,a,',',2x,a,',',2x,2(i2,1x))
c
          if(rtol.lt.ftol)then
            swap=yval(1)
            yval(1)=yval(ilo)
            do 14 n=1,Nvar        !Nterms
              swap=p(1,n)
              p(1,n)=p(ilo,n)
              p(ilo,n)=swap
 14         continue
            go to 751
          endif
c
          if(iter.gt.maxiter)go to 751   ! escape
          iter=iter+2
          if(iter.lt.10)then
            write(*,234)iter
          endif
          if((iter.ge.10).and.(iter.lt.100))then
            write(*,235)iter
          endif
          if((iter.ge.100).and.(iter.lt.1000))then
            write(*,236)iter
          endif
c
 234      format('iter = ',i1)
 235      format('iter = ',i2)
 236      format('iter = ',i3)
c
c   Start a new iteration.  This block is basically the function amotry
c
          fac=-1.0d0
c
          fac1=(1.0d0-fac)/dble(Nterms)
          fac2=fac1-fac
          do 111 j=1,Nvar      !Nterms
            ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
            var(j)=ptry(j)
 111      continue
c
          ichilabel=iter
          ibest=0
          call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
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
     @      isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,
     @      ibest,ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
     @      small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     @      icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
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
          ytry=chi1
          if(chi1.lt.small)then
            small=chi1
            do mmm=1,Nvmax
              sss(mmm)=var(mmm)
            enddo
          endif
c
          i16=i16+1
          if(ytry.lt.yval(ihi))then
            yval(ihi)=ytry
            do 112 j=1,Nvar!     Nterms
              psum(j)=psum(j)-p(ihi,j)+ptry(j)
              p(ihi,j)=ptry(j)
 112        continue
          endif
c
c   This is the end of the function block
c
          if(ytry.le.yval(ilo))then
            fac=2.0d0
c
c   Repeat the above block.
c
            fac1=(1.0d0-fac)/dble(Nterms)
            fac2=fac1-fac
            do 211 j=1,Nvar      !Nterms
              ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
              var(j)=ptry(j)
 211        continue
c
            ichilabel=iter
            ibest=0
            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,ymodR,
     &         ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,RV1,
     @         RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,fracs1,
     %         fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
     @         chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,
     %         chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,
     @         zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,
     @         xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,ydataR,errR,
     @         zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,resI,NdataJ,
     @         xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,
     @         zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,
     @         xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,
     @         Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,
     %         saveasini,savxdataU,savydataU,saverrU,savxdataB,
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
            ynewtry=chi1
            if(chi1.lt.small)then
              small=chi1
              do mmm=1,Nvmax
                sss(mmm)=var(mmm)
              enddo
            endif
c
            i16=i16+1
            if(ynewtry.lt.yval(ihi))then
              yval(ihi)=ynewtry
              do 212 j=1,Nvar    !Nterms
                psum(j)=psum(j)-p(ihi,j)+ptry(j)
                p(ihi,j)=ptry(j)
 212          continue
            endif
c
c   This is the end of the function block
c
            ytry=ynewtry
          else if(ytry.ge.yval(inhi))then
            ysave=yval(ihi)
            fac=0.5d0
c
c   Do another function block.
c
            fac1=(1.0d0-fac)/dble(Nterms)
            fac2=fac1-fac
            do 311 j=1,Nvar    !Nterms
              ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
              var(j)=ptry(j)
 311        continue
c
            ichilabel=iter
            ibest=0
            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,ymodR,
     &         ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,RV1,
     @         RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,fracs1,
     %         fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
     @         chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,
     %         chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,
     @         zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,
     @         xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,ydataR,errR,
     @         zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,resI,NdataJ,
     @         xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,
     @         zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,
     @         xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,
     @         Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,
     %         saveasini,savxdataU,savydataU,saverrU,savxdataB,
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
            ynewtry=chi1
            if(chi1.lt.small)then
              small=chi1
              do mmm=1,Nvmax
                sss(mmm)=var(mmm)
              enddo
            endif
              i16=i16+1
            if(ynewtry.lt.yval(ihi))then
              yval(ihi)=ynewtry
              do 312 j=1,Nvar   !Nterms
                psum(j)=psum(j)-p(ihi,j)+ptry(j)
                p(ihi,j)=ptry(j)
 312          continue
            endif
c
c   This is the end of the function block
c
            if(ynewtry.ge.ysave)then
              do 16 i=1,Nvar+1    !Nterms+1
                if(i.ne.ilo)then
                  do 15 j=1,Nvar         !Nterms
                    psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                    p(i,j)=psum(j)
                    var(j)=psum(j)
 15               continue
c
                ichilabel=i
                ibest=0
                call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     @            ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @            ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     %            xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     &            fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     &            chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi1,
     @            NdataU,xdataU,ydataU,errU,zeroU,resU,NdataB,xdataB,
     @            ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,errV,
     @            zeroV,resV,NdataR,xdataR,ydataR,errR,zeroR,resR,
     @            NdataI,xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,
     @            ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,
     @            zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,
     @            xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,
     @            ggamma2,Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,
     @            Nvmax,svar,var,saveasini,savxdataU,savydataU,saverrU,
     @            savxdataB,savydataB,saverrB,savxdataV,savydataV,
     @            saverrV,savxdataR,savydataR,saverrR,savxdataI,
     &            savydataI,saverrI,savxdataJ,savydataJ,saverrJ,
     @            savxdataH,savydataH,saverrH,savxdataK,savydataK,
     @            saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,savyRV2,
     @            saverrRV2,ifrac,ilum,i16,isavNU,isavNB,isavNV,isavNR,
     @            isavNI,isavNJ,isavNH,isavNK,isavRV1,isavRV2,isvel1,
     @            isvel2,Ndatamax,ibest,ifixgamma,savesep,ichilabel,
     @            resRV1,resRV2,thresh,
     @            small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     @            icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,
     @            NRV3,parmstring,planetparm,dynparm,line,fill1,fill2,
     @    omega1,
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
                  yval(i)=chi1
                  if(chi1.lt.small)then
                    small=chi1
                    do mmm=1,Nvmax
                      sss(mmm)=var(mmm)
                    enddo
                  endif
                  i16=i16+1
                endif        !if i.ne.ilo
 16           continue
              iter=iter+Nvar      !Nterms
c
              if(iter.lt.10)then
                write(*,234)iter
              endif
              if((iter.ge.10).and.(iter.lt.100))then
                write(*,235)iter
              endif
              if((iter.ge.100).and.(iter.lt.1000))then
                write(*,236)iter
              endif
              go to 1
            endif   ! ynewtry.ge.ysave
          else
            iter=iter-1
c
            if(iter.lt.10)then
              write(*,234)iter
            endif
            if((iter.ge.10).and.(iter.lt.100))then
              write(*,235)iter
            endif
            if((iter.ge.100).and.(iter.lt.1000))then
              write(*,236)iter
            endif
c
c            write(*,*)'iter = ',iter
          endif     ! ytry.le.yval(ilo)
c
          go to 2
c
 751      continue                  ! come here when delta chi is small
c
c   reset the variables at their best values and print the chi^2
c
          do mmm=1,Nvmax
            var(mmm)=sss(mmm)
          enddo
c    
          ichilabel=1
          ibest=99
          call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
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
     @      isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,
     @      ibest,ifixgamma,savesep,ichilabel,resRV1,resRV2,
     @      thresh,small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     @      icnarray,
     @      RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
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
c
          chiall=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+chisqK
     &               +chisqRV1+chisqRV2+ochi)
          chiall=chi1
          write(*,*)' '
c
              i16=i16+1
c
          do 3749 iq=1,Nvar !Nterms
            write(*,1001)var(iq),vstep(iq),svar(iq)(1:15)
            p(1,iq)=var(iq)     
            vstart(iq)=var(iq)
 3749     continue
c
          call printS(chiall)
c
          if((isw30.ge.3).and.(isw23.ge.1))then
            call writeeclipse(Ncycle,Ttimes,Tseps,isw30,Nmaxeclipse,
     @           Tdur1,Tdur2)
          endif
c         
          if(icnU.ne.430) call wlinmod(Nphase,xmod,ymodU,'modelU.linear'
     @        ,isw7)
          if(icnB.ne.430) call wlinmod(Nphase,xmod,ymodB,'modelB.linear'
     @        ,isw7)
          if(icnV.ne.430) call wlinmod(Nphase,xmod,ymodV,'modelV.linear'
     @        ,isw7)
          if(icnR.ne.430) call wlinmod(Nphase,xmod,ymodR,'modelR.linear'
     @        ,isw7)
          if(icnI.ne.430) call wlinmod(Nphase,xmod,ymodI,'modelI.linear'
     @        ,isw7)
          if(icnJ.ne.430) call wlinmod(Nphase,xmod,ymodJ,'modelJ.linear'
     @        ,isw7)
          if(icnH.ne.430) call wlinmod(Nphase,xmod,ymodH,'modelH.linear'
     @        ,isw7)
          if(icnK.ne.430) call wlinmod(Nphase,xmod,ymodK,'modelK.linear'
     @        ,isw7)
c
          if(icnU.ne.430)call wmagmod(Nphase,xmod,ymodU,'modelU.mag',
     @        isw7,zeroU)
          if(icnB.ne.430)call wmagmod(Nphase,xmod,ymodB,'modelB.mag',
     @        isw7,zeroB)
          if(icnV.ne.430)call wmagmod(Nphase,xmod,ymodV,'modelV.mag',
     @        isw7,zeroV)
          if(icnR.ne.430)call wmagmod(Nphase,xmod,ymodR,'modelR.mag',
     @        isw7,zeroR)
          if(icnI.ne.430)call wmagmod(Nphase,xmod,ymodI,'modelI.mag',
     @        isw7,zeroI)
          if(icnJ.ne.430)call wmagmod(Nphase,xmod,ymodJ,'modelJ.mag',
     @        isw7,zeroJ)
          if(icnH.ne.430)call wmagmod(Nphase,xmod,ymodH,'modelH.mag',
     @        isw7,zeroH)
          if(icnK.ne.430)call wmagmod(Nphase,xmod,ymodK,'modelK.mag',
     @        isw7,zeroK)

          call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
          call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
          call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
          if(iidint.ge.1)call wlinmod(Nphase,xmod,ymodd,'lcdisk.linear',
     @          isw7)
          call wRVmod(NRVphase,xRVmod,RV1,'star1.RV',isw7,ggamma1)
          call wRVmod(NRVphase,xRVmod,RV2,'star2.RV',isw7,ggamma2)
          if(isw30.ge.3)call wRVmod(NRVphase,xRVmod,RV3,'star3.RV',
     @          isw7,ggamma1)
          call wlinmod(NRVphase,xRVmod,dRV1,'star1.delRV',isw7)
          call wlinmod(NRVphase,xRVmod,dRV2,'star2.delRV',isw7)
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
     %            errU,'ELCresidualsU.fold')
          if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,resB,
     %            errB,'ELCresidualsB.fold')
          if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,resV,
     %            errV,'ELCresidualsV.fold')
          if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,resR,
     %            errR,'ELCresidualsR.fold')
          if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,resI,
     %            errI,'ELCresidualsI.fold')
          if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,resJ,
     %            errJ,'ELCresidualsJ.fold')
          if(icnH.ne.430)call wELCdata(Ndatamax,NdataH,xdataH,resH,
     %            errH,'ELCresidualsH.fold')
          if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,resK,
     %            errK,'ELCresidualsK.fold')
          if(icnRV1.ne.430) call wELCdata(Ndatamax,NRV1,xRV1,resRV1,
     %            errRV1,'ELCresidualsRV1.fold')
          if(icnRV2.ne.430) call wELCdata(Ndatamax,NRV2,xRV2,resRV2,
     %              errRV2,'ELCresidualsRV2.fold')
          if(icnRV3.ne.430) call wELCdata(Ndatamax,NRV3,xRV3,resRV3,
     %              errRV3,'ELCresidualsRV3.fold')
c
          icount=icount+1
          if(icount.le.Nstep(2))go to 749
c          
 69       continue   !close(2)              ! close the output file
c
c   Finally, make a file similar to ELC.inp with the current parameters.
c
          call writegridout(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,
     @       omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,
     @       betarim,rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,
     @       Nref,rLx,Period,fm,separ,gamma1,t3,g3,SA3,density,sw1,sw2,
     @       sw3,T0,idraw,iecheck,iidint,iatm,ism1,icnU,icnB,icnV,icnR,
     &       icnI,icnJ,icnH,icnK,iRVfilt,isw1,isw2,isw3,isw4,ilaw,wave,
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
          if(isw30.gt.0)call writebody3grid(Nalph3,Nbet3,tertperiod,
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

          call recordloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,
     @      Nvmax,Nvar,svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv,
     @      vstep1)
c
 1001     format(1x,2(f15.7,2x),a15)
c
          close(45)
          close(46)
          close(47)
          if(isw30.gt.0)close(48)
          if(isw30.ge.3)close(49)
          close(55)
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
          include 'dynamicssubs.for'
