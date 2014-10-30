          program gridELC
c
c   December 24, 1999
c
c   This is an optimizer code based on the 'grid search' technique
c   discussed in Bevington (1969).  This program will call the
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
          parameter(Nmaxphase=750000,Ndatamax=200000)
          parameter(Nvmax=60)
          parameter(Nmaxeclipse=1000)
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm,obv,eobv,sobv to 9.
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
     @      xRV3(Ndatamax),yRV3(Ndatamax),
     @      errRV3(Ndatamax)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,3),
     @         dwavey(8,3)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     %      sigvar(Nvmax),stepsave(Nvmax),
     $      drv1(Nmaxphase),drv2(Nmaxphase)
          dimension xRVmod(Nmaxphase)
          dimension resU(Ndatamax),resB(Ndatamax),resV(Ndatamax)
          dimension resR(Ndatamax),resI(Ndatamax),resJ(Ndatamax)
          dimension resH(Ndatamax),resK(Ndatamax),resRV1(Ndatamax)
          dimension resRV2(Ndatamax),resRV3(Ndatamax)
c
          character*40 Udatafile,svar(Nvmax),Hdatafile,Kdatafile,RV1file,
     @                 RV2file
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
c
          character*1000 parmstring
          character*2000 planetparm
          character*1700 dynparm,line
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
c   Define variables for the spots
c
          dimension spotdparm(2,4),spot1parm(2,4),spot2parm(2,4)
c
c   New variables are needed to allow for the fitting of period and T0.
c
          dimension savxdataU(Ndatamax),savydataU(Ndatamax),
     @      saverrU(Ndatamax),savydataB(Ndatamax),saverrB(Ndatamax),
     @      savydataV(Ndatamax),saverrV(Ndatamax),savydataR(Ndatamax),
     @      saverrR(Ndatamax),savydataI(Ndatamax),saverrI(Ndatamax),
     &      savydataJ(Ndatamax),saverrJ(Ndatamax),savydataH(Ndatamax),
     &      saverrH(Ndatamax),savydataK(Ndatamax),saverrK(Ndatamax),
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
          dimension fracs1(Nmaxphase,4),fracs2(Nmaxphase,4)
          dimension fracs3(Nmaxphase,4),fracs4(Nmaxphase,4)
          dimension fracs5(Nmaxphase,4),fracs6(Nmaxphase,4)
          dimension fracs7(Nmaxphase,4),fracs8(Nmaxphase,4)
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
          dimension Ncycle(34),Ttimes(34,Nmaxeclipse)
          dimension Tseps(34,Nmaxeclipse),Tdur1(34,Nmaxeclipse)
          dimension Nobscycle(34),obsTtimes(34,Nmaxeclipse)
          dimension obsTerr(34,Nmaxeclipse),Tdur2(34,Nmaxeclipse)
          dimension icnarray(34)
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
     @      atmint4,atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin

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
          call getinput(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,omega1,
     $       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     $       rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,Nref,rLx,
     @       Period,fm,separ,gamma,t3,g3,SA3,density,sw1,sw2,sw3,T0,
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
c  disable threshold for gridELC
c
          thresh=0.0d0
c
c  If the flag iatm>0, then load the model atmosphere table.
c
          if(iatm.ge.1)then
            call loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &      atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,atmint8,
     &      Tmax,Tmin,gmax,gmin)
          endif
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
c
c   Do we force the gamma to be the input value?
c
          ifixgamma=isw4
c
c   Initialize the Sobel sequence here.
c
          nnn=-11
          call sobseq(nnn,xsob)
c
          gamma1=gamma
c
c   Disable the draw option
c
          idraw=0
c
c   Fix the inner disk radius at fill2.
c
          if(teff2.gt.0.0d0)rinner=fill2
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
c   Get the parameters for the grid search.
c
          call getloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,
     @      Nvmax,Nvar,svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)

c   UPDATE FEBRUARY 5, 2005
c
c   Add also the ability to include the luminosity ratio in the total chi^2
c
          ifrac=-99
          ilum=-99
          do 111 i=1,Nobv
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
 111      continue
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
     @        NRV3,xRV3,yRV3,errRV3,Ndatamax,icnRV3,Nmaxeclipse)
          endif
c
          isvel1=0
          isvel2=0
          if(icnRV1.ne.430)isvel1=299
          if(icnRV2.ne.430)isvel2=299
c
          gamma1=gamma
          gamma2=gamma
c
c   Define the variables whose values will be printed on the screen at the
c   end of the iteration.
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
c   Start the grid search here.  This code is more or less from Bevington
c   (1969).
c
          chi1=0.0
          chi2=0.0
          chi3=0.0
          chilimb=0.0d0
          small=1.0d33
          Ndattot=NdataU+NdataB+NdataV+NdataR+NdataI+NdataJ+NdataH+
     @       NdataK+NRV1+NRV2+Nobv+NRV3
c
          Nterms=0
c
c   UPDATE May 30, 2002
c
c   Change the upper loop limit to max(1,Nvar)
c
          do 699 i7=1,max(1,Nvar)
            var(i7)=vstart(i7)     !initialize the variables
            kkk=icnvrt(svar(i7)(1:2))
            stepsave(i7)=vstep(i7)
            if(kkk.ne.430)then
              Nterms=Nterms+1
            endif
 699      continue
c
c   UPDATE May 30, 2002
c
c   Change the if-then statement to Nvar.eq.0
c
 749      if(Nvar.eq.0)go to 69   !no valid variables specified
c
          write(*,*)'Nterms, Nvar, Ndattot',Nterms,Nvar,Ndattot
          i16=0
          do 750 ijk=1,Nstep(1)
            do 700 i7=1,Nvar     !Nterms    
              kkk=icnvrt(svar(i7)(1:2))
              if(kkk.eq.430)go to 700
c
              ibest=0
              ichilabel=1
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
     @           saverrRV1,savxRV2,savyRV2,saverrRV2,ifrac,ilum,i16,
     @           isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,isavNH,
     @           isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,ibest,
     @           ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
     @           small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     &           icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,
     @           NRV3,
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
              call printgriditer(i16+1,'model number = ',ijk,
     @          'iteration number = ',svar(i7),isw29)
c
              if(Ndattot-Nterms.ge.1)then
                chi1=chi1/dabs(dble(Ndattot-Nterms))
              endif

              if(chi1.lt.small)then
                small=chi1
                do mmm=1,Nvmax
                  sss(mmm)=var(mmm)
                enddo
              endif

              i16=i16+1
              fn=0.0
              delta=vstep(i7) 
              var(i7)=var(i7)+delta
c
              ibest=0
              ichilabel=2
              call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     @           ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @           ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @           xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     &           fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi2,
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
     @           saverrRV1,savxRV2,savyRV2,saverrRV2,ifrac,ilum,i16,
     @           isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,isavNH,
     @           isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,ibest,
     @           ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
     @           small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     &           icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,
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
              if(Ndattot-Nterms.ge.1)then
                chi2=chi2/dabs(dble(Ndattot-Nterms))
              endif

              call printgriditer(i16+1,'model number = ',ijk,
     @          'iteration number = ',svar(i7),isw29)
c
              if(chi2.lt.small)then
                small=chi2
                do mmm=1,Nvmax
                  sss(mmm)=var(mmm)
                enddo
              endif
c
              i16=i16+1
c
              diff=chi1-chi2
              if(diff.gt.0.0)go to 5061
c
              delta=-delta
              var(i7)=var(i7)+delta
c  
              csave=chi1
              chi1=chi2
              chi2=csave
c 
c   incriment the variable var(i7) until chi square increases
c
 5061         fn=fn+1.0
c 
              var(i7)=var(i7)+delta
c 
              ibest=0
              ichilabel=3
              call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     @           ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @           ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @           xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     &           fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi3,
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
     @           saverrRV1,savxRV2,savyRV2,saverrRV2,ifrac,ilum,i16,
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
              if(Ndattot-Nterms.ge.1)then
                chi3=chi3/dabs(dble(Ndattot-Nterms))
              endif
              call printgriditer(i16+1,'model number = ',ijk,
     @        'iteration number = ',svar(i7),isw29)
c
              if(chi3.lt.small)then
                small=chi3
                do mmm=1,Nvmax
                  sss(mmm)=var(mmm)
                enddo
              endif
c
              i16=i16+1
c
              diff23=chi3-chi2
              if(diff23.lt.0)then
                chi1=chi2
                chi2=chi3
                go to 5061
              endif          
c
c  Find the minimum of parabola defined by the last three points.
c
              if(chi3-chi2.eq.0.0)go to 700
              delta=delta*(1.0/(1.0+(chi1-chi2)/(chi3-chi2))+0.5)
              var(i7)=var(i7)-delta
              free=dabs(dble(Ndattot-Nterms))
              sigvar(i7)=vstep(i7)*sqrt(2.0/(free*(chi3-2.0*chi2+chi1)))
              vstep(i7)=vstep(i7)*fn/3.0
c
              if((Nterms.eq.1).or.(i7.eq.Nvar))then   ! (i7.eq.Nterms))then
c 
                ibest=0
                ichilabel=4
                call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     @           ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @           ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @           xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     &           fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     &           chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi4,
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
     @           saverrRV1,savxRV2,savyRV2,saverrRV2,ifrac,ilum,i16,
     @           isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,isavNH,
     @           isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,ibest,
     @           ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
     @           small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     &           icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,
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
                if(Ndattot-Nterms.ge.1)then
                  chi4=chi4/dabs(dble(Ndattot-Nterms))
                endif

                call printgriditer(i16+1,'model number = ',ijk,
     @           'iteration number = ',svar(i7),isw29)
c
                if(chi4.lt.small)then
                  small=chi4
                  do mmm=1,Nvmax
                    sss(mmm)=var(mmm)
                  enddo
                endif
c
                i16=i16+1
              endif

 700        continue  !loop over the 7 terms
c
            write(*,*)' '
c
            do 7749 iq=1,Nvar      !Nterms
              write(*,1001)var(iq),sigvar(iq),svar(iq)(1:15)
 7749       continue
c
            call printS(small)  

            if((isw30.ge.3).and.(isw23.ge.1))then
              call writeeclipse(Ncycle,Ttimes,Tseps,isw30,Nmaxeclipse,
     @           Tdur1,Tdur2)
            endif
c         
            if(icnU.ne.430)call wlinmod(Nphase,xmod,ymodU,
     @          'modelU.linear',isw7)
            if(icnB.ne.430)call wlinmod(Nphase,xmod,ymodB,
     @          'modelB.linear',isw7)
            if(icnV.ne.430)call wlinmod(Nphase,xmod,ymodV,
     @          'modelV.linear',isw7)
            if(icnR.ne.430)call wlinmod(Nphase,xmod,ymodR,
     @          'modelR.linear',isw7)
            if(icnI.ne.430)call wlinmod(Nphase,xmod,ymodI,
     @          'modelI.linear',isw7)
            if(icnJ.ne.430)call wlinmod(Nphase,xmod,ymodJ,
     @          'modelJ.linear',isw7)
            if(icnH.ne.430)call wlinmod(Nphase,xmod,ymodH,
     @          'modelH.linear',isw7)
            if(icnK.ne.430)call wlinmod(Nphase,xmod,ymodK,
     @          'modelK.linear',isw7)
            call wRVmod(NRVphase,xRVmod,RV1,'star1.RV',isw7,ggamma1)
            call wRVmod(NRVphase,xRVmod,RV2,'star2.RV',isw7,ggamma2)
            call wRVmod(NRVphase,xRVmod,RV3,'star3.RV',isw7,ggamma1)
c
            if(icnU.ne.430)call wmagmod(Nphase,xmod,ymodU,'modelU.mag',
     @         isw7,zeroU)
            if(icnB.ne.430)call wmagmod(Nphase,xmod,ymodB,'modelB.mag',
     @         isw7,zeroB)
            if(icnV.ne.430)call wmagmod(Nphase,xmod,ymodV,'modelV.mag',
     @         isw7,zeroV)
            if(icnR.ne.430)call wmagmod(Nphase,xmod,ymodR,'modelR.mag',
     @         isw7,zeroR)
            if(icnI.ne.430)call wmagmod(Nphase,xmod,ymodI,'modelI.mag',
     @         isw7,zeroI)
            if(icnJ.ne.430)call wmagmod(Nphase,xmod,ymodJ,'modelJ.mag',
     @         isw7,zeroJ)
            if(icnH.ne.430)call wmagmod(Nphase,xmod,ymodH,'modelH.mag',
     @         isw7,zeroH)
            if(icnK.ne.430)call wmagmod(Nphase,xmod,ymodK,'modelK.mag',
     @         isw7,zeroK)
c
            if(itime.gt.0)then
              if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,
     @           ydataU,errU,'ELCdataU.fold')
              if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,
     @           ydataB,errB,'ELCdataB.fold')
              if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,
     @           ydataV,errV,'ELCdataV.fold')
              if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,
     @           ydataR,errR,'ELCdataR.fold')
              if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,
     @           ydataI,errI,'ELCdataI.fold')
              if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,
     @           ydataJ,errJ,'ELCdataJ.fold')
              if(icnH.ne.430) call wELCdata(Ndatamax,NdataH,xdataH,
     @           ydataH,errH,'ELCdataH.fold')
              if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,
     @           ydataK,errK,'ELCdataK.fold')
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
c   Reset the variables at their best values.
c
          do mmm=1,Nvmax
            var(mmm)=sss(mmm)
          enddo
c
 750      continue  ! loop over grid searches
c
          continue                  ! come here when delta chi is small
c
c   reset the variables at their best values and print the chi^2
c
          do mmm=1,Nvmax
            var(mmm)=sss(mmm)
          enddo
c
 69       ibest=0
          ichilabel=1
          call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     @       ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @       ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @       xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     &       fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     &       chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chiall,
     @       NdataU,xdataU,ydataU,errU,zeroU,resU,NdataB,xdataB,
     @       ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,errV,zeroV,
     @       resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
     @       xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,
     @       errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,
     @       NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,
     @       errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,sobv,
     @       obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,saveasini,
     @       savxdataU,savydataU,saverrU,savxdataB,savydataB,
     @       saverrB,savxdataV,savydataV,saverrV,savxdataR,
     @       savydataR,saverrR,savxdataI,savydataI,saverrI,
     @       savxdataJ,savydataJ,saverrJ,savxdataH,savydataH,
     @       saverrH,savxdataK,savydataK,saverrK,savxRV1,savyRV1,
     @       saverrRV1,savxRV2,savyRV2,saverrRV2,ifrac,ilum,i16,
     @       isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,isavNH,
     @       isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,ibest,
     @       ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
     @       small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     &       icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
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
          call printgriditer(i16+1,'model number = ',ijk,
     @        'iteration number = ',svar(i7),isw29)
c
          write(*,*)' '
          do 3749 iq=1,Nvar      !Nterms
            write(*,1001)var(iq),sigvar(iq),svar(iq)(1:15)
 3749     continue
c
          call printS(chiall)
c
          i16=i16+1

          if((isw30.ge.3).and.(isw23.ge.1))then
            call writeeclipse(Ncycle,Ttimes,Tseps,isw30,Nmaxeclipse,
     @           Tdur1,Tdur2)
          endif
c         
          if(icnU.ne.430)call wlinmod(Nphase,xmod,ymodU,'modelU.linear',
     @       isw7)
          if(icnB.ne.430)call wlinmod(Nphase,xmod,ymodB,'modelB.linear',
     @       isw7)
          if(icnV.ne.430)call wlinmod(Nphase,xmod,ymodV,'modelV.linear',
     @       isw7)
          if(icnR.ne.430)call wlinmod(Nphase,xmod,ymodR,'modelR.linear',
     @       isw7)
          if(icnI.ne.430)call wlinmod(Nphase,xmod,ymodI,'modelI.linear',
     @       isw7)
          if(icnJ.ne.430)call wlinmod(Nphase,xmod,ymodJ,'modelJ.linear',
     @       isw7)
          if(icnH.ne.430)call wlinmod(Nphase,xmod,ymodH,'modelH.linear',
     @       isw7)
          if(icnK.ne.430)call wlinmod(Nphase,xmod,ymodK,'modelK.linear',
     @       isw7)
c
          if(icnU.ne.430)call wmagmod(Nphase,xmod,ymodU,'modelU.mag',
     @       isw7,zeroU)
          if(icnB.ne.430)call wmagmod(Nphase,xmod,ymodB,'modelB.mag',
     @       isw7,zeroB)
          if(icnV.ne.430)call wmagmod(Nphase,xmod,ymodV,'modelV.mag',
     @       isw7,zeroV)
          if(icnR.ne.430)call wmagmod(Nphase,xmod,ymodR,'modelR.mag',
     @       isw7,zeroR)
          if(icnI.ne.430)call wmagmod(Nphase,xmod,ymodI,'modelI.mag',
     @       isw7,zeroI)
          if(icnJ.ne.430)call wmagmod(Nphase,xmod,ymodJ,'modelJ.mag',
     @       isw7,zeroJ)
          if(icnH.ne.430)call wmagmod(Nphase,xmod,ymodH,'modelH.mag',
     @       isw7,zeroH)
          if(icnK.ne.430)call wmagmod(Nphase,xmod,ymodK,'modelK.mag',
     @       isw7,zeroK)
c
          call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
          call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
          call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
          if(iidint.ge.1)call wlinmod(Nphase,xmod,ymodd,'lcdisk.linear',
     @       isw7)
c
          call wRVmod(NRVphase,xRVmod,RV1,'star1.RV',isw7,ggamma1)
          call wRVmod(NRVphase,xRVmod,RV2,'star2.RV',isw7,ggamma2)
          if(isw30.ge.3)call wRVmod(NRVphase,xRVmod,RV3,'star3.RV',
     @          isw7,ggamma1)
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
c   Finally, make a file similar to ELC.inp with the current parameters.
c
          call writegridout(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,
     @      omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,
     %      betarim,rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,
     &      Nref,rLx,Period,fm,separ,gamma1,t3,g3,SA3,density,sw1,sw2,
     @      sw3,T0,idraw,iecheck,iidint,iatm,ism1,icnU,icnB,icnV,icnR,
     &      icnI,icnJ,icnH,icnK,iRVfilt,isw1,isw2,isw3,isw4,ilaw,wave,
     @      dbolx,dboly,dwavex,dwavey,ecc,argper,pshift,sw5,sw6,sw7,
     @      sw8,sw9,ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,
     @      spot2parm,spotdparm,primmass,primK,primrad,ratrad,frac1,
     @      frac2,ecosw,temprat,idark1,idark2,isw12,isw13,isw21,isw22,
     @      isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,
     @      sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,isw25,isw26,
     @      isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,ocose,osine,
     @      omegadot,contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49)
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
c
          call recordloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,
     $      Nvmax,Nvar,svar,var,vstart,stepsave,Nstep,Nobv,sobv,obv,
     @      eobv,vstep)
c
 1001     format(1x,2(f15.7,2x),a15)

          close(45)
          close(46)
          if(isw24.gt.0)close(47)
          if(isw30.gt.0)close(48)
          if(isw30.ge.3)close(49)

          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
          include 'dynamicssubs.for'
