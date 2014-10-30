          program markovELC
c
c   November 21, 2008
c
c   This program will read the file gridloop.opt and optimize
c   the model using a Monte Carlo Markov Chain.
c
c   The meaning of the triplets of numbers is as follows
c   
c   central value     variance for jump      N_sigma_limit
c
c   The parameter is bounded between the range 
c   (mean)-(N_sigma_limit*(variance))  to (mean)+(N_sigma_limit*(variance))
c
          implicit double precision (a-h,o-z)
c
c          parameter(Nmaxphase=670000,Ndatamax=100000)
c          parameter(Nvmax=60)


          parameter(Nmaxphase=770000,Ndatamax=200000) 
c          parameter(Nmaxphase=202000,Ndatamax=44000)
          parameter(Nvmax=60,Nchainmax=60000)
          parameter(Nmaxeclipse=1000)
c
          dimension jumpcount(Nvmax),jumptotal(Nvmax)
          dimension obsparm(19),obv(19),eobv(19)
          character*40 sobv(19),command
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
     @      dwavey(8,3)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     %      savestep(Nvmax),Nstepnew(Nvmax)
          character*40 Udatafile,svar(nvmax),Hdatafile,Kdatafile,
     &             RV1file,RV2file
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          dimension parmarray(Nvmax,Nchainmax),chiarr(Nchainmax)
          dimension iparm(Nvmax),chiparm(Nvmax,Nchainmax),
     @        saveparm(Nvmax,Nchainmax)
          dimension yscr(Nchainmax),zscr(Nchainmax)
          dimension xRVmod(Nmaxphase)
c
          character*7 extension
          character*2 chainstring
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
c   This array is needed for the sub-random Sobel sequence, used in
c   place of ran9
c
          dimension xsob(2)
c
          dimension powercoeff(8,9)
c
          character*1000 parmstring
          character*2000 planetparm
          character*1700 dynparm,line
          character*1000 pstringparm(100000)
          character*30 outstring,outstring1
          dimension pstringchi(100000)
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
          dimension resU(Ndatamax),resB(Ndatamax),resV(Ndatamax)
          dimension resR(Ndatamax),resI(Ndatamax),resJ(Ndatamax)
          dimension resH(Ndatamax),resK(Ndatamax),resRV1(Ndatamax)
          dimension resRV2(Ndatamax),resRV3(Ndatamax)
c
          dimension compfracs(8,3),ochidisk(8)
c
          dimension fracs1(Nmaxphase,4),fracs2(Nmaxphase,4)
          dimension fracs3(Nmaxphase,4),fracs4(Nmaxphase,4)
          dimension fracs5(Nmaxphase,4),fracs6(Nmaxphase,4)
          dimension fracs7(Nmaxphase,4),fracs8(Nmaxphase,4)
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
c            
c          common /stringblock/ parmstring,planetparm,dynparm,line
c
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
     &     atmint4,atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin
c
          common /intatm/  Nlines,Nmu,Nalph3,Nbet3,itconj,it1,it2,it3,
     @      it4
c
          common /medblock/ rmed 
c
          common /ranblock/ idum
c
c   initialize variables
c
          ilimbcheck=0
          ipstringcount=1
c
c   Open the parameter file and read all of the parameters. 
c   Pass the parameters to the light curve routines
c   via a common block.  
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
     @       bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,sw29,
     @       sw30,contam,Tconj,beam1,beam2,isw25,isw26,isw27,isw28,
     @       isw29,isw30,isw31,isw32,isw33,isw34,ocose,osine,omegadot,
     @       contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49)
c
           savesep=separ
c
c   Set the threshold to record models in files
c
          thresh=sw48
c
c  If the flag iatm>0, then load the model atmosphere table.
c
          gmin=0.0d0
          gmax=0.0d0
          Tmin=0.0d0
          Tmax=0.0d0
          if(iatm.ge.1)then
            call loadtable(maxlines,maxmu,Nlines,atmT,atmg,atmmu,Nmu,
     &         atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,
     &         atmint8,Tmax,Tmin,gmax,gmin)
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
c   Add a new variable to enable median fitting.
c   if sw6 > 1.0, then use median fitting
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
          gamma1=gamma
          gamma2=gamma
c
c   Fix the inner disk radius at fill2.
c
          if(teff2.gt.0.0)rinner=fill2
c
c   UPDATE August 13, 2001
c
c   If the flag isw9 (=ielite) is more than 0, insert the parameters
c   set in the ELC.inp file into the population.
c
          mmm=-11
          call sobseq(mmm,xsob)
c
c   open up optional file markovELC.inp and read
c   lengthchain
c   Nchain
c   ireport
c
          lengthchain=1000
          Nchain=50
          ireport=20         
          ichain=0
c
          open(unit=44,file='markovELC.inp',status='old',err=2)
c
          read(44,*,end=6)lengthchain
          read(44,*,end=6)Nchain
          read(44,*,end=6)ireport
c
6         close(44)
          go to 3
 
2         write(*,*)'markovELC.inp not found.  Using default values'
          write(*,*)'of lengthchain=1000, Nchain=50, and ireport=20'
c
c   Initialize the steps.
c
3         do 1 i=1,nvmax
            Nstep(i)=1
            var(i)=0.0d0
            svar(i)='none'
            jumpcount(i)=0
            jumptotal(i)=0
            sss(i)=0.0d0
 1        continue
c
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
c   Get the parameters for the grid search.
c
          call getloopopt(Udatafile,Bdatafile,Vdatafile,Rdatafile,
     $      Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,RV2file,
     @      Nvmax,Nvar,svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,eobv)


c
c   UPDATE December 3, 2009
c
c   Check to see that the lower bound is less than the upper bound
c
           do 24 kkk=1,Nvar
             if(vstart(kkk).ge.vstep(kkk))then
               write(*,4321)kkk,svar(kkk)
               write(*,4322)vstart(kkk),vstep(kkk)
               stop
              endif
24          continue
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
          if(icnU.ne.430)call kopydata(Ndatamax,NdataU,xdataU,ydataU,
     @       errU,isavNU,savxdataU,savydataU,saverrU)
          if(icnB.ne.430)call kopydata(Ndatamax,NdataB,xdataB,ydataB,
     @       errB,isavNB,savxdataB,savydataB,saverrB)
          if(icnV.ne.430)call kopydata(Ndatamax,NdataV,xdataV,ydataV,
     @       errV,isavNV,savxdataV,savydataV,saverrV)
          if(icnI.ne.430)call kopydata(Ndatamax,NdataI,xdataI,ydataI,
     @       errI,isavNI,savxdataI,savydataI,saverrI)
          if(icnJ.ne.430)call kopydata(Ndatamax,NdataJ,xdataJ,ydataJ,
     @       errJ,isavNJ,savxdataJ,savydataJ,saverrJ)
          if(icnH.ne.430)call kopydata(Ndatamax,NdataH,xdataH,ydataH,
     @       errH,isavNH,savxdataH,savydataH,saverrH)
          if(icnK.ne.430)call kopydata(Ndatamax,NdataK,xdataK,ydataK,
     @       errK,isavNK,savxdataK,savydataK,saverrK)
          if(icnR.ne.430)call kopydata(Ndatamax,NdataR,xdataR,ydataR,
     @       errR,isavNR,savxdataR,savydataR,saverrR)
          if(icnRV1.ne.430)call kopydata(Ndatamax,NRV1,xRV1,yRV1,errRV1,
     @       isavRV1,savxRV1,savyRV1,saverrRV1)
          if(icnRV2.ne.430)call kopydata(Ndatamax,NRV2,xRV2,yRV2,errRV2,
     @       isavRV2,savxRV2,savyRV2,saverrRV2)
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
     @       NRV3,xRV3,yRV3,errRV3,Ndatamax,icnRV3,Nmaxeclipse)
          endif
c
          isvel1=0
          isvel2=0
          if(icnRV1.ne.430)isvel1=299
          if(icnRV2.ne.430)isvel2=299
c
          if(isw32.lt.0)then
            idum=isw32
          else
            idum=-1234567
          endif
c
c   if sw25 > 0, use that as an error for asini.
c
          saveasini=sw5
c
c   Start the looping here.  
c
          do 9999 ichain=1,Nchain
            imod=0
            chi1=0.0d0
            small=9.0d33
            ismall=0
            Ndattot=0
            Ndattot=NdataU+NdataB+NdataV+NdataR+NdataI+NdataJ+NdataH+
     @        NdataK+NRV1+NRV2+NRV3
            Nterms=0
            do 699 i7=1,Nvar
              var(i7)=ran1(idum)*(vstep(i7)-vstart(i7))+vstart(i7)
              if(ichain.eq.1)then
                savestep(i7)=dabs(vstart(i7)-vstep(i7))/
     @             (2.0d0*dble(Nstep(i7)))
              endif
              kkk=icnvrt(svar(i7)(1:2))
              if(kkk.eq.430)then
                go to 749
              else
                Nterms=Nterms+1
              endif
 699        continue
c
            write(*,*)' '
            write(*,*)'Nterms, Nvar, Ndattot',Nterms,Nvar,Ndattot
          
 749        continue
c
c   Open an output file for the fitting statistics.
c
            write(extension,3333)ichain+1000000-1
c
            open(unit=45,file='generation.'//extension,status='unknown')
            if(isw24.ge.1)open(unit=47,file='ELCratio.'//extension,
     @           status='unknown')
            open(unit=55,file='chi.'//extension,status='unknown')
            open(unit=46,file='ELCparm.'//extension,status='unknown')
            if(isw30.ge.1)open(unit=48,file='ELCbody3parm.'//extension,
     @        status='unknown')
            if(isw30.ge.3)open(unit=49,file='ELCdynparm.'//extension,
     @        status='unknown')
c   
            write(*,*)' '
            do 99 j=1,Nterms
              write(*,96)j,svar(j),savestep(j)
              parmarray(j,1)=var(j)
              iparm(j)=1
              ipstringcount=1
 99         continue
c
            i16=1                          !counter
            if(isw22.ge.1)ifastflag=1
            if(i16.eq.1)ifastflag=0
            ifastflag=0
c
            if((isw9.ge.1))then
              call varassign(Nvmax,svar,var,fill1,fill2,omega1,omega2,Q,
     @          finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,rLx,
     @          separ,gamma,t3,g3,sa3,ecc,argper,pshift,spot1parm,
     @          spot2parm,spotdparm,period,T0,alb1,alb2,dwavex,dwavey,
     @          primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @          bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2,
     @          contam,ocose,osine,isw29,tertperiod,tertt0,tertecos,
     @          tertesin,tertincl,tertOmega,tertQ,Tgrav1,Tgrav2,
     @          tertconj,omegadot,contamS0,contamS1,contamS2,contamS3,
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
     @       P8ratrad,sw72,sw73)
c
c   if ilaw=4, then make sure the sum of the x and y coefficients
c   does not exceed 1
c
              if(ilaw.eq.4)then
                do jjj=1,8
                  rlimbsum=dwavex(jjj,1)+dwavey(jjj,1)
                  if(rlimbsum.gt.1.0d0)then
                    write(*,*)'the sum of x and y exceeds 1.0 ',
     $               'for index ',jjj
                    stop
                  endif
                  rlimbsum=dwavex(jjj,2)+dwavey(jjj,2)
                  if(rlimbsum.gt.1.0d0)then
                    write(*,*)'the sum of x and y exceeds ',
     $               '1.0 for index ',jjj
                    stop
                  endif
                enddo
              endif
              ilimbcheck=0
              do 799 j=1,Nterms
                kk=icnvrt(svar(j)(1:2))
                if(kk.ge.1752.and.kk.le.1759)ilimbcheck=99
                if(kk.ge.1784.and.kk.le.1791)ilimbcheck=99
                if(kk.ge.1816.and.kk.le.1823)ilimbcheck=99
                if(kk.ge.1720.and.kk.le.1727)ilimbcheck=99
                varx=var(j)
                parmarray(j,1)=varx
                if(varx.lt.vstart(j))then
                  write(*,4323)svar(j)
                  write(*,4324)varx,vstart(j),vstep(j)
                  stop
                endif
                if(varx.gt.vstep(j))then
                  write(*,4323)svar(j)
                  write(*,4324)varx,vstart(j),vstep(j)
                  stop
                endif
 799          continue
c
            endif
c
            chilimb=0.0d0
            if(ilimbcheck.gt.1)then
              if(ilaw.eq.4)then
                call getchilimb(Nvmax,Nterms,svar,dwavex,dwavey,chilimb)
              endif
              if(ilaw.eq.14)then
                call getchilimb(Nvmax,Nterms,svar,dwavex,dwavey,chilimb)
              endif
            endif
c
            ifastflag=0
            ibest=0
            ichilabel=1
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
            imod=imod+1
            chiarr(i16)=chi1
            nnn=0
c
            call printiter(imod,'model number = ',ichain,
     @          'chain number = ')
c
            if(chi1.lt.small)then
              ismall=ismall+1
              small=chi1
              do mmm=1,Nvmax
                sss(mmm)=var(mmm)
              enddo
            endif
c
c   Main loop
c
            do 10000 ig=2,lengthchain
  
              do 1121 j=1,Nterms
  
                if((jumpcount(j).gt.50).or.(jumptotal(j).gt.80))then
                  fracjump=dble(jumpcount(j))/dble(jumptotal(j))
                  if(fracjump.lt.0.4d0)then
                    stepnew=savestep(j)*(fracjump-0.4d0+1.0d0)
                    if(stepnew.gt.2.0d0*dabs(vstart(j)-vstep(j)))then
                      stepnew=0.5d0*dabs(vstart(j)-vstep(j))
                    endif
                    write(*,1120)fracjump,j,savestep(j),stepnew
                    savestep(j)=stepnew
                    jumpcount(j)=0
                    jumptotal(j)=0
                    go to 1121
                  endif    
                  if(fracjump.ge.0.4d0)then
                    stepnew=savestep(j)*(fracjump-0.4d0+1.0d0)
                    if(stepnew.gt.2.0d0*dabs(vstart(j)-vstep(j)))then
                      stepnew=0.5d0*dabs(vstart(j)-vstep(j))
                    endif
                    write(*,1120)fracjump,j,savestep(j),stepnew
                    savestep(j)=stepnew
                    jumpcount(j)=0
                    jumptotal(j)=0
                    go to 1121
                  endif
                endif
c
1121          continue
              i16=ig
c
c   pick a random parameter and offset it by a random gaussian deviate
c
1234          j=1+idnint(dble(Nterms)*ran1(idum))
              if(j.lt.1)j=1
              if(j.gt.Nterms)j=Nterms
              jumpindex=j
              vx1=parmarray(j,ig-1)
              vx2=savestep(j)
 997          rjump=vx2*gasdev(idum)
              varx=vx1+rjump
              if(varx.lt.vstart(j))go to 997
              if(varx.gt.vstep(j))go to 997
              var(j)=varx
c
              write(*,*)' '
              if(rmed.ge.1.0d0)then
                call printsmallmed(small)
              else
                call printsmall(small)
              endif
c
              call pstring(' jump',8,rjump,outstring,lll)
              if(j.lt.10)then
                write(*,996,advance='no')j,svar(j),outstring(1:lll)
              endif
              if((j.ge.10).and.(j.lt.100))then
                write(*,1996,advance='no')j,svar(j),outstring(1:lll)
              endif
c
c   This call to assignvar and the subsequent check is needed to
c   make sure the sum of the limb darkening coefficients don't 
c   exceed 1.  A chi^2 penalty will be added inside the monster
c   subroutine.  However, it seems that the routine can get
c   stuck, even with the chi^2 penalty.  This call to assignvar
c   will make it impossible to violate the coniditon inside the monster
c   routine.
c
              call assignvar(Nvmax,svar,var,fill1,fill2,omega1,omega2,Q,
     @          finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,rLx,
     @          separ,gamma,t3,g3,sa3,ecc,argper,pshift,spot1parm,
     @          spot2parm,spotdparm,period,T0,alb1,alb2,dwavex,dwavey,
     @          primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @          bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2,
     @          contam,ocose,osine,isw29,tertperiod,tertt0,tertecos,
     @          tertesin,tertincl,tertOmega,tertQ,Tgrav1,Tgrav2,
     @          tertconj,omegadot,contamS0,contamS1,contamS2,contamS3,
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
     @       P8ratrad,sw72,sw73)
c
            if((ilaw.eq.4).or.(ilaw.eq.14))then
              if(ilimbcheck.gt.1)then
                kk=icnvrt(svar(j)(1:2))
                if(kk.ge.1752.and.kk.le.1759)then
                  do jjj=1,8
c                    if(j.eq.jjj)then
                      rlimbsum=dwavex(jjj,1)+dwavey(jjj,1)
                      if(rlimbsum.gt.1.0d0)go to 997
                      rlimbsum=dwavex(jjj,2)+dwavey(jjj,2)
                      if(rlimbsum.gt.1.0d0)go to 997
c                    endif
                  enddo
                endif
                if(kk.ge.1784.and.kk.le.1791)then
                  do jjj=1,8
c                    if(j.eq.jjj)then
                      rlimbsum=dwavex(jjj,1)+dwavey(jjj,1)
                      if(rlimbsum.gt.1.0d0)go to 997
                      rlimbsum=dwavex(jjj,2)+dwavey(jjj,2)
                      if(rlimbsum.gt.1.0d0)go to 997
c                    endif
                  enddo
                endif
                if(kk.ge.1816.and.kk.le.1823)then
                  do jjj=1,8
c                    if(j.eq.jjj)then
                      rlimbsum=dwavex(jjj,1)+dwavey(jjj,1)
                      if(rlimbsum.gt.1.0d0)go to 997
                      rlimbsum=dwavex(jjj,2)+dwavey(jjj,2)
                      if(rlimbsum.gt.1.0d0)go to 997
c                    endif
                  enddo
                endif
                if(kk.ge.1720.and.kk.le.1727)then
                  do jjj=1,8
c                    if(j.eq.jjj)then
                      rlimbsum=dwavex(jjj,1)+dwavey(jjj,1)
                      if(rlimbsum.gt.1.0d0)go to 997
                      rlimbsum=dwavex(jjj,2)+dwavey(jjj,2)
                      if(rlimbsum.gt.1.0d0)go to 997
c                    endif
                  enddo
                endif
              endif
             endif

              if(i16.eq.1)ifastflag=0
c
              chilimb=0.0d0
              if(ilimbcheck.gt.1)then
                if(ilaw.eq.4)then
                  call getchilimb(Nvmax,Nterms,svar,dwavex,dwavey,chilimb)
                endif
                if(ilaw.eq.14)then
                  call getchilimb(Nvmax,Nterms,svar,dwavex,dwavey,chilimb)
                endif
              endif
c
              ifastflag=0
              ibest=0
              ichilabel=1
              call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     @           ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @           ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @           xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     @           fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     @           chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi1,
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
              imod=imod+1
              chiarr(i16)=chi1
              call printiter(imod,'model number = ',ichain,
     @          'chain number = ')
c
              if(chi1.lt.small)then
                ismall=ismall+1
                small=chi1
                do mmm=1,Nvmax
                  sss(mmm)=var(mmm)
                enddo
              endif
c 
              jumptotal(jumpindex)=jumptotal(jumpindex)+1
              if(chiarr(i16).lt.chiarr(i16-1))then
c
c   The chi^2 is lower.  Put the parameter into the chain and
c   continue the loop.
c
                do 876 j=1,Nterms
                  parmarray(j,ig)=var(j)
876             continue
                chiparm(jumpindex,iparm(jumpindex))=chiarr(i16)
                saveparm(jumpindex,iparm(jumpindex))=var(jumpindex)
                iparm(jumpindex)=iparm(jumpindex)+1
                jumpcount(jumpindex)=jumpcount(jumpindex)+1
                pstringparm(ipstringcount)=parmstring
                pstringchi(ipstringcount)=chiarr(i16)
                ipstringcount=ipstringcount+1
  
                call pstring('Delta_chi',5,chiarr(i16)-chiarr(i16-1),
     $              outstring,lll)
                if(jumpindex.lt.10)then
                  write(*,95)outstring(1:lll),jumpindex
                else
                  write(*,1195)outstring(1:lll),jumpindex
                endif
c
                mmm=jumpcount(jumpindex)
                nnn=jumptotal(jumpindex)
                call writejump(mmm,nnn,jumpindex)
                go to 8888
              else
                deltachi=dabs(chiarr(i16)-chiarr(i16-1))
                prob1=1.0d0
c  
                call chistring('chi_start',chiarr(i16-1),outstring,lll1)
                call chistring(' chi_end',chiarr(i16),outstring1,lll2)
c  
                write(*,2991)outstring(1:lll1),outstring1(1:lll2)
c
                if(prob1.lt.1.0d-100)then
                  prob=dexp(-deltachi*0.5d0)
                else
                  prob=dexp(-deltachi*0.5d0)
                endif
                uran=ran1(idum)
                if(uran.lt.prob)then
                  do 877 j=1,Nterms
                    parmarray(j,ig)=var(j)
877               continue
                  chiparm(jumpindex,iparm(jumpindex))=chiarr(i16)
                  saveparm(jumpindex,iparm(jumpindex))=var(jumpindex)
                  iparm(jumpindex)=iparm(jumpindex)+1
                  jumpcount(jumpindex)=jumpcount(jumpindex)+1
                  pstringparm(ipstringcount)=parmstring
                  pstringchi(ipstringcount)=chiarr(i16)
                  ipstringcount=ipstringcount+1
  
                  call pstring('Delta_chi',5,chiarr(i16)-chiarr(i16-1),
     @              outstring,lll)
                  write(*,93)outstring(1:lll),uran,prob,jumpindex
                  mmm=jumpcount(jumpindex)
                  nnn=jumptotal(jumpindex)
                  call writejump(mmm,nnn,jumpindex)
c
                  go to 8888
                endif
                var(j)=parmarray(j,i16-1)
                call pstring('Delta_chi',5,chiarr(i16)-chiarr(i16-1),
     @             outstring,lll)
                if(jumpindex.lt.10)then
                  write(*,92)outstring(1:lll),uran,prob,jumpindex
                endif
                if((jumpindex.ge.10).and.(jumpindex.lt.100))then
                  write(*,292)outstring(1:lll),uran,prob,jumpindex
                endif
                go to 1234 
              endif
c
c   reset the variables at their best values and print the chi^2
c
8888          if((ismall.gt.ireport).or.(ig.eq.lengthchain))then
                ismall=0
                do mmm=1,Nvmax
                  var(mmm)=sss(mmm)
                enddo
c
                chilimb=0.0d0
                if(ilimbcheck.gt.1)then
                  if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @              dwavey,chilimb)
                  if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @              dwavey,chilimb)
                endif
c
                ifastflag=0
                ibest=0
                ichilabel=1
                call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     @           ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @           ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @           xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     @           fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     @           chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi1,
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
     %           icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,
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
                chiarr(i16)=chi1
                chiall=chi1
                write(*,*)' '
                call printS(chiall)
c  
                if((isw30.ge.3).and.(isw23.ge.1))then
                  call writeeclipse(Ncycle,Ttimes,Tseps,isw30,
     @               Nmaxeclipse,Tdur1,Tdur2)
                endif
c         
                if(icnU.ne.430)call wlinmod(Nphase,xmod,ymodU,
     @            'modelU.linear',isw7)
                if(icnB.ne.430)call wlinmod(Nphase,xmod,ymodB,
     @            'modelB.linear',isw7)
                if(icnV.ne.430)call wlinmod(Nphase,xmod,ymodV,
     @            'modelV.linear',isw7)
                if(icnR.ne.430)call wlinmod(Nphase,xmod,ymodR,
     @            'modelR.linear',isw7)
                if(icnI.ne.430)call wlinmod(Nphase,xmod,ymodI,
     @            'modelI.linear',isw7)
                if(icnJ.ne.430)call wlinmod(Nphase,xmod,ymodJ,
     @            'modelJ.linear',isw7)
                if(icnH.ne.430)call wlinmod(Nphase,xmod,ymodH,
     @            'modelH.linear',isw7)
                if(icnK.ne.430)call wlinmod(Nphase,xmod,ymodK,
     @            'modelK.linear',isw7)
c
                if(icnU.ne.430)call wmagmod(Nphase,xmod,ymodU,
     @            'modelU.mag',isw7,zeroU)
                if(icnB.ne.430)call wmagmod(Nphase,xmod,ymodB,
     @            'modelB.mag',isw7,zeroB)
                if(icnV.ne.430)call wmagmod(Nphase,xmod,ymodV,
     @            'modelV.mag',isw7,zeroV)
                if(icnR.ne.430)call wmagmod(Nphase,xmod,ymodR,
     @            'modelR.mag',isw7,zeroR)
                if(icnI.ne.430)call wmagmod(Nphase,xmod,ymodI,
     @            'modelI.mag',isw7,zeroI)
                if(icnJ.ne.430)call wmagmod(Nphase,xmod,ymodJ,
     @            'modelJ.mag',isw7,zeroJ)
                if(icnH.ne.430)call wmagmod(Nphase,xmod,ymodH,
     @            'modelH.mag',isw7,zeroH)
                if(icnK.ne.430)call wmagmod(Nphase,xmod,ymodK,
     @            'modelK.mag',isw7,zeroK)

                call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
                call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
                call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
c
                if(iidint.ge.1)call wlinmod(Nphase,xmod,ymodd,
     @            'lcdisk.linear',isw7)
                call wRVmod(NRVphase,xRVmod,RV1,'star1.RV',isw7,ggamma1)
                call wRVmod(NRVphase,xRVmod,RV2,'star2.RV',isw7,ggamma2)
                if(isw30.ge.3)call wRVmod(NRVphase,xRVmod,RV3,
     @            'star3.RV',isw7,ggamma1)
c
                
                if(itime.gt.0)then
                  if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,
     @                ydataU,errU,'ELCdataU.fold')
                  if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,
     @                ydataB,errB,'ELCdataB.fold')
                  if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,
     @                ydataV,errV,'ELCdataV.fold')
                  if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,
     @                ydataR,errR,'ELCdataR.fold')
                  if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,
     @                ydataI,errI,'ELCdataI.fold')
                  if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,
     @                ydataJ,errJ,'ELCdataJ.fold')
                  if(icnH.ne.430)call wELCdata(Ndatamax,NdataH,xdataH,
     @                ydataH,errH,'ELCdataH.fold')
                  if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,
     @                ydataK,errK,'ELCdataK.fold')
                  if(icnRV1.ne.430) call wELCdata(Ndatamax,NRV1,xRV1,
     @                yRV1,errRV1,'ELCdataRV1.fold')
                  if(icnRV2.ne.430) call wELCdata(Ndatamax,NRV2,xRV2,
     @                yRV2,errRV2,'ELCdataRV2.fold')
                  if(icnRV3.ne.430) call wELCdata(Ndatamax,NRV3,xRV3,
     @                yRV3,errRV3,'ELCdataRV3.fold')
c
                endif

                if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,
     @              resU,errU,'ELCresidualsU.fold')
                if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,
     @              resB,errB,'ELCresidualsB.fold')
                if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,
     @              resV,errV,'ELCresidualsV.fold')
                if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,
     @              resR,errR,'ELCresidualsR.fold')
                if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,
     @              resI,errI,'ELCresidualsI.fold')
                if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,
     @              resJ,errJ,'ELCresidualsJ.fold')
                if(icnH.ne.430)call wELCdata(Ndatamax,NdataH,xdataH,
     @              resH,errH,'ELCresidualsH.fold')
                if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,
     @              resK,errK,'ELCresidualsK.fold')
                if(icnRV1.ne.430) call wELCdata(Ndatamax,NRV1,xRV1,
     @              resRV1,errRV1,'ELCresidualsRV1.fold')
                if(icnRV2.ne.430) call wELCdata(Ndatamax,NRV2,xRV2,
     @              resRV2,errRV2,'ELCresidualsRV2.fold')
                if(icnRV3.ne.430) call wELCdata(Ndatamax,NRV3,xRV3,
     @              resRV3,errRV3,'ELCresidualsRV3.fold')
c
                call writegridout(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,
     &           omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,
     @           betarim,rinner,router,tdisk,xi,Ntheta,Nradius,alb1,
     @           alb2,Nref,rLx,Period,fm,separ,gamma1,t3,g3,SA3,density,
     @           sw1,sw2,sw3,T0,idraw,iecheck,iidint,iatm,ism1,icnU,
     @           icnB,icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,isw2,
     @           isw3,isw4,ilaw,wave,dbolx,dboly,dwavex,dwavey,ecc,
     @           argper,pshift,sw5,sw6,sw7,sw8,sw9,ikeep,isynch,isw5,
     @           isw6,isw7,isw8,isw9,spot1parm,spot2parm,spotdparm,
     @           primmass,primK,primrad,ratrad,frac1,frac2,ecosw,
     @           temprat,idark1,idark2,isw12,isw13,isw21,isw22,isw23,
     @           isw24,bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,
     @           sw28,sw29,sw30,contam,Tconj,beam1,beam2,isw25,isw26,
     @           isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,ocose,
     @           osine,omegadot,contamS0,contamS1,contamS2,contamS3,
     @           sw47,sw48,sw49)
c
                if(isw30.gt.0)call writebody3grid(Nalph3,Nbet3,
     @             tertperiod,tertt0,tertecos,tertesin,tertincl,
     @             tertOmega,tertQ,dwavex,dwavey,itconj,it1,it2,it3,it4,
     @             tertconj,tertratrad,hh,sw72,sw73,
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

                do 754 jjj=1,Nterms
                  Nstepnew(jjj)=idnint(dabs(vstart(jjj)-vstep(jjj))/
     %              (2.0d0*savestep(jjj)))
                  if(Nstepnew(jjj).eq.0)nstepnew(jjj)=1
754             continue

                call recordloopopt(Udatafile,Bdatafile,Vdatafile,
     @            Rdatafile,Idatafile,Jdatafile,Hdatafile,Kdatafile,
     @            RV1file,RV2file,Nvmax,Nvar,svar,var,vstart,savestep,
     @            Nstepnew,Nobv,sobv,obv,eobv,vstep)
c
                write(command,101)ichain+1000-1
                call system(command)
                write(command,102)ichain+1000-1
                call system(command)
                if(isw30.gt.0)write(command,103)ichain+1000-1
                call system(command)
              endif     !end if ismall > ireport
c
10000       continue
c
            close(45)
            close(46)
            if(isw24.gt.0)close(47)
            if(isw30.gt.0)close(48)
            if(isw30.ge.3)close(49)
            close(55)
c
c   output stuff about each chain
c
            do 9998 j=1,Nterms
              if(j.lt.10)write(chainstring,9997)j
              if(j.ge.10)write(chainstring,9996)j
              do 9995 ii=1,iparm(j)-1
                yscr(ii)=chiparm(j,ii)
                zscr(ii)=dble(ii)
 9995         continue

c
c  UPDATE JULY 5, 2011
c
c  Check to see that iparm(j) is more than 2
c
              if(iparm(j).ge.3)call sort2(iparm(j)-1,yscr,zscr)
c
              rrrmed=yscr((iparm(j)-1)/2)
c
              open(unit=55,file='markovchain'//chainstring//'.'//
     @          extension,status='unknown')
c
              iflag=0
              do 9994 ii=1,iparm(j)-1
                if(chiparm(j,ii).lt.rrrmed)iflag=iflag+1
                if(iflag.ge.10)write(55,9993)saveparm(j,ii)
9994          continue
              close(55)
9998        continue
c
            do 8995 ii=1,ipstringcount-1
              yscr(ii)=pstringchi(ii)
              zscr(ii)=dble(ii)
8995        continue

            if(ipstringcount.ge.3)call sort2(ipstringcount-1,yscr,zscr)
  
            rrrmed=yscr((ipstringcount-1)/2)

            open(unit=56,file='markovchainparm.'//extension,
     $         status='unknown')
            iflag=0
            do 8994 ii=1,ipstringcount-1
              if(pstringchi(ii).lt.rrrmed)iflag=iflag+1
              if(iflag.ge.10)write(56,8993)pstringparm(ii)
8994        continue
            close(56)
c
9999      continue   !continue loop over chains
c
3333      format(i7)
96        format('variable #',i2,' = ',a2,', step size = ',f15.8)
4323      format('error:  initial value of variable ',a2, 
     $        ' is out of range')
4324      format(3(f16.8,3x))
1120      format('frac  ',f7.5,', old step for parameter #',
     %             i2' = ',f15.9,' new step = ',f15.9)
996       format('variable #',i1,1x,'(',a2,'),'1x,a$)
1996      format('variable #',i2,1x,'(',a2,'),'1x,a$)
95        format(/a,', jump taken for parameter #',i1) 
1195      format(/a,', jump taken for parameter #',i2) 
2991      format(/a,',',a)
4321      format('error on variable ',i2,' (',a2,
     &       '):  lower bound is not less ','than the upper bound')
4322      format(24x,2(f16.9,2x))
 101      format('cp gridELC.opt gridELC.',i4)
 102      format('cp gridELC.inp ELC.',i4)
 103      format('cp gridELCbody3.inp ELCbody3.',i4)
8993      format(a)
9993      format(f18.9)
9997      format('0',i1)
9996      format(i2)
92        format(a,', ran=',f8.6,', probability=',f8.6,',',
     $            /'     jump not taken for parameter #',i1)
292       format(a,', ran=',f8.6,', probability=',f8.6,',',
     $            /'     jump not taken for parameter #',i2)
93        format(a,', ran=',f7.5,', probability=',f7.5,',',
     @             /' jump taken for parameter #',i2)
c
          end
c
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
          subroutine rnkpop(n,arr,indx,irank)
c
          integer n,indx(n),m,nstack,irank(n)
          real*8 arr(n)
          parameter(m=7,nstack=50)
          integer i
c
          call indexx(n,arr,indx)
c
          do 1 i=1,n
            irank(indx(i))=n-i+1
 1        continue
          return
          end
c
c    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine select(np,jfit,fdif,idad)
c
          implicit double precision(a-h,o-z)
c
          dimension jfit(np)
c
          common /ranblock/ idum
c
          np1=np+1
          urand=ran1(idum)
          dice=urand*dble(np*np1)
          rtfit=0.0d0
          do 1 i=1,np
            rtfit=rtfit+dble(np1)+fdif*dble(np1-2*jfit(i))
            if(rtfit.ge.dice)then
              idad=i
              go to 2
            endif
 1        continue
 2        return
          end
c
c  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
          subroutine encode(n,nd,ph,ign)
c
          double precision ph(n),z
          dimension ign(nd*n)
c
          z=10.0d0**nd
          ii=0
c
          do 1 i=1,n
            ip=idint(ph(i)*z)
            do 2 j=nd,1,-1
              ign(ii+j)=mod(ip,10)
              ip=ip/10
 2          continue
            ii=ii+nd
 1        continue
          return
          end
c
c   *************************************
c
          subroutine cross(n,nd,pcross,ign1,ign2)
c
c          double precision pcross,urand
          implicit double precision (a-h,o-z)

          dimension ign1(nd*n),ign2(nd*n)
c
          common /ranblock/ idum
c
          urand=ran1(idum)
          if(urand.lt.pcross)then
            urand=ran1(idum)
            ispl=idint(urand*dble(n*nd))+1
            do 10 i=ispl,n*nd
              it=ign2(i)
              ign2(i)=ign1(i)
              ign1(i)=it
 10         continue
          endif
          return
          end
c
c   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
          subroutine mutate(n,nd,pmut,ign)
c
c          double precision pmut,urand
          implicit double precision(a-h,o-z)
          dimension ign(n*nd)
c
          common /ranblock/ idum

          do 10 i=1,n*nd
            urand=ran1(idum)
            if(urand.lt.pmut)then
              urand=ran1(idum)
              ign(i)=int(dint(urand*10.0d0))
            endif
 10       continue
          return
          end
c
c   ##################################
c
          subroutine decode(n,nd,ign,ph)
c
          dimension ign(n*nd)
          double precision z,ph(n)
c
          z=10.0d0**(-nd)
          ii=0
          do 1 i=1,n
            ip=0
            do 2 j=1,nd
              ip=10*ip+ign(ii+j)
 2          continue
            ph(i)=dble(ip)*z
            ii=ii+nd
 1        continue
          return
          end
c
c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
          subroutine genrep(ndim,n,np,ip,ph,rnewph)
c
          double precision ph(ndim,2),rnewph(ndim,np)
c
          i1=2*ip-1
          i2=i1+1
          do 1 k=1,n
            rnewph(k,i1)=ph(k,1)
            rnewph(k,i2)=ph(k,2)
 1        continue
          return
          end
c
c   &&&&&&&&&&&&&&&&&&&&
c 
          subroutine adjmut(np,fitns,ifit,pmutmn,pmutmx,pmut)
c
          implicit double precision (a-h,o-z)
          double precision pmutmn,pmutmx,pmut,fitns,rat,rprod
c
          dimension fitns(np),ifit(np)
c
          rdiflo=0.05d0
          rdifhi=0.25d0
          delta=1.5d0
c
          rdif=dabs(fitns(ifit(np))-fitns(ifit(np/2)))/
     $      dabs(fitns(ifit(np))+fitns(ifit(np/2)))
c
          if(rdif.le.rdiflo)then
            rprod=pmut*delta
c            pmut=min(pmutmx,prod)
            if(pmutmx.lt.rprod)pmut=pmutmx
            if(rprod.le.pmutmx)pmut=rprod 
          endif
          if(rdif.ge.rdifhi)then
            rat=pmut/delta 
c            pmut=max(pmutmn,rat)
            if(pmutmn.gt.rat)pmut=pmutmn
            if(rat.ge.pmutmn)pmut=rat
          endif
c
          return
          end
c
c #!$!%(@^%#%*(&(#@#*$^@^$##*#@&(%)
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
          include 'dynamicssubs.for'
