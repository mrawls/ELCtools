         program mqELC
c
c   December 24, 1999
c
c   This is an optimizer code based on the Levenberg-Marquardt method.
c   The code is from Bevington (1969).
c
c 
c   UPDATE January 15, 2002
c
c   Update the code so that alb1 and alb2 can be adjusted (albedos of
c   star 1 and star 2).  Update subroutines assignvar, varassign,
c   newwritevar, and writevar
c
          implicit double precision (a-h,o-z)
          parameter(Nmaxphase=500000,Ndatamax=200000)
          parameter(Nvmax=60)
          parameter(Nmaxeclipse=1000)
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
          dimension ymodU1(Nmaxphase),ymodB1(Nmaxphase),
     $      ymodV1(Nmaxphase),ymodR1(Nmaxphase),ymodI1(Nmaxphase),
     $      ymodJ1(Nmaxphase),ymodH1(Nmaxphase),ymodK1(Nmaxphase),
     &      RV11(Nmaxphase),RV21(Nmaxphase)
          dimension ymodU2(Nmaxphase),ymodB2(Nmaxphase),
     $      ymodV2(Nmaxphase),ymodR2(Nmaxphase),ymodI2(Nmaxphase),
     $      ymodJ2(Nmaxphase),ymodH2(Nmaxphase),ymodK2(Nmaxphase),
     &      RV12(Nmaxphase),RV22(Nmaxphase)
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
     @      xRV3(Ndatamax),yRV3(Ndatamax),errRV3(Ndatamax)
         dimension riden(Nvmax,Nvmax),beta(Nvmax),alpha(Nvmax,Nvmax),
     $      b(Nvmax),sigmaa(Nvmax),array(Nvmax,Nvmax),
     @      deriv(Nmaxphase,Nvmax),drv1(Nmaxphase),drv2(Nmaxphase)
          dimension wave(8),dbolx(8,2),dboly(8,2),
     @      dwavex(8,3),dwavey(8,3)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax)
          dimension xRVmod(Nmaxphase),vstep1(Nvmax),resRV2(Ndatamax),
     @         resRV3(Ndatamax)
          dimension resU(Ndatamax),resB(Ndatamax),resV(Ndatamax)
          dimension resR(Ndatamax),resI(Ndatamax),resJ(Ndatamax)
          dimension resH(Ndatamax),resK(Ndatamax),resRV1(Ndatamax)
          character*40 Udatafile,svar(Nvmax),Hdatafile,Kdatafile,RV1file
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          character*40 RV2file  
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
          character*1000 parmstring
          character*2000 planetparm
          character*1700 dynparm,line
c                                       
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
     @      saverrH(Ndatamax),savydataK(Ndatamax),saverrK(Ndatamax),
     @      savyRV1(Ndatamax),saverrRV1(Ndatamax),savyRV2(Ndatamax),
     @      saverrRV2(Ndatamax),savxdataB(Ndatamax),savxdataV(Ndatamax),
     &      savxdataR(Ndatamax),savxdataI(Ndatamax),savxdataJ(Ndatamax),
     &      savxdataH(Ndatamax),savxdataK(Ndatamax),savxRV1(Ndatamax),
     @      savxRV2(Ndatamax)
c
c   UPDATE JULY 7, 2004
c
c   This array is needed for the sub-random Sobel sequence, used in
c   place of ran9
c
          dimension xsob(2)
          dimension gaplow(9999),gaphigh(9999)
c
c   UPDATE May 8, 2006
c
c   Add this array for power-law limb darkening coefficients.  
c   8=filter index, 9=coefficient c1,c2,c3 ...
c
          dimension powercoeff(8,9)
c
          dimension compfracs(8,3),ochidisk(8),sss(Nvmax)
          dimension fracs1(Nmaxphase,4),fracs2(Nmaxphase,4)
          dimension fracs3(Nmaxphase,4),fracs4(Nmaxphase,4)
          dimension fracs5(Nmaxphase,4),fracs6(Nmaxphase,4)
          dimension fracs7(Nmaxphase,4),fracs8(Nmaxphase,4)
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
c     @     Teff1,Teff2,Tgrav1,Tgrav2,betarim,rinner,router,tdisk,xi,
c     &      alb1,alb2,rLx,Period,fm,separ,gamma,wave,dbolx,dboly,dwavex,
c     @      dwavey,t3,g3,SA3,density,sw1,sw2,sw3,T0,ecc,argper,pshift,
c     @      sw5,sw6,sw7,sw8,sw9,primmass,primK,primrad,ratrad,frac1,
c     @      frac2,ecosw,temprat,bigI,bigbeta,sw23,sw24,powercoeff,sw25,
c     @      sw26,sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,ocose,
c     @      osine,sw42,contamS0,contamS1,contamS2,contamS3,sw47,sw48,
c     @      sw49,gaplow,gaphigh,
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
     @      atmint4,atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin
c
          common /intatm/  Nlines,Nmu,Nalph3,Nbet3,itconj,it1,it2,it3,
     @      it4
c
          common /medblock/ rmed 
c
c
c   Open the parameter file and read all of the parameters. 
c   Pass the parameters to the light curve routines
c   via a common block.  
c
        ifastflag=0   !disable fast genetic mode
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
     @       isw29,isw30,isw31,isw32,isw33,isw34,ocose,osine,sw42,
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
     &         atmint1,atmint2,atmint3,atmint4,atmint5,atmint6,atmint7,
     &         atmint8,Tmax,Tmin,gmax,gmin)
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
c   August 2, 2001
c
c   Use these flags when fitting for the period and/or T0
c
          itime=isw7
c
c   Disable median fitting.
c
          rmed=0.0d0
c
c   Disable the draw option
c
          idraw=0
c
          ifixgamma=isw4
c
c   Fix the inner disk radius at fill2.
c
          if(teff2.gt.0.0d0)rinner=fill2
c
c   If the value of savesep is negative, then the value of separ
c   should be computed from fm, finc, Q, and the period.  If savesep
c   is negative, then set separ=savesep so that the separation will
c   be computed correctly inside the lightcurve subroutine.
c
              savesep=separ
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
c   UPDATE JULY 7, 2004
c
c   Initialize the Sobel sequence here.
c
          nnn=-11
          call sobseq(nnn,xsob)
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
c   August 2, 2001
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
     @        NRV3,xRV3,yRV3,errRV3,Ndatamax,icnRV3,Nmaxeclipse)
          endif
c
          isvel1=0
          isvel2=0
          if(icnRV1.ne.430)isvel1=299
          if(icnRV2.ne.430)isvel2=299
          gamma1=gamma
          gamma2=gamma
c
          open(unit=45,file='generation.0000',status='unknown')
          open(unit=46,file='ELCparm.0000',status='unknown')
          if(isw24.ge.1)open(unit=47,file='ELCratio.0000',
     @      status='unknown')
          if(isw30.ge.1)open(unit=48,file='ELCbody3parm.0000',
     @      status='unknown')
          if(isw30.ge.3)open(unit=49,file='ELCdynparm.0000',
     @      status='unknown')
          open(unit=55,file='chi.0000',status='unknown')
c
c   Start the grid search here.  This code is more or less from Bevington
c   (1969).
c
          chi1=0.0d0
          chi2=0.0d0
          chi3=0.0d0
          flamda=0.001d0
          small=1.0d30
          imod=0
          Ndattot=NdataU+NdataB+NdataV+NdataR+NdataI+NdataJ+NdataH+
     @       NdataK+NRV1+NRV2+NRV3
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
 749      if(Nterms.eq.0)go to 69   !no valid variables specified
c
          write(*,*)'Nterms, Nvar, Ndata ',Nterms,Nvar,Ndattot
          do 750 ijk=1,Nstep(1)
c
            do 34 j=1,nterms
              beta(j)=0.0d0
              do 34 k=1,j
 34         alpha(j,k)=0.0d0
c
c  Find the chi square at first point.
c
            chilimb=0.0d0
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
     @         thresh,small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,
     @         obsTerr,icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,
     @         ggamma3,NRV3,
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
            imod=imod+1
            chi1=chi1/dabs(dble(Ndattot-Nterms))
            call printiter(imod,'model number = ',ijk,
     @        'iteration number = ')
c
            if(chi1.lt.small)then
              small=chi1
              do mmm=1,Nvmax
                sss(mmm)=var(mmm)
              enddo
            endif
c 
            do 50 j=1,nterms
              delta=vstep(j) 
              oldvar=var(j)
              var(j)=var(j)+delta
c
              ifastflag=0
              ibest=0
              ichilabel=2
              call monster(Nphase,Nmaxphase,xmod,ymodU1,ymodB1,ymodV1,
     @           ymodR1,ymodI1,ymodJ1,ymodH1,ymodK1,ymods1,ymods2,ymods3,
     @           ymodd,RV11,RV21,drv1,drv2,obsparm,ifastflag,NRVphase,
     @           xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     @           fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     @           chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi2,
     @           NdataU,xdataU,ydataU,errU,zeroU,resU,NdataB,xdataB,
     @           ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,errV,zeroV,
     @           resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
     @           xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,
     &           errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,
     %           NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,
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
     @           ifixgamma,savesep,ichilabel,thresh,small,
     @           Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     @           icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,
     @           NRV3,parmstring,planetparm,dynparm,line,
     @    fill1,fill2,omega1,
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

              chi2=chi2/dabs(dble(Ndattot-Nterms))
              call printiter(imod,'model number = ',ijk,
     @          'iteration number = ')
c
              if(chi2.lt.small)then
                small=chi2
                do mmm=1,Nvmax
                  sss(mmm)=var(mmm)
                enddo
              endif
c
              var(j)=oldvar
              var(j)=var(j)-delta
c
c
              ifastflag=0
              ibest=0
              ichilabel=3
              call monster(Nphase,Nmaxphase,xmod,ymodU2,ymodB2,ymodV2,
     @           ymodR2,ymodI2,ymodJ2,ymodH2,ymodK2,ymods1,ymods2,ymods3,
     @           ymodd,RV12,RV22,drv1,drv2,obsparm,ifastflag,NRVphase,
     @           xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     @           fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     @           chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi3,
     @           NdataU,xdataU,ydataU,errU,zeroU,resU,NdataB,xdataB,
     @           ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,errV,zeroV,
     @           resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
     @           xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,
     &           errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,
     %           NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,
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
     @           ifixgamma,savesep,ichilabel,thresh,small,
     @           Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
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
              imod=imod+1
              chi3=chi3/dabs(dble(Ndattot-Nterms))
              call printiter(imod,'model number = ',ijk,
     @          'iteration number = ')
c
              if(chi3.lt.small)then
                small=chi3
                do mmm=1,Nvmax
                  sss(mmm)=var(mmm)
                enddo
              endif
c
c   Get the derivatives.
c
              jcount=0
              if(icnU.ne.430)then 
                do 6677 m=1,NdataU
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodU,NdataU,xdataU,
     @                      m,yout,zeroU)
                  call evalfit(9,Nphase,xmod,ymodU1,NdataU,xdataU,
     @                      m,yout1,zeroU1)
                  call evalfit(9,Nphase,xmod,ymodU2,NdataU,xdataU,
     @                      m,yout2,zeroU2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     %                (ydataU(m)-yout)*deriv(jcount,j)/(errU(m)*errU(m))
                  do 6646 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %               deriv(jcount,j)*deriv(jcount,k)/(errU(m)*errU(m))
 6646             continue
 6677           continue
              endif
c
              if(icnB.ne.430)then 
                do 6678 m=1,NdataB
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodB,NdataB,xdataB,
     @                       m,yout,zeroB)
                  call evalfit(9,Nphase,xmod,ymodB1,NdataB,xdataB,
     @                       m,yout1,zeroB1)
                  call evalfit(9,Nphase,xmod,ymodB2,NdataB,xdataB,
     @                       m,yout2,zeroB2)
c
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataB(m)-yout)*deriv(jcount,j)/(errB(m)*errB(m))
                  do 6647 k=1,j
                    alpha(j,k)=alpha(j,k)+
     $                deriv(jcount,j)*deriv(jcount,k)/(errB(m)*errB(m))
 6647             continue
 6678           continue
              endif
c
              if(icnV.ne.430)then 
                do 6679 m=1,NdataV
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodV,NdataV,xdataV,
     $                   m,yout,zeroV)
                  call evalfit(9,Nphase,xmod,ymodV1,NdataV,xdataV,
     &                   m,yout1,zeroV1)
                  call evalfit(9,Nphase,xmod,ymodV2,NdataV,xdataV,
     @                   m,yout2,zeroV2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     %            (ydataV(m)-yout)*deriv(jcount,j)/(errV(m)*errV(m))
                  do 6648 k=1,j
                    alpha(j,k)=alpha(j,k)+
     $              deriv(jcount,j)*deriv(jcount,k)/(errV(m)*errV(m))
 6648             continue
 6679           continue
              endif
c
              if(icnR.ne.430)then 
                do 6680 m=1,NdataR
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodR,NdataR,xdataR,
     &                  m,yout,zeroR)
                  call evalfit(9,Nphase,xmod,ymodR1,NdataR,xdataR,
     @                  m,yout1,zeroR1)
                  call evalfit(9,Nphase,xmod,ymodR2,NdataR,xdataR,
     @                 m,yout2,zeroR2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataR(m)-yout)*deriv(jcount,j)/(errR(m)*errR(m))
                  do 6649 k=1,j
                    alpha(j,k)=alpha(j,k)+
     $               deriv(jcount,j)*deriv(jcount,k)/(errR(m)*errR(m))
 6649             continue
 6680           continue
              endif
c
              if(icnI.ne.430)then 
                do 6681 m=1,NdataI
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodI,NdataI,xdataI,
     @                       m,yout,zeroI)
                  call evalfit(9,Nphase,xmod,ymodI1,NdataI,xdataI,
     @                       m,yout1,zeroI1)
                  call evalfit(9,Nphase,xmod,ymodI2,NdataI,xdataI,
     @                       m,yout2,zeroI2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataI(m)-yout)*deriv(jcount,j)/(errI(m)*errI(m))
                  do 6650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %               deriv(jcount,j)*deriv(jcount,k)/(errI(m)*errI(m))
6650              continue
6681            continue
              endif
c
              if(icnJ.ne.430)then 
                do 16681 m=1,NdataJ
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodJ,NdataJ,xdataJ,
     @                      m,yout,zeroJ)
                  call evalfit(9,Nphase,xmod,ymodJ1,NdataJ,xdataJ,
     @                      m,yout1,zeroJ1)
                  call evalfit(9,Nphase,xmod,ymodJ2,NdataJ,xdataJ,
     @                      m,yout2,zeroJ2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataJ(m)-yout)*deriv(jcount,j)/(errJ(m)*errJ(m))
                  do 16650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %               deriv(jcount,j)*deriv(jcount,k)/(errJ(m)*errJ(m))
16650           continue
16681         continue
              endif
c
              if(icnK.ne.430)then 
                do 26681 m=1,NdataK
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodK,NdataK,xdataK,
     @                   m,yout,zeroK)
                  call evalfit(9,Nphase,xmod,ymodK1,NdataK,xdataK,
     @                   m,yout1,zeroK1)
                  call evalfit(9,Nphase,xmod,ymodK2,NdataK,xdataK,
     @                   m,yout2,zeroK2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataK(m)-yout)*deriv(jcount,j)/(errK(m)*errK(m))
                  do 26650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %               deriv(jcount,j)*deriv(jcount,k)/(errK(m)*errK(m))
26650           continue
26681         continue
              endif
c
              if(icnH.ne.430)then 
                do 36681 m=1,NdataH
                  jcount=jcount+1
                  call evalfit(9,Nphase,xmod,ymodH,NdataH,xdataH,
     @                   m,yout,zeroH)
                  call evalfit(9,Nphase,xmod,ymodH1,NdataH,xdataH,
     @                   m,yout1,zeroH1)
                  call evalfit(9,Nphase,xmod,ymodH2,NdataH,xdataH,
     @                   m,yout2,zeroH2)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $             (ydataH(m)-yout)*deriv(jcount,j)/(errH(m)*errH(m))
                  do 36650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %               deriv(jcount,j)*deriv(jcount,k)/(errH(m)*errH(m))
36650           continue
36681         continue
              endif
c 
              if(icnRV1.ne.430)then 
                do 46681 m=1,NRV1
                  jcount=jcount+1
                  call evalfit(0,Nphase,xmod,RV1,NRV1,xRV1,
     @                    m,yout,gamma1)
                  call evalfit(0,Nphase,xmod,RV11,NRV1,xRV1,
     @                    m,yout1,gamma11)
                  call evalfit(0,Nphase,xmod,RV12,NRV1,xRV1,
     @                    m,yout2,gamma12)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $              (xRV1(m)-yout)*deriv(jcount,j)/(errRV1(m)*
     #              errRV1(m))
                  do 46650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %                 deriv(jcount,j)*deriv(jcount,k)/(errRV1(m)*
     %                 errRV1(m))
46650             continue
46681           continue
              endif
c
              if(icnRV2.ne.430)then 
                do 56681 m=1,NRV2
                  jcount=jcount+1
                  call evalfit(0,Nphase,xmod,RV2,NRV2,xRV2,
     @                  m,yout,gamma2)
                  call evalfit(0,Nphase,xmod,RV21,NRV2,xRV2,m,
     @                  yout1,gamma21)
                  call evalfit(0,Nphase,xmod,RV22,NRV2,xRV2,
     @                  m,yout2,gamma22)
                  d1=(yout1-yout)/delta
                  d2=(yout-yout2)/delta
                  d3=(yout1-yout2)/(2.0*delta)
                  deriv(jcount,j)=(d1+d2+d3)/3.0
                  beta(j)=beta(j)+
     $               (xRV2(m)-yout)*deriv(jcount,j)/(errRV2(m)*
     %               errRV2(m))
                  do 56650 k=1,j
                    alpha(j,k)=alpha(j,k)+
     %                 deriv(jcount,j)*deriv(jcount,k)/(errRV2(m)*
     &                 errRV2(m))
56650             continue
56681           continue
              endif
c
              var(j)=oldvar
 50         continue               !continue the loop over Nterms
c
c   we have filled up the beta and alpha matrices
c
            do 53 j=1,nterms
              do 53 k=1,j
 53         alpha(k,j)=alpha(j,k)
c
c   invert the matrices
c
 71         do 74 j=1,nterms
              do 73 k=1,nterms
                riden(j,k)=0.0d0
 73           array(j,k)=alpha(j,k)/sqrt(alpha(j,j)*alpha(k,k))
 74         array(j,j)=(1.0d0+flamda) 
            do 652 j=1,nterms
              riden(j,j)=1.0d0
 652        continue
            call gaussj(array,Nterms,Nvmax,riden,Nterms)
c
            do 84 j=1,nterms
              b(j)=var(j)
              do 84 k=1,nterms
 84         b(j)=b(j)+beta(k)*array(j,k)/sqrt(alpha(j,j)*alpha(k,k))
c
c   Evaluate the function and chi^2.  Note the array 'b' is used instead of
c   'var'.
c 
            ifastflag=0
            ibest=0
            ichilabel=4
            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     @         ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @         ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @         xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,
     @         fracs7,fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,
     @         chisqJ,chisqH,chisqK,chisqRV1,chisqRV2,chilimb,chi4,
     @         NdataU,xdataU,ydataU,errU,zeroU,resU,NdataB,xdataB,
     @         ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,errV,
     @         zeroV,resV,NdataR,xdataR,ydataR,errR,zeroR,resR,
     @         NdataI,xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,
     @         ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,
     %         zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,
     @         xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,
     @         ggamma2,Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,
     &         Nvmax,svar,var,saveasini,savxdataU,savydataU,saverrU,
     &         savxdataB,savydataB,saverrB,savxdataV,savydataV,
     @         saverrV,savxdataR,savydataR,saverrR,savxdataI,
     &         savydataI,saverrI,savxdataJ,savydataJ,saverrJ,
     @         savxdataH,savydataH,saverrH,savxdataK,savydataK,
     &         saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,savyRV2,
     @         saverrRV2,ifrac,ilum,i16,isavNU,isavNB,isavNV,isavNR,
     @         isavNI,isavNJ,isavNH,isavNK,isavRV1,isavRV2,isvel1,
     &         isvel2,Ndatamax,ibest,ifixgamma,savesep,ichilabel,thresh,
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
            chi4=chi4/dabs(dble(Ndattot-Nterms))
            call printiter(imod,'model number = ',ijk,
     @          'iteration number = ')
c
            if(chi4.lt.small)then
              small=chi4
              do mmm=1,Nvmax
                sss(mmm)=b(mmm)
              enddo
            endif
c
            diff=chi1-chi4
c
            if(diff) 95, 101, 101
 95         flamda=10.0d0*flamda
            write(*,*)'flamda = ',flamda
            if(flamda.gt.1.0d12)stop
            go to 71
c
 101        write(*,*)' '
c
            do 103 j=1,nterms
              var(j)=b(j)
              var(j)=sss(j)
              vstep(j)=0.9*vstep(j)
              sigmaa(j)=sqrt(array(j,j)/alpha(j,j))
              write(*,1001)var(j),sigmaa(j),svar(j)(1:15)
 103        continue
c
            flamda=flamda/10.0
            write(*,*)'flamda = ',flamda
c
            call printS(small)
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
            call wRVmod(NRVphase,xRVmod,RV1,'star1.RV',isw7,ggamma1)
            call wRVmod(NRVphase,xRVmod,RV2,'star2.RV',isw7,ggamma2)
            if(isw30.ge.3)call wRVmod(NRVphase,xRVmod,RV3,'star3.RV',
     @          isw7,ggamma1)
            call wlinmod(NRVphase,xRVmod,dRV1,'star1.delRV',isw7)
            call wlinmod(NRVphase,xRVmod,dRV2,'star2.delRV',isw7)
c
c   If we are fitting for the period and/or the T0, write the current
c   folded light curves
c
            if(itime.gt.0)then
              if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,
     @              ydataU,errU,'ELCdataU.fold')
              if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,
     @              ydataB,errB,'ELCdataB.fold')
              if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,
     @              ydataV,errV,'ELCdataV.fold')
              if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,
     @              ydataR,errR,'ELCdataR.fold')
              if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,
     @              ydataI,errI,'ELCdataI.fold')
              if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,
     @              ydataJ,errJ,'ELCdataJ.fold')
              if(icnH.ne.430) call wELCdata(Ndatamax,NdataH,xdataH,
     @              ydataH,errH,'ELCdataH.fold')
              if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,
     @              ydataK,errK,'ELCdataK.fold')
              if(icnRV1.ne.430) call wELCdata(Ndatamax,NRV1,xRV1,yRV1,
     @              errRV1,'ELCdataRV1.fold')
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
 750      continue  ! loop over iterations
c
c   reset the variables at their best values and print the chi^2
c
          do mmm=1,Nvmax
            var(mmm)=sss(mmm)
          enddo
c 
          ifastflag=0
          ibest=0
          ichilabel=1
          call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,ymodR,
     @       ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,RV1,
     @       RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,fracs1,
     @       fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
     @       chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,
     @       chisqRV1,chisqRV2,chilimb,chiall,NdataU,xdataU,ydataU,errU,
     @       zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,
     &       xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,ydataR,errR,
     %       zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,resI,NdataJ,
     $       xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,
     @       zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,
     @       xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,
     @       Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,
     @       saveasini,savxdataU,savydataU,saverrU,savxdataB,
     @       savydataB,saverrB,savxdataV,savydataV,saverrV,savxdataR,
     @       savydataR,saverrR,savxdataI,savydataI,saverrI,savxdataJ,
     @       savydataJ,saverrJ,savxdataH,savydataH,saverrH,savxdataK,
     @       savydataK,saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,
     @       savyRV2,saverrRV2,ifrac,ilum,i16,isavNU,isavNB,isavNV,
     @       isavNR,isavNI,isavNJ,isavNH,isavNK,isavRV1,isavRV2,
     @       isvel1,isvel2,Ndatamax,ibest,ifixgamma,savesep,ichilabel,
     @       thresh,small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,
     @       obsTerr,icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,
     @       ggamma3,NRV3,parmstring,planetparm,dynparm,line,
     @    fill1,fill2,omega1,
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
           call printiter(imod,'model number = ',ijk,
     @        'iteration number = ')
c
          write(*,*)' '
c
          do 3749 iq=1,Nterms
            write(*,1001)var(iq),sigmaa(iq),svar(iq)(1:15)
 3749     continue
c
          call printS(chiall)
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
          if(iidint.ge.1)call wlinmod(Nphase,xmod,ymodd,
     @        'lcdisk.linear',isw7)
c
          call wRVmod(NRVphase,xRVmod,RV1,'star1.RV',isw7,ggamma1)
          call wRVmod(NRVphase,xRVmod,RV2,'star2.RV',isw7,ggamma2)
          if(isw30.ge.3)call wRVmod(NRVphase,xRVmod,RV3,'star3.RV',
     @         isw7,ggamma1)
          call wlinmod(NRVphase,xRVmod,dRV1,'star1.delRV',isw7)
          call wlinmod(NRVphase,xRVmod,dRV2,'star2.delRV',isw7)
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
     @       omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,
     @       betarim,rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,
     @       Nref,rLx,Period,fm,separ,gamma1,t3,g3,SA3,density,sw1,sw2,
     @       sw3,T0,idraw,iecheck,iidint,iatm,ism1,icnU,icnB,icnV,icnR,
     @       icnI,icnJ,icnH,icnK,iRVfilt,isw1,isw2,isw3,isw4,ilaw,wave,
     @       dbolx,dboly,dwavex,dwavey,ecc,argper,pshift,sw5,sw6,sw7,
     @       sw8,sw9,ikeep,isynch,isw5,isw6,isw7,isw8,isw9,spot1parm,
     @       spot2parm,spotdparm,primmass,primK,primrad,ratrad,frac1,
     @       frac2,ecosw,temprat,idark1,idark2,isw12,isw13,isw21,isw22,
     @       isw23,isw24,bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,
     @       sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,isw25,isw26,
     @       isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,ocose,
     @       osine,
     @       sw42,contamS0,contamS1,contamS2,contamS3,sw47,sw47,sw49)
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
     @      Nvmax,Nvar,svar,var,vstart,vstep,Nstep,Nobv,sobv,obv,
     @      eobv,vstep1)
c
 1001     format(1x,2(f15.7,2x),a15)

 69       close(45)
          close(46)
          if(isw24.gt.0)close(47)
          if(isw30.gt.0)close(48)
          if(isw30.ge.3)close(49)
          close(55)
c
          end
c
c   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
          include 'dynamicssubs.for'
