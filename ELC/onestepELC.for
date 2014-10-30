          program stepELC
c
c   June 29, 2014
c
c   better grid search, can be parallel
c
          implicit double precision (a-h,o-z)
c
          parameter(Nmaxphase=750000,Ndatamax=200000)
          parameter(Nvmax=60)
          parameter(Nmaxeclipse=600)
c
          dimension obsparm(19),obv(19),eobv(19),xxx(3),yyy(3)
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
     $      drv1(Nmaxphase),drv2(Nmaxphase),
     @      xRV3(Ndatamax),yRV3(Ndatamax),
     @      errRV3(Ndatamax)
          dimension wave(8),dbolx(8,2),dboly(8,2),dwavex(8,3),
     @       dwavey(8,3)
          dimension var(Nvmax),vstart(Nvmax),vstep(Nvmax),Nstep(Nvmax),
     %      stepsave(Nvmax)
          character*40 Udatafile,svar(nvmax),Hdatafile,Kdatafile,
     &             RV1file,RV2file
          character*40 Bdatafile,Vdatafile,Rdatafile,Idatafile,Jdatafile
          character*40 command
          dimension parmarray(Nvmax,10000),chiarr(10000),indxchi(10000)
          dimension dummy(Nvmax),oldchiarr(10000)
          dimension xRVmod(Nmaxphase),oldparmarray(Nvmax,10000)
          dimension resU(Ndatamax),resB(Ndatamax),resV(Ndatamax)
          dimension resR(Ndatamax),resI(Ndatamax),resJ(Ndatamax)
          dimension resH(Ndatamax),resK(Ndatamax),resRV1(Ndatamax)
          dimension resRV2(Ndatamax),resRV3(Ndatamax)
c
          character*7 extension
c
          parameter (maxlines=1300,maxmu=115)   !was 1100
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
c   UPDATE May 8, 2006
c
c   Add this array for power-law limb darkening coefficients.  
c   8=filter index, 9=coefficient c1,c2,c3 ...
c
          dimension powercoeff(8,9)
c
          character*1000 parmstring
          character*2000 planetparm      
          character*1700 dynparm,line
          character*2000 bigstring(1000,4)
          
          character*30 outstring,outstring1,outstring2
c
c   New variables are needed to allow for the fitting of period and T0.
c
          dimension savxdataU(Ndatamax),savydataU(Ndatamax),
     @      saverrU(Ndatamax),
     &      savydataB(Ndatamax),saverrB(Ndatamax),savydataV(Ndatamax),
     $      saverrV(Ndatamax),
     &      savydataR(Ndatamax),saverrR(Ndatamax),savydataI(Ndatamax),
     $      saverrI(Ndatamax),
     &      savydataJ(Ndatamax),saverrJ(Ndatamax),savydataH(Ndatamax),
     $      saverrH(Ndatamax),
     $      savydataK(Ndatamax),saverrK(Ndatamax),savyRV1(Ndatamax),
     $      saverrRV1(Ndatamax),
     %      savyRV2(Ndatamax),saverrRV2(Ndatamax),savxdataB(Ndatamax),
     $      savxdataV(Ndatamax),
     &      savxdataR(Ndatamax),savxdataI(Ndatamax),savxdataJ(Ndatamax),
     &      savxdataH(Ndatamax),savxdataK(Ndatamax),savxRV1(Ndatamax),
     $      savxRV2(Ndatamax)
c
c
c   UPDATE OCTOBER 10, 2007
c
c   Add this array to keep track of the light curves of the individual
c   components in all band passes.
c
          dimension compfracs(8,3),ochidisk(8)
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
c
          dimension gaplow(9999),gaphigh(9999)
c
          dimension sss(Nvmax),sssloop(Nvmax),ssssave(Nvmax)
c
          dimension Ncycle(34),Ttimes(34,Nmaxeclipse)
          dimension Tseps(34,Nmaxeclipse),Tdur1(34,Nmaxeclipse)
          dimension Nobscycle(34),obsTtimes(34,Nmaxeclipse)
          dimension obsTerr(34,Nmaxeclipse),Tdur2(34,Nmaxeclipse)
          dimension icnarray(34)
          dimension xSC(9999),ySC(9999),probarr(9999)
          dimension idxarr1(9999),idxarr2(9999),zzzarr(9999)

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
c     @     Nref,idraw,iecheck,iidint,iatm,ism1,ilaw,icnU,icnB,icnV,icnR,
c     @     icnI,icnJ,icnH,icnK,icnRV1,icnRV2,iRVfilt,isw1,isw2,isw3,
c     @     isw4,ikeep,isynch,isw5,isw6,isw7,isw8,isw9,idark1,idark2,
c     @     isw12,isw13,isw21,isw22,isw23,isw24,isw25,isw26,isw27,isw28,
c     @     isw29,isw30,isw31,isw32,isw33,isw34,NSC
cc
c          common /fracblock/ compfracs          
cc
c          common /thirdblock/ tertperiod,tertt0,tertecos,tertesin,
c     @        tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72,sw73
c
          common /realatm/ atmT,atmg,atmmu,atmint1,atmint2,atmint3,
     &      atmint4,atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin

          common /intatm/  Nlines,Nmu,Nalph3,Nbet3,itconj,it1,it2,it3,
     @      it4
c
          common /ranblock/ idum
c
          common /medblock/ rmed 
c
c  /realatm/, /intatm/, /intatm/, /medblock/

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
c   UPDATE May 8, 2006
c
c   Add isw21-isw24, sw21-sw24, powercoeff to list above.
c
c
c   UPDATE October 28, 2002
c
c   Define a new variable called sepsave, which is equal to the
c   value of separ given in ELC.inp
c
           savesep=separ
c
            if(ilaw.eq.4)then
              do jjj=1,8
                rlimbsum1=dwavex(jjj,1)+dwavey(jjj,1)
                if(rlimbsum1.gt.1.0d0)then
                  write(*,*)'the sum of x and y exceeds 1.0 for index ',
     @                jjj
                  stop
                endif
                rlimbsum2=dwavex(jjj,2)+dwavey(jjj,2)
                if(rlimbsum2.gt.1.0d0)then
                  write(*,*)'the sum of x and y exceeds 1.0 for index ',
     @                jjj
                  stop
                endif
              enddo
            endif
c
c   Load the atmosphere table here, rather than in the subroutine
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
c   UPDATE JULY 30, 2004
c
c   set Nbin = 0
c
c          Nbin=0
c
c   UPDATE May 27, 2002
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
c   UPDATE April 15, 2002
c
c   Define the values of gamma1 and gamma2
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
          isw9save=isw9
c
c   UPDATE JULY 7, 2004
c
c   Initialize the Sobel sequence here.
c
          mmm=-11
          call sobseq(mmm,xsob)
c
c   Initialize the steps.
c
          do 1 i=1,nvmax
            Nstep(i)=1
            var(i)=0.0d0
            svar(i)='none'
            sss(i)=0.0d0
            sssloop(i)=0.0d0
            ssssave(i)=0.0d0
 1        continue
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
4321      format('error on variable ',i2,' (',a2,
     &       '):  lower bound is not less ','than the upper bound')
4322      format(24x,2(f16.9,2x))
c
c   Add also the ability to include the luminosity ratio in the total chi^2
c
          ifrac=-99
          ilum=-99
          do 111 i=1,Nobv
            stepsave(i)=vstep(i)
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
          ilimbcheck=-55
          do j=1,Nvar
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
          enddo
c
          if(ilaw.ne.4)ilimbcheck=-99
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
     @       NRV3,xRV3,yRV3,errRV3,Ndatamax,icnRV3,Nmaxeclipse)
          endif
c
          isvel1=0
          isvel2=0

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
c   UPDATE November 6, 2008
c 
c   if sw25 > 0, use that as an error for asini.
c
          saveasini=sw5
c
c   Start the looping here.  
c
          chi1=0.0d0
          small=9.0d33
          smallloop=small
          Ndattot=0
          Ndattot=NdataU+NdataB+NdataV+NdataR+NdataI+NdataJ+NdataH+
     $        NdataK+NRV1+NRV2+NRV3
c
c
c   Open an output file for the fitting statistics.
c
          open(unit=45,file='generation.0000',status='unknown')
          if(isw24.ge.1)open(unit=47,file='ELCratio.0000',
     $         status='unknown')
          open(unit=55,file='chi.0000',status='unknown')
          open(unit=46,file='ELCparm.0000',status='unknown')
          if(isw30.ge.1)open(unit=48,file='ELCbody3parm.0000',
     @      status='unknown')
          if(isw30.ge.3)open(unit=49,file='ELCdynparm.0000',
     @      status='unknown')

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
          write(*,*)'Nterms, Nvar, Ndattot',Nterms,Nvar,Ndattot
c
c   
c   Define the random parameter sets.  Nstep(1) is the number of random
c   sets to define.
c

          Niter=10
          Nloop=7

c
c   open up optional file stepELC.inp and read
c   Niter
c   Nloop
c
c
          open(unit=44,file='stepELC.inp',status='old',err=2)
c
          read(44,*,end=6)Niter
          read(44,*,end=6)Nloop
6         close(44)
          go to 3
 
2         write(*,*)'stepELC.inp not found.  Using default values'
          write(*,*)'of Niter=10, Nloop=5'
c
c   Initialize the steps.
c


c          np=Nstep(1)     ! number of members in the population 
c          ngen=Nstep(2)   ! number of generations
c
c          do 99 j=1,Nterms
c            vx1=vstart(j)
c            vx2=vstep(j)
c            vxmult=(vx2-vx1)
c            vxadd=vx1
c            do 98 i=1,np                 
c              urand=ran1(idum)
c              varx=urand*vxmult+vxadd
c              parmarray(j,i)=varx
c 98         continue
c 99       continue
c
c   Evaluate the fitness of the initial population.
c
c
c            if(isw22.ge.1)ifastflag=1
c            if(i16.eq.1)ifastflag=0
c            do 760 j=1,Nterms
c              var(j)=parmarray(j,i16)
c 760        continue
c
c   UPDATE August 16, 2001
c
c   If the isw9 flag is 1 or greater, then alter the first population member
c   and put the parameters inside ELC.inp into the var array
c
c
 3        call varassign(Nvmax,svar,var,fill1,fill2,omega1,omega2,Q,
     @          finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,rLx,
     @          separ,gamma,t3,g3,sa3,ecc,argper,pshift,spot1parm,
     @          spot2parm,spotdparm,period,T0,alb1,alb2,dwavex,dwavey,
     @          primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @          bigI,bigbeta,powercoeff,density,Tconj,beam1,beam2,
     @          contam,ocose,osine,isw29,tertperiod,tertt0,tertecos,
     @          tertesin,tertincl,tertOmega,tertQ,Tgrav1,Tgrav2,
     @          tertconj,omegadot,
     @          contamS0,contamS1,contamS2,contamS3,
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
     @       P8ratrad)

          do j=1,Nterms
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
          enddo
c
4323      format('error:  initial value of variable ',a2, 
     $        ' is out of range')
4324      format(3(f16.8,3x))

          do jj=1,Nterms
            sss(jj)=var(jj)
            sssloop(jj)=var(jj)
            ssssave(jj)=var(jj)
          enddo
c
          imodel=0
c
c   loop over Nterms
c
          do 900 iloop=1,Nloop
          do 850 iterm=1,Nterms
            indexhi=vstep(iterm)
            indexlow=vstart(iterm)
            rup=vstep(iterm)
            rlow=vstart(iterm)
            rbest=var(iterm)
c
c   save the variables here
c
            do kk=1,Nterms
              sssloop(kk)=var(kk)
            enddo
            do 800 iout=1,Nstep(iterm)
            do ll=1,Niter/2
              rstep=((rbest-rlow)/dble(Niter/2-1))
              do j=1,Nterms
                varx=var(j)
                if(j.eq.iterm)then
                  varx=(dble(ll-1))*rstep+rlow
                endif
                parmarray(j,ll)=varx
              enddo
            enddo
            do ll=1,Niter/2
              rstep=((rup-rbest)/dble(Niter/2))
              do j=1,Nterms
                varx=var(j)
                if(j.eq.iterm)then
                  varx=(dble(ll))*rstep+rbest
                endif
                parmarray(j,ll+Niter/2)=varx
              enddo
            enddo
c
            if(iout.eq.Nstep(iterm))then
              do ll=1,Niter
                urand=ran1(idum)
                do j=1,Nterms
                  varx=var(j)
                  if(j.eq.iterm)then
                    varx=urand*(rup-rlow)+rlow
                  endif
                  parmarray(j,ll)=varx
                enddo
              enddo
            endif   
c
!$omp  parallel firstprivate(Nphase,xmod,ymodU,ymodB,ymodV,
!$omp&      ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,
!$omp@      RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,fracs1,
!$omp@      fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
!$omp@      chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,chisqRV1,
!$omp&      chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,zeroU,resU,
!$omp@      NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,
!$omp@      errV,zeroV,resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
!$omp@      xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,errJ,
!$omp%      zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,NdataK,
!$omp$      xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,errRV1,NRV2,
!$omp@      xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,sobv,obv,eobv,ochi,
!$omp&      ochidisk,ochilr,svar,var,saveasini,
!$omp@      savxdataU,savydataU,saverrU,savxdataB,savydataB,
!$omp@      saverrB,savxdataV,savydataV,saverrV,
!$omp$      savxdataR,savydataR,saverrR,savxdataI,savydataI,
!$omp&      saverrI,savxdataJ,savydataJ,saverrJ,
!$omp&      savxdataH,savydataH,saverrH,savxdataK,savydataK,
!$omp@      saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,savyRV2,saverrRV2,
!$omp@      ifrac,ilum,i16,isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,
!$omp@      isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,
!$omp@      ibest,ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
!$omp@      Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
!$omp@      icnarray,Tdur1,Tdur2,
!$omp@      RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
!$omp@      fill1,fill2,omega1,
!$omp@    omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
!$omp@    rinner,router,tdisk,xi,alb1,alb2,rLx,Period,fm,separ,
!$omp@    gamma,wave,dbolx,dboly,dwavex,dwavey,t3,g3,SA3,density,
!$omp@    sw1,sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
!$omp@    primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
!$omp@    bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,
!$omp@    sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,omegadot,
!$omp@    contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49,gaplow,
!$omp@    gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
!$omp@    P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
!$omp@    P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
!$omp@    P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
!$omp@    P5esin,P5incl,P5Omega,P5Q,P5ratrad,P6tconj,P6period,P6T0,
!$omp@    P6ecos,P6esin,P6incl,P6Omega,P6Q,P6ratrad,P7tconj,P7period,
!$omp@    P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,
!$omp@    P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,
!$omp@    xSC,ySC,spot1parm,spot2parm,spotdparm,Nalph1,Nbet1,Nalph2,
!$omp@    Nbet2,Ntheta,Nradius,Nref,idraw,iecheck,iidint,iatm,ism1,
!$omp@    ilaw,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
!$omp@    iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,isw5,isw6,isw7,
!$omp@    isw8,isw9,idark1,idark2,isw12,isw13,isw21,isw22,isw23,isw24,
!$omp@    isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,
!$omp@    NSC,compfracs,tertperiod,tertt0,tertecos,tertesin,
!$omp@    tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72,sw73)
!$omp@  shared(small,sss,sssloop,bigstring,chiarr,oldchiarr,nnn,
!$omp@  smallloop,indexsmall,imodel,ssssave)
!$omp@  private(parmstring,planetparm,dynparm,line)
!$omp@  shared(/realatm/, /intatm/, /medblock/)

cc!$omp@ shared(small,sss,bigstring,chiarr,oldchiarr,nnn)

          smallloop=1.0d33
!$omp do
          do 750 i16=1,Niter   !np

!$omp critical            
            if(isw22.ge.1)ifastflag=1
            if(i16.eq.1)ifastflag=0
            do  j=1,Nterms
              var(j)=parmarray(j,i16)
            enddo
!$omp end critical            

          chilimb=0.0d0
          ifastflag=0
          ibest=0
          ichilabel=i16
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
     @      icnarray,
     @      RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
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
!$omp critical
c            bigstring(i16,1)=line
c            bigstring(i16,2)=parmstring
c            bigstring(i16,3)=planetparm
c            bigstring(i16,4)=dynparm
            if(ifastflag.ge.1)then
              if(chi1.lt.3.0d0*small*dabs(dble(Ndattot-Nterms)))then
                ifastflag=0
              endif
            endif
c
            imodel=imodel+1
            chiarr(i16)=chi1
            oldchiarr(i16)=chi1
            if(chi1.lt.small)then
              small=chi1
              do mmm=1,Nvmax
                sss(mmm)=var(mmm)
              enddo
            endif
            if(chi1.lt.smallloop)then
              smallloop=chi1
              indexsmall=i16
              rbest=var(iterm)
c              do mmm=1,Nvmax
c                sssloop(mmm)=var(mmm)
c              enddo
            endif
!$omp end critical
            if(rmed.ge.1.0d0)then
              call printsmallmed(small)
            else
              call printsmall(small)
            endif
            nnn=0
            call printiter(i16,'step number = ',iout,
     @      'iteration number = ')
            call printiter(iterm,'variable number = ',iloop,
     @      'outer loop number = ')
            call printiter1(imodel,'model number = ')
c

 750      continue     !continue loop over generation members
!$omp enddo
!$omp end parallel

c          do jj=1,Nterms
c            var(jj)=sssloop(jj)
c          enddo
          if(indexsmall.gt.1)then
            rlow=parmarray(iterm,indexsmall-1)
            ylow=chiarr(indexsmall-1)
          else
            rlow=parmarray(iterm,1)
            ylow=chiarr(1)
          endif
          if(indexsmall.lt.Niter)then
            rup=parmarray(iterm,indexsmall+1)
            yup=chiarr(indexsmall+1)
          else
            rup=parmarray(iterm,Niter)
            yup=chiarr(Niter)
          endif
c
          xxx(1)=rlow
          xxx(2)=rbest
          xxx(3)=rup
          yyy(1)=ylow
          yyy(2)=smallloop
          yyy(3)=yup
          if(((xxx(1).lt.xxx(2)).and.(xxx(2).lt.xxx(3))).or.
     @         ((xxx(1).gt.xxx(2)).and.(xxx(2).gt.xxx(3))))then
            call quadreg(xxx,yyy,xzero)
            rbest=xzero
            if(rbest.lt.vstart(iterm))rbest=vstart(iterm)
            if(rbest.gt.vstep(iterm))rbest=vstep(iterm)
          endif
c
 800      continue
 850      continue
          do kk=1,Nterms
            var(kk)=sss(kk)
          enddo
 900      continue
c
c
c  We want to record the parameters here.
c
c          open(unit=66,file='demcmc_fitparm.1000000',status='unknown')
c          open(unit=67,file='demcmc_starparm.1000000',status='unknown')
c          open(unit=68,file='demcmc_planetparm.1000000',
c     @           status='unknown')
c          open(unit=69,file='demcmc_dynparm.1000000',status='unknown')
c          open(unit=70,file='demcmc_chi.1000000',status='unknown')
c          do jj=1,np
c            lll=lnblnk(bigstring(jj,1))
c            write(66,3344)bigstring(jj,1)(26:lll)
c            lll=lnblnk(bigstring(jj,2))
c            write(67,3344)bigstring(jj,2)(1:lll)
c            lll=lnblnk(bigstring(jj,3))
c            write(68,3344)bigstring(jj,3)(1:lll)
c            lll=lnblnk(bigstring(jj,4))
c            write(69,3344)bigstring(jj,4)(1:lll)
c            tempchi=chiarr(jj)
c            if(tempchi.gt.1.d13)tempchi=999999999999.9999d0
c            write(70,3345)tempchi
c          enddo
 3344     format(a)
 3345     format(f30.14)
c
c          close(66)
c          close(67)
c          close(68)
c          close(69) 
c          close(70)
c
c 
c   reset the variables at their best values and print the chi^2
c
          do mmm=1,Nvmax
            var(mmm)=sss(mmm)
            do nnn=1,Nterms
              oldparmarray(mmm,nnn)=parmarray(mmm,nnn)
            enddo
          enddo
c
          chilimb=0.0d0
          if(ilimbcheck.gt.1)then
            if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @           dwavey,chilimb)
            if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @           dwavey,chilimb)
          endif
c
          ifastflag=0
          ibest=99
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
     @      small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     @      icnarray,
     @      RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
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
          imodel=imodel+1
          chiall=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+
     @        chisqK+chisqRV1+chisqRV2+ochi+chilimb)

          chiall=chi1
          write(*,*)' '
          call printS(chiall)
c
          if((isw30.ge.3).and.(isw23.ge.1))then
            call writeeclipse(Ncycle,Ttimes,Tseps,isw30,Nmaxeclipse,
     @           Tdur1,Tdur2)
          endif
c         
          if(icnU.ne.430)call wlinmod(Nphase,xmod,ymodU,'modelU.linear',
     @        isw7)
          if(icnB.ne.430)call wlinmod(Nphase,xmod,ymodB,'modelB.linear',
     @        isw7)
          if(icnV.ne.430)call wlinmod(Nphase,xmod,ymodV,'modelV.linear',
     @        isw7)
          if(icnR.ne.430)call wlinmod(Nphase,xmod,ymodR,'modelR.linear',
     @        isw7)
          if(icnI.ne.430)call wlinmod(Nphase,xmod,ymodI,'modelI.linear',
     @        isw7)
          if(icnJ.ne.430)call wlinmod(Nphase,xmod,ymodJ,'modelJ.linear',
     @        isw7)
          if(icnH.ne.430)call wlinmod(Nphase,xmod,ymodH,'modelH.linear',
     @        isw7)
          if(icnK.ne.430)call wlinmod(Nphase,xmod,ymodK,'modelK.linear',
     @        isw7)
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
            if(icnH.ne.430)call wELCdata(Ndatamax,NdataH,xdataH,ydataH,
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
            call writegridout(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,
     &         omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,
     @         betarim,rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,
     @         Nref,rLx,Period,fm,separ,gamma1,t3,g3,SA3,density,sw1,
     @         sw2,sw3,T0,idraw,iecheck,iidint,iatm,ism1,icnU,icnB,icnV,
     @         icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,isw2,isw3,isw4,
     @         ilaw,wave,dbolx,dboly,dwavex,dwavey,ecc,argper,pshift,
     @         sw5,sw6,sw7,sw8,sw9,ikeep,isynch,isw5,isw6,isw7,isw8,
     @         isw9,spot1parm,spot2parm,spotdparm,primmass,primK,
     @         primrad,ratrad,frac1,frac2,ecosw,temprat,idark1,idark2,
     @         isw12,isw13,isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,
     @         sw24,powercoeff,sw25,sw26,sw27,sw28,sw29,sw30,contam,
     @         Tconj,beam1,beam2,isw25,isw26,isw27,isw28,isw29,isw30,
     @         isw31,isw32,isw33,isw34,ocose,osine,omegadot,contamS0,
     @         contamS1,
     @         contamS2,contamS3,sw47,sw48,sw49)
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
     $        Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,
     %        RV2file,Nvmax,Nvar,
     $        svar,var,vstart,stepsave,Nstep,Nobv,sobv,obv,eobv,vstep)
c
            write(command,101)1000000
            call system(command)
            write(command,102)1000000
            call system(command)
            write(command,103)1000000
            if(isw30.gt.0)call system(command)

          close(45)
          close(46)
          if(isw24.gt.0)close(47)
          if(isw30.gt.0)close(48)
          close(55)
c
c   set a flag to signal the need to update the model
c
          iupdate=0
c
c   Main generation  loop
c
c          scalegam=2.38d0/sqrt(dble(2*Nterms))
c          siggam=0.1d0
c          do 10000 ig=1,ngen
cc
cc   Main population loop
cc
c            write(extension,3333)ig+1000000
c 3333       format(i7)
cc
c            if((ig.ge.50).and.(mod(ig,50).eq.0))scalegam=1.0d0
cc            if((ig.ge.50).and.(mod(ig,100).eq.0))scalegam=0.00001d0
cc
c            open(unit=45,file='generation.'//extension,status='unknown')
c            if(isw24.ge.1)open(unit=47,file='ELCratio.'//extension,
c     @         status='unknown')
c            open(unit=55,file='chi.'//extension,status='unknown')
c            open(unit=46,file='ELCparm.'//extension,status='unknown')
c            if(isw30.ge.1)open(unit=48,file='ELCbody3parm.'//extension,
c     @        status='unknown')            
c            if(isw30.ge.3)open(unit=49,file='ELCdynparm.'//extension,
c     @        status='unknown')            
c
c
c   loop over the population
c
c            njump=0
c            do 7750 i16=1,np
c              chiarr(i16)=oldchiarr(i16)
cc
c              ifastflag=0   !disabled January 1, 2013
cc
cc   pick two random integers between 1 and np    
cc
cc
cc
c 6969         i1=int(ran1(idum)*dble(np))+1
c              if(i1.lt.1)i1=1
c              if(i1.gt.np)i1=np
c              if(i1.eq.i16)go to 6969
c 6979         i2=int(ran1(idum)*dble(np))+1
c              if(i2.lt.1)i2=1
c              if(i2.gt.np)i2=np
c              if(i2.eq.i16)go to 6979
c              if(i2.eq.i1)go to 6979
cc
c              idxarr1(i16)=i1
c              idxarr2(i16)=i2
c              zzzarr(i16)=gasdev(idum)
c              probarr(i16)=ran1(idum)
c 7750       continue
cc
cc
c!$omp parallel firstprivate(Nphase,xmod,ymodU,ymodB,ymodV,
c!$omp&      ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,
c!$omp@      RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,fracs1,
c!$omp@      fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
c!$omp@      chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,chisqRV1,
c!$omp&      chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,zeroU,resU,
c!$omp@      NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,
c!$omp@      errV,zeroV,resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
c!$omp@      xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,errJ,
c!$omp%      zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,NdataK,
c!$omp$      xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,errRV1,NRV2,
c!$omp@      xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,sobv,obv,eobv,ochi,
c!$omp&      ochidisk,ochilr,svar,var,saveasini,
c!$omp@      savxdataU,savydataU,saverrU,savxdataB,savydataB,
c!$omp@      saverrB,savxdataV,savydataV,saverrV,
c!$omp$      savxdataR,savydataR,saverrR,savxdataI,savydataI,
c!$omp&      saverrI,savxdataJ,savydataJ,saverrJ,
c!$omp&      savxdataH,savydataH,saverrH,savxdataK,savydataK,
c!$omp@      saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,savyRV2,saverrRV2,
c!$omp@      ifrac,ilum,i16,isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,
c!$omp@      isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,
c!$omp@      ibest,ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
c!$omp@      Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
c!$omp@      icnarray,Tdur1,Tdur2,
c!$omp@      RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
c!$omp@      fill1,fill2,omega1,
c!$omp@    omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
c!$omp@    rinner,router,tdisk,xi,alb1,alb2,rLx,Period,fm,separ,
c!$omp@    gamma,wave,dbolx,dboly,dwavex,dwavey,t3,g3,SA3,density,
c!$omp@    sw1,sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
c!$omp@    primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
c!$omp@    bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,
c!$omp@    sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,omegadot,
c!$omp@    contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49,gaplow,
c!$omp@    gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
c!$omp@    P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
c!$omp@    P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
c!$omp@    P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
c!$omp@    P5esin,P5incl,P5Omega,P5Q,P5ratrad,P6tconj,P6period,P6T0,
c!$omp@    P6ecos,P6esin,P6incl,P6Omega,P6Q,P6ratrad,P7tconj,P7period,
c!$omp@    P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,
c!$omp@    P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,
c!$omp@    xSC,ySC,spot1parm,spot2parm,spotdparm,Nalph1,Nbet1,Nalph2,
c!$omp@    Nbet2,Ntheta,Nradius,Nref,idraw,iecheck,iidint,iatm,ism1,
c!$omp@    ilaw,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
c!$omp@    iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,isw5,isw6,isw7,
c!$omp@    isw8,isw9,idark1,idark2,isw12,isw13,isw21,isw22,isw23,isw24,
c!$omp@    isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,
c!$omp@    NSC,compfracs,tertperiod,tertt0,tertecos,tertesin,
c!$omp@    tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72,sw73)
c!$omp@ shared(small,sss,bigstring,chiarr,nnn,njump,iupdate,
c!$omp@ idxarr1,idxarr2,parmarray,zzzarr,zzz,probarr,oldchiarr)
c!$omp@ firstprivate(Deltachi,outstring,prob,uran,idum,i1,i2,
c!$omp@ lll,lll1,lll2,lll3,siggam,scalegam,oldparmarray)
c!$omp@ private(parmstring,planetparm,dynparm,line)
c!$omp do
c
c            do 8751 i16=1,np
c!$omp critical
c              i1=idxarr1(i16)
c              i2=idxarr2(i16)
c              call istring('index',i16,outstring,lll)
c              call istring('i1',i1,outstring1,lll1)
c              call istring('i2',i2,outstring2,lll2)
c
c              write(*,*)' '
c              write(*,6990)outstring(1:lll),outstring1(1:lll1),
c     @          outstring2(1:lll2)
cc
c              call chistring('chi_index',oldchiarr(i16),
c     $              outstring,lll)
c              call chistring('chi_i1',oldchiarr(i1),
c     $              outstring1,lll2)
c              call chistring('chi_i2',oldchiarr(i2),
c     $              outstring2,lll3)
c              write(*,6990)outstring(1:lll),outstring1(1:lll2),
c     @           outstring2(1:lll3)
c 6990         format(a,', ',a,', ',a)
cc
cc              zzz=siggam*gasdev(idum)  
cc
c              zzz=siggam*zzzarr(i16)
c              call pstring('z',8,zzz,outstring,lll)
c              call pstring('gamma_0',8,scalegam,outstring1,
c     @              lll1)
c              call pstring('gamma',8,scalegam*(1.0d0+zzz),outstring2,
c     @              lll2)
c              write(*,6991,advance='no')outstring(1:lll),
c     @             outstring1(1:lll1),outstring2(1:lll2)
c 6991         format(a,', ',a,', ',a)
c
c              do 7760 j=1,Nterms
c                var(j)=oldparmarray(j,i16)      +
c     @            (oldparmarray(j,i1)-oldparmarray(j,i2))*scalegam*(1.0d0+zzz)
c                if(var(j).lt.vstart(j))var(j)=vstart(j)
c                if(var(j).gt.vstep(j))var(j)=vstep(j)
c 7760         continue
cc
c!$omp end critical
c              chilimb=0.0d0
c              if(ilimbcheck.gt.1)then
c                if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,dwavex,
c     @              dwavey,chilimb)
c                if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,dwavex,
c     @              dwavey,chilimb)
c              endif
cc
c              ibest=0
c          ichilabel=1
c          call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
c     &      ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,
c     @      RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,fracs1,
c     @      fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
c     @      chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,chisqRV1,
c     &      chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,zeroU,resU,
c     @      NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,
c     @      errV,zeroV,resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
c     @      xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,errJ,
c     %      zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,NdataK,
c     $      xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,errRV1,NRV2,
c     @      xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,sobv,obv,eobv,ochi,
c     &      ochidisk,ochilr,Nvmax,svar,var,saveasini,
c     @      savxdataU,savydataU,saverrU,savxdataB,savydataB,
c     @      saverrB,savxdataV,savydataV,saverrV,
c     $      savxdataR,savydataR,saverrR,savxdataI,savydataI,
c     &      saverrI,savxdataJ,savydataJ,saverrJ,
c     &      savxdataH,savydataH,saverrH,savxdataK,savydataK,
c     @      saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,savyRV2,saverrRV2,
c     @      ifrac,ilum,i16,isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,
c     @      isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,
c     @      ibest,ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
c     @      small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
c     $      icnarray,
c     @      RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
c     @    parmstring,planetparm,dynparm,line,fill1,fill2,omega1,
c     @    omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
c     @    rinner,router,tdisk,xi,alb1,alb2,rLx,Period,fm,separ,
c     @    gamma,wave,dbolx,dboly,dwavex,dwavey,t3,g3,SA3,density,
c     @    sw1,sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
c     @    primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
c     @    bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,
c     @    sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,omegadot,
c     @    contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49,gaplow,
c     @    gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
c     @    P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
c     @    P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
c     @    P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
c     @    P5esin,P5incl,P5Omega,P5Q,P5ratrad,P6tconj,P6period,P6T0,
c     @    P6ecos,P6esin,P6incl,P6Omega,P6Q,P6ratrad,P7tconj,P7period,
c     @    P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,
c     @    P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,
c     @    xSC,ySC,spot1parm,spot2parm,spotdparm,Nalph1,Nbet1,Nalph2,
c     @    Nbet2,Ntheta,Nradius,Nref,idraw,iecheck,iidint,iatm,ism1,
c     @    ilaw,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
c     @    iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,isw5,isw6,isw7,
c     @    isw8,isw9,idark1,idark2,isw12,isw13,isw21,isw22,isw23,isw24,
c     @    isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,isw33,isw34,
c     @    NSC,compfracs,tertperiod,tertt0,tertecos,tertesin,
c     @    tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72,sw73,
c     @    Nmaxeclipse,Tdur1,Tdur2)
cc
c!$omp critical
c
c              Deltachi=chi1-oldchiarr(i16)   !was chiarr
c              if(Deltachi.lt.0.0d0)then   ! take jump and swap
c                do j=1,Nterms
c                   parmarray(j,i16)=var(j)
c                enddo
c                chiarr(i16)=chi1
c                call pstring('Delta_chi',5,Deltachi,
c     $              outstring,lll)
c                write(*,95)outstring(1:lll)
c                njump=njump+1
c                bigstring(i16,1)=line
c                bigstring(i16,2)=parmstring
c                bigstring(i16,3)=planetparm
c                bigstring(i16,4)=dynparm
cc
c              else
c                prob=exp(-abs(Deltachi)*0.50d0)
cc
cc                uran=ran1(idum)
cc
c                uran=probarr(i16)
c                if(uran.lt.prob)then
c                  call chistring('chi_start',oldchiarr(i16),outstring,lll1)
c                  call chistring(' chi_end',chi1,outstring1,lll2)
c                  write(*,2991)outstring(1:lll1),outstring1(1:lll2)
c                  call pstring('Delta_chi',5,Deltachi,
c     @              outstring,lll)
c                  write(*,93)outstring(1:lll),uran,prob
c                  njump=njump+1
cc
c                  do j=1,Nterms
c                    parmarray(j,i16)=var(j)
c                  enddo
c                  chiarr(i16)=chi1
c                  bigstring(i16,1)=line
c                  bigstring(i16,2)=parmstring
c                  bigstring(i16,3)=planetparm
c                  bigstring(i16,4)=dynparm
c                endif
c              endif
cc
c              if(chi1.lt.small)then
c                small=chi1
c                iupdate=999
c                do mmm=1,Nvmax
c                  sss(mmm)=var(mmm)
c                enddo
c              endif
cc
c!$omp end critical
c              if(rmed.ge.1.0d0)then
c                call printsmallmed(small)
c              else
c                call printsmall(small)
c              endif
c
c              call printiter(i16,'iteration number = ',ig,
c     @        'generation number = ')
c
c
c 8751       continue             !continue loop over chain
cc
c!$omp enddo
c!$omp end parallel
c
c            write(*,*)' '
c
c            call istring('number of jumps',njump,outstring,lll)
c            call istring('Npop',np,outstring1,lll1)
c
c            write(*,7751)outstring(1:lll),outstring1(1:lll1)
c 7751       format(a,', ',a)
c
c            fracjump=dble(njump)/dble(np)
c
c            if(fracjump.lt.0.2d0)scalegam=scalegam*0.9d0
c            if(fracjump.gt.0.31d0)scalegam=scalegam*1.1d0
c            if((fracjump.ge.0.20d0).and.(fracjump.le.0.31d0))then
c              scalegam=scalegam*sqrt(fracjump/0.25d0)
c            endif
c
c            call pstring('jump fraction',8,fracjump,
c     @              outstring,lll)
c            write(*,7742)outstring(1:lll)
c
c            call pstring('new gamma',8,scalegam,
c     @              outstring,lll)
c            write(*,7742)outstring(1:lll)
c
c 7742       format(a)
c
c 93         format(a,', ran=',f7.5,', probability=',f7.5,',',
c     @         ' jump taken')
c 95         format(a,', jump taken')
c 2991       format(a,',',a)
cc
cc
c
c
cc
cc
cc  We want to record the parameters here.
cc
c          open(unit=66,file='demcmc_fitparm.'//extension,status='unknown')
c          open(unit=67,file='demcmc_starparm.'//extension,
c     @           status='unknown')
c          open(unit=68,file='demcmc_planetparm.'//extension,
c     @           status='unknown')
c          open(unit=69,file='demcmc_dynparm.'//extension,status='unknown')
c          open(unit=70,file='demcmc_chi.'//extension,status='unknown')
c          do jj=1,np
c            lll=lnblnk(bigstring(jj,1))
c            write(66,3344)bigstring(jj,1)(26:lll)
c            lll=lnblnk(bigstring(jj,2))
c            write(67,3344)bigstring(jj,2)(1:lll)
c            lll=lnblnk(bigstring(jj,3))
c            write(68,3344)bigstring(jj,3)(1:lll)
c            lll=lnblnk(bigstring(jj,4))
c            write(69,3344)bigstring(jj,4)(1:lll)
c            tempchi=chiarr(jj)
c            if(tempchi.gt.1.d13)tempchi=999999999999.9999d0
c            write(70,3345)tempchi
c          enddo
cc
cc
cc   reset the variables at their best values and print the chi^2 if needed
cc
          write(*,*)'iupdate = ' ,iupdate
          if(iupdate.ne.0)then
            do mmm=1,Nvmax
              var(mmm)=sss(mmm)
              do nnn=1,np
                oldparmarray(mmm,nnn)=parmarray(mmm,nnn)
              enddo
            enddo
c
            chilimb=0.0d0
            if(ilimbcheck.gt.1)then
              if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @           dwavey,chilimb)
              if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @           dwavey,chilimb)
            endif
c
            ibest=99
            ichilabel=0
            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     &       ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,
     @       RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,fracs1,
     @       fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,chisqU,
     @       chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,chisqRV1,
     &       chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,zeroU,resU,
     @       NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,xdataV,ydataV,
     @       errV,zeroV,resV,NdataR,xdataR,ydataR,errR,zeroR,resR,NdataI,
     @       xdataI,ydataI,errI,zeroI,resI,NdataJ,xdataJ,ydataJ,errJ,
     %       zeroJ,resJ,NdataH,xdataH,ydataH,errH,zeroH,resH,NdataK,
     $       xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,yRV1,errRV1,NRV2,
     @       xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,sobv,obv,eobv,ochi,
     &       ochidisk,ochilr,Nvmax,svar,var,saveasini,
     @       savxdataU,savydataU,saverrU,savxdataB,savydataB,
     @       saverrB,savxdataV,savydataV,saverrV,
     $       savxdataR,savydataR,saverrR,savxdataI,savydataI,
     &       saverrI,savxdataJ,savydataJ,saverrJ,
     &       savxdataH,savydataH,saverrH,savxdataK,savydataK,
     @       saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,
     @       savyRV2,saverrRV2,
     @       ifrac,ilum,i16,isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,
     @       isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,Ndatamax,
     @       ibest,ifixgamma,savesep,ichilabel,resRV1,resRV2,thresh,
     @       small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,obsTerr,
     @       icnarray,
     @       RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,
     @       parmstring,planetparm,dynparm,line,fill1,fill2,omega1,
     @       omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,
     @       rinner,router,tdisk,xi,alb1,alb2,rLx,Period,fm,separ,
     @       gamma,wave,dbolx,dboly,dwavex,dwavey,t3,g3,SA3,density,
     @       sw1,sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     @       primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @       bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,
     @       sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,omegadot,
     @       contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49,gaplow,
     @       gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,
     @       P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @       P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,
     @       P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @       P5esin,P5incl,P5Omega,P5Q,P5ratrad,P6tconj,P6period,P6T0,
     @       P6ecos,P6esin,P6incl,P6Omega,P6Q,P6ratrad,P7tconj,P7period,
     @       P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,
     @       P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,
     @       xSC,ySC,spot1parm,spot2parm,spotdparm,Nalph1,Nbet1,Nalph2,
     @       Nbet2,Ntheta,Nradius,Nref,idraw,iecheck,iidint,iatm,ism1,
     @       ilaw,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
     @       iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,isw5,isw6,isw7,
     @       isw8,isw9,idark1,idark2,isw12,isw13,isw21,
     @       isw22,isw23,isw24,
     @       isw25,isw26,isw27,isw28,isw29,isw30,isw31,
     @       isw32,isw33,isw34,
     @       NSC,compfracs,tertperiod,tertt0,tertecos,tertesin,
     @       tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,sw72,sw73,
     @       Nmaxeclipse,Tdur1,Tdur2)
c
              chiall=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+
     @          chisqK+chisqRV1+chisqRV2+ochi+chilimb)
              write(*,*)' '
              chiall=chi1
              call printS(chiall)
c
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
c
              if(icnU.ne.430)call wmagmod(Nphase,xmod,ymodU,
     @          'modelU.mag',isw7,zeroU)
              if(icnB.ne.430)call wmagmod(Nphase,xmod,ymodB,
     @          'modelB.mag',isw7,zeroB)
              if(icnV.ne.430)call wmagmod(Nphase,xmod,ymodV,
     @          'modelV.mag',isw7,zeroV)
              if(icnR.ne.430)call wmagmod(Nphase,xmod,ymodR,
     @          'modelR.mag',isw7,zeroR)
              if(icnI.ne.430)call wmagmod(Nphase,xmod,ymodI,
     @          'modelI.mag',isw7,zeroI)
              if(icnJ.ne.430)call wmagmod(Nphase,xmod,ymodJ,
     @          'modelJ.mag',isw7,zeroJ)
              if(icnH.ne.430)call wmagmod(Nphase,xmod,ymodH,
     @          'modelH.mag',isw7,zeroH)
              if(icnK.ne.430)call wmagmod(Nphase,xmod,ymodK,
     @          'modelK.mag',isw7,zeroK)
  
              call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
              call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
              call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
              if(iidint.ge.1)call wlinmod(Nphase,xmod,ymodd,
     @          'lcdisk.linear',isw7)
c
              call wRVmod(NRVphase,xRVmod,RV1,'star1.RV',isw7,ggamma1)
              call wRVmod(NRVphase,xRVmod,RV2,'star2.RV',isw7,ggamma2)
              if(isw30.ge.3)call wRVmod(NRVphase,xRVmod,RV3,
     @         'star3.RV',isw7,ggamma1)
c
              if(itime.gt.0)then
               if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,
     @             ydataU,errU,'ELCdataU.fold')
               if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,
     @             ydataB,errB,'ELCdataB.fold')
               if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,
     @             ydataV,errV,'ELCdataV.fold')
               if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,
     @             ydataR,errR,'ELCdataR.fold')
               if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,
     @             ydataI,errI,'ELCdataI.fold')
               if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,
     @             ydataJ,errJ,'ELCdataJ.fold')
               if(icnH.ne.430)call wELCdata(Ndatamax,NdataH,xdataH,
     @             ydataH,errH,'ELCdataH.fold')
               if(icnK.ne.430)call wELCdata(Ndatamax,NdataK,xdataK,
     @             ydataK,errK,'ELCdataK.fold')
               if(icnRV1.ne.430) call wELCdata(Ndatamax,NRV1,xRV1,yRV1,
     %              errRV1,'ELCdataRV1.fold')
               if(icnRV2.ne.430) call wELCdata(Ndatamax,NRV2,xRV2,yRV2,
     %              errRV2,'ELCdataRV2.fold')
               if(icnRV3.ne.430) call wELCdata(Ndatamax,NRV3,xRV3,yRV3,
     %              errRV3,'ELCdataRV3.fold')
              endif
c
              if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,resU,
     %             errU,'ELCresidualsU.fold')
              if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,resB,
     %            errB,'ELCresidualsB.fold')
              if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,resV,
     %            errV,'ELCresidualsV.fold')
              if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,resR,
     %            errR,'ELCresidualsR.fold')
              if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,resI,
     %            errI,'ELCresidualsI.fold')
              if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,resJ,
     %             errJ,'ELCresidualsJ.fold')
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
              call writegridout(Nalph1,Nbet1,Nalph2,Nbet2,fill1,fill2,
     &         omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,
     @         betarim,rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,
     @         Nref,rLx,Period,fm,separ,gamma1,t3,g3,SA3,density,sw1,
     @         sw2,sw3,T0,idraw,iecheck,iidint,iatm,ism1,icnU,icnB,icnV,
     @         icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,isw2,isw3,isw4,
     @         ilaw,wave,dbolx,dboly,dwavex,dwavey,ecc,argper,pshift,
     @         sw5,sw6,sw7,sw8,sw9,ikeep,isynch,isw5,isw6,isw7,isw8,
     @         isw9,spot1parm,spot2parm,spotdparm,primmass,primK,
     @         primrad,ratrad,frac1,frac2,ecosw,temprat,idark1,idark2,
     @         isw12,isw13,isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,
     @         sw24,powercoeff,sw25,sw26,sw27,sw28,sw29,sw30,contam,
     @         Tconj,beam1,beam2,isw25,isw26,isw27,isw28,isw29,isw30,
     @         isw31,isw32,isw33,isw34,ocose,osine,omegadot,contamS0,
     @         contamS1,
     @         contamS2,contamS3,sw47,sw48,sw49)
c
              if(isw30.gt.0)call writebody3grid(Nalph3,Nbet3,tertperiod,
     @         tertt0,tertecos,tertesin,tertincl,tertOmega,tertQ,dwavex,
     @         dwavey,itconj,it1,it2,it3,it4,tertconj,
     @         tertratrad,hh,sw72,sw73,
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
     $         Idatafile,Jdatafile,Hdatafile,Kdatafile,RV1file,
     %         RV2file,Nvmax,Nvar,
     $         svar,var,vstart,stepsave,Nstep,Nobv,sobv,obv,eobv,vstep)
c
              write(command,101)ig+1000000
              call system(command)
              write(command,102)ig+1000000
              call system(command)
              write(command,103)ig+1000000
              if(isw30.gt.0)call system(command)
c
              iupdate=0
            endif
            close(45)
            close(46)
            close(66)
            close(67)
            close(68)
            close(69)
            close(70)
            if(isw24.gt.0)close(47)
            if(isw30.gt.0)close(48)
            if(isw30.ge.3)close(49)
            close(55)
c
c   copy the current chi^2 array into old one
c
            do jkl=1,Nterms
              oldchiarr(jkl)=chiarr(jkl)
            enddo

10000     continue   !continue over generations
c
c
          close(38)
c
 101      format('cp gridELC.opt gridELC.',i7)
 102      format('cp gridELC.inp ELC.',i7)
 103      format('cp gridELCbody3.inp ELCbody3.',i7)
c
          stop
          end
c
c     &&&&&&&&&&&&&&&&&&&&&&&&
c
          subroutine indexx(n,arr,indx)
c
          implicit none
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
          implicit none
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
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
c&&&&&&&&&&&&&&&&&&&&&!!!!!@@@@@@@@^^^^^&&&&&&&  
c
c
          subroutine quadreg(x,y,xzero)
c
          implicit double precision(a-h,o-z)
c
          dimension x(3),y(3)
c
          x1=x(1)
          x2=x(2)
          x3=x(3)
          y1=y(1)
          y2=y(2)
          y3=y(3)
c
          denom=(x1-x2)*(x1-x3)*(x2-x3)
          a=(x3*(-y1+y2)+x2*(y1-y3)+x1*(-y2+y3))/denom
          b=(x3**2*(y1-y2)+x1**2*(y2-y3)+x2**2*(-y1+y3))/denom
          c=(x3*(x2*(x2-x3)*y1+x1*(-x1+x3)*y2)+x1*(x1-x2)*x2*y3)/denom
          xzero=-b/(2.0d0*a)
c
          return
          end
c
c &&&&&&&&&&&&&%%%%%%%%%%%%%%$$$$$$$$$$$$$
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
          include 'dynamicssubs.for'
