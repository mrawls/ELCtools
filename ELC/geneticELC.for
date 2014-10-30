          program geneticELC
c
c   May 11, 2001
c
c   This program will read in the light curves and parameters specified
c   in the 'gridloop.opt' file and optimize the fits using a genetic
c   algorithm (P. Charbonneau, 1995, ApJS, 101, 309).
c
c   
c   July 10, 2001
c
c   Program modified to include "black sheep" (A. Bobinger, 2000,
c   A&A 357, 1170).  After each population replacement and ranking,
c   replace the worst members with copies of the best member where
c   the genes are varied by small random amounts.  Specifically, the
c   number of black sheep = int(Npop*[1-pcross]/2.0), and 
c   parray(j,i)=parray(j,i)+0.05*gasdev(idum)
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
c   UPDATE August 13, 2001
c
c   Change the number of maximum variables to 25 (was 19).  Change
c   the 19 to 25 in the dimension of parmarray,dummy,fakevar and svar
c
          parameter(Nmaxphase=760000,Ndatamax=200000)
          parameter(Nvmax=60)
          parameter(Nmaxeclipse=1000)
c
c   UPDATE September 11, 2001
c
c   Change the dimension of obsparm,obv,eobv,sobv to 9.
c
c
c   UPDATE March 19, 2002
c
c   Add the parameter 'ifixgamma' to the end of the argument list
c   of checklcfit.  It is missing from some calls, but not others.
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
          dimension dummy(Nvmax)
          dimension xRVmod(Nmaxphase)
          dimension resU(Ndatamax),resB(Ndatamax),resV(Ndatamax)
          dimension resR(Ndatamax),resI(Ndatamax),resJ(Ndatamax)
          dimension resH(Ndatamax),resK(Ndatamax),resRV1(Ndatamax)
          dimension resRV2(Ndatamax),resRV3(Ndatamax)
c
c   UPDATE August 13, 2001
c
c   Make the variable 'line' 200 characters
c
c   UPDATE October 28, 2002
c
c   Make the length of line 300 (was 200)
c
c
c          character*900 line        ! was 132
          character*7 extension
c
c   Dimension the variables needed for the model atmospheres here.
c
          parameter (maxlines=1300,maxmu=115)   !was 1100
          dimension atmT(maxlines),atmg(maxlines),atmmu(maxlines,maxmu),
     %       Nmu(maxlines)
          dimension atmint1(maxlines,maxmu),atmint2(maxlines,maxmu)
          dimension atmint3(maxlines,maxmu),atmint4(maxlines,maxmu)
          dimension atmint5(maxlines,maxmu),atmint6(maxlines,maxmu)
          dimension atmint7(maxlines,maxmu),atmint8(maxlines,maxmu)
c
          dimension spotdparm(2,4),spot1parm(2,4),spot2parm(2,4)
c
c   Here are some variables needed for the genetic part.
c
          parameter(nd=6,ipmax=8192,nmax=64,idmax=8,
     %       pmutmn=0.0005d0,pmutmx=0.25d0)

          dimension ign1(nmax*idmax),ign2(nmax*idmax),ifit(ipmax),
     @        jfit(ipmax)
          dimension fitns(ipmax),ph(nmax,2),oldph(nmax,ipmax),
     &        rnewph(nmax,ipmax)
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
c   Add a new character string and common block for a 'parameter string'
c   This string of parameters will be fed to the genetic code to make it
c   easier to compute uncertainties on the physical quantities like mass
c   and radius.
c
c
c   UPDATE November 28, 2001
c
c   Change parmstring from character*199 to character*201
c
c   UPDATE January 16, 2002
c
c   parmstring was character*201, now should be character*227
c
c   UPDATE June 7, 2002
c
c   Make the length of parmstring character*237.
c
c   UPDATE October 28, 2002
c
c   Make the length of parmstring character*249.
c
c   UPDATE October 22, 20028
c
c   Make the length of parmstring character*259.
c

          character*1000 parmstring
          character*2000 planetparm      
          character*1700 dynparm,line
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
cc
          common /realatm/ atmT,atmg,atmmu,atmint1,atmint2,atmint3,
     &      atmint4,atmint5,atmint6,atmint7,atmint8,Tmax,Tmin,gmax,gmin
c
          common /intatm/  Nlines,Nmu,Nalph3,Nbet3,itconj,it1,it2,it3,
     @      it4
c
          common /ranblock/ idum
c
          common /medblock/ rmed 
c
c
c   Open the parameter file and read all of the parameters. 
c   Pass the parameters to the light curve routines
c   via a common block.  
c
c   Add the spot parameters to getinput and recordparm.
c
c   UPDATE August 10, 2004
c
c   Add the 8 real and 4 integer variables to the list.
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
c
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
c   Here is a new block to allow for the inclusion of the disk fraction
c   in the total chi^2
c
c   UPDATE FEBRUARY 5, 2005
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
4323      format('error:  initial value of variable ',a2, 
     $            ' is out of range')
4324      format(3(f16.8,3x))

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
          Ndattot=0
          Ndattot=NdataU+NdataB+NdataV+NdataR+NdataI+NdataJ+NdataH+
     $        NdataK+NRV1+NRV2+NRV3
c
c   Open an output file for the fitting statistics.
c
          open(unit=38,file='geneticELC.out',status='unknown')
          open(unit=45,file='generation.1000000',status='unknown')
          if(isw24.ge.1)open(unit=47,file='ELCratio.1000000',
     $         status='unknown')
          open(unit=55,file='chi.1000000',status='unknown')
          open(unit=46,file='ELCparm.1000000',status='unknown')
          if(isw30.ge.1)open(unit=48,file='ELCbody3parm.1000000',
     @      status='unknown')
          if(isw30.ge.3)open(unit=49,file='ELCdynparm.1000000',
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
          if(isw9save.gt.Nstep(1))then
            write(*,*)'ielete is too large'
            stop
          endif
c
          write(*,*)'Nterms, Nvar, Ndattot',Nterms,Nvar,Ndattot
c
c   Here are some variables needed for the genetic part.
c
          pcross=0.85d0 
          pmut=0.005d0
          fdif=1.0d0
c   
c   Define the random parameter sets.  Nstep(1) is the number of random
c   sets to define.
c

          np=Nstep(1)     ! number of members in the population 
          ngen=Nstep(2)   ! number of generations

          if(np.gt.8192)then
            write(*,*)'Error:  Maximum population allowed is 8192'
            stop
          endif
          do 99 j=1,Nterms
            vx1=vstart(j)
            vx2=vstep(j)
            vxmult=(vx2-vx1)
            vxadd=vx1
            do 98 i=1,np                 
              urand=ran1(idum)
              varx=urand*vxmult+vxadd
              parmarray(j,i)=varx
              oldph(j,i)=urand
 98         continue
 99       continue
c
c   Evaluate the fitness of the initial population.
c
          do 7750 i16=1,np
            if(isw22.ge.1)ifastflag=1
            if(i16.eq.1)ifastflag=0
            do  j=1,Nterms
              var(j)=parmarray(j,i16)
            enddo
c
            if((isw9save.ge.1).and.(i16.eq.1))then
c
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
     @          P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     @          P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,
     @          P4esin,P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,
     @          P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,P5ratrad,P6tconj,
     $          P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,P6ratrad,
     @          P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,
     @          P7ratrad,P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,
     @          P8Omega,P8Q,P8ratrad,sw72,sw73)
c
              do j=1,Nterms
                vx1=vstart(j)
                vx2=vstep(j)
                vxmult=(vx2-vx1)
                vxadd=vx1
                varx=var(j)
                urand=(varx-vxadd)/vxmult
                if(urand.gt.1.0d0)urand=1.0d0
                if(urand.lt.0.0d0)urand=0.0d0
                parmarray(j,1)=varx
                oldph(j,1)=urand
              enddo
c
            endif
c
c   UPDATE December 20, 2012
c
c   if ielete (isw9) is more than 1, add additional preset
c   parameters to the mix
c
            if((isw9save.gt.1).and.(i16.gt.1).and.(i16.le.isw9))then
c
              kkk=i16
c
              call getgridinput(kkk-1,iNalph1,iNbet1,iNalph2,iNbet2,
     @          fill1,fill2,omega1,omega2,rdphase,Q,finc,Teff1,Teff2,
     @          Tgrav1,Tgrav2,betarim,rinner,router,tdisk,xi,iNtheta,
     @          iNradius,alb1,alb2,iNref,rLx,Period,fm,separ,gamma,t3,
     @          g3,SA3,density,sw1,sw2,sw3,T0,iidraw,iiecheck,iiidint,
     @          iiatm,iism1,iicnU,iicnB,iicnV,iicnR,iicnI,iicnJ,iicnH,
     @          iicnK,iiRVfilt,iisw1,iisw2,iisw3,iisw4,iilaw,wave,dbolx,
     @          dboly,dwavex,dwavey,ecc,argper,pshift,sw5,rsw6,rsw7,
     @          rsw8,rsw9,iikeep,iisynch,iisw5,iisw6,iisw7,iisw8,iisw9,
     @          spot1parm,spot2parm,spotdparm,primmass,primK,primrad,
     @          ratrad,frac1,frac2,ecosw,temprat,iidark1,iidark2,iisw12,
     @          iisw13,iisw21,iisw22,iisw23,iisw24,bigI,bigbeta,rsw23,
     @          rsw24,powercoeff,sw25,rsw26,rsw27,rsw28,rsw29,rsw30,
     @          contam,Tconj,beam1,beam2,iisw25,iisw26,iisw27,iisw28,
     @          iisw29,iisw30,iisw31,iisw32,iisw33,iisw34,ocose,osine,
     @          omegadot,contamS0,contamS1,contamS2,contamS3,sw47,sw48,
     @          sw49)
c
              if((isw30.ge.1).and.(isw7.ge.2))then
                  call getgridbody3(kkk-1,iNalph3,iNbet3,tertperiod,
     @              tertt0,tertecos,tertesin,tertincl,tertOmega,tertQ,
     @              dwavex,dwavey,itconj,it1,it2,it3,it4,tertconj,
     @              tertratrad,hh,sw72,sw73,P2tconj,P2period,P2T0,
     @              P2ecos,P2esin,P2incl,P2Omega,P2Q,P2ratrad,P3tconj,
     @              P3period,P3T0,P3ecos,P3esin,P3incl,P3Omega,P3Q,
     @              P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,P4incl,
     @              P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @              P5esin,P5incl,P5Omega,P5Q,P5ratrad,P6tconj,P6period,
     @              P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,P6ratrad,
     @              P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,
     @              P7Q,P7ratrad,P8tconj,P8period,P8T0,P8ecos,P8esin,
     @              P8incl,P8Omega,P8Q,P8ratrad)
c                             
              endif
c
              call varassign(Nvmax,svar,var,fill1,fill2,omega1,omega2,
     @            Q,finc,Teff1,Teff2,betarim,rinner,router,tdisk,xi,rLx,
     @            separ,gamma,t3,g3,sa3,ecc,argper,pshift,spot1parm,
     @            spot2parm,spotdparm,period,T0,alb1,alb2,dwavex,dwavey,
     @            primmass,primK,primrad,ratrad,frac1,frac2,ecosw,
     @            temprat,bigI,bigbeta,powercoeff,density,Tconj,beam1,
     @            beam2,contam,ocose,osine,isw29,tertperiod,tertt0,
     @            tertecos,tertesin,tertincl,tertOmega,tertQ,Tgrav1,
     @            Tgrav2,tertconj,omegadot,contamS0,contamS1,contamS2,
     @            contamS3,P2tconj,P2period,P2T0,P2ecos,P2esin,P2incl,
     @            P2Omega,P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,
     @            P3esin,P3incl,P3Omega,P3Q,P3ratrad,P4tconj,P4period,
     @            P4T0,P4ecos,P4esin,P4incl,P4Omega,P4Q,P4ratrad,
     @            P5tconj,P5period,P5T0,P5ecos,P5esin,P5incl,P5Omega,
     @            P5Q,P5ratrad,P6tconj,P6period,P6T0,P6ecos,P6esin,
     @            P6incl,P6Omega,P6Q,P6ratrad,P7tconj,P7period,P7T0,
     @            P7ecos,P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,
     @            P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @            P8ratrad,sw72,sw73)
c
              do j=1,Nterms
                vx1=vstart(j)
                vx2=vstep(j)
                vxmult=(vx2-vx1)
                vxadd=vx1
                varx=var(j)
                urand=(varx-vxadd)/vxmult
                if(urand.gt.1.0d0)urand=1.0d0
                if(urand.lt.0.0d0)urand=0.0d0
                parmarray(j,kkk)=varx
                oldph(j,kkk)=urand
              enddo
c
            endif
 7750     continue
c
          iupdate=0
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
!$omp@  shared(small,sss,bigstring,chiarr,oldchiarr,nnn,fitns)
!$omp@  private(parmstring,planetparm,dynparm,line)
!$omp@  shared(/realatm/, /intatm/, /medblock/)
!$omp do

          do 750 i16=1,np
c
!$omp critical
            do j=1,Nterms
              var(j)=parmarray(j,i16)
            enddo
!$omp end critical
            chilimb=0.0d0
            if(ilimbcheck.gt.1)then
              if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,
     @             dwavex,dwavey,chilimb)
              if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,
     @             dwavex,dwavey,chilimb)
            endif
c
            ifastflag=0
            ibest=0
            ichilabel=1
            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     &        ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,
     @        RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     @        fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,
     @        chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,
     @        chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,
     @        zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,
     @        xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,ydataR,errR,
     @        zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,resI,NdataJ,
     @        xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,
     @        zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,
     @        yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,
     @        sobv,obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,
     @        saveasini,savxdataU,savydataU,saverrU,savxdataB,savydataB,
     @        saverrB,savxdataV,savydataV,saverrV,savxdataR,savydataR,
     @        saverrR,savxdataI,savydataI,saverrI,savxdataJ,savydataJ,
     @        saverrJ,savxdataH,savydataH,saverrH,savxdataK,savydataK,
     @        saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,savyRV2,
     @        saverrRV2,ifrac,ilum,i16,isavNU,isavNB,isavNV,isavNR,
     @        isavNI,isavNJ,isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,
     @        Ndatamax,ibest,ifixgamma,savesep,ichilabel,resRV1,resRV2,
     @        thresh,small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,
     @        obsTerr,icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,
     @        ggamma3,NRV3,parmstring,planetparm,dynparm,line,fill1,
     @        fill2,omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,
     @        Tgrav2,betarim,rinner,router,tdisk,xi,alb1,alb2,rLx,
     @        Period,fm,separ,gamma,wave,dbolx,dboly,dwavex,dwavey,t3,
     @        g3,SA3,density,sw1,sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,
     @        sw7,sw8,sw9,primmass,primK,primrad,ratrad,frac1,frac2,
     @        ecosw,temprat,bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,
     @        sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,
     @        omegadot,contamS0,contamS1,contamS2,contamS3,sw47,sw48,
     @        sw49,gaplow,gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,
     @        P2incl,P2Omega,P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,
     @        P3esin,P3incl,P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,
     @        P4ecos,P4esin,P4incl,P4Omega,P4Q,P4ratrad,P5tconj,
     @        P5period,P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @        P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @        P6ratrad,P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,
     @        P7Omega,P7Q,P7ratrad,P8tconj,P8period,P8T0,P8ecos,P8esin,
     @        P8incl,P8Omega,P8Q,P8ratrad,xSC,ySC,spot1parm,spot2parm,
     @        spotdparm,Nalph1,Nbet1,Nalph2,Nbet2,Ntheta,Nradius,Nref,
     @        idraw,iecheck,iidint,iatm,ism1,ilaw,icnU,icnB,icnV,icnR,
     @        icnI,icnJ,icnH,icnK,icnRV1,icnRV2,iRVfilt,isw1,isw2,isw3,
     @        isw4,ikeep,isynch,isw5,isw6,isw7,isw8,isw9,idark1,idark2,
     @        isw12,isw13,isw21,isw22,isw23,isw24,isw25,isw26,isw27,
     @        isw28,isw29,isw30,isw31,isw32,isw33,isw34,NSC,compfracs,
     @        tertperiod,tertt0,tertecos,tertesin,tertincl,tertOmega,
     @        tertQ,tertconj,tertratrad,hh,sw72,sw73,Nmaxeclipse,Tdur1,
     @        Tdur2)
c
            if(ifastflag.ge.1)then
              if(chi1.lt.3.0d0*small*dabs(dble(Ndattot-Nterms)))then
                ifastflag=0
              endif
            endif
c
!$omp critical
            chiarr(i16)=chi1
            fitns(i16)=1.0d0/chi1    ! the fitness is  prop.to. chi^2
            nnn=0
            call printiter(i16,'iteration number = ',nnn,
     @      'generation number = ')
c
            if(chi1.lt.small)then
              small=chi1
              do mmm=1,Nvmax
                sss(mmm)=var(mmm)
                iupdate=99
              enddo
            endif
            if(rmed.ge.1.0d0)then
              call printsmallmed(small)
            else
              call printsmall(small)
            endif
!$omp end critical
 750      continue     !continue loop over generation members
!$omp enddo
!$omp end parallel
c 
c   Sort the chiarr and print out parameters for the best Nstep(2) sets.
c
          call rnkpop(np,fitns,ifit,jfit)
          call indexx(Np,chiarr,indxchi)
          write(38,8850)nnn,chiarr(indxchi(1)),chiarr(indxchi(np/2)),
     @            pmut
c
c   reset the variables at their best values and print the chi^2
c
          if(iupdate.gt.0)then
            do mmm=1,Nvmax
              var(mmm)=sss(mmm)
            enddo
c
            chilimb=0.0d0
            if(ilimbcheck.gt.1)then
              if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @           dwavey,chilimb)
              if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @           dwavey,chilimb)
            endif
            ifastflag=0
            ibest=99
            ichilabel=1
            call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     &        ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,ymodd,
     @        RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,xRVmod,
     @        fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,fracs8,
     @        chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,chisqK,
     @        chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,ydataU,errU,
     @        zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,resB,NdataV,
     @        xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,ydataR,errR,
     @        zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,resI,NdataJ,
     @        xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,ydataH,errH,
     @        zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,resK,NRV1,xRV1,
     @        yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,ggamma1,ggamma2,Nobv,
     @        sobv,obv,eobv,ochi,ochidisk,ochilr,Nvmax,svar,var,
     @        saveasini,savxdataU,savydataU,saverrU,savxdataB,savydataB,
     @        saverrB,savxdataV,savydataV,saverrV,savxdataR,savydataR,
     @        saverrR,savxdataI,savydataI,saverrI,savxdataJ,savydataJ,
     @        saverrJ,savxdataH,savydataH,saverrH,savxdataK,savydataK,
     @        saverrK,savxRV1,savyRV1,saverrRV1,savxRV2,savyRV2,
     @        saverrRV2,ifrac,ilum,i16,isavNU,isavNB,isavNV,isavNR,
     @        isavNI,isavNJ,isavNH,isavNK,isavRV1,isavRV2,isvel1,isvel2,
     @        Ndatamax,ibest,ifixgamma,savesep,ichilabel,resRV1,resRV2,
     @        thresh,small,Ncycle,Ttimes,Tseps,Nobscycle,obsTtimes,
     &        obsTerr,icnarray,RV3,xRV3,yRV3,errRV3,icnRV3,resRV3,
     @        ggamma3,NRV3,parmstring,planetparm,dynparm,line,fill1,
     &        fill2,omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,
     @        Tgrav2,betarim,rinner,router,tdisk,xi,alb1,alb2,rLx,
     @        Period,fm,separ,gamma,wave,dbolx,dboly,dwavex,dwavey,t3,
     @        g3,SA3,density,sw1,sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,
     @        sw7,sw8,sw9,primmass,primK,primrad,ratrad,frac1,frac2,
     @        ecosw,temprat,bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,
     @        sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,
     @        omegadot,contamS0,contamS1,contamS2,contamS3,sw47,sw48,
     @        sw49,gaplow,gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,
     @        P2incl,P2Omega,P2Q,P2ratrad,P3tconj,P3period,P3T0,P3ecos,
     @        P3esin,P3incl,P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,
     @        P4ecos,P4esin,P4incl,P4Omega,P4Q,P4ratrad,P5tconj,
     @        P5period,P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @        P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @        P6ratrad,P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,
     @        P7Omega,P7Q,P7ratrad,P8tconj,P8period,P8T0,P8ecos,P8esin,
     @        P8incl,P8Omega,P8Q,P8ratrad,xSC,ySC,spot1parm,spot2parm,
     @        spotdparm,Nalph1,Nbet1,Nalph2,Nbet2,Ntheta,Nradius,Nref,
     @        idraw,iecheck,iidint,iatm,ism1,ilaw,icnU,icnB,icnV,icnR,
     @        icnI,icnJ,icnH,icnK,icnRV1,icnRV2,iRVfilt,isw1,isw2,isw3,
     @        isw4,ikeep,isynch,isw5,isw6,isw7,isw8,isw9,idark1,idark2,
     @        isw12,isw13,isw21,isw22,isw23,isw24,isw25,isw26,isw27,
     @        isw28,isw29,isw30,isw31,isw32,isw33,isw34,NSC,compfracs,
     @        tertperiod,tertt0,tertecos,tertesin,tertincl,tertOmega,
     @        tertQ,tertconj,tertratrad,hh,sw72,sw73,Nmaxeclipse,Tdur1,
     @        Tdur2)
c
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
c
            call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
            call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
            call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
            if(iidint.ge.1)call wlinmod(Nphase,xmod,ymodd,
     @        'lcdisk.linear', isw7)
c
            call wRVmod(NRVphase,xRVmod,RV1,'star1.RV',isw7,ggamma1)
            call wRVmod(NRVphase,xRVmod,RV2,'star2.RV',isw7,ggamma2)
            if(isw30.ge.3)call wRVmod(NRVphase,xRVmod,RV3,'star3.RV',
     @          isw7,ggamma1)
c
            if(itime.gt.0)then
              if(icnU.ne.430)call wELCdata(Ndatamax,NdataU,xdataU,
     %           ydataU,errU,'ELCdataU.fold')
              if(icnB.ne.430)call wELCdata(Ndatamax,NdataB,xdataB,
     &           ydataB,errB,'ELCdataB.fold')
              if(icnV.ne.430)call wELCdata(Ndatamax,NdataV,xdataV,
     @           ydataV,errV,'ELCdataV.fold')
              if(icnR.ne.430)call wELCdata(Ndatamax,NdataR,xdataR,
     @           ydataR,errR,'ELCdataR.fold')
              if(icnI.ne.430)call wELCdata(Ndatamax,NdataI,xdataI,
     @           ydataI,errI,'ELCdataI.fold')
              if(icnJ.ne.430)call wELCdata(Ndatamax,NdataJ,xdataJ,
     @           ydataJ,errJ,'ELCdataJ.fold')
              if(icnH.ne.430)call wELCdata(Ndatamax,NdataH,xdataH,
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
            iupdate=0
c
          endif   !if iupdate > 0
c
          close(45)
          close(46)
          if(isw24.gt.0)close(47)
          if(isw30.gt.0)close(48)
          close(55)
c
c   Main generation  loop
c
          do 10000 ig=1,ngen
c
c   Main population loop
c
            write(extension,3333)ig+1000000
 3333       format(i7)

            open(unit=45,file='generation.'//extension,status='unknown')
            if(isw24.ge.1)open(unit=47,file='ELCratio.'//extension,
     @         status='unknown')
            open(unit=55,file='chi.'//extension,status='unknown')
            open(unit=46,file='ELCparm.'//extension,status='unknown')
            if(isw30.ge.1)open(unit=48,file='ELCbody3parm.'//extension,
     @        status='unknown')            
            if(isw30.ge.3)open(unit=49,file='ELCdynparm.'//extension,
     @        status='unknown')            
c
            do 20000 ip=1,np/2
c
              n=Nterms
c
c   1.  Pick two parents.
c
              call select(np,jfit,fdif,ip1)
21000         call select(np,jfit,fdif,ip2)
              if(ip1.eq.ip2)go to 21000
c
c   2.  Encode parent phenotypes.
c
              do 19000 kkk=1,Nterms
                dummy(kkk)=oldph(kkk,ip1)
19000         continue
              call encode(n,nd,dummy,ign1)

              do 19001 kkk=1,Nterms
                dummy(kkk)=oldph(kkk,ip2)
19001         continue
              call encode(n,nd,dummy,ign2)
c
c   3.  Breed.
c
              call cross(n,nd,pcross,ign1,ign2)
              call mutate(n,nd,pmut,ign1)
              call mutate(n,nd,pmut,ign2)
c
c   4.  Decode offspring genotypes.
c
              call decode(n,nd,ign1,dummy)
              do 19002 kkk=1,Nterms
                ph(kkk,1)=dummy(kkk)
19002         continue

              call decode(n,nd,ign2,dummy)
              do 19003 kkk=1,Nterms
                ph(kkk,2)=dummy(kkk)
19003         continue

c             
c   5.  Insert into population.
c
              call genrep(nmax,n,np,ip,ph,rnewph)
c
20000       continue
c
c   We now must replace the old population by the new, and evaluate
c   and rank the fitness.  First, put the fittest member of the old
c   population into the new one.
c
            do 44000 k=1,Nterms
              rnewph(k,1)=oldph(k,ifit(np))
c
c   update to keep the top two best models
c
              rnewph(k,2)=oldph(k,ifit(np-1))
44000       continue
c
            do 22000 i=1,np
              do 33000 k=1,Nterms
                oldph(k,i)=rnewph(k,i)
c                parmarray(k,i)=rnewph(k,i)
33000         continue
22000       continue
c
            do 999 j=1,Nterms
              vx1=vstart(j)
              vx2=vstep(j)
              vxmult=(vx2-vx1)
              vxadd=vx1
              do 998 i=1,np    !!was 1                 
                urand=oldph(j,i)
                varx=urand*vxmult+vxadd
                parmarray(j,i)=varx
 998          continue
 999        continue
c
c   Get the new fitness.   make parallel here
c
!$omp parallel firstprivate(Nphase,xmod,ymodU,ymodB,ymodV,
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
!$omp@ shared(small,sss,bigstring,chiarr,nnn,njump,iupdate,
!$omp@ idxarr1,idxarr2,parmarray,zzzarr,zzz,probarr,oldchiarr,fitns)
!$omp@ firstprivate(Deltachi,outstring,prob,uran,idum,i1,i2,
!$omp@ lll,lll1,lll2,lll3,siggam,scalegam,oldparmarray)
!$omp@ private(parmstring,planetparm,dynparm,line)
!$omp do

            do 7755 i16=1,np
c
              if(isw22.ge.1)ifastflag=1
              if(i16.eq.1)ifastflag=0
              ifastflag=0   !disabled January 1, 2013
c
!$omp critical  
              do  j=1,Nterms
                var(j)=parmarray(j,i16)
              enddo
!$omp end critical
              chilimb=0.0d0
              if(ilimbcheck.gt.1)then
                if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @              dwavey,chilimb)
                if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @              dwavey,chilimb)
              endif
c
              ibest=0
              ichilabel=1
              call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     &          ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @          ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     &          xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,
     @          fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,
     @          chisqK,chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,
     &          ydataU,errU,zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,
     @          resB,NdataV,xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,
     $          ydataR,errR,zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,
     @          resI,NdataJ,xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,
     @          ydataH,errH,zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,
     @          resK,NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     @          ggamma1,ggamma2,Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,
     &          Nvmax,svar,var,saveasini,savxdataU,savydataU,saverrU,
     %          savxdataB,savydataB,saverrB,savxdataV,savydataV,saverrV,
     $          savxdataR,savydataR,saverrR,savxdataI,savydataI,saverrI,
     @          savxdataJ,savydataJ,saverrJ,savxdataH,savydataH,saverrH,
     @          savxdataK,savydataK,saverrK,savxRV1,savyRV1,saverrRV1,
     &          savxRV2,savyRV2,saverrRV2,ifrac,ilum,i16,isavNU,isavNB,
     @          isavNV,isavNR,isavNI,isavNJ,isavNH,isavNK,isavRV1,
     @          isavRV2,isvel1,isvel2,Ndatamax,ibest,ifixgamma,savesep,
     &          ichilabel,resRV1,resRV2,thresh,small,Ncycle,Ttimes,
     @          Tseps,Nobscycle,obsTtimes,obsTerr,icnarray,RV3,xRV3,
     @          yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,parmstring,
     %          planetparm,dynparm,line,fill1,fill2,omega1,omega2,
     @          dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,rinner,
     @          router,tdisk,xi,alb1,alb2,rLx,Period,fm,separ,gamma,
     @          wave,dbolx,dboly,dwavex,dwavey,t3,g3,SA3,density,sw1,
     @          sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     @          primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @          bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,
     @          sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,omegadot,
     @          contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49,
     @          gaplow,gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,
     @          P2incl,P2Omega,P2Q,P2ratrad,P3tconj,P3period,P3T0,
     @          P3ecos,P3esin,P3incl,P3Omega,P3Q,P3ratrad,P4tconj,
     @          P4period,P4T0,P4ecos,P4esin,P4incl,P4Omega,P4Q,P4ratrad,
     @          P5tconj,P5period,P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,
     &          P5ratrad,P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,
     &          P6Omega,P6Q,P6ratrad,P7tconj,P7period,P7T0,P7ecos,
     @          P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,P8period,
     @          P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,xSC,ySC,
     @          spot1parm,spot2parm,spotdparm,Nalph1,Nbet1,Nalph2,Nbet2,
     &          Ntheta,Nradius,Nref,idraw,iecheck,iidint,iatm,ism1,ilaw,
     @          icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
     @          iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,isw5,isw6,isw7,
     @          isw8,isw9,idark1,idark2,isw12,isw13,isw21,isw22,isw23,
     @          isw24,isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,
     @          isw33,isw34,NSC,compfracs,tertperiod,tertt0,tertecos,
     @          tertesin,tertincl,tertOmega,tertQ,tertconj,tertratrad,
     @          hh,sw72,sw73,Nmaxeclipse,Tdur1,Tdur2)
c
!$omp critical
c
              chiarr(i16)=chi1
              fitns(i16)=1.0d0/chi1    ! the fitness is prop.to. chi^2
              call printiter(i16,'iteration number = ',ig,
     @          'generation number = ')
c
              if(chi1.lt.small)then
                small=chi1
                do mmm=1,Nvmax
                  sss(mmm)=var(mmm)
                  iupdate=99
                enddo
              endif
              if(rmed.ge.1.0d0)then
                call printsmallmed(small)
              else
                call printsmall(small)
              endif
c
!$omp end critical
c
 7755      continue
!$omp enddo
!$omp end parallel
c 
c   Sort the chiarr and print out parameters for the best Nstep(2) sets.
c
            call rnkpop(np,fitns,ifit,jfit)
            call indexx(Np,chiarr,indxchi)
c
c   July 10, 2001
c
c   Add the "black sheep stuff here.  Take the worst individuals and
c   replace them with mutated copies of the best guy
c
            Nsheep=idnint(0.6d0*dble(np)*(1.0d0-pcross))+2
c
c   Increase the number of black sheep once the number of generations 
c   is larger than 75% of the population number.
c
            lll=idnint(0.75d0*dble(np))
            if(ig.gt.lll)Nsheep=idnint(0.7d0*dble(np)*(1.0d0-pcross))+2
            if(ig.gt.2*lll)Nsheep=
     @          idnint(0.8d0*dble(np)*(1.0d0-pcross))+2
            if(ig.gt.3*lll)Nsheep=
     @          idnint(1.0d0*dble(np)*(1.0d0-pcross))+2
            if(Nsheep.ge.np)Nsheep=np-1
c
            i16=i16-1               ! adjust the counter
c
            do  isheep=np-Nsheep,np
c
c   Add the ability to insert a population member here
c
              iopened=0
              if(isheep.eq.np)then
c
                call insgridinput(iopened,iNalph1,iNbet1,iNalph2,iNbet2,
     @            fill1,fill2,omega1,omega2,rdphase,Q,finc,Teff1,Teff2,
     @            Tgrav1,Tgrav2,betarim,rinner,router,tdisk,xi,iNtheta,
     @            iNradius,alb1,alb2,iNref,rLx,Period,fm,separ,gamma,t3,
     @            g3,SA3,density,sw1,sw2,sw3,T0,iidraw,iiecheck,iiidint,
     @            iiatm,iism1,iicnU,iicnB,iicnV,iicnR,iicnI,iicnJ,iicnH,
     @            iicnK,iiRVfilt,iisw1,iisw2,iisw3,iisw4,iilaw,wave,
     @            dbolx,dboly,dwavex,dwavey,ecc,argper,pshift,sw5,rsw6,
     &            rsw7,rsw8,rsw9,iikeep,iisynch,iisw5,iisw6,iisw7,iisw8,
     %            iisw9,spot1parm,spot2parm,spotdparm,primmass,primK,
     $            primrad,ratrad,frac1,frac2,ecosw,temprat,iidark1,
     @            iidark2,iisw12,iisw13,iisw21,iisw22,iisw23,iisw24,
     @            bigI,bigbeta,rsw23,rsw24,powercoeff,sw25,rsw26,rsw27,
     @            rsw28,rsw29,rsw30,contam,Tconj,beam1,beam2,iisw25,
     @            iisw26,iisw27,iisw28,iisw29,iisw30,iisw31,iisw32,
     @            iisw33,iisw34,ocose,osine,omegadot,contamS0,contamS1,
     &            contamS2,contamS3,sw47,sw48,sw49)
c
                if((isw30.ge.1).and.(isw7.ge.2).and.(iopened.gt.0))then
                  call insgridbody3(iopened,iNalph3,iNbet3,tertperiod,
     @              tertt0,tertecos,tertesin,tertincl,tertOmega,tertQ,
     @              dwavex,dwavey,itconj,it1,it2,it3,it4,tertconj,
     @              tertratrad,hh,sw72,sw73,P2tconj,P2period,P2T0,
     &              P2ecos,P2esin,P2incl,P2Omega,P2Q,P2ratrad,P3tconj,
     @              P3period,P3T0,P3ecos,P3esin,P3incl,P3Omega,P3Q,
     &              P3ratrad,P4tconj,P4period,P4T0,P4ecos,P4esin,P4incl,
     &              P4Omega,P4Q,P4ratrad,P5tconj,P5period,P5T0,P5ecos,
     @              P5esin,P5incl,P5Omega,P5Q,P5ratrad,P6tconj,P6period,
     @              P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,P6ratrad,
     @              P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,
     @              P7Q,P7ratrad,P8tconj,P8period,P8T0,P8ecos,P8esin,
     @              P8incl,P8Omega,P8Q,P8ratrad)
                endif
c
                if(iopened.gt.0)then
                  call varassign(Nvmax,svar,var,fill1,fill2,omega1,
     @              omega2,Q,finc,Teff1,Teff2,betarim,rinner,router,
     @              tdisk,xi,rLx,separ,gamma,t3,g3,sa3,ecc,argper,
     %              pshift,spot1parm,spot2parm,spotdparm,period,T0,alb1,
     @              alb2,dwavex,dwavey,primmass,primK,primrad,ratrad,
     &              frac1,frac2,ecosw,temprat,bigI,bigbeta,powercoeff,
     @              density,Tconj,beam1,beam2,contam,ocose,osine,isw29,
     @              tertperiod,tertt0,tertecos,tertesin,tertincl,
     @              tertOmega,tertQ,Tgrav1,Tgrav2,tertconj,omegadot,
     @              contamS0,contamS1,contamS2,contamS3,P2tconj,
     &              P2period,P2T0,P2ecos,P2esin,P2incl,P2Omega,P2Q,
     @              P2ratrad,P3tconj,P3period,P3T0,P3ecos,P3esin,P3incl,
     &              P3Omega,P3Q,P3ratrad,P4tconj,P4period,P4T0,P4ecos,
     @              P4esin,P4incl,P4Omega,P4Q,P4ratrad,P5tconj,P5period,
     @              P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     @              P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,
     @              P6Q,P6ratrad,P7tconj,P7period,P7T0,P7ecos,P7esin,
     @              P7incl,P7Omega,P7Q,P7ratrad,P8tconj,P8period,P8T0,
     @              P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,sw72,sw73)
c
                endif
              endif
c
              i1=int(ran1(idum)*2.0d0)+1
              if(i1.lt.1)i1=1
              if(i1.gt.2)i1=2
              do j=1,Nterms
                vx1=vstart(j)
                vx2=vstep(j)
                vxmult=(vx2-vx1)
                vxadd=vx1
c                urand=oldph(j,indxchi(1))+0.05d0*gasdev(idum)
                urand=oldph(j,indxchi(i1))+0.05d0*gasdev(idum)
                if(urand.gt.1.0d0)urand=1.0d0
                if(urand.lt.0.0d0)urand=0.0d0
                varx=urand*vxmult+vxadd
                if(iopened.gt.0)then
                  parmarray(j,indxchi(isheep))=var(j)
                  urand=(varx-vxadd)/vxmult
                  if(urand.gt.1.0d0)urand=1.0d0
                  if(urand.lt.0.0d0)urand=0.0d0
                else
                  parmarray(j,indxchi(isheep))=varx
                  var(j)=varx
                endif
                oldph(j,indxchi(isheep))=urand
              enddo
            enddo
c
c   We now have created a mutated copy of the fittest individual for the
c   member at the bottom.  Evaluate the fitness, etc.
c
!$omp parallel firstprivate(Nphase,xmod,ymodU,ymodB,ymodV,
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
!$omp@      ifrac,ilum,isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,
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
!$omp@ shared(small,sss,bigstring,chiarr,nnn,njump,iupdate,indxchi,
!$omp@ idxarr1,idxarr2,parmarray,zzzarr,zzz,probarr,oldchiarr,i16)
!$omp@ firstprivate(Deltachi,outstring,prob,uran,idum,i1,i2,
!$omp@ lll,lll1,lll2,lll3,siggam,scalegam,oldparmarray)
!$omp@ private(parmstring,planetparm,dynparm,line,iii)
!$omp do
            do 9950 isheep=np-Nsheep,np
              if(isw22.ge.1)ifastflag=1
c
!$omp critical
              i16=i16+1
c              write(*,*)i16,isheep+Nsheep
              iii=indxchi(isheep)
              do j=1,Nterms
                var(j)=parmarray(j,iii)
              enddo
!$omp end critical
c
              chilimb=0.0d0
              if(ilimbcheck.gt.1)then
                if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @            dwavey,chilimb)
                if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @            dwavey,chilimb)
              endif
c
              ibest=0
              ichilabel=1
              call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     &          ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @          ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @          xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,
     &          fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,
     @          chisqK,chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,
     &          ydataU,errU,zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,
     @          resB,NdataV,xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,
     @          ydataR,errR,zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,
     @          resI,NdataJ,xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,
     @          ydataH,errH,zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,
     @          resK,NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     @          ggamma1,ggamma2,Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,
     @          Nvmax,svar,var,saveasini,savxdataU,savydataU,saverrU,
     &          savxdataB,savydataB,saverrB,savxdataV,savydataV,saverrV,
     $          savxdataR,savydataR,saverrR,savxdataI,savydataI,saverrI,
     @          savxdataJ,savydataJ,saverrJ,savxdataH,savydataH,saverrH,
     @          savxdataK,savydataK,saverrK,savxRV1,savyRV1,saverrRV1,
     @          savxRV2,savyRV2,saverrRV2,ifrac,ilum,i16,isavNU,isavNB,
     @          isavNV,isavNR,isavNI,isavNJ,isavNH,isavNK,isavRV1,
     @          isavRV2,isvel1,isvel2,Ndatamax,ibest,ifixgamma,savesep,
     @          ichilabel,resRV1,resRV2,thresh,small,Ncycle,Ttimes,
     @          Tseps,Nobscycle,obsTtimes,obsTerr,icnarray,RV3,xRV3,
     @          yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,parmstring,
     @          planetparm,dynparm,line,fill1,fill2,omega1,omega2,
     @          dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,rinner,
     @          router,tdisk,xi,alb1,alb2,rLx,Period,fm,separ,gamma,
     @          wave,dbolx,dboly,dwavex,dwavey,t3,g3,SA3,density,sw1,
     @          sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     @          primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @          bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,
     @          sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,omegadot,
     @          contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49,
     &          gaplow,gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,
     &          P2incl,P2Omega,P2Q,P2ratrad,P3tconj,P3period,P3T0,
     @          P3ecos,P3esin,P3incl,P3Omega,P3Q,P3ratrad,P4tconj,
     @          P4period,P4T0,P4ecos,P4esin,P4incl,P4Omega,P4Q,P4ratrad,
     &          P5tconj,P5period,P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,
     &          P5ratrad,P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,
     @          P6Omega,P6Q,P6ratrad,P7tconj,P7period,P7T0,P7ecos,
     @          P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,P8period,
     @          P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,xSC,ySC,
     @          spot1parm,spot2parm,spotdparm,Nalph1,Nbet1,Nalph2,
     @          Nbet2,Ntheta,Nradius,Nref,idraw,iecheck,iidint,iatm,
     &          ism1,ilaw,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,
     &          icnRV1,icnRV2,iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,
     @          isw5,isw6,isw7,isw8,isw9,idark1,idark2,isw12,isw13,
     @          isw21,isw22,isw23,isw24,isw25,isw26,isw27,isw28,isw29,
     @          isw30,isw31,isw32,isw33,isw34,NSC,compfracs,tertperiod,
     &          tertt0,tertecos,tertesin,tertincl,tertOmega,tertQ,
     %          tertconj,tertratrad,hh,sw72,sw73,Nmaxeclipse,Tdur1,
     @          Tdur2)
c
!$omp critical
              chiarr(iii)=chi1
              fitns(iii)=1.0d0/chi1    ! the fitness is prop.to. chi^2
              call printiter(isheep+Nsheep+1,'iteration number = ',ig,
     @          'generation number = ')
c
              if(chi1.lt.small)then
                small=chi1
                do mmm=1,Nvmax
                  sss(mmm)=var(mmm)
                  iupdate=99
                enddo
              endif
              if(rmed.ge.1.0d0)then
                call printsmallmed(small)
              else
                call printsmall(small)
              endif
!$omp end critical
c
 9950       continue     ! do 9950 isheep=np,np-Nsheep,-1
!$omp end do
!$omp end parallel
c
c   Now rerank the population, and continue with the normal genetic 
c   procedure.
c
            call rnkpop(np,fitns,ifit,jfit)
            call indexx(Np,chiarr,indxchi)
c
c   UPDATE November 1, 2002
c
c   Add a new loop to fine tune the variables.  If we have had more than
c   50 generations then
c   step through them ONE
c   AT A TIME and tweak by 0.005*gasdev(idum).
c
            if(ig.lt.50)go to 9952  
  
            Nsheepsave=Nsheep
            Nsheep=Nterms
            ikount=0
c
            do isheep=np-Nsheep,np
              ikount=ikount+1
              if(isw22.ge.1)ifastflag=1
              i5=int(ran1(idum)*5.0d0)+1
              if(i5.lt.1)i5=1
              if(i5.gt.5)i5=5
              do j=1,Nterms
                vx1=vstart(j)
                vx2=vstep(j)
                vxmult=(vx2-vx1)
                vxadd=vx1
c                urand=oldph(j,indxchi(1))
c                if(j.eq.ikount)urand=oldph(j,indxchi(1))+
c     @               0.005d0*gasdev(idum)
                urand=oldph(j,indxchi(i5))
                if(j.eq.ikount)urand=oldph(j,indxchi(i5))+
     @               0.005d0*gasdev(idum)
                if(urand.gt.1.0d0)urand=1.0d0
                if(urand.lt.0.0d0)urand=0.0d0
                varx=urand*vxmult+vxadd
                parmarray(j,indxchi(isheep))=varx
                var(j)=varx
                oldph(j,indxchi(isheep))=urand
                urand=ran1(idum)
                if(urand.lt.0.5d0)then
                  i1=int(ran1(idum)*dble(Nterms))+1
                  if(i1.lt.1)i1=1
                  if(i1.gt.Nterms)i1=Nterms
                  if(i1.ne.j)then
c                     urand=oldph(i1,indxchi(1))+
c     @                  0.005d0*gasdev(idum)
                     urand=oldph(i1,indxchi(i5))+
     @                  0.005d0*gasdev(idum)
                    if(urand.gt.1.0d0)urand=1.0d0
                    if(urand.lt.0.0d0)urand=0.0d0
                    vx1=vstart(i1)
                    vx2=vstep(i1)
                    vxmult=(vx2-vx1)
                    vxadd=vx1
                    varx=urand*vxmult+vxadd
                    parmarray(i1,indxchi(isheep))=varx
                    var(i1)=varx
                    oldph(i1,indxchi(isheep))=urand
                  endif
                endif
              enddo
            enddo
c
c   We now have created a mutated copy of the fittest individual for the
c   member at the bottom.  Evaluate the fitness, etc.
c
!$omp parallel firstprivate(Nphase,xmod,ymodU,ymodB,ymodV,
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
!$omp@      ifrac,ilum,isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,
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
!$omp@ shared(small,sss,bigstring,chiarr,nnn,njump,iupdate,i16,
!$omp@ idxarr1,idxarr2,parmarray,zzzarr,zzz,probarr,oldchiarr)
!$omp@ firstprivate(Deltachi,outstring,prob,uran,idum,i1,i2,
!$omp@ lll,lll1,lll2,lll3,siggam,scalegam,oldparmarray)
!$omp@ private(parmstring,planetparm,dynparm,line,iii)
!$omp do
            do 9951 isheep=np-Nsheep,np
!$omp critical
              i16=i16+1
              iii=indxchi(isheep)
c
c              write(*,*)i16,isheep+Nsheep+Nsheepsave+2
              do j=1,Nterms
                var(j)=parmarray(j,iii)
              enddo
!$omp end critical
c
              chilimb=0.0d0
              if(ilimbcheck.gt.1)then
                if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @              dwavey,chilimb)
                if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @              dwavey,chilimb)
              endif
c 
              ibest=0
              ichilabel=1
              call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     &          ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @          ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @          xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,
     @          fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,
     %          chisqK,chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,
     %          ydataU,errU,zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,
     &          resB,NdataV,xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,
     &          ydataR,errR,zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,
     @          resI,NdataJ,xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,
     @          ydataH,errH,zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,
     @          resK,NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     @          ggamma1,ggamma2,Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,
     @          Nvmax,svar,var,saveasini,savxdataU,savydataU,saverrU,
     &          savxdataB,savydataB,saverrB,savxdataV,savydataV,saverrV,
     $          savxdataR,savydataR,saverrR,savxdataI,savydataI,saverrI,
     @          savxdataJ,savydataJ,saverrJ,savxdataH,savydataH,saverrH,
     &          savxdataK,savydataK,saverrK,savxRV1,savyRV1,saverrRV1,
     @          savxRV2,savyRV2,saverrRV2,ifrac,ilum,i16,isavNU,isavNB,
     &          isavNV,isavNR,isavNI,isavNJ,isavNH,isavNK,isavRV1,
     &          isavRV2,isvel1,isvel2,Ndatamax,ibest,ifixgamma,savesep,
     @          ichilabel,resRV1,resRV2,thresh,small,Ncycle,Ttimes,
     &          Tseps,Nobscycle,obsTtimes,obsTerr,icnarray,RV3,xRV3,
     @          yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,parmstring,
     @          planetparm,dynparm,line,fill1,fill2,omega1,omega2,
     @          dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,rinner,
     &          router,tdisk,xi,alb1,alb2,rLx,Period,fm,separ,gamma,
     &          wave,dbolx,dboly,dwavex,dwavey,t3,g3,SA3,density,sw1,
     $          sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     @          primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @          bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,
     @          sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,omegadot,
     @          contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49,
     @          gaplow,gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,
     @          P2incl,P2Omega,P2Q,P2ratrad,P3tconj,P3period,P3T0,
     @          P3ecos,P3esin,P3incl,P3Omega,P3Q,P3ratrad,P4tconj,
     @          P4period,P4T0,P4ecos,P4esin,P4incl,P4Omega,P4Q,P4ratrad,
     &          P5tconj,P5period,P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,
     &          P5ratrad,P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,
     &          P6Omega,P6Q,P6ratrad,P7tconj,P7period,P7T0,P7ecos,
     @          P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,P8period,
     &          P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,xSC,ySC,
     &          spot1parm,spot2parm,spotdparm,Nalph1,Nbet1,Nalph2,Nbet2,
     &          Ntheta,Nradius,Nref,idraw,iecheck,iidint,iatm,ism1,ilaw,
     @          icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
     @          iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,isw5,isw6,isw7,
     @          isw8,isw9,idark1,idark2,isw12,isw13,isw21,isw22,isw23,
     &          isw24,isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,
     &          isw33,isw34,NSC,compfracs,tertperiod,tertt0,tertecos,
     &          tertesin,tertincl,tertOmega,tertQ,tertconj,tertratrad,
     &          hh,sw72,sw73,Nmaxeclipse,Tdur1,Tdur2)
c
!$omp critical
              chiarr(iii)=chi1
              fitns(iii)=1.0d0/chi1    ! the fitness is prop.to. chi^2
              call printiter(isheep+Nsheep+Nsheepsave+2,
     @          'iteration number = ',ig,'generation number = ')
c
              if(chi1.lt.small)then
                small=chi1
                do mmm=1,Nvmax
                  sss(mmm)=var(mmm)
                  iupdate=99
                enddo
              endif
              if(rmed.ge.1.0d0)then
                call printsmallmed(small)
              else
                call printsmall(small)
              endif
!$omp end critical
c
 9951       continue             ! do 9950 isheep=np,np-Nterms,-1
!$omp enddo
!$omp end parallel
c
c   Now rerank the population, and continue with the normal genetic 
c   procedure.
c
            call rnkpop(np,fitns,ifit,jfit)
            call indexx(Np,chiarr,indxchi)
c
 9952       write(38,8850)ig,chiarr(indxchi(1)),chiarr(indxchi(np/2)),
     @          pmut
c
 8850       format(i4,2x,2(f13.7,2x),2x,f6.4)
            call adjmut(np,fitns,ifit,pmutmn,pmutmx,pmut)
c
c   reset the variables at their best values and print the chi^2
c
            if(iupdate.gt.0)then
              do mmm=1,Nvmax
                var(mmm)=sss(mmm)
              enddo
c
              chilimb=0.0d0
              if(ilimbcheck.gt.1)then
                if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @             dwavey,chilimb)
                if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @             dwavey,chilimb)
              endif
c
              ibest=99
              ichilabel=1
              call monster(Nphase,Nmaxphase,xmod,ymodU,ymodB,ymodV,
     &          ymodR,ymodI,ymodJ,ymodH,ymodK,ymods1,ymods2,ymods3,
     @          ymodd,RV1,RV2,drv1,drv2,obsparm,ifastflag,NRVphase,
     @          xRVmod,fracs1,fracs2,fracs3,fracs4,fracs5,fracs6,fracs7,
     &          fracs8,chisqU,chisqB,chisqV,chisqR,chisqI,chisqJ,chisqH,
     @          chisqK,chisqRV1,chisqRV2,chilimb,chi1,NdataU,xdataU,
     &          ydataU,errU,zeroU,resU,NdataB,xdataB,ydataB,errB,zeroB,
     @          resB,NdataV,xdataV,ydataV,errV,zeroV,resV,NdataR,xdataR,
     @          ydataR,errR,zeroR,resR,NdataI,xdataI,ydataI,errI,zeroI,
     &          resI,NdataJ,xdataJ,ydataJ,errJ,zeroJ,resJ,NdataH,xdataH,
     @          ydataH,errH,zeroH,resH,NdataK,xdataK,ydataK,errK,zeroK,
     @          resK,NRV1,xRV1,yRV1,errRV1,NRV2,xRV2,yRV2,errRV2,
     @          ggamma1,ggamma2,Nobv,sobv,obv,eobv,ochi,ochidisk,ochilr,
     &          Nvmax,svar,var,saveasini,savxdataU,savydataU,saverrU,
     &          savxdataB,savydataB,saverrB,savxdataV,savydataV,saverrV,
     $          savxdataR,savydataR,saverrR,savxdataI,savydataI,saverrI,
     @          savxdataJ,savydataJ,saverrJ,savxdataH,savydataH,saverrH,
     @          savxdataK,savydataK,saverrK,savxRV1,savyRV1,saverrRV1,
     @          savxRV2,savyRV2,saverrRV2,ifrac,ilum,i16,isavNU,isavNB,
     @          isavNV,isavNR,isavNI,isavNJ,isavNH,isavNK,isavRV1,
     @          isavRV2,isvel1,isvel2,Ndatamax,ibest,ifixgamma,savesep,
     @          ichilabel,resRV1,resRV2,thresh,small,Ncycle,Ttimes,
     @          Tseps,Nobscycle,obsTtimes,obsTerr,icnarray,RV3,xRV3,
     @          yRV3,errRV3,icnRV3,resRV3,ggamma3,NRV3,parmstring,
     @          planetparm,dynparm,line,fill1,fill2,omega1,omega2,
     @          dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,betarim,rinner,
     @          router,tdisk,xi,alb1,alb2,rLx,Period,fm,separ,gamma,
     @          wave,dbolx,dboly,dwavex,dwavey,t3,g3,SA3,density,sw1,
     @          sw2,sw3,T0,ecc,argper,pshift,sw5,sw6,sw7,sw8,sw9,
     @          primmass,primK,primrad,ratrad,frac1,frac2,ecosw,temprat,
     @          bigI,bigbeta,sw23,sw24,powercoeff,sw25,sw26,sw27,sw28,
     @          sw29,sw30,contam,Tconj,beam1,beam2,ocose,osine,omegadot,
     @          contamS0,contamS1,contamS2,contamS3,sw47,sw48,sw49,
     @          gaplow,gaphigh,P2tconj,P2period,P2T0,P2ecos,P2esin,
     @          P2incl,P2Omega,P2Q,P2ratrad,P3tconj,P3period,P3T0,
     @          P3ecos,P3esin,P3incl,P3Omega,P3Q,P3ratrad,P4tconj,
     @          P4period,P4T0,P4ecos,P4esin,P4incl,P4Omega,P4Q,P4ratrad,
     @          P5tconj,P5period,P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,
     @          P5ratrad,P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,
     @          P6Omega,P6Q,P6ratrad,P7tconj,P7period,P7T0,P7ecos,
     @          P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,P8period,
     @          P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,P8ratrad,xSC,ySC,
     @          spot1parm,spot2parm,spotdparm,Nalph1,Nbet1,Nalph2,Nbet2,
     @          Ntheta,Nradius,Nref,idraw,iecheck,iidint,iatm,ism1,ilaw,
     @          icnU,icnB,icnV,icnR,icnI,icnJ,icnH,icnK,icnRV1,icnRV2,
     @          iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,isw5,isw6,isw7,
     @          isw8,isw9,idark1,idark2,isw12,isw13,isw21,isw22,isw23,
     @          isw24,isw25,isw26,isw27,isw28,isw29,isw30,isw31,isw32,
     @          isw33,isw34,NSC,compfracs,tertperiod,tertt0,tertecos,
     @          tertesin,tertincl,tertOmega,tertQ,tertconj,tertratrad,
     @          hh,sw72,sw73,Nmaxeclipse,Tdur1,Tdur2)
c
              chiall=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+
     @           chisqK+chisqRV1+chisqRV2+ochi+chilimb)
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
     @           'modelU.linear',isw7)
              if(icnB.ne.430)call wlinmod(Nphase,xmod,ymodB,
     @           'modelB.linear',isw7)
              if(icnV.ne.430)call wlinmod(Nphase,xmod,ymodV,
     @           'modelV.linear',isw7)
              if(icnR.ne.430)call wlinmod(Nphase,xmod,ymodR,
     @           'modelR.linear',isw7)
              if(icnI.ne.430)call wlinmod(Nphase,xmod,ymodI,
     @           'modelI.linear',isw7)
              if(icnJ.ne.430)call wlinmod(Nphase,xmod,ymodJ,
     @           'modelJ.linear',isw7)
              if(icnH.ne.430)call wlinmod(Nphase,xmod,ymodH,
     @           'modelH.linear',isw7)
              if(icnK.ne.430)call wlinmod(Nphase,xmod,ymodK,
     @           'modelK.linear',isw7)
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
     &          omega1,omega2,dphase,Q,finc,Teff1,Teff2,Tgrav1,Tgrav2,
     @          betarim,rinner,router,tdisk,xi,Ntheta,Nradius,alb1,alb2,
     @          Nref,rLx,Period,fm,separ,gamma1,t3,g3,SA3,density,sw1,
     @          sw2,sw3,T0,idraw,iecheck,iidint,iatm,ism1,icnU,icnB,
     @          icnV,icnR,icnI,icnJ,icnH,icnK,iRVfilt,isw1,isw2,isw3,
     @          isw4,ilaw,wave,dbolx,dboly,dwavex,dwavey,ecc,argper,
     @          pshift,sw5,sw6,sw7,sw8,sw9,ikeep,isynch,isw5,isw6,isw7,
     @          isw8,isw9,spot1parm,spot2parm,spotdparm,primmass,primK,
     @          primrad,ratrad,frac1,frac2,ecosw,temprat,idark1,idark2,
     @          isw12,isw13,isw21,isw22,isw23,isw24,bigI,bigbeta,sw23,
     @          sw24,powercoeff,sw25,sw26,sw27,sw28,sw29,sw30,contam,
     @          Tconj,beam1,beam2,isw25,isw26,isw27,isw28,isw29,isw30,
     @          isw31,isw32,isw33,isw34,ocose,osine,omegadot,contamS0,
     @          contamS1,contamS2,contamS3,sw47,sw48,sw49)
c
              if(isw30.gt.0)call writebody3grid(Nalph3,Nbet3,tertperiod,
     @           tertt0,tertecos,tertesin,tertincl,tertOmega,tertQ,
     @           dwavex,dwavey,itconj,it1,it2,it3,it4,tertconj,
     @           tertratrad,hh,sw72,sw73,P2tconj,P2period,P2T0,P2ecos,
     @           P2esin,P2incl,P2Omega,P2Q,P2ratrad,P3tconj,P3period,
     @           P3T0,P3ecos,P3esin,P3incl,P3Omega,P3Q,P3ratrad,P4tconj,
     @           P4period,P4T0,P4ecos,P4esin,P4incl,P4Omega,P4Q,
     @           P4ratrad,P5tconj,P5period,P5T0,P5ecos,P5esin,P5incl,
     @           P5Omega,P5Q,P5ratrad,P6tconj,P6period,P6T0,P6ecos,
     @           P6esin,P6incl,P6Omega,P6Q,P6ratrad,P7tconj,P7period,
     @           P7T0,P7ecos,P7esin,P7incl,P7Omega,P7Q,P7ratrad,P8tconj,
     @           P8period,P8T0,P8ecos,P8esin,P8incl,P8Omega,P8Q,
     @           P8ratrad)
c
              call recordloopopt(Udatafile,Bdatafile,Vdatafile,
     &           Rdatafile,Idatafile,Jdatafile,Hdatafile,Kdatafile,
     @           RV1file,RV2file,Nvmax,Nvar,svar,var,vstart,stepsave,
     @           Nstep,Nobv,sobv,obv,eobv,vstep)
c
              write(command,101)ig+1000000
              call system(command)
              write(command,102)ig+1000000
              call system(command)
              write(command,103)ig+1000000
              if(isw30.gt.0)call system(command)
              iupdate=0
            endif  !if iupdate > 0
c
            close(45)
            close(46)
            if(isw24.gt.0)close(47)
            if(isw30.gt.0)close(48)
            if(isw30.ge.3)close(49)
            close(55)

10000     continue
c
          call indexx(Np,chiarr,indxchi)
c
c   reset the variables at their best values and print the chi^2
c
          do mmm=1,Nvmax
            var(mmm)=sss(mmm)
          enddo
c
          if(savesep.lt.0.0d0)separ=savesep
c
          ifastflag=0
          chilimb=0.0d0
          if(ilimbcheck.gt.1)then
            if(ilaw.eq.4)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @         dwavey,chilimb)
            if(ilaw.eq.14)call getchilimb(Nvmax,Nterms,svar,dwavex,
     @         dwavey,chilimb)
          endif
c
          ibest=99
          ichilabel=1
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
     &      ochidisk,ochilr,Nvmax,svar,var,saveasini,savxdataU,
     @      savydataU,saverrU,savxdataB,savydataB,saverrB,savxdataV,
     @      savydataV,saverrV,savxdataR,savydataR,saverrR,savxdataI,
     @      savydataI,saverrI,savxdataJ,savydataJ,saverrJ,savxdataH,
     @      savydataH,saverrH,savxdataK,savydataK,saverrK,savxRV1,
     @      savyRV1,saverrRV1,savxRV2,savyRV2,saverrRV2,ifrac,ilum,i16,
     @      isavNU,isavNB,isavNV,isavNR,isavNI,isavNJ,isavNH,isavNK,
     @      isavRV1,isavRV2,isvel1,isvel2,Ndatamax,ibest,ifixgamma,
     @      savesep,ichilabel,resRV1,resRV2,thresh,small,Ncycle,Ttimes,
     @      Tseps,Nobscycle,obsTtimes,obsTerr,icnarray,RV3,xRV3,yRV3,
     @      errRV3,icnRV3,resRV3,ggamma3,NRV3,parmstring,planetparm,
     @      dynparm,line,fill1,fill2,omega1,omega2,dphase,Q,finc,Teff1,
     @      Teff2,Tgrav1,Tgrav2,betarim,rinner,router,tdisk,xi,alb1,
     @      alb2,rLx,Period,fm,separ,gamma,wave,dbolx,dboly,dwavex,
     @      dwavey,t3,g3,SA3,density,sw1,sw2,sw3,T0,ecc,argper,pshift,
     @      sw5,sw6,sw7,sw8,sw9,primmass,primK,primrad,ratrad,frac1,
     @      frac2,ecosw,temprat,bigI,bigbeta,sw23,sw24,powercoeff,sw25,
     @      sw26,sw27,sw28,sw29,sw30,contam,Tconj,beam1,beam2,ocose,
     &      osine,omegadot,contamS0,contamS1,contamS2,contamS3,sw47,
     @      sw48,sw49,gaplow,gaphigh,P2tconj,P2period,P2T0,P2ecos,
     @      P2esin,P2incl,P2Omega,P2Q,P2ratrad,P3tconj,P3period,P3T0,
     @      P3ecos,P3esin,P3incl,P3Omega,P3Q,P3ratrad,P4tconj,P4period,
     @      P4T0,P4ecos,P4esin,P4incl,P4Omega,P4Q,P4ratrad,P5tconj,
     @      P5period,P5T0,P5ecos,P5esin,P5incl,P5Omega,P5Q,P5ratrad,
     &      P6tconj,P6period,P6T0,P6ecos,P6esin,P6incl,P6Omega,P6Q,
     @      P6ratrad,P7tconj,P7period,P7T0,P7ecos,P7esin,P7incl,P7Omega,
     @      P7Q,P7ratrad,P8tconj,P8period,P8T0,P8ecos,P8esin,P8incl,
     @      P8Omega,P8Q,P8ratrad,xSC,ySC,spot1parm,spot2parm,spotdparm,
     @      Nalph1,Nbet1,Nalph2,Nbet2,Ntheta,Nradius,Nref,idraw,iecheck,
     @      iidint,iatm,ism1,ilaw,icnU,icnB,icnV,icnR,icnI,icnJ,icnH,
     @      icnK,icnRV1,icnRV2,iRVfilt,isw1,isw2,isw3,isw4,ikeep,isynch,
     @      isw5,isw6,isw7,isw8,isw9,idark1,idark2,isw12,isw13,isw21,
     @      isw22,isw23,isw24,isw25,isw26,isw27,isw28,isw29,isw30,isw31,
     @      isw32,isw33,isw34,NSC,compfracs,tertperiod,tertt0,tertecos,
     @      tertesin,tertincl,tertOmega,tertQ,tertconj,tertratrad,hh,
     @      sw72,sw73,Nmaxeclipse,Tdur1,Tdur2)
c
          chiall=(chisqU+chisqB+chisqV+chisqR+chisqI+chisqJ+chisqH+
     @         chisqK+chisqRV1+chisqRV2+ochi+chilimb)
          write(*,*)' '
          chiall=chi1
          call printS(chiall)
c
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

          call wlinmod(Nphase,xmod,ymods1,'lcstar1.linear',isw7)
          call wlinmod(Nphase,xmod,ymods2,'lcstar2.linear',isw7)
          call wlinmod(Nphase,xmod,ymods3,'lcstar3.linear',isw7)
          if(iidint.ge.1)call wlinmod(Nphase,xmod,ymodd,'lcdisk.linear',
     @          isw7)
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
c  
c
c &&&&&&&&&&&&&
c
          include 'lcsubs.for'
          include 'optimizesubs.for'
          include 'dynamicssubs.for'
