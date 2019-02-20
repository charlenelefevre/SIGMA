      SUBROUTINE GRAIN_DIST_WD01(INDEX,DTYPE,A,DNDA)
      IMPLICIT NONE

!History:  Feb 12, 2002:  Corrected bug--previous version did not multiply
!          very small carbonaceous grain size dist by b_C.

!from http://physics.gmu.edu/~joe/files/grain_dist.f
!Originally from http://cita.utoronto.ca/~weingart

!Output:   DNDA = (1/n_H) dn_gr/da, where n_gr(a) is the number density
!                 of grains with radius < a and n_H is the H nucleus
!                 number density   (cm^-1)
!Input:   DTYPE = 1 for silicate dust
!                 2 for carbonaceous dust
!             A = grain radius (cm)

!          INDEX  is an integer indicating which set of conditions to
!          adopt

!        INDEX    R_V    10^5 bc case

!          1      3.1      0.0       A
!          2      3.1      1.0       A
!          3      3.1      2.0       A
!          4      3.1      3.0       A
!          5      3.1      4.0       A
!          6      3.1      5.0       A
!          7      3.1      6.0       A
!          8      4.0      0.0       A
!          9      4.0      1.0       A
!          10     4.0      2.0       A
!          11     4.0      3.0       A
!          12     4.0      4.0       A
!          13     5.5      0.0       A
!          14     5.5      1.0       A
!          15     5.5      2.0       A
!          16     5.5      3.0       A
!          17     4.0      0.0       B
!          18     4.0      1.0       B
!          19     4.0      2.0       B
!          20     4.0      3.0       B
!          21     4.0      4.0       B
!          22     5.5      0.0       B
!          23     5.5      1.0       B
!          24     5.5      2.0       B
!          25     5.5      3.0       B

      INTEGER INDEX,DTYPE
      REAL A,DNDA

      REAL ALPHAGARR(25),BETAGARR(25),ATGARR(25),ACGARR(25),CGARR(25),
     &     ALPHASARR(25),BETASARR(25),ATSARR(25),CSARR(25),ACS,
     &     BC5ARR(25),ALPHAG,BETAG,ATG,ACG,CG,ALPHAS,BETAS,ATS,CS,BC5,
     &     DNDAVSG,A01,A02,SIG,B1,B2

      DATA ALPHAGARR/
     &-2.25,-2.17,-2.04,-1.91,-1.84,-1.72,-1.54,-2.26,-2.16,-2.01,
     &-1.83,-1.64,-2.35,-2.12,-1.94,-1.61,-2.62,-2.52,-2.36,-2.09,
     &-1.96,-2.80,-2.67,-2.45,-1.90/
      DATA BETAGARR/
     &-0.0648,-0.0382,-0.111,-0.125,-0.132,-0.322,-0.165,-0.199,-0.0862,
     &-0.0973,-0.175,-0.247,-0.668,-0.67,-0.853,-0.722,-0.0144,-0.0541,
     &-0.0957,-0.193,-0.813,0.0356,0.0129,-0.00132,-0.0517/
      DATA ATGARR/
     &0.00745,0.00373,0.00828,0.00837,0.00898,0.0254,.0107,0.0241,
     &0.00867,0.00811,0.0117,0.0152,0.148,0.0686,0.0786,0.0418,0.0187,
     &0.0366,0.0305,0.0199,0.0693,0.0203,0.0134,0.0275,0.012/
      DATA ACGARR/
     &0.606,0.586,0.543,0.499,0.489,0.438,0.428,0.861,0.803,0.696,
     &0.604,0.536,1.96,1.35,0.921,0.72,5.74,6.65,6.44,4.6,3.48,3.43,
     &3.44,5.14,7.28/
      DATA CGARR/
     &9.94e-11,3.79e-10,5.57e-11,4.15e-11,2.90e-11,3.20e-12,9.99e-12,
     &5.47e-12,4.58e-11,3.96e-11,1.42e-11,5.83e-12,4.82e-14,3.65e-13,
     &2.57e-13,7.58e-13,6.46e-12,1.08e-12,1.62e-12,4.21e-12,2.95e-13,
     &2.74e-12,7.25e-12,8.79e-13,2.86e-12/
      DATA ALPHASARR/
     &-1.48,-1.46,-1.43,-1.41,-2.1,-2.1,-2.21,-2.03,-2.05,-2.06,-2.08,
     &-2.09,-1.57,-1.57,-1.55,-1.59,-2.01,-2.11,-2.05,-2.1,-2.11,-1.09,
     &-1.14,-1.08,-1.13/
      DATA BETASARR/
     &-9.34,-10.3,-11.7,-11.5,-0.114,-0.0407,0.3,0.668,0.832,0.995,
     &1.29,1.58,1.1,1.25,1.33,2.12,0.894,1.58,1.19,1.64,2.1,-0.37,
     &-0.195,-0.336,-0.109/
      DATA ATSARR/
     &0.172,0.174,0.173,0.171,0.169,0.166,0.164,0.189,0.188,0.185,
     &0.184,0.183,0.198,0.197,0.195,0.193,0.198,0.197,0.197,0.198,0.198,
     &0.218,0.216,0.216,0.211/
      DATA CSARR/
     &1.02e-12,1.09e-12,1.27e-12,1.33e-12,1.26e-13,1.27e-13,1.e-13,
     &5.2e-14,4.81e-14,4.7e-14,4.26e-14,3.94e-14,4.24e-14,4.e-14,
     &4.05e-14,3.2e-14,4.95e-14,3.69e-14,4.37e-14,3.63e-14,3.13e-14,
     &1.17e-13,1.05e-13,1.17e-13,1.04e-13/
      DATA BC5ARR/
     &0.,1.,2.,3.,4.,5.,6.,0.,1.,2.,3.,4.,0.,1.,2.,3.,0.,1.,2.,3.,4.,
     &0.,1.,2.,3./

      ALPHAG=ALPHAGARR(INDEX)
      BETAG=BETAGARR(INDEX)
      ATG=ATGARR(INDEX)*1.E-4
      ACG=ACGARR(INDEX)*1.E-4
      CG=CGARR(INDEX)
      ALPHAS=ALPHASARR(INDEX)
      BETAS=BETASARR(INDEX)
      ATS=ATSARR(INDEX)*1.E-4
      ACS=1.E-5
      CS=CSARR(INDEX)
      BC5=BC5ARR(INDEX)
      IF (DTYPE .EQ. 1) THEN
       DNDA=(CS/A)*(A/ATS)**ALPHAS
       IF (BETAS .GE. 0.) THEN
        DNDA=DNDA*(1.+BETAS*A/ATS)
       ELSE
        DNDA=DNDA/(1.-BETAS*A/ATS)
       ENDIF
       IF (A .GT. ATS) DNDA=DNDA*EXP(((ATS-A)/ACS)**3)
      ENDIF
      IF (DTYPE .EQ. 2) THEN
       DNDA=(CG/A)*(A/ATG)**ALPHAG
       IF (BETAG .GE. 0.) THEN
        DNDA=DNDA*(1.+BETAG*A/ATG)
       ELSE
        DNDA=DNDA/(1.-BETAG*A/ATG)
       ENDIF
       IF (A .GT. ATG) DNDA=DNDA*EXP(((ATG-A)/ACG)**3)
       A01=3.5E-8
       A02=3.E-7
       SIG=0.4
       B1=2.0496E-7
       B2=9.6005E-11
       DNDAVSG=(B1/A)*EXP(-0.5*(LOG(A/A01)/SIG)**2)+
     &         (B2/A)*EXP(-0.5*(LOG(A/A02)/SIG)**2)
       IF (DNDAVSG .GE. 0.0001*DNDA) DNDA=DNDA+BC5*DNDAVSG
      ENDIF

      RETURN
      END
