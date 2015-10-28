cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 1st order propagation for:
c 1. Dispersion (Sellmeier form + resonance in Lorentz shape)
c 2. Instantaneous chi3
c 
c Main structure based on uv2ph.f 2014-04-09 version
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      MODULE Global
      IMPLICIT NONE
      REAL,PARAMETER::
     & printsmall=1E-50,
     & cFT=0.39894,c0=137.036,aCoulm=8.47837E-30,aCoul=1.60218E-19,
     & aHz=6.57969E15,aPHz=6.57969,aGHz=6.57969E6,
     & aeV=27.211384,afs=0.024188843,a0s=2.418843E-17,lam=45.563353,
     & In0cyc=3.51E16,                                                  ! cycle-averaged
     & In0abs=7.02E16,                                                  ! absolute field
     & In0fun=6.44E15,                                                  ! fundamental units combined
     & aum=5.2917721E-5,acm=5.2917721E-9,am=5.2917721E-11,
     & acm3=6.7483346E24,torr=9.6564746E18,
     & aJ=4.35975E-18,aVm=5.14221E11,
     & u11=2.4048256,a11=8.6584065,Ieff0=0.56553644,Iave0=0.26951411    ! Iave0: I_ave / I_center
      COMPLEX,PARAMETER::i=(0.0,1.0),clx=(1.0,0.0)
      CHARACTER(LEN=40)::title
      INTEGER::
     & iw,Nw,Lw,it,Nt,Lt,iz,Lz,Ncr,Ntee,KTcr,KWcr,NsR
      REAL::
     & pi,eps,w1,w2,dw,t1,t2,dt,fbrL,fbrR,fbrA,Aeff,z,dz,dzl,dzh,den,
     & U0J,Td,w0,CEP,tolz0,tolzlim,
     & kt,temp,pbar,cput0,cput1,Tcr,T0cr,Wcr,W0cr,ng,dtt,tolz,
     & dznext,tolR,Om,wR,KR,n2,coeff2
      REAL,ALLOCATABLE::
     & t(:),w(:),k(:),Cr(:),Ww(:),GW(:),GB(:)
      COMPLEX,ALLOCATABLE::
     & ewt(:,:),e0t(:),bg(:),chi1(:),E(:)
      CONTAINS

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Input
      INTEGER::ios
      REAL::t1fs,t2fs,Tdfs,Lm,Dm,Tcrfs,T0crfs,dtfs,dttfs,
     & dzlm,dzhm,TRfs,OmGHz
c Wcr & Tcr = width of margin
c W0cr & T0cr = width of supergaussian
c KWcr & KTcr = order of supergaussian
      NAMELIST/NUMR/w1,w2,dw,Wcr,W0cr,KWcr,
     &              dtfs,t1fs,t2fs,Tcrfs,T0crfs,KTcr,
     &              dzlm,dzhm,tolz0,tolzlim
      NAMELIST/SYST/pbar,temp,Lm,Dm,OmGHz,wR,TRfs,n2
      NAMELIST/PUMP/w0,U0J,Tdfs,CEP
      NAMELIST/SHOW/dttfs,Lz,Lw
      READ(*,'(A48)') title
      READ(*,NML=NUMR)
      READ(*,NML=SYST)
      READ(*,NML=PUMP)
      READ(*,NML=SHOW)
c numerical
      Tcr=Tcrfs/afs
      T0cr=T0crfs/afs
      dzl=dzlm/am
      dzh=dzhm/am
      t1=t1fs/afs
      t2=t2fs/afs
      dt=dtfs/afs
      Ncr=NINT(Tcr/dt)
c system
      !den=pbar*750.06/temp*torr/acm3
      !WRITE(*,*) 'density(atu)=',den
      fbrL=Lm/am
      fbrR=Dm/am/2.0
      fbrA=pi*fbrR*fbrR
      kt=u11/fbrR
      Aeff=a11/kt/kt
      KR=afs/TRfs/2.0
      Om=OmGHz/aGHz
c pump
      Td=Tdfs/afs
c display
      dtt=dttfs/afs
      Ntee=4000
      WRITE(*,*) 'Input done'
 15   FORMAT(10(A14,ES14.5E3))
      RETURN
      END SUBROUTINE Input

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initial values
      SUBROUTINE Dimen
c time
      Nt=NINT((t2-t1)/dt)
      ALLOCATE(t(Nt),Cr(Nt),E(Nt))
      DO it=1,Nt
        t(it)=t1+dt*REAL(it-1)
      END DO
      WRITE(*,*) 'Nt =',Nt
      IF (dtt<=dt) THEN
        Lt=1
      ELSE IF (dtt>dt) THEN
        Lt=NINT(dtt/dt)
      END IF
      IF ((Nt-1)/Lt+1>5000) WRITE(*,*) 'not enough rows in matrix'
c frequency
      Nw=NINT((w2-w1)/dw)
      ALLOCATE(w(Nw),k(Nw),Ww(Nw),bg(Nw),chi1(Nw))
      DO iw=1,nw
        w(iw)=w1+dw*REAL(iw-1)
        k(iw)=w(iw)/c0
      END DO
      WRITE(*,*) 'Nw =',Nw
c Full field: exp(-iEt), exp(iwt), exp(-iw0t)
      ALLOCATE(ewt(Nw,Nt),e0t(Nt))
      FORALL (iw=1:Nw,it=1:Nt)
        ewt(iw,it)=EXP(i*w(iw)*t(it))
      END FORALL
      FORALL (it=1:Nt)
        e0t(it)=EXP(-i*w0*t(it))
      END FORALL
c prepare supporting functions for propagation
      CALL CrpFn
      CALL Dispersion
      CALL CnvFn
      coeff2=2.0*n2/c0*pbar*temp/273.0
      WRITE(*,*) 'Dimen done'
 10   FORMAT(10ES14.5E3)
 15   FORMAT(A14,10ES14.5E3)
      RETURN
      END SUBROUTINE Dimen

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Dispersion
      INTEGER::ios,iomeg,Nomeg
      REAL::beta1,beta2,xx
      REAL,ALLOCATABLE::omeg(:),beta(:)
      ALLOCATE(beta(Nw))
c chi1: Sellmeier and resonance
      DO iw=1,Nw
        xx=(w(iw)*1000.0/lam)**2 ! wavelength for Sellmeier
        chi1(iw)=pbar*273.0/temp*2.0*0.012055*(0.26783/(46.301-xx)
     &                                        +0.29481/(50.578-xx)
     &                                        +5.0333/(112.74-xx))
     &          +Om*Om/(wR*wR-w(iw)*w(iw)-i*w(iw)*2.0*KR)
      END DO
c beta, beta1, beta2, ng
      beta=SQRT(k*k*(1.0+REAL(chi1))-kt*kt)
      OPEN(40,FILE='X1.dat',STATUS='UNKNOWN')
      WRITE(40,5) 'f(Hz)','wlen(nm)','Re[X1]','Im[X1]'
      WRITE(41,5) 'f(Hz)','wlen(nm)','beta(1/m)','beta1(fs/m)',
     & 'beta2(fs2/m)'
      DO iw=2,Nw-1
        beta1=(beta(iw+1)-beta(iw-1))/2.0/dw
        IF (ABS(w(iw)-w0)<0.5*dw) ng=beta1*c0
        beta2=(beta(iw+1)+beta(iw-1)-2.0*beta(iw))/dw/dw
        WRITE(40,10) w(iw)*aHz,lam/w(iw),chi1(iw)
        WRITE(41,10) w(iw)*aHz,lam/w(iw),beta(iw)/am,beta1*afs/am,
     &               beta2*afs*afs/am
      END DO
      WRITE(40,*) '# ng =',ng
      CLOSE(40)
      CLOSE(41)
 5    FORMAT(10A14)
 10   FORMAT(10ES14.5E3)
      RETURN
      END SUBROUTINE Dispersion

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c supergaussion function
      FUNCTION SGau(s,L,x)
      INTEGER::L
      REAL::SGau,s,x
      SGau=EXP(-0.5*(ABS(x)/s)**L)
      RETURN
      END FUNCTION SGau

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c supergaussion boundary processor
      SUBROUTINE SGbound(x,y,Nx,xlo,xhi,s,L)
      INTEGER::L,ix,Nx
      REAL::xlo,xhi,s
      REAL,ALLOCATABLE::x(:)
      COMPLEX,ALLOCATABLE::y(:)
      DO ix=1,Nx
        IF (x(ix)<xlo) THEN
          y(ix)=y(ix)*SGau(s,L,x(ix)-xlo)
        ELSE IF (x(ix)>xhi) THEN
          y(ix)=y(ix)*SGau(s,L,x(ix)-xhi)
        END IF
      END DO
      RETURN
      END SUBROUTINE SGbound

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c create supergaussian cropping functions in time and in frequency
      SUBROUTINE CrpFn
      COMPLEX,ALLOCATABLE::Crx(:),Wwx(:)
c supergaussian cropping in time
      OPEN(213,FILE='Cr.dat',STATUS='UNKNOWN')
      ALLOCATE(Crx(Nt))
      Crx=(1.0,0.0)
      CALL SGbound(t,Crx,Nt,t(1)+Tcr,t(Nt)-Tcr,T0cr,KTcr)
      Cr=REAL(Crx)
      DEALLOCATE(Crx)
      DO it=1,Nt
        WRITE(213,10) t(it)*afs,Cr(it)
      END DO
      CLOSE(213)
c supergaussian cropping in frequency
      OPEN(215,FILE='Ww.dat',STATUS='UNKNOWN')
      ALLOCATE(Wwx(Nw))
      Wwx=(1.0,0.0)
      CALL SGbound(w,Wwx,Nw,w(1)+Wcr,w(Nw)-Wcr,W0cr,Kwcr)
      Ww=REAL(Wwx)
      DEALLOCATE(Wwx)
      DO iw=1,Nw
        WRITE(215,10) w(iw)*aPHz,Ww(iw)
      END DO
      CLOSE(215)
      WRITE(*,*) 'Cr(t), Ww(w) prepared'
 10   FORMAT(10ES14.5E3)
 15   FORMAT(A14,10ES14.5E3)
      RETURN
      END SUBROUTINE CrpFn

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Convolution functions
      SUBROUTINE CnvFn
      INTEGER::jt
      REAL::ttt,a1,a3
      bg=Ww*k/2.0/ng*(1.0+chi1-kt*kt/k/k-ng*ng)
      OPEN(218,FILE='bg.dat',STATUS='UNKNOWN')
      WRITE(218,5) '#w(Hz)','wlen(nm)','Re_bg','Im_bg'
      DO iw=1,Nw
        WRITE(218,10) w(iw)*aHz,lam/w(iw),bg(iw)
      END DO
      CLOSE(218)
      ALLOCATE(GW(-Ncr:Ncr),GB(-Ncr:Ncr))
      OPEN(220,FILE='Ct.dat',STATUS='UNKNOWN')
      WRITE(220,5) '#t(fs)','G','S','B'
      DO it=-Ncr,Ncr
        ttt=dt*REAL(it)
        a1=0.0
        a3=0.0
        DO iw=1,Nw
          a1=a1+dw*COS(w(iw)*ttt)*Ww(iw)
          a3=a3-dw*AIMAG(EXP(-i*w(iw)*ttt)*bg(iw))
        END DO
        GW(it)=a1/pi                                                    ! GW: conv. func. for transmission window Ww
        GB(it)=a3/pi                                                    ! GB: conv. func. for fiber+gas dispersion
        WRITE(220,10) ttt*afs,GW(it),GB(it)
      END DO
      CLOSE(220)
      WRITE(*,*) 'G(tau), S(tau), B(tau) prepared'
 5    FORMAT(10A14)
 10   FORMAT(10ES14.5E3)
      RETURN
      END SUBROUTINE CnvFn

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c garr: real envelope function
c farr: complex input function to be processed
      FUNCTION Convo(garr,farr)
      INTEGER::jt
      REAL,ALLOCATABLE::garr(:)
      COMPLEX,ALLOCATABLE::Convo(:),farr(:),farr1(:),farr2(:) 
      ALLOCATE(Convo(Nt),farr1(Nt),farr2(Nt))
      farr1=farr
      farr2=(0.0,0.0)
      DO it=1,Nt
        DO jt=it-Ncr,it+Ncr
          IF (jt<1.OR.jt>Nt) CYCLE
          farr2(it)=farr2(it)+dt*garr(it-jt)*farr1(jt)
        END DO
      END DO
      Convo=farr2*Cr
      RETURN
      END FUNCTION Convo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Pulse(dura,ene_J,freq,phase)
      REAL::dura,ene_J,freq,phase,pow_W,amp,devi,intens,aaa
      devi=dura/SQRT(2.0*LOG(2.0))                                      ! pulse duration deviation (atu)
      pow_W=ene_J/(devi*a0s)*SQRT(2.0/pi)                               ! power at pulse peak (W)
      intens=pow_W/(fbrA*acm*acm)
      amp=SQRT(intens/In0cyc)                                           ! field amplitude (atu)
      aaa=0.0
      OPEN(17,FILE='Ein.dat',STATUS='UNKNOWN')
      WRITE(17,5) '#t(fs)','E_in(V/m)'
      DO it=1,Nt
        E(it)=amp*EXP(-(t(it)/devi)**2)*COS(freq*t(it)+phase)
        WRITE(17,10) t(it)*afs,REAL(E(it))*aVm
        aaa=aaa+dt*ABS(E(it))**2
      END DO
      CLOSE(17)
      WRITE(*,15) 'pulse energy (J)',aaa*In0abs*(fbrA*acm*acm)*a0s,
     &            'pk power (W)',pow_W,
     &            'fiber area (cm2)',fbrA*acm*acm,
     &            'ave pk int (W/cm2)',intens,
     &            'centr pk int (W/cm2)',intens/Iave0,
     &            'field amp (atu)',amp
      WRITE(*,*) 'Pump pulse prepared'
 5    FORMAT(20A14)
 10   FORMAT(20ES14.5E3)
 15   FORMAT(20(A26,ES14.5E3,/))
      RETURN
      END SUBROUTINE Pulse

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Fourier transform fw(w) = int{ exp(iwt) f(t) dt }
      FUNCTION FT(arrayt)
      COMPLEX,ALLOCATABLE::FT(:),arrayt(:)
      ALLOCATE(FT(Nw))
      DO iw=1,Nw
        FT(iw)=0.5*dt*(ewt(iw,1)*arrayt(1)+ewt(iw,Nt)*arrayt(Nt))
        DO it=2,Nt-1
          FT(iw)=FT(iw)+dt*ewt(iw,it)*arrayt(it)
        END DO
      END DO
      FT=cFT*FT
      RETURN
      END FUNCTION FT

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Inverse FT f(t) = int{ exp(-iwt) fw(w) dw }
      FUNCTION IFT(arrayw)
      COMPLEX,ALLOCATABLE::IFT(:),arrayw(:)
      ALLOCATE(IFT(Nt))
      DO it=1,Nt
        IFT(it)=0.5*dw*( CONJG(ewt(iw,1))*arrayw(1)
     &         +CONJG(ewt(iw,Nt))*arrayw(Nw))
        DO iw=2,Nw-1
          IFT(it)=IFT(it)+dw*CONJG(ewt(iw,it))*arrayw(iw)
        END DO
      END DO
      IFT=cFT*IFT
      RETURN
      END FUNCTION IFT

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Propagation in the space
      SUBROUTINE MacroStep
      REAL::Amp,dz0
      COMPLEX,ALLOCATABLE::E00(:),E1(:),E2(:),DelE(:)
      ALLOCATE(E00(Nt),E1(Nt),E2(Nt),DelE(Nt))
      E00=E
      Amp=SQRT(SUM(ABS(E00)**2))
 10   E=E00
      CALL MacroRK(dz)
      E1=E
      E=E00
      CALL MacroRK(dz/2.0)
      CALL MacroRK(dz/2.0)
      E2=E
      DelE=E2-E1
      tolz=SQRT(SUM(ABS(DelE)**2))/Amp
      dz0=dz*ABS(tolz0/tolz)**0.2
      IF (dzl==dzh) THEN ! pass
        dznext=dz
      ELSE IF (dz<dzl) THEN                                             ! dz too low
        dz=dzl
        GO TO 10
      ELSE IF (dz>dzh) THEN                                             ! dz too high
        dz=dzh        
        GO TO 10
      ELSE IF (tolz>tolzlim) THEN                                      ! error too high
        IF (dz==dzl) THEN                                               ! lower limit => keep dz !!! pass
          dznext=dz
        ELSE IF (dz0>=dzl) THEN                                        ! decrease dz to goal and recalculate
          dz=dz0
          GO TO 10
        ELSE IF (dz0<dzl) THEN                                         ! decrease dz to dzl
          dz=dzl
          GO TO 10
        END IF
      ELSE IF (tolz<tolzlim) THEN                                      ! error within limit !!! pass
        IF (dz==dzh) THEN                                               ! upper limit => keep dz
          dznext=dz
        ELSE IF (dz0<=dzh) THEN                                        ! increase next dz to goal
          dznext=dz0
        ELSE IF (dz0>dzh) THEN                                         ! increase next dz to dzh
          dznext=dzh
        END IF
      END IF
      E=E2+DelE/15.0
      RETURN
      END SUBROUTINE MacroStep
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 4th order Runge-Kutta process in z
      SUBROUTINE MacroRK(dzRK)
      REAL::dzRK
      COMPLEX,ALLOCATABLE::k1(:),k2(:),k3(:),k4(:),kk(:),Ec(:)
      ALLOCATE(k1(Nt),k2(Nt),k3(Nt),k4(Nt),kk(Nt),Ec(Nt))
      Ec=E
      CALL PPG(k1)
      E=Ec+0.5*dzRK*k1
      CALL PPG(k2)
      E=Ec+0.5*dzRK*k2
      CALL PPG(k3)
      E=Ec+dzRK*k3
      CALL PPG(k4)
      kk=(k1+2.0*k2+2.0*k3+k4)/6.0
      E=Ec+dzRK*kk
      RETURN
      END SUBROUTINE MacroRK

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Propagator: calculate k for known E
      SUBROUTINE PPG(arr)
      COMPLEX,ALLOCATABLE::arr(:)
      arr=Convo(GB,E)+coeff2*E*E*Deriv(E)
      RETURN
      END SUBROUTINE PPG

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION Deriv(func)
      COMPLEX,ALLOCATABLE::Deriv(:),func(:)
      ALLOCATE(Deriv(Nt))
      Deriv(1)=(func(2)-func(1))/dt
      FORALL (it=2:Nt-1)
        Deriv(it)=(func(it+1)-func(it-1))/2.0/dt
      END FORALL
      Deriv(Nt)=(func(Nt)-func(Nt-1))/dt
      RETURN
      END FUNCTION Deriv

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Macro
      REAL::U,I0,yy,TE,TR,Tstep
      REAL,ALLOCATABLE::pco(:),F(:),FF(:),EE(:),EwEw(:)
      COMPLEX,ALLOCATABLE::Ew(:)
      ALLOCATE(pco(Nt),F(Nt),FF(Nt),EE(Nt),EwEw(Nw),Ew(Nw))
      CALL Pulse(Td,U0J,w0,CEP)
      OPEN(75,FILE='delay_mat.dat',STATUS='UNKNOWN')
      OPEN(76,FILE='wlen_mat.dat',STATUS='UNKNOWN')
      OPEN(77,FILE='freq_mat.dat',STATUS='UNKNOWN')
      DO it=1,Nt,Lt
        WRITE(75,10) t(it)*afs
      END DO
      DO iw=1,Nw,Lw
        WRITE(76,10) lam/w(iw)
      END DO
      DO iw=1,Nw,Lw
        WRITE(77,10) w(iw)*aPHz
      END DO
      CLOSE(75)
      CLOSE(76)
      CLOSE(77)
      WRITE(*,*) 'propagtion starts'
      OPEN(51,FILE='Et.dat',STATUS='UNKNOWN')
      OPEN(52,FILE='It.dat',STATUS='UNKNOWN')
      OPEN(53,FILE='If.dat',STATUS='UNKNOWN')
      OPEN(54,FILE='Il.dat',STATUS='UNKNOWN')
      OPEN(62,FILE='U.dat',STATUS='UNKNOWN')
      OPEN(65,FILE='dz.dat',STATUS='UNKNOWN')
      OPEN(71,FILE='Et_mat.dat',STATUS='UNKNOWN')
      OPEN(72,FILE='Ft_mat.dat',STATUS='UNKNOWN')
      OPEN(43,FILE='If_mat.dat',STATUS='UNKNOWN')
      OPEN(44,FILE='Il_mat.dat',STATUS='UNKNOWN')
      OPEN(74,FILE='z_mat.dat',STATUS='UNKNOWN')
      WRITE(51,5) '#z(cm)','t(fs)','E','F'
      WRITE(52,5) '#z(cm)','t(fs)','I_E','I_F'
      WRITE(53,5) '#z(cm)','freq(PHz)','I_E'
      WRITE(54,5) '#z(cm)','wlen(nm)','I_E'
      WRITE(62,5) '#z(cm)','U(uJ)','I0(W/cm2)'
      WRITE(65,5) '#cpu(min)','z(cm)','dz(um)','tolz'
      iz=1
      z=0.0
      dz=dzh
      dznext=dz
      DO ! z loop
c pulse energy
        yy=0.0
        EE=ABS(E)**2  ! instantaneous intensity
        DO it=1,Nt
          yy=yy+dt*EE(it)
        END DO
        U=yy*In0abs*(fbrA*acm*acm)*a0s
        I0=MAXVAL(EE)*In0abs/2.0/Iave0
        WRITE(62,10) z*acm,U*1E6,I0
c E and P report
        IF (MOD(iz-1,Lz)==0) THEN
          Ew=FT(E)
          EwEw=ABS(Ew)**2
          F=2.0*ABS(CONJG(e0t)*IFT(Ew)) ! profile defined by max field
          FF=F*F
          DO it=1,Nt,Lt
            WRITE(51,10) z*acm,t(it)*afs,REAL(E(it)),F(it)
            WRITE(52,10) z*acm,t(it)*afs,EE(it),FF(it)
          END DO
          WRITE(51,*)
          WRITE(52,*)
          DO iw=1,Nw,Lw
            WRITE(53,10) z*acm,w(iw)*aPHz,EwEw(iw)
            WRITE(54,10) z*acm,lam/w(iw),w(iw)**2/2.0/pi/c0*EwEw(iw)
          END DO
          WRITE(53,*)
          WRITE(54,*)
          WRITE(71,50) (REAL(E(it)),it=1,Nt,Lt)
          WRITE(72,50) (F(it),it=1,Nt,Lt)
          !write(*,*) Nw,Lw,Nt,Lt
          WRITE(43,50) (EwEw(iw),iw=1,Nw,Lw)
          WRITE(44,50) (w(iw)**2/2.0/pi/c0*EwEw(iw),iw=1,Nw,Lw)
          WRITE(74,10) z*acm
        END IF
c end of propagation
        IF (z>=fbrL) EXIT
c propagation
        dz=dznext
        CALL MacroStep
        CALL CPU_TIME(cput1)
        TE=(cput1-cput0)/60.0
        WRITE(65,10) TE,z*acm,dz*aum,tolz
        IF (MOD(iz-1,Lz)==0) THEN
          TR=TE*(fbrL/z-1.0)
          WRITE(*,15) title,'z(m)=',z*am,'dz(m)=',dz*am,
     &                'U(J)=',U,'I0=',I0,
     &                'ET(min)=',TE,'RT(min)=',TR
        END IF
        iz=iz+1
        z=z+dz
      END DO ! z loop
      CLOSE(51)
      CLOSE(53)
      CLOSE(62)
      CLOSE(65)
      CLOSE(66)
      CLOSE(71)
      CLOSE(72)
      CLOSE(43)
      CLOSE(44)
      CLOSE(74)
 5    FORMAT(20A14)
 10   FORMAT(20ES14.5E3)
 12   FORMAT(F14.5,10ES14.5E3)
 15   FORMAT(A15,4(A8,ES11.3,2X),2(A10,F8.2))
 21   FORMAT(ES14.5E3,I14,5ES14.5E3)
 50   FORMAT(5000ES14.5E3)
      RETURN
      END SUBROUTINE Macro


      END MODULE Global

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM Main
      USE Global
      IMPLICIT NONE
      CALL CPU_TIME(cput0)
      pi=ACOS(-1.0)
      eps=1.0/4.0/pi
      CALL Input
      CALL Dimen
      CALL Macro
      CALL CPU_TIME(cput1)
      OPEN(500,FILE='cput.dat',STATUS='UNKNOWN')
      WRITE(500,*) 'CPU time (min) =',(cput1-cput0)/60.0
      CLOSE(500)
      STOP
      END PROGRAM Main

