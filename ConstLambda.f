      program constlambda

c***********************************************************************
c   Isothermes Modell mit konstantem Ausdehnungsparameter
c   Lambda := 4 pi e ne l / B.
c   Einlesen der Parameter Sichtwinkel theta     in grad
c                          Magnetfeld  bmg       in MG
c                          Temperatur  tkev      in keV
c                          log(Lambda) loglambda ohne Einheiten
c                          Geometrie:  0 kein 1/cos Faktor fuer Lambda
c                                      1 mit  1/cos Faktor fuer Lambda
c   aus file 'ConstLambda_Ein'.
c   Ausgegeben wird Intensitaet_nue in erg/s/cm**2/sr/Hz
c                                   fuer ordentlichen       Strahl i_0
c                                   fuer ausserordentlichen Strahl i_1
c                   Gesamtintensitaet_nue i = i_0 + i_1
c   gegen           nue in Hz
c   auf Standardausgabe.
c***********************************************************************

      implicit none

      integer nnue

      parameter(nnue=10000)

      integer i,geomfak
      real*8 sigma,pi,h,c,bk,xme,e,xmp,ryd,ab
      real*8 thetag,bmg,tkev,loglambda,t,omc
      real*8 omega(nnue),phi0(nnue),phi1(nnue),bb(nnue)
      real*8 tau0,tau1,i0,i1
      
      common /konsta/ sigma,pi,h,c,bk,xme,e,xmp,ryd,ab
      common /om/ omega

      call kondat

c   thetag    Sichtwinkel/grad
c   bmg       Magnetfeld/MG
c   tkev      Temperatur/keV
c   loglambda log(Lambda)
c   geomfak   Geometriefaktor
      
      open(10,file='ConstLambda_Ein',status='old')
      read(10,99) thetag
      read(10,99) bmg
      read(10,99) tkev
      read(10,99) loglambda
      read(10,100) geomfak
      close(10)

      omc = e*bmg/xme/c *1d6

      if (geomfak .eq. 1) loglambda = loglambda -
     &                                log10(dcos(thetag*pi/1.8d2))

      call cyccalc(thetag,tkev,phi0,phi1)

      t=tkev/8.616975442d-8
      call planck(omc,t,bb)

      write(*,'(a,f7.3)') "#Sichtwinkel      [grad] :",thetag
      write(*,'(a,f7.3)') "#Magnetfeld       [MG]   :",bmg
      write(*,'(a,f7.3)') "#Temperatur       [keV]  :",tkev
      write(*,'(a,f7.3)') "#log(Lambda)      [1]    :",loglambda
      write(*,'(a,i7)')   "#Geometrie=0 oder 1      :",geomfak
      write(*,'(a)')      "#   nue/Hz       i_0         i_1          i"
      
      do 10 i=1,nnue
        tau0  = 1d1**loglambda * phi0(i)
        if (tau0 .ge. 1d2) then
          i0 = bb(i)
        else
          i0 = bb(i) * (1d0 - 1d0/dexp(tau0))
        endif
        
        tau1  = 1d1**loglambda * phi1(i)
        if (tau1 .ge. 1d2) then
          i1 = bb(i)
        else
          i1 = bb(i) * (1d0 - 1d0/dexp(tau1))
        endif

c     intensitaet_nue in erg/s/cm**2/sr/Hz gegen nue/Hz
        if(i0.eq.0d0) i0=1d-90
        if(i1.eq.0d0) i1=1d-90
10      write(*,'(5(1pd12.4))') omega(i)*omc/2d0/pi,i0,i1,i0+i1

99    format(t27,f12.5)
100   format(t27,i3)

      end

c***********************************************************************
c***********************************************************************
	subroutine cyccalc(thetag,tkev,phi0,phi1)

c***********************************************************************
c   Programm zur Berechnung von Zyklotronabsorptionskoeffizienten
c   des ordentlichen (0) und ausserordentlichen (1) Strahls fuer eine 
c   isotherme Schicht der Elektronendichte xne mit einem konstanten 
c   Sichtwinkel zum konstanten Magnetfeld durch die Schicht.
c   Berechnung nach der abgewandelten Gl(9) in Thompson & Cawthorne('87)
c   MNRAS 224,425. Aenderungen unter Beruecksichtigung der Gln in 
c   Chanmugam & Dulk('81) ApJ 244,569. 
c   Die Parameter werden aus dem file CCalc_Ein gelesen.
c   Die AbsKoeff werden in die files CYCPHI0 bzw CYCPHI1 geschrieben.
c   Dabei ist die Abszisse omega/omegacyklotron (im folgenden 'omc'),
c   die Ordinate sind die AbsKoeff in Einheiten omegaplasma^2/omc/c, das
c   entspricht den bei TC('87) in Gl(9) benutzten PHI's.
c***********************************************************************

	implicit real*8 (a-h,o-z),integer (i-n)

c***********************************************************************
c   nnue : Anzahl der Frequenzstuetzstellen
c   n : Anzahl der Tiefenstuetzstellen
c***********************************************************************

	parameter(nnue=10000)
        parameter(n=50)

        dimension omega(nnue)
        dimension phi0(nnue),phi1(nnue)
        
        common /konsta/ sigma,pi,h,c,bk,xme,e,xmp,ryd,ab
c        common /wz/ xmwz,rwz,gwz,b,xmpunkt,omc,xlambda,xl,tmax
c        common /beob/ thetag,theta,x
        common /om/ omega
	common /ifreq/ nbmin,nbmax,npoi

	external sfint0,sfint1

c        call kondat
	call setdat0
	call setdat1
	call moddat

c***********************************************************************
c   theta,b,tkev,xne  sind Beob.winkel,Magnetfeld,Temperatur (in keV)
c   und Dichte der Absorptionskoeffizienten.
c   nbmin und nbmax sind die Nummern der kleinsten bzw. groessten
c   aufgeloesten Harmonischen.
c***********************************************************************

c	open(1,file='CCalc_Ein',status='old')
c	read(1,1000) thetag
c	read(1,1000) b
c	read(1,1000) tkev
c	read(1,1000) xne
c	read(1,1001) nbmin
c	read(1,1001) nbmax
c	close(1)
	
        b=1d7               !werte egal wg normierung
        xne=1.6579d15       !werte egal wg normierung
        nbmin=1
        nbmax=25
	t=tkev/8.616975442d-8
	theta=thetag*pi/1.8d2
        omc=e*b/xme/c

	call aequisetomega(omega)
c	call setomega(omega)

 	call calcyc(t,theta,b,xne,omega,sfint0,phi0,0)
 	
 	call calcyc(t,theta,b,xne,omega,sfint1,phi1,1)
	
c	open(10,file='CYCPHI0',status='new')
c	open(20,file='CYCPHI1',status='new')

c	write(10,'(a)') '# ordentl. Absorptionskoeffizient fuer'
c	write(10,999)   '#  Theta/grad=',thetag,', B/G=',b,', Temp/K=',t
c        write(10,999)   '#  Temp/keV=',tkev,', Dichte/g/cm^3=',xne*xme
c	write(20,'(a)') '# ausserordentl. Absorptionskoeffizient fuer'
c	write(20,999)   '#  Theta/grad=',thetag,', B/G=',b,', Temp/K=',t
c        write(20,999)   '#  Temp/keV=',tkev,', Dichte/g/cm^3=',xne*xme 

c	do 10 i=1,nnue
c	  if (phi0(i).ne.0d0) then
c	    write(10,1002) omega(i),phi0(i)
c	  else
c	    write(10,1002) omega(i),1d-40
c	  endif
c	  if (phi1(i).ne.0d0) then
c	    write(20,1002) omega(i),phi1(i)
c	  else
cc	    write(20,1002) omega(i),1d-40
c	  endif 
c10	continue	

c 999    format (a,f6.3,4(a,1pg11.4e2))
c1000	format (t19,f12.5)
c1001	format (t19,i3)
c1002	format (f9.5,1pe17.8e3,' ')
	
c	close(10)
c	close(20)

c	write(*,*) 'PROGRAMM BEENDET'
		
	end
c***********************************************************************
c***********************************************************************

        subroutine aequisetomega(omega)

c***********************************************************************
c   Die subroutine aequisetomega erzeugt eine aequidistante
c   Frequenzverteilung in ZyklotronHarmNummer omega/omega_c
c***********************************************************************
        
	implicit real*8(a-h,o-z),integer(i-n)

	common /ifreq/ nbmin,nbmax,npoi
	parameter (nnue=10000)
	dimension omega(nnue)

        fmin=dble(nbmin)-5d-1
        fmax=dble(nbmax)+1d0
	nha=nbmin
	
	do 5 n=1,nbmax-nbmin+1
	  omega(n)=dble(nha)
5	  nha=nha+1     
	
        do 10 nc=1,nnue-(nbmax-nbmin+1)
          omega(nc+nbmax-nbmin+1)=fmin+(fmax-fmin)/dble(nnue)*dble(nc-1)
10	continue

c***********************************************************************
c   Heapsort
c***********************************************************************
	l=nnue/2+1
	ir=nnue
11	continue
	if (l.gt.1) then
	  l=l-1
	  rra=omega(l)
	else
	  rra=omega(ir)
	  omega(ir)=omega(1)
	  ir=ir-1
	  if (ir.eq.1) then
	    omega(1)=rra
	    goto 100
	  endif
	endif
	i=l
	j=l+l
21	if (j.le.ir) then
	  if (j.lt.ir) then
	    if (omega(j).lt.omega(j+1)) j=j+1
	  endif
	  if (rra.lt.omega(j)) then
	    omega(i)=omega(j)
	    i=j
	    j=j+j
	  else
	    j=ir+1
	  endif
	  goto21
	endif
	omega(i)=rra
	goto 11
c***********************************************************************	
	
100	omega(1)=fmin
	omega(nnue)=fmax

	end
c***********************************************************************
c***********************************************************************

        subroutine setomega(omega)

c***********************************************************************
c   Die subroutine setomega erzeugt eine Frequenzverteilung
c    in ZyklotronHarmNummer omega/omega_c
c***********************************************************************
       
	implicit real*8(a-h,o-z),integer(i-n)

	common /ifreq/ nbmin,nbmax,npoi

c***********************************************************************
c   nnue ist die Anzahl der Frequenzpunkte
c***********************************************************************

        parameter (nnue=10000)

c***********************************************************************
c   omega enthaelt die Frequenzstuetzpunkte
c***********************************************************************

        dimension osp(nnue),omega(nnue)

c***********************************************************************
c   nbmin = niedrigste aufgeloeste Harmonische
c   nbmax = hoechste  aufgeloeste Harmonische
c   npoi  = Anzahl Frequenzpunkte pro Linienfluegel. Anzahl
c           der Punkte pro Linie ist damit 2*npoi+1, ( zwei
c           Fluegel plus ein Punkt im Linienzentrum ).
c   nc    = Arbeitsinteger
c   nst   = Anzahl der Frequenzpunkte in den Linien insgesamt
c***********************************************************************

c#####  nbmin=nmin
c#####  nbmax=nmax

        npoi=6
        nc=1
        nst=(nbmax-nbmin+1)*(2*npoi+1)

c***********************************************************************
c   Die Frequenzpunkte werden in folgender Weise gesetzt:
c
c
c   <____________________ cyclotronlinie _________________>
c
c Frequenz- .    .   .  .  .  . . .  .  .  .   .    . punkte
c               
c   <...... npoi ............>   <........ npoi ..........>
c                              1
c                         Linienzentrum
c***********************************************************************

c************* Berechnung von npoi und nst *****************************
c	xnnue=dble(nnue)
c	xnst=8d-1*xnnue/2d0
c	nst=idnint(xnst)-1
c	q=dble(nbmax-nbmin+1)
c	nq=idnint(q/2d0)
c	if (dble(nq).eq.(q/2d0)) then
c	  nst=2*nst
c	else
c	  nst=2*nst-1
c	endif
c	xnst=dble(nst)
c	xnpoi=(xnst/q-1)/2d0
c	npoi=idnint(xnpoi)
c	nst=(nbmax-nbmin+1)*(2*npoi+1)
c***********************************************************************

	nc=1

        do 111 n1=nbmin,nbmax
          omega(nc)=dble(n1)
          nc=nc+1

          do 112 n2=1,npoi
            domega=1d1**(-dble(n2)*4d-1)

            do 113 n3=1,2
              domega=-domega
              omega(nc)=dble(n1)+domega
              nc=nc+1
  113       continue

 112      continue

111     continue 



c***********************************************************************
c   fmin und fmax geben das Frequenzminimum bzw. -maximum zur
c   die Berechnung des Kontinuums an. Die noch zur Verfuegung
c   stehenden (nnue-nst) Frequenzstuetzpunkte werden aequidistant
c   zwischen fmin und fmax verteilt.
c***********************************************************************

        fmin=dble(nbmin)-5d-1
        fmax=dble(nbmax)+1d0

        do 115 nc=nst+1,nnue
115       omega(nc)=fmin+(fmax-fmin)/(nnue-nst-1)*(nc-nst-1)

c***********************************************************************
c   Die gesetzten Frequenzstuetzpunkte werden jetzt in aufstei-
c   gender Reihenfolge geordnet ( d.h.: omega(n+1) > omega(n) ).
c***********************************************************************

        do 120 j=1,nnue
          omin=1d2

          do 121 i=1,nnue
            if (omega(i).lt.omin) then
              omin=omega(i)
              nmin=i
            endif
 121      continue

          osp(j)=omega(nmin)
          omega(nmin)=1d2
120     continue

c	open(1,file='freq',status='new')

        do 122 i=1,nnue
122       omega(i)=osp(i)

c 122   write(1,*) omega(i)
c	close(1)

        end        

c***********************************************************************
c***********************************************************************

      subroutine planck(omc,t,bb)

c***********************************************************************
c*** die subroutine planck berechnet die kirchhoff-planck-intensitaet
c*** (b_nue) pro polarisationsmode zu den nnue frequenzen nue
c*** und zur temperatur t
c***********************************************************************

      implicit real*8(a-h,o-z), integer(i-n)

      parameter (nnue=10000)
      
      real*8 nue,omega(nnue),bb(nnue)

      common /konsta/ sigma,pi,h,c,bk,xme,e,xmp,ryd,ab       
      common /om/ omega

      do 1 i=1,nnue
        nue   = omega(i)*omc/2d0/pi
        bb1   = h*nue**3/c**2
        bb2   = dexp(h*nue/bk/t)-1d0
1       bb(i) = bb1/bb2

      end
      
c***********************************************************************
c***********************************************************************

        subroutine calcyc(t,theta,b,xne,omega,sfint,phi,imod)

c***********************************************************************
c   Die Subroutine CALCYC berechnet die Cyclotronabsorbtionsko-
c   effizienten fuer den ausserordentlichen/ordentlichen Strahl nach 
c   dem Artikel von Chanmugam & Dulk (CD) (Ap.J.(1987) 244, S.569ff) 
c   und den Gesamtabs.koeff nach dem Artikel von Thomson & Cawthorne
c   (TC) (MNRAS (1987) 224, S.425ff)
c***********************************************************************

        implicit real*8(a-h,o-z),integer(i-n)

        parameter (nnue=10000)
        dimension xlphi(nnue,25),phi(nnue),omega(nnue),
     &      xni(nnue),ai(nnue),nnh(25),cycmax(25),ncmax(25)
 
        common /konsta/  sigma,pi,h,c,xkb,xme,e,xmp,ryd,ab
        common /ifreq/ nbmin,nbmax,npoi

        external sfint

	data nnh /8,8,7,7,6,6,5,5,4,4,3,3,2,2,1,1,1,1,1,1,1,1,1,1,1/

c####   nbmin=nmin
c####   nbmax=nmax

        do 10 nc=1,nnue

          do 11 nh=1,nbmax
            xlphi(nc,nh)=0d0
 11       continue

10      continue

        xmue=xme*c*c/xkb/t
           
c***********************************************************************
c   Unter- und Obergrenze des zu berechnenden Integrals sind
c   ug,og. Integrationsvariable ist beta.
c***********************************************************************
	
        if (t.lt.1.5d6) then
          s=dsin(theta)
          co=dcos(theta)
          a=xmue/2d0
          nharm=1
          oharm=1d0
          nh2=2*nharm-2
          nh1=2*nharm-1
          fakul=1d0
  
          do 930 nc=1,nnue
            if (omega(nc).gt.7.5d0) then
              phi(nc)=0d0
              goto 930
            endif
            if (omega(nc).lt.oharm+5d-1) then
	      omp=dsqrt(4d0*pi*xne*e*e/xme)
              omc=e*b/xme/c
              op=omp/omega(nc)/omc
              oc=1d0/omega(nc)
              if (imod.eq.1) then    
                ath=(2d0*(op*op-1d0)*dcos(theta))/(oc*dsin(theta)**2
     &		     +dsqrt(oc*oc*dsin(theta)**4+4d0*(op*op-1d0)**2
     &		     *dcos(theta)**2))
              else
                ath=(2d0*(op*op-1d0)*dcos(theta))/(oc*dsin(theta)**2
     &		     -dsqrt(oc*oc*dsin(theta)**4+4d0*(op*op-1d0)**2
     &		     *dcos(theta)**2))
 	      endif
              got=(ath*(co-1d0/co+oharm/omega(nc))-oharm/omega(nc))**2
     &                                                            /co 
              profil=dexp(-a/(omega(nc)*co)**2*(omega(nc)-oharm)**2)
              phi(nc)=dsqrt(pi)*got*profil*s**nh2*(omega(nc)/2)
     &              **nh1/fakul*a**(1.5d0-oharm)/(1d0+ath)
            else
              nharm=nharm+1
              oharm=oharm+1d0
              nh2=2*nharm-2
              nh1=2*nharm-1
              fakul=fakul*oharm
            endif
930       continue
	
        return
        endif

        do 130 nh=nbmin,nbmax
          cycmax(nh)=0d0

          do 65 io=1,nnue
 65         if (omega(io).gt.dble(nh)) goto 66

          goto 222
 66       continue
          io=io-1
	  
          do 140 ni=1,2
            if (ni.eq.1) nc=io+1
            if (ni.eq.2) nc=io

            do 150 nlauf=1,nnue
              fmem=0d0
              if (ni.eq.1) nc=nc-1
              if (ni.eq.2) nc=nc+1
              if (nc.gt.nnue .or. nc.eq.0) goto 140
              xn=(omega(nc)/dble(nh))**2
              an=1d0+xn*dcos(theta)*dcos(theta)
              ap=xn*dcos(theta)/an
              aq=(xn-1d0)/an
              det=ap*ap-aq
              if (det.le.0d0) goto 150
              og=ap+dsqrt(det)
              ug=ap-dsqrt(det)

c***********************************************************************
c   Beginn der Integration (Simpson)
c***********************************************************************

              ug1=ug
              og1=og
              dxt=1d-5
  570         continue
              xt=ug1+dxt

              call sfint(nh,theta,b,xne,omega(nc),xt,xmue
     &                     ,fint,xni(nc),ai(nc))

              f1=fint
              if (xt.gt.og .and. cycmax(nh).ne.0d0) goto 140
              if (xt.gt.og .and. cycmax(nh).eq.0d0) goto 150
c	      if (nh.le.10 .and. f1.lt.1d-10 .and. fmem.gt.f1) goto150
              if (f1.lt.1d-26 .and. fmem.gt.f1 ) goto 150
              fmem=f1
              if (f1.lt.1d-26 ) then
                dxt=5.*dxt
              else
                goto 580
              endif
              goto 570
  580         continue
              xmem=xt
              simps=0d0

              do 610 intr=1,2
                xt=xmem
                dx=dabs(xmem)*1d-2
  555           continue
                if (intr.eq.1) then
                  if (xt+dx.gt.og) dx=og-xt

                  call sfint(nh,theta,b,xne,omega(nc),xt,xmue
     &                         ,fint,xni(nc),ai(nc))

                  f1=fint

                  call sfint(nh,theta,b,xne,omega(nc),xt+dx/2d0,xmue
     &                         ,fint,xni(nc),ai(nc))

                  f2=fint

                  call sfint(nh,theta,b,xne,omega(nc),xt+dx,xmue
     &                         ,fint,xni(nc),ai(nc))

                  f3=fint
                else
                  if (xt-dx.lt.ug) dx=xt-ug

                  call sfint(nh,theta,b,xne,omega(nc),xt-dx,xmue
     &                         ,fint,xni(nc),ai(nc))

                  f1=fint

                  call sfint(nh,theta,b,xne,omega(nc),xt-dx/2d0,xmue
     &                         ,fint,xni(nc),ai(nc))

                  f2=fint

                  call sfint(nh,theta,b,xne,omega(nc),xt,xmue
     &                         ,fint,xni(nc),ai(nc))

                  f3=fint
                endif
                st1=dabs(-3d0*f1+4d0*f2-f3)/dx
                st2=dabs(f1-2d0*f2+f3)/(dx/2d0)**2
                fmit=(f1+4d0*f2+f3)/6d0
                simps=simps+fmit*dx
                if (fmit*dx.lt.1d-3*simps) goto 666
                if (intr.eq.1) then
                   xt=xt+dx
                else
                   xt=xt-dx
                endif
                if (xt.ge.og .or. xt.le.ug) goto 666
                if (st2.ne.0d0) then
                   dxnew=st1/st2/1d0
                   if (dxnew.lt.2d0*dx) then
                      dx=dxnew
                   else
                      dx=2d0*dx
                   endif
                endif
                goto 555            
 666            continue
 610         continue

             cyc=simps
 600         continue

c***********************************************************************
c   Ende der Integration.
c***********************************************************************

              xlphi(nc,nh)=cyc
              if (cyc.gt.cycmax(nh)) then
                 cycmax(nh)=cyc
                 ncmax(nh)=nc
              endif

c***********************************************************************
c   Interpolation der 'Ausreisser'        
c	   
c           if(                 nh .gt. nbmin                  .and.
c     &                         ni .eq. 1                      .and.
c     &        xlphi(nc+npoi+2,nh) .ge. xlphi(nc+npoi+2,nh-1)  .and.
c     &        (cyc                .lt. xlphi(nc,nh-1)*1d-2  .or.
c     &         cyc                .lt. xlphi(nc+1,nh)*1d-1      )  ) 
c     &                     xlphi(nc,nh)=xlphi(nc+1,nh)/2d0
c***********************************************************************

c***********************************************************************
c   dyn Abbruch der "linken Seite" rel zum Abskoeff der Frequenz in 
c   der vorherigen Harmonischen
c***********************************************************************

              if (                 nh.gt.nbmin                 .and.
     &                             ni.eq.1                     .and.
     &            xlphi(nc+npoi+2,nh).lt.xlphi(nc+npoi+2,nh-1) .and.
     &                             nc.lt.ncmax(nh-1)           .and.
     &            cyc.lt.xlphi(nc,nh-1)*1d-4                  ) goto 140

c***********************************************************************
c   dyn Abbruch der "rechten Seite" rel zum Abskoeff der vorherigen
c   Frequenz der selben Harmonischen
c***********************************************************************

              if ( nc.gt.1                        .and.
     &             ni.eq.2                        .and.
     &            cyc.lt.xlphi(nc-1,nh)*1d1**(-nnh(nh))) goto 140
	      if (nc.eq.1 .and. cyc.lt.1d-33) goto 140

c***********************************************************************
c             if(ni.eq.2              .and.
c     &          nh.gt.10             .and.
c     &         cyc.lt.1d-23               ) goto 140
c***********************************************************************
  150       continue

 140      continue

c***********************************************************************
c   Interpolation der 'Ausreisser'
c*******         'linke Seite' :
  
          do 128 ip=1,nnue
            if (          ip.ne.1              .and.
     &                    ip.ne.nnue-2         .and.
     &                    ip.lt.ncmax(nh)      .and.
     &          xlphi(ip,nh).lt.xlphi(ip-1,nh)      ) then
              xlphi(ip,nh)=dsqrt(xlphi(ip-1,nh)*xlphi(ip+1,nh))
              mist=ip+2
  13          continue                          
              if (xlphi(ip,nh).lt.xlphi(ip-1,nh)) then
                xlphi(ip,nh)=dsqrt(xlphi(ip-1,nh)*xlphi(mist,nh))
                mist=mist+1
                if (mist.eq.nnue) goto 128
                goto 13
              endif             
            endif                          
 128	  continue

c*******         'rechte Seite' :

          do 129  ip=nnue,1,-1
            if (          ip.ne.3              .and.
     &                    ip.ne.nnue           .and.
     &                    ip.gt.ncmax(nh)      .and.
     &          xlphi(ip,nh).lt.xlphi(ip+1,nh)      ) then
              xlphi(ip,nh)=dsqrt(xlphi(ip+1,nh)*xlphi(ip-1,nh))
              mist=ip-2
  23          continue                          
              if (xlphi(ip,nh).lt.xlphi(ip+1,nh)) then
                xlphi(ip,nh)=dsqrt(xlphi(ip+1,nh)*xlphi(mist,nh))
                mist=mist-1
                if (mist.eq.1) goto 129
                goto 23
              endif
            endif
 129	  continue

c***********************************************************************
c   dyn Abb der Berechnung der Abskoeff
c***********************************************************************
          if (        nh.eq.nbmin           .and. 
     &        cycmax(nh).lt.1d-33                ) goto 125
130     continue
	
125     continue
222     continue

        do 300 nc=1,nnue
          phi(nc)=0d0

          do 310 nh=nbmin,nbmax
            phi(nc)=phi(nc)+xlphi(nc,nh)
 310      continue

300     continue

c***********************************************************************
c   Fuer t>=5keV werden die hohen Harmonischen durch eine Formel
c   von Robinson und Melrose genaehert.
c***********************************************************************

c        thetad=theta/3.14*180.
c        tkev=t*.8616975442e-7
c        if(tkev.ge.5.)then

c        do 400 i=ncmax,nnue
c          call robmel (omega(i), theta, xmue,  fip, fim)
c          phi(i)=fip+fim
c 400    continue
c 401    continue

c        dif=log(phi(ncmax-1))-log(phi(ncmax))

c        do 410 i=ncmax,nnue
c        phi(i)=exp(log(phi(i))+dif)
c 410    continue

c        endif

        end

c***********************************************************************
c***********************************************************************
	
        subroutine sfintg(nb,theta,un,ge,omega,beta,xmue,fint,
     &                    brau,cht)

c***********************************************************************
c        fint ist der Integrand aus TC.
c***********************************************************************

        implicit real*8(a-h,o-z),integer(i-n)
     
        pi=3.14d0
        xnb=dble(nb)
        x1=1d0-beta*beta
        x2=1d0-beta*dcos(theta)         
        x3=omega/xnb      
        x4=x1-(x3*x2)**2
        if (x4.le.0d0) then
          fint=0d0
          goto 101
        endif
        betas=dsqrt(x4)
        z=xnb*dsin(theta)*betas/x2
        y1=((dcos(theta)-beta)/dsin(theta)*bes0(nb,z))**2
        y2=(betas*bes1(nb,z))**2
        y3=1d0/x3/x2
        fint=(y1+y2)/xnb/(x2*x3)**4
        if (xmue.lt.5d0)
     &     fint=fint*dexp(-xmue*y3)*pi*xmue**4/2d0/xmbes2(xmue)
        if (xmue.ge.5d0 .and. xmue.lt.2d1)
     &     fint=fint*dexp(xmue*(1d0-y3))*pi*xmue**2/2d0/xmbes2(xmue)
        if (xmue.ge.2d1)
     &     fint=fint*dexp(xmue*(1d0-y3))*xmue**2.5d0/2d0/xmbes2(xmue)
101     continue

        end

c***********************************************************************
c***********************************************************************

        subroutine sfint1(nb,theta,b,xne,omega,beta,xmue,fint,xn1,a1)

c***********************************************************************
c   fint ist der Integrand aus CD.
c***********************************************************************
       
        implicit real*8(a-h,o-z),integer(i-n)
     
        pi=3.141592654d0
	xme=9.11d-28
	e=4.8d-10
	c=2.99792456d10

c***********************************************************************
c   Der Brechungsindex xn1 kann wahlweise berechnet oder auf 1.
c   gesetzt werden. Hierfuer braucht man nur die unerwuenschte
c   Berechnung auszukommentieren.
c***********************************************************************

	omp=dsqrt(4d0*pi*xne*e*e/xme)
        omc=e*b/xme/c
        op=omp/omega/omc
        oc=1d0/omega
c	xn1q=1d0-(2d0*op*op*(op*op-1d0))/(dsqrt(oc**4
c     &      *dsin(theta)**4+4d0*oc*oc*(op*op-1d0)**2
c     &      *dcos(theta)**2)+2d0*(op*op-1d0)+oc*oc
c     &      *dsin(theta)**2)
c	xn1=dsqrt(xn1q)
        xn1=1d0
        a1=(2d0*(op*op-1d0)*dcos(theta))/(oc*dsin(theta)**2
     &	    +dsqrt(oc*oc*dsin(theta)**4+4d0*(op*op-1d0)**2
     &	    *dcos(theta)**2))
        xnb=dble(nb)
        x1=1d0-beta*beta
        x2=1d0-xn1*beta*dcos(theta)         
        x3=omega/xnb      
        x4=x1-(x3*x2)**2
        if (x4.le.0d0) then
          fint=0d0
          goto 101
        endif
        betas=dsqrt(x4)
        z=xn1*xnb*dsin(theta)*betas/x2
        y1=a1*(((dcos(theta)/xn1)-beta)/dsin(theta))*bes0(nb,z)
        y2=betas*bes1(nb,z)
        y3=1d0/x3/x2
	y4=xn1/(1d0+a1*a1)
        fint=y4*((y1-y2)**2)/xnb/(x2*x3)**4
        if (xmue.lt.5d0)
     &     fint=fint*dexp(-xmue*y3)*pi*xmue**4/2d0/xmbes2(xmue)
        if (xmue.ge.5d0 .and. xmue.lt.2d1)
     &     fint=fint*dexp(xmue*(1d0-y3))*pi*xmue**2/2d0/xmbes2(xmue)
        if (xmue.ge.2d1)
     &     fint=fint*dexp(xmue*(1d0-y3))*xmue**2.5d0/2d0/xmbes2(xmue)
101     continue

        end

c***********************************************************************
c***********************************************************************
	
        subroutine sfint0(nb,theta,b,xne,omega,beta,xmue,fint,xn0,a0)

c***********************************************************************
c   fint ist der Integrand aus CD.
c***********************************************************************

        implicit real*8(a-h,o-z),integer(i-n)
     
        pi=3.141592654d0
	xme=9.11d-28
	e=4.8d-10
	c=2.99792456d10

c***********************************************************************
c   Der Brechungsindex xn0 kann wahlweise berechnet oder auf
c   1. gesetzt werden. Hierfuer braucht man nur die ungewuenschte
c   Berechnung auszukommentieren.
c***********************************************************************
		       
	omp=dsqrt(4d0*pi*xne*e*e/xme)
        omc=e*b/xme/c
        op=omp/omega/omc
        oc=1d0/omega

c	xn0=dsqrt(1d0+(2d0*op*op*(op*op-1d0))/(dsqrt(oc**4
c     &     *dsin(theta)**4+4d0*oc*oc*(op*op-1d0)**2
c     &	    *dcos(theta)**2)-2d0*(op*op-1d0)-oc*oc
c     &	    *dsin(theta)**2))
        xn0=1d0
        a0=(2d0*(op*op-1d0)*dcos(theta))/(oc*dsin(theta)**2
     &	    -dsqrt(oc*oc*dsin(theta)**4+4d0*(op*op-1d0)**2
     &	    *dcos(theta)**2))

        xnb=dble(nb)
        x1=1d0-beta*beta
        x2=1d0-xn0*beta*dcos(theta)         
        x3=omega/xnb      
        x4=x1-(x3*x2)**2
        if (x4.le.0d0) then
          fint=0d0
          goto 101
        endif
        betas=dsqrt(x4)
        z=xn0*xnb*dsin(theta)*betas/x2
        y1=a0*(((dcos(theta)/xn0)-beta)/dsin(theta))*bes0(nb,z)
        y2=betas*bes1(nb,z)
        y3=1d0/x3/x2
	y4=xn0/(1d0+a0*a0)
        fint=y4*((y1-y2)**2)/xnb/(x2*x3)**4
        if (xmue.lt.5d0)
     &     fint=fint*dexp(-xmue*y3)*pi*xmue**4/2d0/xmbes2(xmue)
        if (xmue.ge.5d0 .and. xmue.lt.2d1)
     &     fint=fint*dexp(xmue*(1d0-y3))*pi*xmue**2/2d0/xmbes2(xmue)
        if (xmue.ge.2d1)
     &     fint=fint*dexp(xmue*(1d0-y3))*xmue**2.5d0/2d0/xmbes2(xmue)
101     continue

        end

c***********************************************************************
c***********************************************************************

        function bes0(nb,arg)

c***********************************************************************
c   bes0 ist die Besselfunktion(arg) der Ordnung nb.
c***********************************************************************
       
        implicit real*8(a-h,o-z),integer(i-n)

        dimension x(3),y(3)

        common /bes/ bess0(25,250),bess1(25,250),xbe(250)

        if ((nb.le.10) .and. (arg.gt.9.9d0)) print*,'warnung!',arg
	if ((nb.eq.11) .and. (arg.gt.10.9d0)) print*,'warnung!',arg
	if ((nb.eq.12) .and. (arg.gt.11.9d0)) print*,'warnung!',arg
	if ((nb.eq.13) .and. (arg.gt.12.9d0)) print*,'warnung!',arg
	if ((nb.eq.14) .and. (arg.gt.13.9d0)) print*,'warnung!',arg
	if ((nb.eq.15) .and. (arg.gt.14.9d0)) print*,'warnung!',arg
        if ((nb.le.16) .and. (arg.gt.15.9d0)) print*,'warnung!',arg
	if ((nb.eq.17) .and. (arg.gt.16.9d0)) print*,'warnung!',arg
	if ((nb.eq.18) .and. (arg.gt.17.9d0)) print*,'warnung!',arg
	if ((nb.eq.19) .and. (arg.gt.18.9d0)) print*,'warnung!',arg
	if ((nb.eq.20) .and. (arg.gt.19.9d0)) print*,'warnung!',arg
	if ((nb.eq.21) .and. (arg.gt.20.9d0)) print*,'warnung!',arg
	if ((nb.eq.22) .and. (arg.gt.21.9d0)) print*,'warnung!',arg
	if ((nb.eq.23) .and. (arg.gt.22.9d0)) print*,'warnung!',arg
	if ((nb.eq.24) .and. (arg.gt.23.9d0)) print*,'warnung!',arg
	if ((nb.eq.25) .and. (arg.gt.24.9d0)) print*,'warnung!',arg
        if (arg.lt.2d0) then
          no=int(3.3d0*dlog(arg)/dlog(1d1)+4.2d0)+1
          if (no.lt.0) no=0
          fnk=1d0

          do 52 nf=1,nb
52          fnk=fnk*dble(nf)
          bes0=1d0/fnk
          fk=1d0

          do 58 k=1,no
            fk=fk*dble(k)
            fnk=fnk*dble(nb+k)
            bes0=bes0+(-(arg**2)/4d0)**k/fk/fnk
58        continue
	
          bes0=bes0*(arg/2d0)**nb
          goto 55
        endif
        inum=int(1d1*(arg+1d-1))
        x(1)=xbe(inum)
        x(2)=xbe(inum+1)
        x(3)=xbe(inum+2)
        y(1)=bess0(nb,inum)
        y(2)=bess0(nb,inum+1)
        y(3)=bess0(nb,inum+2)
        x0=arg 
        y0=0d0

        do 501 j=1,3
          prod=1d0

          do 500 i=1,3      
            if (j.ne.i) prod=prod*(x(i)-x0)/(x(i)-x(j))
 500      continue

          y0=y0+y(j)*prod
501     continue 

          bes0=y0
55      continue

        end


c***********************************************************************
c***********************************************************************

        function bes1(nb,arg)

c***********************************************************************
c   bes1 ist die erste Ableitung der Besselfkt der Ordnung nb.
c***********************************************************************
       
        implicit real*8(a-h,o-z),integer(i-n)

        dimension x(3),y(3)

        common /bes/ bess0(25,250),bess1(25,250),xbe(250)

        if ((nb.le.10) .and. (arg.gt.9.9d0)) print*,'warnung!',arg
	if ((nb.eq.11) .and. (arg.gt.10.9d0)) print*,'warnung!',arg
	if ((nb.eq.12) .and. (arg.gt.11.9d0)) print*,'warnung!',arg
	if ((nb.eq.13) .and. (arg.gt.12.9d0)) print*,'warnung!',arg
	if ((nb.eq.14) .and. (arg.gt.13.9d0)) print*,'warnung!',arg 
	if ((nb.eq.15) .and. (arg.gt.14.9d0)) print*,'warnung!',arg
        if ((nb.eq.16) .and. (arg.gt.15.9d0)) print*,'warnung!',arg
	if ((nb.eq.17) .and. (arg.gt.16.9d0)) print*,'warnung!',arg
	if ((nb.eq.18) .and. (arg.gt.17.9d0)) print*,'warnung!',arg
	if ((nb.eq.19) .and. (arg.gt.18.9d0)) print*,'warnung!',arg
	if ((nb.eq.20) .and. (arg.gt.19.9d0)) print*,'warnung!',arg 
	if ((nb.eq.21) .and. (arg.gt.20.9d0)) print*,'warnung!',arg
	if ((nb.eq.22) .and. (arg.gt.21.9d0)) print*,'warnung!',arg 
	if ((nb.eq.23) .and. (arg.gt.22.9d0)) print*,'warnung!',arg
	if ((nb.eq.24) .and. (arg.gt.23.9d0)) print*,'warnung!',arg 
	if ((nb.eq.25) .and. (arg.gt.24.9d0)) print*,'warnung!',arg
        if (arg.lt.2d0) then
          no=int(3.3d0*dlog(arg)/dlog(1d1)+4.2d0)+1
          if (no.lt.0) no=0
          fnk=1d0

          do 52 nf=1,nb
52          fnk=fnk*dble(nf)

          bes1=dble(nb)/fnk
          fk=1d0

          do 58 k=1,no
            fk=fk*dble(k)
            fnk=fnk*dble(nb+k)
            bes1=bes1+(-(arg**2)/4d0)**k/fk/fnk*dble(2*k+nb)
58        continue

          bes1=bes1*(arg/2d0)**(nb-1)*5d-1

          goto 55
        endif
        inum=int(1d1*(arg+1d-1))
        x(1)=xbe(inum)
        x(2)=xbe(inum+1)
        x(3)=xbe(inum+2)
        y(1)=bess1(nb,inum)
        y(2)=bess1(nb,inum+1)
        y(3)=bess1(nb,inum+2)
        x0=arg 
        y0=0d0

        do 501 j=1,3
          prod=1d0

          do 500 i=1,3      
            if (j.ne.i) prod=prod*(x(i)-x0)/(x(i)-x(j))
 500      continue

          y0=y0+y(j)*prod
501     continue 

        bes1=y0
55      continue

        end

c***********************************************************************
c***********************************************************************

        function xmbes2(x0)

c***********************************************************************
c   xmbes2 ist die modifizierte Besselfkt. der Ordnung 2.
c***********************************************************************
	
        implicit real*8(a-h,o-z),integer(i-n)

        dimension x(3),y(3)
         
        common /mbes/ xmdbes(3,100),arg(3,100)

        if (x0.gt.5.2d2) then
          y0=4d-1
          goto 1000
        endif
        if (x0.lt.5d0) iarg=1
        if (x0.ge.5d0 .and. x0.lt.2d1) iarg=2
        if (x0.ge.2d1) iarg=3

        do 100 inum=1,100
100       if (x0.lt.arg(iarg,inum)) goto 110
110     continue
         
        if (inum.le.2) then
          x(1)=arg(iarg,inum-1)
          x(2)=arg(iarg,inum)
          x(3)=arg(iarg,inum+1)
          y(1)=xmdbes(iarg,inum-1)
          y(2)=xmdbes(iarg,inum)
          y(3)=xmdbes(iarg,inum+1)
        else
          x(1)=arg(iarg,inum-2)
          x(2)=arg(iarg,inum-1)
          x(3)=arg(iarg,inum)
          y(1)=xmdbes(iarg,inum-2)
          y(2)=xmdbes(iarg,inum-1)
          y(3)=xmdbes(iarg,inum)
        endif
        y0=0d0

        do 501 j=1,3
          prod=1d0

          do 500 i=1,3      
            if(j.ne.i)prod=prod*(x(i)-x0)/(x(i)-x(j))
 500      continue

          y0=y0+y(j)*prod
501     continue 

1000    continue
        xmbes2=y0

        end

c***********************************************************************
c***********************************************************************

        subroutine moddat

        implicit real*8(a-h,o-z),integer(i-n)
         
        common /mbes/ xmdbes(3,100),arg(3,100)

        dimension xmod1(100),xmod2(100),xmod3(100)

        do 100 i=1,100
100       arg(1,i)=dble(i-1)*5d-2
        do 110 i=1,100
110       arg(2,i)=5d0+dble(i-1)*1.5d-1
        do 120 i=1,100
120       arg(3,i)=2d1+dble(i-1)*5d0

        data xmod1   /
     &   2.000d+00, 1.999d+00, 1.995d+00, 1.989d+00, 1.980d+00,
     &   1.970d+00, 1.957d+00, 1.942d+00, 1.926d+00, 1.907d+00,
     &   1.888d+00, 1.866d+00, 1.843d+00, 1.819d+00, 1.794d+00,
     &   1.768d+00, 1.741d+00, 1.713d+00, 1.684d+00, 1.655d+00,
     &   1.625d+00, 1.595d+00, 1.564d+00, 1.533d+00, 1.502d+00,
     &   1.470d+00, 1.439d+00, 1.407d+00, 1.376d+00, 1.345d+00,
     &   1.313d+00, 1.282d+00, 1.251d+00, 1.221d+00, 1.190d+00,
     &   1.160d+00, 1.130d+00, 1.101d+00, 1.072d+00, 1.043d+00,
     &   1.015d+00, 9.873d-01, 9.600d-01, 9.332d-01, 9.068d-01,
     &   8.809d-01, 8.556d-01, 8.307d-01, 8.063d-01, 7.825d-01,
     &   7.591d-01, 7.363d-01, 7.140d-01, 6.922d-01, 6.709d-01,
     &   6.501d-01, 6.298d-01, 6.100d-01, 5.907d-01, 5.719d-01,
     &   5.536d-01, 5.358d-01, 5.184d-01, 5.015d-01, 4.851d-01,
     &   4.691d-01, 4.536d-01, 4.385d-01, 4.238d-01, 4.096d-01,
     &   3.958d-01, 3.823d-01, 3.693d-01, 3.567d-01, 3.444d-01,
     &   3.325d-01, 3.210d-01, 3.099d-01, 2.991d-01, 2.886d-01,
     &   2.784d-01, 2.686d-01, 2.591d-01, 2.499d-01, 2.410d-01,
     &   2.323d-01, 2.240d-01, 2.159d-01, 2.081d-01, 2.006d-01,
     &   1.933d-01, 1.863d-01, 1.794d-01, 1.729d-01, 1.665d-01,
     &   1.604d-01, 1.545d-01, 1.487d-01, 1.432d-01, 1.379d-01 /

        data xmod2  /
     &   7.879d-01, 7.694d-01, 7.519d-01, 7.355d-01, 7.200d-01,
     &   7.053d-01, 6.915d-01, 6.783d-01, 6.657d-01, 6.538d-01,
     &   6.425d-01, 6.316d-01, 6.213d-01, 6.114d-01, 6.019d-01,
     &   5.928d-01, 5.841d-01, 5.757d-01, 5.676d-01, 5.599d-01,
     &   5.524d-01, 5.452d-01, 5.383d-01, 5.316d-01, 5.251d-01,
     &   5.189d-01, 5.128d-01, 5.069d-01, 5.013d-01, 4.958d-01,
     &   4.904d-01, 4.853d-01, 4.803d-01, 4.754d-01, 4.706d-01,
     &   4.660d-01, 4.616d-01, 4.572d-01, 4.529d-01, 4.488d-01,
     &   4.448d-01, 4.409d-01, 4.370d-01, 4.333d-01, 4.296d-01,
     &   4.261d-01, 4.226d-01, 4.192d-01, 4.159d-01, 4.126d-01,
     &   4.095d-01, 4.064d-01, 4.033d-01, 4.004d-01, 3.975d-01,
     &   3.946d-01, 3.918d-01, 3.891d-01, 3.864d-01, 3.838d-01,
     &   3.812d-01, 3.787d-01, 3.762d-01, 3.738d-01, 3.714d-01,
     &   3.690d-01, 3.667d-01, 3.645d-01, 3.622d-01, 3.601d-01,
     &   3.579d-01, 3.558d-01, 3.537d-01, 3.517d-01, 3.497d-01,
     &   3.477d-01, 3.458d-01, 3.439d-01, 3.420d-01, 3.402d-01,
     &   3.383d-01, 3.366d-01, 3.348d-01, 3.331d-01, 3.313d-01,
     &   3.297d-01, 3.280d-01, 3.264d-01, 3.248d-01, 3.232d-01,
     &   3.216d-01, 3.201d-01, 3.185d-01, 3.170d-01, 3.156d-01,
     &   3.141d-01, 3.127d-01, 3.112d-01, 3.098d-01, 3.085d-01 /

        data xmod3  /
     &   4.371d-01, 4.294d-01, 4.242d-01, 4.206d-01, 4.178d-01,
     &   4.157d-01, 4.140d-01, 4.127d-01, 4.115d-01, 4.105d-01,
     &   4.097d-01, 4.090d-01, 4.083d-01, 4.078d-01, 4.073d-01,
     &   4.069d-01, 4.065d-01, 4.061d-01, 4.058d-01, 4.055d-01,
     &   4.052d-01, 4.049d-01, 4.047d-01, 4.045d-01, 4.043d-01,
     &   4.041d-01, 4.039d-01, 4.038d-01, 4.036d-01, 4.035d-01,
     &   4.034d-01, 4.032d-01, 4.031d-01, 4.030d-01, 4.029d-01,
     &   4.028d-01, 4.027d-01, 4.026d-01, 4.025d-01, 4.024d-01,
     &   4.023d-01, 4.023d-01, 4.022d-01, 4.021d-01, 4.021d-01,
     &   4.020d-01, 4.019d-01, 4.019d-01, 4.018d-01, 4.018d-01,
     &   4.017d-01, 4.017d-01, 4.016d-01, 4.016d-01, 4.015d-01,
     &   4.015d-01, 4.014d-01, 4.014d-01, 4.014d-01, 4.013d-01,
     &   4.013d-01, 4.012d-01, 4.012d-01, 4.012d-01, 4.011d-01,
     &   4.011d-01, 4.011d-01, 4.011d-01, 4.010d-01, 4.010d-01,
     &   4.010d-01, 4.009d-01, 4.009d-01, 4.009d-01, 4.009d-01,
     &   4.008d-01, 4.008d-01, 4.008d-01, 4.008d-01, 4.007d-01,
     &   4.007d-01, 4.007d-01, 4.007d-01, 4.007d-01, 4.006d-01,
     &   4.006d-01, 4.006d-01, 4.006d-01, 4.006d-01, 4.006d-01,
     &   4.005d-01, 4.005d-01, 4.005d-01, 4.005d-01, 4.005d-01,
     &   4.005d-01, 4.004d-01, 4.004d-01, 4.004d-01, 4.004d-01 /

        do 200 i=1,100
200       xmdbes(1,i)=xmod1(i)
        do 210 i=1,100
210       xmdbes(2,i)=xmod2(i)
        do 220 i=1,100
220       xmdbes(3,i)=xmod3(i)

        end

c***********************************************************************
c***********************************************************************

	subroutine robmel (s, theta, xmue, fip, fim)

*======================================================================
*
*	cyclotron absorption coefficient
*
*	analytix formula derived by robinson and melrose
*	in the form given by dulk (1985)
*
*	valid for parameter range: theta .gt. 30 deg
*				   tkev  .gt. 5-10 kev
*				   s     .gt. 10
*
*	input : xmue  = m c**2 / k-boltz * tkev (common-block)
*	        s     : dimensionless frequency
*	        theta : polar angle (rad)
*	output:	fip   : dimensionless absorption coefficient for
*		fim   : ordinary and extraordinary modes,
*			the abs. coeff. in both modes can be calculated
*			by multiplication with 4 pi e- Ne / B
*
*	because of nue << nue-plasma  n+- approx. 1. and no value
* 	for the electron density Ne is required
*
        implicit real*8(a-h,o-z),integer(i-n)

        common /konsta/  sigma,pi,h,c,xkb,xme,e,xmp,ryd,ab

	xnpm = 1.
	x = s * (dsin(theta))**2 / xmue
	g0 = sqrt(1.+2.*(s/xmue)*(1.+4.5*x)**(-1./3.))
	b0 = sqrt(1.-1./g0**2)
	xx = 1.-(xnpm*b0*dcos(theta))**2
	b1 = xnpm*b0*dsin(theta)/sqrt(xx)
	s0 = g0*s*xx
	xi0 = 1./sqrt(1.-b1**2)
	a  = dtan(theta)*dsin(theta)/(2.*s)
	ap = -(a+sqrt(1.+a**2))*dcos(theta) / abs(dcos(theta))
	am = -1./ap
	sc = 1.5 * xi0**3
	yy = 1.-(xnpm*b0)**2
	c2p = ap*dcos(theta)*yy
	c2m = am*dcos(theta)*yy
	z1 = 1./xi0
	z = b1*exp(z1)/(1.+z1)
	
	ww1p = (c2p*(1.+0.85*sc/s0)**(-1./3.)+sqrt(yy*xx))**2
	ww1m = (c2m*(1.+0.85*sc/s0)**(-1./3.)+sqrt(yy*xx))**2
	ww2 = exp(-xmue*(g0-1.)) *
     &		z**(2.*s0)*(1.+4.5297*sc/s0)**(1./6.)*xx
	ww3 = 2.67e-9*xmue**2*(1.-15./8./xmue)*g0**1.5*sqrt(g0**2-1.) *
     &		xi0**2*(xi0**2-1.) /
     &		(xnpm*dsin(theta)**3*sqrt(x)*s0**1.5)
	ww4 = xnpm**2*b0**2*xi0*(dsin(theta)**4) / (2.*(s0+sc))

	fip = ww3*(ww1p+ww4*ap**2)*ww2 / (1.+ap**2)
	fip = fip / (4.*pi*e)
	fim = ww3*(ww1m+ww4*am**2)*ww2 / (1.+am**2)
	fim = fim / (4.*pi*e)

	return

	end

c***********************************************************************
c***********************************************************************

        subroutine kondat
  
        implicit real*8(a-h,o-z),integer(i-n)

        common /konsta/  sigma,pi,h,c,bk,xme,e,xmp,ryd,ab
        
        h     = 6.626196d-27
        c     = 2.99792456d+10
        bk    = 1.380622d-16
        sigma = 5.669610d-05
        pi    = 3.141592654d0
        xme   = 9.1094d-28
        e     = 4.8032d-10
        xmp   = 1.6726d-24
        ryd   = 2.17d-11
        ab    = 5.2918d-9

        end

c***********************************************************************
c***********************************************************************

        subroutine setdat0

        implicit real*8(a-h,o-z),integer(i-n)

        common /bes/ bess0(25,250),bess1(25,250),arg(250)

        dimension bes1(250),bes2(250),bes3(250),bes4(250),bes5(250)
        dimension bes6(250),bes7(250),bes8(250),bes9(250),bes10(250)
        dimension bes11(250),bes12(250),bes13(250),bes14(250),bes15(250)
        dimension bes16(250),bes17(250),bes18(250),bes19(250),bes20(250)
        dimension bes21(250),bes22(250),bes23(250),bes24(250),bes25(250)

	do 10 j=1,250

	  do 20 i=1,25
	    bess0(i,j)=0d0
	    bess1(i,j)=0d0
 20	  continue

10	continue

        do 100 i=1,250
100	  arg(i)=dble(i-1)/1d1

        data bes1 /
     &   0.0000d+00,  0.4994d-01,  0.9950d-01,  0.1483d+00,  0.1960d+00, 
     &   0.2423d+00,  0.2867d+00,  0.3290d+00,  0.3688d+00,  0.4059d+00, 
     &   0.4401d+00,  0.4709d+00,  0.4983d+00,  0.5220d+00,  0.5419d+00, 
     &   0.5579d+00,  0.5699d+00,  0.5778d+00,  0.5815d+00,  0.5812d+00, 
     &   0.5767d+00,  0.5683d+00,  0.5560d+00,  0.5399d+00,  0.5202d+00, 
     &   0.4971d+00,  0.4708d+00,  0.4416d+00,  0.4097d+00,  0.3754d+00, 
     &   0.3391d+00,  0.3009d+00,  0.2613d+00,  0.2207d+00,  0.1792d+00, 
     &   0.1374d+00,  0.9547d-01,  0.5383d-01,  0.1282d-01, -0.2724d-01, 
     &  -0.6604d-01, -0.1033d+00, -0.1386d+00, -0.1719d+00, -0.2028d+00, 
     &  -0.2311d+00, -0.2566d+00, -0.2791d+00, -0.2985d+00, -0.3147d+00, 
     &  -0.3276d+00, -0.3371d+00, -0.3432d+00, -0.3460d+00, -0.3453d+00, 
     &  -0.3414d+00, -0.3343d+00, -0.3241d+00, -0.3110d+00, -0.2951d+00, 
     &  -0.2767d+00, -0.2559d+00, -0.2329d+00, -0.2081d+00, -0.1816d+00, 
     &  -0.1538d+00, -0.1250d+00, -0.9534d-01, -0.6522d-01, -0.3490d-01, 
     &  -0.4683d-02,  0.2515d-01,  0.5433d-01,  0.8257d-01,  0.1096d+00, 
     &   0.1352d+00,  0.1592d+00,  0.1813d+00,  0.2014d+00,  0.2192d+00, 
     &   0.2346d+00,  0.2476d+00,  0.2580d+00,  0.2657d+00,  0.2708d+00, 
     &   0.2731d+00,  0.2728d+00,  0.2697d+00,  0.2641d+00,  0.2559d+00, 
     &   0.2453d+00,  0.2324d+00,  0.2174d+00,  0.2004d+00,  0.1816d+00, 
     &   0.1613d+00,  0.1395d+00,  0.1166d+00,  0.9284d-01,  0.6837d-01,
     &   150*0d0/

        data bes2 /
     &   0.0000d+00,  0.1249d-02,  0.4983d-02,  0.1117d-01,  0.1973d-01, 
     &   0.3060d-01,  0.4367d-01,  0.5879d-01,  0.7582d-01,  0.9459d-01, 
     &   0.1149d+00,  0.1366d+00,  0.1593d+00,  0.1830d+00,  0.2074d+00, 
     &   0.2321d+00,  0.2570d+00,  0.2817d+00,  0.3061d+00,  0.3299d+00, 
     &   0.3528d+00,  0.3746d+00,  0.3951d+00,  0.4139d+00,  0.4310d+00, 
     &   0.4461d+00,  0.4590d+00,  0.4696d+00,  0.4777d+00,  0.4832d+00, 
     &   0.4861d+00,  0.4862d+00,  0.4835d+00,  0.4780d+00,  0.4697d+00, 
     &   0.4586d+00,  0.4448d+00,  0.4283d+00,  0.4093d+00,  0.3879d+00, 
     &   0.3641d+00,  0.3383d+00,  0.3105d+00,  0.2811d+00,  0.2501d+00, 
     &   0.2178d+00,  0.1846d+00,  0.1506d+00,  0.1161d+00,  0.8129d-01, 
     &   0.4657d-01,  0.1214d-01, -0.2172d-01, -0.5475d-01, -0.8670d-01, 
     &  -0.1173d+00, -0.1464d+00, -0.1737d+00, -0.1990d+00, -0.2221d+00, 
     &  -0.2429d+00, -0.2612d+00, -0.2769d+00, -0.2899d+00, -0.3001d+00, 
     &  -0.3074d+00, -0.3119d+00, -0.3135d+00, -0.3123d+00, -0.3082d+00, 
     &  -0.3014d+00, -0.2920d+00, -0.2800d+00, -0.2656d+00, -0.2490d+00, 
     &  -0.2303d+00, -0.2097d+00, -0.1875d+00, -0.1638d+00, -0.1389d+00, 
     &  -0.1130d+00, -0.8638d-01, -0.5929d-01, -0.3197d-01, -0.4684d-02, 
     &   0.2232d-01,  0.4881d-01,  0.7453d-01,  0.9925d-01,  0.1228d+00, 
     &   0.1448d+00,  0.1653d+00,  0.1840d+00,  0.2008d+00,  0.2154d+00, 
     &   0.2279d+00,  0.2380d+00,  0.2458d+00,  0.2512d+00,  0.2542d+00,
     &   150*0d0/

        data bes3 /
     &   0.0000d+00,  0.2082d-04,  0.1663d-03,  0.5593d-03,  0.1320d-02, 
     &   0.2564d-02,  0.4400d-02,  0.6930d-02,  0.1025d-01,  0.1443d-01, 
     &   0.1956d-01,  0.2569d-01,  0.3287d-01,  0.4114d-01,  0.5050d-01, 
     &   0.6096d-01,  0.7252d-01,  0.8515d-01,  0.9880d-01,  0.1134d+00, 
     &   0.1289d+00,  0.1453d+00,  0.1623d+00,  0.1800d+00,  0.1981d+00, 
     &   0.2166d+00,  0.2353d+00,  0.2540d+00,  0.2727d+00,  0.2911d+00, 
     &   0.3091d+00,  0.3264d+00,  0.3431d+00,  0.3588d+00,  0.3734d+00, 
     &   0.3868d+00,  0.3988d+00,  0.4092d+00,  0.4180d+00,  0.4250d+00, 
     &   0.4302d+00,  0.4333d+00,  0.4344d+00,  0.4333d+00,  0.4301d+00, 
     &   0.4247d+00,  0.4171d+00,  0.4072d+00,  0.3952d+00,  0.3811d+00, 
     &   0.3648d+00,  0.3466d+00,  0.3265d+00,  0.3046d+00,  0.2811d+00, 
     &   0.2561d+00,  0.2298d+00,  0.2023d+00,  0.1738d+00,  0.1446d+00, 
     &   0.1148d+00,  0.8460d-01,  0.5428d-01,  0.2404d-01, -0.5908d-02, 
     &  -0.3535d-01, -0.6406d-01, -0.9184d-01, -0.1185d+00, -0.1438d+00, 
     &  -0.1676d+00, -0.1896d+00, -0.2099d+00, -0.2281d+00, -0.2442d+00, 
     &  -0.2581d+00, -0.2696d+00, -0.2787d+00, -0.2853d+00, -0.2895d+00, 
     &  -0.2911d+00, -0.2903d+00, -0.2869d+00, -0.2811d+00, -0.2730d+00, 
     &  -0.2626d+00, -0.2501d+00, -0.2355d+00, -0.2190d+00, -0.2007d+00, 
     &  -0.1809d+00, -0.1598d+00, -0.1374d+00, -0.1141d+00, -0.8997d-01, 
     &  -0.6532d-01, -0.4034d-01, -0.1526d-01,  0.9700d-02,  0.3432d-01,
     &   150*0d0/

        data bes4 /
     &   0.0000d+00,  0.2603d-06,  0.4158d-05,  0.2100d-04,  0.6614d-04, 
     &   0.1607d-03,  0.3315d-03,  0.6101d-03,  0.1033d-02,  0.1641d-02, 
     &   0.2477d-02,  0.3588d-02,  0.5023d-02,  0.6831d-02,  0.9063d-02, 
     &   0.1177d-01,  0.1500d-01,  0.1879d-01,  0.2320d-01,  0.2825d-01, 
     &   0.3400d-01,  0.4045d-01,  0.4765d-01,  0.5560d-01,  0.6431d-01, 
     &   0.7378d-01,  0.8401d-01,  0.9498d-01,  0.1067d+00,  0.1190d+00, 
     &   0.1320d+00,  0.1456d+00,  0.1597d+00,  0.1743d+00,  0.1892d+00, 
     &   0.2044d+00,  0.2198d+00,  0.2353d+00,  0.2507d+00,  0.2661d+00, 
     &   0.2811d+00,  0.2958d+00,  0.3100d+00,  0.3236d+00,  0.3365d+00, 
     &   0.3484d+00,  0.3594d+00,  0.3693d+00,  0.3780d+00,  0.3853d+00, 
     &   0.3912d+00,  0.3956d+00,  0.3985d+00,  0.3996d+00,  0.3991d+00, 
     &   0.3967d+00,  0.3926d+00,  0.3866d+00,  0.3788d+00,  0.3691d+00, 
     &   0.3576d+00,  0.3444d+00,  0.3294d+00,  0.3128d+00,  0.2945d+00, 
     &   0.2748d+00,  0.2537d+00,  0.2313d+00,  0.2077d+00,  0.1832d+00, 
     &   0.1578d+00,  0.1317d+00,  0.1051d+00,  0.7811d-01,  0.5097d-01, 
     &   0.2382d-01, -0.3126d-02, -0.2970d-01, -0.5572d-01, -0.8100d-01, 
     &  -0.1054d+00, -0.1286d+00, -0.1507d+00, -0.1713d+00, -0.1903d+00, 
     &  -0.2077d+00, -0.2233d+00, -0.2369d+00, -0.2485d+00, -0.2581d+00, 
     &  -0.2655d+00, -0.2707d+00, -0.2736d+00, -0.2743d+00, -0.2728d+00, 
     &  -0.2691d+00, -0.2633d+00, -0.2553d+00, -0.2453d+00, -0.2334d+00,
     &   150*0d0/

        data bes5 /
     &   0.0000d+00,  0.2603d-08,  0.8319d-07,  0.6304d-06,  0.2649d-05, 
     &   0.8054d-05,  0.1995d-04,  0.4288d-04,  0.8308d-04,  0.1487d-03, 
     &   0.2498d-03,  0.3987d-03,  0.6101d-03,  0.9008d-03,  0.1290d-02, 
     &   0.1799d-02,  0.2452d-02,  0.3275d-02,  0.4294d-02,  0.5538d-02, 
     &   0.7040d-02,  0.8828d-02,  0.1094d-01,  0.1340d-01,  0.1624d-01, 
     &   0.1950d-01,  0.2321d-01,  0.2739d-01,  0.3207d-01,  0.3728d-01, 
     &   0.4303d-01,  0.4934d-01,  0.5624d-01,  0.6372d-01,  0.7179d-01, 
     &   0.8044d-01,  0.8968d-01,  0.9949d-01,  0.1098d+00,  0.1207d+00, 
     &   0.1321d+00,  0.1439d+00,  0.1561d+00,  0.1687d+00,  0.1816d+00, 
     &   0.1947d+00,  0.2080d+00,  0.2214d+00,  0.2347d+00,  0.2480d+00, 
     &   0.2611d+00,  0.2740d+00,  0.2865d+00,  0.2986d+00,  0.3101d+00, 
     &   0.3209d+00,  0.3310d+00,  0.3403d+00,  0.3486d+00,  0.3559d+00, 
     &   0.3621d+00,  0.3671d+00,  0.3708d+00,  0.3731d+00,  0.3741d+00, 
     &   0.3736d+00,  0.3716d+00,  0.3680d+00,  0.3629d+00,  0.3562d+00, 
     &   0.3479d+00,  0.3380d+00,  0.3266d+00,  0.3137d+00,  0.2993d+00, 
     &   0.2835d+00,  0.2663d+00,  0.2478d+00,  0.2282d+00,  0.2075d+00, 
     &   0.1858d+00,  0.1632d+00,  0.1399d+00,  0.1161d+00,  0.9175d-01, 
     &   0.6713d-01,  0.4237d-01,  0.1761d-01, -0.6987d-02, -0.3125d-01, 
     &  -0.5504d-01, -0.7818d-01, -0.1005d+00, -0.1219d+00, -0.1422d+00, 
     &  -0.1613d+00, -0.1790d+00, -0.1953d+00, -0.2099d+00, -0.2229d+00,
     &   150*0d0/

        data bes6 /
     &   0.0000d+00,  0.2169d-10,  0.1387d-08,  0.1577d-07,  0.8838d-07, 
     &   0.3361d-06,  0.9996d-06,  0.2509d-05,  0.5560d-05,  0.1120d-04, 
     &   0.2094d-04,  0.3682d-04,  0.6154d-04,  0.9859d-04,  0.1523d-03, 
     &   0.2280d-03,  0.3321d-03,  0.4721d-03,  0.6569d-03,  0.8965d-03, 
     &   0.1202d-02,  0.1587d-02,  0.2066d-02,  0.2653d-02,  0.3367d-02, 
     &   0.4225d-02,  0.5246d-02,  0.6452d-02,  0.7863d-02,  0.9503d-02, 
     &   0.1139d-01,  0.1356d-01,  0.1602d-01,  0.1881d-01,  0.2193d-01, 
     &   0.2543d-01,  0.2931d-01,  0.3360d-01,  0.3832d-01,  0.4347d-01, 
     &   0.4909d-01,  0.5517d-01,  0.6172d-01,  0.6876d-01,  0.7628d-01, 
     &   0.8428d-01,  0.9275d-01,  0.1017d+00,  0.1111d+00,  0.1209d+00, 
     &   0.1310d+00,  0.1416d+00,  0.1525d+00,  0.1637d+00,  0.1751d+00, 
     &   0.1868d+00,  0.1986d+00,  0.2104d+00,  0.2223d+00,  0.2341d+00, 
     &   0.2458d+00,  0.2574d+00,  0.2686d+00,  0.2795d+00,  0.2900d+00, 
     &   0.2999d+00,  0.3093d+00,  0.3180d+00,  0.3259d+00,  0.3330d+00, 
     &   0.3392d+00,  0.3444d+00,  0.3486d+00,  0.3516d+00,  0.3535d+00, 
     &   0.3541d+00,  0.3535d+00,  0.3516d+00,  0.3483d+00,  0.3436d+00, 
     &   0.3376d+00,  0.3301d+00,  0.3213d+00,  0.3111d+00,  0.2996d+00, 
     &   0.2867d+00,  0.2725d+00,  0.2571d+00,  0.2406d+00,  0.2230d+00, 
     &   0.2043d+00,  0.1847d+00,  0.1644d+00,  0.1432d+00,  0.1215d+00, 
     &   0.9932d-01,  0.7675d-01,  0.5396d-01,  0.3107d-01,  0.8216d-02,
     &   150*0d0/

        data bes7 /
     &   0.0000d+00,  0.1550d-12,  0.1982d-10,  0.3381d-09,  0.2527d-08, 
     &   0.1202d-07,  0.4291d-07,  0.1257d-06,  0.3186d-06,  0.7229d-06, 
     &   0.1502d-05,  0.2908d-05,  0.5309d-05,  0.9225d-05,  0.1537d-04, 
     &   0.2468d-04,  0.3840d-04,  0.5809d-04,  0.8571d-04,  0.1237d-03, 
     &   0.1749d-03,  0.2430d-03,  0.3319d-03,  0.4467d-03,  0.5927d-03, 
     &   0.7766d-03,  0.1005d-02,  0.1287d-02,  0.1631d-02,  0.2048d-02, 
     &   0.2547d-02,  0.3142d-02,  0.3845d-02,  0.4669d-02,  0.5630d-02, 
     &   0.6743d-02,  0.8024d-02,  0.9490d-02,  0.1116d-01,  0.1305d-01, 
     &   0.1518d-01,  0.1756d-01,  0.2022d-01,  0.2317d-01,  0.2643d-01, 
     &   0.3002d-01,  0.3395d-01,  0.3824d-01,  0.4290d-01,  0.4794d-01, 
     &   0.5338d-01,  0.5921d-01,  0.6545d-01,  0.7209d-01,  0.7914d-01, 
     &   0.8660d-01,  0.9445d-01,  0.1027d+00,  0.1113d+00,  0.1203d+00, 
     &   0.1296d+00,  0.1392d+00,  0.1491d+00,  0.1592d+00,  0.1696d+00, 
     &   0.1801d+00,  0.1908d+00,  0.2015d+00,  0.2122d+00,  0.2230d+00, 
     &   0.2336d+00,  0.2441d+00,  0.2543d+00,  0.2643d+00,  0.2739d+00, 
     &   0.2832d+00,  0.2919d+00,  0.3001d+00,  0.3076d+00,  0.3145d+00, 
     &   0.3206d+00,  0.3259d+00,  0.3303d+00,  0.3337d+00,  0.3362d+00, 
     &   0.3376d+00,  0.3379d+00,  0.3371d+00,  0.3351d+00,  0.3319d+00, 
     &   0.3275d+00,  0.3218d+00,  0.3149d+00,  0.3068d+00,  0.2974d+00, 
     &   0.2868d+00,  0.2750d+00,  0.2620d+00,  0.2480d+00,  0.2328d+00,
     &   150*0d0/

        data bes8 /
     &   0.0000d+00,  0.9685d-15,  0.2477d-12,  0.6341d-11,  0.6321d-10, 
     &   0.3758d-09,  0.1611d-08,  0.5509d-08,  0.1597d-07,  0.4078d-07, 
     &   0.9422d-07,  0.2008d-06,  0.4002d-06,  0.7540d-06,  0.1354d-05, 
     &   0.2332d-05,  0.3874d-05,  0.6235d-05,  0.9753d-05,  0.1488d-04, 
     &   0.2218d-04,  0.3239d-04,  0.4643d-04,  0.6543d-04,  0.9076d-04, 
     &   0.1241d-03,  0.1674d-03,  0.2230d-03,  0.2937d-03,  0.3826d-03, 
     &   0.4934d-03,  0.6304d-03,  0.7982d-03,  0.1002d-02,  0.1248d-02, 
     &   0.1543d-02,  0.1894d-02,  0.2309d-02,  0.2797d-02,  0.3366d-02, 
     &   0.4029d-02,  0.4794d-02,  0.5674d-02,  0.6681d-02,  0.7827d-02, 
     &   0.9126d-02,  0.1059d-01,  0.1224d-01,  0.1408d-01,  0.1613d-01, 
     &   0.1841d-01,  0.2092d-01,  0.2369d-01,  0.2673d-01,  0.3004d-01, 
     &   0.3366d-01,  0.3758d-01,  0.4182d-01,  0.4638d-01,  0.5128d-01, 
     &   0.5653d-01,  0.6213d-01,  0.6808d-01,  0.7438d-01,  0.8103d-01, 
     &   0.8804d-01,  0.9539d-01,  0.1031d+00,  0.1111d+00,  0.1194d+00, 
     &   0.1280d+00,  0.1368d+00,  0.1459d+00,  0.1553d+00,  0.1648d+00, 
     &   0.1744d+00,  0.1842d+00,  0.1940d+00,  0.2039d+00,  0.2137d+00, 
     &   0.2235d+00,  0.2331d+00,  0.2426d+00,  0.2518d+00,  0.2608d+00, 
     &   0.2694d+00,  0.2775d+00,  0.2853d+00,  0.2925d+00,  0.2991d+00, 
     &   0.3051d+00,  0.3103d+00,  0.3148d+00,  0.3185d+00,  0.3214d+00, 
     &   0.3233d+00,  0.3243d+00,  0.3242d+00,  0.3232d+00,  0.3211d+00,
     &   150*0d0/

        data bes9 /
     &   0.0000d+00,  0.5381d-17,  0.2753d-14,  0.1057d-12,  0.1405d-11, 
     &   0.1045d-10,  0.5375d-10,  0.2145d-09,  0.7109d-09,  0.2043d-08, 
     &   0.5249d-08,  0.1231d-07,  0.2679d-07,  0.5471d-07,  0.1059d-06, 
     &   0.1956d-06,  0.3469d-06,  0.5936d-06,  0.9843d-06,  0.1586d-05, 
     &   0.2492d-05,  0.3827d-05,  0.5753d-05,  0.8487d-05,  0.1230d-04, 
     &   0.1754d-04,  0.2465d-04,  0.3415d-04,  0.4672d-04,  0.6315d-04, 
     &   0.8440d-04,  0.1116d-03,  0.1462d-03,  0.1896d-03,  0.2438d-03, 
     &   0.3109d-03,  0.3934d-03,  0.4940d-03,  0.6160d-03,  0.7628d-03, 
     &   0.9386d-03,  0.1148d-02,  0.1395d-02,  0.1686d-02,  0.2027d-02, 
     &   0.2425d-02,  0.2885d-02,  0.3417d-02,  0.4027d-02,  0.4725d-02, 
     &   0.5520d-02,  0.6422d-02,  0.7441d-02,  0.8588d-02,  0.9873d-02, 
     &   0.1131d-01,  0.1291d-01,  0.1468d-01,  0.1664d-01,  0.1880d-01, 
     &   0.2117d-01,  0.2376d-01,  0.2658d-01,  0.2966d-01,  0.3299d-01, 
     &   0.3659d-01,  0.4047d-01,  0.4463d-01,  0.4909d-01,  0.5385d-01, 
     &   0.5892d-01,  0.6430d-01,  0.6999d-01,  0.7599d-01,  0.8230d-01, 
     &   0.8892d-01,  0.9584d-01,  0.1031d+00,  0.1105d+00,  0.1183d+00, 
     &   0.1263d+00,  0.1346d+00,  0.1430d+00,  0.1517d+00,  0.1605d+00, 
     &   0.1694d+00,  0.1785d+00,  0.1876d+00,  0.1967d+00,  0.2058d+00, 
     &   0.2149d+00,  0.2238d+00,  0.2327d+00,  0.2413d+00,  0.2496d+00, 
     &   0.2577d+00,  0.2655d+00,  0.2728d+00,  0.2797d+00,  0.2860d+00,
     &   150*0d0/ 

        data bes10 /
     &   0.0000d+00,  0.2691d-19,  0.2753d-16,  0.1586d-14,  0.2812d-13, 
     &   0.2613d-12,  0.1614d-11,  0.7518d-11,  0.2848d-10,  0.9212d-10, 
     &   0.2631d-09,  0.6791d-09,  0.1613d-08,  0.3570d-08,  0.7444d-08, 
     &   0.1474d-07,  0.2791d-07,  0.5080d-07,  0.8924d-07,  0.1520d-06, 
     &   0.2515d-06,  0.4059d-06,  0.6400d-06,  0.9880d-06,  0.1496d-05, 
     &   0.2225d-05,  0.3255d-05,  0.4689d-05,  0.6661d-05,  0.9338d-05, 
     &   0.1293d-04,  0.1769d-04,  0.2395d-04,  0.3210d-04,  0.4259d-04, 
     &   0.5601d-04,  0.7302d-04,  0.9441d-04,  0.1211d-03,  0.1542d-03, 
     &   0.1950d-03,  0.2450d-03,  0.3057d-03,  0.3791d-03,  0.4674d-03, 
     &   0.5730d-03,  0.6986d-03,  0.8474d-03,  0.1023d-02,  0.1228d-02, 
     &   0.1468d-02,  0.1747d-02,  0.2069d-02,  0.2441d-02,  0.2868d-02, 
     &   0.3356d-02,  0.3912d-02,  0.4543d-02,  0.5256d-02,  0.6061d-02, 
     &   0.6964d-02,  0.7975d-02,  0.9104d-02,  0.1036d-01,  0.1175d-01, 
     &   0.1329d-01,  0.1498d-01,  0.1685d-01,  0.1888d-01,  0.2111d-01, 
     &   0.2354d-01,  0.2617d-01,  0.2903d-01,  0.3211d-01,  0.3543d-01, 
     &   0.3900d-01,  0.4282d-01,  0.4690d-01,  0.5125d-01,  0.5587d-01, 
     &   0.6077d-01,  0.6594d-01,  0.7140d-01,  0.7713d-01,  0.8315d-01, 
     &   0.8943d-01,  0.9599d-01,  0.1028d+00,  0.1099d+00,  0.1172d+00, 
     &   0.1247d+00,  0.1324d+00,  0.1404d+00,  0.1484d+00,  0.1567d+00, 
     &   0.1650d+00,  0.1735d+00,  0.1820d+00,  0.1905d+00,  0.1990d+00,
     &   150*0d0/

        data bes11 /
     &   0.0000d+00,  0.1223d-21,  0.2503d-18,  0.2163d-16,  0.5114d-15, 
     &   0.5942d-14,  0.4405d-13,  0.2394d-12,  0.1037d-11,  0.3774d-11, 
     &   0.1198d-10,  0.3403d-10,  0.8820d-10,  0.2116d-09,  0.4755d-09, 
     &   0.1010d-08,  0.2040d-08,  0.3947d-08,  0.7347d-08,  0.1321d-07, 
     &   0.2304d-07,  0.3907d-07,  0.6460d-07,  0.1043d-06,  0.1650d-06, 
     &   0.2559d-06,  0.3897d-06,  0.5837d-06,  0.8607d-06,  0.1251d-05, 
     &   0.1794d-05,  0.2540d-05,  0.3554d-05,  0.4918d-05,  0.6733d-05, 
     &   0.9127d-05,  0.1226d-04,  0.1631d-04,  0.2152d-04,  0.2818d-04, 
     &   0.3660d-04,  0.4720d-04,  0.6044d-04,  0.7688d-04,  0.9716d-04, 
     &   0.1220d-03,  0.1524d-03,  0.1892d-03,  0.2337d-03,  0.2871d-03, 
     &   0.3509d-03,  0.4269d-03,  0.5168d-03,  0.6228d-03,  0.7472d-03, 
     &   0.8928d-03,  0.1062d-02,  0.1259d-02,  0.1486d-02,  0.1748d-02, 
     &   0.2048d-02,  0.2391d-02,  0.2782d-02,  0.3226d-02,  0.3729d-02, 
     &   0.4297d-02,  0.4935d-02,  0.5651d-02,  0.6451d-02,  0.7343d-02, 
     &   0.8335d-02,  0.9434d-02,  0.1065d-01,  0.1199d-01,  0.1346d-01, 
     &   0.1508d-01,  0.1684d-01,  0.1877d-01,  0.2086d-01,  0.2314d-01, 
     &   0.2560d-01,  0.2825d-01,  0.3111d-01,  0.3419d-01,  0.3748d-01, 
     &   0.4100d-01,  0.4476d-01,  0.4875d-01,  0.5299d-01,  0.5748d-01, 
     &   0.6222d-01,  0.6721d-01,  0.7245d-01,  0.7795d-01,  0.8370d-01, 
     &   0.8970d-01,  0.9593d-01,  0.1024d+00,  0.1091d+00,  0.1160d+00, 
     &   0.1231d+00,  0.1304d+00,  0.1379d+00,  0.1455d+00,  0.1532d+00, 
     &   0.1611d+00,  0.1690d+00,  0.1770d+00,  0.1850d+00,  0.1930d+00,
     &   140*0d0/

        data bes12 /
     &   0.0000d+00,  0.5096d-24,  0.2086d-20,  0.2704d-18,  0.8525d-17, 
     &   0.1238d-15,  0.1102d-14,  0.6989d-14,  0.3460d-13,  0.1417d-12, 
     &   0.5000d-12,  0.1563d-11,  0.4420d-11,  0.1149d-10,  0.2783d-10, 
     &   0.6333d-10,  0.1366d-09,  0.2809d-09,  0.5539d-09,  0.1052d-08, 
     &   0.1933d-08,  0.3443d-08,  0.5968d-08,  0.1009d-07,  0.1665d-07, 
     &   0.2693d-07,  0.4268d-07,  0.6645d-07,  0.1017d-06,  0.1533d-06, 
     &   0.2276d-06,  0.3333d-06,  0.4819d-06,  0.6884d-06,  0.9721d-06, 
     &   0.1358d-05,  0.1878d-05,  0.2572d-05,  0.3490d-05,  0.4696d-05, 
     &   0.6264d-05,  0.8292d-05,  0.1089d-04,  0.1421d-04,  0.1840d-04, 
     &   0.2368d-04,  0.3027d-04,  0.3847d-04,  0.4860d-04,  0.6105d-04, 
     &   0.7628d-04,  0.9481d-04,  0.1172d-03,  0.1443d-03,  0.1767d-03, 
     &   0.2155d-03,  0.2616d-03,  0.3163d-03,  0.3807d-03,  0.4565d-03, 
     &   0.5452d-03,  0.6486d-03,  0.7689d-03,  0.9082d-03,  0.1069d-02, 
     &   0.1254d-02,  0.1466d-02,  0.1709d-02,  0.1986d-02,  0.2300d-02, 
     &   0.2656d-02,  0.3058d-02,  0.3510d-02,  0.4019d-02,  0.4589d-02, 
     &   0.5225d-02,  0.5934d-02,  0.6722d-02,  0.7595d-02,  0.8560d-02, 
     &   0.9624d-02,  0.1079d-01,  0.1208d-01,  0.1348d-01,  0.1502d-01, 
     &   0.1669d-01,  0.1851d-01,  0.2048d-01,  0.2261d-01,  0.2491d-01, 
     &   0.2739d-01,  0.3006d-01,  0.3291d-01,  0.3596d-01,  0.3922d-01, 
     &   0.4269d-01,  0.4638d-01,  0.5029d-01,  0.5442d-01,  0.5878d-01, 
     &   0.6337d-01,  0.6819d-01,  0.7325d-01,  0.7853d-01,  0.8405d-01, 
     &   0.8978d-01,  0.9574d-01,  0.1019d+00,  0.1083d+00,  0.1149d+00, 
     &   0.1216d+00,  0.1285d+00,  0.1356d+00,  0.1428d+00,  0.1501d+00, 
     &   0.1575d+00,  0.1650d+00,  0.1726d+00,  0.1802d+00,  0.1877d+00,
     &   130*0d0/

        data bes13 /
     &   0.0000d+00,  0.1960d-26,  0.1605d-22,  0.3120d-20,  0.1312d-18, 
     &   0.2382d-17,  0.2544d-16,  0.1883d-15,  0.1065d-14,  0.4911d-14, 
     &   0.1926d-13,  0.6623d-13,  0.2044d-12,  0.5761d-12,  0.1502d-11, 
     &   0.3665d-11,  0.8433d-11,  0.1844d-10,  0.3852d-10,  0.7728d-10, 
     &   0.1495d-09,  0.2798d-09,  0.5084d-09,  0.8987d-09,  0.1550d-08, 
     &   0.2612d-08,  0.4309d-08,  0.6971d-08,  0.1107d-07,  0.1730d-07, 
     &   0.2659d-07,  0.4028d-07,  0.6017d-07,  0.8872d-07,  0.1292d-06, 
     &   0.1860d-06,  0.2648d-06,  0.3732d-06,  0.5207d-06,  0.7196d-06, 
     &   0.9859d-06,  0.1339d-05,  0.1804d-05,  0.2412d-05,  0.3201d-05, 
     &   0.4218d-05,  0.5520d-05,  0.7177d-05,  0.9274d-05,  0.1191d-04, 
     &   0.1521d-04,  0.1931d-04,  0.2439d-04,  0.3064d-04,  0.3830d-04, 
     &   0.4765d-04,  0.5899d-04,  0.7271d-04,  0.8923d-04,  0.1090d-03, 
     &   0.1327d-03,  0.1608d-03,  0.1941d-03,  0.2335d-03,  0.2797d-03, 
     &   0.3340d-03,  0.3974d-03,  0.4712d-03,  0.5568d-03,  0.6559d-03, 
     &   0.7702d-03,  0.9016d-03,  0.1052d-02,  0.1225d-02,  0.1421d-02, 
     &   0.1644d-02,  0.1897d-02,  0.2183d-02,  0.2505d-02,  0.2868d-02, 
     &   0.3275d-02,  0.3730d-02,  0.4238d-02,  0.4804d-02,  0.5432d-02, 
     &   0.6128d-02,  0.6898d-02,  0.7747d-02,  0.8681d-02,  0.9707d-02, 
     &   0.1083d-01,  0.1206d-01,  0.1340d-01,  0.1486d-01,  0.1644d-01, 
     &   0.1816d-01,  0.2001d-01,  0.2201d-01,  0.2417d-01,  0.2649d-01, 
     &   0.2897d-01,  0.3163d-01,  0.3447d-01,  0.3750d-01,  0.4071d-01, 
     &   0.4413d-01,  0.4775d-01,  0.5157d-01,  0.5560d-01,  0.5984d-01, 
     &   0.6429d-01,  0.6896d-01,  0.7384d-01,  0.7894d-01,  0.8424d-01, 
     &   0.8975d-01,  0.9546d-01,  0.1014d+00,  0.1074d+00,  0.1137d+00, 
     &   0.1201d+00,  0.1267d+00,  0.1335d+00,  0.1403d+00,  0.1473d+00, 
     &   0.1543d+00,  0.1614d+00,  0.1686d+00,  0.1758d+00,  0.1830d+00,
     &   120*0d0/

        data bes14 /
     &   0.0000d+00,  0.7000d-29,  0.1463d-24,  0.3344d-22,  0.1874d-20, 
     &   0.4255d-19,  0.5454d-18,  0.4710d-17,  0.3046d-16,  0.1580d-15, 
     &   0.6885d-15,  0.2606d-14,  0.8776d-14,  0.2680d-13,  0.7529d-13, 
     &   0.1969d-12,  0.4834d-12,  0.1123d-11,  0.2486d-11,  0.5267d-11, 
     &   0.1073d-10,  0.2110d-10,  0.4018d-10,  0.7430d-10,  0.1338d-09, 
     &   0.2349d-09,  0.4034d-09,  0.6781d-09,  0.1118d-08,  0.1810d-08, 
     &   0.2880d-08,  0.4512d-08,  0.6962d-08,  0.1059d-07,  0.1591d-07, 
     &   0.2360d-07,  0.3459d-07,  0.5014d-07,  0.7192d-07,  0.1021d-06, 
     &   0.1436d-06,  0.2002d-06,  0.2766d-06,  0.3789d-06,  0.5151d-06, 
     &   0.6950d-06,  0.9309d-06,  0.1238d-05,  0.1636d-05,  0.2147d-05, 
     &   0.2801d-05,  0.3633d-05,  0.4684d-05,  0.6007d-05,  0.7661d-05, 
     &   0.9721d-05,  0.1227d-04,  0.1542d-04,  0.1928d-04,  0.2401d-04, 
     &   0.2976d-04,  0.3673d-04,  0.4514d-04,  0.5526d-04,  0.6738d-04, 
     &   0.8185d-04,  0.9906d-04,  0.1195d-03,  0.1436d-03,  0.1719d-03, 
     &   0.2052d-03,  0.2441d-03,  0.2895d-03,  0.3423d-03,  0.4035d-03, 
     &   0.4742d-03,  0.5557d-03,  0.6494d-03,  0.7567d-03,  0.8794d-03, 
     &   0.1019d-02,  0.1178d-02,  0.1359d-02,  0.1563d-02,  0.1793d-02, 
     &   0.2052d-02,  0.2344d-02,  0.2670d-02,  0.3035d-02,  0.3442d-02, 
     &   0.3895d-02,  0.4398d-02,  0.4955d-02,  0.5571d-02,  0.6250d-02, 
     &   0.6999d-02,  0.7821d-02,  0.8722d-02,  0.9708d-02,  0.1078d-01, 
     &   0.1196d-01,  0.1323d-01,  0.1462d-01,  0.1612d-01,  0.1774d-01, 
     &   0.1949d-01,  0.2137d-01,  0.2339d-01,  0.2556d-01,  0.2789d-01, 
     &   0.3037d-01,  0.3302d-01,  0.3584d-01,  0.3883d-01,  0.4200d-01, 
     &   0.4536d-01,  0.4891d-01,  0.5265d-01,  0.5658d-01,  0.6071d-01, 
     &   0.6504d-01,  0.6957d-01,  0.7429d-01,  0.7921d-01,  0.8432d-01, 
     &   0.8962d-01,  0.9511d-01,  0.1008d+00,  0.1066d+00,  0.1126d+00, 
     &   0.1188d+00,  0.1251d+00,  0.1315d+00,  0.1380d+00,  0.1447d+00, 
     &   0.1514d+00,  0.1582d+00,  0.1650d+00,  0.1718d+00,  0.1787d+00,
     &   110*0d0/

        data bes15  /
     &   .0000d+00,   .2333d-31,   .7642d-27,   .3344d-24,   .2500d-22,
     &   .7094d-21,   .1091d-19,   .1100d-18,   .8129d-18,   .4744d-17,
     &   .2298d-16,   .9566d-16,   .3516d-15,   .1163d-14,   .3521d-14,
     &   .9866d-14,   .2585d-13,   .6385d-13,   .1497d-12,   .3348d-12,
     &   .7183d-12,   .1484d-11,   .2961d-11,   .5728d-11,   .1077d-10,
     &   .1971d-10,   .3521d-10,   .6150d-10,   .1052d-09,   .1765d-09,
     &   .2908d-09,   .4709d-09,   .7507d-09,   .1179d-08,   .1825d-08,
     &   .2789d-08,   .4208d-08,   .6275d-08,   .9250d-08,   .1349d-07,
     &   .1948d-07,   .2785d-07,   .3945d-07,   .5540d-07,   .7714d-07,
     &   .1065d-06,   .1460d-06,   .1986d-06,   .2683d-06,   .3599d-06,
     &   .4797d-06,   .6352d-06,   .8361d-06,   .1094d-05,   .1423d-05,
     &   .1842d-05,   .2371d-05,   .3036d-05,   .3868d-05,   .4905d-05,
     &   .6192d-05,   .7780d-05,   .9733d-05,   .1212d-04,   .1504d-04,
     &   .1859d-04,   .2288d-04,   .2805d-04,   .3427d-04,   .4171d-04,
     &   .5059d-04,   .6115d-04,   .7368d-04,   .8848d-04,   .1059d-03,
     &   .1264d-03,   .1504d-03,   .1784d-03,   .2109d-03,   .2488d-03,
     &   .2926d-03,   .3432d-03,   .4015d-03,   .4684d-03,   .5451d-03,
     &   .6328d-03,   .7327d-03,   .8464d-03,   .9754d-03,   .1121d-02,
     &   .1286d-02,   .1472d-02,   .1681d-02,   .1916d-02,   .2178d-02,
     &   .2472d-02,   .2799d-02,   .3162d-02,   .3566d-02,   .4013d-02,
     &   .4508d-02,   .5054d-02,   .5655d-02,   .6315d-02,   .7040d-02,
     &   .7834d-02,   .8701d-02,   .9647d-02,   .1068d-01,   .1180d-01,
     &   .1301d-01,   .1432d-01,   .1574d-01,   .1728d-01,   .1892d-01,
     &   .2070d-01,   .2260d-01,   .2464d-01,   .2681d-01,   .2914d-01,
     &   .3161d-01,   .3424d-01,   .3704d-01,   .4000d-01,   .4312d-01,
     &   .4643d-01,   .4991d-01,   .5357d-01,   .5741d-01,   .6143d-01,
     &   .6564d-01,   .7004d-01,   .7462d-01,   .7937d-01,   .8431d-01,
     &   .8943d-01,   .9471d-01,   .1002d+00,   .1058d+00,   .1115d+00,
     &   .1174d+00,   .1235d+00,   .1296d+00,   .1359d+00,   .1422d+00,
     &   .1487d+00,   .1551d+00,   .1617d+00,   .1682d+00,   .1748d+00,
     &   .1813d+00,   .1878d+00,   .1942d+00,   .2005d+00,   .2068d+00,
     &   .2128d+00,   .2187d+00,   .2244d+00,   .2298d+00,   .2350d+00,
     &   .2399d+00,   .2445d+00,   .2488d+00,   .2526d+00,   .2561d+00,
     &   .2591d+00,   .2616d+00,   .2637d+00,   .2652d+00,   .2662d+00,
     &   .2666d+00,   .2664d+00,   .2656d+00,   .2642d+00,   .2621d+00,
     &   .2594d+00,   .2560d+00,   .2519d+00,   .2472d+00,   .2417d+00,
     &   .2356d+00,   .2288d+00,   .2213d+00,   .2131d+00,   .2043d+00,
     &   .1949d+00,   .1848d+00,   .1741d+00,   .1629d+00,   .1512d+00,
     &   .1389d+00,   .1263d+00,   .1131d+00,   .9965d-01,   .8584d-01,
     &   .7176d-01,   .5746d-01,   .4299d-01,   .2842d-01,   .1379d-01,
     &  -.8121d-03,  -.1534d-01,  -.2973d-01,  -.4392d-01,  -.5784d-01,
     &  -.7142d-01,  -.8461d-01,  -.9733d-01,  -.1095d+00,  -.1212d+00,
     &  -.1321d+00,  -.1424d+00,  -.1520d+00,  -.1607d+00,  -.1686d+00,
     &  -.1756d+00,  -.1816d+00,  -.1867d+00,  -.1908d+00,  -.1939d+00,
     &  -.1959d+00,  -.1969d+00,  -.1968d+00,  -.1957d+00,  -.1934d+00,
     &  -.1902d+00,  -.1859d+00,  -.1806d+00,  -.1742d+00,  -.1670d+00,
     &  -.1588d+00,  -.1497d+00,  -.1398d+00,  -.1292d+00,  -.1178d+00,
     &  -.1058d+00,  -.9319d-01,  -.8010d-01,  -.6659d-01,  -.5274d-01,
     &  -.3863d-01,  -.2433d-01,  -.9954d-02,   .4427d-02,   .1872d-01,
     &   .3284d-01,   .4670d-01,   .6021d-01,   .7329d-01,   .8585d-01/

        data bes16  /
     &   .0000d+00,   .7292d-34,   .4777d-29,   .3135d-26,   .3125d-24,
     &   .1109d-22,   .2047d-21,   .2406d-20,   .2034d-19,   .1335d-18,
     &   .7186d-18,   .3292d-17,   .1320d-16,   .4734d-16,   .1543d-15,
     &   .4634d-15,   .1296d-14,   .3401d-14,   .8444d-14,   .1995d-13,
     &   .4506d-13,   .9777d-13,   .2045d-12,   .4137d-12,   .8117d-12,
     &   .1549d-11,   .2879d-11,   .5224d-11,   .9272d-11,   .1612d-10,
     &   .2749d-10,   .4603d-10,   .7579d-10,   .1228d-09,   .1961d-09,
     &   .3086d-09,   .4792d-09,   .7348d-09,   .1113d-08,   .1668d-08,
     &   .2472d-08,   .3625d-08,   .5265d-08,   .7575d-08,   .1080d-07,
     &   .1527d-07,   .2141d-07,   .2979d-07,   .4113d-07,   .5638d-07,
     &   .7675d-07,   .1038d-06,   .1394d-06,   .1861d-06,   .2470d-06,
     &   .3259d-06,   .4275d-06,   .5579d-06,   .7242d-06,   .9353d-06,
     &   .1202d-05,   .1537d-05,   .1957d-05,   .2481d-05,   .3130d-05,
     &   .3934d-05,   .4923d-05,   .6136d-05,   .7618d-05,   .9423d-05,
     &   .1161d-04,   .1426d-04,   .1745d-04,   .2128d-04,   .2586d-04,
     &   .3132d-04,   .3782d-04,   .4553d-04,   .5464d-04,   .6538d-04,
     &   .7801d-04,   .9280d-04,   .1101d-03,   .1303d-03,   .1537d-03,
     &   .1809d-03,   .2123d-03,   .2486d-03,   .2904d-03,   .3384d-03,
     &   .3933d-03,   .4561d-03,   .5277d-03,   .6092d-03,   .7017d-03,
     &   .8064d-03,   .9248d-03,   .1058d-02,   .1209d-02,   .1378d-02,
     &   .1567d-02,   .1778d-02,   .2015d-02,   .2278d-02,   .2571d-02,
     &   .2896d-02,   .3256d-02,   .3654d-02,   .4094d-02,   .4578d-02,
     &   .5110d-02,   .5694d-02,   .6334d-02,   .7034d-02,   .7799d-02,
     &   .8632d-02,   .9538d-02,   .1052d-01,   .1159d-01,   .1274d-01,
     &   .1399d-01,   .1534d-01,   .1679d-01,   .1834d-01,   .2001d-01,
     &   .2181d-01,   .2372d-01,   .2576d-01,   .2794d-01,   .3026d-01,
     &   .3272d-01,   .3534d-01,   .3810d-01,   .4102d-01,   .4411d-01,
     &   .4735d-01,   .5077d-01,   .5435d-01,   .5811d-01,   .6203d-01,
     &   .6613d-01,   .7040d-01,   .7485d-01,   .7946d-01,   .8424d-01,
     &   .8919d-01,   .9429d-01,   .9955d-01,   .1050d+00,   .1105d+00,
     &   .1162d+00,   .1220d+00,   .1279d+00,   .1339d+00,   .1400d+00,
     &   .1461d+00,   .1524d+00,   .1586d+00,   .1649d+00,   .1712d+00,
     &   .1775d+00,   .1837d+00,   .1898d+00,   .1959d+00,   .2019d+00,
     &   .2077d+00,   .2134d+00,   .2189d+00,   .2242d+00,   .2292d+00,
     &   .2340d+00,   .2385d+00,   .2426d+00,   .2465d+00,   .2499d+00,
     &   .2529d+00,   .2555d+00,   .2576d+00,   .2593d+00,   .2605d+00,
     &   .2611d+00,   .2612d+00,   .2607d+00,   .2596d+00,   .2579d+00,
     &   .2556d+00,   .2526d+00,   .2491d+00,   .2449d+00,   .2400d+00,
     &   .2345d+00,   .2283d+00,   .2215d+00,   .2140d+00,   .2059d+00,
     &   .1972d+00,   .1879d+00,   .1780d+00,   .1676d+00,   .1566d+00,
     &   .1452d+00,   .1333d+00,   .1209d+00,   .1082d+00,   .9509d-01,
     &   .8170d-01,   .6807d-01,   .5423d-01,   .4024d-01,   .2615d-01,
     &   .1202d-01,  -.2093d-02,  -.1613d-01,  -.3003d-01,  -.4373d-01,
     &  -.5718d-01,  -.7031d-01,  -.8307d-01,  -.9538d-01,  -.1072d+00,
     &  -.1185d+00,  -.1291d+00,  -.1391d+00,  -.1484d+00,  -.1570d+00,
     &  -.1647d+00,  -.1716d+00,  -.1776d+00,  -.1826d+00,  -.1868d+00,
     &  -.1899d+00,  -.1921d+00,  -.1932d+00,  -.1934d+00,  -.1925d+00,
     &  -.1906d+00,  -.1877d+00,  -.1838d+00,  -.1789d+00,  -.1731d+00,
     &  -.1663d+00,  -.1586d+00,  -.1501d+00,  -.1408d+00,  -.1307d+00,
     &  -.1199d+00,  -.1085d+00,  -.9650d-01,  -.8398d-01,  -.7103d-01/

        data bes17  /
     &   .0000d+00,   .2145d-36,   .2810d-31,   .2767d-28,   .3677d-26,
     &   .1631d-24,   .3613d-23,   .4956d-22,   .4787d-21,   .3537d-20,
     &   .2115d-19,   .1066d-18,   .4665d-18,   .1812d-17,   .6365d-17,
     &   .2048d-16,   .6109d-16,   .1705d-15,   .4482d-15,   .1118d-14,
     &   .2659d-14,   .6060d-14,   .1329d-13,   .2811d-13,   .5757d-13,
     &   .1145d-12,   .2214d-12,   .4174d-12,   .7686d-12,   .1385d-11,
     &   .2444d-11,   .4231d-11,   .7194d-11,   .1203d-10,   .1979d-10,
     &   .3209d-10,   .5129d-10,   .8088d-10,   .1259d-09,   .1938d-09,
     &   .2947d-09,   .4433d-09,   .6600d-09,   .9729d-09,   .1421d-08,
     &   .2056d-08,   .2949d-08,   .4195d-08,   .5920d-08,   .8291d-08,
     &   .1153d-07,   .1591d-07,   .2181d-07,   .2971d-07,   .4021d-07,
     &   .5408d-07,   .7231d-07,   .9614d-07,   .1271d-06,   .1672d-06,
     &   .2187d-06,   .2847d-06,   .3688d-06,   .4756d-06,   .6103d-06,
     &   .7798d-06,   .9921d-06,   .1257d-05,   .1586d-05,   .1993d-05,
     &   .2494d-05,   .3111d-05,   .3865d-05,   .4785d-05,   .5903d-05,
     &   .7258d-05,   .8893d-05,   .1086d-04,   .1322d-04,   .1605d-04,
     &   .1942d-04,   .2343d-04,   .2819d-04,   .3381d-04,   .4044d-04,
     &   .4824d-04,   .5739d-04,   .6809d-04,   .8058d-04,   .9513d-04,
     &   .1120d-03,   .1316d-03,   .1542d-03,   .1803d-03,   .2103d-03,
     &   .2447d-03,   .2842d-03,   .3293d-03,   .3807d-03,   .4392d-03,
     &   .5056d-03,   .5810d-03,   .6661d-03,   .7623d-03,   .8705d-03,
     &   .9923d-03,   .1129d-02,   .1282d-02,   .1453d-02,   .1644d-02,
     &   .1856d-02,   .2093d-02,   .2355d-02,   .2646d-02,   .2967d-02,
     &   .3321d-02,   .3712d-02,   .4142d-02,   .4614d-02,   .5131d-02,
     &   .5698d-02,   .6317d-02,   .6992d-02,   .7727d-02,   .8526d-02,
     &   .9393d-02,   .1033d-01,   .1135d-01,   .1245d-01,   .1363d-01,
     &   .1491d-01,   .1628d-01,   .1775d-01,   .1933d-01,   .2102d-01,
     &   .2282d-01,   .2474d-01,   .2679d-01,   .2897d-01,   .3128d-01,
     &   .3372d-01,   .3631d-01,   .3905d-01,   .4193d-01,   .4497d-01,
     &   .4817d-01,   .5152d-01,   .5503d-01,   .5870d-01,   .6253d-01,
     &   .6653d-01,   .7069d-01,   .7501d-01,   .7948d-01,   .8412d-01,
     &   .8891d-01,   .9384d-01,   .9893d-01,   .1041d+00,   .1095d+00,
     &   .1150d+00,   .1206d+00,   .1262d+00,   .1320d+00,   .1379d+00,
     &   .1438d+00,   .1498d+00,   .1558d+00,   .1618d+00,   .1679d+00,
     &   .1739d+00,   .1799d+00,   .1858d+00,   .1917d+00,   .1974d+00,
     &   .2031d+00,   .2086d+00,   .2139d+00,   .2190d+00,   .2239d+00,
     &   .2286d+00,   .2329d+00,   .2370d+00,   .2408d+00,   .2442d+00,
     &   .2472d+00,   .2499d+00,   .2521d+00,   .2538d+00,   .2551d+00,
     &   .2559d+00,   .2562d+00,   .2560d+00,   .2552d+00,   .2538d+00,
     &   .2519d+00,   .2493d+00,   .2462d+00,   .2425d+00,   .2381d+00,
     &   .2331d+00,   .2275d+00,   .2213d+00,   .2144d+00,   .2070d+00,
     &   .1990d+00,   .1903d+00,   .1812d+00,   .1714d+00,   .1612d+00,
     &   .1505d+00,   .1392d+00,   .1276d+00,   .1156d+00,   .1032d+00,
     &   .9044d-01,   .7744d-01,   .6421d-01,   .5079d-01,   .3723d-01,
     &   .2358d-01,   .9905d-02,  -.3751d-02,  -.1733d-01,  -.3078d-01,
     &  -.4403d-01,  -.5704d-01,  -.6974d-01,  -.8208d-01,  -.9401d-01,
     &  -.1055d+00,  -.1164d+00,  -.1267d+00,  -.1364d+00,  -.1455d+00,
     &  -.1538d+00,  -.1613d+00,  -.1681d+00,  -.1740d+00,  -.1790d+00,
     &  -.1831d+00,  -.1863d+00,  -.1886d+00,  -.1898d+00,  -.1902d+00,
     &  -.1895d+00,  -.1878d+00,  -.1852d+00,  -.1817d+00,  -.1771d+00/

        data bes18  /
     &   .0000d+00,   .5957d-39,   .1561d-33,   .2306d-30,   .4086d-28,
     &   .2265d-26,   .6023d-25,   .9640d-24,   .1064d-22,   .8848d-22,
     &   .5880d-21,   .3260d-20,   .1557d-19,   .6553d-19,   .2479d-18,
     &   .8549d-18,   .2720d-17,   .8066d-17,   .2246d-16,   .5916d-16,
     &   .1482d-15,   .3547d-15,   .8148d-15,   .1803d-14,   .3854d-14,
     &   .7985d-14,   .1607d-13,   .3147d-13,   .6012d-13,   .1122d-12,
     &   .2050d-12,   .3669d-12,   .6443d-12,   .1112d-11,   .1886d-11,
     &   .3148d-11,   .5178d-11,   .8398d-11,   .1344d-10,   .2123d-10,
     &   .3313d-10,   .5112d-10,   .7802d-10,   .1178d-09,   .1762d-09,
     &   .2609d-09,   .3828d-09,   .5568d-09,   .8031d-09,   .1149d-08,
     &   .1631d-08,   .2299d-08,   .3216d-08,   .4467d-08,   .6165d-08,
     &   .8453d-08,   .1152d-07,   .1560d-07,   .2101d-07,   .2813d-07,
     &   .3746d-07,   .4963d-07,   .6541d-07,   .8578d-07,   .1120d-06,
     &   .1454d-06,   .1881d-06,   .2421d-06,   .3104d-06,   .3962d-06,
     &   .5037d-06,   .6378d-06,   .8046d-06,   .1011d-05,   .1266d-05,
     &   .1580d-05,   .1964d-05,   .2433d-05,   .3005d-05,   .3699d-05,
     &   .4538d-05,   .5551d-05,   .6769d-05,   .8230d-05,   .9977d-05,
     &   .1206d-04,   .1454d-04,   .1748d-04,   .2095d-04,   .2505d-04,
     &   .2988d-04,   .3555d-04,   .4218d-04,   .4993d-04,   .5897d-04,
     &   .6947d-04,   .8166d-04,   .9577d-04,   .1121d-03,   .1308d-03,
     &   .1524d-03,   .1772d-03,   .2056d-03,   .2380d-03,   .2750d-03,
     &   .3171d-03,   .3650d-03,   .4192d-03,   .4806d-03,   .5499d-03,
     &   .6280d-03,   .7160d-03,   .8147d-03,   .9255d-03,   .1049d-02,
     &   .1188d-02,   .1342d-02,   .1514d-02,   .1705d-02,   .1917d-02,
     &   .2152d-02,   .2412d-02,   .2699d-02,   .3015d-02,   .3362d-02,
     &   .3745d-02,   .4164d-02,   .4623d-02,   .5124d-02,   .5672d-02,
     &   .6269d-02,   .6919d-02,   .7625d-02,   .8391d-02,   .9221d-02,
     &   .1012d-01,   .1109d-01,   .1213d-01,   .1326d-01,   .1447d-01,
     &   .1577d-01,   .1716d-01,   .1865d-01,   .2024d-01,   .2194d-01,
     &   .2375d-01,   .2568d-01,   .2773d-01,   .2990d-01,   .3220d-01,
     &   .3463d-01,   .3719d-01,   .3990d-01,   .4275d-01,   .4574d-01,
     &   .4888d-01,   .5217d-01,   .5561d-01,   .5920d-01,   .6295d-01,
     &   .6685d-01,   .7090d-01,   .7510d-01,   .7946d-01,   .8396d-01,
     &   .8860d-01,   .9339d-01,   .9831d-01,   .1034d+00,   .1085d+00,
     &   .1138d+00,   .1192d+00,   .1247d+00,   .1303d+00,   .1359d+00,
     &   .1416d+00,   .1474d+00,   .1532d+00,   .1590d+00,   .1648d+00,
     &   .1706d+00,   .1764d+00,   .1821d+00,   .1878d+00,   .1933d+00,
     &   .1988d+00,   .2041d+00,   .2092d+00,   .2142d+00,   .2190d+00,
     &   .2235d+00,   .2278d+00,   .2318d+00,   .2355d+00,   .2389d+00,
     &   .2420d+00,   .2446d+00,   .2469d+00,   .2487d+00,   .2501d+00,
     &   .2511d+00,   .2516d+00,   .2515d+00,   .2510d+00,   .2499d+00,
     &   .2483d+00,   .2461d+00,   .2433d+00,   .2400d+00,   .2361d+00,
     &   .2316d+00,   .2265d+00,   .2208d+00,   .2145d+00,   .2076d+00,
     &   .2002d+00,   .1922d+00,   .1837d+00,   .1746d+00,   .1650d+00,
     &   .1549d+00,   .1444d+00,   .1334d+00,   .1220d+00,   .1102d+00,
     &   .9814d-01,   .8575d-01,   .7309d-01,   .6022d-01,   .4718d-01,
     &   .3402d-01,   .2078d-01,   .7520d-02,  -.5714d-02,  -.1887d-01,
     &  -.3189d-01,  -.4472d-01,  -.5732d-01,  -.6962d-01,  -.8156d-01,
     &  -.9311d-01,  -.1042d+00,  -.1148d+00,  -.1248d+00,  -.1343d+00,
     &  -.1430d+00,  -.1511d+00,  -.1585d+00,  -.1651d+00,  -.1708d+00/

        data bes19  /
     &   .0000d+00,   .1568d-41,   .8217d-36,   .1820d-32,   .4301d-30,
     &   .2981d-28,   .9512d-27,   .1776d-25,   .2242d-24,   .2097d-23,
     &   .1548d-22,   .9446d-22,   .4920d-21,   .2244d-20,   .9144d-20,
     &   .3379d-19,   .1147d-18,   .3615d-18,   .1066d-17,   .2965d-17,
     &   .7819d-17,   .1966d-16,   .4732d-16,   .1095d-15,   .2444d-15,
     &   .5275d-15,   .1104d-14,   .2247d-14,   .4453d-14,   .8612d-14,
     &   .1628d-13,   .3012d-13,   .5463d-13,   .9723d-13,   .1700d-12,
     &   .2923d-12,   .4948d-12,   .8252d-12,   .1357d-11,   .2201d-11,
     &   .3525d-11,   .5578d-11,   .8726d-11,   .1350d-10,   .2066d-10,
     &   .3132d-10,   .4700d-10,   .6990d-10,   .1030d-09,   .1506d-09,
     &   .2183d-09,   .3139d-09,   .4481d-09,   .6350d-09,   .8935d-09,
     &   .1249d-08,   .1734d-08,   .2392d-08,   .3281d-08,   .4472d-08,
     &   .6062d-08,   .8171d-08,   .1096d-07,   .1461d-07,   .1939d-07,
     &   .2561d-07,   .3365d-07,   .4402d-07,   .5733d-07,   .7434d-07,
     &   .9598d-07,   .1234d-06,   .1580d-06,   .2016d-06,   .2561d-06,
     &   .3242d-06,   .4089d-06,   .5139d-06,   .6435d-06,   .8032d-06,
     &   .9992d-06,   .1239d-05,   .1531d-05,   .1887d-05,   .2318d-05,
     &   .2839d-05,   .3467d-05,   .4222d-05,   .5126d-05,   .6208d-05,
     &   .7497d-05,   .9031d-05,   .1085d-04,   .1300d-04,   .1554d-04,
     &   .1854d-04,   .2205d-04,   .2617d-04,   .3099d-04,   .3660d-04,
     &   .4315d-04,   .5075d-04,   .5955d-04,   .6974d-04,   .8150d-04,
     &   .9504d-04,   .1106d-03,   .1285d-03,   .1489d-03,   .1723d-03,
     &   .1990d-03,   .2293d-03,   .2638d-03,   .3029d-03,   .3472d-03,
     &   .3973d-03,   .4538d-03,   .5174d-03,   .5889d-03,   .6691d-03,
     &   .7590d-03,   .8595d-03,   .9718d-03,   .1097d-02,   .1236d-02,
     &   .1391d-02,   .1563d-02,   .1753d-02,   .1963d-02,   .2196d-02,
     &   .2452d-02,   .2734d-02,   .3043d-02,   .3383d-02,   .3756d-02,
     &   .4163d-02,   .4609d-02,   .5094d-02,   .5623d-02,   .6199d-02,
     &   .6824d-02,   .7501d-02,   .8235d-02,   .9029d-02,   .9886d-02,
     &   .1081d-01,   .1181d-01,   .1288d-01,   .1402d-01,   .1526d-01,
     &   .1657d-01,   .1798d-01,   .1949d-01,   .2109d-01,   .2280d-01,
     &   .2462d-01,   .2654d-01,   .2858d-01,   .3075d-01,   .3303d-01,
     &   .3544d-01,   .3798d-01,   .4066d-01,   .4347d-01,   .4642d-01,
     &   .4951d-01,   .5274d-01,   .5611d-01,   .5963d-01,   .6329d-01,
     &   .6710d-01,   .7106d-01,   .7515d-01,   .7939d-01,   .8377d-01,
     &   .8828d-01,   .9292d-01,   .9769d-01,   .1026d+00,   .1076d+00,
     &   .1127d+00,   .1179d+00,   .1232d+00,   .1286d+00,   .1341d+00,
     &   .1396d+00,   .1452d+00,   .1507d+00,   .1564d+00,   .1620d+00,
     &   .1676d+00,   .1732d+00,   .1787d+00,   .1842d+00,   .1895d+00,
     &   .1948d+00,   .2000d+00,   .2050d+00,   .2098d+00,   .2144d+00,
     &   .2189d+00,   .2231d+00,   .2270d+00,   .2307d+00,   .2340d+00,
     &   .2370d+00,   .2397d+00,   .2420d+00,   .2439d+00,   .2455d+00,
     &   .2465d+00,   .2472d+00,   .2473d+00,   .2470d+00,   .2461d+00,
     &   .2448d+00,   .2429d+00,   .2405d+00,   .2375d+00,   .2340d+00,
     &   .2299d+00,   .2253d+00,   .2201d+00,   .2143d+00,   .2080d+00,
     &   .2011d+00,   .1936d+00,   .1857d+00,   .1772d+00,   .1682d+00,
     &   .1587d+00,   .1488d+00,   .1384d+00,   .1276d+00,   .1165d+00,
     &   .1049d+00,   .9312d-01,   .8102d-01,   .6868d-01,   .5614d-01,
     &   .4345d-01,   .3066d-01,   .1779d-01,   .4919d-02,  -.7923d-02,
     &  -.2068d-01,  -.3330d-01,  -.4574d-01,  -.5794d-01,  -.6986d-01/

        data bes20  /
     &   .0000d+00,   .3919d-44,   .4108d-38,   .1365d-34,   .4302d-32,
     &   .3727d-30,   .1427d-28,   .3110d-27,   .4485d-26,   .4720d-25,
     &   .3873d-24,   .2599d-23,   .1477d-22,   .7301d-22,   .3204d-21,
     &   .1269d-20,   .4597d-20,   .1539d-19,   .4808d-19,   .1411d-18,
     &   .3919d-18,   .1035d-17,   .2610d-17,   .6316d-17,   .1471d-16,
     &   .3309d-16,   .7207d-16,   .1523d-15,   .3132d-15,   .6276d-15,
     &   .1228d-14,   .2348d-14,   .4397d-14,   .8074d-14,   .1455d-13,
     &   .2577d-13,   .4488d-13,   .7696d-13,   .1300d-12,   .2166d-12,
     &   .3560d-12,   .5776d-12,   .9260d-12,   .1467d-11,   .2300d-11,
     &   .3567d-11,   .5475d-11,   .8324d-11,   .1254d-10,   .1872d-10,
     &   .2770d-10,   .4067d-10,   .5922d-10,   .8559d-10,   .1228d-09,
     &   .1749d-09,   .2475d-09,   .3478d-09,   .4856d-09,   .6739d-09,
     &   .9296d-09,   .1275d-08,   .1739d-08,   .2358d-08,   .3182d-08,
     &   .4271d-08,   .5704d-08,   .7582d-08,   .1003d-07,   .1321d-07,
     &   .1731d-07,   .2260d-07,   .2938d-07,   .3803d-07,   .4903d-07,
     &   .6296d-07,   .8055d-07,   .1027d-06,   .1304d-06,   .1650d-06,
     &   .2081d-06,   .2615d-06,   .3275d-06,   .4090d-06,   .5090d-06,
     &   .6316d-06,   .7812d-06,   .9635d-06,   .1185d-05,   .1453d-05,
     &   .1777d-05,   .2167d-05,   .2635d-05,   .3196d-05,   .3867d-05,
     &   .4666d-05,   .5617d-05,   .6745d-05,   .8080d-05,   .9656d-05,
     &   .1151d-04,   .1370d-04,   .1626d-04,   .1925d-04,   .2275d-04,
     &   .2683d-04,   .3157d-04,   .3707d-04,   .4345d-04,   .5081d-04,
     &   .5931d-04,   .6909d-04,   .8034d-04,   .9323d-04,   .1080d-03,
     &   .1249d-03,   .1441d-03,   .1660d-03,   .1909d-03,   .2192d-03,
     &   .2512d-03,   .2874d-03,   .3283d-03,   .3743d-03,   .4261d-03,
     &   .4843d-03,   .5496d-03,   .6227d-03,   .7044d-03,   .7955d-03,
     &   .8971d-03,   .1010d-02,   .1136d-02,   .1275d-02,   .1429d-02,
     &   .1600d-02,   .1788d-02,   .1996d-02,   .2225d-02,   .2477d-02,
     &   .2753d-02,   .3055d-02,   .3387d-02,   .3749d-02,   .4145d-02,
     &   .4576d-02,   .5045d-02,   .5556d-02,   .6110d-02,   .6710d-02,
     &   .7360d-02,   .8063d-02,   .8822d-02,   .9640d-02,   .1052d-01,
     &   .1147d-01,   .1249d-01,   .1358d-01,   .1474d-01,   .1599d-01,
     &   .1733d-01,   .1875d-01,   .2027d-01,   .2188d-01,   .2359d-01,
     &   .2541d-01,   .2734d-01,   .2937d-01,   .3152d-01,   .3379d-01,
     &   .3619d-01,   .3870d-01,   .4134d-01,   .4412d-01,   .4702d-01,
     &   .5006d-01,   .5324d-01,   .5655d-01,   .6000d-01,   .6358d-01,
     &   .6731d-01,   .7117d-01,   .7516d-01,   .7929d-01,   .8355d-01,
     &   .8794d-01,   .9245d-01,   .9709d-01,   .1018d+00,   .1067d+00,
     &   .1116d+00,   .1167d+00,   .1218d+00,   .1271d+00,   .1323d+00,
     &   .1377d+00,   .1431d+00,   .1485d+00,   .1539d+00,   .1593d+00,
     &   .1647d+00,   .1701d+00,   .1755d+00,   .1808d+00,   .1860d+00,
     &   .1911d+00,   .1961d+00,   .2010d+00,   .2057d+00,   .2102d+00,
     &   .2145d+00,   .2186d+00,   .2225d+00,   .2261d+00,   .2294d+00,
     &   .2324d+00,   .2351d+00,   .2375d+00,   .2395d+00,   .2410d+00,
     &   .2422d+00,   .2430d+00,   .2433d+00,   .2432d+00,   .2425d+00,
     &   .2414d+00,   .2398d+00,   .2377d+00,   .2351d+00,   .2319d+00,
     &   .2282d+00,   .2240d+00,   .2192d+00,   .2138d+00,   .2080d+00,
     &   .2016d+00,   .1947d+00,   .1872d+00,   .1793d+00,   .1708d+00,
     &   .1619d+00,   .1525d+00,   .1427d+00,   .1325d+00,   .1219d+00,
     &   .1110d+00,   .9967d-01,   .8810d-01,   .7627d-01,   .6422d-01/

        data bes21  /
     &   .0000d+00,   .9332d-47,   .1956d-40,   .9753d-37,   .4097d-34,
     &   .4438d-32,   .2039d-30,   .5184d-29,   .8546d-28,   .1012d-26,
     &   .9228d-26,   .6812d-25,   .4224d-24,   .2262d-23,   .1069d-22,
     &   .4538d-22,   .1753d-21,   .6240d-21,   .2064d-20,   .6398d-20,
     &   .1870d-19,   .5186d-19,   .1371d-18,   .3469d-18,   .8433d-18,
     &   .1976d-17,   .4478d-17,   .9832d-17,   .2097d-16,   .4353d-16,
     &   .8812d-16,   .1742d-15,   .3369d-15,   .6382d-15,   .1185d-14,
     &   .2162d-14,   .3874d-14,   .6831d-14,   .1186d-13,   .2028d-13,
     &   .3420d-13,   .5691d-13,   .9350d-13,   .1518d-12,   .2435d-12,
     &   .3864d-12,   .6067d-12,   .9429d-12,   .1451d-11,   .2213d-11,
     &   .3344d-11,   .5010d-11,   .7443d-11,   .1097d-10,   .1604d-10,
     &   .2329d-10,   .3357d-10,   .4806d-10,   .6833d-10,   .9652d-10,
     &   .1355d-09,   .1891d-09,   .2622d-09,   .3617d-09,   .4961d-09,
     &   .6768d-09,   .9185d-09,   .1240d-08,   .1667d-08,   .2229d-08,
     &   .2966d-08,   .3931d-08,   .5185d-08,   .6811d-08,   .8910d-08,
     &   .1161d-07,   .1506d-07,   .1947d-07,   .2506d-07,   .3215d-07,
     &   .4110d-07,   .5235d-07,   .6645d-07,   .8407d-07,   .1060d-06,
     &   .1332d-06,   .1669d-06,   .2085d-06,   .2596d-06,   .3223d-06,
     &   .3990d-06,   .4925d-06,   .6062d-06,   .7442d-06,   .9111d-06,
     &   .1112d-05,   .1355d-05,   .1646d-05,   .1994d-05,   .2411d-05,
     &   .2907d-05,   .3497d-05,   .4198d-05,   .5027d-05,   .6006d-05,
     &   .7160d-05,   .8518d-05,   .1011d-04,   .1198d-04,   .1416d-04,
     &   .1670d-04,   .1966d-04,   .2310d-04,   .2709d-04,   .3171d-04,
     &   .3704d-04,   .4319d-04,   .5026d-04,   .5839d-04,   .6772d-04,
     &   .7839d-04,   .9058d-04,   .1045d-03,   .1203d-03,   .1384d-03,
     &   .1588d-03,   .1820d-03,   .2082d-03,   .2378d-03,   .2712d-03,
     &   .3087d-03,   .3510d-03,   .3984d-03,   .4516d-03,   .5110d-03,
     &   .5775d-03,   .6516d-03,   .7342d-03,   .8261d-03,   .9281d-03,
     &   .1041d-02,   .1167d-02,   .1305d-02,   .1458d-02,   .1627d-02,
     &   .1813d-02,   .2018d-02,   .2242d-02,   .2489d-02,   .2759d-02,
     &   .3054d-02,   .3376d-02,   .3728d-02,   .4111d-02,   .4528d-02,
     &   .4981d-02,   .5473d-02,   .6006d-02,   .6583d-02,   .7206d-02,
     &   .7879d-02,   .8605d-02,   .9386d-02,   .1023d-01,   .1113d-01,
     &   .1210d-01,   .1313d-01,   .1424d-01,   .1543d-01,   .1669d-01,
     &   .1804d-01,   .1947d-01,   .2100d-01,   .2262d-01,   .2433d-01,
     &   .2615d-01,   .2807d-01,   .3010d-01,   .3224d-01,   .3449d-01,
     &   .3686d-01,   .3935d-01,   .4196d-01,   .4470d-01,   .4756d-01,
     &   .5056d-01,   .5368d-01,   .5693d-01,   .6031d-01,   .6382d-01,
     &   .6746d-01,   .7124d-01,   .7514d-01,   .7917d-01,   .8332d-01,
     &   .8759d-01,   .9199d-01,   .9649d-01,   .1011d+00,   .1058d+00,
     &   .1106d+00,   .1155d+00,   .1205d+00,   .1256d+00,   .1307d+00,
     &   .1359d+00,   .1411d+00,   .1463d+00,   .1516d+00,   .1568d+00,
     &   .1621d+00,   .1673d+00,   .1725d+00,   .1776d+00,   .1827d+00,
     &   .1877d+00,   .1925d+00,   .1972d+00,   .2018d+00,   .2062d+00,
     &   .2105d+00,   .2145d+00,   .2183d+00,   .2219d+00,   .2251d+00,
     &   .2281d+00,   .2308d+00,   .2332d+00,   .2352d+00,   .2369d+00,
     &   .2381d+00,   .2390d+00,   .2395d+00,   .2395d+00,   .2391d+00,
     &   .2382d+00,   .2368d+00,   .2350d+00,   .2326d+00,   .2298d+00,
     &   .2264d+00,   .2225d+00,   .2181d+00,   .2132d+00,   .2078d+00,
     &   .2018d+00,   .1954d+00,   .1884d+00,   .1810d+00,   .1730d+00/

        data bes22  /
     &   .0000d+00,   .2121d-49,   .8893d-43,   .6650d-39,   .3725d-36,
     &   .5044d-34,   .2781d-32,   .8249d-31,   .1554d-29,   .2071d-28,
     &   .2098d-27,   .1704d-26,   .1153d-25,   .6689d-25,   .3405d-24,
     &   .1549d-23,   .6384d-23,   .2414d-22,   .8458d-22,   .2768d-21,
     &   .8518d-21,   .2481d-20,   .6871d-20,   .1818d-19,   .4613d-19,
     &   .1126d-18,   .2655d-18,   .6055d-18,   .1340d-17,   .2881d-17,
     &   .6035d-17,   .1233d-16,   .2463d-16,   .4812d-16,   .9213d-16,
     &   .1730d-15,   .3191d-15,   .5783d-15,   .1031d-14,   .1811d-14,
     &   .3134d-14,   .5348d-14,   .9004d-14,   .1497d-13,   .2459d-13,
     &   .3992d-13,   .6410d-13,   .1018d-12,   .1602d-12,   .2494d-12,
     &   .3848d-12,   .5883d-12,   .8917d-12,   .1340d-11,   .1998d-11,
     &   .2956d-11,   .4341d-11,   .6328d-11,   .9161d-11,   .1317d-10,
     &   .1882d-10,   .2671d-10,   .3768d-10,   .5284d-10,   .7368d-10,
     &   .1022d-09,   .1409d-09,   .1932d-09,   .2637d-09,   .3581d-09,
     &   .4839d-09,   .6509d-09,   .8714d-09,   .1161d-08,   .1541d-08,
     &   .2037d-08,   .2680d-08,   .3512d-08,   .4585d-08,   .5963d-08,
     &   .7725d-08,   .9971d-08,   .1282d-07,   .1644d-07,   .2099d-07,
     &   .2672d-07,   .3391d-07,   .4288d-07,   .5406d-07,   .6795d-07,
     &   .8515d-07,   .1064d-06,   .1325d-06,   .1646d-06,   .2039d-06,
     &   .2519d-06,   .3104d-06,   .3814d-06,   .4674d-06,   .5715d-06,
     &   .6969d-06,   .8477d-06,   .1029d-05,   .1246d-05,   .1504d-05,
     &   .1813d-05,   .2180d-05,   .2615d-05,   .3131d-05,   .3740d-05,
     &   .4458d-05,   .5303d-05,   .6296d-05,   .7459d-05,   .8819d-05,
     &   .1041d-04,   .1226d-04,   .1441d-04,   .1691d-04,   .1980d-04,
     &   .2315d-04,   .2701d-04,   .3147d-04,   .3659d-04,   .4248d-04,
     &   .4923d-04,   .5695d-04,   .6578d-04,   .7585d-04,   .8732d-04,
     &   .1004d-03,   .1152d-03,   .1320d-03,   .1509d-03,   .1724d-03,
     &   .1966d-03,   .2239d-03,   .2546d-03,   .2891d-03,   .3278d-03,
     &   .3711d-03,   .4196d-03,   .4737d-03,   .5341d-03,   .6013d-03,
     &   .6760d-03,   .7590d-03,   .8511d-03,   .9530d-03,   .1066d-02,
     &   .1190d-02,   .1328d-02,   .1479d-02,   .1646d-02,   .1829d-02,
     &   .2029d-02,   .2249d-02,   .2490d-02,   .2753d-02,   .3040d-02,
     &   .3354d-02,   .3695d-02,   .4065d-02,   .4468d-02,   .4905d-02,
     &   .5379d-02,   .5891d-02,   .6445d-02,   .7042d-02,   .7686d-02,
     &   .8380d-02,   .9126d-02,   .9927d-02,   .1079d-01,   .1171d-01,
     &   .1269d-01,   .1375d-01,   .1487d-01,   .1607d-01,   .1735d-01,
     &   .1871d-01,   .2015d-01,   .2168d-01,   .2330d-01,   .2502d-01,
     &   .2683d-01,   .2875d-01,   .3077d-01,   .3290d-01,   .3513d-01,
     &   .3748d-01,   .3995d-01,   .4253d-01,   .4523d-01,   .4805d-01,
     &   .5099d-01,   .5406d-01,   .5726d-01,   .6057d-01,   .6402d-01,
     &   .6758d-01,   .7127d-01,   .7509d-01,   .7902d-01,   .8307d-01,
     &   .8724d-01,   .9152d-01,   .9591d-01,   .1004d+00,   .1050d+00,
     &   .1097d+00,   .1144d+00,   .1193d+00,   .1242d+00,   .1291d+00,
     &   .1342d+00,   .1392d+00,   .1443d+00,   .1494d+00,   .1545d+00,
     &   .1596d+00,   .1647d+00,   .1697d+00,   .1747d+00,   .1796d+00,
     &   .1844d+00,   .1892d+00,   .1938d+00,   .1982d+00,   .2026d+00,
     &   .2067d+00,   .2106d+00,   .2144d+00,   .2179d+00,   .2211d+00,
     &   .2241d+00,   .2268d+00,   .2292d+00,   .2312d+00,   .2329d+00,
     &   .2343d+00,   .2353d+00,   .2358d+00,   .2360d+00,   .2357d+00,
     &   .2350d+00,   .2339d+00,   .2323d+00,   .2302d+00,   .2276d+00/

        data bes23  /
     &   .0000d+00,   .4611d-52,   .3867d-45,   .4337d-41,   .3239d-38,
     &   .5483d-36,   .3628d-34,   .1256d-32,   .2704d-31,   .4053d-30,
     &   .4563d-29,   .4077d-28,   .3009d-27,   .1892d-26,   .1037d-25,
     &   .5055d-25,   .2223d-24,   .8934d-24,   .3314d-23,   .1145d-22,
     &   .3710d-22,   .1135d-21,   .3293d-21,   .9112d-21,   .2413d-20,
     &   .6139d-20,   .1505d-19,   .3566d-19,   .8183d-19,   .1823d-18,
     &   .3952d-18,   .8348d-18,   .1721d-17,   .3469d-17,   .6846d-17,
     &   .1324d-16,   .2512d-16,   .4681d-16,   .8577d-16,   .1546d-15,
     &   .2745d-15,   .4803d-15,   .8288d-15,   .1411d-14,   .2373d-14,
     &   .3942d-14,   .6473d-14,   .1051d-13,   .1689d-13,   .2686d-13,
     &   .4231d-13,   .6601d-13,   .1021d-12,   .1564d-12,   .2378d-12,
     &   .3585d-12,   .5362d-12,   .7961d-12,   .1173d-11,   .1717d-11,
     &   .2496d-11,   .3604d-11,   .5170d-11,   .7372d-11,   .1045d-10,
     &   .1472d-10,   .2063d-10,   .2874d-10,   .3984d-10,   .5492d-10,
     &   .7535d-10,   .1029d-09,   .1397d-09,   .1890d-09,   .2544d-09,
     &   .3409d-09,   .4550d-09,   .6046d-09,   .8001d-09,   .1055d-08,
     &   .1385d-08,   .1811d-08,   .2360d-08,   .3064d-08,   .3964d-08,
     &   .5111d-08,   .6566d-08,   .8408d-08,   .1073d-07,   .1365d-07,
     &   .1732d-07,   .2190d-07,   .2760d-07,   .3469d-07,   .4348d-07,
     &   .5434d-07,   .6772d-07,   .8416d-07,   .1043d-06,   .1290d-06,
     &   .1590d-06,   .1956d-06,   .2400d-06,   .2937d-06,   .3586d-06,
     &   .4367d-06,   .5307d-06,   .6435d-06,   .7784d-06,   .9395d-06,
     &   .1132d-05,   .1360d-05,   .1631d-05,   .1952d-05,   .2331d-05,
     &   .2778d-05,   .3305d-05,   .3924d-05,   .4649d-05,   .5499d-05,
     &   .6491d-05,   .7648d-05,   .8995d-05,   .1056d-04,   .1238d-04,
     &   .1448d-04,   .1691d-04,   .1971d-04,   .2294d-04,   .2666d-04,
     &   .3092d-04,   .3581d-04,   .4141d-04,   .4780d-04,   .5510d-04,
     &   .6341d-04,   .7287d-04,   .8360d-04,   .9578d-04,   .1096d-03,
     &   .1251d-03,   .1427d-03,   .1626d-03,   .1849d-03,   .2100d-03,
     &   .2382d-03,   .2698d-03,   .3051d-03,   .3446d-03,   .3888d-03,
     &   .4379d-03,   .4927d-03,   .5536d-03,   .6212d-03,   .6962d-03,
     &   .7792d-03,   .8710d-03,   .9725d-03,   .1084d-02,   .1208d-02,
     &   .1343d-02,   .1493d-02,   .1656d-02,   .1836d-02,   .2032d-02,
     &   .2247d-02,   .2482d-02,   .2738d-02,   .3017d-02,   .3321d-02,
     &   .3651d-02,   .4010d-02,   .4398d-02,   .4820d-02,   .5275d-02,
     &   .5767d-02,   .6299d-02,   .6871d-02,   .7488d-02,   .8151d-02,
     &   .8864d-02,   .9628d-02,   .1045d-01,   .1132d-01,   .1226d-01,
     &   .1326d-01,   .1433d-01,   .1547d-01,   .1668d-01,   .1797d-01,
     &   .1934d-01,   .2078d-01,   .2232d-01,   .2394d-01,   .2566d-01,
     &   .2747d-01,   .2938d-01,   .3139d-01,   .3350d-01,   .3572d-01,
     &   .3805d-01,   .4049d-01,   .4304d-01,   .4570d-01,   .4849d-01,
     &   .5139d-01,   .5440d-01,   .5754d-01,   .6080d-01,   .6417d-01,
     &   .6767d-01,   .7128d-01,   .7501d-01,   .7886d-01,   .8282d-01,
     &   .8688d-01,   .9106d-01,   .9533d-01,   .9971d-01,   .1042d+00,
     &   .1087d+00,   .1134d+00,   .1181d+00,   .1228d+00,   .1277d+00,
     &   .1325d+00,   .1374d+00,   .1424d+00,   .1473d+00,   .1523d+00,
     &   .1573d+00,   .1622d+00,   .1671d+00,   .1719d+00,   .1767d+00,
     &   .1814d+00,   .1860d+00,   .1905d+00,   .1949d+00,   .1991d+00,
     &   .2031d+00,   .2070d+00,   .2107d+00,   .2141d+00,   .2173d+00,
     &   .2203d+00,   .2230d+00,   .2254d+00,   .2275d+00,   .2292d+00/

        data bes24  /
     &   .0000d+00,   .9606d-55,   .1611d-47,   .2711d-43,   .2700d-40,
     &   .5712d-38,   .4536d-36,   .1831d-34,   .4508d-33,   .7601d-32,
     &   .9511d-31,   .9349d-30,   .7528d-29,   .5127d-28,   .3028d-27,
     &   .1581d-26,   .7419d-26,   .3168d-25,   .1245d-24,   .4539d-24,
     &   .1548d-23,   .4974d-23,   .1512d-22,   .4376d-22,   .1209d-21,
     &   .3206d-21,   .8176d-21,   .2012d-20,   .4789d-20,   .1105d-19,
     &   .2479d-19,   .5413d-19,   .1152d-18,   .2396d-18,   .4873d-18,
     &   .9702d-18,   .1894d-17,   .3629d-17,   .6831d-17,   .1264d-16,
     &   .2303d-16,   .4132d-16,   .7306d-16,   .1274d-15,   .2193d-15,
     &   .3727d-15,   .6259d-15,   .1039d-14,   .1706d-14,   .2770d-14,
     &   .4454d-14,   .7091d-14,   .1118d-13,   .1748d-13,   .2708d-13,
     &   .4160d-13,   .6340d-13,   .9585d-13,   .1438d-12,   .2142d-12,
     &   .3168d-12,   .4653d-12,   .6789d-12,   .9841d-12,   .1418d-11,
     &   .2030d-11,   .2890d-11,   .4090d-11,   .5756d-11,   .8058d-11,
     &   .1122d-10,   .1555d-10,   .2143d-10,   .2941d-10,   .4016d-10,
     &   .5458d-10,   .7386d-10,   .9950d-10,   .1335d-09,   .1783d-09,
     &   .2373d-09,   .3145d-09,   .4151d-09,   .5460d-09,   .7154d-09,
     &   .9340d-09,   .1215d-08,   .1575d-08,   .2035d-08,   .2621d-08,
     &   .3364d-08,   .4305d-08,   .5491d-08,   .6983d-08,   .8853d-08,
     &   .1119d-07,   .1411d-07,   .1773d-07,   .2222d-07,   .2778d-07,
     &   .3463d-07,   .4306d-07,   .5341d-07,   .6607d-07,   .8153d-07,
     &   .1004d-06,   .1232d-06,   .1510d-06,   .1845d-06,   .2251d-06,
     &   .2738d-06,   .3324d-06,   .4027d-06,   .4868d-06,   .5872d-06,
     &   .7068d-06,   .8491d-06,   .1018d-05,   .1218d-05,   .1454d-05,
     &   .1733d-05,   .2062d-05,   .2448d-05,   .2901d-05,   .3432d-05,
     &   .4052d-05,   .4776d-05,   .5620d-05,   .6601d-05,   .7740d-05,
     &   .9060d-05,   .1059d-04,   .1235d-04,   .1439d-04,   .1673d-04,
     &   .1943d-04,   .2253d-04,   .2607d-04,   .3013d-04,   .3477d-04,
     &   .4006d-04,   .4609d-04,   .5295d-04,   .6075d-04,   .6959d-04,
     &   .7960d-04,   .9092d-04,   .1037d-03,   .1181d-03,   .1344d-03,
     &   .1527d-03,   .1732d-03,   .1962d-03,   .2220d-03,   .2509d-03,
     &   .2832d-03,   .3192d-03,   .3593d-03,   .4040d-03,   .4536d-03,
     &   .5087d-03,   .5699d-03,   .6376d-03,   .7124d-03,   .7951d-03,
     &   .8864d-03,   .9870d-03,   .1098d-02,   .1219d-02,   .1353d-02,
     &   .1500d-02,   .1660d-02,   .1836d-02,   .2028d-02,   .2238d-02,
     &   .2466d-02,   .2715d-02,   .2986d-02,   .3280d-02,   .3600d-02,
     &   .3946d-02,   .4321d-02,   .4726d-02,   .5164d-02,   .5637d-02,
     &   .6147d-02,   .6696d-02,   .7286d-02,   .7920d-02,   .8601d-02,
     &   .9331d-02,   .1011d-01,   .1095d-01,   .1184d-01,   .1279d-01,
     &   .1381d-01,   .1489d-01,   .1604d-01,   .1726d-01,   .1856d-01,
     &   .1993d-01,   .2138d-01,   .2292d-01,   .2455d-01,   .2626d-01,
     &   .2807d-01,   .2997d-01,   .3196d-01,   .3406d-01,   .3626d-01,
     &   .3857d-01,   .4098d-01,   .4350d-01,   .4614d-01,   .4888d-01,
     &   .5174d-01,   .5471d-01,   .5779d-01,   .6099d-01,   .6430d-01,
     &   .6773d-01,   .7127d-01,   .7492d-01,   .7868d-01,   .8255d-01,
     &   .8652d-01,   .9060d-01,   .9477d-01,   .9903d-01,   .1034d+00,
     &   .1078d+00,   .1123d+00,   .1169d+00,   .1216d+00,   .1263d+00,
     &   .1310d+00,   .1358d+00,   .1406d+00,   .1454d+00,   .1502d+00,
     &   .1550d+00,   .1598d+00,   .1646d+00,   .1693d+00,   .1740d+00,
     &   .1786d+00,   .1830d+00,   .1874d+00,   .1917d+00,   .1958d+00/

        data bes25  /
     &   .0000d+00,   .1921d-57,   .6444d-50,   .1627d-45,   .2160d-42,
     &   .5712d-40,   .5444d-38,   .2565d-36,   .7214d-35,   .1369d-33,
     &   .1903d-32,   .2058d-31,   .1808d-30,   .1334d-29,   .8484d-29,
     &   .4748d-28,   .2376d-27,   .1078d-26,   .4486d-26,   .1727d-25,
     &   .6204d-25,   .2092d-24,   .6667d-24,   .2017d-23,   .5818d-23,
     &   .1607d-22,   .4263d-22,   .1089d-21,   .2690d-21,   .6433d-21,
     &   .1493d-20,   .3369d-20,   .7405d-20,   .1588d-19,   .3328d-19,
     &   .6824d-19,   .1371d-18,   .2700d-18,   .5221d-18,   .9920d-18,
     &   .1854d-17,   .3410d-17,   .6179d-17,   .1104d-16,   .1944d-16,
     &   .3381d-16,   .5806d-16,   .9850d-16,   .1652d-15,   .2740d-15,
     &   .4498d-15,   .7307d-15,   .1176d-14,   .1873d-14,   .2958d-14,
     &   .4631d-14,   .7188d-14,   .1107d-13,   .1690d-13,   .2562d-13,
     &   .3855d-13,   .5760d-13,   .8546d-13,   .1259d-12,   .1844d-12,
     &   .2683d-12,   .3880d-12,   .5578d-12,   .7973d-12,   .1133d-11,
     &   .1602d-11,   .2252d-11,   .3151d-11,   .4385d-11,   .6073d-11,
     &   .8372d-11,   .1149d-10,   .1569d-10,   .2133d-10,   .2889d-10,
     &   .3895d-10,   .5229d-10,   .6993d-10,   .9316d-10,   .1236d-09,
     &   .1634d-09,   .2153d-09,   .2826d-09,   .3695d-09,   .4817d-09,
     &   .6257d-09,   .8101d-09,   .1045d-08,   .1345d-08,   .1725d-08,
     &   .2205d-08,   .2811d-08,   .3574d-08,   .4529d-08,   .5724d-08,
     &   .7215d-08,   .9068d-08,   .1137d-07,   .1421d-07,   .1773d-07,
     &   .2205d-07,   .2736d-07,   .3387d-07,   .4182d-07,   .5152d-07,
     &   .6333d-07,   .7766d-07,   .9502d-07,   .1160d-06,   .1413d-06,
     &   .1718d-06,   .2083d-06,   .2522d-06,   .3046d-06,   .3672d-06,
     &   .4418d-06,   .5306d-06,   .6358d-06,   .7606d-06,   .9080d-06,
     &   .1082d-05,   .1287d-05,   .1528d-05,   .1812d-05,   .2143d-05,
     &   .2532d-05,   .2985d-05,   .3514d-05,   .4129d-05,   .4844d-05,
     &   .5674d-05,   .6635d-05,   .7746d-05,   .9030d-05,   .1051d-04,
     &   .1221d-04,   .1417d-04,   .1642d-04,   .1899d-04,   .2194d-04,
     &   .2531d-04,   .2915d-04,   .3353d-04,   .3851d-04,   .4417d-04,
     &   .5060d-04,   .5788d-04,   .6611d-04,   .7542d-04,   .8592d-04,
     &   .9775d-04,   .1111d-03,   .1260d-03,   .1429d-03,   .1617d-03,
     &   .1828d-03,   .2064d-03,   .2328d-03,   .2622d-03,   .2949d-03,
     &   .3314d-03,   .3719d-03,   .4169d-03,   .4668d-03,   .5220d-03,
     &   .5831d-03,   .6506d-03,   .7251d-03,   .8073d-03,   .8977d-03,
     &   .9971d-03,   .1106d-02,   .1226d-02,   .1357d-02,   .1501d-02,
     &   .1658d-02,   .1830d-02,   .2017d-02,   .2221d-02,   .2444d-02,
     &   .2685d-02,   .2948d-02,   .3233d-02,   .3541d-02,   .3875d-02,
     &   .4237d-02,   .4627d-02,   .5049d-02,   .5503d-02,   .5992d-02,
     &   .6518d-02,   .7083d-02,   .7689d-02,   .8340d-02,   .9036d-02,
     &   .9781d-02,   .1058d-01,   .1143d-01,   .1233d-01,   .1330d-01,
     &   .1433d-01,   .1542d-01,   .1658d-01,   .1781d-01,   .1911d-01,
     &   .2049d-01,   .2195d-01,   .2349d-01,   .2511d-01,   .2682d-01,
     &   .2862d-01,   .3051d-01,   .3250d-01,   .3458d-01,   .3677d-01,
     &   .3905d-01,   .4144d-01,   .4393d-01,   .4653d-01,   .4923d-01,
     &   .5205d-01,   .5497d-01,   .5801d-01,   .6115d-01,   .6440d-01,
     &   .6777d-01,   .7124d-01,   .7481d-01,   .7849d-01,   .8228d-01,
     &   .8616d-01,   .9014d-01,   .9422d-01,   .9838d-01,   .1026d+00,
     &   .1070d+00,   .1114d+00,   .1158d+00,   .1203d+00,   .1249d+00,
     &   .1295d+00,   .1342d+00,   .1389d+00,   .1436d+00,   .1483d+00/

        do 200 i=1,250
200       bess0(1,i)=bes1(i)
        do 210 i=1,250
210       bess0(2,i)=bes2(i)
        do 220 i=1,250
220       bess0(3,i)=bes3(i)
        do 230 i=1,250
230       bess0(4,i)=bes4(i)
        do 240 i=1,250
240       bess0(5,i)=bes5(i)
        do 250 i=1,250
250       bess0(6,i)=bes6(i)
        do 260 i=1,250
260       bess0(7,i)=bes7(i)
        do 270 i=1,250
270       bess0(8,i)=bes8(i)
        do 280 i=1,250
280       bess0(9,i)=bes9(i)
        do 290 i=1,250
290       bess0(10,i)=bes10(i)
        do 300 i=1,250
300       bess0(11,i)=bes11(i)
        do 310 i=1,250
310       bess0(12,i)=bes12(i)
        do 320 i=1,250
320       bess0(13,i)=bes13(i)
        do 330 i=1,250
330       bess0(14,i)=bes14(i)
        do 340 i=1,250
340       bess0(15,i)=bes15(i)
        do 350 i=1,250
350       bess0(16,i)=bes16(i)
        do 360 i=1,250
360       bess0(17,i)=bes17(i)
        do 370 i=1,250
370       bess0(18,i)=bes18(i)
        do 380 i=1,250
380       bess0(19,i)=bes19(i)
        do 400 i=1,250
400       bess0(20,i)=bes20(i)
        do 410 i=1,250
410       bess0(21,i)=bes21(i)
        do 420 i=1,250
420       bess0(22,i)=bes22(i)
        do 430 i=1,250
430       bess0(23,i)=bes23(i)
        do 440 i=1,250
440       bess0(24,i)=bes24(i)
        do 450 i=1,250
450       bess0(25,i)=bes25(i)

        end

************************************************************************
************************************************************************
	
        subroutine setdat1

        implicit real*8(a-h,o-z),integer(i-n)

        common /bes/ bess0(25,250),bess1(25,250),arg(250)

        dimension bes1(250),bes2(250),bes3(250),bes4(250),bes5(250)
        dimension bes6(250),bes7(250),bes8(250),bes9(250),bes10(250)
        dimension bes11(250),bes12(250),bes13(250),bes14(250),bes15(250)
        dimension bes16(250),bes17(250),bes18(250),bes19(250),bes20(250)
        dimension bes21(250),bes22(250),bes23(250),bes24(250),bes25(250)

      	do 100 i=1,250
100	  arg(i)=dble(i-1)/1d1         

        data bes1 /
     &   0.5000d+00,  0.4981d+00,  0.4925d+00,  0.4832d+00,  0.4703d+00, 
     &   0.4539d+00,  0.4342d+00,  0.4112d+00,  0.3852d+00,  0.3565d+00, 
     &   0.3251d+00,  0.2915d+00,  0.2559d+00,  0.2185d+00,  0.1797d+00, 
     &   0.1399d+00,  0.9922d-01,  0.5812d-01,  0.1692d-01, -0.2405d-01, 
     &  -0.6447d-01, -0.1040d+00, -0.1423d+00, -0.1792d+00, -0.2142d+00, 
     &  -0.2472d+00, -0.2779d+00, -0.3060d+00, -0.3314d+00, -0.3538d+00, 
     &  -0.3731d+00, -0.3891d+00, -0.4019d+00, -0.4112d+00, -0.4170d+00, 
     &  -0.4194d+00, -0.4183d+00, -0.4138d+00, -0.4059d+00, -0.3948d+00, 
     &  -0.3806d+00, -0.3635d+00, -0.3435d+00, -0.3210d+00, -0.2962d+00, 
     &  -0.2692d+00, -0.2404d+00, -0.2100d+00, -0.1782d+00, -0.1455d+00, 
     &  -0.1121d+00, -0.7824d-01, -0.4429d-01, -0.1053d-01,  0.2274d-01, 
     &   0.5524d-01,  0.8667d-01,  0.1168d+00,  0.1453d+00,  0.1721d+00, 
     &   0.1968d+00,  0.2192d+00,  0.2393d+00,  0.2568d+00,  0.2717d+00, 
     &   0.2838d+00,  0.2930d+00,  0.2993d+00,  0.3027d+00,  0.3032d+00, 
     &   0.3007d+00,  0.2955d+00,  0.2875d+00,  0.2769d+00,  0.2638d+00, 
     &   0.2483d+00,  0.2307d+00,  0.2110d+00,  0.1896d+00,  0.1666d+00, 
     &   0.1423d+00,  0.1169d+00,  0.9075d-01,  0.6399d-01,  0.3692d-01, 
     &   0.9807d-02, -0.1709d-01, -0.4352d-01, -0.6924d-01, -0.9401d-01, 
     &  -0.1176d+00, -0.1398d+00, -0.1604d+00, -0.1792d+00, -0.1961d+00, 
     &  -0.2109d+00, -0.2235d+00, -0.2338d+00, -0.2417d+00, -0.2472d+00,
     &   150*0d0/

        data bes2 /
     &   0.0000d+00,  0.2496d-01,  0.4967d-01,  0.7388d-01,  0.9735d-01, 
     &   0.1199d+00,  0.1412d+00,  0.1610d+00,  0.1793d+00,  0.1958d+00, 
     &   0.2102d+00,  0.2226d+00,  0.2327d+00,  0.2404d+00,  0.2457d+00, 
     &   0.2485d+00,  0.2487d+00,  0.2463d+00,  0.2414d+00,  0.2339d+00, 
     &   0.2239d+00,  0.2115d+00,  0.1968d+00,  0.1799d+00,  0.1610d+00, 
     &   0.1402d+00,  0.1178d+00,  0.9378d-01,  0.6851d-01,  0.4217d-01, 
     &   0.1500d-01, -0.1276d-01, -0.4086d-01, -0.6905d-01, -0.9708d-01, 
     &  -0.1247d+00, -0.1516d+00, -0.1777d+00, -0.2026d+00, -0.2261d+00, 
     &  -0.2481d+00, -0.2683d+00, -0.2865d+00, -0.3026d+00, -0.3165d+00, 
     &  -0.3279d+00, -0.3368d+00, -0.3432d+00, -0.3469d+00, -0.3479d+00, 
     &  -0.3462d+00, -0.3419d+00, -0.3349d+00, -0.3253d+00, -0.3132d+00, 
     &  -0.2988d+00, -0.2821d+00, -0.2632d+00, -0.2424d+00, -0.2199d+00, 
     &  -0.1957d+00, -0.1702d+00, -0.1436d+00, -0.1161d+00, -0.8786d-01, 
     &  -0.5925d-01, -0.3046d-01, -0.1753d-02,  0.2663d-01,  0.5444d-01, 
     &   0.8144d-01,  0.1074d+00,  0.1321d+00,  0.1553d+00,  0.1769d+00, 
     &   0.1967d+00,  0.2144d+00,  0.2300d+00,  0.2434d+00,  0.2543d+00, 
     &   0.2629d+00,  0.2689d+00,  0.2725d+00,  0.2734d+00,  0.2719d+00, 
     &   0.2679d+00,  0.2614d+00,  0.2526d+00,  0.2415d+00,  0.2283d+00, 
     &   0.2131d+00,  0.1961d+00,  0.1774d+00,  0.1572d+00,  0.1358d+00, 
     &   0.1133d+00,  0.8993d-01,  0.6595d-01,  0.4157d-01,  0.1703d-01,
     &   150*0d0/

        data bes3 /
     &   0.0000d+00,  0.6243d-03,  0.2490d-02,  0.5572d-02,  0.9834d-02, 
     &   0.1522d-01,  0.2167d-01,  0.2909d-01,  0.3739d-01,  0.4647d-01, 
     &   0.5621d-01,  0.6649d-01,  0.7716d-01,  0.8810d-01,  0.9915d-01, 
     &   0.1102d+00,  0.1210d+00,  0.1315d+00,  0.1415d+00,  0.1508d+00, 
     &   0.1594d+00,  0.1671d+00,  0.1737d+00,  0.1792d+00,  0.1833d+00, 
     &   0.1861d+00,  0.1875d+00,  0.1873d+00,  0.1855d+00,  0.1821d+00, 
     &   0.1770d+00,  0.1703d+00,  0.1619d+00,  0.1519d+00,  0.1403d+00, 
     &   0.1271d+00,  0.1125d+00,  0.9653d-01,  0.7928d-01,  0.6090d-01, 
     &   0.4150d-01,  0.2123d-01,  0.2531d-03, -0.2128d-01, -0.4318d-01, 
     &  -0.6529d-01, -0.8741d-01, -0.1094d+00, -0.1310d+00, -0.1520d+00, 
     &  -0.1723d+00, -0.1918d+00, -0.2101d+00, -0.2272d+00, -0.2429d+00, 
     &  -0.2570d+00, -0.2695d+00, -0.2801d+00, -0.2889d+00, -0.2956d+00, 
     &  -0.3003d+00, -0.3028d+00, -0.3031d+00, -0.3013d+00, -0.2973d+00, 
     &  -0.2911d+00, -0.2828d+00, -0.2724d+00, -0.2600d+00, -0.2457d+00, 
     &  -0.2296d+00, -0.2118d+00, -0.1925d+00, -0.1719d+00, -0.1500d+00, 
     &  -0.1270d+00, -0.1033d+00, -0.7888d-01, -0.5403d-01, -0.2894d-01, 
     &  -0.3817d-02,  0.2113d-01,  0.4568d-01,  0.6965d-01,  0.9282d-01, 
     &   0.1150d+00,  0.1360d+00,  0.1557d+00,  0.1739d+00,  0.1904d+00, 
     &   0.2052d+00,  0.2180d+00,  0.2288d+00,  0.2376d+00,  0.2441d+00, 
     &   0.2485d+00,  0.2507d+00,  0.2506d+00,  0.2483d+00,  0.2438d+00,
     &   150*0d0/

        data bes4 /
     &   0.0000d+00,  0.1041d-04,  0.8308d-04,  0.2794d-03,  0.6587d-03, 
     &   0.1278d-02,  0.2190d-02,  0.3443d-02,  0.5082d-02,  0.7143d-02, 
     &   0.9657d-02,  0.1265d-01,  0.1613d-01,  0.2012d-01,  0.2460d-01, 
     &   0.2958d-01,  0.3504d-01,  0.4094d-01,  0.4725d-01,  0.5394d-01, 
     &   0.6095d-01,  0.6822d-01,  0.7569d-01,  0.8329d-01,  0.9094d-01, 
     &   0.9855d-01,  0.1060d+00,  0.1133d+00,  0.1203d+00,  0.1269d+00, 
     &   0.1330d+00,  0.1385d+00,  0.1434d+00,  0.1475d+00,  0.1508d+00, 
     &   0.1532d+00,  0.1545d+00,  0.1549d+00,  0.1541d+00,  0.1522d+00, 
     &   0.1490d+00,  0.1447d+00,  0.1391d+00,  0.1323d+00,  0.1243d+00, 
     &   0.1150d+00,  0.1045d+00,  0.9294d-01,  0.8024d-01,  0.6652d-01, 
     &   0.5185d-01,  0.3631d-01,  0.2000d-01,  0.3037d-02, -0.1447d-01, 
     &  -0.3240d-01, -0.5063d-01, -0.6900d-01, -0.8740d-01, -0.1057d+00, 
     &  -0.1237d+00, -0.1412d+00, -0.1582d+00, -0.1745d+00, -0.1900d+00, 
     &  -0.2045d+00, -0.2178d+00, -0.2299d+00, -0.2407d+00, -0.2500d+00, 
     &  -0.2577d+00, -0.2638d+00, -0.2683d+00, -0.2709d+00, -0.2718d+00, 
     &  -0.2708d+00, -0.2679d+00, -0.2633d+00, -0.2568d+00, -0.2485d+00, 
     &  -0.2385d+00, -0.2267d+00, -0.2134d+00, -0.1986d+00, -0.1824d+00, 
     &  -0.1649d+00, -0.1462d+00, -0.1265d+00, -0.1060d+00, -0.8474d-01, 
     &  -0.6295d-01, -0.4079d-01, -0.1844d-01,  0.3931d-02,  0.2614d-01, 
     &   0.4800d-01,  0.6935d-01,  0.9001d-01,  0.1098d+00,  0.1286d+00,
     &   150*0d0/

        data bes5 /
     &   0.0000d+00,  0.1301d-06,  0.2078d-05,  0.1049d-04,  0.3302d-04, 
     &   0.8020d-04,  0.1652d-03,  0.3038d-03,  0.5137d-03,  0.8147d-03, 
     &   0.1228d-02,  0.1776d-02,  0.2481d-02,  0.3366d-02,  0.4455d-02, 
     &   0.5770d-02,  0.7332d-02,  0.9159d-02,  0.1127d-01,  0.1368d-01, 
     &   0.1640d-01,  0.1943d-01,  0.2279d-01,  0.2647d-01,  0.3047d-01, 
     &   0.3478d-01,  0.3938d-01,  0.4427d-01,  0.4940d-01,  0.5477d-01, 
     &   0.6032d-01,  0.6603d-01,  0.7185d-01,  0.7773d-01,  0.8363d-01, 
     &   0.8949d-01,  0.9524d-01,  0.1008d+00,  0.1062d+00,  0.1113d+00, 
     &   0.1160d+00,  0.1203d+00,  0.1242d+00,  0.1274d+00,  0.1301d+00, 
     &   0.1321d+00,  0.1333d+00,  0.1338d+00,  0.1335d+00,  0.1322d+00, 
     &   0.1301d+00,  0.1270d+00,  0.1230d+00,  0.1180d+00,  0.1120d+00, 
     &   0.1050d+00,  0.9700d-01,  0.8808d-01,  0.7823d-01,  0.6749d-01, 
     &   0.5590d-01,  0.4352d-01,  0.3041d-01,  0.1664d-01,  0.2288d-02, 
     &  -0.1256d-01, -0.2780d-01, -0.4334d-01, -0.5908d-01, -0.7490d-01, 
     &  -0.9070d-01, -0.1064d+00, -0.1217d+00, -0.1368d+00, -0.1513d+00, 
     &  -0.1652d+00, -0.1783d+00, -0.1906d+00, -0.2020d+00, -0.2123d+00, 
     &  -0.2215d+00, -0.2294d+00, -0.2360d+00, -0.2412d+00, -0.2449d+00, 
     &  -0.2472d+00, -0.2479d+00, -0.2470d+00, -0.2446d+00, -0.2405d+00, 
     &  -0.2349d+00, -0.2277d+00, -0.2190d+00, -0.2088d+00, -0.1972d+00, 
     &  -0.1842d+00, -0.1700d+00, -0.1546d+00, -0.1382d+00, -0.1208d+00,
     &   150*0d0/

        data bes6 /
     &   0.0000d+00,  0.1301d-08,  0.4159d-07,  0.3151d-06,  0.1323d-05, 
     &   0.4021d-05,  0.9953d-05,  0.2138d-04,  0.4138d-04,  0.7397d-04, 
     &   0.1241d-03,  0.1979d-03,  0.3024d-03,  0.4458d-03,  0.6374d-03, 
     &   0.8874d-03,  0.1207d-02,  0.1608d-02,  0.2104d-02,  0.2707d-02, 
     &   0.3432d-02,  0.4293d-02,  0.5302d-02,  0.6475d-02,  0.7824d-02, 
     &   0.9363d-02,  0.1110d-01,  0.1305d-01,  0.1522d-01,  0.1761d-01, 
     &   0.2024d-01,  0.2310d-01,  0.2620d-01,  0.2952d-01,  0.3308d-01, 
     &   0.3685d-01,  0.4083d-01,  0.4500d-01,  0.4934d-01,  0.5383d-01, 
     &   0.5846d-01,  0.6317d-01,  0.6796d-01,  0.7277d-01,  0.7758d-01, 
     &   0.8235d-01,  0.8702d-01,  0.9156d-01,  0.9591d-01,  0.1000d+00, 
     &   0.1039d+00,  0.1074d+00,  0.1105d+00,  0.1132d+00,  0.1155d+00, 
     &   0.1172d+00,  0.1183d+00,  0.1188d+00,  0.1187d+00,  0.1178d+00, 
     &   0.1163d+00,  0.1139d+00,  0.1108d+00,  0.1069d+00,  0.1022d+00, 
     &   0.9672d-01,  0.9039d-01,  0.8325d-01,  0.7532d-01,  0.6661d-01, 
     &   0.5716d-01,  0.4699d-01,  0.3616d-01,  0.2470d-01,  0.1268d-01, 
     &   0.1615d-03, -0.1280d-01, -0.2611d-01, -0.3971d-01, -0.5350d-01, 
     &  -0.6741d-01, -0.8133d-01, -0.9517d-01, -0.1088d+00, -0.1222d+00, 
     &  -0.1352d+00, -0.1478d+00, -0.1597d+00, -0.1710d+00, -0.1816d+00, 
     &  -0.1912d+00, -0.2000d+00, -0.2077d+00, -0.2143d+00, -0.2198d+00, 
     &  -0.2240d+00, -0.2270d+00, -0.2287d+00, -0.2290d+00, -0.2279d+00,
     &   150*0d0/

        data bes7 /
     &   0.0000d+00,  0.1085d-10,  0.6933d-09,  0.7882d-08,  0.4416d-07, 
     &   0.1678d-06,  0.4990d-06,  0.1252d-05,  0.2772d-05,  0.5581d-05, 
     &   0.1042d-04,  0.1831d-04,  0.3057d-04,  0.4892d-04,  0.7548d-04, 
     &   0.1128d-03,  0.1641d-03,  0.2329d-03,  0.3236d-03,  0.4408d-03, 
     &   0.5901d-03,  0.7776d-03,  0.1010d-02,  0.1294d-02,  0.1638d-02, 
     &   0.2050d-02,  0.2539d-02,  0.3114d-02,  0.3785d-02,  0.4560d-02, 
     &   0.5450d-02,  0.6464d-02,  0.7612d-02,  0.8902d-02,  0.1034d-01, 
     &   0.1194d-01,  0.1371d-01,  0.1565d-01,  0.1776d-01,  0.2005d-01, 
     &   0.2253d-01,  0.2519d-01,  0.2803d-01,  0.3104d-01,  0.3423d-01, 
     &   0.3758d-01,  0.4108d-01,  0.4472d-01,  0.4849d-01,  0.5236d-01, 
     &   0.5632d-01,  0.6035d-01,  0.6441d-01,  0.6849d-01,  0.7255d-01, 
     &   0.7656d-01,  0.8049d-01,  0.8430d-01,  0.8796d-01,  0.9142d-01, 
     &   0.9465d-01,  0.9761d-01,  0.1003d+00,  0.1026d+00,  0.1045d+00, 
     &   0.1059d+00,  0.1069d+00,  0.1074d+00,  0.1074d+00,  0.1068d+00, 
     &   0.1056d+00,  0.1038d+00,  0.1013d+00,  0.9818d-01,  0.9437d-01, 
     &   0.8987d-01,  0.8467d-01,  0.7879d-01,  0.7221d-01,  0.6497d-01, 
     &   0.5706d-01,  0.4852d-01,  0.3937d-01,  0.2965d-01,  0.1940d-01, 
     &   0.8663d-02, -0.2512d-02, -0.1407d-01, -0.2594d-01, -0.3807d-01, 
     &  -0.5038d-01, -0.6279d-01, -0.7525d-01, -0.8765d-01, -0.9993d-01, 
     &  -0.1120d+00, -0.1238d+00, -0.1351d+00, -0.1461d+00, -0.1564d+00,
     &   150*0d0/

        data bes8 /
     &   0.0000d+00,  0.7748d-13,  0.9907d-11,  0.1690d-09,  0.1263d-08, 
     &   0.6003d-08,  0.2143d-07,  0.6275d-07,  0.1590d-06,  0.3604d-06, 
     &   0.7485d-06,  0.1448d-05,  0.2641d-05,  0.4585d-05,  0.7630d-05, 
     &   0.1224d-04,  0.1903d-04,  0.2875d-04,  0.4236d-04,  0.6105d-04, 
     &   0.8623d-04,  0.1196d-03,  0.1631d-03,  0.2191d-03,  0.2902d-03, 
     &   0.3795d-03,  0.4904d-03,  0.6266d-03,  0.7924d-03,  0.9923d-03, 
     &   0.1231d-02,  0.1515d-02,  0.1849d-02,  0.2240d-02,  0.2693d-02, 
     &   0.3216d-02,  0.3815d-02,  0.4498d-02,  0.5272d-02,  0.6143d-02, 
     &   0.7119d-02,  0.8206d-02,  0.9412d-02,  0.1074d-01,  0.1220d-01, 
     &   0.1380d-01,  0.1553d-01,  0.1741d-01,  0.1944d-01,  0.2161d-01, 
     &   0.2393d-01,  0.2639d-01,  0.2900d-01,  0.3175d-01,  0.3464d-01, 
     &   0.3765d-01,  0.4077d-01,  0.4401d-01,  0.4734d-01,  0.5074d-01, 
     &   0.5421d-01,  0.5772d-01,  0.6126d-01,  0.6479d-01,  0.6830d-01, 
     &   0.7177d-01,  0.7515d-01,  0.7843d-01,  0.8157d-01,  0.8455d-01, 
     &   0.8733d-01,  0.8988d-01,  0.9217d-01,  0.9416d-01,  0.9582d-01, 
     &   0.9712d-01,  0.9802d-01,  0.9851d-01,  0.9854d-01,  0.9809d-01, 
     &   0.9713d-01,  0.9565d-01,  0.9362d-01,  0.9103d-01,  0.8785d-01, 
     &   0.8408d-01,  0.7972d-01,  0.7475d-01,  0.6919d-01,  0.6303d-01, 
     &   0.5629d-01,  0.4898d-01,  0.4112d-01,  0.3274d-01,  0.2386d-01, 
     &   0.1452d-01,  0.4764d-02, -0.5376d-02, -0.1585d-01, -0.2660d-01,
     &   150*0d0/

        data bes9 /
     &   0.0000d+00,  0.4843d-15,  0.1239d-12,  0.3169d-11,  0.3159d-10, 
     &   0.1878d-09,  0.8047d-09,  0.2751d-08,  0.7969d-08,  0.2034d-07, 
     &   0.4698d-07,  0.1001d-06,  0.1993d-06,  0.3752d-06,  0.6732d-06, 
     &   0.1159d-05,  0.1923d-05,  0.3092d-05,  0.4832d-05,  0.7362d-05, 
     &   0.1096d-04,  0.1599d-04,  0.2290d-04,  0.3222d-04,  0.4463d-04, 
     &   0.6093d-04,  0.8206d-04,  0.1092d-03,  0.1435d-03,  0.1866d-03, 
     &   0.2403d-03,  0.3063d-03,  0.3871d-03,  0.4850d-03,  0.6028d-03, 
     &   0.7435d-03,  0.9105d-03,  0.1107d-02,  0.1338d-02,  0.1606d-02, 
     &   0.1917d-02,  0.2274d-02,  0.2684d-02,  0.3151d-02,  0.3680d-02, 
     &   0.4276d-02,  0.4946d-02,  0.5695d-02,  0.6528d-02,  0.7451d-02, 
     &   0.8469d-02,  0.9587d-02,  0.1081d-01,  0.1214d-01,  0.1359d-01, 
     &   0.1515d-01,  0.1683d-01,  0.1864d-01,  0.2056d-01,  0.2261d-01, 
     &   0.2478d-01,  0.2708d-01,  0.2949d-01,  0.3201d-01,  0.3464d-01, 
     &   0.3738d-01,  0.4020d-01,  0.4311d-01,  0.4609d-01,  0.4913d-01, 
     &   0.5222d-01,  0.5533d-01,  0.5845d-01,  0.6157d-01,  0.6466d-01, 
     &   0.6770d-01,  0.7067d-01,  0.7355d-01,  0.7630d-01,  0.7891d-01, 
     &   0.8134d-01,  0.8358d-01,  0.8558d-01,  0.8734d-01,  0.8880d-01, 
     &   0.8996d-01,  0.9078d-01,  0.9124d-01,  0.9131d-01,  0.9096d-01, 
     &   0.9019d-01,  0.8895d-01,  0.8725d-01,  0.8505d-01,  0.8235d-01, 
     &   0.7914d-01,  0.7540d-01,  0.7114d-01,  0.6634d-01,  0.6102d-01,
     &   150*0d0/

        data bes10 /
     &   0.0000d+00,  0.2690d-17,  0.1376d-14,  0.5284d-13,  0.7024d-12, 
     &   0.5220d-11,  0.2686d-10,  0.1072d-09,  0.3549d-09,  0.1020d-08, 
     &   0.2619d-08,  0.6139d-08,  0.1335d-07,  0.2725d-07,  0.5270d-07, 
     &   0.9728d-07,  0.1724d-06,  0.2948d-06,  0.4885d-06,  0.7866d-06, 
     &   0.1235d-05,  0.1894d-05,  0.2844d-05,  0.4191d-05,  0.6068d-05, 
     &   0.8643d-05,  0.1213d-04,  0.1678d-04,  0.2293d-04,  0.3095d-04, 
     &   0.4130d-04,  0.5454d-04,  0.7130d-04,  0.9234d-04,  0.1185d-03, 
     &   0.1509d-03,  0.1906d-03,  0.2389d-03,  0.2972d-03,  0.3673d-03, 
     &   0.4510d-03,  0.5503d-03,  0.6674d-03,  0.8048d-03,  0.9651d-03, 
     &   0.1151d-02,  0.1366d-02,  0.1614d-02,  0.1897d-02,  0.2219d-02, 
     &   0.2585d-02,  0.2998d-02,  0.3462d-02,  0.3983d-02,  0.4563d-02, 
     &   0.5208d-02,  0.5923d-02,  0.6711d-02,  0.7576d-02,  0.8524d-02, 
     &   0.9559d-02,  0.1068d-01,  0.1190d-01,  0.1322d-01,  0.1463d-01, 
     &   0.1615d-01,  0.1777d-01,  0.1949d-01,  0.2132d-01,  0.2326d-01, 
     &   0.2529d-01,  0.2743d-01,  0.2967d-01,  0.3200d-01,  0.3442d-01, 
     &   0.3692d-01,  0.3950d-01,  0.4214d-01,  0.4484d-01,  0.4759d-01, 
     &   0.5036d-01,  0.5316d-01,  0.5596d-01,  0.5874d-01,  0.6150d-01, 
     &   0.6421d-01,  0.6685d-01,  0.6941d-01,  0.7186d-01,  0.7417d-01, 
     &   0.7633d-01,  0.7832d-01,  0.8010d-01,  0.8166d-01,  0.8297d-01, 
     &   0.8402d-01,  0.8476d-01,  0.8519d-01,  0.8529d-01,  0.8502d-01,
     &   150*0d0/

        data bes11 /
     &   0.0000d+00,  0.1345d-19,  0.1377d-16,  0.7928d-15,  0.1405d-13, 
     &   0.1306d-12,  0.8064d-12,  0.3755d-11,  0.1422d-10,  0.4599d-10, 
     &   0.1313d-09,  0.3387d-09,  0.8041d-09,  0.1779d-08,  0.3708d-08, 
     &   0.7340d-08,  0.1389d-07,  0.2526d-07,  0.4435d-07,  0.7545d-07, 
     &   0.1248d-06,  0.2012d-06,  0.3170d-06,  0.4889d-06,  0.7396d-06, 
     &   0.1099d-05,  0.1606d-05,  0.2311d-05,  0.3280d-05,  0.4592d-05, 
     &   0.6350d-05,  0.8680d-05,  0.1174d-04,  0.1570d-04,  0.2081d-04, 
     &   0.2733d-04,  0.3557d-04,  0.4592d-04,  0.5882d-04,  0.7478d-04, 
     &   0.9439d-04,  0.1183d-03,  0.1474d-03,  0.1825d-03,  0.2245d-03, 
     &   0.2747d-03,  0.3342d-03,  0.4044d-03,  0.4870d-03,  0.5835d-03, 
     &   0.6958d-03,  0.8258d-03,  0.9759d-03,  0.1148d-02,  0.1345d-02, 
     &   0.1570d-02,  0.1825d-02,  0.2113d-02,  0.2438d-02,  0.2802d-02, 
     &   0.3209d-02,  0.3663d-02,  0.4167d-02,  0.4725d-02,  0.5341d-02, 
     &   0.6017d-02,  0.6758d-02,  0.7568d-02,  0.8450d-02,  0.9407d-02, 
     &   0.1044d-01,  0.1156d-01,  0.1276d-01,  0.1405d-01,  0.1542d-01, 
     &   0.1689d-01,  0.1844d-01,  0.2009d-01,  0.2183d-01,  0.2366d-01, 
     &   0.2557d-01,  0.2757d-01,  0.2966d-01,  0.3182d-01,  0.3406d-01, 
     &   0.3637d-01,  0.3874d-01,  0.4116d-01,  0.4363d-01,  0.4613d-01, 
     &   0.4865d-01,  0.5119d-01,  0.5372d-01,  0.5624d-01,  0.5873d-01, 
     &   0.6117d-01,  0.6355d-01,  0.6584d-01,  0.6804d-01,  0.7012d-01, 
     &   0.7206d-01,  0.7384d-01,  0.7544d-01,  0.7685d-01,  0.7803d-01, 
     &   0.7898d-01,  0.7967d-01,  0.8007d-01,  0.8018d-01,  0.7996d-01,
     &   140*0d0/

        data bes12 /
     &   0.0000d+00,  0.6115d-22,  0.1251d-18,  0.1081d-16,  0.2556d-15, 
     &   0.2970d-14,  0.2201d-13,  0.1196d-12,  0.5179d-12,  0.1885d-11, 
     &   0.5980d-11,  0.1698d-10,  0.4400d-10,  0.1055d-09,  0.2370d-09, 
     &   0.5029d-09,  0.1016d-08,  0.1964d-08,  0.3654d-08,  0.6569d-08, 
     &   0.1145d-07,  0.1940d-07,  0.3204d-07,  0.5172d-07,  0.8172d-07, 
     &   0.1266d-06,  0.1927d-06,  0.2883d-06,  0.4248d-06,  0.6169d-06, 
     &   0.8837d-06,  0.1250d-05,  0.1747d-05,  0.2414d-05,  0.3302d-05, 
     &   0.4470d-05,  0.5995d-05,  0.7969d-05,  0.1050d-04,  0.1373d-04, 
     &   0.1781d-04,  0.2293d-04,  0.2932d-04,  0.3723d-04,  0.4698d-04, 
     &   0.5892d-04,  0.7344d-04,  0.9103d-04,  0.1122d-03,  0.1376d-03, 
     &   0.1679d-03,  0.2038d-03,  0.2462d-03,  0.2961d-03,  0.3545d-03, 
     &   0.4226d-03,  0.5016d-03,  0.5931d-03,  0.6984d-03,  0.8194d-03, 
     &   0.9576d-03,  0.1115d-02,  0.1294d-02,  0.1496d-02,  0.1725d-02, 
     &   0.1981d-02,  0.2269d-02,  0.2590d-02,  0.2947d-02,  0.3344d-02, 
     &   0.3782d-02,  0.4266d-02,  0.4798d-02,  0.5382d-02,  0.6020d-02, 
     &   0.6716d-02,  0.7472d-02,  0.8292d-02,  0.9179d-02,  0.1013d-01, 
     &   0.1116d-01,  0.1226d-01,  0.1344d-01,  0.1469d-01,  0.1602d-01, 
     &   0.1744d-01,  0.1893d-01,  0.2050d-01,  0.2216d-01,  0.2389d-01, 
     &   0.2569d-01,  0.2758d-01,  0.2953d-01,  0.3155d-01,  0.3363d-01, 
     &   0.3577d-01,  0.3796d-01,  0.4019d-01,  0.4246d-01,  0.4476d-01, 
     &   0.4707d-01,  0.4939d-01,  0.5170d-01,  0.5400d-01,  0.5626d-01, 
     &   0.5848d-01,  0.6064d-01,  0.6273d-01,  0.6472d-01,  0.6660d-01, 
     &   0.6836d-01,  0.6998d-01,  0.7143d-01,  0.7271d-01,  0.7379d-01, 
     &   0.7465d-01,  0.7528d-01,  0.7566d-01,  0.7577d-01,  0.7560d-01,
     &   130*0d0/

        data bes13 /
     &   0.0000d+00,  0.2548d-24,  0.1043d-20,  0.1352d-18,  0.4261d-17, 
     &   0.6190d-16,  0.5506d-15,  0.3492d-14,  0.1728d-13,  0.7078d-13, 
     &   0.2496d-12,  0.7801d-12,  0.2206d-11,  0.5734d-11,  0.1388d-10, 
     &   0.3156d-10,  0.6804d-10,  0.1399d-09,  0.2757d-09,  0.5235d-09, 
     &   0.9610d-09,  0.1711d-08,  0.2964d-08,  0.5006d-08,  0.8261d-08, 
     &   0.1335d-07,  0.2114d-07,  0.3288d-07,  0.5030d-07,  0.7573d-07, 
     &   0.1123d-06,  0.1644d-06,  0.2375d-06,  0.3389d-06,  0.4781d-06, 
     &   0.6672d-06,  0.9217d-06,  0.1261d-05,  0.1709d-05,  0.2297d-05, 
     &   0.3060d-05,  0.4046d-05,  0.5308d-05,  0.6914d-05,  0.8942d-05, 
     &   0.1149d-04,  0.1467d-04,  0.1861d-04,  0.2348d-04,  0.2945d-04, 
     &   0.3674d-04,  0.4559d-04,  0.5628d-04,  0.6914d-04,  0.8454d-04, 
     &   0.1029d-03,  0.1247d-03,  0.1504d-03,  0.1807d-03,  0.2162d-03, 
     &   0.2577d-03,  0.3060d-03,  0.3619d-03,  0.4265d-03,  0.5008d-03, 
     &   0.5861d-03,  0.6836d-03,  0.7948d-03,  0.9210d-03,  0.1064d-02, 
     &   0.1225d-02,  0.1407d-02,  0.1610d-02,  0.1838d-02,  0.2093d-02, 
     &   0.2375d-02,  0.2689d-02,  0.3036d-02,  0.3419d-02,  0.3840d-02, 
     &   0.4302d-02,  0.4808d-02,  0.5360d-02,  0.5961d-02,  0.6613d-02, 
     &   0.7320d-02,  0.8083d-02,  0.8905d-02,  0.9789d-02,  0.1074d-01, 
     &   0.1175d-01,  0.1283d-01,  0.1398d-01,  0.1520d-01,  0.1649d-01, 
     &   0.1785d-01,  0.1928d-01,  0.2078d-01,  0.2236d-01,  0.2400d-01, 
     &   0.2571d-01,  0.2748d-01,  0.2932d-01,  0.3121d-01,  0.3315d-01, 
     &   0.3515d-01,  0.3719d-01,  0.3926d-01,  0.4136d-01,  0.4348d-01, 
     &   0.4562d-01,  0.4775d-01,  0.4988d-01,  0.5198d-01,  0.5406d-01, 
     &   0.5609d-01,  0.5807d-01,  0.5997d-01,  0.6179d-01,  0.6351d-01, 
     &   0.6512d-01,  0.6660d-01,  0.6792d-01,  0.6909d-01,  0.7008d-01, 
     &   0.7088d-01,  0.7146d-01,  0.7182d-01,  0.7194d-01,  0.7180d-01,
     &   120*0d0/

        data bes14 /
     &   0.0000d+00,  0.9800d-27,  0.8025d-23,  0.1560d-20,  0.6558d-19, 
     &   0.1191d-17,  0.1271d-16,  0.9408d-16,  0.5323d-15,  0.2453d-14, 
     &   0.9617d-14,  0.3307d-13,  0.1020d-12,  0.2875d-12,  0.7494d-12, 
     &   0.1827d-11,  0.4204d-11,  0.9187d-11,  0.1919d-10,  0.3847d-10, 
     &   0.7439d-10,  0.1392d-09,  0.2527d-09,  0.4465d-09,  0.7695d-09, 
     &   0.1296d-08,  0.2137d-08,  0.3455d-08,  0.5484d-08,  0.8560d-08, 
     &   0.1315d-07,  0.1990d-07,  0.2971d-07,  0.4377d-07,  0.6369d-07, 
     &   0.9160d-07,  0.1303d-06,  0.1834d-06,  0.2557d-06,  0.3531d-06, 
     &   0.4832d-06,  0.6556d-06,  0.8824d-06,  0.1178d-05,  0.1562d-05, 
     &   0.2056d-05,  0.2687d-05,  0.3489d-05,  0.4503d-05,  0.5775d-05, 
     &   0.7364d-05,  0.9337d-05,  0.1178d-04,  0.1477d-04,  0.1844d-04, 
     &   0.2290d-04,  0.2831d-04,  0.3484d-04,  0.4268d-04,  0.5206d-04, 
     &   0.6324d-04,  0.7651d-04,  0.9219d-04,  0.1107d-03,  0.1323d-03, 
     &   0.1577d-03,  0.1872d-03,  0.2216d-03,  0.2613d-03,  0.3071d-03, 
     &   0.3598d-03,  0.4202d-03,  0.4893d-03,  0.5680d-03,  0.6575d-03, 
     &   0.7588d-03,  0.8733d-03,  0.1002d-02,  0.1147d-02,  0.1310d-02, 
     &   0.1491d-02,  0.1693d-02,  0.1918d-02,  0.2168d-02,  0.2443d-02, 
     &   0.2748d-02,  0.3082d-02,  0.3450d-02,  0.3853d-02,  0.4293d-02, 
     &   0.4772d-02,  0.5293d-02,  0.5859d-02,  0.6470d-02,  0.7131d-02, 
     &   0.7842d-02,  0.8607d-02,  0.9426d-02,  0.1030d-01,  0.1124d-01, 
     &   0.1223d-01,  0.1329d-01,  0.1441d-01,  0.1559d-01,  0.1684d-01, 
     &   0.1815d-01,  0.1952d-01,  0.2096d-01,  0.2246d-01,  0.2402d-01, 
     &   0.2564d-01,  0.2732d-01,  0.2905d-01,  0.3083d-01,  0.3266d-01, 
     &   0.3453d-01,  0.3643d-01,  0.3836d-01,  0.4032d-01,  0.4229d-01, 
     &   0.4427d-01,  0.4625d-01,  0.4821d-01,  0.5016d-01,  0.5207d-01, 
     &   0.5395d-01,  0.5577d-01,  0.5752d-01,  0.5919d-01,  0.6078d-01, 
     &   0.6225d-01,  0.6361d-01,  0.6483d-01,  0.6591d-01,  0.6682d-01, 
     &   0.6755d-01,  0.6810d-01,  0.6843d-01,  0.6855d-01,  0.6844d-01,
     &   110*0d0/

        data bes15  /
     &   .0000d+00,   .3500d-29,   .5731d-25,   .1672d-22,   .9370d-21,
     &   .2127d-19,   .2726d-18,   .2354d-17,   .1522d-16,   .7894d-16,
     &   .3439d-15,   .1301d-14,   .4381d-14,   .1338d-13,   .3757d-13,
     &   .9820d-13,   .2410d-12,   .5600d-12,   .1239d-11,   .2623d-11,
     &   .5342d-11,   .1050d-10,   .1999d-10,   .3694d-10,   .6647d-10,
     &   .1167d-09,   .2002d-09,   .3364d-09,   .5543d-09,   .8968d-09,
     &   .1426d-08,   .2233d-08,   .3443d-08,   .5236d-08,   .7857d-08,
     &   .1164d-07,   .1706d-07,   .2470d-07,   .3540d-07,   .5022d-07,
     &   .7057d-07,   .9827d-07,   .1356d-06,   .1857d-06,   .2522d-06,
     &   .3399d-06,   .4547d-06,   .6042d-06,   .7973d-06,   .1045d-05,
     &   .1362d-05,   .1765d-05,   .2272d-05,   .2910d-05,   .3707d-05,
     &   .4697d-05,   .5923d-05,   .7431d-05,   .9280d-05,   .1154d-04,
     &   .1428d-04,   .1759d-04,   .2159d-04,   .2639d-04,   .3212d-04,
     &   .3896d-04,   .4707d-04,   .5666d-04,   .6797d-04,   .8125d-04,
     &   .9680d-04,   .1149d-03,   .1360d-03,   .1605d-03,   .1888d-03,
     &   .2214d-03,   .2589d-03,   .3019d-03,   .3510d-03,   .4070d-03,
     &   .4706d-03,   .5427d-03,   .6243d-03,   .7163d-03,   .8197d-03,
     &   .9357d-03,   .1066d-02,   .1211d-02,   .1372d-02,   .1552d-02,
     &   .1751d-02,   .1971d-02,   .2214d-02,   .2481d-02,   .2774d-02,
     &   .3096d-02,   .3448d-02,   .3832d-02,   .4250d-02,   .4703d-02,
     &   .5195d-02,   .5727d-02,   .6301d-02,   .6919d-02,   .7583d-02,
     &   .8295d-02,   .9056d-02,   .9869d-02,   .1073d-01,   .1165d-01,
     &   .1263d-01,   .1366d-01,   .1475d-01,   .1590d-01,   .1710d-01,
     &   .1837d-01,   .1969d-01,   .2106d-01,   .2250d-01,   .2398d-01,
     &   .2552d-01,   .2711d-01,   .2875d-01,   .3043d-01,   .3215d-01,
     &   .3391d-01,   .3569d-01,   .3750d-01,   .3933d-01,   .4117d-01,
     &   .4302d-01,   .4486d-01,   .4669d-01,   .4850d-01,   .5027d-01,
     &   .5201d-01,   .5369d-01,   .5532d-01,   .5687d-01,   .5833d-01,
     &   .5969d-01,   .6095d-01,   .6208d-01,   .6307d-01,   .6392d-01,
     &   .6460d-01,   .6511d-01,   .6543d-01,   .6555d-01,   .6545d-01,
     &   .6513d-01,   .6458d-01,   .6378d-01,   .6272d-01,   .6140d-01,
     &   .5980d-01,   .5792d-01,   .5576d-01,   .5330d-01,   .5054d-01,
     &   .4749d-01,   .4414d-01,   .4049d-01,   .3654d-01,   .3231d-01,
     &   .2778d-01,   .2298d-01,   .1791d-01,   .1258d-01,   .7007d-02,
     &   .1206d-02,  -.4808d-02,  -.1101d-01,  -.1739d-01,  -.2392d-01,
     &  -.3057d-01,  -.3733d-01,  -.4415d-01,  -.5102d-01,  -.5790d-01,
     &  -.6475d-01,  -.7156d-01,  -.7828d-01,  -.8489d-01,  -.9133d-01,
     &  -.9759d-01,  -.1036d+00,  -.1094d+00,  -.1149d+00,  -.1200d+00,
     &  -.1248d+00,  -.1291d+00,  -.1331d+00,  -.1366d+00,  -.1395d+00,
     &  -.1420d+00,  -.1439d+00,  -.1453d+00,  -.1461d+00,  -.1462d+00,
     &  -.1458d+00,  -.1447d+00,  -.1430d+00,  -.1406d+00,  -.1376d+00,
     &  -.1340d+00,  -.1297d+00,  -.1248d+00,  -.1192d+00,  -.1131d+00,
     &  -.1064d+00,  -.9915d-01,  -.9139d-01,  -.8313d-01,  -.7442d-01,
     &  -.6530d-01,  -.5581d-01,  -.4599d-01,  -.3590d-01,  -.2559d-01,
     &  -.1510d-01,  -.4499d-02,   .6159d-02,   .1682d-01,   .2742d-01,
     &   .3789d-01,   .4819d-01,   .5824d-01,   .6800d-01,   .7739d-01,
     &   .8637d-01,   .9487d-01,   .1028d+00,   .1102d+00,   .1170d+00,
     &   .1231d+00,   .1285d+00,   .1331d+00,   .1370d+00,   .1400d+00,
     &   .1422d+00,   .1435d+00,   .1440d+00,   .1435d+00,   .1422d+00,
     &   .1400d+00,   .1370d+00,   .1331d+00,   .1283d+00,   .1227d+00/

        data bes16  /
     &   .0000d+00,   .1167d-31,   .3821d-27,   .1672d-24,   .1250d-22,
     &   .3546d-21,   .5454d-20,   .5495d-19,   .4062d-18,   .2370d-17,
     &   .1148d-16,   .4778d-16,   .1755d-15,   .5808d-15,   .1757d-14,
     &   .4923d-14,   .1289d-13,   .3184d-13,   .7461d-13,   .1669d-12,
     &   .3578d-12,   .7388d-12,   .1474d-11,   .2850d-11,   .5354d-11,
     &   .9796d-11,   .1749d-10,   .3054d-10,   .5222d-10,   .8755d-10,
     &   .1442d-09,   .2334d-09,   .3718d-09,   .5835d-09,   .9028d-09,
     &   .1379d-08,   .2079d-08,   .3097d-08,   .4562d-08,   .6649d-08,
     &   .9592d-08,   .1370d-07,   .1940d-07,   .2721d-07,   .3786d-07,
     &   .5224d-07,   .7153d-07,   .9722d-07,   .1312d-06,   .1758d-06,
     &   .2341d-06,   .3097d-06,   .4071d-06,   .5322d-06,   .6916d-06,
     &   .8939d-06,   .1149d-05,   .1470d-05,   .1871d-05,   .2369d-05,
     &   .2986d-05,   .3748d-05,   .4682d-05,   .5825d-05,   .7216d-05,
     &   .8903d-05,   .1094d-04,   .1340d-04,   .1634d-04,   .1986d-04,
     &   .2405d-04,   .2902d-04,   .3491d-04,   .4185d-04,   .5000d-04,
     &   .5956d-04,   .7073d-04,   .8375d-04,   .9886d-04,   .1164d-03,
     &   .1366d-03,   .1599d-03,   .1866d-03,   .2173d-03,   .2523d-03,
     &   .2923d-03,   .3377d-03,   .3892d-03,   .4474d-03,   .5132d-03,
     &   .5872d-03,   .6703d-03,   .7635d-03,   .8678d-03,   .9840d-03,
     &   .1113d-02,   .1257d-02,   .1417d-02,   .1593d-02,   .1787d-02,
     &   .2001d-02,   .2236d-02,   .2494d-02,   .2777d-02,   .3085d-02,
     &   .3421d-02,   .3786d-02,   .4182d-02,   .4612d-02,   .5076d-02,
     &   .5576d-02,   .6115d-02,   .6694d-02,   .7315d-02,   .7979d-02,
     &   .8688d-02,   .9444d-02,   .1025d-01,   .1110d-01,   .1200d-01,
     &   .1296d-01,   .1396d-01,   .1502d-01,   .1613d-01,   .1730d-01,
     &   .1852d-01,   .1979d-01,   .2111d-01,   .2248d-01,   .2390d-01,
     &   .2537d-01,   .2688d-01,   .2843d-01,   .3002d-01,   .3165d-01,
     &   .3330d-01,   .3499d-01,   .3669d-01,   .3840d-01,   .4013d-01,
     &   .4186d-01,   .4358d-01,   .4529d-01,   .4697d-01,   .4863d-01,
     &   .5025d-01,   .5181d-01,   .5332d-01,   .5476d-01,   .5612d-01,
     &   .5739d-01,   .5855d-01,   .5960d-01,   .6053d-01,   .6132d-01,
     &   .6195d-01,   .6243d-01,   .6273d-01,   .6285d-01,   .6277d-01,
     &   .6249d-01,   .6199d-01,   .6126d-01,   .6030d-01,   .5909d-01,
     &   .5763d-01,   .5591d-01,   .5392d-01,   .5167d-01,   .4914d-01,
     &   .4633d-01,   .4325d-01,   .3989d-01,   .3625d-01,   .3234d-01,
     &   .2816d-01,   .2372d-01,   .1903d-01,   .1409d-01,   .8913d-02,
     &   .3520d-02,  -.2080d-02,  -.7869d-02,  -.1383d-01,  -.1994d-01,
     &  -.2618d-01,  -.3253d-01,  -.3896d-01,  -.4545d-01,  -.5197d-01,
     &  -.5849d-01,  -.6498d-01,  -.7142d-01,  -.7776d-01,  -.8399d-01,
     &  -.9006d-01,  -.9594d-01,  -.1016d+00,  -.1070d+00,  -.1121d+00,
     &  -.1170d+00,  -.1214d+00,  -.1255d+00,  -.1292d+00,  -.1324d+00,
     &  -.1352d+00,  -.1375d+00,  -.1392d+00,  -.1405d+00,  -.1412d+00,
     &  -.1413d+00,  -.1408d+00,  -.1398d+00,  -.1381d+00,  -.1359d+00,
     &  -.1330d+00,  -.1295d+00,  -.1255d+00,  -.1208d+00,  -.1156d+00,
     &  -.1097d+00,  -.1034d+00,  -.9653d-01,  -.8916d-01,  -.8134d-01,
     &  -.7308d-01,  -.6442d-01,  -.5540d-01,  -.4608d-01,  -.3648d-01,
     &  -.2666d-01,  -.1666d-01,  -.6543d-02,   .3645d-02,   .1385d-01,
     &   .2401d-01,   .3408d-01,   .4399d-01,   .5369d-01,   .6313d-01,
     &   .7225d-01,   .8099d-01,   .8930d-01,   .9714d-01,   .1044d+00,
     &   .1112d+00,   .1173d+00,   .1227d+00,   .1275d+00,   .1315d+00/

        data bes17  /
     &   .0000d+00,   .3646d-34,   .2388d-29,   .1567d-26,   .1562d-24,
     &   .5542d-23,   .1023d-21,   .1203d-20,   .1016d-19,   .6672d-19,
     &   .3590d-18,   .1644d-17,   .6593d-17,   .2364d-16,   .7704d-16,
     &   .2313d-15,   .6464d-15,   .1696d-14,   .4211d-14,   .9944d-14,
     &   .2246d-13,   .4871d-13,   .1018d-12,   .2059d-12,   .4039d-12,
     &   .7703d-12,   .1431d-11,   .2596d-11,   .4606d-11,   .8004d-11,
     &   .1364d-10,   .2283d-10,   .3757d-10,   .6085d-10,   .9709d-10,
     &   .1527d-09,   .2370d-09,   .3632d-09,   .5500d-09,   .8234d-09,
     &   .1219d-08,   .1787d-08,   .2594d-08,   .3729d-08,   .5313d-08,
     &   .7505d-08,   .1052d-07,   .1462d-07,   .2017d-07,   .2762d-07,
     &   .3756d-07,   .5074d-07,   .6810d-07,   .9083d-07,   .1204d-06,
     &   .1587d-06,   .2080d-06,   .2711d-06,   .3516d-06,   .4536d-06,
     &   .5822d-06,   .7439d-06,   .9459d-06,   .1197d-05,   .1509d-05,
     &   .1894d-05,   .2367d-05,   .2947d-05,   .3654d-05,   .4513d-05,
     &   .5554d-05,   .6810d-05,   .8321d-05,   .1013d-04,   .1230d-04,
     &   .1487d-04,   .1793d-04,   .2155d-04,   .2582d-04,   .3084d-04,
     &   .3673d-04,   .4363d-04,   .5166d-04,   .6102d-04,   .7186d-04,
     &   .8442d-04,   .9890d-04,   .1156d-03,   .1347d-03,   .1567d-03,
     &   .1817d-03,   .2103d-03,   .2428d-03,   .2796d-03,   .3213d-03,
     &   .3685d-03,   .4216d-03,   .4813d-03,   .5483d-03,   .6233d-03,
     &   .7072d-03,   .8006d-03,   .9046d-03,   .1020d-02,   .1148d-02,
     &   .1289d-02,   .1446d-02,   .1618d-02,   .1807d-02,   .2014d-02,
     &   .2241d-02,   .2489d-02,   .2760d-02,   .3055d-02,   .3375d-02,
     &   .3722d-02,   .4098d-02,   .4504d-02,   .4942d-02,   .5413d-02,
     &   .5920d-02,   .6463d-02,   .7044d-02,   .7664d-02,   .8326d-02,
     &   .9030d-02,   .9778d-02,   .1057d-01,   .1141d-01,   .1229d-01,
     &   .1323d-01,   .1421d-01,   .1524d-01,   .1632d-01,   .1744d-01,
     &   .1862d-01,   .1984d-01,   .2111d-01,   .2242d-01,   .2378d-01,
     &   .2518d-01,   .2662d-01,   .2810d-01,   .2961d-01,   .3115d-01,
     &   .3272d-01,   .3430d-01,   .3591d-01,   .3753d-01,   .3915d-01,
     &   .4077d-01,   .4239d-01,   .4399d-01,   .4557d-01,   .4712d-01,
     &   .4863d-01,   .5010d-01,   .5151d-01,   .5285d-01,   .5412d-01,
     &   .5530d-01,   .5639d-01,   .5737d-01,   .5823d-01,   .5897d-01,
     &   .5956d-01,   .6001d-01,   .6030d-01,   .6042d-01,   .6035d-01,
     &   .6010d-01,   .5964d-01,   .5898d-01,   .5809d-01,   .5698d-01,
     &   .5564d-01,   .5406d-01,   .5223d-01,   .5015d-01,   .4782d-01,
     &   .4523d-01,   .4238d-01,   .3927d-01,   .3590d-01,   .3228d-01,
     &   .2840d-01,   .2428d-01,   .1992d-01,   .1532d-01,   .1050d-01,
     &   .5468d-02,   .2378d-03,  -.5178d-02,  -.1076d-01,  -.1650d-01,
     &  -.2237d-01,  -.2835d-01,  -.3442d-01,  -.4057d-01,  -.4675d-01,
     &  -.5296d-01,  -.5915d-01,  -.6531d-01,  -.7140d-01,  -.7741d-01,
     &  -.8329d-01,  -.8901d-01,  -.9455d-01,  -.9988d-01,  -.1050d+00,
     &  -.1098d+00,  -.1143d+00,  -.1185d+00,  -.1223d+00,  -.1257d+00,
     &  -.1287d+00,  -.1313d+00,  -.1334d+00,  -.1350d+00,  -.1361d+00,
     &  -.1367d+00,  -.1368d+00,  -.1363d+00,  -.1352d+00,  -.1336d+00,
     &  -.1314d+00,  -.1287d+00,  -.1253d+00,  -.1214d+00,  -.1170d+00,
     &  -.1120d+00,  -.1064d+00,  -.1004d+00,  -.9384d-01,  -.8682d-01,
     &  -.7937d-01,  -.7150d-01,  -.6325d-01,  -.5465d-01,  -.4576d-01,
     &  -.3660d-01,  -.2722d-01,  -.1767d-01,  -.7988d-02,   .1768d-02,
     &   .1155d-01,   .2131d-01,   .3099d-01,   .4054d-01,   .4990d-01/

        data bes18  /
     &   .0000d+00,   .1072d-36,   .1405d-31,   .1383d-28,   .1838d-26,
     &   .8153d-25,   .1806d-23,   .2477d-22,   .2393d-21,   .1768d-20,
     &   .1057d-19,   .5326d-19,   .2330d-18,   .9051d-18,   .3178d-17,
     &   .1022d-16,   .3049d-16,   .8505d-16,   .2236d-15,   .5575d-15,
     &   .1326d-14,   .3020d-14,   .6619d-14,   .1400d-13,   .2866d-13,
     &   .5696d-13,   .1101d-12,   .2076d-12,   .3821d-12,   .6880d-12,
     &   .1214d-11,   .2100d-11,   .3570d-11,   .5966d-11,   .9812d-11,
     &   .1590d-10,   .2540d-10,   .4003d-10,   .6229d-10,   .9578d-10,
     &   .1456d-09,   .2189d-09,   .3256d-09,   .4797d-09,   .7000d-09,
     &   .1012d-08,   .1451d-08,   .2063d-08,   .2909d-08,   .4070d-08,
     &   .5654d-08,   .7799d-08,   .1068d-07,   .1454d-07,   .1966d-07,
     &   .2642d-07,   .3529d-07,   .4688d-07,   .6192d-07,   .8135d-07,
     &   .1063d-06,   .1383d-06,   .1789d-06,   .2305d-06,   .2955d-06,
     &   .3771d-06,   .4792d-06,   .6064d-06,   .7642d-06,   .9592d-06,
     &   .1199d-05,   .1494d-05,   .1854d-05,   .2292d-05,   .2824d-05,
     &   .3467d-05,   .4242d-05,   .5174d-05,   .6291d-05,   .7624d-05,
     &   .9212d-05,   .1110d-04,   .1333d-04,   .1596d-04,   .1906d-04,
     &   .2270d-04,   .2696d-04,   .3193d-04,   .3773d-04,   .4446d-04,
     &   .5226d-04,   .6128d-04,   .7168d-04,   .8364d-04,   .9737d-04,
     &   .1131d-03,   .1311d-03,   .1516d-03,   .1749d-03,   .2013d-03,
     &   .2313d-03,   .2651d-03,   .3033d-03,   .3463d-03,   .3945d-03,
     &   .4486d-03,   .5092d-03,   .5767d-03,   .6520d-03,   .7358d-03,
     &   .8287d-03,   .9317d-03,   .1046d-02,   .1171d-02,   .1310d-02,
     &   .1462d-02,   .1629d-02,   .1812d-02,   .2013d-02,   .2231d-02,
     &   .2469d-02,   .2729d-02,   .3010d-02,   .3315d-02,   .3645d-02,
     &   .4001d-02,   .4385d-02,   .4799d-02,   .5243d-02,   .5719d-02,
     &   .6229d-02,   .6774d-02,   .7355d-02,   .7973d-02,   .8631d-02,
     &   .9329d-02,   .1007d-01,   .1085d-01,   .1167d-01,   .1254d-01,
     &   .1345d-01,   .1441d-01,   .1541d-01,   .1645d-01,   .1754d-01,
     &   .1868d-01,   .1986d-01,   .2108d-01,   .2234d-01,   .2364d-01,
     &   .2498d-01,   .2635d-01,   .2776d-01,   .2920d-01,   .3066d-01,
     &   .3215d-01,   .3365d-01,   .3517d-01,   .3670d-01,   .3823d-01,
     &   .3976d-01,   .4128d-01,   .4279d-01,   .4428d-01,   .4573d-01,
     &   .4715d-01,   .4853d-01,   .4985d-01,   .5111d-01,   .5229d-01,
     &   .5340d-01,   .5442d-01,   .5534d-01,   .5615d-01,   .5684d-01,
     &   .5740d-01,   .5782d-01,   .5809d-01,   .5821d-01,   .5815d-01,
     &   .5792d-01,   .5751d-01,   .5690d-01,   .5608d-01,   .5506d-01,
     &   .5382d-01,   .5236d-01,   .5066d-01,   .4874d-01,   .4658d-01,
     &   .4417d-01,   .4153d-01,   .3864d-01,   .3551d-01,   .3214d-01,
     &   .2853d-01,   .2469d-01,   .2062d-01,   .1633d-01,   .1183d-01,
     &   .7119d-02,   .2218d-02,  -.2862d-02,  -.8109d-02,  -.1351d-01,
     &  -.1904d-01,  -.2469d-01,  -.3043d-01,  -.3625d-01,  -.4213d-01,
     &  -.4804d-01,  -.5395d-01,  -.5985d-01,  -.6570d-01,  -.7149d-01,
     &  -.7717d-01,  -.8274d-01,  -.8814d-01,  -.9337d-01,  -.9839d-01,
     &  -.1032d+00,  -.1077d+00,  -.1119d+00,  -.1158d+00,  -.1194d+00,
     &  -.1225d+00,  -.1253d+00,  -.1277d+00,  -.1296d+00,  -.1311d+00,
     &  -.1321d+00,  -.1326d+00,  -.1326d+00,  -.1320d+00,  -.1310d+00,
     &  -.1294d+00,  -.1272d+00,  -.1245d+00,  -.1213d+00,  -.1176d+00,
     &  -.1133d+00,  -.1085d+00,  -.1032d+00,  -.9738d-01,  -.9112d-01,
     &  -.8440d-01,  -.7727d-01,  -.6974d-01,  -.6185d-01,  -.5363d-01/

        data bes19  /
     &   .0000d+00,   .2979d-39,   .7805d-34,   .1153d-30,   .2043d-28,
     &   .1133d-26,   .3011d-25,   .4819d-24,   .5319d-23,   .4422d-22,
     &   .2938d-21,   .1629d-20,   .7775d-20,   .3273d-19,   .1238d-18,
     &   .4268d-18,   .1358d-17,   .4025d-17,   .1121d-16,   .2951d-16,
     &   .7389d-16,   .1768d-15,   .4061d-15,   .8982d-15,   .1920d-14,
     &   .3976d-14,   .7997d-14,   .1566d-13,   .2991d-13,   .5580d-13,
     &   .1019d-12,   .1823d-12,   .3200d-12,   .5517d-12,   .9355d-12,
     &   .1561d-11,   .2567d-11,   .4160d-11,   .6654d-11,   .1051d-10,
     &   .1639d-10,   .2527d-10,   .3855d-10,   .5818d-10,   .8693d-10,
     &   .1287d-09,   .1887d-09,   .2742d-09,   .3953d-09,   .5652d-09,
     &   .8018d-09,   .1129d-08,   .1578d-08,   .2191d-08,   .3021d-08,
     &   .4139d-08,   .5636d-08,   .7627d-08,   .1026d-07,   .1373d-07,
     &   .1827d-07,   .2418d-07,   .3184d-07,   .4171d-07,   .5438d-07,
     &   .7058d-07,   .9118d-07,   .1173d-06,   .1502d-06,   .1915d-06,
     &   .2432d-06,   .3076d-06,   .3876d-06,   .4866d-06,   .6085d-06,
     &   .7583d-06,   .9416d-06,   .1165d-05,   .1437d-05,   .1767d-05,
     &   .2165d-05,   .2645d-05,   .3221d-05,   .3910d-05,   .4734d-05,
     &   .5714d-05,   .6878d-05,   .8256d-05,   .9883d-05,   .1180d-04,
     &   .1405d-04,   .1669d-04,   .1977d-04,   .2337d-04,   .2755d-04,
     &   .3240d-04,   .3802d-04,   .4451d-04,   .5199d-04,   .6059d-04,
     &   .7046d-04,   .8176d-04,   .9468d-04,   .1094d-03,   .1261d-03,
     &   .1452d-03,   .1667d-03,   .1911d-03,   .2186d-03,   .2495d-03,
     &   .2844d-03,   .3234d-03,   .3672d-03,   .4161d-03,   .4707d-03,
     &   .5315d-03,   .5991d-03,   .6741d-03,   .7572d-03,   .8491d-03,
     &   .9505d-03,   .1062d-02,   .1185d-02,   .1320d-02,   .1468d-02,
     &   .1630d-02,   .1807d-02,   .2000d-02,   .2210d-02,   .2438d-02,
     &   .2686d-02,   .2954d-02,   .3245d-02,   .3558d-02,   .3896d-02,
     &   .4259d-02,   .4650d-02,   .5069d-02,   .5517d-02,   .5996d-02,
     &   .6508d-02,   .7053d-02,   .7632d-02,   .8247d-02,   .8899d-02,
     &   .9589d-02,   .1032d-01,   .1109d-01,   .1189d-01,   .1274d-01,
     &   .1363d-01,   .1456d-01,   .1554d-01,   .1655d-01,   .1761d-01,
     &   .1870d-01,   .1984d-01,   .2102d-01,   .2223d-01,   .2348d-01,
     &   .2476d-01,   .2607d-01,   .2742d-01,   .2879d-01,   .3018d-01,
     &   .3160d-01,   .3303d-01,   .3447d-01,   .3592d-01,   .3737d-01,
     &   .3881d-01,   .4025d-01,   .4167d-01,   .4307d-01,   .4445d-01,
     &   .4578d-01,   .4708d-01,   .4832d-01,   .4950d-01,   .5062d-01,
     &   .5166d-01,   .5262d-01,   .5348d-01,   .5424d-01,   .5489d-01,
     &   .5542d-01,   .5582d-01,   .5608d-01,   .5619d-01,   .5615d-01,
     &   .5594d-01,   .5556d-01,   .5499d-01,   .5424d-01,   .5329d-01,
     &   .5214d-01,   .5078d-01,   .4921d-01,   .4742d-01,   .4541d-01,
     &   .4317d-01,   .4071d-01,   .3802d-01,   .3510d-01,   .3195d-01,
     &   .2858d-01,   .2499d-01,   .2118d-01,   .1716d-01,   .1294d-01,
     &   .8523d-02,   .3919d-02,  -.8595d-03,  -.5801d-02,  -.1089d-01,
     &  -.1612d-01,  -.2146d-01,  -.2690d-01,  -.3243d-01,  -.3802d-01,
     &  -.4365d-01,  -.4930d-01,  -.5495d-01,  -.6057d-01,  -.6614d-01,
     &  -.7164d-01,  -.7704d-01,  -.8231d-01,  -.8742d-01,  -.9236d-01,
     &  -.9709d-01,  -.1016d+00,  -.1058d+00,  -.1098d+00,  -.1134d+00,
     &  -.1167d+00,  -.1197d+00,  -.1223d+00,  -.1244d+00,  -.1262d+00,
     &  -.1275d+00,  -.1284d+00,  -.1288d+00,  -.1287d+00,  -.1281d+00,
     &  -.1270d+00,  -.1254d+00,  -.1233d+00,  -.1207d+00,  -.1175d+00/

        data bes20  /
     &   .0000d+00,   .7839d-42,   .4108d-36,   .9101d-33,   .2150d-30,
     &   .1490d-28,   .4755d-27,   .8879d-26,   .1120d-24,   .1048d-23,
     &   .7738d-23,   .4719d-22,   .2458d-21,   .1121d-20,   .4567d-20,
     &   .1687d-19,   .5728d-19,   .1805d-18,   .5322d-18,   .1479d-17,
     &   .3900d-17,   .9803d-17,   .2359d-16,   .5457d-16,   .1218d-15,
     &   .2628d-15,   .5499d-15,   .1119d-14,   .2216d-14,   .4284d-14,
     &   .8096d-14,   .1497d-13,   .2715d-13,   .4829d-13,   .8441d-13,
     &   .1451d-12,   .2455d-12,   .4092d-12,   .6724d-12,   .1090d-11,
     &   .1746d-11,   .2761d-11,   .4316d-11,   .6674d-11,   .1021d-10,
     &   .1546d-10,   .2320d-10,   .3448d-10,   .5079d-10,   .7418d-10,
     &   .1075d-09,   .1545d-09,   .2203d-09,   .3120d-09,   .4387d-09,
     &   .6128d-09,   .8502d-09,   .1172d-08,   .1606d-08,   .2188d-08,
     &   .2963d-08,   .3991d-08,   .5347d-08,   .7125d-08,   .9447d-08,
     &   .1246d-07,   .1637d-07,   .2139d-07,   .2783d-07,   .3605d-07,
     &   .4650d-07,   .5973d-07,   .7642d-07,   .9737d-07,   .1236d-06,
     &   .1563d-06,   .1969d-06,   .2472d-06,   .3092d-06,   .3855d-06,
     &   .4790d-06,   .5933d-06,   .7324d-06,   .9014d-06,   .1106d-05,
     &   .1353d-05,   .1650d-05,   .2007d-05,   .2433d-05,   .2943d-05,
     &   .3549d-05,   .4269d-05,   .5122d-05,   .6129d-05,   .7316d-05,
     &   .8711d-05,   .1035d-04,   .1226d-04,   .1450d-04,   .1710d-04,
     &   .2012d-04,   .2362d-04,   .2768d-04,   .3236d-04,   .3775d-04,
     &   .4394d-04,   .5105d-04,   .5918d-04,   .6848d-04,   .7908d-04,
     &   .9113d-04,   .1048d-03,   .1204d-03,   .1379d-03,   .1578d-03,
     &   .1801d-03,   .2053d-03,   .2335d-03,   .2652d-03,   .3007d-03,
     &   .3403d-03,   .3845d-03,   .4337d-03,   .4883d-03,   .5490d-03,
     &   .6161d-03,   .6904d-03,   .7724d-03,   .8628d-03,   .9622d-03,
     &   .1071d-02,   .1191d-02,   .1322d-02,   .1466d-02,   .1622d-02,
     &   .1793d-02,   .1978d-02,   .2180d-02,   .2399d-02,   .2635d-02,
     &   .2891d-02,   .3167d-02,   .3465d-02,   .3785d-02,   .4129d-02,
     &   .4498d-02,   .4894d-02,   .5316d-02,   .5767d-02,   .6248d-02,
     &   .6760d-02,   .7303d-02,   .7880d-02,   .8490d-02,   .9136d-02,
     &   .9817d-02,   .1053d-01,   .1129d-01,   .1208d-01,   .1291d-01,
     &   .1378d-01,   .1469d-01,   .1564d-01,   .1662d-01,   .1764d-01,
     &   .1871d-01,   .1980d-01,   .2094d-01,   .2210d-01,   .2330d-01,
     &   .2453d-01,   .2579d-01,   .2708d-01,   .2839d-01,   .2972d-01,
     &   .3107d-01,   .3243d-01,   .3380d-01,   .3517d-01,   .3655d-01,
     &   .3792d-01,   .3928d-01,   .4063d-01,   .4196d-01,   .4325d-01,
     &   .4452d-01,   .4574d-01,   .4691d-01,   .4803d-01,   .4908d-01,
     &   .5006d-01,   .5096d-01,   .5178d-01,   .5249d-01,   .5311d-01,
     &   .5361d-01,   .5398d-01,   .5423d-01,   .5434d-01,   .5431d-01,
     &   .5411d-01,   .5376d-01,   .5324d-01,   .5254d-01,   .5166d-01,
     &   .5059d-01,   .4932d-01,   .4786d-01,   .4619d-01,   .4431d-01,
     &   .4222d-01,   .3992d-01,   .3740d-01,   .3467d-01,   .3172d-01,
     &   .2856d-01,   .2520d-01,   .2162d-01,   .1785d-01,   .1388d-01,
     &   .9723d-02,   .5387d-02,   .8812d-03,  -.3783d-02,  -.8595d-02,
     &  -.1354d-01,  -.1860d-01,  -.2377d-01,  -.2902d-01,  -.3435d-01,
     &  -.3972d-01,  -.4512d-01,  -.5054d-01,  -.5594d-01,  -.6131d-01,
     &  -.6662d-01,  -.7185d-01,  -.7697d-01,  -.8197d-01,  -.8681d-01,
     &  -.9147d-01,  -.9593d-01,  -.1002d+00,  -.1041d+00,  -.1079d+00,
     &  -.1113d+00,  -.1143d+00,  -.1171d+00,  -.1195d+00,  -.1214d+00/

        data bes21  /
     &   .0000d+00,   .1960d-44,   .2054d-38,   .6826d-35,   .2151d-32,
     &   .1863d-30,   .7134d-29,   .1554d-27,   .2242d-26,   .2359d-25,
     &   .1936d-24,   .1299d-23,   .7380d-23,   .3647d-22,   .1600d-21,
     &   .6337d-21,   .2295d-20,   .7684d-20,   .2400d-19,   .7043d-19,
     &   .1955d-18,   .5161d-18,   .1302d-17,   .3149d-17,   .7333d-17,
     &   .1649d-16,   .3590d-16,   .7586d-16,   .1559d-15,   .3123d-15,
     &   .6108d-15,   .1168d-14,   .2186d-14,   .4013d-14,   .7230d-14,
     &   .1280d-13,   .2228d-13,   .3819d-13,   .6449d-13,   .1074d-12,
     &   .1764d-12,   .2861d-12,   .4585d-12,   .7262d-12,   .1138d-11,
     &   .1763d-11,   .2706d-11,   .4111d-11,   .6189d-11,   .9234d-11,
     &   .1366d-10,   .2004d-10,   .2917d-10,   .4213d-10,   .6040d-10,
     &   .8598d-10,   .1216d-09,   .1707d-09,   .2382d-09,   .3304d-09,
     &   .4554d-09,   .6241d-09,   .8505d-09,   .1153d-08,   .1554d-08,
     &   .2084d-08,   .2782d-08,   .3694d-08,   .4883d-08,   .6425d-08,
     &   .8415d-08,   .1098d-07,   .1425d-07,   .1843d-07,   .2374d-07,
     &   .3046d-07,   .3893d-07,   .4958d-07,   .6290d-07,   .7951d-07,
     &   .1002d-06,   .1258d-06,   .1574d-06,   .1963d-06,   .2440d-06,
     &   .3024d-06,   .3737d-06,   .4603d-06,   .5654d-06,   .6925d-06,
     &   .8458d-06,   .1030d-05,   .1251d-05,   .1516d-05,   .1831d-05,
     &   .2207d-05,   .2653d-05,   .3182d-05,   .3806d-05,   .4542d-05,
     &   .5408d-05,   .6424d-05,   .7614d-05,   .9004d-05,   .1062d-04,
     &   .1251d-04,   .1470d-04,   .1723d-04,   .2016d-04,   .2354d-04,
     &   .2743d-04,   .3190d-04,   .3702d-04,   .4289d-04,   .4959d-04,
     &   .5723d-04,   .6593d-04,   .7581d-04,   .8702d-04,   .9970d-04,
     &   .1140d-03,   .1302d-03,   .1484d-03,   .1689d-03,   .1918d-03,
     &   .2176d-03,   .2463d-03,   .2785d-03,   .3143d-03,   .3541d-03,
     &   .3984d-03,   .4475d-03,   .5019d-03,   .5621d-03,   .6285d-03,
     &   .7017d-03,   .7823d-03,   .8708d-03,   .9680d-03,   .1074d-02,
     &   .1191d-02,   .1318d-02,   .1457d-02,   .1608d-02,   .1772d-02,
     &   .1950d-02,   .2143d-02,   .2352d-02,   .2578d-02,   .2822d-02,
     &   .3085d-02,   .3368d-02,   .3671d-02,   .3997d-02,   .4346d-02,
     &   .4720d-02,   .5118d-02,   .5543d-02,   .5996d-02,   .6477d-02,
     &   .6988d-02,   .7529d-02,   .8101d-02,   .8706d-02,   .9344d-02,
     &   .1002d-01,   .1072d-01,   .1146d-01,   .1224d-01,   .1305d-01,
     &   .1390d-01,   .1479d-01,   .1571d-01,   .1667d-01,   .1766d-01,
     &   .1868d-01,   .1974d-01,   .2084d-01,   .2196d-01,   .2312d-01,
     &   .2430d-01,   .2551d-01,   .2674d-01,   .2800d-01,   .2927d-01,
     &   .3055d-01,   .3185d-01,   .3316d-01,   .3447d-01,   .3578d-01,
     &   .3708d-01,   .3838d-01,   .3966d-01,   .4091d-01,   .4214d-01,
     &   .4334d-01,   .4449d-01,   .4560d-01,   .4666d-01,   .4765d-01,
     &   .4858d-01,   .4943d-01,   .5020d-01,   .5088d-01,   .5146d-01,
     &   .5193d-01,   .5229d-01,   .5253d-01,   .5264d-01,   .5261d-01,
     &   .5243d-01,   .5211d-01,   .5162d-01,   .5097d-01,   .5015d-01,
     &   .4915d-01,   .4796d-01,   .4659d-01,   .4503d-01,   .4327d-01,
     &   .4131d-01,   .3915d-01,   .3679d-01,   .3423d-01,   .3146d-01,
     &   .2850d-01,   .2533d-01,   .2197d-01,   .1841d-01,   .1467d-01,
     &   .1075d-01,   .6658d-02,   .2400d-02,  -.2012d-02,  -.6568d-02,
     &  -.1126d-01,  -.1606d-01,  -.2098d-01,  -.2598d-01,  -.3106d-01,
     &  -.3619d-01,  -.4136d-01,  -.4655d-01,  -.5174d-01,  -.5691d-01,
     &  -.6204d-01,  -.6711d-01,  -.7209d-01,  -.7696d-01,  -.8170d-01/

        data bes22  /
     &   .0000d+00,   .4666d-47,   .9782d-41,   .4876d-37,   .2048d-34,
     &   .2219d-32,   .1019d-30,   .2591d-29,   .4272d-28,   .5057d-27,
     &   .4612d-26,   .3404d-25,   .2110d-24,   .1130d-23,   .5341d-23,
     &   .2266d-22,   .8756d-22,   .3116d-21,   .1030d-20,   .3193d-20,
     &   .9333d-20,   .2587d-19,   .6838d-19,   .1730d-18,   .4205d-18,
     &   .9851d-18,   .2231d-17,   .4898d-17,   .1044d-16,   .2167d-16,
     &   .4386d-16,   .8669d-16,   .1676d-15,   .3173d-15,   .5893d-15,
     &   .1074d-14,   .1925d-14,   .3392d-14,   .5885d-14,   .1006d-13,
     &   .1696d-13,   .2821d-13,   .4634d-13,   .7518d-13,   .1206d-12,
     &   .1912d-12,   .3001d-12,   .4662d-12,   .7172d-12,   .1093d-11,
     &   .1651d-11,   .2472d-11,   .3670d-11,   .5407d-11,   .7903d-11,
     &   .1147d-10,   .1652d-10,   .2363d-10,   .3358d-10,   .4740d-10,
     &   .6650d-10,   .9273d-10,   .1285d-09,   .1772d-09,   .2428d-09,
     &   .3310d-09,   .4490d-09,   .6058d-09,   .8134d-09,   .1087d-08,
     &   .1446d-08,   .1914d-08,   .2523d-08,   .3311d-08,   .4328d-08,
     &   .5633d-08,   .7303d-08,   .9431d-08,   .1213d-07,   .1555d-07,
     &   .1986d-07,   .2527d-07,   .3205d-07,   .4050d-07,   .5102d-07,
     &   .6406d-07,   .8017d-07,   .1000d-06,   .1244d-06,   .1543d-06,
     &   .1908d-06,   .2353d-06,   .2893d-06,   .3547d-06,   .4338d-06,
     &   .5291d-06,   .6436d-06,   .7808d-06,   .9450d-06,   .1141d-05,
     &   .1374d-05,   .1651d-05,   .1979d-05,   .2367d-05,   .2824d-05,
     &   .3362d-05,   .3994d-05,   .4734d-05,   .5599d-05,   .6609d-05,
     &   .7785d-05,   .9151d-05,   .1074d-04,   .1257d-04,   .1469d-04,
     &   .1713d-04,   .1994d-04,   .2317d-04,   .2687d-04,   .3111d-04,
     &   .3595d-04,   .4147d-04,   .4775d-04,   .5489d-04,   .6299d-04,
     &   .7216d-04,   .8253d-04,   .9423d-04,   .1074d-03,   .1223d-03,
     &   .1389d-03,   .1576d-03,   .1785d-03,   .2019d-03,   .2280d-03,
     &   .2570d-03,   .2894d-03,   .3253d-03,   .3652d-03,   .4093d-03,
     &   .4581d-03,   .5119d-03,   .5713d-03,   .6367d-03,   .7086d-03,
     &   .7875d-03,   .8740d-03,   .9686d-03,   .1072d-02,   .1185d-02,
     &   .1308d-02,   .1442d-02,   .1587d-02,   .1745d-02,   .1916d-02,
     &   .2101d-02,   .2301d-02,   .2517d-02,   .2749d-02,   .2999d-02,
     &   .3268d-02,   .3556d-02,   .3865d-02,   .4195d-02,   .4548d-02,
     &   .4924d-02,   .5325d-02,   .5752d-02,   .6205d-02,   .6685d-02,
     &   .7194d-02,   .7731d-02,   .8299d-02,   .8898d-02,   .9528d-02,
     &   .1019d-01,   .1089d-01,   .1161d-01,   .1237d-01,   .1317d-01,
     &   .1400d-01,   .1486d-01,   .1576d-01,   .1669d-01,   .1765d-01,
     &   .1865d-01,   .1967d-01,   .2073d-01,   .2181d-01,   .2293d-01,
     &   .2406d-01,   .2523d-01,   .2641d-01,   .2761d-01,   .2883d-01,
     &   .3006d-01,   .3130d-01,   .3255d-01,   .3380d-01,   .3505d-01,
     &   .3629d-01,   .3752d-01,   .3874d-01,   .3993d-01,   .4110d-01,
     &   .4224d-01,   .4333d-01,   .4439d-01,   .4539d-01,   .4633d-01,
     &   .4721d-01,   .4802d-01,   .4875d-01,   .4939d-01,   .4994d-01,
     &   .5039d-01,   .5073d-01,   .5096d-01,   .5106d-01,   .5104d-01,
     &   .5088d-01,   .5057d-01,   .5012d-01,   .4951d-01,   .4874d-01,
     &   .4780d-01,   .4669d-01,   .4541d-01,   .4394d-01,   .4229d-01,
     &   .4045d-01,   .3842d-01,   .3620d-01,   .3378d-01,   .3118d-01,
     &   .2839d-01,   .2540d-01,   .2223d-01,   .1888d-01,   .1534d-01,
     &   .1164d-01,   .7763d-02,   .3731d-02,  -.4508d-03,  -.4774d-02,
     &  -.9227d-02,  -.1380d-01,  -.1848d-01,  -.2325d-01,  -.2809d-01/

        data bes23  /
     &   .0000d+00,   .1060d-49,   .4446d-43,   .3325d-39,   .1862d-36,
     &   .2521d-34,   .1390d-32,   .4124d-31,   .7769d-30,   .1035d-28,
     &   .1049d-27,   .8516d-27,   .5760d-26,   .3342d-25,   .1701d-24,
     &   .7735d-24,   .3189d-23,   .1206d-22,   .4223d-22,   .1382d-21,
     &   .4251d-21,   .1238d-20,   .3428d-20,   .9068d-20,   .2301d-19,
     &   .5616d-19,   .1323d-18,   .3017d-18,   .6674d-18,   .1435d-17,
     &   .3005d-17,   .6139d-17,   .1226d-16,   .2394d-16,   .4582d-16,
     &   .8602d-16,   .1586d-15,   .2874d-15,   .5123d-15,   .8992d-15,
     &   .1555d-14,   .2653d-14,   .4466d-14,   .7421d-14,   .1218d-13,
     &   .1977d-13,   .3174d-13,   .5040d-13,   .7923d-13,   .1233d-12,
     &   .1902d-12,   .2906d-12,   .4402d-12,   .6614d-12,   .9856d-12,
     &   .1457d-11,   .2139d-11,   .3116d-11,   .4509d-11,   .6479d-11,
     &   .9250d-11,   .1312d-10,   .1850d-10,   .2593d-10,   .3613d-10,
     &   .5006d-10,   .6899d-10,   .9458d-10,   .1290d-09,   .1750d-09,
     &   .2364d-09,   .3177d-09,   .4250d-09,   .5660d-09,   .7505d-09,
     &   .9910d-09,   .1303d-08,   .1706d-08,   .2226d-08,   .2892d-08,
     &   .3744d-08,   .4828d-08,   .6205d-08,   .7945d-08,   .1014d-07,
     &   .1289d-07,   .1635d-07,   .2065d-07,   .2601d-07,   .3267d-07,
     &   .4089d-07,   .5104d-07,   .6352d-07,   .7882d-07,   .9754d-07,
     &   .1204d-06,   .1481d-06,   .1818d-06,   .2226d-06,   .2718d-06,
     &   .3311d-06,   .4023d-06,   .4877d-06,   .5897d-06,   .7114d-06,
     &   .8563d-06,   .1028d-05,   .1232d-05,   .1473d-05,   .1757d-05,
     &   .2092d-05,   .2485d-05,   .2946d-05,   .3486d-05,   .4116d-05,
     &   .4850d-05,   .5704d-05,   .6696d-05,   .7845d-05,   .9174d-05,
     &   .1071d-04,   .1248d-04,   .1451d-04,   .1685d-04,   .1952d-04,
     &   .2259d-04,   .2609d-04,   .3008d-04,   .3462d-04,   .3979d-04,
     &   .4565d-04,   .5229d-04,   .5980d-04,   .6828d-04,   .7784d-04,
     &   .8861d-04,   .1007d-03,   .1143d-03,   .1295d-03,   .1465d-03,
     &   .1655d-03,   .1868d-03,   .2104d-03,   .2367d-03,   .2658d-03,
     &   .2982d-03,   .3341d-03,   .3737d-03,   .4175d-03,   .4657d-03,
     &   .5188d-03,   .5773d-03,   .6414d-03,   .7118d-03,   .7888d-03,
     &   .8731d-03,   .9650d-03,   .1065d-02,   .1175d-02,   .1293d-02,
     &   .1422d-02,   .1562d-02,   .1714d-02,   .1878d-02,   .2055d-02,
     &   .2246d-02,   .2452d-02,   .2673d-02,   .2911d-02,   .3167d-02,
     &   .3440d-02,   .3733d-02,   .4046d-02,   .4379d-02,   .4735d-02,
     &   .5114d-02,   .5516d-02,   .5943d-02,   .6395d-02,   .6874d-02,
     &   .7380d-02,   .7914d-02,   .8477d-02,   .9069d-02,   .9691d-02,
     &   .1034d-01,   .1103d-01,   .1174d-01,   .1249d-01,   .1327d-01,
     &   .1408d-01,   .1492d-01,   .1579d-01,   .1669d-01,   .1763d-01,
     &   .1859d-01,   .1959d-01,   .2061d-01,   .2166d-01,   .2273d-01,
     &   .2383d-01,   .2495d-01,   .2608d-01,   .2724d-01,   .2841d-01,
     &   .2959d-01,   .3078d-01,   .3197d-01,   .3317d-01,   .3436d-01,
     &   .3554d-01,   .3672d-01,   .3788d-01,   .3902d-01,   .4013d-01,
     &   .4121d-01,   .4225d-01,   .4325d-01,   .4420d-01,   .4510d-01,
     &   .4594d-01,   .4670d-01,   .4740d-01,   .4801d-01,   .4853d-01,
     &   .4896d-01,   .4928d-01,   .4950d-01,   .4960d-01,   .4958d-01,
     &   .4943d-01,   .4915d-01,   .4873d-01,   .4815d-01,   .4743d-01,
     &   .4655d-01,   .4551d-01,   .4429d-01,   .4291d-01,   .4136d-01,
     &   .3962d-01,   .3771d-01,   .3562d-01,   .3334d-01,   .3088d-01,
     &   .2824d-01,   .2542d-01,   .2243d-01,   .1925d-01,   .1591d-01/

        data bes24  /
     &   .0000d+00,   .2305d-52,   .1933d-45,   .2168d-41,   .1620d-38,
     &   .2741d-36,   .1814d-34,   .6277d-33,   .1352d-31,   .2026d-30,
     &   .2281d-29,   .2038d-28,   .1504d-27,   .9452d-27,   .5182d-26,
     &   .2525d-25,   .1110d-24,   .4462d-24,   .1655d-23,   .5716d-23,
     &   .1852d-22,   .5663d-22,   .1643d-21,   .4546d-21,   .1204d-20,
     &   .3062d-20,   .7504d-20,   .1777d-19,   .4078d-19,   .9084d-19,
     &   .1969d-18,   .4157d-18,   .8569d-18,   .1727d-17,   .3406d-17,
     &   .6585d-17,   .1249d-16,   .2327d-16,   .4262d-16,   .7682d-16,
     &   .1363d-15,   .2385d-15,   .4113d-15,   .7001d-15,   .1177d-14,
     &   .1954d-14,   .3207d-14,   .5207d-14,   .8362d-14,   .1329d-13,
     &   .2093d-13,   .3264d-13,   .5044d-13,   .7728d-13,   .1174d-12,
     &   .1769d-12,   .2645d-12,   .3925d-12,   .5782d-12,   .8457d-12,
     &   .1229d-11,   .1773d-11,   .2542d-11,   .3623d-11,   .5132d-11,
     &   .7227d-11,   .1012d-10,   .1409d-10,   .1952d-10,   .2690d-10,
     &   .3687d-10,   .5030d-10,   .6830d-10,   .9230d-10,   .1242d-09,
     &   .1663d-09,   .2217d-09,   .2944d-09,   .3894d-09,   .5129d-09,
     &   .6729d-09,   .8794d-09,   .1145d-08,   .1486d-08,   .1920d-08,
     &   .2474d-08,   .3176d-08,   .4063d-08,   .5181d-08,   .6586d-08,
     &   .8346d-08,   .1054d-07,   .1328d-07,   .1667d-07,   .2088d-07,
     &   .2607d-07,   .3245d-07,   .4029d-07,   .4990d-07,   .6162d-07,
     &   .7590d-07,   .9326d-07,   .1143d-06,   .1397d-06,   .1704d-06,
     &   .2073d-06,   .2517d-06,   .3048d-06,   .3683d-06,   .4440d-06,
     &   .5341d-06,   .6411d-06,   .7679d-06,   .9179d-06,   .1095d-05,
     &   .1303d-05,   .1548d-05,   .1836d-05,   .2172d-05,   .2566d-05,
     &   .3025d-05,   .3559d-05,   .4180d-05,   .4900d-05,   .5734d-05,
     &   .6698d-05,   .7810d-05,   .9092d-05,   .1057d-04,   .1226d-04,
     &   .1420d-04,   .1641d-04,   .1895d-04,   .2184d-04,   .2513d-04,
     &   .2887d-04,   .3312d-04,   .3793d-04,   .4337d-04,   .4953d-04,
     &   .5647d-04,   .6429d-04,   .7308d-04,   .8296d-04,   .9403d-04,
     &   .1064d-03,   .1203d-03,   .1358d-03,   .1531d-03,   .1723d-03,
     &   .1937d-03,   .2174d-03,   .2437d-03,   .2729d-03,   .3051d-03,
     &   .3407d-03,   .3800d-03,   .4232d-03,   .4708d-03,   .5230d-03,
     &   .5803d-03,   .6431d-03,   .7118d-03,   .7868d-03,   .8687d-03,
     &   .9579d-03,   .1055d-02,   .1161d-02,   .1275d-02,   .1399d-02,
     &   .1534d-02,   .1680d-02,   .1837d-02,   .2006d-02,   .2189d-02,
     &   .2385d-02,   .2596d-02,   .2823d-02,   .3065d-02,   .3325d-02,
     &   .3603d-02,   .3899d-02,   .4215d-02,   .4552d-02,   .4909d-02,
     &   .5289d-02,   .5692d-02,   .6119d-02,   .6570d-02,   .7047d-02,
     &   .7549d-02,   .8079d-02,   .8636d-02,   .9221d-02,   .9834d-02,
     &   .1048d-01,   .1115d-01,   .1185d-01,   .1258d-01,   .1334d-01,
     &   .1413d-01,   .1495d-01,   .1581d-01,   .1669d-01,   .1759d-01,
     &   .1853d-01,   .1949d-01,   .2048d-01,   .2150d-01,   .2253d-01,
     &   .2359d-01,   .2467d-01,   .2576d-01,   .2687d-01,   .2800d-01,
     &   .2913d-01,   .3027d-01,   .3142d-01,   .3256d-01,   .3370d-01,
     &   .3484d-01,   .3596d-01,   .3707d-01,   .3815d-01,   .3921d-01,
     &   .4024d-01,   .4124d-01,   .4219d-01,   .4310d-01,   .4395d-01,
     &   .4475d-01,   .4548d-01,   .4613d-01,   .4672d-01,   .4721d-01,
     &   .4762d-01,   .4793d-01,   .4814d-01,   .4824d-01,   .4823d-01,
     &   .4809d-01,   .4782d-01,   .4742d-01,   .4689d-01,   .4621d-01,
     &   .4537d-01,   .4439d-01,   .4325d-01,   .4194d-01,   .4048d-01/

        data bes25  /
     &   .0000d+00,   .4803d-55,   .8055d-48,   .1355d-43,   .1350d-40,
     &   .2856d-38,   .2268d-36,   .9156d-35,   .2253d-33,   .3799d-32,
     &   .4754d-31,   .4672d-30,   .3762d-29,   .2562d-28,   .1513d-27,
     &   .7899d-27,   .3706d-26,   .1582d-25,   .6215d-25,   .2266d-24,
     &   .7731d-24,   .2483d-23,   .7548d-23,   .2183d-22,   .6034d-22,
     &   .1599d-21,   .4077d-21,   .1003d-20,   .2387d-20,   .5509d-20,
     &   .1235d-19,   .2696d-19,   .5739d-19,   .1193d-18,   .2425d-18,
     &   .4828d-18,   .9423d-18,   .1805d-17,   .3396d-17,   .6284d-17,
     &   .1144d-16,   .2052d-16,   .3628d-16,   .6325d-16,   .1088d-15,
     &   .1849d-15,   .3104d-15,   .5150d-15,   .8451d-15,   .1372d-14,
     &   .2205d-14,   .3509d-14,   .5533d-14,   .8643d-14,   .1339d-13,
     &   .2055d-13,   .3131d-13,   .4731d-13,   .7095d-13,   .1056d-12,
     &   .1561d-12,   .2292d-12,   .3343d-12,   .4843d-12,   .6973d-12,
     &   .9979d-12,   .1420d-11,   .2008d-11,   .2825d-11,   .3953d-11,
     &   .5501d-11,   .7617d-11,   .1049d-10,   .1439d-10,   .1964d-10,
     &   .2667d-10,   .3607d-10,   .4856d-10,   .6511d-10,   .8692d-10,
     &   .1156d-09,   .1531d-09,   .2019d-09,   .2654d-09,   .3474d-09,
     &   .4533d-09,   .5892d-09,   .7633d-09,   .9855d-09,   .1268d-08,
     &   .1626d-08,   .2079d-08,   .2650d-08,   .3367d-08,   .4265d-08,
     &   .5387d-08,   .6784d-08,   .8520d-08,   .1067d-07,   .1332d-07,
     &   .1660d-07,   .2062d-07,   .2554d-07,   .3157d-07,   .3892d-07,
     &   .4786d-07,   .5871d-07,   .7185d-07,   .8773d-07,   .1069d-06,
     &   .1299d-06,   .1575d-06,   .1906d-06,   .2302d-06,   .2773d-06,
     &   .3334d-06,   .4001d-06,   .4791d-06,   .5725d-06,   .6828d-06,
     &   .8127d-06,   .9655d-06,   .1145d-05,   .1355d-05,   .1601d-05,
     &   .1888d-05,   .2222d-05,   .2611d-05,   .3063d-05,   .3586d-05,
     &   .4192d-05,   .4892d-05,   .5699d-05,   .6629d-05,   .7698d-05,
     &   .8924d-05,   .1033d-04,   .1194d-04,   .1377d-04,   .1587d-04,
     &   .1825d-04,   .2097d-04,   .2405d-04,   .2754d-04,   .3149d-04,
     &   .3596d-04,   .4100d-04,   .4669d-04,   .5308d-04,   .6027d-04,
     &   .6834d-04,   .7738d-04,   .8751d-04,   .9882d-04,   .1114d-03,
     &   .1255d-03,   .1412d-03,   .1586d-03,   .1779d-03,   .1994d-03,
     &   .2231d-03,   .2494d-03,   .2783d-03,   .3103d-03,   .3455d-03,
     &   .3843d-03,   .4269d-03,   .4736d-03,   .5248d-03,   .5808d-03,
     &   .6421d-03,   .7090d-03,   .7820d-03,   .8614d-03,   .9479d-03,
     &   .1042d-02,   .1144d-02,   .1254d-02,   .1374d-02,   .1503d-02,
     &   .1643d-02,   .1793d-02,   .1955d-02,   .2130d-02,   .2317d-02,
     &   .2518d-02,   .2734d-02,   .2965d-02,   .3211d-02,   .3475d-02,
     &   .3756d-02,   .4055d-02,   .4374d-02,   .4712d-02,   .5071d-02,
     &   .5452d-02,   .5855d-02,   .6280d-02,   .6730d-02,   .7204d-02,
     &   .7703d-02,   .8227d-02,   .8778d-02,   .9356d-02,   .9961d-02,
     &   .1059d-01,   .1125d-01,   .1194d-01,   .1266d-01,   .1340d-01,
     &   .1418d-01,   .1498d-01,   .1581d-01,   .1666d-01,   .1755d-01,
     &   .1846d-01,   .1939d-01,   .2035d-01,   .2133d-01,   .2233d-01,
     &   .2335d-01,   .2439d-01,   .2545d-01,   .2652d-01,   .2760d-01,
     &   .2869d-01,   .2979d-01,   .3089d-01,   .3198d-01,   .3308d-01,
     &   .3416d-01,   .3524d-01,   .3630d-01,   .3734d-01,   .3835d-01,
     &   .3933d-01,   .4028d-01,   .4119d-01,   .4206d-01,   .4287d-01,
     &   .4363d-01,   .4433d-01,   .4496d-01,   .4551d-01,   .4599d-01,
     &   .4638d-01,   .4667d-01,   .4687d-01,   .4697d-01,   .4696d-01/

        do 200 i=1,250
200       bess1(1,i)=bes1(i)
        do 210 i=1,250
210       bess1(2,i)=bes2(i)
        do 220 i=1,250
220       bess1(3,i)=bes3(i)
        do 230 i=1,250
230       bess1(4,i)=bes4(i)
        do 240 i=1,250
240       bess1(5,i)=bes5(i)
        do 250 i=1,250
250       bess1(6,i)=bes6(i)
        do 260 i=1,250
260       bess1(7,i)=bes7(i)
        do 270 i=1,250
270       bess1(8,i)=bes8(i)
        do 280 i=1,250
280       bess1(9,i)=bes9(i)
        do 290 i=1,250
290       bess1(10,i)=bes10(i)
        do 300 i=1,250
300       bess1(11,i)=bes11(i)
        do 310 i=1,250
310       bess1(12,i)=bes12(i)
        do 320 i=1,250
320       bess1(13,i)=bes13(i)
        do 330 i=1,250
330       bess1(14,i)=bes14(i)
 	do 340 i=1,250
340       bess1(15,i)=bes15(i)
        do 350 i=1,250
350       bess1(16,i)=bes16(i)
        do 360 i=1,250
360       bess1(17,i)=bes17(i)
        do 370 i=1,250
370       bess1(18,i)=bes18(i)
        do 380 i=1,250
380       bess1(19,i)=bes19(i)
        do 400 i=1,250
400       bess1(20,i)=bes20(i)
        do 410 i=1,250
410       bess1(21,i)=bes21(i)
        do 420 i=1,250
420       bess1(22,i)=bes22(i)
        do 430 i=1,250
430       bess1(23,i)=bes23(i)
        do 440 i=1,250
440       bess1(24,i)=bes24(i)
        do 450 i=1,250
450       bess1(25,i)=bes25(i)
	
        end
