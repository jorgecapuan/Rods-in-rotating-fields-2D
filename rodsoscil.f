	Program Nanorods_2D_Jorge
	implicit none
	include 'constants.f'
!	include 'mpif.h'
	integer*4 seed(3)
        integer step,k,l,i,j,n,countrat(nmol),Nbonds(nmol),LL(nmol)
!	integer nmolaux,ndimaux,NAaux
	integer rbin,h,sizeHistRdf,countRdf,countpol,countt,countpol2
	
        Real*8 r(ndim,nmol),y_box,s(ndim,nmol),Polim2
        Real*8 vmagli,mominert,vmagrot,u(ndim,nmol),ddx,ddy
	Real*8 rsx(NA,nmol),rsy(NA,nmol),D(NA),usum
	Real*8 ra(ndim,nmol),kienergy,totenergy,potenergy
	Real*8 skienergy,sig,MomLi,tempinst,pi,toten,kinen,poten
	Real*8 FSX(NA,nmol),FSY(NA,nmol),g(ndim,nmol),sumaveg(7)
	Real*8 sm(ndim,nmol),sa(ndim,nmol),Pol2sqr,NemOrdsqr
	Real*8 rv(ndim,nmol),Pol,Pol2,NemOrd,NemOrd2,B(ndim,nmol)
	Real*8 histRdf(100000),histRdf2(100000),eta,sumveg2(7)	
	Real*8 ratkvid(nmol),sckvid(nmol),Polimkvi,nclukvi 
	Real*8 histCom(100000),histCom2(100000),kk,tt,nbondcor(60000)
	Real*8 histang(100000),histang2(100000),Gtst,Gtst2
	Real*8 totenn,potenn,kinenn,NemOrd2sqr,Polsqr
	integer nTimeWin1,TimeWindow1(10000),cDiff1(10000),nstep
	integer nTimeWin2,TimeWindow2(10000),cDiff2(10000)
	Real*8 xp(10000,nmol),yp(10000,nmol),Diffx(10000)
	Real*8 Diffy(10000),rcross(ndim,nmol),magcor(60000),Magtotal
        Real*8 spx(60000,nmol),spy(60000,nmol),corrmi(60000)
	Real*8 magcorr(10000),nbondaver,kk1,jj1,kk2,jj2
        REal*8 gamape,gamapa,gamarot,corpe,corpa,corrot,facrot
	integer switchDiff
	integer kvid,countRdf2,m,nsamples, countcoran, countcoran2
    
        integer switch,nsamp,ntel,ntime(60000),t0,tt0,time0(60000)
        integer t,t0max,tmax
        Real*8 dtime,dt,sx0(nmol,60000),Difftheta(10000)
        Real*8 delt,sy0(nmol,60000),facpe,facpa,VXI,VYI,vacf(60000)
        Real*8 thetap(10000,nmol),Diffpa(10000),Diffpe(10000)
        

	character (len=13) nome_arq_1
        character (len=2) nome_arq_2
        character (len=15) nome_arq_total
        character (len=18) nome_arq_3
        character (len=20) nome_arq_total2
	COMMON /aleat/seed
	 Common /link1/MomLi
!	Common /link2/ra
	Common /link2/step
	common /angle/ countcoran,countcoran2
        seed(1) = 96
	seed(2) = 9687
	seed(3) = 4587999
	
	open(2,file='Initialcoord.dat')
	open(8,file='meanvalues.dat')
	open(4,file='energy.dat')
	open(44,file='blender.xyz')
       open(632,file='nclustaver.dat')

	open(1510,file='nbondaver.dat')
	
	y_box = x_box
!	pi = 4.d0*atan(1.d0)	
!	nmol = dint((4.d0*eta*x_box*y_box)/(pi*NA*sigma*sigma)) + 1
	

!########## Counters ##################################
	kvid = 0 
	countRdf = 0
	countRdf2 = 0
	countpol = 0
	countpol2 = 0
	countt = 0
	countcoran = 0
	countcoran2 = 0
	toten =0
	kinen = 0
	poten = 0
	Polimkvi = 0.d0
	Polim2 = 0.d0
	Gtst = 0.d0
	Gtst2 = 0.d0
	nclukvi = 0 
        t0 = 0
        switchDiff = 1
        switch = 0
	nsamples = final_step-stepaverage
 !       tempreq = 1.d-1
        nsamp = 50
!        Bext = 0.d0
        
        
	Do k=0,6

	sumaveg(k) = 0.d0
	sumveg2(k) = 0.d0
	
	enddo

	Do K = 1,nmol
	sckvid(K) = 0	
	ratkvid(K) = 0
	countrat(k) = 0
	enddo
        
!#####################################################
	
	vmagli = dsqrt(ndim*temp) !This is for 2D 

!######### Moment of Inertia #########################	
	mominert = 0.d0 
	
	Do I = 1,NA
	 mominert = mominert + 0.25d0*(-(NA-1.d0)+2.d0*(I-1.d0))**2	
	enddo
!	write(*,*)mominert
	mominert = mominert + 0.125d0*NA	! Actually this is the adimensional moment of inertia for a rod which is I/m*sigma²
!	write(*,*)mominert	
	mominert = mominert/NA
!	write(*,*)mominert
!	vmagrot=dsqrt(mass*vmag*(sigma)**2/(2*mominert))   !this is the general expression for moment of inertia before the adimensionalization
!	mominert = (1./12.)*NA**2		
	vmagrot = vmagli/dsqrt(2.d0*mominert) !this is the efective moment of inertia after the adimesionalization
!####################################################################
!	skienergy = 0.d0
 	
	Call PosiInit(y_box,u,r,s,D,rsx,rsy)
	write(*,*) 'Number of rods=',nmol
	write(*,*) 'Number of beads=',NA
	write(*,*) 'Simulation box = ',x_box
	write(*,*) 'B =',B0
	write(*,*) 'omega=',omega
         Call ConvertPosi (D,r,s,rsx,rsy)

        pi = 4.d0*atan(1.d0)
!         Call Polarization(nmol,NA,ndim,countpol,countpol2,s)


	Call Initangvel(vmagrot,mominert,u,g)
        Call InitLinearVels(vmagli,rv) 

!################## Center of Mass calculation ########################
	
	Do i = 1,nmol
	r(1,i) = rsx(1,i) + dble(NA-1)*sigma*0.5d0*s(1,i)
        r(2,i) = rsy(1,i) + dble(NA-1)*sigma*0.5d0*s(2,i)
       
	enddo
         Call  ApplyPBCR(switchDiff,y_box,r,rcross)
!#####################################################################
        
	 Do i=1,nmol
	  Do j=1,NA		
	Write(2,*) rsx(j,i),rsy(j,i)
      	  enddo
	 enddo
!	write(3,'(I5,e20.8)') l,dsqrt(rv(1,l)**2+rv(2,l)**2.)
        step = 0
	Call Force (usum,rsx,rsy,FSX,FSY,B,s)

	Call ConverForTor (mominert,D,r,ra,s,sa,FSX,FSY,u,B,rv,1)
        corpe = 0.839d0 + 0.185d0/dble(NA) + 0.233d0/dble(NA)**2 
        corpa = - 0.207d0 + 0.987d0/dble(NA) - 0.133d0/dble(NA)**2
       corrot = - 0.622d0 + 0.917d0/dble(NA) - 0.050d0/dble(NA)**2 
        gamape = mult*4.d0*pi*NA/(log(dble(NA))+corpe)
        gamapa = mult*2.d0*pi*NA/(log(dble(NA))+corpa)
       gamarot = mult*NA*NA*NA*pi/(3.d0*log(dble(NA))+corrot)
        facpe = 1.d0/(1.d0+gamape*deltaT*0.5d0/dble(NA))
         facpa = 1.d0/(1.d0+gamapa*deltaT*0.5d0/dble(NA))
         facrot = 1.d0/(1.d0+gamarot*deltaT*0.5d0/mominert)

	Do step =1,final_step
 !       If(step.eq.1) then
  !      tempreq = temp
  !      Bext = B0
  !      vmagli = dsqrt(ndim*tempreq)
  !      vmagrot = vmagli/dsqrt(2.d0*mominert)	
!	open(1602,file='randomconfig.dat')
!	Do i=1,nmol
!	  Do j=1,NA		
!	Write(1602,'(10E20.8)') rv(1,i),rv(2,i),u(1,i),u(2,i),s(1,i),
!     :s(2,i),rsx(j,i),rsy(j,i),atan2(s(2,i),s(1,i)),mi/dble(NA)
!      	  enddo
!	 enddo
 !       endif
	If (step.GT.stepaverage)countt = countt + 1

	
	
	  Call LeapFrogstep (r,rv,ra,s,sm,u,g,sa,kienergy,usum,
     :totenergy,potenergy,mominert,tempinst)
		
	  Call ApplyPBCR(switchDiff,y_box,r,rcross)
	  Call ConvertPosi (D,r,s,rsx,rsy)
	  Call ApplyPBC(y_box,rsx,rsy)

          Call Force (usum,rsx,rsy,FSX,FSY,B,s)
         
	  Call ConverForTor (mominert,D,r,ra,s,sa,FSX,FSY,u,B,rv,2)
	  kienergy = 0.d0
      
        Do i = 1,nmol
            rv(1,i) = facpe*(rv(1,i) + deltaT*ra(1,i)*0.5d0)
           rv(2,i) = facpa*(rv(2,i) + deltaT*ra(2,i)*0.5d0)
           u(1,i) = facrot*(u(1,i) + deltaT*sa(1,i)*0.5d0)
           u(2,i) = facrot*(u(2,i) + deltaT*sa(2,i)*0.5d0)
!###########Turning the velocities into the standart coordinates#####
           VXI =  rv(1,i)*s(2,i) + rv(2,i)*s(1,i) 
           VYI =  rv(2,i)*s(2,i) - rv(1,i)*s(1,i)
!####################################################################        
           rv(1,i) = VXI
           rv(2,i) = VYI

           kienergy = kienergy + (rv(1,i)*rv(1,i) + rv(2,i)*rv(2,i))
     :+ mominert*(u(1,i)*u(1,i) + u(2,i)*u(2,i))
           enddo
        kienergy = 0.5d0*kienergy/dble(nmol)
	potenergy = usum/dble(nmol*NA)
	totenergy = kienergy + potenergy
	tempinst = 2.d0*kienergy/3.d0
        
	 If(switchDiff.eq.2) then
	call  Diffusion(step,cDiff1,nTimeWin1,TimeWindow1,xp,yp,
     :cDiff2,nTimeWin2,TimeWindow2,r,rcross,Diffx,Diffy,s,thetap,
     :Difftheta,Diffpa,Diffpe)
	Endif
!	 Call AdjustDipole (s)  
	
       
	
	Call CLustering_analysis(step,countt,countrat,LL,sckvid,
     :ratkvid,nclukvi,Polimkvi,Polim2,histang,histang2,r,rsx,rsy,s)
        
!	
          	
!	  

	  
      if(mod(step,2000).eq.0.d0)Write(*,'(I10,4E20.8)')step,
     :totenergy,potenergy,kienergy,tempinst
      if(mod(step,2000).eq.0.d0)Write(4,'(I10,4E20.8)')step,totenergy
     :,potenergy,kienergy,tempinst

!	If (step.le.50000) then
	IF (P.eq.1) then
!###################################################################
	if (mod(step,configtest).eq.0) then

!!	kvid=kvid+1
!!        nome_arq_1 = 'config_'						! arquivo de video
!!       write(nome_arq_2,'(I5.5)')kvid
!!      nome_arq_total = nome_arq_1 // nome_arq_2 
!!       open(unit=18,file =  nome_arq_total //'.dat')
!!	Do i=1,nmol
!!	  Do j=1,NA		
!!	Write(18,'(3E20.8)') rsx(j,i),rsy(j,i),atan2(s(2,i),s(1,i))
!!    	  enddo
!!      enddo
	
	
		write(44,*) nmol*NA
		write(44,*) 'Jorge-Lennard-Jones'
		
		Do n = 1,nmol
		  Do j=1,NA	
			write(44,*)"Default",rsx(j,n),rsy(j,n),0 
		enddo
	       enddo
		

	endif
!###################################################################
	 endif
!	endif
       If (step.GT.stepaverage) then
    
!############### Diffusion  #####################################	
        If (mod(step,nsamp).eq.0) then	
        Call dif(switch,nsamp,ntime,ntel,t0,tt0,time0,
     :s,sx0,sy0,dtime,corrmi,Magtotal,magcor,rsx,rsy,nbondaver,
     :nbondcor)
       
	switch = 1
        If(mod(step,1000).eq.0) switch = 2 
	endif
	If(switchDiff.eq.1) then

       Call DifInit(cDiff1,nTimeWin1,TimeWindow1,Diffx,Diffy,xp
     :,cDiff2,nTimeWin2,TimeWindow2,yp,r,s,thetap,Difftheta,
     :Diffpa,Diffpe)
     
	switchDiff = 2
	endif
!#################################################################
    

!       Call Polarization(countpol,countpol2,Pol,NemOrd,
!     :Polsqr,NemOrdsqr,Pol2,NemOrd2,Pol2sqr,NemOrd2sqr,s)
	

	
	
	

      If(mod(step,1000).eq.0)Call EvalRDF(LL,countRdf,countRdf2,r,
     :y_box,histRdf,histRdf2,histCom,histCom2,s)
	
	toten = toten + totenergy
	poten = poten + potenergy
	kinen = kinen + kienergy       
	
	Endif	
	If (nstep.le.10) then
	If(mod(step-1,10000).eq.0.d0) then
         nstep = nstep + 1	
           nome_arq_3 = 'parcialcond_nstep_'						
           write(nome_arq_2,'(I2)')nstep
          nome_arq_total2 = nome_arq_3 // nome_arq_2
         open(988,file =  nome_arq_total2 //'.dat')
        Do i=1,nmol
	  Do j=1,NA
         Write(988,'(2E20.8)') rsx(j,i),rsy(j,i)
        enddo
        enddo
        endif
        close(988)
        endif

	If(mod(step,10000).eq.0.d0) then	
	open(98,file='finalconditions.dat')
	Do i=1,nmol
	  Do j=1,NA		
	Write(98,'(10E20.8)') rv(1,i),rv(2,i),u(1,i),u(2,i),s(1,i),
     :s(2,i),rsx(j,i),rsy(j,i),atan2(s(2,i),s(1,i)),mi!/dble(NA)
      	  enddo
	 enddo
	totenn = toten/dble(countt)
	potenn = poten/dble(countt)
	kinenn = kinen/dble(countt)
	write(8,'(3E20.10)') totenn,potenn,kinenn	
	close(98)
	close(8)
	endif
	

       
	Enddo ! Steps
	
	
    
		
!	Do i=1,nmol
!	  Do j=1,NA		
!	Write(69,'(10E20.8)') rv(1,i),rv(2,i),u(1,i),u(2,i),s(1,i),
!     :s(2,i),rsx(j,i),rsy(j,i),atan2(s(2,i),s(1,i)),mi!/dble(NA)
!      	  enddo
!	 enddo
     
      end Program
      
      
 !####################################################################     
       Subroutine PosiInit(y_box,u,r,s,D,rsx,rsy)
       Implicit None
       include 'constants.f'
      Integer*4 seed(3)
      Integer n1,i,j,initcell(ndim),k,aspc(nmol),D(NA),m,n,B1,B2
      Real*8 gap(ndim),lcell(ndim),r(ndim,nmol),s(ndim,nmol),pi,
     :y_box,raux(ndim,nmol+1000),saux(ndim,nmol+1000),c(ndim),TunThe      
      Real*8 sm(ndim,nmol),u(ndim,nmol),rand48,mit(ndim,nmol)
      Real*8 rsx(NA,nmol),rsy(NA,nmol),dx,dy,dr,aux(nmol),theta
   
       	COMMON /aleat/seed
      
          pi = 4.d0*atan(1.d0)
             	n1 = -1     
500       Do m=1,nmol
	      if(aux(m).eq.0) then
 	      r(1,m) = (x_box*0.5d0)*(2.d0*rand48(seed)-1.d0)
 	      r(2,m) = (y_box*0.5d0)*(2.d0*rand48(seed)-1.d0)
               theta = 2.d0*pi*rand48(seed)
	      s(1,m) = dcos(theta)
	      s(2,m) = dsin(theta) 
	      aux(m) = 1
	      endif
          enddo
	n1 = n1 + 1
	write(*,*)'Number of overlaps',n1,dr,i,j

	Call ConvertPosi (D,r,s,rsx,rsy)
  

!      stop
          do i=1,nmol-1
	  do j=i+1,nmol
		Do B1 = 1,NA
		 Do B2 = 1,NA
		dx = rsx(B1,i)-rsx(B2,j)
	        dy = rsy(B1,i)-rsy(B2,j)
                 
	    IF (dabs(dx).gt.(0.5d0*x_box)) dx = dx - dsign(x_box,dx)
	    IF (dabs(dy).gt.(0.5d0*y_box)) dy = dy - dsign(y_box,dy) 
		dr = dsqrt(dx*dx + dy*dy)		
				
		if(dr.lt.sigma) then
		aux(j) = 0
		go to 500
		endif	
		enddo
	      enddo
	     enddo
	   enddo
		
	
 	   write(*,*)'Initial conditions have been settled!'


       return
       
       End
           
           
 !###################################################################
 	Subroutine InitLinearVels(vmagli,rv)
	implicit none
	include 'constants.f'
	integer*4 seed(3)
	integer k,i,l
	Real*8 rand48,pi,vsum(ndim),rv(ndim,nmol),vmagli,ang

!	Common /link1/rv
	COMMON /aleat/seed
	
	pi = 4.d0*atan(1.d0)
	
! 		write(*,*) ndim,nmol
	
	Do i=1,ndim
	vsum(i)= 0.d0
	enddo

	Do i=1,nmol
	 ang = 2.d0*pi*rand48(seed)
	 rv(1,i)=vmagli*dcos(ang)
	 rv(2,i)=vmagli*dsin(ang)
	 
	 Do k=1,ndim
	  vsum(k)= vsum(k)+rv(k,i)
	 enddo
	enddo

	Do k=1,ndim
	 vsum(k)= vsum(k)/nmol
	enddo

	Do i=1,nmol
	 Do k=1,ndim
	 rv(k,i)= rv(k,i) - vsum(k) 		
	enddo
	enddo

       
	return
	
	end	
		
!####################################################################

	Subroutine Initangvel (vmagrot,mominert,u,g) 
	Implicit none
	include 'constants.f'
	Integer k,n
	real*8 ang,pi,vmagrot,mominert
	real*8 u(ndim,nmol),g(ndim,nmol)
	Integer*4 seed(3)
	Real*8 rand48
	COMMON /aleat/seed
!	COMMON /sa/sa1,sa2#
	
	
	pi= 4.d0*atan(1.d0) 

	
	Do n=1,nmol
	
	
	    ang=2.d0*pi*rand48(seed)
	
	    u(1,n)=vmagrot*cos(ang)
	    u(2,n)=vmagrot*sin(ang)
	
	 do k=1,ndim
	    g(k,n)=0.d0
	
	 end do
	end do
	

	end



         
!####################################################################
       Subroutine Force (usum,rsx,rsy,FSX,FSY,B,s)
      implicit none
	include 'constants.f'
      integer i,j,B1,B2,step,Nbonds(nmol),sumbond(6)
      integer countt
      Real*8 FSX(NA,nmol),FSY(NA,nmol),modr,sri,srj,B(ndim,nmol)
      Real*8 f1x,f1y,a1,a2,a3,rr,rri,ss,ffx,ffy,l2,Lx,usum,ecut,ss1
      Real*8 dx,dy,s(ndim,nmol),rsx(NA,nmol),rsy(NA,nmol),y_box,sFac
      Real*8 rcuti,modrcut,a1cut,a2cut,a4,rrrcuti,bxi,byi,bxj,byj
      Real*8 sumaveg(7),sumvegpol(7),sumveg2(7),nbondaver,
     :sumbondreal(7)
	 Common /link2/step
!	   open(222,file ='ddd.dat')
        


	y_box = x_box
	l2 = (0.5d0*80.d0)**2+(0.5d0*80.d0)**2
         bxi = 0.d0
	byi = 0.d0
	bxj = 0.d0
	byj = 0.d0
	a1 = 0.d0
	a2 = 0.d0
	a3 = 0.d0
	ss= 0.d0
	sri = 0.d0
	srj = 0.d0
	ffx = 0.d0
	ffy= 0.d0
	f1x = 0.d0
	f1y = 0.d0
!!	If (step.eq.stepaverage) then
!!	Do i = 0,6
!!	sumvegpol(i) = 0.d0
!!	sumveg2(i) = 0.d0
!!	sumaveg(i) = 0.d0
!!	enddo
!!	endif
       
	call AdjustDipole(s)

!	If(step.eq.536) then
!	Do i = 1,nmol
 !       write(*,'(2E30.22,2I4)')s(1,i),s(2,i),i
!	enddo
!	pause	
!	endif

	Do i=1,nmol 

	      B(1,i) = 0.d0
	      B(2,i) = 0.d0

		DO j=1,NA
	      FSX(j,i) = 0.d0
	      FSY(j,i) = 0.d0
	       enddo
	enddo

	usum = 0.d0
	
        Do i =1,nmol-1
		
	  Do j= i+1,nmol
!	  if(j.ne.i) then
!	 write(222,*)i,j
		Do B1 = 1,NA
	       	 Do B2 = 1,NA
	       dx = rsx(B1,i)-rsx(B2,j)
	       dy = rsy(B1,i)-rsy(B2,j)
                 
	    IF (dabs(dx).gt.(0.5d0*x_box)) dx = dx - dsign(x_box,dx)
	    IF (dabs(dy).gt.(0.5d0*y_box)) dy = dy - dsign(y_box,dy)         
     

	      
	     modr = dsqrt(dx*dx+dy*dy)
	     rr = dx*dx+dy*dy
	     rri = 1.d0/dble(rr)     
		
         if(rr.lt.l2) then
		 
!	  sFac = 1.d0/(dsqrt(s(1,i)*s(1,i)+s(2,i)*s(2,i)))
!		s(1,i) = sFac*s(1,i)
!		s(2,i) = sFac*s(2,i)

!	  sFac = 1.d0/(dsqrt(s(1,j)*s(1,j)+s(2,j)*s(2,j)))	
!		s(1,j)=sFac*s(1,j)
!		s(2,j)=sFac*s(2,j)

!	 Call AdjustDipole (s)  
!####################################################################

	 a1 = (mi)*(mi)*rri/dble(modr)
	 a2 = 3.d0*a1*rri
	 a3 = 5.d0*a2*rri
	 ss = s(1,i)*s(1,j) + s(2,i)*s(2,j)
	 sri = s(1,i)*dx + s(2,i)*dy
	 srj = s(1,j)*dx + s(2,j)*dy
	 f1x = (48.d0)*dx*(rri)**7
	 f1y = (48.d0)*dy*(rri)**7
	
      ffx =f1x +(a2*ss - a3*sri*srj)*dx + a2*(s(1,i)*srj+s(1,j)*sri)
      ffy =f1y +(a2*ss - a3*sri*srj)*dy + a2*(s(2,i)*srj+s(2,j)*sri)
!	write(*,'(2e20.12,2I8)') ffx,ffy,i,j
	
	
	 
	bxi = -a1*s(1,j) + a2*srj*dx
	byi = -a1*s(2,j) + a2*srj*dy
	bxj = -a1*s(1,i) + a2*sri*dx
	byj = -a1*s(2,i) + a2*sri*dy
	
	ss1 = s(1,i)*s(1,i)+s(2,i)*s(2,i)

	
	
	
!	If(i.eq.10.or.j.eq.10) write(*,'(4E30.22)')sri,srj,i,j
!	if(i.eq.10.or.j.eq.10 ) write(*,'(2E30.22,2I8)')bxi,bxj,i,j	
	B(1,i) = B(1,i)+bxi
	B(2,i) = B(2,i)+byi
	B(1,j) = B(1,j)+bxj
	B(2,j) = B(2,j)+byj

	FSX(B1,i) = FSX(B1,i)+ffx
	FSY(B1,i) = FSY(B1,i)+ffy

	FSX(B2,j) = FSX(B2,j)-ffx
	FSY(B2,j) = FSY(B2,j)-ffy
	
	a1cut = (mi)*(mi)/(dsqrt(l2))**3
	a2cut = 3.d0*a1cut/l2
	
	ecut = (4.d0)*(1.d0/l2)**6 + a1cut*ss - a2cut*sri*srj 
!	usum = usum + 4.d0*a1*(a1-1.d0) - ecut
    
	usum =usum + a1*ss-a2*sri*srj+(4.d0)*(rri)**6 - ecut
	
	
!	write(*,'(A2,e30.22)')'B',B(1,4)	
       endif
	       enddo
	     enddo	
!		endif	
	  enddo

	enddo

	
       




        Do j = 1,nmol
        usum =usum-mi*NA*B0*(s(1,j)*dcos(omega*deltaT*(step))+
     :s(2,j)*dsin(omega*deltaT*(step)))
        !obs (mi/NA) x NA x B0 = mi x B0 
        
        enddo
     
      
!	pause
	return
	end
!####################################################################
	!######################################################################     
      Subroutine MagneticTorq(mominert,r,B,s,u,rsx,rsy)
      implicit none
      include 'constants.f'
      integer i,j,B1,B2
      Real*8 modr,B(ndim,nmol)
      Real*8 f1x,f1y,a1,a2,a3,rr,rri,l2,usum,ecut
      Real*8 dx,dy,r(ndim,nmol),y_box
      Real*8 rcuti,modrcut,ffx,ffy,ss,sri,srj,u(ndim,nmol)
      Real*8 s(ndim,nmol),mominert,sa(ndim,nmol),bxi,byi
      Real*8 Bperx,Bpery,SB,SU,rsx(NA,nmol),rsy(NA,nmol)
      

	   
	y_box = x_box
	   l2 = (0.5d0*x_box)**2+(0.5d0*y_box)**2
         
	Do i=1,nmol
	
	      B(1,i) = 0.d0
	      B(2,i) = 0.d0

	enddo

	
	
        Do i =1,nmol
		
	  Do j=1,nmol
		IF(j.ne.i) then
	       Do B1 = 1,NA
	        Do B2 = 1,NA
	       dx = rsx(B1,i)-rsx(B2,j)
	       dy = rsy(B1,i)-rsy(B2,j)
                 
	    IF (dabs(dx).gt.(0.5d0*x_box)) dx = dx - dsign(x_box,dx)
	    IF (dabs(dy).gt.(0.5d0*y_box)) dy = dy - dsign(y_box,dy)         
     

	      
	     modr = dsqrt(dx*dx+dy*dy)
	     rr = dx*dx+dy*dy
	     
		            
           
		
         if(rr.lt.l2) then 
	  rri = 1.d0/rr     

!########################################################################
	 a1 = mi*mi*rri/modr
	 a2 = 3.d0*a1*rri
	 
	 srj = s(1,j)*dx + s(2,j)*dy
	
	 
	bxi = -a1*s(1,j) + a2*srj*dx
	byi = -a1*s(2,j) + a2*srj*dy

	B(1,i) = B(1,i)+bxi
	B(2,i) = B(2,i)+byi
	 


       endif
	enddo
	  enddo
	    endif
	  enddo
!	SB = s(1,i)*B(1,i) + s(2,i)*B(2,i)
!	SU = 2.d0*(s(1,i)*u(1,i)+s(2,i)*u(2,i))/deltaT
!	Bperx = B(1,i) - SB*s(1,i)
!	Bpery = B(2,i) - SB*s(2,i)
!	sa(1,i) = Bperx/mominert - SU*s(1,i)
!	sa(2,i) = Bpery/mominert - SU*s(2,i)   
	enddo

	
	
	return
	end
!###################################################################
	Subroutine LeapFrogstep (r,rv,ra,s,sm,u,g,sa,kienergy,usum,
     :totenergy,potenergy,mominert,tempinst)	
	implicit none
	include 'constants.f'
	integer i,j,switchve,step	
	Real*8 r(ndim,nmol),rv(ndim,nmol),ra(ndim,nmol)
        Real*8 UXI,UYI,SXI,SYI,DOT,UX(nmol),UY(nmol),s(ndim,nmol)
	Real*8 mominert,kienergy,usum,potenergy,totenergy,tempinst
	Real*8 VXI,VYI,MomLi,u(ndim,nmol),g(ndim,nmol),VPA,VPE
	Real*8 sm(ndim,nmol),sx,sy,sa(ndim,nmol),fac,facr
        
	  Common /link1/MomLi
!	  Common /link2/ra	
	 Common /link2/step	
	kienergy = 0.d0
	MomLi = 0.d0
      
	 open(24,file ='vel.dat')	
	 Do i=1,nmol
		 
                VXI = rv(1,i)
               VYI = rv(2,i)
!########### Parallel and Perpendicular components to s #############
            VPE =  VXI*s(2,i) - VYI*s(1,i)
            VPA =  VXI*s(1,i) + VYI*s(2,i)
!####################################################################
           rv(1,i) = VPE + deltaT*ra(1,i)*0.5d0
           rv(2,i) = VPA + deltaT*ra(2,i)*0.5d0

!########## Turning the velocities into the standart coordinates#####
           VXI = rv(2,i)*s(1,i) + rv(1,i)*s(2,i)
           VYI = -rv(1,i)*s(1,i) + rv(2,i)*s(2,i)
!#################################################################### 
    

           r(1,i) = r(1,i) + deltaT*VXI
           r(2,i) = r(2,i) + deltaT*VYI
		
		UXI = u(1,i)
		UYI = u(2,i)
	   u(1,i) = UXI + deltaT*sa(1,i)*0.5d0
           u(2,i) = UYI + deltaT*sa(2,i)*0.5d0

	
	 s(1,i) = s(1,i) + deltaT*u(1,i)
	 s(2,i) = s(2,i) + deltaT*u(2,i)

	 enddo
       Call AdjustDipole (s)
	return
	end
!####################################################################
	Subroutine ApplyPBC(y_box,rsx,rsy)
        implicit none
	include 'constants.f'
	integer A,I,j	
	Real*8 region(ndim),y_box,rsx(NA,nmol),rsy(NA,nmol)
	
         region(1) = 0.5d0*x_box
	 region(2) = 0.5d0*y_box
           
	Do I = 1,nmol
	 	
          Do A=1,NA
	   
	
 
      IF(rsx(A,I).gt.region(1)) rsx(A,I)=rsx(A,I)-2.d0*region(1)
      IF(rsx(A,I).lt.(-region(1)))rsx(A,I)=rsx(A,I)+2.d0*region(1)
	
      IF(rsy(A,I).gt.region(2)) rsy(A,I)=rsy(A,I)-2.d0*region(2)
      IF(rsy(A,I).lt.(-region(2)))rsy(A,I)=rsy(A,I)+2.d0*region(2)

	
	   	
         enddo
	enddo


	   
	return
	end		
! ####################################################################
	Subroutine Evalprops(rv,ra,u,sa,totenergy,kienergy,potenergy,
     :usum,mominert,tempinst)
	implicit none
	include 'constants.f'
	integer i,m	
	real*8 vv,v,totenergy,kienergy,potenergy,usum,tempinst
	real*8 rv(ndim,nmol),u(ndim,nmol),MomLi,mominert,UXI,UYI
	real*8 ra(ndim,nmol),sa(ndim,nmol),VXI,VYI

        
	kienergy = 0.d0
	MomLi = 0.d0
	Do i = 1,nmol
		VXI = rv(1,i)
		VYI = rv(2,i)
           VXI = (VXI+rv(1,i))*0.5d0
	   VYI = (VYI+rv(2,i))*0.5d0
		
		UXI = u(1,i)
		UYI = u(2,i)		

		UXI = (UXI+u(1,i))*0.5d0
		UYI = (UYI+u(2,i))*0.5d0	 
 	  kienergy = kienergy + mominert*(UXI*UXI + UYI*UYI) +
     :(VXI*VXI + VYI*VYI)
	 MomLi = MomLi + UXI + UYI
	
	 
	enddo
	kienergy = 0.5*kienergy/nmol
	potenergy = usum/(nmol*NA)
	totenergy = kienergy + potenergy 
	tempinst = 2.d0*kienergy/3.d0	
	return
	end
!####################################################################
	Subroutine ApplyPBCR(switchDiff,y_box,r,rcross)
        implicit none
	include 'constants.f'
	integer A,I,J,switchDiff
	Real*8 region(ndim),y_box
	REal*8 r(ndim,nmol),rcross(ndim,nmol)
         region(1) = 0.5d0*x_box
	 region(2) = 0.5d0*y_box
           
	Do I = 1,nmol
	 Do J= 1,ndim

	  IF(r(J,I).gt.region(J)) then
	   r(J,I)=r(J,I)-2.d0*region(J)
	   If(switchDiff.eq.2)rcross(J,I)=rcross(J,I)+1.d0
	  Endif

	  IF(r(J,I).lt.(-region(J))) then
	   r(J,I)=r(J,I)+2.d0*region(J)
           If(switchDiff.eq.2)rcross(J,I)=rcross(J,I)-1.d0
	  Endif
	
	 enddo		
	enddo
	return
        end
!####################################################################
        SUBROUTINE ConvertPosi (D,r,s,rsx,rsy)
	implicit none
	include 'constants.f'
C *******************************************************************
C ** CONVERTS C-O-M COORDINATES AND BOND VECTOR TO SITE POSITIONS. **
C **                                                               **
C ** THE POSITION OF EACH ATOM IN THE MOLECULE IS DEFINED IN TERMS **
C ** OF THE UNIT BOND VECTOR EX(I),EY(I),EZ(I) AND THE ATOM        **
C ** POSITION VARIABLE D(A): RSX(I,A) = RX(I) + D(A)*EX(I) ETC.    **
C **                                                               **
C ** PRINCIPAL VARIABLES:                                          **
C **                                                               **
C ** INTEGER nmol                         NUMBER OF MOLECULES     **
C ** INTEGER NA                            NUMBER OF ATOMS PER MOL **
C ** REAL    r(2,nmol)                    POSITIONS AT TIME T     **
C ** REAL    s(2,nmol)                    UNIT BOND VEC AT TIME T **
C ** REAL    D(NA)                         ATOM POSITIONS IN MOLEC **
C ** REAL    rs(2,nmol,NA)                ATOM POSITIONS          **
C *******************************************************************

        INTEGER I, A
!	integer step    
        REAL*8 SXI,SYI,s(ndim,nmol),rsx(NA,nmol),rsy(NA,nmol),D(NA)  
        REAL*8 r(ndim,nmol),y_box
	y_box = x_box
  
!	 Common /link2/step
        DO I = 1, nmol
	
!	if (step.ne.0.d0) Write(3,*) r(1,I),r(2,I)
           SXI = s(1,I)
           SYI = s(2,I)


		
           DO A = 1, NA
	         D(A) = ((1.d0-NA)/2.d0+(A-1.d0))
            rsx(A,I) = r(1,I) + D(A)*SXI
            rsy(A,I) = r(2,I) + D(A)*SYI
	        
!!            IF(rsx(A,I).gt.x_box/2.d0) rsx(A,I)=rsx(A,I)-x_box
!!	    IF(rsx(A,I).lt.(-x_box/2.d0)) rsx(A,I)=rsx(A,I)+x_box
!!	    IF(rsy(A,I).gt.(y_box/2.d0)) rsy(A,I)=rsy(A,I)-y_box	
!!	    IF(rsy(A,I).lt.(-y_box/2.d0)) rsy(A,I)=rsy(A,I)+y_box	

	  enddo

        enddo
	 

        RETURN
        END
!####################################################################

   
      SUBROUTINE ConverForTor (mominert,D,r,ra,s,sa,FSX,FSY,u,B,rv,
     :switchve)

C *******************************************************************
C ** CONVERT ATOM FORCES TO TOTAL FORCES AND AUXILIARY TORQUES.    **
C **                                                               **
C ** PRINCIPAL VARIABLES:                                          **
C **                                                               **
C ** INTEGER N                             NUMBER OF MOLECULES     **
C ** INTEGER NA                            NUMBER OF ATOMS PER MOL **
C ** REAL    RX(N),RY(N),RZ(N)             POSITIONS AT TIME T     **
C ** REAL    FX(N),FY(N),FZ(N)             C-O-M FORCES            **
C ** REAL    EX(N),EY(N),EZ(N)             UNIT BOND VEC AT TIME T **
C ** REAL    GX(N),GY(N),GZ(N)             AUXILIARY TORQUE AT T   **
C ** REAL    D(NA)                         ATOM POSITIONS IN MOLEC **
C ** REAL    RSX(N,NA),RSY(N,NA),RSZ(N,NA) ATOM POSITIONS          **
C ** REAL    FSX(N,NA),FSY(N,NA),FSZ(N,NA) FORCES ON EACH ATOM     **
C *******************************************************************
	implicit none
	include 'constants.f'
       INTEGER   I, A,step,switchve
       integer*2 seed(3)
       Real*8 FXI,FYI,GXI,GYI,SXI,SYI,FSXIA,FSYIA,sa(ndim,nmol)
       Real*8 D(NA),FSX(NA,nmol),FSY(NA,nmol),gx,gy,ra(ndim,nmol)
       Real*8 DOT,RXI,RYI,r(ndim,nmol),s(ndim,nmol),mominert,OT
       Real*8 u(ndim,nmol),UXI,UYI,B(ndim,nmol),TXI,TYI,dti,pi
       Real*8 ang2,rand48,rv(ndim,nmol),u1,u2,u3,mag,magr,ang
       Real*8 corpa,corpe,corrot,gamape,gamapa,gamarot,magpe,magpa
     :,magrot,FPA,FPE,VPA,VPE,DOTg,DOTb
!	Common /link2/ra
	Common /link2/step
        COMMON /aleat/seed
        pi = 4.d0*atan(1.d0)
	
!           Do i=1,nmol
!		Do l=1,NA	   
!		write(*,*) s(1,i),s(2,i)
!		enddo
!            enddo


C  *******************************************************************
	 corpe = 0.839d0 + 0.185d0/dble(NA) + 0.233d0/dble(NA)**2 
        corpa = - 0.207d0 + 0.987d0/dble(NA) - 0.133d0/dble(NA)**2
       corrot = - 0.622d0 + 0.917d0/dble(NA) - 0.050d0/dble(NA)**2 

        gamape = mult*4.d0*NA*pi/(log(dble(NA))+corpe)
        gamapa = mult*2.d0*NA*pi/(log(dble(NA))+corpa)
       gamarot = mult*NA*NA*NA*pi/(3.d0*log(dble(NA))+corrot)
   
        DO I = 1, nmol
        u1 = rand48(seed)
	u2 = rand48(seed)
	u3 = rand48(seed)
        if (log(u1).lt.-10.d5) u1 = 10.d-5
	magpe = dsqrt((-4.0d0*gamape*temp*log(u1)))/dsqrt(deltaT)
        magpa = dsqrt((-4.0d0*gamapa*temp*log(u1)))/dsqrt(deltaT)
       magrot = dsqrt((-4.0d0*gamarot*temp*log(u1)))/dsqrt(deltaT)
        
!	mag = min(mag,1.d0)
	ang = 2.0d0*pi*u2
	ang2 = 2.0d0*pi*u3	
	 
           FXI = 0.d0
           FYI = 0.d0
           GXI = 0.d0
           GYI = 0.d0
	   

	   UXI = u(1,I)
	   UYI = u(2,I)
           SXI = s(1,I)
           SYI = s(2,I)
		
		
           DO A = 1, NA
	   D(A) = ((1.d0-NA)/2.d0+(A-1.d0))
            FSXIA = FSX(A,I)
            FSYIA = FSY(A,I)
              FXI = FXI + FSXIA
              FYI = FYI + FSYIA
              
              GXI = GXI + D(A)*FSXIA
              GYI = GYI + D(A)*FSYIA
             

           enddo
        IF(switchve.eq.1) then
	DOTg = GXI*SXI+GYI*SYI
        DOTb = (B(1,I)+mi*NA*B0*dcos(omega*deltaT*(step)))*SXI + 
     :(B(2,I)+mi*NA*B0*dsin(omega*deltaT*(step)))*SYI
	
        
          gx =  GXI + (B(1,I)+mi*NA*B0*dcos(omega*deltaT*(step)))
     :- DOTg*SXI - DOTb*SXI
          gy =  GYI + (B(2,I)+mi*NA*B0*dsin(omega*deltaT*(step)))
     :- DOTg*SYI - DOTb*SYI  
          
  !####################################         IF(switchve.eq.1) then
           OT = (UXI*UXI + UYI*UYI)
	   dti = deltaT
      sa(1,I) = gx/dble(mominert) -OT*SXI -gamarot*UXI/dble(mominert) 
     : - magrot*(dcos(ang2)*SYI)/dble(mominert)

      sa(2,I) = gy/dble(mominert) -OT*SYI -gamarot*UYI/dble(mominert)
     : + magrot*(dcos(ang2)*SXI)/dble(mominert)

!########### Parallel and Perpendicular components to s #############
        FPE = FXI*SYI - FYI*SXI
        FPA = FXI*SXI + FYI*SYI 
        VPE = rv(1,I)*SYI - rv(2,I)*SXI
        VPA = rv(1,I)*SXI + rv(2,I)*SYI 
   
!####################################################################

      ra(1,I) = FPE/dble(NA) - gamape*VPE/dble(NA) +
     :(magpe*(dcos(ang)*dsin(ang2)*SYI+ 
     :dsin(ang)*dsin(ang2)*SXI))/dble(NA)
      ra(2,I) = FPA/dble(NA) - gamapa*VPA/dble(NA) + 
     :(magpa*(-dcos(ang)*dsin(ang2)*SXI + 
     :dsin(ang)*dsin(ang2)*SYI))/dble(NA)	
	  
            Else

        
        DOTg = GXI*SXI+GYI*SYI
        DOTb = (B(1,I)+mi*NA*B0*dcos(omega*deltaT*(step)))*SXI + 
     :(B(2,I)+mi*NA*B0*dsin(omega*deltaT*(step)))*SYI        
          gx =  GXI + (B(1,I)+mi*NA*B0*dcos(omega*deltaT*(step)))
     :- DOTg*SXI - DOTb*SXI
          gy =  GYI + (B(2,I)+mi*NA*B0*dsin(omega*deltaT*(step)))
     :- DOTg*SYI - DOTb*SYI  

             dti = 1.d0/dble(deltaT)
            OT =  2.d0*(UXI*SXI + UYI*SYI)*dti
         sa(1,I) = gx/dble(mominert) - OT*SXI
     :- magrot*(dcos(ang2)*SYI)/dble(mominert)
	 sa(2,I) = gy/dble(mominert) - OT*SYI 
     :+ magrot*(dcos(ang2)*SXI)/dble(mominert)
!########### Parallel and Perpendicular components to s #############
      
        FPE = FXI*SYI - FYI*SXI
        FPA = FXI*SXI + FYI*SYI 
!####################################################################
       ra(1,I) = FPE/dble(NA) +
     :(magpe*(dcos(ang)*dsin(ang2)*SYI+ 
     :dsin(ang)*dsin(ang2)*SXI))/dble(NA)

      ra(2,I) = FPA/dble(NA) + 
     :(magpa*(-dcos(ang)*dsin(ang2)*SXI + 
     :dsin(ang)*dsin(ang2)*SYI))/dble(NA)
        Endif
	 enddo
       
        RETURN
        END

!####################################################################
	Subroutine AdjustDipole(s)
	 implicit none
	include 'constants.f'
	integer i,j
	Real*8 s(ndim,nmol),sFac	
	
	Do i=1,nmol
		sFac = 1.d0/(dsqrt(s(1,i)*s(1,i)+s(2,i)*s(2,i)))	
	 Do j =1,ndim
		s(j,i)=sFac*s(j,i)
	 enddo
	enddo

	return
	end	
!####################################################################

      Subroutine AdjustTempLine(vmagli,rv)
	implicit none
	include 'constants.f'
	Integer	n,k
	Real*8	vvsum,rv(ndim,nmol),vmagli,vFac,tal
	tal = 1.1d0*deltaT
         vvsum = 0.d0
	 Do n=1,nmol
	  Do k=1,ndim
	     vvsum = vvsum + rv(k,n)*rv(k,n)
	      enddo
	   
	  enddo
	   vFac = vmagli*vmagli/(vvsum/dble(nmol))
	    vFAc = dsqrt(1 + (deltaT/tal)*(vFac-1))	
	  Do n = 1,nmol
	   Do k=1,ndim
		rv(k,n) = rv(k,n)*vFac
	   enddo
	  enddo		
	return
	end
!####################################################################
	 Subroutine AdjustTempang(mominert,vmagli,u)
	implicit none
	include 'constants.f'
	Integer	n,k
	Real*8	vvsum,u(ndim,nmol),vmagli,vFac,mominert,tal
	tal = 1.1d0*deltaT
         vvsum = 0.d0
	 Do n=1,nmol
	  Do k=1,ndim
	     vvsum = vvsum + mominert*u(k,n)*u(k,n)
	      enddo
	   
	  enddo
	  vFac = vmagli*vmagli/(2.d0*vvsum/dble(nmol))
	 vFAc = dsqrt(1 + (deltaT/tal)*(vFac-1))	
	  Do n = 1,nmol
	   Do k=1,ndim
		u(k,n) = u(k,n)*vFac
	   enddo
	  enddo		
	return
	end
!#####################################################################
	Subroutine ApplyThermostat(mominert,rv,ra,sa,u)
	   implicit none
	include 'constants.f'
	 integer n,k
	 Real*8 s1,s2,rv(ndim,nmol),ra(ndim,nmol),vFFac,sa(ndim,nmol)
	 Real*8 u(ndim,nmol),mominert
	
	 s1 = 0.d0
	 s2 = 0.d0
	 Do n = 1,nmol
	  Do k = 1,ndim
            s1 = s1 + rv(k,n)*ra(k,n)
	    s2 = s2 + rv(k,n)*rv(k,n) 
	   enddo
	 enddo

	   Do n =1,nmol
             Do k=1,ndim	 
	   s1 = s1 + mominert*u(k,n)*sa(k,n)
	   s2 = s2 + mominert*u(k,n)*u(k,n)
	     Enddo
            Enddo
	
	vFFac = -s1/s2
          
	 Do n=1, nmol
	  Do k=1,ndim
   	    ra(k,n) = ra(k,n)+ vFFac*rv(k,n)
	 enddo
	 enddo
	
	Do n=1, nmol
	  Do k=1,ndim
   	    sa(k,n) = sa(k,n)+ vFFac*u(k,n)
	 enddo
	 enddo
  
	return
	end
!####################################################################	 
c       *************************************************************
c       **************Funcao de numeros aleatorios*******************


        FUNCTION rand48(seed)
c
c                                Fortran 77 function that implements
c                                the random number generator used
c                                by the C utility erand48. It produces
c                                the same values as erand48 for the
c                                three Integer*4s seed(3) and a sequence
c                                of random numbers that agree in all but
c                                the least significant digits. 
c                                (Very small discrepancies in the least 
c                                significant digits are produced by the 
c                                different truncations used by Fortran
c                                and C.)   
c
c                                        Claudio Rebbi - Boston University
c                                                        May 1992 
c
      Integer*4 seed(3)
      Integer*4 i1,i2,i3,i11,i21,i31,i12,i22,i13
      Integer*4 E66D,DEEC,FFFF
      PARAMETER(E66D=58989, DEEC=57068, FFFF=65535)
      Real*8 rand48

      i1=seed(1)
      IF(i1.LT.0) i1=i1+65536
      i2=seed(2)
      IF(i2.LT.0) i2=i2+65536
      i3=seed(3)
      IF(i3.LT.0) i3=i3+65536

        i11=i1*E66D
        i21=i2*E66D
        i31=i3*E66D
        i12=i1*DEEC
        i22=i2*DEEC
        i13=i1*5
        i1=IAND(i11,FFFF)+11
        i11=ISHFT(i11,-16)+ISHFT(i1,-16)
        i1=IAND(i1,FFFF)
        i11=i11+IAND(i21,FFFF)+IAND(i12,FFFF)
        i2=IAND(i11,FFFF)
        i3=ISHFT(i11,-16)+ISHFT(i21,-16)+ISHFT(i12,-16)+
     &     IAND(i31,FFFF)+IAND(i22,FFFF)+IAND(i13,FFFF)
        i3=IAND(i3,FFFF)
c
c       rand48=i3*2**(-16)+i2*2**(-32)
c
        rand48=i3*1.52587890625D-05+i2*2.328306D-10

      IF(i1.GE.32768) i1=i1-65536
      seed(1)=i1
      IF(i2.GE.32768) i2=i2-65536
      seed(2)=i2
      IF(i3.GE.32768) i3=i3-65536
      seed(3)=i3

      RETURN
      END

!####################################################################
!                          PROPERTIES                               #
!#################################################################### 


 

! ******************************************************************

! **************Radial distribution function************************

! ******************************************************************

      Subroutine EvalRDF(LL,countRdf,countRdf2,r,y_box,histRdf,
     :histRdf2,histCom,histCom2,s)

 	Implicit none
	include 'constants.f'	
	integer countRdf,n,i,j,j1,j2,limitRdf,k,sizeHistRdf
     :,countRdf2,LL(nmol)       
	Real*8 region(ndim),y_box,rangeRdf,rrRange,deltaR,rr,
     :normFac,dr(ndim+1),pi,regionH(ndim),r(ndim,nmol),rbin
     	Real*8 histRdf(100000),dens,histRdf2(100000),histCom(100000)
	Real*8 histCom2(100000),s(ndim,nmol),ss,p2cos

	  pi = 4.d0*atan(1.d0)
	  
	 limitRdf = 1000
	 deltaR = 0.05d0
          

	 countRdf = countRdf +1
	
	 
	region(1) = x_box	 
	regionH(1) = 0.5d0*region(1)
	region(2) = x_box	 
	regionH(2) = 0.5d0*region(2)
	rangeRdf = regionH(1)
        sizeHistRdf = int(rangeRdf/deltaR)      
	dens = nmol/(x_box*x_box)

	 If (countRdf.eq.1.) then
		Do n = 1,sizeHistRdf
		  histRdf(n) = 0.d0
	          histCom(n) = 0.d0
		enddo

	 endif
		If(countRdf2.eq.0.d0) then
		
		Do n = 1,sizeHistRdf
		  histRdf2(n) = 0.d0
	          histCom2(n) = 0.d0
		enddo
		endif

	 rrRange = rangeRdf*rangeRdf
	 
	 
	  Do j1=1,nmol-1
	    Do j2=j1+1,nmol
	     Do k=1,ndim

	     	dr(k)= r(k,j1)-r(k,j2)

	             IF (dabs(dr(k)).gt.(regionH(k)))
     :dr(k) = dr(k) - dsign(x_box,dr(k))
	
	      enddo

	     rr = dr(1)*dr(1) + dr(2)*dr(2)
	     	
	     if(rr.LT.rrRange) then
!		ss = s(1,j1)*s(1,j2) + s(2,j1)*s(2,j2)
!	     p2cos = 3.d0*(ss*ss) - 1.d0 
		
	     n = dint(dsqrt(rr)/deltaR) +1
	     histRdf(n) = histRdf(n)+1! p2cos
	!     histRdf(n) = histRdf(n)*p2cos
	!	write(*,*)p2cos,histRdf(n)
!		Write(*,*) LL(j1),LL(j2)
	     If(LL(j1).eq.LL(j2)) then
		
	     histCom(n) = histCom(n)+ 1.d0
	!	histCom(n) = histCom(n)	     
	   
	     endif
		 
			
	     endif	

	     enddo
	    enddo 
	  
	   If (countRdf.eq.limitRdf) then
               countRdf2 = countRdf2 +1
	     normFac = (1.0d0/dble(dens*pi*nmol*countRdf))!*5.d0/2.d0 ! 5/2 term is due to h220 corelation function
	     
	
              open(9,file='radial.dat')	
	      open(12,file='radialtotal.dat')
	      open(999,file='Gcom.dat')	
	      open(122,file='Gcomtotal.dat')

   	       Do i=1,sizeHistRdf

 	        rbin = (i-0.5d0)*deltaR
                histRdf(i)=histRdf(i)*normFac/(rbin*deltaR)
	        histRdf2(i) = histRdf2(i) + histRdf(i)

		
		
		histCom(i) = histCom(i)*(normFac)/(rbin*deltaR)
	        histCom2(i) = histCom2(i) + histCom(i)
				
		
                write(9,*)rbin,histRdf(i)
                write(12,*)rbin,histRdf2(i)/dble(countRdf2)
		write(999,*)rbin,histCom(i)
                write(122,*)rbin,histCom2(i)/dble(countRdf2)
	      enddo

	   
           close(9)
	  close(12)
	  close(999)
	  close(122)
         countRdf = 0.d0	 

	endif
	 
	  
	
	return
	end
!#################################################################### 
 !#################################################################### 
  	 Subroutine Polarization(countpol,countpol2,Pol,NemOrd,
     :Polsqr,NemOrdsqr,Pol2,NemOrd2,Pol2sqr,NemOrd2sqr,s)
	implicit none	
	include 'constants.f'
	 
	 Integer n,countpol,limitpol,countpol2
	 Real*8 Qxx,Qyy,Qxy,Lambda1,Lambda2,Lambda,vx,vy,NemOrd
	 Real*8 Modv,vvx,vvy,Pol,s(ndim,nmol),aux,NemOrd2,Pol2
	 Real*8 theta(nmol),Pol2sqr,NemOrdsqr,spol,nempol
	 Real*8 NemOrd2sqr, Polsqr
	  
	 Qxx = 0.d0
	 Qyy = 0.d0
	 Qxy = 0.d0
	 limitpol = 1000
         countpol = countpol + 1
	
	 If (countpol.eq.1) then
	    NemOrd = 0.d0
            Pol = 0.d0
            Polsqr = 0.d0
	    NemOrdsqr = 0.d0
	 endif
   	 
 	     If (countpol2.eq.0) then
	    
               NemOrd2 = 0.d0
               Pol2 = 0.d0
	    	Pol2sqr = 0.d0
	   	NemOrd2sqr= 0.d0 
	    endif
	
	  
	
        Do n=1,nmol
	 theta(n)= atan2(s(2,n),s(1,n))
	 enddo

	 Do n=1,nmol
	      Qxx = Qxx + 3.d0*Cos(theta(n))*Cos(theta(n))-1.d0
	      Qyy = Qyy + 3.d0*Sin(theta(n))*Sin(theta(n))-1.d0
	      Qxy = Qxy + 3.d0*Cos(theta(n))*Sin(theta(n))
	 enddo
	 
	  Qxx = Qxx/(2*nmol)
	  Qyy = Qyy/(2*nmol)
	  Qxy = Qxy/(2*nmol)
	 
	 Lambda1 = ((Qxx+Qyy)+sqrt((Qxx+Qyy)**2-
     :4.d0*(-Qxy**2+Qxx*Qyy)))/2.d0
	 Lambda2 = ((Qxx+Qyy)-sqrt((Qxx+Qyy)**2-
     :4.d0*(-Qxy**2+Qxx*Qyy)))/2.d0
	 
	 If (Lambda1.gt.Lambda2) Lambda=Lambda1
	 If (Lambda1.lt.Lambda2) Lambda=Lambda2
	 
	 vx = 1.d0
	 vy = vx*(Lambda-Qxx)/Qxy
	 
	 Modv = dsqrt(vx*vx+vy*vy)
	 vvx=vx/Modv
	 vvy=vy/Modv

	 aux=0.d0
	 Do n=1,nmol
!	 	 Pol = Pol + (cos(theta(n))*vvx+sin(theta(n))*vvy)/nmol
  		aux=aux+(dcos(theta(n))*vvx+dsin(theta(n))*vvy)/nmol
                
	 enddo 

	   Pol = Pol + abs(aux)
           Polsqr = Polsqr + abs(aux)*abs(aux)
               
	   NemOrd = NemOrd + Lambda
	   NemOrdsqr = NemOrdsqr + Lambda*Lambda

	   If(countpol.eq.limitPol) then
            countpol2 = countpol2 + 1
                	  
	    Pol2 = Pol2 + Pol/dble(countpol)
	    NemOrd2 = NemOrd2 + NemOrd/dble(countpol)

            Pol2sqr = Pol2sqr + Polsqr/dble(countpol)
	    NemOrd2sqr = NemOrd2sqr + NemOrdsqr/dble(countpol)
  
      spol = dsqrt(Pol2sqr/dble(countpol2)-(Pol2/dble(countpol2))**2)
        spol = spol/dsqrt(dble(countpol*countpol2))
      nempol = dsqrt(NemOrd2sqr/dble(countpol2)-
     :(NemOrd2/dble(countpol2))**2)  
       nempol = nempol/dsqrt(dble(countpol*countpol2))
	    open(10,file='polarization.dat')
	    open(11,file='nematicorder.dat')
	    open(13,file='nematicordertotal.dat')
	    open(14,file='polarizationtotal.dat')

	    write(10,'(I4,E21.14)')NA,Pol2
	    write(11,'(I4,E21.14)')NA,NemOrd2
      write(14,'(I4,2E21.14)')NA,Pol2/dble(countpol2),spol
      write(13,'(I4,2E21.14)')NA,NemOrd2/dble(countpol2),nempol
!	    close(10)
!	    close(11)
	    close(13)
	    close(14)
	    countpol = 0.d0
	   endif
	   
	 return
         
	 
	 end     	
	                                                
!###################################################################
       Subroutine EvalVeldist(countVel,countVel2,rv,histVel,histVel2)
	implicit none
	include 'constants.f'
	integer j,k,countVel,countVel2,sizeHistvel
	Real*8 deltaV,histSum,vv,histVel(100000),histVel2(100000)
	Real*8 rbin,rangeVel,rv(ndim,nmol)
	integer n,limitVel

	rangeVel = 0.1d0
	sizeHistvel = 1000
	limitVel = 100	
	deltaV = rangeVel/sizeHistvel	
	countVel = countVel + 1	
	
	if(countVel.eq.1) then

	Do j=1,sizeHistvel
	histVel(j) = 0.d0
	
	enddo
	endif

	if(countVel2.eq.1) then

	Do j=1,sizeHistvel
	histVel2(j) = 0.d0
	
	enddo
	endif

	
!	deltaV = rangeVel/sizeHistVel
	
	Do n=1,nmol
	 vv = rv(1,n)*rv(1,n) + rv(2,n)*rv(2,n)
	 j = dint(dsqrt(vv)/deltaV) + 1
	!	write(*,*)(dsqrt(vv)/deltaV)
	 !  pause
	If (j.gt.sizeHistvel) j = sizeHistvel
	 histVel(j) = histVel(j) + 1.d0  
	enddo	

	If(countVel.eq.limitVel) then
	  countVel2 = countVel2 + 1	
		histSum = 0.d0
		open(90,file='veldis.dat')	
	        open(120,file='veldistotal.dat')
	
        	Do j=1,sizeHistvel
		histSum = histSum + histVel(j)
		enddo
	        
		Do k =1,sizeHistvel
	        rbin = (k-0.5)*deltaV
		histVel(k) = histVel(k)/dble(histSum) 
	 	histVel2(k) = histVel2(k) + histVel(k)
		write(90,*)rbin,histVel(k)
                write(120,*)rbin,histVel2(k)/dble(countVel2)
		enddo
	close(120)
        close(90)
	countVel = 0.d0
	endif
	
	return
	end
!###################################################################


!####################################################################
! NA = number of particles (or patchy) that forms a molecule 
! RCL = Limit of distance which a particle belongs on a given cluster
! rs(1,J,A) x-component of position of a atom 'A' in molecule 'J'   
! rs(2,J,A) same above to y-component
! 
!####################################################################
      Subroutine CLustering_analysis(step,countt,countrat,LL,sckvid,
     :ratkvid,nclukvi,Polimkvi,Polim2,histang,histang2,r,rsx,rsy,s)
    
	Implicit none
	include 'constants.f'
	Integer countcoran,countcoran2
	Integer n,I,LK,L(nmol),LL(nmol),J,A,B,K,countt,step
	Integer aux(nmol),K1,K2,PP,sizecluster(nmol),ncluster,u
	Real*8 dr,RCLSQ,RCL,RXJ(NA),RYJ(NA),y_box,Rg2y(nmol)
	REAL*8 dx,dy,mm,summ(nmol),rsx(NA,nmol),rsy(NA,nmol)
	Integer aux2(nmol),aux4(nmol),aux3(nmol,nmol)
	Real*8 Rcmx(nmol),r(ndim,nmol),Rcmy(nmol),raux_x(nmol,nmol)
	Real*8 raux_y(nmol,nmol),Rat_Gira(nmol),ddx,ddy,ddr

	Real*8 gg,ff,tt(nmol),Polim2
	Integer kvid,aux6(nmol),aux7(nmol),countrat(nmol)
	
        Real*8 Polimkvi,Polim,sckvid(nmol),ratkvid(nmol),Polimspol
	Real*8 Rat_ave(nmol),nclukvi,Rtest,Rg2x(nmol) 
	Integer  HH,countrat3,Gtstclus
	
	Integer histsum,sizehist,ig,limitang
	Real*8 histang(10000),deltacos,ss,histang2(10000)
	Real*8 s(ndim,nmol),sbin,Gtst,Gtst2,Gstpol
	character (len=7) nome_arq_1
        character (len=5) nome_arq_2
        character (len=12) nome_arq_total
	common /angle/ countcoran,countcoran2

            
		
	  countcoran = countcoran + 1
	   RCL = 1.6d0
	   y_box = x_box
!!	   nmol = 0
!!	   sigma = 1.d0 
!!           NA = 1
!!	   x_box = 70
!!	   y_box = 70
	 
!!	 Polimkvi = 0 
!!	 nclukvi = 0 
!!	Do K = 1,10000
!!	sckvid(K) = 0	
!!	ratkvid(K) = 0
!!	countrat(k) = 0
!!	countrat2(k) = 0
!!	enddo
	
!!!        Do kvid= 1,10000
!!!	n = 0

!!!	if(mod(kvid,500).eq.0)write(*,*)kvid
!!!        nome_arq_1 = 'config_'						! arquivo de video
!!!        write(nome_arq_2,'(I5.5)')kvid

!!!       nome_arq_total = nome_arq_1 // nome_arq_2 
!!!       open(unit=18,file = '/home/capuan/Documents/Nanorodskratos/
!!!     :/eta0_3/Naeq1/Files_Config/' // nome_arq_total //'.dat')
!!	open(18,file='config_01078.dat')
	
!!!15    n=n+1
!!!      Do j=1,NA	
!!!      read(18,'(2e20.8)',END=500) rs(1,n,j),rs(2,n,j)      !Ler arquivo das configurações
!!!      enddo	
!      s(1,n)=cos(tt(n))
!      s(2,n)=sin(tt(n))	
!!!      goto 15
!!!    mi = NA*mm
	
!!500    continue
      
!!      nmol=n-1 !kvid 


!!	close(18)

!!	if(mod(kvid,500).eq.0)write(*,*)'nmol=',nmol,NA
!!	close(5)
!!	close(6)
		deltacos = 0.005d0
		sizehist = int(2.d0/deltacos)+1
				
		
		
		limitang = 1000
		
		If(countcoran.eq.1)then
		Do i = 1,sizehist
		histang(i) = 0.d0	
		enddo	
		histsum = 0		
		endif
		If(countcoran2.eq.0)then
		Do i = 1,sizehist
		histang2(i) = 0.d0	
		enddo	
		endif
	RCLSQ = RCL*RCL
	
	Do 10  I =1,nmol
	L(I) = I
	LL(I)= I
10	Enddo

	Do 50 I = 1, nmol-1
!	 write(*,*)'aaa',I,L(I)	
	IF (I.EQ.L(I)) then
	  
		J=I
	  Do A=1,NA 
		
		RXJ(A) = rsx(A,J)
		RYJ(A) = rsy(A,J)
	  Enddo
		Do u=1,nmol
		 aux2(u) = 0		
		enddo
	  
	Do K = I+1,nmol !20
!	   write(*,*) K	
	   LK = L(K)
	   
	  IF(LK.EQ.K) then
	    Do A=1,NA
	    Do B =1,NA		
	   dx = RXJ(A) - rsx(B,K)
	   dy = RYJ(A) - rsy(B,K)
	  IF (dabs(dx).gt.(0.5d0*x_box)) dx = dx - dsign(x_box,dx)
	  IF (dabs(dy).gt.(0.5d0*y_box)) dy = dy - dsign(y_box,dy) 
	  
!	   dx = dx - ANINT(dx)
!	   dy = dy - ANINT(dy)
	   dr = dx*dx + dy*dy
	  IF ((dr.LE.RCLSQ).AND.(aux2(K).eq.0)) THEN
!!	   Call eneergy(dx,dy,s,mi,Epot,J,K)
!!	  If (Epot.LT.0.d0) then
	  !###########################################
		! Angle Corelation :) 		
!		nrodbond = nrodbond + 1
		
        	ss = s(1,I)*s(1,K) + s(2,I)*s(2,K) 
		
		ig = dint((ss+1)/deltacos) +1 !I had to add 'sizehist' because the minimum value is -sizehist+1			
		histang(ig) = histang(ig) + 1.

	  !###########################################
			
	    L(K) = L(J)
	    LL(K) = LL(J)
	    L(J) = LK
	    aux2(K)=1
!!	 endif
	  ENDIF
		ENDDO   !B
		ENDDO	!A
	 ENDIF

20	ENDDO 

	     J=L(J)
	 Do A=1,NA 
		
		RXJ(A) = rsx(A,J)
		RYJ(A) = rsy(A,J)
	  Enddo

		Do u=1,nmol
		 aux2(u) = 0		
		enddo
	 

30	IF(J.NE.I) THEN
	Do K = I+1,nmol   !40
		LK = L(K)
!		write(*,*) I,LK,K	
	  IF(LK.EQ.K) then
	    Do A=1,NA
	    Do B =1,NA		
	   dx = RXJ(A) - rsx(B,K)
	   dy = RYJ(A) - rsy(B,K)

	  IF (dabs(dx).gt.(0.5d0*x_box)) dx = dx - dsign(x_box,dx)
	  IF (dabs(dy).gt.(0.5d0*y_box)) dy = dy - dsign(y_box,dy)
	   dr = dx*dx + dy*dy
!	write(*,*)'aqui1'
	  IF ((dr.LE.RCLSQ).AND.(aux2(K).eq.0)) THEN
!!		 Call eneergy(dx,dy,s,mi,Epot,J,K)
!!	  If (Epot.LT.0.d0) then
!		write(*,*)'aqui2'
	!	nrodbond = nrodbond + 1
		
        	ss = s(1,I)*s(1,K) + s(2,I)*s(2,K) 
		ig = dint((ss+1)/deltacos) +1 !I had to add 'sizehist' because the minimum value is -sizehist+1		
		histang(ig) = histang(ig) + 1.	   
		 L(K) = L(J)
	       LL(K) = LL(J)
	        L(J) = LK
	       aux2(K)=1
!!	  endif
	  ENDIF
		ENDDO   !B
		ENDDO	!A
	 ENDIF

40	ENDDO 

			J=L(J)
	Do A=1,NA 
		
		RXJ(A) = rsx(A,J)
		RYJ(A) = rsy(A,J)
	  Enddo
	GO TO 30

	ENDIF
        ENDIF
	
		
50    ENDDO
	

!################## Center of Mass calculation ########################
	
!!	Do i = 1,nmol
!!	r(1,i)=0.d0
!!	r(2,i)=0.d0
!!	DO j=1,NA
!!	 ddx = rs(1,i,1)-rs(1,i,j)
!!	 ddy = rs(2,i,1)-rs(2,i,j)
!!        If(abs(ddx).GT.(0.5d0*x_box))rs(1,i,j)=
!!     :rs(1,i,j)-dsign(x_box,rs(1,i,j))
!!	If(abs(ddy).GT.(0.5d0*y_box)) rs(2,i,j)= 
!!     : rs(2,i,j)-dsign(y_box,rs(2,i,j))	 
	
!!	r(1,i) = r(1,i) + rs(1,i,j)/NA
!!	r(2,i) = r(2,i)+ rs(2,i,j)/NA
!!	Do k = 1,nmol
	 
!!!	  IF(r(1,k).gt.0.5d0*x_box) r(1,k)=r(1,k)-x_box
!!!	  IF(r(1,k).lt.(-0.5d0*x_box))r(1,k)=r(1,k)+x_box
!!!	  IF(r(2,k).gt.0.5d0*y_box) r(2,k)=r(2,k)-x_box
!!!	  IF(r(2,k).lt.(-0.5d0*y_box))r(2,k)=r(2,k)+x_box
	 		
!!!	enddo
		
!!!	Enddo
!!	open(122,file='/home/capuan/Desktop/center.dat')
!!	write(122,*) r(1,i),r(2,i)
!!!	enddo
!#####################################################################




         Gtstclus = 1
!############### Analysis of clusters ###############################	
	Do I = 1,nmol
	aux(I) = 0
	sizecluster(I) = 1
	
	aux4(I) = 0
	aux6(I) = 0
	aux7(I) = 0
	Rcmx(I) = 0 
	Rcmy(I) = 0
	Rg2x(I) = 0
	Rg2y(I) = 0
		
	Rat_Gira(I) = 0
	Rat_ave(I) = 0
!	jj(I) = 0
!	Q(I) = 0
!	aux5(I) = 0
	Do K=1,nmol	
	aux3(I,K) = 0
	enddo
	enddo
	ncluster = 0

	Do K1 = 1,nmol-1
	 Do K2 = K1+1,nmol

	 If(LL(K1).eq.LL(K2)) then
	

          PP = LL(K1)   
!!! 		ddx =  r(1,P) - r(1,K2)
!!!	        ddy =  r(2,P) - r(2,K2)
		
!!!	 If(abs(ddx).GT.(0.5d0*x_box))  then 
!!!	r(1,K2) = r(1,K2)-dsign(x_box,r(1,K2))
!!!	endif

!!!	 If(abs(ddy).GT.(0.5d0*y_box)) then
!!!	r(2,K2)=r(2,K2)-dsign(y_box,r(2,K2))
!!!	endif 

!	 ddx =  r(1,K1) - r(1,K2)
!	 ddy =  r(2,K1) - r(2,K2)
!	 ddr = ddx*ddx+ddy*ddy
!	Rtest = (RCL+sigma*NA)*(RCL+sigma*NA)	

!!	If(ddr.GT.Rtest) then
!!	r(1,K2)= r(1,K2)-dsign(x_box,r(1,K2))
!!	r(2,K2)=r(2,K2)-dsign(y_box,r(2,K2))
!!	endif	
	


	  If(aux(PP).eq.0) then ! this is to avoid extra counting

	  ncluster = ncluster + 1 !number of clusters 
	  aux(PP) = 1
!!!	  Rcmx(ncluster) = r(1,K1)
!!!	  Rcmy(ncluster) = r(2,K1)
	  raux_x(ncluster,K1) = r(1,K1)
	  raux_y(ncluster,K1) = r(2,K1)
	  
	  Endif
	  If((aux4(k2).eq.0)) then	
	   sizecluster(ncluster) = sizecluster(ncluster)+1 !number of particles belonging to a given 'ncluster' cluster
!!!	   Rcmx(ncluster) = Rcmx(ncluster) + r(1,K2)
!!!	   Rcmy(ncluster) = Rcmy(ncluster) + r(2,K2)
	   
	   raux_x(ncluster,K2) = r(1,K2)
	   raux_y(ncluster,K2) = r(2,K2)
	   aux3(ncluster,K1) = 1
	   aux3(ncluster,K2) = 1 
	   aux4(k2) = 1

	 Endif
	   
	Endif
	Enddo
	Enddo


!!	 Do I=1,nmol
!!		open(123,file='/home/capuan/Desktop/center2.dat')
!!	write(123,*) r(1,i),r(2,i)
!!	enddo
	Do k = 1,ncluster
	Do j = 1,nmol
	Do i = 1,nmol
	If(aux3(k,j).eq.1.AND.aux3(k,i).eq.1) then
	 	ddx =   raux_x(k,j) -  raux_x(k,i)
	        ddy =  raux_y(k,j) - raux_y(k,i)
	 
	IF (dabs(ddx).gt.(0.5d0*x_box)) ddx = ddx - dsign(x_box,ddx)
	IF (dabs(ddy).gt.(0.5d0*y_box)) ddy = ddy - dsign(y_box,ddy)	
	Rg2x(k) = Rg2x(k) + ddx*ddx
	Rg2y(k) = Rg2y(k) + ddy*ddy
	endif
  
	enddo
	enddo
	enddo



	  Polim=0
	  K1 = 1 	
	DO K=1,NCLUSTER
	
!	 WRITE(*,*) K,SIZECLUSTER(K)
	if (sizecluster(k).ne.0) then
!!!	 Rcmx(K) = Rcmx(K)/real(sizecluster(K))
!!!	 Rcmy(K) = Rcmy(K)/real(sizecluster(K))
		
	    Gtstclus = sizecluster(K1)
	
	  if(sizecluster(K).gt.Gtstclus) then
	   Gtstclus = sizecluster(K)
	    K1 = K
	  Endif

	  Polim = Polim + dble(sizecluster(K))/dble(nmol)
	  Rat_Gira(K) = Rg2x(k)+ Rg2y(k)
		
	  Do J=1,nmol
	  IF(aux3(K,J).eq.1) then
!!	  nome_arq_1 = 'confis_'						! arquivo de video
!!        write(nome_arq_2,'(I5.5)')K
	
!!        nome_arq_total = nome_arq_1 // nome_arq_2 

!!	 open(unit=102, file ='/home/jorge/Desktop/'
!!     ://nome_arq_total//'.dat')
!!	   write(102,*) raux_x(k,J),raux_y(K,J)
!!!	  Rat_Gira(K) = Rat_Gira(K)+(raux_x(K,J)-Rcmx(K))**2+ 
!!!     :(raux_y(K,J)-Rcmy(K))**2
	  
!!!		Q(K) = Q(K)+1
!!	Rat_Gira(K) = Rat_Gira(K)+(dsqrt(raux_x(K,J)**2+raux_y(K,J)**2)-
!!     :dsqrt(Rcmx(K)**2+Rcmy(K)**2))**2
	 endif
!	 write(*,*)'contado=',I,sizecluster(K)
	  enddo
!!	  close(102)832  
        Rat_Gira(K) = dsqrt(0.5d0*Rat_Gira(K))
	Rat_Gira(K) = Rat_Gira(K)/DBLE(sizecluster(K)) 
!	  write(*,*)sizecluster(K),Rat_Gira(K),Rcmx(K),Rcmy(K)
!	 WRITE(*,*) 'PARTÍCULAS EM CLUSTERS',I
	  endif	
	 ENDDO
	 Gtst = Gtst + Real(Gtstclus)/dble(nmol)

	 Gtst2 = Gtst2+Real(Gtstclus*Gtstclus)/dble(nmol*nmol)
!	 Write(*,*) 'media',Gtst/countt
	 
	Do k1 = 1,ncluster
	hh = sizecluster(k1)  
	Rat_ave(hh) = Rat_ave(hh)+Rat_Gira(k1)
	enddo
!!	IF(ncluster.gt.1) then
!!	Do k1= 1,ncluster-1
!!	 Do k2 = k1,ncluster
!!	 if(sizecluster(K1).eq.sizecluster(K2)) then
!!	  HH = sizecluster(K1)
!!	  if(aux6(HH).eq.0)Rat_ave(HH) = Rat_Gira(k1)
!!	 if(aux5(k2).eq.0) then
!!	 Rat_ave(HH)= Rat_ave(HH)+Rat_Gira(k2)
!!	  aux5(k2) = 1
!!	  jj(HH)= jj(HH) + 1
!	 write(*,*) jj(HH)
!!	  aux6(HH)=1
!!	 endif 
!!	  endif
	 
!!	 enddo
!!	 enddo
!!	Else
!!	 Rat_ave(sizecluster(ncluster)) = Rat_Gira(ncluster)
!!	endif
	
 	 
!!	write(*,*)'clusteres=',ncluster
	
		
	
	
	 	 
!########### How many clusters have 'n' particles? ##################
	Do I=1,nmol
        summ(I)= 0
	Enddo

	Do I=1,nmol
	DO J=1,ncluster
	 If(sizecluster(J).eq.I) summ(I)= summ(I)+1
	enddo
	enddo  
	Do I=1,nmol
	if(summ(I).ne.0) Rat_ave(I) = Rat_ave(I)/dble(summ(I))
  !      if(mod(kvid,500).eq.0.and.summ(I).ne.0)write(*,*)I,Rat_ave(I)
!	if(summ(I).ne.0) write(*,*) 'Existem', summ(I),'com',I,'particulas'   
	enddo
 	
	     

!####################################################################
!	Write(*,*) 'Numero de clusters=',ncluster
!	Do I=1,nmol
!	if(summ(I).ne.0) write(*,*) 'Existem', summ(I),'com',I,'particulas'   
!	if(summ(I).ne.0) Write(3,*) I,summ(I)	
!	Enddo
	
!################## average countings ###############################	
	Do k = 1,ncluster
	hh = sizecluster(k)
	if(aux7(hh).eq.0)countrat(hh) = countrat(hh) + 1 	
	aux7(hh) = 1
	
	enddo
!#####################################
	If(step.gt.stepaverage) then
	Polim2 = Polim2 + Polim*Polim
	
	Polimkvi = Polimkvi + Polim

	Polimspol= dsqrt(Polim2/dble(countt)- 
     :(Polimkvi/dble(countt))**2)
        Polimspol = Polimspol/dsqrt(dble(countt))
	nclukvi = nclukvi + ncluster
	endif

	Do K = 1,nmol
	sckvid(K) = sckvid(K)+ summ(K)	
	
	enddo
	Do K = 1,nmol
!!	hh = sizecluster(K)
	ratkvid(k) = ratkvid(k) + Rat_ave(k)
	enddo
		
!     Take a look at picture I took.




!!	Enddo ! Configurations

	if (mod(step,500).eq.0) write(632,*)step*deltaT,ncluster
	

	
	
	Gstpol= dsqrt(Gtst2/dble(countt)-(Gtst/dble(countt))**2)
        Gstpol = Gstpol/dsqrt(dble(countt))
	If(mod(step,10000).eq.0) then	
	open(354,file='greatestcluster.dat')
	open(35,file='cluster_histogramet.dat')
	open(59,file='polimerization.dat')
	open(60,file='ncluster.dat')
	
	open(51,file='Ratiogyration.dat')

		
	write(354,*) NA,Gtst/(countt),Gstpol
	write(59,'(I4,2E21.14)')NA,Polimkvi/dble(countt),Polimspol
	write(60,*)NA,nclukvi/dble(countt)

	Do k=1,nmol
	if(countrat(k).ne.0)write(35,*)k,
     :real(sckvid(k))/dble(countt)
	if(countrat(k).ne.0) write(51,*)k,ratkvid(k)/countrat(k)
	
	enddo
	
	close(354)
	close(35)
	close(59)
	close(60)
	close(51)	

	Endif
	
	IF(countcoran.eq.limitang) then
	countcoran2 = countcoran2 + 1
	Do i = 1,sizehist
	histsum = histsum + histang(i)
	enddo
	open(121,file='anglecor.dat')
	open(122,file='anglecortotal.dat')
	Do i = sizehist,1,-1
	histang(i) = histang(i)/dble(histsum)
	histang2(i) = histang2(i) + histang(i)
	sbin = 1.d0-dble(i)/dble(sizehist)	
	write(121,*)sbin,histang(i)
	write(122,*)sbin,histang2(i)/dble(countcoran2)
	enddo
	close(121)
	close(122)
	
	countcoran = 0
	endif
	
        END 
	
!##########################################################################################
!!!!!!!!!!!!!!Dynamical Properties !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!##########################################################################################


       Subroutine DifInit(cDiff1,nTimeWin1,TimeWindow1,Diffx,Diffy,xp
     :,cDiff2,nTimeWin2,TimeWindow2,yp,r,s,thetap,Difftheta,
     :Diffpa,Diffpe)
	implicit none
	include 'constants.f'
	integer n,k,nTimeWin1,TimeWindow1(10000),cDiff1(10000)
	integer nTimeWin2,TimeWindow2(10000),cDiff2(10000)
	Real*8 xp(10000,nmol),yp(10000,nmol),Diffx(10000),s(2,nmol)
	Real*8 Diffy(10000),r(ndim,nmol),magcorr(10000)
        Real*8 spy(60000,nmol),spx(60000,nmol),corrmi(60000)
        Real*8 thetap(10000,nmol),Difftheta(10000)
        Real*8 Diffpa(10000), Diffpe(10000),pi
        pi = 4.d0*atan(1.d0)
	write(*,*)"Initializing calculation of diffusion..."
	nTimeWin1 = int (log(dble(final_step-stepaverage) - 
     :2.d0*log(10.d0))/log(1.2d0))
	nTimeWin2 = int ((final_step-stepaverage)*0.0025d0) 
	do n=1,nTimeWin1
		TimeWindow1(n) = 10*int(10*1.2d0**n)
		
		Diffx(n) = 0.0d0
		Diffy(n) = 0.0d0
		cDiff1(n) = 0
                Difftheta(n) = 0.d0
                Diffpa(n) = 0.d0
                Diffpe(n) = 0.d0
                
		do k=1,nmol
			xp(n,k) = r(1,k)
			yp(n,k) = r(2,k)
                        thetap(n,k) = atan2(s(2,k),s(1,k))
         if(thetap(n,k).lt.0.d0) thetap(n,k) = thetap(n,k) + 2.d0*pi
                enddo
	enddo
        return
	end

!####################################################################
	Subroutine Diffusion(step,cDiff1,nTimeWin1,TimeWindow1,xp,yp,
     :cDiff2,nTimeWin2,TimeWindow2,r,rcross,Diffx,Diffy,s,thetap,
     :Difftheta,Diffpa,Diffpe)
	implicit none
	include 'constants.f'
      integer nTimeWin1,TimeWindow1(10000),n,k,cDiff1(10000),step
      integer nTimeWin2,TimeWindow2(10000),cDiff2(10000)
      Real*8 xp(10000,nmol),yp(10000,nmol),Diffx(10000),s(ndim,nmol)
      Real*8 Diffy(10000),xr,yr,rcross(ndim,nmol),y_box,r(ndim,nmol)
      Real*8 spy(60000,nmol),spx(60000,nmol),corrmi(60000)
      Real*8 magxx(10000),magyy(10000),magy(10000),magx(10000)
      Real*8 magcorr(10000),Difftheta(10000),tr,thetap(10000,nmol)
      Real*8 dpa,dpe,Diffpa(10000),Diffpe(10000),dxr,dyr,pi
        pi = 4.d0*atan(1.d0)
	y_box = x_box
!############## MSD #################################################
	if (mod(step,10).eq.0) then
	do n=1,nTimeWin1
	 if (mod(step,TimeWindow1(n)).eq.0) then
	  do k=1,nmol
		xr = rcross(1,k)*x_box + r(1,k)
		yr = rcross(2,k)*y_box + r(2,k)
                tr = atan2(s(2,k),s(1,k))
                if(tr.lt.0.d0) tr = tr + 2.d0*pi
                dxr = xr - xp(n,k)
                dyr = yr - yp(n,k)
		Diffx(n) = Diffx(n) + (dxr)**2.0d0
		Diffy(n) = Diffy(n) + (dyr)**2.0d0
                dpa = dxr*s(1,k) + dyr*s(2,k)
                dpe = dxr*s(2,k) - dyr*s(1,k)
                Diffpa(n) = Diffpa(n) + dpa**2.d0
                Diffpe(n) = Diffpe(n) + dpe**2.d0
        
            Difftheta(n) =  Difftheta(n) + (tr - thetap(n,k))**2.0d0
       		xp(n,k) = xr
		yp(n,k) = yr
               thetap(n,k) = tr
       	  enddo
	cDiff1(n) = cDiff1(n) + 1
         endif
        enddo
        endif


!####################################################################

!####################Correlation Functions###########################

	
!####################################################################
        if (mod(step,1000).eq.0) then

	  open(105,file="msd.dat")
      


	  do n=1,nTimeWin1
	  if (cDiff1(n).gt.0) then
	  write(105,'(7E15.6)')TimeWindow1(n)*deltaT,
     :Diffx(n)/nmol/cDiff1(n),Diffy(n)/nmol/cDiff1(n),
     :(Diffx(n)+Diffy(n))/nmol/cDiff1(n),Difftheta(n)/nmol/cDiff1(n)
     :,Diffpa(n)/nmol/cDiff1(n),Diffpe(n)/nmol/cDiff1(n)
	  endif
	  enddo



	 close(105)

         close(2151)

	endif	
        return
	end

        !############################################################
        Subroutine dif(switch,nsamp,ntime,ntel,t0,tt0,time0,
     :s,sx0,sy0,dtime,corrmi,Magtotal,magcor,rsx,rsy,nbondaver,
     :nbondcor)
        implicit none
        include 'constants.f'
        integer switch,nsamp,ntel,ntime(60000),t0,tt0,time0(60000)
        integer t,t0max,tmax,it0,i, delt,delt1,step
        Real*8 dtime,corrmi(60000),sx0(nmol,60000),rv(2,nmol)
        Real*8 sy0(nmol,60000),time,s(2,nmol),corrminew(60000)
        Real*8 magcor(60000), Magtotal,nbondcor(60000),nbondaver
        Real*8 rsx(NA,nmol), rsy(NA,nmol)
        Common /link2/step
        t0max = 450
        tmax = 4500
        it0 = 100
   
        ! t0max*it0 should not be smaller than tmax
 !       write(*,*)ntel,ntime(1),t0
!        pause
  20    If (switch.eq.0) then
        ntel = 0 !time counter
        dtime = deltaT*nsamp !time between two samples
    
              do 10 i =1,tmax  !tmax = total number of time step 
                ntime(i) = 0   !number of samples for time i
!                vacf(i) = 0
                corrmi(i) = 0.d0
   !             r2t(i) = 0
                magcor(i) = 0.d0
                nbondcor(i) = 0.d0
        
   10         enddo

          else if(switch.eq.1) then    !sample
                ntel = ntel + 1
          If(mod(ntel,it0).eq.0) then !decide to take a new t = 0
                t0 = t0 + 1             !update number of t = 0
                tt0 = mod(t0-1,t0max) + 1 !see note 1
                time0(tt0) = ntel    !store the time of t = 0
                  do 40 i =1,nmol
   !               x0(i,tt0) = x(i) 
                 sx0(i,tt0)= s(1,i) !store orientation for given t = 0
                 sy0(i,tt0)= s(2,i)
                 
   40             enddo
      !       call Net_Magnetization(Magtotal0,step)
      !          mag0(tt0) = Magtotal0
        Endif  
         Do 50 t=1,min(t0,t0max)   !update vacf for t=0
                delt=ntel-time0(t)+1 !actual time minus t = 0
                if(delt.lt.tmax) then 
                  ntime(delt)=ntime(delt)+1
                       !  write(*,*)delt 
                     
                  do i =1,nmol
!                  vacf(delt) = vacf(delt)+rv(1,i)*vx0(i,t)+  !update velocity correlation
!     :rv(2,i)*vy0(i,t) 
          corrmi(delt)=corrmi(delt)+ sx0(i,t)*s(1,i)+sy0(i,t)*s(2,i)
          
                  enddo
           call Net_Magnetization(Magtotal,s,step)
                magcor(delt) = magcor(delt) + Magtotal
                
           call nbonds(rsx,rsy,nbondaver)
                nbondcor(delt) = nbondcor(delt) + nbondaver 
                endif
       
  50    enddo
        Else if(switch.eq.2) then    !determine results
        open(106,file="corrmi.dat")
           Do i = 1,tmax
            time=dtime*(i+0.5)     ! time
        !    vacf(i)=vacf(i)/(nmol*ntime(i))  ! volume velocity autocorr.
             write(106,'(4E15.6)')time,corrmi(i)/(nmol*ntime(i)),
     :magcor(i)/(ntime(i)),nbondcor(i)/ntime(i)
             
           enddo
        close (106)
        Endif
        Return
        end
!####################################################################        
        Subroutine Net_Magnetization(Magtotal,s,step)
        implicit none
        include 'constants.f'
        REal*8 Mag,Magtotal,Magtotal2,Magerr,mit(ndim,nmol),s(2,nmol)
       
        integer countmag,i,step
        Magtotal = 0.d0
        Mag = 0.d0
    
!        countmag = countmag +1
        Do i = 1,nmol
        Mag = Mag + s(1,i)*dcos(omega*deltaT*(step)) +
     :s(2,i)*dsin(omega*deltaT*(step))
      
        enddo
        Mag = Mag/dble(nmol)
        
        Magtotal = Magtotal + Mag
        
!        Magnn(1) = Magnn(1) + Magn(1)/dble(nmol)
!        Magnn(2) = Magnn(2) + Magn(2)/dble(nmol)
       
!       Magtotal2 = Magtotal2 + Mag*Mag
        
!      Magerr = dsqrt(Magtotal2/countmag-(Magtotal/countmag)**2)
!     :/dsqrt(dble(countmag))
        
!        if (mod(countmag,100).eq.0) then
!        open(2306,file='magnetization.dat')
!        write(2306,*)B0, Magtotal/countmag,Magerr
!        close(2306)     
!        endif
        Return
        end
!####################################################################
       Subroutine nbonds(rsx,rsy,nbondaver)
       implicit none
       include 'constants.f'
       integer i,j,B1,B2,Nbond(nmol),aa(nmol,nmol)
       real*8 dx,dy,rsx(NA,nmol),rsy(NA,nmol),modr,rr,nbondaver
       do j = 1,nmol
       do i = 1,nmol
        aa(i,j) = 0
       enddo
      enddo
        nbondaver = 0.d0
        Do i =1,nmol
	Nbond(i) = 0
	enddo
        Do i =1,nmol-1
	  Do j= i+1,nmol
               Do B1 = 1,NA
	       Do B2 = 1,NA
	       dx = rsx(B1,i)-rsx(B2,j)
	       dy = rsy(B1,i)-rsy(B2,j)
                 
	    IF (dabs(dx).gt.(0.5d0*x_box)) dx = dx - dsign(x_box,dx)
	    IF (dabs(dy).gt.(0.5d0*x_box)) dy = dy - dsign(x_box,dy)         
     

	      
	     modr = dsqrt(dx*dx+dy*dy)
	     rr = dx*dx+dy*dy
   
         IF((aa(i,j).eq.0).and.(rr.lt.(1.6d0*1.6d0))) then
	    Nbond(i) = Nbond(i) + 1
	    Nbond(j) = Nbond(j) + 1
          aa(i,j) = 1  

	   Endif
        enddo
        enddo
        enddo
        enddo
        !If(step.gt.stepaverage) then
	Do i=1,nmol
        nbondaver = nbondaver + Nbond(i)/dble(nmol)
        enddo
	
!	endif
        Return
        end
