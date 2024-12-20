	integer nmol
         Parameter (nmol = 1500) !Number of rods
	integer NA
         Parameter (NA = 2 ) !Number of atoms per molecule
        Real*8 sigma	
	Parameter (sigma = 1.d0) ! particle diameter
	integer ndim
	Parameter (ndim = 2) !Number of dimensions
	Real*8 piap
	Parameter (piap = 3.14159265359)!~Pi
	Real*8 rodareatot
	Parameter (rodareatot = NA*nmol*piap*(sigma*0.5)**2) !Number of dimensions
	integer final_step
	Parameter (final_step = 7.d6 ) ! Number os steps
	integer stepaverage
	Parameter (stepaverage = 1.d6) ! Steps of Relaxation
	Real*8 areafrac
	Parameter (areafrac = 1.d-1 )!area fraction
	Real*8 x_box
	Parameter (x_box=dsqrt(rodareatot/areafrac)) ! length of the box {-x_box/2, x_box/2} 
	Real*8 omega
	Parameter (omega = 30.d0)
	Real*8 B0
	Parameter (B0 = 10.d0)
	Real*8 temp
	Parameter (temp= 1.d0) ! temperature
	Real*8 deltaT
	Parameter (deltaT = 5.d-4) ! Time interval
	integer steptemp
	Parameter (steptemp=10 ) ! Number of times of steps which temperature is adjusted in 
	Real*8 mi
         Parameter (mi = 660.d0/100.d0)  ! Dipole moment
	integer P
	Parameter (P = 2)  ! video? Yes = 1, Otherwise = any value
	Integer configtest
	Parameter (configtest = 100)
	Real*8 mult
	Parameter (mult = 2.d0)
