	subroutine acRBsinf(E,y,l1,j1,s1,l2,j2,s2,green,eps00)
		use f90_kind
		use constants
		use parametros
    !
    !     Calculate the gf props between 2-atoms on a semi-infinite ribbon with armchair edges
    ! 

    !        INPUT (integer) :: l1,j1,s1,l2,j2,s2 -- impurity locations
    !        INPUT (complex) :: E+iy -- complex energy
    !        OUTPUT (complex) :: green -- required propagator

    !     All sites described by triple (l.a1, j.a2,s)
    !     ========================================
    !     a1, a2 are unit vectors of graphene
    !     s=1,2 is the intracell index

    !     First impurity position described by (l1,l2,s1)
    !     Second impurity position described by (l2, j2, s2)

    !     Semi-infinite ribbon with armchair border, quantised along the direction (l,-l) i.e. (l,-l,1)   
    !     and has vanishing density of states in directions (l,l,s) and (N+l,-N+l,s) for s=1,2
    !     LDOS also vanishes along the ring (l,-l,1)

    !     (0,1,2) ______(1,1,1)
    !            /      \
    !           /        \
    !   (0,1,1)/          \______
    !          \          /      \
    !           \        /        \
    !     (0,0,2)\______/(1,0,1)   \
    !            /      \          /
    !           /        \        /
    !   (0,0,1)/          \______/
    !          \          /
    !           \        /    
    !            \______/
    implicit none

    !     declare inputs
    !     ---------------
    ! impurity locations
	integer, intent(in) :: l1, j1, l2, j2, s1, s2

	integer	:: nn
	integer	:: ds, dl, dj, sl, sj

!	complex energy E+i y
	real(double), intent(in)	:: E,y,eps00
    complex(double)	:: cw
!	required propagator	
    complex(double), intent(out)	:: green
    complex(double)	:: numerator, denom, denom2, sinn1, sinn2, k_n

	real(double)	:: arg
    complex(double)	:: cos_pole,sin_pole,exp_pole,exp_minus_pole,cosn,exp_2phi

	cw = cmplx(E,y,double)

!	impurity separation
    dl = l1-l2
    dj = j1-j2
    ds = s1-s2 
	sl = l1+l2
	sj = j1+j2

    green = zero

	do nn=1,2*(P+1)
		if ((nn.eq.(P+1))) then
			go to 30
		end if
		arg  = nn*pi/(2*(P+1)-1)
		cosn = cos(arg)

		cos_pole = (((cw-eps00)/t)**2 - 1)/(4*cosn) - cosn

!		k_n = acos(cos_pole) - imaginary arcosine, to check which pole is inside the countour
		k_n	= -zi*log(cos_pole+zi*sqrt(zum-(cos_pole**2)))

		if (aimag(k_n).lt.0.d0) then
			sin_pole	= -sqrt(zum-(cos_pole**2))
		else
			sin_pole	= sqrt(zum-(cos_pole**2))
		end if

		denom		= cosn*sin_pole
		sinn1		= sin(arg*(l1-j1))
		sinn2		= sin(arg*(l2-j2))
		exp_pole	= cos_pole+zi*sin_pole

		if((s1.eq.1).and.(s2.eq.1)) then
			green = green+sinn1*sinn2*((exp_pole**abs(dl+dj))-(exp_pole**abs(sl+sj)))/denom
		else if((s1.eq.2).and.(s2.eq.2)) then
			exp_minus_pole = (cos_pole-zi*sin_pole)
			exp_2phi = (1.d0+2.d0*cosn*exp_pole)/(1.d0+2.d0*cosn*exp_minus_pole)

			green = green+sinn1*sinn2*(exp_pole**abs(dl+dj)-exp_2phi*(exp_pole**abs(sl+sj)))/denom
		else
			numerator = (exp_pole**abs(dl+dj))-(exp_pole**abs(sl+sj))+(2.d0*cosn*((exp_pole**abs(dl+dj+ds))-(exp_pole**abs(sl+sj+1))))
			green = green + sinn1*sinn2*numerator/denom
		end if

30	continue
	end do

	if(ds.eq.0)then
	   green = zi*(cw-eps00)*green/(4.d0*(P+1)*t**2)
	else
	   green = zi*green/(4.d0*(P+1)*t)
	end if

	denom2	= ((cw-eps00)**2 - t**2)*(P+1)
	sinn1	= sin(pi*(l1-j1)/2.d0)
	sinn2	= sin(pi*(l2-j2)/2.d0)
	if((dl+dj.eq.0).and.(sl+sj.ne.0)) then
		if(ds.eq.0) then
			green	= green + (sinn1*sinn2*(cw-eps00)/denom2)
		else
			green	= green + (sinn1*sinn2*t/denom2)
		end if
	end if

	return
	end subroutine acRBsinf