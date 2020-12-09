import numpy as np
import math as m
import collections

#input file:
#x, y, z, R, dist, u, v, w, v_r, v_phi, v_x, v_y

# output file:
#e, r_max, r_min, R_max, R_min, Z_max, E, L_z, (L_x^2+L_y^2)^1/2, I_3, Z_max
#e : signed e (>0 for prograde, <0 for retrograde)
#r_max, r_min : maximum and minimum galactocentric radii (kpc)
#R_max, R_min : maximum and minimum distances along the galactic plane (kpc)
#Z_max : maximum distance away from the galactic plane (kpc)

#For unbound stars, the value '-9.999' is assigned.

#E : total energy ((km/s)^2) (<0 for bound stars)
#L_i : i-component of angular momentum (kpc km/s)
#I_3 : 3rd integral of motion ((kpc km/s)^2)

def run(req_dict, calc_dict):
	#nmdat = 200000

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global sah,c,sah2,c2,qh,rho_0,rt #halo
	global x,y,vx,vy,r,z,vr,vz,vp,lambda_,nu,vlam,vnu #phas
	global energy,I_2,I_3 #inte
	global eccentricity,elsecc,cradecc,abun #feoh
	global rmin,rmax,crmin,crmax,zeta,zmax #apoc

	global gmb, b, vesc #els

	#set basic parameters
	set_params()

	#set up tau array
	mesh()

	#set file names
	#read file

	#read data, derive eccentricities and write results

	unbound_count = 0

	#data = np.genfromtxt(input_file,delimiter=',',names=True,case_sensitive="lower",dtype=None)

	#nmdat = len(data)
	nmdat = 1
	
	
	x = calc_dict['x_gc']
	y = calc_dict['y_gc']
	z = calc_dict['z_gc']
	r = calc_dict['R']
	dist = req_dict['dist']
	vx = calc_dict['vx_gc']
	vy = calc_dict['vy_gc']
	vz = calc_dict['vz_gc']
	vr = calc_dict['vRg']
	vp = calc_dict['vTg']

	global is_bound
	is_bound = 1
	#if distance is negative, kick out
	if(dist <= 0.):
		is_bound = -1
		getecc()
		#all values will be "-9.999"

	#get orbital parameters from Staeckel potential

	#determine (lambda,nu) from (R,z)
	position()

	#determine (v_lambda,v_nu) from (vr,vz,lambda_,nu,z)
	transform_velocities()

	#determine (E,I_2,I_3)
	get_EI2I3()

	#get boundaries of orbits
	gbound()
	
	vtot=0.
	
	if(is_bound == 0):
		unbound_count = unbound_count + 1
		vtot = m.sqrt( vr**2 + vp**2 + vz**2 )
		#o6.write("unbound: {}\nvtot = {}\n".format(unbound_count,vtot))
		#print "unbound: {}\nvtot = {}\n".format(unbound_count,vtot)

	#get eccentricities & write 'ecc.out'
	getecc()

	if(vp >= 0.):
		rsignecc = eccentricity
	if(vp < 0.):
		rsignecc = eccentricity
	if(eccentricity < -9.):
		rsignecc = -9.999

	#o20.write("{},{},{},{},{},{}\n".format(rsignecc,rmax,rmin,crmax,crmin,zmax))

	#angular momentum
	rlz = r * vp
	rlx = y * vz - z * vy
	rly = z * vx - x * vz
	rlt = m.sqrt( rlx**2 + rly**2 )

	true_energy = -energy

	#o40.write("{},{},{},{},{}\n".format(true_energy,rlz,rlt,I_3,zmax))

	#o6.write("the number of stars employed={}\n".format(i-1))
	#o6.write("the number of unbound stars={}\n".format(unbound_count))
	
	bound = False
	if(is_bound==1):
		bound = True

	return_dict = collections.OrderedDict([('Energy',true_energy),('L_z',rlz),('L_p',rlt),('I_3',I_3),('Z_max',zmax),
	('ecc',rsignecc),('r_apo',rmax),('r_peri',rmin),('R_apo_P',crmax),('R_peri_P',crmin),('bound',bound)])
	return return_dict

#set the ELS potential
#see Chiba & Yoshii 1998 appendix
def set_params():

	global gmb, b, vesc #els
	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global sah,c,sah2,c2,qh,rho_0,rt #halo

	#ELS potential at solar dist. = -GM/(b + bq)
	#where q = m.sqrt( (rsun/b)**2 + 1. )
	#M is disk mass
	#b is scale length
	
	#can evaluate b & q using definitions of vesc and vsun

	vesc = 500.
	
	rsun = 8.2 #Bland-Hawthorne 2016
	vsun = 236. #Kawata 2019
	
	q = 1. / ( 1. - (m.sqrt(2.)*vsun/vesc)**2 )
	#c = m.sqrt( q**2 - 1.)
	b = rsun / m.sqrt( q**2 - 1.)
	
	#now can evaluate for gm/b
	gmb = (vsun**2) * q / (q - 1.)
	
	#gmb = vsun * (1. + q) * m.sqrt(q)/m.sqrt( q**2 - 1.)
	#gmb = gmb**2
	
#set Staeckel potential
	rt = 200.
	
	#disk parameters (perfect oblate disk)
	
	gamma_dummy  = 0.125
	del2 = 4.0**2
	neg_gamma = gamma_dummy**2
	neg_alpha = del2 + neg_gamma
	neg_alpha_sqrt  = m.sqrt( neg_alpha )
	q   = gamma_dummy / neg_alpha_sqrt
      
    #halo parameters (s=2 model)
    
	disk_mass = 9.0 * (10**10)
	rho_0  = 2.45 * (10**7)
	c  = 6.0
	c2 = c**2
	sah2 = del2 + c**2
	sah  = m.sqrt( sah2 )
	qh   = c / sah

#COMMENT	
def mesh():

	nm = 10001
	
	global tau_array #arry
	tau_array = [None] * nm
	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	
	tmax = m.log10(40000.)
	tmin = m.log10(neg_gamma + (1. * m.pow(10.,-4.)) * neg_gamma )
	dt = (tmax - tmin) / float(nm)
	
	#CHECK ALL nm RANGES BEFORE FINAL VER
	for i in range (0,nm):
		dum = tmin + dt * (i)
		tau_array[i] = 10.**dum

#get e_ELS, R_max, and R_min from ELS potential
# def gels(elsecc,elsrmax,elsrmin):
# 
# 	global x,y,vx,vy,r,z,vr,vz,vp,lambda_,nu,vlam,vnu #phas
# 	global gmb, b, vesc #els
# 
# 	#(1) energy and angular momentum
# 	
# 	elsene = felsene(vr,vp,r)
# 	elsang = felsang(vp,r)
# 	
# 	if(elsene >= 0.):
# 		#unbounded
# 		elsecc = -9.999
# 		elsrmax = -9.999
# 		return
# 
# 	#(2) R_min, R_max, and e_ELS
# 	
# 	#neg_alpha, beta, and gamma
# 	ralp2 = 1. + 4.*gmb * (b/elsang)**2 
# 	rgam  = (gmb - 2.*elsene) * (b/elsang)**2 / ralp2
# 	rbet2 = rgam**2 + 2.*elsene * (b/elsang)**2 / ralp2
# 	ralp = m.sqrt(ralp2)
# 	rbet = m.sqrt(rbet2)
# 	
# 	#R_min, R_max, and e_ELS
# 	r1 = m.sqrt( 1. - 2.*(rgam-rbet) ) * b / (rgam - rbet)
# 
# 	if(( rgam+rbet ) > 0.5):
# 		r2 = 0.
# 	else:
# 		r2 = m.sqrt( 1. - 2.*(rgam+rbet) ) * b / (rgam + rbet)
# 
# 	elsecc  = ( r1 - r2 ) / ( r1 + r2 )
# 	elsrmax = r1
# 	elsrmin = r2
# 	
# def felsene(vr,vp,r):
# 
# 	global gmb, b, vesc #els
# 	felsene = 0.5 * ( vr**2 + vp**2 ) - gmb / ( 1. + m.sqrt( r**2 / b**2 + 1. ) )
#      
# def felsang(vp,r):
# 
# 	felsang = r * vp
	
#G(tau) for disk	
def G_disk(t):

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global sah,c,sah2,c2,qh,rho_0,rt #halo

	G_grav = 4.3013 * (10**(-6))
	eps = 1. * (10**(-11))
	
	#perfect oblate disk
	
	dum = t - neg_gamma
	
	if(dum > eps):
		#from equ. (7) in Sommer-Larsen & Zhen
		G_disk = ((2.*G_grav*disk_mass) / (m.pi*m.sqrt(t-neg_gamma))) * m.atan(m.sqrt((t-neg_gamma)/neg_gamma))
	elif(dum <= eps):
		#small angle approx
		#arctan(theta)~theta, so can simplify equation
		G_disk = 2.*G_grav*disk_mass/m.pi/m.sqrt(neg_gamma)
		
	return G_disk

#G(tau) for halo (0 at origin)	
def fgh(t):

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global sah,c,sah2,c2,qh,rho_0,rt #halo
	
	G_grav = 4.3013 * (10**(-6))
	eps = 1. * (10**(-11))

	#s=2 model halo of de Zeeuw et al. (1986)
	
	dum = t - neg_gamma
	
	#possible sign error below, check later
	
	if(dum > eps):
	
		b = del2 + c2 - neg_alpha
		
		#from second line of equ. (9) from Sommer-Larsen & Zhen
		dum1 = m.log( (del2+neg_gamma+b)/(neg_gamma+b) ) - (t-neg_gamma+del2)/2./(t-neg_gamma) * m.log( (t+b)/(neg_gamma+b) )
		
		#from third & fourth lines of equ. (9) from Sommer-Larsen & Zhen
		dum2 = 1./m.sqrt(t-neg_gamma) * m.atan( m.sqrt( (t-neg_gamma)/(neg_gamma+b) ) ) - 1./m.sqrt(del2) * m.atan( m.sqrt(del2)/m.sqrt(neg_gamma+b) )
		
		#combine the above two parts...
		fgh = dum1 + (del2-neg_gamma-b)/m.sqrt(neg_gamma+b) * dum2
		
		#...and multiply by the terms in the first line of equ. (9) from Sommer-Larsen & Zhen
		fgh = - 4.*m.pi*G_grav*rho_0 * (-neg_gamma-b) * fgh
	
	elif(dum <= eps):
		#use small angle approx & Taylor expansion of ln to simplify fgh
		b = del2 + c2 - neg_alpha
		dum1 = m.log( (del2+neg_gamma+b)/(neg_gamma+b) ) - del2/2./(neg_gamma+b)
		dum2 = 1./m.sqrt(neg_gamma+b) - 1./m.sqrt(del2) * m.atan( m.sqrt(del2)/m.sqrt(neg_gamma+b) )
		fgh = dum1 + (del2-neg_gamma-b)/m.sqrt(neg_gamma+b) * dum2
		fgh = - 4.*m.pi*G_grav*rho_0 * (-neg_gamma-b) * fgh
		
	return fgh

#G(tau) for halo (0 at the cutoff radius)		
def G_halo(t):

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global sah,c,sah2,c2,qh,rho_0,rt #halo
	
	lambda_d = rt**2 + neg_alpha
	
	cons = fgh(lambda_d) + G_disk(lambda_d)

	G_halo = fgh(t) - cons
	
	return G_halo

#G(tau) for disk + halo
#from equ. (6) in Sommer-Larsen & Zhen	
def G(t):

	G = G_disk(t) + G_halo(t)
	
	return G

#dG(tau)/dtau for potential	
def dG_dt(t):

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global sah,c,sah2,c2,qh,rho_0,rt #halo
	#perfect oblate disk + s=2 model halo of de Zeeuw et al. (1986)
	
	G_grav = 4.3013 * m.pow(10.,-6.)
	
	#this is probably a mistake
	#dG_disk_dt = m.sqrt( t - neg_gamma ) / gamma_dummy * ( 1./(1.+dG_disk_dt**2) - atan(dG_disk_dt) / dG_disk_dt )


	#debugging
	#print "t = {}".format(t)
	#print "gamma_dummy = {}".format(gamma_dummy)
	
	dG_disk_dt = m.sqrt(t-neg_gamma)/gamma_dummy
	dG_disk_dt = (G_grav*disk_mass/(m.pi*(gamma_dummy**3)*(dG_disk_dt**2))) * ( (1./(1.+(dG_disk_dt**2))) - (m.atan(dG_disk_dt)/dG_disk_dt) )

	b = sah2 - neg_alpha
	p = m.sqrt( (t-neg_gamma) / (neg_gamma+b) )
	dum1 = del2/(t-neg_gamma) * m.log( (t+b)/(neg_gamma+b) ) - (t-neg_gamma+del2)/(t+b)
	dum2 = (del2-neg_gamma-b)/(neg_gamma+b) * ( 1./(1.+p**2) - m.atan(p) / p )
	dG_halo_dt = 2.*m.pi*G_grav*rho_0 * (neg_gamma+b)/(t-neg_gamma) * ( dum1 + dum2 )


	dG_dt = dG_disk_dt + dG_halo_dt
	
	return dG_dt

#Phi(tau) for potential
#from equ. (5) in Sommer-Larsen & Zhen
def phi(lambda_,nu):

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk

	phi = -((lambda_-neg_gamma) * G(lambda_) - (nu-neg_gamma) * G(nu)) / (lambda_ - nu)
	
	return phi

#for B(tau)
#from equ. (13) in Sommer-Larsen & Zhen?      
def B(t):

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global energy,I_2,I_3 #inte
	
	#debugging
# 	print "t={}".format(t)
# 	print "neg_alpha={}".format(neg_alpha)
# 	print "E={}".format(energy)
# 	print "neg_gamma={}".format(neg_gamma)
# 	print "I2={}".format(I_2)
# 	print "I3={}".format(I_3)
# 	print "G(t)={}".format(G(t))
	
	B = - (t-neg_alpha)*(t-neg_gamma) * energy - (t-neg_gamma) * I_2 - (t-neg_alpha) * I_3 + (t-neg_alpha)*(t-neg_gamma) * G(t)

	return B

#for dB(tau)/dtau
def dB_dt(t):

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global energy,I_2,I_3 #inte
	
	dB_dt = - (2.*t-neg_alpha-neg_gamma) * energy - I_2 - I_3 + (2.*t-neg_alpha-neg_gamma) * G(t) + (t-neg_alpha)*(t-neg_gamma) * dG_dt(t)

	return dB_dt

#determine (lambda_,nu) from (R,z) (*cylindrical* coords)
#start w/ 1 = R**2/(tau+alpha) + z**2(tau+gamma)
#lambda and nu are the roots of this equation
#rearrange as A*tau**2 + B*tau + C = 0
#solve quadratic equation to get lambda and nu
def position():

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global x,y,vx,vy,r,z,vr,vz,vp,lambda_,nu,vlam,vnu #phas

	#A = 1, unneeded
	neg_B = neg_alpha + neg_gamma + r**2 + z**2
	C = neg_alpha*neg_gamma + neg_gamma*r**2 + neg_alpha*z**2
	lambda_ = 0.5 * ( neg_B + m.sqrt( neg_B**2 - 4.*C ) )
	nu  = 0.5 * ( neg_B - m.sqrt( neg_B**2 - 4.*C ) )

#determine (v_lam,v_nu) for given (vr,vz,lambda_,nu,z)
#literal black magic	
def transform_velocities():

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global x,y,vx,vy,r,z,vr,vz,vp,lambda_,nu,vlam,vnu #phas

	cos_theta = (lambda_-neg_gamma)*(nu-neg_alpha)/(neg_gamma-neg_alpha)/(lambda_-nu)
	cos_theta = m.sqrt( cos_theta )
	sin_theta = (lambda_-neg_alpha)*(nu-neg_gamma)/(neg_alpha-neg_gamma)/(lambda_-nu)
	sin_theta = m.sqrt( sin_theta )
	
	#this is how sign funct. works in fortran
	def sign(a,b):
		c=None
		if(b>=0.):
			c = abs(a)
		elif(b<0.):
			c = -abs(a)
		return c

	vlam =   cos_theta * vr + sin_theta * sign(1.0,z) * vz
	vnu  = - sin_theta * vr + cos_theta * sign(1.0,z) * vz

#determine (E,I_2,I_3) from (v_lam,v_p,v_nu,lambda_,nu) and (v_x,v_y,v_z,x,y,z,R)	
def get_EI2I3():

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global x,y,vx,vy,r,z,vr,vz,vp,lambda_,nu,vlam,vnu #phas
	global energy,I_2,I_3 #inte

	#negative of the Hamiltonian, from equ. (10) & equ. (12) in Sommer-Larsen & Zhen
	energy = 0.5 * (vlam**2 + vp**2 + vnu**2) + phi(lambda_,nu)
	energy = - energy

	#from equ. (14) and equ. (12) in Sommer-Larsen & Zhen
	I_2  = 0.5 * ( r * vp )**2

	L_x = y * vz - z * vy
	L_y = z * vx - x * vz
	
	#from equ. (15) in Sommer-Larsen & Zhen
	I_3  = 0.5 * (L_x**2 + L_y**2) + del2 * ( 0.5 * vz**2 - z**2 * ( G(lambda_) - G(nu) )/(lambda_ - nu) )
	
def wrfb():

	nm = 10001
	
	global tau_array
	
	for i in range(0,nm):
		tau = tau_array[i]
		
		bt  = B(tau)
		bt  = bt / (10**6)
		o20.write("{},{}".format(tau,bt))


#get boundaries of orbits		
def gbound():

	nm = 10001
	
	global tau_array
	global x,y,vx,vy,r,z,vr,vz,vp,lambda_,nu,vlam,vnu #phas
	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global energy,I_2,I_3 #inte	
	global atsol #boun
	
	global is_bound
	
	#at_solution?
	atsol = [None] * 3
	
	B_array = [None] * nm
	t_initial_array = [None] * 3
	
	#small numbers?
	epsilon_array = [10**(-5),10**(-3),10**(-3)]
	
	#(0) skip for unbounded orbits
	
	if(energy < 0.):
		is_bound = 0
		return
	else:
		is_bound = 1
	
	#(1) search three nearest solutions: t_initial_array(3)
	
	for i in range(0,nm):
		tau = tau_array[i]
		B_array[i] = B(tau)
		
	j_initial = 0
	
	#CHECK THIS RANGE
	for i in range(0,nm-1):

		#debugging
		#print "loop {}/{}".format(i,nm)
		
		dum = B_array[i] * B_array[i+1]
		
		#if B transitioning from + to - or from - to +
		if(dum <= 0.):
			j_initial = j_initial + 1
			if((j_initial-1)==0):
				t_initial_array[j_initial-1] = tau_array[i]
			elif((j_initial-1)==1):
				t_initial_array[j_initial-1] = tau_array[i+1]
			else:
				t_initial_array[j_initial-1] = tau_array[i+1]
			#only look for 3 transitions/"zero points"
			if((j_initial-1)==2):
				break
	
	#if B doesn't transition 3 times, star is unbound			
	if(j_initial < 3):
		#print("failed for 3 nearest solutions\n")
		#print("	R = {}, z = {}".format(r,z))
		is_bound = 0
		return
	
	#(2) search 3 exact solutions by Newton method
	
	for j in [0,1,2]:
				
		t_initial = t_initial_array[j]
		
		#CHECK THIS RANGE
		for k in range(1,200):
			t_next = t_initial - B(t_initial) / dB_dt(t_initial)
			dum = abs(t_next - t_initial) / t_next
			if(dum < epsilon_array[j]):
				atsol[j] = t_next
				break
			else:
				if(k < 200):
					t_initial = t_next
				else:
					o6.write("failed for exact solutions, for j=".format(j))
					is_bound = 0
	
#get eccentricities in the Staeckel potential			
def getecc():

	global neg_alpha_sqrt,gamma_dummy,neg_alpha,neg_gamma,q,del2,disk_mass #disk
	global eccentricity,elsecc,cradecc,abun #feoh
	global rmin,rmax,crmin,crmax,zeta,zmax #apoc
	global atsol #boun
	
	global is_bound

	#unbound
	if( is_bound <= 0. ):
		eccentricity  = -9.999
		elsecc  = -9.999
		cradecc = -9.999
		rmax = -9.999
		crmax = -9.999
		rmin = -9.999
		crmin = -9.999
		zeta = -9.999
		zmax = -9.999
		return
		
	nu0 = atsol[0]
	lambda_1 = atsol[1]
	lambda_2 = atsol[2]
	
	if(nu0 <= neg_gamma):
		#print("nu0 <= neg_gamma")
		pass
		
	if(nu0 >= lambda_1):
		#print("nu0 >= lambda_1")
		eccentricity  = -9.999
		elsecc  = -9.999
		cradecc = -9.999
		rmax = -9.999
		crmax = -9.999
		rmin = -9.999
		crmin = -9.999
		zeta = -9.999
		zmax = -9.999
		return
		
	#(1) eccentricity in r
	rmax = m.sqrt(( lambda_2 - neg_alpha ) + ( nu0 - neg_gamma ))

	try:
		rmin = m.sqrt(lambda_1 - neg_alpha)
	except:
		rmin = -9.999

	eccentricity = ( rmax - rmin ) / ( rmax + rmin )
	
	#(2) eccentricity in R
	crmax = m.sqrt(lambda_2 - neg_alpha)

	try:
		crmin = m.sqrt(lambda_1 - neg_alpha)
	except:
		crmin = -9.999

	cradecc = ( crmax - crmin ) / ( crmax + crmin )
	
	#(3) width in the nu direction
	zeta = m.sqrt( nu0 - neg_gamma )
	
	#this is probably a mistake
	#zmax = m.sqrt(( lambda_2 - neg_gamma ) * ( nu0 - neg_gamma ))
	
	zmax = m.sqrt((( lambda_2 - neg_gamma ) * ( nu0 - neg_gamma )) / (neg_alpha - neg_gamma))
