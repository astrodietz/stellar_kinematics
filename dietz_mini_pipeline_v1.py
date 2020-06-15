import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from galpy.util import bovy_coords as b_c
import galpy.potential as pot
from astropy import units as u
import numpy as np
from scipy.stats import norm
import math as m
import os
import datetime
import collections
import staeckel_orbit
import sys
from sys import argv, stdout

input_file = argv[1]
main_title = argv[2]
dir_name = argv[3]

final_dir_name = "{}/{}_kin".format(dir_name,main_title)

if not os.path.exists(final_dir_name):
    os.makedirs(final_dir_name)
    
final_dir_name = final_dir_name + "/"

output_file = "{}{}_kin.csv".format(final_dir_name,main_title)

f = open(input_file)
o = open(output_file, 'w')

#create array from csv file
a1 = np.genfromtxt(f, autostrip=True, delimiter=",", names=True, case_sensitive="lower", dtype=None)

#if using pipeline with condor
if(len(sys.argv)==6):
	start = int(argv[4])
	stop = int(argv[5])
else:
	start = 0
	stop = len(a1)

#condor option
if(stop > len(a1)):
	stop = len(a1)

loop_length = stop - start

global Xsun
Xsun = 8.
global Zsun
Zsun = 0.
global vcirc
vcirc = 220.
global vsun
vsun1 = [-9.0, 12.0 + vcirc, 7.0]

global vsun2
vsun2=[-9.,12.,7.]

global N_orbits
N_orbits = 1000
print "\n"

#returns (heliocentric & galactocentric) x,y,z
def xyz(req_dict):

	#get heliocentric position
	x_hc, y_hc, z_hc = b_c.lbd_to_XYZ(req_dict['l'], req_dict['b'], req_dict['dist'], degree=True)

	#get galactocentric position
	x_gc, y_gc, z_gc = b_c.XYZ_to_galcenrect(x_hc, y_hc, z_hc, Xsun, Zsun)
			
	return float(x_gc), y_gc, float(z_gc)

#returns (heliocentric & galactocentric) cartesian velocities
def uvw(req_dict):

	#get vx, vy, vz (heliocentric) in R-handed coord. syst.
	u, v, w = b_c.vrpmllpmbb_to_vxvyvz(req_dict['rv'], req_dict['pml'], req_dict['pmb'], req_dict['l'], req_dict['b'], req_dict['dist'], XYZ=False, degree=True)
	
	#get vx, vy, vz (galactocentric) in L-handed coord syst (corrected for sun & disk motion)
	vx_gc, vy_gc, vz_gc = b_c.vxvyvz_to_galcenrect(u, v, w, vsun1, Xsun, Zsun)
		
	return u, v, w, float(vx_gc), vy_gc, float(vz_gc)

#returns (L-handed) galactocentric cylindrical coords	
def cyl_coords(calc_dict):

	R, phi, z = b_c.rect_to_cyl(calc_dict['x_gc'],calc_dict['y_gc'],calc_dict['z_gc'])
	
	return R, phi

#returns (L-handed) galactocentric cylindrical velocities
def cylindrical_vs(calc_dict):
		
	vRg, vTg, vZg = b_c.rect_to_cyl_vec(calc_dict['vx_gc'],calc_dict['vy_gc'],calc_dict['vz_gc'],calc_dict['x_gc'],calc_dict['y_gc'],calc_dict['z_gc'])
	
	return vRg, vTg

#returns uncertainties in u,v,w
def uvw_unc(req_dict, opt_dict):

	cov_radec=np.zeros((2,2))
	cov_radec[0,0]=opt_dict['epmra']**2
	cov_radec[1,1]=opt_dict['epmdec']**2

	covar_pmllbb=b_c.cov_pmrapmdec_to_pmllpmbb(cov_radec, req_dict['ra'], req_dict['dec'], degree=True, epoch=2000.0)

	cov_vxvyvz=b_c.cov_dvrpmllbb_to_vxyz(req_dict['dist'], opt_dict['edist'], opt_dict['erv'], req_dict['pml'], req_dict['pmb'], covar_pmllbb, req_dict['l'], req_dict['b'] , plx=False, degree=True)

	du=vx_e=m.sqrt(cov_vxvyvz[0,0])
	dv=vy_e=m.sqrt(cov_vxvyvz[1,1])
	dw=vz_e=m.sqrt(cov_vxvyvz[2,2])
		
	return du, dv, dw

#returns a dict of orbital parameters from staeckel_orbit.py
def orbit_1_staeckel(req_dict, calc_dict):

	returned_orb = staeckel_orbit.run(req_dict, calc_dict)
	
	if(-9.999 in returned_orb):
		returned_orb['bound']=False
	
	return returned_orb
	
#returns orbital parameter uncertainties
def orbit_N(method,req_dict, opt_dict, calc_dict):

	count_orbits=0.

	orbit_params={'Energy':[], 'L_z':[], 'L_p':[], 'I_3':[], 'Z_max':[],
	'ecc':[], 'r_apo':[], 'r_peri':[], 'R_apo_P':[], 'R_peri_P':[]}

	orbit_unc={'Energy_unc':-9.999, 'L_z_unc':-9.999, 'L_p_unc':-9.999, 'I_3_unc':-9.999, 'Z_max_unc':-9.999,
	'ecc_unc':-9.999, 'r_apo_unc':-9.999, 'r_peri_unc':-9.999, 'R_apo_P_unc':-9.999, 'R_peri_P_unc':-9.999, 'N_orbits':-9.999}
			
	#loop over N orbits
	for i in range(0,N_orbits):
	
		#generate input parameters by randomly selecting from a normal distribution with spread = input uncertainty
		rand_dict = req_params.copy()
	
		rand_dict['dist'] = np.random.normal(loc=req_dict["dist"], scale=opt_dict["edist"])
		rand_dict['rv'] = np.random.normal(loc=req_dict["rv"], scale=opt_dict["erv"])
		
		rand_dict['pmra'] = np.random.normal(loc=req_dict["pmra"], scale=opt_dict["epmra"])
		rand_dict['pmdec'] = np.random.normal(loc=req_dict["pmdec"], scale=opt_dict["epmdec"])
		
		rand_dict['pml'], rand_dict['pmb'] = b_c.pmrapmdec_to_pmllpmbb(rand_dict['pmra'], rand_dict['pmdec'], req_dict['ra'], req_dict['dec'], degree=True, epoch=2000.0)
		
		#use randomly generated input parameters to calculate positions, velocities
		rand_x, rand_y, rand_z = xyz(rand_dict)
		rand_u, rand_v, rand_w, rand_vx_gc, rand_vy_gc, rand_vz_gc = uvw(rand_dict)
		
		rand_calc_dict=collections.OrderedDict([('x_gc',rand_x),('y_gc',rand_y),('z_gc',rand_z),
		('u',rand_u),('v',rand_v),('w',rand_w),
		('vx_gc',rand_vx_gc),('vy_gc',rand_vy_gc),('vz_gc',rand_vz_gc),
		('R',-9.999),('phi',-9.999),('vRg',-9.999),('vTg',-9.999)])
		
		rand_calc_dict['R'], rand_calc_dict['phi'] = cyl_coords(rand_calc_dict)
		rand_calc_dict['vRg'], rand_calc_dict['vTg'] = cylindrical_vs(rand_calc_dict)
		
		returned_orb = orbit_1_staeckel(rand_dict,rand_calc_dict)
			
		for key in returned_orb:
			if(returned_orb['bound']==True and key in orbit_params):
				orbit_params[key].append(returned_orb[key])

		if(returned_orb['bound']==True):
			count_orbits+=1.
			
	orbit_unc['N_orbits'] = count_orbits
	
	fits_dir_name = "{}{}_uncertainty_fits".format(final_dir_name,method)
	
	if not os.path.exists(fits_dir_name):
		os.makedirs(fits_dir_name)
    
	fits_dir_name = fits_dir_name + "/"
	
	f, axarr = plt.subplots(5,2, figsize=(12, 10))
	
	row_num = 0
	col_num = 0
	
	for key in orbit_params:
	
		data = orbit_params[key]

		#fit a gaussian to the randomly generated orbital parameters
		mu, std = norm.fit(data)
		
		#uncertainty on parameter = standard dev of fit
		orbit_unc['{}_unc'.format(key)]=std

		#plot the histogram
		axarr[row_num,col_num].hist(data, bins=25, normed=True, facecolor='none', edgecolor="black", histtype="step")

		#plot the fit
		xmin, xmax = axarr[row_num,col_num].get_xlim()
		x = np.linspace(xmin, xmax, 100)
		p = norm.pdf(x, mu, std)
		axarr[row_num,col_num].plot(x, p, 'k', linewidth=2)
		
		#plot the calculated parameter value
		axarr[row_num,col_num].axvline(x=calc_dict[key], color='b', ls ='--', lw=2)
		
		#plot the mean of the fit
		axarr[row_num,col_num].axvline(x=mu, color='k', lw =2)
		
		axarr[row_num,col_num].set_title("{} (mu = {}, std = {})".format(key, round(mu,2), round(std,2)))
		
		if(col_num==0):
			col_num+=1
		else:
			col_num=0
			row_num+=1		
			
	plt.tight_layout()
	plt.savefig("{}{}_{}_{}_fits.png".format(fits_dir_name,main_title,req_dict['name'],method))
	plt.close()

	return orbit_unc
	
def calculate():

	global req_params
	global opt_params
	global calc_params

	bad_input = False
	
	#check for invalid required params
	for key in req_params:
		if(key != "name"):
			if((m.isnan(float(req_params[key]))==True) or (req_params[key]==-9.999) or (req_params[key]==-9999.9)):
				bad_input = True
			if(isinstance(req_params[key], str)==True):
				bad_input = True
	
	#valid required params, continue to calculations
	if(bad_input==False):		
		#get heliocentric & galactocentric x,y,z
		calc_params['x_gc'], calc_params['y_gc'], calc_params['z_gc'] = xyz(req_params)

		#calculate galactocentric distance
		calc_params['R_gal']=np.sqrt(calc_params['x_gc']**2+calc_params['y_gc']**2+calc_params['z_gc']**2)

		#get heliocentric & galactocentric cartesian velocities	
		calc_params['u'], calc_params['v'], calc_params['w'], calc_params['vx_gc'], calc_params['vy_gc'], calc_params['vz_gc'] = uvw(req_params)
	
		#get cylindrical coordinates
		calc_params['R'], calc_params['phi'] = cyl_coords(calc_params)
		
		#get cylindrical velocities
		calc_params['vRg'], calc_params['vTg'] = cylindrical_vs(calc_params)
		
		staeckel_orb = orbit_1_staeckel(req_params,calc_params)
		for key in calc_params:
			if(key in staeckel_orb):
				calc_params[key]=staeckel_orb[key]
							
		uncert_list = [opt_params['epmra'],opt_params['epmdec'],opt_params['erv'],opt_params['edist']]
		
		calc_uncert = True
		
		#check for invalid uncertainty params
		for err in uncert_list:
			if((m.isnan(float(err))==True) or (err==-9.999) or (err==-9999.9)):
				calc_uncert = False
			if(isinstance(err, str)==True):
				calc_uncert = False
	
		#valid uncertainty params, continue to calculations
		if(calc_uncert==True):
			staeckel_orb_unc = orbit_N("staeckel",req_params,opt_params,calc_params)
			for key in calc_params:
				if key in staeckel_orb_unc:
					calc_params[key]=staeckel_orb_unc[key]
	
	#invalid required params	
	if(bad_input==True):
		no_xyz = False
		for key in ['l','b','dist']:
			if((m.isnan(float(req_params[key]))==True) or (req_params[key]==-9.999)):
				no_xyz = True
			if(isinstance(req_params[key], str)==True):
				no_xyz = True
					
		if(no_xyz==False):
			#get galactocentric x,y,z
			calc_params['x_gc'], calc_params['y_gc'], calc_params['z_gc'] = xyz(req_params)
			#calculate galactocentric distance
			calc_params['R_gal']=np.sqrt(calc_params['x_gc']**2+calc_params['y_gc']**2+calc_params['z_gc']**2)
			#get cylindrical coordinates
			calc_params['R'], calc_params['phi'] = cyl_coords(calc_params)
	
	if(-9.999 not in [opt_params['edist'],opt_params['epmra'],opt_params['epmdec']]):
		#get uncertainties in u,v,w
		calc_params['du'], calc_params['dv'], calc_params['dw'] = uvw_unc(req_params, opt_params)		
		
	for key in req_params:
		o.write("{},".format(req_params[key]))
	for key in opt_params:
		o.write("{},".format(opt_params[key]))
	for key in calc_params:
		o.write("{},".format(calc_params[key]))
		
	#get rid of that last comma & add a new line
	o.seek(-1, os.SEEK_END)
	o.truncate()
	o.write("\n")
			
#function calls all conversion sub-functions
def convert():

	global req_params
	
	#if position given in ra & dec, convert to l & b
	if(req_params.get('l')==-9.999):
		req_params['l'], req_params['b'] = b_c.radec_to_lb(req_params['ra'], req_params['dec'], degree=True, epoch=2000.0)
	
	#if position given in l & b, convert to ra & dec
	if(req_params.get('ra')==-9.999):
		req_params['ra'], req_params['dec'] = b_c.lb_to_radec(req_params['l'], req_params['b'], degree=True, epoch=2000.0)
	
	#if proper motion given in pmra & pmdec, convert to pml & pmb
	if(req_params.get('pml')==-9.999):
		req_params['pml'], req_params['pmb'] = b_c.pmrapmdec_to_pmllpmbb(req_params['pmra'], req_params['pmdec'], req_params['ra'], req_params['dec'], degree=True, epoch=2000.0)

#cycle through rows of data, perform conversions & calculations
def loop():

	global req_params
	global opt_params
	global calc_params
	
	#write output file header
	for key in req_params:
		o.write("{},".format(key))
	for key in opt_params:
		o.write("{},".format(key))
	for key in calc_params:
		o.write("{},".format(key))
		
	#get rid of that last comma & add a new line
	o.seek(-1, os.SEEK_END)
	o.truncate()
	o.write("\n")
	
	#function that lets you rename dict keys
	def change_key(dict, old_key, new_key):
		dict[new_key] = dict.pop(old_key)

#---start looping through stars
	
	#for row in range(0,loop_length):
	for row in range(start,stop):	

		#print a progress statement
		stdout.write("\r{}/{}: {}".format(row+1,loop_length,a1[row]['name']))
		stdout.flush()
		
		reset()
	
#-------start filling param dicts	

		for param_key in req_params:
			try:
				#if star has the required param, store it
				req_params[param_key] = a1[row][param_key]
			except:
				#if not, move on
				pass
			
		for param_key in opt_params:
			try:
				#if star has the optional param, store it
				opt_params[param_key] = a1[row][param_key]
			except:
				#if not, move on
				pass
					
#-------end filling param dicts
		
		req_params['pmra'] = float(req_params['pmra'])
		opt_params['epmra'] = float(opt_params['epmra'])
		req_params['pmdec'] = float(req_params['pmdec'])
		opt_params['epmdec'] = float(opt_params['epmdec'])
		req_params['pml'] = float(req_params['pml'])
		req_params['pmb'] = float(req_params['pmb'])

 		convert()
 		calculate()
						
#---end looping through stars

print("\noriginal number of stars: {}\n".format(loop_length))

#param dicts
global req_params
global opt_params
global calc_params

#reset param dict vals
def reset():

	global req_params
	global opt_params
	global calc_params

	#params that must be given
	req_params=collections.OrderedDict([('name','None'),
	('ra',-9.999), ('dec',-9.999),
	('pmra',-9.999), ('pmdec',-9.999),
	('l',-9.999), ('b',-9.999),
	('pml',-9.999), ('pmb',-9.999),
	('rv',-9.999), ('dist',-9.999)])
	
	#extra params
	opt_params=collections.OrderedDict([('epmra',-9.999), ('epmdec',-9.999),
	('erv',-9.999), ('edist',-9.999),
	('kin_in',None)])

	#params that will be calculated
	calc_params=collections.OrderedDict([('x_gc',-9.999),('y_gc',-9.999),('z_gc',-9.999),('R_gal',-9.999),
	('u',-9.999),('du',-9.999),('v',-9.999),('dv',-9.999),('w',-9.999),('dw',-9.999),
	('vx_gc',-9.999),('vy_gc',-9.999),('vz_gc',-9.999),
	('R',-9.999),('phi',-9.999),('vRg',-9.999),('vTg',-9.999),
	('N_orbits',-9.999),('bound',False),
	('Energy',-9.999),('Energy_unc',-9.999),('L_z',-9.999),('L_z_unc',-9.999),('L_p',-9.999),('L_p_unc',-9.999),('I_3',-9.999),('I_3_unc',-9.999),('Z_max',-9.999),('Z_max_unc',-9.999),
	('ecc',-9.999),('ecc_unc',-9.999),('r_apo',-9.999),('r_apo_unc',-9.999),('r_peri',-9.999),('r_peri_unc',-9.999),('R_apo_P',-9.999),('R_apo_P_unc',-9.999),('R_peri_P',-9.999),('R_peri_P_unc',-9.999)])

reset()
loop()

f.close()
o.close()