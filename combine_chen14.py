import math as m
import numpy as np
import os.path

main_dir="chen14"
main_name="chen14_unc"

#how many total files (num jobs ran)
last_file=7862

dir_lim=int(m.ceil(last_file/1000.))

#final file to write to
o=open("{}/{}_kin.csv".format(main_dir,main_name),"a")

first_file=True

start=0
stop=0

for run_dir in range(0,dir_lim):
#each directory as 1000 sub-directories...
	start=stop
	stop=start+1000

	if(run_dir==dir_lim-1):
		stop=last_file	

	for run in range(start,stop):
		file_exists=True
		file_empty=False
		kin_dir="{}/runs_{}/{}_run_{}_kin".format(main_dir,run_dir,main_name,run)
		kin_file="{}/{}_run_{}_kin.csv".format(kin_dir,main_name,run)
		
		if(first_file==True):
			for line in open(kin_file):
				o.write(line)
				first_file=False
		#don't want to copy the header for each subsequent file
		else:
			try:
				f=open(kin_file)
			except:
				file_exists=False
			if(file_exists==True):
				try:
					f.next() #skip header
				except:
					file_empty=True
				if(file_empty==False):
					for line in f:
						o.write(line)
					f.close()
				else:
					print "file empty: {}".format(kin_file)
			else:
				print "file doesn't exist: {}".format(kin_file)

o.close()

