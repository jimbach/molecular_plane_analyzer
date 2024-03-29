#!/usr/bin/python

##############################
# molecular_plane_analyzer
# 2020, Jim Bachmann
#
##############################

import os,sys,numpy as np

#######   28.07.2013 - analysis_per_frame   ########
####### script to analyse .gro trajectories ########

if 1==1:																									#check the input
	SCRIPT=sys.argv[0].split('/')
	NARGS=len(sys.argv)-1
	THELP=['-h','-help']

	if NARGS<2 or (NARGS>0 and sys.argv[1] in THELP):
		print""
		print"parameterexpected:",SCRIPT[len(SCRIPT)-1],
		print"[input.xyz] [control_file.fcf]"
		print""
		sys.exit(0)

	if (sys.argv[1][-4:] != ".xyz") or (sys.argv[2][-4:] != ".fcf"):
		print""
		print" wrong input format"
		print" expected: [input.xyz] [control_file.fcf]"
		print""
		sys.exit(0)

if 1==1:												#define output file
	STDOUT=sys.stdout

if 1==1:												#initialize variables
													#these variables are meant to be read only
	number_of_atoms = 0								#contains the number of atoms			
	inputDATA = []										#use adressing: inputDATA[atomnumber][0=atomtype 1=x 2=y 3=z]
	ipoint =  0
	reset = 0
#	count_atoms = 0
	timestep = 0.0
	atomdistkey, moldistkey, comdistkey, angleskey, dipbi_vec_key, C5_P3HT_ang_key, com_key, masses_data, read_done = 0,0, 0, 0, 0, 0, 0, 0, 0
	vector_correlation_key, c5_planaritaet_key, p3_sidechain_key = 0, 0, 0
	ATOMD_list, MOLD_list, COMD_list, ANGLES_list, masses_list = [],[],[],[],[]								#for reading data
	no_timestep = 0
	no_end_timestep = 0
	readcenters=0
	readreference=0
	necessary=0
	if NARGS<4: no_timestep = 1
	if NARGS<5: no_end_timestep = 1

coords=open(sys.argv[1])										#for headlines
format_control_file=open(sys.argv[2])	

num_of_vectors=1
centerlist=[]
referencelist=[]

######## In this part the control file is read
while 1:
	cf_line = format_control_file.readline()

	if not cf_line:	break										#stop when reaching end of file
	if len(cf_line)<=1:	continue								#skip empty lines
	cf_keyword = cf_line.split()[0]																					
	if cf_keyword[:1]=="#": continue								#skip comments

	if cf_keyword=="NORMALVEC":		normal_vec_key = 1


	if normal_vec_key == 1:
		if readcenters==1:
			if len(cf_line.split())>1:
				centerlist.append(map(lambda y:int(y), cf_line.split()))
				print cf_line
			elif cf_keyword == 'END_DEFINECENTER':
				if readcenters==1:
					readcenters = 2
			else:
				print"Error reading centers! After 'DEFINECENTER' the atom index line for centers should follow"
				print""
				sys.exit(0)
				
		if readreference==1:
			if len(cf_line.split())>1:
				referencelist.append(map(lambda y:int(y), cf_line.split()))
				print cf_line
			elif cf_keyword == 'END_DEFINEREFERENCE':
				if readreference==1:
					readreference = 2
			else:
				print"Error reading centers! After 'DEFINEREFERENCE' the atom index line for centers should follow"
				print""
				sys.exit(0)
				
		if cf_keyword == 'DEFINECENTER':
			if readcenters==2:
				num_of_vectors=num_of_vectors+1
			readcenters = 1
			necessary=necessary+1
			print 'reading center:'

		if cf_keyword == 'DEFINEREFERENCE':
			readreference = 1
			necessary=necessary+1
			print 'reading reference:'
			
		if cf_keyword == 'ATOMSPERMOLECULE':
			atoms_per_molecule=int(cf_line.split()[1])
			necessary=necessary+1
			print str(atoms_per_molecule)+' atoms per molecule'
			
		if cf_keyword == 'BOXVECTORS':
			boxX,boxY,boxZ=float(cf_line.split()[1]),float(cf_line.split()[2]),float(cf_line.split()[3])
			necessary=necessary+1
			print str(boxX)+' X boxvector, '+str(boxY)+' Y boxvector, '+str(boxZ)+' Z boxvector '
			
		if cf_keyword == 'NORMALVECTORSTEP':
			normal_vec_step=int(cf_line.split()[1])
			necessary=necessary+1
			print 'calculating normal vector each: '+str(normal_vec_step)+' step/s'
			
		if cf_keyword == 'END_NORMALVEC':
			normal_vec_key = 2

#check if all necessary keywords were supplied
if necessary < 5:
	print "found " + str(necessary) + ' out of 5 necessary inputs'
	print"Error! Please specify all input parameters! DEFINECENTER, DEFINEREFERENCE, NORMALVECTORSTEP, BOXVECTORS and ATOMSPERMOLECULE need to be set."
	print""
	sys.exit(0)
	
######## Define your scripts here

######## Calculate normal-Vectors of molecule surface
def normal_vectors(ATOM_list,  frame, number_of_atoms, number_of_molecules, number_of_vectors, centerlist, referencelist, boxX, boxY, boxZ):
	"""calculates the normal Vectors of a molecule-surface """	
	NORMAL_VECTORS_out=open("normalvectors.xyz",'a')
	
	NORMAL_VECTORS_check=open("check_normalvectors.xyz",'a')
	
	ANGLEDIST_VECTORS_out=open("angle_dist_vectors.xyz",'a')
	
	
	Vorzugsrichtung_out=open("sum_over_all_normalvectors.dat",'a')
	
	
	
	NORMAL_VECTORS_check.write(str(  (number_of_vectors*(2 + int(len(referencelist))) )*number_of_molecules )+'\n')
	NORMAL_VECTORS_check.write(str(frame)+'\n')
	
	NORMAL_VECTORS_out.write(str(   (number_of_vectors*2)*number_of_molecules )+'\n')
	NORMAL_VECTORS_out.write(str(frame)+'\n')

	ANGLEDIST_VECTORS_out.write(str(   (number_of_vectors*2)*number_of_molecules )+'\n')
	ANGLEDIST_VECTORS_out.write(str(frame)+'\n')
	
	vorzugsvektor = np.asarray([0,0,0])
	reference=0
	
	
	for l in range(0,number_of_molecules):
		for vec_count in range(0,number_of_vectors):
			r = []
			#get the 'middle' of the DPI molecule
			
		
			#get geometric center
			center=np.asarray([0,0,0])
			#append all vectors defining the center
			for i in range(len(centerlist[vec_count])):
				# x y z
				r.append( np.asarray( [ ATOM_list[centerlist[vec_count][i] + reference][1],ATOM_list[centerlist[vec_count][i]  + reference][2],ATOM_list[centerlist[vec_count][i]  + reference][3] ] ) ) 
				center = center + r[i]
			center = center / len(centerlist[vec_count])
		
			#calculate all reference vectors
			r=[]
			for i in range(len(referencelist[vec_count])):
				r.append(np.asarray( [ ATOM_list[referencelist[vec_count][i]  + reference][1],ATOM_list[referencelist[vec_count][i]  + reference][2],ATOM_list[referencelist[vec_count][i]  + reference][3] ] ) )
		
			#get the cross products
			normal_vec= np.asarray([0,0,0])
			for i in range(len(r)):
				if (i<(len(r)-1)):
					
					#minimum image convenction
					if ( abs(r[i][0]-center[0]) < abs(r[i][0]-center[0]+boxX)  ) and ( abs(r[i][0]-center[0]) < abs(r[i][0]-center[0]-boxX)  ):
						dx = r[i][0]-center[0]
					elif (abs(r[i][0]-center[0]+boxX) < abs(r[i][0]-center[0]-boxX) ):
						dx = r[i][0]-center[0]+boxX
					else:
						dx = r[i][0]-center[0]-boxX
					
					if ( abs(r[i][1]-center[1]) < abs(r[i][1]-center[1]+boxY)  ) and ( abs(r[i][1]-center[1]) < abs(r[i][1]-center[1]-boxY)  ):
						dy = r[i][1]-center[1]
					elif (abs(r[i][1]-center[1]+boxY) < abs(r[i][1]-center[1]-boxY) ):
						dy = r[i][1]-center[1]+boxY
					else:
						dy = r[i][1]-center[1]-boxY
						
					if ( abs(r[i][2]-center[2]) < abs(r[i][2]-center[2]+boxZ)  ) and ( abs(r[i][2]-center[2]) < abs(r[i][2]-center[2]-boxZ)  ):
						dz = r[i][2]-center[2]
					elif (abs(r[i][2]-center[2]+boxZ) < abs(r[i][2]-center[2]-boxZ) ):
						dz = r[i][2]-center[2]+boxZ
					else:
						dz = r[i][2]-center[2]-boxZ
					
					
					#minimum image convenction
					if ( abs(r[i+1][0]-center[0]) < abs(r[i+1][0]-center[0]+boxX)  ) and ( abs(r[i+1][0]-center[0]) < abs(r[i+1][0]-center[0]-boxX)  ):
						dx2 = r[i+1][0]-center[0]
					elif (abs(r[i+1][0]-center[0]+boxX) < abs(r[i+1][0]-center[0]-boxX) ):
						dx2 = r[i+1][0]-center[0]+boxX
					else:
						dx2 = r[i+1][0]-center[0]-boxX
					
					if ( abs(r[i+1][1]-center[1]) < abs(r[i+1][1]-center[1]+boxY)  ) and ( abs(r[i+1][1]-center[1]) < abs(r[i+1][1]-center[1]-boxY)  ):
						dy2 = r[i+1][1]-center[1]
					elif (abs(r[i+1][1]-center[1]+boxY) < abs(r[i+1][1]-center[1]-boxY) ):
						dy2 = r[i+1][1]-center[1]+boxY
					else:
						dy2 = r[i+1][1]-center[1]-boxY
						
					if ( abs(r[i+1][2]-center[2]) < abs(r[i+1][2]-center[2]+boxZ)  ) and ( abs(r[i+1][2]-center[2]) < abs(r[i+1][2]-center[2]-boxZ)  ):
						dz2 = r[i+1][2]-center[2]
					elif (abs(r[i+1][2]-center[2]+boxZ) < abs(r[i+1][2]-center[2]-boxZ) ):
						dz2 = r[i+1][2]-center[2]+boxZ
					else:
						dz2 = r[i+1][2]-center[2]-boxZ
					
					vec1=np.asarray([dx,dy,dz])
					vec2=np.asarray([dx2,dy2,dz2])
					#dy = r[i][1]-center[1]
					#dz = r[i][2]-center[2]
					#normal_vec = normal_vec + np.cross(r[i]-center, r[i+1]-center)
					normal_vec = normal_vec + np.cross(vec1, vec2)/np.linalg.norm( np.cross(vec1, vec2) )
				elif ( i==len(r) ):
					
					#minimum image convenction
					dx = min( abs(r[i][0]-center[0]), abs(r[i][0]-center[0]+boxX), abs(r[i][0]-center[0]-boxX) )
					dy = min( abs(r[i][1]-center[1]), abs(r[i][1]-center[1]+boxX), abs(r[i][1]-center[1]-boxX) )
					dz = min( abs(r[i][2]-center[2]), abs(r[i][2]-center[2]+boxX), abs(r[i][2]-center[2]-boxX) )
					
					dx2 = min( abs(r[0][0]-center[0]), abs(r[0][0]-center[0]+boxX), abs(r[0][0]-center[0]-boxX) )
					dy2 = min( abs(r[0][1]-center[1]), abs(r[0][1]-center[1]+boxX), abs(r[0][1]-center[1]-boxX) )
					dz2 = min( abs(r[0][2]-center[2]), abs(r[0][2]-center[2]+boxX), abs(r[0][2]-center[2]-boxX) )
					
					vec1=np.asarray([dx,dy,dz])
					vec2=np.asarray([dx2,dy2,dz2])
					#dy = r[i][1]-center[1]
					#dz = r[i][2]-center[2]
					#normal_vec = normal_vec + np.cross(r[i]-center, r[i+1]-center)
					normal_vec = normal_vec + np.cross(vec1, vec2)/np.linalg.norm(np.cross(vec1, vec2))
					#normal_vec = normal_vec + np.cross(r[i]-center, r[0]-center) 
		
			#normalize
			normal_vec = normal_vec/np.linalg.norm(normal_vec)
		
			#	#direction vectors
			#	dir_vec1 = np.asarray([r[6][0],r[6][1],r[6][2]])
			#	dir_vec2 = np.asarray([r[7][0],r[7][1],r[7][2]])
			#	dir_vec3 = np.asarray([r[8][0],r[8][1],r[8][2]])
			#	dir_vec4 = np.asarray([r[9][0],r[9][1],r[9][2]])
			#	dir_vec = (dir_vec1 + dir_vec2 + dir_vec3 + dir_vec4)/4.0
			#	dir_vec = mittelpunkt_dpi - dir_vec			

			#change direction according to bending
			#	if normal_vec[0]!=0 and dir_vec[0]!=0:
			#		if normal_vec[0]/abs(normal_vec[0])!=dir_vec[0]/abs(dir_vec[0]):
			#			normal_vec[0] = -normal_vec[0]
			#	if normal_vec[1]!=0 and dir_vec[1]!=0:
		#			if normal_vec[1]/abs(normal_vec[1])!=dir_vec[1]/abs(dir_vec[1]):
		#				normal_vec[1] = -normal_vec[1]
		#		if normal_vec[2]!=0 and dir_vec[2]!=0:
		#			if normal_vec[2]/abs(normal_vec[2])!=dir_vec[2]/abs(dir_vec[2]):
		#				normal_vec[2] = -normal_vec[2]

			
			#check for a vorzugsvektor
			vorzugsvektor = vorzugsvektor + normal_vec
			NORMAL_VECTORS_out.write('C' +'	'+ str((center)).replace('[', '').replace(']','')+ '\n')
			NORMAL_VECTORS_out.write('O' +'	'+ str((center+normal_vec)).replace('[', '').replace(']','')+ '\n')
			ANGLEDIST_VECTORS_out.write('C' +'	'+ str((center)).replace('[', '').replace(']','')+ '\n')
			ANGLEDIST_VECTORS_out.write('H' +'	'+ str((normal_vec)).replace('[', '').replace(']','')+ '\n')
			NORMAL_VECTORS_check.write('C' +'	'+ str((center)).replace('[', '').replace(']','')+ '\n')
			NORMAL_VECTORS_check.write('O' +'	'+ str((center+normal_vec)).replace('[', '').replace(']','')+ '\n')
			NORMAL_VECTORS_check.write('H' +'	'+ str((normal_vec)).replace('[', '').replace(']','')+ '\n')
			for i in range(len(referencelist[vec_count])):
				NORMAL_VECTORS_check.write('N' +'	'+ str(ATOM_list[referencelist[vec_count][i] + reference][1])+'	'+ str(ATOM_list[referencelist[vec_count][i] + reference][2])+'	'+ str(ATOM_list[referencelist[vec_count][i] + reference][3])+ '\n')

		#always skip 1 entire molecule forward
		reference=reference + number_of_atoms/number_of_molecules
		
	#calculate a vorzugsvektor as sum over DIPBI normal vectors divided by number of DPBI vectors
	vorzugsvektor = vorzugsvektor/(number_of_vectors*frames*number_of_molecules)
	Vorzugsrichtung_out.write(str(frame) + '   ' + str((vorzugsvektor)).replace('[', '').replace(']','')+'\n')
	NORMAL_VECTORS_out.close()
	NORMAL_VECTORS_check.close()
	return
######## End of Calculate normal-Vectors of molecule surface
	



######## In this part the data is read and the scripts are applied
# variables to check if we are beginning to read a new frame with string "header"
begin=1
header=''
frames=0

while 1:
	line1 = coords.readline()		#to avoid smaug killing

	if ipoint>1 and (ipoint == number_of_atoms+2):
		frames=frames+1
		if normal_vec_key == 2:
			if frames%normal_vec_step == 0:
				normal_vectors(inputDATA,frames,number_of_atoms,number_of_molecules,num_of_vectors,centerlist,referencelist, boxX, boxY, boxZ)
		##### end of analyze part. please dont modify the following code
		inputDATA = []
		ipoint = 0

	if not line1:
		break

	line=line1.split()			#manage keywords
	
	#check if we are beginning to read a new frame
	if (begin==1):     
		number_of_atoms=int(line1)
		number_of_molecules=number_of_atoms/atoms_per_molecule
		print 'found ', number_of_molecules, ' molecules'
		begin=0
		
	if ipoint>1 and ipoint<(number_of_atoms+2):						#+2 to skip top lines in .xyz file
		#building .xyz frames as a list [atomname, x, y, z]
		inputDATA.append([str(line[0]),float(line[1]),float(line[2]),float(line[3])])
	
	ipoint = ipoint+1
coords.close()
