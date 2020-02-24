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
		print"[vectors1.xyz] [vectors2.xyz] boxX boxY boxZ cutoff"
		print""
		sys.exit(0)

	if (sys.argv[1][-4:] != ".xyz") or (sys.argv[2][-4:] != ".xyz") or (NARGS>6):
		print""
		print" wrong input format"
		print" expected: [vectors1.xyz] [vectors2.xyz] boxX boxY boxZ cutoff"
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
	vectors1, vectors2 = [],[]						#for reading data
	no_timestep = 0
	no_end_timestep = 0
	readcenters=0
	readreference=0
	necessary=0
	if NARGS<4: no_timestep = 1
	if NARGS<5: no_end_timestep = 1

vectorfile1=open(sys.argv[1])										#for headlines
vectorfile2=open(sys.argv[2])	
boxX=float(sys.argv[3])
boxY=float(sys.argv[4])
boxZ=float(sys.argv[5])
cutoff=float(sys.argv[6])


num_of_vectors=1
centerlist=[]
referencelist=[]

	
######## Define your scripts here

	
######## Calculate angle between nearby DIPBI of C5 rings
def neighbour_direction(fromvectors, tovectors, boxX, boxY, boxZ, cutoff):
	
	#sort the coordinates into a grid for fast neighbour searching
	box_grid_size_x = cutoff
	box_grid_size_y = cutoff
	box_grid_size_z = cutoff
	
	#we cannot encounter distances bigger than this within this theory
	maxdist=np.sqrt(cutoff**3)
	
	box_max_x =int(boxX/box_grid_size_x)
	box_max_y =int(boxY/box_grid_size_y)
	box_max_z =int(boxZ/box_grid_size_z)
	
	print 'grid cells x: '+str(box_max_x)+' y: '+str(box_max_y)+' z: '+str(box_max_z) 
	
	vec_corr=open("distance_to_angle.dat",'a')
	grid_fromvectors_array = [[[[] for col in range(box_max_x)] for row in range(box_max_y)] for depth in range(box_max_z) ]
	grid_tovectors_array = [[[[] for col in range(box_max_x)] for row in range(box_max_y)] for depth in range(box_max_z) ]
	#generate neighbour list for the from vectors
	for counter in range(0,len(fromvectors)/2):	#/2 because normal vectors and positions are saved
		#calculate grid index
		x = int(fromvectors[2*counter][0]/box_grid_size_x)
		y = int(fromvectors[2*counter][1]/box_grid_size_y)
		z = int(fromvectors[2*counter][2]/box_grid_size_z)
		#push it round the box
		if x > box_grid_size_x: x = x - box_max_x
		if y > box_grid_size_y: y = y - box_max_y
		if z > box_grid_size_z: z = z - box_max_z
		if x < 0 : x = x + box_max_x
		if y < 0 : y = y + box_max_y
		if z < 0 : z = z + box_max_z
		#print x,y,z,fromvectors[2*counter]
		#write to array
		#append the indices of the normal vectors found in this grid cell into an array
		grid_fromvectors_array[x][y][z].append(2*counter)
	#generate neighbour list for the to vectors
	for counter in range(0,len(tovectors)/2):	#/2 because normal vectors and positions are saved
		#calculate grid index
		x = int(tovectors[2*counter][0]/box_grid_size_x)
		y = int(tovectors[2*counter][1]/box_grid_size_y)
		z = int(tovectors[2*counter][2]/box_grid_size_z)
		#push it round the box
		if x > box_grid_size_x: x = x - box_max_x
		if y > box_grid_size_y: y = y - box_max_y
		if z > box_grid_size_z: z = z - box_max_z
		if x < 0 : x = x + box_max_x
		if y < 0 : y = y + box_max_y
		if z < 0 : z = z + box_max_z
		#write to array
		grid_tovectors_array[x][y][z].append(2*counter)

	for counter_x in range(0,box_max_x):
		print str(float(counter_x)/box_max_x)+'% complete'
		for counter_y in range(0,box_max_y):
			for counter_z in range(0,box_max_z):
				if len(grid_tovectors_array[counter_x][counter_y][counter_z])>0:
					#sit on P3HT and watch the DIPBI nearby
					for tovector_count in range(0,len(grid_tovectors_array[counter_x][counter_y][counter_z])):
						#in your own grid cell and the 26 neighbouring ones
						for i_x in range(0,3):
							for i_y in range(0,3):
								for i_z in range(0,3):
									#counter including neighbour cells
									i_help_x = counter_x + i_x -1 #(0 1 2) -1 -> (-1 0 1)
									i_help_y = counter_y + i_y -1
									i_help_z = counter_z + i_z -1
									#manage pbc
									if i_help_x==-1: i_help_x = box_max_x #python counting
									if i_help_x==box_max_x: i_help_x = 0
									if i_help_y==-1: i_help_y = box_max_y
									if i_help_y==box_max_y: i_help_y = 0
									if i_help_z==-1: i_help_z = box_max_z
									if i_help_z==box_max_z: i_help_z = 0
									#print str(i_help_x)+'	'+str(i_help_y)+'	'+str(i_help_z)
									for fromvector_count in range(0,len(grid_fromvectors_array[i_help_x][i_help_y][i_help_z])-1):
										ortsvec1 = fromvectors[grid_fromvectors_array[i_help_x][i_help_y][i_help_z][fromvector_count]]
										ortsvec2 = tovectors[grid_tovectors_array[counter_x][counter_y][counter_z][tovector_count]]
										normalvec1 =fromvectors[grid_fromvectors_array[i_help_x][i_help_y][i_help_z][fromvector_count]+1]
										normalvec2 =tovectors[grid_tovectors_array[counter_x][counter_y][counter_z][tovector_count]+1]
										#print ortsvec1,ortsvec2,normalvec1,normalvec2
										#minimum image convention
										dx=min(abs(ortsvec1[0] - ortsvec2[0]),abs(ortsvec1[0] - ortsvec2[0]-boxX),abs(ortsvec1[0] - ortsvec2[0]+boxX))
										dy=min(abs(ortsvec1[1] - ortsvec2[1]),abs(ortsvec1[1] - ortsvec2[1]-boxY),abs(ortsvec1[1] - ortsvec2[1]+boxY))
										dz=min(abs(ortsvec1[2] - ortsvec2[2]),abs(ortsvec1[2] - ortsvec2[2]-boxZ),abs(ortsvec1[2] - ortsvec2[2]+boxZ))
										distance = np.sqrt(dx**2+dy**2+dz**2)
										#omit self computation
										if distance!=0:
											if distance>maxdist:
												print 'Error, too large distance computed, distance > sqrt(cutoff**3) not possible!'
											angle_corr = np.arccos(np.dot(normalvec1,normalvec2))*360/(2*np.pi)
											vec_corr.write(str(grid_fromvectors_array[i_help_x][i_help_y][i_help_z][fromvector_count]) +'	'+str(grid_tovectors_array[counter_x][counter_y][counter_z][tovector_count]) +'	' + str(distance)+'	' + str(angle_corr)+  '\n')
					
								
	vec_corr.close()	
	#control_out.close()			
	return



######## In this part the data is read and the scripts are applied
# variables to check if we are beginning to read a new frame with string "header"
begin=1
header=''
frames=0
centers1=[]
centers2=[]
readcenter=1
readcenter2=1
number_of_vecs_from=0
number_of_vecs_to=0
count2=0

#read vectors from which to calc
while 1:
	line1 = vectorfile1.readline()		#to avoid smaug killing

	#check if we are beginning to read a new frame
	if (begin==1):     
		number_of_atoms=int(line1)
		begin=0
		
	if ipoint>1 and (ipoint == number_of_atoms+2):
		frames=frames+1
		
		#read vectors to which to calc
		readcenters2=1
		while 1:
			line2=vectorfile2.readline()
			
			if not line2:
				break
			
			if count2==0:
				number_of_tovectors=int(line2)
			
			splitline=line2.split()			#manage keywords
			if len(splitline)>1:
				vectors2.append(np.asarray( [ float(splitline[1]), float(splitline[2]), float(splitline[3]) ] ) )
			
			#framewise analysis
			if count2==number_of_tovectors+1: break
			
			count2=count2+1
			
		neighbour_direction(vectors1, vectors2, boxX, boxY, boxZ, cutoff)
		##### end of analyze part. please dont modify the following code
		vectors1,vectors2 = [], []
		ipoint = 0
		count2 = 0
		number_of_vecs_from=0

	if not line1:
		break

	line=line1.split()			#manage keywords
	
		
	if ipoint>1 and ipoint<(number_of_atoms+2):						
		#build a vector list of vector 1
		vectors1.append(np.asarray( [ float(line[1]), float(line[2]), float(line[3]) ] ) )
		
	ipoint = ipoint+1
vectorfile1.close()
vectorfile2.close()
