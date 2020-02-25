#!/usr/bin/python

import os,sys


if 1==1:																									#check the input
	SCRIPT=sys.argv[0].split('/')
	NARGS=len(sys.argv)-1
	THELP=['-h','-help']

	if NARGS<>10 or (NARGS>0 and sys.argv[1] in THELP):
		print""
		print"parameterexpected:",SCRIPT[len(SCRIPT)-1],
		print"[input.dat] [output.dat] [x_range] [y_range] [col_1] [col_2] [grid size x (int)] [grid size y (int)] [x shift] [y shift]"
		print""
		sys.exit(0)
inP = open(sys.argv[1])
outP = open(sys.argv[2], 'w')

#initialize grid array - it will contain the number of angles within this grid
grid_array = [[0 for col in range(int(sys.argv[8]))] for row in range(int(sys.argv[7]))]
grid_size_x = float(sys.argv[3])/float(sys.argv[7])
grid_size_y = float(sys.argv[4])/float(sys.argv[8])
i_abs = 0
omitted=0

while 1:
	inPline = inP.readline()
	if not inPline: break
	inPList = inPline.split()
	try:	
		x = int((float(inPList[int(sys.argv[5])-1].replace('[','').replace(']','').replace(',',''))+float(sys.argv[9]) )/grid_size_x)
		y = int((float(inPList[int(sys.argv[6])-1].replace('[','').replace(']','').replace(',',''))+float(sys.argv[10])) /grid_size_y)
		grid_array[x][y] = grid_array[x][y] + 1
		i_abs = i_abs + 1
	except IndexError:
		#print 'error in inputfile, line: ' + str(i_abs + 1) + ' not in range or not a valid number!'
		omitted=omitted+1

print 'found '+str(i_abs)+' data points'
if omitted>0:
	print 'omitted '+str(omitted)+' data points because they were not in range'
	print 'used '+str(i_abs-omitted)+' data points within range'
	print 'some data points were not in range and were omitted. this is normal behaviour. increase xrange/yrange if you dont want to see this message'

#for i_x in range(100):
#	for i_y in range(100):
#		outP.write(str(grid_size_x*i_x ) + '	' +str(grid_size_y*i_y) +'	' + str(float(grid_array[i_x][i_y])/i_abs) + '\n')

for i_x in range(int(sys.argv[7])):
	for i_y in range(int(sys.argv[8])):
		outP.write(str(grid_size_x*i_x ) + '	' +str(grid_size_y*i_y) +'	' + str(float(grid_array[i_x][i_y])/i_abs) + '\n')

inP.close()
outP.close()
