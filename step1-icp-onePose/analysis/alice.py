import numpy as np
import matplotlib.pyplot as plt
import math

def ReadFile(path, fileName):
	print("Data File Path : {}".format(path))
	print("File Name : {}".format(fileName))

	# read
	f = open(path+fileName)
	lines = f.readlines()

	e_axis	= []
	e_angle = []
	for line in lines:
		line = line.strip().split()
		e_axis.	append(float(line[0]))
		e_angle.append(float(line[1]))

	return e_axis, e_angle 

if __name__ == '__main__':
	print("hello")

	path = './'
	fileName = 'ICP_Errors.txt'
	e_axis, e_angle = ReadFile(path, fileName)

	hist, bin_edges = np.histogram(e_angle)
	print(hist)
	print(bin_edges)

	plt.bar(hist, bin_edges[:-1],   label = 'Angle')

	plt.legend(frameon=True)
	plt.xlabel('Error of Angle',fontdict={'family' : 'Times New Roman', 'size': 12})
	plt.ylabel('Counts',fontdict={'family' : 'Times New Roman', 'size': 12})
	plt.title('Errors of the angle')
	plt.tight_layout()
	plt.savefig('figure-errorOfAngle.png',dpi=300)

	plt.show()


