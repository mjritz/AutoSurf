import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import stats



f1 = open(sys.argv[1], 'r')
tail = f1.readlines()
f1.close()
#end = f2.readlines()                                                                                    
#f2.close() 
#water = f3.readlines()                                                                                    
#f3.close() 

# initialize some variable to be lists:
tail_data = np.loadtxt(sys.argv[1],skiprows=4)
data1 = np.delete(tail_data, 0, 1)
z_tail_density = np.delete(data1, 1, 1)
tail_z_density_half = z_tail_density[1200:]
tail_z_max = np.argmax(tail_z_density_half[:,1])
tail_half_area = np.trapz(tail_z_density_half[:,1],tail_z_density_half[:,0] )
print 'tail half area'
print tail_half_area
max_tail_row = tail_z_density_half[tail_z_max]

head_data = np.loadtxt(sys.argv[2],skiprows=4)
data2 = np.delete(head_data, 0,1)
z_head_density = np.delete(data2, 1, 1)
head_z_density_half = z_head_density[1200:]
head_z_max = np.argmax(head_z_density_half[:,1])
head_half_area = np.trapz(head_z_density_half[:,1],head_z_density_half[:,0] )
print 'head half area'
print head_half_area
max_head_row = head_z_density_half[head_z_max] 

#print max_tail_row[0]
print 'max_tail_z - max_head_z'
print max_tail_row[0]-max_head_row[0] 

water_data = np.loadtxt(sys.argv[3],skiprows=4)
data3 = np.delete(water_data, 0,1)
z_water_density = np.delete(data3, 1, 1)
water_z_density_half = z_water_density[600:]
gradient = np.gradient(water_z_density_half[:, 1])
max_slope_z = water_z_density_half[np.argmin(gradient),0]

print 'max_tail_z - max_gradient_water'
print max_tail_row[0]-max_slope_z

print tail_half_area, head_half_area,max_tail_row[0]-max_slope_z,max_tail_row[0]-max_head_row[0]

#print max_head_row[0]
#head_z_max = np.argmax(head_z_density_half[:,1])
#head_half_area = np.trapz(head_z_density_half[:,1],head_z_density_half[:,0] )
#print head_half_area
#max_head_row = head_z_density_half[head_z_max]

#print np.amax(data[:3])

# scan the rows of the file stored in lines, and put the values into some
# variables
#for line in tail[4:]:
#  p = line.split()
#  x.append(float(p[1]))
#  all_y.append(float(p[3]))

#xv = np.array(x)
#all_yv = np.array(all_yv)
#xyv = np.append(xv[:],all_yv[:])
#print xyv
#for line in end[4:]:                                                                                    
#for line in end: 
#single_peak = all_y[1200:]  
#max_line=np.argmax(single_peak)
#half_single_peak = single_peak[max_line:]
#print half_single_peak[:]

#print yv2
#yv3 = np.array(y3)

# Not sure...
#plt.pause(.1)

# Plot the data    
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(xv,yv1, 'g--', xv,yv1, 'g-')  
#ax2 = ax.twin()
#ax.plot(xv,yv2, 'r--',xv,yv2, 'r-')
#ax = ax.twin()
#ax.plot(xv1,yv3, 'b--',xv1,yv3, 'b-')
#ax.plot(xv,yv, 'r--', xv,yv, 'ro', xt,yt, 'b--', xt, yt, 'bo')

# Plot the data
# This would work: plt.plot(xv, yv), but cannot set x and y labels
#ax.invert_xaxis()

# Add axis labels
#ax.set_xlim(-75,75)
#ax.set_xlabel(r'$Box\ Length\ (z)$')
#ax.set_ylabel(r'Density')
#ax3.set_axis_off()
#ax2.set_axis_off()
#plt.show()
#fig.savefig(sys.argv[4], dpi=300, figsize=(8, 6))  

#Integrate on positive side of tail plot from (curve max, +infinity)
#area = np.trapz(half_single_peak)
#print area
