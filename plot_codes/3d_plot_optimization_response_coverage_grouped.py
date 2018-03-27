import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import itertools
#####################################################################
#visualization tweaks from https://stackoverflow.com/questions/18602660/matplotlib-bar3d-clipping-problems
def sph2cart(r, theta, phi):
    '''spherical to Cartesian transformation.'''
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

def sphview(ax):
    '''returns the camera position for 3D axes in spherical coordinates'''
    r = np.square(np.max([ax.get_xlim(), ax.get_ylim()], 1)).sum()
    theta, phi = np.radians((90-ax.elev, ax.azim))
    return r, theta, phi
#
# end of apodemus's code

def getDistances(view,xpos,ypos,dz1):
    distances  = []
    a = np.array((xpos, ypos, dz1))
    for i in range(len(xpos)):
	#print i, a[0,i], a[1,i], a[2,i]
        distance = (a[0, i] - view[0])**2 + (a[1, i] - view[1])**2 + (a[2, i] - view[2])**2
	distances.append(np.sqrt(distance))
    return distances

########################################################################
def plot_3d_figure(df, ax1, title):
    
    bin_widths = [1,3, 2,2,2,3,3]*8
    
    ax1.view_init(elev=41., azim=-18)
    age_group_xpos = df.age_group_xpos.unique()
    age_group_xtick = df.age_group_xtick.unique()
    age_groups = df.age_group.unique()
    scenario_groups = df.scenario.unique()
    
    age_group_dict = {key:value for (key,value) in zip(age_group_xtick, age_groups)}
    
    xpos,ypos = np.meshgrid(age_group_xpos,scenario_groups)
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    x,y = np.meshgrid(age_groups,scenario_groups)
    #x1,y1 = np.meshgrid(age_group_xpos,scenario_groups)
    lx = x.shape[0]
    ly = x.shape[1]
    n = lx*ly
    zpos1 = np.zeros(n)
    dx = np.array([(num-0.5) for num in bin_widths])
    #dy = dx.copy()
    dy = 0.01*np.ones_like(zpos1)
    dz1 = np.ones(x.shape)
    dz2 = np.ones(x.shape)
    
    
    for num1 in xrange(x.shape[0]):
	for num2 in xrange(x.shape[1]):
	    xnum = x[num1][num2]
	    ynum = y[num1][num2]
	    dz1[num1][num2] = df[(df['age_group'] == xnum) & (df['scenario'] == ynum) & (df["type"] == "low_efficacy")]["percentage_vaccinated"]
	    dz2[num1][num2] = df[(df['age_group'] == xnum) & (df['scenario'] == ynum) & (df["type"] == "high_efficacy")]["percentage_vaccinated"]
    
     
    dz1 = dz1.flatten()
    dz2 = dz2.flatten()
    zpos2 = dz1
    
    # Calculate the distance of each bar from the camera.
    x2, y2, z2 = sph2cart(*sphview(ax1))
    camera = np.array((x2,y2,0))
    z_order = getDistances(camera,xpos,ypos,dz1)
    max1 = max(z_order)
    opacity = 1
    
    
    for i in range(n):
	pl = ax1.bar3d(xpos[i], ypos[i], zpos1[i], dx[i], dy[i], dz1[i],
		color = "#999999", alpha=opacity, zsort='max')
	pl = ax1.bar3d(xpos[i], ypos[i], zpos2[i], dx[i], dy[i], dz2[i],
		color = "#377eb8", alpha=opacity, zsort='max')
	pl._sort_zpos = max1- z_order[i]
	
	
    ax1.set_title(title, fontsize=22)
    ax1.set_xlabel('Age group', labelpad =10, fontsize=18)
    ax1.set_ylabel('Vaccine distribution', labelpad =10, fontsize=18)
    ax1.set_xlim(1,18)
    ax1.set_xticks(age_group_xtick)
    ax1.set_xticklabels([age_group_dict[num] for num in age_group_xtick], rotation = 0,va='bottom',ha = "right")
    ax1.xaxis.labelpad = 15
    ax1.set_ylim(5,90)
    ax1.set_zlim(-1,120)
    ax1.set_zticks(np.arange(0,110,20))
    ax1.set_zlabel('Coverage (%)', fontsize=18)
    plt.autoscale(enable=True, axis='both', tight=True)

#######################################################################




fig = plt.figure(figsize=(14,10))

dt1 = pd.read_csv("../data/optimization_results/cleaned_optimization_infection_response_coverage_100iter.csv")
s1 = fig.add_subplot(221, projection='3d')
title = "Minimizing incidence"
plot_3d_figure(dt1, s1, title)

dt2 = pd.read_csv("../data/optimization_results/cleaned_optimization_hospitalization_response_coverage_100iter.csv")
print ("hospitalization")
s2 = fig.add_subplot(222, projection='3d')
title = "Minimizing hospitalization"
plot_3d_figure(dt2, s2, title)

dt3 = pd.read_csv("../data/optimization_results/cleaned_optimization_death_response_coverage_100iter.csv")
print ("death")
s3 = fig.add_subplot(223, projection='3d')
title = "Minimizing mortality"
plot_3d_figure(dt3, s3, title)

dt4 = pd.read_csv("../data/optimization_results/cleaned_optimization_DALY_response_coverage_100iter.csv")
s4 = fig.add_subplot(224, projection='3d')
title = "Minimizing DALY"
plot_3d_figure(dt4, s4, title)
#print ("check!"), len(scenario_groups),np.arange(1,8,1), [scenario_groups[num] for num in np.arange(0,8,1)]


plt.savefig("../plots/minimizing_burden_response_coverage.png",  pad_inches=0)
#plt.show()
#ax1.view_init(elev=58., azim=-24)

