import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import itertools
from mayavi import mlab

#######################################################################
df = pd.read_csv("../data/optimization_results/cleaned_optimization_DALY_response_coverage_10iter.csv")

bin_widths = [1,3, 2,2,2,3,3]*8


age_group_xpos = df.age_group_xpos.unique()
age_group_xtick = df.age_group_xtick.unique()
age_groups = df.age_group.unique()
scenario_groups = df.scenario.unique()

age_group_dict = {key:value for (key,value) in zip(age_group_xtick, age_groups)}

print age_group_dict
 
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
dy = 0.05*np.ones_like(zpos1)
dz1 = np.ones(x.shape)
dz2 = np.ones(x.shape)



for num1 in xrange(x.shape[0]):
    for num2 in xrange(x.shape[1]):
        xnum = x[num1][num2]
        ynum = y[num1][num2]
        dz1[num1][num2] = df[(df['age_group'] == xnum) & (df['scenario'] == ynum) & (df["type"] == "low_efficacy")]["percentage_vaccinated"]
	dz2[num1][num2] = df[(df['age_group'] == xnum) & (df['scenario'] == ynum) & (df["type"] == "high_efficacy")]["percentage_vaccinated"]



zpos2 = dz1
my_distinct_colors = ["#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe" ,"#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000","#aaffc3","#808000"]

color_list = list(itertools.chain.from_iterable(itertools.repeat(x, 7) for x in my_distinct_colors))



from numpy import pi, sin, cos, mgrid
dphi, dtheta = pi/250.0, pi/250.0
[phi,theta] = mgrid[0:pi+dphi*1.5:dphi,0:2*pi+dtheta*1.5:dtheta]
m0 = 4; m1 = 3; m2 = 2; m3 = 3; m4 = 6; m5 = 2; m6 = 6; m7 = 4;
r = sin(m0*phi)**m1 + cos(m2*phi)**m3 + sin(m4*theta)**m5 + cos(m6*theta)**m7
x = r*sin(phi)*cos(theta)
y = r*cos(phi)
z = r*sin(phi)*sin(theta)

# View it.
from mayavi import mlab
mlab.options.offscreen = True
mlab.test_contour3d()
mlab.savefig('example.png')


"""
mlab.init_notebook()
s = np.abs(np.random.random((3, 3)))
x, y = np.indices(s.shape)
bar1 = mlab.barchart(s)
mlab.show()
#bar1 = mlab.barchart(x, y, dz1)
#bar2 = mlab.barchart(xpos, ypos, dz2)

"""