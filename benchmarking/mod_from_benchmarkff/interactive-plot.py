import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import plotmol

from bokeh.palettes import Spectral4
from bokeh.plotting import Figure, output_file, save, show
from plotmol.plotmol import default_tooltip_template
from scipy import stats

import numpy as np

p = None
with open("metrics.pickle" , "rb") as pkf:
    p = pickle.load(pkf)

en_xtb, en_am1 = p[0]
rm_xtb, rm_am1 = p[1]
sm_xtb, sm_am1 = p[3]

en_xtb = np.hstack(en_xtb)
en_am1 = np.hstack(en_am1)
rm_xtb = np.hstack(rm_xtb)
rm_am1 = np.hstack(rm_am1)
sm_xtb = np.hstack(sm_xtb)

print("original size:", len(en_xtb))
dup_remov_dict = {}
for i,k in enumerate(sm_xtb):
    dup_remov_dict[k] = [en_xtb[i], en_am1[i], rm_xtb[i], rm_am1[i]]

nodup = []
for k,v in dup_remov_dict.items():
    nodup.append( [k, *v] )

print(nodup[0])

sm_xtb, en_xtb, en_am1, rm_xtb, rm_am1 = zip(*nodup)
print("without duplicates:", len(en_xtb))


sm_xtb = list(sm_xtb)

for i, sm in enumerate(sm_xtb):
    conf = sm.rfind("-")
    if conf >= 0:
        sm_xtb[i] = sm_xtb[i][:conf]


to_sort_en = list(zip(en_xtb, en_am1, sm_xtb))
sortt_en = sorted(to_sort_en, key = lambda x: x[0]**2 + x[1]**2)


en_xtb, en_am1, sm_en = zip(*sortt_en)


to_sort_rm = list(zip(rm_xtb, rm_am1, sm_xtb))

sortt_rm = sorted(to_sort_rm, key = lambda x: x[0]**2 + x[1]**2)
rm_xtb, rm_am1, sm_rm = zip(*sortt_rm)

print(max(rm_xtb))
print(max(rm_am1))
print(min(rm_xtb))
print(min(rm_am1))

length = len(en_xtb)




# the top 25% of molecules in terms of energy "distance"
figure = Figure(
        tooltips = default_tooltip_template(),
        title = f"qcarchive energy ref",
        x_axis_label = "xtb",
        y_axis_label = "am1",
        plot_width = 800,
        plot_height = 800,
        x_range = [-70,40],
        y_range = [-90,50]
        )

figure.line(x = [-90, 50],
            y = [-90, 50])
plotmol.scatter(figure,
        x = en_xtb[-length//4:],
        y = en_am1[-length//4:],
        smiles = sm_en[-length//4:],
                marker = "o",
                marker_size = 10,
                marker_color = "black",
                fill_color = "blue",
                legend_label = f"ddE, largest 25%"
                )
                
#show(figure)
output_file("energy_top25.html")
save(figure)



# energy of all molecules
figure = Figure(
        tooltips = default_tooltip_template(),
        title = f"qcarchive energy ref",
        x_axis_label = "xtb",
        y_axis_label = "am1",
        plot_width = 800,
        plot_height = 800,
        x_range = [-70,40],
        y_range = [-90,50]
        )

figure.line(x = [-90, 50],
            y = [-90, 50])
plotmol.scatter(figure,
        x = en_xtb,
        y = en_am1,
        smiles = sm_en,
                marker = "o",
                marker_size = 10,
                marker_color = "black",
                fill_color = "blue",
                legend_label = f"ddE"
                )
                
#show(figure)
output_file("energy.html")
save(figure)



# the top 25% of molecules in terms of rmsd "distance"
figure = Figure(
        tooltips = default_tooltip_template(),
        title = f"qcarchive rmsd ref",
        x_axis_label = "xtb",
        y_axis_label = "am1",
        plot_width = 800,
        plot_height = 800,
        x_range = [1,7],
        y_range = [0,5]
        )
figure.line(x = [0,7],
            y = [0,7])
plotmol.scatter(figure,
        x = rm_xtb[-length//4:],
        y = rm_am1[-length//4:],
        smiles = sm_rm[-length//4:],
                marker = "o",
                marker_size = 10,
                marker_color = "black",
                fill_color = "blue",
                legend_label = f"rmsd, top 25% magnitude"
                )
                
#show(figure)
output_file("rmsd_upper25.html")
save(figure)



# rmsd for all molecules
figure = Figure(
        tooltips = default_tooltip_template(),
        title = f"qcarchive rmsd ref",
        x_axis_label = "xtb",
        y_axis_label = "am1",
        plot_width = 800,
        plot_height = 800,
        x_range = [1,7],
        y_range = [0,5]
        )

figure.line(x = [0,7],
            y = [0,7])
            
plotmol.scatter(figure,
        x = rm_xtb[:],
        y = rm_am1[:],
        smiles = sm_rm[:],
                marker = "o",
                marker_size = 10,
                marker_color = "black",
                fill_color = "blue",
                legend_label = f"rmsd"
                )
                
#show(figure)
output_file("rmsd.html")
save(figure)
