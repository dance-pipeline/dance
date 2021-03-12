import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
import plotmol

from bokeh.palettes import Spectral4
from bokeh.plotting import Figure, output_file, save, show
from plotmol.plotmol import default_tooltip_template


make_scatter = True
FILL_ALPHA = 0.05

p = None
with open("metrics.pickle" , "rb") as pkf:
    p = pickle.load(pkf)

en_xtb, en_am1 = p[0]
rm_xtb, rm_am1 = p[1]
sm_xtb, sm_am1 = p[3]
print("enxtb0:", en_xtb[0])
print("enam10:", en_am1[0])

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


if make_scatter:

    # energy of all molecules
    figure = Figure(
            tooltips = default_tooltip_template(),
            title = f"Relative conformer energy error compared to reference QCA energy",
            x_axis_label = "xtb kcal/mol",
            y_axis_label = "am1 kcal/mol",
            plot_width = 800,
            plot_height = 800,
            x_range = [-90,50],
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



    # density energy of all molecules
    figure = Figure(
            tooltips = default_tooltip_template(),
            title = f"Relative conformer energy error compared to reference QCA energy",
            x_axis_label = "xtb kcal/mol",
            y_axis_label = "am1 kcal/mol",
            plot_width = 800,
            plot_height = 800,
            x_range = [-90,50],
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
            line_alpha=0,
            fill_color = "blue",
            fill_alpha = FILL_ALPHA,
            legend_label = f"ddE"
            )
                    
    #show(figure)
    output_file("energy_density.html")
    save(figure)

    # rmsd density for all molecules
    figure = Figure(
            tooltips = default_tooltip_template(),
            title = f"RMSD to reference QCA conformation from a given method",
            x_axis_label = "xtb Å",
            y_axis_label = "am1 Å",
            plot_width = 800,
            plot_height = 800,
            x_range = [0,7],
            y_range = [0,7]
            )

    figure.line(x = [0,7],
                y = [0,7])
                
    plotmol.scatter(figure,
            x = rm_xtb,
            y = rm_am1,
            smiles = sm_rm,
            marker = "o",
            marker_size = 10,
            line_alpha=0,
            fill_color = "blue",
            fill_alpha = FILL_ALPHA,
            legend_label = f"rmsd"
            )
                    
    #show(figure)
    output_file("rmsd_density.html")
    save(figure)



    # rmsd for all molecules
    figure = Figure(
            tooltips = default_tooltip_template(),
            title = f"RMSD to reference QCA conformation from a given method",
            x_axis_label = "xtb Å",
            y_axis_label = "am1 Å",
            plot_width = 800,
            plot_height = 800,
            x_range = [0,7],
            y_range = [0,7]
            )

    figure.line(x = [0,7],
                y = [0,7])
                
    plotmol.scatter(figure,
            x = rm_xtb,
            y = rm_am1,
            smiles = sm_rm,
            marker = "o",
            marker_size = 10,
            marker_color = "black",
            fill_color = "blue",
            legend_label = f"rmsd"
            )
                    
    #show(figure)
    output_file("rmsd.html")
    save(figure)


plt.rcParams["figure.figsize"] = [12,12]

en_vals = np.arange(0,90)
xtb_dens = []
am1_dens = []
for i,val in enumerate(en_vals):
    if i==0: continue
    xtb_dens.append( sum(1 for x in en_xtb if en_vals[i-1] < abs(x) < val) )
    am1_dens.append( sum(1 for a in en_am1 if en_vals[i-1] < abs(a) < val) )
xtb_dens = np.array(xtb_dens)/len(en_xtb)
am1_dens = np.array(am1_dens)/len(en_am1)
x = plt.plot(en_vals[:-1], xtb_dens)
a = plt.plot(en_vals[:-1], am1_dens)
plt.title("Distribution of relative conformer energy error compared to reference QCA energy")
plt.legend((x[0],a[0]),("xtb","am1"))
plt.xlabel("ddE kcal/mol")
plt.ylabel("% of molecules with |ddE| less than x")
plt.savefig("energy_dist.png")
plt.clf()

plt.rcParams["figure.figsize"] = [12,12]
rm_vals = np.arange(0,7,0.1)
xtb_rdens = []
am1_rdens = []
for i,val in enumerate(rm_vals):
    if i==0: continue
    xtb_rdens.append( sum(1 for x in rm_xtb if rm_vals[i-1] <= abs(x) < val) )
    am1_rdens.append( sum(1 for a in rm_am1 if rm_vals[i-1] <= abs(a) < val) )
xtb_rdens = np.array(xtb_rdens)/len(rm_xtb)
am1_rdens = np.array(am1_rdens)/len(rm_am1)
x = plt.plot(rm_vals[:-1], xtb_rdens)
a = plt.plot(rm_vals[:-1], am1_rdens)
plt.title("Distribution of RMSD to reference QCA conformation from a given method")
plt.legend((x[0],a[0]),("xtb","am1"))
plt.xlabel("RMSD Å")
plt.ylabel("molecules with RMSD less than x-value, %")
plt.savefig("rmsd_dist.png")
plt.clf()


