import LFPy
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib.collections import PolyCollection
#import scipy
#from scipy.signal import butter, lfilter
#import matplotlib.animation as animation
from neuron import h
import os
import math

def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def plotstuff(cell, electrode):
    '''plotting'''
    fig = plt.figure()
    plt.plot(electrode.x, electrode.y,'.',  marker='o', markersize=3, color='r',zorder=0)
#    rotation = {'x' : 0, 'y' : math.pi, 'z' : 0} #-math.pi/9 # Mainen
#    cell.set_rotation(**rotation) 
# Plot neuron morphology
    zips=[]
    for x, y in cell.get_idx_polygons(projection=('x', 'y')):
        #### PATCH Radu do not plot big boxes
#        tmp=list(zip(x,y))
#        if not 250 in tmp[1]:
        if PolyArea(x,y)<10000:
            zips.append(list(zip(x, y)))
        ##### END PATCH
    polycol = PolyCollection(zips,edgecolors='k',facecolors='k')
    ax = fig.add_subplot(111)
    ax.add_collection(polycol)
    ax.axis(ax.axis('equal'))
    plt.xlabel('Distance $\mu$m - (Ox)')
    plt.ylabel('Distance $\mu$m - (0y) ')
    plt.grid(True)
    plt.title(r'$Neuron$ $Morphology$')   
    return fig




# =============================================================================
# ================================= MORPHOLOGY ================================
# =============================================================================


LA="1000"
DA="2"
LD="50"
DD="2"

st=1/1000

filename="./BSR_LA"+LA+"_DA"+DA+"_LD"+LD+"_DD"+DD+"_demo.hoc"
cell_parameters = {    
        'morphology' : filename, 
        'v_init' : -65,
        'passive' : True,
        'passive_parameters' : {'g_pas' : 1./30000, 'e_pas' : -65},
        'cm' : 1.0,
        'Ra' : 150,
        'dt' : st,
        'tstart' : 0.,
        'tstop' : 20.,
        'nsegs_method' : 'lambda_f', # spatial discretization method
        'lambda_f' : 100., # frequency where length constants are computed
                }
cell = LFPy.Cell(**cell_parameters)
                

# =============================================================================
# ================================= stimulation parameters================================
# =============================================================================
                
amp=0.2               
dur=10
delay=1
#                if cell_parameters["passive"]==False:
#                    amp=1.95**(int(ida)/1.95)/30+1.6*int(ila)/10000+int(ild)/5000*int(idd)*1.5
#                else:
#                    amp=1.65**(int(ida)/1.95)/30+0.805*int(ila)/10000+int(ild)/5000*int(idd)*1.1 #-0.4 for Mainen equivalent

stim = {
        'idx' : cell.get_closest_idx(x=0,y=0,z=0),
        'record_current' : True,
        'pptype' : 'IClamp',
        'amp' : amp,
        'dur' : dur, #0.01
        'delay' : delay, #5
        }
                

stimulus=LFPy.StimIntElectrode(cell,**stim)
                
                
# =============================================================================
# ================================= SIMULATION  ===============================
# =============================================================================
                

cell.simulate(rec_imem=True)
#cell.imem[np.isnan(cell.imem)]=0.0


                
# =============================================================================
# ================================= electrodes ================================
# =============================================================================

hstep = np.array(range(-250,1251,125))
vstep = np.array(range(250,49,-50))
N = vstep.shape
Ny = hstep.shape

x_elec = np.tile(hstep, N)
y_elec = np.sort(np.tile(vstep, Ny))[::-1]
z_elec = np.zeros(x_elec.shape)


meshgrid = {
        'sigma':0.33,
        'x': x_elec,
        'y': y_elec, 
        'z': z_elec,
        }

# Electrodes initialization
meshgrid_electrodes = LFPy.RecExtElectrode(cell, **meshgrid)  
meshgrid_electrodes.calc_lfp()
                
                
                
# =============================================================================
# ================================= plot and save  ===============================
# =============================================================================
                
timeind = (cell.tvec > np.argmax(cell.somav)*st-3) & (cell.tvec <= np.argmax(cell.somav)*st+5)


fileout="Vlfpy_BS"+"_LA"+LA+"_DA"+DA+"_LD"+LD+"_DD"+DD+"demo.txt"
np.savetxt(fileout,meshgrid_electrodes.LFP.T[timeind])

fileout2="Vm_BS"+"_LA"+LA+"_DA"+DA+"_LD"+LD+"_DD"+DD+"demo.txt"
np.savetxt(fileout2,cell.somav[timeind])

fileout3="Im_BS"+"_LA"+LA+"_DA"+DA+"_LD"+LD+"_DD"+DD+"demo.txt"
np.savetxt(fileout3,cell.imem[0,timeind])

elpos=(x_elec,y_elec,z_elec)
np.savetxt("elpos_demo",elpos)                


fig = plt.figure()
# ================= STIMULATION PLOT ==========================================
plt.subplot(411)
pl.plot(cell.tvec[timeind], stimulus.i[timeind])
plt.axis('tight')
plt.ylabel(r'$I_s$ (nA)', va='center')
plt.grid(True)
# ================= SOMA PLOT =================================================
plt.subplot(412)
plt.plot(cell.tvec[timeind], cell.somav[timeind],lw=1)
plt.axis('tight')
plt.ylabel(r'$V_m$ (mV)', va='center')
plt.grid(True)
                
## ================= Imem  ============================================
plt.subplot(413)
plt.plot(cell.tvec[timeind], cell.imem[0,timeind],lw=1)#,label='x='+str(elec_parameters["x"])+' - y='+str(elec_parameters["y"]))
plt.axis('tight')
plt.ylabel(r'$I_m$ (nA)', va='center')
plt.grid(True)
                
## ================= ELECTRODE LFP  ============================================
plt.subplot(414)
plt.plot(cell.tvec[timeind], meshgrid_electrodes.LFP.T[timeind],lw=1)#,label='x='+str(elec_parameters["x"])+' - y='+str(elec_parameters["y"]))
plt.ylabel(r'$V_e$ (mV)', va='center')
plt.xlabel('time (ms)')
plt.grid(True)
                
plt.show()


plotstuff(cell,meshgrid_electrodes)
