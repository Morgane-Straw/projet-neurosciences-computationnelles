####### Inhibitory network of Exponential LIF neurons with random connectivity #######
#### Without a specific structure ####

from brian2 import *
import numpy
from numpy import random
import numpy.linalg
import scipy.io
import pickle


Nbinh = 1000 # number of inhibitory neurons

# Parameters
C = 281*pF
gL = 10.6*nS
EL = -70.6*mV
VT = -50.4*mV 
Vcut= -30*mV
DeltaT = 2*mV
tauw = 40*ms
b = 0*nA  #0.0805*nA
taum = C/gL
Vr = EL
a = 3*gL

#For connections and synapses
taui = 10.*ms # inhibitory time constant
wi = 0.02*nS # inhibitory synaptic weight
Ei = -75*mV # Synaptic reversal potential
pCon = 0.2 # connection probability


#---------------
###Clock
Tsim = 1200*ms # simulation duration
defaultclock.dt = 0.02*ms

#---------------

### Model
eqsInh="""
    dv/dt = (gL*(EL-v) + gL*DeltaT*exp((v-VT)/DeltaT) + I - w + gi*(Ei-v))/C : volt
    dw/dt = (a*(v-EL)-w)/tauw : amp
    dgi/dt = -gi*(1./taui) : siemens
    I : amp
"""

### Create the neuron/group of neurons
Ninh=NeuronGroup(Nbinh, model=eqsInh, threshold = 'v>Vcut', reset="v=Vr;w+=b", method='euler')

### Create the connections
Ci = Synapses(Ninh, Ninh, on_pre='gi+=wi')
Ci.connect(p=pCon)

#########################
print ("connections done")
# print ("waiting for plots")
#########################

### Initialization
#Constant initialisation
Ninh.v=-70*mV+(randn(len(Ninh))*2-1)*mV#70.6*mV
Ninh.gi = 3*nS
Ninh.I = 0.82*nA


#### Monitors #####
Minh = SpikeMonitor(Ninh)
Nbspike = Minh.count / Tsim

# A population rate monitor
Pinh = PopulationRateMonitor(Ninh)
    
# To record a few traces
trace=StateMonitor(Ninh,'v',record=[20, 55, 56]) # records voltage for neurons 20, 55, and 56
###################


### Run the simulation
run(Tsim)

print( "The number of spikes in the inh population is ", size(Minh.t))

### Plots
figure(1)
plot(Minh.t/ms, Minh.i, '.k')
xlabel('Time (ms)')
ylabel('Neuron index')
title('Raster plot Adex 1000 neurons')
savefig('RasterPlotInh.png')

figure(2, figsize=(15, 11))
plot(trace.t/ms,trace[20].v/mV, 'g')
xlabel('Time (ms)')
ylabel('Vm neuron 20 (mV)')
    
figure(3, figsize=(15, 11))
plot(trace.t/ms,trace[55].v/mV, 'b')
xlabel('Time (ms)')
ylabel('Vm neuron 55 (mV)')

figure(4, figsize=(15, 11))
plot(trace.t/ms,trace[56].v/mV, 'k')
xlabel('Time (ms)')
ylabel('Vm neuron 56 (mV)')
savefig('VmTraces.png')

# figure(5, figsize=(15, 11))
# plot(Pinh.t/ms,Pinh.smooth_rate(50*ms)) #smoothed curve of firing rates

### Record spike times into a file
output = open('SpikeTimes_inh.pkl', 'wb') #in a pickle file
SpikeTimes= array(Minh.t)
pickle.dump(SpikeTimes,output,-1)
output.close()

mystring1 = "%s%d%s" % ("SpikesTimes_inh_a",a,".txt") #in a text file
FMI2 = open(mystring1,'w')
FMI2.flush()
savetxt(FMI2,Minh.t)
FMI2.close

print ("done")
show()
