from neuron import h
import matplotlib.pyplot as plt 

soma = h.Section(name = 'soma')

soma.psection()
{'point_processes' : {},
 'density_mechs' : {},
 'ions' : {},
 'morphology' : {'L' : 100.0,
    'diam' : [500.0],
    'pts3D' : [],
    'parent' : None,
    'trueparent' : None},
 'nseg' : 1,
 'Ra' : 35.4,
 'cm' : [1.0],
 'regions' : set(),
 'species' : set(),
 'name' : 'soma',
 'hoc_internal_name' : '__nrnsec_0x104lef000',
 'cell' : None}
                    
soma.L = 20
soma.diam = 20      

soma.insert('hh')

iclamp = h.IClamp(soma(0.5))
iclamp.delay = 2
iclamp.dur = 0.1
iclamp.amp = 0.9

v = h.Vector().record(soma(0.5)._ref_v) # Membrane Potential Vector
t = h.Vector().record(h._ref_t) # Time Stamp Vector

h.load_file('stdrun.hoc')

plt.figure()
plt.plot(t,v)
plt.xlabel('t (ms)')
plt.ylabel('soma_v (mV)')
plt.show()

           