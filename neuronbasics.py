from neuron import h
import matplotlib.pyplot as plt
soma = h.Section(name = 'soma')
h.topology()

soma.psection()['morphology']['L']

soma.L = 20
soma.diam = 20

dir(soma)

import textwrap
print(textwrap.fill(','.join(dir(h))))

help(soma.connect)

soma.insert('hh')

print("type(soma) = {}".format(type(soma)))
print("type(soma(0.5)) = {}".format(type(soma(0.5))))

mech = soma(0.5).hh
print(dir(mech))

print(mech.gkbar)
print(soma(0.5).hh.gkbar)

iclamp = h.IClamp(soma(0.5))
print([item for item in dir(iclamp) if not item.startswith('__')])

iclamp.delay = 2
iclamp.dur = 0.1
iclamp.amp = 0.9

v = h.Vector().record(soma(0.5)._ref_v)             # Membrane potential Vector
t = h.Vector().record(h._ref_t)                     # Time stamp vector

h.load_file('stdrun.hoc')
h.finitialize(-65)
h.continuerun(40)
plt.figure()
plt.plot(t, v)
plt.xlabel('t (ms)')
plt.ylabel('v (mV)')
plt.show()

try:
    from neuron.units import ms, mV
except ModuleNotFoundError:
    ms = 1
    mV = 1
    
try:
    h.PlotShape(False).plot(plt)
    s = h.PlotShape(True)
except AttributeError:
    s = h.PlotShape(True)