"""Script to test APO Exposure timer calc"""

# execfile('APOinputclasses.py')
from Python.APOinputclasses import Sky, Target, Instrument, Observation

import matplotlib.pyplot as plt

sky = Sky(lunar_phase=0, seeing=1)
star = Target(24.807, 'VEGAMAG', [5000, 6000], temp=6000)
inst1 = Instrument('Arces')

# inst2 = Instrument('Arctic')

ob1 = Observation(star, sky, inst1)
# ob2 = Observation(star, sky, inst1)

sn1 = ob1.SNfromTime(100)
# sn2 = ob2.SNfromTime(1E11)

t1 = ob1.TimefromSN(50)

#plt.plot(t1[0], t1[1])
# plt.plot(times[0], times[1])

#plt.plot(sn1[2][0], sn1[2][1])
# plt.plot(times[0], times[1])


#print(t1)
