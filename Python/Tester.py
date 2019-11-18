"""Script to test APO Exposure timer calc"""

# execfile('APOinputclasses.py')
from Python.APOinputclasses import Sky, Target, Instrument, Observation

import matplotlib.pyplot as plt

sky = Sky(lunar_phase=0, seeing=1.5)
star = Target(6.5, 'VEGAMAG', [5000, 6000], temp=6000)
inst1 = Instrument('Arces')

# inst2 = Instrument('Arctic')

ob1 = Observation(star, sky, inst1)
# ob2 = Observation(star, sky, inst1)

sn1 = ob1.SNfromTime(1)
# sn2 = ob2.SNfromTime(1E11)

t1 = ob1.TimefromSN(300)

plt.plot(t1[0], t1[1])
# plt.plot(times[0], times[1])

# plt.plot(sn1[0], sn1[1])
# plt.plot(times[0], times[1])

3
print(t1)
