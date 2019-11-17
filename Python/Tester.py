"""Script to test APO Exposure timer calc"""

# execfile('APOinputclasses.py')
from Python.APOinputclasses import Sky, Target, Instrument, Observation
import matplotlib.pyplot as plt

sky = Sky(lunar_phase=0)
star = Target(20, 'VEGAMAG', [5000, 7000], temp=6000)
inst = Instrument('Arces')
ob = Observation(star, sky, inst)
inst1 = Instrument('Arces')
inst2 = Instrument('Arctic')

ob1 = Observation(star, sky, inst1)
ob2 = Observation(star, sky, inst2)

sn1 = ob1.SNfromTime(100)
sn2 = ob2.SNfromTime(100)


plt.plot(sn1[0], sn1[1])
# plt.plot(times[0], times[1])

plt.plot(sn[0], sn[1])
#plt.plot(times[0], times[1])


print(sn2)
