"""Script to test APO Exposure timer calc"""

# execfile('APOinputclasses.py')
from Python.APOinputclasses import Sky, Target, Instrument, Observation
import matplotlib.pyplot as plt

sky = Sky(lunar_phase=0)
star = Target(20, 'VEGAMAG', [5000, 7000], temp=6000)
inst = Instrument('Arctic')
ob = Observation(star, sky, inst)

times = ob.TimefromSN(10)

sn = ob.SNfromTime(1000)


#plt.plot(sn[0], sn[1])
#plt.plot(times[0], times[1])


print(sn)
