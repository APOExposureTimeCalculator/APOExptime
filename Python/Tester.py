"""Script to test APO Exposure timer calc"""
import matplotlib.pyplot as plt
# execfile('APOinputclasses.py')
from Python.APOinputclasses import Sky, Target, Instrument, Observation


sky = Sky(lunar_phase=0, seeing=2)
star = Target(19, 'VEGAMAG', [5000, 6000], temp=6000)
inst1 = Instrument('NICFPS')


ob1 = Observation(star, sky, inst1)


sn1 = ob1.SNfromTime(100)
print(sn1)

t1 = ob1.TimefromSN(50)

#plt.plot(t1[0], t1[1])
# plt.plot(times[0], times[1])

#plt.plot(sn1[0][0], sn1[0][1])
# plt.plot(times[0], times[1])

