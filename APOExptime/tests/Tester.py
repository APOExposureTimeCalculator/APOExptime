"""Script to test APO Exposure timer calc"""

# execfile('APOinputclasses.py')
from Python.APOinputclasses import Sky, Target, Instrument, Observation

sky = Sky(lunar_phase=0)
star = Target(20, 'VEGAMAG', [5000, 7000], temp=6000)
inst = Instrument('Arctic')
ob = Observation(star, sky, inst)

times = ob.TimefromSN(100)

sn = ob.SNfromTime(100)

print(times)
print(sn)
