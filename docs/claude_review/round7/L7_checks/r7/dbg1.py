import numpy as np, warnings
warnings.filterwarnings('ignore')
from kapcurve import curve
from kapcheck import summary, analyse
p = [3, -2, 5/3, 7/2]
c, res, tv = curve(p)
print('resid', res)
print('contacts tv', tv)
s = summary(c)
print(s)
for k, v in analyse(c).items():
    if k[0] != 'hyp': print(k, v)
