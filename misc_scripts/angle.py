import numpy as np
from sys import argv

nf = int(argv[1])
for i in range(nf):
    print(f"{i:03}", 2*np.pi/nf * i)
