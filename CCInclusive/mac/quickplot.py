import sys
filename = sys.argv[1]
contents = open(filename,'r').read().split()
myx = [x.split(',')[0] for x in contents]
myy = [x.split(',')[1] for x in contents]
import matplotlib.pyplot as plt
import numpy as np
myx = np.array(myx)
myy = np.array(myy)
plt.plot(myx,myy,'ro--')
plt.grid(True)
plt.xlabel('E_mu Value Plugged In [MEV]',fontsize=16)
plt.ylabel('abs("x") Value',fontsize=16)
plt.title('Minimization Visualization: True Energy 629 MeV, Minimized E 1000 MeV',fontsize=16)
plt.show()
