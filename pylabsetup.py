import matplotlib as mpl
import pylab as plt
from pylab import rc
import numpy as np


def myhisteps(x,y):
    x=np.array(x)
    y=np.array(y)
    if len(x)==len(y):
        dx=np.diff(x)*0.5
        x=np.concatenate([x[:1]-dx[0],x[:-1]+dx,x[-1:]+dx[-1]])
    elif not len(x)==len(y)+1:
        print ('''either pass two arrays of equal length (to make the steps surrounding the x values, len(x)=len(y))
              or the x  change-points, a-la' bayesian block (ending one block and starting the next as
change-points, len(x) len(y)+1 )''')
        return -1,-1
    return np.array(zip(x,x)).flatten()[1:-1],np.array(zip(y,y)).flatten()

rc('axes', linewidth=1.2)
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.size'] = 26.
mpl.rcParams['font.family'] = 'serif'
#mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = [  'Times New Roman', 'Times','Palatino', 'Charter', 'serif']
mpl.rcParams['font.sans-serif'] = ['Helvetica']
mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] = 22.
mpl.rcParams['ytick.labelsize'] = 22.
mpl.rcParams['xtick.major.size']= 15.
mpl.rcParams['xtick.minor.size']= 10.
mpl.rcParams['ytick.major.size']= 15.
mpl.rcParams['ytick.minor.size']= 10.
#mpl.rcParams['figure.autolayout']= True

#fontsize=26
#mpl.rc('axes',  titlesize=fontsize)
#mpl.rc('axes',  labelsize=fontsize)
#mpl.rc('xtick', labelsize=fontsize)
#mpl.rc('ytick', labelsize=fontsize)
#mpl.rc('font', size=fontsize, family='serif', serif='Utopia',
#              style='normal', variant='normal',
#              stretch='normal', weight='normal')
#mpl.rc('font',**{'family':'serif','serif':[ 'Times New Roman', 'Times', 'serif'],
#                 'sans-serif':['Helvetica'], 'size':19, 
#                 'weight':'normal'})
mpl.rc('axes',**{'labelweight':'normal', 'linewidth':1})
mpl.rc('axes',**{'labelweight':'normal', 'linewidth':1})
mpl.rc('ytick',**{'major.pad':5, 'color':'k'})
mpl.rc('xtick',**{'major.pad':5, 'color':'k'})
params = {'legend.fontsize': 24,
          #'legend.linewidth': 1,
          'legend.numpoints':1,
          'legend.handletextpad':1
      }

plt.rcParams.update(params)   
plt.minorticks_on()



