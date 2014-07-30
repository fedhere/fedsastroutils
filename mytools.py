import numpy as np


def maskarray(myarray, value):
    mymask=np.asarray([False]*len(myarray))
    mymask[np.where(myarray == value)]=True
    return mymask
