'''
This file is part of CameraModel.
Copyright 2012 David W. Hogg (NYU) <http://cosmo.nyu.edu/hogg/>.

Look at a diffractive shadow from a circular opening.
'''

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':12})
    rc('text', usetex=True)
import numpy as np
import pylab as plt
from CameraModel import *

class ShadowBox(Camera):
    '''
    ## class `ShadowBox`:

    A particular kind of `Camera` that is a circular aperture in a
    screen!  It is illuminated by a perfect plane wave (for now).

    # initialization input:

    * `fratio`: F ratio of the primary optics (default 4).
    * `distorted`: If `True`, apply random distortions.
    '''
    def __init__(self, fratio=1., distorted=False):
        self.stages = []
        D = 8. # camera aperture diameter in microns
        res = 0.2 # resolution of grid in microns
        transmitter = OpticalSurface(D, res)
        f = fratio * D # distance from screen to detector
        if distorted:
            transmitter.distort_randomly()
        detector = OpticalSurface(2. * D, res, square=True)
        detector.shift(f)
        self.stages.append(CameraStage(transmitter, detector))
        return None

if __name__ == '__main__':
    np.random.seed(42)
    lam = 0.55 # microns
    for dlam in (None, 0.15): # microns
        for distorted in (False, ): # only do undistorted
            sb = ShadowBox(distorted=distorted)
            image = sb.take_one_image(lam, dlam=dlam, nlam=32)
            Ny = np.round(np.sqrt(len(image))).astype(int)
            logI = np.log(image.reshape((Ny, Ny)))
            vmax = np.max(logI)
            plt.clf()
            plt.gray()
            plt.imshow(logI, vmin=vmax-8., vmax=vmax)
            plt.axis('equal')
            if distorted:
                suffix = 'distorted'
            else:
                suffix = 'undistorted'
            dlamlabel = dlam
            if dlamlabel is None:
                dlamlabel = 0
            plt.savefig('shadowbox-%3.2f-%3.2f-%s.png' % (lam, dlamlabel, suffix))
