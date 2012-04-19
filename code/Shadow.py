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
    def __init__(self, fratio=0.5, distorted=False):
        self.stages = []
        D = 12. # camera aperture diameter in microns
        res = 0.1 # resolution of grid in microns (1/4 wave-ish)
        transmitter = OpticalSurface(D, res)
        f = fratio * D # distance from screen to detector
        if distorted:
            transmitter.distort_randomly()
        detector = OpticalSurface(0.75 * D, 2. * res, square=True)
        detector.shift(f) # z-axis offset
        detector.x[:,0] += 0.5 * detector.D
        detector.x[:,1] += 0.5 * detector.D
        self.stages.append(CameraStage(transmitter, detector))
        return None

if __name__ == '__main__':
    np.random.seed(42)
    lam = 0.55 # microns
    for distorted in (False, ): # only do undistorted
        plt.figure(figsize=(10,10))
        plt.clf()
        k = 0
        for dlam in (None, lam / 1.5): # microns
            sb = ShadowBox(distorted=distorted)
            image = sb.take_one_image(lam, dlam=dlam, nlam=128)
            Ny = np.round(np.sqrt(len(image))).astype(int)
            I = image.reshape((Ny, Ny))
            vmax = np.mean(I)
            k += 1
            plt.subplot(2,2,k)
            plt.gray()
            plt.imshow(I, vmin=0., vmax=vmax, interpolation='nearest')
            plt.axis('equal')
            k += 1
            plt.subplot(2,2,k)
            plt.plot(I[0,:])
            plt.axhline(vmax)
            if distorted:
                suffix = 'distorted'
            else:
                suffix = 'undistorted'
            plt.savefig('shadowbox-%3.2f-%s.png' % (lam, suffix))
