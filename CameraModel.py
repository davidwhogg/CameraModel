'''
This file is part of CameraModel.
Copyright 2012 David W. Hogg (NYU) <http://cosmo.nyu.edu/hogg/>.
'''

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':12})
    rc('text', usetex=True)
import numpy as np
import pylab as plt

class OpticalSurface():
    '''
    A surface of a lens or a mirror or a wavefront.

    Input at initialization: the diameter D of the (circular) surface
    and the resolution xres (really the grid spacing).  At
    initialization, the surface is centered on the origin and normal
    to the z (x[:,2]) axis.
    
    Optional input at initialization: If square, then make the surface
    square not circular.
    '''
    def __init__(self, D, xres, square=False):
        self.D = D
        self.xres = xres
        self.square = square
        x0, x1 = np.mgrid[0:D / xres + 2,0:D / xres + 2]
        x0 = (x0.flatten() - np.median(x0)) * xres
        x1 = (x1.flatten() - np.median(x1)) * xres
        if not square:
            I = (x0**2 + x1**2) < (0.25 * D**2)
            x0 = x0[I]
            x1 = x1[I]
        K = len(x0)
        self.x = np.zeros((K, 3)).astype(float)
        self.x[:,0] = x0
        self.x[:,1] = x1
        return None

    def copy(self):
        other = OpticalSurface(0., 1.)
        other.D = self.D
        other.xres = self.xres
        other.square = self.square
        other.x = self.x.copy()
        return other

    def shift(self, deltaz):
        '''
        Move the optical surface in the z (x[:,2]) direction by
        distance deltaz.
        '''
        self.x[:,2] += deltaz
        return None

    def distort_parabolically(self, f):
        '''
        Distort the optical surface from a flat plane into a mirror of
        focal length f.  This routine assumes that the surface is
        centered and aligned along the z (x[:,2]) axis.  Positive and
        negative f correspond to the same focal length but different
        orientations; RTFSC.
        '''
        self.x[:,2] += 0.25 * (self.x[:,0]**2 + self.x[:,1]**2) / f
        return None

    def distort_randomly(self, amp=0.5, order=10):
        '''
        Shift the surface in the z direction by random amounts
        according to a hard-coded and insane set of rules.
        '''
        for o in range(order + 1):
            for yo in range(o + 1):
                xo = o - yo
                self.x[:,2] += amp * np.random.normal() * (self.x[:,0] / self.D)**xo * (self.x[:,1] / self.D)**yo
        return None

    def drill_hole(self, r):
        '''
        Remove surface points within a circular region of radius r at
        the center.
        '''
        I = (self.x[:,0]**2 + self.x[:,1]**2) > (r**2)
        self.x = self.x[I]
        return None

    def field(self, lam, a, transmitter):
        '''
        Given a transmitter OpticalSurface emitting radiation at
        wavelength lam with complex amplitude map a[:], compute and
        return complex amplitudes map at the receiver surface (self)
        assuming no reflections or absorptions anywhere (that is,
        vacuum).
        '''
        dist = np.sqrt(np.sum((transmitter.x[:,np.newaxis,:] - self.x[np.newaxis,:,:])**2, axis=2))
        return np.sum((a[:,np.newaxis] / dist) * np.exp(2. * np.pi * 1j * dist / lam), axis=0)

class CameraStage():
    '''
    A stage in a CameraModel made up of a transmitter OpticalSurface
    and a receiver OpticalSurface.  See classes OpticalSurface() and
    Camera() for more context and some detail.

    Input at initialization: A transmitter surface (of class
    OpticalSurface) and a receiver surface (same).
    '''
    def __init__(self, transmitter, receiver):
        self.transmitter = transmitter.copy()
        self.receiver = receiver.copy()
        return None

    def receive(self, lam, amplitude):
        '''
        Input: the wavelength lambda at which we are working, and the
        specification of the complex amplitude map at the transmitter
        surface

        Output: the complex amplitude map at the receiver surface
        '''
        return self.receiver.field(lam, amplitude, self.transmitter)

class Camera():
    '''
    The camera model after which this project is named is made up of a
    set of CameraStages.  In each CameraStage, there is a transmitter
    OpticalSurface that is broadcasting wavefronts and a receiver
    OpticalSurface that is detecting them.  In a Camera with multiple
    CameraStages, the properties of the transmitter OpticalSurface are
    set by the receiver OpticalSurface from the previous CameraStage,
    and the information obtained by the receiver OpticalSurface is
    used to set the transmitter OpticalSurface of the next
    CameraStage.  This might be an idiotic methodology.

    Input at initialization: a list of stages (iterable of
    CameraStage elements)
    '''
    def __init__(self, stages):
        self.stages = stages
        return None

    def take_one_image(self, lam, amplitude=None):
        '''
        Input: the wavelength lam at which the image is to be taken
        and the specification of the complex amplitude map at the
        transmitter surface of the zeroth stage

        Output: the intensity at the receiver surface of the final
        stage
        '''
        if amplitude is None:
            amp = np.ones(len(self.stages[0].transmitter.x))
        else:
            amp = amplitude
        for stage in self.stages:
            print amp
            amp = stage.receive(lam, amp)
        print amp
        return np.real(amp)**2 + np.imag(amp)**2

class Coronograph(Camera):
    '''
    A particular kind of Camera that has the properties of a
    coronograph.  The zeroth CameraStage has a parabolic incoming
    wavefront transmitter, converging on a focal plane receiver with a
    hole drilled in it.  The first CameraStage makes that focal plane
    the transmitter back to an identical parabolic reflector receiver.
    The third has that reflector transmit back to the final focal
    plane.
    '''
    def __init__(self, fratio=4, distorted=True):
        self.stages = []
        D = 25. # camera aperture in microns
        res = 0.5 # resolution in microns
        transmitter = OpticalSurface(D, res)
        f = fratio * D # focal length
        transmitter.distort_parabolically(0.5 * f)
        if distorted:
            transmitter.distort_randomly()
        receiver = OpticalSurface(f, res)
        holeradius = 4.5 * fratio # microns
        receiver.drill_hole(holeradius)
        receiver.shift(f)
        self.stages.append(CameraStage(transmitter, receiver))
        self.stages.append(CameraStage(receiver, transmitter))
        detector = OpticalSurface(f, 2. * res, square=True)
        detector.shift(f)
        self.stages.append(CameraStage(transmitter, detector))
        return None

def main():
    cg = Coronograph()
    for lam in np.arange(1.0, 2.0, 0.02):
        image = cg.take_one_image(lam)
        Ny = np.round(np.sqrt(len(image))).astype(int)
        logI = np.log(image.reshape((Ny, Ny)))
        vmax = np.max(logI)
        plt.clf()
        plt.gray()
        plt.imshow(logI, vmin=vmax-8., vmax=vmax)
        plt.axis('equal')
        plt.savefig('cm%3.2f.png' % lam)
    return None

if __name__ == '__main__':
    main()
