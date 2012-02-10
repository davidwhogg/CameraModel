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

class OpticalSurface(ndarray):
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
        self.x = np.zeros((K, 3))
        self.x[:,0] = x0
        self.x[:,1] = x1
        return None

    def distort_parabolically(self, f):
        '''
        Distort the optical surface from a flat plane into a mirror of
        focal length f.  This routine assumes that the surface is
        centered and aligned along the z (x[:,2]) axis.
        
        If random, then do the same but not on a regular square grid.
        '''
        self.x[:,2] += 0.25 * (self.x[:,0]**2 + self.x[:,1]**2) / f # all on a parabola
        return None

    def distort_randomly(self, amp=1., order=10):
        '''
        Shift the surface in the z direction by random amounts
        according to a hard-coded and insane set of rules.
        '''
        for o in range(order + 1):
            for yo in range(o + 1):
                xo = o - yo
                self.x[:,2] += amp * np.random.normal() * (self.x[:,0] / D)**xo * (self.x[:,1] / self.D)**yo
        return None

    def field(self, lam, a, x):
        '''
        Given a set of K monochromatic sources emitting radiation at
        wavelength lam from K three-dimensional positions x[:] with K
        complex amplitudes a[:], compute N complex amplitudes b[:] at
        N three-dimensional positions self.x[:], assuming no
        reflections or absorptions anywhere (that is, vacuum).
        '''
        dist = np.sqrt(np.sum((x[:,np.newaxis,:] - self.x[np.newaxis,:,:])**2, axis=2))
        return np.sum((a[np.newaxis,:] / dist) * np.exp(2. * np.pi * 1j * dist / lam), axis=0)

class CameraStage():
    '''
    A stage in a CameraModel made up of a transmitter OpticalSurface
    and a receiver OpticalSurface.  See class Camera() for more
    context and some detail.

    Input at initialization: A transmitter surface (of class
    OpticalSurface) and a receiver surface (same).
    '''
    def __init__(self, transmitter, receiver):
        self.transmitter = transmitter
        self.reciever = receiver
        return None

    def receive(self, lam, amplitude):
        '''
        Input: the wavelength lambda at which we are working, and the
        specification of the complex amplitude map at the transmitter
        surface

        Output: the complex amplitude map at the receiver surface
        '''
        return receiver.field(lam, amplitude, self.transmitter)

class Camera():
    '''
    The Camera Model after which this project is named is made up of a
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

    def take_one_image(self, lam, amplitude):
        '''
        Input: the wavelength lam at which the image is to be taken
        and the specification of the complex amplitude map at the
        transmitter surface of the zeroth stage

        Output: the intensity at the receiver surface of the final
        stage
        '''
        amp = amplitude
        for stage in self.stages:
            amp = stage.receive(lam, amp)
        return np.real(amp)**2 + np.imag(amp)**2

    def plot_lens_and_intensity(x,y,I):
        '''
        Basic plotting, with everything (stupidly) hard-coded.

        THIS DOES NOT CURRENTLY WORK
        '''
        Ny = np.round(np.sqrt(len(y))).astype(int)
        plt.subplot(2,2,1)
        plt.plot(x[:,2], x[:,1], 'k.', alpha=0.25)
        plt.axis('equal')
        plt.subplot(2,2,2)
        plt.plot(x[:,0], x[:,1], 'k.', alpha=0.25)
        plt.axis('equal')
        plt.subplot(2,2,3)
        plt.plot(x[:,0], x[:,2], 'k.', alpha=0.25)
        plt.axis('equal')
        plt.subplot(2,2,4)
        plt.gray()
        logI = np.log(I.reshape((Ny, Ny)))
        vmax = np.max(logI)
        plt.imshow(logI, vmin=vmax-8., vmax=vmax)
        plt.axis('equal')
        return None

def main():
    '''
    THIS DOES NOT CURRENTLY WORK
    '''
    lam = 1.2 # wavelength in microns
    D = 20. * lam # camera aperture in microns
    f = 4. * D # focal length in microns
    xp = parabolic_mirror(D, 0.5 * f, 0.5 * lam, random=True)
    x = distort_points(xp)
    K = len(x)
    a = np.ones(K).astype('complex') # all identical phases to start
    Ny = 64
    y = focal_plane_array(Ny, 0.2 * f * lam / D)
    y[:,2] = f
    I = intensity(lam, a, x, y)
    plt.clf()
    plot_lens_and_intensity(x,y,I)
    plt.savefig('cm.png')
    return None

if __name__ == '__main__':
    main()
