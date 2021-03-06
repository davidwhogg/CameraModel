'''
This file is part of CameraModel.
Copyright 2012 David W. Hogg (NYU) <http://cosmo.nyu.edu/hogg/>.

This implements the basic camera model as a set of stages, each of
which has optical surfaces.  It is slow, because it doesn't yet use
the FFT tricks Greengard suggests.  Issues include:

* Needs to implement FFT tricks.
* Needs to permit incoherent superposition of many wavelengths!
'''

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':12})
    rc('text', usetex=True)
import numpy as np
import numpy.polynomial.chebyshev as cheby
import pylab as plt

class OpticalSurface():
    '''
    ## class `OpticalSurface`:

    A surface of a lens or a mirror or a wavefront.  At
    initialization, the surface is flat, circular, centered on the
    origin, and normal to the z (`x[:,2]`) axis.

    # initialization input:

    * `D`: The diameter of the surface.
    * `xres`: The spacing of the control-point grid.
    * `square`: If `True`, make the surface square.
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
        '''
        ## `copy()`:

        Duh.
        '''
        other = OpticalSurface(0., 1.)
        other.D = self.D
        other.xres = self.xres
        other.square = self.square
        other.x = self.x.copy()
        return other

    def shift(self, deltaz):
        '''
        ## `shift()`:

        Move the optical surface in the z (`x[:,2]`) direction by
        distance `deltaz`.
        '''
        self.x[:,2] += deltaz
        return None

    def distort_parabolically(self, f):
        '''
        ## `distort_parabolically()`:

        Distort the optical surface from a flat plane into a mirror of
        focal length `f`.  This routine assumes that the surface is
        centered and aligned along the z (`x[:,2]`) axis.  Positive
        and negative `f` correspond to the same focal length but
        different orientations; RTFSC.
        '''
        self.x[:,2] += 0.25 * (self.x[:,0]**2 + self.x[:,1]**2) / f
        return None

    def distort_randomly(self, amp=0.06, startorder=10, endorder=15):
        '''
        ## `distort_randomly()`:

        Shift the surface in the z (`x[:,2]`) direction by random
        amounts according to a hard-coded and insane set of rules.

        # bugs:

        There is almost certainly a chebyshev-based one-liner that can
        be done here, but I can't figure out the chebyshev coefficient
        formatting for two-d chebyshev functions.  And the numpy
        documentation I see online is inconsistent with what I get
        locally on my computer; not sure what's up.
        '''
        xx = 2. * self.x[:,0] / self.D # -1 < xx < 1
        yy = 2. * self.x[:,1] / self.D # -1 < yy < 1
        for o in range(startorder, endorder + 1):
            for yo in range(o + 1):
                xo = o - yo
                cx = cheby.Chebyshev(np.append(np.zeros(xo), 1.))
                cy = cheby.Chebyshev(np.append(np.zeros(yo), 1.))
                self.x[:,2] += (amp / o) * np.random.normal() * cx(xx) * cy(yy)
        return None

    def drill_hole(self, r):
        '''
        ## `drill_hole()`:

        Remove surface points within a circular region of radius `r`
        at the center.
        '''
        I = (self.x[:,0]**2 + self.x[:,1]**2) > (r**2)
        self.x = self.x[I]
        return None

    def field(self, lam, a, transmitter):
        '''
        ## `field()`:

        Given a `transmitter`, which is an `OpticalSurface`, emitting
        radiation at wavelength `lam` with complex amplitude map
        `a[:]`, compute and return complex amplitudes map at the
        receiver surface (self) assuming no reflections or absorptions
        anywhere (that is, vacuum).

        This function does the heavy lifting, which is the convolution
        by the Green function.  It is way stupid because it does a
        huge quadratic operation.
        '''
        dist = np.sqrt(np.sum((transmitter.x[:,np.newaxis,:] - self.x[np.newaxis,:,:])**2, axis=2))
        return np.sum((a[:,np.newaxis] / dist) * np.exp(2. * np.pi * 1j * dist / lam), axis=0)

class CameraStage():
    '''
    ## class `CameraStage`:

    A stage in a `CameraModel` made up of a transmitter
    `OpticalSurface` and a receiver `OpticalSurface`.  See classes
    `OpticalSurface()` and `Camera()` for more context and some
    detail.

    # initialization input:

    * `transmitter`: An `OpticalSurface`.
    * `receiver`: An `OpticalSurface`.
    '''
    def __init__(self, transmitter, receiver):
        self.transmitter = transmitter.copy()
        self.receiver = receiver.copy()
        return None

    def receive(self, lam, amplitude):
        '''
        ## `receive()`:

        Perform the Green function convolution of the complex
        `amplitude` on the `transmitter` to get the complex amplitude
        on the `receiver`.

        # input:

        * `lam`: the wavelength lambda at which we are working.
        * `amplitude`: the complex amplitude map at the transmitter
          surface.

        # output:

        The complex amplitude map at the receiver surface.
        '''
        return self.receiver.field(lam, amplitude, self.transmitter)

class Camera():
    '''
    ## `Camera`:

    The camera model after which this project is named is made up of a
    set of `CameraStage`s.  In each `CameraStage`, there is a
    transmitter `OpticalSurface` that is broadcasting wavefronts and a
    receiver `OpticalSurface` that is detecting them (via `receive()`
    functions, which are convolutions).  In a `Camera` with multiple
    `CameraStages`, the properties of the transmitter `OpticalSurface`
    are set by the receiver `OpticalSurface` from the previous
    `CameraStage`, and the information obtained by the receiver
    `OpticalSurface` is used to set the transmitter `OpticalSurface`
    of the next `CameraStage`.  This might be an idiotic methodology.

    # initialization input:

    * `stages`: A list (iterable) of `CameraStage` elements.
    '''
    def __init__(self, stages):
        self.stages = stages
        return None

    def take_one_image(self, lam, dlam=None, nlam=16):
        '''
        ## `take_one_image()`:

        Transmit an input complex amplitude map through the whole
        `Camera`.  In the future, this needs to permit something other
        than a normal plane wave!  And something other than a "white"
        spectrum (in whatever units we are working).

        # input:

        * `lam`: The wavelength lam at which the image is to be taken.
        * `dlam`: If not `None`, the range of wavelengths around the
          central wavelength.
        * `nlam`: The number of wavelengths to use (not used if `dlam`
          is `None`)

        # output:

        The intensity (not amplitude!) at the receiver surface of the
        final `CameraStage`.
        '''
        if dlam is None:
            lams = [lam]
        else:
            ddlam = dlam / float(nlam)
            lams = (lam - 0.5 * dlam + 0.5 * ddlam) + ddlam * np.arange(nlam)
        receivedamp = 0.
        for lam in lams:
            phase = np.random.uniform(0., 2. * np.pi)
            print lam, phase
            amp = np.ones(len(self.stages[0].transmitter.x)) * np.exp(1j * phase)
            for stage in self.stages:
                amp = stage.receive(lam, amp)
            receivedamp += amp
        return np.real(receivedamp)**2 + np.imag(receivedamp)**2

class Coronograph(Camera):
    '''
    ## class `Coronograph`:

    A particular kind of `Camera` that has the properties of a
    coronograph.  The zeroth `CameraStage` has a parabolic incoming
    wavefront transmitter, converging on a focal plane receiver with a
    hole drilled in it.  The first `CameraStage` makes that focal
    plane the transmitter back to an identical parabolic reflector
    receiver.  The third has that reflector transmit back to the final
    focal plane.

    There are many dumb choices here, including magic numbers chosen
    without analysis.  RTFSC.

    # initialization input:

    * `fratio`: F ratio of the primary optics (default 4).
    * `distorted`: If `False`, don't apply random distortions.
    '''
    def __init__(self, fratio=4, distorted=True):
        self.stages = []
        D = 30. # camera aperture diameter in microns
        res = 0.5 # resolution in microns
        transmitter = OpticalSurface(D, res)
        f = fratio * D # focal length
        transmitter.distort_parabolically(0.5 * f)
        if distorted:
            transmitter.distort_randomly()
        receiver = OpticalSurface(f, 0.5 * fratio * res)
        holeradius = 4. * fratio # microns
        receiver.drill_hole(holeradius)
        receiver.shift(f)
        self.stages.append(CameraStage(transmitter, receiver))
        self.stages.append(CameraStage(receiver, transmitter))
        detector = OpticalSurface(f, 2. * res, square=True)
        detector.shift(f)
        self.stages.append(CameraStage(transmitter, detector))
        return None

def main():
    '''
    For command-line operation, to show stuff to Leslie Greengard
    (NYU).
    '''
    np.random.seed(42)
    lam = 1.5
    for distorted in (False, True):
        cg = Coronograph(distorted=distorted)
        image = cg.take_one_image(lam)
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
        plt.savefig('coronograph-%3.2f-%s.png' % (lam, suffix))
    return None

if __name__ == '__main__':
    main()
