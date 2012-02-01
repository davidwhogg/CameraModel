if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':12})
    rc('text', usetex=True)
import numpy as np
import pylab as plt

def field(lam,a,x,y):
    '''
    Given a set of K monochromatic sources emitting radiation at
    wavelength lam from K three-dimensional positions x[:] with K
    complex amplitudes a[:], compute N complex amplitudes b[:] at N
    three-dimensional positions y[:], assuming no reflections or
    absorptions anywhere (that is, vacuum).
    '''
    dist = np.sqrt(np.sum((x[np.newaxis,:,:] - y[:,np.newaxis,:])**2, axis=2))
    return np.sum((a[np.newaxis,:] / dist) * np.exp(2. * np.pi * 1j * dist / lam), axis=1) 

def intensity(lam,a,x,y):
    '''
    Same as field(l,a,x,y) but returning N intensities I[:] at the N
    three-dimensional positions y[:].
    '''
    b = field(lam,a,x,y)
    return np.real(b)**2 + np.imag(b)**2

def parabolic_mirror(D, f, xres, random=False):
    '''
    Return a set of three-dimensional points x[:] that constitute a
    mirror of radius D and focal length f centered at (0,0,0) and
    centered and aligned along the z (x[:,2]) axis.  Set the (regular,
    square-grid) spacing between points in the x-y plane to xres.

    If random, then do the same but not on a regular square grid.
    '''
    if random:
        Ntmp = np.ceil(D / xres).astype(int)**2
        x0 = np.random.uniform(-0.5 * D, 0.5 * D, Ntmp)
        x1 = np.random.uniform(-0.5 * D, 0.5 * D, Ntmp)
    else:
        x0, x1 = np.mgrid[0:D / xres + 2,0:D / xres + 2]
        x0 = (x0.flatten() - np.median(x0)) * xres
        x1 = (x1.flatten() - np.median(x1)) * xres
    I = (x0**2 + x1**2) < (0.25 * D**2)
    x0 = x0[I]
    x1 = x1[I]
    K = len(x0)
    x = np.zeros((K, 3))
    x[:,0] = x0
    x[:,1] = x1
    x[:,2] = 0.5 * (x0**2 + x1**2) / f # all on a parabola
    return x

def distort_mirror(x):
    '''
    Take N input 3-dimensional points x[:] and shift them in the z
    direction by random amounts according to a hard-coded and insane
    set of rules.
    '''
    amp = 1.0
    order = 10
    dx = np.zeros_like(x)
    D = max(x[:,0]) - min(x[:,0])
    for o in range(order + 1):
        for yo in range(o + 1):
            xo = o - yo
            dx[:,2] += amp * np.random.normal() * (x[:,0] / D)**xo * (x[:,1] / D)**yo
    return x + dx

def focal_plane_array(Ny, yres):
    '''
    Return a set of Ny**2 points y[:] that correspond to a square
    focal-plane detector array centered on and aligned normal to the z
    (y[:,2]) axis.
    '''
    y0, y1 = np.mgrid[0:Ny, 0:Ny]
    y0 = (y0.flatten() - np.median(y0)) * yres
    y1 = (y1.flatten() - np.median(y1)) * yres
    N = len(y0)
    y = np.zeros((N, 3))
    y[:,0] = y0
    y[:,1] = y1
    return y

def plot_lens_and_intensity(x,y,I):
    '''
    Basic plotting, with everything (stupidly) hard-coded.
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
    lam = 1.2 # wavelength in microns
    D = 20. * lam # camera aperture in microns
    f = 4. * D # focal length in microns
    xp = parabolic_mirror(D, f, 0.5 * lam, random=True)
    x = distort_mirror(xp)
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
