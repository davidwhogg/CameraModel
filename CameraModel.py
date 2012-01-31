if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':12})
    rc('text', usetex=True)
    import pylab as plt
import numpy as np

def intensity(l,a,x,y):
    '''
    Given a set of K coherently emitting monochromatic sources
    emitting radiation at wavelength l from K three-dimensional
    positions x[:] with K complex amplitudes a[:], compute N
    intensities I[:] at N three-dimensional positions y[:], assuming
    no reflections or absorptions anywhere.
    '''
    dist = np.sqrt(np.sum((x[np.newaxis,:,:] - y[:,np.newaxis,:])**2, axis=2))
    field = np.sum((a[np.newaxis,:] / dist) * np.exp(2. * np.pi * 1j * dist / l), axis=1) 
    return np.real(field)**2 + np.imag(field)**2

def main():
    lam = 1.2 # wavelength in microns
    xres = 0.5 * lam # resolution of antenna grid
    D = 20. * lam # camera aperture in microns
    f = 4. * D # focal length in microns

    # set up perfect parabolic lens or mirror
    x0, x1 = np.mgrid[0:2 * D / lam + 1,0:2 * D / lam + 1]
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

    # distort lens
    order = 3
    for o in range(order + 1):
        for yo in range(o + 1):
            xo = o - yo
            x[:,2] += np.random.normal() * (x[:,0] / D)**xo * (x[:,1] / D)**yo

    # set up antenna phases
    a = np.ones(K).astype('complex') # all identical

    # set up focal plane and put it in focus
    yres = 0.2 * f * lam / D
    Ny = 65
    y0, y1 = np.mgrid[0:Ny, 0:Ny]
    y0 = (y0.flatten() - np.median(y0)) * yres
    y1 = (y1.flatten() - np.median(y1)) * yres
    N = len(y0)
    y = np.zeros((N, 3))
    y[:,0] = y0
    y[:,1] = y1
    y[:,2] = f

    # compute and plot intensity
    I = intensity(lam, a, x, y)
    plt.clf()
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
    plt.imshow(I.reshape((Ny, Ny)), vmin=0., vmax=1.)
    plt.axis('equal')
    plt.savefig('cm.png')
    return None

if __name__ == '__main__':
    main()
