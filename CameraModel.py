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
    D = 100. * lam # camera aperture in microns
    f = 500. * lam # focal length in microns
    K = 200 # number of emitting points

    # set up parabolic lens or mirror
    x = np.zeros((K,3))
    x[:,0] = (2. * np.arange(K) - K + 1.) * D / K
    x[:,1] = 0.5 * x[:,0]**2 / f # all on a parabola

    # set up antenna phases
    a = np.ones(K).astype('complex') # all identical
    # tilt the incoming plane wave by a small angle
    a *= np.exp(1j * np.linspace(0.,5.,K))

    # set up focal plane and put it in focus
    N = 3000
    y = np.zeros((N,3))
    y[:,0] = (np.arange(N) - 0.5 * N + 0.5) * 0.01 * lam
    y[:,1] = f

    # compute and plot intensity
    I = intensity(lam, a, x, y)
    plt.clf()
    plt.subplot(2,1,1)
    for xi, ai in zip(x,a):
        plt.plot([xi[0], xi[0] - np.imag(ai)], [xi[1], xi[1] + np.real(ai)], 'r-')
        plt.plot(xi[0], xi[1], 'k.')
    plt.subplot(2,1,2)
    plt.plot(y[:,0], I, 'k-')
    plt.savefig('cm.png')
    return None

if __name__ == '__main__':
    main()
