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
    lam = 1.2
    K = 11
    deltax = 2.7 * lam
    a = np.ones(K).astype('complex')
    x = np.zeros((K,3))
    x[:,0] = (np.arange(K) - 0.5 * K + 0.5) * deltax
    N = 3000
    y = np.zeros((N,3))
    y[:,0] = (np.arange(N) - 0.5 * N + 0.5) * 0.01 * lam
    y[:,1] = 50. * lam
    I = intensity(lam, a, x, y)
    plt.clf()
    plt.plot(y[:,0], I, 'k-')
    plt.savefig('cm.png')
    return None

if __name__ == '__main__':
    main()
