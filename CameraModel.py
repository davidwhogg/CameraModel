import numpy as np

def intensity(l,a,x,y):
    """
    Given a set of K coherently emitting monochromatic sources
    emitting radiation at wavelength l from K three-dimensional
    positions x[:] with K complex amplitudes a[:], compute N
    intensities I[:] at N three-dimensional positions y[:], assuming
    no reflections or absorptions anywhere.
    """
    dist = np.sqrt(np.sum((x[:,np.newaxis(),:] - y[np.newaxis(),:,:])^2, axis=2))
    return np.sum((a[:,np.newaxis()] / dist) * np.exp(2. * np.pi * 1j * dist / l), axis=1)^2

def main():
    lam = 1.2
    K = 5
    a = np.random.normal(size=K) + np.random.normal(size=K)j
    print "a", a
    x = np.zeros((K,3))
    x[:,0] = (np.arange(K) - 0.5 * K + 0.5) * 2.
    print "x", x
    N = 20
    y = np.zeros((N,3))
    y[:,0] = (np.arange(N) - 0.5 * N + 0.5) * 0.2 * lam
    y[:,1] = 100. * lam
    print "y", y
    return None

if __name__ == '__main__':
    main()
