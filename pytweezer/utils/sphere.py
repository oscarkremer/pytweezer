import numpy as np

def sphere(r=1.0, disc=21):
    phis = np.linspace(-np.pi, np.pi, disc)
    thetas = np.linspace(np.pi, 0, disc)
    X = np.zeros((disc, disc))
    Y = np.zeros((disc, disc))
    Z = np.zeros((disc, disc))
    for i, theta in enumerate(thetas):
        for j, phi in enumerate(phis):
            X[i][j] = r*np.cos(phi)*np.sin(theta)
            Y[i][j] = r*np.sin(phi)*np.sin(theta)
            Z[i][j] = r*np.cos(theta)
    return X, Y, Z