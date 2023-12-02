    import numpy as np
    import matplotlib.pyplot as plt
    #Creating a 1D array that has numbers from -5 to 5
    #with 0.1 differences i.e 100 equally spaced points
    points = np.arange(-5, 5, .01)
    xs, ys = np.meshgrid(points, points)
    print(ys)
    v = np.sqrt(xs ** 2 + ys ** 2)
    print(v)
    plt.imshow(v, cmap = plt.cm.gray); plt.colorbar
