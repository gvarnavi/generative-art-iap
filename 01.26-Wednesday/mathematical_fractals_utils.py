import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
from matplotlib.colors import ListedColormap

sns.set(font_scale=1.5, style='white')
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['mathtext.default'] = 'regular'

def cobweb_plot(fmap, points, *args):
    ''' Make a cobweb plot of the function func.
        Inputs:
            fmap - function of the iterative map
            points - points of the trajectory
            args - additional arguments taken by func
        Outputs:
            f, ax - figure and axis objects of resulting plot
    '''
    # set up figure and axes
    f = plt.figure(figsize=(7,7))
    ax = f.add_subplot()
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xlabel('$x_n$')
    ax.set_ylabel('$x_{n+1}$')

    # plot the map func
    x = np.linspace(0,1,100)
    y = fmap(x, *args)
    ax.plot(x, y, lw=1.5, color='k')

    # plot dashed line y = x
    ax.plot(x, x, ls='dashed', color='gray')

    # plot iterated path initialized at x0 over N iterations
    ax.plot(points[:,0], points[:,1], lw=1.5, color='darkslateblue')
    return f, ax


def orbit_diagram(x, p):
    ''' Plot of the stable orbits x over a range of parameter values p.
        Inputs:
            x - array of orbits; the second dimension of x must match the length of p
            p - parameter values
        Outputs:
            f, ax - figure and axis objects of resulting plot
    '''
    # set up figure and axes
    f = plt.figure(figsize=(9,7))
    ax = f.add_subplot()
    ax.set_xlabel('r')
    ax.set_ylabel('x')
    ax.set_xlim([p.min(), p.max()])

    # plot orbit diagram
    ax.plot(p, x.T, ',', color='darkslateblue', alpha=0.5)
    return f, ax

def get_color_image(image, palette='default'):
    ''' Get color image from two-dimensional array.
        Inputs:
          image - two-dimensional image array
          palette - colormap name
        Outputs:
          image - three-dimensional color image array
    '''
    # get colormap
    if palette=='default':
        colors = plt.get_cmap('twilight_shifted').colors
        d = len(colors)
        palette = ListedColormap(np.concatenate([colors[:int(0.5*d)], colors[:int(0.5*d)][::-1]]))
    else: palette = plt.get_cmap(palette)

    # convert image to rgba color values
    image = (image - image.min())/(image.max() - image.min())
    image = palette(image)
    return image

def fractal_plot(image, palette='default'):
    ''' Plot the fractal image.
        Inputs:
          image - two-dimensional image array
          palette - colormap name
        Outputs:
          f, ax - figure and axis objects of resulting plot
    '''
    # set up figure and axes
    f = plt.figure(figsize=(9,9))
    ax = f.add_subplot()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # convert to color image
    image = get_color_image(image, palette=palette)

    # display image
    ax.imshow(np.flipud(image))
    return f, ax

def fractal_movie(images, palette='default'):
    ''' Create an animation of fractal images.
        Inputs:
          images - three-dimensional array of images
          palette - colormap name
        Outputs:
          ani - animation
    '''
    # set up figure and axes
    f = plt.figure(figsize=(9,9))
    ax = f.add_subplot()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # convert to color image
    images = get_color_image(images, palette=palette)

    # initialize animation
    image = np.zeros((images.shape[0], images.shape[1], images.shape[3]), dtype=np.int64)
    img = ax.imshow(image)
    def init():
        img.set_data(images[:,:,0])
        return img,

    # animation function
    def animate(i):
        img.set_data(images[:,:,i+1])
        return img,

    # create animation
    ani = animation.FuncAnimation(f, animate, init_func=init, frames=images.shape[2]-1, interval=100, blit=True)
    plt.close(f)
    return ani
