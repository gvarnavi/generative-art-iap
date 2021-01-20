import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

sns.set(font_scale=1.5, style='white')
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['mathtext.default'] = 'regular'

def plot_attractor(image, palette='inferno'):
    ''' Plot the attractor image.
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

    # create a border by padding with some zeros
    npad = int(0.1*image.shape[0])
    image = np.pad(image, ((npad,npad), (npad,npad)), mode='constant')

    # display image
    ax.imshow(np.flipud(image), cmap=palette)
    return f, ax

def plot_attractor_quilt(image, palette='inferno', tile=3):
    ''' Plot the attractor image.
        Inputs:
          image - two-dimensional image array
          palette - colormap name
          tile - how many times to tile the image on each side
        Outputs:
          f, ax - figure and axis objects of resulting plot
    '''
    # set up figure and axes
    f = plt.figure(figsize=(9,9))
    ax = f.add_subplot()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # tile the image
    image = np.tile(image, (tile, tile))

    # create a border by padding with some zeros
    npad = int(0.1*image.shape[0])
    image = np.pad(image, ((npad,npad), (npad,npad)), mode='constant')

    # display image
    ax.imshow(np.flipud(image), cmap=palette)
    return f, ax

def get_segments3d(s):
    ''' Convert a trajectory of points in 3d into individual line segments.
        Inputs:
          s - array storing 3d-coordinates of the full trajectory
        Outputs:
          segments - array of line segments
    '''
    points = s.T.reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

def colorline3d(s, palette='CMRmap'):
    ''' Create a color gradient along a trajectory.
        Inputs:
          s - array storing 3d-coordinates of the full trajectory
          palette - colormap name
        Outputs:
          lc - line collection object for plotting
    '''
    N = s.shape[1]
    segments = get_segments3d(s)
    lc = Line3DCollection(segments, array=np.linspace(0,1,N), cmap=palette)
    return lc

def set_limits3d(ax, s):
    ''' Set axis limits of a three-dimensional plot.
        Inputs:
          ax - axis handle
          s - array storing 3d-coordinates of plot
    '''
    ax.set_xlim([s[0,:].min(), s[0,:].max()])
    ax.set_ylim([s[1,:].min(), s[1,:].max()])
    ax.set_zlim([s[2,:].min(), s[2,:].max()])

def plot_trajectory(s, palette='CMRmap'):
    ''' Plot the trajectory s.
        Inputs:
          s - array storing 3d-coordinates of the full trajectory
          palette - colormap name
        Outputs:
          f, ax - figure and axis objects of resulting plot
    '''
    # set up figure and axes
    f = plt.figure(figsize=(9,9))
    ax = f.add_subplot(projection='3d')
    set_limits3d(ax, s)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # plot trajectory
    ax.add_collection3d(colorline3d(s, palette))
    return f, ax

def animate_trajectory(S, palette='rainbow'):
    ''' Create a movie of the trajectory s.
        Inputs:
          S - list of arrays storing 3d-coordinates of the full trajectory
          palette - colormap name
        Outputs:
          ani - animation
    '''
    # set up figure and axes
    f = plt.figure(figsize=(9,9))
    ax = f.add_subplot(projection='3d')
    set_limits3d(ax, np.hstack(S))
    ax.axis('off')
    ax.view_init(30,0)
        
    # choose a different color for each trajectory
    colors = plt.get_cmap(palette)(np.linspace(0, 1,len(S)))

    # initialize lines and points artists for each trajectory
    lines = sum([ax.plot([], [], [], '-', c=c) for c in colors], [])
    markers = sum([ax.plot([], [], [], 'o', c=c) for c in colors], [])
    def init():
        for l, m in zip(lines, markers):
            l.set_data([], [])
            m.set_data([], [])
        return lines + markers

    # animation
    spf = 2             # time steps per frame
    def animate(i):
        i = (spf*i)%S[0].shape[1]
        for l, m, s in zip(lines, markers, S):
            x, y, z = s[:,:i]
            l.set_data(x, y)
            l.set_3d_properties(z)
            m.set_data(x[-1:], y[-1:])
            m.set_3d_properties(z[-1:])

        ax.view_init(30, 0.6*i/spf)
        return lines + markers

    ani = animation.FuncAnimation(f, animate, init_func=init, frames=S[0].shape[1]//spf, interval=30, blit=True)
    plt.close(f)
    return ani
