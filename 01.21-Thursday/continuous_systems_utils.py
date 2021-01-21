import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import Normalize
import matplotlib.animation as animation
from IPython.display import clear_output

# calculate RGB colors from a matplotlib colormap.
def get_colors(theta, colormap='hsv'):
    '''
    Map an array of periodic values (0 to 2*pi) to RGBA color channels of a matplotlib colormap.
    theta - the array of values.
    colormap - the name of the colormap to use.
    '''
    theta = np.mod(theta, 2*np.pi)                            # shift all theta to lie on [0,2*pi], since values are cyclic.
    norm = Normalize(vmin=0.,vmax=2*np.pi)                    # set the color range from [0,2*pi]
    cmap = plt.get_cmap(colormap)                             # get the colormap object
    colors = cmap(norm(theta))                                # evaluate the colormap to get color channels (RGBA)
    colors[..., -1] = 0.8                                     # change the opacity
    return colors

# Plot a single snapshot of the Kuramoto model.
def plot_kuramoto(theta, colormap='hsv'):
    '''
    Plot a lattice of oscillators (phases) as a color image.
    theta - the array of oscillator phases.
    colormap - an optional specification for a colormap to use.
    '''
    plt.figure(figsize=(6,6))                                 # initialize figure
    colors = get_colors(theta, colormap)                      # get colors associated with each theta value
    plt.imshow(colors, extent=[-1, 1, -1, 1])                 # plot the colors
    plt.xlim(-1,1); plt.ylim(-1,1)
    plt.axis('off')
    plt.show()
    clear_output(wait=True)

# Animate Kuramoto model.
def animate_kuramoto(thetas):
    frames = len(thetas)
    
    # initialize a figure.
    fig, ax =plt.subplots(1,1,figsize=(6,6))
    colors = get_colors(thetas[0,:,:])
    im = ax.imshow(colors, extent=[-1, 1, -1, 1])
    ax.set_xlim(-1,1); ax.set_ylim(-1,1)
    ax.axis('off'); fig.tight_layout()

    def animate(i):
        colors = get_colors(thetas[i,:,:])
        im.set_data(colors)
        return im,

    ani = animation.FuncAnimation(fig, animate, frames=frames, interval=50, blit=True)
    plt.close(fig)
    return ani

# define a custom plotting function which plots swarmalator positions x and y, and colors agents by their phase.
def plot_swarm(x, y, theta, colormap='hsv'):
    '''
    Plot swarmalators in 2D with phase represented by color.
    x - the x-coordinates of the swarmalators.
    y - the y-coordinates of the swarmalators.
    theta - the phases of the swarmalators.
    colormap - an optional specification for a colormap to use.
    '''
    plt.figure(figsize=(7,7))
    plt.scatter(x, y, color=get_colors(theta, colormap))
    plt.tick_params(labelsize=14)
    plt.xlim(-2,2); plt.ylim(-2,2)
    plt.grid(False)
    plt.show()
    clear_output(wait=True)

def animate_swarm(q):
    frames = len(q)
    
    # initialize a figure.
    fig, ax = plt.subplots(1,1,figsize=(6,6))
    scat = ax.scatter(q[0,:,0], q[0,:,1], color=get_colors(q[0,:,2]))
    ax.tick_params(labelsize=14)
    ax.set_xlim(-2,2); ax.set_ylim(-2,2)
    ax.grid(False); fig.tight_layout()
    
    def animate(i):
        scat.set_offsets(q[i,:,:2])
        colors = get_colors(q[i,:,2])
        scat._facecolors = colors
        scat._edgecolors = colors
        return scat,

    ani = animation.FuncAnimation(fig, animate, frames=frames, interval=50, blit=True)
    plt.close(fig)
    return ani
