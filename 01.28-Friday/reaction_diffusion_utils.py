import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns

sns.set(font_scale=1.5, style='white')
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['mathtext.default'] = 'regular'

def plot_pattern(u, colormap='twilight'):
    '''Plot the concentration as a color image.
    u - grid of concentrations of one species.
    colormap - an optional specification for a colormap to use.
    '''
    fig, ax = plt.subplots(1,1,figsize=(6,6))
    im = ax.imshow(u, cmap=plt.get_cmap(colormap),
                   vmin = 0.3, vmax = 1,
                   interpolation='bicubic',
                   extent=[-1, 1, -1, 1], alpha=0.8)
    ax.tick_params(labelsize=14)
    ax.set_xlim(-1,1); ax.set_ylim(-1,1)
    ax.axis('off')
    ax.grid(False)
    
    # return handles to figure, axis, and image
    return fig, ax, im

def animate_pattern(out, colormap='twilight'):
    '''
    Animate the concentration over time.
    out - 3d array of frames.
    colormap - an optional specification for a colormap to use.
    '''
    frames = out.shape[-1]
    fig, ax, im = plot_pattern(out[:,:,0], colormap)

    def animate(i):
        '''Plot updates for animation.'''
        im.set_array(out[:,:,i])
        return im,

    ani = animation.FuncAnimation(fig, animate, frames=frames, interval=50, blit=True)
    plt.close(fig)
    return ani