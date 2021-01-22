import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm

# Helper function to read in a data text file, sort data by state,
# and write the sorted numeric data to a binary file.
def load_data(feature):
    state_names = []
    feature_values = []
    f = open("data/"+feature+".txt","r")
    contents = f.readlines()
    f.close()
    header = contents[0]
    for i,line in enumerate(contents[1:51]):
        a = line.split()
        if len(a) == 2:
            state_names += [a[0]]
        elif len(a) == 3:
            state_names += [a[0]+" "+a[1]]
        else:
            print("Encountered invalid line (%d): %s"%(i+1,a))
            break
        feature_values += [float(a[-1])]
    sort_indices = [i for (s,i) in sorted((s,i) for (i,s) in enumerate(state_names))]
    feature_values = np.array(feature_values)[sort_indices]
    feature_values.astype('float').tofile("data/"+feature+".bin")
    print("Saved numeric data to 'data/%s.bin'"%feature)
    return header

# Helper function to make a truncated colormap.
def make_colormap(colormap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=colormap.name, a=minval, b=maxval),
        colormap(np.linspace(minval, maxval, n)))
    return new_cmap

# Add a colorbar.
def add_cbar(data, fig, ax, colormap):
    cnorm = colors.Normalize(vmin=np.min(data), vmax=np.max(data))
    sm = cm.ScalarMappable(cmap=colormap, norm=cnorm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.05, shrink=0.7)
    return cbar

# Helper function to make a plot of a map with associated colorbar.
def plot_map(data, header, colormap, white_edges, feature_transform=None):
    # Shift and rescale data.
    norm_data = data - np.min(data)
    norm_data /= np.max(norm_data)
    
    # Apply colormap to create a color image.
    im = (255.*colormap(norm_data)[:,:,:3]).astype(np.uint8)
    im[data==0]=255.
    m,n,p = im.shape

    # detect edges.
    dx, dy = np.gradient(norm_data)
    dr = np.sqrt(dx*dx+dy*dy)
    if white_edges:
        im[dr>0]=255.
    else:
        im[dr>0]=0.
        
    if feature_transform:
        f = np.fromfile("out/"+feature_transform+"_fields.bin")
        u = f[:m*n].reshape((m,n))
        X = f[m*n:2*m*n].reshape((m,n))
        Y = f[2*m*n:].reshape((m,n))
        X = np.stack([X,Y], axis=-1)
        cart = np.zeros((m,n,p),dtype=np.uint8)
        for i in range(m):
            for j in range(n):
                i2=int(X[i,j,0]+0.5)
                j2=int(X[i,j,1]+0.5)
                if i2<0: i2=0
                elif i2>m-1: i2=m-1
                if j2<0: j2=0
                elif j2>n-1: j2=n-1
                cart[i,j,:]=im[i2,j2,:]
        im = cart
        
    # Visualize the map with edges.
    fig, ax = plt.subplots(1,1,figsize=(24,16))
    ax.imshow(im)
    cbar = add_cbar(data, fig, ax, colormap) # add a corresponding colorbar
    ax.axis('off')
    ax.set_title(header, size=24)
    cbar.ax.tick_params(labelsize=20)
    fig.tight_layout()
    plt.show()
    return fig

# Helper function to make a map given a data size and feature.
def populate_map(size, feature):
    if size == 'large':
        m=816; n=1216
    elif size == 'med':
        m=408; n=608
    else:
        m=204; n=304
    feature_values = np.fromfile("data/"+feature+".bin")
    state_labels = np.fromfile("maps/"+size+".bin", dtype=np.uint8).reshape((m,n))
    
    # Set values in each state.
    feature_map = np.zeros((m,n))
    for i in range(m):
        for j in range(n):
            if state_labels[i,j]>0:
                feature_map[i,j] = feature_values[state_labels[i,j]-1]
    return feature_map # return a map with feature values populated in each state.

