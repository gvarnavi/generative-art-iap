import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi

# factals returning points
# plot a line segment given one endpoint, an angle, and length:
def draw_line(x, y, theta, r, ax):
    ax.plot([x, x + r*np.cos(theta)],
            [y, y + r*np.sin(theta)],
            color=[0.6,0.4,0.4], lw=1)
    
def fractal_tree(k, x, y, theta, r, pts, ax):
    # If the recursion depth k = 1, terminate and add a point.
    if k == 1:
        pts += [[x + r*np.cos(theta),y + r*np.sin(theta)]]
        if ax is not None:
            draw_line(x,y,theta,r,ax)
        
    # otherwise, draw the current branch, then create 2 new segments +/- phi away from theta,
    # and a factor alpha times the length.
    else:
        pts += [[x + r*np.cos(theta),y + r*np.sin(theta)]]
        if ax is not None:
            draw_line(x,y,theta,r,ax)
        # the coordinates of the branching vertex
        xj = x + r*np.cos(theta)
        yj = y + r*np.sin(theta)
        phi1 = np.pi/3 + (0.1*np.random.random()-0.2)
        phi2 = np.pi/3 + (0.1*np.random.random()-0.2)
        fractal_tree(k-1, xj, yj, theta + phi1, 0.6*r, pts, ax)
        fractal_tree(k-1, xj, yj, theta - phi2, 0.6*r, pts, ax)
        
def dragon_curve(k, x, y, theta, r, pts, ax):
    # If the recursion depth k = 1, terminate and draw a single line.
    if k == 1:
        pts += [[x + r*np.cos(theta),y + r*np.sin(theta)]]
        if ax is not None:
            draw_line(x, y, theta, r, ax)
        
    # otherwise, subdivide the line into 2 new segments:
    else:
        x1 = x
        y1 = y
        dragon_curve(k-1, x1, y1, theta + np.pi/4, r/np.sqrt(2), pts, ax)
        x2 = x + r*np.cos(theta)
        y2 = y + r*np.sin(theta)
        dragon_curve(k-1, x2, y2, theta + 3*np.pi/4, r/np.sqrt(2), pts, ax)
        
def make_fractal(fun, n, theta, ax=None):
    if n < 1:
        raise ValueError('Iteration must be greater than zero.')
    pts = [[0,0]]
    fun(k=n, x=0, y=0, theta=theta, r=1, pts=pts, ax=ax)
    return np.array(pts)

def get_finite_polygons(vor):
    # return only finite regions, distinguished by not having any vertices of negative index.
    vertices = np.array(vor.vertices)
    regions = []
    for i, region in enumerate(vor.point_region):
        vert = vor.regions[region]
        if all(v >= 0 for v in vert): # look for all positive indices.
            regions.append(vert)
    return regions, vertices

def get_finite_polygons_with_base(vor, pts, base=False):
    vertices = np.asarray(vor.vertices.tolist())
    regions = []
    centers = []
    
    # add in first point, which will have incomplete base.
    if base:
        region = vor.point_region[1]
        vert = vor.regions[region]
        vert_coords = vertices[vert]
        ind = np.argsort(vert_coords[:,1])
        b1 = [1.5*vert_coords[ind[0],0],0]
        b2 = [1.5*vert_coords[ind[1],0],0]
        vertices = np.vstack([vertices,b1,b2])
        regions.append([vert[ind[0]],vert[ind[1]],len(vertices)-1,len(vertices)-2])
        centers.append(pts[0])
    for i, region in enumerate(vor.point_region):
        vert = vor.regions[region]
        if all(v >= 0 for v in vert): # finite region
            regions.append(vert)
            centers.append(pts[i])
    return regions, vertices, np.asarray(centers)

def cactus_plot(pts, ax, cmap=plt.cm.viridis):
    vor = Voronoi(pts)
    regions, vertices, centers = get_finite_polygons_with_base(vor, pts, True)
    for i, (point, region) in enumerate(zip(centers,regions)):
        polygon = vertices[region]
        centroid = polygon_centroid(polygon[:,0],polygon[:,1])
        area = polygon_area(polygon[:,0],polygon[:,1])
        if i==0:
            color = cmap((0.2*np.random.random()+0.8))
            ax.fill(*zip(*polygon), alpha=0.2, color=color)
        else:
            if i==1:
                max_area = area
            if (np.linalg.norm(centroid - point))/np.sqrt(area) < 0.4:
                color = cmap((0.4*np.random.random()+0.6*area/max_area))
                ax.fill(*zip(*polygon), alpha=0.2, color=color)
    ax.scatter(pts[:-4,0], pts[:-4,1], color=[0.9,0.4,0.6], zorder=1000, s=10, marker='*')
    ax.axis('equal');ax.axis('off')

def dragon_plot(pts, ax, cmap=plt.cm.cubehelix):
    vor = Voronoi(pts)
    regions, vertices, centers = get_finite_polygons_with_base(vor, pts)
    for i, (point, region) in enumerate(zip(centers,regions)):
        polygon = vertices[region]
        centroid = polygon_centroid(polygon[:,0],polygon[:,1])
        area = polygon_area(polygon[:,0],polygon[:,1])
        if i==0:
            max_area = area
        if (np.linalg.norm(centroid - point))/np.sqrt(area) < 0.5:
            color = cmap((0.4*np.random.random()+0.6*area/max_area))
            ax.fill(*zip(*polygon), alpha=0.2, color=color)
    ax.scatter(pts[:-4,0], pts[:-4,1], color=[0.9,0.4,0.6], zorder=1000, s=10, marker='*')
    ax.axis('equal');ax.axis('off')


def polygon_area(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def polygon_centroid(x,y):
    n = len(x)
    return np.array([np.sum(x)/n, np.sum(y)/n])


