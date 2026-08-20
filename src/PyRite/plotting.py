import pathlib
import matplotlib.pyplot as plt
import PyRite
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FixedLocator
import matplotlib.colors as mcolors
from matplotlib import colormaps
import numpy as np

def set_style(style="default", axisRatio = None):
    style_path = (
        pathlib.Path(PyRite.__file__).parent
        / "styles"
        / f"{style}.mplstyle"
    )
    plt.style.use(style_path)

    if axisRatio is not None:
        width, _ = plt.rcParams["figure.figsize"]
        plt.rcParams["figure.figsize"] = (width, width * axisRatio)


def updateAxisColor(c, lw = None):
    
    plt.rcParams.update({
        "axes.edgecolor": c,
        "text.color" : c,
        "axes.labelcolor": c,
        "axes.titlecolor": c,
        "xtick.color": c,
        "ytick.color":c,
    })

    if lw != None:
        plt.rcParams.update({
        "axes.linewidth" : lw
    })
        
def updateFaceColor(c):
     plt.rcParams.update({
        "figure.facecolor": c,
        "axes.facecolor" : c,
        "savefig.facecolor": c,
        "savefig.edgecolor": c,
    })

def remakeTicks(x, y, ax=None, xoff = 0, yoff = 0, xminor = True, yminor = True):
    """
    :param x: distance between ticks in x
    :param y: distance between ticks in y
    :param ax: axis
    :param xoff: offset in x, default 0
    :param yoff: offset in y, default 0
    :return: 
    """
    if ax == None:
        ax = plt.gca()
    ax.xaxis.set_major_locator(MultipleLocator(x, offset = xoff))
    ax.yaxis.set_major_locator(MultipleLocator(y, offset = yoff))

    if xminor == True:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    #elif isinstance(xminor, float) or isinstance(xminor, int):
    if yminor == True:
        ax.yaxis.set_minor_locator(AutoMinorLocator())

def makeSubplots(fig, nx, ny, pads, procentagex = [], procentagey = [],
                 visibility = [], xspace = 0, yspace = 0, colorbar = []):

    """
    :param fig: pyplot figure object
    :param nx: how many subplots in x axis
    :param ny: how many subplots in y axis
    :param pads: how much free space is needed from: left, bottom, right, top (in canvas, so 0-1)
    :param procentagex: what should be the procetage taken by each element on x axis; if empty, all equal (0, 100), lenght of nx
    :param procentagey: what should be the procetage taken by each element on y axis; if empty, all equal (0, 100), lenght of ny
    :param visibility: which axis should be visible; lenght of nx*ny; if empty, all elements will be visible
    :param xspace: space between axis in x (canvas); default = 0
    :param yspace: space between axis in y (canvas); default = 0
    :param colorbar: 3 element list, direction (str, x/y), width and label pad; if empty, no colorbar built
    :return: list of all built axis, starting from left bottom going right, then to the next column etc; colorbar is the last one
    """             
    x0, y0, xk, yk = pads

    dx = 1-x0-xk - xspace*(nx-1)
    dy = 1-y0-yk - yspace*(ny-1)
    
    if len(colorbar) == 0:
        pass
    else:
        if colorbar[0] == 'x':
            dy -= colorbar[1]
            dy -= colorbar[2]
        elif colorbar[0] == 'y':
            dx -= colorbar[1]
            dx -= colorbar[2]
        else:
            print("Wrong colorbar (x or y)")

    axes = []

    y0Here = y0

    naxis = 0
    for y in range(ny):
        if len(procentagey) != ny:
            if len(procentagey) != 0:
                print("Wrong number of variables in procentagey")
            dyHere = dy/ny
        else:
            dyHere = dy*procentagey[y]/100
        x0Here = x0
        for x in range(nx):
            if len(procentagex) != nx:
                if len(procentagex) != 0:
                    print("Wrong number of variables in procentagey")
                dxHere = dx/nx
            else:
                dxHere = dx*procentagex[x]/100

            if len(visibility) != nx*ny:
                if len(visibility) != 0:
                    print("Wrong number of variables in visibility")
                visibleHere = True
            else:
                visibleHere = visibility[naxis]
            
            ax = fig.add_axes([x0Here, y0Here,dxHere, dyHere], visible = visibleHere)
            axes.append(ax)
            x0Here += dxHere + xspace
            naxis += 1
        y0Here += dyHere + yspace

    if len(colorbar) == 0:
        pass
    else:
        if colorbar[0] == 'x':
            x0 = x0
            dx = 1-x0-xk
            y0 = y0Here
            dy = colorbar[1]
        elif colorbar[0] == 'y':
            x0 = x0Here
            dx = colorbar[1]
            y0 = y0
            dy = 1-y0-yk
        axes.append(fig.add_axes([x0, y0, dx, dy]))
    if nx*ny == 1:
        return axes[0]
    return axes

def truncate_colormap(cmap_name, minval=0.0, maxval=1.0, n=256):
    cmap = colormaps.get_cmap(cmap_name)
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        f'trunc({cmap.name},{minval:.2f},{maxval:.2f})',
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap

def showMask(ax, maskFile, maskColor = 'dodgerblue', h = ''):
    """
    :param ax: axis object on which you want to show the mask
    :param maskFile: if str, opens it with astropy.io.fits, otherwise checks if 2D array
    :param maskColor: pyplot color
    :return: None
    """
    if isinstance(maskFile, str):
        maskToShow = fits.open(maskFile)
        if h == '':
            try:
                maskToShow = maskToShow[0].data
            except:
                maskToShow = maskToShow[1].data
        else:
            maskToShow = maskToShow[h].data
            
    else:
        maskToShow = maskFile

    xTick = np.arange(maskToShow.shape[0])
    yTick = np.arange(maskToShow.shape[1])



    YC, XC = np.meshgrid(xTick + (xTick[1] - xTick[0]) / 2+1,
                         yTick + (yTick[1] - yTick[0]) / 2+1)
    if len(np.unique(maskFile)) <2:
        levels = [0.5]
        if maskColor == 'random':
            maskColor = np.random.random(3)
        colors = [maskColor]
        lws = [1]
    else:
        levels = np.arange(0, len(np.unique(maskFile)))+0.5
        if maskColor == 'random':
            colors = []
            for i in levels:
                maskColor = np.random.random(3)
                colors.append(maskColor)
        else:
            colors = [maskColor]*len(levels)
        lws = [1]*len(levels)
    ax.contour(XC, YC, maskToShow.T, levels=levels, colors=colors,
               linewidths = lws)


    maskToShow = np.log10(maskToShow)
    maskToShow[maskToShow>-np.inf] = 1
    Y, X = np.meshgrid(xTick+1, yTick+1)
    ax.pcolormesh(X+0.5, Y+0.5, maskToShow.T,
                  alpha=0.45, zorder=100,
                  cmap='Blues_r')