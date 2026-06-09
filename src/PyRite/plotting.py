import pathlib
import matplotlib.pyplot as plt
import PyRite
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FixedLocator

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


def updateAxisColor(c):
    
    plt.rcParams.update({
        "axes.edgecolor": c,
        "text.color" : c,
        "axes.labelcolor": c,
        "axes.titlecolor": c,
        "xtick.color": c,
        "ytick.color":c,
    })

def remakeTicks(x, y, ax=None, xoff = 0, yoff = 0):
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
                visibleHere = visibility[y*(ny-1) + x]
            
            ax = fig.add_axes([x0Here, y0Here,dxHere, dyHere], visible = visibleHere)
            axes.append(ax)
            x0Here += dxHere + xspace
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
    return axes