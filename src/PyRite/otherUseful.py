import numpy as np
from uncertainties import unumpy
from astropy.stats import sigma_clipped_stats

def is_pixel_in_ellipse(image_size, center, a, b, theta, scale = 1, ap = 0):
    # Generate the grid of coordinates
    a = a*scale +ap
    b = b*scale +ap
    y, x = np.indices(image_size)
    theta = np.deg2rad(theta)

    # Shift the coordinates to the ellipse center
    x_shifted = x - center[0]
    y_shifted = y - center[1]

    # Apply the rotation matrix to the shifted coordinates
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    x_rot = x_shifted * cos_theta + y_shifted * sin_theta
    y_rot = -x_shifted * sin_theta + y_shifted * cos_theta

    # Ellipse equation (x_rot / a)^2 + (y_rot / b)^2 <= 1
    ellipse_mask = ((x_rot / a) ** 2 + (y_rot / b) ** 2) <= 1
    return ellipse_mask

def mag_to_jy(mag, mag_err=None):
    if mag_err== None:
        flux = 10**((23.9-mag)/2.5)*1e-6
    else:
        x = unumpy.uarray(mag, mag_err)
        flux = 10**((23.9-x)/2.5)*1e-6
    return flux

def jy_to_mag(flux, fluxerr = None):
    if fluxerr== None:
        mag = -2.5*np.log10(flux/1e-6)+23.9
    else:
        x = unumpy.uarray(flux, fluxerr)
        mag = -2.5*np.log10(x)+23.9
    return mag

def change_image(image):
    """
    Display a galaxy image with automatic astronomical contrast scaling.
    Only requires the image array.
    """

    # Robust background estimate
    mean, median, std = sigma_clipped_stats(image, sigma=3)

    # Shift background to zero
    img = image - median

    # Clip negative/background-dominated pixels
    img = np.clip(img, 0, None)

    # Robust upper scale
    vmax = np.percentile(img, 98)

    # Asinh stretch
    stretched = np.arcsinh(10 * img / vmax)
    stretched /= stretched.max()
    return stretched