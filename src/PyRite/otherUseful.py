import numpy as np
from uncertainties import unumpy as unp
from astropy.io import fits
from astropy import units as u
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
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
        flux_err = -99.9
    else:
        x = unp.uarray(mag, mag_err)
        flux = 10**((23.9-x)/2.5)*1e-6
        flux, flux_err = unp.nominal_values(flux), unp.std_devs(flux)
    return flux, flux_err

def jy_to_mag(flux, fluxerr = None):
    if fluxerr== None:
        mag = -2.5*np.log10(flux/1e-6)+23.9
        mag_err = -99.9
    else:
        x = unp.uarray(flux, fluxerr)
        mag = -2.5*np.log10(x)+23.9
        mag, mag_err = unp.nominal_values(mag), unp.std_devs(mag) 
    return mag, mag_err

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

def make_cutout(big_image, coord, width_arcsec, output_name, ext = 1):
    hdu = fits.open(big_image)  # loading the fits file
    data = hdu[ext].data  # data from the fits file
    try:
        header = hdu[1].header  # header information from the fits file (this contains a lot of info like the pixel size, number of pixels, convert pixel values to sky coordinates etc)
    except:
        header = hdu[0].header
    header['EXPTIME'] = 1
    size = u.Quantity([width_arcsec, width_arcsec], u.arcsec)  # defining the width of cutout in astropy
    wcs = WCS(header)  # information on the pixel to sky coordinate conversion from the header
    cutout = Cutout2D(data, coord, size, wcs,
                      fill_value=0)  # Making the cutout from the given data, around the specified coordinates and size
    cutout_header = cutout.wcs.to_header()  # A new header information after making the cutout (this is because after the cutout, the total number of pixels changed and this changes the pixel to sky coordinate conversion)
    fits.writeto(output_name, cutout.data, cutout_header,
                 overwrite=True)  # writing the cutout data to a new fits file with the output filename given above
