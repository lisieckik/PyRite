import numpy as np
from uncertainties import unumpy as unp
from astropy.io import fits
from astropy import units as u
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
import json
from importlib.resources import files
import subprocess

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

def make_cutout(big_image, ra, dec, width_arcsec, output_name, ext=1, verbose = False, removeSIP = False):
    if isinstance(big_image, str):
        hdu = fits.open(big_image)
    else:
        hdu = big_image

    data = hdu[ext].data
    header = hdu[ext].header
    header['EXPTIME'] = 1
    size = u.Quantity([width_arcsec, width_arcsec], u.arcsec)
    if removeSIP:
        for key in list(header):
            if (
                key.startswith('A_') or
                key.startswith('B_') or
                key.startswith('AP_') or
                key.startswith('BP_')
                ):
                del header[key]
    wcs = WCS(header, relax=True)
    coord = SkyCoord(
        ra=ra* u.deg,
        dec=dec * u.deg,
    )
    if verbose:
        position = wcs.world_to_pixel(coord)
        print("Image shape:", data.shape)
        print("Object pixel position:", position)

    cutout = Cutout2D(
        data,
        coord,
        size,
        wcs,
        fill_value=0)
    cutout_header = cutout.wcs.to_header()
    fits.writeto(
        output_name,
        cutout.data,
        cutout_header,
        overwrite=True
    )
    hdu.close()

def cigale_filters(v=25):
    path = files(__package__) / f"CIGALE{v}_filters_parsed.json"
    with path.open("r") as f:
        loaded = json.load(f)
    return loaded

def prepareKernel(kernelFWHM = 2., kernelSize = 5, output = ''):
    """
    :param kernelFWHM: size of gaussian (pix)
    :param kernelSize: size of the kernel file (pix)
    :param output: name of the output, if default, file will be named kernel_%.2f_%i.fits with both sizes
    :return: output
    """
    if output == '':
        output = 'kernel_%.2f_%i.fits'%(kernelFWHM, kernelSize)
    result = subprocess.run(
        [
            "astmkprof",
            '--kernel=gaussian,%.2f,%.2f'%(kernelFWHM, kernelSize),
            '--oversample=1',
            "--output=%s" % output
        ],
        capture_output=True,
        text=True
    )
    print('File created:', output)
    return output