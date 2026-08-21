import os
import subprocess
import numpy as np
from astropy.io import fits
from astropy.table import Table
from .otherUseful import prepareKernel

os.environ["PATH"] = (
    "/home/lisieckik/local/bin:"
    + os.environ["PATH"]
)

def runNoiseChisel(inputImage,
                   snquant = 0.99,
                   qthresh = 0.3,
                   dt = 0.,
                   erode = 2,
                   erodengb = 4,
                   snminarea=10,
                   minskyfrac=0.7,
                   outliernumngb = 15,
                   interpnumngb = 15,
                   largeTileSize = 150,
                   tileSize = 30,
                   output = 'detected.fits',
                   kernel = '',
                   parameters = '',
                   h = 1,
                   verbose = False):
    """
    :param inputImage: path to the input map
    :param snquant: GnuAstro docs
    :param qthresh: GnuAstro docs
    :param dt: GnuAstro docs
    :param erode: GnuAstro docs
    :param erodengb: GnuAstro docs
    :param snminarea: GnuAstro docs
    :param minskyfrac: GnuAstro docs
    :param outliernumngb: GnuAstro docs
    :param output: output path
    :param kernel: path to kernel
    :param parameters: some pretested sets of parameters, available options: JADES
    :param h: hdu index
    :param verbose: if False, no messages appear
    :return: output path
    """
    if parameters == 'JADES':
        snquant = 0.99,
        qthresh = 0.33,
        dt = 0.05,
        erode = 3,
        erodengb = 4,

    if kernel =='':
        kernel = 'kernel_%.2f_%i.fits' % (2, 5)
        if not os.path.exists(kernel):
            prepareKernel()


    result = subprocess.run(
        [
            "astnoisechisel",
            inputImage,
            '--hdu=%s'%h,
            '--dthresh=%f' % dt,
            '--snquant=%f' % snquant,
            '--qthresh=%f' % qthresh,
            '--erode=%i' % erode,
            '--erodengb=%i' % erodengb,
            "--output=%s" % output,
            '--kernel=%s'%kernel,
            '--snminarea=%i'%snminarea,
            '--minskyfrac=%f'%minskyfrac,
            '--outliernumngb=%i'%outliernumngb,
            '--interpnumngb=%i'%interpnumngb,
            '--largetilesize=%i,%i'%(largeTileSize,largeTileSize),
            '--tilesize=%i,%i'%(tileSize,tileSize)
        ],
        capture_output=True,
        text=True
    )
    if verbose:
        print(result.stdout)
        print(result.stderr)
    print('File created:', output)
    return output
def prepareUnsharpedImage(inputImage,
                    kernel='',
                    output = 'unsharped_diff_image.fits',
                    keepTempFiles = False,
                    h = 1,
                    verbose = False):
    if kernel =='':
        kernel = 'kernel_%.2f_%i.fits' % (2, 5)
        if not os.path.exists(kernel):
            prepareKernel()


    result = subprocess.run(
        [
            "astconvolve",
            inputImage,
            '-h%s'%h,
            '--kernel=%s'%kernel,
            "--domain=spatial",
            "--output=temp_unsharped.fits"

        ],
        capture_output=True,
        text=True
    )
    if verbose:
        print(result.stdout)
        print(result.stderr)
    print('File created:', 'temp_unsharped.fits')
    result = subprocess.run(
        [
            "astarithmetic",
            inputImage,
            '-h%s'%h,
            'temp_unsharped.fits',
            '-h1',
            '-',
            "--output=%s"%output

        ],
        capture_output=True,
        text=True
    )
    print('File created:', output)
    if verbose:
        print(result.stdout)
        print(result.stderr)
    if not keepTempFiles:
        os.remove('temp_unsharped.fits')
        print('File deleted:', 'temp_unsharped.fits')
    return output
def segmentNoiseChiselResults(detectedImage,
                              kernel='',
                              gthresh=0.5,
                              minriverlength=15,
                              snminarea=15,
                              minskyfrac=0.6,
                              output = 'segmented.fits',
                              verbose = False
                              ):
    if kernel =='':
        kernel = 'kernel_%.2f_%i.fits' % (2, 5)
        if not os.path.exists(kernel):
            prepareKernel()

    result = subprocess.run(
        [
            "astsegment",
            detectedImage,
            '--kernel=%s'%kernel,
            '--gthresh=%.f'%gthresh,
            '--minriverlength=%i'%minriverlength,
            '--output=%s'%output,
            '--snminarea=%i'%snminarea,
            '--minskyfrac=%f'%minskyfrac
        ],
        capture_output=True,
        text=True
    )
    if verbose:
        print(result.stdout)
        print(result.stderr)
    print('File created:', output)
    return output
def makeCatalogue(nchiselImage,
                  zp = 27.99959,
                  output='catalog.fits',
                  verbose = False,onlySB = False):
    if onlySB:
            result = subprocess.run(
        [
            "astmkcatalog",
            nchiselImage,
            '--ids',
            '--x',
            '--y',
            '--area',
            '--sum',
            '--sum-error',
            '--semi-major',
            '--semi-minor',
            '--position-angle',
            '--output=%s'%output
        ],
        capture_output=True,
        text=True
    )
    else:
        result = subprocess.run(
            [
                "astmkcatalog",
                nchiselImage,
                '--zeropoint=%f'%zp,
                '--ids',
                '--x',
                '--y',
                '--ra',
                '--dec',
                '--area',
                '--magnitude',
                '--semi-major',
                '--semi-minor',
                '--position-angle',
                '--output=%s'%output
            ],
            capture_output=True,
            text=True
        )
    if verbose:
        print(result.stdout)
        print(result.stderr)