import numpy as np
import os, shutil
from astropy.io import fits
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
import sep
from matplotlib.patches import Ellipse
from .otherUseful import *
from .cosmologyRelated import cosmo

class AutoGalfitClass:
    def __init__(self, imagePath, zeropoint, fitObjects = [], psfFile = 'PSF.fits',
                 ext = 1, dx = -1,
                 outputFile = "model.fits", sigmaFile = "none",
                 maskFile = "none",consFile = "none", imageRegion = [-1],
                 sizeConvolution = [-1], optimize = True, circuralizedRe = True,
                 redshift = -99):
        """
        :param imagePath: path to the image (in the same directory)
        :param zeropoint: value of zeropoint
        :param fitObjects: list of FitObjects to fit
        :param psfFile: path to psf file
        :param ext: extension of image to be used default 0
        :param dx: pixelSize
        :param outputFile: how to call the output model, default model.fits
        :param sigmaFile: path to sigma image, default none
        :param maskFile: path to mask, default none
        :param consFile: path to constrains file, default non
        :param imageRegion: list of points definding the region to fit [x0, x1, y0, y1]
              if x0 == -1, whole image is used, default [-1]
        :param sizeConvolution: list of points definding convolution size [x,y]
              if x == -1, the size of imageRegion is used, default [-1]
        :param optimize: if True, run galfit
        :param circuralizedRe: if True, circuralized radius will be saved
        """
        self.galfitPath = "/home/lisieckik/Galfit/galfit"
        self.imagePath = imagePath
        self.ext = ext
        self.header = fits.open(imagePath)
        self.image = self.header[ext].data
        self.header = self.header[ext].header
        self.fitObjects = fitObjects
        self.zeropoint = zeropoint
        self.dxName = 'Arcsec'
        self.__workingDirectory = ""
        self.cRe = circuralizedRe
        if dx == -1:
            try:
                self.dx = abs(self.header['CD1_1'])*3600
            except:
                raise Exception("You need to specify pixelSize!")
        else:
            self.dx = dx
        
        self.redshift = redshift
        if redshift >0:
            self.D_A = cosmo.arcsec_per_kpc_proper(redshift)
        
        self.outputFile = outputFile
        self.sigmaFile = sigmaFile
        self.psfFile = psfFile
        self.maskFile = maskFile
        self.consFile = consFile
        if imageRegion[0] == -1:
            self.imageRegion = [0, self.image.shape[0], 0, self.image.shape[1]]
        else:
            self.imageRegion = imageRegion
        if sizeConvolution[0] == -1:
            self.sizeConvolution = [self.imageRegion[1]-self.imageRegion[0], self.imageRegion[3]-self.imageRegion[2]]
        else:
            self.sizeConvolution = sizeConvolution
        self.optimize = optimize
    def showImage(self, scale = '', maskOn = False, imageSave = False, reverseMask = False, imageName = "image"):
        """
        :param scale: if log, image is in log
        :param maskOn: if True, mask is shown
        :return:
        """
        xRange = (self.imageRegion[1]-self.imageRegion[0]) * self.dx
        yRange = (self.imageRegion[3]-self.imageRegion[2]) * self.dx
        xTick = np.linspace(-xRange / 2, xRange / 2, (self.imageRegion[1]-self.imageRegion[0]) + 1)
        yTick = np.linspace(-yRange / 2, yRange / 2, (self.imageRegion[3]-self.imageRegion[2]) + 1)
        Y, X = np.meshgrid(xTick, yTick)

        if scale == 'log':
            toShow = self.image[self.imageRegion[0]:self.imageRegion[1],
                                self.imageRegion[2]:self.imageRegion[3]] - np.min(self.image[self.imageRegion[0]:self.imageRegion[1],
                                                                                             self.imageRegion[2]:self.imageRegion[3]])
            aSort = np.argsort(toShow)
            toShow = toShow + toShow[aSort[0,2], aSort[1,2]]
            toShow = np.log10(toShow)
            vM = np.percentile(toShow, [10, 99.])

        else:
            toShow = self.image[self.imageRegion[0]:self.imageRegion[1],
                                self.imageRegion[2]:self.imageRegion[3]]
            vM = np.percentile(toShow, [10, 95.])

        if maskOn:
            if self.maskFile == 'none':
                print("There is no mask to show!")
            else:
                maskToShow = np.zeros(self.image.shape)
                maskPoints = np.loadtxt(self.maskFile).astype(int)
                for i in maskPoints:
                    if reverseMask:
                        maskToShow[i[1], i[0]] = 1
                    else:
                        maskToShow[i[0], i[1]] = 1
                maskToShow = maskToShow[self.imageRegion[0]:self.imageRegion[1],
                                    self.imageRegion[2]:self.imageRegion[3]]
                plt.pcolormesh(X,Y,maskToShow.T, alpha = 0.1, zorder =100, cmap = 'Grays_r')

                xTick = np.linspace(-xRange / 2, xRange / 2, (self.imageRegion[1]-self.imageRegion[0]))
                yTick = np.linspace(-yRange / 2, yRange / 2, (self.imageRegion[3]-self.imageRegion[2]))
                YC, XC = np.meshgrid(xTick + (xTick[1]-xTick[0])/2, yTick+ (yTick[1]-yTick[0])/2)
                plt.contour(XC,YC, maskToShow.T, levels = [0.5], colors = ['r'], lw = '3')

        plt.pcolormesh(X,Y, toShow.T, vmin = vM[0], vmax = vM[1])
        plt.xlabel(self.dxName)
        plt.ylabel(self.dxName)
        if imageSave:
            plt.savefig('%s%s.png'%(self.__workingDirectory, imageName))
            plt.close()
        else:
            plt.show()
    def setWorkingDir(self, workingDirPath):
        """
        :param workingDirPath: path to copy all files and run galfit there
                                allows to easy split run into threads
        :return:
        """
        self.__workingDirectory = workingDirPath
        if self.__workingDirectory != "" and not os.path.exists(workingDirPath):
            os.mkdir(workingDirPath)
        if workingDirPath != "":
            shutil.copy2(self.imagePath, "%s%s"%(self.__workingDirectory, self.imagePath.split("/")[-1]))
            self.imagePath = self.imagePath.split("/")[-1]

            if self.sigmaFile != 'none':
                shutil.copy2(self.sigmaFile, "%s%s"%(self.__workingDirectory,self.sigmaFile.split("/")[-1]))
                self.sigmaFile = self.sigmaFile.split("/")[-1]

            shutil.copy2(self.psfFile, "%s%s"%(self.__workingDirectory,self.psfFile.split("/")[-1]))
            self.psfFile = self.psfFile.split("/")[-1]
            if self.maskFile != 'none':
                shutil.copy2(self.maskFile, "%s%s"%(self.__workingDirectory,self.maskFile.split("/")[-1]))
                self.maskFile = self.maskFile.split("/")[-1]
            if self.consFile != 'none':
                shutil.copy2(self.consFile, "%s%s"%(self.__workingDirectory,self.consFile))
        if workingDirPath[-1] != "/" and workingDirPath[-1] != "\\":
            self.__workingDirectory += "/"
    def prepFeedmeFile(self):
        """
        Prepares the feedme file with use of fitObjects
        :return:
        """
        self.__inputFile__ =  "A) %s[%i]  # Input data image (FITS file)\n"%(self.imagePath, self.ext)
        self.__inputFile__ += "B) %s  # Output data image block\n"%(self.outputFile)
        self.__inputFile__ += "C) %s  # Sigma image name (made from data if blank or \"none\")\n"%(self.sigmaFile)
        self.__inputFile__ += "D) %s  # Input PSF image and (optional) diffusion kernel\n"%(self.psfFile)
        self.__inputFile__ += "E) 1  # PSF fine sampling factor relative to data\n"
        self.__inputFile__ += "F) %s  # Bad pixel mask (FITS image or ASCII coord list)\n"%(self.maskFile)
        self.__inputFile__ += "G) %s  # File with parameter constraints (ASCII file)\n"%(self.consFile)
        self.__inputFile__ += "H) %i %i %i %i # Image region to fit (xmin xmax ymin ymax)\n"%(self.imageRegion[0],self.imageRegion[1],self.imageRegion[2],self.imageRegion[3],)
        self.__inputFile__ += "I) %i %i # Size of the convolution box (x y)\n"%(self.sizeConvolution[0],self.sizeConvolution[1])
        self.__inputFile__ += "J) %f  # Magnitude photometric zeropoint\n"%(self.zeropoint)
        self.__inputFile__ += "K) %f %f # Plate scale (dx dy)    [arcsec per pixel]\n"%(self.dx, self.dx)
        self.__inputFile__ += "O) regular  # Display type (regular, curses, both)\n"

        if self.optimize:
            self.__inputFile__ += "P) 0  # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n1\n\n"
        else:
            self.__inputFile__ += "P) 1  # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n1\n\n"

        n0 = 1
        for o in self.fitObjects:
            o.checkPos(self.image)
            o.setText()
            self.__inputFile__ += o.combineText(n0)
            self.__inputFile__ += '\n\n\n'
            n0+=1
        try:
            resFile = open("%s/galfit.feedme"%(self.__workingDirectory), 'w')
        except:
            resFile = open("galfit.feedme", 'w')
        resFile.write(self.__inputFile__)
        resFile.close()
    def runGalfit(self, parameterOutputFile = 'bestFit.txt', runID = -99, saveImage = True, percentagesContrast = [0,98]):
        """
        Runs galfit once
        :param parameterOutputFile: File to save the results, str
        :param runID: How to call output in results, int
        :param saveImage: Do you want to save image? Bool
        :return:
        """
        self.prepFeedmeFile()
        proc = self.galfitPath + " galfit.feedme"
        if self.__workingDirectory != "":
            process = Popen(proc, stdout=PIPE, stderr=PIPE, shell=True, cwd='%s' %self.__workingDirectory)
        else:
            process = Popen(proc, stdout=PIPE, stderr=PIPE, shell=True)

        stdout, stderr = process.communicate()

        chi2, res = ReadResults("%sgalfit.01"%self.__workingDirectory)
        __, inp = ReadResults("%sgalfit.feedme"%self.__workingDirectory)


        if self.__workingDirectory != "":
            resultFile = '%s%s' %(self.__workingDirectory, parameterOutputFile)
        else:
            resultFile = '%s' %(parameterOutputFile)
        
        if self.redshift >0:
            self.D_A = cosmo.arcsec_per_kpc_proper(self.redshift)

        if not os.path.exists(resultFile) or runID == -99:
            resultFile = open(resultFile, 'w')

            resultFile.write("# runID chi2")
            for k in inp.keys():
                resultFile.write(" x_%i_in y_%i_in" % (k, k))
                resultFile.write(" mag_%i_in" % (k))
                if inp[k].objectType == 'sersic':
                    resultFile.write(" apix_%i_in n_%i_in" % (k, k))
                    resultFile.write(" ba_%i_in PA_%i_in" % (k, k))

            for k in inp.keys():
                resultFile.write(" x_%i_out y_%i_out" % (k, k))
                resultFile.write(" mag_%i_out" % (k))
                if inp[k].objectType == 'sersic':
                    resultFile.write(" apix_%i_out n_%i_out" % (k, k))
                    resultFile.write(" ba_%i_out PA_%i_out" % (k, k))
                    if self.cRe:
                        resultFile.write(" circuralizedRe_%i_out"%(k))
                        if self.redshift>0:
                            resultFile.write(" circuralizedRe_kpc_%i_out"%(k))

            resultFile.write('\n')
        else:
            resultFile = open(resultFile, 'a')



        resultFile.write('%i %f '%(runID, chi2))
        for k in inp.keys():
            resultFile.write(" %.3f %.3f" % (inp[k].X[0], inp[k].Y[0]))
            resultFile.write(" %.4f" % (inp[k].mag[0]))
            if inp[k].objectType == 'sersic':
                resultFile.write(" %.4f %.4f" % (inp[k].apix[0], inp[k].n[0]))
                resultFile.write(" %.4f %.4f" % (inp[k].ba[0], inp[k].PA[0]))

        for k in inp.keys():
            if chi2>0:
                resultFile.write(" %.3f %.3f" % (res[k].X[0], res[k].Y[0]))
                resultFile.write(" %.4f" % (res[k].mag[0]))
                if inp[k].objectType == 'sersic':
                    resultFile.write(" %.4f %.4f" % (res[k].apix[0], res[k].n[0]))
                    resultFile.write(" %.4f %.4f" % (res[k].ba[0], res[k].PA[0]))
                    if self.cRe:
                        cirRe = res[k].apix[0] * np.sqrt(res[k].ba[0]) * self.dx
                        resultFile.write(" %.4f"%(cirRe))
                        if self.redshift>0:
                            resultFile.write(" %.4f"%((cirRe/self.D_A).value))

            else:
                resultFile.write(" %.3f %.3f" % (-1,-1))
                resultFile.write(" %.4f" % (-1))
                if inp[k].objectType == 'sersic':
                    resultFile.write(" %.4f %.4f" % (-1,-1))
                    resultFile.write(" %.4f %.4f" % (-1,-1))
                    if self.cRe:
                        cirRe = -1
                        resultFile.write(" -1")
                        if self.redshift>0:
                            resultFile.write(" -1")
        resultFile.write('\n')
        resultFile.close()

        if saveImage:
            try:
                fig = plt.figure(figsize=(9, 3))
                data = fits.open('%smodel.fits'%(self.__workingDirectory))

                image0 = data[1].data
                model = data[2].data
                residual = data[3].data

                mV = np.min(residual)
                image0 = image0 - mV
                model = model - mV
                residual = residual - mV

                uV = np.sort(np.unique(residual))[1]*0.1
                image0 = image0 + uV
                model = model + uV
                residual = residual + uV

                v = np.percentile(np.log10(image0), percentagesContrast)
                plt.subplot(131)
                plt.imshow(np.log10(image0), vmin = v[0], vmax=v[1])


                plt.subplot(132)
                plt.imshow(np.log10(model), vmin = v[0], vmax=v[1])


                plt.subplot(133)
                plt.imshow(np.log10(residual), vmin = v[0], vmax=v[1])

                plt.tight_layout()
                plt.savefig('%s%iresult.png'%(self.__workingDirectory, abs(runID)))
                data.close()
                plt.close()
            except:
                print("Image could not be produce!")
    def galfitBootstrap(self, n = 100, saveModels = False, parameterOutputFile = 'bestFit.txt', percentagesContrast = [0,98]):
        """
        Allows for multiple galfit run with random parameters defined in FitObjects
        :param n: Number of runs, int
        :param saveModels: Do you want to save all the models? bool
        :param parameterOutputFile: How to call output file, str
        :return:
        """
        try:
            os.remove('%s%s' %(self.__workingDirectory, parameterOutputFile))
        except:
            pass

        if saveModels:
            if self.__workingDirectory == "":
                if not os.path.exists('%s' %('models')):
                    os.mkdir('%s' %('models'))
            else:
                if not os.path.exists('%s%s' %(self.__workingDirectory, 'models')):
                    os.mkdir('%s%s' %(self.__workingDirectory, 'models'))
        for i in range(n):
            if np.round(i/n,4)*100 % 10 == 0:
                print(np.round(i/n,2)*100, "% Done")
            self.prepFeedmeFile()
            self.runGalfit(runID=i, saveImage=saveModels, parameterOutputFile = parameterOutputFile, percentagesContrast = percentagesContrast)
            try:
                os.remove("%sgalfit.01"%self.__workingDirectory)
            except:
                pass
            if saveModels:
                shutil.copy2('%s%s' %(self.__workingDirectory, self.outputFile), 
                             '%smodels/%i_%s' %(self.__workingDirectory, i, self.outputFile))
                shutil.copy2('%s%iresult.png'%(self.__workingDirectory, abs(i)), 
                             '%smodels/%iresult.png'%(self.__workingDirectory, abs(i)))
                os.remove('%s%iresult.png'%(self.__workingDirectory, abs(i)))
    def prepareAuxFiles(self, newName='', ap = 3, substractBKG = True):
        """
        This functions prepares sigma image, mask
        and background substracted image with use of sep
        :param ap: number of additional pixels for masking, int
        :return:
        """
        sex = QuickSextractor(self.image, self.header)
        if newName == '':
            sex.saveSigmaImage("%ssigmaImageQS.fits"%(self.__workingDirectory))
            self.sigmaFile = "%ssigmaImageQS.fits"%(self.__workingDirectory)
        else:
            sex.saveSigmaImage("%ssigmaImageQS_%s.fits"%(self.__workingDirectory, newName))
            self.sigmaFile = "%ssigmaImageQS_%s.fits"%(self.__workingDirectory, newName)

        if substractBKG:
            if newName == '':
                self.image = sex.saveBKGsubImage("%sbkgsub_ImageQS.fits"%(self.__workingDirectory))
                self.imagePath = "%sbkgsub_ImageQS.fits"%(self.__workingDirectory)
            else:
                self.image = sex.saveBKGsubImage("%sbkgsub_ImageQS_%s.fits"%(self.__workingDirectory, newName))
                self.imagePath = "%sbkgsub_ImageQS_%s.fits"%(self.__workingDirectory, newName)


        greymask = np.zeros(self.image.shape)
        for o in sex.objects:
            x = o['x']
            y = o['y']
            a = o['a']
            b = o['b']
            t = np.rad2deg(o['theta'])
            r = np.sqrt((x - self.image.shape[0] / 2 + 1) ** 2 + (y - self.image.shape[1] / 2 + 1) ** 2)
            if r>5:
                nIter = 0
                centerEmpty = False
                if r/a <7:
                    mask = is_pixel_in_ellipse(self.image.shape,
                            (x - 1, y - 1), a, b, t, scale=r/a*0.7, ap=ap)
                else:

                    mask = is_pixel_in_ellipse(self.image.shape,
                            (x - 1, y - 1), a, b, t, scale=5, ap=ap)

                greymask[mask] = 1

        if newName == '':
           mask = open('%smaskQS.txt'%(self.__workingDirectory), 'w')
        else:
           mask = open('%smaskQS_%s.txt'%(self.__workingDirectory, newName), 'w')
        

        for i in range(greymask.shape[0]):
            for j in range(greymask.shape[1]):
                if greymask[i, j] == 1:
                    mask.write('%i %i\n' % (i, j))
        mask.close()
        if newName == '':
            self.maskFile = '%smaskQS.txt'%(self.__workingDirectory)
        else:
            self.maskFile = '%smaskQS_%s.txt'%(self.__workingDirectory, newName)

        if self.__workingDirectory != "":
            self.setWorkingDir(self.__workingDirectory)
        self.ext = 0
    def prepareConsFile(self, consName='consQS.txt'):
        """
        This method prepares the the constrains file. Removes the previous one if exists!
        :param consName: Name of the constrains file, deflaut: consQS.txt
        :return:
        """
        self.consFile = consName
        if os.path.exists("%s%s"%(self.__workingDirectory, consName)):
            os.remove("%s%s"%(self.__workingDirectory, consName))
            self.consFile = "%s%s"%(self.__workingDirectory, consName)
        f = open(consName, 'w')
        f.close()
    def addCon(self, con):
        """
        This method adds one line of constrain to the automatic file.
        :param con: Line defining the constrain as galfit needs it, str
        eg:
            2-1 x -2 2     meaning x parameter of 2nd object has to be within -2,2 pixels from 1st
        """
        self.prepareConsFile()
        f = open(self.consFile, 'a')
        f.write(con + '\n')
class FitObject:
    def __init__(self, objectType, posX=-1, posY=-1, freeX = True, freeY = True,
                 mag = 22, freeMag = True, apix = 3, freeApix = True,
                 n = 2, freeN = True, ba = 0.5, freeBA = True, PA = 0, freePA = True, addModel = True,
                 randomX = 0, randomY = 0, randomMag = 0, randomApix = 0, randomN = 0, randomBA = 0,
                 randomPA = 0):
        """
        This class holds are information for galfit to fit.
        :param objectType: What kind of object you want to fit? Possible values: sersic, psf
        :param posX: Position on x [pixels] of the obejct, default center of the image, int
        :param posY: Position on y [pixels] of the obejct, default center of the image, int
        :param freeX: Keep x a free parameter, bool, default: True
        :param freeY: Keep y a free parameter, bool, default: True
        :param mag: Value of magnitude for start guess, float, default: 22
        :param freeMag: Keep mag a free parameter, bool, default: True
        :param apix: Value of semi-major axis in pix for start guess, float, default: 3
        :param freeApix: Keep apix a free parameter, bool, default: True
        :param n: Value of sersic index, float, default: 2
        :param freeN: Keep apix n free parameter, bool, default: True
        :param ba: Value of ellipticity, float, default: 0.5
        :param freeBA: Keep apix ba free parameter, bool, default: True
        :param PA: Value of position angle, float, default: 0
        :param freePA: Keep apix PA free parameter, bool, default: True
        :param addModel: Get the model into residual, default: True
        :param randomX: Range for random values, if 1, value will be found randomly as x-1, x+1, defualt: 0
        :param randomY: Range for random values, if 1, value will be found randomly as y-1, y+1, defualt: 0
        :param randomMag: Range for random values, if 1, value will be found randomly as mag-1, mag+1, defualt: 0
        :param randomApix: Range for random values, if 1, value will be found randomly as apix-1, apix+1, defualt: 0
        :param randomN: Range for random values, if 1, value will be found randomly as n-1, n+1, defualt: 0
        :param randomBA: Range for random values, if 1, value will be found randomly as ba-1, ba+1, defualt: 0
        :param randomPA: Range for random values, if 1, value will be found randomly as PA-1, PA+1, defualt: 0
        """
        self.X = [posX, freeX]
        self.Y = [posY, freeY]
        self.mag = [mag, freeMag]
        self.apix = [apix, freeApix]
        self.n = [n, freeN]
        self.ba = [ba, freeBA]
        self.PA = [PA, freePA]

        self.addModel = addModel
        self.objectType = objectType
        self.randomX = randomX
        self.randomY = randomY
        self.randomMag = randomMag
        self.randomApix = randomApix
        self.randomN = randomN
        self.randomBA = randomBA
        self.randomPA = randomPA

        self.setText()
    def setText(self):
        """
        Prepares text to save in feedme
        :return:
        """
        if self.randomX >0:
            xHere = self.X[0] + np.random.uniform(-self.randomX,self.randomX)
        else:
            xHere = self.X[0]

        if self.randomY > 0:
            yHere = self.Y[0] + np.random.uniform(-self.randomY, self.randomY)
        else:
            yHere = self.Y[0]
        if self.randomMag > 0:
            magHere = self.mag[0] + np.random.uniform(-self.randomMag, self.randomMag)
        else:
            magHere = self.mag[0]
        if self.randomApix > 0:
            apixHere = self.apix[0] + np.random.uniform(-self.randomApix, self.randomApix)
        else:
            apixHere = self.apix[0]

        if self.randomN > 0:
            nHere = self.n[0] + np.random.uniform(-self.randomN, self.randomN)
        else:
            nHere = self.n[0]

        if self.randomBA > 0:
            baHere = self.ba[0] + np.random.uniform(-self.randomBA, self.randomBA)
        else:
            baHere = self.ba[0]

        if self.randomPA > 0:
            PAHere = self.PA[0] + np.random.uniform(-self.randomPA, self.randomPA)
        else:
            PAHere = self.PA[0]


        self.text =     ["0) %s             #  object type"%self.objectType]
        if self.X[1] and self.Y[1]:
            self.text.append("1) %.2f %.2f   1 1    #  position x, y" %(xHere, yHere))
        elif freeX:
            self.text.append("1) %.2f %.2f   1 0    #  position x, y" %(xHere, yHere))
        else:
            self.text.append("1) %.2f %.2f   0 0    #  position x, y" %(xHere, yHere))

        if self.mag[1]:
            self.text.append("3) %.2f        1      #  Integrated magnitude" %(magHere))
        else:
            self.text.append("3) %.2f        0      #  Integrated magnitude" %(magHere))

        if self.objectType == 'sersic':
            if self.apix[1]:
                self.text.append("4) %.2f        1      #  R_e (half-light radius)   [pix]" %(apixHere))
            else:
                self.text.append("4) %.2f        0      #  R_e (half-light radius)   [pix]" %(apixHere))

            if self.n[1]:
                self.text.append("5) %.2f        1      #  Sersic index n (de Vaucouleurs n=4)" %(nHere))
            else:
                self.text.append("5) %.2f        0      #  Sersic index n (de Vaucouleurs n=4)" %(nHere))

            if self.ba[1]:
                self.text.append("9) %.2f        1      #  axis ratio (b/a)" %(baHere))
            else:
                self.text.append("9) %.2f        0      #  axis ratio (b/a)" %(baHere))

            if self.PA[1]:
                self.text.append("10) %i         1      #  position angle (PA) [deg: Up=0, Left=90]" %(PAHere))
            else:
                self.text.append("10) %i         0      #  position angle (PA) [deg: Up=0, Left=90]" %(PAHere))

        if self.addModel:
            self.text.append("Z) 0                  #  output option (0 = resid., 1 = Don't subtract) )")
        else:
            self.text.append("Z) 1                  #  output option (0 = resid., 1 = Don't subtract) )")
    def checkPos(self, image):
        """
        Assigns central position to the obejct.
        :param image: Array.
        :return:
        """
        if self.X[0] == -1 or self.Y[0] == -1:
            self.text[1] = "1) %.2f %.2f   1 1    #  position x, y" % (image.shape[0]/2, image.shape[1]/2)
            self.X = [image.shape[0]/2, self.X[1]]
            self.Y = [image.shape[1]/2, self.Y[1]]
    def combineText(self, objectNumber):
        """
        Allows access to the text needed to build feedme.
        :param objectNumber: id of the object
        :return:
        """
        T = "# Object number: %i\n"%objectNumber
        for i in self.text:
            T += i
            T += "\n"
        return T
    def printParams(self, objectNumber = -1):
        print(self.combineText(objectNumber))
class QuickSextractor:
    def __init__(self, image, header, thresh = 1.3, deblend_nthresh = 32, deblend_cont = 0.0001):
        """
        This class runs sextractor with python and allows user to build necessary files for galfit.
        :param image: Array
        :param header: Header of the image
        :param thresh: Sextractor parameter.
        :param deblend_nthresh: Sextractor parameter.
        :param deblend_cont: Sextractor parameter.
        """
        self.image = image
        self.image = self.image.byteswap().newbyteorder()
        self.image = self.image.astype(np.float32)

        self.mean, self.std = np.mean(self.image), np.std(self.image)
        bkg = sep.Background(self.image)
        self.bkg_image = bkg.back()
        self.bkg_rms = bkg.rms()
        self.header = header
        self.header['EXPTIME'] = 1
        self.data_sub = self.image - self.bkg_image

        self.objects = sep.extract(self.data_sub, thresh,
                                   err=bkg.globalrms,
                                   deblend_nthresh=deblend_nthresh,
                                   deblend_cont=deblend_cont)
    def saveSigmaImage(self, outputFile):
        fits.writeto(outputFile, self.bkg_rms, header=self.header, overwrite=True)
    def saveBKGsubImage(self, outputFile):
        fits.writeto(outputFile, self.data_sub, header=self.header, overwrite=True)
        return self.data_sub
    def showDetections(self):
        """
        This method shows the detected objects on the image.
        """
        toShow = self.image - np.min(self.image)
        aSort = np.argsort(toShow)
        toShow = toShow + toShow[aSort[0, 2], aSort[1, 2]]
        toShow = np.log10(toShow)
        vM = np.percentile(toShow, [10, 99.])
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        for o in self.objects:
            ell = Ellipse(xy = (o['x'], o['y']), width=o['a']*6, height=o['b']*6,
                          angle=np.rad2deg(o['theta']), edgecolor='r', fc = 'None')
            ax.add_patch(ell)
        ax.imshow(toShow, vmin=vM[0], vmax=vM[1])
        ax.scatter(self.objects['x'], self.objects['y'], marker = 'x', c = 'r',zorder = 20, s =20)
        plt.show()
        plt.close()
    def makeFitObject(self, objectType, zeropoint,
                      randomMag = 0.2, randomApix=1,n = 3, randomN=1, randomBA=0.2, randomPA=45, randomX=5, randomY=5):
        dx_pos = self.objects['x'] - self.image.shape[0]
        dy_pos = self.objects['y'] - self.image.shape[1]
        r = np.sqrt(dx_pos**2 + dy_pos**2)
        r0 = np.argsort(r)[0]
        print(self.objects[r0])

        FO = FitObject(objectType,
                    mag=zeropoint -2.5*np.log10(self.objects['flux'][r0]), randomMag=randomMag,
                    apix=self.objects['a'][r0]*0.4, randomApix=randomApix,
                    n=n, randomN=randomN,
                    ba=self.objects['b'][r0]/self.objects['a'][r0], randomBA=randomBA, 
                    PA=self.objects['theta'][r0] + 90,
                    randomPA=randomPA, randomX=randomX, randomY=randomY)
        return FO

def ReadResults(file):
    """
    This function reads the results from galfit file and saves it into FitObject class.
    :param file: galfit result file eg. galfit.01
    :return:
    """
    chi2 = -1
    try:
        file = open(file, 'r').readlines()
        results = {}
        comp = False
        for l in file:
            lsplit = l.split()
            if len(lsplit)>1:
                if lsplit[0] == "#" and lsplit[1] == "Chi^2/nu":
                    chi2 = float(lsplit[3][0:-1])

                if lsplit[0] == "#" and (lsplit[1] == "Component" or lsplit[1] == "Object"):
                    comp = True
                    cNumber = int(lsplit[-1])
                if comp:
                    if lsplit[0] == "0)":
                        cHere = FitObject(lsplit[1])
                    elif lsplit[0] == "1)":
                        cHere.X = [float(lsplit[1]), int(lsplit[3])]
                        cHere.Y = [float(lsplit[2]), int(lsplit[4])]
                    elif lsplit[0] == "3)":
                        cHere.mag = [float(lsplit[1]), int(lsplit[2])]
                    elif lsplit[0] == "4)":
                        cHere.apix = [float(lsplit[1]), int(lsplit[2])]
                    elif lsplit[0] == "5)":
                        cHere.n = [float(lsplit[1]), int(lsplit[2])]
                    elif lsplit[0] == "9)":
                        cHere.ba = [float(lsplit[1]), int(lsplit[2])]
                    elif lsplit[0] == "10)":
                        cHere.PA = [float(lsplit[1]), int(lsplit[2])]
                    elif lsplit[0] == "Z)":
                        comp = False
                        cHere.setText()
                        results[cNumber] = cHere
    except:
        results = {}
    return chi2, results