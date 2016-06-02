#!/usr/bin/env python
import sys, os
import argparse
import subprocess
import ppgplot
import scipy.optimize
import numpy
import timeClasses
import math
import spectrumClasses

pathToCode = "/Users/rashley/code/cyclotronFitting"
# pathToCode = "."

global iteration, mainPlotWindow, currentPlotWindow

def radians(angle):
	return angle/180. * math.pi

def replaceExpChar(value):
	if 'D' in value:
		value = value.replace('D', 'E')
	
	return value
	
def getModelSpectrum(angle, field, temperature, log_lambda, geometry, intensityparameter):
	modelParams = {}
	modelParams['angle'] = float(angle)
	modelParams['field'] = float(field)
	modelParams['temperature'] = float(temperature)
	modelParams['log_lambda'] = float(log_lambda)
	modelParams['geometry'] = geometry
	print "Creating a model for: ", modelParams
	filename = "/tmp/%3.2f_%3.2f_%3.2f_%3.2f_%3.2f.dat"%(angle, field, temperature, log_lambda, geometry)
	if not os.path.exists(filename):
		parameterFile = open("ConstLambda_Ein", "w")
		parameterFile.write("Sichtwinkel      [grad] : %f\n"%angle)
		parameterFile.write("Magnetfeld       [MG]   : %f\n"%field)
		parameterFile.write("Temperatur       [keV]  : %f\n"%temperature)
		parameterFile.write("log(Lambda)      [1]    : %f\n"%log_lambda)
		parameterFile.write("Geometrie= 0 oder 1     : 0\n")
		parameterFile.close()

		modelCommand = [pathToCode + "/ConstLambda.bin"]
		outfile = open(filename, "w")
		subprocess.call(modelCommand, stdout = outfile)
		outfile.close()
					
	# Load the computed model
	print "Loading model:", filename		
	inputfile = open(filename, 'r')
	freqs = []
	i_0s = []
	i_1s = []
	i_s = []
	for line in inputfile:
		if line[0] == '#':
			continue
		if " warnung! " in line:
			continue
		data = line.split()
		replaceExpChar(data[0])
		freqs.append(float(replaceExpChar(data[0])))
		i_0s.append(float(replaceExpChar(data[1])))
		i_1s.append(float(replaceExpChar(data[2])))
		i_s.append(float(replaceExpChar(data[3])))		
			
	
	intensity = i_s
	if intensityparameter=="i_0": intensity = i_0s
	if intensityparameter=="i_1": intensity = i_1s
		
	c = 3.E8
	angstroms = 1E-10
	wavelengths = [c / f for f in freqs]
	flambdas = []
	for fnu,l in zip(intensity, wavelengths):
		flambda = fnu * c / (l*l)
		flambdas.append(flambda) 
		# print l, flambda
	wavelengths = [w / angstroms for w in wavelengths] 
	lowerWavelength = min(wavelengths)
	upperWavelength = max(wavelengths)
	lowerFlambda = min(flambdas)
	upperFlambda = max(flambdas)
	
	spectrum = spectrumClasses.spectrumObject()
	spectrum.setData(wavelengths, flambdas)
	spectrum.name = "model: angle=%f, temp=%f, field=%f"%(angle, temperature, field)
	spectrum.sortData()
	
	return spectrum
	
def getSampledModel(wavelengths, angle, field, temperature, log_lambda):
	global colour
	geometry = 0
	model = getModelSpectrum(angle, field, temperature, log_lambda, geometry)
	model.resample(wavelengths)
	modelArea = model.integrate()
	model.divide(modelArea)
	model.divide(1/observedArea)
	return model.flux
	
def quadratic(x, a0, a1, a2):
	y = a0*x*x + a1 *x + a2
	return y
	
def computeChiSq(spectrum, model):
	chi = 0
	for f1, f2 in zip(spectrum.flux, model):
		chi+= (f1-f2)**2 
	return chi
	
def getChiSqByParameters(params, *args):
	global iteration, mainPlotWindow, currentPlotWindow, colour, inclination, phase
	print "Params:", params
	beta = params[0]
	log_lambda = params[1]
	scale_factor = params[2]
	linear_offset = params[3]
	print "Args:", args
	temperature = args[0]
	field = args[1]
	
	# cos(theta) = cos(i)cos(beta) - sin(i)sin(beta)cos(phi + pi/2)
	cosTheta = math.cos(radians(inclination)) * math.cos(radians(beta)) - math.sin(radians(inclination)) * math.sin(radians(beta))*math.cos(phase + math.pi/2.)
	angle = math.acos(cosTheta) / math.pi * 180
	print "Angle: %f [deg], Field: %f [MG], Temperature:%f [keV], log_lambda: %f, scale: %f, offset: %f"%(angle, field, temperature, log_lambda, scale_factor, linear_offset)
	model = getSampledModel(observedSpectrum.wavelengths, angle, field, temperature, log_lambda)
	model = [m * scale_factor + linear_offset for m in model]
	chi = computeChiSq(observedSpectrum, model)
	allChiSqs.append(chi)
	print "Chi-squared:", chi
	startWavelength = min(observedSpectrum.wavelengths)
	endWavelength = max(observedSpectrum.wavelengths)
	
	# Draw the most recent iteration
	ppgplot.pgslct(currentPlotWindow)
	ppgplot.pgsci(1)
	ppgplot.pgenv(startWavelength, endWavelength, lowerFlux, upperFlux, 0, 0)
	ppgplot.pgline(observedSpectrum.wavelengths, observedSpectrum.flux)
	ppgplot.pgsci(4)
	ppgplot.pgline(observedSpectrum.wavelengths, model)
	ppgplot.pgsci(1)
	ppgplot.pglab("wavelength", "flux", "Current fit: %d"%iteration)
	
	# Overplot the iteration on the original diagram
	print "overplotting"
	ppgplot.pgslct(mainPlotWindow)
	ppgplot.pgsci(colour)
	ppgplot.pgline(observedSpectrum.wavelengths, model)
	colour += 1
	if colour>15: colour = 1
	ppgplot.pgsci(1)
	
	# Re-generate the Chi-Squared plot
	ppgplot.pgslct(chiSqPlotWindow)
	
	if iteration > 9:
		ppgplot.pgenv(0, iteration+1, 0, max(allChiSqs), 0, 0)
	else:
		ppgplot.pgenv(0, 10, 0, max(allChiSqs), 0, 0)
	iterations = range(iteration+1)
	ppgplot.pgpt(iterations, allChiSqs, 2)
	minCh = min(allChiSqs)
	medCh = numpy.median(allChiSqs)
	maxCh = max(allChiSqs)
	ppgplot.pglab("Iteration [n]", "Chi-squared", "Chi-squared values [%.2f, %.2f, %.2f]"%(minCh, medCh, maxCh))
	
	
	iteration += 1
	return chi
	
def computeViewingAngle(phase, inclination, beta, phi):
	phaseAngle = phase * 360
	cosTheta = math.cos(radians(inclination)) * math.cos(radians(beta)) - math.sin(radians(inclination)) * math.sin(radians(beta))*math.cos(radians(phaseAngle - phi) + math.pi/2.)
	theta = math.acos(cosTheta) / math.pi * 180
	return theta
	
if __name__ == "__main__":
	iteration = 0
	parser = argparse.ArgumentParser(description='Runs the FORTRAN code to make a model spectrum and then plot it.')
	 
	arg = parser.parse_args()
	print arg
	
	beta = 40        	# beta (angle of the magnetic pole to spin axis)
	field = 30   		# Field in MGauss
	temperature = 20 	# Temperature in keV
	log_lambda = 1 		# Optical depth parameter
	geometry = 0		# Mysterious geometry parameter
	
	inclination =81    	# Inclination of orbit
	phi = 20				# Angle offset of magnetic axis to orbital axis (not used yet)
	
	phaseArray = []
	angleArray = []
	# Plot theta as a function of orbital phase for these parameters
	for phase in numpy.arange(0.5, 1.51, 0.01):
		theta = computeViewingAngle(phase, inclination, beta, phi)
		phaseArray.append(phase)
		angleArray.append(theta)
		
	mainPlotWindow = ppgplot.pgopen("/xs")	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgsci(1)
	ppgplot.pgenv(min(phaseArray), max(phaseArray), 0, 180, 0, 0)
	ppgplot.pgline(phaseArray, angleArray)
	ppgplot.pgsls(2)
	ppgplot.pgline([0.5, 1.5], [90, 90])
	ppgplot.pglab("orbital phase", "viewing angle", "Viewing angle \gh as a function of orbital phase.")
	ppgplot.pgtext(0.6, 150, "i:%2.0f \gb:%2.0f \gf:%2.0f"%(inclination, beta, phi))
		
	
	modelPlotWindow = ppgplot.pgopen("models_i_81_b_40.ps/ps")	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgask(False)
	mainFluxMax = 0
	mainFluxMin = 1E99
	
	for phase in numpy.arange(0.5, 1.6, 0.1):
		
		ppgplot.pgsci(1)
		theta = computeViewingAngle(phase, inclination, beta, phi)
		intensities = ["i", "i_0", "i_1"]
		lineStyles = [1, 2, 3]
	
		index = 0
		for i, l in zip(intensities, lineStyles): 
			modelSpectrum = getModelSpectrum(theta, field, temperature, log_lambda, geometry, i)
			requiredWavelengths = numpy.arange(5700, 8900, 1)
			modelSpectrum.resample(requiredWavelengths)
			fluxArray = modelSpectrum.getFlux()
			wavelengthArray = modelSpectrum.getWavelengths()	
			minWavelength = numpy.min(wavelengthArray)
			maxWavelength = numpy.max(wavelengthArray)
			minFlux = numpy.min(fluxArray)
			maxFlux = numpy.max(fluxArray)
			
			if i == "i":
				if mainFluxMax<maxFlux: mainFluxMax = maxFlux
				if mainFluxMin>minFlux: mainFluxMin = minFlux
			if index == 0: 
				ppgplot.pgenv(minWavelength, maxWavelength, 0, maxFlux, 0, 0)
				yRange = maxFlux
			ppgplot.pgsls(l)
			ppgplot.pgline(wavelengthArray, fluxArray)
			textX = 7000
			textY = modelSpectrum.getNearestFlux(textX) + (yRange *.02)
			ppgplot.pgtext(textX, textY, i)
		
			index+=1
	
		ppgplot.pgsls(1)
		ppgplot.pglab("wavelength \A", "i", "Sight angle \gh: %2.1f, Phase: %2.2f, i: %2.1f, beta: %2.1f, field: %3.1f, temp: %3.1f"%(theta, phase, inclination, beta, field, temperature))
	
	print "Flux range:", mainFluxMax, mainFluxMin
	# Now plot i(total) across many phases
	modelPlotWindow = ppgplot.pgopen("/xs")	
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgask(False)
	ppgplot.pgenv(minWavelength, maxWavelength, mainFluxMin, mainFluxMax, 0, 0)
	
	for phase in numpy.arange(0.5, 1.6, 0.02):
		print "Phase:", phase
		theta = computeViewingAngle(phase, inclination, beta, phi)
		
		modelSpectrum = getModelSpectrum(theta, field, temperature, log_lambda, geometry, "i")
		requiredWavelengths = numpy.arange(5700, 8900, 1)
		modelSpectrum.resample(requiredWavelengths)
		fluxArray = modelSpectrum.getFlux()
		wavelengthArray = modelSpectrum.getWavelengths()	
		ppgplot.pgsls(1)
		ppgplot.pgline(wavelengthArray, fluxArray)
		textX = 7000
		textY = modelSpectrum.getNearestFlux(textX) + (yRange *.02)
		ppgplot.pgtext(textX, textY, "%2.1f"%phase)
		
	sys.exit()
	
	# Snip out Halpha 
	spectrum.snipWavelengthRange(6550, 6570)
	
	lowerWavelength = min(spectrum.wavelengths)
	upperWavelength = max(spectrum.wavelengths)
	observedSpectrumRange = (lowerWavelength, upperWavelength)
	lowerFlux = min(spectrum.flux)
	upperFlux = max(spectrum.flux)
	lowerFlux = 0
	
	mainPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pgsci(1)
	ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
	ppgplot.pgline(spectrum.wavelengths, spectrum.flux)
	ppgplot.pglab("wavelength", "flux", spectrum.objectName)
	
	observedArea = spectrum.integrate()
	print "Wavelength range of observations:", observedSpectrumRange
		
	# modelPlotWindow = ppgplot.pgopen(arg.device)	
	# pgPlotTransform = [0, 1, 0, 0, 0, 1]	
	# ppgplot.pgask(False)
	
	currentPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]	
	ppgplot.pgslct(currentPlotWindow)
	ppgplot.pgenv(lowerWavelength, upperWavelength, lowerFlux, upperFlux, 0, 0)
	ppgplot.pgline(spectrum.wavelengths, spectrum.flux)
	ppgplot.pglab("wavelength", "flux", "Current fit")
	
	chiSqPlotWindow = ppgplot.pgopen(arg.device)	
	ppgplot.pgask(False)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]	
	ppgplot.pgslct(chiSqPlotWindow)
	ppgplot.pgenv(0, 10, 0, 100, 0, 0)
	ppgplot.pglab("Iteration [n]", "Chi-squared", "Chi-squared values")
	
	
	colour = 2
	allChiSqs = []
	
	observedSpectrum = spectrum
	#         Angle      Temp     Scale factor  Linear offset
	bounds = [(60, 90), (10, 40), (0.5, 1.5),   (0.5, 1.5)]
	#        Field strength
	fixed = (36, 33)

	scale = 1.0
	offset = 0.1
	guess = [angle, log_lambda, scale, offset]
	fixed = (temperature, field)
	iteration = 0
	results = scipy.optimize.minimize(getChiSqByParameters, guess, args = fixed, method = 'Nelder-Mead')
	
	print results
	
