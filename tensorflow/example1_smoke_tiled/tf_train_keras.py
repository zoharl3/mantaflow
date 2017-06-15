#******************************************************************************
#
# MantaFlow fluid solver framework
# Copyright 2017 Daniel Hook, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL) 
# http://www.gnu.org/licenses
#
# Manta & tensor flow example with tiles using Keras
# (note, this is a slightly simplified version of the tf_train script, with a deeper network and better regularization)
#
#******************************************************************************

# NT_DEBUG todo
# test load/save
# test output mode
# test vel...

import time
import os
import shutil
import sys
import math

import tensorflow as tf
import numpy as np

# load manta tools
sys.path.append("../tools")
import tilecreator as tiCr
import uniio
import paramhelpers as ph
from convautoenc import ConvolutionalAutoEncoder

import keras
#import tf.keras as keras  # TODO, later for TF v1.2 

# path to sim data, trained models and output are also saved here
basePath = '../data/'

# main mode switch:
outputOnly = True  # apply model, or run full training?
outputInputs = False							 

simSizeLow   = 64
tileSizeLow  = 16
upRes        = 4

# dont use for training! for applying model, add overlap here if necessary (i.e., cropOverlap>0) 
# note:  cropTileSizeLow + (cropOverlap * 2) = tileSizeLow
cropOverlap     = 0
cropTileSizeLow = tileSizeLow - 2*cropOverlap

emptyTileValue  = 0.01
learningRate    = 0.00005
trainingEpochs  = 10000 # for large values, stop manualy with ctrl-c...
dropout         = 0.9   # slight...
batchSize       = 96
testInterval    = 200
saveInterval    = 1000
#saveInterval    = 50 # increase for long runs...
fromSim = toSim = -1
keepAll         = False
numTests        = 10      # evaluate on 10 data points from test data
randSeed        = 1
fileFormat      = "npz"
brightenOutput  = -1 # multiplied with output to brighten it up
outputDataName  = '' # name of data to be regressed; by default, does nothing (density), e.g. if output data is pressure set to "pressure"
bWidth          = -1 # boundaryWidth to be cut away. 0 means 1 cell, 1 means two cells. In line with "bWidth" in manta scene files																																  

# run this many iterations per keras fit call
kerasChunk = 200

# optional, add velocity (x,y,z) as additional channels to input?
useVelocities   = 0


# ---------------------------------------------

# load an existing model when load_ values > -1
# when training , manually abort when it's good enough
# then enter test_XXXX id and model checkpoint ID below to load

loadModelTest = -1
loadModelNo   = -1
testPathStartNo = 1

# command line params, explanations mostly above with variables
outputOnly      = int(ph.getParam( "out",             outputOnly ))>0
trainingEpochs  = int(ph.getParam( "trainingEpochs",  trainingEpochs ))
loadModelTest   = int(ph.getParam( "loadModelTest",   loadModelTest))
loadModelNo     = int(ph.getParam( "loadModelNo",     loadModelNo))
basePath        =     ph.getParam( "basePath",        basePath        )
useVelocities   = int(ph.getParam( "useVelocities",   useVelocities  ))
testPathStartNo = int(ph.getParam( "testPathStartNo", testPathStartNo  ))
fromSim         = int(ph.getParam( "fromSim",         fromSim  )) 
toSim           = int(ph.getParam( "toSim",           toSim  ))
alwaysSave      = int(ph.getParam( "alwaysSave",      False  )) # by default, only save checkpoint when cost is lower, can be turned off here
randSeed        = int(ph.getParam( "randSeed",        randSeed )) 
simSizeLow      = int(ph.getParam( "simSizeLow",      simSizeLow )) 
upRes           = int(ph.getParam( "upRes",           upRes ))
#fileFormat     =     ph.getParam( "fileFormat",      fileFormat) # create pngs for inputs
outputInputs    = int(ph.getParam( "outInputs",       outputInputs)) 
brightenOutput  = int(ph.getParam( "brightenOutput",  brightenOutput)) 
outputDataName  =    (ph.getParam( "outName",         outputDataName))
ph.checkUnusedParams()

# initialize
simSizeHigh  = simSizeLow   * upRes
tileSizeHigh = tileSizeLow  * upRes
if outputOnly: # dont discard
	emptyTileValue = -1.

if toSim==-1:
	toSim = fromSim
tiCr.setBasePath(basePath)

np.random.seed(randSeed)
tf.set_random_seed(randSeed)

if not outputOnly:
	# run train!
	loadModelTest = -1
	if fromSim==-1:
		fromSim = toSim  = 1000 # short default, use single sim

	if cropOverlap>0:
		print("Error - dont use cropOverlap != 0 for training...")
		exit(1)

else:
	keepAll = True

# ---------------------------------------------

# create output dir
def next_test_path(folder_no = 1):
	test_path_addition = 'test_%04d/' % folder_no
	while os.path.exists(basePath + test_path_addition):
		folder_no += 1
		test_path_addition = 'test_%04d/' % folder_no 
		test_folder_no = folder_no
	test_path = basePath + test_path_addition
	print("Using test dir '%s'" % test_path)
	os.makedirs(test_path)
	return (test_path, folder_no)

# create model loading path
if not loadModelTest == -1:
	if not os.path.exists(basePath + 'test_%04d/' % loadModelTest):
		print('ERROR: Test to load does not exist.')
	# search for newest model if no loadModelNo is given
	if loadModelNo == -1:
		for currModel in range(0, 999):
			if os.path.isfile(basePath + 'test_%04d/model_%04d.kkpt' % (loadModelTest, currModel)):
				loadModelNo = currModel
		if loadModelNo == -1:
			print('ERROR: Model with id below 200 does not exist. Please specify model id as "loadModelNo".')
			exit()
		# print('Latest model: %d.' % loadModelNo)

load_path = basePath + 'test_%04d/model_%04d.kkpt' % (loadModelTest, loadModelNo)

(test_path,test_folder_no) = next_test_path(testPathStartNo)
if not outputOnly: uniio.backupFile(__file__, test_path)

# custom Logger to write Log to file
class Logger(object):
	def __init__(self):
		self.terminal = sys.stdout
		self.log = open(test_path + "logfile.log", "a")

	def write(self, message):
		self.terminal.write(message)
		self.log.write(message)

	def flush(self): 
		# to avoid errormsg, " AttributeError: 'Logger' object has no attribute 'flush' "
		pass
sys.stdout = Logger()

print("Call: " + str(" ".join(sys.argv) ) )

# print Variables to log
def print_variables():
	print('\nUsing variables:')
	print('fromSim: {}'.format(fromSim))
	print('toSim: {}'.format(toSim))
	print('simSizeLow: {}'.format(simSizeLow))
	print('tileSizeLow: {}'.format(tileSizeLow))
	print('cropOverlap: {}'.format(cropOverlap))
	print('cropTileSizeLow: {}'.format(cropTileSizeLow))
	print('upRes: {}'.format(upRes))
	print('emptyTileValue: {}'.format(emptyTileValue))
	print('learningRate: {}'.format(learningRate))
	print('trainingEpochs: {}'.format(trainingEpochs))
	print('dropout: {}'.format(dropout))
	print('batchSize: {}'.format(batchSize))
	print('\n')

print_variables()

n_input  = tileSizeLow  ** 2 
n_output = tileSizeHigh ** 2
n_inputChannels = 1 
if useVelocities:
	n_inputChannels = 4
n_input *= n_inputChannels

clFMs = int(8 / n_inputChannels)
#print( "inputs " + format(len( tiCr.tile_data['inputs_train']) ))

# ---------------------------------------------
# model

model = keras.models.Sequential()
model.add( keras.layers.Conv2D(clFMs/2, (2,2), activation='relu', strides=(2,2), input_shape=(tileSizeLow,tileSizeLow,n_inputChannels), padding='same' ) )

model.add( keras.layers.Conv2D(clFMs  , (2,2), activation='relu', strides=(2,2), padding='same' ) )
model.add( keras.layers.BatchNormalization() )  
model.add( keras.layers.Flatten() ) # not really needed

model.add( keras.layers.Dense(4*4*clFMs, activation='relu') )
model.add( keras.layers.BatchNormalization() )  
model.add( keras.layers.Dropout(0.25) )  
model.add( keras.layers.Reshape( (4,4,clFMs) ) ) # for flatten, not really needed

model.add( keras.layers.convolutional.Conv2DTranspose(clFMs,   (2,2), activation='relu', strides=(2,2), padding='same' ) )
model.add( keras.layers.convolutional.Conv2DTranspose(clFMs/2, (2,2), activation='relu', strides=(2,2), padding='same' ) )
model.add( keras.layers.convolutional.Conv2DTranspose(clFMs/4, (2,2), activation='relu', strides=(2,2), padding='same' ) )
model.add( keras.layers.convolutional.Conv2DTranspose(1,       (4,4), activation='relu', strides=(2,2), padding='same' ) )

model.compile( loss='mse', optimizer='adam') 

if 1: # count DOFs?
	DOFs = 0
	for l in model.layers:
		for i in l.get_weights():
			DOFs += i.size
	print("Total DOFs: " + format(DOFs) )

if not loadModelTest == -1:
	model.load_weights( load_path )
	print("Model restored from %s." % load_path)

# load test data, note no split into train & test/validation set necessary here, done by keras during fit
tiCr.loadTestDataNpz(fromSim, toSim, emptyTileValue, cropTileSizeLow, cropOverlap, 1.0, 0.0, load_vel=useVelocities, low_res_size=simSizeLow, upres=upRes, keepAll=keepAll)

# ---------------------------------------------
# training
if not outputOnly:
	# manually reshape data to remove third spatial dimension
	dataSize = tiCr.tile_data['outputs_train'].shape[0]
	#print( "Data sizes "+ format(dataSize) +", inputs " + format( tiCr.tile_data['inputs_train'].shape) + ", outputs " + format( tiCr.tile_data['outputs_train'].shape) ) # 16x16, 64x64 
	tiCr.tile_data['inputs_train']  = np.reshape( tiCr.tile_data['inputs_train'],  [dataSize,tileSizeLow,tileSizeLow,  n_inputChannels] )
	tiCr.tile_data['outputs_train'] = np.reshape( tiCr.tile_data['outputs_train'], [dataSize,tileSizeHigh,tileSizeHigh,1 ] )

	startTime = time.time()
	if 1:
		lastSave   = 1
		save_no    = 0
		trainingEpochs = int(trainingEpochs/kerasChunk)
		for epoch in range(trainingEpochs):
			batch_xs, batch_ys = tiCr.selectRandomTiles(batchSize*kerasChunk) 
			batch_xs  = np.reshape( batch_xs,  [batchSize*kerasChunk,tileSizeLow,tileSizeLow,  n_inputChannels] )
			batch_ys  = np.reshape( batch_ys,  [batchSize*kerasChunk,tileSizeHigh,tileSizeHigh,  1] )

			# ...go!
			hist = model.fit( batch_xs, batch_ys, batch_size=batchSize, epochs=1, validation_split=0.05 )

			#_, cost, summary = sess.run([optimizer, costFunc, lossTrain], feed_dict={x: batch_xs, y_true: batch_ys, keep_prob: dropout})
			#cost = lastCost - 1.0 # NT_DEBUG fixme, read out cost...

			# save model
			doSave = False
			if (lastSave >= saveInterval) or epoch == (trainingEpochs-1): doSave = True 

			if doSave:
				mpath = test_path + 'model_%04d.kkpt' % save_no
				model.save_weights(mpath)
				print('Saved weights ' + mpath )
				save_no += 1
				lastSave = 1
			else:
				lastSave += 1
			print('\nScript Epoch {:04d}/{:04d}'.format((epoch + 1), trainingEpochs))

	print('\n*****TRAINING %d FINISHED*****' % test_folder_no)
	training_duration = (time.time() - startTime) / 60.0
	#print('Training needed %.02f minutes.' % (training_duration))
	print('To apply the trained model, set "outputOnly" to True, and set id for "loadModelTest". E.g. "out 1 loadModelTest %d".' % test_folder_no)
	
else:
	# ---------------------------------------------
	# outputOnly: apply to a full data set, and re-create full outputs from tiles

	batch_xs, batch_ys = tiCr.tile_inputs_all_complete, tiCr.tile_outputs_all_complete # old, inputs_test, outputs_test
	print('Creating %d tile outputs...' % (len(batch_xs)) )

	tdataSize = len(batch_xs)
	# for prediction, we have to get rid of the third spatial dim
	batch_xs = np.reshape( batch_xs, [tdataSize,tileSizeLow,tileSizeLow,  n_inputChannels] )
	batch_ys = np.reshape( batch_ys, [tdataSize,tileSizeHigh,tileSizeHigh,1 ] )

	resultTiles = model.predict( batch_xs )

	# now restore it...
	batch_xs = np.reshape( batch_xs, [tdataSize,1, tileSizeLow,tileSizeLow,  n_inputChannels] )
	batch_ys = np.reshape( batch_ys, [tdataSize,1, tileSizeHigh,tileSizeHigh,1 ] )
	resultTiles = np.reshape( resultTiles, [tdataSize,1, tileSizeHigh,tileSizeHigh,1 ] )

	# simply concat tiles into images...
	tileSizeHiCrop = upRes * cropTileSizeLow
	tilesPerImg = (simSizeHigh // tileSizeHiCrop) ** 2
	img_count = len(tiCr.tile_inputs_all_complete) / tilesPerImg
	tiCr.debugOutputPngsCrop(resultTiles, tileSizeHigh, simSizeHigh, test_path, \
		imageCounter=1, cut_output_to=tileSizeHiCrop, tiles_in_image=tilesPerImg)

	if outputInputs:
		tiCr.debugOutputPngsSingle(batch_xs,         tileSizeLow, simSizeLow, test_path, imageCounter=0, name='input', channel=0)
		if useVelocities: 
			tiCr.debugOutputPngsSingle(batch_xs,         tileSizeLow, simSizeLow, test_path, imageCounter=0, name='in_vel_x', channel=1)
			tiCr.debugOutputPngsSingle(batch_xs,         tileSizeLow, simSizeLow, test_path, imageCounter=0, name='in_vel_y', channel=2)

	print('Output finished, %d pngs written to %s.' % (img_count, test_path) )



