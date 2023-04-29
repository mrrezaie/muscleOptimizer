import opensim as osim
import numpy as np
from time import time
from itertools import product
import matplotlib.pyplot as plt
from scipy.optimize import nnls
import json, os, sys


interval = 10
oneLegOnly = True # if using identical scale factors for both legs
nameRef = 'models/Rajagopal2015_passiveCal_hipAbdMoved.osim'
nameTarg = 'models/s01boot_calibrated.osim'

console = sys.stdout
sys.stdout = open(f'{nameTarg[:-5]}_N{interval}.log', 'w') # print to log file
refJson = nameRef[:-5]+f'_N{interval}.json'
osim.Logger.setLevelString('Off')

t0 = time()

modelRef = osim.Model(nameRef)
modelTarg = osim.Model(nameTarg)
state = modelRef.initSystem()
nameMuscles = [i.getName() for i in modelRef.getMuscles()]

t1 = time()

coordinateMuscle = dict()
for i in modelRef.getCoordinateSet():

	coordinateMuscle[i.getName()] = list()
	if  i.get_locked()==False and i.getMotionType()!=3:
		length0 = [muscle.getLength(state) for muscle in modelRef.getMuscles()]
		r0 = i.getDefaultValue()
		r1 = i.getRangeMin()
		r2 = i.getRangeMax()
		r3 = (r1+r2)/2
		length = list()
		for j in [r1,r2,r3]:
			i.setValue(state, j)
			modelRef.realizePosition(state)
			length.append([muscle.getLength(state) for muscle in modelRef.getMuscles()])

		dl = 1000 * (np.array(length) - length0) # changes in muscle length (mm)
		ok = np.sum(np.abs(dl), axis=0)>1e-1 # sum of absolute difference
		for muscle in np.array(nameMuscles)[ok].tolist(): # list of muscles
			coordinateMuscle[i.getName()].append(muscle)
		
		i.setValue(state, r0) # back to default 
 # e.g 'knee_angle_r': ['bflh_r', 'bfsh_r', 'gaslat_r', 'gasmed_r', 'grac_r', 'recfem_r', 'sart_r', 
 #                      'semimem_r', 'semiten_r', 'tfl_r', 'vasint_r', 'vaslat_r', 'vasmed_r']

t2 = time()
# print(t2-t1)

muscleCoordinate = dict()
for coordinate,ii in coordinateMuscle.items():
	if len(ii)>0:
		for muscle in ii: # each muscle
			if muscle not in muscleCoordinate.keys():
				muscleCoordinate[muscle] = list()
			muscleCoordinate[muscle].append(coordinate)
# e.g. 'gaslat_l': ['knee_angle_l', 'ankle_angle_l', 'subtalar_angle_l']

if oneLegOnly:
	for i in nameMuscles:
		if i.endswith('_l'):
			muscleCoordinate.pop(i, None)

sharedCoordinates = dict()
for i,ii in muscleCoordinate.items():
	ii = ';'.join(ii)
	if ii not in sharedCoordinates.keys():
		sharedCoordinates[ii] = list()
	sharedCoordinates[ii].append(i)
# e.g. 'knee_angle_l;ankle_angle_l;subtalar_angle_l': ['gaslat_l', 'gasmed_l']

len(muscleCoordinate.keys())
len(sharedCoordinates.keys())

coordinateInterval = dict()
for i,ii in coordinateMuscle.items():
	if len(ii) > 0:
		coordinate = modelRef.getCoordinateSet().get(i)
		coordinateInterval[i] = np.linspace(coordinate.getRangeMin(), coordinate.getRangeMax(), interval).tolist()
# e.g. 'knee_angle_r': [0.0, 0.5236, 1.0472, 1.5708, 2.0944]

def getMuscleQuantities(modelFile):
	model = osim.Model(modelFile)
	state = model.initSystem()

	muscleQuantities = dict()
	for i,ii in muscleCoordinate.items():
		nPose = interval**len(ii)
		muscleQuantities[i] = np.empty((nPose,4))

	for i,ii in sharedCoordinates.items():
		coordinates = i.split(';')
		# e.g. ['hip_flexion_r', 'hip_adduction_r', 'hip_rotation_r']
		print('\tCoordinates:',', '.join(coordinates))
		print('\t\tN muscles:', len(ii))

		temp = [coordinateInterval[j] for j in coordinates]
		combine = list(product(*temp)) # all combinations
		for jj,j in enumerate(combine):
			for kk,k in enumerate(coordinates):
				model.updCoordinateSet().get(k).setValue(state, j[kk])

			model.realizePosition(state)
			# model.equilibrateMuscles(state)

			for nameMuscle in ii: # muscles sharing those coordinates
				# print(muscle)
				muscle = model.getMuscles().get(nameMuscle)
				muscle.setActivation(state, 1)
				muscle.computeEquilibrium(state)
				quantities = np.array([ muscle.getLength(state), \
										muscle.getNormalizedFiberLength(state), \
										muscle.getTendonLength(state), \
										muscle.getPennationAngle(state)])
				muscleQuantities[nameMuscle][jj,:] = quantities

		for k in coordinates: # back to defaults
			default = model.updCoordinateSet().get(k).getDefaultValue()
			model.updCoordinateSet().get(k).setValue(state, default)
		# model.realizePosition(state)
	return muscleQuantities


if os.path.isfile(refJson):
	print(f'Pre-calculated muscle quantities exists:\n\tload {refJson}')
	ref = json.load(open(refJson, mode='r'))
	printJson = False
else:
	print(f'\nRef model: {nameRef}')
	ref = getMuscleQuantities(modelRef)
	printJson = True

print(f'\nTarg model: {nameTarg}')
Targ = getMuscleQuantities(modelTarg)


optimized = dict()
for i in muscleCoordinate.keys():
	print('\n',i)
	muscleRef = modelRef.getMuscles().get(i)
	OFL = muscleRef.getOptimalFiberLength()
	TSL = muscleRef.getTendonSlackLength()
	PAO = muscleRef.getPennationAngleAtOptimalFiberLength()

	ref[i] = np.array(ref[i])
	row = ref[i].shape[0]
	# What is this?
	limit = np.sin(PAO) / np.sin(np.arccos(0.1))
	if limit < 0.5: limit = 0.5
	ok  = ref[i][:,1] > limit

	MTL = ref[i][:,0][ok]   # muscle-tendon length
	NFL = ref[i][:,1][ok]   # normalized filber length
	TL  = ref[i][:,2][ok]   # tendon length
	NTL = TL / TSL          # normalized tendon length
	PA  = ref[i][:,3][ok]   # pennation angle
	NFLT= NFL * np.cos(PA)  # normalized fiber length along tendon

	MTL2 = Targ[i][:,0][ok] # muscle-tendon length of scaled model

	# Compute least-squares solution to equation Ax = b
	A = np.vstack((NFLT,NTL)).T
	b = MTL2#.reshape((-1,1))
	x = nnls(A,b)[0]
	if min(x)<=0:
		if max(TL) - min(TL) < 1e-3:
			print('\tWARNING: no change in tendon length;\n\tRECOMPUTE ...')
		fraction = TL / MTL
		proportion = fraction * MTL2 # tendon-muscletendon length proportion

		A1 = np.vstack((NFLT, 0*NTL)).T # normalized fiber length
		b1 = MTL2 - proportion          # actual fiber length
		x1 = nnls(A1,b1)[0][0]          # optimal fiber length

		A2 = np.vstack((0*NFLT, NTL)).T # normalized tendon length
		b2 = MTL2 - NFLT*x1             # actual tendon length
		x2 = nnls(A2,b2)[0][1]          # tendon slack length

		x = np.array([x1,x2])

	error = np.sum((A.dot(x)-b)**2) # sum of squared error

	print(f'\tOFL   : {OFL:.6f};\t{(100*(OFL-x[0])/OFL):.4f} %')
	print(f'\tTSL   : {TSL:.6f};\t{(100*(TSL-x[1])/TSL):.4f} %')
	print(f'\teval  : {sum(ok)} / {row}')
	print(f'\terror : {(error):.4e}')

	optimized[i] = x.tolist()

for i,ii in optimized.items():
	muscleTarg = modelTarg.getMuscles().get(i)
	muscleTarg.setOptimalFiberLength(ii[0])
	muscleTarg.setTendonSlackLength( ii[1])

if oneLegOnly:
	for i,ii in optimized.items():
		other = i[:-2]+'_l'
		muscleTarg = modelTarg.getMuscles().get(other)
		muscleTarg.setOptimalFiberLength(ii[0])
		muscleTarg.setTendonSlackLength( ii[1])
		ref[other] = ref[i]

if printJson:
	class myEncoder(json.JSONEncoder):
		'''Support numpy
		round numbers'''
		def default(self, obj):
			if isinstance(obj, np.ndarray):
				return obj.round(6).tolist()
			return json.JSONEncoder.default(self, obj)

	json.dump(ref, open(refJson,'w'), cls=myEncoder, separators=(',', ':'))

f.close()
modelTarg.printToXML(f'{nameTarg[:-5]}_N{interval}_pyOpt.osim') # nameTarg[:-5]+
t3 = time()

print(f'\n{t3-t0}')
