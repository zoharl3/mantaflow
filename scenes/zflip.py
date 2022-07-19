#
# Very simple flip without level set
# and without any particle resampling
# 
import os, sys
import keyboard

from manta import *

sys.path.append( r'c:\prj\python\\' )
from text_color import *

out = r'c:/prj-external-libs/mantaflow/out/'

os.system( 'rm %s*.png' % out )
os.system( 'rm %s*.txt' % out )

# solver params
dim = 2
particleNumber = 3
res = 25
gs = vec3(res,res,res)
gs.z=1
s = Solver( name='main', gridSize = gs, dim=dim )
s.timestep = .2 # .2, .5, 1(easier to debug)
gravity = -10 * 1e-2; # 1e-2, 1e-3; adaptive

print( '(unscaled) gravity:', gravity )
print( 'timestep:', s.timestep )

# prepare grids and particles
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
pp       = s.create(BasicParticleSystem) 
# add velocity data to particles
pVel     = pp.create(PdataVec3) 

# scene setup
flags.initDomain(boundaryWidth=0) 

t = vec3(0.15, 0.15,0)
#fluidbox = Box( parent=s, p0=gs*( t + vec3(0,0,0) ), p1=gs*( t + vec3(0.4,0.7,1) ) ) # my dam

fluidbox = Box( parent=s, p0=gs*( vec3(0,0,0) ), p1=gs*( vec3(0.4,0.8,1) ) ) # my dam

t = vec3(0.2, 0.4,0)
#fluidbox = Box( parent=s, p0=gs*( t+vec3(0,0,0) ), p1=gs*( t+vec3(0.2,0.2,1) ) ) # square

#fluidbox = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) # manta dam

phiInit = fluidbox.computeLevelset()
flags.updateFromLevelset(phiInit)
# phiInit is not needed from now on!

# note, there's no resamplig here, so we need _LOTS_ of particles...
sampleFlagsWithParticles( flags=flags, parts=pp, discretization=particleNumber, randomness=0.2 ) # .2
    
if (GUI):
    gui = Gui()
    gui.setRealGrid( 0 )
    gui.show()
    gui.pause()
    
#main loop
for t in range( 1, int( 1e3 +1) ): # 2500
    emphasize( '- t=%d' % t );
    mantaMsg('\n(Frame %i), simulation time %f' % (s.frame, s.timeTotal))
    
    print( 'mapPartsToMAC' )
    mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
    
    #vel.printGrid()
    
    markFluidCells( parts=pp, flags=flags )

    # adaptive
    print( 'forces' )
    addGravity(flags=flags, vel=vel, gravity=(0,gravity,0))

    #vel.printGrid()

    # set solid
    setWallBcs(flags=flags, vel=vel) # clears velocity from solid
    print( 'setWallBcs' )
    
    # pressure solve
    print( 'pressure' )
    solvePressure(flags=flags, vel=vel, pressure=pressure)

    #vel.printGrid()

    # we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
    extrapolateMACSimple( flags=flags, vel=vel )
    
    # FLIP velocity update
    print( 'FLIP velocity update' )
    flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=1- 0 )
    
    #vel.printGrid()
    
    # advect particles 
    print( 'advectInGrid' )
    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) # IntEuler, IntRK4

    # position solver
    if 0:
        copyFlagsToFlags(flags, flagsPos)
        mapMassToGrid(flags=flagsPos, density=density, parts=pp, source=pMass, deltaX=deltaX, phiObs=phiObs, dt=s.timestep, particleMass=mass, noDensityClamping =  resampleParticles)          
        
        # resample particles
        if 0:
            gridParticleIndex(parts=pp, indexSys=pindex, flags=flags, index=gpi, counter=gCnt)
            apicMapPartsToMAC(flags=flags, vel=vel, parts=pp, partVel=pVel, cpx=apic_pCx, cpy=apic_pCy, cpz=apic_pCz, mass=apic_mass)
            resampeOverfullCells(vel=vel, density=density, index=gpi, indexSys=pindex, part=pp, pVel=pVel, dt=s.timestep)
    
        # position solver
        solvePressureSystem(rhs=density, vel=vel, pressure=Lambda, flags=flagsPos, cgAccuracy = 1e-3)
        computeDeltaX(deltaX=deltaX, Lambda=Lambda, flags=flagsPos)
        mapMACToPartPositions(flags=flagsPos, deltaX=deltaX, parts=pp, dt=s.timestep)
    
    if 0:
        flags.printGrid()
        vel.printGrid()
    
        #pp.printParts()
        pp.writeParticlesText( out + 'flipt_%04d.txt' % t )
    
    #gui.screenshot( out + 'flipt_%04d.png' % t ); # slow
    
    s.step()
        
print( 'press enter...' )
#keyboard.read_key()
input()

