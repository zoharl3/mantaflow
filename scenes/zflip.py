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

os.system( 'rm %s*' % out )

# solver params
dim = 2
particleNumber = 2
res = 5
gs = vec3(res,res,res)
gs.z=1
s = Solver( name='main', gridSize = gs, dim=dim )
s.timestep = .5
gravity = -10 * 1e-2; # adaptive

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

t = vec3(0,.4,0)
#fluidbox = Box( parent=s, p0=gs*( t + vec3(0,0,0) ), p1=gs*( t + vec3(0.4,0.6,1) ) ) # breaking dam
fluidbox = Box( parent=s, p0=gs*( t+vec3(0,0,0) ), p1=gs*( t+vec3(0.4,0.2,1) ) )

phiInit = fluidbox.computeLevelset()
flags.updateFromLevelset(phiInit)
# phiInit is not needed from now on!

# note, there's no resamplig here, so we need _LOTS_ of particles...
#sampleFlagsWithParticles( flags=flags, parts=pp, discretization=particleNumber, randomness=0.2 ) # .2
sampleFlagsWithParticles( flags=flags, parts=pp, discretization=particleNumber, randomness=0. )

    
if (GUI):
    gui = Gui()
    gui.show()
    #gui.pause()
    
#main loop
for t in range(3): # 2500
    emphasize( ' - t=%d' % t );
    mantaMsg('\n(Frame %i), simulation time %f' % (s.frame, s.timeTotal))
    
    pp.writeParticlesText( out + 'flipt_%04d.txt' % t )
    
    # FLIP 
    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntEuler, deleteInObstacle=False ) 
    
    mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
    
    #extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) 
    #extrapolateMACSimple( flags=flags, vel=vel )
    
    markFluidCells( parts=pp, flags=flags )

    # adaptive
    print( 'forces' )
    addGravity(flags=flags, vel=vel, gravity=(0,gravity,0))

    print( 'vel (x=v, y=u):' )
    vel.printGrid()

    # pressure solve
    print( 'pressure' )
    setWallBcs(flags=flags, vel=vel)    
    #solvePressure(flags=flags, vel=vel, pressure=pressure)

    print( 'vel (x=v, y=u):' )
    vel.printGrid()

    # we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
    #extrapolateMACSimple( flags=flags, vel=vel )
    
    # FLIP velocity update
    print( 'FLIP velocity update' )
    flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0 )
    
    print( 'vel (x=v, y=u):' )
    vel.printGrid()
    
    gui.screenshot( out + 'flipt_%04d.png' % t ); # slow
    
    s.step() 

print( 'press enter...' )
keyboard.read_key()

