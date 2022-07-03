#
# Very simple flip without level set
# and without any particle resampling
# 
from manta import *

# solver params
dim = 2
particleNumber = 3
sres=1;
res = 16*sres
gs = vec3(res,res,res)
gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.5

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

t = vec3(0,.2,0)
fluidbox = Box( parent=s, p0=gs*( t + vec3(0,0,0) ), p1=gs*( t + vec3(0.4,0.6,1) ) ) # breaking dam

phiInit = fluidbox.computeLevelset()
flags.updateFromLevelset(phiInit)
# phiInit is not needed from now on!

# note, there's no resamplig here, so we need _LOTS_ of particles...
sampleFlagsWithParticles( flags=flags, parts=pp, discretization=particleNumber, randomness=0.2 )

    
if (GUI):
    gui = Gui()
    gui.show()
    gui.pause()
    
#main loop
for t in range(2500):
    mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))
    
    # FLIP 
    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntEuler, deleteInObstacle=False ) 
    
    mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
    
    #extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) 
    #extrapolateMACSimple( flags=flags, vel=vel )
    
    markFluidCells( parts=pp, flags=flags )

    # resolution dependent; the grid size isn't normalized
    addGravity(flags=flags, vel=vel, gravity=(0,-0.002*sres,0))

    # pressure solve
    setWallBcs(flags=flags, vel=vel)    
    solvePressure(flags=flags, vel=vel, pressure=pressure)

    # we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
    #extrapolateMACSimple( flags=flags, vel=vel )
    
    # FLIP velocity update
    flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0 )
    
    #gui.screenshot( r'c:\prj-external-libs\mantaflow\out\flipt_%04d.png' % t ); # slow
    
    s.step() 

