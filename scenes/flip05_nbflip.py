##############################################################
# Flip scene with narrow-band FLIP
#
# Reference:
#   [1] Narrow Band FLIP for Liquid Simulations
#       F. Ferstl, R. Ando, R. Westermann, C. Wojtan, N. Thuerey
#       Computer Graphics Forum (Proc. Eurographics 2016)
############################################################## 

import os, sys
from manta import *

out = r'c:/prj-external-libs/mantaflow/out/'

#os.system( 'rm %s*.png' % out )
os.system( 'rm %s*.vdb' % out )

# Toggle between regular FLIP and NB-FLIP
narrowBand = 1

# Choose dimension (2 or 3) and resolution
dim = 3
res = 128

# Configuration 
narrowBandWidth  = 4;                   # Narrow band width in cells  (= R in [1])
combineBandWidth = narrowBandWidth - 1; # Combine band width in cells (= r in [1])

# Solver params
gs = vec3(res,res,res)
if dim==2: gs.z = 1

s = Solver(name='main', gridSize = gs, dim=dim)
mantaMsg("Narrow band: %i" % narrowBand)
mantaMsg('Solver grid resolution: %i x %i x %i' % (gs.x, gs.y, gs.z))

# Adaptive time stepping
s.frameLength = 1.0   # length of one frame (in "world time")
s.timestep    = 1.0
s.timestepMin = 0.5   # time step range
s.timestepMax = 1.0
s.cfl         = 5.0   # maximal velocity per cell, 0 to use fixed timesteps

gravity = (0,-0.003,0)

minParticles = pow(2,dim)

# Prepare grids and particles
flags    = s.create(FlagGrid)

phiParts = s.create(LevelsetGrid) 
phi      = s.create(LevelsetGrid) 
pressure = s.create(RealGrid)

vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
velParts = s.create(MACGrid)
mapWeights  = s.create(MACGrid)

pp        = s.create(BasicParticleSystem) 
pVel      = pp.create(PdataVec3)
mesh      = s.create(Mesh)

# Acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem)
gpi    = s.create(IntGrid)

# Geometry in world units (to be converted to grid space upon init)
flags.initDomain(boundaryWidth=0)
phi.initFromFlags(flags)

# A simple breaking dam
#fluidBasin = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.15,1.0))
#phi.join( fluidBasin.computeLevelset() )

fluidDam = Box( parent=s, p0=gs*( vec3(0,0,0.3) ), p1=gs*( vec3(0.4,0.8,.7) ) ) # my dam
#fluidDam = Box( parent=s, p0=gs*vec3(0,0.15,0), p1=gs*vec3(0.4,0.5,0.8) )
phi.join( fluidDam.computeLevelset() )
    
flags.updateFromLevelset(phi)

sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.1 ) # .1
mapGridToPartsVec3(source=vel, parts=pp, target=pVel )

if 1 and (GUI):
    gui = Gui()
    gui.setRealGridDisplay( 0 )
    gui.setVec3GridDisplay( 0 )
    gui.show( dim==3 )
    #gui.pause()
        
# Main loop
step = -1

if 1:
    # note: when saving pdata fields, they must be accompanied by and listed before their parent pp
    #objects = [flags, phiParts, phi, pressure, vel, pVel, pp]
    objects = [ phi, pp ]
    save( name=out + 'fluid_data_%04d.vdb' % s.frame, objects=objects )

while s.frame < 20000:
    step = step + 1
    
    maxVel = vel.getMax()
    s.adaptTimestep( maxVel )
    mantaMsg( '\nFrame %i, step %i, time-step size %f' % (s.frame, step, s.timestep) )

    # Make sure we have velocities throughout the liquid region
    if narrowBand:
        # Combine particles velocities with advected grid velocities
        mapPartsToMAC(vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights)
        extrapolateMACFromWeight( vel=velParts , distance=2, weight=mapWeights )
        combineGridVel(vel=velParts, weight=mapWeights , combineVel=vel, phi=phi, narrowBand=combineBandWidth, thresh=0)
        velOld.copyFrom(vel)
    else:
        # Map particle velocities to grid
        mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights)
        extrapolateMACFromWeight( vel=vel , distance=2, weight=mapWeights ) 

    # Forces & pressure solve
    #vel.printGrid()
    #flags.printGrid()
    addGravity(flags=flags, vel=vel, gravity=gravity)
    #vel.printGrid()
    setWallBcs(flags=flags, vel=vel)
    solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
    setWallBcs(flags=flags, vel=vel)

    # Extrapolate velocities
    extrapolateMACSimple( flags=flags, vel=vel, distance=(int(maxVel*1.25 + 2.)) )
    
    # Update particle velocities
    flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.95 )

    # Advect particles and grid phi
    # Note: Grid velocities are extrapolated at the end of each step
    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) 
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=1)
    flags.updateFromLevelset(phi)
    
    # Advect grid velocity
    if narrowBand:
        advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)

    # Create level set of particles
    gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
    unionParticleLevelset( pp, pindex, flags, gpi, phiParts )

    if narrowBand:
        # Combine level set of particles with grid level set
        phi.addConst(1.); # shrink slightly
        phi.join( phiParts )
        extrapolateLsSimple(phi=phi, distance=narrowBandWidth+2, inside=True )
    else:
        # Overwrite grid level set with level set of particles
        phi.copyFrom( phiParts )
        extrapolateLsSimple(phi=phi, distance=4, inside=True )

    #extrapolateLsSimple( phi=phi, distance=3 ) # zl why again?
    flags.updateFromLevelset(phi)

    if dim==3:
        phi.createMesh(mesh)
        
    # Resample particles
    pVel.setSource( vel, isMAC=True ) # Set source grids for resampling, used in adjustNumber!
    if narrowBand:
        phi.setBoundNeumann(0) # make sure no particles are placed at outer boundary
        adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, narrowBand=narrowBandWidth ) 
    elif 0:
        adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi ) 

    s.step()

    # optionally save some of the simulation objects to an OpenVDB file (requires compilation with -DOPENVDB=1)
    if 0:
        # note: when saving pdata fields, they must be accompanied by and listed before their parent pp
        objects = [flags, phiParts, phi, pressure, vel, pVel, pp]
        save( name=out + 'fluid_data_%04d.vdb' % s.frame, objects=objects )
