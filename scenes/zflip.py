
# flip5 in the comments refers to the scene script flip05_nbflip.py

import os, sys, math
import keyboard

from manta import *

sys.path.append( r'c:\prj\python\\' )
from text_color import *
from tictoc import *

# prints line number
import logging
logging.basicConfig(
    format="%(pathname)s line%(lineno)s: %(message)s",
    level=logging.INFO
)
#logging.info('') # example

###

# auto-flush
sys.stdout.reconfigure( line_buffering=True )

out = r'c:/prj-external-libs/mantaflow/out/'

os.system( 'rm %s*.*' % out )
os.system( 'cp %s../video.bat %s' % (out, out) )

# (debug) for consistent result; for large res, the step() hangs?
if 0:
    limit_to_one_core()

# flags
bMesh       = 1
bSaveParts  = 0 # .vdb
bSaveUni    = 0 # .uni
if bSaveParts or bSaveUni:
    bMesh = 1

bScreenShot = 1

# solver params
dim = 3 # 2, 3
part_per_cell_1d = 2 # 3, 2(default), 1
it_max = 1400 # 300, 500, 1200, 1500
res = 150 # 32, 48, 64(default), 96, 128(large), 256(, 512 is too large)

b_fixed_vol = 0
b_correct21 = 1
narrowBand = bool( 0 )
narrowBandWidth = 5 # 32:5, 64:6, 96:6, 128:8

###

if b_correct21:
    b_fixed_vol = 0
    narrowBand = bool( 0 )

combineBandWidth = narrowBandWidth - 1

dt = .2 # .2(default), .5, 1(flip5, easier to debug)
gs = vec3( res, res, res )
if dim == 2:
    gs.z = 1
    bMesh = 0
    bSaveParts = 0
if 0: # sample 1D; requires changing sampleLevelsetWithParticles()
    ppc = part_per_cell_1d
else:
    ppc = part_per_cell_1d**dim

boundary_width = 0

s = Solver( name='main', gridSize=gs, dim=dim )
gravity = -0.1 # -0.1
gravity *= math.sqrt( res )
#gravity = -0.003 # flip5

print()
print( 'dim:', dim, ', res:', res, ', ppc:', ppc )
print( 'narrowBand:', narrowBand, ', narrowBandWidth:', narrowBandWidth )
print( 'b_fixed_vol:', b_fixed_vol )
print( 'gravity: %0.02f' % gravity )
print( 'timestep:', dt )

# adaptive time stepping; from flip5
b_adaptive_time_step = 0
if b_adaptive_time_step:
    s.frameLength = 1.0   # length of one frame (in "world time")
    s.timestep    = 1.0
    s.timestepMin = 0.5   # time step range
    s.timestepMax = 1.0
    s.cfl         = 5.0   # maximal velocity per cell, 0 to use fixed timesteps

# size of particles 
radiusFactor = 1.0

# prepare grids and particles
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
vel2     = s.create(MACGrid)
velOld   = s.create(MACGrid)
velParts = s.create(MACGrid)
pressure = s.create(RealGrid)
mapWeights = s.create(MACGrid)
pp       = s.create(BasicParticleSystem)

# add velocity data to particles
pVel     = pp.create(PdataVec3)
pVel2     = pp.create(PdataVec3)
phiParts = s.create(LevelsetGrid)
phiMesh  = s.create(LevelsetGrid)
phiObs   = s.create(LevelsetGrid)
phiObsInit = s.create(LevelsetGrid)
obsVel  = s.create(MACGrid)
fractions = s.create(MACGrid) # sticks to walls and becomes slow without this?
mesh     = s.create(Mesh)
bfs      = s.create(IntGrid) # discrete distance from surface

# Acceleration data for particle
pindex = s.create(ParticleIndexSystem)
gpi    = s.create(IntGrid)

#position solver stuff
usePositionSolver = True
density = s.create(RealGrid)
Lambda = s.create(RealGrid)
deltaX = s.create(MACGrid)
flagsPos = s.create(FlagGrid)
pMass = pp.create(PdataReal)
mass = 1.0 / part_per_cell_1d**3
if dim == 2:
    mass = 1.0 / part_per_cell_1d**2

resampleParticles = False # must be a boolean type

if resampleParticles:
    gCnt = s.create(IntGrid)
    
# scene setup
flags.initDomain( boundaryWidth=boundary_width, phiWalls=phiObs ) 

if 0: # dam
    # my dam
    #fluidbox = Box( parent=s, p0=gs*( vec3(0, 0, 0.3) ), p1=gs*( vec3(0.4, 0.8, .7) ) )
    fluidbox = Box( parent=s, p0=gs*( vec3(0, 0, 0.35) ), p1=gs*( vec3(0.3, 0.6, .65) ) ) # new dam (smaller, less crazy)

    # flip05_nbflip.py
    #fluidbox = Box( parent=s, p0=gs*vec3(0, 0.15, 0), p1=gs*vec3(0.4, 0.5, 0.8) )

    # square
    if 0:
        t1 = 0.4 # 0.15, 0.3, .4
        sz1 = .1 # .2, .4
        t = vec3(t1, t1, 0)
        sz = vec3(sz1, sz1, 1)
        fluidbox = Box( parent=s, p0=gs*( t + vec3(0,0,0) ), p1=gs*( t + sz ) )

    # manta dam
    #fluidbox = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) 

    # phi
    phi = fluidbox.computeLevelset()
    flags.updateFromLevelset( phi )

elif 0: # falling drop
    fluidBasin = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.1,1.0)) # basin
    dropCenter = vec3(0.5,0.3,0.5)
    dropRadius = 0.1
    fluidDrop  = Sphere( parent=s , center=gs*dropCenter, radius=res*dropRadius)
    fluidVel   = Sphere( parent=s , center=gs*dropCenter, radius=res*(dropRadius+0.05) )
    fluidSetVel= vec3(0,-1,0)
    phi = fluidBasin.computeLevelset()
    phi.join( fluidDrop.computeLevelset() ) # add drop
    flags.updateFromLevelset( phi )

elif 0: # basin
    # water
    fluidbox = Box( parent=s, p0=gs*( vec3(0, 0.5, 0) ), p1=gs*( vec3(1, 0.9, 1) ) )
    phi = fluidbox.computeLevelset()
    flags.updateFromLevelset( phi )

    # obstacle
    if 1:
        #mesh2 = s.create(Mesh) # it renders only one mesh (mLocalMesh)?

        #mesh.load( r'c:\prj\mantaflow_mod\resources\cube1.obj' )
        #mesh.scale( Vec3(1) )

        mesh.load( r'c:\prj\mantaflow_mod\resources\funnel.obj' )
        mesh.scale( Vec3(res) ) # the scale needs to be in all axes (i.e. can't use gs)

        mesh.offset( gs * Vec3(0.5, 0.3, 0.5) )
        meshObs = s.create( LevelsetGrid )
        mesh.computeLevelset( meshObs, 2 ) # uses normals, thus a smooth mesh is better
        #meshObs.printGrid()
        #phiObs.printGrid()
        phiObs.join( meshObs )
        phi.subtract( phiObs )

else: # low, full box
    # water
    h = 0.25 # 0.3, 0.9
    fluidbox = Box( parent=s, p0=gs*( vec3(0, 0., 0) ), p1=gs*( vec3(1, h, 1) ) )
    phi = fluidbox.computeLevelset()
    flags.updateFromLevelset( phi )

    # moving obstacle
    if 1:
        obs_rad = .1
        obs_center = gs*vec3( 0.5, 0.9 - obs_rad, 0.5 )
        sphere = Sphere( parent=s, center=obs_center, radius=res*obs_rad )
        phiObsInit.copyFrom( phiObs )
        phiObs.join( sphere.computeLevelset() )

        obsVelVec = vec3( 0, -0, 0 )
        obsVel.setConst( obsVelVec )
        obsVel.setBound( value=Vec3(0.), boundaryWidth=boundary_width+1 )

#phi.printGrid()
#phiObs.printGrid()

sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=part_per_cell_1d, randomness=0.1 ) # 0.05, 0.1, 0.2

# also sets boundary flags for phiObs (and shrinks it)
if 0:
    setObstacleFlags( flags=flags, phiObs=phiObs )
else:
    updateFractions( flags=flags, phiObs=phiObs, fractions=fractions, boundaryWidth=boundary_width )
    setObstacleFlags( flags=flags, phiObs=phiObs, fractions=fractions )

# phi is influenced by the walls for some reason
# create a level set from particles
gridParticleIndex( parts=pp, flags=flags, indexSys=pindex, index=gpi )
unionParticleLevelset( pp, pindex, flags, gpi, phiParts, radiusFactor )
phi.copyFrom( phiParts )
if narrowBand:
    extrapolateLsSimple( phi=phi, distance=narrowBandWidth+2, inside=True )
else:
    extrapolateLsSimple( phi=phi, distance=4, inside=True ) # 4
#phi.printGrid()

copyFlagsToFlags( flags, flagsPos )

print( '# particles:', pp.pySize() )
V0 = float(pp.pySize()) / ppc

if 1 and GUI:
    gui = Gui()
    for i in range( 0 ):
        gui.nextMeshDisplay() # 0:full, 1:invisible, 2:x-ray
    gui.setRealGridDisplay( 0 )
    gui.setVec3GridDisplay( 0 )
    if 1 and dim == 3: # camera
        gui.setCamPos( 0, 0, -2.2 ) # drop
        gui.setCamRot( 35, -30, 0 )
    if bMesh:
        gui.toggleHideGrids()
    gui.show()
    #gui.pause()
else:
    bScreenShot = 0

it = 0
it2 = 0

if bScreenShot:
    gui.screenshot( out + 'frame_%04d.png' % it ); # slow

if bSaveParts:
    if bSaveUni:
        pressure.save( out + 'ref_parts_0000.uni' )
        pp.save( out + 'parts_%04d.uni' % it )

    objects = [ flags, phi, pp ] # need the 3 of them for volumetric .vdb
    #objects = [ pp ]
    fname = out + 'fluid_data_%04d.vdb' % it
    print( fname )
    save( name=fname, objects=objects ) # error in debug mode "string too long?"

# measure
f_measure = open( out + '_measure.txt', 'w' )

# loop
ret = 0
while 1:
    emphasize( '\n-----------------\n- time: %g(/%d; it2=%d)' % ( it, it_max, it2 ) )
    print( '- n=%d' % pp.pySize() )

    if 1 and ret != 0:
        error( f'Error: ret={ret}' )
        ret = 0
        break

    if not it < it_max:
        break

    tic()

    maxVel = vel.getMaxAbs()
    if not b_adaptive_time_step:
        s.frameLength = dt
        s.timestep = ( 1 - it % 1 ) * dt
    else: # flip5
        s.adaptTimestep( maxVel ) # enable init above

    # map particle velocities to grid
    print( '- mapPartsToMAC' )
    if 1 and narrowBand:
        # combine particles velocities with advected grid velocities; stores mapWeights; saves velOld
        mapPartsToMAC( vel=velParts, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights )
        # extrapolate velocities throughout the liquid region
        extrapolateMACFromWeight( vel=velParts , distance=2, weight=mapWeights )
        #vel.printGrid()
        #velParts.printGrid()
        combineGridVel( vel=velParts, weight=mapWeights, combineVel=vel, phi=phi, narrowBand=combineBandWidth, thresh=0 )
        #limit_velocity( vel, pVel, 256 ) # 64:15, 128:?
        velOld.copyFrom( vel ) # save before forces and pressure; like the end of mapPartsToMAC()
        #print( '>> combine' )
        #vel.printGrid()
    elif 1:
        # map particle velocities to grid
        mapPartsToMAC( vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=mapWeights )
        extrapolateMACFromWeight( vel=vel , distance=2, weight=mapWeights )

    # moving obstacle
    if 1 and it < 700:
        #flags.printGrid()
        if obs_center.y - res*obs_rad > 1: # move center
            print( '- move obstacle' )
            obsVelVec += dt*vec3(0, gravity, 0) # acceleration
            obsVel.setConst( obsVelVec )
            obsVel.setBound( value=Vec3(0.), boundaryWidth=boundary_width+1 )
            obs_center += dt*obsVelVec
        else:
            print( '- obstacle stopped' )
            obsVel.setConst(Vec3(0.))
        #obsVel.printGrid()

        sphere = Sphere( parent=s, center=obs_center, radius=res*obs_rad )
        phiObs.copyFrom( phiObsInit )
        phiObs.join( sphere.computeLevelset() )
        if 0:
            setObstacleFlags( flags=flags, phiObs=phiObs )
        else:
            updateFractions( flags=flags, phiObs=phiObs, fractions=fractions, boundaryWidth=boundary_width )
            setObstacleFlags( flags=flags, phiObs=phiObs, fractions=fractions )
        #flags.printGrid()

    # update flags
    if 1:
        print( '- markFluidCells (update flags)' )
        markFluidCells( parts=pp, flags=flags ) # needed for a moving obstacle
        #markFluidCells( parts=pp, flags=flags, phiObs=phiObs )
        if narrowBand and ( not b_fixed_vol or it == 0 ):
            update_fluid_from_phi( flags=flags, phi=phi, band_width=narrowBandWidth )
        #flags.printGrid()

    #vel.printGrid()
    #flags.printGrid()
    # forces
    if 1:
        print( '- forces' )
        
        # gravity
        if 1:
            bscale = 0 # 1:adaptive to grid size; flip5
            g = gravity*dt/1 # see addGravity(); assuming s.mDt=1
            #g = gravity
            addGravity( flags=flags, vel=vel, gravity=(0, g, 0), scale=bool(bscale) )

        # vortex
        if 0:
            #c = gs/2
            c = Vec3( res/2, 0.3*res, res/2 )
            vortex( pp=pp, dt=dt, c=c, rad=0.1*res, h=0.9*res, pVel=pVel2 )
            mapPartsToMAC( vel=vel2, flags=flags, velOld=vel2, parts=pp, partVel=pVel2 )
            vel.add( vel2 )

    # set solid (walls)
    print( '- setWallBcs' )
    #setWallBcs( flags=flags, vel=vel, fractions=fractions, phiObs=phiObs, obvel=obsVel ) # calls KnSetWallBcsFrac
    setWallBcs( flags=flags, vel=vel, phiObs=phiObs, obvel=obsVel ) # calls KnSetWallBcs
    #setWallBcs( flags=flags, vel=vel, phiObs=phiObs ) # clear velocity from solid
    #vel.printGrid()
    #flags.printGrid()

    # pressure solve
    if 1:
        print( '- pressure' )
        tic()
        solvePressure( flags=flags, vel=vel, pressure=pressure, phi=phi )
        print( '  (pressure) ', end='' )
        toc()
        #vel.printGrid()

    dist = min( int(maxVel*1.25 + 2), 8 ) # res
    print( '- extrapolate MAC Simple (dist=%0.1f)' % dist )
    extrapolateMACSimple( flags=flags, vel=vel, distance=dist, intoObs=True )
    #setWallBcs( flags=flags, vel=vel, fractions=fractions, phiObs=phiObs, obvel=obsVel )
    setWallBcs( flags=flags, vel=vel, phiObs=phiObs, obvel=obsVel )
    #setWallBcs( flags=flags, vel=vel, phiObs=phiObs )
    #flags.printGrid()
    #vel.printGrid()

    print( '- set particles\' pos0' )
    set_particles_pos0( pp=pp )

    # FLIP velocity update
    print( '- FLIP velocity update' )
    alpha = .1 # 0
    flipVelocityUpdate( vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=1 - alpha )
    #vel.printGrid()
    
    # advect
    print( '- advect' )
    # advect particles
    pp.advectInGrid( flags=flags, vel=vel, integrationMode=IntEuler, deleteInObstacle=False, stopInObstacle=False ) # IntEuler, IntRK2, IntRK4
    #pushOutofObs( parts=pp, flags=flags, phiObs=phiObs ) # creates issues for correct21 and fixedVol
    # advect phi; why? the particles should determine phi, which should flow on its own; disabling this creates artifacts in flip5; it makes it worse for fixed_vol
    if not b_fixed_vol:
        advectSemiLagrange( flags=flags, vel=vel, grid=phi, order=1 )
        flags.updateFromLevelset( phi ) 
    # advect grid velocity
    if narrowBand:
        advectSemiLagrange( flags=flags, vel=vel, grid=vel, order=2 )

    # fixed volume (my scheme)
    include_walls = false
    if b_fixed_vol:
        phi.setBoundNeumann( 0 ) # make sure no new particles are placed at outer boundary
        #phi.printGrid()

        pVel.setSource( vel, isMAC=True ) # set source grid for resampling, used in insertBufferedParticles()

        dt_bound = 0
        #dt_bound = dt/4
        #dt_bound = s.timestep/4
        #dt_bound = max( dt_bound, dt/4 )

        s.timestep = fixed_volume_advection( pp=pp, pVel=pVel, flags=flags, dt=s.timestep, dt_bound=dt_bound, dim=dim, ppc=ppc, phi=phi, it=it2, use_band=narrowBand, band_width=narrowBandWidth, bfs=bfs )
        if s.timestep < 0:
            ret = -1
            s.timestep *= -1

        # if using band
        if 0 and narrowBand:
            include_walls = true

    # correct21 (position solver, Thuerey21)
    # The band requires fixing, probably identifying non-band fluid cells as full. In the paper, it's listed as future work.
    if b_correct21:
        print( '- position solver' )
        copyFlagsToFlags(flags, flagsPos)
        mapMassToGrid(flags=flagsPos, density=density, parts=pp, source=pMass, deltaX=deltaX, phiObs=phiObs, dt=s.timestep, particleMass=mass, noDensityClamping=resampleParticles)
        #gui.pause()
        
        # resample particles
        if resampleParticles:
            print( '    - resample particles' )
            gridParticleIndex(parts=pp, indexSys=pindex, flags=flags, index=gpi, counter=gCnt)
            #apicMapPartsToMAC(flags=flags, vel=vel, parts=pp, partVel=pVel, cpx=apic_pCx, cpy=apic_pCy, cpz=apic_pCz, mass=apic_mass)
            resampeOverfullCells(vel=vel, density=density, index=gpi, indexSys=pindex, part=pp, pVel=pVel, dt=s.timestep)
    
        # position solver
        print( '    - solve pressure due to density' )
        solvePressureSystem(rhs=density, vel=vel, pressure=Lambda, flags=flagsPos, cgAccuracy = 1e-3)
        computeDeltaX(deltaX=deltaX, Lambda=Lambda, flags=flagsPos)
        mapMACToPartPositions(flags=flagsPos, deltaX=deltaX, parts=pp, dt=s.timestep)
        
        # print
        if 0:
            #flags.printGrid()
            #flagsPos.printGrid()
            density.printGrid()
            Lambda.printGrid()
            deltaX.printGrid()
            #gui.pause()
    
    # create level set from particles
    if 1:
        gridParticleIndex( parts=pp, flags=flags, indexSys=pindex, index=gpi )
        unionParticleLevelset( pp, pindex, flags, gpi, phiParts, radiusFactor ) 
        if narrowBand:
            # combine level set of particles with grid level set
            phi.addConst( 1. ); # shrink slightly
            phi.join( phiParts )
            extrapolateLsSimple( phi=phi, distance=narrowBandWidth+2, inside=True, include_walls=include_walls, intoObs=True )
        else:
            # overwrite grid level set with level set of particles
            phi.copyFrom( phiParts )
            extrapolateLsSimple( phi=phi, distance=4, inside=True, include_walls=include_walls ) # 4

    # resample particles
    if not b_fixed_vol:
        pVel.setSource( vel, isMAC=True ) # set source grids for resampling, used in adjustNumber
        minParticles = ppc
        maxParticles = 2*minParticles # 2, 1(exacerbates artifact in flip5 dam 128?)
        if narrowBand:
            phi.setBoundNeumann( 0 ) # make sure no particles are placed at outer boundary
            #phi.printGrid()
            # vel is used only to get the parent
            adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=minParticles, maxParticles=maxParticles, phi=phi, narrowBand=narrowBandWidth, exclude=phiObs ) 
        elif 0:
            adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=minParticles, maxParticles=maxParticles, phi=phi, exclude=phiObs ) 

    # mark int for measure
    if not b_fixed_vol:
        flags.mark_surface()

    # measure
    if 1:
        m = measure( pp, pVel, flags, gravity, ppc, V0 )
        f_measure.write( f'{m[0]}\n' )
        f_measure.flush()

    # mesh
    if bMesh:
        phiMesh.copyFrom( phi )
        improvedParticleLevelset( pp, pindex, flags, gpi, phiMesh, radiusFactor, 1, 1 , 0.4, 3.5 ) # creates artifacts in dam flip05 128

        # mesh
        phiMesh.setBound( value=0., boundaryWidth=1 )
        phiMesh.createMesh( mesh )

    # print/write
    if 0:
        #flags.printGrid()
        vel.printGrid()
    
        pp.printParts()
        #pp.writeParticlesText( out + 'flipt_%04d.txt' % it )
    
    print( '(iteration) ', end='' )
    toc()

    # step
    print( '- step (%d)' % it )
    s.step()
    #print( 'after step' )

    it += s.timestep / dt
    it2 += 1
    if 0 or abs( it - round(it) ) < 1e-7:
        it = round( it )

        # screenshot
        if bScreenShot:
            gui.screenshot( out + 'frame_%04d.png' % it ) # slow

        # save particle data
        if bSaveParts:
            if bSaveUni:
                # save particle data for flip03_gen.py surface generation scene
                pp.save( out + 'parts_%04d.uni' % it )

            # pdata fields must be before pp
            objects = [ flags, phi, pp ]
            #objects = [ pp ]
            save( name=out + 'fluid_data_%04d.vdb' % it, objects=objects )
        
# video
if 1:
    os.system( r'c:\prj-external-libs\mantaflow\out\video.bat > nul 2>&1' )

# code not reached if quitting manta (with esc); pausing in run.py instead
# pause
if 0:
    print( '(zflip.py) press a key...' )
    keyboard.read_key()
elif 0:
    print( '(zflip.py) press enter...' )
    input()

