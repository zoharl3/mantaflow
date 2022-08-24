#
# Very simple flip without level set
# and without any particle resampling
# 
import os, sys, math
import keyboard, copy

from tictoc import *

from manta import *

sys.path.append( r'c:\prj\python\\' )
from text_color import *

out = r'c:/prj-external-libs/mantaflow/out/'

os.system( 'rm %s*.png' % out )
os.system( 'rm %s*.txt' % out )
os.system( 'rm %s*.uni' % out )
os.system( 'rm %s*.vdb' % out )

# flags
bSaveParts  = 1 # needed from drawing the surface
bSaveUni    = 0

bScreenShot = 1

# solver params
dim = 3 # 2, 3
it_max = 280 # 300, 500, 1200, 1500
part_per_cell_1d = 2 # 3, 2(default), 1
res = 64 # 17(min band), 32, 48, 64(default), 128(large)

dt = .2 # .2, .5, 1(easier to debug)
gs = vec3(res, res, res)
if dim == 2:
    gs.z = 1
    bSaveParts = 0
    
s = Solver( name='main', gridSize=gs, dim=dim )
gravity = -0.1
gravity *= math.sqrt( res )

print( 'gravity:', gravity )
print( 'timestep:', dt )

# size of particles 
radiusFactor = 1.0

# prepare grids and particles
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
pp       = s.create(BasicParticleSystem) 
# add velocity data to particles
pVel     = pp.create(PdataVec3) 
phiObs   = s.create(LevelsetGrid, name='phiObs')
mesh     = s.create(Mesh)

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
mass = 1.0 / (part_per_cell_1d * part_per_cell_1d * part_per_cell_1d) 
if dim == 2:
    mass = 1.0 / (part_per_cell_1d * part_per_cell_1d) 

resampleParticles = False # must be a boolean type

if resampleParticles:
    pindex = s.create(ParticleIndexSystem) 
    gpi = s.create(IntGrid)
    gCnt = s.create(IntGrid)
    
# scene setup
flags.initDomain(boundaryWidth=0) 

if 1: # breaking dam
    # my dam
    fluidbox = Box( parent=s, p0=gs*( vec3(0,0,0.3) ), p1=gs*( vec3(0.4,0.8,.7) ) ) 

    # square
    if 0:
        t1 = 0.4 # 0.15, 0.3, .4
        sz1 = .4 # .2, .4
        t = vec3(t1, t1, 0)
        sz = vec3(sz1, sz1, 1)
        fluidbox = Box( parent=s, p0=gs*( t + vec3(0,0,0) ), p1=gs*( t + sz ) )

    # manta dam
    #fluidbox = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) 

    # phi
    phi = fluidbox.computeLevelset()

else: # falling drop
    fluidBasin = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.1,1.0)) # basin
    dropCenter = vec3(0.5,0.3,0.5)
    dropRadius = 0.1
    fluidDrop  = Sphere( parent=s , center=gs*dropCenter, radius=res*dropRadius)
    fluidVel   = Sphere( parent=s , center=gs*dropCenter, radius=res*(dropRadius+0.05) )
    fluidSetVel= vec3(0,-1,0)
    phi = fluidBasin.computeLevelset()
    phi.join( fluidDrop.computeLevelset() )

flags.updateFromLevelset( phi )

sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=part_per_cell_1d, randomness=0.05 ) # 0.05, 0.2
    
copyFlagsToFlags( flags, flagsPos )
flags.initDomain(boundaryWidth=0, phiWalls=phiObs)

np = pp.pySize()
print( '# particles:', np )
pos1 = s.create(PdataVec3)
pos1.pyResize( np )

if 1 and GUI:
    gui = Gui()
    gui.setRealGridDisplay( 0 )
    gui.setVec3GridDisplay( 0 )
    gui.show()
    #gui.pause()
else:
    bScreenShot = 0

it = 0

if bScreenShot:
    gui.screenshot( out + 'frame_%04d.png' % it ); # slow

if bSaveParts:
    if bSaveUni:
        pressure.save( out + 'ref_parts_0000.uni' )
        pp.save( out + 'parts_%04d.uni' % it )

    #objects = [ flags, phi, pp ]
    objects = [ pp ]
    fname = out + 'fluid_data_%04d.vdb' % it
    print( fname )
    save( name=fname, objects=objects )

# loop
while it < it_max:
    emphasize( '\n-----------------\n- time: %g(/%d)' % ( it, it_max ) )
    print( 'n=%d' % pp.pySize() )

    s.timestep = ( 1 - it % 1 ) * dt

    print( '- mapPartsToMAC' )
    mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
    
    #vel.printGrid()
    
    markFluidCells( parts=pp, flags=flags )

    # forces
    if 1:
        print( '- forces' )
        addGravityNoScale( flags=flags, vel=vel, gravity=(0,gravity,0) )
        #addGravity( flags=flags, vel=vel, gravity=(0,gravity,0) ) # adaptive to grid size

    #vel.printGrid()

    # set solid
    setWallBcs(flags=flags, vel=vel) # clears velocity from solid
    print( '- setWallBcs' )
    
    # pressure solve
    if 1:
        print( '- pressure' )
        solvePressure( flags=flags, vel=vel, pressure=pressure, phi=phi )

    #vel.printGrid()

    # we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
    extrapolateMACSimple( flags=flags, vel=vel, distance=res ) # 4
    
    # save position
    pp.getPosPdata( target=pos1 )
    
    # FLIP velocity update
    print( '- FLIP velocity update' )
    flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=1- 0 )
    
    #vel.printGrid()
    
    # advect particles 
    print( '- advectInGrid' )
    pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntEuler, deleteInObstacle=False ) # IntEuler, IntRK2, IntRK4

    # fixed vol
    include_walls = false
    if 1:
        flags.mark_interface()
        tic()
        s.timestep = fixed_volume_advection( pp=pp, x0=pos1, flags=flags, dt=s.timestep, dim=dim, part_per_cell_1d=part_per_cell_1d, state=0, phi=phi, it=it )
        print( '      ', end='' )
        toc()

        # if using band
        if 0:
            include_walls = true
            bSaveParts = 0

    # position solver, Thuerey21
    if 0:
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
    gridParticleIndex( parts=pp, flags=flags, indexSys=pindex, index=gpi )
    unionParticleLevelset( pp, pindex, flags, gpi, phi, radiusFactor ) 
    extrapolateLsSimple( phi=phi, distance=4, inside=True, include_walls=include_walls ) # 4

    # level set and mesh
    if bSaveParts:
        improvedParticleLevelset( pp, pindex, flags, gpi, phi, radiusFactor, 1, 1 , 0.4, 3.5 )

        # mesh
        phi.setBound( value=0., boundaryWidth=1 )
        phi.createMesh( mesh )

    # print/write
    if 0:
        flags.printGrid()
        vel.printGrid()
    
        #pp.printParts()
        pp.writeParticlesText( out + 'flipt_%04d.txt' % t )
    
    # step
    print( '- step' )
    s.step()

    it = it + s.timestep / dt
    if abs( it - round(it) ) < 1e-7:
        it = round( it )

        if bScreenShot:
            gui.screenshot( out + 'frame_%04d.png' % it ); # slow

        # save particle data
        if bSaveParts:
            if bSaveUni:
                # save particle data for flip03_gen.py surface generation scene
                pp.save( out + 'parts_%04d.uni' % it )

            #objects = [ flags, phi, pp ]
            objects = [ pp ]
            save( name=out + 'fluid_data_%04d.vdb' % it, objects=objects )
        
# pause
if 1:
    print( 'press enter...' )
    #keyboard.read_key()
    input()

