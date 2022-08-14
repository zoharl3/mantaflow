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

bSaveParts = 0
bScreenShot = 1

# solver params
dim = 3 # 2, 3
it_max = 1500 # 1500
part_per_cell_1d = 1 # 3, 2
res = 48 # 32, 48, 64, 128

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

# my dam
fluidbox = Box( parent=s, p0=gs*( vec3(0,0,0) ), p1=gs*( vec3(0.4,0.8,1) ) ) 

# square
#t = vec3(0.15, 0.15,0)
#t = vec3(0.3, 0.3, 0)
#fluidbox = Box( parent=s, p0=gs*( t + vec3(0,0,0) ), p1=gs*( t + vec3(0.4,0.4,1) ) )

# square
#t = vec3(0.4, 0.4,0)
#fluidbox = Box( parent=s, p0=gs*( t+vec3(0,0,0) ), p1=gs*( t+vec3(0.1,0.2,1) ) ) 

# manta dam
#fluidbox = Box( parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(0.4,0.6,1)) 

# phi
phi = fluidbox.computeLevelset()
flags.updateFromLevelset( phi )

sampleFlagsWithParticles( flags=flags, parts=pp, discretization=part_per_cell_1d, randomness=0.2 ) # 0.2
    
copyFlagsToFlags(flags, flagsPos)
flags.initDomain(boundaryWidth=0, phiWalls=phiObs)

np = pp.pySize()
print( '# particles:', np )
pos1 = s.create(PdataVec3)
pos1.pyResize( np )

if GUI:
    gui = Gui()
    gui.setRealGridDisplay( 0 )
    gui.setVec3GridDisplay( 0 )
    gui.show()
    gui.pause()
    
if bScreenShot:
    gui.screenshot( out + 'frame_%04d.png' % 0 ); # slow

if bSaveParts:
    pressure.save( out + 'ref_parts_0000.uni' );

# loop
it = 0
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
    if 1:
        flags.mark_interface()
        tic()
        s.timestep = fixed_volume_advection( pp=pp, x0=pos1, flags=flags, dt=s.timestep, dim=dim, part_per_cell_1d=part_per_cell_1d, state=0 )
        print( '      ', end='' )
        toc()

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
    
    # print/write
    if 0:
        flags.printGrid()
        vel.printGrid()
    
        #pp.printParts()
        pp.writeParticlesText( out + 'flipt_%04d.txt' % t )
    
    # mesh
    if 0 and dim == 3:
        phi.createMesh( mesh )

    # step
    print( '- step' )
    s.step()

    it = it + s.timestep / dt
    if abs( it - round(it) ) < 1e-7:
        it = round( it )

        if bScreenShot:
            gui.screenshot( out + 'frame_%04d.png' % it ); # slow

        # save particle data for flip03_gen.py surface generation scene
        if bSaveParts:
            pp.save( out + 'parts_%04d.uni' % it )

            # note: when saving pdata fields, they must be accompanied by and listed before their parent pp
            objects = [flags, phi, pressure, vel, pVel, pp]
            save( name=out + 'fluid_data_%04d.vdb' % it, objects=objects )
        
if 0:
    print( 'press enter...' )
    #keyboard.read_key()
    input()

