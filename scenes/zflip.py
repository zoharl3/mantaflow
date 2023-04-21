
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

# scene setup
def scene_setup():
    #flags.initDomain( boundaryWidth=boundary_width ) 
    flags.initDomain( boundaryWidth=boundary_width, phiWalls=phiObs ) 

    if 0: # dam
        # my dam
        #fluidbox = Box( parent=s, p0=gs*( Vec3(0, 0, 0.3) ), p1=gs*( Vec3(0.4, 0.8, .7) ) )
        fluidbox = Box( parent=s, p0=gs*( Vec3(0, 0, 0.35) ), p1=gs*( Vec3(0.3, 0.6, .65) ) ) # new dam (smaller, less crazy)

        # flip05_nbflip.py
        #fluidbox = Box( parent=s, p0=gs*Vec3(0, 0.15, 0), p1=gs*Vec3(0.4, 0.5, 0.8) )

        # square
        if 0:
            t1 = 0.4 # 0.15, 0.3, .4
            sz1 = .1 # .2, .4
            t = Vec3(t1, t1, 0)
            sz = Vec3(sz1, sz1, 1)
            fluidbox = Box( parent=s, p0=gs*( t + Vec3(0,0,0) ), p1=gs*( t + sz ) )

        # manta dam
        #fluidbox = Box( parent=s, p0=gs*Vec3(0,0,0), p1=gs*Vec3(0.4,0.6,1)) 

        # phi
        phi = fluidbox.computeLevelset()
        flags.updateFromLevelset( phi )

    elif 0: # falling drop
        fluidBasin = Box( parent=s, p0=gs*Vec3(0,0,0), p1=gs*Vec3(1.0,0.1,1.0)) # basin
        dropCenter = Vec3(0.5,0.3,0.5)
        dropRadius = 0.1
        fluidDrop  = Sphere( parent=s , center=gs*dropCenter, radius=res*dropRadius)
        fluidVel   = Sphere( parent=s , center=gs*dropCenter, radius=res*(dropRadius+0.05) )
        fluidSetVel= Vec3(0,-1,0)
        phi = fluidBasin.computeLevelset()
        phi.join( fluidDrop.computeLevelset() ) # add drop
        flags.updateFromLevelset( phi )

    elif 0: # basin
        # water
        fluidbox = Box( parent=s, p0=gs*( Vec3(0, 0.5, 0) ), p1=gs*( Vec3(1, 0.9, 1) ) )
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

    else: # a low, full box with an obstacle
        # water
        h = 0.25 # 0.3, 0.9
        fluidbox = Box( parent=s, p0=gs*( Vec3(0, 0., 0) ), p1=gs*( Vec3(1, h, 1) ) )
        phi = fluidbox.computeLevelset()
        flags.updateFromLevelset( phi )

        # moving obstacle
        bObs = 1
        if bObs:
            obs_rad = .05*res # .05, .1, .3
            obs_center = gs*Vec3( 0.5, 0.95 - obs_rad/res, 0.5 ) # y:0.5, 0.9
            shape = Box( parent=s, p0=obs_center - Vec3(obs_rad), p1=obs_center + Vec3(obs_rad) )
            #shape = Sphere( parent=s, center=obs_center, radius=obs_rad )
            phiObsInit.copyFrom( phiObs )
            phiObs.join( shape.computeLevelset() )

            obs_vel_vec = Vec3( 0, -0, 0 )
            obsVel.setConst( obs_vel_vec )
            obsVel.setBound( value=Vec3(0.), boundaryWidth=boundary_width+1 )

    #phi.printGrid()
    #phiObs.printGrid()

    sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=part_per_cell_1d, randomness=0.1 ) # 0.05, 0.1, 0.2

    # also sets boundary flags for phiObs (and shrinks it)
    if bObs:
        if 1:
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

# correct21 (position solver, Thuerey21)
# The band requires fixing, probably identifying non-band fluid cells as full. In the paper, it's listed as future work.
class Correct21:
    def __init__( self, dim, s, part_per_cell_1d, pp ):
        self.density = s.create(RealGrid)
        self.Lambda = s.create(RealGrid)
        self.deltaX = s.create(MACGrid)
        self.flagsPos = s.create(FlagGrid)
        self.pMass = pp.create(PdataReal)
        self.mass = 1.0 / part_per_cell_1d**dim
        self.resampleParticles = False # must be a boolean type since passing to cpp later
        if self.resampleParticles:
            gCnt = s.create(IntGrid)

    def main( self ):
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

class Simulation:
    # flags
    bMesh       = 1
    bSaveParts  = 0 # .vdb
    bSaveUni    = 0 # .uni
    if bSaveParts or bSaveUni:
        bMesh = 1

    bScreenShot = 1

    # solver params
    dim = 2 # 2, 3
    part_per_cell_1d = 2 # 3, 2(default), 1
    it_max = 1400 # 300, 500, 1200, 1400
    res = 32 # 32, 48, 64(default), 96, 128(large), 256(, 512 is too large)

    b_fixed_vol = 1
    b_correct21 = 0

    narrowBand = bool( 0 )
    narrowBandWidth = 5 # 32:5, 64:6, 96:6, 128:8

    ###

    #gs = Vec3( res, res, 5 ) # debug thin 3D; at least z=5 if with obstacle (otherwise, it has 0 velocity?)
    gs = Vec3( res, res, res )

    if b_correct21:
        b_fixed_vol = 1
        narrowBand = bool( 0 )

    combineBandWidth = narrowBandWidth - 1
    dt = .2 # .2(default), .5, 1(flip5, easier to debug)
    if dim == 2:
        gs.z = 1
        bMesh = 0
        bSaveParts = 0
        set_print_2D( True )
    if 0: # sample 1D; requires changing sampleLevelsetWithParticles()
        ppc = part_per_cell_1d
    else:
        ppc = part_per_cell_1d**dim

    boundary_width = 0

    s = Solver( name='main', gridSize=gs, dim=dim )
    gravity = -0.1 # -0.1
    gravity *= math.sqrt( res )
    #gravity = -0.003 # flip5

    # grids and particles
    flags    = s.create(FlagGrid)
    flags2   = s.create(FlagGrid)
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
    bObs = 0
    obs_center = Vec3( 0 )
    obs_rad = 0
    obs_vel_vec = Vec3( 0 )
    fractions = s.create(MACGrid) # Sticks to walls and becomes slow without this? Not anymore. Shrinks an obstacle box and draws it inaccurately.
    mesh     = s.create(Mesh)
    bfs      = s.create(IntGrid) # discrete distance from surface

    # acceleration data for particle
    pindex = s.create(ParticleIndexSystem)
    gpi    = s.create(IntGrid)

    correct21 = Correct21( dim, s, part_per_cell_1d, pp )

    def main( self ):
        print()
        print( 'dim:', self.dim, ', res:', self.res, ', ppc:', self.ppc )
        print( 'narrowBand:', self.narrowBand, ', narrowBandWidth:', self.narrowBandWidth )
        print( 'b_fixed_vol:', self.b_fixed_vol )
        print( 'gravity: %0.02f' % self.gravity )
        print( 'timestep:', self.dt )

        # adaptive time stepping; from flip5
        b_adaptive_time_step = 0
        if b_adaptive_time_step:
            self.s.frameLength = 1.0   # length of one frame (in "world time")
            self.s.timestep    = 1.0
            self.s.timestepMin = 0.5   # time step range
            self.s.timestepMax = 1.0
            self.s.cfl         = 5.0   # maximal velocity per cell, 0 to use fixed timesteps

        # size of particles 
        radiusFactor = 1.0

        print( '# particles:', self.pp.pySize() )
        V0 = float( self.pp.pySize() ) / self.ppc

        if 1 and GUI:
            gui = Gui()
            for i in range( 2 ):
                gui.nextMeshDisplay() # 0:full, 1:hide, 2:x-ray
            gui.setRealGridDisplay( 0 )
            gui.setVec3GridDisplay( 1 )
            if 0 and dim == 3: # camera
                gui.setCamPos( 0, 0, -2.2 ) # drop
                gui.setCamRot( 35, -30, 0 )
            if 0 and bMesh:
                gui.toggleHideGrids()
            gui.show()
            #gui.pause()
        else:
            bScreenShot = 0

        it = 0
        it2 = 0

        if self.bScreenShot:
            gui.screenshot( out + 'frame_%04d.png' % it ); # slow

        if self.bSaveParts:
            if self.bSaveUni:
                self.pressure.save( out + 'ref_parts_0000.uni' )
                self.pp.save( out + 'parts_%04d.uni' % it )

            objects = [ flags, phi, pp ] # need the 3 of them for volumetric .vdb
            #objects = [ pp ]
            fname = out + 'fluid_data_%04d.vdb' % it
            print( fname )
            save( name=fname, objects=objects ) # error in debug mode "string too long?"

        # measure
        f_measure = open( out + '_measure.txt', 'w' )

        # loop
        ret = 0
        n_obs_skip = 0
        while 1:
            emphasize( '\n-----------------\n- time: %g(/%d; it2=%d)' % ( it, self.it_max, it2 ) )
            print( '- n=%d' % self.pp.pySize() )

            if 1 and ret != 0:
                error( f'Error: ret={ret}' )
                ret = 0
                break

            if not it < self.it_max:
                break

            tic()

            maxVel = self.vel.getMaxAbs()
            if not b_adaptive_time_step:
                self.s.frameLength = self.dt
                self.s.timestep = ( 1 - it % 1 ) * self.dt
            else: # adaptive for flip5
                self.s.adaptTimestep( maxVel ) # enable init above

            # map particle velocities to grid
            print( '- mapPartsToMAC' )
            if 1 and self.narrowBand:
                # combine particles velocities with advected grid velocities; stores mapWeights; saves velOld
                mapPartsToMAC( vel=self.velParts, flags=self.flags, velOld=self.velOld, parts=self.pp, partVel=self.pVel, weight=self.mapWeights )
                # extrapolate velocities throughout the liquid region
                extrapolateMACFromWeight( vel=self.velParts , distance=2, weight=self.mapWeights )
                #vel.printGrid()
                #velParts.printGrid()
                combineGridVel( vel=self.velParts, weight=self.mapWeights, combineVel=self.vel, phi=self.phi, narrowBand=self.combineBandWidth, thresh=0 )
                #limit_velocity( vel, pVel, 256 ) # 64:15, 128:?
                self.velOld.copyFrom( self.vel ) # save before forces and pressure; like the end of mapPartsToMAC()
                #print( '>> combine' )
                #vel.printGrid()
            elif 1:
                # map particle velocities to grid
                mapPartsToMAC( vel=self.vel, flags=self.flags, velOld=self.velOld, parts=self.pp, partVel=self.pVel, weight=self.mapWeights )
                extrapolateMACFromWeight( vel=self.vel , distance=2, weight=self.mapWeights )

            # moving obstacle
            if bObs:
                #flags.printGrid()
                dv = self.s.timestep * Vec3( 0, 1*self.gravity, 0 )
                if obs_center.y - obs_rad > 2: # move
                    print( '- move obstacle' )
                    obs_vel_vec += dv
                else: # stay
                    print( '- obstacle stopped' )
                    obs_vel_vec = Vec3( 0. )

                # obsVel
                if 1:
                    obs_vel_vec2 = obs_vel_vec + dv # add some velocity in case it stopped--to remove remaining particles from the bottom
                    obsVel.setConst( obs_vel_vec2 )
                    obsVel.setBound( value=Vec3( 0. ), boundaryWidth=boundary_width + 1 )
                    #obsVel.printGrid()

                # phiObs
                print( f'  - obs_center={obs_center}, obs_rad={obs_rad}' )
                p0 = obs_center - Vec3( obs_rad )
                p1 = obs_center + Vec3( obs_rad )
                if dim == 2:
                    p0.z = p1.z = 0.5
                shape = Box( parent=s, p0=p0, p1=p1 )
                #shape = Sphere( parent=s, center=obs_center, radius=obs_rad )
                phiObs.copyFrom( phiObsInit )
                phiObs.join( shape.computeLevelset() )
                #phiObs.printGrid()

                if 1:
                    mark_obstacle_box( flags=flags, p0=p0, p1=p1 )
                elif 1:
                    setObstacleFlags( flags=flags, phiObs=phiObs )
                else: # more precise
                    updateFractions( flags=flags, phiObs=phiObs, fractions=fractions, boundaryWidth=boundary_width )
                    setObstacleFlags( flags=flags, phiObs=phiObs, fractions=fractions )
                #flags.printGrid()

            # update flags; there's also flags.updateFromLevelset()
            if 1:
                print( '- markFluidCells (update flags)' )
                markFluidCells( parts=pp, flags=flags ) # better for a moving obstacle?
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
                    g = gravity*s.timestep/1 # see addGravity(); assuming s.mDt=1
                    #g = gravity
                    addGravity( flags=flags, vel=vel, gravity=(0, g, 0), scale=bool(bscale) )

                # vortex
                if 0:
                    #c = gs/2
                    c = Vec3( res/2, 0.3*res, res/2 )
                    vortex( pp=pp, dt=s.timestep, c=c, rad=0.1*res, h=0.9*res, pVel=pVel2 )
                    mapPartsToMAC( vel=vel2, flags=flags, velOld=vel2, parts=pp, partVel=pVel2 )
                    vel.add( vel2 )

            # set velocity for obstacles
            print( '- setWallBcs' )
            #setWallBcs( flags=flags, vel=vel )
            #setWallBcs( flags=flags, vel=vel, fractions=fractions )
            #setWallBcs( flags=flags, vel=vel, fractions=fractions, phiObs=phiObs, obvel=obsVel ) # calls KnSetWallBcsFrac, which doesn't work?
            #obsVel.printGrid()
            #vel.printGrid()
            setWallBcs( flags=flags, vel=vel, obvel=obsVel ) # calls KnSetWallBcs
            #vel.printGrid()
            #flags.printGrid()

            # pressure solve
            if 1:
                print( '- pressure' )
                tic()
                solvePressure( flags=flags, vel=vel, pressure=pressure, phi=phi )
                print( '  (pressure) ', end='' )
                toc()

                #setWallBcs( flags=flags, vel=vel )
                #setWallBcs( flags=flags, vel=vel, fractions=fractions )
                #setWallBcs( flags=flags, vel=vel, fractions=fractions, phiObs=phiObs, obvel=obsVel )
                setWallBcs( flags=flags, vel=vel, obvel=obsVel )
                #vel.printGrid()

            dist = min( int( maxVel*1.25 + 2 ), 8 ) # res
            print( '- extrapolate MAC Simple (dist=%0.1f)' % dist )
            extrapolateMACSimple( flags=flags, vel=vel, distance=dist, intoObs=False )
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
            if not b_fixed_vol and not b_correct21:
                pushOutofObs( parts=pp, flags=flags, phiObs=phiObs ) # creates issues for correct21 and fixedVol
            # advect phi; why? the particles should determine phi, which should flow on its own; disabling this creates artifacts in flip5; it makes it worse for fixed_vol
            if 1 and not b_fixed_vol:
                advectSemiLagrange( flags=flags, vel=vel, grid=phi, order=1 )
                if 0:
                    flags.updateFromLevelset( phi ) # creates in 3D an extra layer of fluid without particles
            # advect grid velocity
            if narrowBand:
                advectSemiLagrange( flags=flags, vel=vel, grid=vel, order=2 )

            # fixed volume (my scheme)
            #flags.printGrid()
            include_walls = false
            obs_naive = 0
            obs_stop = 0
            if b_fixed_vol:
                phi.setBoundNeumann( 0 ) # make sure no new particles are placed at outer boundary
                #phi.printGrid()

                pVel.setSource( vel, isMAC=True ) # set source grid for resampling, used in insertBufferedParticles()

                dt_bound = 0
                #dt_bound = dt/4
                #dt_bound = s.timestep/4
                #dt_bound = max( dt_bound, dt/4 )

                # obs_vel: it modifies it to either one cell distance or zero, staying in place and losing velocity (unlike particles)
                obs_vel_vec3 = obs_vel_vec
                if obs_naive:
                    obs_vel_vec3 = Vec3(0)
                ret2 = fixed_volume_advection( pp=pp, pVel=pVel, flags=flags, dt=s.timestep, dt_bound=dt_bound, dim=dim, ppc=ppc, phi=phi, it=it2, use_band=narrowBand, band_width=narrowBandWidth, bfs=bfs, obs_center=obs_center, obs_rad=obs_rad, obs_vel=obs_vel_vec3 )
                if not ret2:
                    ret = -1
                    s.timestep *= -1
                else:
                    s.timestep = ret2[0]
                    if not obs_naive:
                        #obs_vel_vec = Vec3( ret2[1], ret2[2], ret2[3] )
                        obs_stop = ret2[1]

                # if using band
                if 0 and narrowBand:
                    include_walls = true
            #flags.printGrid()

            # update obstacle
            if bObs:
                # test obstacle position
                print( '  - obs_stop=%d' % obs_stop )
                if 1 and not obs_stop:
                    obs_center2 = obs_center + s.timestep * obs_vel_vec
                    p0 = obs_center2 - Vec3( obs_rad )
                    p1 = obs_center2 + Vec3( obs_rad )
                    if dim == 2:
                        p0.z = p1.z = 0.5
                    flags2.copyFrom( flags )
                    if not mark_obstacle_box( flags=flags2, p0=p0, p1=p1 ):
                        emphasize( '  - obstacle position is invalid; stopping the obstacle' )
                        assert( not b_fixed_vol or obs_naive )
                        obs_stop = 1

                print( f'  - obs_vel_vec={obs_vel_vec}, dt={dt}, obs_center={obs_center}' )
                obs_center2 = obs_center + s.timestep * obs_vel_vec
                # slow down by skipping grid progress 
                if int( obs_center2.y - obs_rad ) == int( obs_center.y - obs_rad ):
                    obs_center = obs_center2
                else:
                    if not obs_stop:
                        n_obs_skip += 1
                        if n_obs_skip > 0: # 0, 2; how many steps to skip
                            n_obs_skip = 0
                            obs_center = obs_center2
                    else:
                        obs_vel_vec = Vec3(0)

            # correct21
            if b_correct21:
                correct21.main( flags, phiObs )
            
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

# main
if __name__ == '__main__':
    # auto-flush
    sys.stdout.reconfigure( line_buffering=True )

    out = r'c:/prj-external-libs/mantaflow/out/'

    os.system( 'rm %s*.*' % out )
    os.system( 'cp %s../video.bat %s' % (out, out) )

    # (debug) for consistent result; for large res, the step() hangs?
    if 0:
        limit_to_one_core()

    sim = Simulation()
    sim.main()
