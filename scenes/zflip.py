
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
            self.gCnt = s.create(IntGrid)

    def main( self, sol, flags, pp, vel, pindex, gpi, phiObs ):
        print( '- position solver' )
        copyFlagsToFlags(flags, self.flagsPos)
        mapMassToGrid(flags=self.flagsPos, density=self.density, parts=pp, source=self.pMass, deltaX=self.deltaX, phiObs=phiObs, dt=sol.timestep, particleMass=self.mass, noDensityClamping=self.resampleParticles)
        #gui.pause()
        
        # resample particles
        if self.resampleParticles:
            print( '    - resample particles' )
            gridParticleIndex(parts=pp, indexSys=pindex, flags=flags, index=gpi, counter=self.gCnt)
            #apicMapPartsToMAC(flags=flags, vel=vel, parts=pp, partVel=pVel, cpx=apic_pCx, cpy=apic_pCy, cpz=apic_pCz, mass=apic_mass)
            resampeOverfullCells(vel=vel, density=self.density, index=gpi, indexSys=pindex, part=pp, pVel=pVel, dt=s.timestep)
    
        # position solver
        print( '    - solve pressure due to density' )
        solvePressureSystem(rhs=self.density, vel=vel, pressure=self.Lambda, flags=self.flagsPos, cgAccuracy=1e-3)
        computeDeltaX(deltaX=self.deltaX, Lambda=self.Lambda, flags=self.flagsPos)
        mapMACToPartPositions(flags=self.flagsPos, deltaX=self.deltaX, parts=pp, dt=sol.timestep)
        
        # print
        if 0:
            #flags.printGrid()
            #self.flagsPos.printGrid()
            self.density.printGrid()
            self.Lambda.printGrid()
            self.deltaX.printGrid()

class moving_obstacle:
    def __init__( self, sol ):
        self.exists = 0
        self.vel = sol.create(MACGrid)
        self.center = Vec3( 0 )
        self.rad = 0
        self.vel_vec = Vec3( 0 )
        self.phi_init = sol.create(LevelsetGrid)

class Simulation:
    # flags
    bMesh       = 1
    bSaveParts  = 0 # .vdb
    bSaveUni    = 0 # .uni
    if bSaveParts or bSaveUni:
        bMesh = 1

    bScreenShot = 1

    # params
    dim = 2 # 2, 3
    part_per_cell_1d = 2 # 3, 2(default), 1
    it_max = 1400 # 300, 500, 1200, 1400
    res = 32 # 32, 48, 64(default), 96, 128(large), 256(, 512 is too large)

    b_fixed_vol = 1
    b_correct21 = 0

    narrowBand = bool( 1 )
    narrowBandWidth = 5 # 32:5, 64:6, 96:6, 128:8

    ###

    #gs = Vec3( res, res, 5 ) # debug thin 3D; at least z=5 if with obstacle (otherwise, it has 0 velocity?)
    gs = Vec3( res, res, res )

    if dim == 2:
        gs.z = 1
        bMesh = 0
        bSaveParts = 0

    dt = .2 # .2(default), .5, 1(flip5, easier to debug)

    gravity = -0.1 # -0.1
    gravity *= math.sqrt( res )
    #gravity = -0.003 # flip5

    sol = Solver( name='main', gridSize=gs, dim=dim )

    flags = sol.create(FlagGrid)
    vel = sol.create(MACGrid)
    pp = sol.create(BasicParticleSystem)
    phiObs = sol.create(LevelsetGrid)
    phi = []

    obs = moving_obstacle( sol )

    boundary_width = 0

    def setup_scene( self ):
        #self.flags.initDomain( boundaryWidth=self.boundary_width ) 
        self.flags.initDomain( boundaryWidth=self.boundary_width, phiWalls=self.phiObs ) 

        if 1: # dam
            # my dam
            #fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0, 0, 0.3) ), p1=self.gs*( Vec3(0.4, 0.8, .7) ) )
            fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0, 0, 0.35) ), p1=self.gs*( Vec3(0.3, 0.6, .65) ) ) # new dam (smaller, less crazy)

            # flip05_nbflip.py
            #fluidbox = Box( parent=self.sol, p0=self.gs*Vec3(0, 0.15, 0), p1=self.gs*Vec3(0.4, 0.5, 0.8) )

            # square
            if 0:
                t1 = 0.4 # 0.15, 0.3, .4
                sz1 = .1 # .2, .4
                t = Vec3(t1, t1, 0)
                sz = Vec3(sz1, sz1, 1)
                fluidbox = Box( parent=self.sol, p0=self.gs*( t + Vec3(0,0,0) ), p1=self.gs*( t + sz ) )

            # manta dam
            #fluidbox = Box( parent=self.sol, p0=self.gs*Vec3(0,0,0), p1=self.gs*Vec3(0.4,0.6,1)) 

            # phi
            self.phi = fluidbox.computeLevelset()
            self.flags.updateFromLevelset( self.phi )

        elif 0: # falling drop
            fluidBasin = Box( parent=self.sol, p0=self.gs*Vec3(0,0,0), p1=self.gs*Vec3(1.0,0.1,1.0)) # basin
            dropCenter = Vec3(0.5,0.3,0.5)
            dropRadius = 0.1
            fluidDrop  = Sphere( parent=self.sol , center=self.gs*dropCenter, radius=res*dropRadius)
            fluidVel   = Sphere( parent=self.sol , center=self.gs*dropCenter, radius=res*(dropRadius+0.05) )
            fluidSetVel= Vec3(0,-1,0)
            self.phi = fluidBasin.computeLevelset()
            self.phi.join( fluidDrop.computeLevelset() ) # add drop
            self.flags.updateFromLevelset( self.phi )

        elif 0: # basin
            # water
            fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0, 0.5, 0) ), p1=self.gs*( Vec3(1, 0.9, 1) ) )
            self.phi = fluidbox.computeLevelset()
            self.flags.updateFromLevelset( self.phi )

            # obstacle
            if 1:
                #mesh2 = self.sol.create(Mesh) # it renders only one mesh (mLocalMesh)?

                #mesh.load( r'c:\prj\mantaflow_mod\resources\cube1.obj' )
                #mesh.scale( Vec3(1) )

                mesh.load( r'c:\prj\mantaflow_mod\resources\funnel.obj' )
                mesh.scale( Vec3(res) ) # the scale needs to be in all axes

                mesh.offset( self.gs * Vec3(0.5, 0.3, 0.5) )
                meshObs = self.sol.create( LevelsetGrid )
                mesh.computeLevelset( meshObs, 2 ) # uses normals, thus a smooth mesh is better
                #meshObs.printGrid()
                #self.phiObs.printGrid()
                self.phiObs.join( meshObs )
                self.phi.subtract( self.phiObs )

        else: # a low, full box with an obstacle
            # water
            h = 0.25 # 0.3, 0.9
            fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0, 0., 0) ), p1=self.gs*( Vec3(1, h, 1) ) )
            self.phi = fluidbox.computeLevelset()
            self.flags.updateFromLevelset( self.phi )

            # moving obstacle
            self.obs.exists = 1
            if self.obs.exists:
                self.obs.rad = .05*self.res # .05, .1, .3
                self.obs.center = self.gs*Vec3( 0.5, 0.95 - self.obs.rad/self.res, 0.5 ) # y:0.5, 0.9
                shape = Box( parent=self.sol, p0=self.obs.center - Vec3(self.obs.rad), p1=self.obs.center + Vec3(self.obs.rad) )
                #shape = Sphere( parent=self.sol, center=self.obs.center, radius=self.obs.rad )
                self.obs.phi_init.copyFrom( self.phiObs )
                self.phiObs.join( shape.computeLevelset() )

                self.obs.vel_vec = Vec3( 0, -0, 0 )
                self.obs.vel.setConst( self.obs.vel_vec )
                self.obs.vel.setBound( value=Vec3(0.), boundaryWidth=self.boundary_width+1 )

        #self.phi.printGrid()
        #self.phiObs.printGrid()

        sampleLevelsetWithParticles( phi=self.phi, flags=self.flags, parts=self.pp, discretization=self.part_per_cell_1d, randomness=0.1 ) # 0.05, 0.1, 0.2

        # also sets boundary flags for phiObs (and shrinks it)
        if self.obs.exists:
            if 1:
                setObstacleFlags( flags=self.flags, phiObs=self.phiObs )
            else:
                updateFractions( flags=self.flags, phiObs=self.phiObs, fractions=fractions, boundaryWidth=self.boundary_width )
                setObstacleFlags( flags=self.flags, phiObs=self.phiObs, fractions=fractions )

    def main( self ):
        if self.b_correct21:
            self.b_fixed_vol = 0
            self.narrowBand = bool( 0 )

        combineBandWidth = self.narrowBandWidth - 1

        if self.dim == 2:
            set_print_2D( True )

        if 0: # sample 1D; requires changing sampleLevelsetWithParticles()
            ppc = self.part_per_cell_1d
        else:
            ppc = self.part_per_cell_1d**self.dim

        flags2   = self.sol.create(FlagGrid)
        vel2     = self.sol.create(MACGrid)
        velOld   = self.sol.create(MACGrid)
        velParts = self.sol.create(MACGrid)
        pressure = self.sol.create(RealGrid)
        mapWeights = self.sol.create(MACGrid)

        # add velocity data to particles
        pVel     = self.pp.create(PdataVec3)
        pVel2     = self.pp.create(PdataVec3)
        phiParts = self.sol.create(LevelsetGrid)
        phiMesh  = self.sol.create(LevelsetGrid)
        fractions = self.sol.create(MACGrid) # Sticks to walls and becomes slow without this? Not anymore. Shrinks an obstacle box and draws it inaccurately.
        mesh     = self.sol.create(Mesh)
        bfs      = self.sol.create(IntGrid) # discrete distance from surface

        # acceleration data for particle
        pindex = self.sol.create(ParticleIndexSystem)
        gpi    = self.sol.create(IntGrid)

        correct21 = Correct21( self.dim, self.sol, self.part_per_cell_1d, self.pp )

        print()
        print( 'dim:', self.dim, ', res:', self.res, ', ppc:', ppc )
        print( 'narrowBand:', self.narrowBand, ', narrowBandWidth:', self.narrowBandWidth )
        print( 'b_fixed_vol:', self.b_fixed_vol )
        print( 'gravity: %0.02f' % self.gravity )
        print( 'timestep:', self.dt )

        # adaptive time stepping; from flip5
        b_adaptive_time_step = 0
        if b_adaptive_time_step:
            self.sol.frameLength = 1.0   # length of one frame (in "world time")
            self.sol.timestep    = 1.0
            self.sol.timestepMin = 0.5   # time step range
            self.sol.timestepMax = 1.0
            self.sol.cfl         = 5.0   # maximal velocity per cell, 0 to use fixed timesteps

        # setup scene
        self.setup_scene()

        # size of particles 
        radiusFactor = 1.0

        print( '# particles:', self.pp.pySize() )
        V0 = float( self.pp.pySize() ) / ppc

        # phi is influenced by the walls for some reason
        # create a level set from particles
        gridParticleIndex( parts=self.pp, flags=self.flags, indexSys=pindex, index=gpi )
        unionParticleLevelset( self.pp, pindex, self.flags, gpi, phiParts, radiusFactor )
        self.phi.copyFrom( phiParts )
        if self.narrowBand:
            extrapolateLsSimple( phi=self.phi, distance=self.narrowBandWidth+2, inside=True )
        else:
            extrapolateLsSimple( phi=self.phi, distance=4, inside=True ) # 4
        #self.phi.printGrid()

        if 1 and GUI:
            gui = Gui()
            for i in range( 2 ):
                gui.nextMeshDisplay() # 0:full, 1:hide, 2:x-ray
            gui.setRealGridDisplay( 0 )
            gui.setVec3GridDisplay( 1 )
            if 0 and self.dim == 3: # camera
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

            objects = [ self.flags, self.phi, self.pp ] # need the 3 of them for volumetric .vdb
            #objects = [ self.pp ]
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
                self.sol.frameLength = self.dt
                self.sol.timestep = ( 1 - it % 1 ) * self.dt
            else: # adaptive for flip5
                self.sol.adaptTimestep( maxVel ) # enable init above

            # map particle velocities to grid
            print( '- mapPartsToMAC' )
            if 1 and self.narrowBand:
                # combine particles velocities with advected grid velocities; stores mapWeights; saves velOld
                mapPartsToMAC( vel=velParts, flags=self.flags, velOld=velOld, parts=self.pp, partVel=pVel, weight=mapWeights )
                # extrapolate velocities throughout the liquid region
                extrapolateMACFromWeight( vel=velParts , distance=2, weight=mapWeights )
                #vel.printGrid()
                #velParts.printGrid()
                combineGridVel( vel=velParts, weight=mapWeights, combineVel=self.vel, phi=self.phi, narrowBand=combineBandWidth, thresh=0 )
                #limit_velocity( vel, pVel, 256 ) # 64:15, 128:?
                velOld.copyFrom( self.vel ) # save before forces and pressure; like the end of mapPartsToMAC()
                #print( '>> combine' )
                #vel.printGrid()
            elif 1:
                # map particle velocities to grid
                mapPartsToMAC( vel=self.vel, flags=self.flags, velOld=velOld, parts=self.pp, partVel=pVel, weight=mapWeights )
                extrapolateMACFromWeight( vel=self.vel , distance=2, weight=mapWeights )

            # moving obstacle
            if self.obs.exists:
                #self.flags.printGrid()
                dv = self.sol.timestep * Vec3( 0, 1*self.gravity, 0 )
                if self.obs.center.y - self.obs.rad > 2: # move
                    print( '- move obstacle' )
                    self.obs.vel_vec += dv
                else: # stay
                    print( '- obstacle stopped' )
                    self.obs.vel_vec = Vec3( 0. )

                # obsVel
                if 1:
                    obs_vel_vec2 = self.obs.vel_vec + dv # add some velocity in case it stopped--to remove remaining particles from the bottom
                    self.obs.vel.setConst( obs_vel_vec2 )
                    self.obs.vel.setBound( value=Vec3( 0. ), boundaryWidth=self.boundary_width + 1 )
                    #self.obs.vel.printGrid()

                # phiObs
                print( f'  - obs_center={self.obs.center}, obs_rad={self.obs.rad}' )
                p0 = self.obs.center - Vec3( self.obs.rad )
                p1 = self.obs.center + Vec3( self.obs.rad )
                if self.dim == 2:
                    p0.z = p1.z = 0.5
                shape = Box( parent=self.sol, p0=p0, p1=p1 )
                #shape = Sphere( parent=s, center=self.obs.center, radius=self.obs.rad )
                self.phiObs.copyFrom( self.obs.phi_init )
                self.phiObs.join( shape.computeLevelset() )
                #self.phiObs.printGrid()

                if 1:
                    mark_obstacle_box( flags=self.flags, p0=p0, p1=p1 )
                elif 1:
                    setObstacleFlags( flags=self.flags, phiObs=self.phiObs )
                else: # more precise
                    updateFractions( flags=self.flags, phiObs=self.phiObs, fractions=fractions, boundaryWidth=self.boundary_width )
                    setObstacleFlags( flags=self.flags, phiObs=self.phiObs, fractions=fractions )
                #self.flags.printGrid()

            # update flags; there's also flags.updateFromLevelset()
            if 1:
                print( '- markFluidCells (update flags)' )
                markFluidCells( parts=self.pp, flags=self.flags ) # better for a moving obstacle?
                #markFluidCells( parts=self.pp, flags=self.flags, phiObs=self.phiObs )
                if self.narrowBand and ( not self.b_fixed_vol or it == 0 ):
                    update_fluid_from_phi( flags=self.flags, phi=self.phi, band_width=self.narrowBandWidth )
                #self.flags.printGrid()

            #self.vel.printGrid()
            #self.flags.printGrid()
            # forces
            if 1:
                print( '- forces' )
                
                # gravity
                if 1:
                    bscale = 0 # 1:adaptive to grid size; flip5
                    g =  self.gravity*self.sol.timestep/1 # see addGravity(); assuming sol.mDt=1
                    #g =  self.gravity
                    addGravity( flags=self.flags, vel=self.vel, gravity=(0, g, 0), scale=bool(bscale) )

                # vortex
                if 0:
                    #c = gs/2
                    c = Vec3( res/2, 0.3*res, res/2 )
                    vortex( pp=self.pp, dt=s.timestep, c=c, rad=0.1*res, h=0.9*res, pVel=pVel2 )
                    mapPartsToMAC( vel=vel2, flags=self.flags, velOld=vel2, parts=self.pp, partVel=pVel2 )
                    self.vel.add( vel2 )

            # set velocity for obstacles
            print( '- setWallBcs' )
            #setWallBcs( flags=self.flags, vel=self.vel )
            #setWallBcs( flags=self.flags, vel=self.vel, fractions=fractions )
            #setWallBcs( flags=self.flags, vel=self.vel, fractions=fractions, phiObs=self.phiObs, obvel=self.obs.vel ) # calls KnSetWallBcsFrac, which doesn't work?
            #self.obs.vel.printGrid()
            #self.vel.printGrid()
            setWallBcs( flags=self.flags, vel=self.vel, obvel=self.obs.vel ) # calls KnSetWallBcs
            #self.vel.printGrid()
            #self.flags.printGrid()

            # pressure solve
            if 1:
                print( '- pressure' )
                tic()
                solvePressure( flags=self.flags, vel=self.vel, pressure=pressure, phi=self.phi )
                print( '  (pressure) ', end='' )
                toc()

                #setWallBcs( flags=self.flags, vel=self.vel )
                #setWallBcs( flags=self.flags, vel=self.vel, fractions=fractions )
                #setWallBcs( flags=self.flags, vel=self.vel, fractions=fractions, phiObs=self.phiObs, obvel=self.obs.vel )
                setWallBcs( flags=self.flags, vel=self.vel, obvel=self.obs.vel )
                #self.vel.printGrid()

            dist = min( int( maxVel*1.25 + 2 ), 8 ) # res
            print( '- extrapolate MAC Simple (dist=%0.1f)' % dist )
            extrapolateMACSimple( flags=self.flags, vel=self.vel, distance=dist, intoObs=False )
            #self.flags.printGrid()
            #self.vel.printGrid()

            print( '- set particles\' pos0' )
            set_particles_pos0( pp=self.pp )

            # FLIP velocity update
            print( '- FLIP velocity update' )
            alpha = .1 # 0
            flipVelocityUpdate( vel=self.vel, velOld=velOld, flags=self.flags, parts=self.pp, partVel=pVel, flipRatio=1 - alpha )
            #self.vel.printGrid()
            
            # advect
            print( '- advect' )
            # advect particles
            self.pp.advectInGrid( flags=self.flags, vel=self.vel, integrationMode=IntEuler, deleteInObstacle=False, stopInObstacle=False ) # IntEuler, IntRK2, IntRK4
            if not self.b_fixed_vol and not self.b_correct21:
                pushOutofObs( parts=self.pp, flags=self.flags, phiObs=self.phiObs ) # creates issues for correct21 and fixedVol
            # advect phi; why? the particles should determine phi, which should flow on its own; disabling this creates artifacts in flip5; it makes it worse for fixed_vol
            if 1 and not self.b_fixed_vol:
                advectSemiLagrange( flags=self.flags, vel=self.vel, grid=self.phi, order=1 )
                if 0:
                    self.flags.updateFromLevelset( self.phi ) # creates in 3D an extra layer of fluid without particles
            # advect grid velocity
            if self.narrowBand:
                advectSemiLagrange( flags=self.flags, vel=self.vel, grid=self.vel, order=2 )

            # fixed volume (my scheme)
            #self.flags.printGrid()
            include_walls = false
            obs_naive = 0
            obs_stop = 0
            if self.b_fixed_vol:
                self.phi.setBoundNeumann( 0 ) # make sure no new particles are placed at outer boundary
                #self.phi.printGrid()

                pVel.setSource( self.vel, isMAC=True ) # set source grid for resampling, used in insertBufferedParticles()

                dt_bound = 0
                #dt_bound = dt/4
                #dt_bound = self.sol.timestep/4
                #dt_bound = max( dt_bound, dt/4 )

                # obs_vel: it modifies it to either one cell distance or zero, staying in place and losing velocity (unlike particles)
                obs_vel_vec3 = self.obs.vel_vec
                if obs_naive:
                    obs_vel_vec3 = Vec3(0)
                ret2 = fixed_volume_advection( pp=self.pp, pVel=pVel, flags=self.flags, dt=self.sol.timestep, dt_bound=dt_bound, dim=self.dim, ppc=ppc, phi=self.phi, it=it2, use_band=self.narrowBand, band_width=self.narrowBandWidth, bfs=bfs, obs_center=self.obs.center, obs_rad=self.obs.rad, obs_vel=obs_vel_vec3 )
                if not ret2:
                    ret = -1
                    self.sol.timestep *= -1
                else:
                    self.sol.timestep = ret2[0]
                    if not obs_naive:
                        #self.obs.vel_vec = Vec3( ret2[1], ret2[2], ret2[3] )
                        obs_stop = ret2[1]

                # if using band
                if 0 and self.narrowBand:
                    include_walls = true
            #self.flags.printGrid()

            # update obstacle
            if self.obs.exists:
                # test obstacle position
                print( '  - obs_stop=%d' % obs_stop )
                if 1 and not obs_stop:
                    obs_center2 = self.obs.center + self.sol.timestep * self.obs.vel_vec
                    p0 = obs_center2 - Vec3( self.obs.rad )
                    p1 = obs_center2 + Vec3( self.obs.rad )
                    if self.dim == 2:
                        p0.z = p1.z = 0.5
                    flags2.copyFrom( self.flags )
                    if not mark_obstacle_box( flags=flags2, p0=p0, p1=p1 ):
                        emphasize( '  - obstacle position is invalid; stopping the obstacle' )
                        assert( not self.b_fixed_vol or obs_naive )
                        obs_stop = 1

                print( f'  - obs_vel_vec={self.obs.vel_vec}, dt={self.dt}, obs_center={self.obs.center}' )
                obs_center2 = self.obs.center + self.sol.timestep * self.obs.vel_vec
                # slow down by skipping grid progress 
                if int( obs_center2.y - self.obs.rad ) == int( self.obs.center.y - self.obs.rad ):
                    self.obs.center = obs_center2
                else:
                    if not obs_stop:
                        n_obs_skip += 1
                        if n_obs_skip > 0: # 0, 2; how many steps to skip
                            n_obs_skip = 0
                            self.obs.center = obs_center2
                    else:
                        self.obs.vel_vec = Vec3(0)

            # correct21
            if self.b_correct21:
                correct21.main( self.sol, self.flags, self.pp, self.vel, pindex, gpi, self.phiObs )
            
            # create level set from particles
            if 1:
                gridParticleIndex( parts=self.pp, flags=self.flags, indexSys=pindex, index=gpi )
                unionParticleLevelset( self.pp, pindex, self.flags, gpi, phiParts, radiusFactor ) 
                if self.narrowBand:
                    # combine level set of particles with grid level set
                    self.phi.addConst( 1. ); # shrink slightly
                    self.phi.join( phiParts )
                    extrapolateLsSimple( phi=self.phi, distance=self.narrowBandWidth+2, inside=True, include_walls=include_walls )
                else:
                    # overwrite grid level set with level set of particles
                    self.phi.copyFrom( phiParts )
                    extrapolateLsSimple( phi=self.phi, distance=4, inside=True, include_walls=include_walls ) # 4

            # resample particles
            if not self.b_fixed_vol:
                pVel.setSource( self.vel, isMAC=True ) # set source grids for resampling, used in adjustNumber
                minParticles = ppc
                maxParticles = 2*minParticles # 2, 1(exacerbates artifact in flip5 dam 128?)
                if self.narrowBand:
                    self.phi.setBoundNeumann( 0 ) # make sure no particles are placed at outer boundary
                    #self.phi.printGrid()
                    # vel is used only to get the parent
                    adjustNumber( parts=self.pp, vel=self.vel, flags=self.flags, minParticles=minParticles, maxParticles=maxParticles, phi=self.phi, narrowBand=self.narrowBandWidth, exclude=self.phiObs ) 
                elif 0:
                    adjustNumber( parts=self.pp, vel=self.vel, flags=self.flags, minParticles=minParticles, maxParticles=maxParticles, phi=self.phi, exclude=self.phiObs ) 

            # mark int for measure
            if not self.b_fixed_vol:
                self.flags.mark_surface()

            # measure
            if 1:
                m = measure( self.pp, pVel, self.flags,  self.gravity, ppc, V0 )
                f_measure.write( f'{m[0]}\n' )
                f_measure.flush()

            # mesh
            if self.bMesh:
                phiMesh.copyFrom( self.phi )
                improvedParticleLevelset( self.pp, pindex, self.flags, gpi, phiMesh, radiusFactor, 1, 1 , 0.4, 3.5 ) # creates artifacts in dam flip05 128

                # mesh
                phiMesh.setBound( value=0., boundaryWidth=1 )
                phiMesh.createMesh( mesh )

            # print/write
            if 0:
                #self.flags.printGrid()
                self.vel.printGrid()
            
                self.pp.printParts()
                #self.pp.writeParticlesText( out + 'flipt_%04d.txt' % it )
            
            print( '(iteration) ', end='' )
            toc()

            # step
            print( '- step (%d)' % it )
            self.sol.step()
            #print( 'after step' )

            it += self.sol.timestep / self.dt
            it2 += 1
            if 0 or abs( it - round(it) ) < 1e-7:
                it = round( it )

                # screenshot
                if self.bScreenShot:
                    gui.screenshot( out + 'frame_%04d.png' % it ) # slow

                # save particle data
                if self.bSaveParts:
                    if self.bSaveUni:
                        # save particle data for flip03_gen.py surface generation scene
                        self.pp.save( out + 'parts_%04d.uni' % it )

                    # pdata fields must be before pp
                    objects = [ self.flags, self.phi, self.pp ]
                    #objects = [ self.pp ]
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
    if 1:
        limit_to_one_core()

    sim = Simulation()
    sim.main()
