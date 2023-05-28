
# flip5 in the comments refers to the scene script flip05_nbflip.py

import os, sys, math
import keyboard
from pathlib import Path

from manta import *

sys.path.append( r'c:\prj\python\\' )
from text_color import *
from tictoc import *

# print line number
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
        self.density = s.create(RealGrid, name='')
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
        
        # resample particles
        if self.resampleParticles:
            print( '  - resample particles' )
            gridParticleIndex(parts=pp, indexSys=pindex, flags=flags, index=gpi, counter=self.gCnt)
            #apicMapPartsToMAC(flags=flags, vel=vel, parts=pp, partVel=pVel, cpx=apic_pCx, cpy=apic_pCy, cpz=apic_pCz, mass=apic_mass)
            resampeOverfullCells(vel=vel, density=self.density, index=gpi, indexSys=pindex, part=pp, pVel=pVel, dt=s.timestep)
    
        # position solver
        print( '  - solve pressure due to density' )
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
        self.vel = sol.create( MACGrid, name='' )

    def create( self, sol ):
        self.center0 = self.center = Vec3( 0 )
        self.rad = 0
        self.vel_vec = Vec3( 0 )
        self.phi_init = sol.create( LevelsetGrid )
        self.hstart = 0
        self.hstop = 0
        self.skip = 0
        self.skip_last_it = 0 # the iteration where the last skip was made
        self.state = 0
        self.force = Vec3( 0 )
        self.increase_vel = 1
        self.stay = 0
        self.stay_last_it = 0
        self.mesh = sol.create( Mesh, name='mo_mesh' )
        self.file = None # loaded into maya by set_fluid.py
        self.shape = 0 # 0:box, 1:sphere
        self.part = sol.create( obs_particles )

class mesh_generator:
    def __init__( self, dim, gs, sol_main, narrowBand ):
        self.upres = 2 # 1, 2; scale resolution

        self.bScale = self.upres != 1
        self.narrowBand = narrowBand
        self.gs0 = gs

        if 1 and self.bScale:
            self.gs = self.upres*gs
            self.sol = Solver( name='gen_sol', gridSize=self.gs, dim=dim )
        else:
            self.sol = sol_main

        self.flags = self.sol.create(FlagGrid)
        self.phi = self.sol.create(LevelsetGrid)
        self.phiParts = self.sol.create(LevelsetGrid)
        self.pindex = self.sol.create(ParticleIndexSystem) 
        self.gpi = self.sol.create(IntGrid)
        self.mesh = sol_main.create( Mesh, name='mesh' ) # viewing a mesh from a different solver leads to a crash

        self.flags.initDomain( boundaryWidth=0 )

    def update_phi( self, phi ):
        if not self.narrowBand:
            return
        interpolateGrid( self.phi, phi )

    def generate( self, pp ):
        print( '- mesh_generator::generate()' )
        radiusFactor = 2.5 # 1, 2, 2.5

        if self.bScale:
            pp.transformPositions( self.gs0, self.gs )

        self.phi.setBound( value=0., boundaryWidth=1 )
        gridParticleIndex( parts=pp , flags=self.flags, indexSys=self.pindex, index=self.gpi )

        print( '  - union particle level sets' )
        # similar to flip03_gen.py
        #unionParticleLevelset( pp, self.pindex, self.flags, self.gpi, self.phiParts, radiusFactor )
        #averagedParticleLevelset( pp, self.pindex, self.flags, self.gpi, self.phiParts, radiusFactor , 1, 1 )
        improvedParticleLevelset( pp, self.pindex, self.flags, self.gpi, self.phiParts, radiusFactor, 1, 1, 0.4, 3.5 )

        if self.narrowBand:
            self.phi.addConst( 1. )
            self.phi.join( self.phiParts )
        else:
            self.phi.copyFrom( self.phiParts )

        self.phi.setBound( value=0., boundaryWidth=1 )
        self.phi.createMesh( self.mesh )

        if self.bScale:
            pp.transformPositions( self.gs, self.gs0 )
            self.mesh.scale( Vec3(1/self.upres) )

    def save( self, it ):
        fname = out + 'surface_%04d.bobj.gz' % it
        self.mesh.save( fname )

class simulation:
    def __init__( self ):
        # flags
        self.bScreenShot = 1
        self.bMesh       = 1
        self.bSaveMesh   = 1 # .bobj.gz
        self.bSaveVDB    = 0 # .vdb
        self.bSaveUni    = 0 # .uni
        if not self.bMesh:
            self.bSaveMesh = 0

        # params
        self.dim = 3 # 2, 3
        self.part_per_cell_1d = 2 # 3, 2(default), 1
        self.it_max = 2400 # 300, 500, 1200, 1400, 2400
        self.res = 50 # 32, 48/50, 64(default), 96/100, 128(large), 150, 250/256(, 512 is too large)

        self.b_fixed_vol = 1
        self.b_correct21 = 0

        self.narrowBand = bool( 1 )
        self.narrowBandWidth = 6 # 32:5, 64:6, 96:6, 128:8

        #self.gs = Vec3( self.res, self.res, 5 ) # debug thin 3D; at least z=5 if with obstacle (otherwise, it has 0 velocity?)
        self.gs = Vec3( self.res, int(1.5*self.res), self.res ) # tall tank
        #self.gs = Vec3( self.res, self.res, self.res )

        ###

        self.max_gs = max( [self.gs.x, self.gs.y, self.gs.z] )

        self.b2D = self.dim == 2

        if self.dim == 2:
            self.gs.z = 1
            self.bMesh = 0
            self.bSaveVDB = 0
            self.bSaveMesh = 0

        self.dt = .2 # .2(default), .5, 1(flip5, easier to debug)

        self.gravity = -0.1 # -0.1
        self.gravity *= math.sqrt( self.res )
        #self.gravity = -0.003 # flip5

        self.sol = Solver( name='sol', gridSize=self.gs, dim=self.dim )

        # automatic names are given only if the create is called from global context
        self.flags = self.sol.create(FlagGrid, name='flags')
        self.vel = self.sol.create(MACGrid, name='vel')
        self.pp = self.sol.create(BasicParticleSystem, name='pp')
        self.phiObs = self.sol.create(LevelsetGrid, name='')
        self.phi = []

        self.obs = moving_obstacle( self.sol )

        self.boundary_width = 0

        self.scene = {'type':0, 'name':'other'}

    def setup_scene( self ):
        #self.flags.initDomain( boundaryWidth=self.boundary_width ) 
        self.flags.initDomain( boundaryWidth=self.boundary_width, phiWalls=self.phiObs ) 

        if 0: # dam
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

            self.scene['type'] = 0
            self.scene['name'] = 'dam'

        elif 0: # falling drop
            fluidBasin = Box( parent=self.sol, p0=self.gs*Vec3(0,0,0), p1=self.gs*Vec3(1.0,0.1,1.0)) # basin
            dropCenter = Vec3(0.5,0.3,0.5)
            dropRadius = 0.1
            fluidDrop  = Sphere( parent=self.sol , center=self.gs*dropCenter, radius=self.res*dropRadius)
            self.phi = fluidBasin.computeLevelset()
            self.phi.join( fluidDrop.computeLevelset() ) # add drop
            self.flags.updateFromLevelset( self.phi )

            self.scene['type'] = 1
            self.scene['name'] = 'drop'

        elif 0: # basin
            # water
            fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0, 0.5, 0) ), p1=self.gs*( Vec3(1, 0.9, 1) ) )
            self.phi = fluidbox.computeLevelset()
            self.flags.updateFromLevelset( self.phi )

            # obstacle
            if 1:
                mesh = self.sol.create( Mesh, name='omesh' ) # need to switch to it in the gui to view

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
            # scale box height in gs
            # water
            fluid_h = 0.5 # 0.35, 0.5(default)
            fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0, 0., 0) ), p1=self.gs*( Vec3(1, fluid_h, 1) ) )
            print( f'- water level h={fluid_h}*res={fluid_h*self.gs.y}' )
            self.phi = fluidbox.computeLevelset()
            self.flags.updateFromLevelset( self.phi )

            # moving obstacle
            self.obs.exists = 1
            if self.obs.exists:
                self.obs.create( self.sol )
                self.obs.shape = 1 # box:0, sphere:1
                self.obs.rad = .08 if self.obs.shape == 0 else 0.1 # box:.08(default), .3, sphere:.1
                self.obs.rad *= self.res
                # shrink a bit if exactly cell size
                if abs( self.obs.rad - round(self.obs.rad) ) < 1e-7:
                    self.obs.rad *= 0.99
                self.obs.center0 = self.obs.center = self.gs*Vec3( 0.5, 1 - self.obs.rad/self.gs.y, 0.5 ) - Vec3( 0, 1, 0 ) # y:0.6, 0.95(default)

                self.obs.file = open( out + '_obstacle.txt', 'w' )
                self.obs.file.write( f'{self.obs.shape} {self.obs.rad/self.max_gs}\n' )
                self.obs.file.flush()

                span = 0.02 # 0, .01, .05
                fluid_h2 = fluid_h
                self.obs.hstart = fluid_h2*self.gs.y + 2 # must be soon enough to determine the impact speed with obs.vel when switching to state 1
                self.obs.hstop = self.obs.hstart - span*self.gs.y

                p0 = self.obs.center - Vec3(self.obs.rad)
                p1 = self.obs.center + Vec3(self.obs.rad)
                if self.dim == 2:
                    p0.z = p1.z = 0.5
                if self.obs.shape == 0:
                    shape = Box( parent=self.sol, p0=p0, p1=p1 )
                else:
                    shape = Sphere( parent=self.sol, center=self.obs.center, radius=self.obs.rad )
                self.obs.mesh.fromShape( shape )
                self.obs.mesh.save_pos()
                self.obs.mesh.set_color( Vec3( 0.5, 0.2, 0.2 ) )
                self.obs.mesh.set_2D( self.dim == 2 )

                self.obs.phi_init.copyFrom( self.phiObs )
                self.phiObs.join( shape.computeLevelset() )

                self.obs.force = Vec3( 0, self.gravity, 0 )
                self.obs.vel_vec = Vec3( 0, -0, 0 )
                self.obs.vel.setConst( self.obs.vel_vec )
                self.obs.vel.setBound( value=Vec3(0.), boundaryWidth=self.boundary_width+1 )

                tic()
                self.obs.part.create( self.obs.center, self.obs.rad, self.obs.shape, self.gs )
                toc()

                self.scene['type'] = 2 if self.obs.shape == 0 else 3
                self.scene['name'] = 'obs box' if self.obs.shape == 0 else 'obs ball'

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

        if 0 or self.dim == 2:
            set_print_2D( True )

        if 0: # sample 1D; requires changing sampleLevelsetWithParticles()
            ppc = self.part_per_cell_1d
        else:
            ppc = self.part_per_cell_1d**self.dim

        flags2 = self.sol.create(FlagGrid)
        vel2 = self.sol.create(MACGrid)
        velOld = self.sol.create(MACGrid)
        velParts = self.sol.create(MACGrid)
        pressure = self.sol.create(RealGrid)
        mapWeights = self.sol.create(MACGrid)

        volume = self.sol.create( RealGrid, name='volume' ) # volume measure illustration

        # add velocity data to particles
        pVel = self.pp.create(PdataVec3)
        pVel2 = self.pp.create(PdataVec3)
        phiParts = self.sol.create(LevelsetGrid)
        fractions = self.sol.create(MACGrid) # Sticks to walls and becomes slow without this? Not anymore. Shrinks an obstacle box and draws it inaccurately.

        bfs = self.sol.create(IntGrid) # discrete distance from surface

        # acceleration data for particle
        pindex = self.sol.create(ParticleIndexSystem)
        gpi = self.sol.create(IntGrid)

        correct21 = Correct21( self.dim, self.sol, self.part_per_cell_1d, self.pp )

        # info
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

        # setting file
        f_set = open( out + '_settings.txt', 'w' )
        c = self.gs
        f_set.write( '%d %d %d\n' % ( c.x, c.y, c.z ) ) # gs
        f_set.write( f"{self.scene['type']}\n" )
        f_set.flush()

        # create dir
        name = f'{self.res}^{self.dim}'

        if self.b_correct21:
            name += ' cor21'
        elif self.b_fixed_vol:
            if not self.narrowBand:
                name += ' full'
        else:
            name += ' flip'
        
        if self.narrowBand:
            name += f' band{self.narrowBandWidth}'

        name += f", {self.scene['name']}"

        os.mkdir( out + name )

        # size of particles 
        radiusFactor = 1.0

        np = self.pp.pySize()
        V0 = float( np ) / ppc
        #np_max = 2*np
        np_max = ppc * (self.res-2)**self.dim * 0.5
        print( f'# particles: {np}, np_max={np_max}' )

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

        # after initializing particles and before gui
        if self.bMesh:
            mesh_gen = mesh_generator( self.dim, self.gs, self.sol, self.narrowBand )
            mesh_gen.update_phi( self.phi )
            mesh_gen.generate( self.pp )

        if 1 and GUI:
            gui = Gui()
            for i in range( 2 ):
                gui.nextMeshDisplay() # 0:full, 1:hide, 2:x-ray
            gui.setRealGridDisplay( 0 ) # 0:none, 1:volume
            gui.setVec3GridDisplay( 0 ) # 0:none, 1:vel
            if 0 and self.dim == 3: # camera
                gui.setCamPos( 0, 0, -2.2 ) # drop
                gui.setCamRot( 35, -30, 0 )
            #gui.toggleHideGrids()
            gui.show()
            #gui.pause()
        else:
            self.bScreenShot = 0

        it = 0
        it2 = 0

        if self.bSaveMesh:
            mesh_gen.save( it )

        if self.bScreenShot:
            gui.screenshot( out + 'frame_%04d.png' % it ); # slow

        if self.bSaveUni:
            pressure.save( out + 'ref_parts_0000.uni' )
            self.pp.save( out + 'parts_%04d.uni' % it )

        if self.bSaveVDB:
            self.phi.set_name( 'phi' )
            objects = [ self.flags, self.phi, self.pp ] # need the 3 of them for volumetric .vdb and they need to be named (bifrost looks for a channel by name)
            #objects = [ self.pp ]
            fname = out + 'fluid_data_%04d.vdb' % it
            save( name=fname, objects=objects ) # error in debug mode "string too long?"

        # measure
        f_measure = open( out + '_measure.txt', 'w' )

        # loop
        ret = 0
        while 1:
            emphasize( '\n-----------------\n- time: %g(/%d; it2=%d)' % ( it, self.it_max, it2 ) )
            print( '- np=%d, np_max=%d' % ( self.pp.pySize(), np_max ) )

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
                frac = ( 1 - it % 1 )
                if 1 - frac < 1e-5:
                    frac = 1
                self.sol.timestep = frac * self.dt
            else: # adaptive for flip5
                self.sol.adaptTimestep( maxVel )

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
                if self.obs.state != 3:
                    if self.obs.center.y - self.obs.rad + self.obs.vel_vec.y*self.dt > 1.1: # move
                        print( '- obstacle still moves' )
                        if self.obs.increase_vel:
                            dv = self.sol.timestep * self.obs.force
                            self.obs.vel_vec += dv
                            max_y_speed = 7*self.gravity # 7, 10
                            if self.b_fixed_vol and self.obs.vel_vec.y < max_y_speed:
                                print( f'  - limiting speed to {max_y_speed}' )
                                self.obs.vel_vec.y = max_y_speed
                    else: # stay
                        print( '- obstacle reached the bottom' )
                        self.obs.vel_vec = Vec3( 0 )
                        self.obs.force = Vec3( 0 )
                        self.obs.state = 3
                else:
                    print( '- obstacle rests' )

                # obs.vel for boundary conditions
                if 1:
                    #obs_vel_vec2 = self.obs.vel_vec + dv # add some velocity in case it stopped--to remove remaining particles from the bottom
                    obs_vel_vec2 = self.obs.vel_vec + Vec3(0) # force copy
                    if self.obs.state == 1:
                        if self.dim == 3:
                            splash = 9 if self.obs.shape == 0 else 40 # box:9, sphere:40; splash size
                            obs_vel_vec2.y = splash*self.gravity 
                        else:
                            obs_vel_vec2.y = 3*self.gravity
                        print( '  - set obs.vel.y to {obs_vel_vec2.y} due to state 1' )
                    self.obs.vel.setConst( obs_vel_vec2 )
                    self.obs.vel.setBound( value=Vec3( 0 ), boundaryWidth=self.boundary_width + 1 )
                    #self.obs.vel.printGrid()
                elif self.obs.state < 2:
                    self.obs.state = 2

                # phiObs
                print( f'  - obs: center={self.obs.center}, rad={self.obs.rad}, vel_vec={self.obs.vel_vec}, vel.getMaxAbs={self.obs.vel.getMaxAbs()}' )
                p0 = self.obs.center - Vec3( self.obs.rad )
                p1 = self.obs.center + Vec3( self.obs.rad )
                if self.dim == 2:
                    p0.z = p1.z = 0.5
                if self.obs.shape == 0:
                    shape = Box( parent=self.sol, p0=p0, p1=p1 )
                else:
                    shape = Sphere( parent=self.sol, center=self.obs.center, radius=self.obs.rad )
                self.phiObs.copyFrom( self.obs.phi_init )
                self.phiObs.join( shape.computeLevelset() )
                #self.phiObs.printGrid()

                # obs mesh
                self.obs.mesh.load_pos()
                d = self.obs.center - self.obs.center0
                self.obs.mesh.offset( d )

                # obs flags
                if 1:
                    ret2 = mark_obstacle( flags=self.flags, obs=self.obs.part, center=self.obs.center )
                    if self.b_fixed_vol:
                        assert( ret2 )
                elif 1:
                    setObstacleFlags( flags=self.flags, phiObs=self.phiObs )
                else: # more precise
                    updateFractions( flags=self.flags, phiObs=self.phiObs, fractions=fractions, boundaryWidth=self.boundary_width )
                    setObstacleFlags( flags=self.flags, phiObs=self.phiObs, fractions=fractions )
                #self.flags.printGrid()

                # write pos
                c = self.obs.center
                c /= self.gs
                self.obs.file.write( '%g %g %g\n' % ( c.x, c.y, c.z ) )
                self.obs.file.flush()

            # emit
            if 0 and self.pp.pySize() < np_max:
                xi = self.gs * Vec3( 0.5, 0.9, 0.5 )
                v = Vec3( 0, -3.0, 0 ) # -3
                for i in range(-1, 2):
                    for j in range(-1, 2):
                        if self.pp.pySize() >= np_max:
                            break
                        if self.dim == 2:
                            j = 0
                        emit_particles( self.pp, pVel, self.flags, self.part_per_cell_1d, xi + Vec3(i, 0, j), v )
                        if self.dim == 2:
                            break
                V0 = float( self.pp.pySize() ) / ppc # update volume

            # update flags; there's also flags.updateFromLevelset()
            if not self.b_fixed_vol or it == 0:
                #self.flags.printGrid()
                print( '- markFluidCells (update flags)' )
                markFluidCells( parts=self.pp, flags=self.flags ) # marks deep in narrowBand as empty; better for a moving obstacle?
                #markFluidCells( parts=self.pp, flags=self.flags, phiObs=self.phiObs )
                if self.narrowBand:
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
            setWallBcs( flags=self.flags, vel=self.vel, obvel=self.obs.vel ) # calls KnSetWallBcs; default
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
            alpha = .1 # 0, .1
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
            print( f'  - obs_naive={obs_naive}' )
            if self.b_fixed_vol:
                self.phi.setBoundNeumann( 0 ) # make sure no new particles are placed at outer boundary
                #self.phi.printGrid()

                pVel.setSource( self.vel, isMAC=True ) # set source grid for resampling, used in insertBufferedParticles()

                dt_bound = 0
                #dt_bound = dt/4
                #dt_bound = self.sol.timestep/4
                #dt_bound = max( dt_bound, dt/4 )

                obs_vel_vec3 = Vec3(0) if obs_naive else self.obs.vel_vec
                obsp = None
                if self.obs.exists:
                    obsp = self.obs.part
                    obsp.update_center( self.obs.center )

                # obs_vel: it modifies it to either one or zero cell distance, staying in place and losing velocity (unlike particles)
                
                ret2 = fixed_volume_advection( pp=self.pp, pVel=pVel, flags=self.flags, dt=self.sol.timestep, dt_bound=dt_bound, dim=self.dim, ppc=ppc, phi=self.phi, it=it2, use_band=self.narrowBand, band_width=self.narrowBandWidth, bfs=bfs, obs=obsp, obs_vel=obs_vel_vec3 )

                if not ret2:
                    ret = -1
                else:
                    self.sol.timestep = ret2[0]
                    if self.sol.timestep < 0:
                        ret = -1
                        self.sol.timestep = abs( self.sol.timestep )
                    if not obs_naive:
                        #self.obs.vel_vec = Vec3( ret2[1], ret2[2], ret2[3] )
                        obs_stop = ret2[1]

                # if using band
                if 0 and self.narrowBand:
                    include_walls = true
            #self.flags.printGrid()

            # update obstacle
            if self.obs.exists:
                # limit dt to one-cell movement
                while 1:
                    obs_center2 = self.obs.center + self.sol.timestep * self.obs.vel_vec
                    if obs_center2.y < 1.1 + self.obs.rad:
                        obs_center2.y = 1.1 + self.obs.rad
                    if int( self.obs.center.y - self.obs.rad ) - int( obs_center2.y - self.obs.rad ) <= 1:
                        break
                    self.sol.timestep /= 2
                    print( f'  - halving the time step: {self.sol.timestep}' )

                # test obstacle position
                print( '  - obs_stop=%d' % obs_stop )
                if 1 and not obs_stop:
                    flags2.copyFrom( self.flags )
                    if not mark_obstacle( flags=flags2, obs=self.obs.part, center=obs_center2 ):
                        emphasize( '  - obstacle position is invalid; stopping the obstacle' )
                        assert( not self.b_fixed_vol or obs_naive )
                        obs_stop = 1

                self.obs.increase_vel = 1
                print( f'  - obs: .state={self.obs.state}, .vel_vec={self.obs.vel_vec}, dt={self.sol.timestep}, .center={self.obs.center}, .force={self.obs.force}, .skip={self.obs.skip}, .increase_vel={self.obs.increase_vel}, obs_center2={obs_center2}, stay={self.obs.stay}, hstart={self.obs.hstart}, hstop={self.obs.hstop}' )
                if self.obs.state != 3:
                    if self.obs.state == 0 and int( obs_center2.y - self.obs.rad ) == int( self.obs.center.y - self.obs.rad ) and pyNorm( self.obs.vel_vec ) > 0: # state 0 and no grid movement
                        self.obs.center = obs_center2
                    else:
                        if self.obs.state == 1 or not obs_stop:
                            self.obs.stay = 0
                            if self.obs.state < 2 and self.obs.hstop <= self.obs.center.y - self.obs.rad <= self.obs.hstart:
                                if self.obs.state == 0:
                                    self.obs.state = 1
                                    emphasize2( f'  - new obs.state: {self.obs.state}' )
                                # wait by skipping grid progress
                                # count progress by int(it) rather than it2
                                if int(it) > self.obs.skip_last_it:
                                    self.obs.skip_last_it = int(it)
                                    self.obs.skip += 1
                                n_skips = 0 # 0, 1, 2(default), 3, 5; how many steps to skip
                                if 0 and self.dim == 2:
                                    n_skips = min( n_skips, 1 )
                                if self.obs.skip >= n_skips: 
                                    #self.obs.skip = 0
                                    if 1:
                                        if 0:
                                            self.obs.state = 2
                                            emphasize2( f'  - new obs.state: {self.obs.state}' )
                                        if 1:
                                            # speed in water
                                            self.obs.vel_vec.y = self.gravity
                                    if 1:
                                        self.obs.center = obs_center2
                                else:
                                    self.obs.vel_vec.y = 0
                                    emphasize2( f'  - skip {self.obs.skip}/{n_skips}' )
                                self.obs.force.y = 0
                            else:
                                if 1 and self.obs.state == 1:
                                    self.obs.state = 2
                                    emphasize2( f'  - new obs.state: {self.obs.state}' )
                                    if 0:
                                        self.obs.vel_vec.y = self.gravity/2 # 1, 2, 6
                                        #self.obs.vel_vec /= 4
                                        #self.obs.vel_vec = Vec3(0)

                                        self.obs.force = Vec3(0)
                                        #self.obs.force /= 4
                                self.obs.center = obs_center2
                        else:
                            self.obs.increase_vel = 0
                            if int(it) > self.obs.stay_last_it:
                                self.obs.stay_last_it = int(it)
                                self.obs.stay += 1
                            if 1 and self.obs.stay > 50: # 50
                                print( f'  - no obstacle progress (stay={self.obs.stay}); stopping it' )
                                self.obs.vel_vec = Vec3( 0 )
                                self.obs.force = Vec3( 0 )
                                self.obs.state = 3
                                emphasize2( f'  - new obs.state: {self.obs.state}' )

            # correct21
            if self.b_correct21:
                correct21.main( self.sol, self.flags, self.pp, self.vel, pindex, gpi, self.phiObs )
            
            # for narrowBand, before updating phi with the particles
            if self.bMesh:
                mesh_gen.update_phi( self.phi )

            # create level set from particles
            if 1:
                gridParticleIndex( parts=self.pp, flags=self.flags, indexSys=pindex, index=gpi )
                unionParticleLevelset( self.pp, pindex, self.flags, gpi, phiParts, radiusFactor ) 
                if self.narrowBand:
                    # combine level set of particles with grid level set
                    self.phi.addConst( 1. ) # shrink slightly
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

            # update and mark surface for measure
            if not self.b_fixed_vol:
                print( '- markFluidCells (update flags)' )
                markFluidCells( parts=self.pp, flags=self.flags )
                self.flags.mark_surface()

            toc()

            # measure
            if 1:
                m = measure( self.pp, pVel, self.flags, self.gravity, ppc, V0, volume )
                f_measure.write( f'{m[0]}\n' )
                f_measure.flush()
                #self.flags.printGrid()
                #volume.printGrid()

            # mesh
            if self.bMesh:
                tic()
                mesh_gen.generate( self.pp )
                toc()

            # print
            if 0:
                self.flags.printGrid()
                #self.vel.printGrid()
            
                #self.pp.printParts()
                #self.pp.writeParticlesText( out + 'particles_%04d.txt' % it )
            
            print( '(iteration) ', end='' )

            # step; updates gui and when pause takes place
            print( '- step (%g)' % it )
            self.sol.step( int(it) )
            #print( 'after step' )

            # it
            it += self.sol.timestep / self.dt
            it2 += 1

            # save
            if 0 or abs( it - round(it) ) < 1e-7:
                it = round( it )

                # screenshot
                if self.bScreenShot:
                    fname = out + 'frame_%04d.png' % int(it)
                    print( f'- saving: {fname}' )
                    gui.screenshot( fname ) # slow

                # data
                if self.bSaveUni:
                    # save particle data for flip03_gen.py surface generation scene
                    self.pp.save( out + 'parts_%04d.uni' % it )

                if self.bSaveVDB:
                    # pdata fields must be before pp
                    objects = [ self.flags, self.phi, self.pp ]
                    #objects = [ self.pp ]
                    save( name=out + 'fluid_data_%04d.vdb' % it, objects=objects )

                if self.bSaveMesh:
                    mesh_gen.save( it )
                
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

    # del
    for path in Path( out ).glob("*"):
        if path.is_file():
            path.unlink()
        elif path.is_dir():
            path.rmdir()

    # copy
    os.system( 'cp %s../video.bat %s' % (out, out) )

    # (debug) for consistent result; for large res, the step() hangs?
    if 0:
        limit_to_one_core()

    # init matlab
    cmd = "%%close all;\n clear classes; clear java; dbclear all; clear all; pack; jheapcl(0); set(0, 'DefaultFigureWindowState', 'minimized');" + " cd c:/prj/test_data/relative/_tmp;" + " addpath( 'c:/prj/mantaflow_mod' );"
    matlab_eval( cmd )
            
    # simulation
    sim = simulation()
    sim.main()
