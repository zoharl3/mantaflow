
# flip5 in the comments refers to the scene script flip05_nbflip.py

import os, sys, math, shutil
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

def init_matlab():
    cmd = "%%close all;\n clear classes; clear java; dbclear all; clear all; pack; jheapcl(0); set(0, 'DefaultFigureWindowState', 'minimized');" + " cd c:/prj/test_data/relative/_tmp;" + " addpath( 'c:/prj/mantaflow_mod/matlab', 'c:/prj/fluid/matlab' );"
    matlab_eval( cmd )

def toVec3( c ):
    assert( len(c) == 3 )
    return Vec3( c[0], c[1], c[2] )

def test_MAC():
    sol = Solver( gridSize=Vec3(2, 2, 1), dim=2 )
    v = sol.create( MACGrid )

    v.set( 0, 0, 0, Vec3(1, 1, 0) )
    v.set( 0, 1, 0, Vec3(2, 2, 0) )
    v.set( 1, 0, 0, Vec3(3, 3, 0) )
    v.set( 1, 1, 0, Vec3(4, 4, 0) )

    v.printGrid()

    print( f'v.getInterpolated(0.5, 0.5)={v.getInterpolated(Vec3(0.5, 0.5, 0.5))}' )
    print( f'v.getCentered(0, 0)={v.getCentered(Vec3(0, 0, 0))}' ) # the single cell that the 2x2 grid defines
    #print( f'v.getAtMACX(0, 0)={v.getAtMACX(0, 0, 0)}' ) # boundary can't be calculated, asserts in debug
    print( f'v.getAtMACX(1, 0)={v.getAtMACX(1, 0, 0)}' ) # at cell [1, 0], x-component is the one from the cell, y-component is interpolated with previous cell
    #print( f'v.getAtMACX(1, 1)={v.getAtMACX(1, 1, 0)}' ) # boundary
    print( f'v.getAtMACY(0, 1)={v.getAtMACY(0, 1, 0)}' )

# correct19 (position solver, Thuerey19)
# The band requires fixing, probably identifying non-band fluid cells as full. In the paper, it's listed as future work.
class Correct19:
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
        print( '- Correct19.main()' )
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
        solvePressureSystem( rhs=self.density, vel=vel, pressure=self.Lambda, flags=self.flagsPos, cgAccuracy=1e-3 )
        computeDeltaX( deltaX=self.deltaX, Lambda=self.Lambda, flags=self.flagsPos )
        max_deltaX = self.deltaX.getMaxAbs()
        print( f'    - max_deltaX.MaxAbs={max_deltaX}' )
        if 1 and max_deltaX > 10: # 10
            if 0:
                warn( 'deltaX blew up; skipping correction' )
                return
            else:
                warn( 'deltaX blew up; not handling' )
        mapMACToPartPositions( flags=self.flagsPos, deltaX=self.deltaX, parts=pp, dt=sol.timestep )
        
        # print
        if 0:
            #flags.printGrid()
            #self.flagsPos.printGrid()
            self.density.printGrid()
            self.Lambda.printGrid()
            self.deltaX.printGrid()

class moving_obstacle:
    def __init__( self, sim ):
        self.exists = 0
        self.sim = sim
        self.sol = self.sim.sol
        self.vel = self.sol.create( MACGrid, name='' )

    def create( self ):
        self.exists = 1

        self.center0 = self.center = Vec3( 0 )
        self.rad = 0 # radius (at least) in the y-axis
        self.rad3 = Vec3( 0 )
        self.vel_vec = Vec3( 0 )
        self.terminal_speed = 5*self.sim.gravity/3 # 1, 3; terminal velocity

        self.phi_init = self.sol.create( LevelsetGrid )

        self.start_h = 0
        self.stop_h = 0

        self.state = 0
        self.force = Vec3( 0 )
        self.increase_vel = 1
        self.n_stay = 0
        self.stay_last_it = 0

        self.mesh = self.sol.create( Mesh, name='mov_obs_mesh' )
        self.file = None # loaded into maya by set_fluid.py
        self.shape = 0 # 0:box, 1:sphere

        self.part = self.sol.create( obs_particles )

    def init( self, shape ):
        # mesh
        self.mesh.fromShape( shape )
        self.mesh.save_pos()
        self.mesh.set_color( Vec3( 0.5, 0.2, 0.2 ) )
        self.mesh.set_2D( self.sol.is2D() )

        # phi
        self.phi_init.copyFrom( self.sim.phiObs )
        self.sim.phiObs.join( shape.computeLevelset() )

        # vel
        self.vel.setConst( self.vel_vec )
        self.vel.setBound( value=Vec3(0.), boundaryWidth=self.sim.boundary_width+1 )

        # fill obstacle with particles
        tic( 'obs part' )
        self.part.create( self.center, self.rad3, self.shape, self.sim.gs )
        toc()

    def write_pos( self ):
        if self.exists and self.file:
            c = self.center
            c /= self.sim.gs
            self.file.write( '%g %g %g\n' % ( c.x, c.y, c.z ) )
            self.file.flush()

class static_obstacle:
    def __init__( self, sol ):
        self.sol = sol
        self.exists = 0

    def create( self ):
        self.exists = 1
        self.mesh = self.sol.create( Mesh, name='static_obs_mesh' ) # need to switch to it in the gui to view

        self.part = None

        self.vel = self.sol.create( MACGrid )
        self.flags = self.sol.create( FlagGrid )

    # only for this obstacle cells
    def set_wall_bcs( self , flags, vel ):
        if not self.part:
            return
        self.flags.copyFrom( flags )
        self.flags.clear_obstacle()
        mark_obstacle( flags=self.flags, obs=self.part, center=Vec3(0) )
        #setWallBcs( flags=self.flags, vel=vel, obvel=self.vel )
        set_wall_bcs2( flags=self.flags, vel=vel, obvel=self.vel )

class mesh_generator:
    def __init__( self, dim, gs, sol_main, narrowBand, out_dir, upres=2 ):
        self.upres = upres # 1, 2; scale resolution

        self.union_method = 2
        self.bScale = self.upres != 1
        self.narrowBand = narrowBand
        self.gs0 = gs
        self.out_dir = out_dir

        if 1 and self.bScale:
            self.gs = self.upres*gs
            self.sol = Solver( name='gen_sol', gridSize=self.gs, dim=dim )
        else:
            self.sol = sol_main

        self.flags = self.sol.create( FlagGrid )
        self.phi = self.sol.create( LevelsetGrid )
        self.phiParts = self.sol.create( LevelsetGrid )
        self.pindex = self.sol.create( ParticleIndexSystem )
        self.gpi = self.sol.create( IntGrid )
        self.mesh = sol_main.create( Mesh, name='mesh' ) # viewing a mesh from a different solver leads to a crash

        self.flags.initDomain( boundaryWidth=0 )

    def update_phi( self, phi ):
        if not self.narrowBand:
            return
        interpolateGrid( self.phi, phi )

    def generate( self, pp ):
        print( '- mesh_generator::generate()' )
        radiusFactor = 2.5 # 1, 2, (default)2.5

        if self.bScale:
            pp.transformPositions( self.gs0, self.gs )

        self.phi.setBound( value=0., boundaryWidth=1 )
        gridParticleIndex( parts=pp , flags=self.flags, indexSys=self.pindex, index=self.gpi )

        print( f'  - union particle level sets (union_method={self.union_method})' )
        # similar to flip03_gen.py
        if self.union_method == 0:
            unionParticleLevelset( pp, self.pindex, self.flags, self.gpi, self.phiParts, radiusFactor )
        elif self.union_method == 1:
            averagedParticleLevelset( pp, self.pindex, self.flags, self.gpi, self.phiParts, radiusFactor , 1, 1 )
        elif self.union_method == 2:
            improvedParticleLevelset( pp, self.pindex, self.flags, self.gpi, self.phiParts, radiusFactor, 1, 1, 0.4, 3.5 )
        else:
            assert( 0 )

        if self.narrowBand:
            self.phi.addConst( 1. )
            self.phi.join( self.phiParts )
        else:
            self.phi.copyFrom( self.phiParts )

        self.phi.setBound( value=0., boundaryWidth=1 )
        self.phi.createMesh( self.mesh )

        if self.bScale:
            pp.transformPositions( self.gs, self.gs0 )
            self.mesh.scale( Vec3(1.0/self.upres) )

    def save( self, it ):
        fname = self.out_dir + 'surface_%04d.bobj.gz' % it
        self.mesh.save( fname )

# methods
# 0       1          2          3         4
FLIP, FIXED_VOL, CORRECT19, DE_GOES22, MATLAB_FLIP = range( 5 )

class simulation:
    def __init__( self, method ):
        # flags
        self.bScreenShot  = 1
        self.b_fluid_mesh = 1 # generate mesh for fluid
        self.bSaveMesh    = 1 # .bobj.gz
        self.bSaveVDB     = 0 # .vdb
        self.bSaveUni     = 0 # .uni
        if not self.b_fluid_mesh:
            self.bSaveMesh = 0

        # params
        self.part_per_cell_1d = 1 # 1, 2(default), 3
        self.dim = 2 # 2, 3
        self.it_max = 1500 # 300, 500, 1000, 1500, 2500
        self.res = 100 # 32, 48/50, 64, 96/100, 128(large), 150, 250/256(, 512 is too large)

        self.narrowBand = bool( 0 ) # there's an override in main() for some methods
        self.narrowBandWidth = 3 # 3(default,large obs), 6(dam)
        self.inter_control_method = 3 # BAND_INTERFACE_CONTROL_METHOD: full=0, one-sided=1, revert=2, push=3

        self.large_obs = 0
        self.obs_shape = 0 # box:0, sphere:1
        self.b_test_collision_detection = 1 # enable naive test of collision detection for other methods

        if 0: # tall tank
            #self.gs = Vec3( self.res, self.res, 5 ) # debug thin 3D; at least z=5 if with obstacle (otherwise, it has 0 velocity?)
            self.gs = Vec3( self.res, int(1.5*self.res), self.res ) # tall tank
        else: # square tank
            self.gs = Vec3( self.res, self.res, self.res ) # iso

        ###

        self.b2D = self.dim == 2

        # method
        self.method = method

        self.splash = 1
        # disable splash (else it's too fast for flip)
        if self.large_obs:
            self.splash = 0

        if not self.narrowBand:
            self.narrowBandWidth = -1

        self.max_gs = max( [self.gs.x, self.gs.y, self.gs.z] )

        if 0: # sample in 1D; it also requires changing sampleLevelsetWithParticles()
            self.ppc = self.part_per_cell_1d
        else:
            self.ppc = self.part_per_cell_1d**self.dim

        if self.dim == 2:
            self.gs.z = 1
            self.b_fluid_mesh = 0
            self.bSaveVDB = 0
            self.bSaveMesh = 0

        self.dt = .2 # .2(default), .5, 1(flip5, easier to debug)

        self.gravity = -0.02 # -0.02
        self.gravity *= math.sqrt( self.res )
        #self.gravity = -0.003 # flip5

        self.sol = Solver( name='sol', gridSize=self.gs, dim=self.dim )

        # automatic names are given to manta c++ objects only if the create is called from global context
        self.flags = self.sol.create( FlagGrid, name='flags' )
        self.flags2 = self.sol.create( FlagGrid)
        self.vel = self.sol.create( MACGrid, name='vel' )
        self.pp = self.sol.create( BasicParticleSystem, name='pp' )
        self.phiObs = self.sol.create( LevelsetGrid, name='' )
        self.phi = None

        self.obs = moving_obstacle( self )
        self.obs2 = static_obstacle( self.sol )

        self.boundary_width = 0

        self.scene = {'type':0, 'name':'other', 'cam':1}

        self.out_dir = None

    def setup_scene( self ):
        #self.flags.initDomain( boundaryWidth=self.boundary_width ) 
        self.flags.initDomain( boundaryWidth=self.boundary_width, phiWalls=self.phiObs ) 

        if 1: # dam
            # my dam
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
            self.flags.updateFromLevelset( self.phi ) # update fluid flags

            self.scene['type'] = 0
            self.scene['name'] = 'dam'
            self.scene['cam'] = 3

        elif 0: # water drop
            fluidBasin = Box( parent=self.sol, p0=self.gs*Vec3(0,0,0), p1=self.gs*Vec3(1.0,0.1,1.0)) # basin
            dropCenter = Vec3(0.5,0.3,0.5)
            dropRadius = 0.1
            fluidDrop  = Sphere( parent=self.sol , center=self.gs*dropCenter, radius=self.res*dropRadius)
            self.phi = fluidBasin.computeLevelset()
            self.phi.join( fluidDrop.computeLevelset() ) # add drop
            self.flags.updateFromLevelset( self.phi )

            self.scene['type'] = 1
            self.scene['name'] = 'drop'
            self.scene['cam'] = 3

        elif 0: # basin--not selected
            # water
            fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0, 0.6, 0) ), p1=self.gs*( Vec3(1, 0.9, 1) ) )
            self.phi = fluidbox.computeLevelset()
            self.flags.updateFromLevelset( self.phi )

            # obstacle
            if 1:
                self.obs2.create( self.sol )
                self.obs2.mesh.load( r'c:\prj\mantaflow_mod\resources\funnel.obj' )
                s = Vec3(self.res) # the scale needs to be in all axes
                if self.b2D:
                    s.z = 2
                self.obs2.mesh.scale( s )
                self.obs2.mesh.offset( self.gs * Vec3(0.5, 0.3, 0.5) )

                mesh_phi = self.sol.create( LevelsetGrid )
                self.obs2.mesh.computeLevelset( mesh_phi, 2 ) # uses normals, thus a smooth mesh is better
                #mesh_phi.printGrid()
                #self.phiObs.printGrid()
                self.phiObs.join( mesh_phi )
                self.phi.subtract( self.phiObs ) # not to sample particles inside obstacle

        elif 0: # falling obstacle
            # water
            fluid_h = 0.5 # 0.5(default)
            if self.large_obs:
                fluid_h = 0.2 # 0.2(large box)
            fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0, 0, 0) ), p1=self.gs*( Vec3(1, fluid_h, 1) ) )
            print( f'- water level h={fluid_h}*res={fluid_h*self.gs.y}' )
            self.phi = fluidbox.computeLevelset()
            self.flags.updateFromLevelset( self.phi )

            # moving obstacle
            if 1:
                self.obs.create()
                self.obs.shape = self.obs_shape
                
                # rad 
                if self.obs.shape == 0: # box:.08(default), sphere:.1
                    self.obs.rad = 0.08
                else: 
                    self.obs.rad = 0.15 # 0.1, 0.15
                if self.large_obs: # large
                    # margin equals to boundary (1 cell) + side path (between obs and tank)
                    #margin = 2 / self.res # one-cell wide side path; for hi-res results in no progress
                    #margin = 2 / 50 # fixed margin (and obs width), one-cell wide side path for res50; hi-res is faster than res50
                    margin = ( 1 + self.res / 50 ) / self.res # side path changes with res
                    self.obs.rad = 0.5 - margin
                self.obs.rad *= self.res
                # shrink a bit if exactly cell size
                if abs( self.obs.rad - round(self.obs.rad) ) < 1e-7:
                    self.obs.rad *= 0.99
                self.obs.rad3 = Vec3( self.obs.rad )

                # center0
                self.obs.center0 = self.obs.center = self.gs*Vec3( 0.5, 1 - self.obs.rad/self.gs.y, 0.5 ) - Vec3( 0, 1, 0 ) # start from the ceiling
                if 0 and self.b2D:
                    self.obs.center0 = self.obs.center = self.gs*Vec3( 0.5, ( 1 + fluid_h )/2, 0.5 ) # middle of the air
                if 0:
                    self.obs.center0 = self.obs.center = self.gs*Vec3( 0.5, 0.02 + fluid_h + self.obs.rad/self.gs.y, 0.5 ) # near the surface

                fluid_h2 = fluid_h
                self.obs.start_h = fluid_h2*self.gs.y + 2 # must be soon enough to determine the impact speed with obs.vel when switching to state 1
                print( '- start_h={self.obs.start_h}' )

                # shape
                if self.obs.shape == 0:
                    p0 = self.obs.center - Vec3(self.obs.rad)
                    p1 = self.obs.center + Vec3(self.obs.rad)
                    if self.dim == 2:
                        p0.z = p1.z = 0.5
                    shape = Box( parent=self.sol, p0=p0, p1=p1 )
                else:
                    shape = Sphere( parent=self.sol, center=self.obs.center, radius=self.obs.rad )

                # init force and velocity
                if 1:
                    self.obs.force = Vec3( 0, 5*self.gravity, 0 )
                    self.obs.vel_vec = Vec3( 0, -0, 0 )
                else:
                    self.obs.force = Vec3( 0, 0, 0 )
                    self.obs.vel_vec = Vec3( 0, 5*self.gravity, 0 )
                    self.obs.state = 4 # constant speed
                self.obs.init( shape )

                self.scene['type'] = 2 + int( self.large_obs )
                self.scene['name'] = 'obs box' if self.obs.shape == 0 else 'obs ball'
                if self.large_obs:
                    self.scene['name'] = 'large ' + self.scene['name']
                self.scene['cam'] = 2

        elif 0: # fluid on top of obstacle--trying to reproduce the surface extraction bug
            # water
            fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0, 0.7, 0) ), p1=self.gs*( Vec3(1, 0.9, 1) ) )
            self.phi = fluidbox.computeLevelset()
            self.flags.updateFromLevelset( self.phi )

            # obstacle
            if 1:
                self.obs.create()
                self.obs.rad = 0.5 - 1/self.res
                self.obs.rad *= self.res
                self.obs.rad3 = Vec3( self.obs.rad )

                # center0
                self.obs.center0 = self.obs.center = self.gs*Vec3( 0.5, self.obs.rad/self.gs.y, 0.5 ) + Vec3( 0, 1, 0 ) # on the floor

                # shape
                p0 = self.obs.center - Vec3(self.obs.rad)
                p1 = self.obs.center + Vec3(self.obs.rad)
                if self.dim == 2:
                    p0.z = p1.z = 0.5
                shape = Box( parent=self.sol, p0=p0, p1=p1 )
                self.obs.init( shape )

        elif 0: # compress
            # water
            h = 0.4
            fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0.5, 0, 0) ), p1=self.gs*( Vec3(1, h, 1) ) )
            self.phi = fluidbox.computeLevelset()
            self.flags.updateFromLevelset( self.phi )

            # moving obstacle
            if 1:
                self.obs.create()
                self.obs.shape = 0
                self.obs.rad = .1*self.res
                self.obs.rad3 = self.res*Vec3( 0.5, self.obs.rad/self.res, 0.5 ) - Vec3( 1., 0, 0 ) # 1, 1.6
                #self.obs.rad3.x -= 1.01 # space
                self.obs.center0 = self.obs.center = self.gs*Vec3( 0.5, h + (2 + self.obs.rad) /self.res, 0.5 ) - Vec3( 0., 0, 0 ) # 0, 0.6

                p0 = self.obs.center - self.obs.rad3
                p1 = self.obs.center + self.obs.rad3
                if self.dim == 2:
                    p0.z = p1.z = 0.5
                shape = Box( parent=self.sol, p0=p0, p1=p1 )

                # init
                self.obs.force = Vec3( 0, 0, 0 )
                self.obs.vel_vec = Vec3( 0, 5*self.gravity*.5, 0 )
                self.obs.state = 4 # constant speed
                self.obs.init( shape )

            #self.scene['type'] = 4
            self.scene['name'] = 'compress'

        elif 0: # spiral
            # see if speed_limit needs to be changed after solving for pressure
            # water
            fluidbox = Box( parent=self.sol, p0=self.gs*( Vec3(0.5, 0, 0) ), p1=self.gs*( Vec3(1, 0.7, 1) ) ) # tubes:0.6, spiral:0.4
            self.phi = fluidbox.computeLevelset()
            self.flags.updateFromLevelset( self.phi )

            # static obstacle
            if 1:
                self.obs2.create()
                #self.obs2.mesh.set_name( '' )
                self.obs2.mesh.load( r'c:\prj\mantaflow_mod\resources\spiral.obj' )
                s = Vec3( self.res )
                s.x *= 0.7
                if self.b2D:
                    s.z = 4
                self.obs2.mesh.scale( s )
                if self.b2D:
                    self.obs2.mesh.offset( Vec3( 0., 0, -s.z/2 ) )

                mesh_phi = self.sol.create( LevelsetGrid )
                self.obs2.mesh.computeLevelset( mesh_phi, 2 )
                #mesh_phi.printGrid()
                #self.phiObs.printGrid()
                self.phiObs.join( mesh_phi )
                self.phi.subtract( self.phiObs )

                # obs particles
                self.obs2.part = self.sol.create( obs_particles )
                tic( 'obs part' )
                #self.obs2.part.create_from_levelset( mesh_phi )
                self.obs2.part.create_from_mesh( self.obs2.mesh )
                toc()

            # moving obstacle
            if 1:
                self.obs.create()
                self.obs.shape = 0
                self.obs.rad = .9
                left_y = 0.35 - 0.5/self.res
                self.obs.rad3 = self.res*Vec3( (1 - left_y)/2, self.obs.rad/self.res, 0.5 )
                self.obs.rad3 -= Vec3( 1, 0, 0 ) # .5, 2
                self.obs.center0 = self.obs.center = self.gs*Vec3( (1 + left_y)/2, 0.9, 0.5 ) 

                p0 = self.obs.center - self.obs.rad3
                p1 = self.obs.center + self.obs.rad3
                if self.dim == 2:
                    p0.z = p1.z = 0.5
                shape = Box( parent=self.sol, p0=p0, p1=p1 )

                # init
                self.obs.force = Vec3( 0, 0, 0 )
                self.obs.vel_vec = Vec3( 0, 5*self.gravity*1, 0 )
                self.obs.state = 4 # constant speed
                self.obs.init( shape )

                #self.scene['type'] = 5
                self.scene['name'] = 'spiral'

        # common
        #self.phiObs.printGrid()
        #self.flags.printGrid()

        sampleLevelsetWithParticles( phi=self.phi, flags=self.flags, parts=self.pp, discretization=self.part_per_cell_1d, randomness=0.1 ) # 0.05, 0.1, 0.2

        # also sets boundary flags for phiObs (and shrinks it)
        if self.obs.exists or self.obs2.exists:
            if 1: # imprecise
                setObstacleFlags( flags=self.flags, phiObs=self.phiObs )
            else:
                updateFractions( flags=self.flags, phiObs=self.phiObs, fractions=fractions, boundaryWidth=self.boundary_width )
                setObstacleFlags( flags=self.flags, phiObs=self.phiObs, fractions=fractions )
        #self.phi.printGrid()

    # fixed volume (my scheme)
    def fixed_volume( self, pVel, obs_naive, include_walls, ret, it2, bfs, stat ):
        self.phi.setBoundNeumann( 0 ) # make sure no new particles are placed at outer boundary
        #self.phi.printGrid()

        pVel.setSource( self.vel, isMAC=True ) # set source grid for resampling, used in insertBufferedParticles()

        obs_vel_vec3 = Vec3(0) 
        if not obs_naive and self.obs.exists:
            obs_vel_vec3 = self.obs.vel_vec
        obs_part = None
        if self.obs.exists:
            obs_part = self.obs.part
            obs_part.update_center( self.obs.center )

        ret2 = fixed_volume_advection( pp=self.pp, pVel=pVel, flags=self.flags, dt=self.sol.timestep, dim=self.dim, ppc=self.ppc, phi=self.phi, bfs=bfs, it=it2, use_band=self.narrowBand, band_width=self.narrowBandWidth, inter_control_method=self.inter_control_method, obs=obs_part, obs_vel=obs_vel_vec3 )

        obs_stop = 0
        if not ret2:
            ret = -1
        else:
            if ret2[0] < 0:
                ret = -1
            if not obs_naive:
                obs_stop = ret2[1]
            stat['opt_time'] = ret2[2]
            stat['push_time'] = ret2[3]
            stat['avg_num_particles'] = ret2[4]
            stat['avg_nnz'] = ret2[5]

        # if using band
        if 0 and self.narrowBand:
            include_walls = true

        return [ ret, obs_stop, include_walls ]

    def update_obstacle( self, obs_naive, obs_stop, it ):
        print( '- update_obstacle()' )
        # limit dt to one-cell movement
        while 1:
            obs_center2 = self.obs.center + self.sol.timestep * self.obs.vel_vec
            if obs_center2.y < 1.1 + self.obs.rad:
                obs_center2.y = 1.1 + self.obs.rad
            if int( self.obs.center.y - self.obs.rad ) - int( obs_center2.y - self.obs.rad ) <= 1:
                break
            self.sol.timestep /= 2
            print( f'  - halving the time step: {self.sol.timestep}' )

        # collision detection: test obstacle position
        print( '  - obs_stop=%d' % obs_stop )
        if self.b_test_collision_detection and not obs_stop: # (if disabled for flip, then you may want to disable pushOutofObs--it's already disabled by default)
            self.flags2.copyFrom( self.flags )
            self.flags2.clear_obstacle()
            if not mark_obstacle( flags=self.flags2, obs=self.obs.part, center=obs_center2 ):
                emphasize( '  - obstacle position is invalid; stopping the obstacle (obs_stop=1)' )
                assert( self.method != FIXED_VOL or obs_naive )
                obs_stop = 1

        print( f'  - before. obs({it}): .state={self.obs.state}, .vel_vec={self.obs.vel_vec}, .force={self.obs.force}, dt={self.sol.timestep}, .center={self.obs.center}, .increase_vel={self.obs.increase_vel}, obs_center2={obs_center2}, n_stay={self.obs.n_stay}, start_h={self.obs.start_h}' )

        if self.obs.state != 7: # still moving
            if not obs_stop: # not stopped
                self.obs.n_stay = 0

                if self.obs.state < 2 and self.obs.center.y - self.obs.rad <= self.obs.start_h:
                    if self.obs.state == 0: # move to state 2
                        self.obs.state = 2
                        emphasize2( f'  - new obs.state: {self.obs.state}' )
                
                # force due to terminal velocity
                if self.obs.state == 2 and abs( self.obs.force.y ) > 0 and self.obs.vel_vec.y < self.obs.terminal_speed:
                    buoyancy = -5*self.gravity # -1
                    drag = -0.3*self.obs.vel_vec.y # -0.3
                    self.obs.force.y = buoyancy + drag + 5*self.gravity
                    print( f'  - update force due to terminal velocity: buoyancy={buoyancy}, drag={drag}, .force={self.obs.force}' )

                self.obs.center = obs_center2

            elif self.obs.state != 4: 
                # It was stopped and doesn't have a constant speed. It won't be natural if it continues in the same speed. Instead, set to terminal speed since it hit the water already, and the force is 0.
                # Need also to consider persistent trying (the obstacle keeps trying to make a move and doesn't stop completely) when the obstacle tries to clear the way (so particles won't come back).
                # The object just reaches terminal velocity instantly, which is also unnatural. I'm disabling this since it looks bad for flip and like I did it on purpose.
                #self.obs.vel_vec.y = self.obs.terminal_speed

                if int( it ) > self.obs.stay_last_it:
                    self.obs.stay_last_it = int(it)
                    self.obs.n_stay += 1
                if 0 and self.obs.n_stay > 50: # 50
                    print( f'  - no obstacle progress (n_stay={self.obs.n_stay}); stopping it' )
                    self.obs.vel_vec = Vec3( 0 )
                    self.obs.force = Vec3( 0 )
                    self.obs.state = 7
                    emphasize2( f'  - new obs.state: {self.obs.state}' )

        print( f'  - after.  obs({it}): .state={self.obs.state}, .vel_vec={self.obs.vel_vec}, .force={self.obs.force}, dt={self.sol.timestep}, .center={self.obs.center}, .increase_vel={self.obs.increase_vel}, obs_center2={obs_center2}, n_stay={self.obs.n_stay}, start_h={self.obs.start_h}' )

    def move_obstacle( self ):
        print( f'- move_obstacle(), .force={self.obs.force}, .vel_vec={self.obs.vel_vec}' )
        if self.obs.state != 7:
            if self.obs.center.y - self.obs.rad + self.obs.vel_vec.y*self.dt > 1.1: # move
                print( '  - obstacle still moves' )
                if self.obs.increase_vel:
                    dv = self.sol.timestep * self.obs.force
                    self.obs.vel_vec += dv
                    print( f'  - dv={dv}, .vel_vec={self.obs.vel_vec}' )

                    # I apply force due to terminal velocity. It may cause the obs to become too slow or even change direction. I limit the minimal speed after it hits the water.
                    if self.obs.state == 2 and self.obs.vel_vec.y > self.obs.terminal_speed:
                        self.obs.vel_vec.y = self.obs.terminal_speed
                        self.obs.force = Vec3( 0 )
                        print( '    - setting minimal speed to terminal speed:', self.obs.vel_vec )

                    # limit max speed 
                    max_y_speed = 35*self.gravity # 7, 10
                    if self.obs.vel_vec.y < max_y_speed:
                        print( f'    - limiting speed to {max_y_speed}' )
                        self.obs.vel_vec.y = max_y_speed

            else: # stay
                print( '  - obstacle reached the bottom' )
                self.obs.vel_vec = Vec3( 0 )
                self.obs.force = Vec3( 0 )
                self.obs.state = 7
                emphasize2( f'  - new obs.state: {self.obs.state}' )

        else:
            print( '- obstacle rests' )

        # obs.vel for fluid boundary conditions
        if 1:
            obs_vel_vec2 = self.obs.vel_vec + Vec3(0) # force copy

            # splash speed
            if self.splash and self.obs.state > 0:
                self.splash = 0 # instantaneous
                splash_scale = 1
                if self.dim == 3:
                    if self.method == FIXED_VOL:
                        splash_scale = 3
                        if self.obs.shape == 1: 
                            splash_scale *= 1.5
                else: # 2D
                    splash_scale = 1
                obs_vel_vec2.y *= splash_scale
                print( f'  - set obs.vel.y to {obs_vel_vec2.y} due to state 1, splash_scale={splash_scale}' )

            #obs_vel_vec2.x = 1 # debug
            self.obs.vel.setConst( obs_vel_vec2 )

            #self.obs.vel.setBoundMAC( value=Vec3( 0 ), boundaryWidth=self.boundary_width + 1, normalOnly=False )
            #self.obs.vel.setBound( value=Vec3( 0 ), boundaryWidth=self.boundary_width + 1 )
            self.obs.vel.set_bound_MAC2( value=Vec3( 0 ), boundaryWidth=self.boundary_width + 0 )
            
            #self.obs.vel.printGrid()

        else:
            if self.obs.state < 2:
                self.obs.state = 2
                emphasize2( f'  - new obs.state: {self.obs.state}' )
            self.obs.vel.setConst( Vec3(0.) )
            self.obs.vel.setBound( value=Vec3(0.), boundaryWidth=self.boundary_width+1 )

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
            self.flags.clear_obstacle()
            ret2 = mark_obstacle( flags=self.flags, obs=self.obs.part, center=self.obs.center )
            if self.method == FIXED_VOL:
                assert( ret2 )
        elif 0:
            setObstacleFlags( flags=self.flags, phiObs=self.phiObs )
        elif 0: # more precise
            updateFractions( flags=self.flags, phiObs=self.phiObs, fractions=fractions, boundaryWidth=self.boundary_width )
            setObstacleFlags( flags=self.flags, phiObs=self.phiObs, fractions=fractions )
        #self.flags.printGrid()

    def main( self ):
        if self.method != FIXED_VOL:
            if 1 or self.method != FLIP:
                self.narrowBand = bool( 0 )
                self.narrowBandWidth = -1

        combineBandWidth = self.narrowBandWidth - 1

        if 0 or self.dim == 2:
            set_print_2D( True )

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
        #fractions = self.sol.create(MACGrid) # Sticks to walls and becomes slow without this? Not anymore. Shrinks an obstacle box and draws it inaccurately.

        bfs = self.sol.create(IntGrid) # discrete distance from surface

        # acceleration data for particle
        pindex = self.sol.create(ParticleIndexSystem)
        gpi = self.sol.create(IntGrid)

        correct19 = Correct19( self.dim, self.sol, self.part_per_cell_1d, self.pp )

        # info
        print()
        print( 'dim:', self.dim, ', res:', self.res, ', ppc:', self.ppc )
        print( 'narrowBand:', self.narrowBand, ', narrowBandWidth:', self.narrowBandWidth, ', inter_control_method:', self.inter_control_method )
        print( 'method:', self.method )
        print( 'gravity: %0.02f' % self.gravity )
        print( 'timestep:', self.dt )

        # adaptive time stepping; from flip5
        b_adaptive_time_step = 0
        if b_adaptive_time_step:
            self.sol.frameLength = 1.0   # length of one frame (in "world time")
            self.sol.timestep    = 1.0
            self.sol.timestepMin = 0.5   # time step range
            self.sol.timestepMax = 1.0
            self.sol.cfl         = 1.0   # maximal velocity per cell, 0 to use fixed timesteps

        # setup scene
        self.setup_scene()

        # create dir
        name = '[%d,%d,%d]' % ( self.gs.x, self.gs.y, self.gs.z )

        if self.method == FLIP:
            name += ' flip'
        elif self.method == FIXED_VOL:
            if not self.narrowBand:
                name += ' full'
        elif self.method == CORRECT19:
            name += ' idp'
        elif self.method == DE_GOES22:
            name += ' power'
        elif self.method == MATLAB_FLIP:
            name += ' matlab flip'
        else:
            name += ' unknown'
        
        if self.narrowBand:
            name += f' band{self.narrowBandWidth}'

        if self.part_per_cell_1d != 2:
            name += f", {self.part_per_cell_1d}ppc"

        name += f", {self.scene['name']}"
        scene_name = name

        self.out_dir = out_dir_root + name + '/'
        if os.path.exists( self.out_dir ):
            error( f'path already exists: {self.out_dir}' )
            self.out_dir = None
            return
        os.mkdir( self.out_dir )

        # copy
        shutil.copy( out_dir_root + '../video.bat', self.out_dir )
        shutil.copy( out_dir_root + '../tweak_images.py', self.out_dir )

        # files
        f_set = open( self.out_dir + '_settings.txt', 'w' )
        c = self.gs
        f_set.write( '%d %d %d\n' % ( c.x, c.y, c.z ) ) # gs
        f_set.write( f"{self.scene['type']}\n" )
        f_set.write( f"{self.scene['cam']}\n" )
        f_set.flush()

        f_measure = open( self.out_dir + '_measure.txt', 'w' )

        if self.obs.exists:
            self.obs.file = open( self.out_dir + '_obstacle.txt', 'w' )
            self.obs.file.write( f'{self.obs.shape} {self.obs.rad/self.max_gs}\n' )
            self.obs.file.flush()

        # size of particles 
        radiusFactor = 1.0

        np = self.pp.pySize()
        V0 = float( np ) / self.ppc

        # for emitter
        #np_max = 2*np # twice as init
        np_max = self.ppc * (self.res-2)**self.dim * 0.5 # half tank

        print( f'# particles: {np} (np_max={np_max})' )

        # create a level set from particles
        # phi is influenced by the walls for some reason
        gridParticleIndex( parts=self.pp, flags=self.flags, indexSys=pindex, index=gpi )
        unionParticleLevelset( self.pp, pindex, self.flags, gpi, phiParts, radiusFactor )
        self.phi.copyFrom( phiParts )
        if self.narrowBand:
            extrapolateLsSimple( phi=self.phi, distance=self.narrowBandWidth+2, inside=True )
        else:
            extrapolateLsSimple( phi=self.phi, distance=4, inside=True ) # 4
        #self.phi.printGrid()

        # after initializing particles and before gui
        if self.b_fluid_mesh:
            mesh_gen = mesh_generator( self.dim, self.gs, self.sol, self.narrowBand, self.out_dir )
            # to resolve the surface extraction bug where the surface becomes separated on obs; upres=1 is likely the solution
            if self.obs.exists and self.large_obs and self.narrowBand:
                mesh_gen = mesh_generator( self.dim, self.gs, self.sol, self.narrowBand, self.out_dir, 1 )
                mesh_gen.union_method = 0 # might not be needed, just for good measure
                emphasize( f'setting mesh_gen.upres={mesh_gen.upres} (and mesh_gen.union_method={mesh_gen.union_method}) due to large obs' )
            mesh_gen.update_phi( self.phi )
            mesh_gen.generate( self.pp )

        if 1 and GUI:
            gui = Gui()
            gui.set_2D( self.b2D ) # ortho view

            # mesh
            mode = 2
            if 1 and self.b2D and not self.obs.exists:
                mode = 1
            for i in range( mode ):
                gui.nextMeshDisplay() # 0:full, 1:hide, 2:x-ray
            if self.obs2.exists:
                gui.setBackgroundMesh( self.obs2.mesh )
                gui.nextMesh()

            # field
            gui.setRealGridDisplay( 0 ) # 0:none, 1:volume
            gui.setVec3GridDisplay( 0 ) # 0:none, 1:vel
            for i in range( 0 ): # 0:center, 1:wall, 2:color, 3:none
                gui.nextVec3Display()

            # angle cam
            if 0 and self.dim == 3: # camera
                gui.setCamPos( 0, 0, -2.2 )
                gui.setCamRot( 35, -30, 0 )
            
            # hide grid
            if 1 and self.b2D:
                gui.toggleHideGrids()

            gui.show()
            #gui.pause()
        else:
            self.bScreenShot = 0

        it = 0
        it2 = 0

        # measure
        print( '- markFluidCells (update flags) for measure' )
        markFluidCells( parts=self.pp, flags=self.flags )
        self.flags.mark_surface()
        m = measure( self.pp, pVel, self.flags, self.ppc, V0, volume )
        f_measure.write( f'{m[0]}\n' )
        f_measure.flush()

        self.obs.write_pos()

        if self.bSaveMesh:
            mesh_gen.save( it )

        if self.bScreenShot:
            gui.screenshot( self.out_dir + 'frame_%04d.png' % it ); # slow

        if self.bSaveUni:
            pressure.save( self.out_dir + 'ref_parts_0000.uni' )
            self.pp.save( self.out_dir + 'parts_%04d.uni' % it )

        if self.bSaveVDB:
            self.phi.set_name( 'phi' )
            objects = [ self.flags, self.phi, self.pp ] # need the 3 of them for volumetric .vdb and they need to be named (bifrost looks for a channel by name)
            #objects = [ self.pp ]
            fname = self.out_dir + 'fluid_data_%04d.vdb' % it
            save( name=fname, objects=objects ) # error in debug mode "string too long"?

        # stat
        stat = {}
        stat['np0'] = np
        stat['np_avg'] = 0

        measu = []

        ##############################################################
        # loop
        ret = 0
        while 1:
            emphasize( '\n-----------------\n- time: %g(/%d; it2=%d)' % ( it, self.it_max, it2 ) )
            print( f'- scene: {scene_name}' )
            print( '- np=%d (np_max=%d)' % ( self.pp.pySize(), np_max ) )

            if 1 and ret != 0:
                error( f'Error: ret={ret}' )
                ret = 0
                break

            if not it < self.it_max:
                break

            tic( 'loop iteration' )

            # time step
            maxVel = self.vel.getMaxAbs()
            if not b_adaptive_time_step:
                self.sol.frameLength = self.dt
                frac = ( 1 - it % 1 )
                if 1 - frac < 1e-5:
                    frac = 1
                self.sol.timestep = frac * self.dt
            else: # adaptive for flip5
                self.sol.adaptTimestep( maxVel )

            # emit
            if 0 and it > 1e3 and self.pp.pySize() < np_max and ( not measu or measu[0] / self.res**self.dim < 0.9 ): # 1000
                xi = self.gs * Vec3( 0.5, 0.9, 0.5 )
                v = Vec3( 0, -3.0, 0 ) # -3
                n_emitted = 0
                for i in range( -1, 2 ):
                    for j in range( -1, 2 ):
                        if self.pp.pySize() >= np_max:
                            break
                        if self.dim == 2:
                            j = 0
                        emit_particles( self.pp, pVel, self.flags, self.part_per_cell_1d, xi + Vec3( i, 0, j ), v )
                        n_emitted += self.part_per_cell_1d
                        if self.dim == 2:
                            break
                V0 += float( n_emitted ) / self.ppc # update volume
                print( f'- emitted {n_emitted} partices, new V0={V0}, fluid_vol/res^{self.dim}={measu[0] / self.res**self.dim}' )
                
            # matlab fluid
            if self.method in [ DE_GOES22, MATLAB_FLIP ]:
                assert( self.b2D )
                it3 = it2
                # restart matlab
                if 1 and ( it2 + 1 ) % 500 == 0:
                    emphasize( 'restart matlab' )
                    matlab_eval( 'kill_matlab()', False, False )
                    restart_matlab()
                    init_matlab()
                    it3 = 0

                matlab_eval( rf"mlogn( '\n-----------------\n- time: {it}' );" )

                tic( 'matlab_fluid' )
                matlab_fluid( self.method, self.dt, self.res, self.part_per_cell_1d, self.gravity, it3, self.pp, pVel )
                toc()

            # the rest of the methods
            else:
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
                    velOld.copyFrom( self.vel ) # save before forces and pressure; like the end of mapPartsToMAC()
                    #print( '>> combine' )
                    #vel.printGrid()
                elif 1:
                    # map particle velocities to grid
                    mapPartsToMAC( vel=self.vel, flags=self.flags, velOld=velOld, parts=self.pp, partVel=pVel, weight=mapWeights )
                    extrapolateMACFromWeight( vel=self.vel , distance=2, weight=mapWeights )

                # update flags; there's also flags.updateFromLevelset()
                if self.method != FIXED_VOL or it == 0:
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
                        addGravity( flags=self.flags, vel=self.vel, gravity=(0, self.gravity, 0), scale=bool(bscale) )

                    # vortex
                    if 0:
                        #c = gs/2
                        c = Vec3( res/2, 0.3*res, res/2 )
                        vortex( pp=self.pp, dt=s.timestep, c=c, rad=0.1*res, h=0.9*res, pVel=pVel2 )
                        mapPartsToMAC( vel=vel2, flags=self.flags, velOld=vel2, parts=self.pp, partVel=pVel2 )
                        self.vel.add( vel2 )

                # set velocity for obstacles
                # there used to be another setWallBcs after the pressure solve, but it's not necessary: these are boundary conditions (must hold)
                print( '- set wall boundary conditions' )
                #self.obs.vel.printGrid()
                #self.vel.printGrid()

                #setWallBcs( flags=self.flags, vel=self.vel )
                #setWallBcs( flags=self.flags, vel=self.vel, fractions=fractions )
                #setWallBcs( flags=self.flags, vel=self.vel, fractions=fractions, phiObs=self.phiObs, obvel=self.obs.vel ) # calls KnSetWallBcsFrac, which doesn't work?
                #setWallBcs( flags=self.flags, vel=self.vel, obvel=self.obs.vel ) # calls KnSetWallBcs
                set_wall_bcs2( flags=self.flags, vel=self.vel, obvel=self.obs.vel )

                if self.obs2.exists:
                    self.obs2.set_wall_bcs( self.flags, self.vel )

                #self.vel.printGrid()
                #self.flags.printGrid()

                # pressure solve
                speed_limit = 7
                b_bad_pressure = 0
                if 1:
                    print( '- pressure' )
                    tic( 'pressure' )
                    # Solving Poisson eq for the pressure, with Neumann BC on walls (given as velocity, where u = \nabla p) and Dirichlet on empty cells (p=0).
                    # If there's fluid caged in a solid with no empty cells, then there are only Neumann conditions. Then, need to fix one pressure cell (Dirichlet) to eliminate DOF. Setting the flag zeroPressureFixing won't help if there are empty cells in other parts of the domain.
                    solvePressure( flags=self.flags, vel=self.vel, pressure=pressure, phi=self.phi )
                    t = toc()
                    if it2 == 0:
                        stat['pressure'] = t
                    else:
                        stat['pressure'] = ( stat['pressure'] * it2 + t ) / ( it2 + 1 ) # average

                    maxPVel = pVel.getMaxAbs()
                    maxVel = self.vel.getMaxAbs()
                    print( '  - vel.MaxAbs=%0.2f, pVel.MaxAbs=%0.2f, speed_limit=%g' % ( maxVel, maxPVel, speed_limit ) )

                    # limit vel
                    if 1 and maxVel > speed_limit:
                        print( f'  - maxVel is over the speed_limit({speed_limit})' )
                        if 1 and maxVel < 5000: # 40, 200, 5000; may want to use 0 for compressed scenes
                            if 0:
                                print(  f'    - scaling vel to speed_limit={speed_limit}' )
                                self.vel.multConst( Vec3(speed_limit/maxVel) )
                            else:
                                print(  f'    - clamping vel (norm) to speed_limit={speed_limit}' )
                                self.vel.clamp_norm( speed_limit )
                        else:
                            b_bad_pressure = 1
                            if 1:
                                warn( 'velocity blew up; clearing vel' )
                                self.vel.clear()
                            else:
                                warn( 'velocity blew up; reverting to old one' )
                                self.vel.copyFrom( velOld )

                dist = min( int( maxVel*1.25 + 2 ), 8 ) # res
                print( '- extrapolate MAC Simple (dist=%0.1f)' % dist )
                extrapolateMACSimple( flags=self.flags, vel=self.vel, distance=dist, intoObs=False )
                #self.flags.printGrid()
                #self.vel.printGrid()

                print( '- save particle positions in .pos0' )
                set_particles_pos0( pp=self.pp )

                # FLIP velocity update
                print( '- FLIP velocity update' )
                alpha = 0.1 # 0, .1
                flipVelocityUpdate( vel=self.vel, velOld=velOld, flags=self.flags, parts=self.pp, partVel=pVel, flipRatio=1 - alpha )
                #self.vel.printGrid()
                
                # limit pVel
                if 1:
                    limit_particle_velocity( pVel, speed_limit )

                # advect
                print( '- advect' )
                # advect particles
                self.pp.advectInGrid( flags=self.flags, vel=self.vel, integrationMode=IntEuler, deleteInObstacle=False, stopInObstacle=False ) # IntEuler, IntRK2, IntRK4
                if 0 and self.method != FIXED_VOL and self.method != CORRECT19:
                    print( '- pushOutofObs' )
                    pushOutofObs( parts=self.pp, flags=self.flags, phiObs=self.phiObs ) # creates issues for correct19 and fixedVol; can push particles into walls; pushes outside domain
                # advect phi; why? the particles should determine phi, which should flow on its own; disabling this creates artifacts in flip5; it makes it worse for fixed_vol
                if 1 and self.method != FIXED_VOL:
                    advectSemiLagrange( flags=self.flags, vel=self.vel, grid=self.phi, order=1 )
                    if 0:
                        self.flags.updateFromLevelset( self.phi ) # creates in 3D an extra layer of fluid without particles
                # advect grid velocity
                if self.narrowBand:
                    advectSemiLagrange( flags=self.flags, vel=self.vel, grid=self.vel, order=2 )

                # limit step (updates pp)
                self.sol.timestep = limit_time_step_to_one_cell_movement( self.pp, self.sol.timestep )

                # fixed_vol
                #self.flags.printGrid()
                include_walls = false
                obs_naive = 0 # specifically for fixed_vol
                obs_stop = 0
                print( f'- obs_naive={obs_naive}' )
                # fixed volume (my scheme)
                if self.method == FIXED_VOL:
                    [ ret, obs_stop, include_walls ] = self.fixed_volume( pVel, obs_naive, include_walls, ret, it2, bfs, stat )
                #self.flags.printGrid()

                # correct19
                if self.method == CORRECT19:
                    if b_bad_pressure:
                        warn( "skipping correct19 due to bad pressure" )
                    else:
                        correct19.main( self.sol, self.flags, self.pp, self.vel, pindex, gpi, self.phiObs )
                
                # moving obstacle
                if self.obs.exists:
                    self.update_obstacle( obs_naive, obs_stop, it )
                    self.move_obstacle()

                # static obstacle
                if self.obs2.exists and self.obs2.part:
                    mark_obstacle( flags=self.flags, obs=self.obs2.part, center=Vec3(0) )

                # for narrowBand, before updating phi with the particles
                if self.b_fluid_mesh:
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
                if self.method != FIXED_VOL:
                    pVel.setSource( self.vel, isMAC=True ) # set source grids for resampling, used in adjustNumber
                    minParticles = self.ppc
                    maxParticles = 2*minParticles # 2, 1(exacerbates artifact in flip5 dam 128?)
                    if self.narrowBand:
                        self.phi.setBoundNeumann( 0 ) # make sure no particles are placed at outer boundary
                        #self.phi.printGrid()
                        # vel is used only to get the parent
                        adjustNumber( parts=self.pp, vel=self.vel, flags=self.flags, minParticles=minParticles, maxParticles=maxParticles, phi=self.phi, narrowBand=self.narrowBandWidth, exclude=self.phiObs ) 
                    elif 0:
                        adjustNumber( parts=self.pp, vel=self.vel, flags=self.flags, minParticles=minParticles, maxParticles=maxParticles, phi=self.phi, exclude=self.phiObs ) 

            # update and mark surface for measure
            if self.method != FIXED_VOL:
                print( '- markFluidCells (update flags)' )
                if 1 and self.narrowBand:
                    self.flags.updateFromLevelset( self.phi )
                else:
                    markFluidCells( parts=self.pp, flags=self.flags )
                self.flags.mark_surface()

            toc() # iter

            # measure
            measu = measure( self.pp, pVel, self.flags, self.ppc, V0, volume )
            if len( measu ) > 2:
                stat['measure_min'] = measu[1]*100
                stat['measure_max'] = measu[2]*100

            # stat
            if 1:
                if not self.narrowBand or it2 > 0:
                    nIt = it2 if self.narrowBand else it2 + 1
                    stat['np_avg'] = ( ( nIt - 1) * stat['np_avg'] + self.pp.pySize() ) / nIt # average

                f_stat = open( self.out_dir + '_stat.csv', 'w' ) # rewrite

                stat2 = stat.copy()
                stat2['it'] = str( it + 1 )
                stat2['it2'] = str( it2 + 1 )

                # int to string
                for key in stat2.keys():
                    if not isinstance( stat2[key], str ):
                        if isinstance( stat2[key], float ):
                            stat2[key] = '%.1f' % stat2[key]
                        else:
                            stat2[key] = str( stat2[key] )

                #print( stat2 )
                f_stat.write( ','.join( stat2.keys() ) + '\n' )
                f_stat.write( ','.join( stat2.values() ) )
                f_stat.flush()

            # mesh
            if self.b_fluid_mesh:
                tic( 'mesh_gen' )
                mesh_gen.generate( self.pp )
                toc()

            # print
            if 0:
                self.flags.printGrid()
                #self.vel.printGrid()
                self.phiObs.printGrid()

                #self.pp.printParts()
                #self.pp.writeParticlesText( self.out_dir + 'particles_%04d.txt' % it )
            
            # step; updates gui and when pause takes place
            print( '- step (%g)' % it )
            self.sol.step( int(it) )
            #print( 'after step' )

            # it
            it += self.sol.timestep / self.dt
            it2 += 1

            # screenshot
            if self.bScreenShot:
                if abs( it - round(it) ) < 1e-7:
                    fr = round( it )
                else:
                    fr = math.ceil( it )
                fname = self.out_dir + 'frame_%04d.png' % fr
                print( f'- saving: {fname}' )
                gui.screenshot( fname ) # slow

            # integer it
            if 0 or abs( it - round(it) ) < 1e-7:
                it = round( it )

                # write measure
                f_measure.write( f'{ measu[0] }\n' )
                f_measure.flush()

                # write obs pos
                self.obs.write_pos()

                # data
                if self.bSaveUni:
                    # save particle data for flip03_gen.py surface generation scene
                    self.pp.save( self.out_dir + 'parts_%04d.uni' % it )

                if self.bSaveVDB:
                    # pdata fields must be before pp
                    objects = [ self.flags, self.phi, self.pp ]
                    #objects = [ self.pp ]
                    save( name=self.out_dir + 'fluid_data_%04d.vdb' % it, objects=objects )

                if self.bSaveMesh:
                    mesh_gen.save( it )
                
        if 1:
            # crop
            cwd = os.getcwd()
            os.chdir( self.out_dir )
            os.system( 'python tweak_images.py' )
            os.chdir( cwd )

            # video
            os.system( f'"{self.out_dir}/video.bat"' )

        # the following code is not reached if quitting manta (with esc)
        # pause
        if 0:
            print( '(zflip.py) press a key...' )
            keyboard.read_key()
        elif 0:
            print( '(zflip.py) press enter...' )
            input()

# __main__
if 1 and __name__ == '__main__':
    assert( len(sys.argv) >= 2 )
    method = int( sys.argv[1] )

    # auto-flush
    sys.stdout.reconfigure( line_buffering=True )

    out_dir_root = r'c:/prj-external-libs/mantaflow/out/'

    # (debug) for a consistent result; for large res, the step() hangs?
    if 0:
        limit_to_one_core()

    # so kernels won't hog the cpu
    if 1:
        limit_tbb_cores( 4 )

    #setDebugLevel( 10 )

    # init matlab
    if 1:
        init_matlab()

    # test
    #test_MAC()

    # simulation
    sim = simulation( method )
    sim.main()

    print( 'simulation ended' )

    # move matlab log; also in run.py
    path = Path( r'c:\prj\test_data\relative\_tmp\_matlab_out.txt' )
    if path.exists():
        print( f'moving matlab log, {path}' )
        shutil.move( path.as_posix(), sim.out_dir.as_posix() )
    else:
        print( f'no matlab log, {path}' )

    if sim.out_dir:
        os.system( f'copy_log.bat "{sim.out_dir}"' )

    print( 'end of zflip.py' )
    print( '' )
