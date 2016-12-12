"""
Defines L{runSimulation} function which executes compression simulation.
"""
##06.11.2013: doubled checkpointer resolution from 20 to 40 cpts

from esys.lsm import *
from esys.lsm.util import Vec3, BoundingBox, InstallInfo
from time import *
import sys

class Loading(Runnable):
    """
    Objects of this class provide the loading mechanism for a compression
    simulation.
    """
    def __init__(self,lsm,dt,t0,t1,vel,sig_v,sig_h):
        Runnable.__init__(self)
        self.theLSM=lsm
        self.t0=t0 #horizontal and vertical stress at maximum, horizontal stress begins to decline
        self.t1=t1 #time at which horizontal stress has declined to sigma_h
        self.dx=vel*dt
        self.ddx=(vel/(float(t1-t0)))*dt;
        print "dvel", self.ddx
        self.dsig_v=sig_v/float(t0) #stress ramp up increments
        self.dsig_h=(sig_v-sig_h)/float(t1-t0) #stress ramp down increments relative to sigma_n
        self.sig_v=sig_v
        self.sig_h=sig_h
        
        print "dsig_v",self.dsig_v

        self.init_time=time()
        
    def run(self):

        t=self.theLSM.getTimeStep()

        # move sidewalls (deprecated)
        '''
        if(t>self.t0):
            if(t<self.t1):
                dxi=(t-self.t0)*self.ddx
                self.theLSM.moveWallBy("rightWall",Vec3(dxi,0.0,0.0))
                self.theLSM.moveWallBy("leftWall",Vec3(-1.0*dxi,0.0,0.0))
            else:
                self.theLSM.moveWallBy("rightWall",Vec3(self.dx,0.0,0.0))
                self.theLSM.moveWallBy("leftWall",Vec3(-1.0*self.dx,0.0,0.0))
        '''

        #--- apply vertical confining stress ---
        # get wall positions to calculate area for force
        le_wp=self.theLSM.getWallPosition("leftWall")
        ri_wp=self.theLSM.getWallPosition("rightWall")
        fr_wp=self.theLSM.getWallPosition("frontWall")
        ba_wp=self.theLSM.getWallPosition("backWall")
        up_wp=self.theLSM.getWallPosition("upperWall")
        lo_wp=self.theLSM.getWallPosition("lowerWall")
        # calculate length and width
        length=ri_wp[0]-le_wp[0]#for the vertical normal force
        width=ba_wp[2]-fr_wp[2]
        height=up_wp[1]-lo_wp[1]#for the horizontal normal force

        # area
        area_top=length*width#for the vertical normal force
        area_side=width*height#for the horizontal normal force
        # calculate applied force
        if(t<self.t0):
            f_v=t*self.dsig_v*area_top
            f_h=t*self.dsig_v*area_side
        else:
            f_v=self.sig_v*area_top
            if (t<self.t1):
                f_h=(self.sig_v-(float(t-self.t0)*self.dsig_h))*area_side
            else:
                f_h=self.sig_h*area_side

            
        if(f_v>0.0):
            self.theLSM.applyForceToWall("lowerWallInteraction",Vec3(0.0,f_v,0.0));
            self.theLSM.applyForceToWall("upperWallInteraction",Vec3(0.0,-1.0*f_v,0.0));
        if(f_h>0.0):
            self.theLSM.applyForceToWall("leftWallInteraction",Vec3(f_h,0.0,0.0));
            self.theLSM.applyForceToWall("rightWallInteraction",Vec3(-1.0*f_h,0.0,0.0));        

        if((t%10)==0):
            print t, "time steps", " at time: " , time()-self.init_time, "force", f_v, "length", length, "width",width,"area",area_top
            
def runSimulation():
    #    setVerbosity(True)
    # agruments: x-dim, y-dim, z-dim, rmax, nr. of timesteps, cohesion
    mySim=LsmMpi(8,[2,2,2])
    SimID=str(sys.argv[1])
    infile=str(sys.argv[2])
    nt=int(50000)
    xsize=float(120.0)
    ysize=float(50.0)
    zsize=float(30.0)
    
    rmax=1

    t0=int(nt/20) #was:/20
    t1=int(nt/10)#was:/10
    dt=0.01
    dt_save=int(nt/10000)
    dt_snap=int(nt/50)
    vel=(0.075*xsize)/(float(nt)*dt) #does not mean anything with stress controlled walls


   # print "velocity : ", vel
    # normal stress
    sigma_v=float(sys.argv[3])#0.001 ~20MPa @ E=20GPa
    sigma_h=float(sys.argv[4])
    mySim.initVerletModel("RotSphere", rmax*2.2, 0.1)
    mySim.setSpatialDomain(Vec3(-0.3*xsize,0,-0.1*zsize),Vec3(1.3*xsize,ysize,1.1*zsize))
    mySim.setTimeStepSize(dt)
    mySim.readGeometry(infile);
   
    # setup interactions
    # friction for "solid" particles
    mu=0.6
    fip=RotFrictionPrms( "friction", 1.0, mu, mu,1.0,True)
    # dashpot
    dip=LinearDashpotPrms("dashpot",0.1,1.2)
    # elastic for "ductile" particles
    eip1=NRotElasticPrms("elastic1", 1.0,True)
    eip2=NRotElasticPrms("elastic2", 1.0,True)
   
    # bonded

    # bulk bonds
    E_bulk=1.0 # youngs modulus
    gamma_bulk=0.3 # poisson ratio
    C_bulk=float(0.01) # cohesion
    atan_phi_bulk=0.6 # arctan of friction angle
    bip_bulk=BrittleBeamPrms("bonded_bulk",E_bulk,gamma_bulk,C_bulk,atan_phi_bulk,0)
    
    # interface bonds bulk - new glue
    E_newglue_bulk=1.0 # youngs modulus
    gamma_newglue_bulk=0.3 # poisson ratio
    C_newglue_bulk=float(0.01) # cohesion
    atan_phi_newglue_bulk=0.6 # arctan of friction angle
    bip_newglue_bulk=BrittleBeamPrms("bonded_newglue_bulk",E_newglue_bulk,gamma_newglue_bulk,C_newglue_bulk,atan_phi_newglue_bulk,1)

    # interface old glue - new glue
    E_newglue_oldglue=1.0 # youngs modulus
    gamma_newglue_oldglue=0.3 # poisson ratio
    C_newglue_oldglue=float(0.01) # cohesion
    atan_phi_newglue_oldglue=0.6 # arctan of friction angle
    bip_newglue_oldglue=BrittleBeamPrms("bonded_newglue_oldglue",E_newglue_oldglue,gamma_newglue_oldglue,C_newglue_oldglue,atan_phi_newglue_oldglue,2)
 
   # bonds new glue
    E_newglue=1.0 # youngs modulus
    gamma_newglue=0.3 # poisson ratio
    C_newglue=float(0.01)# cohesion
    atan_phi_newglue=0.6 # arctan of friction angle
    bip_newglue=BrittleBeamPrms("bonded_newglue",E_newglue,gamma_newglue,C_newglue,atan_phi_newglue,3)

    #bonds oldglue
    E_oldglue=1.0 # youngs modulus
    gamma_oldglue=0.3 # poisson ratio
    C_oldglue=float(0.01)# cohesion
    atan_phi_oldglue=0.6 # arctan of friction angle
    bip_oldglue=BrittleBeamPrms("bonded_oldglue",E_oldglue,gamma_oldglue,C_oldglue,atan_phi_oldglue,4)

    #interface oldglue - bulk
    E_oldglue_bulk=1.0 # youngs modulus
    gamma_oldglue_bulk=0.3 # poisson ratio
    C_oldglue_bulk=float(0.01)# cohesion
    atan_phi_oldglue_bulk=0.6 # arctan of friction angle
    bip_oldglue_bulk=BrittleBeamPrms("bonded_oldglue_bulk",E_oldglue_bulk,gamma_oldglue_bulk,C_oldglue_bulk,atan_phi_oldglue_bulk,5)

    # create interaction groups
    # elastic interactions ductile - ductile + ductile - brittle
    # assumes _only_ ductile particles have bit 7 (=64) set
    mySim.createInteractionGroupTagged(eip1,64,64,0,0)
    # dashpot interactions ductile - ductile + ductile - brittle
    mySim.createInteractionGroupTagged(dip,64,64,0,0)
    mySim.createInteractionGroup(bip_bulk)
    mySim.createInteractionGroup(bip_newglue_bulk)
    mySim.createInteractionGroup(bip_newglue_oldglue)
    mySim.createInteractionGroup(bip_newglue)
    mySim.createInteractionGroup(bip_oldglue)
    mySim.createInteractionGroup(bip_oldglue_bulk)
    # friction between "solid" particles
    # - assumes all "solid" particles have bit 6 (=32) set to 1
    #   i.e. bulk=32, glue1=33, .... 
    mySim.createInteractionGroupTagged(fip,32,32,32,32)
    mySim.createExclusion("bonded_bulk","friction")    
    mySim.createExclusion("bonded_newglue_bulk","friction")    
    mySim.createExclusion("bonded_newglue_oldglue","friction")    
    mySim.createExclusion("bonded_newglue","friction")    
    mySim.createExclusion("bonded_oldglue","friction")    
    mySim.createExclusion("bonded_oldglue_bulk","friction")    
    # create walls
    # top & bottom
    mySim.createWall("lowerWall",Vec3(0.0,0.0,0.0),Vec3(0.0,1.0,0.0))
    mySim.createWall("upperWall",Vec3(0.0,ysize,0.0),Vec3(0.0,-1.0,0.0))
    wp1=NRotElasticWallPrms("upperWallInteraction","upperWall",1.0)
    wp2=NRotElasticWallPrms("lowerWallInteraction","lowerWall",1.0);
    mySim.createInteractionGroup(wp1)
    mySim.createInteractionGroup(wp2)
    # left & right
    mySim.createWall("leftWall",Vec3(0.0,0.0,0.0),Vec3(1.0,0.0,0.0))
    mySim.createWall("rightWall",Vec3(xsize,ysize,zsize),Vec3(-1.0,0.0,0.0))
    wp3=NRotElasticWallPrms("leftWallInteraction","leftWall",1.0)
    wp4=NRotElasticWallPrms("rightWallInteraction","rightWall",1.0);
    mySim.createInteractionGroup(wp3)
    mySim.createInteractionGroup(wp4)
    # front & back
    mySim.createWall("frontWall",Vec3(0.0,0.0,0.0),Vec3(0.0,0.0,1.0))
    mySim.createWall("backWall",Vec3(xsize,ysize,zsize),Vec3(0.0,0.0,-1.0))
    wp5=NRotElasticWallPrms("frontWallInteraction","frontWall",1.0)
    wp6=NRotElasticWallPrms("backWallInteraction","backWall",1.0);
    mySim.createInteractionGroup(wp5)
    mySim.createInteractionGroup(wp6)
    # setup savers
    pot_file_name="saver/epot_f_"+str(SimID)+".dat"
    ep_prm=InteractionScalarFieldSaverPrms("friction","potential_energy",pot_file_name,"SUM",0,nt,dt_save)
    mySim.createFieldSaver(ep_prm)
    pot_file_name="saver/epot_b_"+str(SimID)+".dat"
    ep_prm=InteractionScalarFieldSaverPrms("bonded_bulk","potential_energy",pot_file_name,"SUM",0,nt,dt_save)
    mySim.createFieldSaver(ep_prm)
    kin_file_name="saver/ekin_"+str(SimID)+".dat"
    ek_prm=ParticleScalarFieldSaverPrms("e_kin",kin_file_name,"SUM",0,nt,dt_save)
    mySim.createFieldSaver(ek_prm)
    wf_file_name="saver/wf_"+str(SimID)+".dat"
    mySim.createFieldSaver(
      WallVectorFieldSaverPrms(
        fileName=wf_file_name,
        fieldName="Force",
        wallName=["lowerWall","upperWall","leftWall","rightWall","frontWall","backWall"],
        fileFormat="RAW_SERIES",
        beginTimeStep=0,
        endTimeStep=nt,
        timeStepIncr=dt_save
      )
    )
    wp_file_name="saver/wp_"+str(SimID)+".dat"
    mySim.createFieldSaver(
      WallVectorFieldSaverPrms(
        fileName=wp_file_name,
        fieldName="Position",
        wallName=["lowerWall","upperWall","leftWall","rightWall","frontWall","backWall"],
        fileFormat="RAW_SERIES",
        beginTimeStep=0,
        endTimeStep=nt,
        timeStepIncr=dt_save
      )
    )
    nb_file_name1="saver/nbonds_bulk"+str(SimID)+".dat"
    nb_file_name2="saver/nbonds_newglue_bulk"+str(SimID)+".dat"
    nb_file_name3="saver/nbonds_oldglue_newglue"+str(SimID)+".dat"
    nb_file_name4="saver/nbonds_newglue"+str(SimID)+".dat"
    nb_file_name5="saver/nbonds_oldglue"+str(SimID)+".dat"
    nb_file_name6="saver/nbonds_oldglue_bulk"+str(SimID)+".dat"
    nb_prm_bulk=InteractionScalarFieldSaverPrms("bonded_bulk","count",nb_file_name1,"SUM",0,nt,dt_save)
    mySim.createFieldSaver(nb_prm_bulk)
    nb_prm_newglue_bulk=InteractionScalarFieldSaverPrms("bonded_newglue_bulk","count",nb_file_name2,"SUM",0,nt,dt_save)
    mySim.createFieldSaver(nb_prm_newglue_bulk)
    nb_prm_newglue_oldglue=InteractionScalarFieldSaverPrms("bonded_newglue_oldglue","count",nb_file_name3,"SUM",0,nt,dt_save)
    mySim.createFieldSaver(nb_prm_newglue_oldglue) 
    nb_prm_newglue=InteractionScalarFieldSaverPrms("bonded_newglue","count",nb_file_name4,"SUM",0,nt,dt_save)
    mySim.createFieldSaver(nb_prm_newglue) 
    nb_prm_oldglue=InteractionScalarFieldSaverPrms("bonded_oldglue","count",nb_file_name5,"SUM",0,nt,dt_save)
    mySim.createFieldSaver(nb_prm_oldglue) 
    nb_prm_oldglue_bulk=InteractionScalarFieldSaverPrms("bonded_oldglue_bulk","count",nb_file_name6,"SUM",0,nt,dt_save)
    mySim.createFieldSaver(nb_prm_oldglue_bulk) 

   # checkpoints
    cp=CheckPointPrms("cpt/cpt",0,nt,dt_snap)
    mySim.createCheckPointer(cp)
    # add loading function
    lf=Loading(mySim,dt,t0,t1,vel,sigma_v,sigma_h)
    mySim.addPreTimeStepRunnable(lf)
    mySim.setNumTimeSteps(nt)
    start_time=time()
    print "######simulaion started:" 
    print strftime("###### %a, %d %b %Y %H:%M:%S gmt", gmtime())
    mySim.run()
    stop_time=time()
    print "runtime: ", stop_time-start_time, " seconds"
    duration=stop_time-start_time
    m, s = divmod(duration, 60)
    h, m = divmod(m, 60)
    print "runtime: %d:%02d:%02d" % (h, m, s)

if (__name__ == "__main__"):
    runSimulation()



