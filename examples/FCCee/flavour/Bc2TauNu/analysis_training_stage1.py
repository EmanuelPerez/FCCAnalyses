import sys
import ROOT
from array import array

print ("Load cxx analyzers ... ",)
ROOT.gSystem.Load("libedm4hep")
ROOT.gSystem.Load("libpodio")
ROOT.gSystem.Load("libawkward")
ROOT.gSystem.Load("libawkward-cpu-kernels")
ROOT.gSystem.Load("libFCCAnalyses")

ROOT.gErrorIgnoreLevel = ROOT.kFatal
_edm  = ROOT.edm4hep.ReconstructedParticleData()
_pod  = ROOT.podio.ObjectID()
_fcc  = ROOT.dummyLoader

print ('edm4hep  ',_edm)
print ('podio    ',_pod)
print ('fccana   ',_fcc)


class analysis():

    #__________________________________________________________
    def __init__(self, inputlist, outname, ncpu):
        self.outname = outname
        if ".root" not in outname:
            self.outname+=".root"

        ROOT.ROOT.EnableImplicitMT(ncpu)
        ROOT.EnableThreadSafety()
        self.df = ROOT.RDataFrame("events", inputlist)
        print (" init done, about to run")
    #__________________________________________________________
    def run(self):
        #df2 = (self.df.Range(1000)        
        df2 = (self.df
               #############################################
               ##          Aliases for # in python        ##
               #############################################
               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               .Alias("Particle0", "Particle#0.index")
               .Alias("Particle1", "Particle#1.index")
              
               #############################################
               ##               Build MC Vertex           ##
               #############################################
               .Define("MCVertexObject", "myUtils::get_MCVertexObject(Particle, Particle0)")
               .Define("MC_Vertex_x",    "myUtils::get_MCVertex_x(MCVertexObject)")
               .Define("MC_Vertex_y",    "myUtils::get_MCVertex_y(MCVertexObject)")
               .Define("MC_Vertex_z",    "myUtils::get_MCVertex_z(MCVertexObject)")
               .Define("MC_Vertex_ind",  "myUtils::get_MCindMCVertex(MCVertexObject)")
               .Define("MC_Vertex_ntrk", "myUtils::get_NTracksMCVertex(MCVertexObject)")
               .Define("MC_Vertex_n",    "int(MC_Vertex_x.size())")
               .Define("MC_Vertex_PDG",  "myUtils::get_MCpdgMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGmother",  "myUtils::get_MCpdgMotherMCVertex(MCVertexObject, Particle)")
               .Define("MC_Vertex_PDGgmother", "myUtils::get_MCpdgGMotherMCVertex(MCVertexObject, Particle)")


               #############################################
               ##              Build Reco Vertex          ##
               #############################################
               .Define("VertexObject", "myUtils::get_VertexObject(MCVertexObject,ReconstructedParticles,EFlowTrack_1,MCRecoAssociations0,MCRecoAssociations1)")

               #############################################
               ##          Build PV var and filter        ##
               #############################################
               .Define("EVT_hasPV",    "myUtils::hasPV(VertexObject)")
               .Define("EVT_NtracksPV", "float(myUtils::get_PV_ntracks(VertexObject))")
               .Define("EVT_NVertex",   "float(VertexObject.size())")
               .Filter("EVT_hasPV==1")


               #############################################
               ##          Build RECO P with PID          ##
               #############################################
               .Define("RecoPartPID" ,"myUtils::PID(ReconstructedParticles, MCRecoAssociations0,MCRecoAssociations1,Particle)")
               
               #############################################
               ##    Build RECO P with PID at vertex      ##
               #############################################
               .Define("RecoPartPIDAtVertex" ,"myUtils::get_RP_atVertex(RecoPartPID, VertexObject)")

               #############################################
               ##         Build vertex variables          ##
               #############################################
               .Define("Vertex_x",        "myUtils::get_Vertex_x(VertexObject)")
               .Define("Vertex_y",        "myUtils::get_Vertex_y(VertexObject)")
               .Define("Vertex_z",        "myUtils::get_Vertex_z(VertexObject)")
               .Define("Vertex_xErr",     "myUtils::get_Vertex_xErr(VertexObject)")
               .Define("Vertex_yErr",     "myUtils::get_Vertex_yErr(VertexObject)")
               .Define("Vertex_zErr",     "myUtils::get_Vertex_zErr(VertexObject)")

               .Define("Vertex_chi2",     "myUtils::get_Vertex_chi2(VertexObject)")
               .Define("Vertex_mcind",    "myUtils::get_Vertex_indMC(VertexObject)")
               .Define("Vertex_ind",      "myUtils::get_Vertex_ind(VertexObject)")
               .Define("Vertex_isPV",     "myUtils::get_Vertex_isPV(VertexObject)")
               .Define("Vertex_ntrk",     "myUtils::get_Vertex_ntracks(VertexObject)")
               .Define("Vertex_n",        "int(Vertex_x.size())")
               .Define("Vertex_mass",     "myUtils::get_Vertex_mass(VertexObject,RecoPartPIDAtVertex)")

               .Define("Vertex_d2PV",     "myUtils::get_Vertex_d2PV(VertexObject,-1)")
               .Define("Vertex_d2PVx",    "myUtils::get_Vertex_d2PV(VertexObject,0)")
               .Define("Vertex_d2PVy",    "myUtils::get_Vertex_d2PV(VertexObject,1)")
               .Define("Vertex_d2PVz",    "myUtils::get_Vertex_d2PV(VertexObject,2)")
               
               .Define("Vertex_d2PVErr",  "myUtils::get_Vertex_d2PVError(VertexObject,-1)")
               .Define("Vertex_d2PVxErr", "myUtils::get_Vertex_d2PVError(VertexObject,0)")
               .Define("Vertex_d2PVyErr", "myUtils::get_Vertex_d2PVError(VertexObject,1)")
               .Define("Vertex_d2PVzErr", "myUtils::get_Vertex_d2PVError(VertexObject,2)")
               
               .Define("Vertex_d2PVSig",  "Vertex_d2PV/Vertex_d2PVErr")
               .Define("Vertex_d2PVxSig", "Vertex_d2PVx/Vertex_d2PVxErr")
               .Define("Vertex_d2PVySig", "Vertex_d2PVy/Vertex_d2PVyErr")
               .Define("Vertex_d2PVzSig", "Vertex_d2PVz/Vertex_d2PVzErr")

               .Define("Vertex_d2MC",     "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,-1)")
               .Define("Vertex_d2MCx",    "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,0)")
               .Define("Vertex_d2MCy",    "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,1)")
               .Define("Vertex_d2MCz",    "myUtils::get_Vertex_d2MC(VertexObject,MCVertexObject,Vertex_mcind,2)")

               .Define("EVT_dPV2DVmin",   "myUtils::get_dPV2DV_min(Vertex_d2PV)")
               .Define("EVT_dPV2DVmax",   "myUtils::get_dPV2DV_max(Vertex_d2PV)")
               .Define("EVT_dPV2DVave",   "myUtils::get_dPV2DV_ave(Vertex_d2PV)")
               
               #############################################
               ##        Build Tau -> 3Pi candidates      ##
               #############################################
               .Define("Tau23PiCandidates",         "myUtils::build_tau23pi(VertexObject,RecoPartPIDAtVertex)")
               
               #############################################
               ##       Filter Tau -> 3Pi candidates      ##
               ############################################# 
               .Define("EVT_NTau23Pi",              "float(myUtils::getFCCAnalysesComposite_N(Tau23PiCandidates))")
               .Filter("EVT_NTau23Pi>0")

               
               #############################################
               ##              Build the thrust           ##
               ############################################# 
               .Define("RP_e",          "ReconstructedParticle::get_e(RecoPartPIDAtVertex)")
               .Define("RP_px",         "ReconstructedParticle::get_px(RecoPartPIDAtVertex)")
               .Define("RP_py",         "ReconstructedParticle::get_py(RecoPartPIDAtVertex)")
               .Define("RP_pz",         "ReconstructedParticle::get_pz(RecoPartPIDAtVertex)")
               .Define("RP_charge",     "ReconstructedParticle::get_charge(RecoPartPIDAtVertex)")
              
               .Define("EVT_thrustNP",      'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
               .Define("RP_thrustangleNP",  'Algorithms::getAxisCosTheta(EVT_thrustNP, RP_px, RP_py, RP_pz)')
               .Define("EVT_thrust",        'Algorithms::getThrustPointing(RP_thrustangleNP, RP_e, EVT_thrustNP, 1.)')
               .Define("RP_thrustangle",    'Algorithms::getAxisCosTheta(EVT_thrust, RP_px, RP_py, RP_pz)')

               
               #############################################
               ##        Get thrust related values        ##
               ############################################# 
               ##hemis0 == negative angle == max energy hemisphere if pointing
               ##hemis1 == positive angle == min energy hemisphere if pointing
               .Define("EVT_thrusthemis0_n",    "Algorithms::getAxisN(0)(RP_thrustangle, RP_charge)")
               .Define("EVT_thrusthemis1_n",    "Algorithms::getAxisN(1)(RP_thrustangle, RP_charge)")
               .Define("EVT_thrusthemis0_e",    "Algorithms::getAxisEnergy(0)(RP_thrustangle, RP_charge, RP_e)")
               .Define("EVT_thrusthemis1_e",    "Algorithms::getAxisEnergy(1)(RP_thrustangle, RP_charge, RP_e)")

               .Define("EVT_ThrustEmax_E",         "EVT_thrusthemis0_e.at(0)")
               .Define("EVT_ThrustEmax_Echarged",  "EVT_thrusthemis0_e.at(1)")
               .Define("EVT_ThrustEmax_Eneutral",  "EVT_thrusthemis0_e.at(2)")
               .Define("EVT_ThrustEmax_Ncharged",  "float(EVT_thrusthemis0_n.at(1))")
               .Define("EVT_ThrustEmax_Nneutral",  "float(EVT_thrusthemis0_n.at(2))")

               .Define("EVT_ThrustEmin_E",         "EVT_thrusthemis1_e.at(0)")
               .Define("EVT_ThrustEmin_Echarged",  "EVT_thrusthemis1_e.at(1)")
               .Define("EVT_ThrustEmin_Eneutral",  "EVT_thrusthemis1_e.at(2)")
               .Define("EVT_ThrustEmin_Ncharged",  "float(EVT_thrusthemis1_n.at(1))")
               .Define("EVT_ThrustEmin_Nneutral",  "float(EVT_thrusthemis1_n.at(2))")


               .Define("Vertex_thrust_angle",   "myUtils::get_Vertex_thrusthemis_angle(VertexObject, RecoPartPIDAtVertex, EVT_thrust)")
               ###0 == negative angle==max energy , 1 == positive angle == min energy
               .Define("Vertex_thrusthemis_emin",    "myUtils::get_Vertex_thrusthemis(Vertex_thrust_angle, 1)")
               .Define("Vertex_thrusthemis_emax",    "myUtils::get_Vertex_thrusthemis(Vertex_thrust_angle, 0)")

               .Define("EVT_ThrustEmin_NDV", "float(myUtils::get_Npos(Vertex_thrust_angle))")
               .Define("EVT_ThrustEmax_NDV", "float(myUtils::get_Nneg(Vertex_thrust_angle))")

               .Define("TrueTau23PiBc_vertex",        "myUtils::get_trueVertex(MCVertexObject,Particle,Particle0, 15, 541)")
               .Define("TrueTau23PiBu_vertex",        "myUtils::get_trueVertex(MCVertexObject,Particle,Particle0, 15, 521)")


               .Define("Tau23PiCandidates_anglethrust", "myUtils::getFCCAnalysesComposite_anglethrust(Tau23PiCandidates, EVT_thrust)")
               .Define("CUT_hasCandEmin",           "myUtils::has_anglethrust_emin(Tau23PiCandidates_anglethrust)")
               .Filter("CUT_hasCandEmin>0")
 

           )
        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [
                

                "EVT_ThrustEmin_E",          "EVT_ThrustEmax_E",
                "EVT_ThrustEmin_Echarged",   "EVT_ThrustEmax_Echarged",
                "EVT_ThrustEmin_Eneutral",   "EVT_ThrustEmax_Eneutral",
                "EVT_ThrustEmin_Ncharged",   "EVT_ThrustEmax_Ncharged",
                "EVT_ThrustEmin_Nneutral",   "EVT_ThrustEmax_Nneutral",
                "EVT_ThrustEmin_NDV",        "EVT_ThrustEmax_NDV",
                "EVT_dPV2DVmin","EVT_dPV2DVmax","EVT_dPV2DVave",
                "EVT_NtracksPV", "EVT_NVertex", "EVT_NTau23Pi",
 
                "TrueTau23PiBc_vertex",
                "TrueTau23PiBu_vertex",                
                
                ]:
            branchList.push_back(branchName)

        #opts = ROOT.RDF.RSnapshotOptions()
        #opts.fCompressionAlgorithm = ROOT.ROOT.kLZ4
        #opts.fCompressionLevel = 3
        #opts.fAutoFlush = -1024*1024*branchList.size()
        #df2.Snapshot("events", self.outname, branchList, opts)
        df2.Snapshot("events", self.outname, branchList)

# example call for standalone file
# python examples/FCCee/flavour/Bc2TauNu/analysis_stage1.py flat_ee_Zbb_Bc2TauNu.root /eos/experiment/fcc/ee/generation/DelphesEvents/fcc_tmp_v03/p8_ee_Zbb_ecm91_EvtGen_Bc2TauNuTAUHADNU/events_003834121.root

# python examples/FCCee/flavour/Bc2TauNu/analysis_stage1.py flat_ee_Zbb_Bu2TauNu.root /eos/experiment/fcc/ee/generation/DelphesEvents/fcc_tmp_v03/p8_ee_Zbb_ecm91_EvtGen_Bu2TauNuTAUHADNU/events_026079857.root

# python examples/FCCee/flavour/Bc2TauNu/analysis_stage1.py  flat_ee_Zbb_Bc2TauNu.root "/eos/experiment/fcc/ee/generation/DelphesEvents/fcc_tmp_v03/p8_ee_Zbb_ecm91_EvtGen_Bc2TauNuTAUHADNU/events_*"

# python examples/FCCee/flavour/Bc2TauNu/analysis_stage1.py flat_ee_Zbb.root  /eos/experiment/fcc/ee/generation/DelphesEvents/fcc_tmp_v03/p8_ee_Zbb_ecm91/events_026734131.root



if __name__ == "__main__":

    if len(sys.argv)<3:
        print ("usage:")
        print ("python ",sys.argv[0]," output.root input.root")
        print ("python ",sys.argv[0]," output.root \"inputdir/*.root\"")
        print ("python ",sys.argv[0]," output.root file1.root file2.root file3.root <nevents>")
        sys.exit(3)


    print ("Create dataframe object from ", )
    fileListRoot = ROOT.vector('string')()
    nevents=0

    if len(sys.argv)==3 and "*" in sys.argv[2]:
        import glob
        filelist = glob.glob(sys.argv[2])
        for fileName in filelist:
            fileListRoot.push_back(fileName)
            print (fileName, " ",)
            print (" ...")


    elif len(sys.argv)>2:
        for i in range(2,len(sys.argv)):
            try:
                nevents=int(sys.argv[i])
                print ("nevents found (will be in the processed events branch in root tree):",nevents)
            except ValueError:
                fileListRoot.push_back(sys.argv[i])
                print (sys.argv[i], " ",)
                print (" ...")

                         
    outfile=sys.argv[1]
    print("output file:  ",outfile)
    if len(outfile.split("/"))>1:
        import os
        os.system("mkdir -p {}".format(outfile.replace(outfile.split("/")[-1],"")))

    if nevents==0:
        for f in fileListRoot:
            tf=ROOT.TFile.Open(str(f),"READ")
            tt=tf.Get("events")
            nevents+=tt.GetEntries()
    print ("nevents ", nevents)
    
    import time
    start_time = time.time()
    ncpus = 8
    analysis = analysis(fileListRoot, outfile, ncpus)
    analysis.run()

    elapsed_time = time.time() - start_time
    print  ("==============================SUMMARY==============================")
    print  ("Elapsed time (H:M:S)     :  ",time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    print  ("Events Processed/Second  :  ",int(nevents/elapsed_time))
    print  ("Total Events Processed   :  ",int(nevents))
    print  ("===================================================================")

    
    outf = ROOT.TFile( outfile, "update" )
    meta = ROOT.TTree( "metadata", "metadata informations" )
    n = array( "i", [ 0 ] )
    meta.Branch( "eventsProcessed", n, "eventsProcessed/I" )
    n[0]=nevents
    meta.Fill()
    p = ROOT.TParameter(int)( "eventsProcessed", n[0])
    p.Write()
    outf.Write()
    outf.Close()

