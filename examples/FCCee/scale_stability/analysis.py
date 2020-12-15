import sys
import ROOT

print ("Load cxx analyzers ... ",)
ROOT.gSystem.Load("libedm4hep")
ROOT.gSystem.Load("libpodio")
ROOT.gSystem.Load("libFCCAnalyses")
#ROOT.gSystem.Load("libmyFCCAnalyses")

ROOT.gErrorIgnoreLevel = ROOT.kFatal
_edm  = ROOT.edm4hep.ReconstructedParticleData()
_pod  = ROOT.podio.ObjectID()
_fcc  = ROOT.getMC_px
_fcc2  = ROOT.getRP2MC_p

print ('edm4hep  ',_edm)
print ('podio    ',_pod)
print ('fccana   ',_fcc)
print ('fccana2  ',_fcc2)
#ROOT.ROOT.EnableThreadSafety()
#ROOT.ROOT.EnableImplicitMT(1)
#ROOT.TTree.SetMaxTreeSize(100000000000)
class analysis():

    #__________________________________________________________
    def __init__(self, inputlist, outname, ncpu):
        self.outname = outname
        if ".root" not in outname:
            self.outname+=".root"

        #ROOT.ROOT.EnableImplicitMT(ncpu)

        self.df = ROOT.RDataFrame("events", inputlist)
        print (" done")
    #__________________________________________________________
    def run(self):
        
        #df2 = (self.df.Range(100)
        df2 = (self.df

               #.Alias("Muon0", "Muon#0.index")
               #.Define("muons",  "getRP(Muon0, ReconstructedParticles)")

               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               .Alias("Particle1", "Particle#1.index")
               .Alias("Particle0", "Particle#0.index")

               # RecoParticles that are associated with a MC muon
               .Define("allmuons",  "selRP_PDG(13, true)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               # Keep only those with E > 2 GeV (otherwise, won't reach the muon detector)
               .Define("muons",   "selRP_E(2.)( allmuons ) ")      
               # make all dimuon combinations from them: 
               .Define("dimuons",         "Pairs(true)(muons, muons)")
               .Define("dimuons_mass",    "getRP_mass(dimuons)")
               .Define("dimuons_charge",    "getRP_charge(dimuons)")

               # for the fakes: first define the charged hadrons
               .Define("ChargedHadrons", "selRP_ChargedHadrons( MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               # Only the ones with  E > 2 GeV for the fake muons :
               .Define("ChargedHadrons_Egt2",  "selRP_E(2.) ( ChargedHadrons )")
               # fake muons :
               .Define("fakeMuons_1em3", "selRP_Fakes( 1e-3, 0.106)(ChargedHadrons_Egt2)" )
               .Define("fakeMuons_1em2", "selRP_Fakes( 1e-2, 0.106)(ChargedHadrons_Egt2)" )
               .Define("fakeMuons_5em2", "selRP_Fakes( 5e-2, 0.106)(ChargedHadrons_Egt2)" )

               # fake muons, with two hypotheses for the fake rate:
               #.Define("fakeMuons_1em3",  "selRP_FakeMuons(1e-3)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               #.Define("fakeMuons_1em2",  "selRP_FakeMuons(1e-2)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               #.Define("fakeMuons_5em2",  "selRP_FakeMuons(5e-2)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")

               .Define("muons_with_fakes_1em3",  "my_mergeParticles( muons, fakeMuons_1em3)")
               .Define("muons_with_fakes_1em2",  "my_mergeParticles( muons, fakeMuons_1em2)")
               .Define("muons_with_fakes_5em2",  "my_mergeParticles( muons, fakeMuons_5em2)")
               # make all dimuon combinations, alsoincluding now the fake muons
               .Define("dimuons_with_fakes_1em3",  "Pairs(true)(muons_with_fakes_1em3, muons_with_fakes_1em3)")
               .Define("dimuons_with_fakes_1em2",  "Pairs(true)(muons_with_fakes_1em2, muons_with_fakes_1em2)")
               .Define("dimuons_with_fakes_5em2",  "Pairs(true)(muons_with_fakes_5em2, muons_with_fakes_5em2)")

               .Define("dimuons_with_fakes_1em3_mass",  "getRP_mass(dimuons_with_fakes_1em3)")
               .Define("dimuons_with_fakes_1em3_charge",  "getRP_charge(dimuons_with_fakes_1em3)")
               .Define("dimuons_with_fakes_1em2_mass",  "getRP_mass(dimuons_with_fakes_1em2)")
               .Define("dimuons_with_fakes_1em2_charge",  "getRP_charge(dimuons_with_fakes_1em2)")
               .Define("dimuons_with_fakes_5em2_mass",  "getRP_mass(dimuons_with_fakes_5em2)")
               .Define("dimuons_with_fakes_5em2_charge",  "getRP_charge(dimuons_with_fakes_5em2)")

               # keep only the combination closest in mass to the JPsi:
               #.Define("dimuons_3p2",   "ResonanceBuilder(23, 3.2)(muons)")
               #.Define("dimuons_3p2_mass",  "getRP_mass(dimuons_3p2)")
               #.Define("dimuons_3p2_charge", "getRP_charge(dimuons_3p2)")

               #.Define("muons_pt", "getRP_pt(muons)")
               #.Define("MCmuons","selMC_PDG(13,  true)(Particle)")
               #.Define("MCmuons_theta","getMC_theta(MCmuons")

               
               # RecoParticles associated with the MC muons that come from JPsi -> mumu
               .Define("muons_from_JPsi", "selMuons_JPsimatch(443, 13, -13)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle,Particle1)")
               .Define("muons_from_JPsi_pt", "getRP_pt( muons_from_JPsi )")
               # The JPsis below are only the combinations for which both muons made a RecoParticle
               .Define("jpsi",         "JPsis()(muons_from_JPsi)")
               .Define("jpsi_mass", "getRP_mass(jpsi)")
               .Define("jpsi_charge", "getRP_charge(jpsi)")

               # select RecoParticles associated with the (K,pi) froma D0 decay
               .Define("KPi_from_D0", "selMuons_JPsimatch( 421, -321, 211)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle,Particle1)")
               .Define("KPi_from_D0bar", "selMuons_JPsimatch( -421, 321, -211)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle,Particle1)")
               .Define("KPi_from_D0_or_D0bar",  "my_mergeParticles( KPi_from_D0, KPi_from_D0bar)")

               #.Define("KPi_from_D0_track_phi",  "getRP2TRK_phi( KPi_from_D0 , EFlowTrack_1 )")
               #.Define("KPi_from_D0_track_tanLambda",   "getRP2TRK_tanLambda( KPi_from_D0 , EFlowTrack_1)")
               #.Define("KPi_from_D0_track_omega",   "getRP2TRK_omega( KPi_from_D0 , EFlowTrack_1)" )
               #.Define("D0",  "D0s()( KPi_from_D0, KPi_from_D0_track_phi, KPi_from_D0_track_tanLambda, KPi_from_D0_track_omega )")
               .Define("D0",  "JPsis()(KPi_from_D0_or_D0bar)")    # yes, the name should be changed..
               .Define("D0_mass",   "getRP_mass(D0)")
               .Define("D0_charge", "getRP_charge(D0)")

               # RecoParticles that are associated with a MC Kaon
               .Define("RP_Kaons",  "selRP_PDG(321, true)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               # RecoParticles that are associated with a MC pi+/-
               .Define("RP_Pions",  "selRP_PDG(211, true)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               # D0 candidates made from these RecoParticles - i.e. assume perfect  PID
               .Define("D0cand",  "Pairs(false)(RP_Kaons, RP_Pions)")
               # Select only the ones of interest  
               .Define("D0cand_in_massWindow", "selRP_mass( 1.5, 2.3) (D0cand)" )
               .Define("D0cand_mass",   "getRP_mass( D0cand_in_massWindow )")
               .Define("D0cand_charge",   "getRP_charge( D0cand_in_massWindow )")
          
               # for D0 candidates with no PID at all :
               # Assign all charged hadrons to kaons :
               .Define("RP_Kaons_noPID",  "selRP_Fakes( 1., 0.494)(ChargedHadrons)" )
               # Assign all charged hadrons to pions :
               .Define("RP_Pions_noPID",  "selRP_Fakes( 1., 0.140)(ChargedHadrons)" )
               # D0 candidates :
               .Define("RP_KPis",  "Pairs(true)(RP_Kaons_noPID, RP_Pions_noPID)" )
               .Define("RP_PiKs",  "Pairs(true)(RP_Pions_noPID, RP_Kaons_noPID)" )
               .Define("D0cand_noPID",  "my_mergeParticles( RP_KPis, RP_PiKs)" )
               .Define("D0cand_noPID_in_massWindow",  "selRP_mass( 1.5, 2.3)(D0cand_noPID)")
               .Define("D0cand_noPID_mass",   "getRP_mass( D0cand_noPID_in_massWindow )")
               .Define("D0cand_noPID_charge",   "getRP_charge( D0cand_noPID_in_massWindow )")

               )

        # select branches for output file
        branchList = ROOT.vector('string')()
        for branchName in [

                #"muons",
                #"dimuons",
                "dimuons_mass",
                "dimuons_charge",
                #"dimuons_with_fakes_1em3",
                #"dimuons_with_fakes_1em2",
                #"dimuons_with_fakes_5em2",
                "dimuons_with_fakes_1em3_mass",
                "dimuons_with_fakes_1em3_charge",
                "dimuons_with_fakes_1em2_mass",
                "dimuons_with_fakes_1em2_charge",
                "dimuons_with_fakes_5em2_mass",
                "dimuons_with_fakes_5em2_charge",
                #"fakeMuons_1em3",
                #"fakeMuons_1em2",
                #"muons_with_fakes_1em2",
                #"muons_with_fakes_1em3",
                #"MCmuons_theta",
                #"muonsIdx_from_JPsi",
                #"muons_pt",
                #"muons_from_JPsi",
                #"muons_from_JPsi_pt",
                #"jpsi",
                "jpsi_mass",
                "jpsi_charge",
	        #
                #"D0",
                "D0_mass",
                "D0_charge",
                #"KPi_from_D0",
                #"KPi_from_D0bar",
                #"D0cand_in_massWindow_mass"
                "D0cand_mass",
                "D0cand_charge",
                "D0cand_noPID_mass",
                "D0cand_noPID_charge"
                
                ]:
            branchList.push_back(branchName)

        opts = ROOT.RDF.RSnapshotOptions()
        opts.fCompressionAlgorithm = ROOT.ROOT.kLZ4
        opts.fCompressionLevel = 3
        opts.fAutoFlush = -1024*1024*branchList.size()
        #df2.Snapshot("events", self.outname, branchList, opts)
        df2.Snapshot("events", self.outname, branchList)

# example call for standalone file
# python FCCeeAnalyses/Z_Zbb_Flavor/dataframe/analysis.py /eos/experiment/fcc/ee/generation/DelphesEvents/fcc_tmp/p8_ee_Ztautau_ecm91/events_012154460.root

if __name__ == "__main__":

    if len(sys.argv)==1:
        print ("usage:")
        print ("python ",sys.argv[0]," file.root")
        sys.exit(3)
    infile = sys.argv[1]
    outDir = 'FCCee/'+sys.argv[0].split('/')[1]+'/'
    import os
    os.system("mkdir -p {}".format(outDir))
    outfile = outDir+infile.split('/')[-1]
    ncpus = 0
    analysis = analysis(infile, outfile, ncpus)
    analysis.run()

    tf = ROOT.TFile(infile)
    entries = tf.events.GetEntries()
    p = ROOT.TParameter(int)( "eventsProcessed", entries)
    outf=ROOT.TFile(outfile,"UPDATE")
    p.Write()
