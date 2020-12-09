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
        
        #df2 = (self.df.Range(1000)
        df2 = (self.df

               #.Alias("Muon0", "Muon#0.index")
               #.Define("muons",  "getRP(Muon0, ReconstructedParticles)")

               .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
               .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
               .Alias("Particle1", "Particle#1.index")
               .Alias("Particle0", "Particle#0.index")

               # RecoParticles that are associated with a MC muon
               .Define("allmuons",  "selRP_PDG(13)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               # Keep only those with E > 2 GeV (otherwise, won't reach the muon detector)
               .Define("muons",   "selRP_E(2.)( allmuons ) ")      
               # make all dimuon combinations from them: 
               .Define("dimuons",         "Dimuons()(muons)")
               .Define("dimuons_mass",    "getRP_mass(dimuons)")
               .Define("dimuons_charge",    "getRP_charge(dimuons)")

               # fake muons, with two hypotheses for the fake rate:
               .Define("fakeMuons_1em3",  "selRP_FakeMuons(1e-3)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               .Define("fakeMuons_1em2",  "selRP_FakeMuons(1e-2)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")
               .Define("fakeMuons_5em2",  "selRP_FakeMuons(5e-2)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle)")

               .Define("muons_with_fakes_1em3",  "my_mergeParticles( muons, fakeMuons_1em3)")
               .Define("muons_with_fakes_1em2",  "my_mergeParticles( muons, fakeMuons_1em2)")
               .Define("muons_with_fakes_5em2",  "my_mergeParticles( muons, fakeMuons_5em2)")
               # make all dimuon combinations, alsoincluding now the fake muons
               .Define("dimuons_with_fakes_1em3",  "Dimuons()(muons_with_fakes_1em3)")
               .Define("dimuons_with_fakes_1em2",  "Dimuons()(muons_with_fakes_1em2)")
               .Define("dimuons_with_fakes_5em2",  "Dimuons()(muons_with_fakes_5em2)")

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
               .Define("muons_from_JPsi", "selMuons_JPsimatch(1)(MCRecoAssociations0,MCRecoAssociations1,ReconstructedParticles,Particle,Particle1)")
               .Define("muons_from_JPsi_pt", "getRP_pt( muons_from_JPsi )")
               # The JPsis below are only the combinations for which both muons made a RecoParticle
               .Define("jpsi",         "JPsis()(muons_from_JPsi)")
               .Define("jpsi_mass", "getRP_mass(jpsi)")
               .Define("jpsi_charge", "getRP_charge(jpsi)")

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
                "jpsi_charge"
                
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
