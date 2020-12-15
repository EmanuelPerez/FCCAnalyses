from config.common_defaults import deffccdicts

#
#
#	STILL UNDER SEVELOPMENT
#
#
#

#python FCCeeAnalyses/ZH_Zmumu/dataframe/finalSel.py 
import sys, os
import ROOT

###Input directory where the files produced at the pre-selection level are
baseDir  = "FCCee/scale_stability/"

###Link to the dictonary that contains all the cross section informations etc...
procDict = os.path.join(os.getenv('FCCDICTSDIR', deffccdicts), '') + "FCCee_procDict_fcc_tmp.json"

#process_list=['p8_ee_Zbb_ecm91','p8_ee_Zuds_ecm91','p8_ee_Zcc_ecm91']
process_list=['p8_ee_Zbb_ecm91' ]

###Optinally Define new variables
#define_list = {"dimuons_mass10": "dimuons_mass*10",
#               "jpsi_mass10":"jpsi_mass*10"}



ROOT.gInterpreter.Declare("""
bool OppositeSign(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> dimuons) {
    if (dimuons.size()<1) return false;
    for (size_t i = 0; i < dimuons.size(); i++) {
       if( dimuons.at(i).charge != 0  )
        return false;
}
    return true;
}
""")

###Dictionnay of the list of cuts. The key is the name of the selection that will be added to the output file
#cut_list = {"sel_dimuOS":"dimuons_charge ==0",
#            "sel_dimuSS":"dimuons_charge  != 0",
#            "sel_jpsiOS":"jpsi_charge  == 0",
#            "sel_jpsiSS":"jpsi_charge  != 0",
#}
cut_list ={ 
            #"all":"1 > 0",
            #"dimuOS":"dimuons_charge.size() ==1 && dimuons_charge[0] == 0",
            #"dimuSS":"dimuons_charge.size() ==1 && dimuons_charge[0] != 0",
            #"jpsiOS":"jpsi_charge.size() == 1 &&  jpsi_charge[0] ==0",
            #"jpsiSS":"jpsi_charge.size() == 1 && jpsi_charge[0] !=0"
            #"sel3":"myFilter(zed_leptonic_m)"
            "dimuOS":"OppositeSign( dimuons_charge )"
            #"dimuOS":"dimuons_3p2_charge.size() == 1 && dimuons_3p2_charge[0]  == 0"
}



###Dictionary for the ouput variable/hitograms. The key is the name of the variable in the output files. "name" is the name of the variable in the input file, "title" is the x-axis label of the histogram, "bin" the number of bins of the histogram, "xmin" the minimum x-axis value and "xmax" the maximum x-axis value.
variables = { 
    #"m_dimu":{"name":"dimuons_mass","title":"m_{#mu#mu} [GeV]","bin":100,"xmin":0,"xmax":10},
    "m_dimu":{"name":"dimuons.mass","title":"m_{#mu#mu} [GeV]","bin":100,"xmin":0,"xmax":10},
    #"m_jpsi":{"name":"jpsi_mass","title":"m_{J#psi} [GeV]","bin":100,"xmin":0,"xmax":10},
    #"m_dimu_zoom":{"name":"dimuons_mass10","title":"m_{#mu#mu} [GeV]","bin":100,"xmin":28,"xmax":35},
    #"m_jpsi_zoom":{"name":"jpsi_mass10","title":"m_{J#psi} [GeV]","bin":100,"xmin":28,"xmax":35},
}

###Number of CPUs to use
NUM_CPUS = 4

###This part is standard to all analyses
import config.runDataFrameFinal as rdf
myana=rdf.runDataFrameFinal(baseDir,procDict,process_list,cut_list,variables)
myana.run(ncpu=NUM_CPUS)

