from config.common_defaults import deffccdicts

 #python FCCeeAnalyses/ZH_Zmumu/dataframe/preSel.py 
import os

basedir=os.path.join(os.getenv('FCCDICTSDIR', deffccdicts), '') + "yaml/FCCee/fcc_tmp/"
outdir="FCCee/scale_stability/"
NUM_CPUS = 15
process_list=['p8_ee_Zbb_ecm91','p8_ee_Zuds_ecm91','p8_ee_Zcc_ecm91']
#process_list=['p8_ee_Zuds_ecm91','p8_ee_Zcc_ecm91']
#process_list=['p8_ee_Zbb_ecm91']
fraction=0.02

import config.runDataFrame as rdf
myana=rdf.runDataFrame(basedir,process_list)
myana.run(ncpu=NUM_CPUS,fraction=fraction,outDir=outdir)


