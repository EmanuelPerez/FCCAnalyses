import ROOT

#
#
#       STILL UNDER SEVELOPMENT
#
#
#


# global parameters
intLumi        = 100.0e+06 #in pb-1
ana_tex        = "e^{+}e^{-} #rightarrow qq"
delphesVersion = "3.4.3-pre06"
energy         = 91.0
collider       = "FCC-ee"
inputDir       = "FCCee/JPsi/"
formats        = ['png','pdf']
yaxis          = ['lin','log']
stacksig       = ['stack','nostack']
outdir         = 'FCCee/JPsi/plots/'

variables = ['m_dimu', 'm_jpsi', 'm_dimu_zoom', 'm_jpsi_zoom' ]

###Dictonnary with the analysis name as a key, and the list of selections to be plotted for this analysis. The name of the selections should be the same than in the final selection
selections = {}
selections['dimuOS'] = ["dimuOS"]
selections['dimuSS'] = ["dimuSS"]

extralabel = {}
extralabel['dimuOS'] = "OS muons"
extralabel['dimuSS'] = "SS muons"


colors = {}
colors['Zbb'] = ROOT.kRed
colors['Zcc'] = ROOT.kBlue+1
colors['Zuds'] = ROOT.kGreen+2

plots = {}
plots['dimuOS'] = {'signal':{'Zbb':['p8_ee_Zbb_ecm91']},
                   'backgrounds':{'Zuds':['p8_ee_Zuds_ecm91']}
                }

plots['dimuSS'] = {'signal':{'Zbb':['p8_ee_Zbb_ecm91']},
                   'backgrounds':{'Zuds':['p8_ee_Zuds_ecm91']}
                }



legend = {}
legend['Zbb'] = 'Zbb'
legend['Zcc'] = 'Zcc'
legend['Zuds'] = 'Zuds'




