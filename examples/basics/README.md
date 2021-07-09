Basic examples
=============
This directory contains a number of examples each showcasing a specific functionality of the FCCAnalyses framework. It serves as a reference guide for how to implement specific common usecases or you can work through the examples one-by-one in order as a tutorial to familiarize yourself with the full functionality of the framework. 

Each example is a stand-alone script for demonstration purposes, and does not make assumptions on a specific physics case. To understand how to write a full analysis with the FCCAnalyses framework please have a look at (insert a link to documentation about code class-structure) - the examples here only illustrate specific technical functionalities. 

By calling `python <example>.py` you can run the specific example over the integrated test file found (add the testdata directory), and it will create a new directory in your current working directory with the name of the example to write the output to. If you prefer to run over your own input file or a different output directory you can run with options:

`python <example>.py -i <path_to_your_inputfile> -o <path_to_your_outputdir>`

Certain examples may have additional options, you can always check what options are available with `python <example>.py -h`. 

Table of contents
=================
  * [Prerequisites](#prerequisites)
    * [RDataFrame](#rdataframe)
    * [EDM4HEP event model](#edm4hep-event-model)
    * [Structure of EDM4HEP files](#structure-of-edm4hep-files)
  * [Reading objects from EDM4HEP](#reading-objects-from-edm4hep)
  * [Writing your own function](#writing-your-own-function)
    * [Inline](#inline)
    * [Using ROOT GInterpreter](#using-root-ginterpreter)
    * [Writing your own class](#writing-your-own-class)
  * [Base collection](#base-collection)

Prerequisites
=================

The FCCAnalyses framework is based on the [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html) interface which allows fast and efficient analysis of [ROOT's TTrees](https://root.cern/doc/master/classTTree.html) and on samples following the [EDM4HEP event data model](https://edm4hep.web.cern.ch/). Some brief explanations and links to further material on the two are given below, a basic understanding of both is necessary for using this framework to write your own analysis code. 

RDataFrame
=================
(to add)

EDM4HEP event model
=================
[Link to EDM4HEP class overview](https://edm4hep.web.cern.ch/namespaceedm4hep.html)

(to add brief intro/pointers)


Structure of EDM4HEP files
==========================

The content of an EDM4HEP file can be seen by opening it in ROOT, and by inspecting the content of the "events" tree with a TBrowser.
Example:
```
root /eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/wzp6_ee_mumuH_ecm240.root
root [0] TBrowser b
```
<img src="figs/browser_events.png" alt="drawing" width="480"/>

As shown in the screenshot above, there are two types of branches:

  - Branches without a pound (#) in their name:  Electron (1), Muon (2), EFlowNeutralHadron (3), Particle (4), Photon (5), ReconstructedParticles (6), EFlowPhoton (7), MCRecoAssociations (8), MissingET (9), ParticleIDs (10), Jet (11), EFlowTrack (12), EFlowTrack\_1 (13). They refer to collections of objects.
  - Branches with a pound in their name:  Each of the object collections listed above, e.g. "Collection", has up to six associated collections of references, 
    i.e. indices that point to another or to the same object collection. They are labeled Collection#i, with i = 0 ... 5. For example, the Muon collection has one single
    associated collection of references, Muon#0. 

To figure out which collection is pointed to by Muon#0 (or by any other collection of references), one can look at the value of Muon#0.collectionID (see screenshot below). 
The collectionID of Muon#0 is the collection number 6, which, in the list of "object collections" above, corresponds to the collection of ReconstructedParticles. 
Indeed, the Muon collection itself contains nothing (see screenshot below): all the information is contained in the ReconstructedParticles. The Muon collection,
together with Muon#0, just provides a convenient way to access, among the ReconstructedParticles, those that were identified as muons.

<img src="figs/browser_Muon0.png" alt="drawing" width="480"/>

The same holds for the Electron and Photon collections. On the other hand, the MissingET collection is already a ReconstructedParticle, as can be seen
by inspecting it in the browser:

<img src="figs/browser_missingET.png" alt="drawing" width="480"/>

The "Particle" collection corresponds to the Monte-Carlo particles. It has two associated collections of references, Particle#0 and Particle#1. As can
be seen by looking at their collectionID, they both point to collection number 4, i.e.  to the Particle collection itself. Particle#0 and
Particle#1 contain, respectively, links to the parents and tp the daughters of the MC particles - as can be seen in the [edm4hep yaml description here](https://github.com/key4hep/EDM4hep/blob/master/edm4hep.yaml#L118-L119).
Examples will be given below, showing how to navigate through the Monte-Carlo record using Particle, Particle#0 and Particle#1.




Reading objects from EDM4HEP
=============
The example read_EDM4HEP.py shows you how to access the different objects such as jets, electrons, muons, missing ET etc. from the EDM4HEP files. Generally a new variable is calculated with a statement inside the `run()` function of the `analysis` class like
`dataframe.Define("<your_variable>", "<accessor_fct (<name_object>)>")`
which creates a column in the RDataFrame named `<your_variable>` and filled with the return value of the `<accessor_fct>` for the given object. 

Here, accessor functions are the functions found in the C++ analyzers code that return a certain variable. Since the analyzers code defines a specific namespace for each module, such as ReconstructedParticle or MCParticle, the full accessor function call looks like `<namespace>::<function_name>(object)`. To access the pT of a reconstructed object you would therefore call `ReconstructedParticle::get_pt(object)` and for a MC-particle the call would be `MCParticle::get_pt(object)`. The namespace corresponds to the file name of the C++ code, making it clear where to look for the source code if you have a question about the internal workings of one such functions. 

Below you can find an overview of the basic, most commonly required functions, to illustrate the naming conventions. This is not an exhaustive list, if you want to find out all functions that are available please take a look in the respective analyzers code itself - [here for reconstructed particles](https://github.com/HEP-FCC/FCCAnalyses/blob/master/analyzers/dataframe/ReconstructedParticle.h) and [here for MC particles](https://github.com/HEP-FCC/FCCAnalyses/blob/master/analyzers/dataframe/MCParticle.h).


| Variable  | Function name | Available for | 
| ------------- | ------------- | ------------- |
| Transverse momentum  | `get_pt(object)`  | `MCParticle`, `ReconstructedParticle` |
| Pseudorapidity  | `get_eta(object)`  | `MCParticle`, `ReconstructedParticle` |
| Energy  | `get_e(object)`  | `MCParticle`, `ReconstructedParticle` |
| Mass  | `get_mass(object)`  | `MCParticle`, `ReconstructedParticle` |
| Charge  | `get_charge(object)`  | `MCParticle`, `ReconstructedParticle` |
| Number (in event)  | `get_n(object)`  | `MCParticle`, `ReconstructedParticle` |
| PDG ID  | `get_pdg(object)`  | `MCParticle` |


If you want to add your own function have a look at the [Writing your own function](#writing-your-own-function) section on this page. 

For the name of the object, in principle the names of the EDM4HEP collections are used - photons, muons and electrons are an exception where a few extra steps are required, as shown in the example here. 

This example also shows how to apply object selection cuts, for example selecting only reconstructed objects with a transverse momentum pT larger than a given threshold by using the `ReconstructedParticle::sel_pt(<threshold>)(<name_object>)` function. 

In the end of the example you can see how the selected variables are written to branches of the output n-tuple, using the `dataframe.Snapshot("<tree_name>", <branch_list> )`, where in all examples here the name of the output-tree is always `events` and the branch_list is defined as a `ROOT.vector('string')` as demonstrated in the example. Note that branches of variables that exist multiple times per event, i.e. anything derived from a collection such as the pT of jets, result in vector branches. This is also true for some quantities that in principle only exist once per event, but are collections in the EDM4HEP format, such as the missing transverse energy. 

<!--
======================================================

GUIDELINES and REQUEST  for new instructions to be added 


FACT: The general newcomer planning to do a  "case study" needs first to have instruction about how to unpack the variable he/she needs from the EDM4HEP into the flat ntuple. 
Very likely he/she has no idea about what the analysis will look like yet.  

More details about how this can be organizes: 

-- introduce links to  the manual for the "data frame" framework from root. 
Since there is a large amount of information and it is easy to get lost please add pointers to specific pages/examples. 

-- have simple examples on how to extract the variables from the edm4hep and add them to the branch in the analysis.py routine: 
    * in the case of a simple basic objects (i.e. jets, electrons, muons...) 
    * in the case of derived objects that need to be calculated: like jet pt
    * where do we find the functions needed already available for a specific object 

-- where in the code do we add NEW functions to obtain new variables/corrections  to be used in analysis.py to create our variables 
    * for instance a recent one: vertexing and new vertices information 
    * invariant masses 
    * etc ...
-- some of this exist but it is buried in the different places, I would suggest to streamlined it in a "fake" example that contains a bit of everything in the same place. 

-- explain more clearly the interplay of the "analysis.py" code and the "FinalSel.py" code in terms of the possibility of having 
calculations done in one or the other. With a couple of concrete examples. 

-- explain which files need to be modified, in case someone wants to use a personal "source file" not in the official repository.  
      * For example, you may need some code that is not generic at all, but very just used by your analysis. 
      So, instead of polluting the existing files MCParticle.cc etc, you may want to
      put this code into some MyStuff.cc.
      * explain how to compile it to use it into the analysis (add it into  the CmakeLists.txt for instance. Not obvious for everyone)

-->