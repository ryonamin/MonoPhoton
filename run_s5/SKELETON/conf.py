#!/usr/bin/python
import os

# samples will be chosen from the follwoing directory.
#===== for l5 sample
#LISTDIR="/hsm/ilc/grid/storm/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/2f_Z_nuNg/ILD_l5_o1_v02/v02-00-01"
#===== for s5 sample
LISTDIR="/hsm/ilc/grid/storm/prod/ilc/mc-opt-3/ild/dst-merged/500-TDR_ws/2f_Z_nuNg/ILD_s5_o1_v02/v02-00-01"

FILESUFFIX="slcio"
# samples will be chosen by searching the follwoing label in the file names.
PROCESSES=["nung"]
#PROCESSES=["bhabhang"]
#PROCESSES=["bhabhang","nung"]

# geometry file
#===== for l5 sample
#GEARFILE="/cvmfs/ilc.desy.de/sw/ILDConfig/v02-00-01/StandardConfig/production/Gear/gear_ILD_l5_v02.xml"
#===== for s5 sample
GEARFILE="/cvmfs/ilc.desy.de/sw/ILDConfig/v02-00-01/StandardConfig/production/Gear/gear_ILD_s5_v02.xml"

# number of input files per one steering file. 
nfilesInOneshot = 5; 

# directory that includes template xml files. 
TEMPLATEXMLORIGDIR = os.environ['MPDIR'] + "/run_s5/XML_TMPLATES"
#===== samples includes MCParticle (l5, s5 sample)
TEMPLATEFILE       = "monophoton.xml"

#===== for l5 sample
#OUTPUT_PREFIX      = "l5_500GeV."
#===== for s5 sample
OUTPUT_PREFIX      = "s5_500GeV."


# no need to change below
# common
LOGDIR      = "log"
XMLDIR      = "generatedXMLs"
TEMPLATEDIR = "templateXMLs"
OUTDIR_ROOT = "root"
SPACE       = "                 "
