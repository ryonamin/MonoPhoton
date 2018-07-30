#!/usr/bin/python
import os

# samples will be chosen from the follwoing directory.
LISTDIR="/hsm/ilc/grid/storm/prod/ilc/mc-dbd/ild/dst-merged/500-TDR_ws"
FILESUFFIX="slcio"
# samples will be chosen by searching the follwoing label in the file names.
PROCESSES=["nung"]
#PROCESSES=["bhabhang"]
#PROCESSES=["bhabhang","nung"]

# geometry file
#===== for DBD sample
GEARFILE="/cvmfs/ilc.desy.de/sw/ILDConfig/v01-16-p10_250/StandardConfig/current/GearOutput.xml"
#===== for l5 sample
#GEARFILE="/cvmfs/ilc.desy.de/sw/ILDConfig/v02-00-01/StandardConfig/production/Gear/gear_ILD_l5_v02.xml"
#===== for s5 sample
#GEARFILE="/cvmfs/ilc.desy.de/sw/ILDConfig/v02-00-01/StandardConfig/production/Gear/gear_ILD_s5_v02.xml"

# number of input files per one steering file. 
nfilesInOneshot = 5; 

# directory that includes template xml files. 
TEMPLATEXMLORIGDIR = os.environ['MPDIR'] + "/run_DBD/XML_TMPLATES"
#===== samples includes MCParticlesSkimmed (DBD sample)
TEMPLATEFILE       = "monophoton_mcskim.xml"
#===== samples includes MCParticle (l5, s5 sample)
#TEMPLATEFILE       = "monophoton.xml"

OUTPUT_PREFIX      = "dbd_500GeV."


# no need to change below
# common
LOGDIR      = "log"
XMLDIR      = "generatedXMLs"
TEMPLATEDIR = "templateXMLs"
OUTDIR_ROOT = "root"
SPACE       = "                 "
