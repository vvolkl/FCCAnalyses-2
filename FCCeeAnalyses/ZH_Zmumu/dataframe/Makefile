

# build the reproducer of the heppy analysis with TDataFrame
# todo: cleaner directory structure, cmake, etc

all: tdataframe_analysis

tdataframe_analysis:  fcc_ana_ZH_Zmumu.cxx
	g++ -I/cvmfs/fcc.cern.ch/sw/views/releases/externals/94.1.0/x86_64-centos7-gcc62-opt/include -g -o fcc_ana_ZH_Zmumu fcc_ana_ZH_Zmumu.cxx `root-config --cflags --glibs` -lROOTVecOps -lROOTDataFrame   


# run the analysis

run:
	./fcc_ana_ZH_Zmumu tree.root root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/fcc_v01/p8_ee_ZH_ecm240/events_067271769.root

# todo: comparison with expected output  (/afs/cern.ch/work/v/vavolkl/public/heppy_analysis_samples/localtest/mgp8_pp_tth01j_5f_hllll_1/heppy.FCChhAnalyses.FCChh.ZH_Zmumu.TreeProducer.TreeProducer_1/tree.root)
