
import sys

import heppy
import FCCeeAnalyses

del sys.modules["heppy"]
del sys.modules["FCCeeAnalyses"]
sys.modules["fccflatprod"] = __import__("fccflatprod")
sys.modules["heppy"] = __import__("fccflatprod")
sys.modules["FCCeeAnalyses"] = __import__("fccflatprod")

import ROOT

#TODO: change library name
ROOT.gSystem.Load("fcc_ana_ZH_Zmumu_cxx")


if __name__ == "__main__":
  execfile(sys.argv[1])

print [a for a in sequence.the_sequence]


from ROOT import select_leptons_add_to_dataframe 
from ROOT import recoil_add_to_dataframe 
from ROOT import initial_dataframe_convert



df = ROOT.RDataFrame("events", comp.files[0])


df = initial_dataframe_convert(df)




df = select_leptons_add_to_dataframe(df, "muons", "muonITags", "selected_muons")
print df

df = recoil_add_to_dataframe(df, "muons",  "recoil")
print df
print df.GetDefinedColumnNames()
