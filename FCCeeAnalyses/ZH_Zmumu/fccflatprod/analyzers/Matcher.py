
import ROOT
from ROOT import noMatchJets_add_to_dataframe

class Matcher:
  def __init__(*args, **kwargs):
    self.delta_r = kwargs["delta_r"]
    self.match_particles = kwargs["match_particles"]
    self.particles = kwargs["particles"]

  def doit(self, dataframe):
    noMatchJets_add_to_dataframe(dataframe, self.particles)
    pass

