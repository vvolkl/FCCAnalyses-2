import ROOT
from ROOT import select_leptons_add_to_dataframe


class Selector:
  def __init__(self, *args, **kwargs):
    self.output = kwargs["output"]
    self.input_objects =  kwargs["input_objects"]
    self.filter_func = kwargs["filter_func"]

  def doit(self, dataframe):
    #select_leptons_add_to_dataframe(dataframe)
    pass
    
    
