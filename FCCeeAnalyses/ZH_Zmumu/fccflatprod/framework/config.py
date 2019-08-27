class Component():
    def __init__(self, name, files, tree_name=None, triggers=None, **kwargs):
      print name

class Analyzer:
  def __init__(self, the_analyzer, *args, **kwargs):
    print args
    print kwargs
    self.the_analyzer = the_analyzer(*args, **kwargs)
    pass

class Sequence():
  def __init__(self, the_sequence):
    self.the_sequence = the_sequence
    pass

class Config():
  def __init__(*args, **kwargs):
    pass

class MCComponent():
  def __init__(self, name , files, **kwargs):
    pass
      

