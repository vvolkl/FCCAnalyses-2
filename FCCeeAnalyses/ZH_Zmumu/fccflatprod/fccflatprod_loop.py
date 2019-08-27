
import sys

import heppy
import FCCeeAnalyses

del sys.modules["heppy"]
del sys.modules["FCCeeAnalyses"]
sys.modules["fccflatprod"] = __import__("fccflatprod")
sys.modules["heppy"] = __import__("fccflatprod")
sys.modules["FCCeeAnalyses"] = __import__("fccflatprod")


if __name__ == "__main__":
  execfile(sys.argv[1])

