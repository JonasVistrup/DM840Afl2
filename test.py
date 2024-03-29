import pydoStuff
from pydoStuff import *

# sanity check for multiple copies of libMØD
modValue = mod.magicLibraryValue()
ourValue = pydoStuff.magicLibraryValue()
if modValue != ourValue:
	print("mod =", modValue)
	print("our =", ourValue)
	raise Exception("Magic values differ! I.e., more than one instance of libmod has been loaded.")
# end of sanity check

# make doStuff a bit more friendly to use
_doStuff_orig = pydoStuff.doStuff
def _doStuff(educts, products, doChemistryCheck=True):
	return mod._unwrap(_doStuff_orig(
		mod._wrap(mod.libpymod._VecGraph, educts),
		mod._wrap(mod.libpymod._VecGraph, products),
		doChemistryCheck))
pydoStuff.doStuff = _doStuff
doStuff = _doStuff
# end of friendlyfier code

g1 = smiles("OCC=O")
g2 = smiles("OC=CO")
res = doStuff([g1], [g2])

printer = GraphPrinter()
#printer.setReactionDefault()
printer.collapseHydrogens = False
printer.withIndex = True

for a in res:
	a.print(printer)
