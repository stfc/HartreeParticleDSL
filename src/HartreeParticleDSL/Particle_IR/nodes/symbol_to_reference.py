'''
This module contains the symbol_to_reference function, which is a mapping
from a Symbol type to the associated Reference, enabling users to create
the correct type of Reference for a given symbol with:
    symbol_to_reference[type(symbol)](symbol)
'''

from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.symbols.arraysymbol import ArraySymbol
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol
from HartreeParticleDSL.Particle_IR.symbols.pointersymbol import PointerSymbol

from HartreeParticleDSL.Particle_IR.nodes.scalar_reference import ScalarReference
from HartreeParticleDSL.Particle_IR.nodes.array_reference import ArrayReference
from HartreeParticleDSL.Particle_IR.nodes.structure_reference import StructureReference
from HartreeParticleDSL.Particle_IR.nodes.pointer_reference import PointerReference

symbol_to_reference = {ScalarTypeSymbol: ScalarReference,
                       ArraySymbol: ArrayReference,
                       StructureSymbol: StructureReference,
                       PointerSymbol: PointerReference}
