'''
This module contains the datatype_to_symbol map, which contains a mapping from
Type to Symbol, enabling creation of the appropriate Symbol for a given Type.
'''

from HartreeParticleDSL.Particle_IR.datatypes.datatype import ScalarType, \
        StructureType, PointerType, ArrayType

from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.symbols.structuresymbol import StructureSymbol
from HartreeParticleDSL.Particle_IR.symbols.pointersymbol import PointerSymbol
from HartreeParticleDSL.Particle_IR.symbols.arraysymbol import ArraySymbol

datatype_to_symbol = {ScalarType: ScalarTypeSymbol,
                      StructureType: StructureSymbol,
                      PointerType: PointerSymbol,
                      ArrayType: ArraySymbol}
