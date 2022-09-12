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
