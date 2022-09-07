import pytest
from HartreeParticleDSL.Particle_IR.symbols.scalartypesymbol import ScalarTypeSymbol
from HartreeParticleDSL.Particle_IR.symbols.symbol import Symbol
from HartreeParticleDSL.Particle_IR.datatypes.datatype import ScalarType
from HartreeParticleDSL.Particle_IR.datatypes.datatype import DOUBLE_TYPE, FLOAT_TYPE
from HartreeParticleDSL.Particle_IR.symbols.symboltable import SymbolTable
from HartreeParticleDSL.Particle_IR.nodes.kern import Kern
from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError

from HartreeParticleDSL.HartreeParticleDSL import _HartreeParticleDSL

def test_symbol_table():
    class test_kern(Kern):
        pass

    with pytest.raises(TypeError) as excinfo:
        a = SymbolTable("123")
    assert ("Argument 'kern' should be of type HartreeParticleDSL or Kern, "
            "but instead got <class 'str'>." in str(excinfo.value))

    tks = test_kern()

    sym_tab = SymbolTable(tks)

    assert sym_tab.is_empty()

    s1 = ScalarTypeSymbol("a", DOUBLE_TYPE)
    s2 = ScalarTypeSymbol("b", DOUBLE_TYPE)

    sym_tab.add(s1)
    sym_tab.add(s2)

    assert sym_tab.lookup("a") is s1
    assert sym_tab.get_symbols()["b"] is s2

    s3 = ScalarTypeSymbol("a", FLOAT_TYPE)
    with pytest.raises(IRGenerationError) as excinfo:
        sym_tab.add(s3)
    assert ("Tried to add a new symbol a "
            "but it was already present in the symbol table."
            in str(excinfo.value))

    with pytest.raises(IRGenerationError) as excinfo:
        sym_tab.add("not sym")
    assert ("Symbol not sym is not a symbol, but <class 'str'>.") in str(excinfo.value)

    s4 = sym_tab.new_symbol("c", FLOAT_TYPE, ScalarTypeSymbol)
    assert sym_tab.lookup("c") is s4

    assert sym_tab.find_or_create("c", FLOAT_TYPE, ScalarTypeSymbol) is s4

    with pytest.raises(IRGenerationError) as excinfo:
        sym_tab.find_or_create("c", DOUBLE_TYPE, ScalarTypeSymbol)
    assert ("Found a symbol with specified name c, but "
            "it had datatype Scalar<FLOAT, SINGLE> which doesn't "
            "match requested datatype Scalar<FLOAT, DOUBLE>.") in str(excinfo.value)

    assert not sym_tab.is_empty()

    z = sym_tab.find_or_create("new", FLOAT_TYPE, ScalarTypeSymbol)
    assert sym_tab.lookup("new") is z

    x = SymbolTable(_HartreeParticleDSL.get_instance())
    assert x._is_global == True

def test_symbol_table_errors():
    class test_kern(Kern):
        pass
    tks = test_kern()

    sym_tab = SymbolTable(tks)

    with pytest.raises(TypeError) as excinfo:
        a = sym_tab.new_symbol("a", FLOAT_TYPE, "b")
    assert ("The symbol_type parameter should be a type class of a subclass of"
            " Symbol, but found <class 'str'>.")
    with pytest.raises(TypeError) as excinfo:
        a = sym_tab.new_symbol(123, FLOAT_TYPE, ScalarTypeSymbol)

    sym_tab._is_global = True
    s4 = sym_tab.new_symbol("c", FLOAT_TYPE, ScalarTypeSymbol)
    assert s4.visibility == Symbol.Visibility.GLOBAL

    with pytest.raises(TypeError) as excinfo:
        sym_tab.lookup(123)
    assert ("Expected the name argument to be a string, but got <class 'int'>."
            in str(excinfo.value))
