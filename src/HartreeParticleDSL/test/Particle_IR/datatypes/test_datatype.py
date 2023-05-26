import pytest

from HartreeParticleDSL.Particle_IR.datatypes.datatype import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError
from psyclone.psyir.symbols import Symbol

def test_notype():
    a = NoType()
    assert str(a) == "NoType"

def test_scalar_type():

    with pytest.raises(TypeError) as excinfo:
        ScalarType(intrinsic="bad input", precision=32)
    assert ("ScalarType expected 'intrinsic' argument to be of "
            "type ScalarType.Intrinsic but found 'str'." in str(excinfo.value))
    with pytest.raises(TypeError) as excinfo:
        ScalarType(intrinsic=ScalarType.Intrinsic.INTEGER, precision="bad input")
    assert ("ScalarType expected 'precision' argument to be of type int, "
            "ScalarType.Precision or DataSymbol, but found 'str'." in str(excinfo.value))

    scalartype = ScalarType(intrinsic=ScalarType.Intrinsic.INTEGER,
                            precision=ScalarType.Precision.DOUBLE)
    scalartype_alt = ScalarType(ScalarType.Intrinsic.REAL, 8)

    assert scalartype.intrinsic == ScalarType.Intrinsic.INTEGER
    assert scalartype.precision == ScalarType.Precision.DOUBLE

    assert str(scalartype) == "Scalar<INTEGER, DOUBLE>"
    assert str(scalartype_alt) == "Scalar<REAL, 8>"

def test_structure_type():
    struct = StructureType()

    scalartype = ScalarType(intrinsic=ScalarType.Intrinsic.INTEGER,
                            precision=ScalarType.Precision.DOUBLE)
    scalartype_alt = ScalarType(ScalarType.Intrinsic.REAL, 8)

    struct.add("myint64", scalartype, Symbol.Visibility.PUBLIC)
    struct.add("myfloat64", scalartype_alt, Symbol.Visibility.PUBLIC)

    with pytest.raises(TypeError) as excinfo:
        struct.add(32, scalartype, Symbol.Visibility.PUBLIC)
    assert ("The name of a component of a StructureType must be a 'str' but got 'int'" in str(excinfo.value))
    with pytest.raises(TypeError) as excinfo:
        struct.add("name", 32, Symbol.Visibility.PUBLIC)
    assert ("The type of a component of a StructureType must be a 'DataType' or 'DataTypeSymbol' but got 'int'") in str(excinfo.value)

    comps = struct.components
    assert comps["myint64"].datatype == scalartype
    assert comps["myfloat64"].datatype == scalartype_alt

    string_rep = "StructureType<>"
    assert string_rep == str(struct)

    assert struct.lookup("myint64").datatype == scalartype
    with pytest.raises(KeyError):
        struct.lookup("notthere")


def test_structure_type_create():
    scalartype = ScalarType(intrinsic=ScalarType.Intrinsic.INTEGER,
                            precision=ScalarType.Precision.DOUBLE)
    scalartype_alt = ScalarType(ScalarType.Intrinsic.REAL, 8)

    creation_list = [("myint64", scalartype, Symbol.Visibility.PUBLIC), ("myfloat64", scalartype_alt, Symbol.Visibility.PUBLIC)]
    bad_list = ["myint64"]

    with pytest.raises(TypeError) as excinfo:
        StructureType.create(bad_list)
    assert ("Each component must be specified using a 3-tuple "
            "of (name, type, visibility) but found a tuple with 7 members: myint64") in str(excinfo.value)

    struct = StructureType.create(creation_list)
#    string_rep = "StructureType<(myint64: Scalar<INTEGER, DOUBLE>), (myfloat64: Scalar<REAL, 8>)>"
    string_rep = "StructureType<>"
    assert string_rep == str(struct)

def test_pointer_type():
    scalartype = ScalarType(intrinsic=ScalarType.Intrinsic.INTEGER,
                            precision=ScalarType.Precision.DOUBLE)
    ptype = PointerType(scalartype)

    with pytest.raises(TypeError) as excinfo:
        PointerType("not type")
    assert ("Attempted to make a PointerType with datatype "
            "of type <class 'str'> instead of a DataType.") in str(excinfo.value)

    assert "PointerType<Scalar<INTEGER, DOUBLE>>" == str(ptype)

def test_array_type():
    scalartype = ScalarType(intrinsic=ScalarType.Intrinsic.INTEGER,
                            precision=ScalarType.Precision.DOUBLE)
    with pytest.raises(TypeError) as excinfo:
        ArrayType("null", [])
    assert ("Attempted to make an ArrayType with datatype "
            "of type <class 'str'> instead of a DataType.") in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        ArrayType(scalartype, "null")
    assert ("shape argument needs to be a nonempty list of int/Extents but found <class 'str'>") in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        ArrayType(scalartype, [])
    assert ("shape argument needs to be a nonempty list of int/Extents but found <class 'list'>") in str(excinfo.value)

    myshape = [32, ArrayType.Extent.DYNAMIC]
    badshape = [32, "string"]

    with pytest.raises(TypeError) as excinfo:
        ArrayType(scalartype, badshape)
    assert( "Each member of the shape argument needs to be an int "
            "or Extent but found <class 'str'>." in str(excinfo.value))

    mytype = ArrayType(scalartype, myshape)
    assert mytype.shape[0] == myshape[0]
    assert mytype.shape[1] == myshape[1]

    str_rep = "ArrayType<Scalar<INTEGER, DOUBLE>: [32, *]>"
    assert str(mytype) == str_rep

def test_inbuilt_types():
    # TODO
    assert str(FLOAT_TYPE) == "Scalar<REAL, SINGLE>"
    assert str(DOUBLE_TYPE) == "Scalar<REAL, DOUBLE>"
    assert str(INT_TYPE) == "Scalar<INTEGER, SINGLE>"
    assert str(LONG_LONG_TYPE) == "Scalar<INTEGER, DOUBLE>"
    assert str(INT32_TYPE) == "Scalar<INTEGER, 32>"
    assert str(INT64_TYPE) == "Scalar<INTEGER, 64>"
    assert str(INT8_TYPE) == "Scalar<INTEGER, 8>"
    assert str(BOOL_TYPE) == "Scalar<BOOLEAN, UNDEFINED>"

    assert str(PARTICLE_POSITION_TYPE) == "ArrayType<Scalar<REAL, DOUBLE>: [3]>"

def test_reset_part():
    type_mapping_str["part"].components["x"] = INT_TYPE
    reset_part_and_config()
    assert type_mapping_str["part"].components.get("x") is None

def test_reset_type_mapping_str():
    type_mapping_str["abcd"] = INT_TYPE
    reset_type_mapping_str()
    assert type_mapping_str.get("abcd") is None
