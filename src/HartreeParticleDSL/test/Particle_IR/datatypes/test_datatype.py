import pytest

from HartreeParticleDSL.Particle_IR.datatypes.datatype import *
from HartreeParticleDSL.HartreeParticleDSLExceptions import IRGenerationError

def test_notype():
    a = NoType()
    assert str(a) == "NoType"

def test_scalar_type():

    with pytest.raises(TypeError) as excinfo:
        ScalarType(intrinsic="bad input", precision=32)
    assert ("Expected 'intrinsic' to be of type ScalarType.Intrinsic "
            "but found <class 'str'>." in str(excinfo.value))
    with pytest.raises(TypeError) as excinfo:
        ScalarType(intrinsic=ScalarType.Intrinsic.INTEGER, precision="bad input")
    assert ("Expected 'precision' to be of type ScalarType.Precision or int "
            "but found <class 'str'>.") in str(excinfo.value)

    scalartype = ScalarType(intrinsic=ScalarType.Intrinsic.INTEGER,
                            precision=ScalarType.Precision.DOUBLE)
    scalartype_alt = ScalarType(ScalarType.Intrinsic.FLOAT, 8)

    assert scalartype.intrinsic == ScalarType.Intrinsic.INTEGER
    assert scalartype.precision == ScalarType.Precision.DOUBLE

    assert str(scalartype) == "Scalar<INTEGER, DOUBLE>"
    assert str(scalartype_alt) == "Scalar<FLOAT, 8>"

def test_structure_type():
    struct = StructureType()

    scalartype = ScalarType(intrinsic=ScalarType.Intrinsic.INTEGER,
                            precision=ScalarType.Precision.DOUBLE)
    scalartype_alt = ScalarType(ScalarType.Intrinsic.FLOAT, 8)

    struct.add("myint64", scalartype)
    struct.add("myfloat64", scalartype_alt)

    with pytest.raises(TypeError) as excinfo:
        struct.add(32, scalartype)
    assert ("Expected name argument to be a string but got <class 'int'>.") in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        struct.add("name", 32)
    assert ("Expected typ argument to be a DataType but got <class 'int'>.") in str(excinfo.value)
    with pytest.raises(IRGenerationError) as excinfo:
        struct.add("myint64", scalartype)
    assert ("Names in a StructureType must be unique but provided duplicate myint64.") in str(excinfo.value)

    comps = struct.components
    assert comps["myint64"] == scalartype
    assert comps["myfloat64"] == scalartype_alt

    string_rep = "StructureType<(myint64: Scalar<INTEGER, DOUBLE>), (myfloat64: Scalar<FLOAT, 8>)>"
    assert string_rep == str(struct)

    assert struct.lookup("myint64") == scalartype
    assert struct.lookup("notthere") is None


def test_structure_type_create():
    scalartype = ScalarType(intrinsic=ScalarType.Intrinsic.INTEGER,
                            precision=ScalarType.Precision.DOUBLE)
    scalartype_alt = ScalarType(ScalarType.Intrinsic.FLOAT, 8)

    creation_list = [("myint64", scalartype), ("myfloat64", scalartype_alt)]
    bad_list = ["myint64"]

    with pytest.raises(TypeError) as excinfo:
        StructureType.create(bad_list)
    assert ("Each component must be specified using a 2-tuple "
            "of (name, type) but found: myint64") in str(excinfo.value)

    struct = StructureType.create(creation_list)
    string_rep = "StructureType<(myint64: Scalar<INTEGER, DOUBLE>), (myfloat64: Scalar<FLOAT, 8>)>"
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
    assert str(FLOAT_TYPE) == "Scalar<FLOAT, SINGLE>"
    assert str(DOUBLE_TYPE) == "Scalar<FLOAT, DOUBLE>"
    assert str(INT_TYPE) == "Scalar<INTEGER, SINGLE>"
    assert str(LONG_LONG_TYPE) == "Scalar<INTEGER, DOUBLE>"
    assert str(INT32_TYPE) == "Scalar<INTEGER, 32>"
    assert str(INT64_TYPE) == "Scalar<INTEGER, 64>"
    assert str(INT8_TYPE) == "Scalar<INTEGER, 8>"
    assert str(BOOL_TYPE) == "Scalar<BOOLEAN, UNDEFINED>"

    assert str(PARTICLE_POSITION_TYPE) == "ArrayType<Scalar<FLOAT, DOUBLE>: [3]>"

def test_reset_part():
    type_mapping_str["part"].components["x"] = INT_TYPE
    reset_part_and_config()
    assert type_mapping_str["part"].components.get("x") is None

def test_reset_type_mapping_str():
    type_mapping_str["abcd"] = INT_TYPE
    reset_type_mapping_str()
    assert type_mapping_str.get("abcd") is None
