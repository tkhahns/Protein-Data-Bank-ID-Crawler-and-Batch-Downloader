"""
This script contains unit tests for testing methods in attributes.py.
Make sure to run from the Phase 2 directory for the correct relative paths.

To run a specific test module, use the command "pytest test/unit/test_something.py".
To run all tests in the test directory, use the command "pytest test/".
Output verbosity can be adjusted by using the relevant flags in the command (e.g. -q, -v, -vv).
"""
import pytest 
from unittest.mock import MagicMock
from attributes import Attributes 

def test_attributes_initialisation():
    test_attributes_pairs = [("id", "VARCHAR"), ("a", "FLOAT")]
    test_primary_keys = ["id"]
    test_foreign_keys = {"id": ("main", "id")}

    test_attributes = Attributes(test_attributes_pairs, test_primary_keys, test_foreign_keys)
    assert test_attributes.attribute_names == ("id", "a")
    assert test_attributes.attribute_types == ("VARCHAR", "FLOAT")
    assert test_attributes.primary_keys == test_primary_keys
    assert test_attributes.foreign_keys == test_foreign_keys
    assert test_attributes.length == 2

def test_attributes_initialisation_default_keys():
    test_attributes_pairs = [("id", "VARCHAR"), ("a", "FLOAT")]

    test_attributes = Attributes(test_attributes_pairs)
    assert test_attributes.attribute_names == ("id", "a")
    assert test_attributes.attribute_types == ("VARCHAR", "FLOAT")
    assert test_attributes.primary_keys == []
    assert test_attributes.foreign_keys == {}
    assert test_attributes.length == 2

def test_table_initialisation_missing_arguments():
    """
    Test that a TypeError is raised when an Attributes object 
    is initialised with no arguments.
    """
    with pytest.raises(TypeError):
        Attributes()

def test_table_initialisation_invalid_attribute_pairs():
    """
    Test that a ValueError is raised when an Attributes object 
    is initialised with an invalid attributes_pairs list.
    """
    test_primary_keys = ["id"]
    test_foreign_keys = {"id": ("main", "id")}
    with pytest.raises(ValueError):
        # zip function requires two values to unpack 
        Attributes([], test_primary_keys, test_foreign_keys)
        Attributes([("id", ), ("a", )], test_primary_keys, test_foreign_keys)

def test_attributes_initialisation_invalid_primary_keys():
    """
    Test that a ValueError is raised when the provided primary keys 
    is not a subset of attributes.
    """
    test_attributes_pairs = [("id", "VARCHAR"), ("a", "FLOAT")]
    test_primary_keys = ["invalid_key"]  # not in attribute_names
    test_foreign_keys = {"id": ("main", "id")}

    with pytest.raises(ValueError, match="Primary keys and foreign keys need to be a subset of attributes"):
        Attributes(test_attributes_pairs, test_primary_keys, test_foreign_keys)

def test_attributes_initialisation_invalid_foreign_keys():
    """
    Test that a ValueError is raised when the provided foreign keys 
    is not a subset of attributes.
    """
    test_attributes_pairs = [("id", "VARCHAR"), ("a", "FLOAT")]
    test_primary_keys = ["id"]  
    test_foreign_keys = {"invalid_key": ("main", "id")}

    with pytest.raises(ValueError, match="Primary keys and foreign keys need to be a subset of attributes"):
        Attributes(test_attributes_pairs, test_primary_keys, test_foreign_keys)


def test_attributes_string(test_attributes):
    expected = "(id VARCHAR, a FLOAT, PRIMARY KEY (id, a), FOREIGN KEY (id) REFERENCES\
                                        main (id))"
    result = str(test_attributes)

    assert expected == result

def test_get_primary_keys(test_attributes):
    expected = ["id", "a"]
    result = test_attributes.get_primary_keys()

    assert expected == result

def test_tuple_to_dict(test_attributes):
    test_values = ("VARCHAR", "FLOAT")
    expected = {"id": "VARCHAR", "a": "FLOAT"}
    result = test_attributes.tuple_to_dict(test_values)

    assert expected == result

def test_tuple_to_dict_invalid_values(test_attributes):
    """
    Test that an IndexError is raised when the number of 
    values provided is fewer than the number of attributes. 
    """
    test_values = ("VARCHAR", )
    with pytest.raises(IndexError):
        test_attributes.tuple_to_dict(test_values)

def test_dict_to_tuple(test_attributes):
    test_values = {"id": "VARCHAR", "a": "FLOAT"}
    expected = ("VARCHAR", "FLOAT")
    result = test_attributes.dict_to_tuple(test_values)

    assert expected == result

def test_dict_to_tuple_invalid_values(test_attributes):
    test_values = {"id": "VARCHAR"}
    expected = ("VARCHAR", None)
    result = test_attributes.dict_to_tuple(test_values)

    assert expected == result

def test_dict_to_tuple_values_not_in_attributes(test_attributes):
    """
    Test that a ValueError is raised when a given key is
    not present in the attribute names.
    """
    test_values = {"id": "VARCHAR", "invalid_key": "FLOAT"}

    with pytest.raises(ValueError):
        test_attributes.dict_to_tuple(test_values)

def test_match_all_columns(test_attributes):
    test_values = ("'1A00'", "1.0")
    expected = "id = '1A00' AND a = 1.0"
    result = test_attributes.match_all_columns(test_values)

    assert expected == result

def test_match_all_columns_invalid_values(test_attributes):
    """
    Test that a ValueError is raised when the number of values 
    given does not match the number of attributes.
    """
    test_values = ("'1A00'", )
    with pytest.raises(ValueError, match="Number of values given does not match number of columns"):
        test_attributes.match_all_columns(test_values)

def test_match_columns(test_attributes):
    test_column_value_pairs = {"id": "'1A00'", "a": "1.0"}
    expected = "id = '1A00' AND a = 1.0"
    result = test_attributes.match_columns(test_column_value_pairs)

    assert expected == result

def test_match_columns_invalid_argument(test_attributes):
    """
    Test that a ValueError is raised when the given dictionary 
    contains columns not part of the attributes.
    """
    test_column_value_pairs = {"id": "'1A00'", "invalid_column": "1.0"}
    with pytest.raises(ValueError, match="Argument contains columns not part of the table"):
        test_attributes.match_columns(test_column_value_pairs)

def test_match_primary_keys(test_attributes):
    test_primary_key_values = ["'1A00'", '1.0']
    expected = "id = '1A00' AND a = 1.0"
    result = test_attributes.match_primary_keys(test_primary_key_values)

    assert expected == result

def test_match_primary_keys_invalid_values(test_attributes):
    """
    Test that a ValueError is raised when the number of values 
    given does not match the number of primary keys.
    """
    test_primary_key_values = ["'1A00'"]
    with pytest.raises(ValueError, match="Number of values given does not match number of keys"):
        test_attributes.match_primary_keys(test_primary_key_values)
        