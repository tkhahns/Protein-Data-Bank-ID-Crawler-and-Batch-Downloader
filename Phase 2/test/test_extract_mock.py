import pytest
from unittest.mock import patch
import gemmi 

import extract 
from extract import ComplexType
import mocks


#@pytest.mark.skip(reason=None)
@pytest.mark.parametrize("complex_type", ComplexType)
def test_get_complex_type(complex_type):
    result = extract.get_complex_type(mocks.MockComplexType(complex_type).get_struct())
    expected = complex_type 
    assert(result == expected)

@patch('polymer_sequence.PolymerSequence', autospec = True)
def test_sheet_table(mock_polymer_sequence):
    mock_sheet_table = mocks.MockSheetTable()
    mock_struct, mock_doc = mock_sheet_table.get_doc()

    expected = [
        ("0A0A", "A", 3, "PPA"), 
        ("0A0A", "B", 3, "PA")
    ]
    result = extract.insert_into_sheet_table(mock_struct, mock_doc, mock_polymer_sequence)

    assert (result == expected)
        
