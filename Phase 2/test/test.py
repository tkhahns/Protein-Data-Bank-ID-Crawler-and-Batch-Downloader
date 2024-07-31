import unittest
import gemmi
import extract

class TestPdbMethods(unittest.TestCase):

    def setUp(self):
        st1 = gemmi.read_structure("")

    def test_sequence(self, control_sequence, sequence):
        self.assertEqual(control_sequence, sequence)

if __name__ == "__main__":
    unittest.main()