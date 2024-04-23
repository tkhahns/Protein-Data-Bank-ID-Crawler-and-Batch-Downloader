from enum import Enum

main_table = "main"
entity_table = "entities"
subchain_table = "subchains"
helix_table = "helices"

ComplexType = Enum('ComplexType', ['Other', 'SingleProtein',
                 'NucleicAcid', 'ProteinNA', 'Saccharide',
                 'ProteinSaccharide', 'SaccharideNA', 'ProteinSaccharideNA',
                 'Proteinmer', 'ComplexProtein'], start=0)