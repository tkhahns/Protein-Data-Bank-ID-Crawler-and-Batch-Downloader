import gemmi 
from gemmi import cif, EntityType, PolymerType

from extract import ComplexType

class MockEntity:
    def __init__(self, entity_type: EntityType, polymer_type: PolymerType = None, subchains: list[str] = []):
        self.ent = gemmi.Entity("mock")
        self.ent.entity_type = entity_type 
        if polymer_type:
            self.ent.polymer_type = polymer_type
        self.ent.subchains = subchains
    
    def get_entity(self) -> gemmi.Entity:
        return self.ent
        
class MockComplexType:
    def __init__(self, complex_type: ComplexType):
        self.mock_struct = gemmi.Structure()
        self.complex_type = complex_type

    def create_mock_entities(self) -> gemmi.EntityList:
        entity_args = {
            ComplexType.Other: [
                (EntityType.Unknown,)
            ],
            ComplexType.NucleicAcid : [
                (EntityType.Polymer, PolymerType.Dna), 
                (EntityType.Polymer, PolymerType.Rna), 
                (EntityType.Polymer, PolymerType.DnaRnaHybrid),
                (EntityType.Polymer, PolymerType.Pna)
            ],
            ComplexType.Saccharide: [
                (EntityType.Branched,),
                (EntityType.Polymer, PolymerType.SaccharideD),
                (EntityType.Polymer, PolymerType.SaccharideL)
            ],
            ComplexType.SingleProtein: [
                (EntityType.Polymer, PolymerType.PeptideL, ["A"])
            ],
            ComplexType.Proteinmer: [
                (EntityType.Polymer, PolymerType.PeptideL, ["A", "B"])
            ], 
            ComplexType.ComplexProtein: [
                (EntityType.Polymer, PolymerType.PeptideD),
                (EntityType.Polymer, PolymerType.PeptideL)
            ],
            ComplexType.ProteinNA: [
                (EntityType.Polymer, PolymerType.PeptideD),
                (EntityType.Polymer, PolymerType.Dna)
            ],
            ComplexType.ProteinSaccharide: [
                (EntityType.Polymer, PolymerType.PeptideD),
                (EntityType.Polymer, PolymerType.SaccharideD)
            ],
            ComplexType.SaccharideNA: [
                (EntityType.Polymer, PolymerType.SaccharideD),
                (EntityType.Polymer, PolymerType.Dna)
            ], 
            ComplexType.ProteinSaccharideNA: [
                (EntityType.Polymer, PolymerType.PeptideD),
                (EntityType.Polymer, PolymerType.SaccharideD),
                (EntityType.Polymer, PolymerType.Dna)
            ]    
        }
        mock_entities = gemmi.EntityList()

        for args in entity_args[self.complex_type]:
            mock_entities.append(MockEntity(*args).get_entity())
            
        return mock_entities

    def get_struct(self) -> gemmi.Structure:
        self.mock_struct.entities = self.create_mock_entities()

        return self.mock_struct

class MockStrand:
    def __init__(self, strand_name: str, sense: int):
        self.strand = gemmi.Sheet.Strand()
        self.strand.name = strand_name
        self.strand.sense = sense

    def get_strand(self) -> gemmi.Sheet.Strand:
        return self.strand

class MockSheet:
    def __init__(self, sheet_name: str, strands: gemmi.Sheet.StrandList):
        self.sheet = gemmi.Sheet("mock")
        self.sheet.name = sheet_name
        self.sheet.strands = strands

    def get_sheet(self) -> gemmi.Sheet:
        return self.sheet

class MockSheetTable:
    def __init__(self):
        self.struct = gemmi.Structure()
        self.doc = cif.Document()
        self.struct.info["_entry.id"] = "0A0A"
        self.mock_sheets = self.create_mock_sheets()

    def create_mock_sheets(self) -> gemmi.SheetList:
        return gemmi.SheetList([
            MockSheet("A", gemmi.Sheet.StrandList([
                MockStrand("1", 1).get_strand(), 
                MockStrand("2", 1).get_strand(), 
                MockStrand("3", -1).get_strand()
            ])).get_sheet(),
            MockSheet("B", gemmi.Sheet.StrandList([
                MockStrand("1", 1).get_strand(), 
                MockStrand("2", 0).get_strand(), 
                MockStrand("3", -1).get_strand()
            ])).get_sheet()
        ])

    def get_doc(self):
        self.struct.sheets = self.mock_sheets
        return self.struct, self.doc
