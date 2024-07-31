from typing import TypeVarTuple, Callable, Generic
from attributes import Attributes
import gemmi
from gemmi import cif

AttributeTypes = TypeVarTuple('AttributeTypes')

class Table(Generic[*AttributeTypes]):
    def __init__(self, name: str, attributes: Attributes[*AttributeTypes],
                 extractor: Callable[[gemmi.Structure, cif.Document], list[tuple[*AttributeTypes]]]):
        self.name = name
        self.attributes = attributes
        self.extractor = extractor

    def attributes_string(self) -> str:
        return f"({','.join(self.attributes)})"

    def create_table(self) -> str:
        return f"CREATE TABLE IF NOT EXISTS {self.name} {str(self.attributes)}"
    
    def retrieve(self, columns=("*",)) -> str:
        return f"SELECT {', '.join(columns)} FROM {self.name}"
    
    def extract_data(self, struct: gemmi.Structure, doc: cif.Document) -> list[Attributes]:
        return self.extractor(struct, doc)
    
    def insert_row(self, data: Attributes):
        args = ', '.join(['?' for i in range(len(data))])
        return f"INSERT INTO {self.name} VALUES({args})"
    
    def update_row(self, data: dict, primary_key_values):
        return f"UPDATE {self.name} SET {self.attributes.match_columns(data), ', '}\
            WHERE {self.attributes.match_primary_keys(primary_key_values)}"
