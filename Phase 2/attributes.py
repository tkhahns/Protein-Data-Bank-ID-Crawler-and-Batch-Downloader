from typing import TypeVarTuple, Generic

AttributeTypes = TypeVarTuple('AttributeTypes')

class Attributes(Generic[*AttributeTypes]):

    def __init__(self, attribute_pairs: list[tuple[str, str]],
                 primary_keys: list[str] = [], foreign_keys: dict[str, tuple[str, str]] = {}) -> None:
        
        self.attribute_names, self.attribute_types = tuple(zip(*attribute_pairs))
        if not set(primary_keys) <= set(self.attribute_names)\
            or not set(foreign_keys) <= set(self.attribute_names):
            raise ValueError("Primary keys and foreign keys need to be a subset of attributes")
        self.primary_keys = primary_keys
        self.foreign_keys = foreign_keys
        self.length = len(self.attribute_names)

    def __str__(self) -> str:
        attribute_pairs = zip(self.attribute_names, self.attribute_types)
        attribute_string =  ', '.join([pair[0] + ' ' + pair[1] for pair in attribute_pairs])
        primary_key_string = f"PRIMARY KEY ({', '.join(self.primary_keys)})"
        foreign_key_string = ', '.join([f"FOREIGN KEY ({key}) REFERENCES\
                                        {self.foreign_keys[key][0]} ({self.foreign_keys[key][1]})" for key in self.foreign_keys])
        string = attribute_string
        if self.primary_keys: string += ', ' + primary_key_string
        if self.foreign_keys: string += ', ' + foreign_key_string
        return f"({string})"
    
    def get_primary_keys(self):
        return self.primary_keys
    
    def tuple_to_dict(self, values: tuple[*AttributeTypes]):
        return {self.attribute_names[i]: values[i] for i in range(self.length)}
    
    def dict_to_tuple(self, values: dict) -> tuple[*AttributeTypes]:
        result = [None] * self.length
        for column in values:
            index = self.attribute_names.index(column)
            result[index] = values[column]
        return tuple(result)
    
    def match_all_columns(self, values: tuple[*AttributeTypes], delim = " AND ") -> str:
        if len(values) != self.length:
            raise ValueError("Number of values given does not match number of columns")
        return delim.join([self.attribute_names[i] + ' = ' + values[i] for i in range(self.length) if values[i] is not None])
    
    def match_columns(self, column_value_pairs: dict, delim = " AND ") -> str:
        if not set(column_value_pairs) <= set(self.attribute_names):
            raise ValueError("Argument contains columns not part of the table")
        return delim.join([column + ' = ' + column_value_pairs[column] for column in column_value_pairs])
    
    def match_primary_keys(self, primary_key_values, delim = " AND ") -> str:
        if len(primary_key_values) != len(self.primary_keys):
            raise ValueError("Number of values given does not match number of keys")
        return delim.join([self.primary_keys[i] + ' = ' + primary_key_values[i] for i in range(len(self.primary_keys))])