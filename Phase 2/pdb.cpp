#include "pdb.hpp"
#include <sstream>
#include <iostream>
#include <gemmi/cif.hpp>

namespace cif = gemmi::cif;

namespace pdb {
    /*
    char convert_amino_3to1(string amino) {
        if (auto search = amino3to1.find(amino); search != amino3to1.end())
            return amino3to1[amino];
        return '?';
    }

    string convert_amino_1to3(char amino) {
        if (auto search = amino1to3.find(amino); search != amino1to3.end())
            return amino1to3[amino];
        return "???";
    }

    string convert_protein_3to1(vector<string> aminos) {
        stringstream converted_string;
        for (string amino : aminos)
            converted_string << convert_amino_3to1(amino);
        return converted_string.str();
    }

    protein::protein(string id, string name) {
        this->id = id;
        this->name = name;
    }
    */

    int main() {
        cif::Document doc = cif::read_file("./1bwh.cif");
        for (cif::Block& block : doc.blocks) {
            cout << *(block.find_value("_entry.id")) << endl;
        }
    }
}