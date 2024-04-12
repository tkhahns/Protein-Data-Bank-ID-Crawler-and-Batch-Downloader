#include "pdb.hpp"
#include "extract.hpp"
#include <gemmi/pdb.hpp>
#include <gemmi/gz.hpp>

namespace pdb {
    void print_relevant_info(string& path) {
        gemmi::Structure protein = gemmi::read_pdb(gemmi::MaybeGzipped(path));
        cout << protein.name << endl;
        cout << protein.models[0].name << endl;
        for (gemmi::Chain chain : protein.models[0].chains)
            cout << chain.name << endl;
        
    }

    string record_search(ifstream& fin, regex& record_type) {
        string line_in;
        while (getline(fin, line_in)) {
            if (regex_match(line_in, record_type))
                return line_in;
        }
        return "";
    }

    string encode_chain(string& to_encode, stringstream* in = new stringstream()) {
        string protein_chain = to_encode.substr(19);
        protein_chain.erase(find_if(protein_chain.rbegin(), protein_chain.rend(),
                [](char ch) {return !isspace(ch);}).base());
        int chain_length = (protein_chain.size() + 1)/4;
        for (int i = 0; i < chain_length; ++i) {
            
        }
    }
}