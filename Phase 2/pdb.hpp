#ifndef PDB_HPP
#define PDB_HPP

#include <map>
#include <vector>

using namespace std;

namespace pdb {
    static map<string, char> amino3to1 {
        {"ala", 'a'},
        {"arg", 'r'},
        {"asn", 'n'},
        {"asp", 'd'},
        {"cys", 'c'},
        {"gln", 'q'},
        {"glu", 'e'},
        {"gly", 'g'},
        {"his", 'h'},
        {"ile", 'i'},
        {"leu", 'l'},
        {"lys", 'k'},
        {"met", 'm'},
        {"phe", 'f'},
        {"pro", 'p'},
        {"ser", 's'},
        {"thr", 't'},
        {"trp", 'w'},
        {"tyr", 'y'},
        {"val", 'v'}
    };

    static map<char, string> amino1to3 {
        {'a', "ala"},
        {'r', "arg"},
        {'n', "asn"},
        {'d', "asp"},
        {'c', "cys"},
        {'q', "gln"},
        {'e', "glu"},
        {'g', "gly"},
        {'h', "his"},
        {'i', "ile"},
        {'l', "leu"},
        {'k', "lys"},
        {'m', "met"},
        {'f', "phe"},
        {'p', "pro"},
        {'s', "ser"},
        {'t', "thr"},
        {'w', "trp"},
        {'y', "tyr"},
        {'v', "val"}
    };

    char convert_amino_3to1(string);
    string convert_amino_1to3(char);

    enum class ChainType {
        A, B, C
    };

    class protein {
        private:
            string id;
            string name;
            int num_chains;
            vector<ChainType> chains;
            vector<string> chain_structures;
            vector<string> helix_structures;
            int protein_atoms;
            int nucleic_acid_atoms;
            int length_a;
            int length_b;
            int length_c;
            int alpha;
            int beta;
            int gamma;
            string space_group;
            int z_value;

        public:
            protein(string, string);
            void add_chain_structures(string, vector<string>);
            void add_crystal_structure(string);
    };
}

#endif
