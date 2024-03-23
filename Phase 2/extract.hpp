#ifndef EXTRACT_HPP
#define EXTRACT_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <regex>
#include <array>

using namespace std;

namespace pdb {
    static array<regex, 8> data = {
        regex("HEADER .*[a-z\d]{4}"),
        regex("COMPND   2 MOLECULE: .*;"),
        regex("COMPND   3 CHAIN: .*;"),
        regex("REMARK   3   PROTEIN ATOMS            : .*"),
        regex("REMARK   3   NUCELIC ACID ATOMS       : .*"),
        regex("SEQRES  [ 0-9]{2} [A-Z] [0-9]{4}  .*"),
        regex("HELIX   [ 0-9]{2} .*"),
        regex("CRYST1")
    };
    string record_search(ifstream, regex);
}

#endif
