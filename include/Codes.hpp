#pragma once

#include "Code.hpp"

#include <cstddef>
#include <string>

class SteaneXCode : public Code {
public:
    SteaneXCode() : Code("./resources/codes/testCode.txt") {
        k = 3; // TODO not 1?
        n = 7;
        d = 3;
    }
};

class SteaneCode : public Code {
public:
    SteaneCode() : Code("./resources/codes/testCode.txt", "./resources/codes/testCode.txt") {
        k = 3; // TODO not 1?
        n = 7;
        d = 3;
    }
};

class HGPcode : public Code {
public:
    HGPcode(const std::string& inFile, const std::size_t size) : Code(inFile) { k = size; };
    HGPcode() : Code("./resources/codes/hgp_(4,7)-[[900,36,10]]_hz.txt") {
        k = 36;
        n = 900;
    }
    HGPcode(const std::string& hxIn, const std::string& hzIn, const std::size_t size) : Code(hxIn, hzIn) { k = size; };
};

// Surface codes

class ToricXCode8 : public Code {
public:
    ToricXCode8() : Code("./resources/codes/toric_(nan,nan)-[[8,2,2]]_hx.txt") {
        k = 2;
        n = 8;
        d = 2;
    }
};

class ToricCode8 : public Code {
public:
    ToricCode8() : Code("./resources/codes/toricCodes/toric_(nan,nan)-[[8,2,2]]_hz.txt") {
        k = 2;
        n = 8;
        d = 2;
    }
};

class ToricXCode18 : public Code {
public:
    ToricXCode18() : Code("./resources/codes/toric_(nan,nan)-[[18,2,3]]_hx.txt") {
        k = 2;
        n = 18;
        d = 3;
    }
};

class ToricCode18 : public Code {
public:
    ToricCode18() : Code("./resources/codes/toricCodes/toric_(nan,nan)-[[18,2,3]]_hz.txt") {
        k = 2;
        n = 18;
        d = 3;
    }
};

class ToricXCode32 : public Code {
public:
    ToricXCode32() : Code("./resources/codes/toric_(nan,nan)-[[32,2,4]]_hx.txt") {
        k = 2;
        n = 32;
        d = 4;
    }
};

class ToricCode32 : public Code {
public:
    ToricCode32() : Code("./resources/codes/toricCodes/toric_(nan,nan)-[[32,2,4]]_hz.txt") {
        k = 2;
        n = 32;
        d = 4;
    }
};

class ToricCode50 : public Code {
public:
    ToricCode50() : Code("./resources/codes/toricCodes/toric_(nan,nan)-[[50,2,5]]_hz.txt") {
        k = 2;
        n = 50;
        d = 5;
    }
};

class ToricCode72 : public Code {
public:
    ToricCode72() : Code("./resources/codes/toricCodes/toric_(nan,nan)-[[72,2,6]]_hz.txt") {
        k = 2;
        n = 72;
        d = 6;
    }
};

class ToricCode98 : public Code {
public:
    ToricCode98() : Code("./resources/codes/toricCodes/toric_(nan,nan)-[[98,2,7]]_hz.txt") {
        k = 2;
        n = 98;
        d = 7;
    }
};

// QLDPC codes from https://arxiv.org/abs/2308.07915

class BivarBikeCode72 : public Code {
public:
    BivarBikeCode72() : Code("./resources/codes/bbCodes/bb_(,)_[[72,12,6]]_hx.txt",
                             "./resources/codes/bbCodes/bb_(,)_[[72,12,6]]_hz.txt") {
        k = 12;
        n = 72;
        d = 6;
    }
};

class BivarBikeCode90 : public Code {
public:
    BivarBikeCode90() : Code("./resources/codes/bbCodes/bb_(,)_[[90,8,10]]_hx.txt",
                             "./resources/codes/bbCodes/bb_(,)_[[90,8,10]]_hz.txt") {
        k = 8;
        n = 90;
        d = 10;
    }
};

class BivarBikeCode144 : public Code {
public:
    BivarBikeCode144() : Code("./resources/codes/bbCodes/bb_(,)_[[144,12,12]]_hx.txt",
                              "./resources/codes/bbCodes/bb_(,)_[[144,12,12]]_hz.txt") {
        k = 12;
        n = 144;
        d = 12;
    }
};

class BivarBikeCode288 : public Code {
public:
    BivarBikeCode288() : Code("./resources/codes/bbCodes/bb_(,)_[[288,12,18]]_hx.txt",
                              "./resources/codes/bbCodes/bb_(,)_[[288,12,18]]_hz.txt") {
        k = 12;
        n = 288;
        d = 18;
    }
};
