/////////////////////////////////////////////////////////////////////////////////
// Program name       : homo.cpp
//
// Version            : 2.0
//
// Author             : Lars S Jermiin
//
// Institutions       : Australian National University
//                      Research School of Biology
//                      Acton, ACT 2601, Australia
//
//                      Univerity College Dublin
//                      School of Biology & Environmental Science
//                      Belfield, Dublin 4, Ireland
//
// Emails             : lars.jermiin [at] anu.edu.au
//                      lars.jermiin [at] ucd.ie
//
// URL                : https://github.com/lsjermiin/Homo2.0
//
// Date begun         : 14 April, 2019
//
// Date modified      : 19 October, 2019
//
// Copyright          : Copyright © 2019 Lars Sommer Jermiin.
//                      All rights reserved.
//
// Responsibility     : The copyright holder takes no legal responsibility for
//                      the correctness of results obtained using this program.
//
// Summary            : Homo is designed to conduct the matched-pairs test of
//                      symmetry for alignments of:
//
//                        1   Nucleotides (A|C|G|T) (4 states)
//                        2   Nucleotides recoded CTR = (C|T|AG) (3 states)
//                        3   Nucleotides recoded AGY = (A|G|CT) (3 states)
//                        4   Nucleotides recoded ATS = (A|T|CG) (3 states)
//                        5   Nucleotides recoded CGW = (C|G|AT) (3 states)
//                        6   Nucleotides recoded ACK = (A|C|GT) (3 states)
//                        7   Nucleotides recoded GTM = (G|T|AC) (3 states)
//                        8   Nucleotides recoded KM = (GT|AC)   (2 states)
//                        9   Nucleotides recoded RY = (AG|CT)   (2 states)
//                       10   Nucleotides recoded SW = (GC|AT)   (2 states)
//                       11   Nucleotides recoded AB = (A|CGT)   (2 states)
//                       12   Nucleotides recoded CD = (C|AGT)   (2 states)
//                       13   Nucleotides recoded GH = (G|ACT)   (2 states)
//                       14   Nucleotides recoded TV = (T|ACG)   (2 states)
//                       15   Di-nucleotides (AA|AC|...|TG|TT)  (16 states)
//                       16   Di-nucleotides 1st position (A|C|G|T) (4 states)
//                       17   Di-nucleotides 2nd position (A|C|G|T) (4 states)
//                       18   Codons (AAA|AAC|...|TTG|TTT) (64 states)
//                       19   Codons 1st + 2nd positions (AA|AC|...|TG|TT)  (16 states)
//                       20   Codons 1st + 3rd positions (AA|AC|...|TG|TT)  (16 states)
//                       21   Codons 2nd + 3rd positions (AA|AC|...|TG|TT)  (16 states)
//                       22   Codons 1st + 2nd positions (A|C|G|T) (4 states)
//                       23   Codons 1st + 3rd positions (A|C|G|T) (4 states)
//                       24   Codons 2nd + 3rd positions (A|C|G|T) (4 states)
//                       25   Codons 1st position (A|C|G|T) (4 states)
//                       26   Codons 2nd position (A|C|G|T) (4 states)
//                       27   Codons 3rd position (A|C|G|T) (4 states)
//                       28   10-state genotype data (A|C|G|T|K|M|R|Y|S|W)  (10 states)
//                       29   14-state genotype data (A|C|G|T|K|M|R|Y|S|W|B|D|H|V) (14 states)
//                       30   Amino acids (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)  (20 states)
//                       31   Recoded amino acids Dayhoff-6 = (AGPST|DENQ|HKR|MIVL|WFY|C)   (6 states)
//
//                      Homo also computes four compositional distances for each
//                      type of data:
//                        *   Bowkers's distance metric
//                        *   Stewart's distance metric
//                        *   Euclidean distance metric
//                      The input is the off-diagonal elements of a divergence
//                      matrix or the marginal sum of the off-diagonal elements
//                      of this matrix.
//
//                      Sequences must be stored in the FASTA format.
//
//                      Characters are converted to integers to speed up the
//                      program.
//
// Nucleotides        : Alphabet: [A,C.G,T/U,-] = [0,1,2,3,4].
//
//                      Ambiguous characters (i.e., ?, N, B, D, H, K, M, R, S,
//                      V, W and Y) are treated as if they were alignment gaps
//                      (-) (i.e., as missing data).
//
//                      10-state genotypes : Alphabet: [A,C,G,K,M,R,S,T/U,W,Y,-] =
//                      [0,1,2,3,4,5,6,7,8,9,10].
//
//                      Ambiguous characters (i.e., ?, N, B, D, G and V) are
//                      treated as if they were alignment gaps (-) (i.e., as
//                      missing data).
//
//                      14-state genotypes : Alphabet: [A,C,G,T/U,K,M,R,S,W,Y,B,D,H,V,-] =
//                      [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14].
//
//                      Ambiguous characters (i.e., ? and N) are treated as if
//                      they were alignment gaps (-) (i.e., as missing data).
//
// Amino acids        : Alphabet: [A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,-] =
//                      [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
//
//                      Ambiguous characters (i.e., ?, X and Z) are treated as
//                      if they were alignment gaps (-) (i.e., as missing data).
//
// Manuscript         : Jermiin et al. (2019).
//                      Software for analysing compositional heterogeneity across
//                      aligned sequence data. Syst. Biol. (in prep.)
//
// References         : Bowker A. H. (1948). A test of symmetry in contingency
//                      tables. Journal of the American Statistical Association
//                      43: 572-574.
//
/////////////////////////////////////////////////////////////////////////////////

#include <cctype>
#include <cmath>
#include <limits>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#define SQR(a) ((a) * (a))

using namespace std;
using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

// The following variables are declared externally because they
// are needed in different functions
const unsigned TWO(2);          // for 2-state alphabet (recoded DNA)
const unsigned THREE(3);        // for 3-state alphabet (recoded DNA)
const unsigned FOUR(4);         // for 4-state alphabet (DNA)
const unsigned SIX(6);          // for 6-state alphabet (recoded amino acids)
const unsigned TEN(10);         // for 10-state alphabet (genotype data)
const unsigned FOURTEEN(14);    // for 14-state alphabet (genotype data)
const unsigned SIXTEEN(16);     // for 16-state alphabet (subset of codons)
const unsigned TWENTY(20);      // for 20-state alphabet (amino acids)
const unsigned SIXTYFOUR(64);   // for 64-state alphabet (codons)
const unsigned max_array(65);   // corresponding to amino acids & gaps
string sites;                  // string controlling whether a site is used or not
vector<string> taxon;           // 2D container for sequence names
vector<vector<int> > alignment; // 2D container for sequence data



// This function translates a string of characters into a vector of integers
vector<int> Translator(unsigned datatype, string seq) {
    int unit; // integer for singlet, duplet or triplet (codon)
    string duplet(""), triplet(""); // strings for dinucleotides and codons
    vector<int> seq_data;
   
    switch (datatype) {
        case 1: // Nucleotides (A|C|G|T)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 2: // Nucleotides (C|T|R)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(2); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'R': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 3: // Nucleotides (A|G|Y)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(2); break;
                    case 'U': seq_data.push_back(2); break;
                    case 'Y': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 4: // Nucleotides (A|T|S)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'C': seq_data.push_back(2); break;
                    case 'S': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 5: // Nucleotides (C|G|W)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(2); break;
                    case 'U': seq_data.push_back(2); break;
                    case 'W': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 6: // Nucleotides (A|C|K)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(2); break;
                    case 'U': seq_data.push_back(2); break;
                    case 'K': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 7: // Nucleotides (G|T|M)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'G': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(2); break;
                    case 'C': seq_data.push_back(2); break;
                    case 'M': seq_data.push_back(2); break;
                    default : seq_data.push_back(3); break; // In case of other characters
                }
            }
            break;
        case 8: // Nucleotides (K|M)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'G': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(0); break;
                    case 'U': seq_data.push_back(0); break;
                    case 'K': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'M': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 9: // Nucleotides (R|Y)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(0); break;
                    case 'R': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'Y': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 10: // Nucleotides (S|W)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(0); break;
                    case 'R': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'Y': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 11: // Nucleotides (A|B)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'B': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 12: // Nucleotides (C|D)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'C': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'D': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 13: // Nucleotides (G|H)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'G': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'T': seq_data.push_back(1); break;
                    case 'U': seq_data.push_back(1); break;
                    case 'H': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 14: // Nucleotides (T|V)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'T': seq_data.push_back(0); break;
                    case 'U': seq_data.push_back(0); break;
                    case 'A': seq_data.push_back(1); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(1); break;
                    case 'V': seq_data.push_back(1); break;
                    default : seq_data.push_back(2); break; // In case of other characters
                }
            }
            break;
        case 15: // Di-nucleotides (AA|AC|..|TG|TT)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 2) {
                for (string::size_type j = i; j != i + 2; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': duplet.push_back('0'); break;
                        case 'C': duplet.push_back('1'); break;
                        case 'G': duplet.push_back('2'); break;
                        case 'T': duplet.push_back('3'); break;
                        case 'U': duplet.push_back('3'); break;
                        default : duplet.push_back('4'); break;
                    }
                }
                unit = stoi(duplet);
                duplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 16: // Di-nucleotides 1st position (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 2) {
                switch (toupper(seq[i])) {
                    case 'A' : seq_data.push_back(0); break; // 1st A
                    case 'C' : seq_data.push_back(1); break; // 1st C
                    case 'G' : seq_data.push_back(2); break; // 1st G
                    case 'T' : seq_data.push_back(3); break; // 1st T
                    default  : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 17: // Di-nucleotides 2nd position (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 2) {
                switch (toupper(seq[i+1])) {
                    case 'A' : seq_data.push_back(0); break; // 2nd A
                    case 'C' : seq_data.push_back(1); break; // 2nd C
                    case 'G' : seq_data.push_back(2); break; // 2nd G
                    case 'T' : seq_data.push_back(3); break; // 2nd T
                    default  : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 18: // Codons (AAA|AAC|...|TTG|TTT)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 000: seq_data.push_back(0);  break; // AAA
                    case 001: seq_data.push_back(1);  break; // AAC
                    case 002: seq_data.push_back(2);  break; // AAG
                    case 003: seq_data.push_back(3);  break; // AAT
                    case 010: seq_data.push_back(4);  break; // ACA
                    case 011: seq_data.push_back(5);  break; // ACC
                    case 012: seq_data.push_back(6);  break; // ACG
                    case 013: seq_data.push_back(7);  break; // ACT
                    case 020: seq_data.push_back(8);  break; // AGA
                    case 021: seq_data.push_back(9);  break; // AGC
                    case 022: seq_data.push_back(10); break; // AGG
                    case 023: seq_data.push_back(11); break; // AGT
                    case 030: seq_data.push_back(12); break; // ATA
                    case 031: seq_data.push_back(13); break; // ATC
                    case 032: seq_data.push_back(14); break; // ATG
                    case 033: seq_data.push_back(15); break; // ATT
                    case 100: seq_data.push_back(16); break; // CAA
                    case 101: seq_data.push_back(17); break; // CAC
                    case 102: seq_data.push_back(18); break; // CAG
                    case 103: seq_data.push_back(19); break; // CAT
                    case 110: seq_data.push_back(20); break; // CCA
                    case 111: seq_data.push_back(21); break; // CCC
                    case 112: seq_data.push_back(22); break; // CCG
                    case 113: seq_data.push_back(23); break; // CCT
                    case 120: seq_data.push_back(24); break; // CGA
                    case 121: seq_data.push_back(25); break; // CGC
                    case 122: seq_data.push_back(26); break; // CGG
                    case 123: seq_data.push_back(27); break; // CGT
                    case 130: seq_data.push_back(28); break; // CTA
                    case 131: seq_data.push_back(29); break; // CTC
                    case 132: seq_data.push_back(30); break; // CTG
                    case 133: seq_data.push_back(31); break; // CTT
                    case 200: seq_data.push_back(32); break; // GAA
                    case 201: seq_data.push_back(33); break; // GAC
                    case 202: seq_data.push_back(34); break; // GAG
                    case 203: seq_data.push_back(35); break; // GAT
                    case 210: seq_data.push_back(36); break; // GCA
                    case 211: seq_data.push_back(37); break; // GCC
                    case 212: seq_data.push_back(38); break; // GCG
                    case 213: seq_data.push_back(39); break; // GCT
                    case 220: seq_data.push_back(40); break; // GGA
                    case 221: seq_data.push_back(41); break; // GGC
                    case 222: seq_data.push_back(42); break; // GGG
                    case 223: seq_data.push_back(43); break; // GGT
                    case 230: seq_data.push_back(44); break; // GTA
                    case 231: seq_data.push_back(45); break; // GTC
                    case 232: seq_data.push_back(46); break; // GTG
                    case 233: seq_data.push_back(47); break; // GTT
                    case 300: seq_data.push_back(48); break; // TAA
                    case 301: seq_data.push_back(49); break; // TAC
                    case 302: seq_data.push_back(50); break; // TAG
                    case 303: seq_data.push_back(51); break; // TAT
                    case 310: seq_data.push_back(52); break; // TCA
                    case 311: seq_data.push_back(53); break; // TCC
                    case 312: seq_data.push_back(54); break; // TCG
                    case 313: seq_data.push_back(55); break; // TCT
                    case 320: seq_data.push_back(56); break; // TGA
                    case 321: seq_data.push_back(57); break; // TGC
                    case 322: seq_data.push_back(58); break; // TGG
                    case 323: seq_data.push_back(59); break; // TGT
                    case 330: seq_data.push_back(60); break; // TTA
                    case 331: seq_data.push_back(61); break; // TTC
                    case 332: seq_data.push_back(62); break; // TTG
                    case 333: seq_data.push_back(63); break; // TTT
                    default:  seq_data.push_back(64); break; // In case of other characters
                }
            }
            break;
        case 19: // Codons 1st + 2nd positions (AA|AC|...|TG|TT)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(2,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 20: // Codons 1st + 3rd positions (AA|AC|...|TG|TT)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(1,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 21: // Codons 2nd + 3rd positions (AA|AC|...|TG|TT)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(0,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case 00: seq_data.push_back(0);  break; // AA
                    case 01: seq_data.push_back(1);  break; // AC
                    case 02: seq_data.push_back(2);  break; // AG
                    case 03: seq_data.push_back(3);  break; // AT
                    case 10: seq_data.push_back(4);  break; // CA
                    case 11: seq_data.push_back(5);  break; // CC
                    case 12: seq_data.push_back(6);  break; // CG
                    case 13: seq_data.push_back(7);  break; // CT
                    case 20: seq_data.push_back(8);  break; // GA
                    case 21: seq_data.push_back(9);  break; // GC
                    case 22: seq_data.push_back(10); break; // GG
                    case 23: seq_data.push_back(11); break; // GT
                    case 30: seq_data.push_back(12); break; // TA
                    case 31: seq_data.push_back(13); break; // TC
                    case 32: seq_data.push_back(14); break; // TG
                    case 33: seq_data.push_back(15); break; // TT
                    default: seq_data.push_back(16); break; // In case of other characters
                }
            }
            break;
        case 22: // Codons 1st + 2nd positions (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    if (j == i || j == i + 1) {
                        switch (toupper(seq[j])) {
                            case 'A': seq_data.push_back(0); break;
                            case 'C': seq_data.push_back(1); break;
                            case 'G': seq_data.push_back(2); break;
                            case 'T': seq_data.push_back(3); break;
                            case 'U': seq_data.push_back(3); break;
                            default : seq_data.push_back(4); break;
                        }
                    }
                }
            }
            break;
        case 23: // Codons 1st + 3rd positions (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    if (j == i || j == i + 2) {
                        switch (toupper(seq[j])) {
                            case 'A': seq_data.push_back(0); break;
                            case 'C': seq_data.push_back(1); break;
                            case 'G': seq_data.push_back(2); break;
                            case 'T': seq_data.push_back(3); break;
                            case 'U': seq_data.push_back(3); break;
                            default : seq_data.push_back(4); break;
                        }
                    }
                }
            }
            break;
        case 24: // Codons 2nd + 3rd positions (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    if (j == i + 1 || j == i + 2) {
                        switch (toupper(seq[j])) {
                            case 'A': seq_data.push_back(0); break;
                            case 'C': seq_data.push_back(1); break;
                            case 'G': seq_data.push_back(2); break;
                            case 'T': seq_data.push_back(3); break;
                            case 'U': seq_data.push_back(3); break;
                            default : seq_data.push_back(4); break;
                        }
                    }
                }
            }
            break;
        case 25: // Codons 1st position (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(2,1);
                triplet.erase(1,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case   0: seq_data.push_back(0); break;
                    case   1: seq_data.push_back(1); break;
                    case   2: seq_data.push_back(2); break;
                    case   3: seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 26: // Codons 2nd position (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(2,1);
                triplet.erase(0,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case   0: seq_data.push_back(0); break;
                    case   1: seq_data.push_back(1); break;
                    case   2: seq_data.push_back(2); break;
                    case   3: seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 27: // Codons 3rd position (A|C|G|T)
            // Check length
            for (string::size_type i = 0; i != seq.size(); i = i + 3) {
                for (string::size_type j = i; j != i + 3; j++) {
                    switch (toupper(seq[j])) {
                        case 'A': triplet.push_back('0'); break;
                        case 'C': triplet.push_back('1'); break;
                        case 'G': triplet.push_back('2'); break;
                        case 'T': triplet.push_back('3'); break;
                        case 'U': triplet.push_back('3'); break;
                        default : triplet.push_back('4'); break;
                    }
                }
                triplet.erase(1,1);
                triplet.erase(0,1);
                unit = stoi(triplet);
                triplet.clear();
                switch (unit) {
                    case   0: seq_data.push_back(0); break;
                    case   1: seq_data.push_back(1); break;
                    case   2: seq_data.push_back(2); break;
                    case   3: seq_data.push_back(3); break;
                    default : seq_data.push_back(4); break; // In case of other characters
                }
            }
            break;
        case 28: // 10-state genotype data
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    case 'K': seq_data.push_back(4); break;
                    case 'M': seq_data.push_back(5); break;
                    case 'R': seq_data.push_back(6); break;
                    case 'Y': seq_data.push_back(7); break;
                    case 'S': seq_data.push_back(8); break;
                    case 'W': seq_data.push_back(9); break;
                    default : seq_data.push_back(10);break; // In case of other characters
                }
            }
            break;
        case 29: // 14-state genotype data
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'G': seq_data.push_back(2); break;
                    case 'T': seq_data.push_back(3); break;
                    case 'U': seq_data.push_back(3); break;
                    case 'K': seq_data.push_back(4); break;
                    case 'M': seq_data.push_back(5); break;
                    case 'R': seq_data.push_back(6); break;
                    case 'Y': seq_data.push_back(7); break;
                    case 'S': seq_data.push_back(8); break;
                    case 'W': seq_data.push_back(9); break;
                    case 'B': seq_data.push_back(10);break;
                    case 'D': seq_data.push_back(11);break;
                    case 'H': seq_data.push_back(12);break;
                    case 'V': seq_data.push_back(13);break;
                    default : seq_data.push_back(14);break; // In case of other characters
                }
            }
            break;
        case 30: // amino acids (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'C': seq_data.push_back(1); break;
                    case 'D': seq_data.push_back(2); break;
                    case 'E': seq_data.push_back(3); break;
                    case 'F': seq_data.push_back(4); break;
                    case 'G': seq_data.push_back(5); break;
                    case 'H': seq_data.push_back(6); break;
                    case 'I': seq_data.push_back(7); break;
                    case 'K': seq_data.push_back(8); break;
                    case 'L': seq_data.push_back(9); break;
                    case 'M': seq_data.push_back(10);break;
                    case 'N': seq_data.push_back(11);break;
                    case 'P': seq_data.push_back(12);break;
                    case 'Q': seq_data.push_back(13);break;
                    case 'R': seq_data.push_back(14);break;
                    case 'S': seq_data.push_back(15);break;
                    case 'T': seq_data.push_back(16);break;
                    case 'V': seq_data.push_back(17);break;
                    case 'W': seq_data.push_back(18);break;
                    case 'Y': seq_data.push_back(19);break;
                    default : seq_data.push_back(20);break; // In case of other characters
                }
            }
            break;
        default: // Dayhoff-6 (AGPST|DENQ|HKR|MIVL|WFY|C)
            for (string::size_type i = 0; i != seq.size(); ++i) {
                switch (toupper(seq[i])) {
                    case 'A': seq_data.push_back(0); break;
                    case 'G': seq_data.push_back(0); break;
                    case 'P': seq_data.push_back(0); break;
                    case 'S': seq_data.push_back(0); break;
                    case 'T': seq_data.push_back(0); break;
                    case 'D': seq_data.push_back(1); break;
                    case 'E': seq_data.push_back(1); break;
                    case 'N': seq_data.push_back(1); break;
                    case 'Q': seq_data.push_back(1); break;
                    case 'H': seq_data.push_back(2); break;
                    case 'K': seq_data.push_back(2); break;
                    case 'R': seq_data.push_back(2); break;
                    case 'M': seq_data.push_back(3); break;
                    case 'I': seq_data.push_back(3); break;
                    case 'L': seq_data.push_back(3); break;
                    case 'V': seq_data.push_back(3); break;
                    case 'F': seq_data.push_back(4); break;
                    case 'W': seq_data.push_back(4); break;
                    case 'Y': seq_data.push_back(4); break;
                    case 'C': seq_data.push_back(5); break;
                    default : seq_data.push_back(6); break; // In case of other characters
                }
            }
            break;
    }
    return(seq_data);
}


// Function that reads input file and stores data in two 2D containers
void Read_Input(string inname, unsigned datatype){
    unsigned long alignment_length(0);
    unsigned long counter(0);
    string seq(""), str(""), tmp(""); // temporary string used to store input
    vector<int> sequence;     // temporary vector used to store input
    ifstream infile;
    
    infile.open(inname.c_str());
    if (!infile) {
        cerr << "Program aborted due to error: input file not found" << endl;
        exit(1);
    }
    while (getline(infile, str)) {
        if (!str.empty()) {
            // remove blank space in string
            tmp.clear();
            for (std::string::size_type i = 0; i != str.size(); ++i) {
                if (!isblank(str[i])) {
                    tmp.push_back(str[i]);
                }
            }
            if (tmp[0] == '>') {
                if (seq.size() > 0) {
                    if (datatype > 14 && datatype < 18) {
                        if (seq.size() % 2 != 0) {
                            std::cerr << "\nERROR: expected sequence of di-nucleotides" << "\n" << std::endl;
                            exit(1);
                        }
                    }
                    if (datatype > 17 && datatype < 28) {
                        if (seq.size() % 3 != 0) {
                            std::cerr << "\nERROR: expected sequence of codons" << "\n" << std::endl;
                            exit(1);
                        }
                    }
                    sequence = Translator(datatype, seq);
                    alignment.push_back(sequence); // stores sequence in vector
                    if (alignment_length == 0)
                        alignment_length = sequence.size();
                    sequence.clear();
                    seq.clear();
                }
                tmp.erase(tmp.begin()); // removes first character from name
                taxon.push_back(tmp); // stores sequence name in vector
            } else {
                seq += tmp;
            }
            str.clear();
        }
    }
    // Store last sequence in vector
    if (seq.size() > 0) {
        if (datatype > 14 && datatype < 18) {
            if (seq.size() % 2 != 0) {
                cerr << "\nProgram aborted due to error: expected sequence of di-nucleotides" << "\n" << endl;
                exit(1);
            }
        }
        if (datatype > 17 && datatype < 28) {
            if (seq.size() % 3 != 0) {
                cerr << "\nProgram aborted due to error: expected sequence of codons" << "\n" << endl;
                exit(1);
            }
        }
       sequence = Translator(datatype, seq);
        alignment.push_back(sequence);
    } else {
        cerr << "Program aborted due to error: last sequence empty" << "\n" << endl;
        exit(1);
    }
    //Check whether the sequence names are unique
    for (vector<string>::const_iterator iter1 = taxon.begin(); iter1 != taxon.end(); ++iter1) {
        for (vector<string>::const_iterator iter2 = iter1 + 1; iter2 != taxon.end(); ++iter2) {
            if (*iter1 == *iter2) {
                cerr << "Program aborted due to error: Sequence name not unique -- look for " << *iter1 << "\n" << endl;
                exit(1);
            }
        }
    }
    // Check whether the sequences have the same length
    for (vector<vector<int> >::const_iterator iter = alignment.begin()+1; iter != alignment.end(); ++iter) {
        ++counter;
        sequence = *iter;
        if (sequence.size() != alignment_length) {
            cerr << "Program aborted due to error: sequences 1 and " << counter << " differ in length!\n" << endl;
            exit(1);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// NOTE PERTAINING TO THE FOLLOWING FOUR FUNCTIONS                            //
//                                                                            //
// Source: gaussian_distribution_tail.c                                       //
// Source: chi-square_distribution_tail.c                                     //
// Author: Dick Horn (mathretprogr@gmail.com)                                 //
// Note: Used with permission from the author (Wednesday, 9 July 2014 4:30)   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


// This function returns the probability that a random variable with a standard
// Normal (Gaussian) distribution has a value greater than "x"
long double xGaussian_Distribution_Tail( long double x ) {
    long double sqrt2 = 0.7071067811865475244008443621048490L;
    return  0.5L * erfcl(sqrt2 * x );
}

// The number of degrees of freedom, nu, is an even integer, nu = 2*n.
// The argument x is chi^2 / 2.
static long double Sum_Poisson_Terms(long double x, int n) {
    int k;
    long double term;
    long double sum;
    
    term = 1.0L;
    sum = 1.0L;
    for (k = 1; k < n; k++) {
        term *= (x / k);
        sum += term;
    }
    return expl(-x)*sum;
}

static long double Sum_Over_Odd_Terms(long double x, int dof) {
    int k;
    int n;
    long double term;
    long double sum;
    long double sqrtx;
    long double twooverpi;
    
    twooverpi = 0.6366197723675813430755350534900574L;
    
    if (dof == 1) return 2.0L * xGaussian_Distribution_Tail( sqrtl(x) );
    n = (dof - 1) / 2;
    sqrtx = sqrtl(x);
    term = sqrtx;
    sum = sqrtx;
    for (k = 2; k <=n; k++) {
        term *= ( x / (k + k - 1) );
        sum += term;
    };
    return 2.0L * xGaussian_Distribution_Tail(sqrtx) + sqrtl(twooverpi) * expl(-x/2.0L)*sum;
}

// This function returns the probability that a random chi-squared-distributed
// variable has a value greater than "x"
long double xChi_Square_Distribution_Tail(long double x, int dof) {
    
    if (dof <= 0) return 0.0L;
    
    if (dof % 2 == 0)
        return Sum_Poisson_Terms(x/2.0L, dof/2);
    else
        return Sum_Over_Odd_Terms(x,dof);
}



//Main program that calls the various functions above
int main(int argc, char** argv){
    unsigned rows_columns(0);
    unsigned dataType(0);
    int df(0);
    unsigned long sum_dm(0), rejectB(0);
    unsigned long total, counter;
    unsigned long dm[max_array][max_array];            // 2D divergence matrix
    double d_efs(0.0); // Euclidean distance -- full symmetry
    double d_ems(0.0); // Euclidean distance -- marginal symmetry
    double d_cfs(0.0); // Compositional distance -- full symmetry
    double row_sum[max_array], col_sum[max_array];
    double min_p(1.0), family_Wise_Error_Rate(0.05);
    double max_dcfs(-numeric_limits<double>::max());
    double min_dcfs(numeric_limits<double>::max());
    long double BTS(0.0), pB(0.0);
    vector<int> sites;
    vector<double> row_of_double;
    vector<vector<double> > mat_dobs, mat_p, mat_dcfs, mat_defs, mat_dems;
    string nature_of_data, survey;
    string inName, outName1, outName2, outName3, outName4, outName5;
    ofstream outfile1, outfile2, outfile3, outfile4, outfile5;
    
    if(argc != 4) {
        cerr << "\nHomo v2.0 Copyright 2019, Lars Jermiin" << endl;
        cerr << " Contact: lars.jermiin [at] anu.edu.au / ucd.ie" << endl;
        cerr << "\nERROR -- use command: homo <infile> <b|f> <1|...|31>\n" << endl;
        cerr << "  infile   Fasta-formatted alignment" << endl;
        cerr << "     b|f   Brief or full report of results" << endl;
        cerr << "       1   Nucleotides; 4 states (A|C|G|T)" << endl;
        cerr << "       2   Nucleotides; 3 states (C|T|AG)" << endl;
        cerr << "       3   Nucleotides; 3 states (A|G|CT)" << endl;
        cerr << "       4   Nucleotides; 3 states (A|T|CG)" << endl;
        cerr << "       5   Nucleotides; 3 states (C|G|AT)" << endl;
        cerr << "       6   Nucleotides; 3 states (A|C|GT)" << endl;
        cerr << "       7   Nucleotides; 3 states (G|T|AC)" << endl;
        cerr << "       8   Nucleotides; 2 states (GT|AC)" << endl;
        cerr << "       9   Nucleotides; 2 states (AG|CT)" << endl;
        cerr << "      10   Nucleotides; 2 states (GC|AT)" << endl;
        cerr << "      11   Nucleotides; 2 states (A|CGT)" << endl;
        cerr << "      12   Nucleotides; 2 states (C|AGT)" << endl;
        cerr << "      13   Nucleotides; 2 states (G|ACT)" << endl;
        cerr << "      14   Nucleotides; 2 states (T|ACG)" << endl;
        cerr << "      15   Di-nucleotides (pos 12); 16 states (AA|AC|...|TG|TT)" << endl;
        cerr << "      16   Di-nucleotides (pos  1);  4 states (A|C|G|T)" << endl;
        cerr << "      17   Di-nucleotides (pos  2);  4 states (A|C|G|T)" << endl;
        cerr << "      18   Codons (pos 123); 64 states (AAA|AAC|...|TTG|TTT)" << endl;
        cerr << "      19   Codons (pos  12); 16 states (AA|AC|...|TG|TT)" << endl;
        cerr << "      20   Codons (pos  13); 16 states (AA|AC|...|TG|TT)" << endl;
        cerr << "      21   Codons (pos  23); 16 states (AA|AC|...|TG|TT)" << endl;
        cerr << "      22   Codons (pos  12);  4 states (A|C|G|T)" << endl;
        cerr << "      23   Codons (pos  13);  4 states (A|C|G|T)" << endl;
        cerr << "      24   Codons (pos  23);  4 states (A|C|G|T)" << endl;
        cerr << "      25   Codons (pos   1);  4 states (A|C|G|T)" << endl;
        cerr << "      26   Codons (pos   2);  4 states (A|C|G|T)" << endl;
        cerr << "      27   Codons (pos   3);  4 states (A|C|G|T)" << endl;
        cerr << "      28   Genotypes; 10 states (A|C|G|T|K|M|R|Y|S|W)" << endl;
        cerr << "      29   Genotypes; 14 states (A|C|G|T|K|M|R|Y|S|W|B|D|H|V)" << endl;
        cerr << "      30   Amino acids; 20 states (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)" << endl;
        cerr << "      31   Amino acids;  6 states (AGPST|DENQ|HKR|MIVL|WFY|C) [D6]" << endl;
        cerr << endl;
        exit(1);
    }
    inName = argv[1];
    survey = argv[2];
    nature_of_data = argv[3];
    dataType = stoi(nature_of_data);
    if (dataType < 1 || dataType > 31) {
        cerr << "\nPROGRAM ABORTED - incorrect choice of data: [1|...|31]\n" << endl;
        exit(1);
    }
    if (toupper(survey[0]) != 'F' && toupper(survey[0]) != 'B') {
        cerr << "\nPROGRAM ABORTED - incorrect choice of output: [b|f]\n" << endl;
        exit(1);
    }
    if (toupper(survey[0]) == 'F') {
        outName1.clear();
        for (string::size_type i = 0; i != inName.size() && inName[i] != '.'; ++i) {
            outName1 += inName[i];
        }
        outName5 = outName1 + "_dems.dis";
        outName4 = outName1 + "_defs.dis";
        outName3 = outName1 + "_dcfs.dis";
        outName2 = outName1 + "_Pvalues.csv";
        outName1 = outName1 + "_Summary.csv";
    }
    switch (dataType) {
        case  1: rows_columns = FOUR; break;
        case  2: rows_columns = THREE; break;
        case  3: rows_columns = THREE; break;
        case  4: rows_columns = THREE; break;
        case  5: rows_columns = THREE; break;
        case  6: rows_columns = THREE; break;
        case  7: rows_columns = THREE; break;
        case  8: rows_columns = TWO; break;
        case  9: rows_columns = TWO; break;
        case 10: rows_columns = TWO; break;
        case 11: rows_columns = TWO; break;
        case 12: rows_columns = TWO; break;
        case 13: rows_columns = TWO; break;
        case 14: rows_columns = TWO; break;
        case 15: rows_columns = SIXTEEN; break;
        case 16: rows_columns = FOUR; break;
        case 17: rows_columns = FOUR; break;
        case 18: rows_columns = SIXTYFOUR; break;
        case 19: rows_columns = SIXTEEN; break;
        case 20: rows_columns = SIXTEEN; break;
        case 21: rows_columns = SIXTEEN; break;
        case 22: rows_columns = FOUR; break;
        case 23: rows_columns = FOUR; break;
        case 24: rows_columns = FOUR; break;
        case 25: rows_columns = FOUR; break;
        case 26: rows_columns = FOUR; break;
        case 27: rows_columns = FOUR; break;
        case 28: rows_columns = TEN; break;
        case 29: rows_columns = FOURTEEN; break;
        case 30: rows_columns = TWENTY; break;
        default: rows_columns = SIX; break;
    }
    Read_Input(inName, dataType);
    sites = alignment[0];
    for (vector<double>::size_type i = 0; i != taxon.size(); i++) {
        row_of_double.push_back(0.0);
    }
    for (vector<vector<double> >::size_type i = 0; i != taxon.size(); i++) {
        mat_dobs.push_back(row_of_double);
    }
    mat_p = mat_dobs;
    mat_defs = mat_dobs;
    mat_dems = mat_dobs;
    mat_dcfs = mat_dobs;
    if (toupper(survey[0]) == 'F') {
        outfile1.open(outName1.c_str());
        outfile1 << "Program, Homo" << endl;
        outfile1 << "Version, 2.0" << endl;
        outfile1 << "Input file," << inName << endl;
        outfile1 << "Characters,";
        switch (dataType) {
            case  1: outfile1 << "Nucleotides (A|C|G|T)" << endl; break;
            case  2: outfile1 << "Nucleotides recoded (C|T|AG)" << endl; break;
            case  3: outfile1 << "Nucleotides recoded (A|G|CT)" << endl; break;
            case  4: outfile1 << "Nucleotides recoded (A|T|CG)" << endl; break;
            case  5: outfile1 << "Nucleotides recoded (C|G|AT)" << endl; break;
            case  6: outfile1 << "Nucleotides recoded (A|C|GT)" << endl; break;
            case  7: outfile1 << "Nucleotides recoded (G|T|AC)" << endl; break;
            case  8: outfile1 << "Nucleotides recoded (GT|AC)" << endl; break;
            case  9: outfile1 << "Nucleotides recoded (AG|CT)" << endl; break;
            case 10: outfile1 << "Nucleotides recoded (GC|AT)" << endl; break;
            case 11: outfile1 << "Nucleotides recoded (A|CGT)" << endl; break;
            case 12: outfile1 << "Nucleotides recoded (C|AGT)" << endl; break;
            case 13: outfile1 << "Nucleotides recoded (G|ACT)" << endl; break;
            case 14: outfile1 << "Nucleotides recoded (T|ACG)" << endl; break;
            case 15: outfile1 << "Di-nucleotides (AA|AC|...|TG|TT)" << endl; break;
            case 16: outfile1 << "Di-nucleotides 1st position (A|C|G|T)" << endl; break;
            case 17: outfile1 << "Di-nucleotides 2nd position (A|C|G|T)" << endl; break;
            case 18: outfile1 << "Codons (AAA|AAC|...|TTG|TTT)" << endl; break;
            case 19: outfile1 << "Codons 1st + 2nd positions (AA|AC|...|TG|TT)" << endl; break;
            case 20: outfile1 << "Codons 1st + 3rd positions (AA|AC|...|TG|TT)" << endl; break;
            case 21: outfile1 << "Codons 2nd + 3rd positions (AA|AC|...|TG|TT)" << endl; break;
            case 22: outfile1 << "Codons 1st + 2nd positions (A|C|G|T)" << endl; break;
            case 23: outfile1 << "Codons 1st + 3rd positions (A|C|G|T)" << endl; break;
            case 24: outfile1 << "Codons 2nd + 3rd positions (A|C|G|T)" << endl; break;
            case 25: outfile1 << "Codons 1st position (A|C|G|T)" << endl; break;
            case 26: outfile1 << "Codons 2nd position (A|C|G|T)" << endl; break;
            case 27: outfile1 << "Codons 3rd position (A|C|G|T)" << endl; break;
            case 28: outfile1 << "Genotypes (A|C|G|T|K|M|R|Y|S|W)" << endl; break;
            case 29: outfile1 << "Genotypes (A|C|G|T|K|M|R|Y|S|W|B|D|H|V)" << endl; break;
            case 30: outfile1 << "Amino acids (A|G|P|S|T|D|E|N|Q|H|K|R|M|I|V|L|W|F|Y|C)" << endl; break;
            default: outfile1 << "Recoded amino acids (AGPST|DENQ|HKR|MIVL|WFY|C) [D6]" << endl; break;
        }
        outfile1 << "Sequences," << taxon.size() << endl << endl;
        outfile1 << "Taxon 1,Taxon 2,Bowker,df,p,defs,dems,dcfs,Sites" << endl;
        cout << endl;
    }
    total = taxon.size() * (taxon.size() - 1)/2;
    counter = total;
    family_Wise_Error_Rate = family_Wise_Error_Rate/total;
    for (vector<vector<int> >::size_type iter1 = 0; iter1 != alignment.size(); ++iter1) {
        for (vector<vector<int> >::size_type iter2 = iter1 + 1; iter2 != alignment.size(); ++iter2) {
            // set all elements in divergence matrix to zero
            for (size_t i = 0; i != max_array; ++i) {
                for (size_t j = 0; j != max_array; ++j) {
                    dm[i][j] = 0;
                }
            }
            // generate divergence matrix
            for (vector<int>::size_type i = 0; i != sites.size(); ++i) {
                dm[alignment[iter1][i]][alignment[iter2][i]]++;
            }
            // generate the sum over all elements in the divergence matrix
            sum_dm = 0;
            for (size_t m = 0; m != rows_columns; ++m) {
                for (size_t n = 0; n != rows_columns; ++n) {
                    sum_dm += dm[m][n];
                }
            }
            // Bowker's matched-pairs test of symmetry
            BTS = 0.0;
            df = 0;
            for (size_t m = 0; m != rows_columns; ++m) {
                for (size_t n = m+1; n != rows_columns; ++n) {
                    if (dm[m][n] + dm[n][m] > 0) {
                        ++df;
                        BTS = BTS + ((long double)(SQR(dm[m][n] - dm[n][m])))/(dm[m][n] + dm[n][m]);
                    }
                }
            }
            if (df == 0) {
                pB = 1.0;
            } else {
                pB = xChi_Square_Distribution_Tail(BTS, df);
            }
            if (pB < min_p) {
                min_p = pB;
            }
            if (pB < family_Wise_Error_Rate) {
                ++rejectB;
            }
            if (toupper(survey[0]) == 'F') {
                mat_p[iter1][iter2] = pB;
                mat_p[iter2][iter1] = pB;
                // Compositional distance -- full symmetry
                if (df == 0) {
                    mat_dcfs[iter1][iter2] = 0.0;
                    mat_dcfs[iter2][iter1] = 0.0;
                } else {
                    d_cfs = sqrt(BTS/df);
                    mat_dcfs[iter1][iter2] = d_cfs;
                    mat_dcfs[iter2][iter1] = d_cfs;
                }
                if (mat_dcfs[iter1][iter2] > max_dcfs) {
                    max_dcfs = mat_dcfs[iter1][iter2];
                }
                if (mat_dcfs[iter1][iter2] < min_dcfs) {
                    min_dcfs = mat_dcfs[iter1][iter2];
                }
                //Euclidean distance based on distance matrix
                d_efs = 0.0;
                for (size_t m = 0; m != rows_columns; ++m) {
                    for (size_t n = m+1; n != rows_columns; ++n) {
                        // dm is defined as unsigned long so this if-else statement is needed
                        if (dm[m][n] > dm[n][m]) {
                            d_efs = d_efs + SQR(((double)(dm[m][n] - dm[n][m]))/sum_dm);
                        } else {
                            d_efs = d_efs + SQR(((double)(dm[n][m] - dm[m][n]))/sum_dm);
                        }
                    }
                }
                d_efs = sqrt(d_efs);
                mat_defs[iter1][iter2] = d_efs;
                mat_defs[iter2][iter1] = d_efs;
                //Euclidean distance based on marginal sum of off-diagonal values of distance matrix
                d_ems = 0.0;
                // set all marginal elements in divergence matrix to zero
                for (size_t i = 0; i != rows_columns; ++i) {
                    row_sum[i] = 0.0;
                    col_sum[i] = 0.0;
                }
                // generate vectors of marginal frequencies
                for (size_t i = 0; i != rows_columns; ++i) {
                    for (size_t j = 0; j != rows_columns; ++j) {
                        if (i != j) {
                            row_sum[i] += dm[i][j];
                        }
                    }
                    row_sum[i] = row_sum[i]/sum_dm;
                }
                for (size_t i = 0; i != rows_columns; ++i) {
                    for (size_t j = 0; j != rows_columns; ++j) {
                        if (i != j) {
                            col_sum[i] += dm[j][i];
                        }
                    }
                    col_sum[i] = col_sum[i]/sum_dm;
                }
                for (size_t i = 0; i != rows_columns; ++i) {
                    d_ems = d_ems + (double)(SQR(row_sum[i] - col_sum[i]));
                }
                d_ems = sqrt(d_ems);
                mat_dems[iter1][iter2] = d_ems;
                mat_dems[iter2][iter1] = d_ems;
                if (iter1 < iter2) {
                    cout << "\rNumber of comparisons left = " << --total;
                    fflush(NULL);
                    outfile1 << taxon[iter1] << ",";
                    outfile1 << taxon[iter2] << ",";
                    outfile1 << BTS << ",";
                    outfile1 << df << ",";
                    outfile1 << pB << ",";
                    outfile1 << d_efs << ",";
                    outfile1 << d_ems << ",";
                    if (df == 0) {
                        outfile1 << "Nan,";
                    } else {
                        outfile1 << d_cfs << ",";
                    }
                    outfile1 << sum_dm << endl;
                }
            }
        }
    }
    if (toupper(survey[0]) == 'F') {
        outfile1.close();
        outfile2.open(outName2.c_str());
        outfile2 << taxon.size() << endl;
        outfile3.open(outName3.c_str());
        outfile3 << taxon.size() << endl;
        outfile4.open(outName4.c_str());
        outfile4 << taxon.size() << endl;
        outfile5.open(outName5.c_str());
        outfile5 << taxon.size() << endl;
        for (vector<vector<int> >::size_type i = 0; i != taxon.size(); i++) {
            outfile2 << taxon[i];
            outfile3 << left << setw(10) << taxon[i];
            outfile4 << left << setw(10) << taxon[i];
            outfile5 << left << setw(10) << taxon[i];
            for (vector<vector<int> >::size_type j = 0; j != taxon.size(); j++) {
                outfile2 << "," << fixed << setprecision(11) << mat_p[i][j];
                outfile3 << "\t" << fixed << mat_dcfs[i][j];
                outfile4 << "\t" << fixed << mat_defs[i][j];
                outfile5 << "\t" << fixed << mat_dems[i][j];
            }
            outfile2 << endl;
            outfile3 << endl;
            outfile4 << endl;
            outfile5 << endl;
        }
        outfile2.close();
        outfile3.close();
        outfile4.close();
        outfile5.close();
        cout << endl;
        cout << endl;
    }
    if (toupper(survey[0]) == 'F') {
        cout << "--------------------------------------------------------------------" << endl;
        cout << "   SUMMARY OF ANALYSIS WITH HOMO 2.0" << endl;
        cout << endl;
        cout << "   Table with all estimates ................... " << outName1 << endl;
        cout << "   Matrix with estimates of p-values .......... " << outName2 << endl;
        cout << "   Matrix with estimates of d_cfs ............. " << outName3 << endl;
        cout << "   Matrix with estimates of d_efs ............. " << outName4 << endl;
        cout << "   Matrix with estimates of d_ems ............. " << outName5 << endl;
        cout << endl;
        cout << "   Positions .................................. " << sites.size() << endl;
        cout << "   Number of tests ............................ " << counter << endl;
        cout << "   Smallest p-value (Bowker) .................. " << scientific << min_p << endl;
        cout << "   Family-wise error rate (0.05/tests) ........ " << scientific << (double) 0.05/counter << endl;
        cout << "   Proportion of rejected tests ............... " << fixed << (double) rejectB/counter << endl;
        cout << "   Min(d_cfs) ................................. " << fixed << min_dcfs << endl;
        cout << "   Max(d_cfs) ................................. " << fixed << max_dcfs << endl;
        if (min_p < (0.05/counter)) {
            cout << endl;
            cout << "WARNING:" << endl << endl;
            cout << "   At least one pair of sequences is unlikely to have evolved" << endl;
            cout << "   under the same Markovian process. For further details, see" << endl;
            cout << "   " << outName1 << endl;
        }
        cout << "--------------------------------------------------------------------" << endl;
        cout << endl;
    } else {
        cout << "File, Taxa, Sites, Min(p), Failed" << endl;
        cout << inName << "," << taxon.size() << "," << sites.size() << "," ;
        cout << scientific << min_p << "," << fixed << (double) rejectB/counter << endl;
    }
    return 0;
}
