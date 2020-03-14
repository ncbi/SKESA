/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

#ifndef _GENETIC_CODE_
#define _GENETIC_CODE_

#include <tuple>
#include <string>
#include <array>
#include <list>
#include <ctype.h>
#include "Model.hpp"
#include "Integer.hpp"

using namespace std;
namespace DeBruijn {
    

    class GeneticCode {
    public:
        typedef tuple<string, string, int, string> TableInfo;  // aa, starts/stops, number, name 
        GeneticCode(int g = 1) : m_tochar{'T','C','A','G'}, m_toint{{'T',0}, {'C',1}, {'A',2}, {'G',3}} {
            Init();
            m_g = -1;
            for(int i = 0; i < (int)m_data.size(); ++i) {
                if(get<2>(m_data[i]) == g) {
                    m_g = i;
                    break;
                }
            }
            if(m_g < 0)
                throw runtime_error("Unknown genetic code");

            for(int i = 0; i < 64; ++i) {
                if(get<1>(m_data[m_g])[i] == 'M') {
                    m_starts.push_back(IntToCodon(i));
                }
            }

            m_skesa_translation.resize(64);
            for(unsigned i = 0; i < 64; ++i) { // skesa codon is in reverse order; bin2NT = {'A','C','T','G'}
                string codon(3, 'A');
                codon[2] = bin2NT[i&3];
                codon[1] = bin2NT[(i>>2)&3];
                codon[0] = bin2NT[i>>4];
                m_skesa_translation[i] = AA(codon);
                if(IsStart(codon))
                    m_skesa_starts.push_back(i);
            }
        }
        char AA(const string& codon) const {
            string ucodon;
            for(char c : codon)
                ucodon.push_back(toupper(c));
            if(ucodon.find_first_not_of("ACGT") == string::npos)
                return get<0>(m_data[m_g])[CodonToInt(ucodon)];
            else
                return 'X';
        }
        string Translate(const string& nuc, bool usealts) const {
            string prot;
            int p = 0;
            if(usealts && IsStart(nuc.substr(p, 3))) {
                prot.push_back('M');
                p += 3;
            }
            for( ; p <= (int)nuc.size()-3; p += 3)
                prot.push_back(AA(nuc.substr(p, 3)));
            return prot;
        }
        char AA(int skesa_codon) const { return m_skesa_translation[skesa_codon]; }
        string Translate(const TKmer& kmer, int kmer_len, bool usealts) {
            string prot;
            int i = kmer_len-3;
            if(usealts && IsStart(kmer.Codon(i))) {
                prot.push_back('M');
                i -= 3;
            }
            for( ; i >= 0; i -= 3)
                prot.push_back(AA(kmer.Codon(i)));
            return prot;
        }
        list<string> Codons(char aa) const {
            list<string> codons;
            for(int i = 0; i < 64; ++i) {
                if(get<0>(m_data[m_g])[i] == aa) {
                    codons.push_back(IntToCodon(i));
                }
            }
            return codons;
        }       
        const list<string>& Starts() const { return m_starts; }
        bool IsStart(const string& codon) const { return find(m_starts.begin(), m_starts.end(), codon) != m_starts.end(); }
        bool IsStart(int skesa_codon) const { return find(m_skesa_starts.begin(), m_skesa_starts.end(), skesa_codon) != m_skesa_starts.end(); }
        list<string> Stops() const {
            list<string> stops;
            for(int i = 0; i < 64; ++i) {
                if(get<1>(m_data[m_g])[i] == '*') {
                    stops.push_back(IntToCodon(i));
                }
            }
            return stops;
        }
        string CodeName() const { return get<3>(m_data[m_g]); }
        const array<TableInfo, 25>& AllInfo() const { return m_data; }


    private:
        int CodonToInt(const string& codon) const { return m_toint.at(codon[0])*16+m_toint.at(codon[1])*4+m_toint.at(codon[2]); }
        string IntToCodon(int p) const {
            string codon(3, 'A');
            codon[0] = m_tochar[p>>4];
            codon[1] = m_tochar[(p>>2)&3];
            codon[2] = m_tochar[p&3];
            return codon;
        }

        int m_g;
        array<char, 4> m_tochar;
        map<char, int> m_toint;
        array<TableInfo, 25> m_data;
        string m_skesa_translation;
        list<string> m_starts;
        list<int> m_skesa_starts;
        
        void Init() {
            m_data = {
                make_tuple("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M---------------M----------------------------", 1,  "Standard Code"),
                make_tuple("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", "----------**--------------------MMMM----------**---M------------", 2,  "Vertebrate Mitochondrial Code"),
                make_tuple("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**----------------------MM---------------M------------", 3,  "Yeast Mitochondrial Code"),
                make_tuple("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------", 4,  "Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code"),
                make_tuple("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "---M------**--------------------MMMM---------------M------------", 5,  "Invertebrate Mitochondrial Code"),
                make_tuple("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------", 6,  "Ciliate, Dasycladacean and Hexamita Nuclear Code"),
                make_tuple("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------", 9,  "Echinoderm and Flatworm Mitochondrial Code"),
                make_tuple("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------", 10, "Euplotid Nuclear Code"),
                make_tuple("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------", 11, "Bacterial, Archaeal and Plant Plastid Code"),
                make_tuple("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------", 12, "Alternative Yeast Nuclear Code"),
                make_tuple("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "---M------**----------------------MM---------------M------------", 13, "Ascidian Mitochondrial Code"),
                make_tuple("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "-----------*-----------------------M----------------------------", 14, "Alternative Flatworm Mitochondrial Code"),
                make_tuple("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------*---*--------------------M----------------------------", 16, "Chlorophycean Mitochondrial Code"),
                make_tuple("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------", 21, "Trematode Mitochondrial Code"),
                make_tuple("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "------*---*---*--------------------M----------------------------", 22, "Scenedesmus obliquus Mitochondrial Code"),
                make_tuple("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--*-------**--*-----------------M--M---------------M------------", 23, "Thraustochytrium Mitochondrial Code"),
                make_tuple("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M------**-------M---------------M---------------M------------", 24, "Pterobranchia Mitochondrial Code"),
                make_tuple("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**-----------------------M---------------M------------", 25, "Candidate Division SR1 and Gracilibacteria Code"),
                make_tuple("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------", 26, "Pachysolen tannophilus Nuclear Code"),
                make_tuple("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------", 27, "Karyorelict Nuclear Code"),
                make_tuple("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*--------------------M----------------------------", 28, "Condylostoma Nuclear Code"),
                make_tuple("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------", 29, "Mesodinium Nuclear Code"),
                make_tuple("FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------", 30, "Peritrich Nuclear Code"),
                make_tuple("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------", 31, "Blastocrithidia Nuclear Code"),
                make_tuple("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M-------*-------M---------------M---------------M------------", 33, "Cephalodiscidae Mitochondrial UAA-Tyr Code")
            };
        }        
    };

}; // namespace
#endif /* _GENETIC_CODE_ */
