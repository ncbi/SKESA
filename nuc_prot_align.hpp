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

#ifndef _NUC_PROT_ALIGN_
#define _NUC_PROT_ALIGN_

#include "genetic_code.hpp"
#include "glb_align.hpp"

using namespace std;
namespace DeBruijn {

class CCigar_NAtoAA : public CCigarBase {
public:
    CCigar_NAtoAA(int qto, int sto) : CCigarBase(qto, 3*sto+2) {}

    TRange SubjectRange() const { return TRange(m_sfrom/3, (m_sto+1)/3-1); }
    TCharAlign ToAlign(const  char* query, const  char* subject) const {
        TCharAlign align;
        query += m_qfrom;
        subject += m_sfrom/3;
        for(auto& element : m_elements) {
            if(element.m_type == 'M') {
                align.first.insert(align.first.end(), query, query+element.m_len);
                query += element.m_len;
                for(int l = 0; l < element.m_len/3; ++l)
                    align.second.insert(align.second.end(), 3, *subject++);
                if(element.m_len%3 != 0)
                    align.second.insert(align.second.end(), element.m_len%3, *subject);
            } else if(element.m_type == 'D') {
                align.first.insert(align.first.end(), element.m_len, '-');
                if(element.m_len%3 != 0)
                    align.second.insert(align.second.end(), element.m_len%3, *subject++);
                for(int l = 0; l < element.m_len/3; ++l)
                    align.second.insert(align.second.end(), 3, *subject++);
            } else {
                align.first.insert(align.first.end(), query, query+element.m_len);
                query += element.m_len;
                align.second.insert(align.second.end(), element.m_len, '-');
            }   
        }
        return align;
    }
    int Score(const  char* query, const  char* subject, int gopen, int gapextend, int fsopen, const char delta[256][256], const GeneticCode& gcode) const {
        int score = 0;
        const char* pstart = nullptr;
        if(m_sfrom == 0 && *subject == 'M')
            pstart = subject;

        query += m_qfrom;
        subject += m_sfrom/3;
        for(auto& element : m_elements) {
            if(element.m_type == 'M') {                
                for(int l = 0; l < element.m_len/3; ++l) { // only whole codons
                    string codon(query, query+3);
                    if(subject == pstart && gcode.IsStart(codon))
                        score += delta[(int)'M'][(int)'M'];
                    else
                        score += delta[(int)gcode.AA(codon)][(int)*subject];
                    query += 3;
                    ++subject;
                }
            } else if(element.m_type == 'D') {
                int gap = element.m_len/3;
                subject += gap;
                int opn = gopen;
                if(element.m_len%3 > 0) {
                    ++gap;
                    ++subject;
                    query += 3-element.m_len%3;
                    opn = fsopen;
                }
                score -= opn+gapextend*gap;
            } else {
                int gap = element.m_len/3;
                query += element.m_len;
                int opn = gopen;
                if(element.m_len%3 > 0) {
                    ++gap;
                    opn = fsopen;
                }
                score -= opn+gapextend*gap;
            }
        }

        return score;
    }
    void PrintAlign(const  char* query, const  char* subject, const char delta[256][256], const GeneticCode& gcode, ostream& os) const {
        bool from_prot_start = (m_sfrom == 0 && *subject == 'M');
        TCharAlign align = ToAlign(query, subject);
        os << align.first << "\n";        
        for(unsigned p = 0; p < align.first.size(); p += 3) {
            for( ; align.second[p] == '-'; ++p) 
                os << " ";            
            string codon = align.first.substr(p, 3);
            if(codon.find('-') == string::npos) {
                char aa = gcode.AA(codon);
                if(aa == align.second[p] || (from_prot_start && p == 0 && gcode.IsStart(codon)))
                    os << "|||";
                else if(delta[(int)aa][(int)align.second[p]] > 0)
                    os << "+++";
                else
                    os << "   ";
            } else {
                os << "   ";
            }
        }
        os << "\n" << align.second << "\n";
    }
};

struct SRawMemoryNAtoAA {
    SRawMemoryNAtoAA(size_t na, size_t nb) {
        for(auto& item : s)
            item = new CScore[nb+1];       // constructor called
        for(auto& item : gapb)
            item = new CScore[nb+1];
        for(auto& item : fsb1)
            item = new CScore[nb+1];
        for(auto& item : fsb2)
            item = new CScore[nb+1];
        mtrx = new uint16_t[(na+1)*(nb+1)]; //not initialised
        for(size_t i = 0; i <= 2*nb; ++i)   //two first raws
            mtrx[i] = CCigar_NAtoAA::Zero;
     }
    ~SRawMemoryNAtoAA() {
        for(auto& item : s)
            delete[] item;
        for(auto& item : gapb)
            delete[] item;
        for(auto& item : fsb1)
            delete[] item;
        for(auto& item : fsb2)
            delete[] item;
        delete[] mtrx;
    }
    void Rotate() {
        rotate(s.begin(), s.end()-1, s.end());
        rotate(gapb.begin(), gapb.end()-1, gapb.end());
        rotate(fsb1.begin(), fsb1.end()-1, fsb1.end());
        rotate(fsb2.begin(), fsb2.end()-1, fsb2.end());
    }

    array<CScore*,4> s;    // best scores in 4 last raws
    array<CScore*,4> gapb; // best gaps (whole codons) in b (insertions in a)
    array<CScore*,4> fsb1; // best frameshifts in b (3*n+1 insertions in a)
    array<CScore*,4> fsb2; // best frameshifts in b (3*n+2 insertions in a)
    uint16_t* mtrx;        // backtracking info (Astart/Bstart gap start, Agap/Bgap best score has gap and should be backtracked to Asrt/Bsart; Zero stop bactracking)
};
CCigar_NAtoAA BackTrackNAtoAA(int ia, int ib, uint16_t* m, int nb) {
    CCigar_NAtoAA track(ia, ib);
    while((ia >= 0 || ib >= 0) && !(*m&CCigar_NAtoAA::Zero)) {
        if(*m&CCigar_NAtoAA::Agap) {
            int len = 3;
            while(!(*m&CCigar_NAtoAA::Astart)) {
                len += 3;
                --m;
                --ib;
            }
            --m;
            --ib;
            track.PushFront(CCigar::SElement(len,'D'));
        } else if(*m&CCigar_NAtoAA::AgapFS1) {
            int len = 1;
            while(!(*(m-nb-1)&CCigar_NAtoAA::AstartFS1)) {
                len += 3;
                --m;
                --ib;
            }
            ia -= 2;
            m -= 2*(nb+1)+1; // shift for 2 extra a bases and 1 b base
            --ib;
            track.PushFront(CCigar::SElement(len,'D'));           
            track.PushFront(CCigar::SElement(2,'M'));  // 2 extra a bases added to left diagonal               
        } else if(*m&CCigar_NAtoAA::AgapFS2) {
            int len = 2;
            while(!(*m&CCigar_NAtoAA::AstartFS2)) {
                len += 3;
                --m;
                --ib;
            }
            --ia;
            m -= nb+2;   // shift for 1 extra a base and 1 b base
            --ib;
            track.PushFront(CCigar::SElement(len,'D'));           
            track.PushFront(CCigar::SElement(1,'M'));   // 1 extra a base added to left diagonal 
        } else if(*m&CCigar_NAtoAA::Bgap) {
            int len = 3;
            while(!(*m&CCigar_NAtoAA::Bstart)) {
                len += 3;
                m -= 3*(nb+1);
            }
            m -= 3*(nb+1);
            ia -= len;
            track.PushFront(CCigar::SElement(len,'I'));
        } else if(*m&CCigar_NAtoAA::BgapFS1) {
            int len =1;
            while(!(*m&CCigar_NAtoAA::BstartFS1)) {
                len += 3;
                m -= 3*(nb+1);                
            }
            m -= nb+1;
            ia -= len;
            track.PushFront(CCigar::SElement(len,'I'));            
        } else if(*m&CCigar_NAtoAA::BgapFS2) {
            int len =2;
            while(!(*(m-nb-1)&CCigar_NAtoAA::BstartFS2)) {
                len += 3;
                m -= 3*(nb+1);                
            }
            m -= 2*(nb+1);
            ia -= len;
            track.PushFront(CCigar::SElement(len,'I'));            
        } else {
            ia -= 3;
            --ib;
            m -= 3*(nb+1)+1; //shift for 3 a bases and 1 b base
            track.PushFront(CCigar::SElement(3,'M'));
        }
    }

    return track;
}

CCigar_NAtoAA LclAlignNAtoAA(const  char* a, int na, const  char*  b, int nb, int rho, int sigma, int rhofs, const char delta[256][256], const GeneticCode& gcode) {
    // rho - new gap penalty (one base gap rho+sigma)
    // sigma - extension penalty
    // rhofs - new frmeshift penalty (one base gap rho+sigma)

    SRawMemoryNAtoAA memory(na, nb);
	CScore*& s = memory.s[0];         // best scores in current a-raw
	CScore*& sm1 = memory.s[1];       // best scores in -1 a-raw
	CScore*& sm2 = memory.s[2];       // best scores in -2 a-raw
	CScore*& sm3 = memory.s[3];       // best scores in -3 a-raw
    CScore*& gapb = memory.gapb[0];   // best score with b-gap (3*n insertion in a) in current a-raw
    CScore*& gapbm3 = memory.gapb[3]; // best score with b-gap (3*n insertion in a) in -3 a-raw
    CScore*& fsb1 = memory.fsb1[0];   // best score with b-frameshift (3*n+1 insertion in a) in current a-raw
    CScore*& fsb1m3 = memory.fsb1[3]; // best score with b-frameshift (3*n+1 insertion in a) in -3 a-raw
    CScore*& fsb2 = memory.fsb2[0];   // best score with b-frameshift (3*n+2 insertion in a) in current a-raw
    CScore*& fsb2m3 = memory.fsb2[3]; // best score with b-frameshift (3*n+2 insertion in a) in -3 a-raw
    uint16_t* mtrx = memory.mtrx;     // backtracking info (Astart/Bstart gap start, Agap/Bgap best score has gap and should be backtracked to Asrt/Bsart; Zero stop bactracking)

    CScore rsa(-rho-sigma, 0);        // new gapa
    CScore rsb(-rho-sigma, 1);        // new gapb  
    CScore rsafs(-rhofs-sigma, 0);    // new a-frameshift
    CScore rsbfs(-rhofs-sigma, 1);    // new b-frameshift

    CScore max_score;
    uint16_t* max_ptr = mtrx;
    uint16_t* m = mtrx+nb+2*(nb+1); // shift to 3rd raw
    for(int i = 2; i < na; ++i) {
		*(++m) = CCigar_NAtoAA::Zero;
        string codon(a+i-2, a+i+1);
        const char* matrix = delta[(int)gcode.AA(codon)];
        CScore gapa; // 3*n deletion from a
        CScore fsa1; // 3*n+1 deletion from a
        CScore fsa2; // 3*n+2 deletion from a
		for(int j = 0; j < nb; ) {
            // b[j] current aa
            // any_score[j] score for previous aa
            // any_score[j+1] score for current aa
			*(++m) = 0;
			CScore ss = sm3[j]+CScore(matrix[(int)b[j]], 1);  // diagonal extension
            if(j == 0 && b[j] == 'M' && gcode.IsStart(codon)) // possible alt start
                ss = sm3[j]+CScore(delta[(int)'M'][(int)'M'], 1);

            gapa += CScore(-sigma, 0);                 // gapa extension (3*n deletion)
            if(s[j]+rsa > gapa) {
                gapa = s[j]+rsa;                       // new gapa
                *m |= CCigar_NAtoAA::Astart;           // the end of last aligned NA codon, first inserted AA
            }
            fsa1 += CScore(-sigma, 0);                 // fsa extension (3*n+1 deletion)
            if(sm2[j]+rsafs > fsa1) {
                fsa1 = sm2[j]+rsafs;
                *(m-nb-1) |= CCigar_NAtoAA::AstartFS2; // first not aligned NA base, first inserted AA
            }
            fsa2 += CScore(-sigma, 0);                 // fsa extension (3*n+2 deletion)
            if(sm1[j]+rsafs > fsa2) {
                fsa2 = sm1[j]+rsafs;
                *m |= CCigar_NAtoAA::AstartFS1;        // first not aligned NA base, first inserted AA
            }
			
			CScore& gapbj = gapb[++j];                 // j increased here
			gapbj = gapbm3[j]+CScore(-sigma, 1);       // gapb extension
			if(sm3[j]+rsb > gapbj) {
				gapbj = sm3[j]+rsb;                    // new gapb
				*m |= CCigar_NAtoAA::Bstart;           // the end of first insertion codon in NA, last aligned AA
			}
            CScore& fsb1j = fsb1[j];
            fsb1j = fsb1m3[j]+CScore(-sigma, 1);       // fsb1 extension
            if(sm1[j]+rsbfs > fsb1j) {
                fsb1j = sm1[j]+rsbfs;             
                *m |= CCigar_NAtoAA::BstartFS1;        // first not aligned NA base, last aligned AA        
            }
            CScore& fsb2j = fsb2[j];
            fsb2j = fsb2m3[j]+CScore(-sigma, 1);       // fsb2 extension
            if(sm2[j]+rsbfs > fsb2j) {
                fsb2j = sm2[j]+rsbfs;             
                *(m-nb-1) |= CCigar_NAtoAA::BstartFS2; // first not aligned NA base, last aligned AA        
            }

            int type = 0;
            CScore score = ss;
            if(gapa > score) {
                score = gapa;
                type = CCigar_NAtoAA::Agap;
            }
            if(fsa1 > score) {
                score = fsa1;
                type = CCigar_NAtoAA::AgapFS1;
            }
            if(fsa2 > score) {
                score = fsa2;
                type = CCigar_NAtoAA::AgapFS2;
            }
            if(gapbj > score) {
                score = gapbj;
                type = CCigar_NAtoAA::Bgap;
            }
            if(fsb1j > score) {
                score = fsb1j;
                type = CCigar_NAtoAA::BgapFS1;
            }
            if(fsb2j > score) {
                score = fsb2j;
                type = CCigar_NAtoAA::BgapFS2;
            }

            s[j] = score;
            *m |= type;

            if(type == 0 && score > max_score) {
                max_score = score;
                max_ptr = m;
            }
            if(s[j].Score() <= 0) {
                s[j] = CScore();
                *m |= CCigar_NAtoAA::Zero;  
            }
		}
        memory.Rotate();		
	}

    int ia = (max_ptr-mtrx)/(nb+1)-1;
    int ib = (max_ptr-mtrx)%(nb+1)-1;
    m = max_ptr;
    return BackTrackNAtoAA(ia, ib, m, nb);
}

} // namespace
#endif  // _NUC_PROT_ALIGN_

