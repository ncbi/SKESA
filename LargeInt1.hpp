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

/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/
/** \file LargeInt<1>.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native u_int64_t type
 */

template<>  class LargeInt<1>
{
public:
    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    LargeInt<1>(const u_int64_t& c=0)  noexcept {  value[0] = c;  }

    LargeInt<1>(const std::string& kmer) noexcept : LargeInt<1>(0) {
        int sizeKmer = kmer.size();
        for (int i = 0; i < sizeKmer; i++) {
            operator<<=(2);
            value[0] += std::find(bin2NT.begin(), bin2NT.end(),  kmer[i]) - bin2NT.begin();
        }
    }  

    template <typename T>
    LargeInt<1>(const T& a, const T& b) noexcept : LargeInt<1>(0) {
        for(T i = a; i < b; ++i) {
            operator<<=(2);
            value[0] += std::find(bin2NT.begin(), bin2NT.end(),  *i) - bin2NT.begin();
        }
    }

    u_int64_t getVal () const  { return *value; }  

    static const char* getName ()  { return "LargeInt<1>"; }

    static size_t getSize ()  { return 8*sizeof(u_int64_t); }

    LargeInt<1> operator+  (const LargeInt<1>& other)   const   {  return value[0] + other.value[0];  }
    LargeInt<1> operator-  (const LargeInt<1>& other)   const   {  return value[0] - other.value[0];  }
    LargeInt<1> operator|  (const LargeInt<1>& other)   const   {  return value[0] | other.value[0];  }
    LargeInt<1> operator*  (const int& coeff)           const   {  return value[0] * coeff;        }
    LargeInt<1> operator/  (const u_int32_t& divisor)   const   {  return value[0] / divisor;      }
    u_int32_t   operator%  (const u_int32_t& divisor)   const   {  return value[0] % divisor;      }
    LargeInt<1> operator^  (const LargeInt<1>& other)   const   {  return value[0] ^ other.value[0];  }
    LargeInt<1> operator&  (const LargeInt<1>& other)   const   {  return value[0] & other.value[0];  }
    LargeInt<1> operator&  (const char& other)          const   {  return value[0] & other;        }
    LargeInt<1> operator~  ()                           const   {  return ~value[0];               }
    LargeInt<1> operator<< (const int& coeff)           const   {  return value[0] << coeff;       }
    LargeInt<1> operator>> (const int& coeff)           const   {  return value[0] >> coeff;       }
    bool        operator!= (const LargeInt<1>& c)       const   {  return value[0] != c.value[0];     }
    bool        operator== (const LargeInt<1>& c)       const   {  return value[0] == c.value[0];     }
    bool        operator<  (const LargeInt<1>& c)       const   {  return value[0] < c.value[0];      }
    bool        operator<= (const LargeInt<1>& c)       const   {  return value[0] <= c.value[0];     }

    LargeInt<1>& operator+=  (const LargeInt<1>& other)    {  value[0] += other.value[0]; return *this; }
    LargeInt<1>& operator^=  (const LargeInt<1>& other)    {  value[0] ^= other.value[0]; return *this; }

    LargeInt<1>& operator<<=  (const int& coeff)  { value[0] <<= coeff; return *this; } 
    LargeInt<1>& operator>>=  (const int& coeff)  { value[0] >>= coeff; return *this; }

    u_int8_t  operator[]  (size_t idx) const   {  return (value[0] >> (2*idx)) & 3; }

    uint8_t Codon(size_t idx) const   { return (value[0] >> (2*idx)) & 63; }
    

    /********************************************************************************/
    friend std::ostream & operator<<(std::ostream & s, const LargeInt<1> & l)
    {
        s << std::hex << l.value[0] << std::dec;  return s;
    }

    /********************************************************************************/
    /** Print corresponding kmer in ASCII
     * \param[in] sizeKmer : kmer size (def=32).
     */
    std::string toString (size_t sizeKmer) const
    {
        int i;
        u_int64_t temp = value[0];

        std::string seq(sizeKmer,'A');
        for (i=sizeKmer-1; i>=0; i--)
        {
            seq[i] = bin2NT[ temp&3 ];
            temp = temp>>2;
        }

        return seq;
    }

    u_int64_t oahash() const {
        return oahash64(value[0]);
    }    

    /********************************************************************************/
    inline static u_int64_t revcomp64 (const u_int64_t& x, size_t sizeKmer)
    {
        u_int64_t res = x;

        // OLD VERSION (with lookup table)
        // unsigned char* kmerrev  = (unsigned char *) (&(res));
        // unsigned char* kmer     = (unsigned char *) (&(x));
        // for (size_t i=0; i<8; ++i)  {  kmerrev[8-1-i] = revcomp_4NT [kmer[i]];  }

        res = ((res>> 2 & 0x3333333333333333) | (res & 0x3333333333333333) <<  2);
        res = ((res>> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) <<  4);
        res = ((res>> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) <<  8);
        res = ((res>>16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
        res = ((res>>32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
        res = res ^ 0xAAAAAAAAAAAAAAAA;
        
        return (res >> (2*(32-sizeKmer))) ;
    }

    /********************************************************************************/
    template<typename Map>
    static LargeInt<1> polynom (const char* data, size_t size, Map fct)
    {
        LargeInt<1> res (0);
        for (size_t i=0; i<size; ++i)  {  res.value[0] = 4 * res.value[0] + fct(data[i]);  }
        return res;
    }
    u_int64_t* getPointer() { return value; }
    const u_int64_t* getPointer() const { return value; }

private:
    friend LargeInt<1> revcomp (const LargeInt<1>& i,   size_t sizeKmer);

    u_int64_t value[1];
};

/********************************************************************************/
inline LargeInt<1> revcomp (const LargeInt<1>& x, size_t sizeKmer)
{
    return LargeInt<1>::revcomp64 (x.value[0], sizeKmer);
}
