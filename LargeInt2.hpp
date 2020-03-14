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

/** \file LargeInt<2>.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Integer class relying on native u_int64_t type
 */

/********************************************************************************/
#if  INT128_FOUND == 1
/********************************************************************************/

u_int64_t revcomp64 (const u_int64_t& x, size_t sizeKmer)
{
    u_int64_t res = x;

    unsigned char* kmerrev  = (unsigned char *) (&(res));
    unsigned char* kmer     = (unsigned char *) (&(x));

    for (size_t i=0; i<8; ++i)  {  kmerrev[8-1-i] = revcomp_4NT [kmer[i]];  }

    return (res >> (2*( 32 - sizeKmer))) ;
}

template<>  class LargeInt<2>
{
public:
    /** Constructor.
     * \param[in] c : initial value of the large integer. */
    LargeInt<2>(const __uint128_t& c=0) noexcept {  value[0] = c;  }

    LargeInt<2>(const std::string& kmer) noexcept : LargeInt<2>(0) {
        int sizeKmer = kmer.size();
        for (int i = 0; i < sizeKmer; i++) {
            operator<<=(2);
            value[0] += std::find(bin2NT.begin(), bin2NT.end(),  kmer[i]) - bin2NT.begin();
        }
    }  

    template <typename T>
    LargeInt<2>(const T& a, const T& b) noexcept : LargeInt<2>(0) {
        for(T i = a; i < b; ++i) {
            operator<<=(2);
            value[0] += std::find(bin2NT.begin(), bin2NT.end(),  *i) - bin2NT.begin();
        }
    }

    u_int64_t getVal () const  { return *value; }

    static const char* getName ()  { return "LargeInt<2>"; }

    static const size_t getSize ()  { return 8*sizeof(__uint128_t); }

    LargeInt<2> operator+  (const LargeInt<2>& other)     const   {  return value[0] + other.value[0];  }
    LargeInt<2> operator-  (const LargeInt<2>& other)     const   {  return value[0] - other.value[0];  }
    LargeInt<2> operator|  (const LargeInt<2>& other)     const   {  return value[0] | other.value[0];  }
    LargeInt<2> operator*  (const int& coeff)              const   {  return value[0] * coeff;        }
    LargeInt<2> operator/  (const u_int32_t& divisor)      const   {  return value[0] / divisor;      }
    u_int32_t    operator%  (const u_int32_t& divisor)      const   {  return value[0] % divisor;      }
    LargeInt<2> operator^  (const LargeInt<2>& other)     const   {  return value[0] ^ other.value[0];  }
    LargeInt<2> operator&  (const LargeInt<2>& other)     const   {  return value[0] & other.value[0];  }
    LargeInt<2> operator&  (const char& other)             const   {  return value[0] & other;        }
    LargeInt<2> operator~  ()                              const   {  return ~value[0];               }
    LargeInt<2> operator<< (const int& coeff)              const   {  return value[0] << coeff;       }
    LargeInt<2> operator>> (const int& coeff)              const   {  return value[0] >> coeff;       }
    bool         operator!= (const LargeInt<2>& c)         const   {  return value[0] != c.value[0];     }
    bool         operator== (const LargeInt<2>& c)         const   {  return value[0] == c.value[0];     }
    bool         operator<  (const LargeInt<2>& c)         const   {  return value[0] < c.value[0];      }
    bool         operator<= (const LargeInt<2>& c)         const   {  return value[0] <= c.value[0];     }

    LargeInt<2>& operator+=  (const LargeInt<2>& other)    {  value[0] += other.value[0]; return *this; }
    LargeInt<2>& operator^=  (const LargeInt<2>& other)    {  value[0] ^= other.value[0]; return *this; }

    LargeInt<2>& operator<<=  (const int& coeff)  { value[0] <<= coeff; return *this; } 
    LargeInt<2>& operator>>=  (const int& coeff)  { value[0] >>= coeff; return *this; }

    u_int8_t  operator[]  (size_t idx) const   {  return (value[0] >> (2*idx)) & 3; }

    uint8_t Codon(size_t idx) const   { return (value[0] >> (2*idx)) & 63; }

    /** Output stream overload. NOTE: for easier process, dump the value in hexadecimal.
     * \param[in] os : the output stream
     * \param[in] in : the integer value to be output.
     * \return the output stream.
     */
    friend std::ostream & operator<<(std::ostream & os, const LargeInt<2> & in)
    {
        __uint128_t x = in.value[0];

        u_int64_t high_nucl = (u_int64_t) (x>>64);
        u_int64_t low_nucl  = (u_int64_t)(x&((((__uint128_t)1)<<64)-1));

        if (high_nucl == 0) {   os << std::hex <<                     low_nucl << std::dec;  }
        else                {   os << std::hex << high_nucl << "." << low_nucl << std::dec;  }
        return os;
    }
    
    /********************************************************************************/
    /** Print corresponding kmer in ASCII
     * \param[in] sizeKmer : kmer size (def=32).
     */
    std::string toString (size_t sizeKmer) const
    {
        std::string seq(sizeKmer,'A');
        for (size_t i=0; i<sizeKmer; i++)  {  seq[sizeKmer-i-1] = bin2NT [(*this)[i]];  }

        return seq;
    }

    u_int64_t oahash() const {
        return oahash64((u_int64_t)(value[0]>>64)) ^ oahash64 ((u_int64_t)(value[0]&((((__uint128_t)1)<<64)-1)));
    }    
    
    /********************************************************************************/
    template<typename Map>
    static LargeInt<2> polynom (const char* data, size_t size, Map fct)
    {
        LargeInt<2> res (0);
        for (size_t i=0; i<size; ++i)  {  res.value[0] = res.value[0] * 4 + fct(data[i]);  }
        return res;
    }
    u_int64_t* getPointer() { return reinterpret_cast<u_int64_t*>(value); }
    const u_int64_t* getPointer() const { return reinterpret_cast<const u_int64_t*>(value); }

private:
    friend LargeInt<2> revcomp (const LargeInt<2>& i,   size_t sizeKmer);

    __uint128_t value[1];
};

/********************************************************************************/
inline LargeInt<2> revcomp (const LargeInt<2>& in, size_t sizeKmer)
{
    //                  ---64bits--   ---64bits--
    // original kmer: [__high_nucl__|__low_nucl___]
    //
    // ex:            [         AC  | .......TG   ]
    //
    //revcomp:        [         CA  | .......GT   ]
    //                 \_low_nucl__/\high_nucl/

    const __uint128_t& x = in.value[0];

    u_int64_t high_nucl = (u_int64_t)(x>>64);
    int nb_high_nucl = sizeKmer>32?sizeKmer - 32:0;

    __uint128_t revcomp_high_nucl = revcomp64 (high_nucl, nb_high_nucl);

    if (sizeKmer<=32) revcomp_high_nucl = 0; // srsly dunno why this is needed. gcc bug? u_int64_t x ---> (x>>64) != 0

    u_int64_t low_nucl = (u_int64_t)(x&((((__uint128_t)1)<<64)-1));
    int nb_low_nucl = sizeKmer>32?32:sizeKmer;

    __uint128_t revcomp_low_nucl = revcomp64 (low_nucl, nb_low_nucl);

    return (revcomp_low_nucl<<(2*nb_high_nucl)) + revcomp_high_nucl;
}

/********************************************************************************/
#endif //INT128_FOUND
/********************************************************************************/
