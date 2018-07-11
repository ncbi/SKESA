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

/** \file LargeInt.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Class that manages large integers
 *
 * arbitrary-precision integer library
 * very limited: only does what minia needs (but not what minia deserves)
 * This file holds interfaces related to the Design Pattern Observer.
 */

#ifndef _GATB_CORE_TOOLS_MATH_LARGEINT_HPP_
#define _GATB_CORE_TOOLS_MATH_LARGEINT_HPP_

/********************************************************************************/

#include <stdint.h>
#include <algorithm>
#include <iostream>
#include <array>
#include <algorithm>

#include "config.hpp"

/********************************************************************************/
namespace DeBruijn {
/********************************************************************************/

extern std::array<const unsigned  char, 256> revcomp_4NT;
extern std::array<const char, 4> bin2NT;

inline static u_int64_t oahash64 (u_int64_t elem)
{
    //    return std::hash<u_int64_t>()(elem); 
    
    u_int64_t code = elem;
    code = code ^ (code >> 14); //supp
    code = (~code) + (code << 18);
    code = code ^ (code >> 31);
    code = code * 21;
    code = code ^ (code >> 11);
    code = code + (code << 6);
    code = code ^ (code >> 22);
		
    return code;
}


/** \brief Large integer class
 *
 * The LargeInt class provides methods for integer calculus. It has a template parameter
 * 'precision' giving the number of bits used the integer representation. For instance:
 *  - LargeInt<1>  : representation of integers up to 2^64
 *  - LargeInt<2>  : representation of integers up to 2^128
 *  - etc
 *
 *  This template class has a specialization for precision=1. In this case, native 64 bits
 *  integers are used.
 *
 *  This template class may have a specialization for precision=2. If the used operating
 *  system allows it, native 128 bits integers are used.
 *
 *  In the other cases, the LargeInt provides a generic integer calculus class. Note that
 *  such an implementation could be optimized in several ways, including direct assembly
 *  code for maximum speed.
 *
 *  The LargeInt class is hugely used throughout the GATB project since it encodes kmers values.
 *
 *  The LargeInt class is mainly used with the IntegerTemplate class, where 4 specializations
 *  of LargeInt are used as template types of IntegerTemplate.
 *
 *  \see IntegerTemplate
 */
template<int precision>  class LargeInt {
public:

    /** Get the name of the class used by the variant (ie. one of the Ti template class parameters)
     * \return the class name.
     */
    static const char* getName ()
    {
        static char buffer[256];
        static bool first = true;
        if (first)  {  first = false;  snprintf (buffer, sizeof(buffer), "LargeInt<%d>", precision);  }
        return buffer;
    }

    /** Get the 64 less significant bits of the LargeInt object as a native integer type.
     * \return (part of) the LargeInt object as a native integer type.
     */
    u_int64_t getVal() const  { return this->value[0]; }

    /** Get the size of an instance of the class
     * \return the size of an object (in bits).
     */
    static const size_t getSize ()  { return 8*sizeof(u_int64_t)*precision; }

    /********************************************************************************/
    /** Constructor.
     * \param[in] val : initial value of the large integer. */
    LargeInt(const u_int64_t& val = 0) noexcept {
        value[0] = val;   
        for (int i = 1; i < precision; i++)  
            value[i] = 0;
    }

    LargeInt(const std::string& kmer) noexcept : LargeInt(0) {
        int sizeKmer = kmer.size();
        for (int i = 0; i < sizeKmer; i++) {
            operator<<=(2);
            value[0] += std::find(bin2NT.begin(), bin2NT.end(),  kmer[i]) - bin2NT.begin();
        }
    }

    template <typename T> 
    LargeInt(const T& a, const T& b) noexcept : LargeInt(0) {
        for(T i = a; i < b; ++i) {
            operator<<=(2);
            value[0] += std::find(bin2NT.begin(), bin2NT.end(),  *i) - bin2NT.begin();
        }
    }

    /********************************************************************************/
    /** Operator +
     * \param[in] other : operand
     * \return sum of object and the operand.
     */
    LargeInt operator+ (const LargeInt& other) const
    {
        LargeInt result;
        int carry = 0;
        for (int i = 0 ; i < precision ; i++)
        {
            result.value[i] = this->value[i] + other.value[i] + carry;
            carry = (result.value[i] < this->value[i]) ? 1 : 0;
        }

        return result;
    }

    /********************************************************************************/
    /** Operator -
     * \param[in] other : operand
     * \return subtraction of object and the operand.
     */
    LargeInt operator- (const LargeInt& other) const
    {
        LargeInt result;
        int carry = 0;
        for (int i = 0 ; i < precision ; i++)
        {
            result.value[i] = this->value[i] - other.value[i] - carry;
            carry = (result.value[i] > this->value[i]) ? 1 : 0;
        }

        return result;
    }

    /********************************************************************************/
    /** Operator /
     * \param[in] divisor : operand
     * \return division of the object by the divisor.
     */
    LargeInt operator/(const uint32_t& divisor) const
    {
        LargeInt result;
        std::fill( result.value, result.value + precision, 0 );

        // inspired by Divide32() from http://subversion.assembla.com/svn/pxcode/RakNet/Source/BigInt.cpp

        u_int64_t r = 0;
        uint32_t mask32bits = ~0;
        for (int i = precision-1; i >= 0; --i)
        {
            for (int j = 1; j >= 0; --j) // [j=1: high-32 bits, j=0: low-32 bits] of array[i]
            {
                u_int64_t n = (r << 32) | ((this->value[i] >> (32*j)) & mask32bits );
                result.value[i] = result.value[i] | (((n / divisor) & mask32bits) << (32*j));
                r = n % divisor;
            }
        }

        return result;
    }

    /********************************************************************************/
    /** Operator %
     * \param[in] divisor : operand
     * \return modulo of the object by the operand.
     */
    uint32_t operator%(const uint32_t& divisor) const
    {
        u_int64_t r = 0;
        uint32_t mask32bits = ~0;
        for (int i = precision-1; i >= 0; --i)
        {
            for (int j = 1; j >= 0; --j) // [j=1: high-32 bits, j=0: low-32 bits] of array[i]
            {
                u_int64_t n = (r << 32) | ((this->value[i] >> (32*j)) & mask32bits );
                r = n % divisor;
            }
        }

        return (uint32_t)r;
    }

    /********************************************************************************/
    /** Operator ^
     * \param[in] other : operand
     * \return operator^ of the object by the operand.
     */
    LargeInt operator^(const LargeInt& other) const
    {
        LargeInt result;
        for (int i=0 ; i < precision ; i++)
            result.value[i] = this->value[i] ^ other.value[i];

        return result;
    }

    /********************************************************************************/
    /** Operator |
     * \param[in] other : operand
     * \return operator| of the object by the operand.
     */
    LargeInt operator|(const LargeInt& other) const
    {
        LargeInt result;
        for (int i=0 ; i < precision ; i++)
            result.value[i] = this->value[i] | other.value[i];
        
        return result;
    }
    
    /********************************************************************************/
    /** Operator &
     * \param[in] other : operand
     * \return operator& of the object by the operand.
     */
    LargeInt operator&(const LargeInt& other) const
    {
        LargeInt result;
        for (int i=0 ; i < precision ; i++)
            result.value[i] = this->value[i] & other.value[i];

        return result;
    }

    /********************************************************************************/
    /** Operator &
     * \param[in] other : operand
     * \return operator& of the object by the operand.
     */
    LargeInt operator&(const char& other) const
    {
        LargeInt result;
        result.value[0] = this->value[0] & other;
        return result;
    }

    /********************************************************************************/
    /** Operator ~
     * \return negation of the object
     */
    LargeInt operator~() const
    {
        LargeInt result;
        for (int i=0 ; i < precision ; i++)
            result.value[i] = ~this->value[i];

        return result;
    }

    /********************************************************************************/
    /** Operator <<. Note: this method is likely to be hugely used when we want to get
     * neighbors of a given kmer encoded as a LargeInt object.
     * \param[in] coeff : operand
     * \return left shift of the object
     */
    LargeInt operator<<(const int& coeff) const
    {
        LargeInt result (0);

        int large_shift = coeff / 64;
        int small_shift = coeff % 64;

        for (int i = large_shift ; i < precision-1; i++)
        {
            result.value[i] = result.value[i] | (this->value[i-large_shift] << small_shift);

            if (small_shift == 0) // gcc "bug".. u_int64_t x; x>>64 == 1<<63, x<<64 == 1
            {
                result.value[i+1] = 0;
            }
            else
            {
                result.value[i+1] = this->value[i-large_shift] >> (64 - small_shift);
            }

        }
        result.value[precision-1] = result.value[precision-1] | (this->value[precision-1-large_shift] << small_shift);

        return result;
    }

    /********************************************************************************/
    /** Operator >>. Note: this method is likely to be hugely used when we want to get
     * neighbors of a given kmer encoded as a LargeInt object.
     * \param[in] coeff : operand
     * \return right shift of the object
     */
    LargeInt operator>>(const int& coeff) const
    {
        LargeInt result (0);

        int large_shift = coeff / 64;
        int small_shift = coeff % 64;

        result.value[0] = (this->value[large_shift] >> small_shift);

        for (int i = 1 ; i < precision - large_shift ; i++)
        {
            result.value[i] = (this->value[i+large_shift] >> small_shift);
            if (small_shift == 0) // gcc "bug".. u_int64_t x; x>>64 == 1<<63, x<<64 == 1
            {
                result.value[i-1] =  result.value[i-1];
            }
            else
            {
                result.value[i-1] =  result.value[i-1] | (this->value[i+large_shift] << (64 - small_shift));
            }
        }

        return result;
    }

    /********************************************************************************/
    /** Operator !=
     * \param[in] c : operand
     * \return inequality
     */
    bool operator!=(const LargeInt& c) const
    {
        for (int i = 0 ; i < precision ; i++)
            if( this->value[i] != c.value[i] )
                return true;
        return false;
    }

    /********************************************************************************/
    /** Operator ==
     * \param[in] c : operand
     * \return equality
     */
    bool operator==(const LargeInt& c) const
    {
        for (int i = 0 ; i < precision ; i++)
            if( this->value[i] != c.value[i] )
                return false;
        return true;
    }

    /********************************************************************************/
    /** Operator <
     * \param[in] c : operand
     */
    bool operator<(const LargeInt& c) const
    {
        for (int i = precision-1 ; i>=0 ; --i)
            if( this->value[i] != c.value[i] )
                return this->value[i] < c.value[i];

        return false;
    }

    /********************************************************************************/
    /** Operator <=
     * \param[in] c : operand
     */
    bool operator<=(const LargeInt& c) const
    {
        return operator==(c) || operator<(c);
    }

    /********************************************************************************/
    /** Operator +=
     * \param[in] other : operand
     * \return addition and affectation
     */
    LargeInt& operator+=  (const LargeInt& other)
    {
        // NOT so easy to optimize because of the carry
        *this = *this + other;
        return *this;
    }

    /********************************************************************************/
    /** Operator ^=
     * \param[in] other : operand
     * \return xor and affectation
     */
    LargeInt& operator^=  (const LargeInt& other)
    {
        for (int i=0 ; i < precision ; i++)  {  this->value[i] ^= other.value[i];  }
        return *this;
    }

    /********************************************************************************/
    /** Operator &=
     * \param[in] other : operand
     * \return and and affectation
     */
    LargeInt& operator&=  (const LargeInt& other)
    {
        for (int i=0 ; i < precision ; i++)  {  this->value[i] &= other.value[i];  }
        return *this;
    }

    /********************************************************************************/
    /** Operator |=
     * \param[in] other : operand
     * \return or and affectation
     */
    LargeInt& operator|=  (const LargeInt& other)
    {
        for (int i=0 ; i < precision ; i++)  {  this->value[i] |= other.value[i];  }
        return *this;
    }

    /********************************************************************************/
    /** Operator <<=
     * \param[in] coeff : operand
     * \return left shift and affectation
     */
    LargeInt& operator<<=  (const int& coeff)
    {
        *(this) = (*this) << coeff;  return *this;
    }

    /********************************************************************************/
    /** Operator >>=
     * \param[in] coeff : operand
     * \return right shift and affectation
     */
    LargeInt& operator>>=  (const int& coeff)
    {
        *(this) = (*this) >> coeff;  return *this;
    }

    /********************************************************************************/
    /** Output stream operator for the IntegerTemplate class
     * \param[in] s : the output stream to be used.
     * \param[in] l : the object to output
     * \return the modified output stream.
     */
    friend std::ostream & operator<<(std::ostream & s, const LargeInt<precision> & l)
    {
        int i=0;

        /** We want to display the number in hexa (easier to do...) */
        s << std::hex;

        /** We skip leading 0. */
        for (i=precision-1; i>=0 && l.value[i]==0; i--)  {}

        /** We dump the different parts of the large integer. */
        for (  ; i>=0 ; i--)  { s << l.value[i];   if (i>=1) { s << ".";  }  }

        /** We go back to decimal format. */
        s << std::dec;

        /** We return the output stream. */
        return s;
    }

    /********************************************************************************/
    /** Computes a kmer value as polynom. We may have conversion from the data buffer to
     * a nucleotide code. This is done through the provided functor.
     * \param[in] data : kmer given as a buffer of nucleotides
     * \param[in] size : size of the kmer
     * \param[in] fct : convert the ith entry in the buffer into a nucleotide code  (A=0, C=1, T=2 and G=3)
     */
    template<typename Map>
    static LargeInt polynom (const char* data, size_t size, Map fct)
    {
        LargeInt res (0);
        for (size_t i=0; i<size; ++i)  {  res = res * 4 + fct(data[i]);  }
        return res;
    }

    /********************************************************************************/
    /** Print corresponding kmer in ASCII
     * \param[in] sizeKmer : kmer size.
     * \return the ASCII string
     */
    std::string toString (size_t sizeKmer) const
    {
        std::string seq(sizeKmer,'A');
        for (size_t i=0; i<sizeKmer; i++)  {  seq[sizeKmer-i-1] = bin2NT [(*this)[i]];  }

        return seq;
    }

    /********************************************************************************/
    /** Operator[] access the ith nucleotide in the given integer. For instance a[4] get the 5th nucleotide of
     * a kmer encoded as an Integer object.
     * \param[in] idx : index of the nucleotide to be retrieved
     * \return the nucleotide value as follow: A=0, C=1, T=2 and G=3
     */
    u_int8_t  operator[]  (size_t idx) const    {  return (this->value[idx/32] >> (2*idx%64)) & 3; }

    u_int64_t oahash() const {
        // hash = XOR_of_series[hash(i-th chunk iof 64 bits)]   
        u_int64_t result = 0, chunk, mask = ~0;
        LargeInt intermediate = *this;
        for (size_t i=0;i<precision;i++) {
            chunk = (intermediate & mask).value[0];
            intermediate = intermediate >> 64;
            result ^= oahash64 (chunk);
        }
        return result;
    }
    u_int64_t* getPointer() { return value; }

private:
    template<int T>  friend LargeInt<T> revcomp (const LargeInt<T>& i, size_t sizeKmer);

    u_int64_t value[precision];
};

/********************************************************************************/
template<int precision>  inline LargeInt<precision> revcomp (const LargeInt<precision>& x, size_t sizeKmer)
{
    const LargeInt<precision> res = x;

    unsigned char* kmerrev  = (unsigned char *) (&(res.value[0]));
    unsigned char* kmer     = (unsigned char *) (&(x.value[0]));

    for (size_t i=0; i<8*precision; ++i)
    {
        kmerrev[8*precision-1-i] = revcomp_4NT [kmer[i]];
    }

    return (res >> (2*( 32*precision - sizeKmer))  ) ;
}

/********************************************************************************/
/********************     SPECIALIZATION FOR precision=1     ********************/
/********************************************************************************/
#include "LargeInt1.hpp"

/********************************************************************************/
/********************     SPECIALIZATION FOR precision=2     ********************/
/********************************************************************************/
#include "LargeInt2.hpp"

/********************************************************************************/
} /* end of namespace */
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_LARGEINT_HPP_ */
