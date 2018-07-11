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

#ifndef _KmerInit_
#define _KmerInit_

#include "LargeInt.hpp"
#include <boost/variant.hpp>
#include <unordered_map>

/******************

This is the only place where we manipulate boost::variant directly. The rest of the code MUST use these definitions and boost::visitor

*******************/

using namespace std;
namespace DeBruijn {

#define MaxPrec 16  // kmer up to 512

    template<template<int, typename...> class BoundedType, typename... V> 
    using BoostVariant = boost::variant<BoundedType<1,V...>,  BoundedType<2,V...>,  BoundedType<3,V...>,  BoundedType<4,V...>,
                                        BoundedType<5,V...>,  BoundedType<6,V...>,  BoundedType<7,V...>,  BoundedType<8,V...>, 
                                        BoundedType<9,V...>,  BoundedType<10,V...>, BoundedType<11,V...>, BoundedType<12,V...>,
                                        BoundedType<13,V...>, BoundedType<14,V...>, BoundedType<15,V...>, BoundedType<16,V...>>;

    // for TKmer
    typedef BoostVariant<LargeInt> TLargeIntN;


    // for TKmerCount
    template<int N> using TLargeIntVec = vector<pair<LargeInt<N>,size_t>>;
    typedef BoostVariant<TLargeIntVec> TKmerCountN;
        
    // for TKmerMap
    struct SKmerHash {
        template<typename T> 
        size_t operator() (const T& kmer) const { return kmer.oahash(); }
    };
    template<int N, class V> using TLargeIntMap = unordered_map<LargeInt<N>,V,SKmerHash>;
    template<class V> using TKmerMapN = BoostVariant<TLargeIntMap, V>;
    
    // This variadic template could be used in construsctors of all boost::variants used in this code
    template<typename Variant, template<int, typename...> class BoundedType, typename... Params>
    Variant CreateVariant(int p) {
        switch(p) {
        case 1 :  return BoundedType<1, Params...>();
        case 2 :  return BoundedType<2, Params...>();
        case 3 :  return BoundedType<3, Params...>();
        case 4 :  return BoundedType<4, Params...>();
        case 5 :  return BoundedType<5, Params...>();
        case 6 :  return BoundedType<6, Params...>();
        case 7 :  return BoundedType<7, Params...>();
        case 8 :  return BoundedType<8, Params...>();
        case 9 :  return BoundedType<9, Params...>();
        case 10 : return BoundedType<10, Params...>();
        case 11 : return BoundedType<11, Params...>();
        case 12 : return BoundedType<12, Params...>();
        case 13 : return BoundedType<13, Params...>();
        case 14 : return BoundedType<14, Params...>();
        case 15 : return BoundedType<15, Params...>();
        case 16 : return BoundedType<16, Params...>();
        default :  throw runtime_error("Not supported kmer length");
        }
    }
            
}; // namespace
#endif /* _KmerInit_ */
