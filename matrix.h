//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//  Copyright (c) 2004 Riccardo Riedl (www.riccardo-riedl.de)
//
//  Permission to use, copy, modify, distribute and sell this software and its
//  documentation for any purpose is hereby granted without fee, provided that
//  the above copyright notice appears in all copies and that both that
//  copyright notice and this permission notice appear in supporting
//  documentation.
//
//  There are no representations made about the suitability of this software
//  for any purpose. It is provided "as is" without expressed or implied
//  warranty.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//  filename:       matrix.h
//  org. location:  http://www.riccardo-riedl.de/content/source/matrix/matrix.h
//  last update:    06-16-2004
//
//  author:         riccardo riedl
//  contact:        http://www.riccardo-riedl.de/mail.php
//
//  description:    this source code defines and declares the c++
//                  template class 'matrix'. this class implements a simple
//                  generic matrix
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//=============================================================================
//  naming convention for members and parameters
//-----------------------------------------------------------------------------
//
//  width:  n   column-index: j
//  height: m   row-index:    i
//  value:  a   linear-index: l
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//  exception class hierarchy
//-----------------------------------------------------------------------------
//
//                                         |--- eBadMatrixAlloc
//                         |--- eMemory ---|
//                         |               |--- eBadBufferAlloc
//  eMatrix::eException ---| 
//                         |               |--- eBadIndex
//                         |               |
//                         |               |--- eSizesMismatch
//                         |---  eSize  ---|    
//                                         |--- eEmpty
//                                         |
//                                         |--- eNoSquare
//
//-----------------------------------------------------------------------------
//  note:   an error message of a catched eException object can
//          be printed to an offstream with the << operator
//-----------------------------------------------------------------------------


//=============================================================================
//  preprocessor directives
//-----------------------------------------------------------------------------

#pragma once

#include <iostream>
#include <vector>

//  defines if an exception shall be thrown when it is tried to resize a matrix
//  with already allocated memory, since this will overwrite the existing
//  values (this does not apply to dynamically increasing or decrising the
//  matrix with AddColumn, AddRow, DeleteColum and DeleteRow)
//
#define RESIZE_EXCEPTION_ON

//  the maximum length of an exception message
//
#define MAX_E_STRING 128

//  by default any index within the matrix starts at 0
//  with #define INDEX_START 1 in the source file before including this header
//  file the start can be set to 1 (only 0 or 1 is allowed and I guess only 0
//  and 1 make sense)
//
#ifndef INDEX_START
    #define INDEX_START 0
#elif (INDEX_START != 0) && (INDEX_START != 1)
    #error INDEX_START must be defined as 0 or 1.
#endif

//  when compiling optimize the class for speed and not for size and set the
//  warning level to 4 (all this will be reset at the end of this file)
//
#pragma optimize("t", on)
#pragma warning(push, 4)

//  definition of index type and the larger type for internal linear addressing
//  (one ulong must be able to hold (ushort * ushort) )
//
typedef unsigned short int ushort;
typedef unsigned long  int ulong;


//  localizable exception messages
//
const char szNoMessage    [] = "No message available for this exception.\n";
const char szBadAlloc     [] = "Failed to allocate %lu bytes.\n";
const char szBadIndex     [] = "Illegal index access with row-index %u and column-index %u in a (%u|%u)-matrix.\n";
const char szSizeMismatch [] = "Illegal operation for (%u|%u)- and (%u|%u)-matrix.\n";
const char szNoSquare     [] = "Illegal operation for (%u|%u)-matrix, where matrix must be square.\n";
const char szEmpty        [] = "Illegal operation. Matrix must not be empty.\n";


namespace eMatrix
{

    using std::ostream;

    //=========================================================================
    //  base class for exceptions
    //-------------------------------------------------------------------------
    
    //  holds a message and can print it to an output stream
    //
    class eException {
    protected:
        char m_str[MAX_E_STRING];
    public:    
        eException(void) { ; }
        eException(const char* str) { strncpy(m_str, str, MAX_E_STRING); }
        friend ostream& operator << (ostream& os, const eException& e) {
            os << ( strlen(e.m_str) ? e.m_str : szNoMessage );
            return os;
        }
    } ;


    //=========================================================================
    //  categories for exception class
    //-------------------------------------------------------------------------
    
    //  eMemory indicates a failure to allocate memory and
    //  is passing the size in bytes of the requested memory
    //
    struct eMemory : protected eException {
        ulong m_memsize;
        eMemory(ulong memsize) : m_memsize(memsize) {
            _snprintf(m_str, MAX_E_STRING, szBadAlloc, m_memsize);
        }
    } ;    
    
    //  eSize indicates an invalid dimension or index for an operation
    //
    struct eSize : protected eException {
        ushort m_m, m_n;
        eSize(ushort m, ushort n) : m_m(m), m_n(n) { ; }
    } ;


    //=========================================================================
    //  memory exception classes
    //-------------------------------------------------------------------------
    
    //  thrown when allocation of memory for the matrix failed
    //
    struct eBadMatrixAlloc : private eMemory { 
        eBadMatrixAlloc(ulong size) : eMemory(size) { ; }
    } ;

    //  thrown when allocation for an internal buffer failed,
    //  which was needed for a certain operation
    //  member m_fMatrixValid will indicated if
    //  the matrix itself will be corrupted as well
    //
    struct eBadBufferAlloc : private eMemory {
        bool m_fValid;
        eBadBufferAlloc(ulong size, bool fValid) : eMemory(size), m_fValid(fValid) { ; }
    } ;


    //=========================================================================
    //  dimension exception classes
    //-------------------------------------------------------------------------

    //  indicates that an operation was tried
    //  for which the matrix must not be empty
    //
    struct eEmpty : protected eSize {
        eEmpty(ushort m, ushort n) : eSize(m, n) {
            _snprintf(m_str, MAX_E_STRING, szEmpty);
        }
    } ;

    //  indicates that an operation was tried for
    //  which the matrix must have equal height and width
    //
    struct eNoSquare : private eSize {
        eNoSquare(ushort m, ushort n) : eSize(m, n) {
            _snprintf(m_str, MAX_E_STRING, szNoSquare, m_m, m_n);
        }
    } ;

    //  indicates that an operation was tried with to
    //  matrix for which they must have the same dimension
    //
    struct eSizesMismatch : private eSize { 
        ushort m_m2, m_n2;
        eSizesMismatch(ushort m1, ushort n1, ushort m2, ushort n2) : eSize(m1, n1), m_m2(m2), m_n2(n2) {
            _snprintf(m_str, MAX_E_STRING, szSizeMismatch, m_m, m_n, m_m2, m_n2);
        }
    } ;
        
    //  indicates attempt to access the matrix which an illegal index
    //
    struct eBadIndex : private eSize {
        ushort m_i, m_j;
        eBadIndex(ushort i, ushort j, ushort m, ushort n) : eSize(m, n), m_i(i), m_j(j) {
            _snprintf(m_str, MAX_E_STRING, szBadIndex, m_i, m_j, m_m, m_n);
        } 
    } ;

}   //  end of namespace eMatrix


using namespace eMatrix;

using std::vector;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;

template <class TYPE> class matrix
{

    //=========================================================================
    //  data members
    //-------------------------------------------------------------------------

    TYPE*  m_pMatrix;   //  internal buffer to hold values stored in matrix
    ushort m_n;         //  current width of the matrix
    ushort m_m;         //  current height of the matrix

public:

    //=========================================================================
    //  construction and destruction
    //-------------------------------------------------------------------------

    matrix (void);
    matrix (const ushort m, const ushort n, const TYPE a = (TYPE) 0.0);
    matrix (const matrix& src);

    virtual ~matrix (void);

    //=========================================================================
    //  operators
    //-------------------------------------------------------------------------

    //  assignment
    //
    matrix& operator = (const matrix& src);

    //  file in- and output
    //
    friend ofstream& operator << (ofstream& os, const matrix& src);
    friend ifstream& operator >> (ifstream& is, matrix& dest);

    //  math operators
    //
    friend matrix       operator + (matrix& m1, matrix& m2);
    friend matrix       operator - (matrix& m1, matrix& m2);
    friend matrix       operator * (matrix& m1, matrix& m2);
    friend matrix       operator * (matrix& m1, TYPE a);
    friend vector<TYPE> operator * (matrix& m,  vector<TYPE>& v);

    //  compare
    //
    bool    operator == (const matrix& src);    //  compares dimension and values
    bool    operator && (const matrix& src);    //  compares dimension only

    //=========================================================================
    //  dimension manipulation and query
    //-------------------------------------------------------------------------
    
    void SetSize  (const ushort m, const ushort n);          //  (m|n)-matrix
    void SetSize  (const ushort n) { return SetSize(n, n); } //  (n|n)-matrix

    void ClearMem (void); //  deletes all memory

    ushort Height (void) const { return m_m; }
    ushort Width  (void) const { return m_n; }

    //=========================================================================
    //  fields manipulation and query
    //-------------------------------------------------------------------------

    TYPE GetAt    (const ushort i, const ushort j);
    void SetAt    (const ushort i, const ushort j, TYPE a);
    void SetAll   (TYPE a);
    void Zero     (void) { SetAll( (TYPE) 0.0 ); }

    //=========================================================================
    //  dynamical resizing by adding or deleting rows or columns
    //-------------------------------------------------------------------------

    ushort AddRow     (void);        //  returns new height
    ushort AddColumn  (void);        //  returns new width

    ushort DeleteRow    (ushort i);   //  returns new height
    ushort DeleteColumn (ushort j);   //  returns new width

    //=========================================================================
    //  debug
    //-------------------------------------------------------------------------

    void Dump (bool wide = true );   //  prints matrix to cerr

private:

    //=========================================================================
    //  macros
    //-------------------------------------------------------------------------

    inline bool EMPTY (void) { return ( ((ulong) m_n) * ((ulong) m_m)) ? false : true; }
    inline bool VALID (ushort i, ushort j) { return ( (i >= 0) && (i < m_m) && (j >= 0) && (j < m_n) ); }
    
    inline ulong LINADDRESS (ushort i, ushort j) { return ( (((ulong)i) * ((ulong)m_n)) + ((ulong)j) ); }
    inline ulong DIMENSION  (void) { return (((ulong)m_n) * ((ulong)m_m)); }
    inline ulong MEMSIZE    (void) { return ((((ulong)m_n) * ((ulong)m_m)) * sizeof(TYPE) ); }
    
} ;


//=============================================================================
//  construction and destruction
//-----------------------------------------------------------------------------

template <class TYPE> matrix<TYPE>::matrix()
    : m_m(0L), m_n(0L), m_pMatrix(0) {
    ;
}

template <class TYPE> matrix<TYPE>::matrix (const ushort m, const ushort n, const TYPE a /* = (TYPE) 0.0*/)
    : m_m(m), m_n(n), m_pMatrix(0) {
    if (!EMPTY()) {
        m_pMatrix = new TYPE[DIMENSION()];
        if (!m_pMatrix) {
            throw eMatrix::eBadMatrixAlloc(MEMSIZE());
        }
        SetAll(a);
    }
}

template <class TYPE> matrix<TYPE>::matrix (const matrix& src) : m_pMatrix(0) {
    (*this) = src;
}

template <class TYPE> matrix<TYPE>::~matrix(void) {
    if (m_pMatrix) {
        delete m_pMatrix;
    }
}


//=============================================================================
//  file in- and output
//-----------------------------------------------------------------------------

template <class TYPE> ofstream& operator << (ofstream& os, const matrix<TYPE>& src) {
    try {
        os.write((char*) &src.m_n,      (streamsize) sizeof(ushort));
        os.write((char*) &src.m_m,      (streamsize) sizeof(ushort));
        os.write((char*) src.m_pMatrix, (streamsize) src.MEMSIZE() );
        return os;	
    }
    catch (...) {
        throw;
    }
}

template <class TYPE> ifstream& operator >> (ifstream& is, matrix<TYPE>& dest) {
    try {
        is.read((char*) &dest.m_n, (streamsize) sizeof(ushort));
        is.read((char*) &dest.m_m, (streamsize) sizeof(ushort));
        if (dest.m_pMatrix) {
            delete dest.m_pMatrix;
        }
        dest.m_pMatrix  = new TYPE[dest.DIMENSION()];
        streamsize size = dest.MEMSIZE();
        is.read((char*) dest.m_pMatrix, size);
        return is;
    }
    catch (...) {
        throw;
    }
}


//=============================================================================
//  math operators
//-----------------------------------------------------------------------------

template <class TYPE> matrix<TYPE> operator + (matrix<TYPE>& m1, matrix<TYPE>& m2) {
    if (m1.EMPTY() || m2.EMPTY()) {
        throw eEmpty(m1.m_m, m1.m_n);
    }
    if (!(m1 && m2)) {  //  compares dimensions
        throw eSizesMismatch(m1.m_m, m1.m_n, m2.m_m, m2.m_n);
    }
    try {
        matrix<TYPE> res;
        res = m1;
        for (ushort i = 0; i < res.m_m; i++) {
            for (ushort j = 0; j < res.m_n; j++) {
                res.m_pMatrix[res.LINADDRESS(i, j)] += + m2.m_pMatrix[res.LINADDRESS(i, j)];
            }
        }
        return res;
    }
    catch(...) {
        throw;
    }
}

template <class TYPE> matrix<TYPE> operator - (matrix<TYPE>& m1, matrix<TYPE>& m2) {
    if (m1.EMPTY() || m2.EMPTY()) {
        throw eEmpty(m1.m_m, m1.m_n);
    }
    if (!(m1 && m2)) {  //  compares dimensions
        throw eSizesMismatch(m1.m_m, m1.m_n, m2.m_m, m2.m_n);
    }
    try {
        matrix<TYPE> res;
        res = m1;
        for (ushort i = 0; i < res.m_m; i++) {
            for (ushort j = 0; j < res.m_n; j++) {
                res.m_pMatrix[res.LINADDRESS(i, j)] -= + m2.m_pMatrix[res.LINADDRESS(i, j)];
            }
        }
        return res;
    }
    catch(...) {
        throw;
    }
}

template <class TYPE> vector<TYPE> operator * (matrix<TYPE>& m, vector<TYPE>& v) {
    if (m.m_n != v.size()) {
        throw eSizesMismatch(m.m_m, m.m_n, v.size(), 1);
    }
    if (!m.DIMENSION()) {
        throw eEmpty(m.m_m, m.m_n);
    }
    vector<TYPE> res(m.m_m);
    int a;
    try {
        for (ushort i = 0; i < m.m_m; i++) {
            a = (TYPE) 0.0;
            for (ushort j = 0; j < m.m_n; j++) {
                a += m.m_pMatrix[m.LINADDRESS(i, j)] * v[j];
            }
            res[i] = a;
        }
        return res;
    }
    catch(...) {
        throw;
    }
}

template <class TYPE> matrix<TYPE> operator * (matrix<TYPE>& m1, matrix<TYPE>& m2) {
    if (!m1.DIMENSION()) {
        throw eEmpty(m1.m_m, m1.m_n);
    }
    if (m1.m_n != m2.m_m) {
        throw eSizesMismatch(m1.m_m, m1.m_n, m2.m_m, m2.m_n);
    }
    try {
        matrix<TYPE> res(m2.m_n, m1.m_m);
        for (ushort i = 0; i < m2.m_m; i++) {
            for (ushort j = 0; j < m1.m_n; j++) {
                TYPE a = (TYPE) 0.0;
                for (ushort x = 0; x < m1.m_n; x++) {
                    a += (m1.m_pMatrix[res.LINADDRESS(i, x)] * m2.m_pMatrix[res.LINADDRESS(x, j)]);
                }    
                res.m_pMatrix[res.LINADDRESS(i, j)] = a;
            }
        }
        return res;
    }
    catch(...) {
        throw;
    }
}

template <class TYPE> matrix<TYPE> operator * (matrix<TYPE>& m1, TYPE a) {
    if (!m1.DIMENSION()) {
        throw eEmpty(m1.m_m, m1.m_n);
    }
    try {
        matrix<TYPE> res;
        res.SetSize(m1.m_m, m1.m_n);
        for (ushort i = 0; i < m1.m_m; i++) {
            for (ushort j = 0; j < m1.m_n; j++) {
                res.m_pMatrix[res.LINADDRESS(i, j)] = (m1.m_pMatrix[m1.LINADDRESS(i, j)] * a);
            }
        }
        return res;
    }
    catch (...) {
        throw;
    }
}


//=============================================================================
//  compare
//-----------------------------------------------------------------------------

template <class TYPE> bool matrix<TYPE>::operator == (const matrix<TYPE>& src) {
    if ( operator &&(src) ) {   //  compares dimensions
        if (m_n && m_m) {
            if (memcmp(m_pMatrix, src.m_pMatrix, MEMSIZE())) {
                return false;
            }
        }
        return true;
    }
    return false;
}

template <class TYPE> bool matrix<TYPE>::operator && (const matrix<TYPE>& src) {
    if ( (m_n == src.m_n) && (m_m == src.m_m) ) {
        return true;
    } else {
        return false;
    }
}

//=============================================================================
//  assignment
//-----------------------------------------------------------------------------

template <class TYPE>
matrix<TYPE>& matrix<TYPE>::operator = (const matrix<TYPE>& src)
{
    if (this != &src) {
        m_n = src.m_n;
        m_m = src.m_m;
        if (m_pMatrix) {
            delete m_pMatrix;
        }
        m_pMatrix = new TYPE[DIMENSION()];
        if (!m_pMatrix) {
            throw eBadMatrixAlloc(MEMSIZE());
        }
        memcpy(m_pMatrix, src.m_pMatrix, MEMSIZE());
    }
    return (*this);
}


//=============================================================================
//  functions to resize the matrix and query its size
//-----------------------------------------------------------------------------

template <class TYPE> void matrix<TYPE>::SetSize(const ushort m, const ushort n) {
    m_m = m; 
    m_n = n;
    if (m_pMatrix) {
#ifdef RESIZE_EXCEPTION_ON
        throw eMatrix::eException("\nTried to resize an already allocated matrix\n");
#else
        delete m_pMatrix;
    }
    if (!EMPTY()) {
        m_pMatrix = new TYPE[DIMENSION()];
        if (!m_pMatrix) {
            throw eMatrix::eBadMatrixAlloc(MEMSIZE());
        }
        Zero();
#endif
    }
}

template <class TYPE> void matrix<TYPE>::ClearMem (void) {
    if (m_pMatrix) {
        delete m_pMatrix;
        m_pMatrix = 0;
    }
    m_m = 0;
    m_n = 0;
}

//=============================================================================
//  get and set values from the matrix
//-----------------------------------------------------------------------------

#define __I(x) (x - INDEX_START)

template <class TYPE> TYPE matrix<TYPE>::GetAt(const ushort i, const ushort j) {
    if (!VALID(__I(i), __I(j))) {
        throw eMatrix::eBadIndex(i, j, m_m, m_n);
    }
    return (m_pMatrix[LINADDRESS(__I(i), __I(j))]);
}

template <class T> void matrix<T>::SetAt(const ushort i, const ushort j, TYPE a) {
    if (!VALID(__I(i), __I(j))) {
        throw eMatrix::eBadIndex(i, j, m_m, m_n);
    }
    m_pMatrix[LINADDRESS(__I(i), __I(j))] = a;
}

template <class TYPE> void matrix<TYPE>::SetAll(TYPE value) {
    for (ulong x = 0; x < DIMENSION(); x++) {
        m_pMatrix[x] = value;
    }
}

//=============================================================================
//  expand matrix by adding and removing columns or rows
//-----------------------------------------------------------------------------

template <class TYPE> ushort matrix<TYPE>::AddColumn(void) {   
    TYPE* pBuffer = 0;
    if (!EMPTY()) {
        pBuffer = new TYPE[DIMENSION()];
        if (!pBuffer) {
            throw eBadBufferAlloc(MEMSIZE(), true); //  matrix is still fine
        }
        memcpy(pBuffer, m_pMatrix, MEMSIZE());
    }
    m_n++;    
    if (!m_pMatrix) { delete m_pMatrix; }
    m_pMatrix = new TYPE[DIMENSION()];
    if (!m_pMatrix) {
        m_n--;
        delete pBuffer;
        throw eBadMatrixAlloc(MEMSIZE());
    }
    for (ushort i = 0; i < m_m; i++) {
        for (ushort j = 0; j < (m_n - 1); j++) {
            m_pMatrix[LINADDRESS(i, j)] = pBuffer[j + (i * (m_n - 1))];
        }
        m_pMatrix[LINADDRESS(i, m_n - 1)] = (TYPE) 0.0;
    }
    if (pBuffer) {
        delete pBuffer;
    }
    return (m_n - 1);
}

template <class TYPE> ushort matrix<TYPE>::AddRow(void) {    
    TYPE* pBuffer = 0;
    if (!EMPTY()) {
        pBuffer = new TYPE[DIMENSION()];
        if (!pBuffer) {
            throw eBadBufferAlloc(MEMSIZE(), true); //  matrix is still fine
        }
        memcpy(pBuffer, m_pMatrix, MEMSIZE());
    }
    m_m++;
    if (!m_pMatrix) {
        delete m_pMatrix;
    }
    m_pMatrix = new TYPE[DIMENSION()];
    if (!m_pMatrix) {
        m_m--;
        delete pBuffer;
        throw eBadMatrixAlloc(MEMSIZE());
    }
    for (ushort i = 0; i < (m_m - 1); i++) {
        for (ushort j = 0; j < m_n; j++) {
            m_pMatrix[LINADDRESS(i, j)] = pBuffer[LINADDRESS(i, j)];
        }
    }
    for (ushort j = 0; j < m_n; j++)
        m_pMatrix[LINADDRESS(m_m - 1, j)] = (TYPE) 0.0;
    if (pBuffer) {
        delete pBuffer;
    }
    return (m_m - 1);
}

template <class TYPE> ushort matrix<TYPE>::DeleteColumn(ushort c)
{
    if (!m_n) {
        throw eException("\nCan not delete one of the columns when there are none.\n");
    }
    if (__I(c) >= m_n) {
        throw eException();
    }
    TYPE* pBuffer = 0;
    if (!EMPTY()) {
        pBuffer = new TYPE[DIMENSION()];
        if (!pBuffer) {
            throw eBadBufferAlloc(MEMSIZE(), true);
        }
        memcpy(pBuffer, m_pMatrix, MEMSIZE()); 
    }
    m_n--;
    if (!EMPTY()) {
        if (m_pMatrix) {
            delete m_pMatrix;
        }
        m_pMatrix = new TYPE[DIMENSION()];
        if (!m_pMatrix) {
            m_n++;
            delete pBuffer;
            throw eBadMatrixAlloc(MEMSIZE());
        }
        short s;    //  shift
        for (ushort i = 0; i < m_m; i++) {
            s = 0;
            for (ushort j = 0; j < m_n; j++) {
                if (j == __I(c)) {
                    s = 1;
                }
                m_pMatrix[LINADDRESS(i, j)] = pBuffer[(j + s) + (i * (m_n + 1))];
            }
        }
    }
    if (pBuffer) {
        delete pBuffer;
    }
    return (m_n - 1);
}

template <class TYPE> ushort matrix<TYPE>::DeleteRow(ushort r) {
    if (!m_m) {
        throw eException("\nCan not delete one of the rows when there are none.\n");
    }
    if (__I(r) >= m_m) {
        throw eException();
    }  
    TYPE* pBuffer = 0;
    if (!EMPTY()) {
        pBuffer = new TYPE[DIMENSION()];
        if (!pBuffer) {
            throw eBadBufferAlloc(MEMSIZE(), true);
        }
        memcpy(pBuffer, m_pMatrix, MEMSIZE());
    }
    m_m--;
    if (!EMPTY()) {
        if (m_pMatrix) {
            delete m_pMatrix;
        }
        m_pMatrix = new TYPE[DIMENSION()];
        if (!m_pMatrix) {
            m_m++;
            delete pBuffer;
            throw eBadMatrixAlloc(MEMSIZE());
        }
        short s;    //  shift
        for (ushort i = 0; i < m_m; i++) {
            if (__I(r) == i) {
                s = 1;
            }
            for (ushort j = 0; j < m_n; j++) {
                m_pMatrix[LINADDRESS(i, j)] = pBuffer[j + ((i + s) * (m_n))];
            }
        }
    }
    if (pBuffer) {
        delete pBuffer;
    }
    return (m_n - 1);
}


//=============================================================================
//  misc
//-----------------------------------------------------------------------------

template <class TYPE> void matrix<TYPE>::Dump(bool wide /* = true */)
{
    cerr << endl << "address: " << (void*) m_pMatrix << endl;
    cerr << "memsize: " << MEMSIZE() << endl;
    cerr << "values:" << endl;
    if (m_pMatrix) {
        for (ushort i = 0; i < m_m; i++) {
            for (ushort j = 0; j < m_n; j++) {
                cerr << m_pMatrix[LINADDRESS(i, j)] << (wide ? '\t' : ' '); 
            }
            cerr << endl;
        }
    }
}


//  reset precompiler to previous state
//
#pragma optimize("", off)
#pragma warning(pop)


//-----------------------------------------------------------------------------
//  end of file
//=============================================================================
