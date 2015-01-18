////////////////////////////////////////////////////////////////////
//
// $Id$
//
// STL string utility
//
// Copyright (c) 2005 by Takashi Kanai. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _STRUTIL_H
#define _STRUTIL_H 1

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

class StrUtil {

public:

  // Å‰‚Ì•¶š‚ğŒ©‚Â‚¯‚é
  void first_word( std::string& str, std::string& fw ) {
    std::istringstream in(str); in >> fw;
  };

  // Å‰‚Ì•¶š‚ğŒ©‚Â‚¯‚é
  void nth_word( std::string& str, int n, std::string& fw ) {
    std::istringstream in(str);
    int i = 0;
    std::string t;
    while ( i < n ) { in >> t; ++i; }
    fw = t;
  };

  // •¶š—ñƒJƒEƒ“ƒg
  int word_count( std::string& str ) {
    std::istringstream in(str);

    int count = 0; 
    std::string sstr;
    while ( in >> sstr ) ++count;

    return count;
  };

  // ®”Œ^‚ğ•¶š—ñŒ^‚É•ÏŠ·
  std::string itos( int n ) {
    std::stringstream str_stream;
    str_stream << n;
    return str_stream.str();
  };

  // ®”Œ^‚ğ•¶š—ñŒ^‚É•ÏŠ·
  std::string ftos( float n ) {
    std::stringstream str_stream;
    str_stream << n;
    return str_stream.str();
  };
};

#endif // _STRUTIL_H

  
