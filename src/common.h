#ifndef _COMMON_H_
#define _COMMON_H_

#include <string>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <stdio.h>
#include <cstdarg>
#include <fstream>
#include <cstdio>

// bringing stand libraray objects in scope
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::make_pair;
using std::set;

// Prints a message (printf style) and flushes the output (useful during debugging)
inline void message(const char* msg, ...) { va_list args; va_start(args, msg); vprintf(msg, args); va_end(args); fflush(stdout); }
// Prints an error message (printf style), flushes the output (useful during debugging) and asserts
inline void error(const char* msg, ...) { va_list args; va_start(args, msg); vprintf(msg, args); va_end(args); assert(false); fflush(stdout); }
// Checks for error, prints an error message (printf style), flushes the output (useful during debugging) and asserts
inline void error_if_not(bool check, const char* msg, ...) { if(check) return; va_list args; va_start(args, msg); vprintf(msg, args); va_end(args); assert(false); fflush(stdout); }
// simple string creation (printf style) --- warning: out should be < 4096
inline string tostring(const char* msg, ...) { char buf[4096]; va_list args; va_start(args, msg); vsprintf(buf, msg, args); va_end(args); return string(buf); }

// Python-style range: iterates from min to max in range-based for loops
struct range {
    int min = 0; // start
    int max = 1; // end
    
    // iterator
    struct iterator {
        int current; // current value
        iterator(int value) : current(value) { } // constructor
        int& operator*() { return current; } // iterator support
        iterator& operator++() { ++current; return *this; } // iterator support
        bool operator!=(const iterator& other) { return current != other.current; } // iterator support
    };
    
    // constructor for range in [0,max)
    range(int max) : max(max) { }
    // constructor for range in [min,max)
    range(int min, int max) : min(min), max(max) { }
    
    // iterator support
    iterator begin() { return iterator(min); }
    // iterator support
    iterator end() { return iterator(max); }
};

// load a text file into a buffer
inline string load_text_file(const char* filename) {
    auto text = string("");
    auto f = fopen(filename,"r");
    error_if_not(f, "could not open file: ", filename);
    char line[4096];
    while (fgets(line, 4096, f)) text += line;
    fclose(f);
    return text;
}

#endif
