#ifndef _JSON_H_
#define _JSON_H_

#include "common.h"

// simple generic value for serialization and command line options
// modeled on the JSON model, but with builtin semantics for fast
// different number formats and arrays of basic types
struct jsonvalue {
    // typedefs
    typedef vector<jsonvalue> array;
    typedef map<string,jsonvalue> object;
    
    // possible types of jsonvalue
    enum _Type { nullt, boolt, doublet, stringt, arrayt, objectt };
    _Type _type = nullt;    // current type
    union {
        bool    _b; // bool value
        double  _d; // number value
        string* _s; // string value
        array*  _a; // generic array value
        object* _o; // object type
    };
    
    // constructor
    jsonvalue() : _type(nullt) { }
    
    // value constructors
    explicit jsonvalue(bool b) : _type(boolt), _b(b) { }
    explicit jsonvalue(int i) : _type(doublet), _d(i) { }
    explicit jsonvalue(double d) : _type(doublet), _d(d) { }
    explicit jsonvalue(const string& s) : _type(stringt), _s(new string(s)) { }
    explicit jsonvalue(const char* s) : _type(stringt), _s(new string(s)) { }
    explicit jsonvalue(const array& a) : _type(arrayt), _a(new vector<jsonvalue>(a)) { }
    explicit jsonvalue(const object& o) : _type(objectt), _o(new map<string,jsonvalue>(o)) { }
    
    // copy constructor
    jsonvalue(const jsonvalue& j) : _type(nullt) { set(j); }
    
    // destuctor
    ~jsonvalue() { _clear(); }
    
    // assignment
    jsonvalue& operator=(const jsonvalue& j) { set(j); return *this; }
    
    // clear
    void _clear() {
        if(_type==stringt) delete _s;
        if(_type==arrayt) delete _a;
        if(_type==objectt) delete _o;
        _type = nullt;
    }
    // set
    void set(const jsonvalue& j) {
        if(_type != nullt) _clear();
        _type = j._type;
        switch(_type) {
            case nullt: break;
            case boolt: _b = j._b; break;
            case doublet: _d = j._d; break;
            case stringt: _s = new string(*j._s); break;
            case arrayt: _a = new vector<jsonvalue>(*j._a); break;
            case objectt: _o = new map<string,jsonvalue>(*j._o); break;
            default: error("wrong type");
        }
    }
    
    // type checking
    bool is_null() const { return _type == nullt; }
    bool is_bool() const { return _type == boolt; }
    bool is_number() const { return _type == doublet; }
    bool is_string() const { return _type == stringt; }
    bool is_array() const { return _type == arrayt; }
    bool is_object() const { return _type == objectt; }
    
    // getters for values
    bool as_bool() const { error_if_not(is_bool(), "wrong type"); return _b; }
    int as_int() const { error_if_not(is_number(), "wrong type"); return (int)_d; }
    float as_float() const { error_if_not(is_number(), "wrong type"); return (float)_d; }
    double as_double() const { error_if_not(is_number(), "wrong type"); return _d; }
    string as_string() const { error_if_not(is_string(), "wrong type"); return *_s; }
    
    // getters for arrays and objects
    const vector<jsonvalue>& as_array_ref() const { error_if_not(is_array(), "wrong type"); return *_a; }
    const map<string,jsonvalue>& as_object_ref() const { error_if_not(is_object(), "wrong type"); return *_o; }

    // proprties of arrays and objects
    int array_size() const { return as_array_ref().size(); }
    const jsonvalue& array_element(int idx) const { error_if_not(idx >= 0 and idx < array_size(), "wrong element index"); return as_array_ref().at(idx); }
    bool object_contains(const string& name) const { return as_object_ref().find(name) != as_object_ref().end(); }
    const jsonvalue& object_element(const string& name) const { error_if_not(object_contains(name), "wrong element name"); return as_object_ref().find(name)->second; }
};

// json loading
jsonvalue load_json(const string& filename);

// command line specification
struct CommandLine {
    // description of command line argument
    struct Argument {
        string                  name;   // complete name
        string                  flag;   // flag string name
        string                  desc;   // description
        string                  type;   // supported: int, float, double, bool (flag), string
        bool                    opt;    // whether it is optional argument
        jsonvalue               def;    // default value
    };
    
    // parameters
    string              progname;   // program name
    string              progdesc;   // program description
    vector<Argument>    options;    // program options (flags)
    vector<Argument>    arguments;  // program arguments
    
    CommandLine(const string& progname, const string& progdesc,
                const vector<Argument>& options, const vector<Argument>& arguments) :
        progname(progname), progdesc(progdesc), options(options), arguments(arguments) { }
};

// parsing command line arguments to json
jsonvalue parse_cmdline(const string& args, const CommandLine& cmd);
jsonvalue parse_cmdline(int argc, char** argv, const CommandLine& cmd);

#endif
