#include "json.h"
#include "picojson.h"

// json value conversion from parser
static jsonvalue _to_jsonvalue(const picojson::value& pjson) {
    if(pjson.is<picojson::null>()) return jsonvalue();
    else if(pjson.is<bool>()) return jsonvalue( pjson.get<bool>() );
    else if(pjson.is<double>()) return jsonvalue( pjson.get<double>() );
    else if(pjson.is<string>()) return jsonvalue( pjson.get<string>() );
    else if(pjson.is<picojson::array>()) {
        auto json = jsonvalue::array();
        for(auto& j : pjson.get<picojson::array>()) json.push_back(_to_jsonvalue(j));
        return jsonvalue(json);
    }
    else if(pjson.is<picojson::object>()) {
        auto json = jsonvalue::object();
        for(auto& nj : pjson.get<picojson::object>()) json[nj.first] = _to_jsonvalue(nj.second);
        return jsonvalue(json);
    } else { error("unknown type"); return jsonvalue(); }
}

// json handling
jsonvalue load_json(const string& filename) {
    // open file
    auto stream = std::ifstream(filename);
    error_if_not(bool(stream), "cannot open file: %s", filename.c_str());
    // read json
    picojson::value pjson;
    stream >> pjson;
    auto err = picojson::get_last_error();
    error_if_not(err.empty(), "json reading error: %s", err.c_str());
    stream.close();
    // conversion
    auto json = _to_jsonvalue(pjson);
    // done
    return json;
}

// print usage information
static void _cmdline_print_usage(const CommandLine& cmd) {
    auto usage = "usage: " + cmd.progname;
    for(auto opt : cmd.options) {
        usage += " ";
        auto optname = (opt.flag == "") ? "--"+opt.name : "-"+opt.flag;
        auto optval = (opt.type != "bool") ? " <"+opt.name+">" : "";
        if(opt.opt) usage += "[" + optname + optval + "]";
        else usage += optname + optval;
    }
    for(auto arg : cmd.arguments) {
        usage += " ";
        if(arg.opt) usage += "[" + arg.name + "]";
        else usage += arg.name;
    }
    usage += "\n\n";
    if(not cmd.options.empty() or not cmd.arguments.empty()) {
        usage += "options:\n";
        for(auto opt : cmd.options) {
            usage += "    " + ( (opt.flag == "")?"":"-"+opt.flag+"/" ) +
                     "--" + opt.name + ( (opt.type=="bool")?"":" "+opt.name) + "\n";
            usage += "        " + opt.desc + "\n";
        }
        for(auto arg : cmd.arguments) {
            usage += "    " + arg.name + ( (arg.type=="bool")?"":" "+arg.name) + "\n";
            usage += "        " + arg.desc + "\n";
        }
        usage += "\n\n";
    }
    message(usage.c_str());
}

// parse error
static void _cmdline_parse_error(const string& msg, const CommandLine& cmd) {
    _cmdline_print_usage(cmd);
    message("%s\n\n",msg.c_str());
    exit(EXIT_FAILURE);
}

// parse value
static jsonvalue _cmdline_parse_value(const string& str, const CommandLine::Argument& opt, const CommandLine& cmd) {
    try {
        if("int" == opt.type) {
            return jsonvalue( std::stoi(str) );
        } else if("float" == opt.type) {
            return jsonvalue( std::stod(str) );
        } else if("double" == opt.type) {
            return jsonvalue( std::stod(str) );
        } else if("string" == opt.type) {
            return jsonvalue( str );
        } else error("unknown type");
    } catch (std::exception& e) {
        _cmdline_parse_error(tostring("invalid value for %s", opt.name.c_str()),cmd);
    }
    return jsonvalue();
}

// parsing values
jsonvalue parse_cmdline(const vector<string>& args, const CommandLine& cmd) {
    auto parsed = jsonvalue::object();
    auto largs = args; // make a copy to change
    set<string> visited;
    for(auto opt : cmd.options) {
        auto pos = -1;
        for(auto i : range(largs.size())) {
            if (largs[i] == "--"+opt.name or largs[i] == "-"+opt.flag) { pos = i; break; }
        }
        if(pos < 0) {
            if(opt.opt) parsed[opt.name] = opt.def;
            else _cmdline_parse_error(tostring("required option -%s",opt.flag.c_str()),cmd);
        } else {
            if("bool" != opt.type) {
                if(pos == largs.size()-1) _cmdline_parse_error(tostring("no value for argument %s",opt.name.c_str()), cmd);
                auto sval = largs[pos+1];
                largs.erase(largs.begin()+pos,largs.begin()+pos+2);
                parsed[opt.name] = _cmdline_parse_value(sval,opt,cmd);
            } else {
                parsed[opt.name] = jsonvalue(true);
                largs.erase(largs.begin()+pos);
            }
        }
    }
    for(auto arg : largs) {
        if(arg[0] == '-') _cmdline_parse_error(tostring("unknown option %s", arg.c_str()),cmd);
    }
    for(auto arg : cmd.arguments) {
        if(largs.empty()) {
            if(arg.opt) parsed[arg.name] = arg.def;
            else _cmdline_parse_error(tostring("missing required argument %s",arg.name.c_str()),cmd);
        } else {
            parsed[arg.name] = _cmdline_parse_value(largs[0], arg, cmd);
            largs.erase(largs.begin());
        }
    }
    if(not largs.empty()) _cmdline_parse_error("too many arguments",cmd);
    return jsonvalue(parsed);
}

// parsing values
jsonvalue parse_cmdline(int argc, char** argv, const CommandLine& cmd) {
    auto args = vector<string>();
    for(auto i : range(1,argc)) args.push_back(argv[i]);
    return parse_cmdline(args,cmd);
}

