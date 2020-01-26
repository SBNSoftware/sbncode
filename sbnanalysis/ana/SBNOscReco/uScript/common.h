#ifndef uscript_common_h
#define uscript_common_h

#include <string_view>

// Enable to print out debug information in VM
// #define DEBUG_TRACE_EXECUTION

// Enable to print out debug information in compiler
// #define DEBUG_PRINT_CODE

// Enable to print timing information at end of program
#define USCRIPT_TIMING


// Get the name of a type
// From: https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c
template <typename T>
constexpr auto type_name()
{
    std::string_view name, prefix, suffix;
#ifdef __clang__
    name = __PRETTY_FUNCTION__;
    prefix = "auto type_name() [T = ";
    suffix = "]";
#elif defined(__GNUC__)
    name = __PRETTY_FUNCTION__;
    prefix = "constexpr auto type_name() [with T = ";
    suffix = "]";
#elif defined(_MSC_VER)
    name = __FUNCSIG__;
    prefix = "auto __cdecl type_name<";
    suffix = ">(void)";
#endif
    name.remove_prefix(prefix.size());
    name.remove_suffix(suffix.size());
    return name;
}  

#ifdef USCRIPT_TIMING
#include <atomic>
#include <iostream>
class uscriptTimer: public std::atomic<unsigned> {
  std::string predicate;
public:
  uscriptTimer(const char *p): predicate(p) {}
  ~uscriptTimer() {
    std::cout << predicate << (load()/1.e6) << std::endl;
  }
};

extern uscriptTimer __uscript_global_compiletime;
extern uscriptTimer __uscript_global_vmtime;
#endif

#endif
