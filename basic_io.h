#include <iostream>

// print funcs ###############################################################
/* My print   works like python */
template <typename T>
void print(T t)
{
    std::cout << t << std::endl ;
}

template<typename T, typename... Args>
void print(T t, Args... args)  // recursive variadic function
{
    std::cout << t ;

    print(args...) ;
}
// ############################################################################

