#ifndef BIOINFORMATIKA_PROJEKT_DEBUG_UTIL_HPP
#define BIOINFORMATIKA_PROJEKT_DEBUG_UTIL_HPP

#include <iostream>

#define DEBUG(x) do { std::cerr << x << std::endl; } while (0)

#define LOG_ERR(x) do { std::cerr << "ERROR: " << x << std::endl; } while(0)

#endif //BIOINFORMATIKA_PROJEKT_DEBUG_UTIL_HPP
