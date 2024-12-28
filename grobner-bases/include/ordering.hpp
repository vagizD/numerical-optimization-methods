#pragma once

#include "monomial.hpp"


namespace ADAAI {

enum TermOrdering {
    None,
    L,
    IL,
    TD,
    TDIL,
    BTDIL
};

template<typename F, size_t N>
bool Lexicographical(const Term<F, N>& t1, const Term<F, N>& t2);

}  // namespace ADAAI


