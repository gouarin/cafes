#include <array>
#include <type_traits>

template<std::size_t Dimensions, std::size_t Ndm=1>
struct context{
    using array1d = std::array<double, Dimensions>;
    using array2d = std::array<std::array<double, Dimensions>, Ndm>;
    std::conditional<Ndm==1, array1d, array2d> h;
};

int main()
{
    context<2> s1{{{.1, .1}}};
    context<2, 2> s2{ {{ {{.1, .1}}, {{.1, .1}} }} };
    return 0;
}