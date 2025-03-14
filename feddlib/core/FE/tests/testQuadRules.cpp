// Test quadrature rules by integrating some functions and checking the result.

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/FE/Helper.hpp"

using namespace FEDD;

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

template <std::size_t dim, typename Lambda> int check_integration(const int degree, std::string FEType, const SC expected_result, const int test, Lambda &f, const char *file, const int line) {
    vec_dbl_ptr_Type w = Teuchos::rcp(new vec_dbl_Type(0));
    vec2D_dbl_ptr_Type p;

    Helper::getQuadratureValues(dim, degree, p, w, FEType);
    SC integral = 0.0;
    for (int i = 0; i < w->size(); i++) {
        if constexpr (dim == 1) {
            integral += w->at(i) * std::invoke(std::forward<Lambda>(f), p->at(i).at(0));
        } else if constexpr (dim == 2) {
            integral += w->at(i) * std::invoke(std::forward<Lambda>(f), p->at(i).at(0), p->at(i).at(1));
        } else if constexpr (dim == 3) {
            integral += w->at(i) * std::invoke(std::forward<Lambda>(f), p->at(i).at(0), p->at(i).at(1), p->at(i).at(2));
        } else {
            std::cout << "Test (" << test << ") " << "Dimension must be 1, 2, or 3. dim = " << dim << std::endl;
            return EXIT_FAILURE;
        }
    }

    SC error_result = fabs(integral - expected_result);
    if (error_result > std::numeric_limits<double>::epsilon() * 100.0) {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(2) << error_result;
        std::cout << "Test (" << test << ") " << "Integral does not match expected result: error = " << oss.str() << std::endl << "    " << file << ":" << line << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    int test = 0;

    /////////////// 1D Tests ///////////////

    // Test 1: 1D, polynomial order 0
    {
        auto f = [](SC x) -> SC { return 1.0; };
        SC r = 1.0; // expected result
        const int dim = 1;
        int degree = 0;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 2: 1D, polynomial order 1
    {
        auto f = [](SC x) -> SC { return 1 - x; };
        SC r = 0.5; // expected result
        const int dim = 1;
        int degree = 1;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 3: 1D, polynomial order 2
    {
        auto f = [](SC x) -> SC { return 1.0 + 2.0 * x - 0.5 * x * x; };
        SC r = 1.0 + 1.0 / 1.2; // expected result
        const int dim = 1;
        int degree = 2;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 4: 1D, polynomial order 3
    {
        auto f = [](SC x) -> SC { return -1.0 + x - 4 * x * x * x; };
        SC r = -1.5; // expected result
        const int dim = 1;
        int degree = 3;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 5: 1D, polynomial order 4
    {
        auto f = [](SC x) -> SC { return -2.0 - x * x + 0.5 * x * x * x * x; };
        SC r = -67.0 / 30.0; // expected result
        const int dim = 1;
        int degree = 4;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 6: 1D, polynomial order 5
    {
        auto f = [](SC x) -> SC { return x - 4 * x * x * x + pow(x, 5); };
        SC r = -1.0 / 3.0; // expected result
        const int dim = 1;
        int degree = 5;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 7: 1D, polynomial order 6
    {
        auto f = [](SC x) -> SC { return x * x * x + 5 * x * x * x * x - 14 * pow(x, 6); };
        SC r = -3.0 / 4.0; // expected result
        const int dim = 1;
        int degree = 6;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 8: 1D, polynomial order 7
    {
        auto f = [](SC x) -> SC { return 0.5 - x + x * x - 4.0 / 3.0 * x * x * x + x * x * x * x - 3.0 * pow(x, 5) - 14.0 * pow(x, 6) + 4.0 * pow(x, 7); };
        SC r = -1.8; // expected result
        const int dim = 1;
        int degree = 7;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }


    /////////////// 2D Tests ///////////////

    // Test 9: 2D, triangle, polynomial order 0
    {
        auto f = [](SC x, SC y) -> SC { return 1.0; };
        SC r = 0.5; // expected result
        const int dim = 2;
        int degree = 0;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 10: 2D, triangle, polynomial order 1
    {
        auto f = [](SC x, SC y) -> SC { return 3.0 + 2.0*x + y; };
        SC r = 2.0; // expected result
        const int dim = 2;
        int degree = 1;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 11: 2D, triangle, polynomial order 2
    {
        auto f = [](SC x, SC y) -> SC { return -4.0 + x - y + 2.0*x*x + 4.0*y*y + 12.0*x*y; };
        SC r = -1.0; // expected result
        const int dim = 2;
        int degree = 2;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 12: 2D, triangle, polynomial order 3
    {
        auto f = [](SC x, SC y) -> SC { return 0.5 + 2.0*x + 3.0*y - 2.0*x*x - 3.0*y*y + 6.0*x*y + x*x*x - 2.0*y*y*y + x*x*y - 2.0*x*y*y; };
        SC r = 0.85; // expected result
        const int dim = 2;
        int degree = 3;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 13: 2D, triangle, polynomial order 4
    {
        auto f = [](SC x, SC y) -> SC { return -5.0 - 0.5*x + y + x*x + 0.5*y*y - x*y - 0.5*x*x*y + x*y*y + 1.5*x*x*x - y*y*y + x*x*x*x - x*x*x*y + 3.0*x*x*y*y - x*y*y*y + 1.5*y*y*y*y; };
        SC r = -2.21-2.0/300.0; // expected result
        const int dim = 2;
        int degree = 4;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 14: 2D, triangle, polynomial order 5
    {
        auto f = [](SC x, SC y) -> SC { return 1.0 - x + 2.0*y - 2.5*x*x + 1.5*y*y + x*y + 4.0*x*x*y + 2.0*x*y*y - 3.0*x*x*x + 0.5*y*y*y - 1.5*x*x*x*x - 2.0*x*x*x*y + 2.0*x*x*y*y + 0.5*x*y*y*y + y*y*y*y - 42.0*x*x*x*x*x - 21.0*x*x*x*x*y + -21.0*x*x*x*y*y + 42.0*x*x*y*y*y - 63.0*x*y*y*y*y + 21.0*y*y*y*y*y; };
        SC r = -5.0/90.0 - 0.2125; // expected result
        const int dim = 2;
        int degree = 5;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

//    // Test 15: 2D, triangle, polynomial order 6
//    {
//        auto f = [](SC x, SC y) -> SC { return ; };
//        SC r = -5.0/90.0 - 0.2125; // expected result
//        const int dim = 2;
//        int degree = 5;
//        std::string FEType = "P";
//        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
//            return EXIT_FAILURE;
//    }
//
//    // Test 16: 2D, triangle, polynomial order 7
//    {
//        auto f = [](SC x, SC y) -> SC { return ; };
//        SC r = -5.0/90.0 - 0.2125; // expected result
//        const int dim = 2;
//        int degree = 5;
//        std::string FEType = "P";
//        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
//            return EXIT_FAILURE;
//    }


    /////////////// 3D Tests ///////////////

    // Test 17: 3D, tetrahedron, polynomial order 0
    {
        auto f = [](SC x, SC y, SC z) -> SC { return 1.0; };
        SC r = 1.0/6.0; // expected result
        const int dim = 3;
        int degree = 0;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 18: 3D, tetrahedron, polynomial order 1
    {
        auto f = [](SC x, SC y, SC z) -> SC { return 0.5 - x + 2.0*y + 3.0*z; };
        SC r = 0.25; // expected result
        const int dim = 3;
        int degree = 1;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    // Test 19: 3D, tetrahedron, polynomial order 2
    {
        auto f = [](SC x, SC y, SC z) -> SC { return -0.5 + 2.0*x + y + z + 5.0*x*x + 4.0*y*y + 3.0*z*z + 2.0*x*y + x*z - y*z; };
        SC r = 0.3; // expected result
        const int dim = 3;
        int degree = 2;
        std::string FEType = "P";
        if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE)
            return EXIT_FAILURE;
    }

    std::cout << "testQuadRules: All " << test << " tests passed." << std::endl;

    return (EXIT_SUCCESS);
}
