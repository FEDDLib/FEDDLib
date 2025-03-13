// Test quadrature rules by integrating some functions and checking the result.

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>


#include "feddlib/core/General/DefaultTypeDefs.hpp"


#include "feddlib/core/FE/FE.hpp"



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
    /*
        {
            auto f = [](SC x, SC y) -> SC { return x*y; };
            SC r = 1.0; // expected result
            const int dim = 2;
            int degree = 0;
            std::string FEType = "P";
            if (check_integration<dim>(degree, FEType, r, ++test, f, __FILE__, __LINE__) == EXIT_FAILURE) return EXIT_FAILURE;
        }
    */

    std::cout << "testQuadRules: All " << test << " tests passed." << std::endl;

    return (EXIT_SUCCESS);
}
