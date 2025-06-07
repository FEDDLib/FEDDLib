#ifndef Helper_hpp
#define Helper_hpp

#include <iostream>
#include <sstream>
#include <iomanip>

//#include "AssembleFE_decl.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/General/SmallMatrix.hpp"

namespace FEDD {

/*! Helper class of static functions that contains rudimental finite element components.
 * It contains basis functions, quadrature rules, transformation functions, and other stuff.
 *
 *
 */
class Helper {

public:
    typedef double SC; // TODO: [JK] Does this not clash with the default definition in DefaultTypeDefs.hpp, which is used elsewhere?

    /// Order of derivative of a function.
    enum VarType {
        Deriv0 = 0, ///< order 0, f(x)
        Deriv1 = 1  ///< order 1, gradient(f(x))
    };

    /// @brief Compute surface normal of corresponding surface
    /// @param[in] dim Dimension
    /// @param[in] pointsRep List of all repeated nodes
    /// @param[in] nodeList Ids of local surface points
    /// @param[out] v_E normal vector
    /// @param[out] norm_v_E Normal vector length
    static void computeSurfaceNormal(int dim,
                                     vec2D_dbl_ptr_Type pointsRep,
                                     vec_int_Type nodeList,
                                     vec_dbl_Type &v_E,
                                     double &norm_v_E);


    /// @brief Build transformation of element to reference element depending on FEType
    /// @param[in] element Finite element
    /// @param[in] pointsRep List of repeated points
    /// @param[out] B Resulting transformation matrix
    /// @param[in] FEType FE Discretization 
    static void buildTransformation(const vec_int_Type& element,
                             vec2D_dbl_ptr_Type pointsRep,
                             SmallMatrix<SC>& B,
                             std::string FEType="P");

    /// @brief Build transformation of element to reference element depending on FEType
    /// @param[in] element Finite element
    /// @param[in] pointsRep List of repeated points
    /// @param[out] B Resulting transformation matrix
    /// @param[in] b Point to transform from
    /// @param[in] FEType FE Discretization 
    static void buildTransformation(const vec_int_Type& element,
                             vec2D_dbl_ptr_Type pointsRep,
                             SmallMatrix<SC>& B,
                             vec_dbl_Type& b,
                             std::string FEType="P");


    /// @brief Transformation of a surface to the reference element
    /// @param[in] element Finite element
    /// @param[in] pointsRep List of repeated points
    /// @param[out] B Resulting transformation matrix
    /// @param[in] b Point to transform from
    /// @param[in] FEType FE Discretization 
    static void buildTransformationSurface(const vec_int_Type& element,
                                    vec2D_dbl_ptr_Type pointsRep,
                                    SmallMatrix<SC>& B,
                                    vec_dbl_Type& b,
                                    std::string FEType="P");

    /// @brief Returning gradient of phi evaluated at the quadrature points
    /// @param[in] Dimension Dimension
    /// @param[in] intFE number corresponding to FE disc.
    /// @param[in] i basisfunction i
    /// @param[in] QuadPts quadpoints
    /// @param[out] value vector including values
    static void gradPhi(int Dimension,
                    int intFE,
                    int i,
                    vec_dbl_Type &QuadPts,
                    vec_dbl_ptr_Type &value);
    
    /// @brief Get quadrature formula.
    ///
    /// A quadrature formula is used to integrate approximate an integral of f(x) over the domain Omega, i.e., 
    ///    integral_Omega f(x) dx,
    /// via the following sum:
    ///    sum_i f(xi)*wi,
    /// where xi is a quadrature point, and wi a quadrature weight.
    ///
    /// Source of the formulas: Most of the quadrature formulas can be found in http://code-aster.org/doc/v11/en/man_r/r3/r3.01.01.pdf 01/2021
    ///
    /// @param[in] Dimension Space dimension d of the domain of f:IR^d-->IR.
    /// @param[in] Degree Order of polynomial that must be integrated exactly in exact arithmetic with the returned formula.
    /// @param[out] QuadPts Quadrature points
    /// @param[out] QuadW Quadrature weights
    /// @param[in] FEType Finite element type, e.g., "P1", "Q2" etc.
    static void getQuadratureValues(int Dimension,
                            int Degree,
                            vec2D_dbl_ptr_Type &QuadPts,
                            vec_dbl_ptr_Type &QuadW,
                            std::string FEType);

    /*!
    \brief Returns quadrature formula on surface element.
    
    Is distinguishes between needing Element or Surface information. !! Input can be improved with just delivering the coordinates of the surface nodes to determine the quad points
    
    Keep in mind that elementwise quadPoints are defined on reference element whereas surface quadPoints at hand are defined on the input surface, which is typically not the reference Element. 
    
    @param[in] dim Space dimension of the underlying domain, e.g., 3 for a triangle in IR^3..
    @param[in] FEType Finite element type of the underlying element (P1, P2, Q1 etc.), e.g., a "P2" tetrahedron for which a quadrature formula on one of its surface elements (triangles) is required.
    @param[out] QuadW Vector to be filled with the quadrature weights
    @param[in] vec_LO_Type surfaceIDs for which you need the quadrature points.
    @param[in] points The repeated(!) points of current problem to identify the surface node ids. 
    @returns Quadrature points
    */
    static vec2D_dbl_Type getQuadratureValuesOnSurface(int dim,
                                                       std::string FEType,
                                                       vec_dbl_Type &QuadW,
                                                       vec_LO_Type surfaceIDs,
                                                       vec2D_dbl_ptr_Type points);
    
    /// @brief Full matrix representation of gradient of a basis function for each quadrature point
    /// @param DPhi grad Phi per quadpoint dim:(quadpoint,i,j)
    /// @param weightsDPhi Quadrature weights
    /// @param Dimension Dimension
    /// @param FEType Finite Element Type
    /// @param Degree Integration degree
    /// @return 
    static int getDPhi(	vec3D_dbl_ptr_Type &DPhi,
                	vec_dbl_ptr_Type &weightsDPhi,
                    int Dimension,
                    std::string FEType,
                    int Degree);

    //  @brief Natalie new function to get viscosity at center of mass
    /// @param DPhi grad Phi p
    /// @param Dimension Dimension
    /// @param FEType Finite Element Type
    static void getDPhiAtCM(vec3D_dbl_ptr_Type &DPhi,
                            int dim,
                            std::string FEType);


    /// @brief Applying the transformation matriX B to the gradient of phi, as is done in when transforming the gradient of phi to the reference element
    /// @param dPhiIn 
    /// @param dPhiOut 
    /// @param Binv 
    static void applyBTinv( vec3D_dbl_ptr_Type& dPhiIn,
                    vec3D_dbl_Type& dPhiOut,
                    const SmallMatrix<SC>& Binv);


    /*!
    \brief Determine polynomial degree of a finite element basis function or its gradient that is required to select the correct quadrature formula for exact integration.

    \details
    This is a wrapper. 
    See Helper::requiredQuadratureDegreeForBasisfunction for details in case the degree is required for a basis function.
    See Helper::requiredQuadratureDegreeForGradientOfBasisfunction for details in case the degree is required for the gradient of a basis function.

    \param[in] dim  Dimension of the domain.
    \param[in] FEType  Finite element type, e.g., "P1", "Q2" etc.
    \param[in] orderOfDerivative  {Deriv0,Deriv1} Order of the derivative: order=0 (type=Deriv0) is the original function, order=1 (type=Deriv1) is the gradient, i.e., the first derivative.

    \return a polynomial degree for a finite element basis function or its gradient
    */
    static UN determineDegree(UN dim,
                              std::string FEType,
                              VarType orderOfDerivative);
                       

    /// @brief Get basisfunction phi per quadrature point
    /// @param Phi Basisfunction phi per quad point with (quadpoint,i)
    /// @param weightsPhi Quadrature weights
    /// @param dim dimension
    /// @param FEType Finite element discretization
    /// @param Degree Integration degree
    /// @param FETypeQuadPoints 
    /// @return 
    static int getPhi(vec2D_dbl_ptr_Type &Phi,
                            vec_dbl_ptr_Type &weightsPhi,
                            int dim,
                            std::string FEType,
                            int Degree,
               			    std::string FETypeQuadPoints="");

    static int getFuncAtQuadNodes(vec_dbl_ptr_Type &funcVals, RhsFunc_Type &rhsFunc, int dim, std::string FEType,
                                  int Degree, std::string FETypeQuadPoints = "");

    /// @brief Get phi i
    /// @param dim
    /// @param intFE
    /// @param i
    /// @param p
    /// @param value
    static void phi(int dim, int intFE, int i, vec_dbl_Type &p, double *value);

    static int getPhiGlobal(vec2D_dbl_ptr_Type &Phi,
                            vec_dbl_ptr_Type &weightsPhi,
                            int dim,
                            std::string FEType,
                            int Degree)
    { TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "getPhiGlobal not implemented yet.");};


private:
	
	Helper(){};

    /*!
    \brief Determine polynomial degree of a finite element basis function that is required to select the correct quadrature formula for exact integration.

    \details
    See Helper::requiredQuadratureDegreeForGradientOfBasisfunction for the analogous function to determine the degree for the gradient of a finite element basis function.

    For a finite element basis function: 
    For edges, triangles, and tetrahedra (i.e., for Pk elements), this function returns the polynomial degree (i.e., the maximum degree in any direction). This can be used to select the correct quadrature formula to integrate the basis function exactly.
    For edges, rectangles, and cuboids (i.e., for Qk elements), this function returns the polynomial degree on the reference element in x, y, or z direction. This can be used to select the correct tensor-product quadrature formula for exact integration of the basis function.

    Examples in 2D:
    a) FEType="P1" (triangle): phi_1(x,y) = 1-x-y is linear. 
       (All other basis functions are similar.)
       1 is returned.
    b) FEType="Q1" (quadrilateral): phi_1(x,y) = x*y is quadratic but linear in each coordinate direction.
       (All other basis functions are similar.)
       1 is returned.

    \param[in] dim  Dimension of the domain.
    \param[in] FEType  Finite element type, e.g., "P1", "Q2" etc.

    \return a polynomial degree for a finite element basis function
    */
    static UN requiredQuadratureDegreeForBasisfunction(UN dim, std::string FEType);

    /*!
    \brief Determine polynomial degree of the gradient of a finite element basis function that is required to select the correct quadrature formula for exact integration.

    \details
    See Helper::requiredQuadratureDegreeForBasisfunction for the analogous function to determine the degree for a finite element basis function.

    For the gradient of a finite element basis function: 
    For edges, triangles, and tetrahedra (i.e., for Pk elements), this function returns the polynomial degree (i.e., the maximum degree in any direction). This can be used to select the correct quadrature formula to integrate the gradient of a basis function exactly.
    For edges, rectangles, and cuboids (i.e., for Qk elements), this function returns the polynomial degree on the reference element in x, y, or z direction. This can be used to select the correct tensor-product quadrature formula for exact integration of the gradient of a basis function.

    Examples in 2D:
    a) FEType="P1" (triangle): phi_1(x,y) = 1-x-y is linear. 
       grad(phi_1(x,y)) = [-1,-1] is constant.
       (All other basis functions are similar.)
       0 is returned.
    b) FEType="Q1" (quadrilateral): phi_1(x,y) = x*y is quadratic but linear in each coordinate direction.
       grad(phi_1(x,y)) = [y,x] is also linear in each coordinate direction.
       (All other basis functions are similar.)
       1 is returned.

    \param[in] dim  Dimension of the domain.
    \param[in] FEType  Finite element type, e.g., "P1", "Q2" etc.

    \return a polynomial degree for the gradient of a finite element basis function
    */
    static UN requiredQuadratureDegreeForGradientOfBasisfunction(UN dim, std::string FEType);

};
}
#endif
