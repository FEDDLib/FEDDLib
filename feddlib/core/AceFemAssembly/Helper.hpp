#ifndef Helper_hpp
#define Helper_hpp

//#include "AssembleFE_decl.hpp"
#include "feddlib/core/FEDDCore.hpp"

namespace FEDD {

class Helper {
  
public:

    enum VarType {Std=0,Grad=1};

	static void gradPhi(	int Dimension,
                    int intFE,
                    int i,
                    vec_dbl_Type &QuadPts,
                    vec_dbl_ptr_Type &value);
    
    /*! Most of the quadrature formulas can be found in http://code-aster.org/doc/v11/en/man_r/r3/r3.01.01.pdf 01/2021  */
    static void getQuadratureValues(int Dimension,
                            int Degree,
                            vec2D_dbl_ptr_Type &QuadPts,
                            vec_dbl_ptr_Type &QuadW,
                            std::string FEType);
    
    static int getDPhi(	vec3D_dbl_ptr_Type &DPhi,
                	vec_dbl_ptr_Type &weightsDPhi,
                    int Dimension,
                    std::string FEType,
                    int Degree);

    static UN determineDegree(UN dim,
                       std::string FEType,
                       UN degFunc);

	static UN determineDegree2(UN dim, 
								std::string FEType1, 		
								std::string FEType2, 
								VarType type1,
								VarType type2, 
								UN extraDeg = 0);

    static int getPhi(vec2D_dbl_ptr_Type &Phi,
                            vec_dbl_ptr_Type &weightsPhi,
                            int dim,
                            std::string FEType,
                            int Degree,
               			    std::string FETypeQuadPoints="");

	static void phi(int dim,
			  int intFE,
			  int i,
			  vec_dbl_Type &p,
			  double* value);


private:
	
	Helper(){};

};
}
#endif
