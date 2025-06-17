#include "Helper.hpp"

#include <string>

namespace FEDD {

UN Helper::requiredQuadratureDegreeForBasisfunction(UN dim, std::string FEType){
    UN deg;
    if (!FEType.compare("P0"))
        deg = 0;
    else if ( !FEType.compare("P1") || !FEType.compare("P1-disc") )
        deg = 1;
    else if (!FEType.compare("P2"))
        deg = 2;
    else if (!FEType.compare("Q1"))
        deg = 1;
    else if (!FEType.compare("Q2"))
        deg = 2;
    else if ( (dim == 3) && (!FEType.compare("Q2-20")) ) // Q2 serendipity
        deg = 2;
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown finite element type: " + FEType + ", dimension: " + std::to_string(dim) + ".");

    // TODO: [JK] 2025/04 Removed P2-CR element, since it is unclear what type of element it really is.
    //                    P2-CR sounds like polynomial 2 Crouzeix-Raviart. 
    //                    3D: I assume that this should have a degree of freedom on each edge of the tetrahedron and on each face, 10 in total.
    //                    This fits to the 10 DOFs mentioned in the Mesh class.
    //                    But this would use standard P2 basis functions on the element itself. 
    //                    It spans the standard P2 space and only uses different nodes, where the basis functions assume the value 1 and 0, respectively.
    //                    Thus, I don't understand why a polynomial degree of 4 is returned below.
    //                    With a little more investigation, we can re-add this element. Probably, the degree below should be 2 for the 
    //                    basis function and 1 for the gradient.
    //if ( (dim == 3) && (!FEType1.compare("P2-CR")) ) {
        // if (type == Helper::Deriv0)
        //    deg = 4; // [JK] Was ist das fuer ein Element, das in einem Tetraeder mit einer P2-Formel Ordnung 4 rausbekommt?!
        // else if (type == Helper::Deriv1)
        //    deg = 3;
    //}

    return deg;
}

UN Helper::requiredQuadratureDegreeForGradientOfBasisfunction(UN dim, std::string FEType){
    // In 1D, Qk finite elements are the same as Pk finite elements.
    UN deg;
    if (!FEType.compare("P0"))
        deg = 0;
    else if ( !FEType.compare("P1") || !FEType.compare("P1-disc") )
        deg = 0;
    else if (!FEType.compare("P2"))
        deg = 1;
    else if (!FEType.compare("Q1"))
        deg = (dim > 1 ? 1 : 0);
        // Example: f(x,y) = x*y.
        // Gradient(f(x,y)) = [y,x]
        // The required degree in x direction of the first argument is 0 but 1 for the second argument. 
        // The same holds for the y direction. Thus, we need degree 1 in each (x and y) direction.
    else if (!FEType.compare("Q2"))
        deg = (dim > 1 ? 2 : 1);
    else if ( (dim == 3) && (!FEType.compare("Q2-20")) ) // Q2 serendipity
        deg = 2;
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown finite element type: " + FEType + ", dimension: " + std::to_string(dim) + ".");

    return deg;
}

UN Helper::determineDegree(UN dim, std::string FEType, VarType orderOfDerivative){
    UN deg;
    if (orderOfDerivative == Deriv0)        // Deriv0  = 0 = no derivative
        deg = requiredQuadratureDegreeForBasisfunction(dim,FEType);
    else if (orderOfDerivative == Deriv1)  // Deriv1 = 1 = first derivative = gradient
        deg = requiredQuadratureDegreeForGradientOfBasisfunction(dim,FEType);
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown order of derivative: " + orderOfDerivative);
    return deg;
}

void Helper::buildTransformationSurface(const vec_int_Type& element,
                                                 vec2D_dbl_ptr_Type pointsRep,
                                                 SmallMatrix<SC>& B,
                                                 vec_dbl_Type& b,
                                                 std::string FEType){
    // small matrix always square
    TEUCHOS_TEST_FOR_EXCEPTION( (B.size()<2 || B.size()>3), std::logic_error, "Initialize SmallMatrix for transformation.");
    UN index;
    UN index0 = element.at(0);
    if (FEType[0]=='P') {
        for (UN j=0; j<B.size()-1; j++) {
            index = element.at(j+1);
            for (UN i=0; i<B.size(); i++) { // dimension
                B[i][j] = pointsRep->at(index).at(i) - pointsRep->at(index0).at(i);
            }
        }
        for (UN i=0; i<B.size(); i++)
            b[i] = pointsRep->at(index0).at(i);
    }
    else if (FEType[0]=='Q'){
        TEUCHOS_TEST_FOR_EXCEPTION( B.size()!=3, std::logic_error, "No Transformation for surface integrals.");
    }
}

void Helper::computeSurfaceNormal(int dim,
                                   vec2D_dbl_ptr_Type pointsRep,
                                   vec_int_Type nodeList,
                                   vec_dbl_Type &v_E,
                                   double &norm_v_E)
{

    vec_dbl_Type p1(dim),p2(dim);

    if(dim==2){
        v_E[0] = pointsRep->at(nodeList[0]).at(1) - pointsRep->at(nodeList[1]).at(1);
        v_E[1] = -(pointsRep->at(nodeList[0]).at(0) - pointsRep->at(nodeList[1]).at(0));
        norm_v_E = std::sqrt(std::pow(v_E[0],2)+std::pow(v_E[1],2));

    }
    else if(dim==3){

        p1[0] = pointsRep->at(nodeList[0]).at(0) - pointsRep->at(nodeList[1]).at(0);
        p1[1] = pointsRep->at(nodeList[0]).at(1) - pointsRep->at(nodeList[1]).at(1);
        p1[2] = pointsRep->at(nodeList[0]).at(2) - pointsRep->at(nodeList[1]).at(2);

        p2[0] = pointsRep->at(nodeList[0]).at(0) - pointsRep->at(nodeList[2]).at(0);
        p2[1] = pointsRep->at(nodeList[0]).at(1) - pointsRep->at(nodeList[2]).at(1);
        p2[2] = pointsRep->at(nodeList[0]).at(2) - pointsRep->at(nodeList[2]).at(2);

        v_E[0] = p1[1]*p2[2] - p1[2]*p2[1];
        v_E[1] = p1[2]*p2[0] - p1[0]*p2[2];
        v_E[2] = p1[0]*p2[1] - p1[1]*p2[0];

        norm_v_E = std::sqrt(std::pow(v_E[0],2)+std::pow(v_E[1],2)+std::pow(v_E[2],2));

    }

}

void Helper::buildTransformation(const vec_int_Type& element,
                                          vec2D_dbl_ptr_Type pointsRep,
                                          SmallMatrix<SC>& B,
                                          std::string FEType){

    TEUCHOS_TEST_FOR_EXCEPTION( (B.size()<2 || B.size()>3), std::logic_error, "Initialize SmallMatrix for transformation.");
    UN index;
    UN index0 = element.at(0);
    if (FEType[0]=='P') {
        for (UN j=0; j<B.size(); j++) {
            index = element.at(j+1);
            for (UN i=0; i<B.size(); i++) {
                B[i][j] = pointsRep->at(index).at(i) - pointsRep->at(index0).at(i);
            }
        }
    }
    else if (FEType[0]=='Q'){
        TEUCHOS_TEST_FOR_EXCEPTION( B.size()!=3, std::logic_error, "Transformation for quadrilateral elements only in 3D.");
        std::vector<int> indexVec(3);
        indexVec[0] = element[1]; indexVec[1] = element[3]; indexVec[2] = element[4];
        for (UN j=0; j<B.size(); j++) {
            for (UN i=0; i<B.size(); i++) {
                B[i][j] = ( pointsRep->at( indexVec[j] ).at(i) - pointsRep->at( index0 ).at(i) ) / 2.;
            }
        }
    }
}
void Helper::buildTransformation(const vec_int_Type& element,
                                          vec2D_dbl_ptr_Type pointsRep,
                                          SmallMatrix<SC>& B,
                                          vec_dbl_Type& b,
                                          std::string FEType){

    TEUCHOS_TEST_FOR_EXCEPTION( (B.size()<2 || B.size()>3), std::logic_error, "Initialize SmallMatrix for transformation.");
    UN index;
    UN index0 = element.at(0);
    if (FEType[0]=='P') {
        for (UN j=0; j<B.size(); j++) {
            index = element.at(j+1);
            for (UN i=0; i<B.size(); i++) {
                B[i][j] = pointsRep->at(index).at(i) - pointsRep->at(index0).at(i);
            }
        }
        for (UN i=0; i<B.size(); i++)
            b[i] = pointsRep->at(index0).at(i);
    }
    else if (FEType[0]=='Q'){
        TEUCHOS_TEST_FOR_EXCEPTION( B.size()!=3, std::logic_error, "Transformation for quadrilateral elements only in 3D.");
        std::vector<int> indexVec(3);
        indexVec[0] = element[1]; indexVec[1] = element[3]; indexVec[2] = element[4];
        for (UN j=0; j<B.size(); j++) {
            for (UN i=0; i<B.size(); i++) {
                B[i][j] = ( pointsRep->at( indexVec[j] ).at(i) - pointsRep->at( index0 ).at(i) ) / 2.;
            }
        }
        for (UN i=0; i<B.size(); i++)
            b[i] = pointsRep->at(index0).at(i);

    }
}


void Helper::applyBTinv( vec3D_dbl_ptr_Type& dPhiIn,
                                    vec3D_dbl_Type& dPhiOut,
                                    const SmallMatrix<SC>& Binv){
    UN dim = Binv.size();
    for (UN w=0; w<dPhiIn->size(); w++){
        for (UN i=0; i < dPhiIn->at(w).size(); i++) {
            for (UN d1=0; d1<dim; d1++) {
                for (UN d2=0; d2<dim; d2++) {
                    dPhiOut[w][i][d1] += dPhiIn->at(w).at(i).at(d2) * Binv[d2][d1];
                }
            }
        }
    }
}

int Helper::getDPhi(vec3D_dbl_ptr_Type &DPhi,
                     vec_dbl_ptr_Type &weightsDPhi,
                     int dim,
		             std::string FEType,
		             int Degree){

    int 			nmbLocElPts;
    int 			intFE;
    vec_dbl_ptr_Type 	value(new vec_dbl_Type(dim,0.0));
    vec2D_dbl_ptr_Type	QuadPts;

    if (dim==2) {
        getQuadratureValues(dim, Degree, QuadPts, weightsDPhi, FEType);
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 3;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 6;
            intFE = 2;
        }

        DPhi.reset(new vec3D_dbl_Type(weightsDPhi->size(),vec2D_dbl_Type(nmbLocElPts,vec_dbl_Type(2,0.0))));

        for (int k=0; k<DPhi->size(); k++ ){
            for (int i=0; i<DPhi->at(0).size(); i++) {
                gradPhi(dim,intFE,i,QuadPts->at(k),value);
                for (int j=0; j<2; j++) {
                    DPhi->at(k).at(i).at(j) = value->at(j);
                }
            }
        }
    }

    else if(dim==3){
    	getQuadratureValues(dim, Degree, QuadPts, weightsDPhi, FEType);
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 4;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 10;
            intFE = 2;
        }
        else if (FEType == "Q1") {
            nmbLocElPts = 8;
            intFE = 3;
        }
        else if (FEType == "Q2") {
            nmbLocElPts = 27;
            intFE = 4;
        }
        else if (FEType == "Q2-20") {
            nmbLocElPts = 20;
            intFE = 5;
        }
        else if (FEType == "P1-disc") {
            nmbLocElPts = 4;
            intFE = 6;
        }
        else if (FEType == "P1-disc-global")
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error   ,"grad of P1-disc-global not implemented yet.");

        DPhi.reset( new vec3D_dbl_Type( weightsDPhi->size(), vec2D_dbl_Type( nmbLocElPts, vec_dbl_Type(3,0.0) ) ) );
        for (int k=0; k<DPhi->size(); k++ ){
            for (int i=0; i<DPhi->at(0).size(); i++) {
                gradPhi(dim,intFE,i,QuadPts->at(k),value);
                for (int j=0; j<3; j++) {
                    DPhi->at(k).at(i).at(j) = value->at(j);
                }
            }
        }
    }

    return intFE;
}


void Helper::gradPhi(int dim,
                int intFE,
                int i,
                vec_dbl_Type &p,
                vec_dbl_ptr_Type &value){
    if (dim==2) {
        switch (intFE) {
            case 0://P0
                switch (i) {
                    case 0:
                        value->at(0)= 0.;
                        value->at(1)= 0.;
                        break;
                }
                break;
            case 1://P1
                switch (i) {
                    case 0:
                        value->at(0)= -1.;
                        value->at(1)= -1.;
                        break;
                    case 1:
                        value->at(0)= 1.;
                        value->at(1)= 0.;
                        break;
                    case 2:
                        value->at(0)= 0.;
                        value->at(1)= 1.;
                        break;
                }
                break;
            case 2://P2
                switch (i) {
                    case 0:
                        value->at(0)= 1. - 4.*(1 - p[0] - p[1]);
                        value->at(1)= 1. - 4.*(1 - p[0] - p[1]);
                        break;
                    case 1:
                        value->at(0)= 4.*p[0] - 1;
                        value->at(1)= 0.;
                        break;
                    case 2:
                        value->at(0)= 0.;
                        value->at(1)= 4.*p[1] - 1;
                        break;
                    case 3:
                        value->at(0)= 4 * (1. - 2*p[0] - p[1]);
                        value->at(1)= -4 * p[0];
                        break;
                    case 4:
                        value->at(0)= 4.*p[1];
                        value->at(1)= 4.*p[0];
                        break;
                    case 5:
                        value->at(0)= - 4.*p[1];
                        value->at(1)= 4 * (1. - p[0] - 2*p[1]);
                        break;
                }
                break;
        }
    }
    else if(dim==3) {
        switch (intFE) {
            case 0://P0
                switch (i) {
                    case 0:
                    value->at(0)= 0.;
                    value->at(1)= 0.;
                    value->at(2)= 0.;
                    break;
                }
                break;
            case 1://P1
                switch (i) {
                    case 0:
                        value->at(0)= -1.;
                        value->at(1)= -1.;
                        value->at(2)= -1.;
                        break;
                    case 1:
                        value->at(0)= 1.;
                        value->at(1)= 0.;
                        value->at(2)= 0.;
                        break;
                    case 2:
                        value->at(0)= 0.;
                        value->at(1)= 1.;
                        value->at(2)= 0.;
                        break;
                    case 3:
                        value->at(0)= 0.;
                        value->at(1)= 0.;
                        value->at(2)= 1.;
                        break;
                }
                break;
            case 2://P2
                switch (i) {
                    case 0:
                        value->at(0)= -3. + 4.*p[0] + 4.*p[1] + 4.*p[2];
                        value->at(1)= -3. + 4.*p[0] + 4.*p[1] + 4.*p[2];
                        value->at(2)= -3. + 4.*p[0] + 4.*p[1] + 4.*p[2];
                        break;
                    case 1:
                        value->at(0)= 4.*p[0] - 1;
                        value->at(1)= 0.;
                        value->at(2)= 0.;
                        break;
                    case 2:
                        value->at(0)= 0.;
                        value->at(1)= 4.*p[1] - 1;
                        value->at(2)= 0.;
                        break;
                    case 3:
                        value->at(0)= 0.;
                        value->at(1)= 0.;
                        value->at(2)= 4.*p[2] - 1;
                        break;
                    case 4:
                        value->at(0)= 4. - 8.*p[0] - 4.*p[1] - 4.*p[2];
                        value->at(1)= - 4.*p[0];
                        value->at(2)= - 4.*p[0];
                        break;
                    case 5:
                        value->at(0)= 4.*p[1];
                        value->at(1)= 4.*p[0];
                        value->at(2)= 0.;
                        break;
                    case 6:
                        value->at(0)= - 4.*p[1];
                        value->at(1)= 4. - 4.*p[0] - 8.*p[1] - 4.*p[2];
                        value->at(2)= - 4.*p[1];
                        break;
                    case 7:
                        value->at(0)= - 4.*p[2];
                        value->at(1)= - 4.*p[2];
                        value->at(2)= 4. - 4.*p[0] - 4.*p[1] - 8.*p[2];
                        break;
                    case 8:
                        value->at(0)= 4.*p[2];
                        value->at(1)= 0.;
                        value->at(2)= 4.*p[0];
                        break;
                    case 9:
                        value->at(0)= 0.;
                        value->at(1)= 4.*p[2];
                        value->at(2)= 4.*p[1];
                        break;  
                }
                break;
            }
        }
}

void Helper::phi(int dim,
                  int intFE,
                  int i,
                  vec_dbl_Type &p,
                  double* value){
    
    if (dim==1) {
        switch (intFE) {
            case 0: //P0
                switch (i) {
                    case 0:
                        *value = 1.;
                        break;
                }
                break;
            case 1://P1
                switch (i) {
                    case 0:
                        *value = ( 1. - p.at(0) );
                        break;
                    case 1:
                        *value = p.at(0);
                        break;
                }
                break;
            case 2://P2
                switch (i) {
                    case 0:
                        *value = ( 1. - 3. * p[0] + 2. * p[0] *  p[0] );
                        break;
                    case 1:
                        *value = ( - p[0] + 2. * p[0] *  p[0] );
                        break;
                    case 2:
                        *value = ( 4. * p[0] - 4. * p[0] *  p[0] );
                        break;
                        
                }
                break;
            default:
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "Only P0,P1,P2 1D basis functions available." );
                break;
        }
    }
    else if (dim==2) {
        switch (intFE) {
            case 0://P0
                switch (i) {
                    case 0:
                        *value = 1.;
                        break;
                }
                break;
            case 1://P1
                switch (i) {
                    case 0:
                        *value = (1. - p.at(0)-p.at(1));
                        break;
                    case 1:
                        *value = p.at(0);
                        break;
                    case 2:
                        *value = p.at(1);
                        break;
                }
                break;
            case 2://P2
                switch (i) {
                    case 0:
                        *value = -(1. - p.at(0)-p.at(1)) * (1 - 2.*(1-p.at(0) - p.at(1)));
                        break;
                    case 1:
                        *value = -p.at(0) *  (1 - 2*p.at(0));
                        break;
                    case 2:
                        *value = -p.at(1) *  (1 - 2*p.at(1));
                        break;
                    case 3:
                        *value = 4*p.at(0) * (1 - p.at(0)-p.at(1));
                        break;
                    case 4:
                        *value = 4*p.at(0)*p.at(1);
                        break;
                    case 5:
                        *value = 4*p.at(1) * (1 - p.at(0)-p.at(1));
                        break;
                }
                break;
        }
    }
    else if(dim==3){
        switch (intFE) {
            case 1://P1
                switch (i) {
                    case 0:
                        *value = (1. - p.at(0)-p.at(1)-p.at(2));
                        break;
                    case 1:
                        *value = p.at(0);
                        break;
                    case 2:
                        *value = p.at(1);
                        break;
                    case 3:
                        *value = p.at(2);
                        break;
                }
                break;
            case 2: //P2
                switch (i) {
                    case 0:
                        *value = (1. - p.at(0)-p.at(1)-p.at(2)) * (1 - 2*p.at(0) - 2*p.at(1) - 2*p.at(2));
                        break;
                    case 1:
                        *value = p.at(0) * (2*p.at(0) - 1);
                        break;
                    case 2:
                        *value = p.at(1) * (2*p.at(1) - 1);
                        break;
                    case 3:
                        *value = p.at(2) * (2*p.at(2) - 1);
                        break;
                    case 4:
                        *value = 4*p.at(0) * (1 - p.at(0)-p.at(1)-p.at(2));
                        break;
                    case 5:
                        *value = 4*p.at(0)*p.at(1);
                        break;
                    case 6:
                        *value = 4*p.at(1) * (1 - p.at(0)-p.at(1)-p.at(2));
                        break;
                    case 7:
                        *value = 4*p.at(2) * (1 - p.at(0)-p.at(1)-p.at(2));
                        break;
                    case 8:
                        *value = 4*p.at(0)*p.at(2);
                        break;
                    case 9:
                        *value = 4*p.at(1)*p.at(2);
                        break;
                }
                break;
                           
              }
                
        }

}

int Helper::getPhi(vec2D_dbl_ptr_Type &Phi,
                            vec_dbl_ptr_Type &weightsPhi,
                            int dim,
                            std::string FEType,
                            int Degree,
                            std::string FETypeQuadPoints){

    int 			nmbLocElPts;
    int 			intFE;
    double  		value;
    vec2D_dbl_ptr_Type	QuadPts;
    if (dim==1) {
        getQuadratureValues(dim, Degree, QuadPts, weightsPhi, FEType);
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 2;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 3;
            intFE = 2;
        }
        Phi.reset( new vec2D_dbl_Type( weightsPhi->size(), vec_dbl_Type( nmbLocElPts, 0.0 ) ) );
        for (int k=0; k<Phi->size(); k++ ){
            for (int i=0; i<Phi->at(0).size(); i++) {
                phi(dim,intFE,i,QuadPts->at(k),&value);
                Phi->at(k).at(i) = value;
            }
        }

    }
    else if (dim==2) {
        getQuadratureValues(dim, Degree, QuadPts, weightsPhi, FEType);
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 3;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 6;
            intFE = 2;
        }

        Phi.reset(new vec2D_dbl_Type(weightsPhi->size(),vec_dbl_Type(nmbLocElPts,0.0)));

        for (int k=0; k<Phi->size(); k++ ){
            for (int i=0; i<Phi->at(0).size(); i++) {
                phi(dim,intFE,i,QuadPts->at(k),&value);
                Phi->at(k).at(i) = value;
            }
        }
    }
    else if(dim==3){
        if (FETypeQuadPoints!="")
            getQuadratureValues(dim, Degree, QuadPts, weightsPhi, FETypeQuadPoints);
        else
            getQuadratureValues(dim, Degree, QuadPts, weightsPhi, FEType);
        
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 4;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 10;
            intFE = 2;
        }
        else if (FEType == "Q1") {
            nmbLocElPts = 8;
            intFE = 3;
        }
        else if (FEType == "Q2") {
            nmbLocElPts = 27;
            intFE = 4;
        }
        else if (FEType == "Q2-20") {
            nmbLocElPts = 20;
            intFE = 5;
        }
        else if (FEType == "P1-disc") {
            nmbLocElPts = 4;
            intFE = 1;
        }
        else if (FEType == "P1-disc-global")
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "P1-disc-global not implemented yet.");
        

        Phi.reset(new vec2D_dbl_Type(weightsPhi->size(),vec_dbl_Type(nmbLocElPts,0.0)));

        for (int k=0; k<Phi->size(); k++ ){
            for (int i=0; i<Phi->at(0).size(); i++) {
                phi(dim,intFE,i,QuadPts->at(k),&value);
                Phi->at(k).at(i) = value;
            }
        }
    }
    return intFE;
}

int Helper::getFuncAtQuadNodes(vec_dbl_ptr_Type &funcVals,
                               RhsFunc_Type &rhsFunc, int dim,
                               std::string FEType, int Degree,
                               std::string FETypeQuadPoints) {

    int nmbLocElPts;
    int intFE;
    double value;
    vec2D_dbl_ptr_Type QuadPts;
    vec_dbl_ptr_Type weightsPhi = Teuchos::rcp(new vec_dbl_Type(0));
    // dummy var for passing to rhs func
    std::vector<double> paras(1);
    if (dim == 1 || dim == 2) {
        getQuadratureValues(dim, Degree, QuadPts, weightsPhi, FEType);

        funcVals.reset(new vec_dbl_Type(weightsPhi->size(), 0.0));
        for (int k = 0; k < funcVals->size(); k++) {
            rhsFunc(QuadPts->at(k).data(), &funcVals->at(k), paras.data());
        }
    } else if (dim == 3) {
        if (FETypeQuadPoints != "") {
            getQuadratureValues(dim, Degree, QuadPts, weightsPhi,
                                FETypeQuadPoints);
        } else {
            getQuadratureValues(dim, Degree, QuadPts, weightsPhi, FEType);
        }
        funcVals.reset(new vec_dbl_Type(weightsPhi->size(), 0.0));

        for (int k = 0; k < funcVals->size(); k++) {
            rhsFunc(QuadPts->at(k).data(), &funcVals->at(k), paras.data());
        }
    }
    return intFE;
}

vec2D_dbl_Type Helper::getQuadratureValuesOnSurface(int dim, std::string FEType, vec_dbl_Type &QuadW, vec_LO_Type surfaceIDs, vec2D_dbl_ptr_Type points){

	vec2D_dbl_Type QuadPts(QuadW.size(), vec_dbl_Type(dim));
	
	if(dim==2){
		double x0 = points->at(surfaceIDs.at(0)).at(0);
		double y0 = points->at(surfaceIDs.at(0)).at(1);
		double x1 = points->at(surfaceIDs.at(1)).at(0);
		double y1 = points->at(surfaceIDs.at(1)).at(1);
		

		if(FEType == "P1"){
			
			QuadPts[0][0] =  (x0+x1)/2.;
			QuadPts[0][1] =  (y0+y1)/2.;

			QuadW[0] = 1.;
		}
		else if(FEType == "P2"){

			QuadPts[0][0] =  x0;
			QuadPts[0][1] =  y0;
			QuadPts[1][0] =  (x0+x1)/2.;
			QuadPts[1][1] =  (y0+y1)/2.;
			QuadPts[2][0] =  x1;
			QuadPts[2][1] =  y1;

			QuadW[0] = 1.;
			QuadW[1] = 4.;
			QuadW[2] = 1.;
		}
		
	}	
	else if(dim==3){
		// Here we choose as quadpoints the midpoints of the triangle sides
		double x0 = points->at(surfaceIDs.at(0)).at(0);
		double y0 = points->at(surfaceIDs.at(0)).at(1);
		double z0 = points->at(surfaceIDs.at(0)).at(2);
		double x1 = points->at(surfaceIDs.at(1)).at(0);
		double y1 = points->at(surfaceIDs.at(1)).at(1);
		double z1 = points->at(surfaceIDs.at(1)).at(2);
		double x2 = points->at(surfaceIDs.at(2)).at(0);
		double y2 = points->at(surfaceIDs.at(2)).at(1);
		double z2 = points->at(surfaceIDs.at(2)).at(2);

		if(FEType == "P1"){
			// In my case: As nabla phi is a constant function, quad points don't really matter in that case ...
			QuadPts[0][0] =   1/3.;
			QuadPts[0][1] =   1/3.;
			QuadPts[0][2] =   1/3.;

			QuadW[0] = 1.;
		}
		else if(FEType == "P2"){
			QuadPts[0][0] =  (x0+x1)/2.;
			QuadPts[0][1] =  (y0+y1)/2.;
			QuadPts[0][2] =  (z0+z1)/2.;
			QuadPts[1][0] =  (x0+x2)/2.;
			QuadPts[1][1] =  (y0+y2)/2.;
			QuadPts[1][2] =  (z0+z2)/2.;
			QuadPts[2][0] =  (x1+x2)/2.;
			QuadPts[2][1] =  (y1+y2)/2.;
			QuadPts[2][2] =  (z1+z2)/2.;

			QuadW[0] = 1/3.;
			QuadW[1] = 1/3.;
			QuadW[2] = 1/3.;
		}
	}

	return QuadPts;	

}

void Helper::getQuadratureValues(int dim,
                                  int Degree,
                                  vec2D_dbl_ptr_Type &QuadPts,
                                  vec_dbl_ptr_Type &QuadW,
                                  std::string FEType){
    // quadrature formulas exact up to a certain polynomial degree

    TEUCHOS_TEST_FOR_EXCEPTION(Degree<0, std::runtime_error, "Quadrature rule: negative degree specified.");
    double a, b, c, P1, P2, volume_ref_element;

    double b1,b2,c1,c2,d,e,f,g,h,i,j;
    if (dim == 1){
        // points are for interval [0,1]
        volume_ref_element = 1.0; // length of unit interval [0,1]

        TEUCHOS_TEST_FOR_EXCEPTION(Degree>7, std::runtime_error, "Quadrature rule in 1D only up to degree 7.");

        if (Degree <= 1) {
            QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(1,0.0)));
            QuadW->resize(1);
            QuadPts->at(0).at(0) = 0.5;
            QuadW->at(0) = 1.0;
        } else if (Degree <= 3) {
            QuadPts.reset(new vec2D_dbl_Type(2,vec_dbl_Type(1,0.0)));
            QuadW->resize(2);
            QuadPts->at(0).at(0) = 0.5 - 0.5/std::sqrt(3.0);
            QuadPts->at(1).at(0) = 0.5 + 0.5/std::sqrt(3.0);
            QuadW->at(0) = 0.5;
            QuadW->at(1) = 0.5;
	} else if (Degree <= 5) {
            QuadPts.reset(new vec2D_dbl_Type(3,vec_dbl_Type(1,0.0)));
            QuadW->resize(3);
            QuadPts->at(0).at(0) = 0.5 - 0.5*std::sqrt(3.0/5.0);
            QuadPts->at(1).at(0) = 0.5;
            QuadPts->at(2).at(0) = 0.5 + 0.5*std::sqrt(3.0/5.0);
            QuadW->at(0) = 5.0 / 18.0;
            QuadW->at(1) = 8.0 / 18.0;
            QuadW->at(2) = 5.0 / 18.0;
        } else if (Degree <= 7) {
            QuadPts.reset(new vec2D_dbl_Type(4,vec_dbl_Type(1,0.0)));
            QuadW->resize(4);
            QuadPts->at(0).at(0) = 0.5 * (1.0 - std::sqrt( 3.0/7.0 + 2.0/7.0*std::sqrt(6.0/5.0) ));
            QuadPts->at(1).at(0) = 0.5 * (1.0 - std::sqrt( 3.0/7.0 - 2.0/7.0*std::sqrt(6.0/5.0) ));
            QuadPts->at(2).at(0) = 0.5 * (1.0 + std::sqrt( 3.0/7.0 - 2.0/7.0*std::sqrt(6.0/5.0) ));
            QuadPts->at(3).at(0) = 0.5 * (1.0 + std::sqrt( 3.0/7.0 + 2.0/7.0*std::sqrt(6.0/5.0) ));
            QuadW->at(0) = 0.5 * (18.0 - std::sqrt(30.0))/36.0;
            QuadW->at(1) = 0.5 * (18.0 + std::sqrt(30.0))/36.0;
            QuadW->at(2) = 0.5 * (18.0 + std::sqrt(30.0))/36.0;
            QuadW->at(3) = 0.5 * (18.0 - std::sqrt(30.0))/36.0;
        }
    } else if (dim == 2) {
        volume_ref_element = 0.5; // area of unit triangle

        TEUCHOS_TEST_FOR_EXCEPTION(Degree>7, std::runtime_error, "Quadrature rules for triangles only up to degree 7.");

        if (Degree <= 1) {
            QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(2,0.0)));
            QuadW->resize(1);
            QuadPts->at(0).at(0) = 1/3.;
            QuadPts->at(0).at(1) = 1/3.;
            QuadW->at(0) = 1/2.;
        } else if (Degree <= 2) {
            QuadPts.reset(new vec2D_dbl_Type(3,vec_dbl_Type(2,0.0)));
            QuadW->resize(3);

            QuadPts->at(0).at(0) = 0.5;
            QuadPts->at(0).at(1) = 0.5;

            QuadPts->at(1).at(0) = 0.0;
            QuadPts->at(1).at(1) = 0.5;

            QuadPts->at(2).at(0) = 0.5;
            QuadPts->at(2).at(1) = 0.0;

            a = 1/6.;
            QuadW->at(0) = a;
            QuadW->at(1) = a;
            QuadW->at(2) = a;
        } else if (Degree <= 3) {
            QuadPts.reset(new vec2D_dbl_Type(4,vec_dbl_Type(2,0.0)));
            QuadW->resize(4);

            QuadPts->at(0).at(0) = 1.0/3.0;
            QuadPts->at(0).at(1) = 1.0/3.0;

            QuadPts->at(1).at(0) = 0.6;
            QuadPts->at(1).at(1) = 0.2;

            QuadPts->at(2).at(0) = 0.2;
            QuadPts->at(2).at(1) = 0.6;

            QuadPts->at(3).at(0) = 0.2;
            QuadPts->at(3).at(1) = 0.2;

            a = 0.2604166666666665;
            QuadW->at(0) = -0.28125;
            QuadW->at(1) = a;
            QuadW->at(2) = a;
            QuadW->at(3) = a;
        } else if (Degree <= 3) {
            // [1] Dunavant, (1985)
            //     High Degree Efficient Symmetrical Gaussian Quadrature Rules For The Triangle
            //
            // [2] Akin, (1994)
            //     Finite elements for analysis and design (p.200)
            QuadPts.reset(new vec2D_dbl_Type(4,vec_dbl_Type(2,0.0)));
            QuadW->resize(4);

            QuadPts->at(0).at(0) = 1.0/3.0;
            QuadPts->at(0).at(1) = 1.0/3.0;

            QuadPts->at(1).at(0) = 0.6;
            QuadPts->at(1).at(1) = 0.2;

            QuadPts->at(2).at(0) = 0.2;
            QuadPts->at(2).at(1) = 0.6;

            QuadPts->at(3).at(0) = 0.2;
            QuadPts->at(3).at(1) = 0.2;

            a = 0.2604166666666665;
            QuadW->at(0) = -0.28125;
            QuadW->at(1) = a;
            QuadW->at(2) = a;
            QuadW->at(3) = a;
        } else if (Degree <= 4) {
            // [1] Dunavant, (1985)
            //     High Degree Efficient Symmetrical Gaussian Quadrature Rules For The Triangle
            //
            // [2] Akin, (1994)
            //     Finite elements for analysis and design (p.200)
            QuadPts.reset(new vec2D_dbl_Type(6,vec_dbl_Type(2,0.0)));
            QuadW->resize(6);

            a = 0.10810301816807;
            b = 0.445948490915965;
            QuadPts->at(0).at(0) = a;
            QuadPts->at(0).at(1) = b;

            QuadPts->at(1).at(0) = b;
            QuadPts->at(1).at(1) = a;

            QuadPts->at(2).at(0) = b;
            QuadPts->at(2).at(1) = b;

            a = 0.816847572980459;
            b = 0.09157621350977101;
            QuadPts->at(3).at(0) = a;
            QuadPts->at(3).at(1) = b;

            QuadPts->at(4).at(0) = b;
            QuadPts->at(4).at(1) = a;

            QuadPts->at(5).at(0) = b;
            QuadPts->at(5).at(1) = b;

            a = 0.1116907948390055;
            b = 0.054975871827661;
            QuadW->at(0) = a;
            QuadW->at(1) = a;
            QuadW->at(2) = a;
            QuadW->at(3) = b;
            QuadW->at(4) = b;
            QuadW->at(5) = b;
        } else if (Degree <= 5) {
            QuadPts.reset(new vec2D_dbl_Type(7,vec_dbl_Type(2,0.0)));
            QuadW->resize(7);

            a = 0.470142064105115;
            b = 0.101286507323456;
            P1 = 0.066197076394253;
            P2 = 0.062969590272413;

            QuadPts->at(0).at(0) = 1/3.;
            QuadPts->at(0).at(1) = 1/3.;

            QuadPts->at(1).at(0) = a;
            QuadPts->at(1).at(1) = a;
				 
            QuadPts->at(2).at(0) = 1-2.*a;
            QuadPts->at(2).at(1) = a;
				 
            QuadPts->at(3).at(0) = a;
            QuadPts->at(3).at(1) = 1-2.*a;
				 
            QuadPts->at(4).at(0) = b;
            QuadPts->at(4).at(1) = b;
				 
            QuadPts->at(5).at(0) = 1-2.*b;
            QuadPts->at(5).at(1) = b;
				 
            QuadPts->at(6).at(0) = b;
            QuadPts->at(6).at(1) = 1-2.*b;

            QuadW->at(0) = 9/80.;
            QuadW->at(1) = P1;
            QuadW->at(2) = P1;
            QuadW->at(3) = P1;
            QuadW->at(4) = P2;
            QuadW->at(5) = P2;
            QuadW->at(6) = P2;
        } else if (Degree <= 6) {
            // [1] Dunavant, (1985)
            //     High Degree Efficient Symmetrical Gaussian Quadrature Rules For The Triangle
            //
            // [2] Akin, (1994)
            //     Finite elements for analysis and design (p.200)
            QuadPts.reset(new vec2D_dbl_Type(12,vec_dbl_Type(2,0.0)));
            QuadW->resize(12);

            SC quadPts_[2][12] =
            {
                {0.501426509658179,0.249286745170910,0.249286745170910,0.873821971016996,0.063089014491502,0.063089014491502,0.053145049844817,0.310352451033784,0.636502499121399,0.636502499121399,0.053145049844817,0.310352451033784},
                {0.249286745170910,0.501426509658179,0.249286745170910,0.063089014491502,0.873821971016996,0.063089014491502,0.310352451033784,0.053145049844817,0.053145049844817,0.310352451033784,0.636502499121399,0.636502499121399}
            };
            for (int i = 0; i < QuadW->size(); i++) {
                QuadPts->at(i).at(0) = quadPts_[0][i];
                QuadPts->at(i).at(1) = quadPts_[1][i];
            }

            SC quadW_[12] = {0.0583931378631895,0.0583931378631895,0.0583931378631895,0.0254224531851035,0.0254224531851035,0.0254224531851035,0.041425537809187,0.041425537809187,0.041425537809187,0.041425537809187,0.041425537809187,0.041425537809187};
            for (int i = 0; i < QuadW->size(); i++) {
                QuadW->at(i) = quadW_[i];
            }
        } else if (Degree <= 7) {
            // [1] Dunavant, (1985)
            //     High Degree Efficient Symmetrical Gaussian Quadrature Rules For The Triangle
            //
            // [2] Akin, (1994)
            //     Finite elements for analysis and design (p.200)
            QuadPts.reset(new vec2D_dbl_Type(13,vec_dbl_Type(2,0.0)));
            QuadW->resize(13);

            SC quadPts_[2][13] =
            {
                {0.333333333333333,0.47930806784192,0.26034596607904,0.26034596607904,0.869739794195568,0.06513010290221601,0.06513010290221601,0.048690315425316,0.312865496004874,0.63844418856981,0.63844418856981,0.048690315425316,0.312865496004874},
                {0.333333333333333,0.26034596607904,0.47930806784192,0.26034596607904,0.06513010290221601,0.869739794195568,0.06513010290221601,0.312865496004874,0.048690315425316,0.048690315425316,0.312865496004874,0.63844418856981,0.63844418856981}
            };
            for (int i = 0; i < QuadW->size(); i++) {
                QuadPts->at(i).at(0) = quadPts_[0][i];
                QuadPts->at(i).at(1) = quadPts_[1][i];
            }

            SC quadW_[13] = {-0.074785022233841,0.087807628716604,0.087807628716604,0.087807628716604,0.026673617804419,0.026673617804419,0.026673617804419,0.0385568804451285,0.0385568804451285,0.0385568804451285,0.0385568804451285,0.0385568804451285,0.0385568804451285};
            for (int i = 0; i < QuadW->size(); i++) {
                QuadW->at(i) = quadW_[i];
            }
        }
    } else if(dim == 3) {
        volume_ref_element = 1.0/6.0; // area of unit tetrahedron

        if (FEType.at(0)=='P') { // TODO: [JK] Why is the element type only queried here and not in 2D?

            TEUCHOS_TEST_FOR_EXCEPTION(Degree>6, std::runtime_error, "Tetrahedron quadrature rules only up to degree 6 available.");
            
            if (Degree <= 1) {
                QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(3,0.0)));
                QuadW->resize(1);
                QuadPts->at(0).at(0) = 0.25;
                QuadPts->at(0).at(1) = 0.25;
                QuadPts->at(0).at(2) = 0.25;
                QuadW->at(0)         = 1/6.;
            } else if (Degree <= 2) {
                // Akin - 1994 - Finite elements for analysis and design (p.198)

                QuadPts.reset(new vec2D_dbl_Type(4,vec_dbl_Type(3,0.0)));
                QuadW->resize(4);

                a = (5.0 + 3.0*std::sqrt(5.0))/20.0;
                b = (5.0 - 1.0*std::sqrt(5.0))/20.0;

                QuadPts->at(0).at(0) = a;
                QuadPts->at(0).at(1) = b;
                QuadPts->at(0).at(2) = b;

                QuadPts->at(1).at(0) = b;
                QuadPts->at(1).at(1) = a;
                QuadPts->at(1).at(2) = b;

                QuadPts->at(2).at(0) = b;
                QuadPts->at(2).at(1) = b;
                QuadPts->at(2).at(2) = a;

                QuadPts->at(3).at(0) = b;
                QuadPts->at(3).at(1) = b;
                QuadPts->at(3).at(2) = b;

                a = 1.0/24.0;
                QuadW->at(0) = a;
                QuadW->at(1) = a;
                QuadW->at(2) = a;
                QuadW->at(3) = a;
            } else if (Degree <= 3) {
                QuadPts.reset(new vec2D_dbl_Type(5,vec_dbl_Type(3,0.0)));
                QuadW->resize(5);
                a = .25;
                b = 1./6.;
                c = .5;
                QuadPts->at(0).at(0) = a;
                QuadPts->at(0).at(1) = a;
                QuadPts->at(0).at(2) = a;
                
                QuadPts->at(1).at(0) = b;
                QuadPts->at(1).at(1) = b;
                QuadPts->at(1).at(2) = b;
                
                QuadPts->at(2).at(0) = b;
                QuadPts->at(2).at(1) = b;
                QuadPts->at(2).at(2) = c;
                
                QuadPts->at(3).at(0) = b;
                QuadPts->at(3).at(1) = c;
                QuadPts->at(3).at(2) = b;
                
                QuadPts->at(4).at(0) = c;
                QuadPts->at(4).at(1) = b;
                QuadPts->at(4).at(2) = b;
                
                QuadW->at(0) = -2./15.;
                QuadW->at(1) = 3./40.;
                QuadW->at(2) = 3./40.;
                QuadW->at(3) = 3./40.;
                QuadW->at(4) = 3./40.;
            } else if (Degree <= 4) {
                QuadPts.reset(new vec2D_dbl_Type(11,vec_dbl_Type(3,0.0)));
                QuadW->resize(11);
                
                a = .785714285714286;
                b = .071428571428571;
                c = .100596423833201;
                d = .399403576166799;
                
                QuadPts->at(0).at(0) 	= .25;
                QuadPts->at(0).at(1)    = .25;
                QuadPts->at(0).at(2)    = .25;
                
                QuadPts->at(1).at(0) 	= a;
                QuadPts->at(1).at(1)    = b;
                QuadPts->at(1).at(2)    = b;
                
                QuadPts->at(2).at(0) 	= b;
                QuadPts->at(2).at(1)    = b;
                QuadPts->at(2).at(2)    = b;
                
                QuadPts->at(3).at(0) 	= b;
                QuadPts->at(3).at(1)    = b;
                QuadPts->at(3).at(2)    = a;
                
                QuadPts->at(4).at(0) 	= b;
                QuadPts->at(4).at(1)    = a;
                QuadPts->at(4).at(2)    = b;
                
                QuadPts->at(5).at(0) 	= c;
                QuadPts->at(5).at(1)    = d;
                QuadPts->at(5).at(2)    = d;
                
                QuadPts->at(6).at(0) 	= d;
                QuadPts->at(6).at(1)    = c;
                QuadPts->at(6).at(2)    = d;
                
                QuadPts->at(7).at(0) 	= d;
                QuadPts->at(7).at(1)    = d;
                QuadPts->at(7).at(2)    = c;
                
                QuadPts->at(8).at(0) 	= d;
                QuadPts->at(8).at(1)    = c;
                QuadPts->at(8).at(2)    = c;
                
                QuadPts->at(9).at(0) 	= c;
                QuadPts->at(9).at(1)    = d;
                QuadPts->at(9).at(2)    = c;
                
                QuadPts->at(10).at(0) 	= c;
                QuadPts->at(10).at(1)   = c;
                QuadPts->at(10).at(2)   = d;
                
                a = -0.078933333333333/6.0;
                b = 0.045733333333333/6.0;
                c = 0.149333333333333/6.0;
                
                QuadW->at(0) = a;
                
                QuadW->at(1) = b;
                QuadW->at(2) = b;
                QuadW->at(3) = b;
                QuadW->at(4) = b;
                
                QuadW->at(5) = c;
                QuadW->at(6) = c;
                QuadW->at(7) = c;
                QuadW->at(8) = c;
                QuadW->at(9) = c;
                QuadW->at(10) = c;
            } else if (Degree <= 5) {
                QuadPts.reset(new vec2D_dbl_Type(15,vec_dbl_Type(3,0.0)));
                QuadW->resize(15);
                a 	= 0.25;
                b1 	= (7.+std::sqrt(15.))/34.;
                b2 	= (7.-std::sqrt(15.))/34.;
                c1 	= (13.-3.*std::sqrt(15.))/34.;
                c2 	= (13.+3.*std::sqrt(15.))/34.;
                d 	= (5.-std::sqrt(15.))/20.;
                e 	= (5.+std::sqrt(15.))/20.;
                
                QuadPts->at(0).at(0) 	= a;
                QuadPts->at(0).at(1)    = a;
                QuadPts->at(0).at(2)    = a;
                
                QuadPts->at(1).at(0) 	= b1;
                QuadPts->at(1).at(1)    = b1;
                QuadPts->at(1).at(2)    = b1;
                
                QuadPts->at(2).at(0) 	= b1;
                QuadPts->at(2).at(1)    = b1;
                QuadPts->at(2).at(2)    = c1;
                
                QuadPts->at(3).at(0) 	= b1;
                QuadPts->at(3).at(1)    = c1;
                QuadPts->at(3).at(2)    = b1;
                
                QuadPts->at(4).at(0) 	= c1;
                QuadPts->at(4).at(1)    = b1;
                QuadPts->at(4).at(2)    = b1;
                
                QuadPts->at(5).at(0) 	= b2;
                QuadPts->at(5).at(1)    = b2;
                QuadPts->at(5).at(2)    = b2;
                
                QuadPts->at(6).at(0) 	= b2;
                QuadPts->at(6).at(1)    = b2;
                QuadPts->at(6).at(2)    = c2;
                
                QuadPts->at(7).at(0) 	= b2;
                QuadPts->at(7).at(1)    = c2;
                QuadPts->at(7).at(2)    = b2;
                
                QuadPts->at(8).at(0) 	= c2;
                QuadPts->at(8).at(1)    = b2;
                QuadPts->at(8).at(2)    = b2;
                
                QuadPts->at(9).at(0) 	= d;
                QuadPts->at(9).at(1)    = d;
                QuadPts->at(9).at(2)    = e;
                
                QuadPts->at(10).at(0) 	= d;
                QuadPts->at(10).at(1)   = e;
                QuadPts->at(10).at(2)   = d;
                
                QuadPts->at(11).at(0) 	= e;
                QuadPts->at(11).at(1)	= d;
                QuadPts->at(11).at(2)	= d;
                
                QuadPts->at(12).at(0) 	= d;
                QuadPts->at(12).at(1)	= e;
                QuadPts->at(12).at(2)	= e;
                
                QuadPts->at(13).at(0) 	= e;
                QuadPts->at(13).at(1)	= d;
                QuadPts->at(13).at(2)	= e;
                
                QuadPts->at(14).at(0) 	= e;
                QuadPts->at(14).at(1)	= e;
                QuadPts->at(14).at(2)	= d;
                
                
                P1 	= (2665.-14.*std::sqrt(15.))/226800.;
                P2 	= (2665.+14.*std::sqrt(15.))/226800.;
                b	= 5./567.;
                
                QuadW->at(0) = 8./405.;
                QuadW->at(1) = P1;
                QuadW->at(2) = P1;
                QuadW->at(3) = P1;
                QuadW->at(4) = P1;

                QuadW->at(5) = P2;
                QuadW->at(6) = P2;
                QuadW->at(7) = P2;
                QuadW->at(8) = P2;

                QuadW->at(9) 			= b;
                QuadW->at(10)           = b;
                QuadW->at(11)           = b;
                QuadW->at(12) 			= b;
                QuadW->at(13)           = b;
                QuadW->at(14)           = b;
            } else if (Degree <= 6) {
                // Keast - 1985 - Moderate-Degree Tetrahedral Quadrature Formulas (p.342)
                QuadPts.reset(new vec2D_dbl_Type(24,vec_dbl_Type(3,0.0)));
                QuadW->resize(24);
                a = .356191386222544953;
                b = .214602871259151684;
                c = .877978124396165982;
                d = .406739585346113397/10.0;
                f = .329863295731730594/10.0;
                g = .322337890142275646;
                h = .269672331458315867;
                i = .636610018750175299/10.0;
                j = .603005664791649076;
                
                QuadPts->at(0).at(0) 	= a;
                QuadPts->at(0).at(1)    = b;
                QuadPts->at(0).at(2)    = b;

                QuadPts->at(1).at(0) 	= b;
                QuadPts->at(1).at(1)    = b;
                QuadPts->at(1).at(2)    = b;

                QuadPts->at(2).at(0) 	= b;
                QuadPts->at(2).at(1)    = b;
                QuadPts->at(2).at(2)    = a;

                QuadPts->at(3).at(0) 	= b;
                QuadPts->at(3).at(1)    = a;
                QuadPts->at(3).at(2)    = b;

                QuadPts->at(4).at(0) 	= c;
                QuadPts->at(4).at(1)    = d;
                QuadPts->at(4).at(2)    = d;
                
                QuadPts->at(5).at(0) 	= d;
                QuadPts->at(5).at(1)    = d;
                QuadPts->at(5).at(2)    = d;

                QuadPts->at(6).at(0) 	= d;
                QuadPts->at(6).at(1)    = d;
                QuadPts->at(6).at(2)    = c;
                
                QuadPts->at(7).at(0) 	= d;
                QuadPts->at(7).at(1)    = c;
                QuadPts->at(7).at(2)    = d;

                QuadPts->at(8).at(0) 	= f;
                QuadPts->at(8).at(1)    = g;
                QuadPts->at(8).at(2)    = g;

                QuadPts->at(9).at(0) 	= g;
                QuadPts->at(9).at(1)    = g;
                QuadPts->at(9).at(2)    = g;

                QuadPts->at(10).at(0) 	= g;
                QuadPts->at(10).at(1)   = g;
                QuadPts->at(10).at(2)   = f;

                QuadPts->at(11).at(0) 	= g;
                QuadPts->at(11).at(1)   = f;
                QuadPts->at(11).at(2)   = g;
                
                QuadPts->at(12).at(0) 	= h;
                QuadPts->at(12).at(1)   = i;
                QuadPts->at(12).at(2)   = i;
                
                QuadPts->at(13).at(0) 	= i;
                QuadPts->at(13).at(1)   = h;
                QuadPts->at(13).at(2)   = i;
                
                QuadPts->at(14).at(0) 	= i;
                QuadPts->at(14).at(1)   = i;
                QuadPts->at(14).at(2)   = h;
                
                QuadPts->at(15).at(0) 	= j;
                QuadPts->at(15).at(1)   = i;
                QuadPts->at(15).at(2)   = i;

                QuadPts->at(16).at(0) 	= i;
                QuadPts->at(16).at(1)   = j;
                QuadPts->at(16).at(2)   = i;

                QuadPts->at(17).at(0) 	= i;
                QuadPts->at(17).at(1)   = i;
                QuadPts->at(17).at(2)   = j;
                
                QuadPts->at(18).at(0) 	= i;
                QuadPts->at(18).at(1)   = h;
                QuadPts->at(18).at(2)   = j;

                QuadPts->at(19).at(0) 	= h;
                QuadPts->at(19).at(1)   = j;
                QuadPts->at(19).at(2)   = i;
                
                QuadPts->at(20).at(0) 	= j;
                QuadPts->at(20).at(1)   = i;
                QuadPts->at(20).at(2)   = h;
                
                QuadPts->at(21).at(0) 	= i;
                QuadPts->at(21).at(1)   = j;
                QuadPts->at(21).at(2)   = h;

                QuadPts->at(22).at(0) 	= h;
                QuadPts->at(22).at(1)   = i;
                QuadPts->at(22).at(2)   = j;
                
                QuadPts->at(23).at(0) 	= j;
                QuadPts->at(23).at(1)   = h;
                QuadPts->at(23).at(2)   = i;
                
                a = 0.665379170969464506/100.0;
                b = 0.167953517588677620/100.0;
                c = 0.922619692394239843/100.0;
                d = 0.803571428571428248/100.0;
                
                QuadW->at(0)  = a;
                QuadW->at(1)  = a;
                QuadW->at(2)  = a;
                QuadW->at(3)  = a;
                QuadW->at(4)  = b;
                QuadW->at(5)  = b;
                QuadW->at(6)  = b;
                QuadW->at(7)  = b;
                QuadW->at(8)  = c;
                QuadW->at(9)  = c;
                QuadW->at(10) = c;
                QuadW->at(11) = c;
                QuadW->at(12) = d;
                QuadW->at(13) = d;
                QuadW->at(14) = d;
                QuadW->at(15) = d;
                QuadW->at(16) = d;
                QuadW->at(17) = d;
                QuadW->at(18) = d;
                QuadW->at(19) = d;
                QuadW->at(20) = d;
                QuadW->at(21) = d;
                QuadW->at(22) = d;
                QuadW->at(23) = d;
            } // quadrature formulas for tetrahedra
        } // Pk finite elements (P1, P2, ...)
    } // 3D quadrature formulas

    {
        // Sanity check: do weights sum to 1 for unit interval [0,1], 1/2 for unit triangle, 1/6 for unit tetrahedron.
        // This does not imply a correct quadrature rule, but it must hold.
        double volume = 0.0;
        for (int i = 0; i < QuadW->size(); i++) {
            volume += QuadW->at(i);
        }
        double error_volume = std::fabs(volume - volume_ref_element);
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(2) << error_volume;
        std::string str = oss.str();
        TEUCHOS_TEST_FOR_EXCEPTION(
            error_volume > std::numeric_limits<double>::epsilon()*10.0, 
            std::runtime_error, 
            "Quadrature weights do not sum to (approximately) " << volume_ref_element << 
            ", the length/area/volume of the reference element. [error:" + str  + "].");
    }
}


void Helper::getDPhiAtCM(vec3D_dbl_ptr_Type &DPhi,
                         int dim,
                         std::string FEType)
{
    int	nmbLocElPts;
    int	intFE;
    vec_dbl_Type CM(dim, 0.0);
    vec_dbl_ptr_Type value(new vec_dbl_Type(dim,0.0));
    TEUCHOS_TEST_FOR_EXCEPTION(dim == 1,std::logic_error, "getDPhiAtCMNot implemented for dim=1");

    if (dim==2)
    {
        
        // As we are in the reference element the center of mass is just:
        CM[0] = 1.0 / 3.0;
        CM[1] = 1.0 / 3.0;
        
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 3;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 6;
            intFE = 2;
        }

        DPhi.reset(new vec3D_dbl_Type(1,vec2D_dbl_Type(nmbLocElPts,vec_dbl_Type(2,0.0))));

        for (int k=0; k<DPhi->size(); k++ )
        {
            for (int i=0; i<DPhi->at(0).size(); i++) 
            {
                gradPhi(dim,intFE,i,CM,value);
                for (int j=0; j<2; j++) 
                {
                    DPhi->at(k).at(i).at(j) = value->at(j);
                }
            }
        }
    }
    else if(dim==3)
    {
        // As we are in the reference element the center of mass is just:
        CM[0] = 1.0 / 4.0;
        CM[1] = 1.0 / 4.0;
        CM[2] = 1.0 / 4.0;

        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 4;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 10;
            intFE = 2;
        }
        DPhi.reset(new vec3D_dbl_Type(1,vec2D_dbl_Type(nmbLocElPts,vec_dbl_Type(3,0.0))));
        for (int k=0; k<DPhi->size(); k++ )
        {
            for (int i=0; i<DPhi->at(0).size(); i++)
            {
                gradPhi(dim,intFE,i,CM,value);
                for (int j=0; j<3; j++) 
                {
                    DPhi->at(k).at(i).at(j) = value->at(j);
                }
            }
        }
    }

}



}


