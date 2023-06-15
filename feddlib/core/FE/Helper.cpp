#include "Helper.hpp"

using namespace std;
namespace FEDD {


UN Helper::determineDegree(UN dim, std::string FEType, UN degFunc){
   
	UN deg;
    if (!FEType.compare("P0"))
        deg = 0;
    else if (!FEType.compare("P1"))
        deg = 1;
    else if (!FEType.compare("P2"))
        deg = 2;
    else if (!FEType.compare("Q2"))
        deg = 2;
    
    deg += degFunc;

    if (deg==0)
        deg = 1;

    return deg;
}

UN Helper::determineDegree(UN dim, std::string FEType, int type){
    UN deg;
    if (!FEType.compare("P0")) {
        deg = 0;
    }
    else if (!FEType.compare("P1")) {
        if (type==Std)
            deg = 1;
        else if (type==Grad)
            deg = 0;
    }
    else if (!FEType.compare("P2")) {
        if (type==Std)
            deg = 2;
        else if (type==Grad)
            deg = 1;
    }
    else if (!FEType.compare("Q2")) {
        if (type==Std)
            deg = 2;
        else if (type==Grad)
            deg = 2;
    }
    
    if (deg==0)
        deg = 1;

    return deg;
}

UN Helper::determineDegree(UN dim, std::string FEType1, std::string FEType2, int type1,int type2, UN extraDeg){

    TEUCHOS_TEST_FOR_EXCEPTION( dim==2 && ( FEType1=="P2-CR" || FEType2=="P2-CR"), std::runtime_error, "P2-CR should be only available in 3D.");
    UN deg1, deg2;
    if (!FEType1.compare("P0")) {
        deg1 = 0;
    }
    else if ( !FEType1.compare("P1") || !FEType1.compare("P1-disc") ) {
        if (type1==Std)
            deg1 = 1;
        else if (type1==Grad)
            deg1 = 0;
    }
    else if (!FEType1.compare("P2")) {
        if (type1==Std)
            deg1 = 2;
        else if (type1==Grad)
            deg1 = 1;
    }
    else if (!FEType1.compare("P2-CR")) {
        if (type1==Std)
            deg1 = 4;
        else if (type1==Grad)
            deg1 = 3;
    }

    else if (!FEType1.compare("Q2")) {
        if (type1==Std)
            deg1 = 2;
        else if (type1==Grad)
            deg1 = 2;
    }
    else if (!FEType1.compare("Q2-20")) {
        if (type1==Std)
            deg1 = 2;
        else if (type1==Grad)
            deg1 = 2;
    }
   
    if (!FEType2.compare("P0")) {
        deg2 = 0;
    }
    else if ( !FEType2.compare("P1") || !FEType2.compare("P1-disc") ) {
        if (type2==Std)
            deg2 = 1;
        else if (type2==Grad)
            deg2 = 0;
    }
    else if (!FEType2.compare("P2")) {
        if (type2==Std)
            deg2 = 2;
        else if (type2==Grad)
            deg2 = 1;
    }
    else if (!FEType2.compare("P2-CR")) {
        if (type2==Std)
            deg2 = 4;
        else if (type2==Grad)
            deg2 = 3;
    }

    else if (!FEType2.compare("Q2")) {
        if (type2==Std)
            deg2 = 2;
        else if (type2==Grad)
            deg2 = 2;
    }
    else if (!FEType2.compare("Q2-20")) {
        if (type2==Std)
            deg2 = 2;
        else if (type2==Grad)
            deg2 = 2;
    }

    UN deg = deg1+deg2+extraDeg;
    if (deg==0)
        deg = 1;
    
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
        norm_v_E = sqrt(pow(v_E[0],2)+pow(v_E[1],2));

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

        norm_v_E = sqrt(pow(v_E[0],2)+pow(v_E[1],2)+pow(v_E[2],2));

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

/*!

\brief Returns neccesary quadrature Values. Is distinguishes between needing Element or Surface information. !! Input can be improved with just delivering the coordinates of the surface nodes to determine the quad points

@param[in] dim Dimension for which the quadrature points are needed.
@param[in] FEType Finite element type for which the quadrature points are needed.
@param[in] QuadW Vector to be filled with the quadrature weights accordingly
@param[in] vec_LO_Type surfaceIDs for which you need the quadrature points.
@param[in] points The repeated(!) points of current problem to identify the surface node ids. 

@param[out] QuadPts Quadrature points
@param[out] QuadW Quadrature weights

\brief Keep in mind that elementwise quadPoints are defined on reference element whereas surface quadPoints at hand are defined on the input surface, which is typically not the reference Element. 

*/

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
    double a, b, c, P1, P2;

    double b1,b2,c1,c2,d,e,f,g,h,i,j;
    if (dim==1){
        // points are for interval [0,1]
        TEUCHOS_TEST_FOR_EXCEPTION(Degree>2, std::runtime_error, "Quadrature rule in 1d only up to degree 3.");
        switch (Degree) {
            case 0:
                QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(1,0.0)));
                QuadW->resize(1);
                QuadPts->at(0).at(0) = 0.5;
                QuadW->at(0) = 1.;
                break;
            case 1:
                QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(1,0.0)));
                QuadW->resize(1);
                QuadPts->at(0).at(0) = 0.5;
                QuadW->at(0) = 1.;
                break;
            case 2:
                QuadPts.reset(new vec2D_dbl_Type(2,vec_dbl_Type(1,0.0)));
                QuadW->resize(2);
                QuadPts->at(0).at(0) = - 0.5/sqrt(3.)+0.5;
                QuadPts->at(1).at(0) = 0.5/sqrt(3.)+0.5;
                QuadW->at(0) = .5;
                QuadW->at(1) = .5;
                break;
            case 3:
                QuadPts.reset(new vec2D_dbl_Type(2,vec_dbl_Type(1,0.0)));
                QuadW->resize(2);
                QuadPts->at(0).at(0) = - 0.5/sqrt(3.)+0.5;
                QuadPts->at(1).at(0) = 0.5/sqrt(3.)+0.5;
                QuadW->at(0) = .5;
                QuadW->at(1) = .5;
                break;
            default:
                break;
        }
    }
    if (dim==2) {

        TEUCHOS_TEST_FOR_EXCEPTION(Degree>7, std::runtime_error, "Quadrature rule in 2d only up to degree 7.");
        if (Degree==3 || Degree==4)
            Degree=5;

        if (Degree==6)
            Degree=7;
        switch (Degree) {
            case 1:

                QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(2,0.0)));
                QuadW->resize(1);
                QuadPts->at(0).at(0) = 1/3.;
                QuadPts->at(0).at(1) = 1/3.;
                QuadW->at(0)	= 1/2.;
                break;

            case 2:

                QuadPts.reset(new vec2D_dbl_Type(3,vec_dbl_Type(2,0.0)));
                QuadW->resize(3);
                a = 1/6.;
                QuadPts->at(0).at(0) 	= 0.5;
                QuadPts->at(0).at(1)    = 0.5;

                QuadPts->at(1).at(0) 	= 0.;
                QuadPts->at(1).at(1) 	= 0.5;

                QuadPts->at(2).at(0) 	= 0.5;
                QuadPts->at(2).at(1) 	= 0.;

                QuadW->at(0) 		= a;
                QuadW->at(1)            = a;
                QuadW->at(2)            = a;

                break;

            case 5:
                QuadPts.reset(new vec2D_dbl_Type(7,vec_dbl_Type(2,0.0)));
                QuadW->resize(7);
                a = 0.470142064105115;
                b = 0.101286507323456;
                P1 = 0.066197076394253;
                P2 = 0.062969590272413;

                QuadPts->at(0).at(0) 	= 1/3.;
                QuadPts->at(0).at(1)    = 1/3.;

                QuadPts->at(1).at(0) 	= a;
                QuadPts->at(1).at(1) 	= a;

                QuadPts->at(2).at(0) 	= 1-2.*a;
                QuadPts->at(2).at(1) 	= a;

                QuadPts->at(3).at(0) 	= a;
                QuadPts->at(3).at(1) 	= 1-2.*a;

                QuadPts->at(4).at(0) 	= b;
                QuadPts->at(4).at(1) 	= b;

                QuadPts->at(5).at(0) 	= 1-2.*b;
                QuadPts->at(5).at(1) 	= b;

                QuadPts->at(6).at(0) 	= b;
                QuadPts->at(6).at(1) 	= 1-2.*b;

                QuadW->at(0) 			= 9/80.;
                QuadW->at(1)            = P1;
                QuadW->at(2)            = P1;
                QuadW->at(3) 			= P1;
                QuadW->at(4)            = P2;
                QuadW->at(5)            = P2;
                QuadW->at(6)            = P2;

                break;
            case 7:
                // 28 Punkte
                
                QuadPts.reset(new vec2D_dbl_Type(28,vec_dbl_Type(2,0.0)));
                QuadW.reset(new vec_dbl_Type(28,0.0));
                
                // x punkt
                QuadPts->at(0).at(0) = 0.777777777777778;
                QuadPts->at(1).at(0) = 0.111111111111111;
                QuadPts->at(2).at(0) = 0.111111111111111;
                QuadPts->at(3).at(0) = 0.666666666666667;
                QuadPts->at(4).at(0) = 0.222222222222222;
                QuadPts->at(5).at(0) = 0.111111111111111;
                QuadPts->at(6).at(0) = 0.222222222222222;
                QuadPts->at(7).at(0) = 0.111111111111111;
                QuadPts->at(8).at(0) = 0.666666666666667;
                QuadPts->at(9).at(0) = 0.555555555555556;
                QuadPts->at(10).at(0) = 0.333333333333333;
                QuadPts->at(11).at(0) = 0.111111111111111;
                QuadPts->at(12).at(0) = 0.333333333333333;
                QuadPts->at(13).at(0) = 0.111111111111111;
                QuadPts->at(14).at(0) = 0.555555555555556;
                QuadPts->at(15).at(0) = 0.555555555555556;
                QuadPts->at(16).at(0) = 0.222222222222222;
                QuadPts->at(17).at(0) = 0.222222222222222;
                QuadPts->at(18).at(0) = 0.444444444444444;
                QuadPts->at(19).at(0) = 0.444444444444444;
                QuadPts->at(20).at(0) = 0.111111111111111;
                QuadPts->at(21).at(0) = 0.444444444444444;
                QuadPts->at(22).at(0) = 0.333333333333333;
                QuadPts->at(23).at(0) = 0.222222222222222;
                QuadPts->at(24).at(0) = 0.333333333333333;
                QuadPts->at(25).at(0) = 0.222222222222222;
                QuadPts->at(26).at(0) = 0.444444444444444;
                QuadPts->at(27).at(0) = 0.333333333333333;
                
                // y punkt
                QuadPts->at(0).at(1) = 0.111111111111111;
                QuadPts->at(1).at(1) = 0.111111111111111;
                QuadPts->at(2).at(1) = 0.777777777777778;
                QuadPts->at(3).at(1) = 0.222222222222222;
                QuadPts->at(4).at(1) = 0.111111111111111;
                QuadPts->at(5).at(1) = 0.666666666666667;
                QuadPts->at(6).at(1) = 0.666666666666667;
                QuadPts->at(7).at(1) = 0.222222222222222;
                QuadPts->at(8).at(1) = 0.111111111111111;
                QuadPts->at(9).at(1) = 0.333333333333333;
                QuadPts->at(10).at(1) = 0.111111111111111;
                QuadPts->at(11).at(1) = 0.555555555555556;
                QuadPts->at(12).at(1) = 0.555555555555556;
                QuadPts->at(13).at(1) = 0.333333333333333;
                QuadPts->at(14).at(1) = 0.111111111111111;
                QuadPts->at(15).at(1) = 0.222222222222222;
                QuadPts->at(16).at(1) = 0.222222222222222;
                QuadPts->at(17).at(1) = 0.555555555555556;
                QuadPts->at(18).at(1) = 0.444444444444444;
                QuadPts->at(19).at(1) = 0.111111111111111;
                QuadPts->at(20).at(1) = 0.444444444444444;
                QuadPts->at(21).at(1) = 0.333333333333333;
                QuadPts->at(22).at(1) = 0.222222222222222;
                QuadPts->at(23).at(1) = 0.444444444444444;
                QuadPts->at(24).at(1) = 0.444444444444444;
                QuadPts->at(25).at(1) = 0.333333333333333;
                QuadPts->at(26).at(1) = 0.222222222222222;
                QuadPts->at(27).at(1) = 0.333333333333333;
                
                // Gewichte
                QuadW->at(0) 			= 0.342410714285714/2.0;
                QuadW->at(1) 			= 0.342410714285714/2.0;
                QuadW->at(2) 			= 0.342410714285714/2.0;
                QuadW->at(3) 			= -0.561160714285714/2.0;
                QuadW->at(4) 			= -0.561160714285714/2.0;
                QuadW->at(5) 			= -0.561160714285714/2.0;
                QuadW->at(6) 			= -0.561160714285714/2.0;
                QuadW->at(7) 			= -0.561160714285714/2.0;
                QuadW->at(8) 			= -0.561160714285714/2.0;
                QuadW->at(9) 			= 1.295089285714286/2.0;
                QuadW->at(10) 			= 1.295089285714286/2.0;
                QuadW->at(11) 			= 1.295089285714286/2.0;
                QuadW->at(12) 			= 1.295089285714286/2.0;
                QuadW->at(13) 			= 1.295089285714286/2.0;
                QuadW->at(14) 			= 1.295089285714286/2.0;
                QuadW->at(15) 			= 0.172767857142857/2.0;
                QuadW->at(16) 			= 0.172767857142857/2.0;
                QuadW->at(17) 			= 0.172767857142857/2.0;
                QuadW->at(18) 			= -1.354910714285714/2.0;
                QuadW->at(19) 			= -1.354910714285714/2.0;
                QuadW->at(20) 			= -1.354910714285714/2.0;
                QuadW->at(21) 			= -0.408482142857143/2.0;
                QuadW->at(22) 			= -0.408482142857143/2.0;
                QuadW->at(23) 			= -0.408482142857143/2.0;
                QuadW->at(24) 			= -0.408482142857143/2.0;
                QuadW->at(25) 			= -0.408482142857143/2.0;
                QuadW->at(26) 			= -0.408482142857143/2.0;
                QuadW->at(27) 			= 1.566517857142857/2.0;
                
                break;
                
            }
    }
    else if(dim==3){
        if (FEType.at(0)=='P') {
            if (Degree==2)
                Degree=3;
            if (Degree==4)
                Degree=5;

            TEUCHOS_TEST_FOR_EXCEPTION(Degree>6, std::runtime_error, "Tetrahedron quadrature rules only up to degree 6 available.");
            
            switch (Degree) {
                case 1:
                    QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(3,0.0)));
                    QuadW->resize(1);
                    QuadPts->at(0).at(0) 	= 0.25;
                    QuadPts->at(0).at(1) 	= 0.25;
                    QuadPts->at(0).at(2) 	= 0.25;
                    QuadW->at(0)			= 1/6.;
                    break;
                    
                case 3:
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
                    
                    QuadW->at(0)		 = -2./15.;
                    QuadW->at(1)		 = 3./40.;
                    QuadW->at(2)		 = 3./40.;
                    QuadW->at(3)		 = 3./40.;
                    QuadW->at(4)		 = 3./40.;
                    break;
                case 4:
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
                    
                    a = -.078933333333333;
                    b = .045733333333333;
                    c= .149333333333333;
                    
                    
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
                    
                case 5:
                    QuadPts.reset(new vec2D_dbl_Type(15,vec_dbl_Type(3,0.0)));
                    QuadW->resize(15);
                    a 	= 0.25;
                    b1 	= (7.+sqrt(15.))/34.;
                    b2 	= (7.-sqrt(15.))/34.;
                    c1 	= (13.-3.*sqrt(15.))/34.;
                    c2 	= (13.+3.*sqrt(15.))/34.;
                    d 	= (5.-sqrt(15.))/20.;
                    e 	= (5.+sqrt(15.))/20.;
                    
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
                    
                    
                    P1 	= (2665.-14.*sqrt(15.))/226800.;
                    P2 	= (2665.+14.*sqrt(15.))/226800.;
                    b	= 5./567.;
                    
                    QuadW->at(0) 			= 8./405.;
                    QuadW->at(1)            = P1;
                    QuadW->at(2)            = P1;
                    QuadW->at(3) 			= P1;
                    QuadW->at(4)            = P1;
                    
                    QuadW->at(5)            = P2;
                    QuadW->at(6)            = P2;
                    QuadW->at(7)            = P2;
                    QuadW->at(8)            = P2;
                    
                    QuadW->at(9) 			= b;
                    QuadW->at(10)           = b;
                    QuadW->at(11)           = b;
                    QuadW->at(12) 			= b;
                    QuadW->at(13)           = b;
                    QuadW->at(14)           = b;
                    
                    break;
                case 6: //Keast
                    QuadPts.reset(new vec2D_dbl_Type(24,vec_dbl_Type(3,0.0)));
                    QuadW->resize(24);
                    a = .356191386222545;
                    b = .214602871259152;
                    c = .877978124396166;
                    d = .040673958534611;
                    f = .032986329573173;
                    g = .322337890142276;
                    h = .269672331458316;
                    i = .063661001875018;
                    j = .603005664791649;
                    
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
                    QuadPts->at(23).at(2)   = j;
                    
                    a = .039922750258168;
                    b = .010077211055321;
                    c = .055357181543654;
                    d = .048214285714286;
                    
                    QuadW->at(0)    = a;
                    QuadW->at(1)    = a;
                    QuadW->at(2)    = a;
                    QuadW->at(3)    = a;
                    QuadW->at(4)    = b;
                    QuadW->at(5)    = b;
                    QuadW->at(6)    = b;
                    QuadW->at(7)    = b;
                    QuadW->at(8)    = c;
                    QuadW->at(9)    = c;
                    QuadW->at(10)   = c;
                    QuadW->at(11)   = c;
                    QuadW->at(12)   = d;
                    QuadW->at(13)   = d;
                    QuadW->at(14)   = d;
                    QuadW->at(15)   = d;
                    QuadW->at(16)   = d;
                    QuadW->at(17)   = d;
                    QuadW->at(18)   = d;
                    QuadW->at(19)   = d;
                    QuadW->at(20)   = d;
                    QuadW->at(21)   = d;
                    QuadW->at(22)   = d;
                    QuadW->at(23)   = d;
            }
        }
    }
    
}
}


