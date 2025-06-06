#ifndef SMALLMATRIX_hpp
#define SMALLMATRIX_hpp

#include "feddlib/core/FEDDCore.hpp"

/*!
 Declaration of SmallMatrix
 
 @brief  SmallMatrix
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
     /*!
    \class SmallMatrix
    \brief This class represents a templated small Matrix of type T. Primarily created for 2x2 and 3x3 matrices. Consequently, most features are limited to a certain matrix size.

    */
template <class T>
class SmallMatrix{
    
public:
    typedef unsigned UN;
    typedef std::vector<T> Vec_T_Type;
    typedef const std::vector<T> Vec_T_Const_Type;
    typedef std::vector<Vec_T_Type> Vec2D_T_Type;
    
    SmallMatrix();
    
    SmallMatrix(int size);
    
    SmallMatrix(int size, T v);

    Vec_T_Type &operator[](int i);
    
    Vec_T_Const_Type &operator[](int i) const;
    
    SmallMatrix<T> &operator=(SmallMatrix<T> &sm);
    
    SmallMatrix<T> operator*(SmallMatrix<T> &other);
    
    SmallMatrix<T> operator+(SmallMatrix<T> &other);
    
    SmallMatrix<T> &operator*=(SmallMatrix<T> &other);
    
    SmallMatrix<T> &operator+=(SmallMatrix<T> &other);
    
    //T &operator[](int i);
    
    int multiply(SmallMatrix<T> &bMat, SmallMatrix<T> &cMat); //this*B=C
    
    int add(SmallMatrix<T> &bMat, SmallMatrix<T> &cMat); //this+B=C
    
    int innerProduct(SmallMatrix<T> &bMat, T &res);
    
    double innerProduct(SmallMatrix<T> &bMat);
    
    void trace(T &res);
    
    void scale(T scalar);
    
    int size() const;
    
    void resize(UN size);
    
    void print();
    
    double computeInverse(SmallMatrix<T> &inverse);
        
    double computeDet();
    
    double computeScaling();
private:
    Vec2D_T_Type     values_;
    int         size_;
};



#include "SmallMatrix.hpp"

template<class T>
SmallMatrix<T>::SmallMatrix():
values_(0){
    
}
template<class T>
SmallMatrix<T>::SmallMatrix(int size):
values_(size,Vec_T_Type(size)),
size_(size){
    
}
template<class T>
SmallMatrix<T>::SmallMatrix(int size, T v):
values_( size, Vec_T_Type( size, v ) ),
size_(size){
    
}
template<class T>
typename SmallMatrix<T>::Vec_T_Type &SmallMatrix<T>::operator[](int i) {
    
    return values_.at(i);
}
template<class T>
typename SmallMatrix<T>::Vec_T_Const_Type &SmallMatrix<T>::operator[](int i) const{
    
    return values_.at(i);
}
template<class T>
int SmallMatrix<T>::multiply(SmallMatrix<T> &bMat, SmallMatrix &cMat){
    
    if (size_==2) {
        cMat[0][0] = values_[0][0]*bMat[0][0] + values_[0][1]*bMat[1][0];
        cMat[0][1] = values_[0][0]*bMat[0][1] + values_[0][1]*bMat[1][1];
        cMat[1][0] = values_[1][0]*bMat[0][0] + values_[1][1]*bMat[1][0];
        cMat[1][1] = values_[1][0]*bMat[0][1] + values_[1][1]*bMat[1][1];
    }
    else if(size_==3){
        cMat[0][0] = values_[0][0]*bMat[0][0] + values_[0][1]*bMat[1][0] + values_[0][2]*bMat[2][0];
        cMat[0][1] = values_[0][0]*bMat[0][1] + values_[0][1]*bMat[1][1] + values_[0][2]*bMat[2][1];
        cMat[0][2] = values_[0][0]*bMat[0][2] + values_[0][1]*bMat[1][2] + values_[0][2]*bMat[2][2];
        
        cMat[1][0] = values_[1][0]*bMat[0][0] + values_[1][1]*bMat[1][0] + values_[1][2]*bMat[2][0];
        cMat[1][1] = values_[1][0]*bMat[0][1] + values_[1][1]*bMat[1][1] + values_[1][2]*bMat[2][1];
        cMat[1][2] = values_[1][0]*bMat[0][2] + values_[1][1]*bMat[1][2] + values_[1][2]*bMat[2][2];
        
        cMat[2][0] = values_[2][0]*bMat[0][0] + values_[2][1]*bMat[1][0] + values_[2][2]*bMat[2][0];
        cMat[2][1] = values_[2][0]*bMat[0][1] + values_[2][1]*bMat[1][1] + values_[2][2]*bMat[2][1];
        cMat[2][2] = values_[2][0]*bMat[0][2] + values_[2][1]*bMat[1][2] + values_[2][2]*bMat[2][2];
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "SmallMatrix wrong size!");
    
    return 0;
}

template<class T>
int SmallMatrix<T>::add(SmallMatrix<T> &bMat, SmallMatrix &cMat){
    
	for(int i=0; i< size_; i++){
		for(int j=0; j< size_;j++){
			cMat[i][j] =  values_[i][j]+bMat[i][j];
		}
	}

    /*if (size_==2) {
        cMat[0][0] = values_[0][0]+bMat[0][0];
        cMat[0][1] = values_[0][1]+bMat[0][1];
        cMat[1][0] = values_[1][0]+bMat[1][0];
        cMat[1][1] = values_[1][1]+bMat[1][1];
    }
    else if(size_==3){
        cMat[0][0] = values_[0][0]+bMat[0][0];
        cMat[0][1] = values_[0][1]+bMat[0][1];
        cMat[0][2] = values_[0][2]+bMat[0][2];
        
        cMat[1][0] = values_[1][0]+bMat[1][0];
        cMat[1][1] = values_[1][1]+bMat[1][1];
        cMat[1][2] = values_[1][2]+bMat[1][2];
        
        cMat[2][0] = values_[2][0]+bMat[2][0];
        cMat[2][1] = values_[2][1]+bMat[2][1];
        cMat[2][2] = values_[2][2]+bMat[2][2];
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "SmallMatrix wrong size!");*/
    
    return 0;
}


template<class T>
int SmallMatrix<T>::innerProduct(SmallMatrix<T> &bMat, T &res){
    res = 0.;
    if (size_==2) {
        res =   values_[0][0]*bMat[0][0] + values_[0][1]*bMat[0][1] +
        values_[1][0]*bMat[1][0] + values_[1][1]*bMat[1][1];
    }
    else if(size_==3){
        res =   values_[0][0]*bMat[0][0] + values_[0][1]*bMat[0][1] + values_[0][2]*bMat[0][2] +
        values_[1][0]*bMat[1][0] + values_[1][1]*bMat[1][1] + values_[1][2]*bMat[1][2] +
        values_[2][0]*bMat[2][0] + values_[2][1]*bMat[2][1] + values_[2][2]*bMat[2][2];
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "SmallMatrix wrong size!");
    
    return 0;
}

template<class T>
double SmallMatrix<T>::innerProduct(SmallMatrix<T> &bMat){
    
    if (size_==2) {
        return   values_[0][0]*bMat[0][0] + values_[0][1]*bMat[0][1] +
        values_[1][0]*bMat[1][0] + values_[1][1]*bMat[1][1];
    }
    else if(size_==3){
        return  values_[0][0]*bMat[0][0] + values_[0][1]*bMat[0][1] + values_[0][2]*bMat[0][2] +
        values_[1][0]*bMat[1][0] + values_[1][1]*bMat[1][1] + values_[1][2]*bMat[1][2] +
        values_[2][0]*bMat[2][0] + values_[2][1]*bMat[2][1] + values_[2][2]*bMat[2][2];
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "SmallMatrix wrong size!");
    
    return 0.;
}

template<class T>
void SmallMatrix<T>::trace(T &res)
{
    // In res steht das Resultat
    if(size_ == 2)
    {
        res = values_[0][0] + values_[1][1];
    }
    else if(size_ == 3)
    {
        res = values_[0][0] + values_[1][1] + values_[2][2];
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "SmallMatrix wrong size!");

    
}
    
template<class T>
void SmallMatrix<T>::scale(T scalar){
    
    for (int i=0; i<size_; i++) {
        for (int j=0; j<size_; j++) {
            values_[i][j] *= scalar;
        }
    }
    
}

template<class T>
int SmallMatrix<T>::size() const{
    return size_;
}

template<class T>
void SmallMatrix<T>::resize(UN size){
    for (int i=0; i<size_; i++)
        values_[i].resize(size);
    
    values_.resize(size,Vec_T_Type(size));
    size_ = size;
}

template<class T>
typename SmallMatrix<T>::SmallMatrix &SmallMatrix<T>::operator=(SmallMatrix<T> &sm){
    
    size_ = sm.size();
    values_.resize(size_);
    
    for (int i=0; i<sm.size(); i++) {
        values_.at(i) = sm[i];
    }
    
    return *this;
}

template<class T>
typename SmallMatrix<T>::SmallMatrix SmallMatrix<T>::operator*(SmallMatrix<T> &other){
    
    size_ = other.size();
    SmallMatrix<T> result(size_);
    
    this->multiply(other, result);
    
    return result;
}

template<class T>
typename SmallMatrix<T>::SmallMatrix SmallMatrix<T>::operator+(SmallMatrix<T> &other){
    
    size_ = other.size();
    SmallMatrix<T> result(size_);
    
    this->add(other, result);
    
    return result;
}

template<class T>
typename SmallMatrix<T>::SmallMatrix &SmallMatrix<T>::operator*=(SmallMatrix<T> &other){
    
    this->multiply(other, *this);
    
    return *this;
}

template<class T>
typename SmallMatrix<T>::SmallMatrix &SmallMatrix<T>::operator+=(SmallMatrix<T> &other){
    
    this->add(other, *this);
    
    return *this;
}

template<class T>
void SmallMatrix<T>::print(){
    
    for (int i=0; i<size_; i++) {
        std::cout << "row "<< i << " values:";
        for (int j=0; j<size_; j++) {
            std::cout << " "<< values_[i][j];
        }
        std::cout << std::endl;
    }
}

template<class T>
double SmallMatrix<T>::computeInverse(SmallMatrix<T> &inverse){
    
    T det = this->computeDet();

    inverse.resize(size_);
    
    if (size_==2) {
        inverse[0][0] = values_[1][1] / det;
        inverse[0][1] = (-values_[0][1]) / det;
        inverse[1][0] = (-values_[1][0]) / det;
        inverse[1][1] = values_[0][0] / det;
    }
    else if (size_==3) {
        inverse[0][0] = (values_[1][1] * values_[2][2] - values_[1][2] * values_[2][1]) / det;
        inverse[0][1] = (values_[0][2] * values_[2][1] - values_[0][1] * values_[2][2]) / det;
        inverse[0][2] = (values_[0][1] * values_[1][2] - values_[0][2] * values_[1][1]) / det;
        
        inverse[1][0] = (values_[1][2] * values_[2][0] - values_[1][0] * values_[2][2]) / det;
        inverse[1][1] = (values_[0][0] * values_[2][2] - values_[0][2] * values_[2][0]) / det;
        inverse[1][2] = (values_[0][2] * values_[1][0] - values_[0][0] * values_[1][2]) / det;
        
        inverse[2][0] = (values_[1][0] * values_[2][1] - values_[1][1] * values_[2][0]) / det;
        inverse[2][1] = (values_[0][1] * values_[2][0] - values_[0][0] * values_[2][1]) / det;
        inverse[2][2] = (values_[0][0] * values_[1][1] - values_[0][1] * values_[1][0]) / det;

    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Matrix size is not 2 or 3.");
    return det;
}

template<class T>
double SmallMatrix<T>::computeDet( ){
    
    double det;
    if (size_==2) {
        det = values_[0][0] * values_[1][1] - values_[1][0] * values_[0][1];
    }
    else if(size_==3){
        det =   values_[0][0] * values_[1][1] * values_[2][2] +
                values_[0][1] * values_[1][2] * values_[2][0] +
                values_[0][2] * values_[1][0] * values_[2][1] -
                values_[2][0] * values_[1][1] * values_[0][2] -
                values_[2][1] * values_[1][2] * values_[0][0] -
                values_[2][2] * values_[1][0] * values_[0][1] ;
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Matrix size is not 2 or 3.");


    return det;
}
    
template<class T>
double SmallMatrix<T>::computeScaling( ){
    
    double scaling;
    if (size_==2)
        scaling = std::sqrt( values_[0][0] * values_[0][0] + values_[1][0] * values_[1][0] );
    else if(size_==3){
        
        double c1 = ( values_[1][0] * values_[2][1] - values_[2][0] * values_[1][1] );
        double c2 = ( values_[2][0] * values_[0][1] - values_[0][0] * values_[2][1] );
        double c3 = ( values_[0][0] * values_[1][1] - values_[1][0] * values_[0][1] );
        
        scaling = std::sqrt( c1*c1 + c2*c2 + c3*c3 );
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Matrix size is not 2 or 3.");
    
    
    return scaling;
}
    
    //Pseudo-Inverse: (A^T A)^-1 * A^T
    
}
#endif
