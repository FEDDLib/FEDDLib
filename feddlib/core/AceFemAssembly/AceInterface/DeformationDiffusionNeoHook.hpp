#ifndef DEFORMATIONDIFFUSIONNEOHOOK_hpp
#define DEFORMATIONDIFFUSIONNEOHOOK_hpp

#ifdef __cplusplus
extern "C"{
#endif// __cplusplus
/*!
 @brief AceGen-Fortran(FEAP)-C headers
 @author Sharan N Ramesh
 @version 1.0
 @copyright SH
 */
  
/*!
 \brief Fortran-C binding for computing the Tangent Stiffnes Matrix and Residuum. All Values must be passed by reference!
 @param[in] v Working vector for AceGen, size fixed to 1060
 @param[in] d [Input] Material Parameters size 2
 @param[in] ul [Input] Solution vector. Size ul(3,10) in Fortran(column major). Must be sent in as a vector in column major.
 @param[in] ul0 Currently unused, send a vector of size 30.
 @param[in] xl [Input] Position of nodes in reference coordinates. Size xl(3,10) in Fortran(column major). Must be sent in as a vector in column major form.
 @param[in] s [Output] Tangent Stiffness Matrix. Size s(30,30) in Fortran(column major). Must be sent in as a vector in column major form.
 @param[in] p [Output] Residuum vector, Size p(30).
 @param[in] ht History parameters at time n. Currently unused.
 @param[in] hp History parameters at time n+1. Currently unused.
 @param[in] deltat Time increment. Double value.
 */
void skr_DDNH(double* v, double* d, double* ul, double* ul0, double* xl, double* s, double* p, double* ht, double* hp, double* deltat);

#ifdef __cplusplus
}
#endif// __cplusplus
#endif// DEFORMATIONDIFFUSIONNEOHOOK_hpp
