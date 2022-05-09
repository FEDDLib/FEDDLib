#ifndef NEOHOOKQUADRATICTETS2_hpp
#define NEOHOOKQUADRATICTETS2_hpp

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
 */
void skr2(double* v, double* d, double* ul, double* ul0, double* xl, double* s, double* p, double* ht, double* hp);

/*!
\brief Fortran-C binding for computing Post-processing values. All Values must be passed by reference!
 @param[in] v Working vector for AceGen, size fixed to 1060
 @param[in] d [Input] Material Parameters size 2
 @param[in] ul [Input] Solution vector. Size ul(3,10) in Fortran(column major). Must be sent in as a vector in column major.
 @param[in] ul0 Currently unused, send a vector of size 30.
 @param[in] xl [Input] Position of nodes in reference coordinates. Size xl(3,10) in Fortran(column major). Must be sent in as a vector in column major form.
 @param[in] s [Output] Tangent Stiffness Matrix. Size s(30,30) in Fortran(column major). Must be sent in as a vector in column major form.
 @param[in] p [Output] Residuum vector, Size p(30).
 @param[in] ht History parameters at time n. Currently unused.
 @param[in] hp History parameters at time n+1. Currently unused.
 @param[in] sg Unused unknown
 @param[in] sg0 Unused unkown
 @param[in] sxd Unused unkown
 @param[in] gpost Gauss Point Post-processing data. Size gpost(64,21).
 @param[in] npost Nodal Post-processing data. Size npost(10,6).
 */
void spp2(double* v, double* d, double* ul, double* ul0, double* xl, double* s, double* p, double* ht, double* hp, double* sg, double* sg0, double* sxd, double* gpost, double* npost);
#ifdef __cplusplus
}
#endif// __cplusplus
#endif// NEOHOOKQUADRATICTETS2_hpp