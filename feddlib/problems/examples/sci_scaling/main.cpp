#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/SCI.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include "feddlib/core/General/AceGenInterfaceCheck.hpp"
#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_StackedTimer.hpp>

void rhsDummy2D(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    return;
}

void rhsDummy(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void zeroBC(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;

    return;
}


void zeroDirichlet(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void reactionFunc(double* x, double* res, double* parameters){
	
    double m = 0.0;	
    res[0] = m * x[0];

}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void inflowChem(double* x, double* res, double t, const double* parameters)
{
	if(t>=parameters[0])
    	res[0] = 1.;
    else	
    	res[0] = 0.;
    return;
}


void rhsX(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    
    res[0] = parameters[1];
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void rhsY(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] =  parameters[1];
    res[2] = 0.;
    return;
}

void rhsZ(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = parameters[1];
    return;
}

// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStep (lambda)
// 3 : LoadStep end time
// 4 : Flag 

void rhsYZ(double* x, double* res, double* parameters){

    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];

  	res[0] =0.;
    res[1] =0.;
    res[2] =0.;

    if(parameters[0]+1.e-12 < TRamp)
        force = (parameters[0]+loadStepSize) * parameters[1] / TRamp ;
    else
        force = parameters[1];


    if(parameters[6] == 4  || parameters[6] == 5){
      	res[0] = force;
        res[1] = force;
        res[2] = force;
    }
    
    return;
}


// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStep (lambda)
// 3 : LoadStep end time
// 4 : Flag 

void rhsArtery(double* x, double* res, double* parameters){

    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];

  	res[0] =0.;
    res[1] =0.;
    res[2] =0.;
    
    if(parameters[0]+1.e-12 < TRamp)
        force = (parameters[0]+loadStepSize) * parameters[1] / TRamp ;
    else
        force = parameters[1];


    if(parameters[6] == 5){
      	res[0] = force;
        res[1] = force;
        res[2] = force;
    }
    
    return;
}

// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStepSize
// 3 : LoadStep end time
// 4 : Flag 

void rhsHeartBeatCube(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double lambda=0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double heartBeatStart = parameters[4];
    
	double a0    = 11.693284502463376;
	double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
		  0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
		  0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
	double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
		  -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
		  0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
		         
    double Q = 0.5*a0;
    

    double t_min = parameters[0] - fmod(parameters[0],1.0); //FlowConditions::t_start_unsteady;
    double t_max = t_min + 0.52; // One heartbeat lasts 1.0 second    
    double y = M_PI * ( 2.0*( parameters[0]-t_min ) / ( t_max - t_min ) -0.87  );
    
    for(int i=0; i< 20; i++)
        Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
    
    
    // Remove initial offset due to FFT
    Q -= 0.026039341343493;
    Q = (Q - 2.85489)/(7.96908-2.85489);
    
    bool Qtrue=false;
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0]+1.0e-10 < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))+1.e-10< 0.5)
    	lambda= 0.75;
    else{
        lambda = 0.75+0.25*Q;//*0.005329; // 0.775+0.125 * cos(4*M_PI*(parameters[0]));
        Qtrue = true; 
    } 
  
    double forceDirection = force/fabs(force);
    if(parameters[6]==5 || parameters[6]==4){
        res[0] =lambda*force;//+forceDirection*Q;
        res[1] =lambda*force;//+forceDirection*Q;
        res[2] =lambda*force;//+forceDirection*Q;        
    } 
    
   /* if(parameters[0]< heartBeatStart){
    	Q = 0.;
    }
    
    if(parameters[0]+1e-12 < TRamp)
        force = force * (parameters[0]+loadStepSize);
    
    if(parameters[6] == 5 || parameters[6] == 4){
     	res[0] = force+Q*0.005329;
        res[1] = force+Q*0.005329;
       	res[2] = force+Q*0.005329;
       	      	
    }*/
      
}
// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStepSize
// 3 : LoadStep end time
// 4 : Flag 
void rhsHeartBeatArtery(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double lambda=0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double heartBeatStart = parameters[4];
    
	double a0    = 11.693284502463376;
	double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
		  0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
		  0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
	double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
		  -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
		  0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
		         
    double Q = 0.5*a0;
    

    double t_min = parameters[0] - fmod(parameters[0],1.0); //FlowConditions::t_start_unsteady;
    double t_max = t_min + 0.5; // One heartbeat lasts 1.0 second    
    double y = M_PI * ( 2.0*( parameters[0]-t_min ) / ( t_max - t_min )-1.);// -0.87  );
    
    for(int i=0; i< 20; i++)
        Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
    
    
    // Remove initial offset due to FFT
    Q -= 0.026039341343493;
    Q = (Q - 2.85489)/(7.96908-2.85489);
    
    bool Qtrue=false;
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0]+1.0e-10 < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))+1.e-10< 0.5)
    	lambda= 0.75;
    else{
        lambda = 0.75+0.25*Q;//*0.005329; // 0.775+0.125 * cos(4*M_PI*(parameters[0]));
        Qtrue = true; 
    } 
  
    double forceDirection = force/fabs(force);
    if(parameters[6]==5){
        res[0] =lambda*force;//+forceDirection*Q;
        res[1] =lambda*force;//+forceDirection*Q;
        res[2] =lambda*force;//+forceDirection*Q;        
    } 
    
   /* if(parameters[0]< heartBeatStart){
    	Q = 0.;
    }
    
    if(parameters[0]+1e-12 < TRamp)
        force = force * (parameters[0]+loadStepSize);
    
    if(parameters[6] == 5 || parameters[6] == 4){
     	res[0] = force+Q*0.005329;
        res[1] = force+Q*0.005329;
       	res[2] = force+Q*0.005329;
       	      	
    }*/
      
}

void rhsHeartBeatArteryPhases(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double lambda=1.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double heartBeatStart = parameters[4];

    double t = parameters[0];
    double heartBeatStart1 = parameters[5];
    double heartBeatEnd1 = parameters[6];
    double heartBeatStart2 = parameters[7];
    double heartBeatEnd2 = parameters[8];

    if(parameters[0]+1e-12 < TRamp)
        lambda = 1*(parameters[0]+loadStepSize)/ TRamp;
    else if(t > heartBeatStart1 && t<heartBeatEnd1)
    {    
        double a0    = 11.693284502463376;
        double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
            0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
            0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
        double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
            -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
            0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
                    
        double Q = 0.5*a0;
        

        double t_min = t - fmod(t,1.0)+heartBeatStart1-std::floor(t)+0.505; ; //FlowConditions::t_start_unsteady;
        double t_max = t_min + 1.0; // One heartbeat lasts 1.0 second    
        double y = M_PI * ( 2.0*( t-t_min ) / ( t_max - t_min ) -1.0)  ;
        
        for(int i=0; i< 20; i++)
            Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
        
        
        // Remove initial offset due to FFT
        Q -= 0.026039341343493;
        Q = (Q - 2.85489)/(7.96908-2.85489);

        if( t < heartBeatStart1 + 0.5)
		    lambda = 0.8 + 0.2*cos(2*M_PI*t);
        else 
    	    lambda= 0.6 + 0.95*Q;
    }
    else if(t>=heartBeatEnd1-1e-08 && t<heartBeatStart2)
    {
        if(t < heartBeatEnd1 + 0.5)
            lambda =  0.6;
        else if(t < heartBeatEnd1 + 1.0)
            lambda =  0.6 + 1.2* 0.5 * ( ( 1 - cos( M_PI*t/0.5) ));
        else
           lambda =  1.8;
    }
    else if(t>=heartBeatStart2-1.e-8 )
    {
        double a0    = 11.693284502463376;
        double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
            0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
            0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
        double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
            -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
            0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
                    
        double Q = 0.5*a0;
        

        double t_min = t - fmod(t,1.0)+heartBeatStart2-std::floor(t)+0.5; ; //FlowConditions::t_start_unsteady;
        double t_max = t_min + 1.0; // One heartbeat lasts 1.0 second    
        double y = M_PI * ( 2.0*( t-t_min ) / ( t_max - t_min ) -1.0)  ;
        
        for(int i=0; i< 20; i++)
            Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
        
        
        // Remove initial offset due to FFT
        Q -= 0.026039341343493;
        Q = (Q - 2.85489)/(7.96908-2.85489);

        if( t < heartBeatStart2 + 0.5)
		    lambda = 1.50 + 0.30*cos(2*M_PI*t);
        else 
    	    lambda= 1.20 + 1.2*Q;
    }
    else
    {
        lambda = 1.0;
    }

    double forceDirection = force/fabs(force);
    if(parameters[10]==5){
        res[0] =lambda*force;//+forceDirection*Q;
        res[1] =lambda*force;//+forceDirection*Q;
        res[2] =lambda*force;//+forceDirection*Q;        
    } 
      
}

void rhsHeartBeatArteryPulse(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double lambda=0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double heartBeatStart = parameters[4];
    
	    
    bool Qtrue=false;
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0]+1.0e-10 < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0]+1.0e-10 < heartBeatStart + 1.0)
		lambda = 0.75;
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))+1.e-10> 0.6)
    	lambda= 0.75;
    else{ // ( parameters[0] >= heartBeatStart + 0.5){ [0,0.6]
        // Within one second the heart beat passes through the artery
        double t= parameters[0] - std::floor(parameters[0]); // [0.5,1] Intervall
        double z= x[2];
        double t_z = t-1./80.*z; //  x  = 0.1 / 8 = 0.0125

        if(t_z > 0. && t_z < 0.5){
            double a0    = 11.693284502463376;
            double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
                0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
                0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
            double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
                -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
                0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
                        
            double Q = 0.5*a0;
        

            double t_min = 0; //parameters[0] - fmod(parameters[0],1.0); //FlowConditions::t_start_unsteady;
            double t_max = t_min + 0.5; // One heartbeat lasts 1.0 second    
            double y = M_PI * ( 2 * ( t_z-t_min ) / ( t_max - t_min ) - 1 );
           
            
            for(int i=0; i< 20; i++)
                Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
            
            
            // Remove initial offset due to FFT
            Q -= 0.026039341343493;
            Q = (Q - 2.85489)/(7.96908-2.85489);
            if(Q < 0 )
                Q = 0.;

            lambda = 0.75+0.25*Q;

        }
        else
            lambda = 0.75;
        
        if(lambda<0.75+1.e-12)
            lambda=0.75;

    } 
  
    double forceDirection = force/fabs(force);
    if(parameters[6]==5){
        res[0] =lambda*force;//+forceDirection*Q;
        res[1] =lambda*force;//+forceDirection*Q;
        res[2] =lambda*force;//+forceDirection*Q;        
    } 
    
}

// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStepSize
// 3 : LoadStep end time
// 4 : Flag 
void rhsArteryPaperPulse(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double lambda=0.;
    double heartBeatStart = parameters[4];
    
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0] < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))< 0.5)
    	lambda= 0.75;
    else{
        double tinc = parameters[0] - std::floor(parameters[0]);
        double Q = -sin(1/16.*M_PI*x[2]-M_PI*(tinc-0.5)*3.0);
        if(Q < 0.+1.e-12){
            Q = 0.;
            lambda=0.75;
        }
        else{
            lambda =0.75+0.25*Q;//0.875 - 0.125
        }
    }
    

    if(parameters[6]==5){
        res[0] =lambda*force;
        res[1] =lambda*force;
        res[2] =lambda*force; 
        
       if(fabs(lambda*force)<0.75*0.016)
        std::cout << " ALARMAAAAA lamba=" << lambda << " force=" << force << " t= " << parameters[0] << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ----------------" << std::endl;
    }
      
}

void rhsArteryPaper(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double lambda=0.;
    double heartBeatStart = parameters[4];
    
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0] < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))< 0.5)
    	lambda= 0.75;
    else
        lambda = 0.875 - 0.125 * cos(4*M_PI*(parameters[0]));
     
 
    if(parameters[6]==5){
        res[0] =lambda*force;
        res[1] =lambda*force;
        res[2] =lambda*force; 
        
       
    }
      
}

void rhsCubePaper(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[2] = 0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double lambda=0.;
    
    double heartBeatStart = parameters[4];
    
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0] < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))< 0.5)
    	lambda= 0.75;
    else
        lambda = 0.875 - 0.125 * cos(4*M_PI*(parameters[0]));
     
     
    if(parameters[6] == 5 || parameters[6] == 4){
        res[0] =lambda*force;
        res[1] =lambda*force;
        res[2] =lambda*force; 
    }
    
    

            
}
void dummyFunc(double* x, double* res, double t, const double* parameters)
{
    return;
}


typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

using namespace FEDD;
using namespace Teuchos;
using namespace std;

int main(int argc, char *argv[])
{


    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef ExporterParaView<SC,LO,GO,NO> ExporterPV_Type;
    typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);
    if (!FEDD::aceGenInterfaceAvailable())
        return EXIT_FAILURE;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    std::string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
   
    std::string xmlProblemFile = "parametersProblemSCI.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");    
    
    std::string xmlProblemStructureFile = "parametersProblemStructure.xml";  
    myCLP.setOption("problemfileStructure",&xmlProblemStructureFile,".xml file with Inputparameters.");    
 
    std::string xmlSolverFileSCI = "parametersSolverSCI.xml"; 
    myCLP.setOption("solverfileSCI",&xmlSolverFileSCI,".xml file with Inputparameters.");
    
    std::string xmlPrecFileStructure = "parametersPrecStructure.xml";
    myCLP.setOption("precfileStructure",&xmlPrecFileStructure,".xml file with Inputparameters.");
    
    std::string xmlPrecFileChem = "parametersPrecChem.xml";
    myCLP.setOption("precfileChem",&xmlPrecFileChem,".xml file with Inputparameters.");
    
 	//string xmlBlockPrecFile = "parametersPrecBlock.xml";
    //myCLP.setOption("blockprecfile",&xmlBlockPrecFile,".xml file with Inputparameters.");
   
    //string xmlPrecFile = "parametersPrec.xml";
    //myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");

    std::string xmlPrecCEFile = "parametersPrecCE.xml";
    myCLP.setOption("precCEfile",&xmlPrecCEFile,".xml file with Inputparameters.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        mpiSession.~GlobalMPISession();
        return 0;
    }
	Teuchos::RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Structure-chemical interaction",true));
    bool verbose (comm->getRank() == 0);
    TimeMonitor::setStackedTimer(stackedTimer);
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
       
        ParameterListPtr_Type parameterListProblemStructure = Teuchos::getParametersFromXmlFile(xmlProblemStructureFile);
        
        ParameterListPtr_Type parameterListSolverSCI = Teuchos::getParametersFromXmlFile(xmlSolverFileSCI);

        ParameterListPtr_Type parameterListPrecStructure = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
        ParameterListPtr_Type parameterListPrecChem = Teuchos::getParametersFromXmlFile(xmlPrecFileChem);
  
        ParameterListPtr_Type parameterListPrec;


 		int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        std::string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","unstructured");
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);
        std::string      discType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");
        std::string precMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        int         n;

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;     

        bool chemistryExplicit_ =    parameterListAll->sublist("Parameter").get("Chemistry Explicit",false);

        if(chemistryExplicit_)
            parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecCEFile);
        else
            parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
       

        parameterListAll->setParameters(*parameterListSolverSCI);
        parameterListAll->setParameters(*parameterListPrec);
                    
        /*else if(!precMethod.compare("Teko"))
            parameterListAll->setParameters(*parameterListPrecTeko);
        else if(precMethod == "Diagonal" || precMethod == "Triangular")
            parameterListAll->setParameters(*parameterListPrecBlock);*/
		parameterListAll->setParameters(*parameterListProblemStructure);
        
        ParameterListPtr_Type parameterListChemAll(new Teuchos::ParameterList(*parameterListPrecChem)) ;
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Chem") );
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter") );
        parameterListChemAll->setParameters(*parameterListSolverSCI);
        parameterListChemAll->setParameters(*parameterListPrecChem);
        parameterListChemAll->setParameters(*parameterListProblemStructure);

        
        ParameterListPtr_Type parameterListStructureAll(new Teuchos::ParameterList(*parameterListPrec));
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );
        parameterListStructureAll->setParameters(*parameterListPrec);
        parameterListStructureAll->setParameters(*parameterListProblem);
        parameterListStructureAll->setParameters(*parameterListProblemStructure);
		
        TimePtr_Type totalTime(TimeMonitor_Type::getNewCounter("FEDD - main - Total Time"));
        TimePtr_Type buildMesh(TimeMonitor_Type::getNewCounter("FEDD - main - Build Mesh"));

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);

        int size = comm->getSize() - numProcsCoarseSolve;

        // #####################
        // Mesh bauen und wahlen
        // #####################
    
        if (verbose)
        {
            std::cout << "###############################################" <<std::endl;
            std::cout << "############ Starting SCI  ... ################" <<std::endl;
            std::cout << "###############################################" <<std::endl;
        }

        DomainPtr_Type domainP1chem;
        DomainPtr_Type domainP1struct;
        DomainPtr_Type domainP2chem;
        DomainPtr_Type domainP2struct;
        
        
        DomainPtr_Type domainChem;
        DomainPtr_Type domainStructure;
        
        std::string bcType = parameterListAll->sublist("Parameter").get("BC Type","Cube");
        
        std::string rhsType = parameterListAll->sublist("Parameter").get("RHS Type","Constant");
    
        int minNumberSubdomains=1;
       
        if (!meshType.compare("structured")) {
		    TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured mesh.");
		    /*if (dim == 2) {
		        n = (int) (std::pow( size/minNumberSubdomains ,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
		        std::vector<double> x(2);
		        x[0]=0.0;    x[1]=0.0;
		        domainStructure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
		        domainChem.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
		    }
		    else if (dim == 3){*/
		        n = (int)(std::pow( size/minNumberSubdomains, 1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
		        std::vector<double> x(3);
		        x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
		        domainStructure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
		        domainChem.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
		    //}
		    // "Square5Element" and not "Square": the 6-element subcube decomposition
		    // gives tetrahedra with a negative Jacobian determinant, which the AceGen
		    // SCI elements cannot assemble on.
		    domainStructure->buildMesh( 3,"Square5Element", dim, discType, n, m, numProcsCoarseSolve);
		    domainChem->buildMesh( 3,"Square5Element", dim, discType, n, m, numProcsCoarseSolve);
		}
        else if (!meshType.compare("unstructured")) {
        
            domainP1chem.reset( new Domain_Type( comm, dim ) );
		    domainP1struct.reset( new Domain_Type( comm, dim ) );
		    domainP2chem.reset( new Domain_Type( comm, dim ) );
		    domainP2struct.reset( new Domain_Type( comm, dim ) );
            
            MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
            domainP1Array[0] = domainP1struct;
            
            ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );                    
            
            pListPartitioner->set("Build Edge List",true);
		    pListPartitioner->set("Build Surface List",true);
		                    
		    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
		    
		    int volumeID=10;
		    if(bcType=="Artery" || bcType == "Artery Full" || bcType == "Artery Plaque")
		    	volumeID = 15;
		    else if(bcType=="Realistic Artery 1" || bcType=="Realistic Artery 2" )
		    	volumeID = 21;
		    	
		    //partitionerP1.readAndPartition(volumeID);
            bool convertMesh = parameterListAll->sublist("Parameter").get("Convert Mesh",false);
            std::string unit = parameterListAll->sublist("Parameter").get("Mesh Unit","cm");

            if(convertMesh)
                partitionerP1.readAndPartition(15,unit , true ); // Convert it from mm to sm
            else
                partitionerP1.readAndPartition(15); 

            domainP1struct->exportElementFlags();
            domainP1struct->exportNodeFlags();
		    
            if (!discType.compare("P2")){
				domainP2chem->buildP2ofP1Domain( domainP1struct );
				domainP2struct->buildP2ofP1Domain( domainP1struct );

				domainChem = domainP2chem;   //domainP2chem;
				domainStructure = domainP2struct;   
			}        
			else{
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only P2 discretization allowed");                               

				domainStructure = domainP1struct;
				domainChem = domainP1struct;
			}
        }
        domainStructure->setDofs(dim);
        domainChem->setDofs(1);
        /*vec2D_dbl_ptr_Type pointsUni_ = domainStructure->getPointsUnique();
    
        for(int i=0; i<domainStructure->getMapUnique()->getNodeNumElements(); i++) {
            cout << " Node = " << domainStructure->getMapUnique()->getGlobalElement(i) << " " << (*pointsUni_)[i][0] << " " << (*pointsUni_)[i][1] << " " << (*pointsUni_)[i][2] << endl; 
        }
        for (int i=0; i<domainStructure->getElementsC()->numberElements(); i++){
            cout << " Element T = " << i << " ";
            for(int j=0; j< domainStructure->getElementsC()->getElement(i).getVectorNodeList().size(); j++)
                cout << domainStructure->getElementsC()->getElement(i).getVectorNodeList().at(j)  << " ";
            cout << endl;
        }
        domainStructure->getElementMap()->print();
        cout << " Num Elements " << domainStructure->getMesh()->getNumElementsGlobal() << endl;*/
        if (parameterListAll->sublist("General").get("ParaView export subdomains",false) ){
		   // ########################
		    // Flags check
		    // ########################

			Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaF(new ExporterParaView<SC,LO,GO,NO>());

			Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(domainStructure->getMapUnique()));
			vec_int_ptr_Type BCFlags = domainStructure->getBCFlagUnique();

			Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);
			for(int i=0; i< entries.size(); i++){
				entries[i] = BCFlags->at(i);
			}

			Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

			exParaF->setup("Flags", domainStructure->getMesh(), discType);

			exParaF->addVariable(exportSolutionConst, "Flags", "Scalar", 1,domainStructure->getMapUnique());

			exParaF->save(0.0);


            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaE(new ExporterParaView<SC,LO,GO,NO>());

			Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolutionE(new MultiVector<SC,LO,GO,NO>(domainStructure->getElementMap()));
			
			Teuchos::ArrayRCP< SC > entriesE  = exportSolutionE->getDataNonConst(0);
			for(int i=0; i<domainStructure->getElementsC()->numberElements(); i++){
				entriesE[i] = domainStructure->getElementsC()->getElement(i).getFlag();
			}

			Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConstE = exportSolutionE;

			exParaE->setup("Flags_Elements", domainStructure->getMesh(), "P0");

			exParaE->addVariable(exportSolutionConstE, "Flags_Elements", "Scalar", 1,domainStructure->getElementMap());

			exParaE->save(0.0);
		
	

            
            if (verbose)
                std::cout << "\t### Exporting subdomains ###\n";

            typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
            typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
            typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
            typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
            typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
            // Same subdomain for solid and chemistry, as they have same domain
            {
               MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domainStructure->getElementMap() ) );
                MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
                vecDecomposition->putScalar(comm->getRank()+1.);
                
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                
                exPara->setup( "subdomains_solid", domainStructure->getMesh(), "P0" );
                
                exPara->addVariable( vecDecompositionConst, "subdomains", "Scalar", 1, domainStructure->getElementMap());
                exPara->save(0.0);
                exPara->closeExporter();
            }
           
        
     
        }

        
    
        domainChem->setReferenceConfiguration();
        domainStructure->setReferenceConfiguration();

        MultiVectorPtr_Type nodes(domainStructure->getNodeListMV());

        // nodes->writeMM("nodes.mm");

        vec2D_dbl_Type diffusionTensor(dim,vec_dbl_Type(3));
        double D0 = parameterListAll->sublist("Parameter Diffusion").get("D0",1.);
        for(int i=0; i<dim; i++){
            diffusionTensor[0][0] =1;
            diffusionTensor[1][1] =1;
            diffusionTensor[2][2] =1;

            if(i>0){
                diffusionTensor[i][i-1] = 0;
                diffusionTensor[i-1][i] = 0;
            }
            else
                diffusionTensor[i][i+1] = 0;				
        }
 
        
        Teuchos::RCP<SmallMatrix<int>> defTS;

        if(chemistryExplicit_){
            defTS.reset( new SmallMatrix<int> (1) );
            // Stucture
            (*defTS)[0][0] = 1;
        }
        else {
            defTS.reset( new SmallMatrix<int> (2) );
            // Stucture
            (*defTS)[0][0] = 1;
            // Chem
            (*defTS)[1][1] = 1;
        }

       
			

        SCI<SC,LO,GO,NO> sci(domainStructure, discType,
                                domainChem, discType, diffusionTensor, reactionFunc,
                                parameterListStructureAll,
                                parameterListChemAll,
                                parameterListAll,
                                defTS);
        
        sci.info();
        
            
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) ); 
            
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryChem( new BCBuilder<SC,LO,GO,NO>( ) ); 
        
    
        // Struktur-RW
        
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStructure( new BCBuilder<SC,LO,GO,NO>( ) );

        if(dim == 2)
        {
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only 3D Test available");                               
        }
        else if(dim == 3 && bcType=="Cube")
        {

            bcFactory->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_X", dim);
            bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Y", dim);
            bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
            
            bcFactory->addBC(zeroDirichlet3D, 0, 0, domainStructure, "Dirichlet", dim);
            bcFactory->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_X_Y", dim);
            bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Y_Z", dim);
            bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_X_Z", dim);
            
            bcFactoryStructure->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_X", dim);
            bcFactoryStructure->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Y", dim);
            bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
            
            bcFactoryStructure->addBC(zeroDirichlet3D, 0, 0, domainStructure, "Dirichlet", dim);
            bcFactoryStructure->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_X_Y", dim);
            bcFactoryStructure->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Y_Z", dim);
            bcFactoryStructure->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_X_Z", dim);

             bcFactoryStructure->addBC(zeroDirichlet3D, 99, 0, domainStructure, "Dirichlet", dim);


        }
        else if(dim==3 && bcType=="Artery"){
        
			bcFactory->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_Y", dim);
			bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_X", dim);
			bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 4, 0, domainStructure, "Dirichlet_Z", dim);
			
			bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Z", dim);

			bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Y_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_X_Z", dim);

			bcFactory->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_X", dim);
			bcFactory->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Y", dim);

			bcFactory->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet_Y_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 12, 0, domainStructure, "Dirichlet_X_Z", dim);


			bcFactoryStructure->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_Y", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_X", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 4, 0, domainStructure, "Dirichlet_Z", dim);
			
			bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Z", dim);


			bcFactoryStructure->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Y_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_X_Z", dim);

			bcFactoryStructure->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_X", dim);

			bcFactoryStructure->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Y", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet_Y_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 12, 0, domainStructure, "Dirichlet_X_Z", dim);
        
        }
        else if(dim==3 && bcType=="Artery Full"){
        
			bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
            bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dim);
			
			bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dim);

			

			bcFactoryStructure->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
			
            bcFactoryStructure->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dim);


			bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dim);

        }
        else if(dim==3 && bcType=="Artery Plaque"){
        
			bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
            bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dim);
            bcFactory->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Z", dim); // Plaque surface
			bcFactory->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet_Z", dim); // plaque surface
			
			bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dim);
			

			bcFactoryStructure->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
			
            bcFactoryStructure->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dim);
             bcFactoryStructure->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet_Z", dim);


			bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dim);

        }
         else if(dim==3 && bcType=="Artery Realistic Plaque"){
        
			bcFactory->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_Z", dim);
            bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dim);
            bcFactory->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Z", dim); // Plaque surface
			
			bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dim);
			

			bcFactoryStructure->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_Z", dim);
            bcFactoryStructure->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Z", dim);
            bcFactoryStructure->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Z", dim);

			bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_X_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Y_Z", dim);

        }
		else if(dim==3 && bcType=="Realistic Artery 1"){
		

            bcFactoryStructure->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 17, 0, domainStructure, "Dirichlet", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 19, 0, domainStructure, "Dirichlet", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 15, 0, domainStructure, "Dirichlet", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 18, 0, domainStructure, "Dirichlet", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 12, 0, domainStructure, "Dirichlet", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 16, 0, domainStructure, "Dirichlet", dim);	
			
			
		    bcFactory->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet", dim);
			bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet", dim);
			bcFactory->addBC(zeroDirichlet3D, 17, 0, domainStructure, "Dirichlet", dim);
			bcFactory->addBC(zeroDirichlet3D, 19, 0, domainStructure, "Dirichlet", dim);
			bcFactory->addBC(zeroDirichlet3D, 15, 0, domainStructure, "Dirichlet", dim);
			bcFactory->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet", dim);
			bcFactory->addBC(zeroDirichlet3D, 18, 0, domainStructure, "Dirichlet", dim);
			bcFactory->addBC(zeroDirichlet3D, 12, 0, domainStructure, "Dirichlet", dim);
			bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet", dim);
			bcFactory->addBC(zeroDirichlet3D, 16, 0, domainStructure, "Dirichlet", dim);	

			
        }
        else if(dim==3 && bcType=="Realistic Artery 2"){
			bcFactoryStructure->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 17, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 19, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 15, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 18, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 12, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 16, 0, domainStructure, "Dirichlet_Z", dim);	
			
			bcFactoryStructure->addBC(zeroDirichlet3D, 77, 0, domainStructure, "Dirichlet_X", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 78, 0, domainStructure, "Dirichlet_Y", dim);
			
		  
			bcFactory->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 17, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 19, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 15, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 18, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 12, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 16, 0, domainStructure, "Dirichlet_Z", dim);	
			
			bcFactory->addBC(zeroDirichlet3D, 77, 0, domainStructure, "Dirichlet_X", dim);
			bcFactory->addBC(zeroDirichlet3D, 78, 0, domainStructure, "Dirichlet_Y", dim);
			
        }
       
        // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
        if (!sci.problemStructure_.is_null())
            sci.problemStructure_->addBoundaries(bcFactoryStructure);
        else
            sci.problemStructureNonLin_->addBoundaries(bcFactoryStructure);
        
        // RHS dummy for structure
        if (dim==2) {
            if (!sci.problemStructure_.is_null())
                sci.problemStructure_->addRhsFunction( rhsX,0 );
            else
                sci.problemStructureNonLin_->addRhsFunction( rhsX,0 );
            
        }
        else if (dim==3) {
            
            if (!sci.problemStructure_.is_null()){
                if(bcType=="Cube"){
					if(rhsType=="Constant")
		    		 	sci.problemStructure_->addRhsFunction( rhsYZ,0 );
		    		if(rhsType=="Paper")
		    		 	sci.problemStructure_->addRhsFunction( rhsCubePaper,0 );
		    		if(rhsType=="Heart Beat")
		        		 	sci.problemStructure_->addRhsFunction( rhsHeartBeatCube,0 );
					
				}
				else if(bcType=="Artery" || bcType=="Realistic Artery 1" || bcType=="Realistic Artery 2" || bcType == "Artery Full"|| bcType == "Artery Plaque"  ){
					if(rhsType=="Constant")
		    		 	sci.problemStructure_->addRhsFunction( rhsArtery,0 );
					if(rhsType=="Paper")
            		 	sci.problemStructure_->addRhsFunction( rhsArteryPaper,0 );
            		if(rhsType=="Heart Beat")
            		 	sci.problemStructure_->addRhsFunction( rhsHeartBeatArtery,0 );
                    if(rhsType=="Heart Beat Phases")
            		 	sci.problemStructure_->addRhsFunction( rhsHeartBeatArteryPhases,0 );
                    if(rhsType=="Paper Pulse")
            		 	sci.problemStructure_->addRhsFunction( rhsArteryPaperPulse,0 );
                    if(rhsType=="Heart Beat Pulse")
            		 	sci.problemStructure_->addRhsFunction(rhsHeartBeatArteryPulse,0 );
				}
                     
                double force = parameterListAll->sublist("Parameter").get("Volume force",1.);
                if(bcType == "Artery Full") // The surface normals are opposite to usual
                    force = force * -1.;
                sci.problemStructure_->addParemeterRhs( force );
                double loadStep = parameterListAll->sublist("Parameter").get("Load Step Size",1.);
                double loadRampEnd= parameterListAll->sublist("Parameter").get("Load Ramp End",1.);
                sci.problemStructure_->addParemeterRhs( loadStep );
                sci.problemStructure_->addParemeterRhs( loadRampEnd );
                double heartBeatStart= parameterListAll->sublist("Parameter").get("Heart Beat Start",70.);
                sci.problemStructure_->addParemeterRhs( heartBeatStart );
                sci.problemStructure_->addParemeterRhs( parameterListProblem->sublist("Parameter").get("Heart Beat Start 1",1.) );
                sci.problemStructure_->addParemeterRhs( parameterListProblem->sublist("Parameter").get("Heart Beat End 1",2.) );
                sci.problemStructure_->addParemeterRhs( parameterListProblem->sublist("Parameter").get("Heart Beat Start 2",3.) );
                sci.problemStructure_->addParemeterRhs( parameterListProblem->sublist("Parameter").get("Heart Beat End 2",4.) );
                sci.problemStructure_->addParemeterRhs( 0. ); // degree of the load function in space: develop's surface integral reads it from the last parameter

            }
            else{             
                if(bcType=="Cube"){
					if(rhsType=="Constant")
		    		 	sci.problemStructureNonLin_->addRhsFunction( rhsYZ,0 );
		    		if(rhsType=="Paper")
		        		sci.problemStructureNonLin_->addRhsFunction( rhsCubePaper,0 );
		    		if(rhsType=="Heart Beat")
		        		sci.problemStructureNonLin_->addRhsFunction( rhsHeartBeatCube,0 );
					
				}
				else if(bcType=="Artery" || bcType=="Realistic Artery 1" || bcType=="Realistic Artery 2" || bcType == "Artery Full" || bcType == "Artery Plaque" || bcType == "Artery Realistic Plaque"  ){
					if(rhsType=="Constant")
		    		 	sci.problemStructureNonLin_->addRhsFunction( rhsArtery,0 );
					if(rhsType=="Paper")
            		 	sci.problemStructureNonLin_->addRhsFunction( rhsArteryPaper,0 );
            		if(rhsType=="Heart Beat")
            		 	sci.problemStructureNonLin_->addRhsFunction( rhsHeartBeatArtery,0 );
                    if(rhsType=="Heart Beat Phases")
            		 	sci.problemStructureNonLin_->addRhsFunction( rhsHeartBeatArteryPhases,0 );
                    if(rhsType=="Paper Pulse")
            		 	sci.problemStructureNonLin_->addRhsFunction(rhsArteryPaperPulse,0 );
                    if(rhsType=="Heart Beat Pulse")
            		 	sci.problemStructureNonLin_->addRhsFunction(rhsHeartBeatArteryPulse,0 );
				}
                
                double force = parameterListAll->sublist("Parameter").get("Volume force",1.);
                if(bcType == "Artery Full")
                    force = force * -1.;
                sci.problemStructureNonLin_->addParemeterRhs( force );
                double loadStep = parameterListAll->sublist("Parameter").get("Load Step Size",1.);
                double loadRampEnd= parameterListAll->sublist("Parameter").get("Load Ramp End",1.);
                sci.problemStructureNonLin_->addParemeterRhs( loadStep );
                sci.problemStructureNonLin_->addParemeterRhs( loadRampEnd );
                double heartBeatStart= parameterListAll->sublist("Parameter").get("Heart Beat Start",70.);
                sci.problemStructureNonLin_->addParemeterRhs( heartBeatStart );
                sci.problemStructureNonLin_->addParemeterRhs( parameterListProblem->sublist("Parameter").get("Heart Beat Start 1",1.) );
                sci.problemStructureNonLin_->addParemeterRhs( parameterListProblem->sublist("Parameter").get("Heart Beat End 1",2.) );
                sci.problemStructureNonLin_->addParemeterRhs( parameterListProblem->sublist("Parameter").get("Heart Beat Start 2",3.) );
                sci.problemStructureNonLin_->addParemeterRhs( parameterListProblem->sublist("Parameter").get("Heart Beat End 2",4.) );
                sci.problemStructureNonLin_->addParemeterRhs( 0. ); // degree of the load function in space: develop's surface integral reads it from the last parameter

            }
            

        }
        if (dim==2)
        {
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only 3D Test available");                               
                            
        }
        else if(dim==3 && bcType=="Cube")
        {

            std::vector<double> parameter_vec(1, parameterListAll->sublist("Parameter").get("Inflow Start Time",0.));
            bcFactory->addBC(inflowChem, 0, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
            bcFactory->addBC(inflowChem, 1, 1, domainChem, "Dirichlet", 1, parameter_vec); // inflow of Chem
            bcFactory->addBC(inflowChem, 7, 1, domainChem, "Dirichlet", 1,parameter_vec);            		
            //bcFactory->addBC(zeroDirichlet, 8, 1, domainChem, "Dirichlet", 1);
            bcFactory->addBC(inflowChem, 9, 1, domainChem, "Dirichlet", 1,parameter_vec);
            /*bcFactory->addBC(zeroDirichlet, 2, 1, domainChem, "Dirichlet", 1);
            bcFactory->addBC(zeroDirichlet, 3, 1, domainChem, "Dirichlet", 1);            
            bcFactory->addBC(zeroDirichlet, 4, 1, domainChem, "Dirichlet", 1);            
            bcFactory->addBC(zeroDirichlet, 5, 1, domainChem, "Dirichlet", 1);            
           // bcFactory->addBC(zeroDirichlet, 6, 1, domainChem, "Dirichlet", 1);            
            */
            
            bcFactoryChem->addBC(inflowChem, 0, 0, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
            bcFactoryChem->addBC(inflowChem, 1, 0, domainChem, "Dirichlet", 1, parameter_vec); // inflow of Chem
            bcFactoryChem->addBC(inflowChem, 7, 0, domainChem, "Dirichlet", 1,parameter_vec);            		
            bcFactoryChem->addBC(inflowChem, 9, 0, domainChem, "Dirichlet", 1,parameter_vec);
           /* bcFactoryChem->addBC(zeroDirichlet, 2, 0, domainChem, "Dirichlet", 1);
            bcFactoryChem->addBC(zeroDirichlet, 3, 0, domainChem, "Dirichlet", 1);            
            bcFactoryChem->addBC(zeroDirichlet, 4, 0, domainChem, "Dirichlet", 1);            
            bcFactoryChem->addBC(zeroDirichlet, 5, 0, domainChem, "Dirichlet", 1);            
            //bcFactoryChem->addBC(zeroDirichlet, 6, 0, domainChem, "Dirichlet", 1);
            bcFactoryChem->addBC(zeroDirichlet, 8, 0, domainChem, "Dirichlet", 1);
            
            */
        }
        else if(dim==3 && bcType=="Artery"){
           std::vector<double> parameter_vec(1, parameterListAll->sublist("Parameter").get("Inflow Start Time",0.));
           bcFactory->addBC(inflowChem, 5, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
		   bcFactory->addBC(inflowChem, 13, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
		   bcFactory->addBC(inflowChem, 14, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
		   bcFactory->addBC(inflowChem, 7, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
		   bcFactory->addBC(inflowChem, 10, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
		   
           bcFactoryChem->addBC(inflowChem, 5, 0, domainChem, "Dirichlet", 1,parameter_vec);
           bcFactoryChem->addBC(inflowChem, 13, 0, domainChem, "Dirichlet", 1,parameter_vec);
           bcFactoryChem->addBC(inflowChem, 14, 0, domainChem, "Dirichlet", 1,parameter_vec);
           bcFactoryChem->addBC(inflowChem, 7, 0, domainChem, "Dirichlet", 1,parameter_vec);
           bcFactoryChem->addBC(inflowChem, 10, 0, domainChem, "Dirichlet", 1,parameter_vec);
        }
        else if(dim==3 && bcType=="Realistic Artery"){
           std::vector<double> parameter_vec(1, parameterListAll->sublist("Parameter").get("Inflow Start Time",0.));
           
           bcFactory->addBC(inflowChem, 5, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
		   bcFactory->addBC(inflowChem, 15, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
		   bcFactory->addBC(inflowChem, 14, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
		   
           bcFactoryChem->addBC(inflowChem, 5, 0, domainChem, "Dirichlet", 1,parameter_vec);
           bcFactoryChem->addBC(inflowChem, 15, 0, domainChem, "Dirichlet", 1,parameter_vec);
           bcFactoryChem->addBC(inflowChem, 14, 0, domainChem, "Dirichlet", 1,parameter_vec);
           
        }
        else if(dim==3 && bcType=="Artery Full"){
           std::vector<double> parameter_vec(1, parameterListAll->sublist("Parameter").get("Inflow Start Time",0.));
           
            bcFactory->addBC(inflowChem, 5, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
            bcFactory->addBC(inflowChem, 8, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
            bcFactory->addBC(inflowChem, 9, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem

		   
           bcFactoryChem->addBC(inflowChem, 5, 0, domainChem, "Dirichlet", 1,parameter_vec);
           bcFactoryChem->addBC(inflowChem, 8, 0, domainChem, "Dirichlet", 1,parameter_vec);
           bcFactoryChem->addBC(inflowChem, 9, 0, domainChem, "Dirichlet", 1,parameter_vec);
         
        }
        else if(dim==3 && bcType=="Artery Realistic Plaque"){
            std::vector<double> parameter_vec(1, parameterListAll->sublist("Parameter").get("Inflow Start Time",0.));
           
            bcFactory->addBC(inflowChem, 11, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
         
            bcFactoryChem->addBC(inflowChem, 11, 0, domainChem, "Dirichlet", 1,parameter_vec);
        } 
        else if(dim==3 && bcType=="Artery Plaque"){
           std::vector<double> parameter_vec(1, parameterListAll->sublist("Parameter").get("Inflow Start Time",0.));
           
            bcFactory->addBC(inflowChem, 6, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem outer wall

           /*bcFactory->addBC(inflowChem, 5, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
            bcFactory->addBC(inflowChem, 8, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
            bcFactory->addBC(inflowChem, 9, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem

		    bcFactory->addBC(zeroDirichlet, 19, 1, domainChem, "Dirichlet", 1); // inflow of Chem
            bcFactory->addBC(zeroDirichlet, 20, 1, domainChem, "Dirichlet", 1); // inflow of Chem*/

           bcFactoryChem->addBC(inflowChem, 6, 0, domainChem, "Dirichlet", 1,parameter_vec);


           /*bcFactoryChem->addBC(inflowChem, 5, 0, domainChem, "Dirichlet", 1,parameter_vec);
           bcFactoryChem->addBC(inflowChem, 8, 0, domainChem, "Dirichlet", 1,parameter_vec);
           bcFactoryChem->addBC(inflowChem, 9, 0, domainChem, "Dirichlet", 1,parameter_vec);

            bcFactoryChem->addBC(zeroDirichlet, 19, 1, domainChem, "Dirichlet", 1); // inflow of Chem
            bcFactoryChem->addBC(zeroDirichlet, 20, 1, domainChem, "Dirichlet", 1); // inflow of Chem*/


         
        }


        // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
        // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
        sci.problemChem_->addBoundaries(bcFactoryChem);
        
          
        // #####################
        // Zeitintegration
        // #####################
        sci.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

        sci.initializeProblem();

        sci.initializeCE();
        // Matrizen assemblieren
        sci.assemble();
                    

                    
        DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

        // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
        // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
        daeTimeSolver.defineTimeStepping(*defTS);

        // Uebergebe das (nicht) lineare Problem
        daeTimeSolver.setProblem(sci);

        // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massmatrizen auf den Zeilen, welche in
        // defTS definiert worden sind.
        daeTimeSolver.setupTimeStepping();

        daeTimeSolver.advanceInTime();
    }
    TimeMonitor_Type::report(std::cout);
    stackedTimer->stop("Structure-chemical interaction");
	StackedTimer::OutputOptions options;
	options.output_fraction = options.output_histogram = options.output_minmax = true;
	stackedTimer->report((std::cout),comm,options);
	
    return(EXIT_SUCCESS);
}
   /*  Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());


            exportSolution.reset(new MultiVector<SC,LO,GO,NO>(domainStructure->getMapVecFieldUnique()));
            exportSolution->putScalar(0.0);
            
            exportSolutionConst = exportSolution;
            exPara->setup("BC COND", domainStructure->getMesh(), discType);
            
            exPara->addVariable(exportSolutionConst, "BC Cond", "Vector", dim, domainStructure->getMapUnique());

            vec_int_ptr_Type flags = domainStructure->getBCFlagUnique();
            vec2D_dbl_ptr_Type nodes = domainStructure->getPointsUnique();

            entries  = exportSolution->getDataNonConst(0);

            double TRamp = 1.;
            double dt = 0.02; //parameterListAll->sublist("Timestepping Parameter").get("dt",1.0);
            double tMax = 5.0; //parameterListAll->sublist("Timestepping Parameter").get("Final time",1.0);
            double force = parameterListProblem->sublist("Parameter").get("Volume force",10.);
            double loadStepSize = 0.01;
            double heartBeatStart = 2.0;
            double r=0.;
            vec_dbl_Type res(3);
            double a = 2.;
            double lambda;
            for(double t=0.; t < tMax ; t= t+dt){

                for(int i=0; i< nodes->size(); i++){

                    if(flags->at(i) == 5){

                        vec_dbl_Type x = nodes->at(i);
                       
                        double lambda=0.;
                        
                        if(t+1e-12 < TRamp)
                            lambda = 0.875*(t+loadStepSize)/ TRamp;
                        else if(t <= TRamp+1.e-12)
                            lambda = 0.875;
                        else if (t < heartBeatStart)
                            lambda = 0.875;
                        else if( t < heartBeatStart + 0.5)
                            lambda = 0.8125+0.0625*cos(2*M_PI*t);
                        else if( t >= heartBeatStart + 0.5 && (t - std::floor(t))< 0.5)
                            lambda= 0.75;
                        else{
                            double tinc = t - std::floor(t);
                            double Q = -sin(1/16.*M_PI*x[2]-M_PI*(tinc-0.5)*3.0);
                            if(Q< 0){
                                Q = 0.;
                                lambda=0.75;
                            }
                            else
                                 lambda =0.75+0.25*Q;//0.875 - 0.125
                        }
                

                        double forceDirection = force/fabs(force);
                        res[0] =lambda*force;
                        res[1] =lambda*force;//+forceDirection*Q;
                        res[2] =lambda*force;//+forceDirection*Q;        
                        
                        for(int d=0; d<dim ; d++)
                            entries[i*dim+d] = res[d];
                    }
                }
                exPara->save(t);

            }
      

    */