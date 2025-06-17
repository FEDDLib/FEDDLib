#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef FEDD_HAVE_ACEGENINTERFACE
#include <aceinterface.h>
#endif

int main()
{

#ifdef FEDD_HAVE_ACEGENINTERFACE
    double positions[] = {0. , 0. ,0., 1., 0., 0., 0., 1., 0., 0., 0., 1., 0.5, 0., 0., 0.5, 0.5, 0., 0., 0.5, 0., 0., 0., 0.5, 0.5, 0., 0.5, 0., 0.5, 0.5};
    //{0.0, 0.5, 0.5, 0.5, 0.5, 1., 0.5, 0., 0.5, 0.5, 1., 0.5, 0.25, 0.5, 0.75, 0.5, 0.25, 0.75, 0.25, 0.25, 0.5, 0.25, 0.75, 0.5, 0.5, 0.75, 0.75, 0.5, 0.5, 0.5};
    double displacements[] = {0.1, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // find number of entries in positions
    int n1 = sizeof(positions) / sizeof(positions[0]);
    // find number of entries in displacements
    int n2 = sizeof(displacements) / sizeof(displacements[0]);
    // printf("n1=%d\n", n1);
    // printf("n2=%d\n", n2);
    char **a = getDataNames();
    // Print all entries in a
    // int i = 0;
    // while (a[i] != NULL)
    // {
    //     printf("%s\n", a[i]);
    //     i++;
    // }
    double concentrations[] = {0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double rates[] = {0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double accelerations[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double defd[] = {30e0, 0.12e1, 0.9e0, 0.3387e-1, -0.3387e-1, 50e0,
                     0.50247e0, 0.18745e0, 0.4e0, 0.2e0, 0.2e0, 0.134e0,
                     0.166e-2, 0.66e-4, 0.14636600000000002e3, 0.10097e-2, 0.9291e1, 0.2668e2,
                     0.15173775e3, 0.27566199999999996e1, 0.1152507e2, 0.127631e1, 0.308798e1, 0.3e0,
                     0.2e0, 0.5e0, 0.6e-4, 0e0, 1000e0, 1e0,
                     0e0};
    double history[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0};
    double deltaT = 1.;
    double time = 0.0;
    double subIterationTolerance = 1e-7;
    int integrationCode = 18;
    double historyUpdated[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0};
           
        
    double **stiffnessMatrixKuu = (double **)malloc(30 * sizeof(double *));
    for (int i = 0; i < 30; i++)
        stiffnessMatrixKuu[i] = (double *)malloc(30 * sizeof(double));
        
    double **stiffnessMatrixKuc = (double **)malloc(30 * sizeof(double *));
    for (int i = 0; i < 30; i++)
        stiffnessMatrixKuc[i] = (double *)malloc(10 * sizeof(double));
        
    double **stiffnessMatrixKcu = (double **)malloc(10 * sizeof(double *));
    for (int i = 0; i < 10; i++)
        stiffnessMatrixKcu[i] = (double *)malloc(30 * sizeof(double));
        
    double **stiffnessMatrixKcc = (double **)malloc(10 * sizeof(double *));
    for (int i = 0; i < 10; i++)
        stiffnessMatrixKcc[i] = (double *)malloc(10 * sizeof(double));

	double **massMatrixMc = (double **)malloc(10 * sizeof(double *));
    for (int i = 0; i < 10; i++)
        massMatrixMc[i] = (double *)malloc(10 * sizeof(double));
        
  //  double **stiffnessMatrixKcc = (double **)malloc(10 * sizeof(double *));
   // for (int i = 0; i < 10; i++)
   //     stiffnessMatrixKcc[i] = (double *)malloc(10 * sizeof(double));
    // find number of entries in concentrations
    int n3 = sizeof(concentrations) / sizeof(concentrations[0]);
    // find number of entries in rates
    int n4 = sizeof(rates) / sizeof(rates[0]);
    // find number of entries in accelerations
    int n5 = sizeof(accelerations) / sizeof(accelerations[0]);
    //printf("n3=%d\n", n3);
    //printf("n4=%d\n", n4);
    //printf("n5=%d\n", n5);

    getStiffnessMatrixKuu(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0],  &defd[0], &history[0],  subIterationTolerance, deltaT, time, integrationCode, &historyUpdated[0], stiffnessMatrixKuu);
    
    //getStiffnessMatrixKuu(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0], &defd[0], &history[0], subIterationTolerance, deltaT, time, integrationCode, &historyUpdated[0], stiffnessMatrixKuu);
    
   getStiffnessMatrixKcc(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0], &defd[0], &history[0], subIterationTolerance, deltaT, time, integrationCode, stiffnessMatrixKcc);


    
    	getStiffnessMatrixKuc(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0], &defd[0], &history[0], subIterationTolerance, deltaT, time, integrationCode,  stiffnessMatrixKuc);

    
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            printf("%f ", stiffnessMatrixKcc[i][j]);
        }
        printf("\n");
    }
    double *residuumRint = (double *)malloc(30 * sizeof(double));

     double *residuumRc = (double *)malloc(10 * sizeof(double));

	getResiduumVectorRc(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0],  &defd[0], &history[0], subIterationTolerance, deltaT, time,integrationCode, &historyUpdated[0], residuumRc);

	getResiduumVectorRint(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0],  &defd[0], &history[0], subIterationTolerance, deltaT, time,integrationCode, &historyUpdated[0], residuumRint);
    printf("\n");
    printf("Residuum");
    for (int j = 0; j < 10; j++)
    {
        printf("%f ", residuumRc[j]);
    }
    printf("\n");
    

    for (int i = 0; i < 30; i++)
    {
        free(stiffnessMatrixKuu[i]);
        free(stiffnessMatrixKuc[i]);
    }

    free(stiffnessMatrixKuu);
    free(stiffnessMatrixKuc);
    
    for (int i = 0; i < 10; i++)
    {
		free(stiffnessMatrixKcc[i]);
		free(stiffnessMatrixKcu[i]);

    }
    free(stiffnessMatrixKcc);
    free(stiffnessMatrixKcu);
#endif	


    return 0;
}

