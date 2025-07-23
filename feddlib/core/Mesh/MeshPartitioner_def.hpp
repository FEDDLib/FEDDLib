#ifndef MeshPartitioner_def_hpp
#define MeshPartitioner_def_hpp

#include "MeshPartitioner_decl.hpp"

/*!
 Definition of MeshPartitioner
 
 @brief  MeshPartitioner
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC, class LO, class GO, class NO>
MeshPartitioner<SC,LO,GO,NO>::MeshPartitioner()
{

}

template <class SC, class LO, class GO, class NO>
MeshPartitioner<SC,LO,GO,NO>::MeshPartitioner( DomainPtrArray_Type domains, ParameterListPtr_Type pL, std::string feType, int dimension )
{
    domains_ = domains;
    pList_ = pL;
    feType_ = feType;
    comm_ = domains_[0]->getComm();
    rankRanges_.resize( domains_.size() );
    dim_ = dimension;
}
    
template <class SC, class LO, class GO, class NO>
MeshPartitioner<SC,LO,GO,NO>::~MeshPartitioner(){

}


    
template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::readAndPartition( int volumeID, std::string meshUnit, bool convertToCM)
{
	if(volumeID != 10 ){
		if(this->comm_->getRank()==0){
			std::cout << " #### WARNING: The volumeID was set manually and is no longer 10. Please make sure your volumeID corresponds to the volumeID in your mesh file. #### " << std::endl;
		}
	}
    //Read
    std::string delimiter = pList_->get( "Delimiter", " " );
    for (int i=0; i<domains_.size(); i++) {
        std::string meshName = pList_->get( "Mesh " + std::to_string(i+1) + " Name", "noName" );
        TEUCHOS_TEST_FOR_EXCEPTION( meshName == "noName", std::runtime_error, "No mesh name given.");
        domains_[i]->initializeUnstructuredMesh( domains_[i]->getDimension(), "P1",volumeID, meshUnit, convertToCM ); //we only allow to read P1 meshes.
        domains_[i]->readMeshSize( meshName, delimiter );
    }
    
    this->determineRanks();

    for (int i=0; i<domains_.size(); i++){
        this->readAndPartitionMesh( i );
        domains_[i]->getMesh()->rankRange_ = rankRanges_[i];
    }
    
}
template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::determineRanks(){
    bool verbose ( comm_->getRank() == 0 );
    vec_int_Type fractions( domains_.size(), 0 );
    bool autoPartition = pList_->get( "Automatic partition", false );
    if ( autoPartition ) {
        //determine sum of elements and fractions based of domain contributions
        GO sumElements = 0;
        for (int i=0; i<domains_.size(); i++)
            sumElements += domains_[i]->getNumElementsGlobal();
            
        for (int i=0; i<fractions.size(); i++)
            fractions[i] = (domains_[i]->getNumElementsGlobal()*100) / sumElements;

        int diff = std::accumulate(fractions.begin(), fractions.end(), 0) - 100;
        auto iterator = fractions.begin();
        while (diff>0) {
            (*iterator)--;
            iterator++;
            diff--;
        }
        iterator = fractions.begin();
        while (diff<0) {
            (*iterator)++;
            iterator++;
            diff++;
        }

        this->determineRanksFromFractions( fractions );
        
        if (verbose) {
            std::cout << "\t --- ---------------- ---" << std::endl;
            std::cout << "\t --- Mesh Partitioner ---" << std::endl;
            std::cout << "\t --- ---------------- ---" << std::endl;
            std::cout << "\t --- Automatic partition for "<< comm_->getSize() <<" ranks" << std::endl;
            for (int i=0; i<domains_.size(); i++) {
                std::cout << "\t --- Fraction mesh "<<std::to_string(i+1) << " : " << fractions[i] <<
                            " of 100; rank range: " << std::get<0>( rankRanges_[i] )<< " to "<< std::get<1>( rankRanges_[i] ) << std::endl;
            }
        }
        
    }
    else if( autoPartition == false && pList_->get("Mesh 1 fraction ranks",-1) >= 0 ){
        for (int i=0; fractions.size(); i++)
            fractions[i] = pList_->get("Mesh " + std::to_string(i+1) + " fraction ranks", -1);
            
        TEUCHOS_TEST_FOR_EXCEPTION( std::accumulate(fractions.begin(), fractions.end(), 0) != 100, std::runtime_error, "Fractions do not sum up to 100!");
        this->determineRanksFromFractions( fractions );
        
        if (verbose) {
            std::cout << "\t --- ---------------- ---" << std::endl;
            std::cout << "\t --- Mesh Partitioner ---" << std::endl;
            std::cout << "\t --- ---------------- ---" << std::endl;
            std::cout << "\t --- Fraction partition for "<< comm_->getSize() <<" ranks" << std::endl;
            for (int i=0; i<domains_.size(); i++) {
                std::cout << "\t --- Fraction mesh "<<std::to_string(i+1) << " : " << fractions[i] <<
                " of 100; rank range: " << std::get<0>( rankRanges_[i] )<< " to "<< std::get<1>( rankRanges_[i] ) << std::endl;
            }
        }
    }
    else if( autoPartition == false && pList_->get("Mesh 1 fraction ranks",-1) < 0 && pList_->get("Mesh 1 number ranks",-1) > 0 ){
        int size = comm_->getSize();
        vec_int_Type numberRanks(domains_.size());
        for (int i=0; i<numberRanks.size(); i++)
            numberRanks[i] = pList_->get("Mesh " + std::to_string(i+1) + " number ranks",0);
            
        TEUCHOS_TEST_FOR_EXCEPTION( std::accumulate(numberRanks.begin(), numberRanks.end(), 0) > size, std::runtime_error, "Too many ranks requested for mesh partition!");
        this->determineRanksFromNumberRanks( numberRanks );
        if (verbose) {
            std::cout << "\t --- ---------------- ---" << std::endl;
            std::cout << "\t --- Mesh Partitioner ---" << std::endl;
            std::cout << "\t --- ---------------- ---" << std::endl;
            std::cout << "\t --- Rank number partition for "<< comm_->getSize() <<" ranks" << std::endl;
            for (int i=0; i<domains_.size(); i++) {
                std::cout << "\t --- Rank range mesh "<<std::to_string(i+1) << " :" << std::get<0>( rankRanges_[i] )<< " to "<< std::get<1>( rankRanges_[i] ) << std::endl;
            }
        }
        
    }
    else{
        for (int i=0; i<domains_.size(); i++)
            rankRanges_[i] = std::make_tuple( 0, comm_->getSize()-1 );
        if (verbose) {
            std::cout << "\t --- ---------------- ---" << std::endl;
            std::cout << "\t --- Mesh Partitioner ---" << std::endl;
            std::cout << "\t --- ---------------- ---" << std::endl;
            std::cout << "\t --- Every mesh on every rank" << std::endl;
        }
    }
}

template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::determineRanksFromFractions( vec_int_Type& fractions ){
   
    int lowerRank = 0; int size = comm_->getSize();
    int upperRank = 0;
    for (int i=0; i<fractions.size(); i++){
        upperRank = lowerRank + fractions[i] / 100. * size - 1;
        if (upperRank<lowerRank)
            upperRank++;
        rankRanges_[i] = std::make_tuple( lowerRank, upperRank );
        if (size>1)
            lowerRank = upperRank + 1;
        else
            lowerRank = upperRank;
    }
    int startLoc = 0;
    while ( upperRank > size-1 ) {
        for (int i=startLoc; i<rankRanges_.size(); i++){
            if (i>0)
                std::get<0>( rankRanges_[i] )--;
            std::get<1>( rankRanges_[i] )--;
        }
        startLoc++;
        upperRank--;
    }
    startLoc = 0;
    while ( upperRank < size-1 ) {
        for (int i=startLoc; i<rankRanges_.size(); i++){
            if (i>0)
                std::get<0>( rankRanges_[i] )++;
            std::get<1>( rankRanges_[i] )++;
        }
        startLoc++;
        upperRank++;
    }
}
    
template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::determineRanksFromNumberRanks( vec_int_Type& numberRanks ){
    
    int lowerRank = 0; int size = comm_->getSize();
    int upperRank = 0;
    for (int i=0; i<numberRanks.size(); i++){
        upperRank = lowerRank + numberRanks[i] - 1;
        rankRanges_[i] = std::make_tuple( lowerRank, upperRank );
        lowerRank = upperRank + 1;
    }
    
    int startLoc = 0;
    while ( upperRank > size-1 ) {
        for (int i=startLoc; i<rankRanges_.size(); i++){
            if (i>0)
                std::get<0>( rankRanges_[i] )--;
            std::get<1>( rankRanges_[i] )--;
        }
        startLoc++;
        upperRank--;
    }
    startLoc = 0;
    while ( upperRank < size-1 ) {
        for (int i=startLoc; i<rankRanges_.size(); i++){
            if (i>0)
                std::get<0>( rankRanges_[i] )++;
            std::get<1>( rankRanges_[i] )++;
        }
        startLoc++;
        upperRank++;
    }
    
}

/// Reading and partioning of the mesh. Input File is .mesh. Reading is serial and at some point the mesh entities are distributed along the processors.

template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::readAndPartitionMesh( int meshNumber ){
            
    typedef Teuchos::OrdinalTraits<GO> OTGO;

    MeshUnstrPtr_Type meshUnstr = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domains_[meshNumber]->getMesh() );
    
	// Reading nodes
    meshUnstr->readMeshEntity("node");
    // We delete the point at this point. We only need the flags to determine surface elements. We will load them again later.
    meshUnstr->pointsRep_.reset();
    // Reading elements
    meshUnstr->readMeshEntity("element");
    // Reading surfaces
    meshUnstr->readMeshEntity("surface");
    // Reading line segments 
    meshUnstr->readMeshEntity("line");


    
    bool verbose ( comm_->getRank() == 0 );
    bool buildEdges = pList_->get("Build Edge List", true);
    bool buildSurfaces = pList_->get("Build Surface List", true);

	// Adding surface as subelement to elements
    if (buildSurfaces)
        this->setSurfacesToElements( meshNumber );
    else
        meshUnstr->deleteSurfaceElements();
    
	// Serially distributed elements
    ElementsPtr_Type elementsMesh = meshUnstr->getElementsC();
    
    // Setup Metis
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // METIS DOKU: https://www.lrz.de/services/software/mathematik/metis/metis_5_0.pdf

    /* 
    options[METIS_OPTION_NUMBERING]
    Used to indicate which numbering scheme is used for the adjacency structure of a graph or the element-
    node structure of a mesh. The possible values are:
     0 C-style numbering is assumed that starts from 0.
     1 Fortran-style numbering is assumed that starts from 1.
    */
    options[METIS_OPTION_NUMBERING] = 0;

    /*
    options[METIS_OPTION_SEED]
    Specifies the seed for the random number generator.
    */
    options[METIS_OPTION_SEED] = 666;

    /*
    ptions[METIS_OPTION_CONTIG]
    Specifies that the partitioning routines should try to produce partitions that are contiguous. Note that if the
    input graph is not connected this option is ignored.
    */
    options[METIS_OPTION_CONTIG] = pList_->get("Contiguous",false); //0: Does not force contiguous partitions; 1: Forces contiguous partitions.
    
    /*
    options[METIS_OPTION_NCUTS]
    Specifies the number of different partitionings that it will compute. The final partitioning is the one that
    achieves the best edgecut or communication volume. Default is 1.
    */
    options[METIS_OPTION_NCUTS] = pList_->get("NCUTS",1);

    /*
    options[METIS_OPTION_MINCONN]
    Specifies that the partitioning routines should try to minimize the maximum degree of the subdomain graph,
    i.e., the graph in which each partition is a node, and edges connect subdomains with a shared interface.
    */
    options[METIS_OPTION_MINCONN] = pList_->get("MINCONN",1); //0; // 1: Explicitly minimize the maximum connectivity.
    
    /* 
     options[METIS_OPTION_OBJTYPE]: Specifies the type of objective. Possible values are:
     METIS_OBJTYPE_CUT Edge-cut minimization.
     METIS_OBJTYPE_VOL Total communication volume minimization.
    */
    idx_t objtype = METIS_OBJTYPE_CUT;
    if( pList_->get("OBJTYPE","METIS_OBJTYPE_CUT") == "METIS_OBJTYPE_CUT")
        objtype =METIS_OBJTYPE_CUT;
    else if(pList_->get("OBJTYPE","METIS_OBJTYPE_CUT") == "METIS_OBJTYPE_VOL")
        objtype = METIS_OBJTYPE_VOL;
    options[METIS_OPTION_OBJTYPE] =objtype; // METIS_OBJTYPE_CUT;// or METIS_OBJTYPE_VOL
    
    /*
     options[METIS_OPTION_RTYPE]
     Determines the algorithm used for refinement. Possible values are:
     METIS_RTYPE_FM FM-based cut refinement.
     METIS_RTYPE_GREEDY Greedy-based cut and volume refinement.
     METIS_RTYPE_SEP2SIDED Two-sided node FM refinement.
     METIS_RTYPE_SEP1SIDED One-sided node FM refinement.
    */
    idx_t rtype = METIS_RTYPE_FM;
    if(pList_->get("RTYPE","METIS_RTYPE_FM")=="METIS_RTYPE_FM")
        rtype = METIS_RTYPE_FM;
    else if(pList_->get("RTYPE","METIS_RTYPE_FM")== "METIS_RTYPE_GREEDY")
         rtype = METIS_RTYPE_GREEDY;
    else if(pList_->get("RTYPE","METIS_RTYPE_FM") == "METIS_RTYPE_SEP2SIDED")
         rtype = METIS_RTYPE_SEP2SIDED;
    else if(pList_->get("RTYPE","METIS_RTYPE_FM") == "METIS_RTYPE_SEP1SIDED")
         rtype = METIS_RTYPE_SEP1SIDED;

    options[METIS_OPTION_RTYPE] = rtype; //pList_->get("RTYPE","METIS_RTYPE_FM"); // METIS_RTYPE_GREEDY;

    /*
    options[METIS_OPTION_NITER]
    Specifies the number of iterations for the refinement algorithms at each stage of the uncoarsening process.
    Default is 10.
    */
    options[METIS_OPTION_NITER] = pList_->get("NITER",50); // 50; // default is 10
    
    /*
     options[METIS_OPTION_CCORDER]
     Specifies if the connected components of the graph should first be identified and ordered separately.
    */
    options[METIS_OPTION_CCORDER] = pList_->get("CCORDER",1);

    /*options[METIS OPTION DBGLVL]
     Specifies the amount of progress/debugging information will be printed during the execution of the algo-
     rithms. The default value is 0 (no debugging/progress information). A non-zero value can be supplied that
     is obtained by a bit-wise OR of the following values
    */
    options[METIS_OPTION_DBGLVL] =  pList_->get("DBGLVL",0);
    
    idx_t ne = meshUnstr->getNumElementsGlobal(); // Global number of elements
    idx_t nn = meshUnstr->getNumGlobalNodes();	// Global number of nodes
    idx_t ned = meshUnstr->getEdgeElements()->numberElements(); // Global number of edges

        
    int dim = meshUnstr->getDimension();
    std::string FEType = domains_[meshNumber]->getFEType();

	// Setup for paritioning with metis
    vec_idx_Type eptr_vec(0); // Vector for local elements ptr (local is still global at this point)
    vec_idx_Type eind_vec(0); // Vector for local elements ids
    
    makeContinuousElements(elementsMesh, eind_vec, eptr_vec);

    idx_t *eptr = &eptr_vec.at(0);
    idx_t *eind = &eind_vec.at(0);
    
    idx_t ncommon;
    int orderSurface;
    if (dim==2) {
        if (FEType=="P1") {
            ncommon = 2;
        }
        else if(FEType=="P2"){
            ncommon = 3;
        }
    }
    else if (dim==3) {
        if (FEType=="P1") {
            ncommon = 3;
        }
        else if(FEType=="P2"){
            ncommon = 6;
        }
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Wrong Dimension.");
    
    idx_t objval = 0;
    vec_idx_Type epart(ne,-1);
    vec_idx_Type npart(nn,-1);

    // Partitioning elements with metis
    if (verbose)
        std::cout << "-- Start partitioning with Metis ... " << std::flush;
    
    {
        FEDD_TIMER_START(partitionTimer," : MeshPartitioner : Partition Elements");
        idx_t nparts = std::get<1>( rankRanges_[meshNumber] ) - std::get<0>( rankRanges_[meshNumber] ) + 1;
        if ( nparts > 1 ) {
            int rank = this->comm_->getRank();
            // upperRange - lowerRange +1
            idx_t returnCode = METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon, &nparts, NULL, options, &objval, &epart[0], &npart[0]);
            if ( verbose )
                std::cout << "\n--\t Metis return code: " << returnCode;
        }
        else{
            for (int i=0; i<ne; i++)
                epart[i] = 0;
        }
    }

    if (verbose){
        std::cout << "\n--\t objval: " << objval << std::endl;
        std::cout << "-- done!" << std::endl;
    }
    
    if (verbose)
        std::cout << "-- Set Elements ... " << std::flush;
    
    vec_GO_Type locepart(0);
    vec_GO_Type pointsRepIndices(0);
    // Global Edge IDs of local elements
    vec_GO_Type locedpart(0);

	// Getting global IDs of element's nodes
    for (int i=0; i<ne; i++) {
        if (epart[i] == comm_->getRank() - std::get<0>( rankRanges_[meshNumber] ) ){
            locepart.push_back(i);
            for (int j=eptr[i]; j<eptr[i+1]; j++)
                pointsRepIndices.push_back( eind[j] ); // Ids of element nodes, globalIDs
        }
    }
    // TODO KHo why erase the vectors here? eind points to the underlying array and is used later.
    eind_vec.erase(eind_vec.begin(), eind_vec.end());
    eptr_vec.erase(eptr_vec.begin(), eptr_vec.end());

    // Sorting ids with global and corresponding local values to create repeated map
    make_unique(pointsRepIndices);
    if (verbose)
        std::cout << "done!" << std::endl;
    
    //  Building repeated node map
    Teuchos::ArrayView<GO> pointsRepGlobMapping =  Teuchos::arrayViewFromVector( pointsRepIndices );
    meshUnstr->mapRepeated_.reset( new Map<LO,GO,NO>(OTGO::invalid(), pointsRepGlobMapping, 0, this->comm_) );
    MapConstPtr_Type mapRepeated = meshUnstr->mapRepeated_;

    // We keep the global elements if we want to build edges later. Otherwise they will be deleted
    ElementsPtr_Type elementsGlobal = Teuchos::rcp( new Elements_Type( *elementsMesh ) );

	// Resetting elements to add the corresponding local IDs instead of global ones
    meshUnstr->elementsC_.reset(new Elements ( FEType, dim ) );
    {
        Teuchos::ArrayView<GO> elementsGlobalMapping =  Teuchos::arrayViewFromVector( locepart );
        // elementsGlobalMapping -> elements per Processor

        meshUnstr->elementMap_.reset(new Map<LO,GO,NO>( (GO) -1, elementsGlobalMapping, 0, this->comm_) );
        
        {
            int localSurfaceCounter = 0;
            for (int i=0; i<locepart.size(); i++) {
                std::vector<int> tmpElement;
                for (int j=eptr[locepart.at(i)]; j<eptr[locepart.at(i)+1]; j++) {
                    //local indices
                    int index = mapRepeated->getLocalElement( (long long) eind[j] );
                    tmpElement.push_back(index);
                }
		        //std::sort(tmpElement.begin(), tmpElement.end());
                FiniteElement fe( tmpElement, elementsGlobal->getElement( locepart.at(i) ).getFlag()  );
                // convert global IDs of (old) globally owned subelements to local IDs
                if (buildSurfaces) {
                    FiniteElement feGlobalIDs = elementsGlobal->getElement( locepart.at(i) );
                    if (feGlobalIDs.subElementsInitialized()){
                        ElementsPtr_Type subEl = feGlobalIDs.getSubElements();
                        subEl->globalToLocalIDs( mapRepeated );
                        fe.setSubElements( subEl );
                    }
                }
                meshUnstr->elementsC_->addElement( fe );
            }
        }
    }

	// Next we distribute the coordinates and flags correctly
    
          
    meshUnstr->readMeshEntity("node"); // We reread the nodes, as they were deleted earlier
    
    if (verbose)
        std::cout << "-- Build Repeated Points Volume ... " << std::flush;
            
    // building the unique map
    meshUnstr->mapUnique_ = meshUnstr->mapRepeated_->buildUniqueMap( rankRanges_[meshNumber] );

    // free(epart);
    if (verbose)
        std::cout << "-- Building unique & repeated points ... " << std::flush;
    {
        vec2D_dbl_Type points = *meshUnstr->getPointsRepeated();
        vec_int_Type flags = *meshUnstr->getBCFlagRepeated();
        meshUnstr->pointsRep_.reset( new std::vector<std::vector<double> >( meshUnstr->mapRepeated_->getNodeNumElements(), std::vector<double>(dim,-1.) ) );
        meshUnstr->bcFlagRep_.reset( new std::vector<int> ( meshUnstr->mapRepeated_->getNodeNumElements(), 0 ) );
        
        int pointIDcont;
        for (int i=0; i<pointsRepIndices.size() ; i++) {
            pointIDcont = pointsRepIndices.at(i);
            for (int j=0; j<dim; j++)
                meshUnstr->pointsRep_->at(i).at(j) = points[pointIDcont][j];
            meshUnstr->bcFlagRep_->at(i) =  flags[pointIDcont];
        }
    }

	// Setting unique points and flags
    meshUnstr->pointsUni_.reset(new std::vector<std::vector<double> >( meshUnstr->mapUnique_->getNodeNumElements(), std::vector<double>(dim,-1.) ) );
    meshUnstr->bcFlagUni_.reset(new std::vector<int> ( meshUnstr->mapUnique_->getNodeNumElements(), 0) );
    GO indexGlobal;
    MapConstPtr_Type map = meshUnstr->getMapRepeated();
    vec2D_dbl_ptr_Type pointsRep = meshUnstr->pointsRep_;
    for (int i=0; i<meshUnstr->mapUnique_->getNodeNumElements() ; i++) {
        indexGlobal = meshUnstr->mapUnique_->getGlobalElement(i);
        for (int j=0; j<dim; j++) {
            meshUnstr->pointsUni_->at(i).at(j) = pointsRep->at( map->getLocalElement( indexGlobal) ).at(j);
        }
        meshUnstr->bcFlagUni_->at(i) = meshUnstr->bcFlagRep_->at( map->getLocalElement( indexGlobal) );
    }

	// Finally we build the edges. As part of the edge build involves nodes and elements,
	// they should be build last to avoid any local and global IDs mix up
	if (!buildEdges)
        elementsGlobal.reset();

    locepart.erase(locepart.begin(),locepart.end());
    if (verbose)
        std::cout << "done!" << std::endl;
        
    if (buildSurfaces){
        this->setEdgesToSurfaces( meshNumber ); // Adding edges as subelements in the 3D case. All dim-1-Subelements were already set
		}
    else
        meshUnstr->deleteSurfaceElements();
        
    if (buildEdges) {
        if (verbose)
            std::cout << "-- Build edge element list ... " << std::endl << std::flush;
                
        buildEdgeListParallel( meshUnstr, elementsGlobal );
        
        if (verbose)
            std::cout << std::endl << " done!"<< std::endl;
        
        MapConstPtr_Type elementMap =  meshUnstr->getElementMap();

        FEDD_TIMER_START(partitionEdgesTimer," : MeshPartitioner : Partition Edges");
        meshUnstr->getEdgeElements()->partitionEdges( elementMap, mapRepeated );
        FEDD_TIMER_STOP(partitionEdgesTimer);

        // edge global indices on different processors
        for( int i=0; i<meshUnstr->getEdgeElements()->numberElements() ; i++){
            locedpart.push_back(meshUnstr->getEdgeElements()->getGlobalID((LO) i));
        }

        // Setup for the EdgeMap
        Teuchos::ArrayView<GO> edgesGlobalMapping =  Teuchos::arrayViewFromVector( locedpart );
        meshUnstr->edgeMap_.reset(new Map<LO,GO,NO>( (GO) -1, edgesGlobalMapping, 0, this->comm_) );
    }

    
    if (verbose)
        std::cout << "done!" << std::endl;
    
    if (verbose)
        std::cout << "-- Partition interface ... " << std::flush;
    meshUnstr->partitionInterface();
    
    if (verbose)
        std::cout << "done!" << std::endl;
    
}


/// Function that (addionally) adds edges as subelements in the 3D case. This is relevant when the .mesh file also contains edge information 
/// in 3D. Otherwise edges are not set as subelements in 3D.
template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::setEdgesToSurfaces(int meshNumber){
    bool verbose ( comm_->getRank() == 0 );
    MeshUnstrPtr_Type meshUnstr = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domains_[meshNumber]->getMesh() );
    ElementsPtr_Type elementsMesh = meshUnstr->getElementsC();
    MapConstPtr_Type mapRepeated = meshUnstr->mapRepeated_;
    if (verbose)
        std::cout << "-- Set edges of surfaces of elements ... " << std::flush;
    
    FEDD_TIMER_START(surfacesTimer," : MeshPartitioner : Set Surfaces of Edge Elements");
    vec2D_int_Type localEdgeIDPermutation;
    setLocalSurfaceEdgeIndices( localEdgeIDPermutation, meshUnstr->getEdgeElementOrder() );
    
    int volumeID = meshUnstr->volumeID_;
    
    ElementsPtr_Type elements = meshUnstr->getElementsC();
    ElementsPtr_Type edgeElements = meshUnstr->getSurfaceEdgeElements();
    
    /* First, we convert the surface edge elements to a 2D array so we can use std::find.*/
    // Can we use/implement find for Elements_Type?
    vec2D_int_Type edgeElements_vec( edgeElements->numberElements() );
    vec_int_Type edgeElementsFlag_vec( edgeElements->numberElements() );
    for (int i=0; i<edgeElements_vec.size(); i++){
        vec_int_Type edge = edgeElements->getElement(i).getVectorNodeListNonConst();
        std::sort( edge.begin(), edge.end() );
        edgeElements_vec.at(i)  = edge;
        edgeElementsFlag_vec.at(i) = edgeElements->getElement(i).getFlag();
    }
    
    vec_int_ptr_Type flags =  meshUnstr->bcFlagRep_;
    int elementEdgeSurfaceCounter;
    for (int i=0; i<elements->numberElements(); i++) {
        elementEdgeSurfaceCounter = 0;
        bool mark = false;
        for (int j=0; j<elements->getElement( i ).size(); j++) {
            if ( flags->at( elements->getElement( i ).getNode( j ) ) < volumeID  )
                elementEdgeSurfaceCounter++;
        }
        if (elementEdgeSurfaceCounter >= meshUnstr->getEdgeElementOrder()){
            //We want to find all surfaces of element i and set the surfaces to the element
            findAndSetSurfaceEdges( edgeElements_vec,  edgeElementsFlag_vec, elements->getElement(i), localEdgeIDPermutation,  mapRepeated );
        }
    }
    if (verbose)
        std::cout << "done!" << std::endl;
}

/// Adding surfaces as subelements to the corresponding elements. This adresses surfaces in 3D or edges in 2D. 
/// If in the 3D case edges are also part of the .mesh file, they will be added as subelements later. 
template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::setSurfacesToElements(int meshNumber){
    
    bool verbose ( comm_->getRank() == 0 );
    MeshUnstrPtr_Type meshUnstr = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domains_[meshNumber]->getMesh() );
    ElementsPtr_Type elementsMesh = meshUnstr->getElementsC(); // Previously read Elements 

    if (verbose)
            std::cout << "-- Set surfaces of elements ... " << std::flush;
    
    FEDD_TIMER_START(surfacesTimer," : MeshPartitioner : Set Surfaces of Elements");
    
    vec2D_int_Type localSurfaceIDPermutation;
    // get permutations
    setLocalSurfaceIndices( localSurfaceIDPermutation, meshUnstr->getSurfaceElementOrder() );
    
    int volumeID = meshUnstr->volumeID_;
    ElementsPtr_Type surfaceElements = meshUnstr->getSurfaceElements();
    
    /* First, we convert the surface Elements to a 2D array so we can use std::find.
     and partition it at the same time. We use a simple linear partition. This is done to reduce
     times spend for std::find. The result is then communicated. We therefore use the unpartitoned elements*/
    // Can we use/implement find for Elements_Type?
    
    int size = this->comm_->getSize();
    
    LO numSurfaceEl = surfaceElements->numberElements();// / size;
    /*LO rest = surfaceElements->numberElements() % size;
    
    vec_GO_Type offsetVec(size);
    for (int i=0; i<size; i++) {
        offsetVec[i] = numSurfaceEl * i;
        if ( i<rest && i == this->comm_->getRank() ) {
            numSurfaceEl++;
            offsetVec[i]+=i;
        }
        else
            offsetVec[i]+=rest;
    }*/

    vec2D_int_Type surfElements_vec( numSurfaceEl );
    vec2D_int_Type surfElements_vec_sorted( numSurfaceEl );

    vec_int_Type surfElementsFlag_vec( numSurfaceEl );
    vec_GO_Type offsetVec(size);
    int offset = offsetVec[this->comm_->getRank()];

    for (int i=0; i<surfElements_vec.size(); i++){
        vec_int_Type surface = surfaceElements->getElement(i ).getVectorNodeListNonConst(); // surfaceElements->getElement(i + offset).getVectorNodeListNonConst();
        surfElements_vec.at(i)  = surface;
        std::sort( surface.begin(), surface.end() ); // We need to maintain a consistent numbering in the surface elements, so we use a sorted and unsorted vector
        surfElements_vec_sorted.at(i) = surface;
        surfElementsFlag_vec.at(i) = surfaceElements->getElement(i).getFlag(); // surfaceElements->getElement(i + offset).getFlag();

    }
 
    // Delete the surface elements. They will be added to the elements in the following loop.
    surfaceElements.reset();
    vec_int_ptr_Type flags =  meshUnstr->bcFlagRep_;
                
    int elementSurfaceCounter;
    int surfaceElOrder = meshUnstr->getSurfaceElementOrder();
    for (int i=0; i<elementsMesh->numberElements(); i++) {
        elementSurfaceCounter = 0;
        for (int j=0; j<elementsMesh->getElement( i ).size(); j++) {
            if ( flags->at( elementsMesh->getElement( i ).getNode( j ) ) < volumeID  )
                elementSurfaceCounter++;           
        }
        
        if (elementSurfaceCounter >= surfaceElOrder){
            FEDD_TIMER_START(findSurfacesTimer," : MeshPartitioner : Find and Set Surfaces");
            //We want to find all surfaces of element i and set the surfaces to the element
            this->findAndSetSurfacesPartitioned( surfElements_vec_sorted,  surfElements_vec, surfElementsFlag_vec, elementsMesh->getElement(i), localSurfaceIDPermutation, offsetVec, i );
        }
    }
    
    if (verbose)
        std::cout << "done!" << std::endl;
}

/// Function that sets edges (2D) or surfaces(3D) to the corresponding element. It determines for the specific element 'element' if it has a surface with a non
/// volumeflag that can be set as subelement 
template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::findAndSetSurfacesPartitioned( vec2D_int_Type& surfElements_vec, vec2D_int_Type& surfElements_vec_unsorted, vec_int_Type& surfElementsFlag_vec, FiniteElement& element, vec2D_int_Type& permutation , vec_GO_Type& linearSurfacePartitionOffset, int globalElID ){
    
    // In general we look through the different permutations the input element 'element' can have and if they correspond to a surface. 
	// The mesh's surface elements 'surfElements_vec' are then used to determine the corresponding surface
	// If found, the nodes are then used to build the subelement and the corresponding surfElementFlag is set.  
	// The Ids are global at this point, as the values are not distributed yet.

    int loc, id1Glob, id2Glob, id3Glob;
    int size = this->comm_->getSize();
    vec_int_Type locAll(size);
    if (dim_ == 2) {
        for (int j=0; j<permutation.size(); j++) {
            id1Glob = element.getNode( permutation.at(j).at(0) ) ;
            id2Glob = element.getNode( permutation.at(j).at(1) ) ;
            
            vec_int_Type tmpSurface(2);
            if (id1Glob > id2Glob){
                tmpSurface[0] = id2Glob;
                tmpSurface[1] = id1Glob;
            }
            else{
                tmpSurface[0] = id1Glob;
                tmpSurface[1] = id2Glob;
            }
            
            loc = searchInSurfaces( surfElements_vec, tmpSurface );
            
            Teuchos::gatherAll<int,int>( *this->comm_, 1, &loc, locAll.size(), &locAll[0] );

            int surfaceRank = -1;
            int counter = 0;
            while (surfaceRank<0 && counter<size) {
                if (locAll[counter] > -1)
                    surfaceRank = counter;
                counter++;
            }
            int surfFlag = -1;
            if (loc>-1)
                surfFlag = surfElementsFlag_vec[loc];
            
            if (surfaceRank>-1) {
                Teuchos::broadcast<int,int>(*this->comm_,surfaceRank,1,&loc);
                Teuchos::broadcast<int,int>(*this->comm_,surfaceRank,1,&surfFlag);
                
                FiniteElement feSurface( tmpSurface, surfFlag);
                if ( !element.subElementsInitialized() )
                    element.initializeSubElements( "P1", 1 ); // only P1 for now
                
                element.addSubElement( feSurface );
            }
        }
    }
    else if (dim_ == 3){
        for (int j=0; j<permutation.size(); j++) {
            
            id1Glob = element.getNode( permutation.at(j).at(0) ) ;
            id2Glob = element.getNode( permutation.at(j).at(1) ) ;
            id3Glob = element.getNode( permutation.at(j).at(2) ) ;
            
            vec_int_Type tmpSurface = {id1Glob , id2Glob , id3Glob};
            sort( tmpSurface.begin(), tmpSurface.end() );
            loc = searchInSurfaces( surfElements_vec, tmpSurface );
            /*Teuchos::gatherAll<int,int>( *this->comm_, 1, &loc, locAll.size(), &locAll[0] );
            
            int surfaceRank = -1;
            int counter = 0;
            while (surfaceRank<0 && counter<size) {
                if (locAll[counter] > -1)
                    surfaceRank = counter;
                counter++;
            }
            int surfFlag = -1;
            if (loc>-1)
                surfFlag = surfElementsFlag_vec[loc];
                        
            if (surfaceRank>-1) {*/
            if(loc > -1 ){
                //Teuchos::broadcast<int,int>(*this->comm_,surfaceRank,1,&loc);
                //Teuchos::broadcast<int,int>(*this->comm_,surfaceRank,1,&surfFlag);

                int surfFlag = surfElementsFlag_vec[loc];
                //cout << " Surfaces set to elements on Proc  " << this->comm_->getRank() << " "  << surfElements_vec_unsorted[loc][0] << " " << surfElements_vec_unsorted[loc][1] << " " << surfElements_vec_unsorted[loc][2] << endl; 
                FiniteElement feSurface( surfElements_vec_unsorted[loc], surfFlag);
                if ( !element.subElementsInitialized() )
                    element.initializeSubElements( "P1", 2 ); // only P1 for now
                
                element.addSubElement( feSurface );
            }
        }
    }
    
}


template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::buildEdgeListParallel( MeshUnstrPtr_Type mesh, ElementsPtr_Type elementsGlobal ){
    FEDD_TIMER_START(edgeListTimer," : MeshReader : Build Edge List");
    ElementsPtr_Type elements = mesh->getElementsC();
    
    TEUCHOS_TEST_FOR_EXCEPTION( elements->getFiniteElementType() != "P1", std::runtime_error ,"Unknown discretization for method buildEdgeList(...).");
    CommConstPtr_Type comm = mesh->getComm();
    bool verbose(comm->getRank()==0);
    
    MapConstPtr_Type repeatedMap =  mesh->getMapRepeated();
    // Building local edges with repeated node list
    vec2D_int_Type localEdgeIndices(0);
    setLocalEdgeIndices( localEdgeIndices );
    EdgeElementsPtr_Type edges = Teuchos::rcp( new EdgeElements_Type() );
    for (int i=0; i<elementsGlobal->numberElements(); i++) {
        for (int j=0; j<localEdgeIndices.size(); j++) {
            
            int id1 = elementsGlobal->getElement( i ).getNode( localEdgeIndices[j][0] );
            int id2 = elementsGlobal->getElement( i ).getNode( localEdgeIndices[j][1] );
            vec_int_Type edgeVec(2);
            if (id1<id2){
                edgeVec[0] = id1;
                edgeVec[1] = id2;
            }
            else{
                edgeVec[0] = id2;
                edgeVec[1] = id1;
            }
            
            FiniteElement edge( edgeVec );
            edges->addEdge( edge, i );
        }
    }
    // we do not need elementsGlobal anymore
    elementsGlobal.reset();
    
    vec2D_GO_Type combinedEdgeElements;
    FEDD_TIMER_START(edgeListUniqueTimer," : MeshReader : Make Edge List Unique");
    edges->sortUniqueAndSetGlobalIDs( combinedEdgeElements );
    FEDD_TIMER_STOP(edgeListUniqueTimer);
    //Next we need to communicate all edge information. This will not scale at all!
    
    edges->setElementsEdges( combinedEdgeElements );
    
    mesh->setEdgeElements( edges );
            
};

template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::setLocalEdgeIndices(vec2D_int_Type &localEdgeIndices ){
    if ( dim_ == 2 ) {
        localEdgeIndices.resize(3,vec_int_Type(2,-1));
        localEdgeIndices.at(0).at(0) = 0;
        localEdgeIndices.at(0).at(1) = 1;
        localEdgeIndices.at(1).at(0) = 0;
        localEdgeIndices.at(1).at(1) = 2;
        localEdgeIndices.at(2).at(0) = 1;
        localEdgeIndices.at(2).at(1) = 2;
    }
    else if( dim_ == 3) {
        localEdgeIndices.resize(6,vec_int_Type(2,-1));
        localEdgeIndices.at(0).at(0) = 0;
        localEdgeIndices.at(0).at(1) = 1;
        localEdgeIndices.at(1).at(0) = 0;
        localEdgeIndices.at(1).at(1) = 2;
        localEdgeIndices.at(2).at(0) = 1;
        localEdgeIndices.at(2).at(1) = 2;
        localEdgeIndices.at(3).at(0) = 0;
        localEdgeIndices.at(3).at(1) = 3;
        localEdgeIndices.at(4).at(0) = 1;
        localEdgeIndices.at(4).at(1) = 3;
        localEdgeIndices.at(5).at(0) = 2;
        localEdgeIndices.at(5).at(1) = 3;
        
    }
}
    
template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::makeContinuousElements(ElementsPtr_Type elements, vec_idx_Type& eind_vec, vec_idx_Type& eptr_vec ){
    
    int nodesPerElement = elements->nodesPerElement();

    int elcounter=0;
    for (int i=0; i<elements->numberElements(); i++) {
        for (int j=0; j<nodesPerElement; j++) {
            eind_vec.push_back( elements->getElement( i ).getNode( j ) );
        }
        eptr_vec.push_back(elcounter);
        elcounter += nodesPerElement;
    }
    eptr_vec.push_back(elcounter);
    
}
 

template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::setLocalSurfaceEdgeIndices( vec2D_int_Type &localSurfaceEdgeIndices, int edgesElementOrder ){
    
    if ( dim_ == 3 ) {
        
        if (edgesElementOrder == 2) { //P1
            localSurfaceEdgeIndices.resize(6,vec_int_Type(2,-1));
            localSurfaceEdgeIndices.at(0).at(0) = 0;
            localSurfaceEdgeIndices.at(0).at(1) = 1;
            localSurfaceEdgeIndices.at(1).at(0) = 0;
            localSurfaceEdgeIndices.at(1).at(1) = 2;
            localSurfaceEdgeIndices.at(2).at(0) = 0;
            localSurfaceEdgeIndices.at(2).at(1) = 3;
            localSurfaceEdgeIndices.at(3).at(0) = 1;
            localSurfaceEdgeIndices.at(3).at(1) = 2;
            localSurfaceEdgeIndices.at(4).at(0) = 1;
            localSurfaceEdgeIndices.at(4).at(1) = 3;
            localSurfaceEdgeIndices.at(5).at(0) = 2;
            localSurfaceEdgeIndices.at(5).at(1) = 3;
        }
    }
}

template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::findAndSetSurfaceEdges( vec2D_int_Type& edgeElements_vec, vec_int_Type& edgeElementsFlag_vec, FiniteElement& element, vec2D_int_Type& permutation, MapConstPtr_Type mapRepeated){
    
    int loc, id1Glob, id2Glob;
    if (dim_ == 3){
        for (int j=0; j<permutation.size(); j++) {           
            id1Glob = mapRepeated->getGlobalElement( element.getNode( permutation.at(j).at(0) ) );
            id2Glob = mapRepeated->getGlobalElement( element.getNode( permutation.at(j).at(1) ) );
            vec_int_Type tmpEdge(0);
            if (id2Glob > id1Glob)
                tmpEdge = {id1Glob , id2Glob};
            else
                tmpEdge = {id2Glob , id1Glob};
            
            loc = searchInSurfaces( edgeElements_vec, tmpEdge );
            
            if (loc>-1) {
                
                int id1 = element.getNode( permutation.at(j).at(0) );
                int id2 = element.getNode( permutation.at(j).at(1) );
                vec_int_Type tmpEdgeLocal(0);
                if (id2>id1)
                    tmpEdgeLocal = { id1 , id2 };
                else
                    tmpEdgeLocal = { id2 , id1 };
                
                // If no partition was performed, all information is still global at this point. We still use the function below and partition the mesh and surfaces later.
                FiniteElement feEdge( tmpEdgeLocal, edgeElementsFlag_vec[loc] );
                // In some cases an edge is the only part of the surface of an Element. In that case there does not exist a triangle subelement. 
                // We then have to initialize the edge as subelement.     
                // In very coarse meshes it is even possible that an interior element has multiple surface edges or nodes connected to the surface, which is why we might even set interior edges as subelements                  
                if ( !element.subElementsInitialized() ){
                    element.initializeSubElements( "P1", 1 ); // only P1 for now                
                    element.addSubElement( feEdge );
                }
                else{
                    ElementsPtr_Type surfaces = element.getSubElements();
                    if(surfaces->getDimension() == 2)   // We set the edge to the corresponding surface element(s)
                        surfaces->setToCorrectElement( feEdge ); // Case we have surface subelements
                    else{ // simply adding the edge as subelement
                        element.addSubElement( feEdge ); // Case we dont have surface subelements
                        // Comment: Theoretically edges as subelement are only truely relevant when we apply a neuman bc on an edge which is currently not the case. 
                        // This is just in case!
                    }
                }                                
            }
        }
    }
    
}

template <class SC, class LO, class GO, class NO>
int MeshPartitioner<SC,LO,GO,NO>::searchInSurfaces( vec2D_int_Type& surfaces, vec_int_Type searchSurface){
    
    int loc = -1;
    
    vec2D_int_Type::iterator it = find(surfaces.begin(),surfaces.end(), searchSurface);
    
    if (it!=surfaces.end())
        loc = distance(surfaces.begin(),it);
    
    return loc;
}
 
template <class SC, class LO, class GO, class NO>
void MeshPartitioner<SC,LO,GO,NO>::setLocalSurfaceIndices(vec2D_int_Type &localSurfaceIndices, int surfaceElementOrder ){
    
    if ( dim_ == 2 ) {
        
        if (surfaceElementOrder == 2) { //P1
            localSurfaceIndices.resize(3,vec_int_Type(3,-1));
            localSurfaceIndices.at(0).at(0) = 0;
            localSurfaceIndices.at(0).at(1) = 1;
            localSurfaceIndices.at(1).at(0) = 0;
            localSurfaceIndices.at(1).at(1) = 2;
            localSurfaceIndices.at(2).at(0) = 1;
            localSurfaceIndices.at(2).at(1) = 2;
        }
        else
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No permutation for this surface yet.");
    }
    else if ( dim_ == 3 ){
        if (surfaceElementOrder == 3) {
            localSurfaceIndices.resize(4,vec_int_Type(3,-1));
            localSurfaceIndices.at(0).at(0) = 0;
            localSurfaceIndices.at(0).at(1) = 1;
            localSurfaceIndices.at(0).at(2) = 2;
            localSurfaceIndices.at(1).at(0) = 0;
            localSurfaceIndices.at(1).at(1) = 1;
            localSurfaceIndices.at(1).at(2) = 3;
            localSurfaceIndices.at(2).at(0) = 1;
            localSurfaceIndices.at(2).at(1) = 2;
            localSurfaceIndices.at(2).at(2) = 3;
            localSurfaceIndices.at(3).at(0) = 0;
            localSurfaceIndices.at(3).at(1) = 2;
            localSurfaceIndices.at(3).at(2) = 3;
        }
        else
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No permutation for this surface yet.");
    }
}

    
}
#endif
