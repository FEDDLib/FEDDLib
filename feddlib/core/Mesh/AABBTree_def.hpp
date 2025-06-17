#ifndef AABBTree_def_hpp
#define AABBTree_def_hpp

#include "AABBTree_decl.hpp"

/*!
Definition of AABBTree

@brief  AABBTree
@author Viktor Grimm
@version 1.0
@copyright VG
*/

namespace FEDD{
    // Default Constructor
    template <class SC, class LO, class GO, class NO>
    AABBTree<SC, LO, GO, NO>::AABBTree():
    empty_(true),
    dim_(),
    numNodes_(),
    minXY_(),
    maxXY_(),
    parentChild_(),
    containedElements_()
    {
    }


    /*
    Creates a new AABBTree from given Element- and node-vector.
    Finds minimal AABB for each element and uses them to call
    AABBTree::createTreeFromRectangles
    */
    template <class SC, class LO, class GO, class NO> \
    void AABBTree<SC, LO, GO, NO>::createTreeFromElements(
        ElementsPtr_Type elements,
        vec2D_dbl_ptr_Type nodes,
        int numberObjects, // = 32,
        double longTol, // = 0.75,
        double volumeTol, // = 0.55,
        bool verbose // = false
    ){

        if (verbose){
            std::cout << "   Starting AABBTree::'createTreeFromElements'" << '\n';
        }

        std::cout << "  Starting to create an AABB-Tree" << '\n';
        int numElems = elements->numberElements();


        vec2D_dbl_ptr_Type triangleMinMax( // point_min_max per triangle
            new vec2D_dbl_Type(
                numElems,
                vec_dbl_Type(4, 0.0)
            )
        );

        vec_int_Type localNodes;
        int currentNode;
        double currentNodeX, currentNodeY;

        for (int elem=0; elem<numElems; elem++){
            localNodes = elements->getElement(elem).getVectorNodeList();
            // Initialize with first point
            currentNode = localNodes[0];
            currentNodeX = nodes->at(currentNode).at(0);
            currentNodeY = nodes->at(currentNode).at(1);

            triangleMinMax->at(elem).at(0) = currentNodeX;
            triangleMinMax->at(elem).at(1) = currentNodeY;

            triangleMinMax->at(elem).at(2) = currentNodeX;
            triangleMinMax->at(elem).at(3) = currentNodeY;
            for ( int node=1; node<3; node++ ){
                currentNode = localNodes[node];
                currentNodeX = nodes->at(currentNode).at(0);
                currentNodeY = nodes->at(currentNode).at(1);
                // Find min_x
                if ( currentNodeX < triangleMinMax->at(elem).at(0) )
                triangleMinMax->at(elem).at(0) = currentNodeX;
                // Find min_y
                if ( currentNodeY < triangleMinMax->at(elem).at(1) )
                triangleMinMax->at(elem).at(1) = currentNodeY;
                // Find max_x
                if ( currentNodeX > triangleMinMax->at(elem).at(2) )
                triangleMinMax->at(elem).at(2) = currentNodeX;
                // Find max_y
                if ( currentNodeY > triangleMinMax->at(elem).at(3) )
                triangleMinMax->at(elem).at(3) = currentNodeY;
            }
        }

        createTreeFromRectangles(
            triangleMinMax,
            numberObjects,
            longTol,
            volumeTol,
            verbose
        );

        if (verbose){
            std::cout << "   Ending AABBTree::'createTreeFromElements'" << '\n';
        }
    }


    /*
    Creates a new AABBTree from the given rectangle list
    Expects each rectangle to conform to one element
    Forms a binary search tree, such that each
    "leaf" node in the tree encloses a maximum number of re-
    ctangles. The tree is formed by recursively subdividing
    the bounding-box of the collection. At each division, a
    simple heuristic is used to determine a splitting axis
    and to position an axis-aligned splitting (hyper-)plane.
    The associated collection of rectangles is partitioned
    between two new child nodes. The dimensions of each node
    in the tree are selected to provide a minimal enclosure
    of the rectangles in its associated sub-tree. Tree nodes
    may overlap as a result.
    */
    template <class SC, class LO, class GO, class NO> \
    void AABBTree<SC, LO, GO, NO>::createTreeFromRectangles(
        vec2D_dbl_ptr_Type initMinMaxXY,
        int numberObjects, // = 32
        double longTol, // = 0.75
        double volumeTol, // = 0.55
        bool verbose // = false
    ){
        // initMinMaxXY:
        //  double-Vector with four columns and numElems rows
        //  containing min_x, min_y, max_x and max_y for each initial rectangle
        // numberObjects:
        //  Desired number of objects per node (bound population control)
        // longTol, volumeTol:
        //  bound tolerance, in respect to side-length along the splitting
        //  axis and volume of the rectangle, that define whether
        //  to split or not


        if (verbose){
            std::cout << "   Starting AABBTree::'createTreeFromRectangles'" << '\n';
        }
        // FIXME: Do some Error-Checking here

        //  Number of Elements
        int numberElements = initMinMaxXY->size();
        std::cout << "   Number of Elements = " << numberElements << '\n';

        //***********************************************************************//
        //
        // Start by defining a rectangle that contains all other rectangles and
        // initializing variables we need later
        //
        //***********************************************************************//

        // Allocate Workspace
        vec2D_dbl_ptr_Type recMin(
            new vec2D_dbl_Type(
                numberElements,
                vec_dbl_Type(2, 0.0)
            )
        );
        vec2D_dbl_ptr_Type recMax(
            new vec2D_dbl_Type(
                numberElements,
                vec_dbl_Type(2, 0.0)
            )
        );

        vec2D_int_ptr_Type parentChild
        (new vec2D_int_Type(
            numberElements,
            vec_int_Type(2, 0)
        )
    );

    std::vector<std::list<int> >  containedElements(
        numberElements,
        std::list<int>()
    );

    vec_int_Type stack(numberElements, 0);

    // find min and max coords
    vec_dbl_Type bigMinMax = initMinMaxXY->at(0);
    for(int elem=1; elem<numberElements; elem++){
        if (initMinMaxXY->at(elem).at(0) < bigMinMax[0])
        bigMinMax[0] = initMinMaxXY->at(elem).at(0);
        if (initMinMaxXY->at(elem).at(1) < bigMinMax[1])
        bigMinMax[1] = initMinMaxXY->at(elem).at(1);
        if (initMinMaxXY->at(elem).at(2) > bigMinMax[2])
        bigMinMax[2] = initMinMaxXY->at(elem).at(2);
        if (initMinMaxXY->at(elem).at(3) > bigMinMax[3])
        bigMinMax[3] = initMinMaxXY->at(elem).at(3);
    }

    // Make rectangles slightly larger
    vec_dbl_Type inflateBy(2, 0.0);
    inflateBy[0] = bigMinMax[2] -  bigMinMax[0];
    inflateBy[1] = bigMinMax[3] -  bigMinMax[1];
    for (int elem=0; elem<numberElements; elem++){
        initMinMaxXY->at(elem).at(0) -= 1e-13 * inflateBy[0];
        initMinMaxXY->at(elem).at(1) -= 1e-13 * inflateBy[1];
        initMinMaxXY->at(elem).at(2) += 1e-13 * inflateBy[0];
        initMinMaxXY->at(elem).at(3) += 1e-13 * inflateBy[1];
    }

    vec2D_dbl_ptr_Type recCenter( // rectangle center
        new vec2D_dbl_Type(
            numberElements,
            vec_dbl_Type(2, 0.0)
        )
    );

    vec2D_dbl_ptr_Type recLength( // rectangle length
        new vec2D_dbl_Type(
            numberElements,
            vec_dbl_Type(2, 0.0)
        )
    );

    // Calculate rectangle center and length
    for (int elem = 0; elem<numberElements; elem++){
        recCenter->at(elem).at(0) = initMinMaxXY->at(elem).at(0) + initMinMaxXY->at(elem).at(2);
        recCenter->at(elem).at(0) *= 0.5;
        recCenter->at(elem).at(1) = initMinMaxXY->at(elem).at(1) + initMinMaxXY->at(elem).at(3);
        recCenter->at(elem).at(1) *= 0.5;

        recLength->at(elem).at(0) = initMinMaxXY->at(elem).at(2) - initMinMaxXY->at(elem).at(0);
        recLength->at(elem).at(0) = initMinMaxXY->at(elem).at(3) - initMinMaxXY->at(elem).at(1);
    }

    // First Rectangle contains all other rectangles

    recMin->at(0).at(0) = bigMinMax[0] - 1e-13 * inflateBy[0];
    recMin->at(0).at(1) = bigMinMax[1] - 1e-13 * inflateBy[1];
    recMax->at(0).at(0) = bigMinMax[2] + 1e-13 * inflateBy[0];
    recMax->at(0).at(1) = bigMinMax[3] + 1e-13 * inflateBy[1];

    for (int elem = 0; elem<numberElements; elem++){
        containedElements[0].push_back(elem);
    }

    if (verbose){
        std::cout << "   First rectangle = " << recMin->at(0).at(0) << ", " << recMin->at(0).at(1)<< ", " << recMax->at(0).at(0) << ", " << recMax->at(0).at(1) << '\n';
    }


    //***********************************************************************//
    //
    // Main Loop: divide nodes untill all constraints are satisfied
    //
    //***********************************************************************//
    stack[0] = 0;
    int nodesOnStack = 0;
    int numberNodes = 1;
    int sonLeft, sonRight, parentNode, axis;
    int countCurrentElements, countShortRectangles, countLongRectangles;
    int countLeft, countRight;
    double lengthAxis, splitPosition, volumeParent, volumeLeft, volumeRight;
    std::list<int> currentElements;
    vec_dbl_Type dd(2, 0.0);
    vec_int_Type ia(2, 0);
    std::vector<bool> isLong(numberElements, false);
    std::vector<bool> isLeft(numberElements, false);
    if (verbose){
        std::cout << "##### Beginning Main Loop #####" << '\n';
    }
    while (nodesOnStack != -1){
        if (verbose){
            std::cout << "   nodesOnStack = " << nodesOnStack << '\n';
        }
        // pop node from stack
        parentNode = stack[nodesOnStack];
        nodesOnStack = nodesOnStack - 1;

        if (verbose){
            std::cout << "   main loop for node = " << parentNode << '\n';
        }

        // push child indexing
        sonLeft = numberNodes;
        sonRight = numberNodes + 1;

        if (verbose){
            std::cout << "   sonLeft = " << sonLeft << '\n';
            std::cout << "   sonRight = " << sonRight << '\n';
        }

        // set of rectangles in parent
        currentElements = containedElements[parentNode];
        countCurrentElements = currentElements.size();
        if (verbose){
            std::cout << "   count currentElements =  " << countCurrentElements << '\n';
        }
        // split plane on longest axis
        dd[0] = recMax->at(parentNode).at(0) - recMin->at(parentNode).at(0);
        dd[1] = recMax->at(parentNode).at(1) - recMin->at(parentNode).at(1);

        // sorting by which axis to split
        if (dd[0] < dd[1]){
            ia[0] = 0;
            ia[1] = 1;
        }
        else{
            ia[1] = 0;
            ia[0] = 1;
        }

        if (verbose){
            std::cout << "   lengthAxis[0] =  " << dd[0] << ", lengthAxis[1] = "<< dd[1] << '\n';
        }

        for (int id=1; id>-1; id--){
            // push rectangles to children
            axis = ia[id];
            lengthAxis = dd[id];

            countShortRectangles = 0;
            countLongRectangles = 0;
            for (auto const& elem : currentElements) {
                if(recLength->at(elem).at(axis) > longTol * lengthAxis){
                    // long rectangles
                    isLong[elem] = true;
                    countLongRectangles += 1;
                }
                else {
                    // short rectangles
                    isLong[elem] = false;
                    countShortRectangles += 1;
                }
            }
            if (countLongRectangles < 0.5*countShortRectangles
                && countLongRectangles < 0.5 * numberObjects){
                    break;
                }
            }
            if (verbose){
                std::cout << "   splitting on axis = " << axis << ", with length = " << lengthAxis << '\n';
                std::cout << "   countShortRectangles = " << countShortRectangles << '\n';
                std::cout << "   countLongRectangles = " << countLongRectangles << '\n';
            }

            if(countShortRectangles == 0){
                // The partition is empty, we are done
                continue;
            }

            // select the split position: take the mean of the set of
            // (non-"long") rectangle centres along axis AX
            splitPosition = 0.0;
            for (auto const& elem : currentElements){
                if (!isLong[elem]){
                    splitPosition = splitPosition + recCenter->at(elem).at(axis);
                }
            }
            splitPosition = splitPosition / countShortRectangles;

            if (verbose){
                std::cout << "   split position = " << splitPosition << '\n';
            }

            // partition elements based on centres
            // if isLeft == true, the elment goes to the left rectangle, if == false
            // the element goes to the right rectangle
            countLeft = 0;
            countRight = 0;
            for (auto const& elem : currentElements){
                if (!isLong[elem]){
                    if (recCenter->at(elem).at(axis) < splitPosition){
                        isLeft[elem] = true;
                        countLeft += 1;
                    }
                    else {
                        isLeft[elem] = false;
                        countRight +=1;
                    }
                }
            }

            if (verbose){
                std::cout << "   countLeft = " << countLeft << '\n';
                std::cout << "   countRight = " << countRight << '\n';
            }
            if (countLeft == 0 || countRight == 0) {
                // The partition is empty, we are done
                continue;
            }

            // Finalise node position and push elemts downwards
            //Initialize with previous rectangle, but in reverse
            recMin->at(sonLeft).at(0) = recMax->at(parentNode).at(0);
            recMin->at(sonLeft).at(1) = recMax->at(parentNode).at(1);
            recMax->at(sonLeft).at(0) = recMin->at(parentNode).at(0);
            recMax->at(sonLeft).at(1) = recMin->at(parentNode).at(1);

            recMin->at(sonRight).at(0) = recMax->at(parentNode).at(0);
            recMin->at(sonRight).at(1) = recMax->at(parentNode).at(0);
            recMax->at(sonRight).at(0) = recMin->at(parentNode).at(0);
            recMax->at(sonRight).at(1) = recMin->at(parentNode).at(0);
            containedElements[parentNode].clear();
            containedElements[sonLeft].clear();
            containedElements[sonRight].clear();
            for(auto elem : currentElements){
                if (isLong[elem]){
                    // long elements will not be pushed down
                    containedElements[parentNode].push_back(elem);
                }
                else{

                    if (isLeft[elem]){
                        // element is in left rectangle, push it down to sonLeft
                        containedElements[sonLeft].push_back(elem);
                        if (initMinMaxXY->at(elem).at(0) < recMin->at(sonLeft).at(0)){
                            recMin->at(sonLeft).at(0) = initMinMaxXY->at(elem).at(0);
                        }
                        if (initMinMaxXY->at(elem).at(1) < recMin->at(sonLeft).at(1)){
                            recMin->at(sonLeft).at(1) = initMinMaxXY->at(elem).at(1);
                        }
                        if (initMinMaxXY->at(elem).at(2) > recMax->at(sonLeft).at(0)){
                            recMax->at(sonLeft).at(0) = initMinMaxXY->at(elem).at(2);
                        }
                        if (initMinMaxXY->at(elem).at(3) > recMax->at(sonLeft).at(1)){
                            recMax->at(sonLeft).at(1) = initMinMaxXY->at(elem).at(3);
                        }
                    }
                    else {
                        // element is in right rectangle, push it down to rightSon
                        containedElements[sonRight].push_back(elem);
                        if (initMinMaxXY->at(elem).at(0) < recMin->at(sonRight).at(0)){
                            recMin->at(sonRight).at(0) = initMinMaxXY->at(elem).at(0);
                        }
                        if (initMinMaxXY->at(elem).at(1) < recMin->at(sonRight).at(1)){
                            recMin->at(sonRight).at(1) = initMinMaxXY->at(elem).at(1);
                        }
                        if (initMinMaxXY->at(elem).at(2) > recMax->at(sonRight).at(0)){
                            recMax->at(sonRight).at(0) = initMinMaxXY->at(elem).at(2);
                        }
                        if (initMinMaxXY->at(elem).at(3) > recMax->at(sonRight).at(1)){
                            recMax->at(sonRight).at(1) = initMinMaxXY->at(elem).at(3);
                        }
                    }
                }
            }

            if (verbose){
                std::cout << "   left rectangle (" <<  recMin->at(sonLeft).at(0) << ", " << recMin->at(sonLeft).at(1)<< "), (" << recMax->at(sonLeft).at(0) << ", " << recMax->at(sonLeft).at(1) << ")"<<'\n';
                std::cout << "   right rectangle (" <<  recMin->at(sonRight).at(0) << ", " << recMin->at(sonRight).at(1)<< "), (" << recMax->at(sonRight).at(0) << ", " << recMax->at(sonRight).at(1) << ")"<<'\n';
            }

            if (countCurrentElements <= numberObjects){
                // if we are already below our desired number of objects we will only
                // split if the volume constraint is not yet satisfied
                if (verbose){
                    std::cout << "   Checking for Volume Constraint (ie objects constraint satisfied)" << '\n';
                }
                volumeParent = \
                (recMax->at(parentNode).at(0) - recMin->at(parentNode).at(0))\
                * (recMax->at(parentNode).at(1) - recMin->at(parentNode).at(1));
                volumeLeft = \
                (recMax->at(sonLeft).at(0) - recMin->at(sonLeft).at(0))\
                * (recMax->at(sonLeft).at(1) - recMin->at(sonLeft).at(1));
                volumeRight = \
                (recMax->at(sonRight).at(0) - recMin->at(sonRight).at(0))\
                * (recMax->at(sonRight).at(1) - recMin->at(sonRight).at(1));
                if (volumeLeft + volumeRight < volumeTol * volumeParent){
                    // parent-child indexing
                    parentChild->at(sonLeft).at(0) = parentNode;
                    parentChild->at(sonRight).at(0) = parentNode;
                    parentChild->at(parentNode).at(1) = sonLeft;

                    stack[nodesOnStack + 1] = sonLeft;
                    stack[nodesOnStack + 2] = sonRight;
                    nodesOnStack = nodesOnStack + 2;
                    numberNodes = numberNodes + 2;
                }
                else{
                    // all Constraints satisfied, no partition needed
                    // pull al elemts back to the parentNode
                    containedElements[parentNode].splice(
                        containedElements[parentNode].begin(),
                        containedElements[sonLeft]
                    );
                    containedElements[parentNode].splice(
                        containedElements[parentNode].begin(),
                        containedElements[sonRight]
                    );
                }
            }
            else {
                // object constraint is not satisfied, so we split
                if (verbose){
                    std::cout <<  "   Objects constraint not satisfied, Setting indices" << '\n';
                }
                parentChild->at(sonLeft).at(0) = parentNode;
                parentChild->at(sonRight).at(0) = parentNode;
                parentChild->at(parentNode).at(1) = sonLeft;
                stack[nodesOnStack + 1] = sonLeft;
                stack[nodesOnStack + 2] = sonRight;
                nodesOnStack = nodesOnStack + 2;
                numberNodes = numberNodes + 2;
            }

            if (verbose){
                std::cout << "   nodesOnStack = " << nodesOnStack << '\n';
            }

            if (verbose){
                std::cout << "   " << '\n';
            }

        } // end MAIN-LOOP
        if (verbose){
            std::cout << "##### End Main Loop #####" << '\n';
        }

        std::cout << "   Number of nodes in Tree =  " << numberNodes << '\n';

        // trim allocation
        recMin->resize(numberNodes);
        recMax->resize(numberNodes);
        parentChild->resize(numberNodes);
        containedElements.resize(numberNodes);

        // Saving data in class-members
        this->dim_ = 2;
        this->numNodes_ = numberNodes;
        this->minXY_ = recMin;
        this->maxXY_ = recMax;
        this->parentChild_ = parentChild;
        this->containedElements_ = containedElements;
        this->empty_ = false;

        if (verbose){
            std::cout << "##### End AABBTree::'createTreeFromRectangles' #####" << '\n';
        }
    }


    /*
    Low level routine that returns the tree-to-item mapping for a collection
    of query points.
    returns:

    The tree-to-item mapping (i : list<int>)  where i corresponds to the
    i-th node and and list<int> contains the number of query points that
    intersect with the i-th node.

    The item-to-tree mapping (i : list<int>) where i corresponds to the
    i-th query point and list<int> contains the number of nodes that
    intersect with the i-th query point.
    */
    template <class SC, class LO, class GO, class NO> \
    std::tuple<std::map<int, std::list<int>>, std::map<int, std::list<int> > > \
    AABBTree<SC, LO, GO, NO>::scanTree(
        vec2D_dbl_ptr_Type queryPoints,
        bool verbose
    ){

        if (verbose){
            std::cout << "   Starting AABBTree::'scanTree'" << '\n';
        }
        // FIXME: Catch some errors where
        if (queryPoints->at(0).size() != this->dim_){
            std::cerr << "Query point dimension =" << queryPoints->at(0).size() << ", but expected dim = " << this->dim_ << '\n';
        }
        int numberQueryPoints = queryPoints->size();

        // Allocate Workspace
        std::map<int, std::list<int>> treeToItem;
        std::map<int, std::list<int>> itemToTree;


        //***********************************************************************//
        //
        // Start by calculating treeToItem
        //
        //***********************************************************************//

        vec_int_Type stack(
            this->numNodes_,
            0
        );
        std::vector<std::list<int> >  pointsInNode(
            this->numNodes_,
            std::list<int>()
        );

        int leftSon, rightSon;
        int nodesOnStack = 0;
        int parentNode, no=0;
        stack[0] = 0;
        for (int point = 0; point<queryPoints->size(); point++){
            pointsInNode[0].push_back(point);
        }
        while (nodesOnStack != -1){
            // pop node from stack
            parentNode = stack[nodesOnStack];
            nodesOnStack = nodesOnStack - 1;

            if (verbose){
                std::cout << "   main loop for node = " << parentNode << '\n';
            }

            if (!this->containedElements_[parentNode].empty()){
                // non-empty node contains items, push onto tree-item mapping
                treeToItem[parentNode] = pointsInNode[parentNode];
                no = no + 1;
                if (verbose){
                    std::cout << "   Node " << parentNode << " is not empty, pushing onto tree-item mapping" << '\n';
                }
            }

            if(this->parentChild_->at(parentNode).at(1) != 0){
                // partition amongst child nodes
                leftSon = parentChild_->at(parentNode).at(1);
                rightSon = parentChild_->at(parentNode).at(1) + 1;
                if(verbose){
                    std::cout << "leftSon = "<< leftSon << ", rightSon = " << rightSon <<'\n';
                }
                // partition of points
                for (auto point : pointsInNode[parentNode]){
                    if(isInRectangle(leftSon, queryPoints->at(point), false)){
                        // point is in leftSon
                        pointsInNode[leftSon].push_back(point);
                    }
                    if(isInRectangle(rightSon, queryPoints->at(point), false)){
                        // point is in rightSon
                        pointsInNode[rightSon].push_back(point);
                    }
                }
                if (!pointsInNode[leftSon].empty()){
                    // leftSon has queryPoints in it, push node to stack
                    nodesOnStack = nodesOnStack + 1;
                    stack[nodesOnStack] = leftSon;
                    if (verbose){
                    }
                }
                if (!pointsInNode[leftSon].empty()){
                    // rightSon has queryPoints in it, push node to stack
                    nodesOnStack = nodesOnStack + 1;
                    stack[nodesOnStack] = rightSon;
                }
            }
        }

        //***********************************************************************//
        //
        // Calculate itemToTree, which is the inverse of treeToItem
        //
        //***********************************************************************//

        for (auto keyValue: treeToItem){
            // iterate through treeToItem, i.e. 'for each node with points in it'
            for (auto value: keyValue.second){
                // iterate through points in node
                itemToTree[value].push_back(keyValue.first);
            }
        }

        if (verbose){
            std::cout << "   End AABBTree::'scanTree'" << '\n';
        }

        return std::make_tuple(treeToItem, itemToTree);
    } // End scanTree

    /*
    Is a given query point in a given rectangle?
    */
    template <class SC, class LO, class GO, class NO> \
    bool \
    AABBTree<SC, LO, GO, NO>::isInRectangle(
        int numberRectangle,
        vec_dbl_Type queryPoint,
        bool verbose
    ){

        bool result = true;
        if (verbose){
            std::cout << "   Start AABBTree::'isInRectangle'" << '\n';
        }

        vec_dbl_Type rectangleMin, rectangleMax;
        rectangleMin = this->minXY_->at(numberRectangle);
        rectangleMax = this->maxXY_->at(numberRectangle);

        for (int dim = 0; dim < this->dim_-1; dim ++){
            if(result \
                && queryPoint[dim] > rectangleMin[dim] \
                && queryPoint[dim] < rectangleMax[dim]){
                    result = true;
                }
                else{
                    result = false;
                }
            }
        if (verbose){
            std::cout << "   End AABBTree::'isInRectangle'" << '\n';
        }
        return result;
    }


    template <class SC, class LO, class GO, class NO> \
    bool AABBTree<SC, LO, GO, NO>::isEmpty(){
        return empty_;
    }

    template <class SC, class LO, class GO, class NO> \
    int AABBTree<SC, LO, GO, NO>::getDim(){
        return dim_;
    }

    template <class SC, class LO, class GO, class NO> \
    int AABBTree<SC, LO, GO, NO>::getNumNodes(){
        return numNodes_;
    }

    template <class SC, class LO, class GO, class NO> \
    std::list<int> AABBTree<SC, LO, GO, NO>::getElements(int node){
        return containedElements_[node];
    }


}


#endif //AABBTree_def_hpp
