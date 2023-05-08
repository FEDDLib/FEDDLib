#ifndef AABBTree_decl_hpp
#define AABBTree_decl_hpp


#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/FE/Elements.hpp"

#include <list>
#include <map>
#include <tuple>
/*!
Declaration of AABBTree

@brief  AABBTree
@author Viktor Grimm
@version 1.0
@copyright VG
*/

namespace FEDD{
    template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
    class AABBTree {
    private:
        // Private Variables

        // Whether the AABBTree is empty
        bool empty_;

        // Dimension of the rectangles. 2 and above are acceptable values
        int dim_;

        // Number of rectangles (i.e. nodes) stored in the tree
        int numNodes_;

        // double-Vector with two columns and numNodes rows
        // containing min_x and min_y for each rectangle
        vec2D_dbl_ptr_Type minXY_;

        // double-Vector with two columns and numNodes rows
        // containing max_x and max_y for each rectangle
        vec2D_dbl_ptr_Type maxXY_;

        // List-Vector with numNodes rows and one column,
        // containing the Elements stored in each node.
        std::vector<std::list<int> > containedElements_;

        // int-Vector with numNodes rows and two columns
        // containing parent/child pointers associated  with each node
        vec2D_int_ptr_Type parentChild_;

        // Private Methods

    public:
        typedef Elements Elements_Type;
        typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;
        // Public Variables
        // Public Methods
        // Default Constructor
        AABBTree();

        void createTreeFromElements(
            ElementsPtr_Type elems,
            vec2D_dbl_ptr_Type nodes,
            int numberObjects = 32,
            double longTol = 0.75,
            double volumeTol = 0.55,
            bool verbose = false
        );

        void createTreeFromRectangles(
            vec2D_dbl_ptr_Type init_min_maxXY,
            int numberObjects = 32,
            double longTol = 0.75,
            double volumeTol = 0.55,
            bool verbose = false
        );

        std::tuple<std::map<int, std::list<int>>, std::map<int, std::list<int> > > scanTree(
            vec2D_dbl_ptr_Type query_points,
            bool verbose
        );


        bool isInRectangle(
            int number_rectangle,
            vec_dbl_Type query_point,
            bool verbose
        );

        bool isEmpty();

        int getDim();

        int getNumNodes();

        // Returns the elements contained in node
        std::list<int> getElements(int node);

    };
}

#endif
