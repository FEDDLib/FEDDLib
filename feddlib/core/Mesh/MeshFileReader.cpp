# include <iostream>
# include <fstream>
# include <vector>
# include "MeshFileReader.hpp"
/*!
 Definition of MeshFileReader
 
 @brief  This class reads data from a given .mesh-file
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
namespace FEDD {
void meshReadSize ( string mesh_filename, int &numNode, int &dim, int &numElement, int &orderElement, int &numSurface, int &orderSurface, int &numEdges, int& orderEdges ){
    
    ifstream file;
    int lineLength;
    string text;
    numNode = 0;
    dim = 0;
    bool allSizesFound = false;
    file.open ( mesh_filename.c_str ( ) );

    if ( ! file )
    {
        cerr << "\n";
        cerr << "  Could not open input file \"" << mesh_filename << "\"\n";
        exit ( 1 );
    }

    for ( ; ; )
    {
        getline ( file, text );
        lineLength = text.length();
        if (lineLength>9) {
            if ( !text.compare(0,9,"Dimension")) {
                dim = stoi(text.substr(10,10));
            }
        }
        if (lineLength>7) {
            if (!text.compare(0,8,"Vertices")) {
                text.erase(0,lineLength);
                getline ( file, text );
                lineLength = text.length();
                numNode = stoi(text);
                text.erase(0,lineLength);
            }
        }
        if (lineLength>4) {
            if (!text.compare(0,5,"Edges")) {
                text.erase(0,lineLength);
                getline ( file, text );
                lineLength = text.length();
                if (dim == 2)
                    numSurface = stoi(text);
                else if(dim == 3)
                    numEdges = stoi(text);
                
            }
        }
        if (lineLength>8) {
            if (!text.compare(0,9,"Triangles")) {
                text.erase(0,lineLength);
                getline ( file, text );
                lineLength = text.length();
                if (dim == 2) {
                    numElement = stoi(text);
                }
                else if(dim == 3){
                    numSurface = stoi(text);
                }
            }
        }
        if (lineLength>9) {
            if (!text.compare(0,10,"Tetrahedra")) {
                text.erase(0,lineLength);
                getline ( file, text );
                lineLength = text.length();
                if (dim == 3) {
                    numElement = stoi(text);
                }
            }
        }
        text.erase(0,lineLength);
    if ( file.eof ( ) || allSizesFound)
      break;

    }
    
    file.close ( );

    if (dim == 2) {
        orderElement = 3;
        orderSurface = 2;
    }
    else if (dim == 3){
        orderElement = 4;
        orderSurface = 3;
        orderEdges = 2;
    }
    return;
}

}
