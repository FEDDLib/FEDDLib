#ifndef MeshFileReader_hpp
#define MeshFileReader_hpp
#include <string>
#include <vector>
#include <fstream>

/*!
 Declaration of MeshFileReader
 
 @brief  MeshFileReader
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
*/

namespace FEDD {

void meshReadSize ( std::string mesh_filename, int &numNode, int &dim, int &numElement, int &orderElement, int &numSurface, int &orderSurface, int &numEdges, int& orderEdges );

template <typename T>
int NextEntry(std::string text, std::string delimiter, int fromPos, int &value, T&random){

    int pos = text.find(delimiter,fromPos);
    value = std::stoi(text.substr(fromPos,pos));
    
    return pos+1;
}

template <typename T>
int NextEntry(std::string text, std::string delimiter, int fromPos, double &value, T&random){

    int pos = text.find(delimiter,fromPos);
    value = stod(text.substr(fromPos,pos));
    
    return pos+1;
}

template <typename T>
void readEntity( std::ifstream& file, std::string& text, std::vector<T>& entity, std::vector<int>& entityFlags, const int& numEntity, const int& orderEntity, const std::string& delimiter, bool isNode = false, int dim = 3 ){
    
    T random; // we use this so NextEntry can deduce the type T
    double dblValue;
    int intValue;
    int flag;
    int pos;
    int lineLength = text.length();
    text.erase(0,lineLength);
    getline ( file, text );
    lineLength = text.length();
    text.erase(0,lineLength);
    if (isNode) {
        for (int i=0; i<numEntity; i++) {
            pos = 0;
            getline ( file, text );
            lineLength = text.length();
            for (int j=0; j<orderEntity; j++) {
                pos = NextEntry(text,delimiter,pos,dblValue,random);
                if (j<dim)
                    entity[dim*i+j] = dblValue;
                
            }
            NextEntry(text,"\n",pos,flag,random);
            entityFlags[i] = flag;
            text.erase(0,lineLength);
        }
    }
    else{
        for (int i=0; i<numEntity; i++) {
            pos = 0;
            getline ( file, text );
            lineLength = text.length();
            for (int j=0; j<orderEntity; j++) {
                pos = NextEntry(text,delimiter,pos,intValue,random);
                entity[orderEntity*i+j] = intValue;
            }
            NextEntry(text,"\n",pos,flag,random);
            entityFlags[i] = flag;
            text.erase(0,lineLength);
        }
    }
}

template <typename T>
void meshReadData ( const std::string& mesh_filename, const std::string& type, const std::string& delimiter, const int &dim, const int& numEntity, const int& orderEntity, std::vector<T>& entity, std::vector<int>& entityFlag ){
    
    std::ifstream file;
    int lineLength;
    int pos;
    std::string text;
    file.open ( mesh_filename.c_str ( ) );
    for ( ; ; )
    {
        getline ( file, text );
        lineLength = text.length();

        if (lineLength>6) {
            if ( !text.compare(0,8,"Vertices") && type == "node" ) {
                readEntity( file, text, entity, entityFlag, numEntity, orderEntity, delimiter, true, dim);
                return;
            }
        }
        if (lineLength>4) {
            if ( !text.compare(0,5,"Edges") &&
                ( (dim == 2 && type == "surface") || (dim == 3 && type == "line") ) ) {
                readEntity( file, text, entity, entityFlag, numEntity, orderEntity, delimiter );
                return;
            }
        }
        if (lineLength>7) {
            if (!text.compare(0,9,"Triangles") &&
                ( (dim == 2 && type == "element") || (dim == 3 && type == "surface") ) ) {
                readEntity( file, text, entity, entityFlag, numEntity, orderEntity, delimiter );
                return;
            }
        }
        if (lineLength>9) {
            if (!text.compare(0,10,"Tetrahedra") &&
                ( dim == 3 && type == "element" ) ) {
                readEntity( file, text, entity, entityFlag, numEntity, orderEntity, delimiter );
                return;
            }
        }
        text.erase(0,lineLength);
        if ( file.eof ( ) )
            break;
    }
    
    file.close ( );
    return;
}

}
#endif
