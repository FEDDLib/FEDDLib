#ifndef MeshFileReader_hpp
#define MeshFileReader_hpp
/*!
 Declaration of MeshFileReader
 
 @brief  MeshFileReader
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
*/
using namespace std;
namespace FEDD {

void meshReadSize ( string mesh_filename, int &numNode, int &dim, int &numElement, int &orderElement, int &numSurface, int &orderSurface, int &numEdges, int& orderEdges );

template <typename T>
int NextEntry(string text, string delimiter, int fromPos, int &value, T&random){

    int pos = text.find(delimiter,fromPos);
    value = stoi(text.substr(fromPos,pos));
    
    return pos+1;
}

template <typename T>
int NextEntry(string text, string delimiter, int fromPos, double &value, T&random){

    int pos = text.find(delimiter,fromPos);
    value = stod(text.substr(fromPos,pos));
    
    return pos+1;
}

template <typename T>
void readEntity( ifstream& file, string& text, vector<T>& entity, vector<int>& entityFlags, const int& numEntity, const int& orderEntity, const string& delimiter, bool isNode = false, int dim = 3 ){
    
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
void meshReadData ( const string& mesh_filename, const string& type, const string& delimiter, const int &dim, const int& numEntity, const int& orderEntity, vector<T>& entity, vector<int>& entityFlag ){
    
    ifstream file;
    int lineLength;
    int pos;
    string text;
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
