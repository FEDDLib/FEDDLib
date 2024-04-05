#include "feddlib/core/IO/IO_ParameterList_decl.hpp"

int main(int argc, char *argv[]) {

    using namespace std;

    FEDD::IO::ParameterList parameterList("main","parameterlist.yaml","yaml");

    parameterList.print(cout);

    return 0;
}