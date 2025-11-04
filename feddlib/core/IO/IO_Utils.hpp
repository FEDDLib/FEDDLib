#ifndef IO_UTILS_hpp
#define IO_UTILS_hpp

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_YamlParameterListCoreHelpers.hpp>

namespace FEDD{
    
    namespace IO {

        /*! 
            @brief Read a parameter list (yaml or xml format) from file.
            @param fileName File name 
            @param fileType File type ("yaml" or "xml") 
        */
        Teuchos::RCP<Teuchos::ParameterList> parseParameterList(const std::string &fileName,
                                                                const std::string &type         = "yaml")
        {
            Teuchos::RCP<Teuchos::ParameterList> parameterList;

            if (!type.compare("yaml"))
            {
                parameterList = Teuchos::getParametersFromYamlFile(fileName);
            }
            else if (!type.compare("xml"))
            {
                parameterList = Teuchos::getParametersFromXmlFile(fileName);
            }
            else
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"File type not supported.");
            }

            return parameterList;
        };


        /*! 
            @brief Read additional parameter list files, specified within the parameter list
            @param parameterList Input parameter list 
            @param nameName Parameter name in the input parameter list which contains an array with the names under which the additional parameter lists will be included 
            @param pathName Parameter name in the input parameter list which contains an array with the paths of the additional parameter lists
            @param typeName Parameter name in the input parameter list which contains an array with the types of the additional parameter lists
        */
        Teuchos::RCP<Teuchos::ParameterList> parseParameterListFiles(Teuchos::RCP<Teuchos::ParameterList>   &parameterList,
                                                                     const std::string                      &nameName       = "ParameterLists File Names",
                                                                     const std::string                      &pathName       = "ParameterLists File Paths",
                                                                     const std::string                      &typeName       = "ParameterLists File Types")
        {                    
            Teuchos::Array<std::string> names(1,""), filePaths(1,""), fileTypes(1,"");
            if (parameterList->isParameter(pathName))
            {
                filePaths = parameterList->get(pathName,filePaths);
                
                if (parameterList->isParameter(nameName))
                {
                    names = parameterList->get(nameName,names);
                    TEUCHOS_TEST_FOR_EXCEPTION(names.size()!=filePaths.size(),std::logic_error,"names.size()!=filePaths.size()");
                } else {
                    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"ParameterLists File Names are not specified.");
                }

                if (parameterList->isParameter(typeName))
                {
                    fileTypes = parameterList->get(typeName,fileTypes);
                    TEUCHOS_TEST_FOR_EXCEPTION(fileTypes.size()!=filePaths.size(),std::logic_error,"fileTypes.size()!=filePaths.size()");
                } else {
                    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"ParameterLists File Types are not specified.");
                }

                int i = 0;
                for (auto it = filePaths.begin(); it != filePaths.end(); it++, i++)
                {
                    Teuchos::RCP<Teuchos::ParameterList> parameterListTmp = parseParameterList(filePaths[i],fileTypes[i]);
                    parameterList->set(names[i],*parameterListTmp);
                }
            }
            return parameterList;
        };

    }

}

#endif
