#ifndef IO_ParameterList_decl_hpp
#define IO_ParameterList_decl_hpp

#include "IO_Utils.hpp"

/*!
 Declaration of ParameterListParser
 
 @brief  BCBuilder
 @author Alexander Heinlein
 @version 1.0
 @copyright AH
 */

namespace FEDD {

    namespace IO {

        /*!
        \class ParameterList
        \brief A wrapper of Teuchos::ParameterList from Trilinos that automatically includes additional parameter lists stored in files 

        This class is a leightweight interface to Teuchos::ParameterList from Trilinos. The main property is that parameters "ParameterLists File Names", "ParameterLists File Paths", and "ParameterLists File Types" may be specified in the parameter list. At construction time, those files are read and the corresponding parameter lists are included into the list itself. 
        */

        class ParameterList {
            
        protected:

            using ParameterList_Type    = Teuchos::ParameterList;
            using ParameterListPtr_Type = Teuchos::RCP<ParameterList_Type>;
            
        public:

            /*! 
            @brief Read FEDD::IO::ParameterList from yaml or xml input file
            @param name Name of the parameter list
            @param fileName File name 
            @param fileType File type ("yaml" or "xml") 
            */
            ParameterList(const std::string &name,
                          const std::string &fileName,
                          const std::string &fileType = "yaml") :
                name_ (name),
                parameterList_ ()
            {
                parameterList_ = parseParameterList(fileName,fileType);
                parameterList_ = parseParameterListFiles(parameterList_,"ParameterLists File Names","ParameterLists File Paths","ParameterLists File Types");
            };

            /*! 
            @brief Construct FEDD::IO::ParameterList from Teuchos::ParameterList
            @param name Name of the parameter list
            @param parameterList Teuchos::RCP of the underlying Teuchos::ParameterList
            */
            ParameterList(const std::string &name,
                          ParameterListPtr_Type &parameterList) :
                name_ (name),
                parameterList_ (parameterList)
            {
                parameterList_ = parseParameterListFiles(parameterList_,"ParameterLists File Names","ParameterLists File Paths","ParameterLists File Types");
            };

            /*! 
            @brief Print parameter list
            */
            std::ostream& print(std::ostream &os, int indent=0, bool showTypes=false, bool showFlags=true, bool showDefault=true)
            {
                os << name_ << std::endl;
                for (size_t i = 0; i < name_.length(); i++)
                {
                    os << "=";
                }
                os << std::endl;
                return parameterList_->print(os,indent,showTypes,showFlags,showDefault);
            };

            /*! 
            @brief Access to the underlying Teuchos::ParamterList
            */
            ParameterListPtr_Type operator() ()
            {
                return parameterList_;
            };

        protected:
            
            const std::string name_; ///< Name of the parameter list
            
            ParameterListPtr_Type parameterList_; ///< Teuchos::RCP of the underlying Teuchos::ParameterList
        
        };

    }

}
#endif
