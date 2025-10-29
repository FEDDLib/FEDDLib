#ifndef HDF5TOOLBOX_hpp
#define HDF5TOOLBOX_hpp

#include "hdf5.h"
#include "H5FDmpio.h"
#include <Tpetra_Core.hpp>
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/LinearAlgebra/Map_decl.hpp"

/*!
 Declaration of HDF5Toolbox
 
 @brief  HDF5Toolbox
 @author Lea Saßmannshausen
 @version 1.0
 @copyright LS
 */

namespace FEDD {
     /*!
    \class HDF5Toolbox
    \brief This class contains the features we use from EpteraExt::HDF5, but we transition it to Tpetra

    */
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class HDF5Toolbox{
    
    public:
        typedef Teuchos::Comm<int> Comm_Type;
        typedef Teuchos::RCP<Comm_Type> CommPtr_Type;
        typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;

        typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
        typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
        typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

        typedef Map<LO,GO,NO> Map_Type;
        typedef Teuchos::RCP<Map_Type> MapPtr_Type;
        typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;

        HDF5Toolbox(CommConstPtr_Type Comm);

        /// @brief Write/export a vector X to the HDF5 file under the group name GroupName
        /// @param GroupName 
        /// @param X 
        void write(const std::string &GroupName, const MultiVectorConstPtr_Type X, bool writeTranspose = false);

        bool isContained(std::string Name, std::string GroupName = "");

        void createGroup(const std::string &GroupName);

        void write(const std::string &GroupName, const std::string &DataSetName, double what);

        void write(const std::string &GroupName, const std::string &DataSetName, int what);

        void write(const std::string &GroupName, const std::string &DataSetName, const std::string &data);

        ~HDF5Toolbox() 
        {
        if (isOpen())
            close();
        }

        void read(const std::string &GroupName, const MapConstPtr_Type Map, MultiVectorPtr_Type X);

        void readIntVectorProperties(const std::string &GroupName, int &GlobalLength);

        void read(const std::string &GroupName, const std::string &DataSetName, GO MySize, int GlobalSize, const hid_t type, void *data);

        void read(const std::string &GroupName, const std::string &DataSetName, int &data);

        void read(const std::string &GroupName, const std::string &DataSetName, std::string &data);

        void tpetraScanSum(const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const GO *sendbuf, GO *recvbuf, int count);

        //! Create a new file.
        void create(const std::string FileName);

        //! Open specified file with given access type.
        void open(const std::string FileName, int AccessType = H5F_ACC_RDWR);

        //! Close the file.
        void close()
        {

            // ssize_t nopen = H5Fget_obj_count(file_id_, H5F_OBJ_ALL);
            
            // if (nopen > 0) {
            //     std::cout << "[Rank " << comm_->getRank()
            //         << "] HDF5 objects still open before H5Fclose: " << nopen << std::endl;
            //     H5Fget_obj_count(file_id_, H5F_OBJ_ALL);
            //     H5Fget_obj_ids(file_id_, H5F_OBJ_ALL, 0, nullptr);
            

            //     ssize_t num = H5Fget_obj_count(file_id_, H5F_OBJ_ALL);
            //     std::vector<hid_t> ids(num);
            //     H5Fget_obj_ids(file_id_, H5F_OBJ_ALL, num, ids.data());

            //     for (auto id : ids) {
            //     unsigned int type = H5Iget_type(id);
            //     const char* tname =
            //         (type == H5I_GROUP)     ? "group" :
            //         (type == H5I_DATASET)   ? "dataset" :
            //         (type == H5I_DATASPACE) ? "dataspace" :
            //         (type == H5I_DATATYPE)  ? "datatype" :
            //         (type == H5I_ATTR)      ? "attribute" : "unknown";

            //     std::cout << "[Rank " << comm_->getRank() << "] Still open: " << tname
            //                 << " (id=" << id << ")\n";
            //     }
            // }

            H5Fclose(file_id_);
            isOpen_ = false;

        }

        //! Flush the content to the file
        void flush()
        {
            H5Fflush(file_id_, H5F_SCOPE_GLOBAL);
        }

        //! Return \c true if a file has already been opened using Open()/Create()
        bool isOpen() const
        {
        return(isOpen_);
        }
            
    private:
        bool isOpen_; // Flag if file is open
        CommConstPtr_Type comm_;
        hid_t       file_id_;
        hid_t	plist_id_;
        herr_t	status;
        std::string FileName_;
};

}

#endif