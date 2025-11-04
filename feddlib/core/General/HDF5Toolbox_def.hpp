#ifndef HDF5TOOLBOX_DEF_hpp
#define HDF5TOOLBOX_DEF_hpp

/*
  

*/

namespace FEDD {

struct FindDataset_t
{
  std::string name;
  bool found;
};

static herr_t FindDataset(hid_t loc_id, const char *name, void *opdata)
{
  std::string& token = ((FindDataset_t*)opdata)->name;
  if (token == name)
    ((FindDataset_t*)opdata)->found = true;

  return(0);
}

template <class SC, class LO, class GO, class NO>
HDF5Toolbox<SC, LO, GO, NO>::HDF5Toolbox(CommConstPtr_Type comm):
  comm_(comm),
  isOpen_(false)
{}

// ----------------------------------------------------------------
/// @brief Write/export a vector X to the HDF5 file under the group name GroupName
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::write(const std::string& GroupName, const MultiVectorConstPtr_Type X, bool writeTranspose)
{

  TEUCHOS_TEST_FOR_EXCEPTION(!isOpen(),std::runtime_error,"HDF5Toolbox:: no file open yet");

  hid_t       group_id, dset_id;
  hid_t       filespace_id, memspace_id;
  herr_t      status;

  // need a linear distribution to use hyperslabs
  MultiVectorPtr_Type linearX;

  MapPtr_Type linearMap = Teuchos::rcp(new Map_Type(X->getMap()->getGlobalNumElements(),X->getMap()->getNodeNumElements(), X->getMap()->getIndexBase(), X->getMap()->getComm()));
  linearX = Teuchos::rcp(new MultiVector_Type(linearMap, X->getNumVectors()));
  linearX->importFromVector(X);


  int NumVectors = X->getNumVectors();
  int GlobalLength = X->getMap()->getGlobalNumElements(); //X.GlobalLength();

  // Whether or not we do writeTranspose or not is
  // handled by one of the components of q_dimsf, offset and count.
  // They are determined by indexT
  int indexT(0);
  if (writeTranspose) indexT = 1; 

  hsize_t q_dimsf[] = {static_cast<hsize_t>(GlobalLength), static_cast<hsize_t>(GlobalLength)};
  q_dimsf[indexT] = NumVectors;

  filespace_id = H5Screate_simple(2, q_dimsf, NULL);

  if (!isContained(GroupName))
    createGroup(GroupName);

  group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  // Create the dataset with default properties and close filespace_id.
  dset_id = H5Dcreate(group_id, "Values", H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create property list for collective dataset write.
  plist_id_ = H5Pcreate(H5P_DATASET_XFER);
#ifdef HAVE_MPI
  H5Pset_dxpl_mpio(plist_id_, H5FD_MPIO_COLLECTIVE);
#endif


  // Select hyperslab in the file.
  hsize_t offset[] = {static_cast<hsize_t>(linearX->getMap()->getGlobalElement(0)-X->getMap()->getIndexBase()),
                      static_cast<hsize_t>(linearX->getMap()->getGlobalElement(0)-X->getMap()->getIndexBase())};
  hsize_t stride[] = {1, 1};
  hsize_t count[] = {static_cast<hsize_t>(linearX->getLocalLength()),
                     static_cast<hsize_t>(linearX->getLocalLength())};
  hsize_t block[] = {1, 1};
  
  for (int n = 0; n < NumVectors; ++n)
  {
    // Select hyperslab in the file
    offset[indexT] = n;
    count [indexT] = 1;

    // Print local info before selecting hyperslab
    // std::cout << "[Rank " << comm_->getRank() << "] "
    //           << "Writing vector " << n
    //           << " | file offset = (" << offset[0] << ", " << offset[1] << ")"
    //           << " | file count = (" << count[0] << ", " << count[1] << ")"
    //           << " | local length = " << linearX->getLocalLength()
    //           << " | global length = " << 1
    //           << std::endl;

    // Get the dataset's file dataspace and select the correct hyperslab
    filespace_id = H5Dget_space(dset_id);

    herr_t selStatus = H5Sselect_hyperslab(
        filespace_id, H5S_SELECT_SET, offset, stride, count, block);

    if (selStatus < 0) {
      std::cerr << "[Rank " << comm_->getRank()
                << "] ERROR selecting hyperslab: "
                << "offset=(" << offset[0] << "," << offset[1] << "), "
                << "count=("  << count[0]  << "," << count[1]  << ")\n";
    }

    // Each process defines dataset in memory (its local buffer)
    hsize_t dimsm[] = { static_cast<hsize_t>(linearX->getLocalLength()) };
    memspace_id = H5Screate_simple(1, dimsm, NULL);

    if (memspace_id < 0) {
      std::cerr << "[Rank " << comm_->getRank()
                << "] ERROR creating memory dataspace of size "
                << dimsm[0] << std::endl;
    }

    // Print detailed info before writing
    // std::cout << "[Rank " << comm_->getRank() << "] "
    //           << "Calling H5Dwrite for vector " << n
    //           << " | memspace dims = " << dimsm[0]
    //           << " | plist_id = " << plist_id_
    //           << std::endl;

    // Perform the write
    herr_t writeStatus = H5Dwrite(
        dset_id, H5T_NATIVE_DOUBLE, memspace_id, filespace_id,
        plist_id_, linearX->getData(n).get());

    if (writeStatus < 0) {
      std::cerr << "[Rank " << comm_->getRank()
                << "] ERROR during H5Dwrite for vector " << n << std::endl;
    }


    // Optional: flush stdout to preserve rank ordering in logs
    std::cout.flush();
    std::cerr.flush();
  }

  hid_t filespace = H5Dget_space(dset_id);
  hsize_t totalElems = H5Sget_simple_extent_npoints(filespace);

  // std::cout << "[Rank " << comm_->getRank()
  //         << "] Dataset total elements = " << totalElems << std::endl;

  hssize_t fileCount = H5Sget_select_npoints(filespace_id);
  hssize_t memCount  = H5Sget_select_npoints(memspace_id);

  // std::cout << "[Rank " << comm_->getRank() << "] "
  //         << "Hyperslab file selection = " << fileCount
  //         << ", memory selection = " << memCount << std::endl;
          
  H5Gclose(group_id);
  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Pclose(plist_id_);
  // Close local memory dataspace
  H5Sclose(memspace_id);
  
  write(GroupName, "GlobalLength", GlobalLength);
  write(GroupName, "NumVectors", NumVectors);
  write(GroupName, "__type__", "Tpetra_MultiVector");
}
// ---------------------------------------------------------
// ==========================================================================
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::read(const std::string& GroupName, const MapConstPtr_Type Map,
                        MultiVectorPtr_Type X)
{
  // gets the length of the std::vector
  int GlobalLength;
  readIntVectorProperties(GroupName, GlobalLength);

  MapPtr_Type linearMap = Teuchos::rcp(new Map_Type(X->getMap()->getGlobalNumElements(),X->getMap()->getNodeNumElements(), X->getMap()->getIndexBase(), X->getMap()->getComm()));

  // we need to first create a linear map, read the std::vector,
  // then import it to the actual nonlinear map
  MultiVectorPtr_Type linearX = Teuchos::rcp(new MultiVector_Type(linearMap, X->getNumVectors()));

  read(GroupName, "Values", linearMap->getNodeNumElements(), linearMap->getGlobalNumElements(),
        H5T_NATIVE_DOUBLE, linearX->getDataNonConst(0).get());

  X->importFromVector(linearX);

}

// ==========================================================================
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::readIntVectorProperties(const std::string& GroupName,
                                           int& GlobalLength)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isContained(GroupName),std::runtime_error, "requested group " <<  GroupName << " not found");

  std::string Label;
  read(GroupName, "__type__", Label);
                    
  read(GroupName, "GlobalLength", GlobalLength);
}

template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::read(const std::string& GroupName, const std::string& DataSetName,
                        GO MySize, int GlobalSize,
                        const hid_t type, void* data)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isOpen(),std::runtime_error,"HDF5Toolbox:: no file open yet");

  // Convert MySize to hsize_t for HDF5 calls
hsize_t MySize_t = MySize;

// Compute the prefix sum across ranks to get the global offset
GO itmp = 0;
tpetraScanSum(comm_, &MySize, &itmp, 1);

// Exclusive prefix (start index for this rank)
hsize_t Offset_t = static_cast<hsize_t>(itmp - MySize);

// === DEBUG OUTPUT ===
// std::cout << "[Rank " << comm_->getRank() << "] "
//           << "Preparing to read dataset '" << DataSetName << "' in group '" << GroupName << "'\n"
//           << "  local size (MySize)     = " << MySize << "\n"
//           << "  computed offset (Offset)= " << Offset_t << "\n"
//           << "  total elements (sum)    = " << itmp << std::endl;

// Open group and dataset
hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
if (group_id < 0) {
  std::cerr << "[Rank " << comm_->getRank() << "] ERROR: H5Gopen failed for group '" 
            << GroupName << "'\n";
}

hid_t dataset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);
if (dataset_id < 0) {
  std::cerr << "[Rank " << comm_->getRank() << "] ERROR: H5Dopen failed for dataset '" 
            << DataSetName << "'\n";
}

// Get the dataset’s file dataspace
hid_t filespace_id = H5Dget_space(dataset_id);
if (filespace_id < 0) {
  std::cerr << "[Rank " << comm_->getRank() << "] ERROR: H5Dget_space failed\n";
}

// Select hyperslab in the file for this process
hsize_t offset[2] = { 0, Offset_t };   // row index, element offset
hsize_t count [2] = { 1, MySize_t };  // one vector, 29 elements

herr_t selStatus = H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

if (selStatus < 0) {
  std::cerr << "[Rank " << comm_->getRank()
            << "] ERROR selecting hyperslab: offset=" << Offset_t
            << " count=" << MySize_t << std::endl;
} else {
  std::cout << "[Rank " << comm_->getRank()
            << "] Selected hyperslab: offset=" << Offset_t
            << " count=" << MySize_t << std::endl;
}

// Create memory dataspace for local buffer
hid_t mem_dataspace = H5Screate_simple(1, &MySize_t, NULL);
if (mem_dataspace < 0) {
  std::cerr << "[Rank " << comm_->getRank()
            << "] ERROR creating memory dataspace (size=" << MySize_t << ")\n";
}

// Print before reading
// std::cout << "[Rank " << comm_->getRank()
//           << "] Reading from file hyperslab offset=" << Offset_t
//           << " count=" << MySize_t
//           << " into local buffer at " << static_cast<void*>(data)
//           << std::endl;

// Perform the read
herr_t status = H5Dread(dataset_id, type, mem_dataspace, filespace_id,
                        H5P_DEFAULT, data);

if (status < 0) {
  std::cerr << "[Rank " << comm_->getRank()
            << "] ERROR during H5Dread (offset=" << Offset_t
            << ", count=" << MySize_t << ")\n";
} else {
  std::cout << "[Rank " << comm_->getRank()
            << "] Successfully read " << MySize_t << " elements from dataset '"
            << DataSetName << "'\n";
}

hid_t filespace = H5Dget_space(dataset_id);
hsize_t totalElems = H5Sget_simple_extent_npoints(filespace);

// std::cout << "[Rank " << comm_->getRank()
//         << "] Dataset total elements = " << totalElems << std::endl;


hssize_t fileCount = H5Sget_select_npoints(filespace_id);
hssize_t memCount  = H5Sget_select_npoints(mem_dataspace);

// std::cout << "[Rank " << comm_->getRank() << "] "
//           << "Hyperslab file selection = " << fileCount
//           << ", memory selection = " << memCount << std::endl;

// Close resources
H5Sclose(mem_dataspace);
H5Sclose(filespace_id);
H5Dclose(dataset_id);
H5Gclose(group_id);

// Optional: ensure ordered output by rank
comm_->barrier();
std::cout.flush();
std::cerr.flush();
//  H5Dclose(filespace_id);
}

// -------------------------------------------------------------------------
/// @brief Create a new file.
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::create(const std::string FileName)
{
  TEUCHOS_TEST_FOR_EXCEPTION(isOpen(),std::runtime_error,"HDF5Toolbox:: an HDF5 is already open, first close the current one - using method Close(), then open/create a new one");

  FileName_ = FileName;

  // Set up file access property list with parallel I/O access
  plist_id_ = H5Pcreate(H5P_FILE_ACCESS);
// #ifdef HAVE_MPI
  {
    // Tell HDF5 what MPI communicator to use for parallel file access
    // for the above property list.
    //
    // HAVE_MPI is defined, so we know that Trilinos was built with
    // MPI.  However, we don't know whether Comm_ wraps an MPI
    // communicator.  Comm_ could very well be a serial communicator.
    MPI_Comm mpiComm = MPI_COMM_NULL; 

    Teuchos::RCP<const Teuchos::MpiComm<int>> mpiWrapper =
    Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm_, false);

    if (!mpiWrapper.is_null()) {
      mpiComm = *(mpiWrapper->getRawMpiComm());
    }
    else {
      // Try Serial communicator next
      Teuchos::RCP<const Teuchos::SerialComm<int>> serialWrapper =
        Teuchos::rcp_dynamic_cast<const Teuchos::SerialComm<int>>(comm_, false);

      if (!serialWrapper.is_null()) {
        // Comm_ is a Teuchos::SerialComm<int>
        // Use MPI_COMM_SELF for serial HDF5 access
        mpiComm = MPI_COMM_SELF;
      }
      else {
        // Unknown communicator subclass
        const char* const errMsg =
          "Tpetra::HDF5::Create: This HDF5 object was created with a "
          "Teuchos::Comm<int> that is neither Teuchos::MpiComm<int> "
          "nor Teuchos::SerialComm<int>. We don't know how to get an "
          "MPI_Comm from it.";
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,errMsg);
      }
    }

    // By this point, mpiComm should be something other than
    // MPI_COMM_NULL.  Otherwise, Comm_ wraps MPI_COMM_NULL.
    if (mpiComm == MPI_COMM_NULL) {
      const char* const errMsg = "TpetraExt::HDF5::Create: The Tpetra_Comm "
        "object with which this HDF5 instance was created wraps MPI_COMM_NULL, "
        "which is an invalid MPI communicator.  HDF5 requires a valid MPI "
        "communicator.";
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,errMsg);
    }

    // Tell HDF5 what MPI communicator to use for parallel file access
    // for the above property list.  For details, see e.g.,
    //
    // http://www.hdfgroup.org/HDF5/doc/UG/08_TheFile.html
    //
    // [last accessed 06 Oct 2011]
    H5Pset_fapl_mpio(plist_id_, mpiComm, MPI_INFO_NULL);
  }
// #endif

#if 0
  unsigned int boh = H5Z_FILTER_MAX;
  H5Pset_filter(plist_id_, H5Z_FILTER_DEFLATE, H5Z_FILTER_MAX, 0, &boh);
#endif

  // create the file collectively and release property list identifier.
  file_id_ = H5Fcreate(FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                      plist_id_);
  H5Pclose(plist_id_);

  isOpen_ = true;
}


// -------------------------------------------------------------------------
/// @brief Checking if the dataset Name is contained in the group GroupName

template <class SC, class LO, class GO, class NO>
bool HDF5Toolbox<SC, LO, GO, NO>::isContained(std::string Name, std::string GroupName)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isOpen(),std::runtime_error,"HDF5Toolbox:: no file open yet");

  FindDataset_t data;
  data.name = Name;
  data.found = false;

  // recursively look for groups
  size_t pos = Name.find("/");
  if (pos != std::string::npos)
  {
    std::string NewGroupName = Name.substr(0, pos);
    if (GroupName != "")
      NewGroupName = GroupName + "/" + NewGroupName;
    std::string NewName = Name.substr(pos + 1);
    return isContained(NewName, NewGroupName);
  }

  GroupName = "/" + GroupName;

  //int idx_f =
  H5Giterate(file_id_, GroupName.c_str(), NULL, FindDataset, (void*)&data);

  return(data.found);
}
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
//! Create group \c GroupName.
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::createGroup(const std::string& GroupName)
{
  hid_t group_id = H5Gcreate(file_id_, GroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);
}
// -------------------------------------------------------------------------



// ==========================================================================
// ==========================================================================

template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::write(const std::string& GroupName, const std::string& DataSetName,
                            double what)
{
  if (!isContained(GroupName))
    createGroup(GroupName);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), H5T_NATIVE_DOUBLE,
                            filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
                           filespace_id, H5P_DEFAULT, &what);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================

template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::write(const std::string& GroupName, const std::string& DataSetName,
                            int what)
{
  if (!isContained(GroupName))
    createGroup(GroupName);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), H5T_NATIVE_INT,
                      filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace_id,
                           H5P_DEFAULT, &what);

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}


// ==========================================================================
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::write(const std::string& GroupName,
                            const std::string& DataSetName,
                            const std::string& data)
{
  if (!isContained(GroupName))
    createGroup(GroupName);

  hsize_t len = 1;

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dataspace_id = H5Screate_simple(1, &len, NULL);

  hid_t atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, data.size() + 1);

  hid_t dataset_id = H5Dcreate(group_id, DataSetName.c_str(), atype,
                               dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Dwrite(dataset_id, atype, H5S_ALL, H5S_ALL,H5P_DEFAULT, data.c_str());
  
  // Close/release resources.
  H5Dclose(dataset_id);
  H5Gclose(group_id);
}
// ==========================================================================
// ==========================================================================

// ==========================================================================
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::read(const std::string& GroupName, const std::string& DataSetName, int& data)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isContained(GroupName),std::runtime_error, "requested group " + GroupName + " not found");

  // Create group in the root group using absolute name.
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  herr_t status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace_id,
                    H5P_DEFAULT, &data);

  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Gclose(group_id);
}

template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::read(const std::string& GroupName,
                           const std::string& DataSetName,
                           std::string& data)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isContained(GroupName),std::runtime_error, "requested group " + GroupName + " not found");

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dataset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  hid_t datatype_id = H5Dget_type(dataset_id);
//  size_t typesize_id = H5Tget_size(datatype_id);
  H5T_class_t typeclass_id = H5Tget_class(datatype_id);

  if(typeclass_id != H5T_STRING)
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error, "requested group " + GroupName + " is not a std::string");

  char data2[160];
  H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL,H5P_DEFAULT, data2) ;                

  data = data2;

  H5Dclose(dataset_id);
  H5Gclose(group_id);
}
// ------------------------------------------
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::open(const std::string FileName, int AccessType)
{
  TEUCHOS_TEST_FOR_EXCEPTION(isOpen(),std::runtime_error,"HDF5Toolbox:: no file open yet");

  FileName_ = FileName;

  // Set up file access property list with parallel I/O access
  plist_id_ = H5Pcreate(H5P_FILE_ACCESS);

#ifdef HAVE_MPI
  MPI_Comm mpiComm = MPI_COMM_WORLD;

  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiWrapper =
    Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm_, false);

  if (!mpiWrapper.is_null()) {
    mpiComm = *(mpiWrapper->getRawMpiComm());
  }

  H5Pset_fapl_mpio(plist_id_, mpiComm, MPI_INFO_NULL);
#endif

  // create the file collectively and release property list identifier.
  file_id_ = H5Fopen(FileName.c_str(), AccessType, plist_id_);
  H5Pclose(plist_id_);

  isOpen_ = true;
}


// ----------------------------------------------
template <class SC, class LO, class GO, class NO>
void HDF5Toolbox<SC, LO, GO, NO>::tpetraScanSum(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                   const GO* sendbuf,
                   GO* recvbuf,
                   int count)
{
#ifdef HAVE_MPI
  // Try to downcast to MPI communicator
  Teuchos::RCP<const Teuchos::MpiComm<int>> mpiComm =
    Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(comm, false);

  if (!mpiComm.is_null()) {
    // Use MPI_Scan directly
    MPI_Comm rawComm = *(mpiComm->getRawMpiComm());
    MPI_Scan(sendbuf, recvbuf, count, MPI_INT, MPI_SUM, rawComm);
    return;
  }
#endif

  // Serial (or non-MPI) fallback
  for (int i = 0; i < count; ++i) {
    recvbuf[i] = sendbuf[i];
  }
}
}
#endif
