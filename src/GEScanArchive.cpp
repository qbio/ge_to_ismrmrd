// Utility class supporting the reading of GE ScanArchive files
#include "GEScanArchive.h"
#include <boost/interprocess/streams/bufferstream.hpp>
#include <experimental/filesystem> 

// construct the helper class from a filename
GEScanArchive::GEScanArchive(const std::string &filename)
{
    m_pH5File = new H5::H5File(filename, H5F_ACC_RDONLY);
    // TODO: Check its actually a scan archive and that we succeeded

    // Read the ArchiveMetaData
    readArchiveMetaData();

}

GEScanArchive::~GEScanArchive()
{
    // assume its a good m_pH5File
    delete m_pH5File;
}

bool GEScanArchive::hasTrajectoryFile() const
{
    return false;
}

bool GEScanArchive::hasNoiseCovarianceFile() const
{
    return false;
}

bool GEScanArchive::hasPureCalibrationFile() const
{
    return false;
}

void GEScanArchive::readArchiveMetaData()
{
    // It is assumed that all GE ScanArchive files do have this header.
    H5::DataSet dataset = m_pH5File->openDataSet( "/ArchiveMetaData/ArchivedFiles.xml" );
    H5::StrType hdf5_strtype = dataset.getStrType();

    /* allocate memory and read the file */
    char *serfile_buffer = new char[hdf5_strtype.getSize()];
    dataset.read(serfile_buffer,H5::StrType(dataset));

    // if the set already contained entries, clear then before reading
    m_ArchivedFileInfoSet.clear();

    boost::interprocess::ibufferstream isb(serfile_buffer,hdf5_strtype.getSize());

    // do the reverse magic of what GE did to create this header
    boost::archive::xml_iarchive ia(isb);
    ia >> BOOST_SERIALIZATION_NVP(m_ArchivedFileInfoSet);

    // don't need the byte buffer anymore
    delete [] serfile_buffer;

    // small debug output routine to see what we got
    std::cout << "Boost read a set size of " << m_ArchivedFileInfoSet.size() << std::endl;
    for (auto it = m_ArchivedFileInfoSet.begin(); it != m_ArchivedFileInfoSet.end(); it++) 
        std::cout << (*it).ArchiveEntry() << " : " << (*it).PathKey() << std::endl; 
}

// static utility functions
bool GEScanArchive::archiveExists(const std::string &filename)
{
    return (std::experimental::filesystem::exists(filename) && GERecon::ScanArchive::IsArchiveFilePath(filename));
}
