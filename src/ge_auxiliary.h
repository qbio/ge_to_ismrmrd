/* Some helper stuff for dealing with GE scan archives 

   These classes are related to classes defined in the GE Orchestra environment, but recreated here to
   facilitate reading of scanarchives with the converters own routines, rather than relying on Orchestra 
   functionality that is either complex, not well understood, or non-existing.
*/
#pragma once

#include <string>

#include <boost/filesystem.hpp>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

namespace GE_AUX
{

//  ArchivedFileInfo is a private nested class lifted from include/recon/System/Archive/DiskArchive.h
// in Orchestra. Basic implementation is provided. There is no claim or guarantee that this class will
// perform at all identically to the class in Orchestra (especially the operator< ), but it fullfills the
// purpose of letting us use boost to deserialize the ArchivedFiles header in the ScanArchive file
class ArchivedFileInfo
{
public:
    /**                                                                                                                                                    
             * Constructor                                                                                                                                         
             */
    ArchivedFileInfo() {}

    /**                                                                                                                                                    
             * Constructor                                                                                                                                         
             * @param pathKey                                                                                                                                      
             * @param fileName                                                                                                                                     
             * @param archiveEntry                                                                                                                                 
             */
    ArchivedFileInfo(const std::string &pathKey, const boost::filesystem::path &fileName, const boost::filesystem::path &archiveEntry)
    {
        myPathKey = pathKey;
        FileName(fileName);
        ArchiveEntry(archiveEntry);
    }

    /**                                                                                                                                                    
             * Destructor                                                                                                                                          
             */
    ~ArchivedFileInfo() {}

    /**                                                                                                                                                    
             * @return const std::string&                                                                                                                          
             */
    const std::string &PathKey() const { return myPathKey; }

    /**                                                                                                                                                    
             * @return const boost::filesystem::path&                                                                                                              
             */
    const boost::filesystem::path &FileName() const { return myFileName; }

    /**                                                                                                                                                    
             * @return const boost::filesystem::path&                                                                                                              
             */
    const boost::filesystem::path &ArchiveEntry() const { return myArchiveEntry; }

    /**                                                                                                                                                    
             * Comparator                                                                                                                                          
             * @param right                                                                                                                                        
             */
    bool operator<(const ArchivedFileInfo &right) const { return myFileName < right.myFileName; }

    /**                                                                                                                                                    
             * Set FileName                                                                                                                                        
             * @param fileName                                                                                                                                     
             */
    void FileName(const boost::filesystem::path &fileName) { myFileName = fileName; }

    /**                                                                                                                                                    
             * Set ArchiveEntry                                                                                                                                    
             * @param archiveEntry                                                                                                                                 
             */
    void ArchiveEntry(const boost::filesystem::path &archiveEntry) { myArchiveEntry = archiveEntry; }

private:
    friend class boost::serialization::access;

    /**                                                                                                                                                    
             * Serialize                                                                                                                                           
             */
    template <class Archive>
    void serialize(Archive &ar, size_t /*version*/)
    {
        ar &boost::serialization::make_nvp("PathKey", myPathKey);
        ar &boost::serialization::make_nvp("FileName", myFileName);
        ar &boost::serialization::make_nvp("ArchiveEntry", myArchiveEntry);
    }

    std::string myPathKey;
    boost::filesystem::path myFileName;
    boost::filesystem::path myArchiveEntry;
};
} // namespace GE_AUX