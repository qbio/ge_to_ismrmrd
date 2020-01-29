// Utility class supporting the reading of GE ScanArchive files

#pragma once
#include "ge_auxiliary.h"

// Orchestra
#include <Orchestra/Common/ArchiveHeader.h>
#include <Orchestra/Common/PrepData.h>
#include <Orchestra/Common/ImageCorners.h>
#include <Orchestra/Common/ScanArchive.h>
#include <Orchestra/Common/SliceInfoTable.h>

#include <Orchestra/Control/ProcessingControl.h>

#include <Orchestra/Legacy/Pfile.h>
#include <Orchestra/Legacy/DicomSeries.h>

#include <Orchestra/Acquisition/ControlPacket.h>
#include <Orchestra/Acquisition/ControlTypes.h>
#include <Orchestra/Acquisition/Core/ArchiveStorage.h>
#include <Orchestra/Acquisition/DataTypes.h>
#include <Orchestra/Acquisition/FrameControl.h>
#include <Orchestra/Acquisition/CartesianFrameCommand.h>

class GEScanArchive
{
private:
    GEScanArchive() {} // default construction is not allowed

public:
    GEScanArchive(const std::string &filename);
    ~GEScanArchive();

    // queries for members that can be found inside GE ScanArchive files
    bool hasTrajectoryFile() const;
    bool hasNoiseCovarianceFile() const;
    bool hasPureCalibrationFile() const;

    // static utility functions
    static bool archiveExists(const std::string &filename);

protected:
    void readArchiveMetaData();

private:
    H5::H5File                        *m_pH5File;  
    std::set<GE_AUX::ArchivedFileInfo> m_ArchivedFileInfoSet;
};