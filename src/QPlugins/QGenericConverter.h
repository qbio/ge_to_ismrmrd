
/** @file QGenericConverter.h */
#ifndef QGENERIC_CONVERTER_H
#define QGENERIC_CONVERTER_H

#include "GenericConverter.h"

// TW GE stuff
#include <imagedb.h>

//namespace PfileToIsmrmrd {

class QGenericConverter: public PfileToIsmrmrd::GenericConverter
{
public:

    virtual std::vector<ISMRMRD::Acquisition> getAcquisitions (GERecon::Legacy::PfilePointer &pfile,
                                                               unsigned int view_num);

    virtual std::vector<ISMRMRD::Acquisition> getAcquisitions (GERecon::ScanArchivePointer &scanArchivePtr,
                                                               unsigned int view_num);

protected:

    virtual std::vector<ISMRMRD::Acquisition> getAcquisitions3D(GERecon::ScanArchivePointer &scanArchivePtr,
                                                               unsigned int view_num);

    virtual std::vector<ISMRMRD::Acquisition> getAcquisitions2D(GERecon::ScanArchivePointer &scanArchivePtr,
                                                               unsigned int view_num);

    virtual int get_view_idx(GERecon::Control::ProcessingControlPointer processingControl,
                             unsigned int view_num, ISMRMRD::EncodingCounters &idx);

    bool isSliceIndexInferred(std::string& psdname) const;
            
};

//} // namespace PfileToIsmrmrd

#endif /* QGENERIC_CONVERTER_H */

