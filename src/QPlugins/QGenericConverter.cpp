
/** @file QGenericConverter.cpp */
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <numeric>
#include <vector>
#include <array>
#include <iostream>
#include <functional>
#include <iterator>

#include <Orchestra/Arc/Arc.h>
#include <Orchestra/Arc/SamplingPattern.h>

#include "QGenericConverter.h"

struct LOADTEST {
   LOADTEST() { std::cerr << __FILE__ << ": shared object loaded"   << std::endl; }
  ~LOADTEST() { std::cerr << __FILE__ << ": shared object unloaded" << std::endl; }
} loadTest;

//namespace PfileToIsmrmrd {

// use the psdname to determine if the converter should infer the slice indices rather than
// use the index stored in the header
// returns:
// true if the psdname matches a sequence where it is known that the packets don't store the slice index
// false if the psdname is not one of the above
bool QGenericConverter::isSliceIndexInferred(std::string& psdname) const
{
   // throw exception if the psdname is empty

   // match strings
   // ssfse/HASTE
   if(psdname.compare("ssfse") == 0)
      return true;

   else
      return false;
   
}


int QGenericConverter::get_view_idx(GERecon::Control::ProcessingControlPointer processingControl,
                                   unsigned int view_num, ISMRMRD::EncodingCounters &idx)
{
    // set all the ones we don't care about to zero
    idx.kspace_encode_step_2 = 0;
    idx.average = 0;
    idx.contrast = 0;
    idx.phase = 0;
    idx.set = 0;
    idx.segment = 0;
    for (int n=0; n<8; n++) {
        idx.user[n] = 0;
    }

    unsigned int nframes   = processingControl->Value<int>("AcquiredYRes");
    unsigned int numSlices = processingControl->Value<int>("NumSlices");

    idx.repetition = view_num / (numSlices * (1 + nframes));

    if (view_num < 1) {
        // this is the mean baseline view return -1
        return -1;
    }

    return 1;
}



std::vector<ISMRMRD::Acquisition> QGenericConverter::getAcquisitions(GERecon::Legacy::PfilePointer &pfile,
                                                                    unsigned int acqMode)
{
    std::vector<ISMRMRD::Acquisition> acqs;

    const GERecon::Control::ProcessingControlPointer processingControl(pfile->CreateOrchestraProcessingControl());
    unsigned int nPhases   = processingControl->Value<int>("AcquiredYRes");
    unsigned int nEchoes   = processingControl->Value<int>("NumEchoes");
    unsigned int nChannels = processingControl->Value<int>("NumChannels");
    unsigned int numSlices = processingControl->Value<int>("NumSlices");

    // Make number of acquisitions to be converted
    acqs.resize(numSlices * nEchoes * nPhases);

    unsigned int acq_num = 0;

    // Orchestra API provides size in bytes.
    // frame_size is the number of complex points in a single channel
    size_t frame_size = processingControl->Value<int>("AcquiredXRes");

    for (int sliceCount = 0 ; sliceCount < numSlices ; sliceCount++)
    {
        for (int echoCount = 0 ; echoCount < nEchoes ; echoCount++)
        {
            for (int phaseCount = 0 ; phaseCount < nPhases ; phaseCount++)
            {
                // Grab a reference to the acquisition
                ISMRMRD::Acquisition& acq = acqs.at(acq_num);

                // Set size of this data frame to receive raw data
                acq.resize(frame_size, nChannels, 0);
                acq.clearAllFlags();

                // Initialize the encoding counters for this acquisition.
                ISMRMRD::EncodingCounters idx;
                get_view_idx(processingControl, 0, idx);

                idx.slice = sliceCount;
                idx.contrast  = echoCount;
                idx.kspace_encode_step_1 = phaseCount;

                acq.idx() = idx;

                // Fill in the rest of the header
                // acq.measurement_uid() = pfile->RunNumber();
                acq.scan_counter() = acq_num;
                acq.acquisition_time_stamp() = time(NULL); // TODO: can we get a timestamp?
                for (int p=0; p<ISMRMRD::ISMRMRD_PHYS_STAMPS; p++) {
                    acq.physiology_time_stamp()[p] = 0;
                }
                acq.available_channels() = nChannels;
                acq.discard_pre() = 0;
                acq.discard_post() = 0;;
                acq.center_sample() = frame_size/2;
                acq.encoding_space_ref() = 0;
                //acq.sample_time_us() = pfile->sample_time * 1e6;

                for (int ch = 0 ; ch < nChannels ; ch++) {
                    acq.setChannelActive(ch);
                }

                // Patient table off-center
                // TODO: fix the patient table position
                acq.patient_table_position()[0] = 0.0;
                acq.patient_table_position()[1] = 0.0;
                acq.patient_table_position()[2] = 0.0;

                // Slice position and orientation
                /* TODO
                static pfile_slice_vectors_t slice_vectors;
                pfile_get_slice_vectors(pfile, idx.slice, &slice_vectors);

                acq.read_dir()[0] = slice_vectors.read_dir.x;
                acq.read_dir()[1] = slice_vectors.read_dir.y;
                acq.read_dir()[2] = slice_vectors.read_dir.z;
                acq.phase_dir()[0] = slice_vectors.phase_dir.x;
                acq.phase_dir()[1] = slice_vectors.phase_dir.y;
                acq.phase_dir()[2] = slice_vectors.phase_dir.z;
                acq.slice_dir()[0] = slice_vectors.slice_dir.x;
                acq.slice_dir()[1] = slice_vectors.slice_dir.y;
                acq.slice_dir()[2] = slice_vectors.slice_dir.z;
                acq.position()[0] = slice_vectors.center.x;
                acq.position()[1] = slice_vectors.center.y;
                acq.position()[2] = slice_vectors.center.z;
                */

                // Set first acquisition flag
                if (idx.kspace_encode_step_1 == 0)
                    acq.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE);

                // Set last acquisition flag
                if (idx.kspace_encode_step_1 == nPhases - 1)
                    acq.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

                // Get data from P-file using KSpaceData object, and copy
                // into ISMRMRD space.
                for (int channelID = 0 ; channelID < nChannels ; channelID++)
                {
                    // VR + JAD - 2016.01.15 - looking at various schemes to stride and read in
                    // K-space data.
                    //
                    // ViewData - will read in "acquisitions", including baselines, starting at
                    //            index 0, going up to slices * echo * (view + baselines)
                    //
                    // KSpaceData (slice, echo, channel, phase = 0) - reads in data, assuming "GE
                    //            native" data order in P-file, gives one slice / image worth of
                    //            K-space data, with baseline views automagically excluded.
                    //
                    // KSpaceData can return different numerical data types.  Picked float to
                    // be consistent with ISMRMRD data type.  This implementation of KSpaceData
                    // is used for data acquired in the "native" GE order.

                    auto kData = pfile->KSpaceData<float>(sliceCount, echoCount, channelID);

                    if (processingControl->Value<bool>("ChopY") == 0) {
                       if (idx.kspace_encode_step_1 % 2 == 1) {
                          kData *= -1.0f;
                       }
                    }

                    for (int i = 0 ; i < frame_size ; i++)
                    {
                       acq.data(i, channelID) = kData(i, phaseCount);
                    }
                }

                acq_num++;
            } // end of phaseCount loop
        } // end of echoCount loop
    } // end of sliceCount loop

    return acqs;
}

/**
*
* @throws std::runtime_error { if plugin encounters data it cannot handle }
*/
std::vector<ISMRMRD::Acquisition> QGenericConverter::getAcquisitions(GERecon::ScanArchivePointer &scanArchivePtr,
                                                                    unsigned int acqMode)
{
   GERecon::Legacy::LxDownloadDataPointer lxData = boost::dynamic_pointer_cast<GERecon::Legacy::LxDownloadData>(scanArchivePtr->LoadDownloadData());

   // Make a decision which function to run based on the 3D flag
   if(lxData->Is3D()) {
      return getAcquisitions3D(scanArchivePtr,acqMode);
   } else {
      std::cout << "XXO 2D stuff XXO" << std::endl;
      return getAcquisitions2D(scanArchivePtr,acqMode);
   }
}

/**
*
* @throws std::runtime_error { if plugin encounters data it cannot handle }
*/

/* TW Note:

   below is going to be some of the saddest code I have ever written. The reason is that GE is stupid in how
   they store their data, so I have to do rather stupid things to read it.

   As a result, there is going to be a ton of stuff this code cannot do, and in the beginning there won't be
   sufficient checks indicating that this code cannot convert a given dataset correctly.
   
   When there is a battery of test data, these checks can be incorporated.

   For sure right now it cannot do:
   - multislab data
   

*/
std::vector<ISMRMRD::Acquisition> QGenericConverter::getAcquisitions3D(GERecon::ScanArchivePointer &scanArchivePtr,
                                                                    unsigned int acqMode)
{
   std::vector<ISMRMRD::Acquisition> acqs;

   GERecon::Acquisition::ArchiveStoragePointer archiveStoragePointer = GERecon::Acquisition::ArchiveStorage::Create(scanArchivePtr);
   GERecon::Legacy::LxDownloadDataPointer lxData = boost::dynamic_pointer_cast<GERecon::Legacy::LxDownloadData>(scanArchivePtr->LoadDownloadData());
   boost::shared_ptr<GERecon::Legacy::LxControlSource> const controlSource = boost::make_shared<GERecon::Legacy::LxControlSource>(lxData);
   GERecon::Control::ProcessingControlPointer processingControl = controlSource->CreateOrchestraProcessingControl();

   int const   packetQuantity = archiveStoragePointer->AvailableControlCount();

   int            packetCount = 0;
   int              dataIndex = 0;
   int                acqType = 0;
   unsigned int       nEchoes = processingControl->Value<int>("NumEchoes");
   unsigned int     nChannels = processingControl->Value<int>("NumChannels");
   unsigned int     numSlices = processingControl->Value<int>("NumSlices");
   const int       frame_size = processingControl->Value<int>("AcquiredXRes");
   const int           nViews = processingControl->Value<int>("AcquiredYRes");
   const int          acqZRes = processingControl->Value<int>("AcquiredZRes");

   bool      lastMeasFlagSet  = false;
   size_t     expectedFrames  = 0;
   size_t     expectedPartitionsPerView = 0;
   size_t     expectedViewsPerPartition = 0;

   size_t     arcCalibrationKYStart = 0; // view index of first calibration ky location
   size_t     arcCalibrationKYEnd = 0;   // view index of last calibration ky location
   size_t     arcCalibrationKZStart = 0; // view index of first calibration kz location
   size_t     arcCalibrationKZEnd = 0;   // view index of last calibration kz location

   // too lazy to initialize these vectors to actual size, make them huge instead
   std::vector<size_t> partitionsPerView(4096,0);
   std::vector<size_t> viewsPerPartition(4096,0);

   // create a sampling pattern for Arc
   boost::scoped_ptr<GERecon::Arc::SamplingPattern> samplingPattern;

   // TODO: Remove these screen prints
   std::cout << "3D Acquisition to deal with !" << std::endl;
   std::cout << "PSDname: " << const_cast<const GERecon::Legacy::LxDownloadData&>(*lxData).ImageHeaderData().psdname << std::endl;
   std::cout << "[AXR,AXY,AXZ] = [" << frame_size << "," << nViews << "," << acqZRes << "]" << std::endl;
  
   // Initialize the trigger counts
   expectedPartitionsPerView = acqZRes;
   if(lxData->IsArc()) {
      const int patternID = processingControl->Value<int>("ArcSamplingPatternID");
      const int kYPeak = processingControl->Value<int>("ArcKYPeak");
      
      std::string spfilename = "/tmp/kacq_yz.txt." + std::to_string(patternID);

      std::cout << "Loading sampling pattern: " << spfilename << std::endl;
      // generate a filename from the patternID
      // this code assumes that the sampling pattern file exists, because it would have been generated by a script
      // before the converter is run
      samplingPattern.reset(new Arc::SamplingPattern(spfilename.c_str(),kYPeak));
      if (!samplingPattern->IsValid()) {
         // Need to throw an exception here, because its game over at this point
         std::cout << "No valid sampking pattern loaded for ID: " << patternID << std::endl;
      } else {
         std::cout << "ARC Acceleration detected: [" << samplingPattern->YAcceleration() << " x " << samplingPattern->ZAcceleration() << "]" << std::endl;
         std::cout << "ARC Sampling Pattern Res: [" << samplingPattern->YRes() << " x " << samplingPattern->ZRes() << "]" << std::endl;

         // Look at the EncodesVector
         std::cout << "Size of encodes vector: [" << samplingPattern->AcquiredLocations().rows() << " x " << samplingPattern->AcquiredLocations().cols() << "]" << std::endl;

         // The following analysis makes only sense for cartesian regular sampling
         std::vector<size_t> ky_locations;
         std::vector<size_t> kz_locations;

         // Don't know how to use an iterator for this type, so use for loop for now
         for(size_t k=0; k<samplingPattern->AcquiredLocations().rows(); k++) {
            ky_locations.push_back(samplingPattern->AcquiredLocations()(k).ky);
            kz_locations.push_back(samplingPattern->AcquiredLocations()(k).kz);
         }

         // first sort the indices
         std::sort(ky_locations.begin(),ky_locations.end());
         std::sort(kz_locations.begin(),kz_locations.end());   

         // then make them unique
         std::vector<size_t>::iterator k_it;
         k_it = std::unique(ky_locations.begin(),ky_locations.end());
         ky_locations.resize(std::distance(ky_locations.begin(),k_it) ); 
         k_it = std::unique(kz_locations.begin(),kz_locations.end());  
         kz_locations.resize(std::distance(kz_locations.begin(),k_it) ); 

         std::cout << "Unique ky locations in the sampling pattern " << ky_locations.size() << std::endl;
         std::cout << "Unique kz locations in the sampling pattern " << kz_locations.size() << std::endl;

         // Find the ARC calibration data

         // init the calibration data assuming no acceleration in either axis
         // will be overwritten by the cases where ky or kz is accelerated
         arcCalibrationKYStart = ky_locations.front();
         arcCalibrationKYEnd = ky_locations.back();
         arcCalibrationKZStart = kz_locations.front();
         arcCalibrationKZEnd = kz_locations.back();

         if(samplingPattern->YAcceleration() > 1) {
            size_t kyCalStartIndex = 0;
            size_t kyCalEndIndex = 0;

            // keep this for debugging
            /*
            for (auto n : ky_locations)
               std::cout << n << ' ';
            std::cout << std::endl;
            */

            std::vector<size_t> d(ky_locations.size(),0);
            std::adjacent_difference(ky_locations.begin(),ky_locations.end(),d.begin());
            // Note: as result of std::adjacent_difference d[0] = 0 ALWAYS!

            // Find the start of the autocalibration range
            for (size_t n = 0; n<d.size()-1; n++) {
               if(d[n+1] == 1) {
                  kyCalStartIndex = n;
                  break;
               } else if (n>0 && d[n+1] != samplingPattern->YAcceleration()) {
                  // TODO: throw exception here
                  std::cout << "ERROR: Sampling pattern is not regular ARC sampling" << std::endl;
                  break;
               } else if (n == d.size()-2) {
                  // TODO: throw exception here
                  std::cout << "ERROR: Now ARC autocalibration data found !" << std::endl;
                  break;
               }
            }

            // find the end of the autocalibration range
            // Note that the first gap after the autocalibration data maybe LESS than Ry in GE data,
            // super wacky
            for (size_t n = kyCalStartIndex; n<d.size()-1; n++) {
               if(d[n+1] > 1) {
                  kyCalEndIndex = n;
                  break;
               } else if (n>0 && d[n+1] != 1) {
                  // TODO: throw exception here
                  std::cout << "ERROR: Sampling pattern is not regular ARC sampling" << std::endl;
                  break;
               } else if (n == d.size()-2) {
                  // if there is no accelerated data in the end, then calibration goes to the end
                  kyCalEndIndex = d.size()-1;
                  break;
               }              
            }

            // check that the sampling past the autocalibration lines is alright
            // note, check starts with +2, because of the "uneven" gap thats possible after the autocalibration data
            for (size_t n = kyCalEndIndex+2; n<d.size(); n++) {
               if(d[n] != samplingPattern->YAcceleration()) {
                  // TODO: throw exception here
                  std::cout << "ERROR: Sampling pattern is not regular ARC sampling" << std::endl;
                  break;
               }
            }
            arcCalibrationKYStart = ky_locations[kyCalStartIndex];
            arcCalibrationKYEnd   = ky_locations[kyCalEndIndex];

            std::cout << "First ACS ky line = " << arcCalibrationKYStart << std::endl;
            std::cout << "Last ACS ky line = " << arcCalibrationKYEnd << std::endl;
         }

         if(samplingPattern->ZAcceleration() > 1) {
            size_t kzCalStartIndex = 0;
            size_t kzCalEndIndex = 0;

            // keep this for debugging
            /*
            for (auto n : kz_locations)
               std::cout << n << ' ';
            std::cout << std::endl;
            */

            std::vector<size_t> d(kz_locations.size(),0);
            std::adjacent_difference(kz_locations.begin(),kz_locations.end(),d.begin());
            // Note: as result of std::adjacent_difference d[0] = 0 ALWAYS!

            // Find the start of the autocalibration range
            for (size_t n = 0; n<d.size()-1; n++) {
               if(d[n+1] == 1) {
                  kzCalStartIndex = n;
                  break;
               } else if (n>0 && d[n+1] != samplingPattern->ZAcceleration()) {
                  // TODO: throw exception here
                  std::cout << "ERROR: Sampling pattern is not regular ARC sampling" << std::endl;
                  break;
               } else if (n == d.size()-2) {
                  // TODO: throw exception here
                  std::cout << "ERROR: Now ARC autocalibration data found !" << std::endl;
                  break;
               }
            }

            // find the end of the autocalibration range
            // Note that the first gap after the autocalibration data maybe LESS than Ry in GE data,
            // super wacky
            for (size_t n = kzCalStartIndex; n<d.size()-1; n++) {
               if(d[n+1] > 1) {
                  kzCalEndIndex = n;
                  break;
               } else if (n>0 && d[n+1] != 1) {
                  // TODO: throw exception here
                  std::cout << "ERROR: Sampling pattern is not regular ARC sampling" << std::endl;
                  break;
               } else if (n == d.size()-2) {
                  // if there is no accelerated data in the end, then calibration goes to the end
                  kzCalEndIndex = d.size()-1;
                  break;
               }              
            }

            // check that the sampling past the autocalibration lines is alright
            // note, check starts with +2, because of the "uneven" gap thats possible after the autocalibration data
            for (size_t n = kzCalEndIndex+2; n<d.size(); n++) {
               if(d[n] != samplingPattern->ZAcceleration()) {
                  // TODO: throw exception here
                  std::cout << "ERROR: Sampling pattern is not regular ARC sampling" << std::endl;
                  break;
               }
            }
            arcCalibrationKZStart = kz_locations[kzCalStartIndex];
            arcCalibrationKZEnd   = kz_locations[kzCalEndIndex];

            std::cout << "First ACS kz line = " << arcCalibrationKZStart << std::endl;
            std::cout << "Last ACS kz line = " << arcCalibrationKZEnd << std::endl;
         }

         // set some of the kspace parameters
         expectedViewsPerPartition = ky_locations.size();
         expectedPartitionsPerView = kz_locations.size();

         expectedFrames = samplingPattern->AcquiredLocations().rows();
      }
   } else {
      expectedFrames = nViews*acqZRes;
      expectedViewsPerPartition = nViews; // ignoring Asset
   }

   // TODO:
   // 
   // Create a loop structure for the partitions and views
   // Start setting flags for FFT conditions
   // Add code to detect partial fourier in 3D direction (check with VIBE scans)
   // Identify the loop order, based on view counter
   //
   // Add code to deal with ARC sampling pattern
   //
   while (packetCount < packetQuantity)
   {
      // encoding IDs to fill ISMRMRD headers.
      int    pe2viewID = 0;
      int    viewID    = 0;

      GERecon::Acquisition::FrameControlPointer const thisPacket = archiveStoragePointer->NextFrameControl();

      // Need to identify opcode(s) here that will mark acquisition / reference / control
      if (thisPacket->Control().Opcode() != GERecon::Acquisition::ScanControlOpcode)
      {

         // this is an issue here, because this limits this to ProgrammableControlPackets
         GERecon::Acquisition::ProgrammableControlPacket const packetContents = thisPacket->Control().Packet().As<GERecon::Acquisition::ProgrammableControlPacket>();

         viewID    = GERecon::Acquisition::GetPacketValue(packetContents.viewNumH,  packetContents.viewNumL);
         pe2viewID = GERecon::Acquisition::GetPacketValue(packetContents.sliceNumH, packetContents.sliceNumL);

         // output something if the opcode is not 1, to make it less noisy
         if(static_cast<uint8_t>(packetContents.opcode) != GERecon::Acquisition::ProgrammableOpcode)
            std::cout << "viewID = " << viewID << " pe2viewID = " << pe2viewID << " opcode = " << (int)(packetContents.opcode) << std::endl << std::flush;
         
         // TW: temporary hack, ignore the ViewCopy opcode
         //     I'm undecided on strategy here, but Packets with no data cause:
         //     - undefined behaviour/segfault on the ->Data() call (Boo GE)
         //     - we have nothing to write into the ISMRMRD file either, when we have no data 
         if(static_cast<uint8_t>(packetContents.opcode) != GERecon::Acquisition::ViewCopyOpcode) {
            if ((viewID < 1) || (viewID > nViews))
            {
               acqType = GERecon::Acquisition::BaselineFrame;
               // nothing else to be done here for basic 2D case
            }    
            else
            {
               acqType = GERecon::Acquisition::ImageFrame;

               acqs.resize(dataIndex + 1);

               // TW build in some segmentation fault protection here, as Orchestra WILL segfault on the ->Data() call
               //    I'm interpreting a 0 rval as meaning that there is no kspace data in the packet as the 
               //    header states that the data size is [FrameSize x NumChannels x NumFrames] 
               if(thisPacket->NumFrames() == 0) {
                  throw std::runtime_error("Attempting to access data of a packet that has no data.");                  
               }
 
               // TW: Now we are seeing some C++11 suddenly some type inference
               auto kData = thisPacket->Data();

               // Grab a reference to the acquisition
               ISMRMRD::Acquisition& acq = acqs.at(dataIndex);

               // Set size of this data frame to receive raw data
               acq.resize(frame_size, nChannels, 0);
               acq.clearAllFlags();

               // Initialize the encoding counters for this acquisition.
               ISMRMRD::EncodingCounters idx;
               get_view_idx(processingControl, 0, idx);

               idx.slice                  = 0; // assuming single slab 3D, this should be always zero
               idx.contrast               = packetContents.echoNum;
               idx.kspace_encode_step_1   = viewID - 1; // TW: Why?
               idx.kspace_encode_step_2   = pe2viewID;

               // increment the counters
               ++partitionsPerView[idx.kspace_encode_step_1];
               ++viewsPerPartition[idx.kspace_encode_step_2];

               acq.idx() = idx;

               // Fill in the rest of the header
               // acq.measurement_uid() = pfile->RunNumber();
               acq.scan_counter() = dataIndex;
               acq.acquisition_time_stamp() = time(NULL);
               for (int p=0; p<ISMRMRD::ISMRMRD_PHYS_STAMPS; p++) {
                  acq.physiology_time_stamp()[p] = 0;
               }
               acq.available_channels()   = nChannels;
               acq.discard_pre()          = 0;
               acq.discard_post()         = 0;
               acq.center_sample()        = frame_size/2;
               acq.encoding_space_ref()   = 0;
               // acq.sample_time_us()       = pfile->sample_time * 1e6;

               for (int ch = 0 ; ch < nChannels ; ch++) {
                  acq.setChannelActive(ch);
               }

               /*
               Flag setting functions, this is a work in progress and will need to be refactored

               */

               // The last in measurement flag is one of the most important ones, as it is an extremely valuable marker
               // that can show if a measurement is complete
               if(dataIndex == expectedFrames-1) {
                  acq.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_MEASUREMENT);
                  lastMeasFlagSet = true;
               }

               if(dataIndex > expectedFrames-1) {
                 // this needs to be replaced with an exception
                 std::cerr << "FATAL ERROR: More data than expected received." << std::endl;
               }

               // set the FFT trigger flags
               if(viewsPerPartition[idx.kspace_encode_step_2] == 1)
                  acq.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1);

               if(viewsPerPartition[idx.kspace_encode_step_2] == expectedViewsPerPartition)
                  acq.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1);

               if(partitionsPerView[idx.kspace_encode_step_1] == 1)
                  acq.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP2);

               if(partitionsPerView[idx.kspace_encode_step_1] == expectedPartitionsPerView)
                  acq.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_ENCODE_STEP2);

               // mark the autocalibration lines, TODO: Add kz box
               if(lxData->IsArc()) {
                  if(idx.kspace_encode_step_1 >= arcCalibrationKYStart 
                     && idx.kspace_encode_step_1 <= arcCalibrationKYEnd
                     && idx.kspace_encode_step_2 >= arcCalibrationKZStart
                     && idx.kspace_encode_step_2 <= arcCalibrationKZEnd)
                     acq.setFlag(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING);
               }

               for (int channelID = 0 ; channelID < nChannels ; channelID++)
               {
                  for (int i = 0 ; i < frame_size ; i++)
                  {
                     // The last dimension here in kData denotes the view
                     // index in the control packet that one must stride
                     // through to get data.  TODO - figure out if this
                     // can be programatically determined, and if so, use
                     // it. Will be needed for cases where multiple lines
                     // of data are contained in a single packet.
		               acq.data(i, channelID) = kData(i, channelID, 0);
                  }
               }

               dataIndex++;
            }
         } // if(static_cast<uint8_t>(packetContents.opcode) == 1)
      }

      packetCount++;
   }

   return acqs;
}

/**
*
* @throws std::runtime_error { if plugin encounters data it cannot handle }
*/
std::vector<ISMRMRD::Acquisition> QGenericConverter::getAcquisitions2D(GERecon::ScanArchivePointer &scanArchivePtr,
                                                                    unsigned int acqMode)
{
   std::vector<ISMRMRD::Acquisition> acqs;

   GERecon::Acquisition::ArchiveStoragePointer archiveStoragePointer = GERecon::Acquisition::ArchiveStorage::Create(scanArchivePtr);
   GERecon::Legacy::LxDownloadDataPointer lxData = boost::dynamic_pointer_cast<GERecon::Legacy::LxDownloadData>(scanArchivePtr->LoadDownloadData());
   boost::shared_ptr<GERecon::Legacy::LxControlSource> const controlSource = boost::make_shared<GERecon::Legacy::LxControlSource>(lxData);
   GERecon::Control::ProcessingControlPointer processingControl = controlSource->CreateOrchestraProcessingControl();

   int const   packetQuantity = archiveStoragePointer->AvailableControlCount();

   int            packetCount = 0;
   int              dataIndex = 0;
   int                acqType = 0;
   unsigned int       nPhases = processingControl->Value<int>("AcquiredYRes");
   unsigned int       nEchoes = processingControl->Value<int>("NumEchoes");
   unsigned int     nChannels = processingControl->Value<int>("NumChannels");
   unsigned int     numSlices = processingControl->Value<int>("NumSlices");
   size_t          frame_size = processingControl->Value<int>("AcquiredXRes");

   std::cout << "PSDname: " << const_cast<const GERecon::Legacy::LxDownloadData&>(*lxData).ImageHeaderData().psdname << std::endl;

   // TW: Some Slice debugging here with lots of noise
   std::cout << "numSlices = " << numSlices << std::endl;


   while (packetCount < packetQuantity)
   {
      // encoding IDs to fill ISMRMRD headers.
      int   sliceID = 0;
      int    viewID = 0;

      GERecon::Acquisition::FrameControlPointer const thisPacket = archiveStoragePointer->NextFrameControl();

      // Need to identify opcode(s) here that will mark acquisition / reference / control
      if (thisPacket->Control().Opcode() != GERecon::Acquisition::ScanControlOpcode)
      {

         // this is an issue here, because this limits this to ProgrammableControlPackets
         GERecon::Acquisition::ProgrammableControlPacket const packetContents = thisPacket->Control().Packet().As<GERecon::Acquisition::ProgrammableControlPacket>();

         viewID  = GERecon::Acquisition::GetPacketValue(packetContents.viewNumH,  packetContents.viewNumL);
         sliceID = GERecon::Acquisition::GetPacketValue(packetContents.sliceNumH, packetContents.sliceNumL);

         // output something if the opcode is not 1, to make it less noisy
         if(static_cast<uint8_t>(packetContents.opcode) != GERecon::Acquisition::ProgrammableOpcode)
            std::cout << "viewID = " << viewID << " sliceID = " << sliceID << " opcode = " << (int)(packetContents.opcode) << std::endl << std::flush;
         

         //std::cout << "sliceID = " << sliceID << " (" << static_cast<int>(packetContents.sliceNumH) << "," << static_cast<int>(packetContents.sliceNumL) << ")" << std::endl;
         if(static_cast<uint8_t>(packetContents.opcode) == GERecon::Acquisition::ProgrammableOpcode) {

            GERecon::Acquisition::CartesianFrameCommand F(packetContents);
            F.Dump(std::cout);
         }
         // TW: temporary hack, ignore the ViewCopy opcode
         //     I'm undecided on strategy here, but Packets with no data cause:
         //     - undefined behaviour/segfault on the ->Data() call (Boo GE)
         //     - we have nothing to write into the ISMRMRD file either, when we have no data 
         if(static_cast<uint8_t>(packetContents.opcode) != GERecon::Acquisition::ViewCopyOpcode) {
            if ((viewID < 1) || (viewID > nPhases))
            {
               acqType = GERecon::Acquisition::BaselineFrame;
               // nothing else to be done here for basic 2D case
            }    
            else
            {
               acqType = GERecon::Acquisition::ImageFrame;

               acqs.resize(dataIndex + 1);

               // TW build in some segmentation fault protection here, as Orchestra WILL segfault on the ->Data() call
               //    I'm interpreting a 0 rval as meaning that there is no kspace data in the packet as the 
               //    header states that the data size is [FrameSize x NumChannels x NumFrames] 
               if(thisPacket->NumFrames() == 0) {
                  throw std::runtime_error("Attempting to access data of a packet that has no data.");                  
               }
 
               // TW: Now we are seeing some C++11 suddenly some type inference
               auto kData = thisPacket->Data();

               // Grab a reference to the acquisition
               ISMRMRD::Acquisition& acq = acqs.at(dataIndex);

               // Set size of this data frame to receive raw data
               acq.resize(frame_size, nChannels, 0);
               acq.clearAllFlags();

               // Initialize the encoding counters for this acquisition.
               ISMRMRD::EncodingCounters idx;
               get_view_idx(processingControl, 0, idx);

               idx.slice                  = sliceID;
               idx.contrast               = packetContents.echoNum;
               idx.kspace_encode_step_1   = viewID - 1;

               acq.idx() = idx;

               // Fill in the rest of the header
               // acq.measurement_uid() = pfile->RunNumber();
               acq.scan_counter() = dataIndex;
               acq.acquisition_time_stamp() = time(NULL);
               for (int p=0; p<ISMRMRD::ISMRMRD_PHYS_STAMPS; p++) {
                  acq.physiology_time_stamp()[p] = 0;
               }
               acq.available_channels()   = nChannels;
               acq.discard_pre()          = 0;
               acq.discard_post()         = 0;
               acq.center_sample()        = frame_size/2;
               acq.encoding_space_ref()   = 0;
               // acq.sample_time_us()       = pfile->sample_time * 1e6;

               for (int ch = 0 ; ch < nChannels ; ch++) {
                  acq.setChannelActive(ch);
               }

               // Set first acquisition flag
               if (idx.kspace_encode_step_1 == 0)
                  acq.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE);

               // Set last acquisition flag
               if (idx.kspace_encode_step_1 == nPhases - 1)
                  acq.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

               if (processingControl->Value<bool>("ChopY") == 0) {
                  if (idx.kspace_encode_step_1 % 2 == 1) {
                     kData *= -1.0f;
                  }
               }

               for (int channelID = 0 ; channelID < nChannels ; channelID++)
               {
                  for (int i = 0 ; i < frame_size ; i++)
                  {
                     // The last dimension here in kData denotes the view
                     // index in the control packet that one must stride
                     // through to get data.  TODO - figure out if this
                     // can be programatically determined, and if so, use
                     // it. Will be needed for cases where multiple lines
                     // of data are contained in a single packet.
		               acq.data(i, channelID) = kData(i, channelID, 0);
                  }
               }

               dataIndex++;
            }
         } // if(static_cast<uint8_t>(packetContents.opcode) == 1)
      }

      packetCount++;
   }

   return acqs;
}

SEQUENCE_CONVERTER_FACTORY_DECLARE(QGenericConverter)

//} // namespace PfileToIsmrmrd

