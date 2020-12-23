//
#ifndef __murat_ana_mva_data_hh__
#define __murat_ana_mva_data_hh__

namespace mu2e { 
  class MVATools;
};

class mva_data {
public: 

  struct data_t {
    const  char *fName;
    const  char *fXmlWeightsFile;
    const  char *fTrackCompHistFile;   // hist file name for "track_comp" job
    double       fMinValue;
    int          fBestID;
  };

  int              fTrainingCode;

  data_t           fData;

  mu2e::MVATools*  fMva;

  mva_data();

  mva_data(const char* Dataset, int MVATrainingCode);

  ~mva_data();

  const char* Name             () { return fData.fName;              }
  const char* XmlWeightsFile   () { return fData.fXmlWeightsFile;    }
  const char* TrackCompHistFile() { return fData.fTrackCompHistFile; }
  double      MinValue         () { return fData.fMinValue;          }
  int         BestID           () { return fData.fBestID;            }
  int         TrainingCode     () { return fTrainingCode;            }

  int         Init             ();
};

#endif
