#ifndef __murat_ana_Index_t_hh__
#define __murat_ana_Index_t_hh__
//-----------------------------------------------------------------------------
// indices used to book histograms and form the names/titles
//-----------------------------------------------------------------------------
  struct EvtIndex_t {
    int sel;
  };
  
  struct TrkIndex_t {
    int sel  {-1};
    int slot {-1};                      // 0-17
    int plane{-1};                      // offline
    int panel{-1};                      // offline
    int pnl12{-1};                      // panel index within the station (0-11)
    int mnid {-1};
    int ch   {-1};                      // straw
  };
  
  struct CrvIndex_t {
    int sel   {-1};
    int sbid  {-1};
    int sipm  {-1};
    int och   {-1};                           // offline channel - 4*sbid+sipm
    int roc   {-1};                           // 1-18 ??? 
    int feb   {-1};                           // 1-24 in offline domain
    int ch    {-1};                           // channel within the FEB (0-63)
    int sector{-1};
    int module{-1};
    int layer {-1};
  };
  
  struct CaloIndex_t {
    int sel   {-1};
    int crate {-1};
    int disk  {-1};
    int cid   {-1};                   //
  };

#endif
