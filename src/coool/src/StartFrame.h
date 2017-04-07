#ifndef __StartFrame__
#define __StartFrame__

#include <TROOT.h>
#include <TApplication.h>
#include <TVirtualX.h>

#include <RQ_OBJECT.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGFileDialog.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TGPicture.h>

#include <vector>
#include <string>

#include "MainFrame.h"


/*! \brief Class corresponding to the start dialog window

    The start button is connected to MainFrame::Init(),
    where the StartFrame window is unmapped, and not deleted.

    \todo Handle root output file name

    \author Colin Bernet
*/

class StartFrame : public TObject {


 private:

  /// Associated window
  TGTransientFrame    *fMain;

  /// Geometry file name
  std::string fGeomFile;

  /// mapping file name
  std::string fMapFile;

  /// data file name. Streams are handled
  std::string fDataFile;

  /// root file name
  std::string fRootFile;

  /// groups description file name
  std::string fGroupFile;

  /// default parameters file name
  std::string fParamFile;

  /// user configuration file name
  std::string fConfigFile;

  /// MainFrame object this belongs to
  MainFrame *fMainframe;

  /// Monitor instance of the MainFrame -2b removed-
  Monitor *fMonitor;

  TGGroupFrame *fCheckFrame,*fTaskFrame,*fInputFrame,*fDecodingFrame;
  TGCompositeFrame *fDecodingFrameDown, *fDecodingFrameDown2;
  TGCompositeFrame *fRawFrame, *fMapFrame, *fGeomFrame, *fDownFrame, *fUpFrame;
  TGLayoutHints *fL1,*fL2, *fL3, *fL4, *fL5;
  TGLabel *fRawLabel, *fMapLabel, *fGeomLabel;
  TGButton *fRawButton, *fMapButton, *fGeomButton, *fStartButton;
  TGPicture *fCompassPic;
  std::vector<TGCheckButton*> fCheckButton;

  TGCheckButton          *fCheckTree;
  TGCheckButton          *fCheckText;
  TGCheckButton          *fCheckTrack;
  TGCheckButton          *fCheckCluster;
  TGCheckButton          *fCheckRef;
  TGCheckButton          *fCheckCalib;
  TGCheckButton          *fCheckReadFailedEvt;
  TGCheckButton          *fCheckExpert;

  TGVerticalFrame        *fMonDetFrame;
  TGVerticalFrame        *fAllNoneFrame;
  TGButton               *fDefButton;
  TGButton               *fAllButton;
  TGButton               *fNoneButton;
  TGButton               *fSetupButton;

  TGTextEntry            *fEvtMax;
  TGLabel                *fEvtLabel;

  TGTextEntry            *fMinSpacingEntry;
  TGLabel                *fMinSpacingLabel;

  TGTextEntry            *fRefEntry;
  TGLabel                *fRefLabel;

  /// file types for mapping files
  static const char* maptype[];

  /// file types for data files
  static const char* datatype[];

  /// file types for geometry files
  static const char* geomtype[];

 public:
  StartFrame(const TGWindow *p,MainFrame *main,
	     UInt_t w, UInt_t h);
  virtual ~StartFrame();

  void Hide() { fMain->UnmapWindow();}


  //Slots

  /// receives signal from fMain->SendCloseMessage()
  void CloseWindow() {
    exit(1);
  }

  /// All basic TGButton's connected to this slot
  void HandleButtons(int id=-1);

  /// TGCheckButton's connected to this slot
  void Checked(int id=-1);

  /// fCheckTree connected to this slot
  void CheckedTree();

  /// fCheckText connected to this slot
  void CheckedText();

  /// fCheckTrack connected to this slot
  void CheckedTrack();

  /// fCheckCluster connected to this slot
  void CheckedCluster();

  /// fCheckRef connected to this slot
  void CheckedRef();

  /// fCheckCalib connected to this slot
  void CheckedCalib();

  /// fCheckReadFailedEvt connected to this slot
  void CheckedReadFailedEvt();

  /// fCheckExpert connected to this slot
  void CheckedExpert();

  /// fDefButton connected to this slot
  void DefClicked();

  /// fAllButton connected to this slot
  void AllClicked();

  /// fNoneButton connected to this slot
  void NoneClicked();

  /// \return root file name
  std::string& GetRootFile() {return fRootFile;}
  /// \return data file name
  std::string& GetDataFile() {return fDataFile;}
  /// \return mapping file name
  std::string& GetMapFile() {return fMapFile;}
  /// \return geometry file name
  std::string& GetGeomFile() {return fGeomFile;}
  /// \return groups description file name
  std::string& GetGroupFile() {return fGroupFile;}
  /// \return parameters file name
  std::string& GetParamFile() {return fParamFile;}
  /// \return configuration file name
  std::string& GetConfigFile() {return fConfigFile;}
  /// \return number of events to be decoded
  int GetNEvent() {
    return atoi(fEvtMax->GetText());
  }
  /// \return minimum spacing between two decoded events in one spill
  int GetMinEventSpacing() {
    return atoi(fMinSpacingEntry->GetText());
  }
  /// \return return reference name entry
  const char* GetRefEntry() {
    return fRefEntry->GetText();
  }

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("StartFrame")
#else
  RQ_OBJECT()
#endif

  ClassDef(StartFrame,0)
};


#endif

