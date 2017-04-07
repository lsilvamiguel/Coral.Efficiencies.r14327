#ifndef __MainFrame__
#define __MainFrame__

#include <TROOT.h>
#include <TApplication.h>
#include <TVirtualX.h>

#include <TGListBox.h>
#include <TGListTree.h>
#include <TGCanvas.h>
#include <TRootEmbeddedCanvas.h>
#include <TGCanvas.h>
#include <TGFrame.h>
#include <TGMenu.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGFileDialog.h>
#include <TTimer.h>
#include <RQ_OBJECT.h>

#include <vector>
#include <list>

#include "Monitor.h"
#include "VariousSettings.h"

class StartFrame;
class TGImprovedPopupMenu;

/*! \brief Class corresponding to the main window of the gui Program.

    The window is not mapped in the constructor, but in MainFrame::Init().
    This forces the user to enter some input parameters via class StartFrame
    before doing anything.
    This class instanciates Monitor, the central monitoring class.

    \todo Open message window when then run number is changing 
    \todo Previous config files opened in file menu

    \author Colin Bernet, Damien Neyret
*/

class MainFrame : public TObject {

  //RQ_OBJECT()

 private:

  enum ETestCommandIdentifiers {
    M_FILE_OPEN,
    M_FILE_SAVE,
    M_FILE_SAVEAS,
    M_GUIFILE_SAVEAS,
    M_PLANEFILE_SAVEAS,
    M_GUIFILE_CREATEPS,
    M_ORDERLY_SHIFT_PLOTS,
    M_FILE_EXIT,
    //M_BUTTON_START,
    M_BUTTON_STOP,
    M_BUTTON_CLEAR,
    M_BUTTON_PANEL,
    M_BUTTON_PRINT,
    M_BUTTON_RESETHISTS,
    M_BUTTON_NEWCANVAS,
    M_BUTTON_PUTSELECTED,
    M_BUTTON_PUTINCANVAS,
    M_BUTTON_PROPERTIES,
    M_BUTTON_SAVE,
    M_BUTTON_RESET,
    M_LISTBOX_HISTO,
    M_CHECKBUTTON_STOPGO,
    M_CHECKBUTTON_RUNCLEAR
  };


  /// Associated window
  TGMainFrame        *fMain;

  // arguments
  int    fArgc;
  char **fArgv;

  /// Monitoring object
  Monitor *fMonitor;

  /// Setting files handling object
  VariousSettings *fVariousSettings;

  /// Start dialog window
  StartFrame *fStartFrame;

  /// Current histogram list -2bremoved-
  std::vector<TH1*>        fCurHists;

  /// Current histogram
  TH1*                fCurHisto;

  /// Current Plane
  Plane              *fCurPlane;

  /// Current Group
  Group              *fCurGroup;

  /// List of Group names
  std::vector<const char*> fGroupListVec;

  //Gui specific
  TGListBox          *fHistListBox;
  TGListTree         *fGroupListTree;
  TGCanvas           *fGroupLTCanvas;
  TTimer             *fTimer;
  TGCompositeFrame   *fCenterFrame,*fLeftFrame,
    *fPlaneFrame,*fHistFrame,*fGroupFrame;
  TGCompositeFrame    *fControlFrame;
  //TGCompositeFrame    *fControlDecoding;
  TGButton           *fClearButton;
  TGButton           *fResetHistsButton;
  //TGButton           *fStartButton;
  TGButton           *fStopButton;
  TGButton           *fPanelButton;
  TGButton           *fPrintButton;
  TGCheckButton      *fRunClearButton;
  TGCheckButton      *fStopGoButton;
  TGTextEntry            *fEvtCounter;
  TGButton	     *fNewCanvasButton;
  TGButton	     *fPutSelectedButton;
  TGButton	     *fPutInCanvasButton;
  TGButton           *fProperties;
  std::map<std::string,TCanvas*>	fCanvasMap;
  TRootEmbeddedCanvas *fCanvasWindow;
  TGMenuBar          *fMenuBar;
  TGPopupMenu        *fMenuFile;
  TGImprovedPopupMenu        *fMenuPrevious;
  TGImprovedPopupMenu        *fMenuPrtPrevious;
  TGLayoutHints      *fMenuBarLayout, *fMenuBarItemLayout;

  /// File types for Save as... menu item
  static const char* filetypes[];

  /// Fill or refill the 'open previous' menu
  void FillPreviousMenu();



 public:
  MainFrame(const TGWindow *p,int argc, char **argv, UInt_t w, UInt_t h);
  virtual ~MainFrame() { CloseWindow(); }

  /// Close this window
  void CloseWindow();

  /// All basic TGButtons connected here
  void HandleButtons(int id=-1);

  /// File Menu connected here
  void HandleMenu(int id);

  /// Previous open files menu connected here
  void HandlePreviousMenu(int id);

  /// Print previous open files menu connected here
  void HandlePrtPreviousMenu(int id);

  /// fHistListBox connected here
  void SelectedHist(int id);

  /// fDetListBox connected here
  //void SelectedDet(int id);

  /// fGroupListBox connected here
  void SelectedGroup(int id);

  /// fGroupListTree connected here
  void ClickedGroup(TGListTreeItem* item, int btn);
  void DoubleClickedGroup(TGListTreeItem* item, int btn);

  /// Builds window and instanciates Monitor fStartFrame Start button connected here
  void Init();

  /// timer for pads updates
  Bool_t HandleTimer(TTimer *t);

  int Argc() const {return fArgc;}

  char* Argv(int argi) {
    if(argi<fArgc)
      return fArgv[argi];
    else return 0;
  }

  /// \return associated root window
  TGMainFrame* GetFrame() {return fMain;}

  /// \return Monitor instance
  Monitor* GetMonitor() const {return fMonitor;}

#if ROOT_VERSION_CODE >= ROOT_30200
  RQ_OBJECT("MainFrame")
#else
  RQ_OBJECT()
#endif

  ClassDef(MainFrame,0)
};



class TGImprovedPopupMenu : public TGPopupMenu  {

  public:

    TGImprovedPopupMenu(const TGWindow* p, UInt_t w = 10, UInt_t h = 10, UInt_t options = 0)
         : TGPopupMenu(p, w, h, options) {};

    ~TGImprovedPopupMenu() { };

    void DeleteAllEntries()
        { TIter next(GetListOfEntries());
          TGMenuEntry *ptr;
          while ((ptr = (TGMenuEntry *) next()))
            DeleteEntry(ptr);
          fCurrent = 0;
          fBorderWidth = 3;
          fHeight      = 6;
          fWidth       = 8;
          fXl          = 16;
        }
};

#endif










