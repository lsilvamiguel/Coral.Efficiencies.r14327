#include "config.h"
using CS::uint32;
using namespace std;
#include <vector>

#include <TROOT.h>
#include <TApplication.h>
#include <TVirtualX.h>

#include <TGListBox.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGMsgBox.h>
#include <TGMenu.h>
#include <TGCanvas.h>
#include <TGComboBox.h>
#include <TGTab.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>
#include <TGFileDialog.h>
#include <TGTextView.h>
#include <TGListView.h>
#include <TGSplitter.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <time.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TEnv.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TGStatusBar.h>
#include <TG3DLine.h>
#include <TRootHelpDialog.h>
#include "CollectErrs.h"

/******************************************************************************************
This class implements a sortable multicollumn list view. The collumn to sort is choosable
by header buttons. ROOT has something similar which drove me mad and did not work.
******************************************************************************************/

class SorTable : public TGCompositeFrame 
{
  public:  
  
  class Entry {
    public:
    enum Display {NONE,INT,DOUBLE,TEXT};
    Entry() {Disp=NONE;}
    Entry(Int_t vi) {IntVal=vi; Disp=INT;}
    Entry(Double_t vd) {DoubleVal=vd; Disp=DOUBLE;}
    Entry(const char *vt) {TextVal=vt; Disp=TEXT;}
    Entry(Int_t val,Int_t norm,bool normit);
    Int_t IntVal;
    Double_t DoubleVal;
    TString TextVal;
    Display Disp;
    
  };
  
  private:
  
  class Row :public vector<Entry> {
    public:
    Row(uint32 size,Int_t id):vector <Entry> (size) {Id=id;}
    Int_t Id;
  };
  class HCollumn : public TGCompositeFrame {
    public:
    HCollumn(TGWindow *p,uint32 id,char *txt,uint32 w,uint32 cw,bool sdir=true,char *tooltip=0);
    ~HCollumn() {delete Button;delete fLH;}
    void Associate(TGWindow *p) {Button->Associate(p);}
    uint32 Width;
    bool sortdir;
    TGLayoutHints *fLH;
    TGTextButton *Button;
    
  };
  struct scs {
    Entry *val;
    uint32 nr;
  };
      
public:   

  SorTable(TGWindow *p,UInt_t wi,UInt_t h,UInt_t opt,Int_t id=0, bool small=false);
  ~SorTable();
  void Associate(TGWindow *w) {MsgWindow=w;}
  void AddCollumn(char *title,uint32 width,bool sortdir=true,char *tooltip=0);
  void LastCollumn();
  void SortBy(uint32 col,bool dir);
  void ClearRows() {Data.clear();}
  void AddRow(Int_t Id=0);
  Entry &GetEntry(uint32 col,uint32 row);
  Int_t GetSelected() {Int_t sel=fLB->GetSelected();if (sel<(Int_t)Data.size() && sel >=0) return Data[sel].Id; else return -1;}
  void Update();
  Bool_t ProcessMessage(Long_t msg,Long_t parm1,Long_t parm2);
private:
  vector<HCollumn *> HCols;
  vector<Row > Data;

  Entry defdum;

  TGWindow *MsgWindow;

  uint32 sortby;
        
  TGCompositeFrame *fCF,*fCF1,*fCF2;
  TGTextButton *fTBFill;
  TGLayoutHints *fLH,*fLH1,*fLH2; 
  TGGC *tggc;
  FontStruct_t fnt;   
  uint32 cwidth;
  TGListBox *fLB;
  TTime lastclick;
  Long_t lastclicked;
  
  ULong_t color,red,orange,yellow,white;
    
  static int c_sort_h(const void *,const void *);

};

/******************************************************************************************
Simple ROOT based help dialog. When closed by user it deletes itself and sets the pointer 
to it pointed to by "me" to zero (to inform the owner).
******************************************************************************************/

class HelpDlg : public TRootHelpDialog {
  HelpDlg **Me;

  public:
  HelpDlg(TGWindow *main,HelpDlg **me,char *title,char *txt);
  ~HelpDlg() {(*Me)=0;}
};

/******************************************************************************************
This class provides the complete display window of MurphyTV.
******************************************************************************************/

class MurphyTV : public TGTransientFrame 
{
private:

/*--------------------------------------------------------
this class provides long and short description of each
error known to MurphyTV. These descriptions are defined
in the file errsdescr.h
--------------------------------------------------------*/
  class ErrZoo {
    struct Err {
      DaqErrorType type;
      const char *title;
      const char *descr;
    };
    enum orrg {NUM_ERRS = 110};
    static Err allerrs[NUM_ERRS];
    
    public:
    
    static uint32 GetNrErrs() {return NUM_ERRS;}
    static DaqErrorType GetType(uint32 err) {if (err<NUM_ERRS) return allerrs[err].type; else return DaqError::MAX_TYPE;}
    static const char *GetTitle(uint32 err) { if (err<NUM_ERRS) return allerrs[err].title; else return "   ";}
    static const char *GetDescription(uint32 err) { if (err<NUM_ERRS) return allerrs[err].descr; else return "   ";}
    
  };


/*--------------------------------------------------------
implementation of the "Catches" display 
--------------------------------------------------------*/
  class CatchList : public TGCompositeFrame 
  {
    enum knar {CB_SPECIAL,LB_CATCHES,LB_PORTS};
  public:
    CatchList(TGWindow *p,UInt_t w,UInt_t h,CollectErrs *ce,Int_t id=0, bool
    calibrationTrigger=false, bool small=false);
    ~CatchList();
    void Associate(TGWindow *w) {MsgWindow=w;}
    void Update(bool redraw,bool redraw2);
    void Update(bool redraw) {Update(redraw,redraw);}
    void SelectSpecial(DaqErrorType type);

   
    Bool_t ProcessMessage(Long_t msg,Long_t gparm1,Long_t parm2);
  private:

    const char *SrcType(uint32 src);
    void SrcInfos(uint32 src);
    TGWindow *MsgWindow;
    Int_t Id;
    CollectErrs *trash;
    DaqErrorType special;
    uint32 asource,dispsrc;
    FontStruct_t fnt;   
    uint32 cwidth;
    
    
    uint32 lastSpillNr;
    long lastSpillChange;
    ULong_t mycolor, mycolor2;
    uint32 count_of_no_data;
 
 
    SorTable *catchtab,*porttab;
    TGCompositeFrame *fCF,*fCF2;
    TGGroupFrame *fGF,*fGF2;

    TGLabel *fLabel,*fLpt;
    TGComboBox *fCombo;
    TGLabel *fLastUpdate;
    TGLabel *fRunNumber;
    TGLabel *fSpillNumber;
    TGCheckButton *fCheck;
    TGLabel *fCalibrationTrigger;
    ULong_t color,red,orange,yellow,white;
    TGTextView *fTVDets;
    TGLayoutHints *fLH,*fLH1,*fLH2,*fLH3,*fLH4, *fLH5;
  };
 
 
/*--------------------------------------------------------
implementation of the "errors" display
--------------------------------------------------------*/
  class ErrsList : public TGCompositeFrame 
  {
    enum knar {LB_ERRS,TE_SRC};
  public:
    ErrsList(TGWindow *p,UInt_t w,UInt_t h,CollectErrs *ce,Int_t id=0, bool small=false);
    ~ErrsList();
    void Associate(TGWindow *w) {MsgWindow=w;}
    void Update(bool redraw=true);
    void ChoseSource(uint32 src);
    
    Bool_t ProcessMessage(Long_t msg,Long_t gparm1,Long_t parm2);
  private:

    TGWindow *MsgWindow;
    Int_t Id;
    CollectErrs *trash;
    uint32 asource;
    Int_t disperr;

    SorTable *errtab;
    TGCompositeFrame *fCF;
    TGGroupFrame *fGF;
    TGTextView *fTVDescr;

    TGLabel *fLabel;
    TGTextEntry *fTESrc;
    TGCheckButton *fCBNorm,*fCBHide;

    TGLayoutHints *fLH,*fLH1,*fLH2,*fLH3;
  };


/*--------------------------------------------------------
implementation of the "monitor" display
--------------------------------------------------------*/

  class HistMonitor : public TGCompositeFrame
  {
    enum utkh {TB_AHA};
    class msg {
      public:
      msg(uint32 ei,float b) {errind=ei;pbad=b;}
      uint32 errind;
      float pbad;
    };
    typedef map <uint32,vector <msg> > RepStat;

  public:
    HistMonitor(TGWindow *p,UInt_t w,UInt_t h,CollectErrs *ce);
    ~HistMonitor();
    Bool_t ProcessMessage(Long_t msg,Long_t gparm1,Long_t parm2);
    void Update();
  private:
    CollectErrs *trash;
    RepStat repstat;
    uint32 DRun;
    bool Ackno;

    TGLayoutHints *fL1,*fL2,*fL3,*fL4;
    TGCompositeFrame *fCFc[3],*fCF2;
    TGGroupFrame *fGF;
    TGTextView *fTVMsg;
    FontStruct_t fnt;
    TGGC *tggc;
    TGTextButton *fTBAha;
    TGLabel *fLShit,*fLMsg;
    TGLabel *fLInfos[3][4];
  };



/*--------------------------------------------------------
implementation of the "Event Sizes" display
--------------------------------------------------------*/

  
  class HistSizes : public TGCompositeFrame
  {
    enum quiek { TE_CATCHSEL1}; 
  public:
    HistSizes(TGWindow *p,UInt_t w,UInt_t h,CollectErrs *ce);
    ~HistSizes();
    void Update();
    Bool_t ProcessMessage(Long_t msg,Long_t gparm1,Long_t parm2);
  private:
    uint32 SourceID1;    
    CollectErrs *trash;    
    TGLayoutHints *fL3,*fL4,*fL5;
    TGCompositeFrame *fF1;
    TRootEmbeddedCanvas *fEc1;
    TGTextEntry *fEdit1;
    TGLabel *fLabel1;
    TH1F *fHist1;              
  };  

/*--------------------------------------------------------
implementation of the "History" display
--------------------------------------------------------*/

  
  class HistHist : public TGCompositeFrame
  {
    enum quak { CB_ERRKAT2};//,CB_ERRKAT2b}; 
  public:
    HistHist(TGWindow *p,UInt_t w,UInt_t h,CollectErrs *ce);
    ~HistHist();
    void Update();
    Bool_t ProcessMessage(Long_t msg,Long_t gparm1,Long_t parm2);
  private:
    int errkat;
    
    CollectErrs *trash;
    
    TGLayoutHints *fL3,*fL4,*fL5;
    TGCompositeFrame *fF2;
    TRootEmbeddedCanvas *fEc2;
    TGComboBox *fCombo2;
    TGLabel *fLabel2;
    TH1F *fHist2;
    TGaxis *fXaxis2,*fYaxis2;
     
    float GetErrKatRate(int bin,int kat); 
    void ChoseTEHist(int kat);            
  };

/*--------------------------------------------------------
implementation of the "Data Flow" display
--------------------------------------------------------*/

  
  class HistHistD : public TGCompositeFrame
  {
  public:
    HistHistD(TGWindow *p,UInt_t w,UInt_t h,CollectErrs *ce);
    ~HistHistD();
    void Update();

  private: 
    CollectErrs *trash;    
    TGLayoutHints *fL3,*fL4,*fL5;
    TRootEmbeddedCanvas *fEc2;
    TH1F *fHist2;
    TGaxis *fXaxis2,*fYaxis2;        
  };  
  
/*----------------------------------------------------------
----------------------------------------------------------*/
    

  enum udfghkd{M_FILE_EXIT,M_HELP_ERRORS,M_HELP_GUI ,TB_CLOSE, CATCHLIST,ERRSLIST};

  private: 
    CatchList *catchlist;
    HistMonitor *monitorhist;
    
    ErrsList *errslist;
    HistHist *history;
    HistSizes *sizeshist;
    HistHistD *dsizeshist;
    
    TGLayoutHints *fL1,*fL2, *fL3, *fL4;
    TGTextButton *fTBClose;
    TGTab *fTab;
  
    TTimer *fTimer;
    int tdiv;
    CollectErrs *trash;
    static bool alive;
    Bool_t HandleTimer(TTimer *);

  public:
    MurphyTV(CollectErrs *ce,TGWindow *main, bool calibrationTrigger=false, bool small=false);
    ~MurphyTV();
    static bool IsOpen() {return alive;}
    Bool_t ProcessMessage(Long_t msg,Long_t gparm1,Long_t parm2);
    void CloseWindow();

};


