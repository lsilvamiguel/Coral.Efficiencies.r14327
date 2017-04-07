// $Id: CheckTracks.h,v 1.15 2008/02/20 13:46:39 rgazda Exp $
#ifndef CheckTracks_h
#define CheckTracks_h

#include <vector>
#include <string>
#include <cstdarg>

/*!
   \file    CheckTracks.h
   \brief   To check coral output alignment tree
   \author  Hugo Pereira
   \version $Revision: 1.15 $
   \date    $Date: 2008/02/20 13:46:39 $
*/

class TChain;
class TPaveText;
class TPaveStats;
class TGraphErrors;
class TVirtualPad;
class Tracks;
class DetFileManager;

#include <TROOT.h>
#include <TObject.h>
#include <TPostScript.h>
#include <TCanvas.h>


/*! \class CheckTracks 
    \brief to check coral output alignment tree

    This class reads detector table and tracks in root tree format
    obtained from coral. It derives from TObject so that it can be used through root.
    It allows to monitor the quality of the alignment through partially automatized commands
    It can read a coral like option file
   
*/
class CheckTracks: public TObject {
  public:
  
  /*! \fn  CheckTracks( 
    const char* trackfileselection="", 
    const char* detectorfile="", 
    bool magnets_on=false );
    \brief constructor.
    \param trackfileselection input tree(s) obtained from coral. Wildcards are accepted.
    \param detectorfile the name of the det.dat file to be loaded
    \param magnets_on controls the expected structure of the tree. \warning It must match the option used in coral
  */
  CheckTracks( 
    const char* trackfileselection="", 
    const char* detectorfile="", 
    bool magnets_on=false );

  
  /*! \fn DetFileManager* LoadDetectorFile( const char* name )
     \brief Load detector.dat file, returns a pointer to the corresponding DetFileManager object
     \param name the name of the file to be loaded.
  */
  DetFileManager* LoadDetectorFile( const char* name );

  //! Returns pointer to the DetFileManager object, if any, 0 otherwise
  DetFileManager* GetDetectorFile( void ) { return df_;}
  
  
  /*! \fn Tracks* LoadTracks(  const char* fileselection, bool magnets_on = false )
     \brief load tracks contained in coral alignment tree output, returns pointer to the corresponding Tracks object.
     \param fileselection input tree(s). wildcards are accepted. Trees are chained.
     \param magnets_on controls the expected structure of the tree. \warning It must match the option used in coral
  */
  Tracks* LoadTracks( const char* fileselection, bool magnets_on = false );
  
  /*! \fn Tracks* AddToTracks( const char* fileselection )
     \brief add tracks contained in coral alignment tree to existing Tracks object, if any. Returns pointer to tracks object.
     \param fileselection input tree(s). wildcards are accepted. Trees are chained.
     \warning added trees must be of same type
  */   
  Tracks* AddToTracks( const char* fileselection );
  
  //! Returns pointer to the Tracks object, if any, 0 otherwise.
  Tracks* GetTracks( void ) { return tracks_; }
  
  //! create root canvas with \param width and \param height geometry, returns pointer to the canvas
  TCanvas* MakeCanvas( int width = 700, int height = 700 );
  
  //! returns pointer to root canvas, if any, 0 otherwise
  TCanvas* GetCanvas( void ) { return cv_; }

  //! Open postscript file of name \param name and type \param type. Returns pointer to the corresponding TPostScript object
  /* \fn TPostScript* OpenPS( const char* name = "checkTracks.ps", int type = 111 )
    \brief open PSFile for writting. Return pointer to TPostScript object.
    \param name the name of the psFile to be open.
    \param the postscript file type.
  */
  TPostScript* OpenPS( const char* name = "checkTracks.ps", int type = 111 );
  
  //! Close the postscript file, delete internal TPostscript object
  void ClosePS( void );
  
  //! return pointers to the TPostScript object if any, 0 otherwise
  TPostScript * GetPS( void ) { return ps_; }
  
  //! Option driven mode. \param file the name of the option file to be loaded example "opt.checkTracks"
  void DrawFromOpt( char *file );
  
  //! Create front page with usefull informations
  void MakeFrontPage( void );   

  /*! \fn void Draw( 
	      const char* var, 
	      const char* cut, 
				const char* opt="", 
				int nEnt=0,
				const char* titleX=0,
				const char* titleY=0 )
     \brief parse the trees loaded int Tracks object, draw according to the paramaters
     \param var expression of the tree parameters to be drawn
     \param cut expression of the tree parameters to be used for selection 
     \param opt root options. 
     \param nEnt the number of entries to be parsed. 0 means all
     \param titleX X axis title
		 \param titleY Y axis title
	*/
  using TObject::Draw; //jj to avoid compilation warning
  void Draw( 
	  const char* var, 
	  const char* cut, 
		const char* opt="", 
		int nEnt=0,
		const char* titleX=0,
		const char* titleY=0 );
  
  /*! \fn void DrawDet( 
	      const char* var, 
				const char* detselection, 
				const char* cut="", 
				const char* opt="", 
				int nEnt=0,
				const char* titleX = 0,
				const char* titleY = 0  )
     \brief parse the trees loaded int Tracks object, draw according to the paramaters. Makes a different plot (one plot per pad) for each detector 
     \param var expression of the tree parameters to be drawn
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut expression of the tree parameters to be used for selection 
     \param opt root options. 
     \param nEnt the number of entries to be parsed. 0 means all
     \param titleX X axis title
		 \param titleY Y axis title
  */
  void DrawDet( 
	  const char* var, 
	  const char* detselection, 
		const char* cut="", 
		const char* opt="", 
		int nEnt=0,
		const char* titleX=0,
		const char* titleY=0 );
  
  /*! \fn void DrawDU_LR( const char* detselection, const char* cut="", int nEnt=0, double OffsetN = 0, double OffsetP = 0 )
     \brief draw two DU histograms on the same plot. One for negative values of T_rVect, one for positives. Usefull for drift like detectors only 
     \param detselection detector selection, using shorten names, character wise wild cards example "DC01 DC0*Y DC01X1__"
     \param cut is used to select entries in the trees
     \param nEnt the max number of entries used.
     \param OffsetN value added to DU for T_RVect<0
     \param OffsetP value added to DU for T_RVect>0 
  */
  void DrawDU_LR( const char* detselection, const char* cut="", int nEnt=0, double OffsetN = 0, double OffsetP = 0 );

  /*! \fn void DrawDUvsU( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 )
     \brief draw DU vs U correlation plots for all detectors matching detselection. Makes one plot/detector.
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut is used to select entries in the trees
     \param opt is used for display options
     \param nEnt the max number of entries used.
  */
  void DrawDUvsU( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );  
  
  /*! \fn void DrawDUvsV( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 )
     \brief draw DU vs V correlation plots for all detectors matching detselection. Makes one plot/detector.
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut is used to select entries in the trees
     \param opt for display options
     \param nEnt the max number of entries used.
  */
  void DrawDUvsV( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );

  /*! \fn void DrawDUvsP( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 )
     \brief draw DU vs P correlation plots for all detectors matching detselection. Makes one plot/detector.
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut is used to select entries in the trees
     \param opt for display options
     \param nEnt the max number of entries used.
  */
  void DrawDUvsP( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );

 /* print main value of residuals vs momentum for selected detectors. It don't plot anything!!! */
 void PrintDUvsP(const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );
 void PrintDUvsU(const char* detselection, float du=100,float theta=10.0,float umin=-1000., float umax=1000., const char* cut="", const char* opt="box", int nEnt=0 );

  /*! \fn void DrawDU( const char* detselection, const char* cut="", int nEnt=0 )
     \brief draw DU, DU vs U and DU vs V correlation plots for all detectors matching detselection. Makes one plot/detector.
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut is used to select entries in the trees
     \param opt for display options
     \param nEnt the max number of entries used.
  */
  void DrawDU( const char* detselection, const char* cut="", int nEnt=0 );

  /*! \fn void DrawDUvsT (theta) const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 )
    \brief draw DU vs T correlation plots for all detectors matching detselection. Makes one plot/detector.
    \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
    \param cut is used to select entries in the trees
    \param opt for display options
    \param nEnt the max number of entries used.
    */
  void DrawDUvsT( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );

  /*! \fn void DrawDUvsT90 (theta+90) const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 )
    \brief draw DU vs T90 correlation plots for all detectors matching detselection. Makes one plot/detector.
    \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
    \param cut is used to select entries in the trees
    \param opt for display options
    \param nEnt the max number of entries used.
    */
  void DrawDUvsT90( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );

  /*! \fn void DrawDUvsZ( const char* detselection, const char* cut="", int nEnt=0 );
     \brief draw DU vs Z plot for all detectors matching detselection. 
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut is used to select entries in the trees
     \param nEnt the max number of entries used.
  */
  void DrawDUvsZ( const char* detselection, const char* cut="", int nEnt=0 );
  
  /*! \fn void DrawTrackMlt( const char* detselection, const char* cut="", const char* opt="", int nEnt=0 );
     \brief draw the number of detectors belonging to detSelection/track. 
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut is used to select entries in the trees
     \param opt for display options
     \int nEnt the max number of entries used.
  */
  void DrawTrackMlt( const char* detselection, const char* cut="", const char* opt="", int nEnt=0 );
  
  /*! \fn void DrawRT( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );
     \brief draw u_track - u_wire versus t_cluster for drift-like detectors
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut is used to select entries in the trees
     \param opt for display options
     \int nEnt the max number of entries used.
     \warning this method works only on drift-like detectors, that is Drift chambers (DC), W45 (DW), Drift Tubes (MB) and Straw tubes (ST)
  */
  void DrawRT( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );

  /*! \fn void FitDUvsU( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 )
     \brief fit DU vs U correlation plots for all detectors matching detselection with a straight line. Makes one plot/detector.
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut is used to select entries in the trees
     \param opt is used for display options
     \param nEnt the max number of entries used.
  */
  void FitDUvsU( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );  

  /*! \fn void FitDUvsV( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 )
     \brief fit DU vs V correlation plots for all detectors matching detselection with a straight line. Makes one plot/detector.
     \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
     \param cut is used to select entries in the trees
     \param opt is used for display options
     \param nEnt the max number of entries used.
  */
  void FitDUvsV( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );  

  /*! \fn void FitDUvsT( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 )
    \brief fit DU vs T (theta)  correlation plots for all detectors matching detselection with a straight line. 
    \Makes one plot/detector.
    \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
    \param cut is used to select entries in the trees
    \param opt is used for display options
    \param nEnt the max number of entries used.
    */
  void FitDUvsT( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );

  /*! \fn void FitDUvsT90( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 )
    \brief fit DU vs T (theta+90)  correlation plots for all detectors matching detselection with a straight line. 
    \Makes one plot/detector.
    \param detselection detector selection, using shorten names, character wise wild cards example "GM01 ST03Y2*b DC01X1__"
    \param cut is used to select entries in the trees
    \param opt is used for display options
    \param nEnt the max number of entries used.
    */
  void FitDUvsT90( const char* detselection, const char* cut="", const char* opt="box", int nEnt=0 );

  bool isBatch_; //!< if true, the input rootfiles are copied to local directory

  private:
  void _UpdateFromCanvas( void );             //!< Update Postcript file from canvas
  void _UpdateFromPad( TVirtualPad* pad );    //!< Update Postcript file from pad
  void _NewPage( void );                      //!< Check ps_ makes a new page
  void _MakeTitle( const char* format, ... ); //!< Add label to the top of each page
  void _PutText( TVirtualPad* pad, double x, double y, double size, const char* format, ... ); //! Add text in virtual pad
  TGraphErrors* _AddTGE( int color = 1 );     //!< add new TGraphErrors for DrawDUvsZ

  std::string glCut_;                            //!< global cut added on all single plot cuts
  int color_;                                    //!< color index for plots with option "same"
  TPostScript *ps_;                              //!< pointer to TPostscript object, if any
  TCanvas *cv_;                                  //!< pointer to TCanvas if any
  DetFileManager *df_;                           //!< pointer to DetFileManagerObject if any
  Tracks *tracks_;                               //!< pointer to Tracks object, if any
  std::vector<std::string> trackFileSelection_;  //!< fileselection entries stored at each call to AddToTracks
  std::vector<std::string> trackFiles_;          //!< files from which tracks are loaded
  std::string detectorFile_;                     //!< name of the detector.dat file
  ClassDef(CheckTracks,1)                        //!< Root Macro
};

#endif
