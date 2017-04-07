#ifndef Higz_h
#define Higz_h

/*! \file higz.h

    C++ interface to higz library.

    \sa http://wwwinfo.cern.ch/asdoc/higz/HIGZMAIN.html
    \author  Alexander Zvyagin
*/

#include "coral_config.h"
#if USE_HIGZ 

#include "cfortran.h"

////////////////////////////////////////////////////////////////////////////////

// Initialization
// 
//                          +----------------------+
//                          | CALL  IGINIT (NWHIGZ) |
//                          +----------------------+
//                                   
// 
// Action: This routine initializes HIGZ. This must be the first function
// to be used in the HIGZ package. 
// 
// Parameter Description: 
// 
// NWHIGZ
//        Minimal ZEBRA dynamic space in memory for the HIGZ division;
//        A value of 0, indicates that allocation will be done automatically.
//        NWHIGZ must be less than NWORDS-5000 where NWORDS is the size of the
//        common block PAWC (see below). 
// 
// The ZEBRA memory allocation must be defined in the application program
// with the common block: 
// 
//       COMMON/PAWC/RPAW(NWORDS)
// 
// If HIGZ is used outside the context of PAW the routine MZPAW must be
// called in the main program in order to initialize the ZEBRA package
// [bib-ZEBRA], before calling IGINIT. Note that packages like
// HBOOK[bib-HBOOK], HPLOT[bib-HPLOT], PAW[bib-PAW] and KUIP[bib-KUIP] call
// MZPAW directly and therefore the user should not issue such a call.
// These packages store dynamic structures in the same common /PAWC/. 
// 
//       CALL MZPAW(NWORDS,'M')
// 

extern "C" void iginit_(const int &NWHIGZ);
inline void IGINIT(int NWHIGZ)
{iginit_(NWHIGZ);}

////////////////////////////////////////////////////////////////////////////////

// Termination
// 
//                              +--------------+
//                              | CALL  IGEND   |
//                              +--------------+
//                                   
// 
// Action: This routine terminates HIGZ. This must be the last call to be issued
// in a HIGZ session. IGEND deactivates and closes all open workstations. It
// also closes the basic graphics package by calling IDAWK, ICLWK, ICLKS. 

extern "C" void igend_ (void);
inline void IGEND(void)
{igend_();}

////////////////////////////////////////////////////////////////////////////////

// The Generic Routine
// 
//               +--------------------------------------------+
//               |CALL  IRQLC (KWKID,LCDNR,ISTAT*,NT*,PX*,PY*) |
//               +--------------------------------------------+
//                                   
// 
// Action: This routine returns the (x,y) position of the cursor in world
// coordinates, and the index the normalization transformation. Its calling
// sequence is compatible with the equivalent GKS routine. Parameter Description:
// 
// KWKID
//        Workstation identifier. 
// LCDNR
//        Locator device. 
//        1
//               Keyboard. 
//        2
//               Graphic tablet. 
//        With the X11 driver LCDNR can have the following values: 
//        10
//               tracking cross 
//        20
//               cross-hair 
//        30
//               rubber circle 
//        40
//               rubber band 
//        50
//               rubber rectangle 
//        99
//               the screen coordinates are taken in XLOC and YLOC. 
//        >0
//               request mode 
//        <0
//               sample mode 
// ISTAT
//        Return status. 
//        0
//               Graphic input has been canceled. 
//        1
//               A point was located and its coordinates are recorded in PX and PY. 
// NT
//        Index of the normalization transformation. 
// PX
//        X coordinate of position of locator 
// PY
//        Y coordinate of position of locator 

extern "C" void irqlc_ (const int &KWKID,const int &LCDNR,int &ISTAT,int &NT,float &PX,float &PY);
inline void IRQLC(int KWKID,int LCDNR,int &ISTAT,int &NT,float &PX,float &PY)
{irqlc_(KWKID,LCDNR,ISTAT,NT,PX,PY);}

////////////////////////////////////////////////////////////////////////////////

// The Two Points Routine
// 
//          +------------------------------------------------------+
//          |CALL  IGLOC2 (KWKID,*NT*,X1*,Y1*,X2*,Y2*,ISTAT*,CHOPT) |
//          +------------------------------------------------------+
//                                   
// 
// Action: This routine returns the graphic cursor position in world coordinates
// space of two points and the corresponding normalization transformation
// number. Rubberbanding is used to visualize the area (box) delimited by the
// two points. Parameter Description: 
// 
// KWKID
//        Workstation identifier 
// NT
//        Index of the normalization transformationsee(CHOPT). 
// X1
//        X coordinate of the cursor position in world coordinates space of the first point. 
// Y1
//        Y coordinate of the cursor position in world coordinates space of the first point. 
// X2
//        X coordinate of the cursor position in world coordinates space of the second point. 
// Y2
//        Y coordinate of the cursor position in world coordinates space of the second point. 
// ISTAT
//        Return status: 
//        0
//               Graphic input has been canceled. 
//        1
//               Two points were located and their coordinates are recorded in X1, Y1, X2, Y2. 
// CHOPT
//        CHARACTER variable specifying the option desired: 
//        ' '
//               NT is an output parameter. 
//        'P'
//               NT is an input and output parameter. In this case, NT contains on input the
//               normalization transformation index with the highest priority.

extern "C" void igloc2_(const int &KWKID,int &NT,float &X1,float &Y1,float &X2,float &Y2,int &ISTAT,const char *CHOPT,int);
inline void IGLOC2(int KWKID,int &NT,float &X1,float &Y1,float &X2,float &Y2,int &ISTAT,const char *CHOPT)
{igloc2_(KWKID,NT,X1,Y1,X2,Y2,ISTAT,CHOPT,strlen(CHOPT));}

////////////////////////////////////////////////////////////////////////////////

// Polymarker colour index.
// 
//                           +--------------------+
//                           |CALL  ISPMCI (ICOLI) |
//                           +--------------------+
//                                   
// 
// Action: This routine sets the polymarker colour index attribute for use by
// future invocations of IPM. The routine IGSET (see section [more info]) can also
// be used with the parameter PMCI. Parameter Description: 
// 
// ICOLI
//        Polymarker colour index. 

extern "C" void ispmci_(const int &ICOLI);
inline void ISPMCI(int ICOLI)
{ispmci_(ICOLI);}

////////////////////////////////////////////////////////////////////////////////

// Marker type
// 
//                            +------------------+
//                            |CALL  ISMK (MTYPE) |
//                            +------------------+
//                                   
// 
// Action: This routine sets the marker type attribute for use by future invocations
// of IPM. All workstations support at least the marker types 1 through 5
// (see below). More marker types may be supported by the underlying graphics
// package. Marker types 20 to 31 are also defined, according to the figure
// [more info], and are independent from the underlying graphics package used.
// If a requested marker type is not supported on a workstation, marker type 1
// (a point) is used when polymarkers are created. The routine IGSET (see section
// [more info]) can also be used with the parameter MTYP. Parameter
// Description: 
// 
// MTYPE
//        Marker type (positive number) 
//        1
//               Point shape (#). 
//        2
//               Plus shape (+). 
//        3
//               Asterisk shape (#). 
//        4
//               Circle shape ( o). 
//        5
//               X shape (x).

extern "C" void ismk_(const int &MTYPE);
inline void ISMK(int MTYPE)
{ismk_(MTYPE);}

////////////////////////////////////////////////////////////////////////////////

// Marker scale factor.
// 
//                           +--------------------+
//                           | CALL  ISMKSC (SSFM) |
//                           +--------------------+
//                                   
// 
// Action: This routine sets the marker scale factor. This scale factor is applied
// on the nominal size of the marker. On all workstation, except PostScript
// files, the marker type 1 is not scalable. The routine IGSET (see section
// [more info]) can also be used with the parameter MSCF. Parameter Description: 
// 
// SSFM
//        Scale factor applied to markers. (>=0.) 
// 
// 
//                                                                                        
//  
// 
// Figure:  HIGZ Marker type (20-31).
// 
//  [MARKER-TYPE]
// 
//                                                                                        
//  
// 
// Figure:  Examples marker scale factor.
// 
//  [MARKER-SIZE]

extern "C" void ismksc_(const float &SSFM);
inline void ISMKSC(float SSFM)
{ismksc_(SSFM);}

////////////////////////////////////////////////////////////////////////////////

// Workstation clear
// 
//                        +--------------------------+
//                        | CALL  ICLRWK (KWKID,KOFL) |
//                        +--------------------------+
//                                   
// 
// Action: This routine clears the output area of a workstation which has been
// previously opened. Parameter Description: 
// 
// KWKID
//        Workstation identifier. On a softcopy device (e.g. a terminal), the output
//        area is cleared. On a hardcopy device, the paper is advanced, so that a
//        fresh area is available for drawing. If KWKID =0 then all active
//        workstations are cleared. 
// KOFL
//        Flag controlling the operation of routine ICLRWK on a workstation for
//        which the output area is already cleared. Possible values are: 
//        0
//               If there has been no output since the previous ICLRWK, nothing happens. 
//        1
//               The output medium is advanced or cleared in any cases. 
// 
// If a change has been requested in the workstation transformation (via ISWKVP or ISWKWN),
// the workstation transformation is recalculated when
// ICLRWK is called. With the GPR, GL, and X11 versions of HIGZ, if the window size has
// changed, the new size will be automatically taken into account
// after a clear workstation. 

extern "C" void iclrwk_(const int &KWKID,const int &KOFL);
inline void ICLRWK(int KWKID,int KOFL)
{iclrwk_(KWKID,KOFL);}

////////////////////////////////////////////////////////////////////////////////

// Clipping
// 
//                           +--------------------+
//                           |CALL  ISCLIP (ICLSW) |
//                           +--------------------+
// 
//                                   
// 
// Action: This routine sets the ``clipping indicator'' for use by future invocations
// of IFA, IPL, IPM and ITX. The clipping indicator specifies where primitives
// should be clipped. Parameter Description: 
// 
// ICLSW
//        Clipping indicator 
//        1
//               Primitives should be clipped at the boundary of the
//               normalization transformation viewport. 
//        0
//               Primitives should be clipped at the edge of the
//               normalized device coordinates space. 

extern "C" void isclip_(const int &ICLSW);
inline void ISCLIP(int ICLSW)
{isclip_(ICLSW);}

////////////////////////////////////////////////////////////////////////////////

// Setting attributes
// 
//                         +------------------------+
//                         |CALL  IGSET (CHNAME,VAL) |
//                         +------------------------+
//                                   
// 
// Action: Routine used to set the value of attributes related to primitives
// and/or macroprimitives. The first parameter is the mnemonic name of the
// parameter, the second is the value to be assigned. Note that all the basic
// primitives attributes can also be set with this routine. 
// 
// CHNAME
//        Character variable specifying the name of the parameter to be set
//        (type CHARACTER*4). This is an UPPERCASE character string. 
// VAL
//        Floating point value of the parameter (must be specified as a REAL number).
//        A value of 0.0 indicates that the parameter value must be reset to its
//        default value. 
// 
// 
// 
// +-CHNAME--+-----------------------------VAL------------------------------##
// +---------+--------------------------------------------------------------##
// | 'FAIS'  | Fill Area Interior Style (0.,1.,2.,3.). See  ISFAIS           ##
// | 'FASI'  | Fill Area Style Index. See  ISFASI                            ##
// | 'LTYP'  | Line TYPe. See  ISLN                                          ##
// | 'BASL'  | BAsic Segment Length. See  ISLN                               ##
// | 'LWID'  | Line WIDth. See  ISLWSC                                       ##
// | 'MTYP'  | Marker TYPe. See  ISMK                                        ##
// |         |                                                              ##
// | 'MSCF'  | Marker SCale Factor. See  ISMKSC                              ##
// | 'PLCI'  | PolyLine Colour Index. See  ISPLCI                            ##
// | 'PMCI'  | PolyMarker Colour Index. See  ISPMCI                          ##
// | 'FACI'  | Fill Area Colour Index. See  ISFACI                           ##
// | 'TXCI'  | TeXt Colour Index. See  ISTXCI                                ##
// | 'TXAL'  | 10*(horizontal alignment) + (vertical alignment). See        ##
// |         |  ISTXAL                                                       ##
// |         |                                                              ##
// | 'CHHE'  | CHaracter HEight. See  ISCHH                                  ##
// | 'TANG'  | Text ANGle (used to calculate the Character up vector). See  ##
// |         |  ISCHUP                                                       ##
// | 'TXFP'  | 10*(TeXt Font) + (TeXt Precision). See  ISTXFP                ##
// | 'TMSI'  | Tick Marks SIze (in world coordinates). See  IGAXIS           ##
// | 'LASI'  | LAbels SIze (in world coordinates). See  IGAXIS               ##
// | 'LAOF'  | LAbels OFfset. See  IGAXIS                                    ##
// |         |                                                              ##
// | 'AWLN'  | Axis Wire LeNght. See  IGAXIS                                 ##
// | 'PASS'  | Text width (given by number of PASSes) of characters drawn   ##
// |         | by  IGTEXT. The width is simulated by shifting the ``pen''    ##
// |         | slightly at each pass.                                       ##
// | 'CSHI'  | Distance between each shifted drawing of the character (in   ##
// |         | percentage of the character height) for characters drawn by  ##
// |         |  IGTEXT                                                       ##
// | 'BORD'  |   0.    The border in  IGBOX,  IGFBOX and  IGARC is not drawn.  ##
// |         |   1.    The border in  IGBOX,  IGFBOX and  IGARC is drawn.      ##
// | 'PICT'  | Starting number for the automatic naming of pictures.        ##
// |         |                                                              ##
// | 'AURZ'  | 1. The last current picture is automatically saved on disk   ##
// |         | when a new picture is created see  IZPICT.                    ##
// |  '*'    | All attributes are set to their default values.              ##
// | 'SHOW'  | The current value and the default of the parameters          ##
// |         | controlled by  IGSET are displayed.                           ##
// | 'BARO'  | Offset of the left edge of the bar with respect to the left  ##
// |         | margin of the bin for a bar chart (expressed as a fraction   ##
// |         | of the bin width). See  IGHIST                                ##
// | 'BARW'  | Width of the bar in a bar chart (expressed as a fraction of  ##
// |         | the bin width). See  IGHIST                                   ##
// | 'NCOL'  | Number of entry in the COLour map.                           ##
// |         |                                                              ##
// +-'CLIP'--+-Clipping-mode:-1.=on-0.=off----------------------------------##
// +-CHNAME--+-----------------VAL-(For-X11-interface-only)-----------------##
// | 'DRMD'  | Drawing mode: 1.=copy 2.=xor 3.=invert                       ##
// | 'SYNC'  | Synchronise the graphics in X11 1.=yes 0.=no                 ##
// | '2BUF'  | 10*(WKID)+(double buffer mode: 1.=on 0.=off)                 ##
// +---------+--------------------------------------------------------------##
// 
// Table:  Overview of  IGSET parameters

extern "C" void igset_(const char *CHNAME,const float &VAL,int);
inline void IGSET(const char *CHNAME,float VAL)
{char s[5]="    "; for( int i=0; i<4 && CHNAME[i]; i++ ) s[i]=CHNAME[i]; igset_(s,VAL,4);}

////////////////////////////////////////////////////////////////////////////////

// Colour representation
// 
// Each colour is defined by an index and percentages of red, green and blue.
// Once a colour is defined it can be used via a reference to its index. If a
// requested colour index is not available on a workstation, colour index 1 is
// used when primitives are created. 
// 
//                     +--------------------------------+
//                     | CALL  ISCR (KWKID,ICI,CR,CG,CB) |
//                     +--------------------------------+
//                                   
// 
// Action: This routine sets the colour representation (red/green/blue) of the
// colour index on a previously opened workstation. On workstations using colour
// tables, this function can change the image immediately. On workstations lacking
// such tables, this new colour definition will be taken into account in the
// next use of this colour. Parameter Description: 
// 
// KWKID
//        Workstation identifier 
// ICI
//        Colour index. 
// CR
//        Intensity of red 0.<=CR< =1. 
// CG
//        Intensity of green 0.<=CG< =1. 
// CB
//        Intensity of blue 0.<=CB< =1. 
// 
// By default the first eight colour indices are defined as follows: 
// 
//    +-+--------+----------------------------##------+--------+-------+-+
//    +-+-Index--+-Colour---------------------##-Red--+-Green--+-Blue--+-+
//    | |   0    | Background colour (White)  ## 1.   |  1.    |  1.   | |
//    | |   1    | Foreground colour (Black)  ## 0.   |  0.    |  0.   | |
//    | |   2    | Red                        ## 1.   |  1.    |  1.   | |
//    | |   3    | Green                      ## 0.   |  1.    |  0.   | |
//    | |   4    | Dark blue                  ## 0.   |  0.    |  1.   | |
//    | |        |                            ##      |        |       | |
//    | |   5    | Yellow                     ## 1.   |  1.    |  0.   | |
//    | |   6    | Magenta (red-purple)       ## 1.   |  0.    |  1.   | |
//    +-+---7----+-Cyan-(light-blue)----------##-0.---+--1.----+--1.---+-+
//                                   
// 
// When a PostScript file is printed on a black and white PostScript printer,
// a grey level simulation of the colours is used according to the figure [more info].

extern "C" void iscr_(const int &KWKID,const int &ICI,const float &CR,const float &CG,const float &CB);
inline void ISCR(int KWKID,int ICI,float CR,float CG,float CB)
{iscr_(KWKID,ICI,CR,CG,CB);}

////////////////////////////////////////////////////////////////////////////////

// Normalization transformation selection
// 
//                            +------------------+
//                            | CALL  ISELNT (NT) |
//                            +------------------+
//                                   
// 
// Action: This routine selects the normalization transformation to be used when
// world coordinates must be mapped to or from normalized device coordinates
// (NDC). These mappings usually take place during invocations of primitives
// ( IFA, IPL, IPM, and ITX) and during graphics input ( IRQLC). Transformation
// 0 always has a window and a viewport that are the unit square (0.-1. by 0.-1.)
// and cannot be changed with ISVP or ISWN. Transformation 0 is selected
// by default. Parameter Description: 
// 
// NT
//        Normalization transformation index (0 

extern "C" void iselnt_(const int &NT);
inline void ISELNT(int NT)
{iselnt_(NT);}

////////////////////////////////////////////////////////////////////////////////

// Fill area colour index.
// 
//                           +--------------------+
//                           |CALL  ISFACI (ICOLI) |
//                           +--------------------+
//                                   
// 
// Action: This routine sets the fill area colour index attribute for use by future
// invocations of IFA. The routine IGSET (see section [more info]) can also be
// used with the parameter FACI. Parameter Description: 
// 
// ICOLI
//        Fill area colour index.

extern "C" void isfaci_(const int &ICOLI);
inline void ISFACI(int ICOLI) 
{isfaci_(ICOLI);}

////////////////////////////////////////////////////////////////////////////////

// Fill area interior style
// 
//                           +--------------------+
//                           | CALL  ISFAIS (INTS) |
//                           +--------------------+
//                                   
// 
// Action: This routine sets the fill area interior style attribute for use by future
// invocations of IFA. The routine IGSET (see section [more info]) can also be
// used with the parameter FAIS. Parameter Description: 
// 
// INTS
//        Fill area interior style. Possible values are: 
//        0
//               Hollow: the perimeter of the filled area, after clipping, is
//               drawn using solid lines. 
//        1
//               Solid: the area is filled solidly. 
//        2
//               Pattern: the area is filled with a dot-dashed pattern. 
//        3
//               Hatched: the area is filled according to the current value
//               of the fill area style index.

extern "C" void isfais_(const int &INTS);
inline void ISFAIS(int INTS)
{isfais_(INTS);}

////////////////////////////////////////////////////////////////////////////////

// Drawing a box
// 
//                        +--------------------------+
//                        | CALL  IGBOX (X1,X2,Y1,Y2) |
//                        +--------------------------+
//                                   
// 
// Action: This routine fills a rectangle according to the ``fill area colour index''
// (see section [more info]), the ``fill area interior style'' (see section [more
// info]), and the ``fill area style index'' (see section [more info]) attributes.
// The border is never drawn unless the interior style is hollow or the routine
// IGSET has been called with 'BORD' and VAL = 1.. As it is shown on the figure [more info],
// the border of the rectangle is drawn according to the values
// of the ``line width scale factor'' (see section [more info]) and the ``polyline colour
// index'' (see section [more info]) attributes, whereas the ``line type''
// is always solid (see section [more info]). Parameter Description: 
// 
// X1
//        X coordinate of 1st corner of the rectangle in WC. 
// X2
//        X coordinate of 2nd corner of the rectangle in WC. 
// Y1
//        Y coordinate of 1st corner of the rectangle in WC. 
// Y2
//        Y coordinate of 2nd corner of the rectangle in WC. 

extern "C" void igbox_(const float &X1,const float &X2,const float &Y1,const float &Y2);
inline void IGBOX(float X1,float X2,float Y1,float Y2)
{igbox_(X1,X2,Y1,Y2);}

////////////////////////////////////////////////////////////////////////////////

// Normalization Transformation window definition
// 
// 
//                   +------------------------------------+
//                   | CALL  ISWN (NT,XMIN,XMAX,YMIN,YMAX) |
//                   +------------------------------------+
//                                   
// 
// Action: This routine sets the boundaries of the window of a normalization
// transformation. The window must be specified in world coordinates. The
// boundaries of the window, together with the boundaries of the viewport
// (which are in normalized device coordinates) determine a transformation from
// world coordinates to normalized device coordinates consisting of separate
// X and Y scale factors and a translation in two dimensions. The normalization
// transformation is selected by using routine ISELNT. Parameter Description: 
// 
// NT
//        Normalization transformation index (0XMIN
//        X coordinate of the lower left hand corner in WC space. 
// XMAX
//        X coordinate of the upper right hand corner in WC space. 
// YMIN
//        Y coordinate of the lower left hand corner in WC space. 
// YMAX
//        Y coordinate of the upper right hand corner in WC space. 
// 
// The last four parameters must satisfy the conditions XMIN < XMAX and YMIN < YMAX. 

extern "C" void iswn_(const int &NT,const float &XMIN,const float &XMAX,const float &YMIN,const float &YMAX);
inline void ISWN(int NT,float XMIN,float XMAX,float YMIN,float YMAX)
{iswn_(NT,XMIN,XMAX,YMIN,YMAX);}

////////////////////////////////////////////////////////////////////////////////

// Normalization Transformation viewport definition
// 
//                   +------------------------------------+
//                   | CALL  ISVP (NT,XMIN,XMAX,YMIN,YMAX) |
//                   +------------------------------------+
//                                   
// 
// Action: This routine sets the boundaries of the viewport of a normalization
// transformation. The viewport must be specified in normalized device
// coordinates. The boundaries of the viewport have two roles: 
// 
//     1.Together with the boundaries of the window (which are in world coordinates)
//     they determine a transformation from world coordinates to
//        normalized device coordinates consisting of separate X and Y scale factors
//        and a translation in two dimensions. 
//     2.When the clipping indicator is 1 (see routine ISCLIP), primitives are
//     clipped to the boundary of the viewport (once the primitives are transformed
//        to normalized device coordinates) 
// 
// The normalization transformation is selected with the routine ISELNT. Parameter Description: 
// 
// NT
//        Normalization transformation index (0XMIN
//        X coordinate of the lower left hand corner in DC space (0.0<=XMIN< =1.0). 
// XMAX
//        X coordinate of the upper right hand corner in DC space (0.0<=XMAX< =1.0). 
// YMIN
//        Y coordinate of the lower left hand corner in DC space (0.0<=YMIN< =1.0). 
// YMAX
//        Y coordinate of the upper right hand corner in DC space (0.0<=YMAX< =1.0). 
// 
// The last four parameters must satisfy the conditions XMIN < XMAX and YMIN < YMAX.

extern "C" void isvp_(const int &NT,const float &XMIN,const float &XMAX,const float &YMIN,const float &YMAX);
inline void ISVP(int NT,float XMIN,float XMAX,float YMIN,float YMAX)
{isvp_(NT,XMIN,XMAX,YMIN,YMAX);}

////////////////////////////////////////////////////////////////////////////////

// Polyline colour index.
// 
//                           +--------------------+
//                           |CALL  ISPLCI (ICOLI) |
//                           +--------------------+
//                                   
// 
// Action: This routine sets the polyline colour index attribute for use by future
// invocations of IPL. The routine IGSET (see section [more info]) can also be
// used with the parameter PLCI. Parameter Description: 
// 
// ICOLI
//        Polyline colour index.

extern "C" void isplci_(const int &ICOLI);
inline void ISPLCI(int ICOLI)
{isplci_(ICOLI);}

////////////////////////////////////////////////////////////////////////////////

// Update workstation
// 
//                         +------------------------+
//                         |CALL  IUWK (KWKID,IRFLG) |
//                         +------------------------+
//                                   
// 
// Action: This routine updates the workstation KWKID. It send all buffered output
// to the screen. In the X11 version of HIGZ, this routine allows to flush
// the X11 buffer. This routine is usually called with the first parameter equal
// to 0 and the second to 1. Parameter Description: 
// 
// KWKID
//        Workstation identifier. KWKID = 0 updates all the current open workstations. 
// IRFLG
//        Regeneration flag: 
//        0
//               postpone update workstation (only when the underlying graphics package is GKS) 
//        1
//               refresh entire display 
//        2
//               update current view

extern "C" void iuwk_(const int &KWKID,const int &IRFLG);
inline void IUWK(int KWKID,int IRFLG)
{iuwk_(KWKID,IRFLG);}

////////////////////////////////////////////////////////////////////////////////

// Line type.
// 
//                            +------------------+
//                            |CALL  ISLN (LTYPE) |
//                            +------------------+
//                                   
// 
// Action: This routine sets the line type attribute for use by future invocations
// of IPL. All workstations support at least line types 1 through 4 (see figure
// [more info]). Other line types may be supported. If a requested line type is not
// supported on a workstation, line type 1 is used when polylines are created.
// The routine IGSET (see section [more info]) can also be used with the parameter
// LTYP. Parameter Description: 
// 
// LTYPE
//        Line type (positive number). 
//        1
//               Solid lines 
//        2
//               Dashed lines 
//        3
//               Dotted lines 
//        4
//               Dashed-dotted lines 
// 
// Note that line type values are dependent upon the underlying graphics package used.
// For the user's convenience, HIGZ defines a number of line types, indicated in the
// figure [more info], which are independent from the basic graphics package used. 

extern "C" void isln_(const int &LTYPE);
inline void ISLN(int LTYPE)
{isln_(LTYPE);}

////////////////////////////////////////////////////////////////////////////////

// Drawing an arc
// 
//                 +----------------------------------------+
//                 | CALL  IGARC (XC,YC,R1,R2,PHIMIN,PHIMAX) |
//                 +----------------------------------------+
//                                   
// 
// Action: This routine draws one or two arcs of a circle. If the two radii are not
// equal the area between the two arcs is filled according to the fill area interior
// style index and the fill area style index. The border is never drawn unless the
// interior style is hollow or the routine IGSET has been called with BORD and
// VAL = 1. If the arc's radii are equal only one arc is drawn. Parameter Description: 
// 
// XC
//        X coordinate of the arc's center in world coordinate space. 
// YC
//        Y coordinate of the arc's center in world coordinate space. 
// R1
//        Radius of first arc. 
// R2
//        Radius of second arc. 
// PHIMIN
//        Starting angle (degrees.) 
// PHIMAX
//        Final angle (degrees.) 

extern "C" void igarc_(const float &XC,const float &YC,const float &R1,const float &R2,const float &PHIMIN,const float &PHIMAX);
inline void IGARC(float XC,float YC,float R1,float R2,float PHIMIN,float PHIMAX)
{igarc_(XC,YC,R1,R2,PHIMIN,PHIMAX);}

////////////////////////////////////////////////////////////////////////////////

// Line width scale factor.
// 
//                           +--------------------+
//                           |CALL  ISLWSC (WIDTH) |
//                           +--------------------+
//                                   
// 
// Action: This routine sets the width of a line for use by future invocations
// of the polyline drawing routine IPL. The actual line width is determined by a
// nominal line width (workstation-dependent) multiplied by the line width scale
// factor. The nominal line width is one pixel on screens. On PostScript
// printers the nominal line width is one ``dot''. Therefore the width of a line
// can vary from a printer to another depending on the printer definition (300 dots
// per inch, 400 dots per inch etc.). The figure [more info] shows some examples of
// various line width. The routine IGSET (see section [more info]) can also
// be used with the parameter LWID. Parameter Description: 
// 
// WIDTH
//        Line width scale factor.

extern "C" void islwsc_(const float &WIDTH);
inline void ISLWSC(float WIDTH)
{islwsc_(WIDTH);}

////////////////////////////////////////////////////////////////////////////////

// Drawing axes
// 
//             +------------------------------------------------+
//             | CALL  IGAXIS (X0,X1,Y0,Y1,WMIN,WMAX,NDIV,CHOPT) |
//             +------------------------------------------------+
//                                   
// 
// Action: This routines allows the user to draw axes on a picture. Parameter Description: 
// 
// X0
//        X coordinate of the origin of the axis in world coordinates space. 
// X1
//        X coordinate of the end of the axis in world coordinates space. 
// Y0
//        Y coordinate of the origin of the axis in world coordinates space. 
// Y1
//        Y coordinate of the end of the axis in world coordinates space. 
// WMIN
//        Lowest value for the tick mark labels written on the axis. 
// WMAX
//        Highest value for the tick mark labels written on the axis. 
// NDIV
//        Number of divisions. calculated according to the following convention: 
//        NDIV = N1 + 100*N2 + 10000*N3
//               where, 
//        N1
//               Number of primary divisions. 
//        N2
//               Number of second order divisions. 
//        N3
//               Number of third order divisions. 
//        Examples: 
//        NDIV=0
//               No tick marks. 
//        NDIV=2
//               produces 2 divisions with one tick mark in the middle of the axis. 
//        Note that, in case numeric labels are requested, N1 indicates the maximum number
//        of primary divisions. An appropriate algorithm calculates a
//        number of primary divisions less or equal to N1, in order to obtain ``reasonable''
//        labels. Option 'N' in CHOPT forces N1 to be used as the exact
//        number of primary divisions. 
// CHOPT
//        Character variable specifying the combinations of options desired. General options 
//        'G'
//               LoGarithmic scale, default is linear. 
//        'B'
//               Blank axis, i.e. the base line constituting the axis is not drawn. However
//               tick marks and labels are drawn. Useful when superimposing two
//               axes. 
//        'A'
//               An arrow is drawn at the end of the axis (position WMAX). 
//        'N'
//               N1 will be used as exact number of divisions. 
//        Orientation of the tick marks on the axis Tick marks are normally drawn on the
//        positive side of the axis. However, if the axis is vertical, i.e. if
//        X0=X1, then they are drawn on the ``negative'' side. Their orientation can
//        be selected by CHOPT. 
//        '+'
//               Tick marks are drawn on the positive side of the axis (default). 
//        '-'
//               Tick marks are drawn on the negative side of the axis. 
//        Specifying '+-' will draw tick marks on both sides of the axis. Orientation
//        of tick marks and labels in the working space Tick marks are normally
//        drawn orthogonal to the axis. However, in case of an oblique axis,
//        they can be drawn vertically. 
//        'V'
//               Tick marks are drawn Vertically (default is perpendicular to axis). 
//        Labeling an axis An axis is normally labeled, unless specified otherwise: 
//        'U'
//               Unlabeled axis (default is labeled). 
//        Position of labels on an axis Labels are normally drawn on the side opposite
//        to the tick marks, unless specified otherwise: 
//        '='
//               Labels are drawn on the same side as the tick marks. 
//        Orientation of labels on an axis. Labels are normally drawn parallel to
//        the axis. However if the axis is vertical, i.e. if X0=X1, then the labels are
//        drawn orthogonally. If the axis is horizontal, i.e. if Y0=Y1, then the labels
//        are Parallel to the axis: 
//        'P'
//               Labels are drawn Parallel to the axis 
//        'O'
//               Labels are drawn Orthogonal to the axis. 
//        Position of labels with respect to the tick marks. Labels are centered on
//        tick marks. However, if the axis is vertical (X0=X1), then they are right
//        adjusted. 
//        'R'
//               Labels are Right adjusted on a tick mark. 
//        'L'
//               Labels are Left adjusted on a tick mark. 
//        'C'
//               Labels are centered on tick a mark. (default) 
//        Direction of labels The default writing direction of labels is from left to right. 
//        'Y'
//               Writing direction is downwards. 
//        Format of labels Training blanks in the label strings are stripped, and then
//        the label is correctly aligned. If the last character of the string is a dot
//        '.', it is also stripped by default. 
//        '.'
//               The dot at the end of a string is mandatory. 
//        Type of labels Labels are by default numeric. 
//        'T'
//               The labels are alphanumeric text strings. In this case 12 default values
//               are provided, namely the 3-character abbreviations of the names of
//               the months: 'JAN', 'FEB', 'MAR',.... These values can be modified by
//               calling the routine IGLBL (see section [more info]). 
//        Optional grid An optional grid (cross-wires) can be drawn as a prolongation
//        of the primary tick marks. 
//        'W'
//               Draw cross-wires at the position of the primary tick marks. The length
//               of the grid can be defined, in world coordinates, with the IGSET
//               parameter AWLN. The current line type is used to draw the grid. 
//        Intrinsic parameters The default values for HIGZ intrinsic parameter settings
//        are shown below expressed as a percentage of the length of the axis
//        (world coordinates): 
//        Primary tick marks:
//               3.0 % 
//        Secondary tick marks:
//               1.5 % 
//        Third order tick marks:
//               .75 % 
//        Length of the arrow:
//               3.0 % 
//        Width of the arrow:
//               .75 % 
//        Characters height for labels:
//               2.0 % 
//        Characters spacing:
//               40% of the character height 
//        Labels offset:
//               4.0 % 
//        The size of the secondary tick marks is always 50% of the primary ones.
//        The size of the third order tick marks is always 50% of the secondary
//        ones. These values can be changed by calls to routine IGSET. The default
//        value is used unless the corresponding option is selected by CHOPT: 
//        'D'
//               The distance between the labels and the axis (the offset) is given
//               by the preceding call to IGSET with the parameter LAOF. 
//        'H'
//               The size (height) of the labels is given by the preceding call to
//               IGSET with the parameter LASI. 
//        'S'
//               The size of the tick marks is given by the preceding call to IGSET
//               with the parameter TMSI.

extern "C" void igaxis_(const float &X0,const float &X1,const float &Y0,const float &Y1,const float &WMIN,const float &WMAX,const int &NDIV,const char *CHOPT,int);
inline void IGAXIS(float X0,float X1,float Y0,float Y1,float WMIN,float WMAX,int NDIV,const char *CHOPT)
{igaxis_(X0,X1,Y0,Y1,WMIN,WMAX,NDIV,CHOPT,strlen(CHOPT));}

////////////////////////////////////////////////////////////////////////////////

// Text colour index.
// 
//                           +--------------------+
//                           |CALL  ISTXCI (ICOLI) |
//                           +--------------------+
//                                   
// 
// Action: This routine sets the text colour index attribute for use by future
// invocations of ITX. The routine IGSET (see section [more info]) can also be
// used with the parameter TXCI. Parameter Description: 
// 
// ICOLI
//        Text colour index. 

extern "C" void istxci_(const int &ICOLI);
inline void ISTXCI(int ICOLI)
{istxci_(ICOLI);}

////////////////////////////////////////////////////////////////////////////////

// Bidimensional matrix drawing
// 
//                  +--------------------------------------+
//                  | CALL  IGTABL (NX,NY,V,NPAR,PAR,CHOPT) |
//                  +--------------------------------------+
//                                   
// 
// Action: This routine draws a 2D matrix (i.e. table) according to the values of CHOPT
// and PAR. The PAR input parameter could be specified to change the
// aspect of the plot (see the description below). The position of the plot on the screen
// is given by the viewport of the current normalization transformation
// currently selected (the window is not used and could be anything). Parameter Description: 
// 
// NX
//        Number of cells in X. 
// NY
//        Number of cells in Y. 
// V(NX,NY)
//        Content of the cells. 
// NPAR
//        Number of parameters in PAR. 
// PAR(NPAR)
//        Array of real parameter. If PAR(i)=0. or NPAR a default value is taken. 
// CHOPT
//        CHARACTER variable specifying the options selected. The possible value
//        of CHOPT and the associate values of PAR are
//        describe below. The default value of CHOPT is 'P'. 
// 
//         Example of MATRIX drawing (see result on figure
//                   [more info] to
//                           [more info])
//                                   
// 
//       program matrix
//       call start_matrix('lego','LA')
// 
//       call start_matrix('lego1','L1A')
//       call start_matrix('lego2','L2')
//       call start_matrix('surf','SA')
//       call start_matrix('surf1','S1A')
//       call start_matrix('surf2','S2A')
//       call start_matrix('surf3','S3A')
//       call start_matrix('surf4','S4A')
//       call start_matrix('surfpol','SPOL')
//       call start_matrix('surfcyl','SCYL')
//       call start_matrix('surfsph','SSPH')
//       call start_matrix('surfpsd','SPSD')
//       end
//       subroutine start_matrix(name,chopt)
//       character*(*) name,chopt
//       parameter (nx=30,ny=30)
//       dimension v(nx,ny)
//       dimension par(29)
// *
// *              Parameters initialisation
// *
//       call vzero(par,29)
//       par(1)=30.
//       par(2)=23.
//       par(3)=-10.
//       par(4)=10.
//       par(5)=-10
//       par(6)=10.
//       par(9)=1030.
//       par(10)=1030.
//       par(11)=510.
//       par(12)=510.
//       par(13)=510.
//       par(14)=1.
//       par(15)=1.
//       par(16)=1.
//       par(20)=0.05
//       par(21)=-61.
//       par(22)=.1
//       par(23)=.1
//       par(24)=.15
//       par(25)=2.
//       par(26)=5.
//       par(27)=7.
//       par(28)=6.
//       par(29)=3.
// *
// *              Matrix filling
// *
//       x=-10.
//       y=-10.
//       s=20./float(nx)
//       do i=1,nx
// 
//          do j=1,ny
//             if(x.ne.0..and.y.ne.0)then
//                v(i,j)=100.*sin(x)/x*sin(y)/y
//             else
//                v(i,j)=100.
//             endif
//             x=x+s
//          enddo
//          y=y+s
//          x=-10.
//       enddo
// *
// *              Matrix drawing
// *
//       call start(NAME,9.,9.)
//       call isfais(0)
//       call igset('BORD',1.)
//       call igset('TXAL',32.)
//       call igset('CHHE',0.25)
//       call igtabl(nx,ny,v,29,par,chopt)
//       call igterm
//       call finish
//       end
// 
// 
// +-+--------------------------------------------------------------------------------------------+-+
// +-+-CHOPT-=-'P'+Polymarker-(scatter-plot)-------------------------------------------+----------+-+
// +-+-PAR-index--+----------------------------PAR-values------------------------------+-default--+-+
// | |     1      | Marker type see  ISMK.                                              |      1.  | |
// | |     2      | Maximum number of random points per cell                           |     50.  | |
// | |     3      | XMIN Lowest X-axis label                                           |   IXMIN  | |
// | |     4      | XMAX Highest Y-axis label                                          |   IXMAX  | |
// | |            |                                                                    |          | |
// | |     5      | YMIN Lowest Y-axis label                                           |   IYMIN  | |
// | |     6      | YMAX Highest Y-axis label                                          |   IYMAX  | |
// | |     7      | ZMIN Lowest Z value                                                |    ZMIN  | |
// | |     8      | ZMAX Highest Z value                                               |    ZMAX  | |
// | |     9      | 1000*IXMIN + IXMAX (Useful for ZOOM)                               |    1-NX  | |
// | |    10      | 1000*IYMIN + IYMAX (Useful for ZOOM)                               |    1-NY  | |
// +-+------------+--------------------------------------------------------------------+----------+-+
// 
//                                                  
//  
// 
// Figure:  Example of the  IGTABL Polymarker option
// 
//  [SCATTER]
// 
// 
// 
// +-+-CHOPT-=-'B'-Boxes--------------------------------------------------------------------------+-+
// +-+------------+--------------------------------------------------------------------+----------+-+
// +-+-PAR-index--+----------------------------PAR-values------------------------------+-default--+-+
// | |     1      | Not used                                                           |          | |
// | |     2      | Not used                                                           |          | |
// | |     3      | XMIN Lowest X-axis label                                           |   IXMIN  | |
// | |     4      | XMAX Highest Y-axis label                                          |   IXMAX  | |
// | |     5      | YMIN Lowest Y-axis label                                           |   IYMIN  | |
// | |            |                                                                    |          | |
// | |     6      | YMAX Highest Y-axis label                                          |   IYMAX  | |
// | |     7      | ZMIN Lowest Z value                                                |    ZMIN  | |
// | |     8      | ZMAX Highest Z value                                               |    ZMAX  | |
// | |     9      | 1000*IXMIN + IXMAX (Useful for ZOOM)                               |    1-NX  | |
// | |    10      | 1000*IYMIN + IYMAX (Useful for ZOOM)                               |    1-NY  | |
// +-+------------+--------------------------------------------------------------------+----------+-+
// 
//                                                  
//  
// 
// Figure:  Example of the  IGTABL Boxes option
// 
//  [BOXES]
// 
// 
// 
// +-+--------------------------------------------------------------------------------------------+-+
// +-+-CHOPT-=-'R'+aRrows--------------------------------------------------------------+----------+-+
// | | PAR index  |                            PAR values                              | default  | |
// +-+------------+--------------------------------------------------------------------+----------+-+
// | |     1      | Not used                                                           |          | |
// | |     2      | Not used                                                           |          | |
// | |     3      | XMIN Lowest X-axis label                                           |   IXMIN  | |
// | |     4      | XMAX Highest Y-axis label                                          |   IXMAX  | |
// | |     5      | YMIN Lowest Y-axis label                                           |   IYMIN  | |
// | |     6      | YMAX Highest Y-axis label                                          |   IYMAX  | |
// | |            |                                                                    |          | |
// | |     7      | ZMIN Lowest Z value                                                |    ZMIN  | |
// | |     8      | ZMAX Highest Z value                                               |    ZMAX  | |
// | |     9      | 1000*IXMIN + IXMAX (Useful for ZOOM)                               |    1-NX  | |
// +-+----10------+-1000*IYMIN-+-IYMAX-(Useful-for-ZOOM)-------------------------------+----1-NY--+-+
// 
//                                                  
//  
// 
// Figure:  Example of the  IGTABL aRrows option
// 
//  [ARROWS]
// 
// 
// 
// 
// +-+-CHOPT-=-'C'-Contour-plot-------------------------------------------------------------------+-+
// +-+------------+--------------------------------------------------------------------+----------+-+
// +-+-PAR-index--+----------------------------PAR-values------------------------------+-default--+-+
// | |     1      | Nlevel (min=2 max=50)                                              |     20.  | |
// | |     2      | 0 use colour to distinguish contours. Line type used is 1.         |      0.  | |
// | |            | 1.XXX use line style to distinguish contours. Colour index used    |          | |
// | |            | is XXX.                                                   |          | |
// | |            | 2.XXX line style and colour are the same for all contours.         |          | |
// | |            | Colour index used is XXX.                                 |          | |
// | |            |                                                                    |          | |
// | |     3      | XMIN Lowest X-axis label                                           |   IXMIN  | |
// | |     4      | XMAX Highest Y-axis label                                          |   IXMAX  | |
// | |     5      | YMIN Lowest Y-axis label                                           |   IYMIN  | |
// | |     6      | YMAX Highest Y-axis label                                          |   IYMAX  | |
// | |     7      | ZMIN Lowest Z value                                                |    ZMIN  | |
// | |     8      | ZMAX Highest Z value                                               |    ZMAX  | |
// | |            |                                                                    |          | |
// | |     9      | 1000*IXMIN + IXMAX (Useful for ZOOM)                               |    1-NX  | |
// +-+----10------+-1000*IYMIN-+-IYMAX-(Useful-for-ZOOM)-------------------------------+----1-NY--+-+
//                                                  
//  
// 
// Figure:  Example of the  IGTABL Contour option
// 
//  [CONTOUR]
// 
// 
// 
// +-+--------------------------------------------------------------------------------------------+-+
// +-+-CHOPT-=-'COL'-COLour-plot-------------------------------------------------------+----------+-+
// +-+-PAR-index--+----------------------------PAR-values------------------------------+-default--+-+
// | |     1      | 0 use the standard 8 colours                                       |      0.  | |
// | |     2      | ...                                                                |          | |
// | |     3      | XMIN Lowest X-axis label                                           |   IXMIN  | |
// | |     4      | XMAX Highest Y-axis label                                          |   IXMAX  | |
// | |            |                                                                    |          | |
// | |     5      | YMIN Lowest Y-axis label                                           |   IYMIN  | |
// | |     6      | YMAX Highest Y-axis label                                          |   IYMAX  | |
// | |     7      | ZMIN Lowest Z value                                                |    ZMIN  | |
// | |     8      | ZMAX Highest Z value                                               |    ZMAX  | |
// | |     9      | 1000*IXMIN + IXMAX (Useful for ZOOM)                               |    1-NX  | |
// | |    10      | 1000*IYMIN + IYMAX (Useful for ZOOM)                               |    1-NY  | |
// +-+------------+--------------------------------------------------------------------+----------+-+
//                                                  
//  
// 
// Figure:  Example of the  IGTABL COLour option
// 
//  [COLOUR]
// 
// 
// 
// 
// +-+--------------------------------------------------------------------------------------------+-+
// +-+-CHOPT-=-'T'+Text----------------------------------------------------------------+----------+-+
// +-+-PAR-index--+----------------------------PAR-values------------------------------+-default--+-+
// | |     1      | Text font                                                          |      1.  | |
// | |     2      | Text Precision                                                     |      0.  | |
// | |     3      | XMIN Lowest X-axis label                                           |   IXMIN  | |
// | |     4      | XMAX Highest Y-axis label                                          |   IXMAX  | |
// | |            |                                                                    |          | |
// | |     5      | YMIN Lowest Y-axis label                                           |   IYMIN  | |
// | |     6      | YMAX Highest Y-axis label                                          |   IYMAX  | |
// | |     7      | ZMIN Lowest Z value                                                |    ZMIN  | |
// | |     8      | ZMAX Highest Z value                                               |    ZMAX  | |
// | |     9      | 1000*IXMIN + IXMAX (Useful for ZOOM)                               |    1-NX  | |
// | |    10      | 1000*IYMIN + IYMAX (Useful for ZOOM)                               |    1-NY  | |
// +-+------------+--------------------------------------------------------------------+----------+-+
//                                                  
//  
// 
// Figure:  Example of the  IGTABL Text option
// 
//  [TABT]
// 
// 
// 
// 
// +-+--------------------------------------------------------------------------------------------+-+
// +-+-CHOPT-=-'K'+character-----------------------------------------------------------+----------+-+
// +-+-PAR-index--+----------------------------PAR-values------------------------------+-default--+-+
// | |     1      | Text font                                                          |      1.  | |
// | |     2      | Text Precision                                                     |      0.  | |
// | |     3      | XMIN Lowest X-axis label                                           |   IXMIN  | |
// | |     4      | XMAX Highest Y-axis label                                          |   IXMAX  | |
// | |            |                                                                    |          | |
// | |     5      | YMIN Lowest Y-axis label                                           |   IYMIN  | |
// | |     6      | YMAX Highest Y-axis label                                          |   IYMAX  | |
// | |     7      | ZMIN Lowest Z value                                                |    ZMIN  | |
// | |     8      | ZMAX Highest Z value                                               |    ZMAX  | |
// | |     9      | 1000*IXMIN + IXMAX (Useful for ZOOM)                               |    1-NX  | |
// | |    10      | 1000*IYMIN + IYMAX (Useful for ZOOM)                               |    1-NY  | |
// +-+------------+--------------------------------------------------------------------+----------+-+
//                                                  
//  
// 
// Figure:  Example of the  IGTABL character K option
// 
// 
//  [TABK]
// 
// 
// 
// 
// +-+-CHOPT-=-'L'-Lego-(mode-0)------------------------------------------------------------------+-+
// | |                                                                                            | |
// | | CHOPT = 'LB' Lego with BARO and BARW                                                       | |
// | | CHOPT = 'L1' Lego with colours (mode 1)                                                    | |
// +-+-CHOPT-=-'L2'-Lego-with-colours-(mode-2)----------------------------------------------------+-+
// | | CHOPT = 'S' Surface (mode 0)                                                               | |
// | | CHOPT = 'S1' Surface with colours (mode 1)                                                 | |
// | | CHOPT = 'S2' Surface with colours (mode 2)                                                 | |
// | |                                                                                            | |
// | | CHOPT = 'S3' Surface with contour plot on top (mode 3)                                     | |
// +-+-CHOPT-=-'S4'-Surface-with-Gouraud-shading-(mode-4)-----------------------------------------+-+
// | | CHOPT = 'CYL' Cylindrical for lego and surface                                             | |
// | | CHOPT = 'SPH' Spherical for lego and surface                                               | |
// | | CHOPT = 'PSD' Pseudo rapidity for lego and surface                                         | |
// +-+-PAR-index--+----------------------------PAR-values------------------------------+-default--+-+
// +-+------------+--------------------------------------------------------------------+----------+-+
// | |     1      | THETA                                                              |     30.  | |
// | |     2      | PHI                                                                |     30.  | |
// | |     3      | XMIN Lowest X-axis label                                           |   IXMIN  | |
// | |     4      | XMAX Highest Y-axis label                                          |   IXMAX  | |
// | |     5      | YMIN Lowest Y-axis label                                           |   IYMIN  | |
// | |     6      | YMAX Highest Y-axis label                                          |   IYMAX  | |
// | |            |                                                                    |          | |
// | |     7      | ZMIN Lowest Z value                                                |    ZMIN  | |
// | |     8      | ZMAX Highest Z value                                               |    ZMAX  | |
// | |     9      | 1000*IXMIN + IXMAX (Useful for ZOOM)                               |    1-NX  | |
// | |    10      | 1000*IYMIN + IYMAX (Useful for ZOOM)                               |    1-NY  | |
// | |    11      | NDVX                                                               |  510.00  | |
// | |    12      | NDVY                                                               |  510.00  | |
// | |            |                                                                    |          | |
// | |    13      | NDVZ                                                               |  510.00  | |
// | |    14      | XCOL                                                               |    1.00  | |
// | |    15      | YCOL                                                               |    1.00  | |
// | |    16      | ZCOL                                                               |    1.00  | |
// | |    17      | XTIC                                                               |    0.02  | |
// | |    18      | YTIC                                                               |    0.02  | |
// | |            |                                                                    |          | |
// | |    19      | ZTIC                                                               |    0.02  | |
// | |    20      | VSIZ                                                               |    0.02  | |
// | |    21      | VFON                                                               |    2.00  | |
// | |    22      | XVAL                                                               |    0.02  | |
// | |    23      | YVAL                                                               |    0.02  | |
// | |    24      | ZVAL                                                               |    0.04  | |
// | |            |                                                                    |          | |
// +-+----25------+-Palette------------------------------------------------------------+----0.04--+-+
// Table:  Values of the  IGTABL Lego and Surface option
// 
//  [tab-IGTABLS]
// 
// 
// 
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Lego option
// 
//  [LEGO]
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Lego L1 option
// 
//  [LEGO1]
// 
// 
// 
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Lego L2 option
// 
//  [LEGO2]
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Surface option
// 
//  [SURF]
// 
// 
// 
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Surface S1 option
// 
//  [SURF1]
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Surface S2 option
// 
//  [SURF2]
// 
// 
// 
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Surface S3 option
// 
//  [SURF3]
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Surface S4 option
// 
//  [SURF4]
// 
// 
// 
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Surface SPOL option
// 
//  [SPOL]
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Surface SCYL option
// 
//  [SCYL]
// 
// 
// 
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Surface SSPH option
// 
//  [SSPH]
// 
//                                                          
//  
// 
// Figure:  Example of the  IGTABL Surface SPSD option
// 
//  [SPSD]
// 
// 
// Note that the options POL, CYL, SPH, and PSD can be used together with
// any lego or surface options. 
// 
// 
// +-+-----------------+--------------------------------------------------+-##
// +-+-CHOPT--+-Description--------------------------------------+-##
// +-+------'H'--------+-Data-are-compacted-as-in-HPLOT.------------------+-##
// | |      'GX'       | loG on X coordinates. A log world                | ##
// +-+-----------------+-coordinates-must-be-defined-before.--------------+-##
// | |      'GY'       | loG on Y coordinates. A log world                | ##
// +-+-----------------+-coordinates-must-be-defined-before.--------------+-##
// | |      'GZ'       | loG on Z coordinates.                            | ##
// +-+------'A'--------+-2nd-vertical-axis-(legos-and-Surfaces-only)------+-##
// | |                 |                                                  | ##
// +-+-----------------+-axis-(for-the-2D-representations).---------------+-##
// +-+------'+'--------+-For-stacked-histograms-(legos).------------------+-##
// +-+------'Z'--------+-Allows-to-display-the-Z-scale.-------------------+-##
// Table:  Other options for  IGTABL
// 
//  [tab-IGTABL]
// 
// 
//     Example of stacked lego plots drawing (see result on figure
//                           [more info])
//                                   
// 
//       program stack
//       parameter (nx=10,ny=10)
//       parameter (npar=25)
//       dimension v1(nx,ny),v2(nx,ny),v3(nx,ny)
//       dimension par(npar)
//       call vzero(par,npar)
//       par(1)  = 30.
//       par(2)  = 23.
//       par(3)  = -10.
//       par(4)  = 10.
//       par(5)  = -10
//       par(6)  = 10.
//       par(9)  = 1000. + nx
//       par(10) = 1000. + ny
//       par(11) = 510.
// 
//       par(12) = 510.
//       par(13) = 510.
//       par(14) = 1.
//       par(15) = 1.
//       par(16) = 1.
//       par(20) = 0.05
//       par(21) = -61.
//       par(22) = .1
//       par(23) = .15
//       par(24) = .1
// *
// *              Matrices filling
// *
//       do i=1,nx
//          do j=1,ny
//             v1(i,j)=float(i)
//             v2(i,j)=float(i+j)
//             v3(i,j)=float(j)
//          enddo
//       enddo
// *
// *              Stack drawing
// *
//       call start('stack',9.,9.)
//       call igset('BARW',0.5)
//       par(25) = 2.
//       call igtabl(nx,ny,v1,npar,par,'+')
//       par(25) = 5.
//       call igtabl(nx,ny,v2,npar,par,'+')
//       par(25) = 3.
//       call igtabl(nx,ny,v3,npar,par,'LB1A')
//       call igterm
//       call finish
// *
//       end

extern "C" void igtabl_(const int &NX,const int &NY,const float* V,const int &NPAR,const float *PAR,const char *CHOPT,int);
inline void IGTABL(int NX,int NY,const float** V,int NPAR,const float *PAR,const char *CHOPT)
{igtabl_(NX,NY,(const float*)V,NPAR,PAR,CHOPT,strlen(CHOPT));}

////////////////////////////////////////////////////////////////////////////////

// Drawing a graph
// 
//                        +--------------------------+
//                        |CALL  IGRAPH (N,X,Y,CHOPT) |
//                        +--------------------------+
//                                   
// 
// Action: This routine draws (in the current normalization transformation) a graph
// with several possible presentations. Parameter Description: 
// 
// N
//        Number of components in the arrays X and Y. 
// X
//        Array of dimension N containing the x coordinates in world coordinates
//        space of the graph to be drawn. 
// Y
//        Array of dimension N containing the y coordinates in world coordinates
//        space of the graph to be drawn. 
// CHOPT
//        CHARACTER variable specifying the options chosen (multiple simultaneous
//        options are possible). 
//        'R'
//               The graph is Rotated, i.e. the values in X are used for the ordinate
//               and the values in Y for the abscissa (default is the contrary). 
//        'L'
//               All points are connected with a straight line. (default) 
//        'F'
//               A Fill area is drawn through the points with the current fill area
//               attributes. The border is never drawn unless the fill area interior style is
//               hollow or the routine IGSET has been called with 'BORD' and VAL = 1.. 
//        'C'
//               The values in Y are plotted in the form of a smooth curve. Spline approximation
//               algorithms are used. This option can be used with option F
//               in order to draw a smooth fill area. 
//        '*'
//               A star is plotted at every point. 
//        'P'
//               A marker is plotted at every point, according to the current polymarker attributes. 
//        'B'
//               The values in Y are plotted in the form of bars. The width of the bar is by
//               default 50% of the interval X(I)-X(I-1). This percentage can be
//               changed by calling IGSET with option BARW. 
//        'A'
//               X and Y axes are drawn on the border of the current normalization transformation. 
//        'GX'
//               Logarithmic scale on the X axis. 
//        'GY'
//               Logarithmic scale on the Y axis. 
// 
//         Example of GRAPH drawing (see result on figure
//                           [more info])
//                                   
// 
//       program  graph
//       character*4 chopt(4)
//       dimension x(9),y(9)
//       parameter (xsize=16.,ysize=20.)
//       data x/0.,.6,.3,.2,-.3,.3,-.2,-.3,-.6/
//       data y/0.,-.2,-.7,-.9,-.2,.2,.9,.7,.2/
//       data chopt/'AL*','AC*','AF*','ACF*'/
// *
//       call start('graph',xsize,ysize)
// *
// *              Viewports definition
// *
//       xnorm  = min(1.,xsize/ysize)
//       xnorm2 = xnorm/2.
//       ynorm  = min(1.,ysize/xsize)
//       ynorm2 = ynorm/2.
//       rmarg  = 0.05
//       rmarg2 = rmarg/2.
//       call isvp(10,rmarg,xnorm2-rmarg2,ynorm2+rmarg2,ynorm-rmarg)
//       call isvp(20,xnorm2+rmarg2,xnorm-rmarg,ynorm2+rmarg2,ynorm-rmarg)
//       call isvp(30,rmarg,xnorm2-rmarg2,rmarg,ynorm2-rmarg2)
//       call isvp(40,xnorm2+rmarg2,xnorm-rmarg,rmarg,ynorm2-rmarg2)
// *
// *              Some attributes setting
// *
//       call isclip(0)
//       call igset('FASI',244.)
//       call igset('BORD',1.)
//       call igset('CHHE',.05)
// *
// *              GRAPH drawing
// *
//       do i=1,4
//          call iswn(10*i,-1.,1.,-1.,1.)
//          call iselnt(10*i)
//          call igset('FAIS',0.)
//          call igbox(-1.,1.,-1.,1.)
//          call itx(.3,.9,'CHOPT = '''//CHOPT(I)//'''')
//          call igset('FAIS',3.)
//          call igraph(9,x,y,chopt(i))
//       enddo
//       call finish
// *
// 
//       end

extern "C" void igraph_(const int &N,const float *X,const float *Y,const char *CHOPT,int);
inline void IGRAPH(int N,const float *X,const float *Y,const char *CHOPT)
{igraph_(N,X,Y,CHOPT,strlen(CHOPT));}

////////////////////////////////////////////////////////////////////////////////

// Menus Input
// 
// +-------------------------------------------------------------------------------------------+
// |CALL  IGMENU (MN,CHTIT,*X1*,*X2*,*Y1*,*Y2*,NBU,CHUSER,N,CHITEM, CHDEF,CHVAL*,ICHOIC*,CHOPT) |
// +-------------------------------------------------------------------------------------------+
//                                   
// 
// Action: This routine displays a menu and returns the user's choice in the variable ICHOIC
// according to the option chosen. This routine works only on one
// menu: the menu management must be performed by the application program but this routine
// provides some facilities to manage several menus
// simultaneously. Parameter Description: 
// 
// MN
//        Menu number. To use segment capabilities of the workstation.
//        If MN=0 the segments are not used. 
// CHTIT
//        Menu title. 
// X1
//        X coordinate of lower left hand corner of menu box 
// Y1
//        Y coordinate of lower left hand corner of menu box 
// X2
//        X coordinate of upper right hand corner of menu box 
// Y2
//        Y coordinate of upper right hand corner of menu box 
// NBU
//        Number of User squares. 
// CHUSER
//        CHARACTER array of length NBU containing the text in the users' squares.
//        The last line of the menu is split into NBU boxes. 
// N
//        Number of items. 
// CHITEM
//        CHARACTER array of length N containing the text for the items. 
// CHDEF
//        CHARACTER array of length N containing the text for the parameters.
//        If CHOPT='P' the menu is split into two columns. The left column contains the
//        items and the right column the default value of the corresponding item.
//        CHDEF(I) (1 is a character string which contains the possible
//        values of the item number I: CHDEF(I)='value1, value2, value3,..., valueN'.
//        If CHDEF(I)=' ' there are no default values. 
// CHVAL*
//        CHARACTER array of length N into which parameter values are written.
//        If CHOPT='P' then CHVAL(I) contains the parameter
//        value for item I. 
// ICHOIC
//        Choice number. The description of the possible values returned in ICHOIC
//        is given in the following table: 
// 
//                   +-+---------+------------------------------------+-+
//                   +-+----0----+-Outside-of-the-menu----------------+-+
//                   +-+---100---+-Title-bar--------------------------+-+
//                   +-+--1,NBU--+-User-keys--------------------------+-+
//                   +-+---1000--+-Right-button-of-the-mouse-clicked--+-+
//                   | |   >0    | Item number                        | |
//                   +-+---------+------------------------------------+-+
//                                          
// 
// CHOPT
//        CHARACTER variable specifying the option(s) selected. 
// 
// The square at the left of the title bar moves and resizes the menu. The square at
// the right of the title bar moves the menu. 
// 
// 
// +-+------+--------------------------------------------------------------------+-+
// | | 'H'  | The picked item is highlighted. The last choice number must be     | |
// +-+------+-given-in-ICHOIC.---------------------------------------------------+-+
// +-+-'D'--+-Display-the-menu.--------------------------------------------------+-+
// +-+-'C'--+-Permit-a-choice-in-the-displayed-menu.-----------------------------+-+
// +-+-'E'--+-Erase-the-menu.----------------------------------------------------+-+
// | | 'P'  | The menu is a menu with parameters.                                | |
// +-+-'R'--+-Return-the-current-position-of-the-menu-in-X1,X2,Y1,Y2.---+-+
// +-+------+--------------------------------------------------------------------+-+
// +-+-'S'--+-Software-characters-are-used-to-draw-the-text-in-the-menu.---------+-+
// | | 'U'  | Update the user text in the user squares with the value in         | |
// | |      | CHUSER. The user square number is given in ICHOIC. The    | |
// | |      | options 'U' and 'H' are incompatible because they used both        | |
// +-+------+-ICHOIC-as-input-parameter.-----------------------------------------+-+
// +-+-'M'--+-Menu-drawn-on-a-Metafile.------------------------------------------+-+
// +-+-'Z'--+-Menu-stored-in-the-ZEBRA-picture.----------------------------------+-+
// | | 'N'  | The last input position is used to find the menu item. With this   | |
// | |      | option choices can be made in several menus at the same time       | |
// | |      | using a DO loop as shown below. NBMENU is the    | |
// | |      | number of menus on the screen.                                     | |
// +-+-'B'--+-A-rubberbanding-box-is-used-for-the-locator.-----------------------+-+
// +-+------+--------------------------------------------------------------------+-+
// | | 'T'  | The title bar is not drawn, then the menu can not be moved         | |
// +-+------+-interactively.-----------------------------------------------------+-+
// +-+-'W'--+-The-menu-is-drawn-with-Width.--------------------------------------+-+
// +-+-'A'--+-The-menu-is-drawn-with-shAdow.-------------------------------------+-+
// +-+-'V'--+-Draw-only-the-vertical-part-of-width-or-shadow.--------------------+-+
// | | 'O'  | Like option 'V' but the width or shadow is aligned on the menu     | |
// | |      | frame.                                                             | |
// +-+-'I'--+-Input-menu.-A-parameter-menu-is-displayed-and- IGMENU-is-entered----+-+
// | |      | directly in request string. This is useful to perform a request    | |
// | |      | string without a very complicated initialization part.             | |
// +-+------+--------------------------------------------------------------------+-+
// +-+-'K'--+-Key-menu.-The-user-keys-are-drawn-as-key.--------------------------+-+
// Table:  Options for  IGMENU
// 
//  [tab-IGMENU]
// 
// 
// Example
// 
// This example program shows how IGMENU can manage several menus at the same time. 
// 
//                     How to manage several menus
//                                   
// 
//       PROGRAM MENU
// *
//       COMMON /PAWC/H(50000)
//       PARAMETER (NBMENU=3)
//       CHARACTER*10 CHU, CHI, CHD, CHV, CHTIT, CHOPT
//       CHARACTER*80 TEXT
//       CHARACTER*16 CHLOC(3)
//       DIMENSION CHU(3),NBU(NBMENU),NBI(NBMENU)
//       DIMENSION CHI(3),CHD(3),CHV(3),CHTIT(NBMENU)
//       DIMENSION X1(NBMENU),X2(NBMENU),Y1(NBMENU),Y2(NBMENU)
// *     Last choice in the menu NB i (useful for HIghligth)
//       DIMENSION ICCH(NBMENU)
//       DATA CHU /'Quit','Exit','GED'/
//       DATA CHI /'Choice 1', '|Choice 2', 'Choice 3'/
// *.______________________________________
// *
// *
// *       Initialize HIGZ
// *
//       CALL MZEBRA(-3)
//       CALL MZPAW(50000,' ')
//       CALL IGINIT(0)
//       CALL IGWKTY(KWKTYP)
//       CALL IGSSE(6,KWKTYP)
//       CALL ISELNT(0)
//       CALL MESSAGE('Example of the IGMENU usage in multiple input')
// *
// *       Initialize and display menu number 1
// *
//   1   ICCH(1)=0
//       X1(1)=0.14
//       X2(1)=0.35
//       Y1(1)=0.1
//       Y2(1)=0.25
//       NBU(1)=2
//       NBI(1)=3
//       CHTIT(1)='MENU 1'
//       CALL IGMENU (0,CHTIT(1),X1(1),X2(1),Y1(1),Y2(1),NBU(1),CHU,
//      +             NBI(1),CHI,CHD,CHV,ICH,'S   D')
// *
// *       Initialize and display menu number 2
// *
// 
//       ICCH(2)=0
//       X1(2)=0.3
//       X2(2)=0.56
//       Y1(2)=0.3
//       Y2(2)=0.45
//       NBU(2)=2
//       NBI(2)=3
//       CHTIT(2)='MENU 2'
//       CALL IGMENU (0,CHTIT(2),X1(2),X2(2),Y1(2),Y2(2),NBU(2),CHU,
//      +             NBI(2),CHI,CHD,CHV,ICH,'S   D')
// *
// *       Initialize and display menu number 3
// *
//       ICCH(3)=0
//       X1(3)=0.05
//       X2(3)=0.95
//       NBU(3)=3
//       NBI(3)=0
//       CHTIT(3)='MENU 3'
//       Y1(3)=0.9
//       Y2(3)=0.935
//       CALL IGMENU (0,CHTIT(1),X1(3),X2(3),Y1(3),Y2(3),NBU(3),CHU,
//      +             NBI(3),CHI,CHD,CHV,ICH,'ST  D')
// *
// *       Initialize the current menu number
// *
//       IMENU=3
// *
// *       Request in the current menu
// *
//    10 CONTINUE
//       IF(IMENU.LT.3)THEN
//          CHOPT='S   CR'
//       ELSE
//          CHOPT='ST  C'
//       ENDIF
//       ICH=ICCH(IMENU)
//       CALL IGMENU (0,CHTIT(IMENU),X1(IMENU),X2(IMENU),
//      +             Y1(IMENU),Y2(IMENU),NBU(IMENU),CHU,
//      +             NBI(IMENU),CHI,CHD,CHV,ICH,CHOPT)
// *
// *       If the choice is outside the menu (ICH=0), we search here
// *       if the input is in an other menu (CHOPT='N')
// *
//       IF(ICH.EQ.0)THEN
//          DO 20  I=1,NBMENU
//             IF(I.LT.3)THEN
//                CHOPT='S CRN'
//             ELSE
//                CHOPT='SCTNKU'
//             ENDIF
//             ICH=ICCH(I)
//             CALL IGMENU (0,CHTIT(I),X1(I),X2(I),Y1(I),Y2(I),
// 
//      +                   NBU(I),CHU,
//      +                   NBI(I),CHI,CHD,CHV,ICH,CHOPT)
//             IF(ICH.NE.0)THEN
//                IMENU=I
//                GOTO 30
//             ENDIF
//    20    CONTINUE
// *
// *       After the DO loop the input is outside all menus
// *
//          CALL MESSAGE('Outside the menus')
//          GOTO 10
//       ENDIF
//       ICCH(IMENU)=ICH
// *
// *       Analyses the result
// *
//    30 CONTINUE
//       IF(ICH.GT.0)THEN
//          WRITE(TEXT,'(''Menu : '',I1,'', choice : '',I1)')IMENU,ICH
//          CALL MESSAGE(TEXT)
//          GOTO 10
//       ENDIF
//       IF(ICH.EQ.-100)THEN
//          WRITE(TEXT,'(''Menu : '',I1,'', title bar'')')IMENU
//          CALL MESSAGE(TEXT)
//          GOTO 10
//       ENDIF
//       IF(ICH.EQ.-1000)THEN
//          CALL MESSAGE('Right button of the mouse')
//          GOTO 10
//       ENDIF
//       IF(ICH.EQ.-1)THEN
//          WRITE(TEXT,'(''QUIT from menu : '',I1)')IMENU
//          CALL MESSAGE(TEXT)
//          CALL IGEND
//          GOTO 999
//       ENDIF
//       IF(ICH.EQ.-2)THEN
//          WRITE(TEXT,'(''EXIT from menu : '',I1)')IMENU
//          CALL MESSAGE(TEXT)
//          CALL IGEND
//          GOTO 999
//       ENDIF
//       IF(ICH.EQ.-3)THEN
//          CALL MESSAGE('Invoke the Graphics Editor')
//          CALL IZPICT('*','S')
//          CALL IZPICT('P1','M')
//          CALL IGRNG(20.,20.)
//          CALL IZGED('P1','S')
//          GOTO 1
//       ENDIF
//       IF(ICH.LT.0)THEN
// 
//         WRITE(TEXT,'(''Menu : '',I1,'', choice : '',I2)')IMENU,ICH
//          CALL MESSAGE(TEXT)
//          GOTO 10
//       ENDIF
// *
//   999 END
//       SUBROUTINE MESSAGE(TEXT)
//       CHARACTER*(*) TEXT
//       CALL IGZSET('G')
//       CALL ISELNT(0)
//       CALL IGSET('FACI',0.)
//       CALL IGSET('FAIS',1.)
//       CALL IGSET('BORD',1.)
//       CALL IGBOX(0.,1.,0.,0.04)
//       CALL IGSET('TXAL',23.)
//       CALL IGSET('CHHE',0.02)
//       CALL IGSET('TXFP',-100.)
//       CALL ITX(0.5,0.02,TEXT)
//       call iuwk(0,0)
//       END

// Due to complications in passing vector of characters I use here cfortran subroutine
PROTOCCALLSFSUB14(IGMENU,igmenu,INT,STRING,PFLOAT,PFLOAT,PFLOAT,PFLOAT,INT,STRINGV,INT,STRINGV,STRINGV,PSTRINGV,PINT,STRING)
#define IGMENU(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14)  CCALLSFSUB14(IGMENU,igmenu,INT,STRING,PFLOAT,PFLOAT,PFLOAT,PFLOAT,INT,STRINGV,INT,STRINGV,STRINGV,PSTRINGV,PINT,STRING,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14)
//void inline IGMENU(MN,CHTIT,X1,X2,Y1,Y2,NBU,CHUSER,N,CHITEM,CHDEF,CHVAL,MN3,MN4)
//{CCALLSFSUB14(IGMENU,igmenu,INT,STRING,PFLOAT,PFLOAT,PFLOAT,PFLOAT,INT,STRINGV,INT,STRINGV,STRINGV,PSTRINGV,PINT,STRING,MN,CHTIT,X1,X2,Y1,Y2,NBU,CHUSER,N,CHITEM,CHDEF,CHVAL,MN3,MN4);}

////////////////////////////////////////////////////////////////////////////////

// Drawing software characters
// 
//                +------------------------------------------+
//                | CALL  IGTEXT (X,Y,CHARS,SIZE,ANGLE,CHOPT) |
//                +------------------------------------------+                      
// 
// Action: This routine draws a software character text, independently from
// underlying graphics package used by HIGZ. IGTEXT can produce over 300
// different graphic signs. The way in which software characters are defined
// is via a string of valid Fortran characters, intermixed by other valid Fortran
// characters, acting as ``escape'' characters (e.g. a change of alphabet, upper
// or lower case). The string is interpreted by IGTEXT and the resulting
// characters are defined according to the figure [more info], which shows the
// list of available software characters. This routine allows the user to mix
// different types of characters (roman, greek, special, upper and lower case,
// sub and superscript). There are a total of 10 control characters. Parameter
// Description: 
// 
// X
//        x coordinate in world coordinatesspace. 
// Y
//        y coordinate in world coordinatesspace. 
// CHARS
//        CHARACTER variable containing the text to be displayed. 
// SIZE
//        Size of the text in world coordinatesspace. 
// ANGLE
//        Inclination angle of the text inclination in degrees. 
// CHOPT
//        CHARACTER variable specifying the text alignment: 
//        'L'
//               The text is Left adjusted starting at the point (X,Y). 
//        'C'
//               The text is Centered around the the point (X,Y). 
//        'R'
//               The text is Right adjusted ending at the point (X,Y). 
//        'S'
//               The text size (length) is returned in ANGLE. 
//        Note that it is not possible to align vertically the text produce by IGTEXT.
//        The way to align vertically software text is to use ITX with the font 0
//        and precision 2 (see ISTXFP). 
// 
// [ESCCHAR] 
// 
// +--------------------------------------------------------------------------------------------------+
// +-+----+--------------------List-of-escape-characters-and-their-meaning--------------------------+-+
// +-+-<--+-go-to-lower-case-----------------------+-+->--+-go-to-upper-case-(default)--------------+-+
// +-+-[--+-go-to-greek-(Roman-=-default)----------+-+-]--+-end-of-greek----------------------------+-+
// +-+-"--+-go-to-special-symbols------------------+-+-#--+-end-of-special-symbols------------------+-+
// | | "  | go to superscript                      | | ?  | go to subscript                         | |
// +-+-!--+-go-to-normal-level-of-script-----------+-+-&--+-backspace-one-character-----------------+-+
// +-+----+----------------------------------------+-+----+-----------------------------------------+-+
// +-+-$--+-termination-character-(optional)-------+-+----+-----------------------------------------+-+
// 
// Note that characters can be also entered directly in lower case or upper case
// instead of using the control characters < and >. The boldface characters may
// be simulated by setting the attributes 'PASS' and 'CSHI' with IGSET. The meaning
// of these attributes is the following: Every stroke used to display the
// character is repeated PASS times, at a distance (in percentage of the character
// height) given by CSHI.

extern "C" void igtext_(const float &X,const float &Y,const char *CHARS,const float &SIZE,const float &ANGLE,const char *CHOPT,int,int);
inline void IGTEXT(float X,float Y,const char *CHARS,float SIZE,float ANGLE,const char *CHOPT)
{igtext_(X,Y,CHARS,SIZE,ANGLE,CHOPT,strlen(CHARS),strlen(CHOPT));}

////////////////////////////////////////////////////////////////////////////////

// Text
// 
//                          +----------------------+
//                          | CALL  ITX (X,Y,CHARS) |
//                          +----------------------+
// 
//                                   
// 
// Action: This routine draws a text string on the currently active workstations
// (there must be at least one). The appearance of the text is controlled by
// attributes set by the current ``text colour index'' (see routine ISTXCI section ),
// the current ``character height'' (see routine ISCHH section ), the current
// ``text orientation'' (see routine ISCHUP section and the option TANG of the
// routine IGSET section ), the current ``text alignment'' (see routine ISTXAL
// section [more info]) and the current ``text font and precision'' (see routine
// ISTXFP section [more info]). Parameter Description: 
// 
// X
//        X coordinate in WC space. 
// Y
//        Y coordinate in WC space. 
// CHARS
//        CHARACTER variable containing the text to be displayed. Only the
//        following characters are allowed to appear in CHARS: 
// 
//             !"#$%&'()*+,-./0123456789:;<=>?
//             @ABCDEFGHIJKLMNOPQRSTUVWXYZ[-]^_
//             abcdefghijklmnopqrstuvwxyz{|}~
//             and the space.
// 
// Software characters (i.e. drawn with lines and not provided by the hardware)
// can be produced with routine IGTEXT. 

extern "C" void itx_(const float &X,const float &Y,const char *CHARS,int);
inline void ITX (float X,float Y,const char *CHARS)
{itx_(X,Y,CHARS,strlen(CHARS));};

////////////////////////////////////////////////////////////////////////////////

// Graphic package control
// 
//                        +--------------------------+
//                        |CALL  IGSSE (IERRF,KWTYPE) |
//                        +--------------------------+
//                                   
// 
// Action: In general, the initialization of the underlaying graphics package
// consists in several calls to different routines, in order to set the environment
// parameters. For user's convenience and for most applications, IGSSE initializes
// the standard graphic package environment. In particular, the default
// primitives attributes and the default window, viewport, workstation window and
// workstation viewport are initialized. Sophisticated applications may need
// to call the specialized basic control routines, namely IOPKS, IOPWK, IACWK,
// ISWKWN and ISWKVP, instead of using IGSSE. IGSSE opens only a
// single workstation. 
// 
// Parameter Description: 
// 
// IERRF
//        Error file logical unit number. 
// KWTYPE
//        Workstation type. See the description of IOPWK section [more info]. 
// 
// IGSSE calls the following routines: 
// 
// IOPKS
//        See section [more info]. 
// IOPWK(1,KONID,KWTYPE)
//        See section [more info]. 
// IACWK(1)
//        See section [more info]. 
// 
// Note that KONID is initialized in IGSSE depending on the underlying graphics
// package used. In general KONID is set to 1. In addition, the workstation
// window and viewport are also initialized in IGSSE as follows: 
// 
//      CALL ISWKWN(1,0.,1.,0.,1.)
//      CALL ISWKVP (1,0.,XMAX,0.,YMAX)
// 
// where XMAX and YMAX are the screen dimensions in pixels. In addition, the
// following primitives attributes (see details below) are initialized: 
// 
// 
//         +-+-Attributes-names----------+-Default-values---------+-+
//         +-+---------------------------+------------------------+-+
//         | | Polyline colour index     | 1                      | |
//         | | Line type                 | 1                      | |
//         | | Line width                | 1.0                    | |
//         | | Polymarker colour index   | 1                      | |
//         | | Marker type               | 1                      | |
//         | | Marker scale factor       | 1.0                    | |
//         | |                           |                        | |
//         | | Fill area colour index    | 1                      | |
//         | | Fill area interior style  | 0                      | |
//         | | Fill area style index     | 1                      | |
//         | | Character height          | 0.01                   | |
//         | | Character up vector       | 0.0,1.0                | |
//         | | Text alignment            | 0,0                    | |
//         | |                           |                        | |
//         | | Text font and precision   | 0,2                    | |
//         | | Text colour index         | 1                      | |
//         | | Clipping indicator        | 1                      | |
//         +-+-GKS-Aspect-source-flag----+-Individual-attributes--+-+
//                                   
// 
// In addition to this initialization role, IGSSE, when it is used in the context
// of the Telnetg program, allows to open the connection between the remote
// machine and the local one even if the X Window System is not available. This
// is done by giving to IGSSE the negative value of the local workstation type.

extern "C" void igsse_(const int &IERRF,const int &KWTYPE);
inline void IGSSE (int IERRF,int KWTYPE)
{igsse_(IERRF,KWTYPE);}

////////////////////////////////////////////////////////////////////////////////

// Get workstation type
// 
//                          +----------------------+
//                          |CALL  IGWKTY (KWTYPE*) |
//                          +----------------------+
//                                   
// 
// Action: This routine gets the workstation type from the standard input.
// Parameter Description: 
// 
// KWTYPE
//        Workstation type. A call to this routine will prompt the user with: 
// 
//        Workstation type (?=HELP) =1
// 
//        Just typing CR will return the default value in KWTYPE. The value of the
//        default depends on the HIGZ installation. Typing ? will give a short help
//        listing on all the different possible workstation types. Any other answer
//        will be interpreted as a new workstation type. Note that with the X11
//        version of HIGZ the routine IGWKTY will accept a workstation type like:
//        n.hostname where n is the line number in the file higzwindows.dat and
//        hostname is the name of the machine on which the graphics will be displayed.
//        In this way it is not necessary to define the variable DISPLAY
//        before using HIGZ. 
//            1.If a workstation type like n.hostname is entered, the hostname is
//            written at the end of the line n in higzwindows.dat. 
//            2.If the workstation type n is entered and if a hostname is present
//            on the line n in higzwindows.dat, the graphics will be redirected to the
//               machine hostname. 
//            3.If the workstation type n is entered and if a hostname is not on the
//            line n in higzwindows.dat, the graphics will be redirected to the machine
//               defined by the variable DISPLAY. 
//            4.If the workstation type n. is entered and if a hostname is present on
//            the line n in higzwindows.dat, the graphics will be redirected to the
//               machine defined by the variable DISPLAY and hostname is removed from
//               the line n in higzwindows.dat. 
// 
// Remark: In the file higzwindows.dat, it is possible to specify the name of the
// window just after the hostname. 

extern "C" void igwkty_(int &KWTYPE);
inline void IGWKTY (int &KWTYPE)
{igwkty_(KWTYPE);}

////////////////////////////////////////////////////////////////////////////////

// Polyline
// 
//                            +------------------+
//                            | CALL  IPL (N,X,Y) |
//                            +------------------+
//                                   
// 
// Action: This routine draws a polyline on the currently active workstations (there must
// be at least one). The polyline connects N points (N>=2) by means
// of N-1 line segments. The X and Y coordinates of the points are in two N-dimensional
// arrays. The appearance of a polyline is controlled by the current
// ``polyline colour index'' (see routine ISPLCI section ), the current ``line type''
// (see routine ISLN section ) and the current ``line width'' (see routine
// ISLWSC section [more info]). Parameter Description: 
// 
// N
//        Number of points. 
// X
//        Array of dimension N containing the x coordinates in WC space. 
// Y
//        Array of dimension N containing the y coordinates in WC space. 

extern "C" void ipl_(const int &N,const float *X,const float *Y);
inline void IPL(int N,const float *X,const float *Y)
{ipl_(N,X,Y);}

////////////////////////////////////////////////////////////////////////////////

// Polymarker
// 
//                            +------------------+
//                            | CALL  IPM (N,X,Y) |
//                            +------------------+
//                                   
// 
// Action: This routine draws a polymarker on the currently active workstations
// (there must be at least one). Markers are placed at N points (N>=1), whose
// x and y coordinates are given in two N-dimensional arrays. The appearance of
// a polymarker is controlled by the current ``polymarker colour index'' (see
// routine ISPMCI section ), the current ``marker type'' (see routine ISMK section
// [more info]) and the current ``marker scale factor'' (see routine
// ISMKSC section [more info]). Parameter Description: 
// 
// N
//        Number of points. 
// X
//        Array of dimension N containing the x coordinates in WC space. 
// Y
//        Array of dimension N containing the y coordinates in WC space.

extern "C" void ipm_(const int &N,const float *X,const float *Y);
inline void IPM(int N,const float *X,const float *Y)
{ipm_(N,X,Y);}

////////////////////////////////////////////////////////////////////////////////

// Character height
// 
//                            +------------------+
//                            | CALL  ISCHH (CHH) |
//                            +------------------+
//                                   
// 
// Action: This routine sets the character height attribute for use by future
// invocations of ITX. The routine IGSET (see section [more info]) can also be used
// with the parameter CHHE. Parameter Description: 
// 
// CHH
//        Character height. The default set by IGSSE is 0.01. The height is given
//        in world coordinates and it must be positive. 

extern "C" void ischh_(const float &);
inline void ISCHH(float CHH)
{ischh_(CHH);}

////////////////////////////////////////////////////////////////////////////////

#endif
#endif
