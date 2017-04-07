// $Id: CsHbookProto.h,v 1.8 2010/06/11 08:55:25 tnagel Exp $
 
/*!
   \file    CsHbookProto.h
   \brief   Prototypes for Hbook.
   \author  Alexander Zvyagin
   \version $Revision: 1.8 $
   \date    $Date: 2010/06/11 08:55:25 $
*/

// http://wwwinfo.cern.ch/asdoc/hbook_html3/hboomain.html

#ifndef CsHbookProto_h
#define CsHbookProto_h

#include <string.h>   // for strlen()

extern "C"
{
  void   hbname_(int const &id, char const *chblok, const void *address, char const *var, int,int);
  void   hbnt_  (int const &id, char const *title, char const *option,int,int);
  void   hbook1_(int const &id, char const *title, int const &nx, const float &xmi, const float & xma, const float &vmx,int);
  void   hbook2_(int const &id, char const *title, int const &nx, const float &xmi, const float & xma,
                                                   int const &ny, const float &ymi, const float & yma, const float &vmx,int);
  void   hcdir_ (char const *const directory,char const *const option, int,int);
  void   hcopy_ (int const &id1, int const &id2, char const *title, int);
  void   hf1_   (int const &id, const float &x,                 const float &w);
  void   hf2_   (int const &id, const float &x, const float &y, const float &w);
  void   hfill_ (int const &id, const float &x, const float &y, const float &w);
  void   hfith_ (int const &id, float (*f)(const float *x), char const *options, int const &n_par,
                 float x[], const float step[], const float x_min[], const float x_max[],
                 float first_deviations[], float &chi_sqr, int );
  void   hfithn_(int const &id, char const *f_name, char const *options, int const &n_par,
                 float x[], const float step[], const float x_min[], const float x_max[],
                 float first_deviations[], float &chi_sqr, int,int );
  void   hfnt_  (int const &id);
  void   hgnt_  (int const &id,int const &nevent,int &error);
  void   hldir_ (char const *directory, char const *option, int,int);
  void   hlimit_(int const &size);
  void   hmdir_ (char const *directory, char const *option, int, int);
  void   hopera_(int const &id1,char const *oper,int const &id2,int const &id3,const float &c1, const float &c2, int);
  void   hpdir_ (char const *directory, char const *option, int, int);
  void   hrend_ (char const *top, int);
  void   hreset_(int const &id,char const *title,int);
  void   hrget_ (int const &id, char const *file_name, char const *options, int,int);
  void   hrin_  (int const &id, int const &icycle, int const &iofset);
  void   hropen_(int const &lun, char const *top, char const *file, char const *option, int const &reclen, int *error, int,int,int);
  void   hrout_ (int const &id, int &cycle, char const *option, int);
  void   hrput_ (int const &id, char const *file_name, char const *options, int, int);
  void   hgive_ (int const &id, char *title, int &nx, float &xmin, float &xax, int &ny, float &ymin, float &ymax,int &nwt, char * &address, int);
  void   hindex_(void);
  float  hi_    (int const &id, int const &i);
  float  hij_   (int const &id, int const &i, int const &j);
  float  hx_    (int const &id, const float &x);
  float  hxy_   (int const &id, const float &x, const float &y);
  float  hxyij_ (int const &id, const float &x, const float &y, int *i, int *j);
  float  hie_   (int const &id, const int &i);
  float  hije_  (int const &id, const int &i, const int &j);
  float  hstati_(int const &id,int const &w,char const *opt,int const &n, int);
  int    hexist_(int const &id);
  void   hdelet_(int const &id);
  void   hnoent_(int const &id, int &n_events);
  void   hlnext_(int &id,char *type,char *title,const char *opt, int,int,int);
};

//******************************************************************************

inline void hbname(int id,const char *chblok,const void *address,char const *var)
{ hbname_(id,chblok,address,var,strlen(chblok),strlen(var)); }

////////////////////////////////////////////////////////////////////////////////
// 
// CALL HBNT (ID,CHTITL,CHOPT)
// 
// Action: Books a CWN. 
// 
// Input parameters: 
// ID Identifier of the Ntuple. 
// CHTITL 
//      Character variable specifying the title associated to the Ntuple. 
// CHOPT 
//      Character variable specifying the desired options. 
//      ' ' for disk resident Ntuples (default). 
//      'D' idem as ' '. 
//      'M' for memory resident Ntuples. 
// 
// The CWN will be stored in the current HBOOK directory. The variables to be stored in the Ntuple will be specified with routine HBNAME or HBNAMC described below. 
// 
// When the CWN will be filled with HFNT, the memory buffers associated with each column will be written to the file and directory corresponding to the current working
// directory when HBNT was called. Remember that when routine HROPEN is called, the current working directory is automatically set to the top directory of that file. It is
// therefore convenient to call HBNT immediately after HROPEN. If this was not the case, routine HCDIR must be called prior to HBNT to set the current working directory.
// When the Ntuple has been filled (via calls to HFNT) the resident buffers in memory as well as the Ntuple header must be written to the file with a call to HROUT. Before calling
// HROUT, the current working directory must be set to the current directory when HBNT was called. 
// 
////////////////////////////////////////////////////////////////////////////////

inline void hbnt(int id,char const *title,char const *option)
{ hbnt_(id,title,option,strlen(title),strlen(option)); }

////////////////////////////////////////////////////////////////////////////////

inline void hbook1(int id,char const *title,int nx,float xmi,float xma,float vmx=0)
{ hbook1_(id,title,nx,xmi,xma,vmx,strlen(title)); }

////////////////////////////////////////////////////////////////////////////////

inline void hbook2(int id,char const *title,int nx,float xmi,float xma,int ny,float ymi,float yma,float vmx=0)
{ hbook2_(id,title,nx,xmi,xma,ny,ymi,yma,vmx,strlen(title)); }

////////////////////////////////////////////////////////////////////////////////

inline void hcdir(char const *const directory,char const *const option)
{
  size_t n=strlen(directory);
  hcdir_(directory,option,n,strlen(option));  // Change directory.

  if( strlen(option)>0 && (*option=='r' || *option=='R') )
    for( int i=n-1; i>=0 && directory[i]==' '; i-- )
      const_cast<char*>(directory)[i]=0;
}

////////////////////////////////////////////////////////////////////////////////

inline void hf1(int id,float x,float w=1)
{ hf1_(id,x,w); }

////////////////////////////////////////////////////////////////////////////////

inline void hf2(int id,float x,float y,float w=1)
{ hf2_(id,x,y,w); }

////////////////////////////////////////////////////////////////////////////////

inline void hfill(int id,float x,float y,float w=1)
{ hfill_(id,x,y,w); }

////////////////////////////////////////////////////////////////////////////////

inline void hfith(int id, float (*f)(const float *x), char const *options, int n_par,
                   float x[], const float step[], const float x_min[], const float x_max[],
                   float first_deviations[], float &chi_sqr)
{ hfith_(id,f,options,n_par,x,step,x_min,x_max,first_deviations,chi_sqr,strlen(options)); }

////////////////////////////////////////////////////////////////////////////////
// 
// CALL HFITHN (ID,CHFUN,CHOPT,NP,*PAR*,STEP,PMIN,PMAX,SIGPAR*,CHI2*)
// 
// Action: Fits the given special function to the contents of a one-dimensional histogram, and optionally superimposes it to the histogram when editing. 
// 
// Input parameters: 
// ID Histogram identifier. 
// CHFUN 
//      Character variable specifying the desired parametric function. Possible keywords are: 
//      G to fit gaussian PAR(1)*exp(-0.5*((x-PAR(2))/PAR(3))**2).
//          The first parameter PAR(1) corresponds to the normalization, the second parameter PAR(2) corresponds to the mean value, while the third parameter
//          PAR(3) corresponds to the half width of the gaussian.     
//      E to fit exponential exp(PAR(1)+PAR(2)*x)     
//      Pn to fit polynomyal PAR(1)+PAR(2)*x+PAR(3)*x**2......+PAR(n+1)*x**n     
//      Any combination of these keywords with the 2 operators + or * is allowed, e.g. 'p4+g', a combination of a 4th degree polynomial and a gaussian, which needs eight
//      parameters or p2*g+g, a second degree polynomyal and 2 gaussians, needing 9 parameters. The order of the parameters in PAR must correspond to the order of the
//      basic functions. For example, in the first case above, PAR(1:5) apply to the polynomial of degree 4 and PAR(6:8) to the gaussian while in the second case
//      PAR(1:3) apply to the polynomial of degree 2, PAR(4:6) to the first gaussian and PAR(7:9) to the second gaussian. Blanks are not allowed in the expression. 
// CHOPT 
//      'B' Some or all parameters are bounded. The arrays STEP, PMIN and PMAX must be specified. By default all parameters vary freely. 
//      'D' The user is assumed to compute derivatives analytically using the routine HDERIV. By default, derivatives are computed numerically. 
//      'E' Perform a detailed error analysis using the MINUIT routines HESSE and MINOS 
//      'F' Force storing of the result of the fit bin by bin with the histogram. 
//      'L' Use the logaritmic Likelihood fitting method. By default the chisquared method is used. 
//      'M' Invoke interactive MINUIT. 
//      'N' The results of the fit are not stored bin by bin with the histogram. By default the function is calculated at the centre of each bin in the specified range. 
//      'Q' Quiet mode. No printing. 
//      'R' Fit a Restricted area of the 1-D histogram. IFTLOW = IQUEST(11) specifies the lower limit of the minimization domain,
//          IFTUP = IQUEST(12) specifies the upper limit of the minimization domain.   
//      'V' Verbose mode. Results are printed after each iteration. By default only final results are printed. 
//      'W' Set event weights to one. By default weights are taken according to statistical errors. 
// NP Number of parameters. 
// PAR Array of dimension NP with initial values for the parameters. 
// STEP Array of dimension NP with initial step sizes for the parameters ('B' option only). 
// PMIN 
//      Array of dimension NP with the lower bounds for the parameters ('B' option only). 
// PMAX 
//      Array of dimension NP with the upper bounds for the parameters ('B' option only). 
// 
// Output parameters: 
// PAR Array of dimension NP with the final fitted values of the parameters. 
// SIGPAR 
//      Array of dimension NP with the standard deviations on the final fitted values of the parameters. 
// CHI2 Chisquared of the fit. 
// 
////////////////////////////////////////////////////////////////////////////////

inline void hfithn(int id, char const *f_name, char const *options, int n_par,
                   float x[], const float step[], const float x_min[], const float x_max[],
                   float first_deviations[], float &chi_sqr)
{
  hfithn_(id,f_name,options,n_par,x,step,x_min,x_max,first_deviations,chi_sqr,
          strlen(f_name),strlen(options));
}

////////////////////////////////////////////////////////////////////////////////

inline void hfnt(int id)
{ hfnt_(id); }

////////////////////////////////////////////////////////////////////////////////

inline void hgnt(int id,int n,int &error)
{ hgnt_(id,n,error); }

////////////////////////////////////////////////////////////////////////////////

inline void hldir(char const *directory,char const *option)
{ hldir_(directory,option,strlen(directory),strlen(option)); }

////////////////////////////////////////////////////////////////////////////////

inline void hlimit(int size)
{ hlimit_(size); }

////////////////////////////////////////////////////////////////////////////////

inline void hmdir(const char *directory,const char *option)
{ hmdir_(directory,option, strlen(directory),strlen(option)); }

////////////////////////////////////////////////////////////////////////////////

inline void hpdir(const char *directory,const char *option)
{ hpdir_(directory,option, strlen(directory),strlen(option)); }

////////////////////////////////////////////////////////////////////////////////

inline void hrend(const char *directory)
{ hrend_(directory, strlen(directory)); }

////////////////////////////////////////////////////////////////////////////////

inline void hreset(int id,const char *title=" ")
{ hreset_(id,title,strlen(title)); }

////////////////////////////////////////////////////////////////////////////////

inline void hrget(int id,const char *name,const char *option)
{ hrget_(id,name,option,strlen(name),strlen(option)); }

////////////////////////////////////////////////////////////////////////////////
// 
// CALL HRIN (ID,ICYCLE,IOFSET)
// 
// Action: Read a histogram from the current directory on the direct access file (or global section) into the current directory in memory. 
// 
// Input parameters: 
// ID Histogram identifier. ID=0 means that all histograms from the current directory on the direct access file (global section) should be read into memory. If a histogram
//      identifier ID already exists in memory a message is printed and it is deleted from memory before reading the new histogram from the file or global section. 
// ICYCLE 
//      Cycle number. If ICYCLE=0 then the lowest cycle is read. To read the highest cycle, use a large number, e.g. 999999. 
// IOFSET 
//      The histogram which is read in memory will have the identifier IDN=ID+IOFSET. Specifying IOFSET different of zero permits to have in memory copies of
//      histograms with the same identifiers ID in different files. This parameter may be very useful when HRIN is called together with routines such as HOPERA or HDIFF.
//      This facility also works for Ntuples. 
// 
////////////////////////////////////////////////////////////////////////////////

inline void hrin(int id,int n,int offset)
{ hrin_(id,n,offset); }

////////////////////////////////////////////////////////////////////////////////
// 
// CALL HROUT (ID,ICYCLE*,CHOPT)
// 
// Action: Write a histogram from the current directory in memory onto the current directory on the direct access file. 
// 
// Input parameters: 
// ID Histogram identifier. ID=0 means write all histograms from the current directory in memory. 
// CHOPT 
//      Character variable specifying the options selected. 
//      'T' Write the whole directory tree hanging from the current directory (if ID=0). 
// Output parameter: 
// ICYCLE 
//      Cycle number. The first time a given histogram with identifier ID is stored on a directory on a direct access file, ICYCLE is set to 1. If the histogram identifier ID
//      already exists on the direct access file, then the call to HROUT will increment the cycle number ICYCLE by one. 
// 
// Experienced users may invoke routines from the ZEBRA RZ  package to purge a directory, i.e. delete all versions of an identifier but the most recent one using routine
// RZPURG. 
// 
////////////////////////////////////////////////////////////////////////////////

inline void hrout(int id,int &cycle,const char *option)
{ hrout_(id,cycle,option,strlen(option)); }

////////////////////////////////////////////////////////////////////////////////

inline void hrput(int const &id, char const *file_name, char const *options)
{ hrput_(id,file_name,options,strlen(file_name),strlen(options)); }

////////////////////////////////////////////////////////////////////////////////
// 
// CALL HGIVE (ID,CHTITL*,NX*,XMI*,XMA*,NY*,YMI*,YMA*,NWT*,LOC*)
// 
// Action: Returns the booking parameters and address of a given histogram. 
// 
// Input parameter: 
// ID Histogram identifier, cannot be zero. 
// Output Parameters: 
// CHTITL 
//      Histogram title (must be declared CHARACTER*80) 
// NX Number of channels in X 
// XMI Lower edge of first X channel 
// XMA Upper edge of last X channel 
// NY Number of channels in Y (zero for a 1-dimensional histogram) 
// YMI Lower edge of first Y channel 
// YMA Upper edge of last Y channel 
// NWT Number of machine words for the title. If there is no title, NWT is returned as 0 . 
// LOC Address of the histogram in the common /PAWC/.    
// 
////////////////////////////////////////////////////////////////////////////////

inline void hgive(int const &id, char *title,
           int &nx, float &xmin, float &xmax,
           int &ny, float &ymin, float &ymax, int &nwt,
           char *&address)
{
  char *tmp;
  int const title_length = 80;
  title[title_length] = 0;
  hgive_(id,title,nx,xmin,xmax,ny,ymin,ymax,nwt,tmp,title_length);
  for( int i=title_length-1; i>=0 && title[i]==' '; title[i--]=0 );
  if( address!=0 )
    address = tmp;
}

inline void hgive(int const &id, char *title,
           int &nx, float &xmin, float &xmax,
           int &ny, float &ymin, float &ymax, int &nwt)
{
  char *addr;
  hgive(id,title,nx,xmin,xmax,ny,ymin,ymax,nwt,addr);
}

////////////////////////////////////////////////////////////////////////////////

inline void hnoent(int const &id, int &n_events)
{ hnoent_(id,n_events); }

////////////////////////////////////////////////////////////////////////////////

inline int hnoent(int const &id)
{ int n_entries; hnoent_(id,n_entries); return n_entries; }

////////////////////////////////////////////////////////////////////////////////

inline void hindex(void)
{ hindex_(); }

////////////////////////////////////////////////////////////////////////////////

inline float hi(int const &id, int const &i)
{ return hi_(id,i); }

////////////////////////////////////////////////////////////////////////////////

inline float hij(int const &id, int const &i, int const &j)
{ return hij_(id,i,j); }

////////////////////////////////////////////////////////////////////////////////

inline float hx(int const &id, const float &x)
{ return hx_(id,x); }

////////////////////////////////////////////////////////////////////////////////

inline float hxy(int const &id, const float &x, const float &y)
{ return hxy_(id,x,y); }

////////////////////////////////////////////////////////////////////////////////

inline float hxyij(int const &id, const float &x, const float &y, int*i, int *j)
{ return hxyij_(id,x,y,i,j); }

////////////////////////////////////////////////////////////////////////////////

inline float hie(int id,int i)
{ return hie_(id,i); }

////////////////////////////////////////////////////////////////////////////////

inline float hije(int id,int i,int j)
{ return hije_(id,i,j); }

////////////////////////////////////////////////////////////////////////////////

inline float hstati(int id,int w,char const *opt,int n)
{ return hstati_(id,w,opt,n,strlen(opt)); }

////////////////////////////////////////////////////////////////////////////////

inline void hdelet(int id)
{ hdelet_(id); }

////////////////////////////////////////////////////////////////////////////////
// 
// CALL HROPEN (LUN,CHTOP,CHFILE,CHOPT,*LREC*,ISTAT*)
// 
// Action: Open a direct access HBOOK file. If several direct access files are opened, they are identified by the top directory only. 
// 
// Input parameter description: 
// 
// LUN Logical unit number associated to the file. 
// CHTOP 
//      Character variable specifying the name of the top directory associated with unit LUN (maximum 8 characters). This is an arbitrary name used to identify the file on unit
//      LUN in subsequent calls to HR.. routines. 
// CHFILE 
//      Character variable specifying the name of the file to be opened. 
// CHOPT 
//      Character variable specifying the options selected 
//      Medium 
//          ' ' Disk (default) 
//          'G' Global Section (see chapter ) 
//      mode 
//          ' ' Existing HBOOK file (default) 
//          'N' Create a new file 
//          'P' Preserve file case 
//          'Q' Override default number of records for new file with contents of IQUEST(10)   
//          'X' The file is/will be in exchange format 
//          'U' Update an existing file 
// LREC 
//      Record length in machine words (recommended value is 1024). If LREC=0 the actual record length is returned on exit. 
// 
// Output parameter description: 
// 
// ISTAT 
//      Return code. ISTAT=0 indicates success. 
// LREC 
//      (Only when (LREC=0) on input) The actual record length of the file on disk. 
// 
// Remarks: 
// 
//      HROPEN performs a Fortran direct access open for the file CHFILE on logical unit LUN. 
//      On VM/CMS no FILEDEF statement should be given.   The filename CHFILE can be given in either of the forms Filename Filetype Filemode or
//      Filename.Filetype.Filemode 
//      On Unix the filename CHFILE will be translated to lowercase.   
//      If LREC=0 on input, HROPEN will automatically determine the record length of existing files. 
//      On MVS systems, the current userid prefix will be automatically added to the front of the file name unless the first character is a dot (.). 
//      The maximum number of records is by default 32000. You can use option Q to change this. 
//      A file declared with HROPEN must be released with HREND. 
//      HROPEN calls HRFILE internally. 
// 
////////////////////////////////////////////////////////////////////////////////

inline void hropen(int lun, char const *top, char const *file,
                   char const *option,int rec_length, int &error)
{
  hropen_(lun,top,file,option,rec_length,&error,
          strlen(top),strlen(file),strlen(option) );
}

inline int hropen(int lun, char const *top, char const *file,
                   char const *option,int rec_length)
{
  int error;
  hropen_(lun,top,file,option,rec_length,&error,
          strlen(top),strlen(file),strlen(option) );
  return error;
}

////////////////////////////////////////////////////////////////////////////////

inline void hcopy( int id1, int id2, char const *title)
{ hcopy_(id1,id2,title,strlen(title)); }

////////////////////////////////////////////////////////////////////////////////

inline void hopera(int id1,char const *oper,int id2,int id3,float c1=1,float c2=1)
{ hopera_(id1,oper,id2,id3,c1,c2,strlen(oper)); }

////////////////////////////////////////////////////////////////////////////////

inline int hexist(int i)
{ return hexist_(i); }

////////////////////////////////////////////////////////////////////////////////
// 
// CALL HLNEXT (*IDH*,CHTYPE*,CHTITL*,CHOPT)
// 
// Action: Scan the contents of the current directory in memory or on an RZ file. 
// 
// Input parameter description: 
// 
// IDH Must be zero for first call 
// CHTYPE 
//      Character variable specifying items to be scanned. 
//      '1' include 1-D histograms 
//      '2' include 2-D histograms 
//      'N' include Ntuples 
//      'D' include subdirectories 
//      ' ' include everything, i.e., equivalent to '12ND'. 
// 
// Output parameter description: 
// 
// IDH On return contains identifier of next histogram. When all histograms are processed, a value of zero is returned. 
// CHTYPE 
//      Character variable specifying type of histogram. 
//      '1' 1-dimensional 
//      '2' 2-dimensional 
//      'N' Ntuple 
//      'D' subdirectory 
//      '?' unknown. 
// CHTYPE 
//      Character variable containing title or subdirectory name. 
// 
// Scan content of current directory
// 
//       IDH=0
//   1   CONTINUE
//       CALL HLNEXT(IDH,CHTYPE,CHTITL,CHOPT)
//       IF(IDH.NE.0) THEN
//          ... process
//          GOTO 1
//       ENDIF
// 
////////////////////////////////////////////////////////////////////////////////

inline void hlnext(int &id,char *type,char *title,const char *opt=" ")
{
  const int l1=5, l2=81;
  type[l1-1] = title[l2-1] = 0;
  hlnext_(id,type,title,opt, l1,l2,strlen(opt));
  for( int i=l1-1; i>=0 && type [i]==' '; type [i--]=0 );
  for( int i=l2-1; i>=0 && title[i]==' '; title[i--]=0 );
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

extern struct hcfitd_struct
{
  double par[24],fit_result;
} hcfitd_;

////////////////////////////////////////////////////////////////////////////////

#endif // CsHbookProto_h
