// $Id$

#ifndef TMtx_h
#define TMtx_h

/*!
  \brief Scalar / Vector / Matrix

  Matrix class. Valid also for scalars (like 1x1 matrix) and vectors (like Mx1 matrix).
*/

class TMtx {

public:
  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  TMtx();                   //!< 1x1 matrix constructor (scalar)
  TMtx(int m);              //!< Mx1 matrix (vector) constructor
  TMtx(int m, int n);       //!< MxN matrix constructor
  TMtx(const TMtx&);     //!< copy constructor

  // Desctructor
  ~TMtx();                  //!< destructor

  //Methods
  void Print(const char* str="") const; //!< Print matrix
  TMtx t();                       //!< Transpose matrix
  TMtx invs(int& ierr);           //!< Invert symmetric matrix
  TMtx i5 (int& ierr);            //!< Invert 5x5 matrix
  TMtx i  (int& ierr);            //!< Invert matrix

  // Overloaded operators
  TMtx& operator = (const TMtx&); //!< "="  operation 
  TMtx  operator * (const TMtx&); //!< "*"  operation
  TMtx  operator + (const TMtx&); //!< "+"  operation
  TMtx  operator - (const TMtx&); //!< "-"  operation
  TMtx& operator +=(const TMtx&); //!< "+=" operation
  TMtx& operator *=(const double &); //!< "*=" operation

  double & operator () (const int i, const int j=1) const; //!< accessor to i,j element 
  operator double();                                 //!< conversion 1x1 matrix to double 
    
private:

  //Matrix dimensions
  int m;
  int n;
  // Base pointer
  double **Mtx;

};


//
// Constructors
//

inline TMtx::TMtx()
{
  NobjCreated++;
  m=1; n=1;
  Mtx = new double* [m];
  for(int i=0; i<m; i++){
    Mtx[i]= new double [n];
  }
};

inline TMtx::TMtx(int i)
{
  NobjCreated++;
  m=i; n=1;
  Mtx = new double* [m];
  for(int i=0; i<m; i++){
    Mtx[i]= new double [n];
  }
};

inline TMtx::TMtx(int i, int j)
{
  NobjCreated++;
  m=i; n=j;
  Mtx = new double* [m];
  for(int i=0; i<m; i++){
    Mtx[i]= new double [n];
  }
};

//
// Copy constructor
//

inline TMtx::TMtx(const TMtx& M)
{
  NobjCreated++;
//  cout<<"Copy constructor for ("<<M.m<<"x"<<M.n<<")"<<endl; 
  m=M.m; n=M.n;
  Mtx = new double* [m];
  for(int i=0; i<m; i++){
    Mtx[i]= new double [n];
    for(int j=0; j<n; j++)
      Mtx[i][j]=M.Mtx[i][j];
  }
}

// Destructor

inline TMtx::~TMtx()
{
  NobjDestructed++;
  for( int i=0; i<m; i++)
    delete [] Mtx[i];
  delete [] Mtx;
};

#endif














