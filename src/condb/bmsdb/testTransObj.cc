/*!
   \file    testTransObj.cc
   \brief   Test Program for Transient object.
   \author  Yoshiyuki Miyachi
   \version $Revision: 1.1 $
   \date    $Date: 2000/06/06 09:48:36 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif
#include "CsSTD.h"
#include "CsBMSconstants.h"

//---- newclass ---------------------------------------------------------
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      test procedure to check class CsBMSconstants
//
int main(void){	
        char header1[80];
	double avsg1[4];
	double evec1[4][4];
	double coff1[4][3][5];
        double p1, p2;
        char filename[] = "bmscoef.dat";
        CsBMSconstants cons;
        cons.readBMSconst(filename);
        cons.prntBMSconst();
        cons.getBMSconst(header1, p1, avsg1, evec1, coff1 );	
        cout << "-------------------------- main print ---------------------"<<endl;    
            cout << resetiosflags(ios::scientific);
            cout << setiosflags(ios::fixed|ios::showpoint);
            cout << header1 << endl;
            cout << "Pbeam = " << p1 << endl;
//
            cout <<  " AVSG16:" << endl;
            cout << setprecision(5); 
            int i, j, ll;
 	    for( i=0; i<4; i++){
            cout << setw(16) << avsg1[i]; avsg1[i]+= 1.;}
            cout << endl;
//
            cout <<  " EVEC16:" << endl;
	    for( i=0; i<4; i++){
	       for( j=0; j<4;j++){
                  cout << setw(16) << evec1[j][i]; evec1[j][i]+= 1.;}
               cout << endl;}
//
            cout <<  " COFF16:" << endl;
	    for( i=0; i<5; i++){
	       for( j=0; j<3;j++){
	         for(  ll=0; ll<4 ;ll++){
                    cout << setw(16) << coff1[ll][j][i]; coff1[ll][j][i]+= 1.;}
                 cout << endl;
                 }
               }
        cout << "---------------------end of main print ---------------------"<<endl;    
        cons.setBMSconst(header1, p1, avsg1, evec1, coff1 );	
        cons.prntBMSconst();
        cout <<"===================== make object cons1 ========================="<< endl;
        CsBMSconstants cons1(header1, p1, avsg1, evec1, coff1 );
        cons1.prntBMSconst();

	exit(0);
}
