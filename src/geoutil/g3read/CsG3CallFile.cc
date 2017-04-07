/*!
  \file		CsG3CallFile.cc
  \brief	A implementation file for a class which convertes geometry information 
			from GEANT 3 call list to CORAL enviroment.
  \author	$Author: tnagel $
  \version	$Revisopn$
  \date		$Date: 2010/01/28 12:51:25 $
*/
#if __GNUG__ >= 2
#  pragma implementation
#endif

#include <cmath>


#include "CsG3CallFile.h"
#include "CsRegistry.h"
#include "CsErrLog.h"
#include "CsOpt.h"

#include <algorithm>

using namespace std;
using namespace CLHEP;

CsG3CallFile* CsG3CallFile::instance_ = NULL ;

CsG3CallFile* CsG3CallFile::Instance(){

  if(instance_ == NULL){
	  instance_ = new CsG3CallFile();

	  CsRegistry csReg_;
	  if( csReg_.EOJRegistration((CsEndOfJob*) instance_) ){
		  CsErrLog::Instance()->mes(elDebugging, 
			  "CsG3CallFile has been registerd successfully.");;
	  } else {
		  CsErrLog::Instance()->mes(elError, 
			  "CsG3CallFile has not been registerd successfully.");;
	  }


  }

  return instance_;

}

bool CsG3CallFile::end(){
		  CsErrLog::Instance()->mes(elDebugging, 
			  "CsG3CallFile is goind to die.");;
  return true;
}


//---- CsG3CallFile ---------------------------------------------------------
CsG3CallFile::CsG3CallFile() :
  g3CallListFile_(""),
  volume_(), 
  position_(),
  medium_(), 
  material_(),
  rotation_(), 
  detector_(),
  HALL_("HALL", 1, "", HepGeom::Vector3D<double>(0.0, 0.0, 0.0), 0, ""),
  HALLVolume_("HALL", "BOX", 1, 3, vector<double>(3, 10000.0) ),
  systemMode_(CsG3CallFile::COMPASS) {

  this->rotation()[0] = HepMatrix(3,3,1);	

  string strMode;
  CsOpt::Instance()->getOpt( "GDB", "system", strMode);

  if(strMode == "COMGEANT"){
	  this->systemMode(CsG3CallFile::COMGEANT);
	  CsErrLog::Instance()->mes(elDebugging, 
		  "All geometery data will be treated as in COMGEANT. (not in CORAL frame system)" );		
  }

  if( this->readCallList() ) {
	  CsErrLog::Instance()->mes(elDebugging, 
		  this->g3CallListFile() + 
		  " has been successfully read." );
  } else {
	  CsErrLog::Instance()->mes(elError, 
		  this->g3CallListFile() + 
		  " has not been successfully read. (ERROR)" );
  }

  HALL_.pVolume( &HALLVolume_ );


}

//CsG3CallFile::CsG3CallFile(const CsG3CallFile& callFile) : 
//	g3CallListFile_(callFile.g3CallListFile()),
//	volume_(callFile.volume()), 
//	position_(callFile.position()),
//	medium_(callFile.medium()), 
//	material_(callFile.material()),
//	rotation_(callFile.rotation()), 
//	detector_(callFile.detector()) {
//}


CsG3CallFile::~CsG3CallFile(){
}

//CsG3CallFile& CsG3CallFile::operator=(const CsG3CallFile& callFile){
//	if(this != &callFile){
//		this->g3CallListFile(	callFile.g3CallListFile());
//		this->volume(			callFile.volume());
//		this->position(			callFile.position());
//		this->medium(			callFile.medium());
//		this->material(			callFile.material());
//		this->rotation(			callFile.rotation());
//		this->detector(			callFile.detector());
//	}
//	return *this;
//}

bool CsG3CallFile::getFileName(){
  if(this->g3CallListFile() == ""){
	  string g3callFile;
	  string mes;

	  mes += "The g3call file ";    
	  if( CsOpt::Instance()->getOpt( "Geom", "g3callFile", g3callFile) ){
		  mes += " was defined in CORAL option file ";
	  } else {
		  mes += " has not been defined in CORAL options file ";
		  CsErrLog::Instance()->mes(elFatal, mes);
		  return false;
	  }
	  mes += g3callFile + " is being used.";
	  CsErrLog::Instance()->mes(elDebugging, mes);

	  this->g3CallListFile(g3callFile);
  }
  return true;
}

bool CsG3CallFile::readCallList(){
  if( ! this->getFileName() ) {
	  CsErrLog::Instance()->mes(elError, 
		  "CsG3CallFile could not get filename.");
	  return false;
  }

  ifstream infile( this->g3CallListFile().c_str()  );

  if(infile.fail()){
	  string mes = "File open failed:\t" + this->g3CallListFile();
	  CsErrLog::Instance()->mes(elError, mes);
	  return false;
  }

  CsErrLog::Instance()->mes(elDebugging, 
	  "G3 call list is opened.");	

  this->readCallList(infile);
  return true;	
}

bool CsG3CallFile::readCallList(istream& is){
  const int max_line_size(1024);
  char line[max_line_size];


  CsErrLog::Instance()->mes(elDebugging, 
	  "G3 call list is now been reading.");	


  while((is.getline(line, max_line_size)) != NULL){

	  istringstream  stream(line);
	  string flag, opt;
	  stream >> flag >> opt;

	  if( flag == "----" ){
		  if(opt == "GSMATE")	this->readGSMATE(stream);
		  if(opt == "GSROTM")	this->readGSROTM(stream);
		  if(opt == "GSMIXT")	this->readGSMIXT(stream);
		  if(opt == "GSVOLU")	this->readGSVOLU(stream);
		  if(opt == "GSPOS" )	this->readGSPOS(stream);
		  if(opt == "GSPOSP")	this->readGSPOS(stream);
		  if(opt == "GSTMED")	this->readGSTMED(stream);
		  if(opt == "GSDET")	this->readGSDET(stream);
		  if(opt == "GSDETH")	this->readGSDETH(stream);
		  if(opt == "GSDETD")	this->readGSDETD(stream);
		  if(opt == "GSDETU")	this->readGSDETU(stream);
	  }
  }

  CsErrLog::Instance()->mes(elDebugging, "Reading Done.");

// link medium to its material
  for(	med_iterator imed = this->medium().begin();
	  imed != this->medium().end();
	  imed++){
	  mat_iterator imat( findMaterial( (*imed).materialId() ) );
	  if( imat != this->material().end() ) (*imed).pMaterial( &(*imat) );
  }


  CsErrLog::Instance()->mes(elDebugging, "Material Information have been linked.");

// link volume to its medium
  for(	vol_iterator iv = this->volume().begin();
	  iv != this->volume().end();
	  iv++){
	  med_iterator imed( findMedium((*iv).materialId()) );
	  if(imed != this->medium().end()) (*iv).pMedium( &(*imed) );
  }

  CsErrLog::Instance()->mes(elDebugging, "Volume Information have been linked to Material information.");

// link its	mother positions
//			rotation matrics
//			detector information

  for(	pos_iterator ip = this->position().begin();
	  ip != this->position().end();
	  ip++){

	  (*ip).pRotation( &(this->rotation()[(*ip).rotation()]) );

	  vol_iterator ivol( this->findVolume( (*ip).name() ) );
	  if( ivol != this->volume().end() )(*ip).pVolume( &(*ivol) );

	  det_iterator idet( findDetector( (*ip).name() ) );
	  if(idet != this->detector().end()) (*ip).pDetector( &(*idet) );
  }


  CsErrLog::Instance()->mes(elDebugging, "Position Information have been linked to Volume Information.");



// find mother volume
  for(	pos_iterator ip = this->position().begin();
	  ip != this->position().end();
	  ip++) {

	  this->findMothers( ip );
  }

  CsErrLog::Instance()->mes(elDebugging, 
	  "Mother Volume Information have been linked.");

// Link Rotation matrix information.
  if( this->systemMode() == CsG3CallFile::COMPASS ){

	  HepMatrix mGtoC(3,3,1), mTtoB(3,3,1);
	  double GeanToCoral[9]  = { 0, 1, 0, 0, 0, 1, 1, 0, 0 };
	  double TubeToBox[9]    = { 0, 0, 1, 1, 0, 0, 0, 1, 0 };

	  for( int i=0; i<9; i++ ) {
		  mGtoC(   i/3+1, i%3+1 ) = GeanToCoral[i];
		  mTtoB(   i/3+1, i%3+1 ) = TubeToBox[i];
	  }

	  for(	pos_iterator ip = this->position().begin();
		  ip != this->position().end();
		  ip++) {

		  if(	(*ip).pVolume()->type() == "TUBE" ||
			  (*ip).pVolume()->type() == "TRD1" ||
			  (*ip).pVolume()->type() == "TRD2" ||
			  (*ip).pVolume()->type() == "CONE" ||
			  (*ip).pVolume()->type() == "SPHE") {

  //			if((*ip).pMother() != NULL)
			  if((*ip).pMother() != &HALL_)
			  if(	(*ip).pMother()->pVolume()->type() != "TUBE" &&
				  (*ip).pMother()->pVolume()->type() != "TRD1" &&
				  (*ip).pMother()->pVolume()->type() != "TRD2" &&
				  (*ip).pMother()->pVolume()->type() != "CONE" &&
				  (*ip).pMother()->pVolume()->type() != "SPHE" ) {

				  (*ip).rM(this->rotation()[(*ip).rotation()] * mGtoC );
			  }
		  }
		  else {
			  (*ip).rM( this->rotation()[(*ip).rotation()]);
		  }

  //		if((*ip).pMother() != NULL)
		  if((*ip).pMother() != &HALL_)
			  if(	(*ip).pMother()->pVolume()->type() == "TUBE" ||
				  (*ip).pMother()->pVolume()->type() == "TRD1" ||
				  (*ip).pMother()->pVolume()->type() == "TRD2" ||
				  (*ip).pMother()->pVolume()->type() == "CONE" ||
				  (*ip).pMother()->pVolume()->type() == "SPHE" ) { 
				  (*ip).center( mTtoB * (*ip).center() ); 
			  }
	  }

  } else {
	  for(	pos_iterator ip = this->position().begin();
		  ip != this->position().end();
		  ip++) {
		  (*ip).rM( this->rotation()[(*ip).rotation()]);
	  }		
  }


  this->HALL().checkDetector();

  CsErrLog::Instance()->mes(elDebugging, "Position structure have been successfully determined.");


// numbering
  int posId(1);
  for(	CsG3CallFile::pos_iterator ip = this->position().begin();
	  ip != this->position().end();
	  ip++){
	  (*ip).id(posId++);
  }	

//#ifdef DEBUG
//	cout << *this << endl;
//#endif //DEBUG

  return true;
}

void CsG3CallFile::findMothers(const pos_iterator& ip){
  string motherName((*ip).mother());

  string mes = "findMothers(" + (*ip).name() + "):\t" +
	  "MOther volume(" + (*ip).mother() + ")";

//	CsErrLog::Instance()->mes(elDebugging, mes);
//	cout << mes << endl;

  if( (*ip).pMother() != NULL  ){
	  return;
  }

  if( (*ip).mother() == "HALL" ){
	  (*ip).pMother( &HALL_ );
	  HALL_.addToDaughter( &(*ip) );
	  return; 
  }

  pos_iterator ipm( this->position().begin() );

  unsigned int number((*ip).n()); // 1

  while(	(ipm = findPosition(  motherName, ipm) )  
			  != this->position().end() ) {

	  this->findMothers( ipm ); 

	  if( (*ipm).n() == 1 && (*ipm).unit() == 1){

		  (*ip).pMother( &(*ipm) ); // defined mother volume in Gspos
		  (*ipm).addToDaughter( &(*ip) );
		  (*ip).motherN( (*ipm).n() );

	  } else {

//			Gspos newPos( (*ip) );

		  this->position().push_back( Gspos( (*ip) ) );
		  pos_iterator newEntry( this->position().end() );
		  newEntry--;

		  (*newEntry).pMother( &(*ipm) );
		  (*newEntry).motherN( (*ipm).n() );
		  (*newEntry).n(number);

		  (*ipm).addToDaughter( &(*newEntry) );

	  }

	  ipm++;
	  number++;
  }
}


ostream& operator<<( ostream& stream, CsG3CallFile& callFile){
  int n(1);


  stream << "================ START ========================" << endl;

  stream << "\n\tMaterial Table Information" << endl;
  for(	CsG3CallFile::mat_iterator imat = callFile.material().begin();
	  imat != callFile.material().end();
	  imat ++) {
	  stream	<<	n++ << "th\t"
			  << (*imat)	<< endl;
  }

  stream << "\n\tGSMED entry" << endl;
  n=1;
  for(	CsG3CallFile::med_iterator imed = callFile.medium().begin();
	  imed != callFile.medium().end();
	  imed++){
	  stream	<< n++ << "th\t" << (*imed) << endl;
  }

  stream << "\n\tGSDET entry" << endl;
  n=1;
  for(	CsG3CallFile::det_iterator idet = callFile.detector().begin();
	  idet != callFile.detector().end();
	  idet++){
	  stream 	<< n++ << '\t'
			  << (*idet).name() << '\t'
			  << (*idet).mother() << '\t'
			  << (*idet).type() << '\t'
			  << endl;
  }


  stream << "\n\tGSVOL entry" << endl;
  for(	CsG3CallFile::vol_iterator ip = callFile.volume().begin();
	  ip != callFile.volume().end();
	  ip++){
	  stream << (*ip) << endl;
  }


  stream << "\n\tGSPOS entry" << endl;
  for(	CsG3CallFile::pos_iterator ip = callFile.position().begin();
	  ip != callFile.position().end();
	  ip++){
	  stream << (*ip) << endl;
  }


  stream << "\n\tVolume Structure" << endl;

  callFile.HALL().dumpDaughters(stream);

  stream << "================ END ========================" << endl;

  return stream;
}

void CsG3CallFile::readGSDET(istream& stream){
  // BEGIN --------------------------------------------------				
  string mother, name, motherA;
  int dummyA, dummyB, dummyC, dummyD;
  int materialId;

  stream	>> mother >> name >> dummyA >> motherA 
		  >> dummyB >> materialId >> dummyC >> dummyD;

  this->detector().push_back(Gsdet(name, mother, materialId));
}

void CsG3CallFile::readGSDETH(istream& stream){
  string mother, name;
  stream >> mother >> name;
  det_iterator idet( findDetector( name ) );
  if(idet != this->detector().end() ){
	  int npara; stream >> npara;
	  for(int i=0;i<npara;i++) {
		  string n; stream >> n;
		  (*idet).addHitName(n);
	  }
	  for(int i=0;i<npara;i++) {
		  int b; stream >> b;
		  (*idet).addHitBit(b);
	  }
	  for(int i=0;i<npara;i++) {
		  double f; stream >> f;
		  (*idet).addHitFactor(f);
	  }
	  for(int i=0;i<npara;i++) {
		  double f; stream >> f;
		  (*idet).addHitOffset(f);
	  }
  }
}

void CsG3CallFile::readGSDETD(istream& stream){
  string mother, name;
  stream >> mother >> name;
  det_iterator idet( findDetector( name ) );
  if(idet != this->detector().end() ){
	  int npara; stream >> npara;
	  for(int i=0;i<npara;i++) {
		  string n; stream >> n;
		  (*idet).addDigitName(n);
	  }
	  for(int i=0;i<npara;i++) {
		  int b; stream >> b;
		  (*idet).addDigitBit(b);
	  }		
  }
}

void CsG3CallFile::readGSDETU(istream& stream){
  string mother, name;
  stream >> mother >> name;
  det_iterator idet( findDetector( name ) );
  if(idet != this->detector().end() ){
	  int npara; stream >> npara;
	  for(int i=0;i<npara;i++) {
		  double para; stream >> para;
		  (*idet).addUser(para);
	  }
  }
}

void CsG3CallFile::readGSROTM(istream& stream){
  unsigned int id;
  HepMatrix rotM(3,3);

  stream >> id;

//	cout << id ;
  for(int i=1; i<4; i++) {
	  double theta, phi; stream >> theta >> phi;

	  double	xTheta(M_PI * theta / 180.0), 
			  xPhi(  M_PI * phi   / 180.0);

	  rotM(1,i) = sin(xTheta) * cos(xPhi);
	  rotM(2,i) = sin(xTheta) * sin(xPhi);
	  rotM(3,i) = cos(xTheta);

  }
  this->rotation()[ id ] = rotM;

}

void CsG3CallFile::readGSVOLU(istream& stream){
  string name, type;
  unsigned int id; 
  int nparam;
  vector<double> param;

  stream >> name >> type >> id >> nparam;

  for(int i=0;i<nparam;i++) {
	  double tmp; stream >> tmp;
	  param.push_back(tmp);
  }
  this->volume().push_back(Gsvolu(name, type, id, nparam, param));
}

void CsG3CallFile::readGSPOS(istream& stream){

  string name;
  stream >> setw(5) >> name;

  string mother, flag;
  unsigned int id;

  int rotation;
  double x, y, z;

  stream >> id >> mother >> x >> y >> z >> rotation >> flag;

//	HepGeom::Vector3D<double> center(10.0*x, 10.0*y, 10.0*z); //(in mm)
  HepGeom::Vector3D<double> center(x, y, z); //(as it is)
  Gspos pos( name, id, mother, center, rotation, flag);

  int npar;
  if(stream >> npar){
	  for(int i=0; i<npar; i++){
		  double par; stream >> par;
		  pos.addParam(par);
	  }
  }

  this->position().push_back(	pos );

}

string CsG3CallFile::readName(istream& stream){
  const int name_length(64);
  char buf[name_length], c;
  stream >> c;
  if(c != '\"') exit(0); 
  stream.getline(buf, name_length, '\"');
  return string(buf);
}

void CsG3CallFile::readGSMATE(istream& stream){
  unsigned int id;
  stream >> id;

  string name(this->readName(stream));
  replace( name.begin(), name.end(), ' ', '-');

  int npar;
  double a, z, dens, radl;
  stream >> a >> z >> dens >> radl >> npar;

  this->material().push_back( 
	  CsMaterial( id, name, dens, radl ) );
}

void CsG3CallFile::readGSMIXT(istream& stream){

  unsigned int id;
  stream >> id;


  string name(this->readName(stream));
  replace( name.begin(), name.end(), ' ', '-');

  int nmat;
  double dens, radl(0.0);
  stream >> dens >> nmat;

  const int npara = abs(nmat);

  vector<double> a(nmat, 0.0), z(nmat, 0.0), w(nmat, 0.0);

  for(int i=0; i<npara; i++) stream >> a[i];
  for(int i=0; i<npara; i++) stream >> z[i];
  for(int i=0; i<npara; i++) stream >> w[i];	

  this->material().push_back( CsMaterial( id, name, dens, nmat, a, z, w) );

}

void CsG3CallFile::readGSTMED(istream& stream){

  unsigned int id;
  stream >> id;

  string name(this->readName(stream));
  replace( name.begin(), name.end(), ' ', '-');

  unsigned int materialId;
  stream >> materialId;

  int isvol, ifield;
  stream >> isvol >> ifield;
  this->medium().push_back( Gstmed(id, name, materialId, isvol, ifield) );
}

CsG3CallFile::vol_iterator CsG3CallFile::findVolume(
	  const string& name){
  return this->findVolume(name, this->volume().begin() );
}

CsG3CallFile::vol_iterator CsG3CallFile::findVolume(
	  const string& name, 
	  const CsG3CallFile::vol_iterator& istart){
  vol_iterator iv;
  for(	iv = istart; iv != this->volume().end(); iv++){
	  if((*iv).name() == name ) break;
  }
  return iv;
}

CsG3CallFile::det_iterator CsG3CallFile::findDetector(
  const string& name){
  return this->findDetector( name, this->detector().begin());
}

CsG3CallFile::det_iterator CsG3CallFile::findDetector(	
	  const string& name, 
	  const CsG3CallFile::det_iterator& istart){
  det_iterator iv;
  for(	iv = istart; iv != this->detector().end(); iv++){
	  if((*iv).name() == name ) break;
  }
  return iv;
}

CsG3CallFile::mat_iterator CsG3CallFile::findMaterial(const unsigned int& id){
  return this->findMaterial(id, this->material().begin());
}

CsG3CallFile::mat_iterator CsG3CallFile::findMaterial(
	  const unsigned int& id, 
	  const CsG3CallFile::mat_iterator& istart){
  mat_iterator im;
  for(	im = istart; im != this->material().end(); im++){
	  if((*im).id() == id ) break;
  }
  return im;
}

CsG3CallFile::med_iterator CsG3CallFile::findMedium(const unsigned int& id){
  return this->findMedium(id, this->medium().begin());
}


CsG3CallFile::med_iterator CsG3CallFile::findMedium(
	  const unsigned int& id, 
	  const CsG3CallFile::med_iterator& istart){
  med_iterator im;
  for(	im = istart; im != this->medium().end(); im++){
	  if((*im).id() == id ) break;
  }
  return im;
}

CsG3CallFile::pos_iterator CsG3CallFile::findPosition(const Gspos* pPos){
  return this->findPosition(pPos, this->position().begin());
}

CsG3CallFile::pos_iterator CsG3CallFile::findPosition(
	  const Gspos* pPos, 
	  const CsG3CallFile::pos_iterator& istart){

  pos_iterator iv;

  for(	iv = istart; iv != this->position().end(); iv++){
	  if(&(*iv) == pPos ) break;
  }

  return iv;
}

CsG3CallFile::pos_iterator CsG3CallFile::findPosition(const string& name){
  return this->findPosition(name, this->position().begin());
}

CsG3CallFile::pos_iterator CsG3CallFile::findPosition(
	  const string& name, 
	  const CsG3CallFile::pos_iterator& istart){

  pos_iterator iv;

  for(	iv = istart; iv != this->position().end(); iv++){
	  if((*iv).name() == name ) break;
  }

  return iv;
}

CsG3CallFile::pos_iterator CsG3CallFile::findPosition(const int& id){
  return this->findPosition(id, this->position().begin());
}

CsG3CallFile::pos_iterator CsG3CallFile::findPosition(
	  const int& id, 
	  const CsG3CallFile::pos_iterator& istart){

  pos_iterator iv;

  for(	iv = istart; iv != this->position().end(); iv++){
	  if((*iv).id() == id ) break;
  }

  return iv;
}

HepGeom::Point3D<double> CsG3CallFile::centerInHALL(const Gspos& pos){
  if( pos.pMother() == NULL) 				return  pos.center();
  if( pos.pMother()->name() ==  "HALL" )	return  pos.center();

  return	this->centerInHALL( *( pos.pMother() ) ) + 
		  this->rotationToHALL( *( pos.pMother() ) ) * 
		  pos.center();

}

HepMatrix CsG3CallFile::rotationToHALL(const Gspos& pos){	
  if( pos.pMother() == NULL)				return  pos.rM();
  if( pos.pMother()->name() ==  "HALL" )	return  pos.rM();

  return this->rotationToHALL( *(pos.pMother()) ) * pos.rM();

}

