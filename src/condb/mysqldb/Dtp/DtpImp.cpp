#include "DtpImp.h"
#include <qcombobox.h>
#include <qpushbutton.h>
#include <qdatetime.h>
#include <qdatetimeedit.h>
#include <qtable.h>
#include <qmessagebox.h>
#include <qfiledialog.h>
#include <qfile.h>
#include <qtextedit.h>
#include <qregexp.h>
#include <qlineedit.h>
#include <qlcdnumber.h>
#include <qcheckbox.h>
#include "registerdialog.h"

//#include </afs/cern.ch/user/c/cbernet/Src/Dtp/FileBrowser/dirview/dirview.h>

#include <regex.h>
#include <unistd.h>
#include <stdlib.h>
#include <fstream.h>

const char* DtpImp::timeFormat_ = "yyyyMMddhhmmss";
const char* DtpImp::timeFormatSasha_ = "yyyy-MM-dd-hh:mm:ss";
const char* DtpImp::persoParams_ = ".dtp";

#include "Pics/window-destroy.xpm"
#include "Pics/empty.xpm"

DtpImp::DtpImp() : 
  Dtp(), db_(0), 
  firstName_(""), lastName_(""), email_("") {

  emptyXpm_ = new QPixmap(( const char** ) magick);
  window_destroyXpm_ = new QPixmap(( const char** ) PenguinCrossed_xpm);

  connectDB();  
  registerUser();

  QDateTime start(QDate(2002,01,01),QTime(0,0));
  startTimeEdit_->setDateTime(start);
  QDateTime end(QDate(2003,01,01),QTime(0,0));
  endTimeEdit_->setDateTime(end);
  dbtable_->setSelectionMode(QTable::MultiRow);
  
}

DtpImp::~DtpImp() {

  if(db_) delete db_;
}

bool DtpImp::registerUser() {
   
  // does the user already have a perso param  file in home directory ?
  string home = getenv("HOME");
  string pfile = home; pfile += "/"; pfile += DtpImp::persoParams_;
  
  ifstream in(pfile.c_str());
  
  if(!in) {
    // no perso file, registration dialog

    QStringList params;
    
    RegisterDialog rd(this, "DB registration", true);
    if(rd.exec() == QDialog::Accepted) {
      rd.parameters(params);

      login_    =params[0].latin1();
      firstName_=params[1].latin1();
      lastName_ =params[2].latin1();
      email_    =params[3].latin1();
      
      
      // creating perso params file
      ofstream out(pfile.c_str());
      if(out) {
	out << (params[0].latin1()) <<endl;
	out << (params[1].latin1()) <<endl;
	out << (params[2].latin1()) <<endl;
	out << (params[3].latin1()) <<endl;
	out.close();
      }
      else {
	cerr<<"DtpImp::registerUser(). cannot create personal parameters file"
	    <<endl;
      }
    }
    else {
      // you don't want to register ? get lost.
      exit(1);
    }
  }
  else {
    in>>login_;
    in>>firstName_;
    in>>lastName_;
    in>>email_;
  }
  loginEdit_->setText(login_.c_str());
  firstNameEdit_->setText(firstName_.c_str());
  lastNameEdit_->setText(lastName_.c_str());
  emailEdit_->setText(email_.c_str()); 
  
  authorID_ = db_->authorID(login_,firstName_,lastName_,email_);
  if(authorID_.empty()) {
    cerr<<"Sorry, cannot register to DB"<<endl;
    return false;
  }
  return true;
}

bool DtpImp::connectDB() {

  if(db_) {
    if( !db_->isConnected() ) delete db_;
    else {
      db_->disconnect();
      connectButton_->setText("Connect");
      return true;
    }
  }
  
  if(servercombo_->currentText() && username_ && password_ &&
     dbcombo_->currentText()) {
    try {
      db_ = new MyInterface(servercombo_->currentText(), 
			    username_, password_,
			    dbcombo_->currentText());
    }
    catch(const char* errmsg) {
      cout<<"DtpImp exception "<<errmsg<<endl;
      delete db_;
      db_=0;
    }
  }  
  if(db_) 
    if(db_->connect()) {
      connectButton_->setText("disconnect");
      locationDB();
      return true;    
    }
  return false;
}


int DtpImp::selectDB() {
  
  if(!db_ || !db_->isConnected()) {
    NotConnected();
    return -1;
  }

  string dir(locationCombo_->currentText());
  string dirid = db_->cdbDirID(dir);
  
  int n=0;
  if(dirid.empty()) {
    dir += " is not in the DB ...\n"; 
    QMessageBox mb("Dtp", 
		   dir.c_str(), 
		   QMessageBox::Information,
		   QMessageBox::Ok  | QMessageBox::Default,
		   QMessageBox::NoButton,
		   QMessageBox::NoButton );  
    mb.exec();  
  }
  else {
    string qry = "select StartTime, EndTime, FileName from runlb.tb_calibration where DetName like '"; 
    qry += namefilter_; qry += "'";
    qry += " and StartTime >= ";
    qry += starttime_.toString(DtpImp::timeFormat_);
    qry += " and EndTime <= ";
    qry += endtime_.toString(DtpImp::timeFormat_);
    qry += " and dirID = ";
    qry += dirid;
    
    if(myFilesCheck_->isChecked()) {
      qry += " and authorID = "; 
      qry += authorID_;
    }

    // cout<<qry<<endl;
    
    db_->query(qry);
    n=db_->getNRows();
    dbtable_->setNumRows(n);
    int i=0;
    do {
      //    cout<<db_->getCol(0)<<"\t"<<db_->getCol(1)<<"\t"<<db_->getCol(2)<<endl;
      dbtable_->setText(i,0,db_->getCol(0));
      dbtable_->setText(i,1,db_->getCol(1));
      dbtable_->setText(i,2,db_->getCol(2));
      i++;
    } while(db_->getNextRow());
    

    db_->endQuery();
    
    dbtable_->adjustColumn(0);
    dbtable_->adjustColumn(1);
    dbtable_->adjustColumn(2);
    
    nSelectedLCD_->display(n);
  }
  return n;
}

void DtpImp::downloadDB() {

  if(!db_ || !db_->isConnected()) {
    NotConnected();
    return;
  }

  string dir(locationCombo_->currentText());  
  int cs = dbtable_->currentSelection();
  if(cs != -1) { 
    QString localdir = QFileDialog::getExistingDirectory("./",
							 this,
							 "get existing directory",
							 "Choose a directory",
							 true );  
    if(!localdir) return;
    
    QTableSelection selection = dbtable_->selection(cs);
    
    int n=0;
    for(int i=selection.topRow(); i<=selection.bottomRow(); i++) {
      string filename = dir; filename += "/";
      filename += dbtable_->text(i,2);
      
      string cp = "cp ";
      cp += filename;
      cp += " ";
      cp += localdir;
      
      if(system(cp.c_str())) {
	string errmsg = "cannot copy :\n"; errmsg += cp;
	errmsg += "\nContinue ?";
	QMessageBox mb("Dtp", 
		       errmsg.c_str(), 
		       QMessageBox::Warning,
		       QMessageBox::Yes | QMessageBox::Default,
		       QMessageBox::No | QMessageBox::Escape,
		       QMessageBox::NoButton );
	if(mb.exec() == QMessageBox::No) break;
      }  
      else n++;
    }
    char msg[100];
    sprintf(msg,"%d files downloaded.", n);
    QMessageBox mb("Dtp", 
		   msg, 
		   QMessageBox::Information,
		   QMessageBox::Ok  | QMessageBox::Default,
		   QMessageBox::NoButton,
		   QMessageBox::NoButton );  
    mb.exec();    
  }
}


void DtpImp::uploadDB() {

  if(!db_ || !db_->isConnected()) {
    NotConnected();
    return;
  }
  
  // testing that endtime > starttime
  if(starttime_ >= endtime_) {
    QMessageBox mb("Dtp", 
		   "Sorry, end time must be after start time", 
		   QMessageBox::Critical,
		   QMessageBox::Ok | QMessageBox::Default,
		   QMessageBox::NoButton,
		   QMessageBox::NoButton );
    
    mb.exec();
    return;
  }

  // testing that TBName is good
  bool batchmode = false;
  string msg;
  bool warn=false;
  if(! db_->checkTBName(namefilter_.latin1())) {
    msg = "TBNames MUST have 8 letters\n";
    msg += namefilter_; 
    msg += " is assumed to be a filter.\n";    
    warn=true;
  }
  if(! applyRangeCheck_->isChecked()) {
    msg += "\nApply Validity Range is NOT checked.\nwill try to guess it from the file name\n";
    warn=true;
  }
  if(warn) {
    QMessageBox mb("Dtp", 
		   msg.c_str(), 
		   QMessageBox::Information,
		   QMessageBox::Ok | QMessageBox::Default,
		   QMessageBox::No | QMessageBox::Escape,
		   QMessageBox::NoButton );
    if(mb.exec() == QMessageBox::No)
      return;
    else 
      batchmode = true;
  }
  
  // check that TBName is already in the DB.

  QStringList files;
  QString localdir;
  if(!batchmode) {
    if(!tbnameInDB(namefilter_)) return;      
    
    // file selection
    QString file = 
      QFileDialog::getOpenFileName("Text (*.dat *.txt)",
				   "./",
				   this,
				   "open files dialog"
				   "Select one file to open");
    
    if(!file) return;
    files += file;
  } else {
    localdir = QFileDialog::getExistingDirectory("./",
						 this,
						 "get existing directory",
						 "Choose a directory",
						 true );  
    if(!localdir) return;
    
    // scanning directory for filenames containing filter
    QString filter(namefilter_); filter.replace(QRegExp("%"),"*");
    QDir seldir(localdir, filter); 
    QStringList matchingfiles = seldir.entryList();    
    files = matchingfiles;
  }
  
  bool prompt = true;
  for(QStringList::iterator i=files.begin(); i!=files.end(); i++) {

    QString file = *i;
    
    // detname is the 1st 8 characters of the file
    QString detname;
    if(batchmode) {
      detname = file; 
      detname.truncate(8);
      file.prepend(localdir);
    }
    else
      detname = namefilter_;
      
    string msg = "You are about to add the following file to the DataBase :\n";
    msg += "\t"; msg += file; msg += "\n";
    msg += "\n";
    msg += "TBName : \t";
    msg += detname;
    msg += "\n";  
    msg += "Start of validity :\t"; 
    msg += starttime_.toString(DtpImp::timeFormatSasha_); 
    msg += "\n";
    msg += "End of validity :\t"; 
    msg += endtime_.toString(DtpImp::timeFormatSasha_);
    msg += "\n";
    
    if(prompt) {
      // cancel is used to copy all files are copied without prompting ;)
      // no is used to copy all files are copied without prompting ;)

      QMessageBox mb("Dtp", 
		     msg.c_str(), 
		     QMessageBox::Information,
		     QMessageBox::Yes | QMessageBox::Default,
		     QMessageBox::No | QMessageBox::Escape,
		     QMessageBox::Cancel);
      mb.setButtonText(QMessageBox::Cancel,"All");
      mb.setButtonText(QMessageBox::No,"Next");
      
      switch (mb.exec()) {
      case QMessageBox::No :
	continue;
      case QMessageBox::Cancel :
	prompt=false;
      }
    }      

    // so far so good, let's do it !
    
    string errmsg;
    if(applyRangeCheck_ -> isChecked())
      errmsg = 
	db_->cdbUpload(
		       starttime_.toString(DtpImp::timeFormatSasha_).latin1(),
		       endtime_.toString(DtpImp::timeFormatSasha_).latin1(),
		       detname.latin1(),
		       file.latin1(),
		       locationCombo_->currentText().latin1(),
		       authorID_);
    else {
      cerr<<file.latin1()<<endl;
      errmsg = 
	db_->cdbUpload(file.latin1(),
		       locationCombo_->currentText().latin1(),
		       authorID_);
    }

    if(!errmsg.empty()) {
      QMessageBox mb("Dtp", 
		     errmsg.c_str(), 
		     QMessageBox::Critical,
		     QMessageBox::Ok | QMessageBox::Default,
		     QMessageBox::NoButton,
		     QMessageBox::NoButton );      
      mb.exec();
    }
  }
  selectDB();
}
  
void DtpImp::removeDB() {
  if(!db_ || !db_->isConnected() ) {
    NotConnected();
    return;
  }
  
  string dir(locationCombo_->currentText());  
  if(dbDirWritable(dir.c_str())) {
    
    int cs = dbtable_->currentSelection();
    if(cs == -1)
      return;

    bool prompt = true;
    QTableSelection selection = dbtable_->selection(cs);
    for(int i=selection.topRow(); i<=selection.bottomRow(); i++) {
 
      string start(dbtable_->text(i,0));
      string end(dbtable_->text(i,1));
      string filename(dbtable_->text(i,2));

      string owner = db_->fileOwner(start,end,filename,authorID_);
      if(owner.empty()) continue;
      if(owner != "myself") { 
	string msg = "Sorry Mr(Mrs) "; msg += lastName_;
	msg += " !\ncontact Mr(Mrs) ";
	msg += owner; msg += " if you want this file to be removed.";
	QMessageBox mb("Dtp", 
		       msg.c_str(), 
		       QMessageBox::Information,
		       QMessageBox::Ok | QMessageBox::Default,
		       QMessageBox::NoButton,
		       QMessageBox::NoButton);
	mb.exec();
	continue;
      }
      
      if(prompt) {
	// cancel is used to copy all files are copied without prompting ;)
	
	string msg = "Do you really want to remove from DB :\n";
	msg += filename; 
	msg += "\nStart of validity :\t"; msg += start; 
	msg += "\nEnd of validity :\t"; msg += end;
	msg += "\n";
	
	QMessageBox mb("Dtp", 
		       msg.c_str(), 
		       QMessageBox::Information,
		       QMessageBox::Yes | QMessageBox::Default,
		       QMessageBox::No | QMessageBox::Escape,
		       QMessageBox::Cancel);
	mb.setButtonText(QMessageBox::Cancel,"All");
	
	switch (mb.exec()) {
	case QMessageBox::No :
	  return ; 
	case QMessageBox::Cancel :
	  prompt=false;
	}
      }

      string qry = "delete from runlb.tb_calibration where StartTime = '";
      qry += start; qry += "' and EndTime = '";
      qry += end; qry += "' and FileName = '";
      qry += filename; qry += "'";
      
      string rm="rm "; rm += dir; rm += "/";
      rm += filename;
      
      if(db_->query(qry)) {
	system(rm.c_str());
	db_->endQuery();
      }      
    }
    selectDB();
  }
}

void  DtpImp::locationDB() const {

  if(!db_ || !db_->isConnected() ) {
    NotConnected();
    return;
  }
  
  string dirname;  
  string qry = "select Directory from runlb.tb_directories where runlb.tb_directories.key = 'calibration'";  
  
  db_->query(qry);
  do {
    locationCombo_->insertItem(db_->getCol(0));
  }  while(db_->getNextRow());
  db_->endQuery();
}

void DtpImp::NotConnected() const {
  QMessageBox mb("Dtp", 
		 "Not connected to DataBase.\nTo connect, use the Options tab.", 
		 QMessageBox::Information,
		 QMessageBox::Ok | QMessageBox::Default,
		 QMessageBox::NoButton,
		 QMessageBox::NoButton );
  mb.exec();
}


bool DtpImp::tbnameInDB(const char*) const {

  if(!db_) {
    NotConnected();
    return false;
  }

  string qry = "select * from runlb.tb_calibration where DetName like'"; 
  qry += namefilter_;
  qry += "'";
  db_->query(qry); 
  if(!db_->getNRows()) {
    string msg(namefilter_); msg += " is not yet registered in the database.\n Are you sure you want to continue ?\n";
    QMessageBox mb("Dtp", 
		   msg.c_str(), 
		   QMessageBox::Warning,
		   QMessageBox::Yes | QMessageBox::Default,
		   QMessageBox::No | QMessageBox::Escape,
		   QMessageBox::NoButton );
    if(mb.exec()==QMessageBox::No) return false;    
  }
  return true;
}
 
bool DtpImp::dbDirWritable(const char* dirname) const {

  string test = dirname; test += "/t";
  string touch = "touch "; touch += test;
  if(system(touch.c_str()))
    return false;
  else {
    string rm = "rm "; rm += test;
    system(rm.c_str());
    return true;
  }
}

bool DtpImp::string_match(const string &str, const string &pattern) {

    int    status;
    regex_t    re;

    if( regcomp(&re, pattern.c_str(), REG_EXTENDED|REG_NOSUB) != 0 )
        return false;

    status = regexec(&re, str.c_str(), (size_t) 0, NULL, 0);
    regfree(&re);

    if (status != 0)
        return false;

    return true;
}




