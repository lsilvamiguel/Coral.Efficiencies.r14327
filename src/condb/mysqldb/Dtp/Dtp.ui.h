/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you wish to add, delete or rename slots use Qt Designer which will
** update this file, preserving your code. Create an init() slot in place of
** a constructor, and a destroy() slot in place of a destructor.
*****************************************************************************/

void Dtp::init() {
    namefilter_ = "%";
    username_ = "toto";
    password_ = "toto"; 
    userEdit_ -> setText(username_);
    passwordEdit_->setText(password_);
    applyRangeCheck_ -> setChecked (true);
}



void Dtp::fileNew()
{

}

void Dtp::fileOpen()
{

}

void Dtp::fileSave()
{

}

void Dtp::fileSaveAs()
{

}

void Dtp::filePrint()
{

}

void Dtp::fileExit()
{

}

void Dtp::editUndo()
{

}

void Dtp::editRedo()
{

}

void Dtp::editCut()
{

}

void Dtp::editPaste()
{

}

void Dtp::editFind()
{

}

void Dtp::helpIndex()
{

}

void Dtp::helpContents()
{

}

void Dtp::helpAbout()
{

}





void Dtp::setNameFilter( const QString & namefilter )
{
    namefilter_ = namefilter;
}

void Dtp::setUsername( const QString & username )
{
    username_ = username;
}

void Dtp::setPassword( const QString & password )
{
    password_ =password;
}

void Dtp::setStartTime( const QDateTime & datetime )
{
    starttime_ = datetime;
}

void Dtp::setEndTime( const QDateTime & datetime )
{
    endtime_=datetime;
}

int Dtp::selectDB()
{
    return 0;
}

void Dtp::downloadDB()
{
    cout<<"download from DB not implemented yet"<<endl;
}

void Dtp::uploadDB()
{
   cout<<"upload to DB not implemented yet"<<endl;
}



void Dtp::removeDB()
{

}

void Dtp::applyDB()
{

}

bool Dtp::connectDB()
{

}

