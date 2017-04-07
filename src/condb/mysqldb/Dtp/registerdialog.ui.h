/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you wish to add, delete or rename slots use Qt Designer which will
** update this file, preserving your code. Create an init() slot in place of
** a constructor, and a destroy() slot in place of a destructor.
*****************************************************************************/


void RegisterDialog::setEmail(const QString & email)
{
    email_ = email;
}

void RegisterDialog::setFirstName(const QString & firstname)
{
    firstName_ = firstname;
}

void RegisterDialog::setLastName(const QString & lastname)
{
    lastName_ = lastname;
}

void RegisterDialog::parameters(QStringList & parameters)
{
    parameters += login_;
    parameters += firstName_;
    parameters += lastName_;
    parameters += email_;   
}

void RegisterDialog::acceptChecked()
{
    if(login_ && firstName_ && lastName_ && email_ && email_.contains("@"))
	accept();
}

void RegisterDialog::setLogin( const QString & login )
{
    login_=login;
}