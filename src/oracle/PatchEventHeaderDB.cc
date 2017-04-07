
#include "PatchEventHeaderDB.h"
#include "CsRegistry.h"
#include "CsRegistrySing.h"
#include "RunInfo.h"

using namespace std;
using namespace oracle::occi;

PatchEventHeaderDB* PatchEventHeaderDB::_instance=0;

PatchEventHeaderDB::PatchEventHeaderDB(void)
{
	_vrbs = ( getenv("PATCH_EVENT_HEADER_VERBOSE") );
	if ( (_bDBOpened = dbOpen()) ) {
		if (!dbGetFileID()) {
			CsErrLog::Instance()->mes( elFatal, "Could not retrieve file ID");
		}
	}
}

PatchEventHeaderDB* PatchEventHeaderDB::Instance()
{
	if( _instance == 0 ) {
		_instance = new PatchEventHeaderDB();
		CsRegistry reg;
		reg.EOJRegistration(_instance);
	}
	return _instance;
}

void PatchEventHeaderDB::DBPatchTriggerMaskForEvent(CsEvent *evt, CsStore *str, CS::DaqEvent *dev) 
{


	// 
	// this piece of code is needed to alter DB (besides evt, str & dev)
	//
	if (_bDBOpened) {
		if (str->getTriggerMask() != dev->GetHeader().GetTrigger()) {
			if ( !applyPatch(dev, str->getTriggerMask()) ) {
				cerr << here(__FUNCTION__) << ". ERROR. Patching triggerMask failed" << endl;
			}
		} else {
			_vrbs && cout << here(__FUNCTION__) << ". No need to update TRIGGER_MASK=" << dev->GetHeader().GetTrigger() << " for RAWEV_FILE_ID=" << _rawFileID << " and EVENT_NUMBER=" << dev->GetEventNumberInRun() << " and BURST=" << dev->GetBurstNumber() << " and EVENT_IN_BURST=" << dev->GetEventNumberInBurst() << endl;
		}
	}
	//
	// END
	//
}

bool PatchEventHeaderDB::dbGetFileID(void)
{
	CsEvent* evt( CsEvent::Instance() );
	CsStore* str( evt->GetStore() );
	
	string sql("select unique file_id from file_maps where run_number=:1 and file_name=:2");
	Statement* stmt(0);
	ResultSet* rs(0);
	try {
		stmt = _oraConn->createStatement();
		stmt->setAutoCommit(false);
		stmt->setSQL(sql);
		stmt->setUInt(1, str->getRun()); 
		stmt->setString(2, CsOraStore::Instance()->GetCurrentChunkName());
		rs = stmt->executeQuery();
	} catch(SQLException& e) {
		cerr << here(__FUNCTION__) << "CoralUserInit/executing \"" 
			  << "CoralUserInit/executing \"" << sql << "\" with run_number=" 
			  << str->getRun() << " and file_name=" 
			  << CsOraStore::Instance()->GetCurrentChunkName() 
			  << ". ERROR " << e.getMessage() << endl;
		return false;
	} catch(...) {
		cerr << here(__FUNCTION__) << ". Unknown exception" << endl;
		return false;
	}
	if ( rs==0 || !rs->next() ) {
		cerr << here(__FUNCTION__) << ". ERROR. No rawFileID found" << endl;
		stmt->closeResultSet(rs);
		_oraConn->terminateStatement(stmt);
		return false;
	} else {
		do {
			_rawFileID = rs->getUInt(1);
		} while(rs->next());
	}
	stmt->closeResultSet(rs);
	_oraConn->terminateStatement(stmt);
	stmt=0;
	rs=0;
	return true;
}

bool PatchEventHeaderDB::dbOpen()
{
	ostringstream user;
	string pwd("cdrprod");
	string database("compr");

	// allow overriding hardcoded pwd and db name
	if ( getenv("ORACLE_HOOKUPDB") ) database = getenv("ORACLE_HOOKUPDB");
	if ( getenv("COMPASS_DB_PWD") ) database = getenv("COMPASS_DB_PWD");	
	if ( getenv("COMMIT_TM_PATCH") ) {
		_vrbs && cout << here(__FUNCTION__) << ". COMMIT of transaction is ENFORCED." << endl;
		_autocommit = ( getenv("TM_PATCH_AUTOCOMMIT") != 0 );
		_vrbs && cout << here(__FUNCTION__) << ". Setting autocommit=" << _autocommit << " for update statements." << endl;
	} else {
		_autocommit = false;
	}
 
	RunInfo runInfo(CsEvent::Instance()->getRunNumber());
	string runPeriod(runInfo.GetRunPeriod());
	int runYear(runInfo.GetRunYear());

	if (runYear != 2006) return false; // only 2006 TM requires the patch

	passwd *etcPasswd;
	if ( (etcPasswd = getpwuid(geteuid())) ) {
      if (etcPasswd->pw_name && 
			 (strncmp(etcPasswd->pw_name,"objsrvv", 7) == 0 ||
			  strncmp(etcPasswd->pw_name,"na58dst", 7) == 0 )
			 ) {
			// OK: only objsrvvy && na58dstX will get the DB open 
			_vrbs && cout << here(__FUNCTION__) << ". Granting access to TriggerMask DB patching to user " << etcPasswd->pw_name << endl;
		} else {
			_vrbs && cout << here(__FUNCTION__) << ". User " << etcPasswd->pw_name << " is not allowed to alter the TriggerMask in the DB" << endl;
			return false; 
		}
	} else {
      cerr << this << here(__FUNCTION__) << ". ERROR. could not determine user name. Cannot allow current program to alter TRIGGER_MASK in the DB" << endl;
		return false;
	}

	// build user's name: COMPASS_YYPPP
	user << "COMPASS_" << setw(2)<<setfill('0') << (runYear>2000 ? runYear-2000 : runYear) << runPeriod;

	_vrbs && cout << here(__FUNCTION__) << ". Opening " << user.str() << "/" << pwd << "@" << database << endl;

	try {
		_oraEnv = Environment::createEnvironment(Environment::OBJECT);
		if(!_oraEnv) {
			cerr << here(__FUNCTION__) << ". Cannot create Oracle environment." << endl;
			return false;
		}
		_oraConn = _oraEnv->createConnection( user.str(), pwd, database );
		if(!_oraConn) {
			cerr << here(__FUNCTION__) << ". Cannot create Oracle connection." << endl;
			return false;
		}
	} catch(SQLException& e) {
      cerr << here(__FUNCTION__) << ". " << e.getMessage() << endl;
		return false;
	} catch(...) {
      cerr << here(__FUNCTION__) << " Unknown exception." << endl;
		return false;
	}

	return true;
}

bool PatchEventHeaderDB::dbClose(bool bWithCommit)
{
	try {
		if (bWithCommit) {
			_vrbs && cout << here(__FUNCTION__) << ". Committing transaction." << endl;
			_oraConn->commit();
		} else {
			_vrbs && cout << here(__FUNCTION__) << ". Rolling-back transaction." << endl;
			_oraConn->rollback();
		}
	} catch(SQLException& e) {
      cerr << here(__FUNCTION__) << ". " << e.getMessage() << endl;
		return false;
	} catch(...) {
      cerr << here(__FUNCTION__) << ". Unknown exception." << endl;
		return false;
	}

	try {
      if(_oraConn && _oraEnv) _oraEnv->terminateConnection(_oraConn); 
      _oraConn = 0;
      if(_oraEnv)	Environment::terminateEnvironment(_oraEnv);
      _oraEnv = 0;
	} catch(SQLException& e) {
		cerr << here(__FUNCTION__) << ". " << e.getMessage() << endl;
		return false;
	} catch(...) {
      cerr << here(__FUNCTION__) << ". Unknown exception." << endl;
		return false;
	}
	_bDBOpened = false;
	return true;
}

bool PatchEventHeaderDB::applyPatch(const CS::DaqEvent *dev, const unsigned int oldTM)
{
	bool retVal(true);
	
	string sql("update event_headers set trigger_mask=:1 where RAWEV_FILE_ID=:2 and EVENT_NUMBER=:3 and BURST=:4 and EVENT_IN_BURST=:5");
	
	_vrbs && cout << here(__FUNCTION__) << "/applying \"" << sql << "\" with TRIGGER_MASK=" << dev->GetHeader().GetTrigger() << " and RAWEV_FILE_ID=" << _rawFileID << " and EVENT_NUMBER=" << dev->GetEventNumberInRun() << " and BURST=" << dev->GetBurstNumber() << " and EVENT_IN_BURST=" << dev->GetEventNumberInBurst() << "; (old TM=" << oldTM << ")" << endl;
	
	Statement* stmt(0);
	ResultSet* rs(0);
	try {
		stmt = _oraConn->createStatement();
		stmt->setAutoCommit(_autocommit);
		stmt->setSQL(sql);
		stmt->setUInt(1, dev->GetHeader().GetTrigger()); // TRIGGER_MASK
		stmt->setUInt(2, _rawFileID); // RAWEV_FILE_ID
		stmt->setUInt(3, dev->GetEventNumberInRun()); // EVENT_NUMBER
		stmt->setUInt(4, dev->GetBurstNumber()); // BURST_NUMBER
		stmt->setUInt(5, dev->GetEventNumberInBurst()); // EVENT_IN_BURST
		rs = stmt->executeQuery();
		_vrbs && cout << here(__FUNCTION__) << ". TriggerMask changed to " << dev->GetHeader().GetTrigger() << " <<<<" << endl;
	} catch(SQLException& e) {
		cerr <<  here(__FUNCTION__) << "/executing \"" << sql << "\" with RAWEV_FILE_ID=" << _rawFileID << " and EVENT_NUMBER=" << dev->GetEventNumberInRun() << " and BURST=" << dev->GetBurstNumber() << " and EVENT_IN_BURST=" << dev->GetEventNumberInBurst() << ". ERROR " << e.getMessage() << endl;
		retVal = false;
	} catch(...) {
		cerr << here(__FUNCTION__) << "/executeQuery: Unknown exception." << endl;
		retVal = false;
	}
	
	stmt->closeResultSet(rs);
	_oraConn->terminateStatement(stmt);
	stmt = 0;
	rs = 0;

	return retVal;	
}

const string PatchEventHeaderDB::here(const char *fnct)
{
	string str(fnct);
	return "PatchEventHeaderDB::" + str;
}

bool PatchEventHeaderDB::end() 
{
	if ( _bDBOpened ) {
		if ( getenv("COMMIT_TM_PATCH") ) {
			dbClose(TRANSACTION_COMMIT);
		} else {
			dbClose(TRANSACTION_ROLLBACK);
		}
	}
	delete _instance;
	_instance = 0;
	return true;
}
