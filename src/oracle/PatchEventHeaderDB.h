#ifndef CS_PatchDB__include
#define CS_PatchDB__include

#include <pwd.h>
#include <iostream>
#include "CsEvent.h"
#include "CsStore.h"
#include "DaqDataDecoding/DaqEvent.h"

#include "occi.h"
#include "CsOraStore.h"

#define TRANSACTION_COMMIT true
#define TRANSACTION_ROLLBACK false

class PatchEventHeaderDB : public CsEndOfJob
{

public:
	static PatchEventHeaderDB* Instance(void);
	void DBPatchTriggerMaskForEvent(CsEvent *evt, 
													  CsStore *str, 
													  CS::DaqEvent *dev);
	bool end();	

protected:
	PatchEventHeaderDB();  // singleton constructor

private:
	//!! ~PatchEventHeaderDB(); // destructor (replaced by end())
	bool dbGetFileID(void);
	bool dbOpen(void);
	bool dbClose(bool bWithCommit);
	bool applyPatch(const CS::DaqEvent *dev, const unsigned int oldTM);
	const std::string here(const char *fnct);

private:
	static PatchEventHeaderDB* _instance; // singleton static attribute 

	oracle::occi::Environment* _oraEnv;
	oracle::occi::Connection*  _oraConn;
	unsigned _rawFileID;
	bool _bDBOpened;
	bool _vrbs;
	bool _autocommit;

};


#endif //CS_PatchDB__include

