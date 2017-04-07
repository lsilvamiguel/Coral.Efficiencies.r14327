#ifndef __RUN_INFO_H__
#define __RUN_INFO_H__

#include <string>

struct RunInfo
{
  int    runNumber;
  int    year;
  std::string period;
  std::string database;
  std::string username;
  std::string password;
public:
  RunInfo(int _runNumber, bool isTestDB = false, bool isDstProduction = false);
  // Parameter isTestDB for compatibility with previous versions
  const char* GetRunPeriod() const { return period.c_str(); }
  int GetRunYear() const { return year; }
  static int GetFirstRunOfPeriod(const std::string& period); // Only for using in production utilities (not for CORAL)!!!
};

#endif
