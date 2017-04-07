

#ifndef __VariousSettings_
#define __VariousSettings_

#include <TROOT.h>
#include <TList.h>
#include <vector>

#include "TPostScript.h"
#include "Monitor.h"
#include <iostream>


class VariousSettings  : public TObject  {

  private:

    /// Config file name
    std::string             fConfigFileName;
    /// working directory
    std::string             fPwd;

    /// list of last config files
    std::vector<std::string>       fConfigFileNameList;

    /// Monitoring object
    Monitor *fMonitor;

    /// filename for personnal options
    static const char* fPersoFile;

    /// maximum number of items in the config files list
    static int fMaxItemInList;


    /// add a name to the last config files list
    void AddNameToList(std::string name);
    void AddNameToList(const char* name)
          { AddNameToList(std::string(name)); };

    /// load one .cfg (or .xml) file
    Bool_t LoadConfigFile(std::string filename, int type_read, Bool_t save_name);

    /// filter function for directory scaning of .cfg files
//    int FilterDotcfg (const struct dirent *);

    /// comparison function to sort directory using inodes numbers
//    int SortByInode (const struct dirent **, const struct dirent **);


  public:

    VariousSettings(Monitor* monitor)
        { fConfigFileName = "";
          fMonitor = monitor;
          fPwd = gSystem->pwd();
        }

    ~VariousSettings()   { }

    /// Load configuration file
    Bool_t LoadConfigSettings(std::string filename, int type_read = 0, Bool_t save_name = true);
    Bool_t LoadConfigSettings(const char* filename, int type_read = 0, Bool_t save_name = true)
        { return LoadConfigSettings(std::string(filename), type_read, save_name); }
    Bool_t LoadConfigSettings()
        { if (fConfigFileName != "") return LoadConfigSettings(fConfigFileName);
            else { std::cerr << "LoadConfigSettings: No default config file name\n"; return false; }
        }
    Bool_t LoadDefaultSettings(std::string filename)
        { return LoadConfigSettings(filename.c_str(), 2, false); }
    // type_read: 0 -> read gui + plane&group settings
    //            1 -> read only gui settings
    //            2 -> read only plane&group settings
    // save_name: true -> fConfigFileName is set with filename, and filename is added in fConfigFileNameList




    /// Save configuration file
    Bool_t SaveConfigSettings(const char* filename, int type_save = 0);
    Bool_t SaveConfigSettings(std::string filename, int type_save = 0)
        { return SaveConfigSettings(filename.c_str(), type_save); }
    Bool_t SaveConfigSettings()
        { if (fConfigFileName != "") return SaveConfigSettings(fConfigFileName);
            else { std::cerr << "SaveConfigSettings: No default config file name\n"; return false; }
        }
    // type_save: 0 -> save gui + plane&group settings (fConfigFileName variable is set with filename)
    //            1 -> save only gui settings
    //            2 -> save only plane&group settings

    /// Fill a TPostScript object with a gui config file
    Bool_t GuiSettingsToPS(const char* filename, TPostScript* psfile);

    /// Create a PS file from a gui settings file
    Bool_t GuiSettingsCreatePS(std::string filename, const TList* filelist = 0);
    Bool_t GuiSettingsCreatePS(const char* filename, const TList* filelist = 0)
         { return GuiSettingsCreatePS(std::string(filename), filelist); }
    Bool_t GuiSettingsCreateOrderlyShiftPlots();

    /// read personal parameters file (.cooolrc)
    bool ReadPersoParam(const char * filename);
    bool ReadPersoParam();

    /// write personal parameters file (.cooolrc)
    bool WritePersoParam(const char * filename);
    bool WritePersoParam()
          { return WritePersoParam(fPersoFile); }


    /// return actual config file name
    const std::string&  GetConfigFileName()  const { return fConfigFileName; }

    /// return actual config file name
    const std::vector<std::string>&  GetConfigFileNameList()  const { return fConfigFileNameList; }

    /// return if fConfigFileName is not empty
    Bool_t HasConfigFileName() { return (fConfigFileName != ""); }


  ClassDef(VariousSettings,0)
};


#endif
