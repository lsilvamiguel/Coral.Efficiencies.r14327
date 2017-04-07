char *MainFrame::HelpGuiTxt="This Document provides a short description of MurphyTV.\n"
"\n"
"1. How to connect to a data source \n"
"2. How to display the Errors\n"
"3. Command line options\n"
"4. How to contact the author\n"
"\n"
"\n"
"1. How to connect to a data source \n"
"\n"
"In the main window (the little one) you can directly enter a source in the\n"
"\"Data Source\" text entry field.\n"
"\n"
"If you want to look at a local file, simply enter the name here\n"
"(or click on [...] for a fileselection dialog).\n"
"For remote files (RFIO) enter Filename@Hostname \n"
"To connect to the online data stream, enter @Hostname:\n"
"After you entered the data source, press ENTER or click on [Connect]\n"
"to start analysing the data.\n"
"The status bar below should now say \"Connected\" (for online data, it can take\n"
"a few seconds)\n"
"\n"
"You can connect at any time to a new data source.\n"
"If you do this, all the informations gathered from the old data are lost.\n"
"\n"
"You can click on [Disconnect] to stop analysing data.\n"
"The informations gathered so far are not lost in this case. \n"
"\n"
"You can enter a mapping file (or directory with mapping files) in the \"Mapping File\"\n"
"text entry field. Some detectors (ECAL ADC, MasterTime) data cannot be \n"
"decoded without mapping file.\n"
"\n"
"\n"
"2. How to display the Errors\n"
"\n"
"In the main window, click on [Display]. Another window pops up.\n"
"(The Display window does not affect the data analysis in any way.\n"
"You can open and close it at any time.)\n"
"\n"
"There are 6 different sections :\n"
"\n"
"2.1 The \"Catches\" section\n"
"\n"
"This is an overview over all SourceIDs in the data, as well as about the Ports of a selected\n"
"SourceID.\n"
"\n"
"The upper table contains for each SourceID the following informations:\n"
"\n"
"-Type: The abbreviations taken from [Compass Notes 8-2001] for the different detectors,\n"
" for which different SourceID ranges are reserved. \n"
" \n"
"-BadEvnt: The number of events which had one or more errors.\n"
"\n"
"-#Errors: The number of errors on this SourceID.  \n"
"\n"
"-#Headers: The number of header words coming from the frontends of the SourceID.\n"
" This does NOT include the SLink Header\n"
" (For GEMs, all words are currently counted as data, so this should be zero)\n"
" \n"
"-#Data: The number of data words coming from the frontends of the SourceID.\n"
"\n"
"-Special: The number of errors of the type you can select in the \"Special\" combo box.\n"
" Use this combo box if you want to know which SourceIDs are producing the error you\n"
" are interested in.\n"
"\n"
"-eg. at: The first trigger number where the error selected with the \"Special\" combo box\n"
" occured. This is useful if you want to see the error with eventDumpAll in the data\n"
" yourself.\n"
"\n"
"-#Spills: The number of spilla affected by this error\n"
"\n"
"Last Spill: The last spill affected by the error\n"
"\n"
"\n"
"The number of errors,header and data words can be normalized to the number of recorded\n"
"events with the \"Normalize\" check box. \n"
"You can chose how to sort the list by clicking on the buttons at the top of the table.\n"
"\n"
"If you double-click a line ,the \"Errors\" Tab is displayed with this sourceID selected.\n"
"Do this to get a list of all error types that occured on a sourceID.\n"
"\n"
"If you click on a SourceID, in the lower table for each Port or GeographicID on\n"
"that SourceID the number of errors,data and header words are displayed.\n"
"\n"
"Some informations are associated only to a port number in the data, others to a\n"
"GeoID, and some to both. Latter case is used to associate also the geoId based\n"
"informations to port numbers, and geoIDs to port numbers.\n"
"For GEMs or the Mastertime the Value of the error word is displayed in the \n"
"GeoID/Value field.\n"
"\n"
"In the \"Attached Detectors\" window at the bottom, the names of the detectors connected\n"
"to the selected sourceID are displayed. This information comes from the mapping file(s).\n"
"If no mapping file(s) are selected, no names are displayed.\n"
"\n"
"\n"
"2.2 The \"Errors\" section\n"
"\n"
"Here you have a table in which for each kind of error it is displayed:\n"
"\n"
"-How often it occured.\n"
"\n"
"-How often it occured on the selected SourceID. You can enter this sourceID in the\n"
" \"SrcID=\" text entry field.\n"
" \n"
"-How many of the events are infected with this error.\n"
"\n"
"-The same for the selected sourceID. \n"
"\n"
"-The number of sourceId's on which the error occured.\n"
"\n"
"The counts are normed to the number of recorded events until you uncheck the\n"
"\"Normalize\" check box.\n"
"\n"
"If you check the \"Show only errors on SrcID\" checkbox, only errors are shown which\n"
"appear on the selected sourceID.\n"
"If you enter a sourceID, this option is automatically selected.\n"
"\n"
"If you click on an error, a short description is displayed at the bottom.\n"
"\n"
"If you double-click on an error, the \"Catches\" section is displayed with the error\n"
"selected as the special interest. Do this to learn on which sourceIDs the error occures. \n"
"\n"
"\n"
"2.3 The \"Monitor\" section\n"
"\n"
"On the top, you see some general informations:\n"
"\n"
"-The run number,as well as the date and time and Trigger Number of the event which \n"
" is currently analysed. \n"
"\n"
"-The number of events analysed.\n"
"\n"
"-The total number of errors found in the analysed events, as well as this number per\n"
" analysed event.\n"
" \n"
"-The number of different SourceIDs (e.g. CATCH boards) found in the data\n"
"\n"
"-The number of SourceIDs which appeared in some events (or are in the list mentioned),\n"
" but not in all events.\n"
" \n"
"-The number of SourceIDs on which not a single error was found.\n"
"\n"
"-The avarage event size (including all kind of headers )\n"
"\n"
"Below , you see a message window. If an error appears on a sourceID ,or if the avarage \n"
"number of events with an error on a sourceID changes dramatically, a short message is \n"
"displayed here. \n"
"\n"
"At the bottom, there is an alarm text. It turns red if new error messages appear.\n"
"If you have read them, click on the button at the right to turn the alarm off till\n"
"there are new error messages again.\n"
"\n"
"To get a description of the different errors, go to the \"Errors\" section.  \n"
"\n"
"\n"
"2.4 The \"History\" section\n"
"\n"
"Here you can see the number of errors plotted versus trigger number.\n"
"Number of errors means the number of recorded errors which fall in one 5000 trigger-\n"
"number bin normalised to the number of recorded events in that bin.\n"
"\n"
"You can display all errors, all errors of a given category or single errors.\n"
"Select this with the combo boxes on the bottom.\n"
"\n"
"\n"
"2.5 The \"Event Sizes\" section\n"
"\n"
"In this histogram the distribution of the total event size (including headers) coming\n"
"from a given SourceID is displayed.\n"
"The SourceID can be selected in the \"SourceID\" text entry.\n"
"\n"
"2.6 The \"Data Flow\" section\n"
"\n"
"Here the total event size per recorded event is plotted versus the trigger number.\n"
"The histogram's bin size is 5000 Trigger numbers.\n"
"\n"
"\n"
"3. Command line options\n"
"\n"
"MurphyTV [-m mapping_file] [-d] [-f] [-C] Data_Source\n"
"\n"
"To open the display window at start use -d\n"
"To connect to a data source at startup use MurphyTV Data_Source \n"
"To specify a mapping file use -m mapping_file\n"
"To use a small display window -f \n"
"To read only events with calibration triggers -C\n"
"Default, i.e. without any additional option, is to read only physics triggers.\n"
"\n"
"\n"
"4. How to contact the author\n"
"\n"
"Send bug reports, critics and suggestions to :\n"
"  \n"
"Christian.Schill@cern.ch\n"
"\n";
