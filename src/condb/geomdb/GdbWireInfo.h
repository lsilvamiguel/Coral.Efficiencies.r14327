//-*-Mode: C++;-*-
#ifndef _GdbWireInfo_h_
#define _GdbWireInfo_h_
/*!
  \file		GdbWireInfo.h
  \brief	Geometry DB wire information object definition file
  \author	$Author: miyachi $
  \version	$Revision: 1.2 $
  \date		$Date: 2000/06/07 08:36:48 $
*/

#if __GNUG__ >= 2
#  pragma interface
#endif

//---- GdbWireInfo -----------------------------------------------------------

class GdbWireInfo {
public:
	GdbWireInfo();

	GdbWireInfo(	const int& id,
					const int& number, 
					const double& pitch,
					const double& angle, 
					const double& firstPosition,
					const int& flag=0);

   GdbWireInfo(	const int& number, 
					const double& pitch,
					const double& angle, 
					const double& firstPosition,
					const int& flag=0);
					
    GdbWireInfo(const GdbWireInfo& wInfo);
    
    ~GdbWireInfo();

	GdbWireInfo& operator=(const GdbWireInfo& wInfo);
	
	int id() const ;
	int number() const ;
	double pitch() const;
	double angle() const;
	double firstPosition() const;
	int flag() const;

	void dumpDetectorsTable(ostream& os);

	friend ostream& operator<<(ostream& os, const GdbWireInfo& wire);
	friend istream& operator>>(istream& is, GdbWireInfo& wire);
	void id(const int& id);

protected:
	void number(const int& number);
	void pitch(const double& pitch);
	void angle(const double& angle);
	void firstPosition(const double& firstPosition);
	void flag(const int& flag);

private:

	int id_;
	int number_;
	double pitch_;
	double angle_;
	double firstPosition_;
	int flag_;
};

#endif
