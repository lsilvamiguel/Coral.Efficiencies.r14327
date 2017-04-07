#ifndef Interval_h
#define Interval_h

class Interval{
public:
  Interval(): since_(0,0),till_(0,0){}
  ~Interval(){}
  Interval(CDB::Time since,CDB::Time till) : since_(since),till_(till){}
  CDB::Time getSince(){return since_;}
  CDB::Time getTill(){return till_;}
protected:
  CDB::Time since_;
  CDB::Time till_;
};
#endif

