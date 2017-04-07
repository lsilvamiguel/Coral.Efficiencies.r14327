#ifndef ObjectsCounter___include
#define ObjectsCounter___include

#include <typeinfo>
#include <map>
#include <iostream>
#include <string>
#include <cstdlib>

class ObjectsCounterMaster
{
  // Methods
  public:
    static void         SetStream               (std::ostream *o);
    static void         Print                   (void);
  private:
    static bool         Init                    (void);
    static void         End                     (void);

  // Attributes
  public:
    static std::map<std::string,int> counter;
  private:
    static std::ostream *stream;
    static bool init;
};


template<class T>
class ObjectsCounter
{
  public:
                        ObjectsCounter          (const ObjectsCounter &) {ObjectsCounterMaster::counter[typeid(T).name()]++;}
                        ObjectsCounter          (void) {ObjectsCounterMaster::counter[typeid(T).name()]++;}
                       ~ObjectsCounter          (void) {ObjectsCounterMaster::counter[typeid(T).name()]--;}
};

#endif // ObjectsCounter___include
