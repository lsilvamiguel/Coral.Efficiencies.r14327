#ifndef StoreDeltaADC___include
#define StoreDeltaADC___include

#include <string>
#include <map>
#include <vector>

class StoreDeltaADC
{
  public:
              StoreDeltaADC ( void  );
              int     InputDeltaADC(const std::string &s);
              int    OutputDeltaADC( std::string &s ) const;
              int    GetDelta ( const std::pair<size_t,size_t> &src_id ) const;
              bool  UpdateDelta ( const std::pair<size_t,size_t> &src_id, size_t delta );

   static   std::pair <int,int>  MiniMaxSADC(  const std::vector<unsigned short>& samples );

  private:
    std::map <  std::pair<size_t,size_t>,  size_t >  srcid_delta_;
};

#endif // StoreDeltaADC___include
