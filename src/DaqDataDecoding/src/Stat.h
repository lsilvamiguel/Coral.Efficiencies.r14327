#ifndef CS__Stat____include
#define CS__Stat____include

#include <set>

namespace CS {

class Stat
{
  public:
  
    virtual            ~Stat                    (void) {}
                        Stat                    (void) : fills_counter(0),buffer_max_size(1000),cut_ratio(0.5) {}
    
    void                Add                     (float v);

    void                SetBufferMaxSize        (unsigned s);
    unsigned            GetBufferMaxSize        (void) const {return buffer_max_size;}

    unsigned            Size                    (void) const {return buffer.size();}

    float               Mean                    (float events_fraction=1) const;
    float               Sigma                   (float events_fraction=1) const;
    
    void                SetCutRatio             (float r) {cut_ratio=r;}
    float               GetCutRatio             (void) const {return cut_ratio;}
    
    const
    std::multiset<float> &   GetBuffer               (void) const {return buffer;}
    std::multiset<float> &   GetBuffer               (void)       {return buffer;}

    unsigned            GetFillsCounter         (void) const {return fills_counter;}

  private:
  
    unsigned            fills_counter;
  
    unsigned            buffer_max_size;
    
    std::multiset<float>buffer;

    // Once buffer size reaches buffer_max_size,
    // the cut_ratio of events will be removed from the buffer
    float               cut_ratio;
};

} // namespace CS

#endif // CS__Stat____include
