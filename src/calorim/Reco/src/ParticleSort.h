#ifndef ParticleSort_h_
#define ParticleSort_h_

#include <vector>
#include <list>

// this function is used to compare particles
// it compares particles by ID, X, Y, Z and E
bool partcmp(const Reco::CalorimeterParticle &a,
             const Reco::CalorimeterParticle &b) {

  if ( a.GetID() == b.GetID() ) {
    if (   a.GetX() == b.GetX() ) {
      if (   a.GetY() == b.GetY() ) {
         if (   a.GetZ() == b.GetZ() ) {
            if (   a.GetE() == b.GetE() ) {
              return false;
            } else
              return ( (a.GetE() - b.GetE() ) < 0 );
         } else
           return ( (a.GetZ() - b.GetZ() ) < 0 );
      } else
        return ( (a.GetY() - b.GetY() ) < 0 );
    } else
      return ( (a.GetX() - b.GetX() ) < 0 );
  } else
    return ( (a.GetID() - b.GetID() ) < 0 );


};

// This struct is there to make a transistion for
// the std::list<> sort function
struct PartCmp {
  bool operator()(const Reco::CalorimeterParticle &a,
                  const Reco::CalorimeterParticle &b) {
    return partcmp( a, b);
  }
};


// SortPart expects a reference to a vector of particles, which
// is sorted afterwards, this implementation needs a lot of malloc
// calls but it respects the structures of the given code pool
void SortPart(std::vector<Reco::CalorimeterParticle> &vect, std::ostream *file = NULL) {

  // Allocate a std::list of particles, which can be sorted
  std::list<Reco::CalorimeterParticle> list;
  // and the compare struct to sort it
  struct PartCmp p;

  // Fill all elements of the vector into the list
  for( std::vector<Reco::CalorimeterParticle>::iterator it = vect.begin();
       it != vect.end(); it++)
    list.push_back(*it);

  // Sort the list
  list.sort(p);

  // Empty the vector
  vect.clear();

  // Fill the vector again
  for( std::list<Reco::CalorimeterParticle>::iterator it = list.begin();
       it != list.end(); it++) {
    vect.push_back(*it);
    if( file )
      file->write((char*)&(*it), sizeof(Reco::CalorimeterParticle) );
  }
}

#endif
