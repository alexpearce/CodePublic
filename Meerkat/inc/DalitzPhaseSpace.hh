#ifndef DALITZ_PHASE_SPACE 
#define DALITZ_PHASE_SPACE

#include "AbsPhaseSpace.hh"

#include "TMath.h"

#include <vector>

/// Class that describes the Dalitz phase space: the phase space of a particle 
/// decaying to three other particles. This is a two-dimensional phase space 
/// defined by the squared invariant masses of the two pairs of final state particles. 

class DalitzPhaseSpace : public AbsPhaseSpace {

  public: 

    //! Constructor
    /*! 
      \param [in] phaseSpaceName name of the phase space
      \param [in] mD mass of the mother particle
      \param [in] mA mass of the 1st daughter
      \param [in] mB mass of the 2nd daughter
      \param [in] mC mass of the 3rd daughter
    */ 
    DalitzPhaseSpace(const char* phaseSpaceName, 
                     Double_t mD, 
                     Double_t mA, 
                     Double_t mB, 
                     Double_t mC); 

    //! Destructor
    virtual ~DalitzPhaseSpace(); 
    
    //! Get dimensionality of the phase space. Always equals to two. 
    /*! 
      \return 2 :) 
    */ 
    UInt_t dimensionality() { return 2; }

    //! Check if the point is within the phase space limits
    /*! 
      \param [in] x point 
      \return true in the point is within the phase space, false otherwise
    */ 
    Bool_t withinLimits(std::vector<Double_t> &x); 

    //! Get lower limit
    /*! 
      \param [in] var variable number (0 or 1 for this class)
      \return lower limit of the phase space for this variable
    */ 
    Double_t lowerLimit(UInt_t var); 
    
    //! Get upper limit
    /*! 
      \param [in] var variable number (0 or 1 for this class)
      \return upper limit of the phase space for this variable
    */ 
    Double_t upperLimit(UInt_t var);

    //! Return limits (lower and upper) for the variable at the certain point of the phase space
    /*! 
      \param [in] var number of the variable
      \param [in] x point in the phase space
      \param [out] lower lower limit
      \param [out] upper upper limit
    */ 
    Bool_t limits(UInt_t var, std::vector<Double_t> &x, Double_t* lower, Double_t* upper); 

  private: 
  
    //! Mass of the particle A
    Double_t m_a;

    //! Mass of the particle B
    Double_t m_b;

    //! Mass of the particle C
    Double_t m_c;

    //! Mass of the mother particle
    Double_t m_d;

    //! Squared mass of the particle A
    Double_t m_a2;

    //! Squared mass of the particle B
    Double_t m_b2;

    //! Squared mass of the particle C
    Double_t m_c2;

    //! Squared mass of the mother particle
    Double_t m_d2;

    //! Sum of squared masses of daughter and mother particles
    Double_t m_SqSum;
    
    //! Lower limit of AB invariant mass values
    Double_t m_MinAB;

    //! Upper limit of AB invariant mass values
    Double_t m_MaxAB;

    //! Lower limit of BC invariant mass values
    Double_t m_MinBC;

    //! Upper limit of BC invariant mass values
    Double_t m_MaxBC;

}; 

#endif
