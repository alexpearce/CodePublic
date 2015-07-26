#include <stdio.h>
#include <vector>

#include "AbsPhaseSpace.hh"

AbsPhaseSpace::AbsPhaseSpace(const char* phaseSpaceName) {
  strncpy(m_name, phaseSpaceName, 255); 
  m_name[255] = 0; 
}

AbsPhaseSpace::~AbsPhaseSpace() {

}
