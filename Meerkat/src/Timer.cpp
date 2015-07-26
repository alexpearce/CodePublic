#include "Timer.hh"

#include <time.h>

static int last_timer; 

void set_timer(void) {
  last_timer = time(0);
}

int timer(int secs) {
  int t = time(0); 
  if ( t - last_timer > secs ) {
    last_timer = t; 
    return 1; 
  } else {
    return 0; 
  }
}
