#ifndef MESSAGES
#define MESSAGES

//! Reset time counter
void set_timer(void); 

//! Return 1 if more than a certain number of seconds have passed since the last call of this function or the last timer reset. 
/*!
  \param [in] secs number of seconds
  \return 0/1
*/ 
int timer(int secs);  

#endif

