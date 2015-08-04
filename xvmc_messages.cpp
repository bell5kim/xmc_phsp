/*****************************************************************************
 * xvmc_messages.cpp: message functions                                      *
 *                                                                           *
 * Copyright (C) 2000    Matthias Fippel                                     *
 *                       Abteilung fuer Medizinische Physik,                 *
 *                       Universitaetsklinikum Tuebingen, Germany            *
 *                                                                           *
 * revisions:                                                                *
 *    initial coding                                      MF 99/12/13        *
 *                                                                           *
 *****************************************************************************/

// ****************************************
// includes
// ****************************************
#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
using namespace std;

#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>

#include "definitions.h"
#include "global.h"
#include "ranmar.h"

// error handler
void xvmc_error(const char *location, const char *reason, const int exit_code)
{
   cerr << "XVMC>" << endl;
   cerr << "XVMC> RUN-TIME ERROR in " << location << endl;
   cerr << "XVMC> Reason: " << reason << "!" << endl;
   cerr << "XVMC> Exiting to system..." << endl;
   cerr << "XVMC>" << endl;
   exit(exit_code);
}

// warning handler (print location and reason)
void xvmc_warning(const char *location, const char *reason, const int mode)
{
   switch (mode)
   {
   case 1:
      cerr << "XVMC>" << endl;
      cerr << "XVMC> RUN-TIME WARNING in " << location << endl;
      cerr << "XVMC> Reason: " << reason << "!" << endl;
      break;
   case 2:
      cerr << "XVMC>" << endl;
      cerr << "XVMC> RUN-TIME WARNING in " << location << endl;
      cerr << "XVMC> Reason: " << reason << "!" << endl;
      cerr << "XVMC>" << endl;
      break;
   default:
      cerr << "XVMC> RUN-TIME WARNING in " << location << endl;
      cerr << "XVMC> Reason: " << reason << "!" << endl;
      break;
   }
   return;
}

// warning handler (print out bool variable, name and value)
void xvmc_warning(const char *name, const bool value, const int mode)
{
   char tf[6];

   if (value) strcpy(tf,"TRUE");
   else       strcpy(tf,"FALSE");

   switch (mode)
   {
   case 1:
      cerr << "XVMC> " << endl;
      cerr << "XVMC> " << "  " << name << ": " << tf << endl;
      break;
   case 2:
      cerr << "XVMC> " << endl;
      cerr << "XVMC> " << "  " << name << ": " << tf << endl;
      cerr << "XVMC> " << endl;
      break;
   default:
      cerr << "XVMC> " << "  " << name << ": " << tf << endl;
      break;
   }
   return;
}

// warning handler (print out integer variable, name and value)
void xvmc_warning(const char *name, const int value, const int mode)
{
   switch (mode)
   {
   case 1:
      cerr << "XVMC> " << endl;
      cerr << "XVMC> " << "  " << name << ": " << value << endl;
      break;
   case 2:
      cerr << "XVMC> " << endl;
      cerr << "XVMC> " << "  " << name << ": " << value << endl;
      cerr << "XVMC> " << endl;
      break;
   default:
      cerr << "XVMC> " << "  " << name << ": " << value << endl;
      break;
   }
   return;
}

// warning handler (print out real variable, name and value)
void xvmc_warning(const char *name, const real value, const int mode)
{
   switch (mode)
   {
   case 1:
      cerr << "XVMC> " << endl;
      cerr << "XVMC> " << "  " << name << ": " << value << endl;
      break;
   case 2:
      cerr << "XVMC> " << endl;
      cerr << "XVMC> " << "  " << name << ": " << value << endl;
      cerr << "XVMC> " << endl;
      break;
   default:
      cerr << "XVMC> " << "  " << name << ": " << value << endl;
      break;
   }
   return;
}

// message handler (1 message)
void xvmc_message(const char *message, const int mode)
{
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message << endl;
      break;
   case 2:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message << endl;
      cout << "XVMC> " << endl;
      break;
   default:
      cout << "XVMC> " << message << endl;
      break;
   }
   return;
}

// message handler (2 messages)
void xvmc_message(const char *message1, const char *message2, const int mode)
{
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << message2 << endl;
      break;
   case 2:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << message2 << endl;
      cout << "XVMC> " << endl;
      break;
   default:
      cout << "XVMC> " << message1 << " " << message2 << endl;
      break;
   }
   return;
}

// message handler (message1, integer number and message2)
void xvmc_message(const char *message1, const int number,
                  const char *message2, const int mode)
{
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      break;
   case 2:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      cout << "XVMC> " << endl;
      break;
   default:
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      break;
   }
   return;
}

// message handler (message1, long integer number and message2)
void xvmc_message(const char *message1, const long number,
                  const char *message2, const int mode)
{
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      break;
   case 2:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      cout << "XVMC> " << endl;
      break;
   default:
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      break;
   }
   return;
}

// message handler (message1, float number and message2)
void xvmc_message(const char *message1, const float number,
                  const char *message2, const int mode)
{
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      break;
   case 2:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      cout << "XVMC> " << endl;
      break;
   default:
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      break;
   }
   return;
}

// message handler (message1, double number and message2)
void xvmc_message(const char *message1, const double number,
                  const char *message2, const int mode)
{
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      break;
   case 2:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      cout << "XVMC> " << endl;
      break;
   default:
      cout << "XVMC> " << message1 << " " << number << " " << message2 << endl;
      break;
   }
   return;
}

// message handler
// (message1, float number1, message2, float number2 and message3)
void xvmc_message(const char *message1, const float number1,
                  const char *message2, const float number2,
                  const char *message3, const int mode)
{
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << endl;
      break;
   case 2:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << endl;
      cout << "XVMC> " << endl;
      break;
   default:
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << endl;
      break;
   }
   return;
}

// message handler
// (message1, double number1, message2, double number2 and message3)
void xvmc_message(const char *message1, const double number1,
                  const char *message2, const double number2,
                  const char *message3, const int mode)
{
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << endl;
      break;
   case 2:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << endl;
      cout << "XVMC> " << endl;
      break;
   default:
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << endl;
      break;
   }
   return;
}

// message handler
// (message1, float number1, message2, float number2,
//  message3, float number3  and message4)
void xvmc_message(const char *message1, const float number1,
                  const char *message2, const float number2,
                  const char *message3, const float number3,
                  const char *message4, const int mode)
{
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << " " << number3 << " "
                       << message4 << endl;
      break;
   case 2:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << " " << number3 << " "
                       << message4 << endl;
      cout << "XVMC> " << endl;
      break;
   default:
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << " " << number3 << " "
                       << message4 << endl;
      break;
   }
   return;
}

// message handler
// (message1, double number1, message2, double number2,
//  message3, double number3  and message4)
void xvmc_message(const char *message1, const double number1,
                  const char *message2, const double number2,
                  const char *message3, const double number3,
                  const char *message4, const int mode)
{
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << " " << number3 << " "
                       << message4 << endl;
      break;
   case 2:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << " " << number3 << " "
                       << message4 << endl;
      cout << "XVMC> " << endl;
      break;
   default:
      cout << "XVMC> " << message1 << " " << number1 << " "
                       << message2 << " " << number2 << " "
                       << message3 << " " << number3 << " "
                       << message4 << endl;
      break;
   }
   return;
}

#ifdef CHECK_ENERGY
// print out results of the energy conservation test
void energy_message(const char *type, const double e,
                    const char *unit, const int mode)
{
   cout.setf( ios::fixed, ios::floatfield );
   switch (mode)
   {
   case 1:
      cout << "XVMC> " << endl;
      cout << "XVMC> " << type << " "
           << setw(20) << setprecision(2) << e << " " << unit << endl;
      break;
   default:
      cout << "XVMC> " << type << " "
           << setw(20) << setprecision(2) << e << " " << unit << endl;
      break;
   }
   cout.unsetf( ios::fixed );
   cout.precision(6);
   return;
}
#endif // CHECK_ENERGY

// CPU time handler, returns process time in seconds
real etime(void)
{
   tms     time_info;
   long    int_time;
   real    cpu_time;

   times(&time_info);
   int_time = long(time_info.tms_utime);
   cpu_time = float(int_time)/float(sysconf(_SC_CLK_TCK));
#ifdef OSF1
   // divide by number of processes (DEC Unix only)
   cpu_time /= float(n_process);
#endif
   return(cpu_time);
}

// print message before starting the first batch,
// (print number of batches, processes and cpu time)
void start_message(int n, int np)
{
   real           etime();

   cout.setf( ios::fixed, ios::floatfield );
   cout << "XVMC> " << endl;
   cout << "XVMC> Batches:  " << setw(3) << n << ","
        << " processes: "     << setw(2) << np << ","
        << " CPU time: "      << setw(8) << setprecision(1) << ZERO << " s,"
        << " random status"   << endl;
   cout << "XVMC> ----------------------------------------------------------"
        << "-------" << endl;
   cout.unsetf( ios::fixed );
   cout.precision(6);

   return;
}

// print message after each batch
// (batch number, process id, cpu time and random generator status)
real batch_message(int n, int np, real cpu_start, ranmar_state state)
{
   real etime();
   real cpu_time = etime()-cpu_start;

   cout.setf( ios::fixed, ios::floatfield );
   cout << "XVMC> finished: " << setw(3) << n << ","
        << " process:   "     << setw(2) << np << ","
        << " CPU time: "      << setw(8) << setprecision(1)
        << cpu_time  << " s,  "
        << setw(4) << state.i1 << setw(4) << state.i2
        << endl;
   cout.unsetf( ios::fixed );
   cout.precision(6);
   return(cpu_time);
}
