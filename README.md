[LEGAL INFO]

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.


[INSTALLATION INSTRUCTIONS UNDER LINUX]

You will need C++11 and preinstalled CAPD 4.0 libraries.
The latter can be downloaded at http://capd.ii.uj.edu.pl/
also see http://mnich.ii.uj.edu.pl/capd/doc4/capd_compilation.html for compilation instructions
and http://mnich.ii.uj.edu.pl/capd/doc4/ for documentation.
Note: it is advisable to run 
  ldconfig
after installation of the libraries.

To compile the source, first open makefile
and change the right hand side of 
  CAPDBINDIR =/usr/local/bin/
to the folder where the capd scripts are installed.
Then call 
  make 
inside the program directory.

To execute the program call
  ./fhn

To show details from the proof at runtime change the verbose variable in fhn.cpp from 0 to 1 and recompile.

Program was tested under GCC 4.9.2 and CAPD SVN revision 568.


[CONTACT INFO]

Aleksander Czechowski
czechows@ii.uj.edu.pl
