tidev
=====

A general framework to calculate the Tidal Evolution of interacting
bodies.

Quick start
-----------

1. Get a copy of tidev from https://github.com/facom/tidev:

   NOTE: You can get an anonymous clone of the project using:

      git clone git://github.com/facom/tidev.git

   Then you will be able to get updates using 'git pull'.

2. Compile utilities:

      make utilbuild

3. Copy an example and prepare a run:

      cp -rf examples/*.config .
      make prepare

4. Launch the run:

      make go

To know more read the doc/MANUAL.txt.

For the contirbutor
-------------------

1. Generate a public key of your account at the server where you will
   develop contributions:

   $ ssh-keygen -t rsa -C "user@email"

2. Upload public key to the github project site
   (https://github.com/facom/tidev).  You will need access to the
   account where the repository was created.

3. Configure git:

   $ git config --global user.name "Your Name"
   $ git config --global user.email "your@email"

4. Get an authorized clone of the master trunk:

   $ git clone git@github.com:facom/tidev.git

License
-------

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.

Copyright (C) 2013 Jorge I. Zuluaga, Mario Melita, Pablo Cuartas,
Bayron Portilla

