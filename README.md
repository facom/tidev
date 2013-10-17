tidev
=====

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
