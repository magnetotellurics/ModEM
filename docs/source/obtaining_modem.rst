Obtaining ModEM
===============

ModEM 2D and 3D MT modeling and inversion codes are currently available for
under the LICENSE agreement (see LICENSE). The latest stable version of ModEM
can always be obtained from the GitHub repository:

https://github.com/magnetotellurics/ModEM

It is advisable that you obtain the code through the GitHub repository, to keep
your version up-to-date with bug fixes and updates.

You can clone the repository by using Git from the command line::

    $ git clone https://github.com/magnetotellurics/ModEM.git
    Cloning into 'ModEM'...
    remote: Enumerating objects: 18608, done.
    remote: Counting objects: 100% (199/199), done.
    remote: Compressing objects: 100% (54/54), done.
    remote: Total 18608 (delta 162), reused 152 (delta 145), pack-reused 18409 (from 1)
    Receiving objects: 100% (18608/18608), 73.43 MiB | 13.34 MiB/s, done.
    Resolving deltas: 100% (14565/14565), done.


ModEM Branches
---------------

As noted in the README, there are several branches on ModEM which have been
converted from ModEM's OSU CEOAS Subversion (SVN) repository. Likewise, there
are some specific branches worth noting.


Main and Develop
^^^^^^^^^^^^^^^^^

The `main branch <https://github.com/magnetotellurics/ModEM>`_ of ModEM holds
the latest release of ModEM and should contain code that has been more used and
tested and should be more stable.

The `develop branch <https://github.com/magnetotellurics/ModEM/tree/develop>`_
on the other hands contains the latest code changes the ModEM. While these
features are new, they may not be as stable as they have not undergone
significant testing.

SVN Branches
^^^^^^^^^^^^^

The complete SVN history has been preserved and converted into a git history.
Some SVN branches, which SVN users my be familar with, have been renamed:

===============   ===============
SVN Branch Name   Git Branch Name
===============   ===============
trunk              trunk
stable-            main 
stable             classic
===============   ===============

