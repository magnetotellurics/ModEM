Contributing
=========================

Users are welcome to contribute to ModEM and we welcome contributions! 

We highly recommend opening a GitHub Issue first so that new contributions can be
discussed with ModEM authors, especially for large pull requests.

Contributions can be made via a `GitHub pull request <_github_pr>`_ at the ModEM
repository: https://github.com/magnetotellurics/ModEM/compare.

.. _github_pr: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests

Steps To Open a Pull Request (PR)
----------------------------------

See: `Contributing to a project`_.

.. _Contributing to a project : https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project

Tips & Testing
---------------

* Please keep pull requests to reasonable sizes (consider breaking up large PR's
  into small ones)
* If you are able, please test your changes using both the GNU gfortran
  compiler and Intel ifort repository (but Gfortran at the minimum if you don't
  have access to the Intel compiler).
* Contribute any relevant documentation to this Documentation (See
  :ref:`building-the-docs`)
* Run the unit tests in the ``unit_testing`` directory (Currently not added)

.. _building-the-docs:

Building this Documentation
----------------------------

This `Sphinx <_sphinx_page>`_ documentation can easily be built if you have Sphinx installed. If you
have Python and Pip installed, you can install Sphinx by:

.. code-block:: bash

    $ pip install -U sphinx

For more information on installing Sphinx see: https://www.sphinx-doc.org/en/master/usage/installation.html.

After making changes to the documentation, you can build locally for testing by doing the following:

.. code-block:: bash

    $ cd ModEM/docs
    $ make html

The HTML files will be built in ``ModEM/docs/build``. You can open these locally
and view them and ensure they look correct. You do not need to (and will not be
able to) commit the built HTML files.

If you change the documentation, please ensure that it builds without errors.

.. _sphinx_page: https://www.sphinx-doc.org/en/master/
