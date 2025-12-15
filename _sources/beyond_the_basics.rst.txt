Beyond the Basics
=====================

Some things that most users will want to do with the inversion require use of
the optional arguments, in multi-step procedures. These include restarting the
inversion, and nesting a forward or inverse calculation within a larger, more
coarsely discretized model domain.

.. _restarting:

Restart the Inversion
---------------------

The inversion can be started from any model parameter file, not just the prior
model.  There are many reasons why a user might wish to do this. The simplest
case would be if the inversion has stopped due to a time limit on a computer
system, or some other reason, and the user wants to continue the run. In this
case, to restart the inversion from last iteration (e.g., iteration # 60) use
the following command line:

.. code-block:: text

    ModEM -I NLCG rFilePrior rFile_Data InvCtrl FwdCtrl CovCtrl root_NLCG_060.prm

One might also want to change some inversion control or covariance parameters.
For example, it is possible that an inversion that has stalled in terms of
decreasing the normalized RMS, might be moved out of a local minimum by
restarting the inversion with a larger value of the damping factor lambda.
Restarting the inversion after removing structures (e.g., setting a conductive
body to the prior resistivity), and possibly freezing the model to some value,
will be useful for hypothesis testing.

The key point that must be understood to start the inversion from something
other than the prior is this: As already noted in the discussion of the ``-C``
option, the inversion search is conducted in the "transformed" model space,
i.e., we directly minimize the L2 norm of :math:`$\tilde{{\bf {m}}}={\bf
{C_{m}}^{-1/2}({\bf m}-{\bf m}_{prior})}$`, and then obtain the actual
(physical) model parameter as :math:`${\bf m}={\bf m}_{prior}+{\bf
C_{m}}^{1/2}\tilde{{\bf m}}$`.  The starting model file has to be given in the
transformed model space, i.e., what is output in the ``{*}{*}{*}.prm`` file. For
a simple restart of the inversion (picking up where the search has gone so far)
you would use the name of the appropriate file from the last iteration completed
(for example ``Example_NLCG_060.prm``) as the ``StartModel`` (last optional
argument) when restarting the inversion with the ``-I`` option. The same procedure
would be used when restarting with a larger value of the damping parameter
lambda.

Some other restarting cases, where either the model is edited, or the covariance
parameters (including frozen sub-domains) involve an additional step. In these
cases you have to first create the appropriate ``{*}{*}{*}.prm`` file. We give
two slightly different examples here; the basic idea is the same in both (and
other cases). Suppose you wish to modify the inversion result, removing some
structure to see if it is required. Let’s call the modified model parameter file
``HypTest.rho`` (containing log resistivity, modified with some tool from the
inversion result). As a first step you could obviously run the forward code,
using ``HypTest.rho`` as the input model. To go further, you might want to
restart the inversion from the modified model. This must first be transformed
(roughened), using the ``-C`` inverse option:

.. code-block:: text

    ModEM -C INV HypoTest.rho HypTest.prm CovCtrl rFilePrior

Then restart the inversion:

.. code-block:: text

    ModEM -I NLCG rFilePrior rFile_Data InvCtrl FwdCtrl CovCtrl HypoTest.prm

Here ``CovCtrl`` is the covariance parameter file, and ``rFilePrior`` is the
prior model file, which in the application just described would be the same as
used for the initial inversion.

Going one step further in the hypothesis test example, one might want to remove
a structure (e.g., a conductive feature) and then freeze the region where
conductivity has been set to background levels (thus verifying whether data can
be fit without the tested structure). To do this the covariance file would have
to be modified to now include the frozen region. This modified covariance file
would then be used in both steps given above (``-C`` and ``-I``).

As another example, one might want to change the covariance (e.g., changing
smoothing length scales, etc.), and restart from results obtained with an
earlier run to (hopefully) accelerate convergence.  This would be done in the
same way, except the input model parameter file on the first (``-C``) step would
be something like ``Example_NLCG_060.rho`` — i.e., convert the final model
parameter file to the transformed (``{*}{*}{*}.prm``) form appropriate to the
new covariance, and then use this for the starting model.


.. _nested-modeling:

Nested Modeling
---------------

Boundary conditions for the ModEM forward model are specified tangential fields
on domain boundaries (top, bottom, and sides). By default, boundary data are
computed from the prior model conductivity (at present using solutions to a
series of 1D MT problems). When this simple default is used it is advisable to
keep boundaries distant from the interior model domain of interest. ModEM also
supports a *nesting* option, where boundary data are taken from a 3D
forward solution in a larger (and typically more coarsely discretized) domain.
This approach requires computing a forward solution on this larger domain, and
saving the electric field solution vector in file ``wFile_EMsoln`` (first
optional argument in ``-F`` option). The model parameter file for this large
domain could be obtained from running the inversion, or from a prior model (for
example containing a realistic ocean over an area too large to use for the final
high resolution inversion run).  Note that there is no option to output the
electric field solution with the ``-I`` option, so if an inverse model from a
large grid is used for nesting, the forward model must be rerun (on the large
grid) using the final inversion result.

Use of the nested option, for either forward (``-F``) or inverse (``-I``)
calculations is triggered by presence of the final (optional seventh) line in
the ``fwdCtrl`` file (see section :ref:`forward-modeling`.  The name of the
electric field solution file (i.e., ``wFile_EMsoln`` used when running
``-F`` on the large grid) should be specified on this line to turn on the
nesting option. 

.. note:: 
    The periods (including order) in the data file used for the large grid
    ``-F`` run must be identical to those used for the nested forward/inverse
    modeling.  The simplest approach is to use the same data file
    (``rFile_Data``) for both large and nested runs, but only the periods (not
    the sites) have to be the same.

Here is a simple example illustrating the use of nested modeling for both
forward modeling (``-F`` option) and inversion (``-I`` option). Assume that you
have constructed a large coarse grid which includes e.g., an ocean. Let us call
the large coarse grid ``Large_Grid.mod`` and the data file ``Data_file.dat``.
In the first run we need to let ModEM compute the electric field solution and
save it in a file. Thus, in the first run we need the ``-F`` option with the
optional output ``wFile_EMsoln``:

.. code-block:: text

    ModEM -F  Large_Grid.mod Data_file.dat temp_Predicted_Large.dat EMsoln_large.soln 

After running this command line, ModEM will produce 2 files: 1) , 

1. ``temp_Predicted_Large.dat`` - the predicted data for the large grid (not important), 
2. ``EMsoln_large.soln``, the electric field solution. To run the nested
    modeling using either forward (``-F``) or inverse (``-I``) options you
    need to first to modify the forward solver control file ``rFile_fwdCtrl`` to
    include the name of the electric field solution obtained from the first run:

.. code-block:: text

    Number of QMR iters per divergence correction : 40
    Maximum number of divergence correction calls : 20
    Maximum number of divergence correction iters : 100
    Misfit tolerance for EM forward solver        : 1.0e-7
    Misfit tolerance for EM adjoint solver        : 1.0e-7
    Misfit tolerance for divergence correction    : 1.0e-5
    Optional EM solution file name for nested BC  : EMsoln_large.soln

For the second run (either forward (``-F``) or inverse (``-I``) ) construct a
smaller and generally finer grid and call it e.g., ``Small_grid.mod``. To run
forward modeling use:

.. code-block:: text

    ModEM -F Small_grid.mod Data_file.dat Predicted_small.dat EMsoln_small.soln FwdCtrl

and for the inversion:

.. code-block:: text

    ModEM -I NLCG Small_grid.mod Data_file.dat InvCtrl FwdCtrl CovCtrl

.. important:: 

    IT IS VERY IMPORTANT TO NOTICE HERE THAT IN THE FIRST AND SECOND RUNS WE
    USED THE SAME DATA FILE ``DATA_FILE.DAT``.  THIS IS BECAUSE ModEM DOES NOT
    INTERPOLATE OVER PERIODS. THUS, TO MAKE SURE THAT THE SAME PERIOD LAYOUT
    SAVED IN ``EMSOLN_LARGE.SOLN`` IS USED IN THE SECOND RUN TO EXTRACT THE
    BOUNDARY VALUES AT THE EDGES OF THE NESTED GRID, USE THE SAME DATA FILE.

.. _flexible-air-layers:

Flexible Air Layers
--------------------
	
Air layers can now be set in the forward solver configuration file. Options can
be printed to screen using ``./Mod3DMT -F`` command without the additional
input arguments. The new capability can be executed by adding two more lines to
the configuration file (lines 8 and 9), after the nested modeling file. If you
are not using the nested modeling option but would like to change the air
layers, you should still include line 7, but leave the name of the EM solution
file in line 7 blank.  Flexible air layers are only implemented for 3D.
	
Three options for the air layers are: ``mirror``, ``fixed height``, ``read
from file``. For backwards compatibility, the default is ``mirror 10 3.
30.``, so that the code should still produce identically the same results to
the previous version, unless the air layers are changed.

Option 1 - Mirrors
^^^^^^^^^^^^^^^^^^^

The air layers based on the uppermost layers of the Earth resistivity model. In
the code, air layers are counted from top to bottom (top is 1, bottom is N). The
logic is as follows,

.. math:: 

	airlayer((N+1)-j) = \alpha^{j-1} * earthlayer(j)

Where :math:`\alpha` and the number of air layers :math:`N` are supplied by the
user, and :math:`j` is an integer index from 1 to :math:`N`. When :math:`j=1`,
the bottom air layer :math:`N` equals the top earth layer 1. When :math:`j=N`,
the top air layer 1 equals :math:`\alpha^{N-1} * earth layer(N)`. Finally, the
third user-supplied parameter `MinTopDz` specifies the minimum width of the top
air layer, in km. If the top layer happens to be thinner than specified by this
parameter, it is adjusted to this value. The parameters are specified in the
following order: :math:`N`, :math:`\alpha`,  `MinTopDz`.

.. code-block:: text

	Air layers mirror|fixed height|read from file : mirror
	Number of air layers and min top dz in km     : 10 3. 30.

This procedure is unnecessarily complex, and has created an unreasonably large
air domain for most grid setups in the past. It is maintained now as the default
for backwards compatibility. However, this option is not recommended.

Option 2 - fixed height
^^^^^^^^^^^^^^^^^^^^^^^

In this setup, the width of the air layers increases logarithmically until it
reaches the user-defined height. Additionally, the user defines the total number
of the air layers, and the logarithmic exponent is calculated automatically from
these two parameters. This is the recommended option.

.. code-block:: text

	Air layers mirror|fixed height|read from file : fixed height	
	Number of air layers and max height in km     : 12 1000.


Option 3: read from file
^^^^^^^^^^^^^^^^^^^^^^^^

In this setup, the user specifies the total number of the air layers, and their
width from top to bottom, in km, directly in the forward solver configuration
file. This provides additional flexibility for advanced users, allowing in
particular to compare modeling results more directly to those produced by
alternative codes or analytical solutions.

.. code-block:: text

    Air layers mirror|fixed height|read from file : read from file
    Number of air layers and dz top to bottom km  : 10 500. 200. 100. ... 2. 1. 0.5

Recommended air layers option is ``fixed height 12 1000``. The users will
likely find out that different air layer setups create substantial differences
in the modeled impedances. We suggest that this issue is carefully explored by
the users before they switch to a different setup.

.. _log10-res-model:

LOG10 Resistivity Model File Option
-----------------------------------

To exercise this new option, create a WS or Mackie format model file with
``LOG10`` resistivity instead of ``LOGE`` or ``LINEAR``. This is
specified in the second line of the resistivity model input file (the first line
is a comment). Internally in the code, we are keeping ``LOGE`` for all
computations; but we convert ``LOG10`` or ``LINEAR`` to ``LOGE``, make
note of the user preference, then convert back for output. This feature is
particularly useful for synthetic models, since it makes the input/output files
more readable. 

Having said that, we should note that we appreciate the inadequacy of an ASCII
input file format for modern high performance 3D modeling situations. An option
utilizing modern binary input/output formats will be added in a future release.

.. _unit-testing:

New Unit Testing Features for Code Developers
----------------------------------------------

This release includes a collection of fully developed and debugged symmetry
tests, gradient test, and the full Jacobian test. The code now passes all tests;
all symmetry tests (including S, i.e., the forward solver) are passed with
arbitrary random perturbations as input. The symmetry test ``[TEST_ADJ]`` has
many sub-options, which can be explored by typing ``./Mod3DMT -A`` command
without the additional input arguments. Options to generate random perturbation
files for model ``[-A m]``, data ``[-A d]``, RHS ``[-A b]``, E-soln ``[-A e]``
are also included.

.. code-block::

    [TEST_GRAD]
    -g  rFile_Model rFile_Data rFile_dModel [rFile_fwdCtrl rFile_EMrhs]
    Runs the ultimate test of the gradient computation based on
    Taylor series approximation.
    [TEST_ADJ]
    -A  J rFile_Model rFile_dModel rFile_Data [wFile_Model wFile_Data]
    Tests the equality d^T J m = m^T J^T d for any model and data.
    Optionally, outputs J m and J^T d.
    [TEST_SENS]
    -S  rFile_Model rFile_dModel rFile_Data wFile_Data [rFile_fwdCtrl]
    Multiplies by the full Jacobian, row by row, to get d = J m.
    Compare to the output of [MULT_BY_J] to test [COMPUTE_J]


