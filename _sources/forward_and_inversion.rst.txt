Forward Modeling and Inversion Options
======================================

.. _forward-modeling:

Forward Modeling
-----------------

.. code-block:: text

    Usage: -F  rFile_Model rFile_Data wFile_Data [wFile_EMsoln rFile_fwdCtrl]

Required files are: 

* **rFile_Model** :  input model grid/resistivity file
* **rFile_Data** : input data file (template: defines data types/locations to return)
* **wFile_Data** : output data file (as defined from input)

The optional files are:

* **wFile_EMsoln** : output file for electric field solution (all periods and
    the 2 polarizations), defined on staggered grid edges
* **rFile_fwdCtrl** : forward solver control file. This allows more control over how the
    iterative solver works. This file is required for nested modeling, as this is
    where you specify the file name for the larger scale regional solution used for
    boundary data (last line below; see :ref:`nested-modeling`). The format for the
    forward solver control file is as follows:


.. code-block::

    Number of QMR iters per divergence correction : 40
    Maximum number of divergence correction calls : 20
    Maximum number of divergence correction iters : 100
    Misfit tolerance for EM forward solver : 1.0e-7
    Misfit tolerance for EM adjoint solver : 1.0e-7
    Misfit tolerance for divergence correction : 1.0e-5
    Optional EM solution file name for nested BC : nested.esoln

The first 6 parameters are numbers and are all required if you use this option.
The last line gives the file for nesting and is optional. If it is present,
electric fields in the file are used to provide boundary data; if it is not
present, default boundary values are used. A hash mark # can be inserted at the
start of the line instead The order of lines is fixed, so to use nested modeling
you have to specify values for the other 6. Default values are given above.
Instead of specifying a file name for forward model control, the user can
specify a single real number. ModEM will interpret this as the value of eps =
misfit tolerance for the forward solver (I.e., the 4th line in the file).

Inversion
---------

The simplest usage specifies only the search algorithm, the prior (and starting)
model file, and the data file. Default values for forward modeling, inversion
control parameters, as well as model covariance, are then used. Additional
optional arguments allow more control over the inversion. The general form for
the command line is:

.. code-block:: text

    Usage: -I NLCG rFile_Model rFile_Data [InvCtrl FwdCtrl CovCtrl StartModel ]

As usual, arguments have to be in the specified order. Thus, to specify
StartModel all other optional arguments must be given.

The required arguments are:

* **NLCG** Inversion search approach (other supported options **DCG**; only
    **NLCG** has been extensively tested and used to date). This argument is
    required by the ``-I`` option, even though most users only use the NLCG option.
* **rFile_Model** : input model grid/resistivity file : This is the prior model, which by
    default is used also as the starting model.
* **rFile_data** :  input data file (template: defines data types/locations to return)

Optional arguments are:

* **InvCtrl** : inversion control file, which has the format

.. code-block:: text

    Model and data output file name : Modular
    Initial damping factor lambda : 1.
    To update lambda divide by : 10.
    Initial search step in model units : 100.
    Restart when rms diff is less than : 2.0e-3
    Exit search when rms is less than : 1.05
    Exit when lambda is less than : 1.0e-4
    Maximum number of iterations : 120

The values given above are the default values. The first line is used to define
the root for output files (model parameter and data) created by the inversion;
these are discussed in more detail below, under :ref:`inversion-outputs`. As
with the forward modeling control file, there is a simpler form for this
optional argument, which modifies only a single (the most commonly modified)
parameter, i.e., lambda = the initial damping parameter for inversion.

* **FwdCtrl** : forward solver control file, described above under the ``-F``
    option.
* **CovCtrl** :  The covariance file. This is used to define the model
    covariance(see section :ref:`model-cov-files`).
* **StartModel** : The starting model file. This must be specified to start from
    sometthing other than the prior. Proper use of this option is discussed below,
    in :ref:`restarting`.

.. _inversion-outputs:

Inversion Ouptuts
------------------

For the NLCG inversion 4 files are produced for each each step in
the iterative search process. The files produced have the name ``root_NLCG_###.rho``,
``root_NLCG_###.prm``, ``root_NLCG_###.dat``, and ``root_NLCG_###.res``,
where ``root`` is the output file root specified in the inversion
control file, and ``###`` is the NLCG iteration number. The output
file root defaults to ``Modular`` if the inversion control file
is not used. The ``{*}{*}{*}.rho`` file is the model parameter file, giving
the resistivity estimate at the current iteration. The ``{*}{*}{*}.prm``
file is also in the model parameter file format, but represents the
transformed model parameter :math:`\tilde{{\bf m}}={\bf C_{m}}^{-1/2}({\bf {m}-{\bf {m}_{prior})}}`.


as discussed below inn :ref:`apply-cov`. The ``{*}{*}{*}.dat`` file contains
predicted data, in the ModEM list format (identical to the input data file), and
the ``{*}{*}{*}.res`` file contains residuals, again in the ModEM list format.
(This file is not useful and will be eliminated in future ModEM releases). In
addition to the suite of model parameter and data files produced, the inversion
provides an ongoing report on results of the various modeling and search steps.
These are displayed to standard output and should be monitored by users. Key
diagnostics of the search process are recorded in the log file
``root_NLCG.log``.

Similar files are produced by the DCG inversion, with the obvious name change.
Note in this case each iteration corresponds to a full Gauss-Newton step, so
there will be many fewer of these output files.  Howwever, this option has seen
little actual use so far.

.. _apply-cov:

Apply Model Covariance
----------------------

This option is used to apply the model covariance (smoothing operator) to a
model parameter, and to invert this operation (roughening). The forward
covariance operator is primarily used for testing; the inverse is required for
some specialized applications discussed below in :ref:`restarting`.  The format
for this option is: 

.. code-block:: text

    -C INV rFile_Model wFile_Model [CovCtrl rFile_Prior]

The required arguments are:

* **INV** : this is the inverse (roughening) operator, most often used. **FWD**
    is the argument used to apply the covariance (smoothing) operator.
* **rFile_Model** and **wFile_Model** are the input and output model.

Optional arguments are:

* **CovCtrl** : the covariance file; if this is not specified default covariance parameters
    are used.

* **rFile_Prior** : the prior model file; in general this is required for the most common
    use of the ``-C`` option.

Explanation: In the case of the INV option the output file contains the
transformed model parameter :math:`\tilde{{\bf m}}={\bf C_{m}}^{-1/2}({\bf{m}-{\bf {m}_{prior})}}`, 
where :math:`{\bf {m}}$, and ${\bf m}_{prior}` are the input and prior model
parameters, provided in files ``rFile_Model`` and ``rFile_prior``, respectively.
If the optional prior model file name is not provided :math:`{\bf {m}_{prior}}`
defaults to 0. The FWD option reverses the transformation, computing 
:math:`{\bf{m}={\bf C_{m}}^{1/2}\tilde{{\bf {m}}}+{\bf {m}_{prior}}}`.

.. warning:: 
    
    The inverse computation may be poorly conditioned if the smoothing
    is too strong, so it is a good idea to proceed with caution, and check your
    model after applying the inverse covariance.  The need for caution is one reason
    why we ask the user to do several steps themselves, instead of packaging this as
    a simpler to use black box.