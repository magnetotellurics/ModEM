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

    Usage: -I NLCG rFile_Model rFile_Data [ InvCtrl FwdCtrl CovCtrl StartModel ]

As usual, arguments have to be in the specified order. Thus, to specify
``StartModel`` all other optional arguments must be given.

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

For the joint distortion inversion (``NLCG_DIST``), the data to be corrected
for galvanic distortion are given with data type ``Full_Impedance_Dist``. The
standard control lines above may be followed by optional distortion lines
(defaults shown):

.. code-block:: text

    Distortion damping factor nu       : 0.1
    Distortion output filename         : distortion.dat
    Reference distortion filename      :
    Distortion type (full|diag)        : full

Unlike the standard lines, these are identified by their label (the text
before the colon) and may be given in any order; any of them may be omitted, 
with the above values as default. 
TODO: I should extend the similar label based syntax to the standard lines at 
some stage, in another commit. 

The distortion matrice C is one real 2x2 matrix per site, independent of
frequency. The distortion matrice C = g * T * S * A, where g is a scalar 
gain, T is a twist matrix, S is a shear matrix, and A is an anisotropic 
matrix. When there is no distortion, C is the identity matrix. C is estmated 
jointly along with the conductivity model, and is regularised to be close to a
reference distorion matrix C^0 (identity if not specified). 

After every forward solve, C is replaced by the exact minimiser of the 
penalty functional for the current EM fields. This is a trivial least-squares 
problem per site and needs no additional EM solves. NLCG therefore minimises 
the reduced functional :math:`\min_C F(\tilde{\bf m}, C)` with the *same* 
line search and lambda strategy as the standard NLCG.

* ``Distortion damping factor nu``: strength of the distortion penalty
  :math:`\tfrac12 \sum_j w_j \|C_j - C^0_j\|_F^2`. The weights :math:`w_j`
  are computed from the observed impedances and errors, so that nu is
  dimensionless: at each site the estimated distortion moves roughly a
  fraction :math:`1/(1+\nu)` of the way from the reference towards the best
  fit allowed by the data. Small values (e.g. 0.1) let C absorb most of the
  site-specific distortion; large values (e.g. 100) hold C at the reference, 
  which is somehow equivalent to ignoring distortion.
* ``Distortion output filename``: the estimated distortion matrices are
  written to ``root_<filename>``.
* ``Reference distortion filename``: optional reference distortion
  :math:`C^0`, in the same format as the output file (sites are matched by
  their index). The identity is used if this is empty or the line is omitted.
  A reference from an independent analysis (e.g. tensor decomposition) can be
  given here; with a large nu the distortion is effectively fixed at it.
* ``Distortion type``: ``full`` estimates the general 2x2 matrix (4
  parameters per site). ``diag`` estimates one gain per horizontal electric
  field component, :math:`C_j = \mathrm{diag}(g_x, g_y)` (2 parameters per
  site); this is the static-shift/splitting part of the distortion and
  reproduces the frequency-independent offset between the xy and yx
  apparent-resistivity curves, while leaving twist and shear to the
  conductivity model. For ``diag`` the reference distortion is projected
  onto its diagonal.

* **FwdCtrl** : forward solver control file, described above under the ``-F``
    option.
* **CovCtrl** :  The covariance file. This is used to define the model
    covariance(see section :ref:`model-cov-files`).
* **StartModel** : The starting model file. This must be specified to start from
    something other than the prior. Proper use of this option is discussed below,
    in :ref:`restarting`.

.. _inversion-outputs:

Inversion Ouptuts
------------------

For the NLCG inversion 4 files are produced for each each step in the iterative
search process. The files produced have the name ``root_NLCG_###.rho``,
``root_NLCG_###.prm``, ``root_NLCG_###.dat``, and ``root_NLCG_###.res``, where
``root`` is the output file root specified in the inversion control file, and
``###`` is the NLCG iteration number. The output file root defaults to
``Modular`` if the inversion control file is not used. The ``{*}{*}{*}.rho``
file is the model parameter file, giving the resistivity estimate at the current
iteration. The ``{*}{*}{*}.prm`` file is also in the model parameter file
format, but represents the transformed model parameter :math:`\tilde{{\bf
m}}={\bf C_{m}}^{-1/2}({\bf {m}-{\bf {m}_{prior})}}`.


As discussed below in :ref:`apply-cov`. The ``{*}{*}{*}.dat`` file contains
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
there will be many fewer of these output files.  However, this option has seen
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
