===============
Radial Velocity
===============

The ``rv.py`` script allows you to perform radial velocity (rv) calculations and create radial velocity plots to aid in planning of radial velocity observations.
The radial velocity calculations are handled by an :ref:`RV <rvclass>` class, while a :ref:`JulianDate <juliandateclass>` class is used to handle date/time conversions.

Parameter file
==============
``rv.py`` requires a parameter file to specify the orbital parameters of the system you wish to analysis. A template is provided in ``data/template_params.txt`` to help you get started. Comment lines starting with ``#`` and in-line comments are ignored.

For a basic rv calculations the standard rv parameters are required, ``k1`` [km/s], ``omega`` [deg], ``eccentricity``, ``tau`` [days], ``period`` [days], as well as the ``name`` parameter.

If the mean system rv offset, ``mean_val`` (usually referred to as gamma), is not provided in the parameter file it is set to 0 km/s. The ``ignore_mean`` keyword in some functions can also be used to use a 0 km/s mean_val.

To include the rv of a companion the parameters ``m_star`` (star mass) and ``msini`` or ``m_true`` (companion mass) are required or ``k2`` the semi-major amplitude of the companion.
If ``k2`` is not provided it is calculated from ``k1`` and the star and companion masses.

.. note::
    A future version could maybe have the option to obtain parameters from planetary databases such as `exoplanet.eu <http://exoplanet.eu/>`_. Although this functionality would be limited to the stars/planets of the databases.


Usage examples
==============
Simple usage cases:

::

    python rv.py data/HD30501_params.txt

.. image:: phase_curve.png
    :height: 400 px
    :width: 600 px
    :scale: 90 %
    :alt: Radial velocity phase curve.
    :align: center

displays the phase curve of HD30501 and its companion.

::

    python rv.py data/HD30501_params.txt -l data/HD30501_obs.txt -m time -d 2013-01-01

Will create a temporal RV curve, marking the locations of the observations provided in the observation list file.

.. image:: time_curve.png
   :height: 400 px
   :width: 600 px
   :scale: 90 %
   :alt: Radial velocity curve, with observations indicated.
   :align: center

More usage examples can be found here:
Modes
======
Different functionally can be accessed from the modes flag. ``-m``, or ``--mode``. The dmoefault mode is ``phase``.

phase
------
Produces the RV phase curve with ``cycle_fraction`` and ``phase_ceter`` configurable.
If the ``k2`` parameter is provided or the mass of the host (``m_host``) and companion (``msini`` or ``m_true``), then the RV for the companion is plotted on the second y-axis.

time
-----
Produces a temporal RV curve of the system. The reference day is today but can be changed using the ``-d``, ``--date`` option.


Marking observations
-------------------
To mark the position of past and/or future observations you the dates can be provided with the ``-o`` (times passed to command line), or  ``-l`` (file with list of observations) options.
The observation times are marked/coloured as ``past`` or ``future``, relative to the reference date.

The temporal space ranges from the maximum extent of the ``cycle_fraction`` after the reference date and any marked observation times.


debug
-----
You can turn on debugging information using the ``--debug`` flag, e.g.::

    python rv.py data/HD30501_params.txt -l data/hd30501_obs.txt --debug


.. _rvclass:

RV class
--------
.. autoclass:: utils.rv_utils.RV
   :members:
   :undoc-members:
   :show-inheritance:


.. _juliandateclass:

JulianDate
----------
Used to convert from datetime objects and strings into julian dates and back again.
This was created because epehm.julain_date() only converted to julian date.

.. autoclass:: utils.rv_utils.JulianDate
   :members:
   :undoc-members:
   :show-inheritance:
