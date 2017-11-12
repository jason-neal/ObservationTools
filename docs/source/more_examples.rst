===============
More Examples
===============

Here are some more examples of rv.py usages.

These use these observation dates.
::

    # obstimes.txt
    2012-02-13
    2012-04-15
    2013-07-21
    2013-09-30
    2016-01-15

phase
-----
Simple usage cases:



displays the phase curve of HD30501 and its companion.

time
----
More time mode examples

::
::

    python rv.py data/HD30501_params.txt

.. image:: phase_curve.png
    :height: 400 px
    :width: 600 px
    :scale: 90 %
    :alt: Radial velocity phase curve.
    :align: center


Show 1.5 phases.
::

    python rv.py data/HD30501_params.txt -c 1.5

.. image:: phase_curve.png
    :height: 400 px
    :width: 600 px
    :scale: 90 %
    :alt: Radial velocity phase curve.
    :align: center

Center phase curve on 0.5

::

    python rv.py data/HD30501_params.txt -p 0.5

.. image:: phase_curve.png
    :height: 400 px
    :width: 600 px
    :scale: 90 %
    :alt: Radial velocity phase curve.
    :align: center

    python rv.py data/HD30501_params.txt -l data/HD30501_obs.txt -m time -d 2013-01-01

.. image:: time_curve.png
   :height: 400 px
   :width: 600 px
   :scale: 90 %
   :alt: Radial velocity curve, with observations indicated.
   :align: center


Using both ``-l`` and ``-o``

::

    python rv.py data/HD30501_params.txt -l data/HD30501_obs.txt -m time -d 2013-01-01 -o 2014-03-12 2014-07-22


.. image:: time_curve.png
   :height: 400 px
   :width: 600 px
   :scale: 90 %
   :alt: Radial velocity curve, with observations indicated.
   :align: center

::

    # All future points
    python rv.py data/HD30501_params.txt -l obstimes.txt -m time -d 2012-01-01


.. image:: time_curve.png
   :height: 400 px
   :width: 600 px
   :scale: 90 %
   :alt: Radial velocity curve, with observations indicated.
   :align: center

::

    python rv.py data/HD30501_params.txt -l obstimes.txt -m time -d 2013-08-01


.. image:: time_curve.png
   :height: 400 px
   :width: 600 px
   :scale: 90 %
   :alt: Radial velocity curve, with observations indicated.
   :align: center

::

    # All past observations
    python rv.py data/HD30501_params.txt -l obstimes.txt -m time -d 2017-01-01


.. image:: time_curve.png
   :height: 400 px
   :width: 600 px
   :scale: 90 %
   :alt: Radial velocity curve, with observations indicated.
   :align: center

And with using the ``cycle_ration``
::

    # All past observations
    python rv.py data/HD30501_params.txt -l obstimes.txt -m time -d 2017-01-01 -c 0.2

.. image:: time_curve.png
   :height: 400 px
   :width: 600 px
   :scale: 90 %
   :alt: Radial velocity curve, with observations indicated.
   :align: center
