# Hydrogen fast filling modeling, optimisation and control

This was my first semester project, which I have done at the Laboratory for Reneweable Energy Science and Engineering (LRESE) at the EPFL.

The main goal consisted in modeling the rapid filling of a gaseous hydrogen tank and to control it in an optimal way so that filling duration is minimised while critical 
evolve in a secure manner. Those are essentially the pressure and the inner tank temperature.

As a first step, a complete mathematical model of the system was derived. It uses thermodynamics tables to ensure the accuracy of the results and has been compared to state-of-the art 
models found in the literature, as well as validated using measurements, which were also found within the literature.

Then, the whole problem was transcripted into an optimisation problem. It was first solved in the aim of finding an optimal control input history, in the same vein as one would proceed for Model-Predictive-Control.

Afterwards, the optimisation problem was tweaked so that it inculdes the full system in closed-loop operation. The objective was to determine optimal coefficients for a PID controller while satisfying several constraints regarding its states.

Finally, an Adaptive PID controller was proposed so that the controller is more robust to model uncertainties. Those were first defined, studied (parametrically) and then the adaptive law was trained and designed.

The results show a nice agreement with the literature while establishing quite novel performance considering that pure optimal control problems had not been really adressed regarding hydrogen fast filling so far.
