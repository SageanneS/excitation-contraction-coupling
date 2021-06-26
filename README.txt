These files contain the base code needed to run the model described in the article:

"Sageanne Senneff, Madeleine M. Lowery,
Effects of extracellular potassium on calcium handling and force generation in a model of excitation-contraction coupling in skeletal muscle,
Journal of Theoretical Biology,
Volume 519,
2021,
110656,
ISSN 0022-5193,
https://doi.org/10.1016/j.jtbi.2021.110656.
(https://www.sciencedirect.com/science/article/pii/S0022519321000783)
Abstract: It is well-established that extracellular potassium (Ko+) accumulation reduces muscle fiber excitability, however the effects of Ko+ on the excitation–contraction coupling (ECC) pathway are less understood. In vivo and in vitro studies following fatiguing stimulation protocols are limited in their ability to capture the effects of Ko+ on force production in combination with other simultaneously changing factors. To address this, a computational model of ECC for slow and fast twitch muscle is presented to explore the relative contributions of excitability-induced and metabolic-induced changes in force generation in response to increasing K+o. The model incorporates mechanisms previously unexplored in modelling studies, including the effects of extracellular calcium on excitability, calcium-dependent inhibition of calcium release, ATP-dependent ionic pumping, and the contribution of ATP hydrolysis to intracellular phosphate accumulation rate. The model was able to capture the frequency-dependent biphasic Force-K+o response observed experimentally. Force potentiation for moderately elevated K+o was driven by increased action potential duration, myoplasmic calcium potentiation, and phosphate accumulation rate, while attenuation of force at higher K+o was due to action potential failure resulting in reduced calcium release. These results suggest that altered calcium release and phosphate accumulation work together with elevated Ko+ to affect force during sustained contractions.
Keywords: Muscle force; Calcium release; Phosphate accumulation; Mathematical model; Sarcoplasmic reticulum"


The model describes excitation-contraction coupling in a single slow twitch skeletal muscle fiber and a single fast twitch skeletal muscle fiber.
The code to run each model is given in "Run_Fast_Twitch_Model.m" and "Run_Slow_Twitch_Model.m", respectively. The file "ics 4.txt" sets the initial
resting membrane potential for each sarcolemma compartment, at each tubular node, for an extracellular potassium concentration of 4 mM.

The data is output as a text file. 

To run the model at steady-state, the current can be introduced after 8 seconds. (Update the bounds on "t" in Lines 379/377 (Fast Twitch/Slow Twitch)).



