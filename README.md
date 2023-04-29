# muscleOptimizer
WIP

Optimize muscle's Optimal Fiber Length and Tendon Slack Length as proposed by [Modenese et al. 2016](https://doi.org/10.1016/j.jbiomech.2015.11.006).

Compared to the [original script](https://github.com/modenaxe/MuscleParamOptimizer), the last three changes lead to significant improvements in computational efficiency:

* The method of identifying joints spanned by muscles is improved; i.e., 3 values within each coordinate's range of motion (min, intermediate, and max) are tested to see if there are any changes in muscles' length with a threshold of 0.1 mm. This works for all models and muscles, e.g., trunk muscles.

* Added an option to optimize the muscles of one leg and set the optimized values for both, in case of using identical scale factors for both legs.

* Added an option to store reference model's muscle quantities as a Json file for later use (other subjects).

* In sampling muscle quantities, the loop changed to iterate over shared coordinates instead of muscles, so that in case of Rajagopal model, it requires only 10 iterations rather than 80 to collect quntities of all muscles (~15 min for each leg):

```
'hip_flexion_r; hip_adduction_r; hip_rotation_r': 
['addbrev_r', 'addlong_r', 'addmagDist_r', 'addmagIsch_r', 'addmagMid_r', 'addmagProx_r', 'glmax1_r', 'glmax2_r', 'glmax3_r', 'glmed1_r', 'glmed2_r', 'glmed3_r', 'glmin1_r', 'glmin2_r', 'glmin3_r', 'iliacus_r', 'piri_r', 'psoas_r'],

'hip_flexion_r; hip_adduction_r; hip_rotation_r; knee_angle_r': 
['bflh_r', 'grac_r', 'recfem_r', 'sart_r', 'semimem_r', 'semiten_r', 'tfl_r'],

'knee_angle_r':
['bfsh_r', 'vasint_r', 'vaslat_r', 'vasmed_r'],

'knee_angle_r; ankle_angle_r; subtalar_angle_r':
['gaslat_r', 'gasmed_r'],

'ankle_angle_r; subtalar_angle_r': 
['edl_r', 'ehl_r', 'fdl_r', 'fhl_r', 'perbrev_r', 'perlong_r', 'soleus_r', 'tibant_r', 'tibpost_r']
```