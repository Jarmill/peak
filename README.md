# Peak Estimation
Estimate the peak of a function p(x) across all trajectories of a dynamical system x'=f(t,x) starting from a set X0 over time \[0, T\]. If possible, recover the trajectories that attain the peak value are also recovered.

An occupation-measure framework is used to find a convergeng sequence of upper bounds in rising degree to the true peak value. The approximate-optimal trajectories may be recovered if moment matrices satisfy rank conditions (up to numerical accuracy).

## Dependencies

- Gloptipoly3: http://homepages.laas.fr/henrion/software/gloptipoly/. A customized version is included in the repository
- YALMIP: https://yalmip.github.io/
- Mosek: https://www.mosek.com/ (or any solver compatible with YALMIP)

Further documentation and directions to come.

## Instructions
The `peak_estimate` routine has two arguments: `p_opt` and `order`. The `order` is the relaxation order involving moments of degree `2*order`. 

`p_opt` is an options structure of type `peak_options` defining properties of the system, with the fields:
```
  var:        Structure of symbolic variables (@mpol)
      t:      time (default empty)
      x:      state
      w:      parametric uncertainty (default empty)

  Tmax:       Maximum time (if var.t is not empty)

  dynamics:   Structure with fields (f, X)
      f:      dynamics x' = f(t,x, w) over the space X.
                  Each entry is a polynomial                
      X:      over what space do the dynamics evolve

  obj:        Functions to maximize along trajectories
      Scalar:     Single function
      List:       Maximize the minimum of all objectives
                     Each entry is a polynomial       
  
  state_supp: Support set of total set X (@supcon)
  state_init: Support set of initial set X0 (@supcon)
  Tmax:       Maximum time to consider (can be infinite if f(t,x) = f(x))
      
  rank_tol:   Rank tolerance for moment matrix to be rank-1
  box:        Box containing valid region X, default to [-1, 1]^n
```

Output is stored in the structure `out`. Trajectories are sampled by the function `switch_sampler` after defining a `sampler()` function to randomly generate a trajectory on X0. The visualizing functions `state_plot, cost_plot, nonneg_plot` illustrate properties of trajectories.

Examine and run `experiments/{pendulum_test, time_var_circ_2_1, sym_attractor_4_1, flow_half_safety}.m` as examples.

## Reference
https://arxiv.org/abs/2009.06120 "Peak Estimation and Recovery with Occupation Measures"

If you find this work useful, please cite:

'''
@misc{miller2020peak,
    title={Peak Estimation and Recovery with Occupation Measures},
    author={Jared Miller and Didier Henrion and Mario Sznaier},
    year={2020},
    eprint={2009.06120},
    archivePrefix={arXiv},
    primaryClass={eess.SY}
}
'''

## Contact
For comments and questions please email [Jared Miller](mailto:miller.jare@northeastern.edu?Subject=CDCS).
