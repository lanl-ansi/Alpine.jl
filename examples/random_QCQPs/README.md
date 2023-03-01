# Random families of non-convex QCQPs 
Scripts for generating a family of random non-convex quadratically-constrained quadratic programs (QCQPs) instances in Section 5 of the paper "Learning to Accelerate the Global Optimization of Quadratically-Constrained Quadratic Programs" by R. Kannan, H. Nagarajan, and D. Deka ([arXiv](https://arxiv.org/abs/2301.00306)). As discussed in the aforementioned paper, these homogeneous families of instances can be very useful for benchmarking machine learning-based models for tuning certain algorithmic parameters of global optimization.

Instances in the above paper were generated using Python 3.8.5 and NumPy 1.19.2.

To generate 1000 random bilinear instances with `N` variables (`N = 10, 20, 50`), run
```
python3 generate_bilinear_instances.py --numVariables N --numInstances 1000
```

To generate 1000 random QCQP instances (including both bilinear and univariate quadratic terms) with `N` variables (`N = 10, 20, 50`), run
```
python3 generate_qcqp_instances.py --numVariables N --numInstances 1000
```

To generate 1000 random pooling instances, run
```
python3 generate_pooling_instances.py --numInstances 1000
```
Additional data files for generating the pooling instances can be found within the `pooling` folder. These data files were generated using the scripts [here](https://github.com/poolinginstances/poolinginstances).


If you find these instances useful in your work, we kindly request that you cite the following paper [arXiv link](https://arxiv.org/abs/2301.00306):
```bibtex
@article{alpine_learning_2022,
  title={Learning to Accelerate the Global Optimization of Quadratically-Constrained Quadratic Programs},
  author={Kannan, Rohit and Nagarajan, Harsha and Deka, Deepjyoti},
  journal={arXiv preprint:2301.00306},
  url={https://arxiv.org/abs/2301.00306},
  year={2022}
}
```