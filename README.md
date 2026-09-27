# MATLAB code for the numerical experiments

This repository contains the MATLAB code for the numerical experiments in the preprint [A mesh- and horizon-robust preconditioner for finite element discretizations of volume-constrained nonlocal diffusion problems](https://arxiv.org/abs/2609.24143) by Jiashu Lu, Yufeng Nie, Haolun Zhang, and Pingrui Zhang.

## Requirements

The code was tested with MATLAB R2025b. No external data or third-party MATLAB package is required. The computations reported in the paper were performed on a MacBook Pro with an Apple M5 Max processor and 128 GB of memory.

The full run of `main_2.m` includes a problem with 17,153 degrees of freedom. Its assembly uses a dense auxiliary accumulator of approximately 2.35 GB, and the spectral computations require additional memory. A high-memory machine is therefore recommended. The `quick` mode is intended only to check the installation and does not reproduce the complete paper tables.

## Experiments

All MATLAB files are stored in the same directory. Start MATLAB in this directory and run:

```matlab
main_1
main_2
main_3
```

- `main_1.m` reproduces the finite element convergence experiment in Section 6.1 and writes Table 6.1 to `results/section_6_1/table_6_1.csv`.
- `main_2.m` reproduces the mesh- and horizon-robustness experiment in Section 6.2. It writes Tables 6.2--6.3 and the detailed PCG and spectral data to `results/section_6_2/`, and generates Figure 6.1 as PDF and PNG files.
- `main_3.m` reproduces the kernel-robustness experiment in Section 6.3. It writes Table 6.4 and the detailed PCG and spectral data to `results/section_6_3/`.

The illustrative plots of the five kernels in Figure 6.2 are not included because they are not numerical results.

For a short installation check, run:

```matlab
main_1('quick')
main_2('quick')
main_3('quick')
```

## Output and restart

The programs create `results/` and `checkpoints/` automatically. Sections 6.2 and 6.3 save one checkpoint for each completed parameter case, so an interrupted full run can be continued by running the same main function again. Delete the corresponding checkpoint directory to recompute an experiment from the beginning.

The unstructured hierarchy and the three test loads use the fixed seed `20260821`. Timings depend on the computer, while the reported spectral quantities and PCG iteration counts are reproducible up to normal floating-point variation.

## Citation

If this code is useful in your work, please cite the preprint:

```bibtex
@article{lu2026mesh,
  title   = {A mesh- and horizon-robust preconditioner for finite element discretizations of volume-constrained nonlocal diffusion problems},
  author  = {Lu, Jiashu and Nie, Yufeng and Zhang, Haolun and Zhang, Pingrui},
  journal = {arXiv preprint arXiv:2609.24143},
  year    = {2026}
}
```
