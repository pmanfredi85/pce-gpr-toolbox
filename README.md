
# PCE-GPR Toolbox

A MATLAB toolbox based on UQLab for **uncertainty quantification (UQ) using a hybrid Polynomial Chaos Expansion (PCE) and Gaussian Process Regression (GPR) method**.  
Designed to simple and user-friendly for research and educational purposes.

## Author

Paolo Manfredi, Politecnico di Torino, Turin (Italy).

Please report any bugs, comments, questions, or suggestions to paolo.manfredi@polito.it

## 📦 Features

- Hybrid UQ method combining PCE and GPR
- Calculation of PCE coefficients up to an arbitrary order
- Includes of confidence information
- Built-in features for training, prediction, calculation of PCE coefficients, and extraction of statistical information
- Includes live script examples
- Modular codebase with user-facing and internal utilities

## 🧭 Repository Structure

```
uq-gpr-toolbox/
├── src/              # MATLAB source code
│   ├── core/         # Main user-facing functions
│   └── utils/        # Auxiliary/internal functions
├── examples/         # Live scripts (.mlx) and their PDF exports
├── docs/             # User guide and documentation
├── LICENSE
├── README.md
└── CITATION.cff
```

## 📂 Getting Started

1. **Clone this repository**:
   ```bash
   git clone https://github.com/pmanfredi85/pce-gpr-toolbox.git
   ```
2. **Add the source folder to your MATLAB path**:
   ```matlab
   addpath(genpath('path_to_repo/src'))
   ```

3. **Explore examples**:
   Open the `.mlx` live scripts in the `examples/` folder to get started.

## 🧪 Requirements

- MATLAB R2023b or newer (recommended)
- UQLab release 2.0.0 (recommended)
- Generalized chi-square distribution toolbox version 2.3.0 (recommended)
- Statistics and Machine Learning Toolbox
- Optimization Toolbox

## 🧩 Third-Party Components

This toolbox includes components from the following external projects:

- [UQLab]: a MATLAB toolbox licensed under the [BSD 3-Clause License](LICENSE)
- [gx2]: a MATLAB toolbox licensed under the [MIT License](3rdparty/MIT_LICENSE.txt)

The relevant files are used in accordance with their license terms. See the full license text in `3rdparty/MIT_LICENSE.txt`.


## 📄 Documentation

The user guide is available as a PDF in the [`docs/`](docs/) folder.

## 🧾 Citation

If you use this toolbox in your research, please cite:

> P. Manfredi, *A hybrid polynomial chaos expansion–Gaussian process regression method for Bayesian uncertainty quantification and sensitivity analysis*, Computer Methods in Applied Mechanics and Engineering, vol. 436, no. 117693 (2025). DOI: https://doi.org/10.1016/j.cma.2024.117693

> P. Manfredi, *Polynomial chaos vs kernel machine learning methods for uncertainty quantification: A comparative study and benchmarking with the hybrid PCE-GPR approach*, Computer Methods in Applied Mechanics and Engineering, vol. 449, part A, no. 118523 (2026). DOI: https://doi.org/10.1016/j.cma.2025.118523

> P. Manfredi, *PCE-GPR: A toolbox for high-dimensional uncertainty quantification*, Zenodo (2025). DOI: https://doi.org/10.5281/zenodo.15348860

Citation metadata is also available in [`CITATION.cff`](CITATION.cff).

## 🔗 Related Resources

- 📘 [UQLab](https://www.uqlab.com): The original uncertainty quantification framework on which this toolbox is based.
- 📘 [gx2](https://github.com/abhranildas/gx2/releases/tag/v2.3.0): The toolbox to support the generalized chi-square distribution.
- 📚 Related publication (theory): https://doi.org/10.1016/j.cma.2024.117693
- 📚 Related publication (benchmarking): https://doi.org/10.1016/j.cma.2025.118523
- 💾 Benchmark datasets: https://doi.org/10.5281/zenodo.15349691

## 🛠 License

This project is licensed under the [BSD 3-Clause License](LICENSE).
