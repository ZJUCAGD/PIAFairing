# [Curve and Surface Fairing]:Progressive-Iterative Fairing of Curves and Surfaces with Localized Control Point Adjustment

> **📢 Important Notice regarding Manuscript Submission**
>
> This repository contains the official source code for the manuscript titled **"[Progressive-Iterative Fairing of Curves and Surfaces with Localized Control Point Adjustment]"**, currently submitted to **_The Visual Computer_**.
>
> This implementation features a novel **Progressive Iterative Approximation (PIA)** algorithm that fair curves and surfaces by adjusting B-spline control points. Developed in **C++**, it leverages the **Eigen** library for linear algebra and **OpenCASCADE 7.8** for geometric modeling.
>
> To ensure transparency and reproducibility, we provide comprehensive documentation, experimental scripts, and datasets. If you use this code in your research, please cite our paper. Details are provided in the [Citation](#citation) section.

[![OCCT Version](https://img.shields.io/badge/OCCT-7.8.0-orange.svg)]()
[![Eigen Version](https://img.shields.io/badge/Eigen-3.4%2B-green.svg)]()

## 📖 Overview

This project implements the curve and surface fairing by PIA proposed in our paper. The method achieves high-quality fairing by minimizing specific energy functionals (e.g., strain energy) while preserving geometric features, utilizing a Progressive Iterative Approximation (PIA) strategy to dynamically adjust B-spline control points.

### Core Algorithm Implementation
The core logic is encapsulated in the `AlgCurFPIAByAdjustingControlPoint` class and `AlgSurfFPIAByAdjustingControlPoint` . For example, the `AlgCurFPIAByAdjustingControlPoint` class:
```
class AlgCurFPIAByAdjustingControlPoint
{
public:
	AlgCurFPIAByAdjustingControlPoint() {}
	~AlgCurFPIAByAdjustingControlPoint();
	bool Init(const Handle(Geom_BSplineCurve)& curve, int energyType);
	void execute(int maxTimes, double error);
	void execute(const std::vector<size_t>& indexList, const std::vector<double>& weights, int maxTimes);
	void setWeights(const std::vector<size_t>& indexList, const std::vector<double>& weights);
	Handle(Geom_BSplineCurve) getResult();
}
```
*   **Initialization (`Init`)**: Accepts a `Handle(Geom_BSplineCurve)` and sets the energy type.
*   **Global Fairing (`execute`)**: Iteratively adjusts all control points until convergence (max iterations or error threshold).
*   **Local Fairing (`execute` with indexList)**: Supports weighted local fairing for specific control point indices, ideal for feature-preserving scenarios.
*   **Weight Configuration (`setWeights`)**: Allows custom constraint weights for individual control points.
*   **Result Retrieval (`getResult`)**: Returns the faired `Handle(Geom_BSplineCurve)` object.
*   **Equation Validation (`equationCheck`)**: Verifies the correctness of the constructed linear systems (for debugging).

The implementation is tightly integrated with **OpenCASCADE 7.8**'s `Geom_BSplineCurve` data structures and uses **Eigen** for efficient solving of large sparse linear systems.

### Experimental Cases & Data Availability
To validate our algorithm and address the reviewer's comments on data transparency, we have encapsulated a comprehensive suite of benchmark tests within the `CurveAndSurfaceTestCases` class. This class covers scenarios ranging from standard mathematical curves to complex industrial models.

Adhering to the FAIR Principles (Findable, Accessible, Interoperable, Reusable), we have processed the experimental data as follows:
- Open Datasets: All non-confidential synthetic curve data and seed files for automated tests have been organized and uploaded to Zenodo for permanent archiving.
- Restricted Industrial Data: Due to copyright/privacy agreements, the original .step or .igs source file of the car model and phone model used in paper cannot be publicly distributed.