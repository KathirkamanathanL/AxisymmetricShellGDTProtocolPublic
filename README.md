# AxisymmetricShellGDTProtocol
Geometric digital twinning protocol for metal shells of revolution and algorithms for assessing the buckling-relevant fabrication tolerances to EN 1993-1-6 (2025).

# Authors
Mr Lijithan Kathirkamanathan and Dr Adam Jan Sadowski  
Department of Civil and Environmental Engineering, Imperial College London, UK 

Dr Marc Seidel  
Siemens Gamesa Renewable Energy, Germany

# Description
This repository presents a standardised protocol for the systematic processing of reality capture surveys of shells of revolution in point cloud format. The protocol performs enhanced registration, can segmentation, systematic and random outlier removal, surface reconstruction, surface inpainting and projection to arbitrary meshes. The protocol facilitates the post-construction quality assessment, damage and lifetime extension evaluations, and advanced structural analysis with real geometries of shells of revolution. Novel algorithms for automatically assessing EN 1993-1-6 (2025) buckling-relevant fabrication tolerances (dimple, out-of-roundness and unintended eccentricity) are provided. The protocol is illustrated on a real high-definition wind turbine support tower (WTST) laser scan and is readily generalisable to other metal shells of revolution. A link to a 5.7 GB augmented sandbox WTST point cloud dataset to be used in conjunction with this repository can be found below.

# Requirements
Both the geometric digital twinning protocol and tolerance algorithms were developed and tested using MATLAB R2023a. The following MATLAB Add-Ons are required to be installed:
- Computer Vision Toolbox
- Global Optimization Toolbox
- Image Processing Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox

# Usage
1) Clone repository
2) Download the WTST dataset from https://doi.org/10.6084/m9.figshare.28903136.v1 into the Input_ReleaseTower directory and unzip.
3) Run protocolRun.m to generate a geometric digital twin of the sample WTST.
4) Run algorithms in the ToleranceAssessment directory to assess buckling-relevant fabrication tolerances of the sample WTST to EN 1993-1-6 (2025).

# Links
- Sample WTST dataset https://doi.org/10.6084/m9.figshare.28903136.v1
