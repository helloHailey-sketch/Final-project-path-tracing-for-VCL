# Introduction to Visual Computing and Interaction

This repository contains my final project for the course **Introduction to Visual Computing and Interaction** at Peking University (PKU).

## Course Homepage

[Introduction to Visual Computing and Interaction - PKU] （https://vcl.pku.edu.cn/course/vci）

## Final Project: Path Tracing

Unlike the Whitted-style ray tracing implemented in Lab 3, path tracing is a global illumination rendering method based on the rendering equation. It accurately simulates light transport by tracing many random paths of light as they bounce and scatter throughout the scene, allowing for effects such as soft shadows, depth of field, caustics, and indirect lighting. Path tracing forms the foundation of modern physically-based rendering (PBR) techniques and is widely used in photorealistic image synthesis.

In this final project, I have successfully reproduced a basic path tracing algorithm. Additionally, I implemented BVH (Bounding Volume Hierarchy) acceleration to improve rendering efficiency. As an extra creative extension, I further reproduced sampling perturbations to achieve a painterly effect, giving the rendered images an artistic, hand-painted style.
