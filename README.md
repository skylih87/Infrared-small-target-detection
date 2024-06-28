# Infrared Small Target Detection Algorithm Using an Augmented Intensity and Density-Based Clustering

My paper has been accepted by `IEEE Transactions on Geosciene and Remote Sensing`. I will release more interesting works and applications on SIRST soon. Please keep following our repository.

![overall](https://github.com/skylih87/Infrared-small-target-detection/assets/133297940/4a9f919e-1284-4ac5-aedd-a6f61e45d4e8)

## Algorithm Introduction

Infrared Small Target Detection Algorithm Using an Augmented Intensity and Density-Based Clustering 
[Paper](https://github.com/user-attachments/files/16023247/Infrared.Small.Target.Detection.Algorithm.Using.an.Augmented.Intensity.and.Density-Based.Clustering.pdf)

I has proposed a novel algorithm with clustering to achieve accurate single-frame infrared small target detection and develop an open-sourced infrared small target dataset (namely, SNU) in this paper. Experiments on both public and our self-developed datasets demonstrate the effectiveness of my method. The contribution of this paper are as follows:

1. A new IR map is created to improve the detection speed and accuracy of the target by using the standard deviation of IR intensity to increase the contrast of the entire image without using patch-based contrast properties using the local features of the image.

2. Our newly proposed detection algorithm uses density-based clustering, which can accurately recognize geometric object forms or very small objects of 2 × 1 size and classify them as one object. To date, there has been no research on detecting targets using clustering, but multitarget detection can be achieved through clustering.

3. Using the characteristics of the size and IR intensity of the small target, we construct a layered window rather than a sliding window. It is possible to quickly extract a small target from among several candidating objects without a complex equation.

4. We develop 300 single-frame IR images by photographing the drone with a thermal imaging camera. Our dataset is composed of numerous target shapes, various target sizes, and diverse backgrounds.

5. We compare the performance of the proposed method with existing algorithms using a self-generated IR imaging dataset and publicly available datasets. Compared to existing methods, our method is more robust to the variations of complex background, target size, target number, and target shape. In addition, it shows a high detection rate and a low false positive rate while having real-time performance.

## Citation

If you find the code useful, please consider citing our paper using the following BibTeX entry.
```
@article{lee2023infrared,
  title={Infrared small target detection algorithm using an augmented intensity and density-based clustering},
  author={Lee, In Ho and Park, Chan Gook},
  journal={IEEE Transactions on Geoscience and Remote Sensing},
  volume={61},
  pages={1--14},
  year={2023},
  publisher={IEEE}
}
```

## Usage

#### On windows:

```
Click on Main.m and run it. 
```
