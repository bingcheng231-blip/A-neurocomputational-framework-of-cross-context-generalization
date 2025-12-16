# Scripts for A neurocomputational framework of cross-context generalization: dynamic representational geometry in high-level visual cortex and dual coding in vmPFC
# Authors:
# Bincheng Wen1,2,3, Chuncheng Zhang1, Changde Du1, Le Chang2,3*, Huiguang He1,3*
# 1 Laboratory of Brain Atlas and Brain-Inspired Intelligence, Key Laboratory of Brain Cognition and Brain-inspired Intelligence Technology, Institute of Automation, Chinese Academy of Sciences, Beijing 100190, China.
# 2 Institute of Neuroscience, Key Laboratory of Brain Cognition and Brain-inspired Intelligence Technology, CAS Center for Excellence in Brain Science and Intelligence Technology, Chinese Academy of Sciences, Shanghai 200031, China.
# 3 University of Chinese Academy of Sciences, Beijing 100049, China.
# *Correspondence: lechang@ion.ac.cn; huiguang.he@ia.ac.cn


In general, if anything across the codes is missing or unlclear, you are most welcome to contact Huiguang He, always happy to chat and explain whatever needed! 


## General description:
The primary analyses and result figures from the paper can be reproduced by running the Main_result.m script.
To run the script correctly, please add the functions in the ./analyses_function folder and the ./analyses_function/utilities folder to the working path

Note: because we did not seed the randomization in the original analyses, some of the results are not 100% identical (due to randomized updampling in the decoding analyses). 

## Data availability
The data provided is intended to facilitate the reproduction of results at different levels within the main analysis.
The data required for the main analysis is located in the https://www.scidb.cn/s/jqiuY3 . 


### How to get all data used in this paper?
The behavior and neuroimage data that support the findings of this study are available upon request from the corresponding author by e-mail after the paper is accepted.


# Software Versions 
## Graphpad 
GraphPad Prism 8.0.2
## MATLAB
R2023a
### Packages for preprocessing
SPM12 (12)
### Packages for fMRI, MEG and behavior task
Psychtoolbox 3
## Python
pythorch 1.13.1
numpy  1.26.4
scipy  1.15.3
matplotlib 3.8.4


## Acknowledgments:
This work was supported in part by the Strategic Priority Research Program of the Chinese Academy of Sciences (grant no. XDB1010202), the National Natural Science Foundation of China under grant nos. 62020106015 and 62206284, the Beijing Natural Science Foundation under grant no. L243016 and the Beijing Nova Program under grant no. 20230484460. We would like to thank E. J. Allen and K. Kay for sharing the NSD fMRI data.

