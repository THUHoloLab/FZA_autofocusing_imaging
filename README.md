## README
FZA_autofocus imaging is a project to achieve the autofocus in the FZA lensless imaging through the sharpness valuation **WTN** and get a clearer imaging by the **ADMM** algorithm.

## Reference
Fangyu Liu, Jiachen Wu, and Liangcai Cao, "Autofocusing of Fresnel zone aperture lensless imaging for QR code recognition," Optics Express 31, 15889-15903 (2023).
https://doi.org/10.1364/OE.489157
https://opg.optica.org/oe/fulltext.cfm?uri=oe-31-10-15889

#### WTN
- Before running this programm, confirm whether the parameters are suitable for your experiment if you have changed the input. 
- The parameters you should pay attention
  - di: the distance between mask and sensor(mm)
  - dp: the size of the pixels
  - r1: the FZA parameter which is the radius of the central circle
  - theta: the weight of the gnorm in the WTN
- Running this programm, you can get the measurement image, the WTN autofocus curve and the BP reconstruction in the estimated focus distance.

#### ADMM
- Pay atttention to uniform the parameters between the WTN and ADMM
- After running the WTN, run the ADMM directly or record the f_distance and original image I then run the ADMM
- The parameters you can change to have different effect
  - lambda
  - mu
  - eta
  - beta
  - N
- Running this program, you can get the dynamic change display through the proess of ADMM and finally have a ADMM reconstruction image in the estimated focus distance.

#### Function
- FZA lensless imaging BP function
  - My...
- sharpness valuation metric 
  - GNORM
  - GRA
  - LAP
  - SMD
  - ToG
  - VAR
- filtering function
  - Gauss_1

#### source
- the experiment orignal measurement image

#### contact
- Fangyu Liu, liu-fy21@mails.tsinghua.edu.cn
- Jiachen Wu, wjc2022@tsinghua.edu.cn
