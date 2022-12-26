# IP4GS
<a href="https://www.r-project.org/" target="_blank"><img src="https://img.shields.io/badge/language-R-orange?style=plastic"></a>
<a href="https://cran.r-project.org/bin/windows/base/old/" target="_blank"><img src="https://img.shields.io/badge/R%20version-%3E%3D%203.6.0-orange?style=plastic"></a>
<a href="https://shiny.rstudio.com/" target="_blank"><img src="https://img.shields.io/badge/Shiny-WebApp-blue?style=plastic"></a>
<a href="https://ngdc.cncb.ac.cn/ip4gs/" target="_blank"><img src="https://img.shields.io/badge/webpage-ready-green?style=plastic"></a>
![](https://img.shields.io/badge/platform-Win%20%7C%20Linux%20%7C%20MacOS-lightgrey?style=plastic)<br/>


[IP4GS](https://ngdc.cncb.ac.cn/ip4gs/), an ***i***nteractive ***p***latform ***for*** ***g***enomic ***s***election, offers a user-friendly interface for performing streamlined GS analysis involving the above-mentioned steps simply through point-and-click actions. IP4GS currently includes seven commonly used models, 11 evaluation metrics, and visualization modules, offering great convenience for plant breeders with limited bioinformatics knowledge to apply GS analysis.

### Quick start
#### 1. Upload data
<div align="center">
<img src="https://github.com/furan2019/IP4GSdata/blob/main/pic/01_upload.jpg" alt="drawing" width="1200"/>
</div>

#### 2. Genotypic data filtration and input data preview
<div align="center">
<img src="https://github.com/furan2019/IP4GSdata/blob/main/pic/02_datapreprocessing.jpg" alt="drawing" width="1200"/>
</div>

#### 3. Data visualization console and dimension reduction console
<div align="center">
<img src="https://github.com/furan2019/IP4GSdata/blob/main/pic/03_datashow.jpg" alt="drawing" width="1200"/>
</div>

#### 4. Preview of DR data
<div align="center">
<img src="https://github.com/furan2019/IP4GSdata/blob/main/pic/04_DR_showDownload.jpg" alt="drawing" width="1200"/>
</div>

#### 5. Genotype-to-phenotype prediction
<div align="center">
<img src="https://github.com/furan2019/IP4GSdata/blob/main/pic/05_G2Pprediction.jpg" alt="drawing" width="1200"/>
</div>

#### 6. Viewing scatter plots of observed and predicted phenotypes
<div align="center">
<img src="https://github.com/furan2019/IP4GSdata/blob/main/pic/06_visualizaiton.jpg" alt="drawing" width="1200"/>
</div>

## Local usage
#### 1. Download UI and server script of IP4GS
#### 2. source UI and sever script
#### 3. install all dependence package in code
#### 4. library all packages and run IP4GS
```R
shinyApp(options = list(launch.browser=T),ui = UI,server = server)
```
