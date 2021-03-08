---
title: "Project Pallette Info"
author: "Aryss Hearne"
date: "3/7/2021"
output: 
  html_document: 
    keep_md: yes
---




```r
library(ggthemes)
library(paletteer)
```

#Pallette Install: One for continious data, one for discrete

```r
#devtools::install_github("johannesbjork/LaCroixColoR")
colors<- LaCroixColoR::lacroix_palette("Pamplemousse", type = "discrete")
colors_cont<-LaCroixColoR::lacroix_palette("Pamplemousse", n = 187, type = "continuous") #as long as the data you're describing is under 187 points long
barplot(rep(1,14), axes=FALSE, col=colors)
```

![](project-palette-info_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
#barplot(rep(1,100), axes=FALSE, col=colors_cont)
```

# Graph Aesthetics:

```r
#plot junk here
#theme_solarized()+
  #scale_color_manual(values=colors)+
  #theme(legend.position="top",
        #axis.text.x=element_text(angle=60, hjust=1))+
#more plot junk here
```

