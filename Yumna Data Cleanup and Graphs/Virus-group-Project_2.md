---
title: "group project excel join"
author: "Aryss Hearne"
date: "2/26/2021"
output: 
  html_document: 
    keep_md: yes
---






```r
library(tidyverse)
```

```
## -- Attaching packages --------------------------------------- tidyverse 1.3.0 --
```

```
## v ggplot2 3.3.3     v purrr   0.3.4
## v tibble  3.0.6     v dplyr   1.0.4
## v tidyr   1.1.2     v stringr 1.4.0
## v readr   1.4.0     v forcats 0.5.1
```

```
## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

```r
library(here)
```

```
## here() starts at C:/Users/yamou/Desktop/Bis15LBuds
```

```r
library(janitor)
```

```
## 
## Attaching package: 'janitor'
```

```
## The following objects are masked from 'package:stats':
## 
##     chisq.test, fisher.test
```

```r
library(naniar)
library(stringr)
library(ggthemes)
library(paletteer)
library(RColorBrewer)
library(ggplot2)
```


```r
virus <- read_csv(here("Yumna Data Cleanup and Graphs", "Data", "viruses.csv")) 
```

```
## 
## -- Column specification --------------------------------------------------------
## cols(
##   `Organism Name` = col_character(),
##   `Organism Groups` = col_character(),
##   BioSample = col_logical(),
##   BioProject = col_character(),
##   Assembly = col_character(),
##   Level = col_character(),
##   `Size(Mb)` = col_double(),
##   `GC%` = col_double(),
##   Host = col_character(),
##   CDS = col_double(),
##   Neighbors = col_double(),
##   `Release Date` = col_datetime(format = ""),
##   `GenBank FTP` = col_character(),
##   `RefSeq FTP` = col_character(),
##   Genes = col_double(),
##   Scaffolds = col_double()
## )
```

```
## Warning: 344 parsing failures.
##  row       col           expected       actual                                                                               file
## 1438 BioSample 1/0/T/F/TRUE/FALSE SAMN02981359 'C:/Users/yamou/Desktop/Bis15LBuds/Yumna Data Cleanup and Graphs/Data/viruses.csv'
## 4401 BioSample 1/0/T/F/TRUE/FALSE SAMN02981224 'C:/Users/yamou/Desktop/Bis15LBuds/Yumna Data Cleanup and Graphs/Data/viruses.csv'
## 8940 BioSample 1/0/T/F/TRUE/FALSE SAMN01137200 'C:/Users/yamou/Desktop/Bis15LBuds/Yumna Data Cleanup and Graphs/Data/viruses.csv'
## 8941 BioSample 1/0/T/F/TRUE/FALSE SAMN01137212 'C:/Users/yamou/Desktop/Bis15LBuds/Yumna Data Cleanup and Graphs/Data/viruses.csv'
## 8944 BioSample 1/0/T/F/TRUE/FALSE SAMN01137140 'C:/Users/yamou/Desktop/Bis15LBuds/Yumna Data Cleanup and Graphs/Data/viruses.csv'
## .... ......... .................. ............ ..................................................................................
## See problems(...) for more details.
```


### Data Tidying starts here ###


```r
virus2<-virus%>%
  separate(`Organism Groups`, into = c("viruss","other","family"),sep=";") %>% 
  select(-"viruss",-"other")
virus2 <- janitor::clean_names(virus2)
```



```r
virus3<-virus2 %>% 
  separate(`host`, into = c("host1","host2","host3","host4", "host5"),sep=",")
```

```
## Warning: Expected 5 pieces. Missing pieces filled with `NA` in 29949 rows [1, 2,
## 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, ...].
```

```r
newv <- virus2 %>% 
  separate_rows(host)
```



```r
virus3<- virus3 %>%
  mutate(strain=word(organism_name, 1L,2, sep=fixed(" ")))

new<- as.data.frame(table(virus3$strain)) %>% 
  #filter(Freq==10) %>% 
  rename(strain=Var1)

virus4 <- full_join(virus3, new, by="strain")

virus4$host_f <- paste(virus4$host1, virus4$host2,virus4$host3,virus4$host4,virus4$host5, sep=",")

virus4 <- virus4 %>% 
  separate_rows(host_f, sep=",") %>% 
  filter(host_f!="NA")

virus4<-virus4 %>%
  filter(host_f!="plants") %>% 
  filter(host_f!="diatom") %>%
  mutate(size_to_gene_ratio= ifelse(genes>0,size_mb/genes,NA )) %>% 
  mutate(number_of_hosts=ifelse(is.na(host1), 0, ifelse(is.na(host2), 1, 
                                       ifelse(is.na(host3), 2,
                                              ifelse(is.na(host4), 3, 
                                                     ifelse(is.na(host5), 4, 5)))))) %>% 
  select(organism_name, family, strain, size_mb, gc_percent, genes, size_to_gene_ratio, number_of_hosts, host_f, host1, host2, host3, host4, host5, Freq)

virus4$host_f <- gsub(" vertebrates", "vertebrates", virus4$host_f)
virus4$host_f <- gsub(" human", "human", virus4$host_f)
virus4$host_f <- gsub("eukaryotic algae", "algae", virus4$host_f)

virus4$number_of_hosts <- as.factor(virus4$number_of_hosts)
```


### Palette Install: One for continious data, two for discrete ###

```r
#devtools::install_github("johannesbjork/LaCroixColoR")

colors<- LaCroixColoR::lacroix_palette("Pamplemousse", type = "discrete")
colors_cont<-LaCroixColoR::lacroix_palette("Pamplemousse", n = 14, type = "continuous")
colors_cont_5<-LaCroixColoR::lacroix_palette("Pamplemousse", n = 5, type = "discrete")
barplot(rep(1,14), axes=FALSE, col=colors)
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
barplot(rep(1,100), axes=FALSE, col=colors_cont)
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

### Graphs plotted that were used in presentation ###


```r
#png('Counts of Number of Hosts.png')
virus4 %>% 
  ggplot(aes(x=number_of_hosts, fill=number_of_hosts))+
  geom_bar()+
  theme_solarized()+
  scale_fill_manual(values=colors)+
  theme(legend.position = "none", axis.text.x=element_text(hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Counts of Number of Hosts",
       x="Number of Hosts",
       y="Count")
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
#dev.off()


#png('Number of Hosts vs Size to Gene Ratio.png')
virus4 %>% 
  ggplot(aes(x=number_of_hosts, y=size_to_gene_ratio, group=number_of_hosts, fill=number_of_hosts))+
  geom_col()+
  theme_solarized()+
  scale_fill_manual(values=colors)+
  theme(legend.position = "none", axis.text.x=element_text(hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Number of Hosts vs Size to Gene Ratio",
       x="Number of Hosts",
       y="Size to Gene Ratio")
```

```
## Warning: Removed 17489 rows containing missing values (position_stack).
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-7-2.png)<!-- -->

```r
#dev.off()

#png('Number of Hosts vs Total Number of Genes.png')
virus4 %>% 
  ggplot(aes(x=number_of_hosts, y=genes, group=number_of_hosts, fill=number_of_hosts))+
  geom_col()+
  theme_solarized()+
  scale_x_discrete()+
  theme_solarized()+
  scale_fill_manual(values=colors_cont_5)+
  theme(legend.position = "none", axis.text.x=element_text(hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Number of Hosts vs Total Number of Genes",
       x="Number of Hosts",
       y="Total Number of Genes")
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-7-3.png)<!-- -->

```r
#dev.off()

#png('Number of Hosts vs Total Size of Genome.png')
virus4 %>% 
  ggplot(aes(x=number_of_hosts, y=size_mb, group=number_of_hosts, fill=number_of_hosts))+
  geom_col()+
  scale_x_discrete()+
  theme_solarized()+
  scale_fill_manual(values=colors_cont_5)+
  theme(legend.position = "none", axis.text.x=element_text(hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Number of Hosts vs Total Size of Genome",
       x="Number of Hosts",
       y="Total Size of Genome")
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-7-4.png)<!-- -->

```r
#dev.off()

#png('Host Identity vs Size to Gene Ratio.png')
virus4 %>% 
  filter(!is.na(host_f)) %>%
  filter(!is.na(size_to_gene_ratio)) %>% 
  ggplot(aes(x=host_f, y=size_to_gene_ratio, group=host_f,fill=host_f))+
  geom_col()+
  theme_solarized()+
  scale_fill_manual(values=colors_cont, name = "Host")+
  theme(legend.position="top",
        axis.text.x=element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Host Identity vs Size to Gene Ratio",
       x="Host",
       y="Size to Gene Ratio")
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-7-5.png)<!-- -->

```r
#dev.off()

#png('Host Identity vs Total Number of Genes.png')
virus4 %>% 
  filter(!is.na(host_f)) %>% 
  filter(genes!=0) %>% 
  ggplot(aes(x=host_f, y=genes, group=host_f,fill=host_f))+
  geom_col()+
  theme_solarized()+
  scale_fill_manual(values=colors_cont, name = "Host")+
  theme(legend.position="top",
        axis.text.x=element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Host Identity vs Total Number of Genes",
       x="Host",
       y="Total Number of Genes")
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-7-6.png)<!-- -->

```r
#dev.off()

#png('Host Identity vs Total Size of Genome.png')
virus4 %>% 
  filter(!is.na(host_f)) %>% 
  ggplot(aes(x=host_f, y=size_mb, group=host_f, fill=host_f))+
  geom_col()+
  theme_solarized()+
  scale_fill_manual(values=colors_cont, name = "Host")+
  theme(legend.position="top",
        axis.text.x=element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Host Identity vs Total Size of Genome",
       x="Host",
       y="Total Size of Genome")
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-7-7.png)<!-- -->

```r
#dev.off()
```


### Some extraneous data exploration that did not make it into the presentation ###


```r
virus4 %>% 
  top_n(1000,size_mb) %>% 
  arrange(desc(size_mb)) %>% 
  ggplot(aes(x=number_of_hosts))+
  geom_bar()+
  theme_solarized()+
  scale_color_manual(values=colors)+
  theme(legend.position="top",
        axis.text.x=element_text(angle=60, hjust=1))
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
virus4 %>% 
  top_n(1000,size_mb) %>% 
  arrange(desc(size_mb)) %>% 
  ggplot(aes(x=host_f))+
  geom_bar()+
  theme_solarized()+
  scale_color_manual(values=colors)+
  theme(legend.position="top",
        axis.text.x=element_text(angle=60, hjust=1))
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```r
virus4 %>% 
  top_n(1000,size_to_gene_ratio) %>% 
  arrange(desc(size_to_gene_ratio)) %>% 
  ggplot(aes(x=number_of_hosts))+
  geom_bar()+
  theme_solarized()+
  scale_color_manual(values=colors)+
  theme(legend.position="top",
        axis.text.x=element_text(angle=60, hjust=1))
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

```r
virus4 %>% 
  top_n(1000,size_to_gene_ratio) %>% 
  arrange(desc(size_to_gene_ratio)) %>% 
  ggplot(aes(x=host_f))+
  geom_bar()+
  theme_solarized()+
  scale_color_manual(values=colors)+
  theme(legend.position="top",
        axis.text.x=element_text(angle=60, hjust=1))
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-8-4.png)<!-- -->

```r
virus4 %>% 
  top_n(1000,genes) %>% 
  arrange(desc(genes)) %>% 
  ggplot(aes(x=number_of_hosts))+
  geom_bar()+
  theme_solarized()+
  scale_color_manual(values=colors)+
  theme(legend.position="top",
        axis.text.x=element_text(angle=60, hjust=1))
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-8-5.png)<!-- -->

```r
virus4 %>% 
  top_n(1000,genes) %>% 
  arrange(desc(genes)) %>% 
  ggplot(aes(x=host_f))+
  geom_bar()+
  theme_solarized()+
  scale_color_manual(values=colors)+
  theme(legend.position="top",
        axis.text.x=element_text(angle=60, hjust=1))
```

![](Virus-group-Project_2_files/figure-html/unnamed-chunk-8-6.png)<!-- -->


