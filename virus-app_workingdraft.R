#install libraries
if (!require("tidyverse")) install.packages('tidyverse')
if (!require("here")) install.packages('here')
if (!require("janitor")) install.packages('janitor')
if (!require("shiny")) install.packages('shiny')
if (!require("shinydashboard")) install.packages('shinydashboard')
if (!require("ggthemes")) install.packages('ggthemes')
if (!require("paletteer")) install.packages('paletteer')
if (!require("bslib")) install.packages('bslib')
#if (!require("johannesbjork/LaCroixColoR")) devtools::install_github("johannesbjork/LaCroixColoR")
# libraries
library(tidyverse)
library(here)
library(janitor)
library(shiny)
library(shinydashboard)
library(ggthemes)
library(paletteer)
library(bslib)
# included data
virus<-readr::read_csv("data/viruses.csv")
#include palette
colors<- LaCroixColoR::lacroix_palette("Pamplemousse", type = "discrete")
moth<-wesanderson::wes_palette("GrandBudapest2",12,"continuous")
# clean data
virus_app<-virus%>%
  separate('Host', into=c("host_type_1", "host_type_2", "host_type_3"), sep=",")%>%
  pivot_longer(host_type_1:host_type_3,
               names_to="host_num",
               values_to="host_type",
               values_drop_na=T)%>%
  separate('Organism Groups', into=c("domain", "kingdom", "subgroup"), sep=";")%>%
  filter(Level=="Complete" & (host_type=="human" | host_type=="land plants"))%>%
  select('host_type', 'Organism Name', 'kingdom','subgroup', 'Size(Mb)', 'GC%', 'Genes','Level')%>%
  rename(organism_name='Organism Name', kingdom=kingdom, subgroup=subgroup, level="Level", size_mb='Size(Mb)', perc_gc='GC%', host=host_type, genes_num='Genes')%>%
  mutate(gene_to_genome_ratio=genes_num/size_mb, na.rm=T)%>%
  arrange(host)




# Introduce App: Fluid Page
ui <- fluidPage(
  #aesthetic
  theme = bs_theme(version = 4, bootswatch = "minty"),
  #title
  titlePanel("Viral Data"),
  #sidebar
  sidebarLayout(
    #sidebar panel 1
    sidebarPanel(
        sliderInput("jitterw", label = "Jitter Width", min=0.5, max=40, value=12, step=5),
        sliderInput("jitterh", label = "Jitter height", min=0, max=.2, value=0, step=.01)
        ),#close panel
    # create a spot for the barplot
    mainPanel(
      tabsetPanel(
        tabPanel("Exploratory Analysis", plotOutput("plot", brush="plot_brush", height=800, width = 1100),
                 verbatimTextOutput("info")
        ),#close tab zero
        tabPanel("Exploratory Analysis Adjusted", plotOutput("plot1", brush="plot_brush", height=800, width = 1100),
        verbatimTextOutput("info1")
        ),#close tab one
        tabPanel("Exploratory Analysis by Kingdom", plotOutput("plot1_5", brush="plot_brush", height=800, width = 1100),
                 verbatimTextOutput("info1_5")
        ),#close tab one point five
        tabPanel("Gene to Genome Ratio Overall", plotOutput("plot2", brush="plot_brush", height=800, width = 1100),
                 verbatimTextOutput("info2")
                 )#closes tab two
    )#closes tabset
      )#close sidebar
    )#end mainPanel  
    )#end fluid Page
    

#server
server <- function(input, output, session) {
  #end clean
  session$onSessionEnded(stopApp)
  #plot info
  output$plot<-renderPlot({
    #plot parameters
    ggplot(virus_app, aes(x=genes_num, y=size_mb, size=perc_gc, color=host, group=host))+
      geom_jitter(alpha=0.25, width=input$jitterw, height=input$jitterh)+
      #aesthetics
      theme_solarized()+
      scale_color_manual(values=colors)+
      theme(legend.position="top",
            axis.text.x=element_text(angle=60, hjust=1))+
      scale_size(range = c(0.1, 10),
                 guide = "none")+
      #labels
      labs(x="Number of Genes",
           y="Genome Size (Mb)",
           size="Percent GC",
           color="Host")+
      #text
      geom_text(aes(x = genes_num, label = organism_name),
                color = "grey50",
                data = filter(virus_app, size_mb > .23 | genes_num >300 | organism_name %in% c("Zika virus", "Human betaherpesvirus", "Horsepox virus")))
  })#end plot0
  #plot info
  output$plot1<-renderPlot({
    #plot parameters
    ggplot(virus_app, aes(x=genes_num, y=size_mb, size=perc_gc, color=host, group=host))+
      geom_jitter(alpha=0.25, width=input$jitterw, height=input$jitterh)+
      scale_x_log10()+
      scale_y_log10()+
      #aesthetics
      theme_solarized()+
      scale_color_manual(values=colors)+
      theme(legend.position="top",
            axis.text.x=element_text(angle=60, hjust=1))+
      scale_size(range = c(0.1, 10),
                 guide = "none")+
      #labels
      labs(x="Number of Genes",
           y="Log of Genome Size (Mb)",
           size="Log of Percent GC",
           color="Host")+
      #text
      geom_text(aes(x = genes_num, label = organism_name),
                color = "grey50",
                data = filter(virus_app, size_mb > .23 | genes_num >300 | organism_name %in% c("Zika virus", "Human betaherpesvirus", "Horsepox virus")))
  })#end plot1
  #plot info 1.5
  output$plot1_5<-renderPlot({
    #plot parameters
    ggplot(virus_app, aes(x=genes_num, y=size_mb, size=perc_gc, color=host, group=host))+
      geom_jitter(alpha=0.25, width=input$jitterw, height=input$jitterh)+
      facet_grid(. ~kingdom)+
      #aesthetics
      theme_solarized()+
      scale_color_manual(values=colors)+
      theme(legend.position="top",
            axis.text.x=element_text(angle=60, hjust=1))+
      scale_size(range = c(0.1, 10),
                 guide = "none")+
      #labels
      labs(x="Number of Genes",
           y="Genome Size (Mb)",
           size="Percent GC",
           color="Host")+
      #text
      geom_text(aes(x = genes_num, label = organism_name),
                color = "grey50",
                data = filter(virus_app, size_mb > .23 | genes_num >300 | organism_name %in% c("Zika virus", "Human betaherpesvirus", "Horsepox virus")))
  })#end plot1
  output$plot2<-renderPlot({
    #plot2 parameters
    virus_app%>%
      ggplot(aes(x=host, y=gene_to_genome_ratio, fill=host, color=kingdom, size=perc_gc))+
      geom_jitter(alpha=0.02, height=input$jitterh)+
      geom_boxplot(varwidth = T, outlier.size = 3.5)+
      facet_grid(. ~kingdom)+
      #aesthetics for plot 2
      theme_solarized()+
      scale_color_manual(values=moth)+
      scale_fill_manual(values=colors)+
      theme(legend.position="top",
            axis.text.x=element_text(angle=60, hjust=1))+
      scale_size(range = c(0.1, 10),
                 guide = "none")+
      #labels for plot 2
      labs(x="Host",
           y="Gene to Genome Ratio",
           color="Host",
           Fill="Kingdom")
  })#end plot 2
  #output point and click
  output$info<-renderPrint({
    brushedPoints(virus_app, input$plot_brush)
  })#end info panel one
  #output point and click
  output$info1<-renderPrint({
    brushedPoints(virus_app, input$plot_brush)
  })#end info panel one
  #output point and click one point five
  output$info1_5<-renderPrint({
    brushedPoints(virus_app, input$plot_brush)
  })#end info panel one point five
  #output point and click two
  output$info2<-renderPrint({
    brushedPoints(virus_app, input$plot_brush)
  })#end info panel two

  }

shinyApp(ui, server)
