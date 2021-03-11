#install libraries
if (!require("tidyverse")) install.packages('tidyverse')
if (!require("here")) install.packages('here')
if (!require("janitor")) install.packages('janitor')
if (!require("shiny")) install.packages('shiny')
if (!require("shinydashboard")) install.packages('shinydashboard')
if (!require("ggthemes")) install.packages('ggthemes')
if (!require("paletteer")) install.packages('paletteer')
if (!require("bslib")) install.packages('bslib')
if (!require("devtools")) install.packages('devtools')
if (!require("shinycssloaders")) install.packages('shinycssloaders')
devtools::install_github("johannesbjork/LaCroixColoR")
# libraries
library(tidyverse)
library(here)
library(janitor)
library(shiny)
library(shinydashboard)
library(ggthemes)
library(paletteer)
library(bslib)
library(shinycssloaders)
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
  #Navbar structure
  navbarPage("Viral Data", theme = bs_theme(version = 4, bootswatch = "minty"),

    # create a spot for the barplot
  navbarMenu("Exploratory Analysis",
             tabPanel("Data Unconstrained", fluid=T,
                      titlePanel("An Exploratory Analysis"),
                      fluidRow(
                        column(4,
                               sliderInput("zoomx", 
                                           label = "Zoom X Axis", 
                                           min=0, 
                                           max=400, 
                                           value=c(0,400))),
                        column(4,
                               sliderInput("zoomy", 
                                           label = "Zoom Y Axis", 
                                           min=0, 
                                           max=.25, 
                                           value=c(0,.25))),
                        fluidRow(
                          column(6,
                                 withSpinner(type =6, plotOutput(outputId = "plot",
                                            brush="plot_brush"
                                            ))),
                          hr(),
                          helpText("words here"),
                          br()
                        ))),
             tabPanel("Data Logged",fluid=T,
                      titlePanel("A Re-examination"),
                      fluidRow(
                          column(6,
                                 withSpinner(type =6, plotOutput(outputId = "plot1",
                                                        brush="plot_brush"
                                 ))),
                          hr(),
                          helpText("words here"),
                          br()
                        )),
             tabPanel("Data by Kingdom", fluid=T,
                      titlePanel("And by Kingdom"),
                      fluidRow(
                        column(4,
                               sliderInput("zoomx", 
                                           label = "Zoom X Axis", 
                                           min=0, 
                                           max=400, 
                                           value=c(0,400))),
                        column(4,
                               sliderInput("zoomy", 
                                           label = "Zoom Y Axis", 
                                           min=0, 
                                           max=.25, 
                                           value=c(0,.25))),
                      fluidRow(
                        column(6,
                               withSpinner(type =6, plotOutput(outputId = "plot1_5",
                                                      brush="plot_brush"
                               ))),
                        hr(),
                        helpText("words here"),
                        br())),
)#close navbar options
  ),#closes navbar Menu Part 1

            tabPanel("Mean of Genomic Data", fluid=T,
                      titlePanel("Mean of Genomic Data"),
                      fluidRow(
                        column(6,
                               withSpinner(type =6, plotOutput(outputId = "plot2",
                                                      brush="plot_brush"
                                                      ))),
                        column(6,
                               dataTableOutput(outputId = "virustable")
                               )))
          #close navbar Menu part 2
)#closes navbar page
)#closes ui/ fluid page

 

#server
server <- function(input, output, session) {
  #end clean
  session$onSessionEnded(stopApp)
  #plot info
  output$plot<-renderPlot({
    #plot parameters
    ggplot(virus_app, aes(x=genes_num, y=size_mb, size=perc_gc, color=host, group=host))+
      geom_point(alpha=0.25)+
      xlim(input$zoomx)+
      ylim(input$zoomy)+
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
      geom_point(alpha=0.25)+
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
      labs(x="Log of Number of Genes",
           y="Log of Genome Size (Mb)",
           size="Percent GC",
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
      geom_point(alpha=0.25)+
      facet_grid(. ~kingdom)+
      xlim(input$zoomx)+
      ylim(input$zoomy)+
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
      geom_point(alpha=0.25)+
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
viralbrush<-reactive({
  user_viralbrush<-input$plot_brush
  brushedPoints(virus_app,user_viralbrush)
})#end brush paramaters

output$virustable<-DT::renderDataTable({
  DT::datatable(unique(viralbrush()[,c("organism_name","host")]),
                colnames=c("Sort"="organism_name", "host"),
                rownames=F
  )
})
  }

shinyApp(ui, server)


#To do: put in zoom slider
#put in citation
#put in hypothesis text