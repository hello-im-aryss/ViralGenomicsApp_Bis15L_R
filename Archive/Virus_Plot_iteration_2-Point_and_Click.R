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
# define palette
colors<- LaCroixColoR::lacroix_palette("Pamplemousse", type = "discrete")
# Clean Data
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
  mutate(gene_to_genome_ratio=genes_num/size_mb)%>%
  arrange(host)

#Introduce App
#ui
ui<-dashboardPage(
  dashboardHeader(title="Viral Data"),
  dashboardSidebar(disable=T),
  dashboardBody(
    fluidRow(
  plotOutput("plot1", brush="plot_brush"),
  verbatimTextOutput("info"),
  #jitter input
  box(title="Plot Options", width=3,
      sliderInput("jitterw", "Jitter Width", min=0.5, max=40, value=12, step=5),
      sliderInput("jitterh", "Jitter height", min=0, max=.2, value=0, step=.01)
      ),#close jitter box
  box(title="Gene to Genome Ratio Overall",
      plotOutput("plot2"))
  )#close fluid row
  )#close dashboard body
)#close dashboard Page

#server
server<-function(input, output, session){
  #end clean
  session$onSessionEnded(stopApp)
  #output input
  output$plot1<-renderPlot({
    #plot parameters
    ggplot(virus_app, aes(x=genes_num, y=size_mb, size=perc_gc, color=host, group=host))+
      geom_jitter(alpha=0.25, width=input$jitterw, height=input$jitterh)+
      facet_grid(. ~kingdom)+
      #aesthetics
      theme_solarized()+
      scale_color_manual(values=colors)+
      theme(legend.position="left",
            axis.text.x=element_text(angle=60, hjust=1))+
      #labels
      labs(x="Number of Genes",
           y="Genome Size (Mb)",
           size="Percent GC",
           color="Host")
  })#end plot1
  output$plot2<-renderPlot({
    #plot2 parameters
    virus_app%>%
      group_by(host)%>%
    ggplot(aes(x=host, y=gene_to_genome_ratio, fill=kingdom))+
      geom_col()+
      coord_flip()+
      #aesthetics for plot 2
      theme_solarized()+
      scale_fill_manual(values=colors)+
      theme(legend.position="left",
            axis.text.x=element_text(angle=60, hjust=1))+
      #labels for plot 2
      labs(x="Host",
           y="Gene to Genome Ratio",
           color="Kingdom")
  })
  #output
  output$info<-renderPrint({
    brushedPoints(virus_app, input$plot_brush)
  })
}


#print
shinyApp(ui, server)


#To dos:
#1. Fix second graph whatever is going on there
#2. Add labels to recognizable viruses in first graph
#3. Look up aesthetic entries on shiny bc the app is ugly as hell
#4. Add citation
#5. Add flavor text?