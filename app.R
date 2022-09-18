#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(gridExtra)
library(ggplot2)
library(deSolve)
library(scales)
library(dplyr)
library(EpiEstim)
library(tidyr)
library(reshape2)
library(cowplot)
library(knitr)
library(forcats)
library(bigsnpr)
library(EnvStats)
library(tidyverse)
library(ggpubr)
library(zoo)
library(scales) # comma formatting
library(rsconnect)
#rsconnect::deployApp('/Users/kj22643/Documents/Documents/Testing_UT/code/shinyapp/testing_strategy_tool_Omicron/code')

source("testing_cost_effectiveness_fxns.R")

PAR_TABLE_PATH<-'par_table.csv'
COST_TABLE_PATH<-'cost_table.csv'
par_table<- read.csv(PAR_TABLE_PATH)
cost_table<-read.csv(COST_TABLE_PATH)
test_freqs<- c(0, 1/30, 1/14, 1/7, 2/7, 3/7, 1)
test_policies<-c('only symptomatic', 'monthly', '2 times per month', 'weekly', '2 times per week', '3 times per week', 'daily')
case_detection_rate = 1/4.2


# Define UI for application (front-end)
ui <- pageWithSidebar(


    # Application title
    titlePanel("University Testing Strategy Tool"),

    # Sidebar with a slider input for number of bins 
    sidebarPanel(
            h3("Campus demographics"),
            sliderInput("N",
                        "Campus size:",
                        min = 0,
                        max = 100000,
                        value = 50000, step = 1000),
            sliderInput("vacc_rate",
                        "Percent optimally immunized:",
                        min = 0,
                        max = 100,
                        value = 60, step = 5, post = "%"), # set as 60 for UT Delta situation
                        # set as 20 for current Omicron situation without a booster mandate
            h6("As additional shots become recommended, this refers to the percent of students 
               with the most up-to-date vaccination status (i.e. boosters)"),
            sliderInput("init_rec",
                        "Percent immune from recent infection:",
                        min = 0,
                        max = 100,
                        value = 32, post ="%", step = 1),
            sliderInput("test_participation_rate",
                        "Percent willing to partipate in testing:",
                        min = 0,
                        max = 100, post = "%", step = 5,
                        value = 75), # Set  at 90 % if considering mandating testing
            sliderInput("pct_pos_isofac",
                        "Percent of students requiring isolation facilities:",
                        min = 0,
                        max = 100, post = "%", step = 5,
                        value = 20), # 20% for UT delta situation, likely to be higher in a smaller
                        # campus with more on-campus housing
            h6("Varies based on degree of on/off campus housing and anticipated student participation"),
            #sliderInput("duration_school_year",
                       # "Days in the semester:",
                       # min = 0,
                       # max = 360,
                       # value = 113)
            h3("Epidemiological variables"),
            # sliderInput("imported_infections",
            #             'Number of imported infections per week',
            #             min = 0,
            #             max = 60,
            #             value = 10, step = 5),
            # h6("Depends on the degree of interaction with the community and community disease prevalence"),
            sliderInput("init_cases",
                        'Initial campus infection prevalence (per 100k)',
                        min = 0,
                        max = 1000,
                        value = 400, step = 100),
            h6("Can be estimated from county-level cases in past 7 days multiplied by reporting rate"),
            
            sliderInput("R0",
                        'Basic reproductive number',
                        min = 2,
                        max = 10,
                        value = 5, step = 0.5), # 5 for Delta, 6.5 for Omicron?
            h6("A new variants increase in transmissibility, increase accordingly"),

            sliderInput("VE_against_inf",
                        "Efficacy of optimal immunization against infection:",
                        min = 0,
                        max = 100,
                        value = 50, step = 1, post = "%"), # for for Delta against infection
                        # for Omicron, with boosters, this might be 60%? But still lots of uncertainty
            # sliderInput("VE_against_symptoms",
            #             "Vaccine efficacy against symptomatic infections:",
            #             min = 0,
            #             max = 100,
            #             value = 64, step = 1, post = "%"),
            # h6("Note: vaccine efficacy against symptoms must be higher than efficacy against infection"),
            # sliderInput("transmissibility_if_vax_inf",
            #             "Relative transmissibility of individuals infected who are boosted:",
            #             min = 0,
            #             max = 100,
            #             value = 50, step = 5, post = "%"),
            h3("Costs"),
            sliderInput("cost_per_rapid_test",
                        "Cost per proactive test:",
                        min = 0,
                        max = 50,
                        value = 6, pre = "$", step = 1),
            sliderInput("cost_per_student_day_in_iso",
                        "Cost per day to isolate positive student:",
                        min = 0,
                        max = 500,
                        value = 300, pre = "$", step = 50),
            sliderInput("cost_per_days_online",
                        "Cost per days of university moving online:",
                        min = 0,
                        max = 200000,
                        value = 100000, pre = "$", step = 1000),
            h3("Test policies and tolerance for disease"),
            selectInput("test", "Population tested:",
                        c("All students" = "u&v",
                          "Half the rate in vaccinated" = "u&0.5v",
                          "Unvaccinated only" = "u")),
            selectInput("risk_tolerance", "Tolerance for symptomatic cases:",
                        c("CDC high (100 cases per 100k in 7 days)" = "CDC red",
                          "1.5x CDC high (150 cases per 100k in 7 days)" = "1.5x CDC red",
                          "2x CDC high (200 cases per 100k in 7 days)" = "2x CDC red",
                          "3x CDC high (300 cases per 100k in 7 days)" = "3x CDC red",
                          "4x CDC high (200 cases per 100k in 7 days)" = "4x CDC red"),
                        selected = "2x CDC red")
        ),

        # Show a plot of the cases over time, the weekly tests, and the costs associated with 
        #each testing frequency
        mainPanel(
           htmlOutput("text"),
           plotOutput("cases_t")
            )
)

#par_table<-read.csv(PAR_TABLE_PATH)
#input = par_table

# Data pre-processing ---



# Define server logic (relationship between inputs and output)
server <- function(input, output) {
    print("here")
  
   out_list= reactive({
       par_table$N = input$N
       par_table$R0 = input$R0
       #par_table$duration_school_year = input$duration_school_year
       par_table$init_vacc = input$vacc_rate/100
       par_table$sigma_v = 1-input$VE_against_inf/100
       #par_table$e_v = input$transmissibility_if_vax_inf/100
       #par_table$sym_red = 1-input$VE_against_symptoms/100
       par_table$init_prev = input$init_cases/100000
       par_table$init_rec = input$init_rec/100
       #par_table$intros_per_week = input$imported_infections
       par_table$w_u = input$test_participation_rate/100
       cost_table$RT_price = input$cost_per_rapid_test
       cost_table$online_price = input$cost_per_days_online
       cost_table$pct_pos_isofac = input$pct_pos_isofac/100
       cost_table$isofac_price = input$cost_per_student_day_in_iso
       
       #par_table$e_v = input$VE_transmission
       #par_table$is_sym = input$sym_det_rate

       
       #risk_tolerance = input$risk_tolerance
       
       close_thres = case_when(
         input$risk_tolerance == "CDC red" ~ 100/100000,
         input$risk_tolerance == "1.5x CDC red" ~ 150/100000,
         input$risk_tolerance == "2x CDC red" ~200/100000,
         input$risk_tolerance == "3x CDC red" ~300/100000,
         input$risk_tolerance == "4x CDC red" ~400/100000)
       risk_tolerance = input$risk_tolerance
       
       population_tested = case_when(
         input$test == "u" ~ 'unvaccinated students only',
         input$test == "u&0.5v" ~'unvaccinated students and half the rate in those vaccinated',
         input$test == "u&v" ~ 'all students'
       )
 
       
       for(k in 1:length(test_freqs)){
         print("in loop")
         testing_freq<-test_freqs[k]
         testing_policy<-test_policies[k]
         par_table$f_u<-test_freqs[k]
         freq_vacc = case_when(
           input$test == "u&v" ~ test_freqs[k],
           input$test == "u&0.5v" ~ 0.5*test_freqs[k],
           input$test == "u" ~0)
         par_table$f_v<-freq_vacc
         test_freq<-rep(test_freqs[k], par_table$duration_school_year+1)
         test_policy<-rep(test_policies[k], par_table$duration_school_year+1)
         out_df<-solve_transmission_eqns(par_table)
         aggr_df_t<-cbind(test_policy, test_freq, out_df)
         
         sum_df<-testing_model(par_table, cost_table, risk_tolerance)
         aggr_df<-cbind(testing_policy, population_tested, sum_df)
         if (k==1){
           df_t<-aggr_df_t
           df<-aggr_df
         }
         else{
           df_t<-rbind(df_t, aggr_df_t)
           df<-rbind(df, aggr_df)
         }
       }
       
       df_t$test_policy<-factor(df_t$test_policy, ordered = TRUE,
                                     stringr::str_wrap(c("only symptomatic", "monthly", "2 times per month", "weekly", "2 times per week", "3 times per week", "daily")))
       
       df$testing_policy<-factor(df$testing_policy, ordered = TRUE,
                                stringr::str_wrap(c("only symptomatic", "monthly", "2 times per month", "weekly", "2 times per week", "3 times per week", "daily")))

       print("here2")
       out_list<-list(df_t, df)

       #out_df
       #out_1<-data.frame(dates = 1:10, Isymdet = rep(b,10))
   
       
                          
       
    })
   
   output$cases_t = renderPlot( height = 800, width = 600, res = 150, {
      out = out_list()
      df_t = out[[1]]
      df = out[[2]]
      N <- df_t$Nvec[1]
      close_thres <- df$close_thres[1]
       gg.cases = ggplot(data = df_t) + 
         geom_hline(yintercept = close_thres*N, color = "red", linetype = "dashed")+
           geom_line(aes(x = time, y = Isymdet, color = test_policy))+
         scale_color_brewer(palette = "Set2") +
         #coord_cartesian(ylim = c(0, 20*N*100/100000))+
          # scale_fill_brewer(palette = "Set2", direction = -1)+
         labs(color = "Test frequency") + 
         theme(text = element_text(size=8),
               axis.text.x = element_text(size = 8),
               axis.text.y = element_text(size = 8))+
           xlab('Time (days)')+ ylab('Symptomatic Detected Cases')+
         ggtitle('Cases in past 7 days at different testing levels')
       
       gg.inf = ggplot(df, aes(x = testing_policy, y = n_inf/1000)) +
         geom_bar(stat = "identity", position = "stack")  +
         xlab('Test frequency') + ylab('Infections (thousands)') +
        theme(axis.text.x = element_blank(), 
              axis.title.x = element_blank(),
              axis.text.y = element_text(size = 8),
              axis.title = element_text(size = 8),
              legend.position = "none")
       
     df_costs_long<- df%>%select(testing_policy, 
                                 pct_inf, cost_PCR,
                                 cost_RT, cost_isofac, cost_sequencing,
                                 cost_contact_tracing, cost_of_online)%>%
       pivot_longer(cols = starts_with("cost_"),names_to = "Source",names_prefix = "cost_", values_to = "cost" )
     df_costs_long$Source<-factor(df_costs_long$Source, ordered = TRUE,
                            stringr::str_wrap(c("of_online", "contact_tracing",
                                                "isofac", 
                                               "PCR", "sequencing", 
                                               "RT")))
     
 
     gg.costs = ggplot(df_costs_long, aes (x = testing_policy, y = cost/N, fill = factor(Source))) +
       geom_bar(stat = "identity", position = "stack") + 
       xlab('Test Frequency') + 
       ylab('Cost to university per student') +
      labs(fill = "Source of cost")+
       theme(text = element_text(size=8),
             axis.text.x = element_text(angle = 30, hjust =1, size = 6),
             legend.text = element_text(size = 6),
             axis.text.y = element_text(size = 8),
             legend.position = "right")+
       #coord_cartesian(ylim = c(0, 700))+
       scale_fill_manual(" Source of Cost",values = c('gray', 'pink', 'coral2',  'orange', 'yellow', 'blue4'),
                         breaks = c('of_online', 'contact_tracing', 'isofac', 
                                    'PCR', "sequencing", "RT"),
                         labels =c('online', 'contact tracing', 'isolation', 
                                   'confirmatory PCR', 'sequencing', 'proactive (rapid) tests')) + 
       ggtitle('Cost breakdown at each testing level')
     gg.costs
     
 
     grid.arrange(gg.cases, gg.costs, nrow=2)
   })
     
   
   
   # Find minimum testing to keep cases below close_thres
output$text = renderUI({
      out = out_list()
      df_t = out[[1]]
      df = out[[2]]
      df_under<-df[df$days_of_online<15,]
      df_min<-df_under[which.min(df_under$testing_freq),]
      str1<-(paste("With the parameters currently selected,", df_min$testing_policy, "testing in", df_min$population_tested,
                   "is needed to stay under the tolerance threshold."))
      str2<-paste("This corresponds to", comma_format()(df_min$n_tests_per_week)," proactive tests per week.")
      str3<-paste0("This would cost about $", round(df_min$cost_testing_per_student, 0), " per student")
                  #round(df_min$cost_per_DO_averted,0), "per days of online averted.")
      HTML(paste(str1, str2,str3, sep = '<br/>'))
})
     
     
  
     
 
    
}

# Run the application 
print("running shiny")
res=shinyApp(ui = ui, server = server)
print("after")

res