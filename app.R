library(shiny)
library(tidyverse)
library(janitor)
library(broom)
library(gt)
library(DT)

# simulation function for power
sim.dataset<-function(n,expb0,expb1,p,t){
    b0<-log(expb0)
    b1<-log(expb1)
    simdata<-data.frame(x=rbinom(n=n,1,p=p))
    simdata$lam=exp(b0+b1*simdata$x)
    simdata$y<-rpois(n=n, lambda=t*simdata$lam)
    return(simdata)}

s.power<-function(sim,tail,nsim,n1,n2,st,expb0,expb1,p,alpha,t){
    set.seed(1)
    N=(n2-n1)/st+1
    b0<-log(expb0)
    b1<-log(expb1)
    coef.power<-data_frame(sig.b0=rep(NA,nsim),
                           sig.b1=rep(NA,nsim))
    power.sam<-data_frame(n=rep(NA,N),
                          power_simu=rep(NA,N))
    
    for (j in 1:N){
        for (dataset_i in 1:nsim){
            simulated_df<-sim.dataset(n1+st*(j-1),expb0,expb1,p,t)
            model<-tidy(glm(y~x,data = simulated_df,family="poisson"))
            coef.power$sig.b0[dataset_i]<-model$p.value[1]<(alpha*(tail==1)+alpha)
            coef.power$sig.b1[dataset_i]<-model$p.value[2]<(alpha*(tail==1)+alpha)
        }
        power.sam$n[j]<-n1+st*(j-1)
        power.sam$power_simu[j]<-summarise(coef.power, 
                                           b0.power=mean(sig.b0),
                                           b1.power=mean(sig.b1))$b1.power}
    if(sim=="Yes"){
        return(power.sam)}else{
            return(data_frame(n=seq(n1,n2,st)))
        }
}

# G-power function for power based on D2007 with variance correction
g.power<-function(tail,n1,n2,st,expb0,expb1,p,alpha,t){  
    b0<-log(expb0)
    b1<-log(expb1)
    N=(n2-n1)/st+1
    power.g<-data_frame(n=rep(NA,N),
                        power_g=rep(NA,N))
    for (j in 1:N){ 
        a<-qnorm(1-alpha/tail)
        b<-b1^2*t
        expb00<-p*exp(b0)*exp(b1)+(1-p)*exp(b0)
        varb1<-(1/(1-p)+1/p/exp(b1))/exp(b0)
        varb00<-1/p/(1-p)/expb00
        zp<-(sqrt((n1+st*(j-1))*b)-a*sqrt(varb00))/sqrt(varb1)
        power.g$power_g[j]<-pnorm(zp)
        power.g$n[j]<-n1+st*(j-1)
    }
    return(power.g)}

# PASS function for power based on S1991
p.power<-function(tail,n1,n2,st,expb0,expb1,p,alpha,t){
    b0<-log(expb0)
    b1<-log(expb1)
    N=(n2-n1)/st+1
    power.p<-data_frame(n=rep(NA,N),
                        power_p=rep(NA,N))
    for (j in 1:N){ 
        a<-exp(b0)*b1*b1*t
        b<-qnorm(1-alpha/tail)*sqrt(1/p/(1-p))
        c<-sqrt(1/(1-p)+1/p/exp(b1))
        zp<-(sqrt((n1+st*(j-1))*a)-b)/c
        power.p$power_p[j]<-pnorm(zp)
        power.p$n[j]<-n1+st*(j-1)
    }
    return(power.p)}

# R package function for power based on D2007 without variance correction
w.power<-function(tail,n1,n2,st,expb0,expb1,p,alpha,t){
    N=(n2-n1)/st+1
    power.w<-data_frame(n=rep(NA,N),
                        power_w=rep(NA,N))
    for (j in 1:N){ 
        a<-qnorm(1-alpha/tail)
        b<-(log(expb1))^2*t
        varb1<-(1/(1-p)+1/p/expb1)/expb0
        zp<-sqrt((n1+st*(j-1))*b/varb1)-a
        power.w$power_w[j]<-pnorm(zp)
        power.w$n[j]<-n1+st*(j-1)
    }
    return(power.w)}


# G-power function for sample size based on D2007 with variance correction
g.size<-function(tail,expb0,expb1,p,alpha,power1,power2,pt,t){  
    b0<-log(expb0)
    b1<-log(expb1)
    M=abs(round((power2-power1)/pt)+1)
    size.g<-data_frame(power=rep(NA,M),
                       size_g=rep(NA,M))
    for (j in 1:M){ 
        a<-qnorm(1-alpha/tail)
        b<-b1^2*t
        c<-qnorm(power1+pt*(j-1))
        expb00<-p*expb0*expb1+(1-p)*expb0
        varb1<-(1/(1-p)+1/p/expb1)/expb0
        varb00<-1/p/(1-p)/expb00
        s<-round((a*sqrt(varb00)+c*sqrt(varb1))^2/b,digits = 0)
        size.g$size_g[j]<-s
        size.g$power[j]<-power1+pt*(j-1)}
    return(size.g)
}

# PASS function for sample size based on S1991
p.size<-function(tail,expb0,expb1,p,alpha,power1,power2,pt,t){
    b0<-log(expb0)
    b1<-log(expb1)
    M=abs(round((power2-power1)/pt)+1)
    size.p<-data_frame(power=rep(NA,M),
                       size_p=rep(NA,M))
    for (j in 1:M){ 
        a<-qnorm(1-alpha/tail)*sqrt(1/p/(1-p))
        b<-qnorm(power1+pt*(j-1))*sqrt(1/(1-p)+1/p/expb1)
        c<-expb0*log(expb1)*log(expb1)*t
        n<-round((a+b)^2/c,digits = 0)
        size.p$size_p[j]<-n
        size.p$power[j]<-power1+pt*(j-1)}
    return(size.p)}

# R package function for sample size based on D2007 without variance correction
w.size<-function(tail,expb0,expb1,p,alpha,power1,power2,pt,t){
    b0<-log(expb0)
    b1<-log(expb1)
    M=abs(round((power2-power1)/pt)+1)
    size.w<-data_frame(power=rep(NA,M),
                       size_w=rep(NA,M))
    for (j in 1:M){ 
        a<-qnorm(1-alpha/tail)+qnorm(power1+pt*(j-1))
        b<-(log(expb1))^2*t
        varb1<-(1/(1-p)+1/p/expb1)/expb0
        n<-round(a^2*varb1/b,digits = 0)
        size.w$size_w[j]<-n
        size.w$power[j]<-power1+pt*(j-1)
    }
    return(size.w)}

# G-power function for effect size based on D2007 with variance correction
g.effect<-function(tail,expb0,p,alpha,n1,n2,st,power1,power2,pt,t){
    
    fn<-function(tail,expb0,b1,p,alpha,n,power,t){
        a<-qnorm(1-alpha/tail)
        b<-qnorm(power)
        expb00<-p*expb0*exp(b1)+(1-p)*expb0
        varb1<-(1/(1-p)+1/p/exp(b1))/expb0
        varb00<-1/p/(1-p)/expb00
        f<-n-(a*sqrt(varb00)+b*sqrt(varb1))^2/b1/b1/t
        return(f)
    }
    
    N=(n2-n1)/st+1
    M=abs(round((power2-power1)/pt)+1)
    
    effect.g<-data_frame(power=rep(NA,M*N),
                         sample_size=rep(NA,M*N),
                         effect_g=rep(NA,M*N))
    
    for (i in 1:N){
        for (j in 1:M){  
            b1s<-uniroot(fn,
                         lower = -1000, 
                         upper = 1000,
                         extendInt = "yes",
                         tail=tail,
                         expb0=expb0,
                         p=p,
                         alpha=alpha,
                         n=n1+st*(i-1),
                         power=power1+pt*(j-1),
                         t=1,
                         tol=0.001)$root
            k=j+M*(i-1)
            effect.g$effect_g[k]<-exp(b1s)
            effect.g$power[k]<-power1+pt*(j-1)
            effect.g$sample_size[k]<-n1+st*(i-1)
        }
    }
    return(effect.g)
}

# PASS function for effect size based on S1991
p.effect<-function(tail,expb0,p,alpha,n1,n2,st,power1,power2,pt,t){
    
    fn<-function(tail,expb0,b1,p,alpha,n,power,t){
        a<-qnorm(1-alpha/tail)
        b<-qnorm(power)
        c<-sqrt(1/p/(1-p))
        d<-sqrt(1/(1-p)+1/p/exp(b1))
        f<-n-(a*c+b*d)^2/b1/b1/expb0/t
        return(f)
    }
    
    N=(n2-n1)/st+1
    M=abs(round((power2-power1)/pt)+1)
    
    effect.p<-data_frame(power=rep(NA,M*N),
                         sample_size=rep(NA,M*N),
                         effect_p=rep(NA,M*N))
    
    for (i in 1:N){
        for (j in 1:M){  
            b1s<-uniroot(fn,
                         lower = -1000, 
                         upper = 1000,
                         extendInt = "yes",
                         tail=tail,
                         expb0=expb0,
                         p=p,
                         alpha=alpha,
                         n=n1+st*(i-1),
                         power=power1+pt*(j-1),
                         t=1,
                         tol=0.001)$root
            k=j+M*(i-1)
            effect.p$effect_p[k]<-exp(b1s)
            effect.p$power[k]<-power1+pt*(j-1)
            effect.p$sample_size[k]<-n1+st*(i-1)
        }
    }
    return(effect.p)
}


# R package function for effect size based on D2007 without variance correction
w.effect<-function(tail,expb0,p,alpha,n1,n2,st,power1,power2,pt,t){
    
    fn<-function(tail,expb0,b1,p,alpha,n,power,t){
        a<-qnorm(1-alpha/tail)
        b<-qnorm(power)
        c<-(1/(1-p)+1/p/exp(b1))/expb0
        f<-n-(a+b)^2*c/b1/b1/t
        return(f)
    }
    
    N=(n2-n1)/st+1
    M=abs(round((power2-power1)/pt)+1)
    
    effect.w<-data_frame(power=rep(NA,M*N),
                         sample_size=rep(NA,M*N),
                         effect_w=rep(NA,M*N))
    
    for (i in 1:N){
        for (j in 1:M){  
            b1s<-uniroot(fn,
                         lower = -1000, 
                         upper = 1000,
                         extendInt = "yes",
                         tail=tail,
                         expb0=expb0,
                         p=p,
                         alpha=alpha,
                         n=n1+st*(i-1),
                         power=power1+pt*(j-1),
                         t=1,
                         tol=0.001)$root
            k=j+M*(i-1)
            effect.w$effect_w[k]<-exp(b1s)
            effect.w$power[k]<-power1+pt*(j-1)
            effect.w$sample_size[k]<-n1+st*(i-1)
        }
    }
    return(effect.w)
}

# R shiny
parameter_tabs <- tabsetPanel(
    id = "params",
    type = "hidden",
    tabPanel("Power",
             numericInput(inputId="expb1",
                          label="expb1",
                          value=1.3,
                          min = 0,
                          max = 1000),
             
             numericInput(inputId="n1",
                          label="Lower bound of sample size",
                          value=5,
                          min = 0,
                          max = 50000),
             
             numericInput(inputId="n2",
                          label="Upper bound of sample size",
                          value=50,
                          min = 0,
                          max = 50000),
             
             numericInput(inputId="st",
                          label="Increment of sample size",
                          value=5,
                          min = 0,
                          max = 10000),
             
             selectInput(inputId = "sim",
                         label = "Want Simulation? Be careful of time!",
                         choices = c("No","Yes")),
             numericInput(inputId="nsim",
                          label="Simulation times",
                          value=10,
                          min = 0,
                          max = 10000),
             actionButton(inputId ="do1",
                          label = "Click me" ), 
             
    ),
    tabPanel("SampleSize", 
             numericInput(inputId="expb1",
                          label="expb1",
                          value=1.3,
                          min = 0,
                          max = 1000),
             
             numericInput(inputId="power1",
                          label="Lower bound of power",
                          value=0.8,
                          min = 0,
                          max = 1),
             
             numericInput(inputId="power2",
                          label="Upper Bound of power",
                          value=0.9,
                          min = 0,
                          max = 1),
             
             numericInput(inputId="pt",
                          label="Increment of power",
                          value=0.01,
                          min = 0,
                          max = 1),
             
             actionButton(inputId ="do2",
                          label = "Click me" ), 
    ),
    
    tabPanel("EffectSize",
             numericInput(inputId="n11",
                          label="Lower bound of sample size",
                          value=5,
                          min = 0,
                          max = 50000),
             
             numericInput(inputId="n21",
                          label="Upper bound of sample size",
                          value=50,
                          min = 0,
                          max = 50000),
             
             numericInput(inputId="st1",
                          label="Increment of sample size",
                          value=5,
                          min = 0,
                          max = 10000),
             
             numericInput(inputId="power11",
                          label="Lower bound of power",
                          value=0.8,
                          min = 0,
                          max = 1),
             
             numericInput(inputId="power21",
                          label="Upper Bound of power",
                          value=0.9,
                          min = 0,
                          max = 1),
             
             numericInput(inputId="pt1",
                          label="Increment of power",
                          value=0.01,
                          min = 0,
                          max = 1),
             
             actionButton(inputId ="do3",
                          label = "Click me" )
             
    )
)
ui <- fluidPage(
    titlePanel("Power Analysis of Poisson Regression with Bernoulli Distribution "),
        tabsetPanel(
            tabPanel("Power anaylysis",
                     sidebarLayout(
                         sidebarPanel(width = 4,
                             selectInput("type", "Type of Power analysis", 
                        choices = c("Power", "SampleSize", "EffectSize")),
                        
                        numericInput(inputId="tail",
                                     label="Alternative (two sided=2 or one side=1)",
                                     value=1,
                                     min = 1,
                                     max = 2),   
                        
                        numericInput(inputId="expb0",
                                     label="expb0",
                                     value=0.85,
                                     min = 0,
                                     max = 1000),
                        
                        numericInput(inputId="alpha",
                                     label="Type I error",
                                     value=0.05,
                                     min = 0,
                                     max = 1),
                        
                        numericInput(inputId="t",
                                     label="Mean exposure",
                                     value=1,
                                     min = 0,
                                     max = 1000),
                        
                        numericInput(inputId="p",
                                     label="Parameter of bernoulli distribution",
                                     value=0.5,
                                     min = 0,
                                     max = 1),
                        parameter_tabs,
                        
                        downloadButton(outputId = "downloadData", 
                                       label = "Download")
                        ),
                     mainPanel(
                         plotOutput("plotp"),
                         DT::dataTableOutput("tablep")
                         )
                     )
                     ),
    tabPanel("Document",
             mainPanel(
                 htmlOutput("inc")
             )
    )
        )
)
server <- function(input, output, session) {
    observeEvent(input$type, {
        updateTabsetPanel(inputId = "params", selected = input$type)
    }) 
    
    table_s<-eventReactive(input$do1,
                           {s.power(input$sim,
                                    input$tail,
                                    input$nsim,
                                    input$n1,
                                    input$n2,
                                    input$st,
                                    input$expb0,
                                    input$expb1,
                                    input$p,
                                    input$alpha,
                                    input$t)}) 
    table_g<-eventReactive(input$do1, 
                           {g.power(input$tail,
                                    input$n1,
                                    input$n2,
                                    input$st,
                                    input$expb0,
                                    input$expb1,
                                    input$p,
                                    input$alpha,
                                    input$t)}) 
    
    table_w <- eventReactive(input$do1, 
                             {w.power(input$tail,
                                      input$n1,
                                      input$n2,
                                      input$st,
                                      input$expb0,
                                      input$expb1,
                                      input$p,
                                      input$alpha,
                                      input$t)})  
    table_p <- eventReactive(input$do1, 
                             {p.power(input$tail,
                                      input$n1,
                                      input$n2,
                                      input$st,
                                      input$expb0,
                                      input$expb1,
                                      input$p,
                                      input$alpha,
                                      input$t)})
    
    table_sg<-eventReactive(input$do2, 
                            {g.size(input$tail,
                                    input$expb0,
                                    input$expb1,
                                    input$p,
                                    input$alpha,
                                    input$power1,
                                    input$power2,
                                    input$pt,
                                    input$t)}) 
    table_sp<-eventReactive(input$do2, 
                            {p.size(input$tail,
                                    input$expb0,
                                    input$expb1,
                                    input$p,
                                    input$alpha,
                                    input$power1,
                                    input$power2,
                                    input$pt,
                                    input$t)}) 
    table_sw<-eventReactive(input$do2, 
                            {w.size(input$tail,
                                    input$expb0,
                                    input$expb1,
                                    input$p,
                                    input$alpha,
                                    input$power1,
                                    input$power2,
                                    input$pt,
                                    input$t)}) 
    table_fg<-eventReactive(input$do3, 
                            {g.effect(input$tail,
                                    input$expb0,
                                    input$p,
                                    input$alpha,
                                    input$n11,
                                    input$n21,
                                    input$st1,
                                    input$power11,
                                    input$power21,
                                    input$pt1,
                                    input$t)}) 
    table_fp<-eventReactive(input$do3, 
                            {p.effect(input$tail,
                                      input$expb0,
                                      input$p,
                                      input$alpha,
                                      input$n11,
                                      input$n21,
                                      input$st1,
                                      input$power11,
                                      input$power21,
                                      input$pt1,
                                      input$t)}) 
    table_fw<-eventReactive(input$do3, 
                            {w.effect(input$tail,
                                      input$expb0,
                                      input$p,
                                      input$alpha,
                                      input$n11,
                                      input$n21,
                                      input$st1,
                                      input$power11,
                                      input$power21,
                                      input$pt1,
                                      input$t)}) 
    
    
    table <- reactive({
        switch(input$type,
               Power = purrr::reduce(list(table_s(),table_g(),table_p(),table_w()),
                                     by="n", 
                                     dplyr::left_join),
               SampleSize = purrr::reduce(list(table_sg(),table_sp(),table_sw()), 
                                         by="power",
                                         dplyr::left_join),
               EffectSize = purrr::reduce(list(table_fg(),table_fp(),table_fw()),
                                         dplyr::left_join)
        )
    })
    
    output$tablep <- DT::renderDataTable(round(table(),digits  = 3),
                                         options = list(scrollX = TRUE),
                                         rownames = FALSE)
    
    
    output$plotp<-renderPlot({
        if (ncol(table())==5 & "power_simu" %in% colnames(table())){
            plotp<-ggplot(data=table(),aes(x=n))+
                geom_line(aes(y=power_simu,color="Simulation"),size=0.5)+
                geom_line(aes(y=power_g,color="G-power"),size=0.5)+
                geom_line(aes(y=power_p,color="PASS"),size=0.5)+  
                geom_line(aes(y=power_w,color="Webpower"),size=0.5)+ 
                scale_colour_manual("", 
                                    breaks = c("Simulation", "G-power", "PASS", "Webpower"),
                                    values = c("red", "green", "blue","black"))+
                ggtitle("Power with different methods by Sample Size")+
                labs(x="Samper Size", y="Power")+theme_bw()
            }else if (colnames(table())[1]=="n"){
                    plotp<-ggplot(data=table(),aes(x=n))+
                        geom_line(aes(y=power_g,color="G-power"),size=0.5)+
                        geom_line(aes(y=power_p,color="PASS"),size=0.5)+  
                        geom_line(aes(y=power_w,color="Webpower"),size=0.5)+ 
                        scale_colour_manual("", 
                                            breaks = c("G-power", "PASS", "Webpower"),
                                            values = c("green", "blue","black"))+
                        ggtitle("Power with different methods by Sample Size")+
                        labs(x="Samper Size", y="Power")+theme_bw()
                    }else if ("size_g" %in% colnames(table()))
                    {plotp<-ggplot(data=table(),aes(x=power))+
                    geom_line(aes(y=size_g,color="G-power"),size=0.5)+
                    geom_line(aes(y=size_p,color="PASS"),size=0.5)+  
                    geom_line(aes(y=size_w,color="Webpower"),size=0.5)+ 
                    scale_colour_manual("", 
                                        breaks = c("G-power", "PASS", "Webpower"),
                                        values = c("green", "blue","black"))+
                    ggtitle("Sample Size with different methods by Power")+
                    labs(x="Power", y="Sample Size") +theme_bw()
                    }else{
                    plotp<-ggplot(data=table(),aes(x=power,group=sample_size))+
                        geom_line(aes(y=effect_g,color="G-power"),size=0.5)+
                        geom_line(aes(y=effect_p,color="PASS"),size=0.5)+  
                        geom_line(aes(y=effect_w,color="Webpower"),size=0.5)+ 
                        scale_colour_manual("", 
                                            breaks = c("G-power", "PASS", "Webpower"),
                                            values = c("green", "blue","black"))+
                        ggtitle("Effect Size with different methods by Power and Sample Size")+
                        labs(x="Power", y="Effect Size") +theme_bw()
                        
                    }
        plotp
    })
    
    getPage<-function() {
        return(includeHTML("document.html"))
    }
    output$inc<-renderUI({getPage()})
    
    output$downloadData <- downloadHandler(
        filename = function() {
            paste(input$type, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(table(), file, row.names = FALSE)
        }
    )
    
}

shinyApp(ui, server)
