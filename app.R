library(shiny)
library(shinythemes)
library(truncnorm)
library(ggplot2)

box.cox.trans <- function(x,lambda=1){
    if (lambda==0){
        return(log(x))
    }else{
        return((x^lambda - 1)/lambda)
    }
}
dens.box.cox.inv <- function(x,mu=5,sigma=1,lambda=1){
    return(dnorm(box.cox.trans(x,lambda=lambda),mean=mu,sd=sigma)*d.box.cox.trans(x,lambda=lambda))
}
d.box.cox.trans <- function(x,lambda=1){
    if (lambda==0){
        return(1/x)
    }else{
        return(x^(lambda - 1))
    }
}

box.cox.inv.trans <- function(x,lambda=1){
    if (lambda==0){
        return(exp(x))
    }else{
        return((x*lambda + 1)^(1/lambda))
    }
}

lod.boxplot.truncation <- function(measured.values,lod,n.lod,right.quantile=0.75){
    
    n.min.values <- 2*n.lod+1
    quantile.factor <- qnorm(0.975)/qnorm(right.quantile)
    
    x <- c(rep(lod-1,n.lod),measured.values)
    if (length(x)<n.min.values){
        return(list(selected.values=NA,upper.truncation=NA,median=NA,lod=lod,n.lod=n.lod))
    }
    
    x.qmr <- quantile(x,probs=c(0.5,right.quantile))
    upper.truncation <- x.qmr[1] + (x.qmr[2]-x.qmr[1])*quantile.factor
    
    x.length.old <- length(x)
    x <- x[x<=upper.truncation]
    if (length(x)<n.min.values){
        return(list(selected.values=NA,upper.truncation=NA,median=NA,lod=lod,n.lod=n.lod))
    }
    if (length(x)==x.length.old){
        return(list(selected.values=x[-(1:n.lod)],upper.truncation=upper.truncation,median=x.qmr[1],lod=lod,n.lod=n.lod))
    }
    
    median.mod <- 0.5/0.975
    right.quantile.mod <- right.quantile/0.975
    
    while(length(x)<x.length.old){
        x.length.old <- length(x)
        x.qmr <- quantile(x,probs=c(median.mod,right.quantile.mod))
        upper.truncation <- x.qmr[1] + (x.qmr[2]-x.qmr[1])*quantile.factor
        x <- x[x<=upper.truncation]
        if (length(x)<n.min.values){
            return(list(selected.values=NA,upper.truncation=NA,median=NA,lod=lod,n.lod=n.lod))
        }
    }  
    
    return(list(selected.values=x[-(1:n.lod)],upper.truncation=upper.truncation,median=x.qmr[1],lod=lod,n.lod=n.lod))
    
}

lodiboxplot95 <- function(measured.values,lod,n.lod,lambda=0,right.quantile=0.75){
    #cat(paste0("Current lambda value is: ", lambda, "\n"))
    transformed.measured.values <- box.cox.trans(measured.values,lambda=lambda)
    transformed.lod <- box.cox.trans(lod,lambda=lambda)

    normal.result <- lod.boxplot.truncation(transformed.measured.values,transformed.lod,n.lod,right.quantile=right.quantile)

    if (is.na(normal.result$upper.truncation)){
        return(list(lower.limit=NA,upper.limit=NA,mu.log=NA,sigma.log=NA,upper.truncation=NA,selected.values=NA,lod=lod,n.lod=n.lod))
    }

    obj.fun <- function(pars){
        pnorm.lod <- ptruncnorm(transformed.lod,b=normal.result$upper.truncation,mean=pars[1],sd=pars[2])

        return(-n.lod*log(pnorm.lod) - length(normal.result$selected.values)*log(1-pnorm.lod) - sum(log(dtruncnorm(normal.result$selected.values,a=transformed.lod,b=normal.result$upper.truncation,mean=pars[1],sd=pars[2]))))
    }

    pars.initial <- c(normal.result$median,(normal.result$upper.truncation-normal.result$median)/qnorm(0.975))
    optim.result <- optim(pars.initial,obj.fun)
    # --> Error in optim: function cannot be evaluated at initial parameters
    #optim.result <- nlminb(pars.initial,obj.fun)
    # --> Warning in nlminb(pars.initial, obj.fun) : NA/NaN function evaluation
    lims <- box.cox.inv.trans(qnorm(c(0.025,0.975),mean=optim.result$par[1],sd=optim.result$par[2]),lambda=lambda)
    lims1perc <- box.cox.inv.trans(qnorm(c(0.01,0.99),mean=optim.result$par[1],sd=optim.result$par[2]),lambda=lambda)
    
    #cat(paste0(" lower.limit is: ", lims[1], "\n"))
    
    return(list(lower.limit=lims[1],upper.limit=lims[2],percentile1=lims1perc[1],percentile99=lims1perc[2],mu.log=optim.result$par[1],sigma.log=optim.result$par[2],upper.truncation=box.cox.inv.trans(normal.result$upper.truncation,lambda=lambda),selected.values=box.cox.inv.trans(normal.result$selected.values,lambda=lambda),minus.log.likelihood=optim.result$value,lod=lod,n.lod=n.lod,lambda=lambda))
}

lod.hist <- function(res.lodiboxplot,xlab="",ylab="Frequency",main="",lwd=2,col.curve=2,breaks=10,col.grid=NULL){
    tmv <- res.lodiboxplot$selected.values
    lambda <- res.lodiboxplot$lambda
    mybreaks <- breaks
    if (is.null(breaks)){
        mybreaks <- 10
    }
    if (length(mybreaks)==1){
        mybreaks <- seq(from=min(tmv),to=max(tmv),length.out=mybreaks+1)
    }else{
        mybreaks <- seq(from=min(c(res.lodiboxplot$lod,min(breaks))),to=max(c(tmv,max(breaks))),length.out=length.breaks)
    }
    
    area <- (length(tmv)+res.lodiboxplot$n.lod)*(mybreaks[2]-mybreaks[1])
    hist(tmv,breaks=mybreaks,xlab=xlab,ylab=ylab,main=main)
    if (!is.null(col.grid)){
        grid(col=col.grid)
    }
    curve(area*dens.box.cox.inv(x,mu=res.lodiboxplot$mu.log,sigma=res.lodiboxplot$sigma.log,lambda=lambda),add=T,lwd=lwd,col=col.curve)
    
}
lod.qqplot <- function(res.lodiboxplot,pch=16,line.col="red",lwd=2,col.grid=NA,xlab="Theoretical Quantiles",ylab="Transformed Sample Quantiles"){
    
    truncated.measured.values <- res.lodiboxplot$selected.values
    n.lod <- res.lodiboxplot$n.lod
    lambda <- res.lodiboxplot$lambda
    
    transformed.measured.values <- box.cox.trans(truncated.measured.values,lambda=lambda)
    
    n.total <- (length(truncated.measured.values) + n.lod)/0.975
    
    quants.norm <- qnorm((n.lod+(1:length(truncated.measured.values)))/(n.total+1))
    
    quants.dat <- transformed.measured.values[order(transformed.measured.values)]
    plot(quants.norm,quants.dat,pch=pch,xlab=xlab,ylab=ylab)
    if (!is.na(col.grid)){
        grid(col=col.grid)
    }
    if (!is.na(line.col)){
        regline <- lm(quants.dat~quants.norm)
        abline(regline,col=line.col,lwd=lwd)
    }
}
compute.r.squared <- function(res.lodiboxplot){
    
    truncated.measured.values <- res.lodiboxplot$selected.values
    n.lod <- res.lodiboxplot$n.lod
    lambda <- res.lodiboxplot$lambda
    transformed.measured.values <- box.cox.trans(truncated.measured.values,lambda=lambda)
    
    n.total <- (length(truncated.measured.values) + n.lod)/0.975
    
    quants.norm <- qnorm((n.lod+(1:length(truncated.measured.values)))/(n.total+1))
    
    quants.dat <- transformed.measured.values[order(transformed.measured.values)]
    
    regline <- lm(quants.dat~quants.norm)
    
    sulm <-summary(regline)
    
    return(list(r.squared=sulm$r.squared,adj.r.squared=sulm$adj.r.squared))
}

plot.r.squared <- function(measured.values,lod,n.lod,lambdas=seq(from=0,to=1,by=0.1),adjusted=F,right.quantile=0.75,lwd=2,ylim=NULL,pch=16,col.grid=NULL){
    rsqs <- rep(-1,length(lambdas))
    strylab <- ""
    if (adjusted){
        strylab <- "Adjusted "	
    }
    for (i in 1:length(lambdas)){
        lbxpi <- lodiboxplot95(measured.values,lod,n.lod,lambda=lambdas[i],right.quantile=right.quantile)
        rsqs2 <- compute.r.squared(lbxpi)
        if (adjusted){
            rsqs[i] <- rsqs2$adj.r.squared
        }else{
            rsqs[i] <- rsqs2$r.squared
        }
    }
    
    plot(lambdas,rsqs,ylim=ylim,xlab=expression(lambda),ylab=paste0(strylab,"R-squared"),pch=pch)
    if (!is.null(col.grid)){
        grid(col=col.grid)
    }
    points(lambdas,rsqs,type="l",lwd=lwd)
    
    r.squared <- rsqs
    res <- data.frame(lambdas,r.squared) 
    colnames(res)[1] <- "lambda"
    
    return(res)
}
plot.reflims <- function(measured.values,lod,n.lod,lambdas=seq(from=0,to=1,by=0.1),right.quantile=0.75,lwd=2,ylim=NULL,ylab="",pch=c(15,16),col=c(2,4),col.grid=NULL){
    lower.limit <- rep(-1,length(lambdas))
    upper.limit <- rep(-1,length(lambdas))
    for (i in 1:length(lambdas)){
        lbxpi <- lodiboxplot95(measured.values,lod,n.lod,lambda=lambdas[i],right.quantile=right.quantile)
        lower.limit[i] <- lbxpi$lower.limit
        upper.limit[i] <- lbxpi$upper.limit
    }
    
    ylims <- ylim
    if (is.null(ylim)){
        ylims <- c(min(lower.limit),max(upper.limit))
    }
    plot(lambdas,lower.limit,ylim=ylims,xlab=expression(lambda),ylab=ylab,pch=pch[1],col=col[1])
    if (!is.null(col.grid)){
        grid(col=col.grid)
    }
    points(lambdas,lower.limit,type="l",lwd=lwd,col=col[1])
    points(lambdas,upper.limit,pch=pch[2],col=col[2])
    points(lambdas,upper.limit,type="l",lwd=lwd,col=col[2])
    
    res <- data.frame(lambdas,lower.limit,upper.limit) 
    colnames(res)[1] <- "lambda"
    
    return(res)
}

ci.lodiboxplot95 <- function(measured.values,lod,n.lod,lambda=0,right.quantile=0.75,conf.level=0.95,n.bootstrap=10){  #n.bootstrap =10 only for test, should be 1000
    #cat(paste0("Current lambda value is: ", lambda, "\n"))
    lower.limits <- rep(NA,n.bootstrap)
    upper.limits <- rep(NA,n.bootstrap)
    percentile1s <- rep(NA,n.bootstrap)
    percentile99s <- rep(NA,n.bootstrap)
    mu.logs <- rep(NA,n.bootstrap)
    sigma.logs <- rep(NA,n.bootstrap)
    upper.truncations <- rep(NA,n.bootstrap)
    
    x <- c(rep(lod/2,n.lod),measured.values)
    for (i in 1:n.bootstrap){
        xi <- sample(x,length(x),replace=T)
        mvi <- subset(xi,xi>=lod)
        n.lodi <- sum(xi<lod)
        
        resi <- lodiboxplot95(mvi,lod,n.lodi,lambda=lambda,right.quantile=right.quantile)
        
        if (!is.na(resi$lower.limit)){
            lower.limits[i] <- resi$lower.limit
            upper.limits[i] <- resi$upper.limit
            percentile1s[i] <- resi$percentile1
            percentile99s[i] <- resi$percentile99
            mu.logs[i] <- resi$mu.log
            sigma.logs[i] <- resi$sigma.log
            upper.truncations[i] <- resi$upper.truncation
        }
    }
    
    conf.lims <- c((1-conf.level)/2,1-(1-conf.level)/2)
    
    lower.limit.ci <- quantile(lower.limits,conf.lims,na.rm=T) 
    upper.limit.ci <- quantile(upper.limits,conf.lims,na.rm=T)      
    percentile1.ci <- quantile(percentile1s,conf.lims,na.rm=T)
    percentile99.ci <- quantile(percentile99s,conf.lims,na.rm=T)
    mu.log.ci <- quantile(mu.logs,conf.lims,na.rm=T)
    sigma.log.ci <- quantile(sigma.logs,conf.lims,na.rm=T)
    upper.truncation.ci <- quantile(upper.truncations,conf.lims,na.rm=T)
    
    #return(list(lower.limit.ci=lower.limit.ci,upper.limit.ci=upper.limit.ci,percentile1.ci=percentile1.ci,percentile99.ci=percentile99.ci,mu.log.ci=mu.log.ci,sigma.log.ci=sigma.log.ci,upper.truncation.ci=upper.truncation.ci))
    cis_df <- data.frame(
        # Lower_Limit_CI = round(lower.limit.ci, 6),
        # Upper_Limit_CI = round(upper.limit.ci, 6),
        # Percentile1_CI = round(percentile1.ci, 6),
        # Percentile99_CI = round(percentile99.ci, 6),
        # Mu_Log_CI = round(mu.log.ci, 6),
        # Sigma_Log_CI = round(sigma.log.ci, 6),
        # Upper_Truncation_CI = round(upper.truncation.ci, 6) not work
        Lower_Limit_CI = format(lower.limit.ci, nsmall = 6),
        Upper_Limit_CI = format(upper.limit.ci, nsmall = 6),
        Percentile1_CI = format(percentile1.ci, nsmall = 6),
        Percentile99_CI = format(percentile99.ci, nsmall = 6),
        Mu_Log_CI = format(mu.log.ci, nsmall = 6),
        Sigma_Log_CI = format(sigma.log.ci, nsmall = 6),
        Upper_Truncation_CI = format(upper.truncation.ci, nsmall = 6)
    )
    return(cis_df)
}
ui <- fluidPage(theme = shinytheme("yeti"),
                navbarPage(
                    "Student Data",
                    tabPanel("Navbar 1"),
                    sidebarLayout(
                        sidebarPanel(
                            fileInput("file", "Upload CSV file"),
                            numericInput("lod", "Enter LOD value:", value = 0, step = 0.1),
                            selectInput("filter_cols", "Select columns to filter",
                                        choices = NULL, multiple = TRUE),
                            uiOutput("filter_ui"),
                            selectInput("select_cols", "Select columns to display",
                                        choices = NULL, multiple = TRUE),
                            
                            numericInput("lambda", "Enter lambda value:", value = 0, step = 0.1),
                            actionButton("hist_btn", "Histogram"),
                            actionButton("qqplot_btn", "QQ Plot"),
                            actionButton("rsquared_btn", "R^2-Plot"),
                            actionButton("reflims_btn", "Plot reflims"),
                            actionButton("conf_level", "Confidence Level")
                        ), # sidebarPanel
                        
                        mainPanel(
                            h1("Output"),
                            h4("Text Output"),
                            verbatimTextOutput('lod'),
                            #verbatimTextOutput('num_na'),
                            verbatimTextOutput('num_string'),
                            verbatimTextOutput('lambda'),

                            verbatimTextOutput('num_x'),
 
                            h4("Output"),
                            textOutput("output_text"),
                            # Only show the selected plot when the button is clicked
                            conditionalPanel(
                                condition = "input.hist_btn > 0",
                                plotOutput("histogram")
                            ),
                            conditionalPanel(
                                condition = "input.qqplot_btn > 0",
                                plotOutput("qqplot")
                            ),
                            conditionalPanel(
                                condition = "input.rsquared_btn > 0",
                                plotOutput("rsquared")
                            ),
                            conditionalPanel(
                                condition = "input.reflims_btn > 0",
                                plotOutput("reflims")
                            ),
                       
                            conditionalPanel(
                                condition = "input.conf_level > 0",
                                #tags$div("11Error: The number of NAs exceeds 50% and the method cannot calculate the reference interval limits.", style = "color:red; font-weight:bold")
                                uiOutput("error_msg")
                                ),
                            conditionalPanel(
                                condition = "input.conf_level > 0",
                                tableOutput("confidence_table")
                            )
                            #tableOutput("table"),
                        )  # mainPanel
                    ),
                    tabPanel("Navbar 2", "This panel is intentionally left blank"),
                    tabPanel("Navbar 3", "This panel is intentionally left blank")
                ) #NavbarPage
)

server <- function(input, output, session) {
    num_string <- 0
    output$lod <- renderText({
        paste("The LOD value entered is:", input$lod)
    })
    
    data <- reactive({
        req(input$file)
        read.csv2(input$file$datapath, header = TRUE, na.strings = c("", "NA"))
    })
    selected_col <- reactiveVal()
    # Update the available columns when a file is uploaded
    observe({
        req(data())
        updateSelectInput(session, "filter_cols", 
                          choices = names(data()), selected = NULL)
        updateSelectInput(session, "select_cols", 
                          choices = names(data()), selected = NULL)
        selected_col(NULL)
    })
    
    # Display the selected columns
    output$ output_text <- renderText({
    #output$table <- renderTable({
        req(input$file)
        req(input$filter_cols)
        req(input$select_cols)
        req(data())
        
        # Set selected column name
        selected_col(input$select_cols[1])
        x <- data()
        for (filter_col in input$filter_cols) {
            if (filter_col %in% c("gender", "sex","geschlecht", "SEX")) {
                # If the selected column is gender or sex, use a dropdown menu
                x <- x[x[[filter_col]] %in% input[[filter_col]], ] #unique(x[[filter_col]]), ]
            } else {
                # Otherwise, use a numeric range
                x <- x[x[[filter_col]] >= input[[paste0(filter_col, "_min")]] & x[[filter_col]] <= input[[paste0(filter_col, "_max")]], ]
            }
        }
        x <- x[, input$select_cols, drop = FALSE]
        if (!is.null(x)) {
            num_na <- sum(apply(x, 1, function(row) any(is.na(row) & names(row) == selected_col())))
            num_string <- sum(apply(x, 1, function(row) any((is.character(row) & !is.na(row) & row == "<LOD") | (is.numeric(row) & row == 0)) & names(row) == selected_col()))
            x <- x[!is.na(x) & !(is.character(x) & x == "<LOD" & names(x) == selected_col() )& x != 0, ]
        }
        if (!is.null(x)) {
            x <- subset(x, !is.na(x) & !((is.character(x) & x == "<LOD" )| (is.numeric(x) & x == 0)))
        }
        output$num_string <- renderText({
            paste("The number of <LOD data is:", num_string)
        })
        output$lambda <- renderText({
            paste("The lambda value entered is:", input$lambda)
        })
        
        num_x <- x
        num_x <- length(num_x)
        output$num_x <- renderText({
            paste("The total number of results is:", num_x)
        })
 
        measured.values.m <- x
        lod <- input$lod
        n.lod.m <-num_string
        tnt.m <- lodiboxplot95(measured.values.m,lod,n.lod.m,lambda=input$lambda)

        output$histogram <- renderPlot({
            if (input$hist_btn > 0) {
                lod.hist(tnt.m, xlab = "TNT male", col.grid = 1)
            }
        })
        output$qqplot <- renderPlot({
            if (input$qqplot_btn > 0) {
                lod.qqplot(tnt.m, col.grid = 1)
            }
        })
        output$rsquared <- renderPlot({
            if (input$rsquared_btn > 0) {
                plot.r.squared(measured.values.m, lod, n.lod.m, col.grid = 1)
            }
        })
        output$reflims <- renderPlot({
            if (input$reflims_btn > 0) {
                plot.reflims(measured.values.m, lod, n.lod.m, col.grid = 1)
            }
        })
        output$error_msg <- renderUI({
            if (input$conf_level > 0) {
                if (num_string >= (num_x/2)) {
                    error_msg <- "Error: The number of NAs exceeds 50% and the method cannot calculate the reference interval limits."
                    return(tags$div(error_msg, style = "color:red; font-weight:bold"))
                }
            
           }
        })
        output$confidence_table <- renderTable({
            if (input$conf_level > 0) {
                if (num_string < (num_x/2)) {
                    cis <- ci.lodiboxplot95(measured.values.m,lod,n.lod.m,lambda=input$lambda)
                    return(cis)
                }
            }
        })
        x
    })
 
    # Create UI for filtering data based on selected columns and range of values
    output$filter_ui <- renderUI({
        req(input$file)
        req(input$filter_cols)
        req(data())
        
        filter_uis <- lapply(input$filter_cols, function(filter_col) {
            if (filter_col %in% c("gender", "sex","geschlecht", "SEX")) {
                choices <- unique(data()[[filter_col]])
                selectInput(filter_col, label = filter_col, choices = choices)

            } else {
                fluidRow(
                    column(width = 4,
                           numericInput(paste0(filter_col, "_min"), paste0(filter_col, " Minimum Value"), value = 0, step = 0.1)),
                    column(width = 4,
                           numericInput(paste0(filter_col, "_max"), paste0(filter_col, " Maximum Value"), value = 10, step = 0.1))
                )
            }
        })
        do.call(tagList, filter_uis)
    })
}

shinyApp(ui = ui, server = server)
