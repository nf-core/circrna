#!/usr/bin/env Rscript

require(aod)
require(plyr)
require(ggplot2)

# The CircTest devs never had containers in mind for their code. It is not published on any R channels, nor are there any stable releases on their github page. I have decided to take their source functions and roll with it instead.

## SUMMMARY

#'@title  Summarizes data
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
#'
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {

    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                    .fun = function(xx, col) {
                        c(N    = length2(xx[[col]], na.rm=na.rm),
                        mean = mean   (xx[[col]], na.rm=na.rm),
                        sd   = sd     (xx[[col]], na.rm=na.rm)
                    )
                    },
                    measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

## RATIO PLOT

Circ.ratioplot <- function(Circ,Linear,CircCoordinates = None,plotrow='1',size=24,ncol=2,groupindicator1=NULL,groupindicator2=NULL,x='Conditions',y='circRNA/(circRNA+Linear)',lab_legend='groupindicator1', circle_description = c(1:3), gene_column = None, y_axis_range = 1, colour_mode = "colour"){

    if( !is.null(groupindicator1) & length(groupindicator1) != ncol(Circ)-length(circle_description) ){
        stop("If provided, the length of groupindicator1 should be equal to the number of samples.")
    }
    if( !is.null(groupindicator2) & length(groupindicator2) != ncol(Circ)-length(circle_description) ){
        stop("If provided, the length of groupindicator2 should be equal to the number of samples.")
    }
    if(is.null(groupindicator1)){
        stop("At least one grouping should be provided through groupindicator1.")
    }
    if(!is.null(groupindicator2)){
        twolevel <- TRUE
    }else{
        twolevel <- FALSE
    }

    rownames.circ <- rownames(Circ)
    Circ <- data.frame(lapply(Circ, as.character), stringsAsFactors=FALSE)
    rownames(Circ) <- rownames.circ

    rownames.linear <- rownames(Linear)
    Linear <- data.frame(lapply(Linear, as.character), stringsAsFactors=FALSE)
    rownames(Linear) <- rownames.linear

    if(!missing(CircCoordinates)){
        rownames.circCoordinates <- rownames(CircCoordinates)
        CircCoordinates <- data.frame(lapply(CircCoordinates, as.character), stringsAsFactors=FALSE)
        rownames(CircCoordinates) <- rownames.circCoordinates
    }else{
        CircCoordinates <- data.frame(Circ[,circle_description])
        rownames(CircCoordinates) <- rownames.circ
        rownames.circCoordinates <- rownames(CircCoordinates)
        CircCoordinates <- data.frame(lapply(CircCoordinates, as.character), stringsAsFactors=FALSE)
        rownames(CircCoordinates) <- rownames.circCoordinates
    }

    groupindicator1 <- factor(groupindicator1,levels=unique(groupindicator1))
    groupindicator2 <- factor(groupindicator2,levels=unique(groupindicator2))

    # Get gene name, if no annotation, output NULL
    if (is.character(plotrow)){
        if ( ! plotrow %in% rownames(CircCoordinates) ){
            stop("Specified 'plotrow' not found.")
        }
    }else{
        if ( is.numeric(plotrow) ){
            if ( ! plotrow %in% 1:nrow(CircCoordinates) ){
                stop("Specified 'plotrow' not found.")
            }
        }else{
            stop("Specified plotrow should be ONE rowname or ONE rownumber.")
        }
    }
    # Choose your own column containing the gene name using gene_column. The genename will be displayed in the plot title if available
    if (missing(gene_column)){
        genename = NULL
    }else{
        genename <- as.character(CircCoordinates[plotrow,gene_column])
        if (genename == '.'){
            genename = NULL
        }
    }
    if(twolevel){
        plotdat <- summarySE( data.frame(Ratio=as.numeric(Circ[plotrow,-circle_description])/(as.numeric(Linear[plotrow,-circle_description])+as.numeric(Circ[plotrow,-circle_description])),
                                                                        groupindicator1,
                                                                        groupindicator2),
                                                measurevar='Ratio',groupvars=c('groupindicator1','groupindicator2') )
    }else{
        plotdat <- summarySE( data.frame(Ratio=as.numeric(Circ[plotrow,-circle_description])/(as.numeric(Linear[plotrow,-circle_description])+as.numeric(Circ[plotrow,-circle_description])),
                                                                        groupindicator1),
                                                                        measurevar='Ratio',groupvars=c('groupindicator1') )
    }
# construct plot
    Q <- ggplot(plotdat, aes(x=groupindicator1, y=Ratio)) +
            geom_boxplot() + theme_classic() +
            theme(axis.text.x = element_blank())+
            theme(axis.text.y = element_text(size=size+4))+
            theme(axis.ticks = element_line(colour = 'black', size = 1)) +
            theme(axis.ticks.x = element_blank())+
            theme(legend.title=element_blank()) +
            theme(text=element_text(size=size+4))+
            theme(legend.text=element_text(size=size)) +
            theme(plot.title = element_text(size=size)) +
            theme(axis.text.y = element_text(margin=margin(5,5,10,5,"pt")))+
            #labs(list(title=paste("Annotation: ", genename, "\nChr ", toString(Circ[plotrow,circle_description]),sep=""),x=x,y=y)) +
            ggtitle(paste("Annotation: ", genename, "\nChr ", toString(Circ[plotrow,circle_description]),sep="")) +
            ylab("circRNA/(circRNA + Linear RNA)") +
            xlab("Sample") +
            geom_errorbar(aes(ymin=Ratio, ymax=Ratio+se), width=.2 , size=2) +
            geom_bar(stat="identity",aes(fill=groupindicator1), color = "black", size=2)

    if (colour_mode == "bw"){
            Q <- Q + scale_fill_grey(start = 0.0, end = 1)
    } else {
            Q <- Q + scale_fill_discrete(name=lab_legend)
    }

            Q <- Q +
            theme(legend.position="bottom") +
            theme(axis.ticks.length = unit(0.5, "cm")) +
            theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.border = element_rect(colour = "black", fill=NA, size=3)) +
                guides(fill=guide_legend(
                                keywidth=0.3,
                                keyheight=0.3,
                                default.unit="inch")
            ) + scale_y_continuous(expand=c(0,0), limits= c(0, y_axis_range))

    if(twolevel){
        Q <- Q + facet_wrap( ~ groupindicator2,ncol=ncol )
    }

    print(Q)
}

## LINEPLOT

Circ.lineplot <- function(Circ,Linear,CircCoordinates = None,plotrow='1',size=18,ncol=2,groupindicator1=NULL,groupindicator2=NULL,x='Conditions',y='Counts', circle_description = c(1:3), gene_column = None){

    require(ggplot2)
    #require(Rmisc)

    if( !is.null(groupindicator1) & length(groupindicator1) != ncol(Circ)-length(circle_description) ){
        stop("If provided, the length of groupindicator1 should be equal to the number of samples.")
    }
    if( !is.null(groupindicator2) & length(groupindicator2) != ncol(Circ)-length(circle_description) ){
        stop("If provided, the length of groupindicator2 should be equal to the number of samples.")
    }
    if(is.null(groupindicator1)){
        stop("At least one grouping should be provided through groupindicator1.")
    }
    if(!is.null(groupindicator2)){
        twolevel <- TRUE
    }else{
        twolevel <- FALSE
    }

    rownames.circ <- rownames(Circ)
    Circ <- data.frame(lapply(Circ, as.character), stringsAsFactors=FALSE)
    rownames(Circ) <- rownames.circ

    rownames.linear <- rownames(Linear)
    Linear <- data.frame(lapply(Linear, as.character), stringsAsFactors=FALSE)
    rownames(Linear) <- rownames.linear

    # if CircCoordinates are available, use them, otherwise get more information from the Circ table, as indicated by the circle_description columns.
    if(!missing(CircCoordinates)){
        rownames.circCoordinates <- rownames(CircCoordinates)
        CircCoordinates <- data.frame(lapply(CircCoordinates, as.character), stringsAsFactors=FALSE)
        rownames(CircCoordinates) <- rownames.circCoordinates
    }else{
        CircCoordinates <- data.frame(Circ[,circle_description])
        rownames(CircCoordinates) <- rownames.circ
        rownames.circCoordinates <- rownames(CircCoordinates)
        CircCoordinates <- data.frame(lapply(CircCoordinates, as.character), stringsAsFactors=FALSE)
        rownames(CircCoordinates) <- rownames.circCoordinates
    }

    groupindicator1 <- factor(groupindicator1,levels=unique(groupindicator1))
    groupindicator2 <- factor(groupindicator2,levels=unique(groupindicator2))

    # Get gene name, if no annotation, output NULL
    if (is.character(plotrow)){
        if ( ! plotrow %in% rownames(CircCoordinates) ){
            stop("Specified 'plotrow' not found.")
        }
    }else{
        if ( is.numeric(plotrow) ){
            if ( ! plotrow %in% 1:nrow(CircCoordinates) ){
                stop("Specified 'plotrow' not found.")
            }
        }else{
            stop("Specified plotrow should be ONE rowname or ONE rownumber.")
        }
    }
    # Choose your own column containing the gene name using gene_column. The genename will be displayed in the plot title if available
    if (missing(gene_column)){
        genename = NULL
    }else{
        genename <- as.character(CircCoordinates[plotrow,gene_column])
        if (genename == '.'){
            genename = NULL
        }
    }

    plot.func <- function(row=plotrow){
        if(twolevel){
            plotdat <- summarySE(data.frame(Counts=c(as.numeric(Circ[row,-circle_description]),as.numeric(Linear[row,-circle_description])),
                                                                            groupindicator1,
                                                                            groupindicator2,
                                                                            Type=c(rep('circRNA',ncol(Circ)-length(circle_description)),rep('linear RNA',ncol(Circ)-length(circle_description)))
            ), measurevar='Counts',groupvars=c('Type','groupindicator1','groupindicator2') )
        }else{
            plotdat <- summarySE(data.frame(Counts=c(as.numeric(Circ[row,-circle_description]),as.numeric(Linear[row,-circle_description])),
                                                                            groupindicator1,
                                                                            Type=c(rep('circRNA',ncol(Circ)-length(circle_description)),rep('linear RNA',ncol(Circ)-length(circle_description)))
            ), measurevar='Counts',groupvars=c('Type','groupindicator1') )
        }

        Q=ggplot(plotdat, aes(x=groupindicator1, y=Counts, group=Type,colour=Type)) +
            theme(text=element_text(size=size))+
            theme_bw()+
            labs( list(title=paste(toString(Circ[row,circle_description]),genename,sep=" "),x=x,y=y) ) +
            ggtitle(paste(toString(Circ[row,circle_description]),genename,sep=" "))+
            geom_errorbar(aes(ymin=Counts-se, ymax=Counts+se), width=.1, position=position_dodge(.1) ) +
            xlab("Condition") +
            geom_line(position=position_dodge(.1)) +
            geom_point(position=position_dodge(.1))
        if (twolevel){
            Q = Q + facet_wrap( ~ groupindicator2,ncol=ceiling(sqrt(length(levels(groupindicator2)))) )
        }

        print(Q)
    }

    return(plot.func(row=plotrow))
}

## FILTER

Circ.filter <- function(circ=circ,linear=linear,Nreplicates=3,filter.sample=4,filter.count=5,percentage=1, circle_description=c(1:3)){
    del_row=c()
    for ( i in 1:nrow(circ) ){
        if ( sum(circ[i,-circle_description]>=filter.count)<filter.sample | circ_max_perc(circ[i,-circle_description], linear[i,-circle_description],Nreplicates=Nreplicates)<percentage)
            del_row = c(del_row,i)
    }

    if (length(del_row) > 0){
        new_dat=circ[-del_row,]
        return(new_dat)
    } else {
        return(circ)
    }
}

circ_max_perc <- function(circ=circ,linear=linear,Nreplicates=3){
    # convert to vector
    circ = as.numeric(circ)
    linear = as.numeric(linear)
    if( length(circ) != length(linear) ){
        stop ('Number of samples in circRNA is not equal to Hostgene.')
    }
    Ngroups = length(circ)/Nreplicates
    # calculate percentage
    circ_sum = unname(tapply(circ, (seq_along(1:length(circ))-1) %/% Nreplicates, sum ))
    linear_sum = unname(tapply(linear, (seq_along(1:length(linear))-1) %/% Nreplicates, sum ))
    perc = max(circ_sum / (circ_sum+linear_sum),na.rm=T)
    return(perc)
}

## CIRC TEST

Circ.test <- function(Circ, Linear, CircCoordinates=None, group, alpha=0.05, plotsig=T, circle_description = c(1:3)){

        # Requre packge
        require(aod)

        # check whether the input matrix are correct
        if ( nrow(Circ)!=nrow(Linear) | ncol(Circ) != ncol(Linear)){
                stop('Circ data and Linear data are not matched, dimention different.')
        }

        # A vector for pvalue and directions indicator
        p.val <- c()
        direction <- c()

        # groups
        if ( length(group) != ncol(Circ)-length(circle_description) ){
                stop("length of 'group' must be equal to the number of samples of 'Circ' and 'Linear'. ")
        }
        group <- factor(group)
        counter <- 0

        ## test
        # construct test matrix for each circRNA

        tmp_df = Circ[,FALSE]

        for (j in seq(1,length(unique(group)))){
                tmp_df[paste("group_",j,"_ratio_mean",sep="")] <- NA
        }

        for ( i in rownames(Circ) ){
                counter <- counter+1

                # total read counts vector
                tot <- round( as.numeric(Linear[i,-circle_description]) + as.numeric(Circ[i,-circle_description]) )

                # circRNA read counts
                circ <- as.numeric(Circ[i,-circle_description])

                # if there is 0 in the total count vector, the model will fail. So permute 0 to 1
                if ( 0 %in% tot ){
                    tot[tot==0]=1
                }

                if (counter %% 1000 == 0){
                        message(paste(counter, "candidates processed"))
                }

                tmp_rations <- data.frame(Ratio=as.numeric(Circ[i,-circle_description])/(as.numeric(Linear[i,-circle_description])+as.numeric(Circ[i,-circle_description])),
                group=group)
                for (rep_group in seq(1,max(as.numeric(levels(group))),1)){
                        tmp_df[i, paste("group_",rep_group,"_ratio_mean",sep="")] <- mean(na.omit(unlist(tmp_rations[tmp_rations$group==rep_group,1])))
                }

                # Constract data frame
                testdat = data.frame(tot,circ,group)

                ## do test
                # Null model
                fitNull <- betabin(cbind(circ,tot-circ) ~ 1, ~ 1, data=testdat)
                # Alternative model
                fitAlt <- betabin(cbind(circ,tot-circ) ~ group, ~ 1, data=testdat)
                # test models
                a <- anova(fitNull,fitAlt)
                p.value <- a@anova.table[,11][2]
                # print(predict(fitAlt,testdat, se.fit=T))
                p.val <- c( p.val, p.value )
                # dir <- 1 # fitAlt@param[2][["group2"]]
                # direction <- c(direction, dir)
        }
        message(paste(counter, "candidates processed in total"))

        Circ$direction <- direction
        #names(Circ$direction ) <- c("direction")
        p.adj <- p.adjust(p.val,n=sum(!is.na(p.val)),'BH')
        # select significant ones
        sig_dat <- Circ[p.adj<=alpha    & !is.na(p.adj),]
        sig_ratios <- tmp_df[p.adj<=alpha    & !is.na(p.adj),]
        sig_p <- p.adj[p.adj<=alpha    & !is.na(p.adj)]
        # direction <- direction[p.adj<=alpha    & !is.na(p.adj)]

        # sort by p-val
        sig_dat <- sig_dat[order(sig_p),]
        sig_ratios <- sig_ratios[order(sig_p),]
        sig_p <- sort(sig_p)

        # A summary table
        if (missing(CircCoordinates)){
                summary_table <- data.frame(sig_dat[,circle_description],sig_p,sig_dat[,circle_description])

                rownames(summary_table) <- rownames(sig_dat)
                names(summary_table) <- c(names(sig_dat)[circle_description],"sig_p",names(sig_ratios)[circle_description])
        } else {
                # summary_table <- cbind(CircCoordinates[rownames(sig_dat),],sig_p,sig_dat$direction)
                # colnames(summary_table) <- c(colnames(CircCoordinates),"sig_p","direction")

                summary_table <- cbind(CircCoordinates[rownames(sig_dat),],sig_p,sig_ratios)
                colnames(summary_table) <- c(colnames(CircCoordinates),"sig_p",colnames(sig_ratios))
        }

        message(paste(nrow(summary_table), "candidates passed the specified thresholds"))

        # return all objects in a list
        return(list(summary_table=summary_table,
                            sig.dat=sig_dat,
                            p.val=p.val,
                            p.adj=p.adj,
                            sig_p=sig_p,
                            ratios=sig_ratios
                            # direction=direction
                        )
                )
}

## MY CODE

args = commandArgs(trailingOnly = TRUE)

circ = read.table(args[1], header=T, sep=",")
linear = read.table(args[2], header=T, sep=",")
pheno = read.table(args[3], header=T, sep=",", row.names = "Sample_ID")


# No need to enforce any filtering for circTest.
# 'filter.sample' - this has been applied to called circs using the 'tool.filter' param
# 'filter.count' - this has been applied to called circs using bsj_filter param
# 'percentage' - set to extremely low value (do not want to discard circRNAs - let the user inspect themselves).

# Need to apply the phenotype csv file correctly to circtest.
n_covars <- ncol(pheno)
if( n_covars == 2){
    covariate_1 <- as.factor(pheno[,1])
    covariate_2 <- as.factor(pheno[,2])
}else{
    covariate_1 <- as.factor(pheno[,1])
}

n_reps <- as.numeric(table(covariate_1)[1])

Circ_filtered <- Circ.filter(circ = circ, linear = linear, Nreplicates = n_reps, filter.sample = 1, filter.count = 1, percentage = 0.00001, circle_description = c(1:4))
Linear_filtered <- linear[rownames(Circ_filtered),]


# groups must be numerically encoded
group = as.numeric(covariate_1)
test <- Circ.test(Circ_filtered, Linear_filtered, group=group, circle_description = c(1:4))
write.table(test$summary_table, "summary_table.txt", row.names=F)


# Apply pheno to output once more..

if( n_covars == 2 ){

    group_indicator1 <- as.character(covariate_1)
    group_indicator2 <- as.character(covariate_2)

    pdf("circ_linear_ratio_plots.pdf", width = 8, height = 10)
    for (i in rownames(test$summary_table))    {
        Circ.ratioplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=groupindicator1, groupindicator2 = group_indicator2,
        circle_description = c(1:4) )
    }
    dev.off()

    pdf("circ_linear_line_plots.pdf", width = 8, height = 10)
    for (i in rownames(test$summary_table))    {
    Circ.lineplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=group_indicator1, groupindicator2 = group_indicator2,
    circle_description = c(1:4) )
    }
    dev.off()

}else{

    group_indicator1 <- as.character(covariate_1)

    pdf("circ_linear_ratio_plots.pdf", width = 8, height = 10)
    for (i in rownames(test$summary_table))    {
        Circ.ratioplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=group_indicator1,
        lab_legend = colnames(pheno)[1],    circle_description = c(1:4) )
    }
    dev.off()

    pdf("circ_linear_line_plots.pdf", width = 8, height = 10)
    for (i in rownames(test$summary_table))    {
    Circ.lineplot(Circ_filtered, Linear_filtered, plotrow=i, groupindicator1=group_indicator1,
    circle_description = c(1:4) )
    }
    dev.off()
}

# include variables, makes life easier in case user wishes to report bugs to workflow.
save.image("circ_test.RData")
