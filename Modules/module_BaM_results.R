require(ggplot2); require(reshape2); require(GGally);




#######################################################################################################
plot.Results_BaM <- function(workspace, 
                             iter, 
                             plot.RC.lin, plot.RC.log, plot.RC.loglog, 
                             grid_RC.ylim, grid_RC.xlim, 
                             grid_RC.xstep, grid_RC.ystep,
                             ticks_RC.y.log, 
                             RC.x.labels, RC.y.labels,
                             u.m.Qgaug, u.m.Hgaug,
                             ncontrol,
                             hP, QP, uQP,
                             h_Gaug, Q_Gaug, uQ_Gaug) {
#######################################################################################################
  
  #-----------------------------------------
  # TUNINGS
  #-----------------------------------------
  #Rfolder="C:\\BEN\\FORGE\\BaM\\trunk\\R"
  file.X        = 'Hgrid.txt'
  file.env.tot  = 'Qrc_TotalU.env'
  file.spag.tot = 'Qrc_TotalU.spag'
  file.env.par  = 'Qrc_ParamU.env'
  file.spag.par = 'Qrc_ParamU.spag'
  file.spag.max = 'Qrc_Maxpost.spag'
  file.mcmc     = 'Results_MCMC_Cooked.txt'
  file.summary  = 'Results_Summary.txt'

  
  
  #--------------------------------------------
  # Preliminaries
  #--------------------------------------------
  setwd(workspace)
  model = ReadModel(ModelFile = "Config_Model.txt")

  
  
  
  #------------------------------------------------------------------
  # MCMC
  #------------------------------------------------------------------
  vertical.length = length(model$par)*4
  mcmc = MCMCplot(doLogPost = T,
                  doPar     = T,
                  doDPar    = T, 
                  MCMCfile  = "Results_MCMC_Cooked.txt" , 
                  type      = "trace",  #="trace", # "histogram","density","scatterplot"
                  xlab      = '',
                  ylab      = '',
                  ncol      = 1, 
                  prior     = NULL,
                  burn      = 0, 
                  slim      = 1,
                  theme_size= 15)
         # ggsave(mcmc, filename =paste(workspace,"/mcmc_it",iter,".png", sep=""), 
         #        bg = "transparent", width = 12, height =vertical.length, dpi = 100)
         
  mcmc2 = MCMCplot(doLogPost = T,
                   doPar     = T,
                   doDPar    = T, 
                   MCMCfile  = "Results_MCMC_Cooked.txt" , 
                   type      = 'density', 
                   prior     = model$par,
                   xlab      = '',
                   ylab      = '',
                   ncol      = 1, 
                   burn      = 0, 
                   slim      = 1,
                   theme_size= 15)  
          # ggsave(mcmc2, filename =paste(workspace,"/mcmc2_it",iter,".png", sep=""), 
          # bg = "transparent", width = 12, height =vertical.length, dpi = 100)
  pdf(paste0(workspace,"/mcmc_it",iter, ".pdf"), 17, vertical.length, useDingbats=F)
  print(plot_grid(mcmc, mcmc2,
                  nrow=1, ncol = 2, 
                  labels = c("Trace plots", "Density plots"),
                  label_size = 20))
  dev.off()
  
  
  
  
  
  
  #-------------------------------------------------------------------
  # CORRELATIONS BETWEEN PARAMETERS
  #-------------------------------------------------------------------
  jpeg(filename  = paste0(workspace,"/correlation_matrix_it",iter,".jpeg"), width = 12, height = 7, units = 'in', res = 100)
  correlation.matrix = correlation.plot(MCMCfile = "Results_MCMC_Cooked.txt" , 
                                        npar     = (model$nPar+3))
  dev.off()
  dev.set(dev.prev())
  

  
  
  
  
  
  
  #------------------------------------------------------------------
  # Residual analysis
  #------------------------------------------------------------------
  # g=Residualplot(nX=model$nX,nY=model$nY)
  # for(i in 1:length(g)){
  #    if(!is.null(g[[i]])){print(g[[i]])}
  # }
  
  
  
  
  
  
  #------------------------------------------------------------------
  # Plot predictions
  #------------------------------------------------------------------
  MAPcurve        <-  "MAP rating curve"
  Ustruct.MAP     <- "95% remnant uncertainty of MAP curve"
  Utot            <- "95% total uncertainty"
  gaugings.period <- "Gaugings of current period"
  gaugings.other  <- "Gaugings of other periods"
  activ.stage     <- "Control activation stages (MAP)"
  col.uncert      <- c("95% total uncertainty"                = "mistyrose1", 
                       "95% remnant uncertainty of MAP curve" = "blue") #"salmon") #"indianred1") 
  cols            <- c("Gaugings of current period"           = "black",  
                       "Gaugings of other periods"            = "gray90",
                       "MAP rating curve"                     = "blue", 
                       "Control activation stages (MAP)"      = "green")
  #cols <- c("Current periods gaugings"= "black", "Other periods gaugings" = "gray", 
  #          "MAP rating curve"= "blue") 
  #al <- c("U_tot"=0.2, "U_par"=0.4, "U_MAP"=0.3)    #"U_par"="indianred1",
  maxpost      <- read.table(file = file.spag.max)
  maxpost.par  <- read.table(file = file.summary)
  Hg           <- read.table(file = file.X)
  mcmc.par     <- read.table(file = file.mcmc, header = TRUE)
  #----------------------------------------------------------------------------------------
  ktransit=NULL;
  for (i in 1:ncontrol) {
    ktransit[[i]] =0
    ktransit[[i]][1] = maxpost.par[16, ncontrol*3+2+i]
    ktransit[[i]][2]  = c(quantile(mcmc.par[,ncontrol*3+3+i], p = c(0.025)))[[1]]
    ktransit[[i]][3] = c(quantile(mcmc.par[,ncontrol*3+3+i], p = c(0.975)))[[1]]
  }
  ktransit.df = data.frame(t(sapply(ktransit, function(x) x[1:max(lengths(ktransit))])))
  names(ktransit.df) = c("MAP", "q2", "q97")

  #----------------------------------------------------------------------------------------
  gamma1    <- maxpost.par[16, ncontrol*3+1];
  str.err   = 0
  if (remnant.err.model == "Linear") {  
    gamma2  <- maxpost.par[16,ncontrol*3+2]
    str.err = gamma1+gamma2*maxpost
  } else {
    str.err = gamma1
  }
  
  #----------------------------------------------------------------------------------------
  env.tot = EnvelopLayer(file  = file.env.tot, 
                         Xfile = file.X, 
                         color = Utot, 
                         alpha = 1)
  env.par = EnvelopLayer(file  = file.env.par, 
                         Xfile = file.X, 
                         color = "U_par", 
                         alpha = 1)
  #spag = SpaghettiLayer(file=file.spag, Xfile=file.X)
  #env.max = EnvelopLayer(file=file.spag.max, Xfile=file.X , color = "black", alpha =1 )
  maxpost.data <- data.frame(x = Hg[which(Hg >= ktransit.df$MAP[1]),1], 
                             y = maxpost[which(Hg >= ktransit.df$MAP[1]),1], 
                             y2 = str.err[which(Hg >= ktransit.df$MAP[1]),1])
  new_line <- element_line(color = "grey", size = 0.1, linetype = 2)
  
  text.ktrans = c()
  for (kkk in 1:length(ktransit.df$MAP)){
      text.ktrans[kkk] = paste0("k",kkk)
  }
  
  
  
  
  
  
  
  
  
  
  
  ################################################################################# 
  #                               Rating Curve plot
  #################################################################################
  #par(family = "LM Roman 10")
  loadfonts()
  plot.Utot = FALSE
  plot.RC.lin          = TRUE             # RC plot in linear scale
  plot.RC.log          = TRUE             # RC plot in logarithmic scale
  plot.RC.loglog       = FALSE            # if "TRUE" you will need to change the x-axis ticks and limits
  
  g = p = NULL
  
  
  if (plot.RC.lin == TRUE) { 
  #LINEAR PLOT:
  #--------------------------------------------------------------------------------
  g = ggplot() +
               coord_cartesian(xlim = grid_RC.xlim, ylim = grid_RC.ylim, expand = c(0))+
               scale_y_continuous(breaks = seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep))+
               scale_x_continuous(breaks = c(seq(grid_RC.xlim[1],  grid_RC.xlim[2],  grid_RC.xstep)), #, ktransit),
                                  labels = c(seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep))) + #,text.ktrans)) + 
               annotate(geom = "rect", ymin = -Inf, ymax= Inf, xmin = unlist(ktransit.df$q2), 
                                           xmax = unlist(ktransit.df$q97), fill= "green", alpha = 0.1)+
               if (plot.Utot == TRUE ) {
                 g = g +  env.tot  # plot total RC uncertainty (at 95%)
               }
               g = g+ 
               geom_ribbon(data = na.omit(maxpost.data), aes(x = x, ymin = (y-y2) ,  #total U of the MAP
                           ymax = (y+y2), fill = Ustruct.MAP), alpha = 0.1) +
               geom_errorbar(aes(x=h_Gaug, ymin =Q_Gaug-2*uQ_Gaug, ymax =Q_Gaug+2*uQ_Gaug, color = gaugings.other), 
                             width=.05, size = 0.2) +
               geom_point(aes(x = h_Gaug, y = Q_Gaug, color = gaugings.other), pch =21, fill ="gray90", size = 1.2) +
               geom_errorbar(aes(x=hP, ymin =QP-2*uQP, ymax =QP+2*uQP, color = gaugings.period), 
                             width=.05, size = 0.2) +
               geom_point(aes(x = hP, y = QP, color =gaugings.period), pch =21, fill = "black", size = 1.3) +
               ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
               xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]"))) + 
               geom_line(data = maxpost.data, aes(x = x, y= y, color = MAPcurve), size = 0.8) +  # MAP curve
               geom_vline(aes(xintercept = ktransit.df$MAP, color = activ.stage),linetype="dashed", size = 0.4) + 
               scale_fill_manual(name=element_blank(), values= col.uncert, breaks=c(Utot,  Ustruct.MAP)) +
               scale_colour_manual(name    = element_blank(), 
                                   values  = cols ,
                                   breaks  = c(gaugings.period, gaugings.other, MAPcurve, activ.stage),
                                   #labels = c(gaugings.period,gaugings.other,MAPcurve,activ.stage),
                                   guide   = guide_legend(override.aes = list(
                                             linetype = c("blank", "blank","solid","solid"),
                                             shape = c(19, 19, NA, NA))) ) +
              theme_light(base_size = 10) +
              theme( #text               = element_text(size=10, family="LM Roman 10")
                    plot.background     = element_rect(fill ="transparent", color = NA)
                    ,panel.grid.major   = element_blank()
                    ,panel.grid.minor   = element_blank()
                    ,panel.background   = element_rect(fill ="transparent")
                    ,axis.ticks         = element_line(colour = "black")
                    ,plot.margin        = unit(c(0.3, 0.5, 0.05, 1.05),"cm")
                    #,axis.text.x        = element_text(colour = c(rep('black', length(seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep)))))
                    ,axis.title.y       = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0))
                    ,legend.position    = "none")
                    # ,axis.line         = element_line(colour = "black")
                    # ,panel.grid.major  = element_line(size=0.4, linetype = "dashed")
                    # ,panel.grid.minor  = element_blank()
                    # ,plot.title        = element_text(hjust = 0.5)
                    # ,panel.grid.minor  = new_line,
                    # ,legend.position   = "bottom", legend.direction="vertical",
                    # ,legend.position   = c(0.3, 0.8), legend.box = "vertical", legend.key.size = unit(0.5, "cm"),
                    # ,legend.text       = element_text(size=12),
                    # ,legend.spacing.x  = unit(0.05, 'cm'),
                    # ,legend.spacing.y  = unit(0.05, 'cm'),
                    # ,legend.key.height = unit(0.5,"line"),
                    # ,legend.key.width  = unit(0.8,"line"),
                    # ,legend.key        = element_rect(colour = "transparent", fill = "transparent"),
                    # ,legend.background = element_rect(colour = "transparent", fill = "transparent"))
               
              ggsave(plot = g, filename =paste(workspace,"/RC_it", iter,".png", sep=""), 
                     bg = "transparent", width = 6, height =5, dpi = 400)
     }         
            
  
  
        if (plot.RC.log == TRUE) {            
        #LOG PLOT
        #---------------------------------------------------------------------------------------
        p <- ggplot() + 
             coord_trans(xlim = grid_RC.xlim, 
                         ylim = grid_RC.ylim.log) +
             scale_y_log10(na.value = -10, 
                           breaks   = ticks_RC.y.log, 
                           labels   = ticks_RC.y.log) +
             scale_x_continuous(breaks = c(seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep)),  #, ktransit),
                                labels = c(seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep))) + # , text.ktrans)) +
             annotate(geom = "rect", 
                      ymin = -Inf, ymax= Inf, xmin = unlist(ktransit.df$q2), xmax = unlist(ktransit.df$q97), 
                      fill= "green", alpha = 0.1)
        
             if (plot.Utot == TRUE ) {
               p = p +  env.tot  # plot total RC uncertainty (at 95%)
             }
             p = p +
             theme_light(base_size = 10)+
             geom_ribbon(data = na.omit(maxpost.data), 
                         aes(x = x, ymin =(y-y2) , 
                         ymax = (y+y2), fill = Ustruct.MAP ), 
                         alpha = 0.1) +
             geom_errorbar(aes(x=h_Gaug, ymin =Q_Gaug-2*uQ_Gaug, ymax =Q_Gaug+2*uQ_Gaug, color =gaugings.other), 
                           width=.02, size = 0.2)+
             geom_point(aes(x = h_Gaug, y = Q_Gaug, color =gaugings.other),    size = 1.2) +
             geom_errorbar(aes(x=hP, ymin =QP-2*uQP, ymax =QP+2*uQP, color =gaugings.period), 
                           width=.02, size = 0.2)+
             geom_point(aes(x = hP, y = QP, color = gaugings.period),  size = 1.2) +
             ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
             xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]"))) + 
             annotation_logticks(base = 10, sides = "l", scaled = TRUE,colour = "black", size = 0.3, linetype = 1)+
             geom_line(data = na.omit(maxpost.data),
                       aes(x = x, y= y, color = MAPcurve),
                       size = 0.8, na.rm = TRUE)+
             geom_vline(aes(xintercept  = ktransit.df$MAP, color= activ.stage), size = 0.4, linetype ="dashed")+   
             scale_fill_manual(name     = element_blank(), 
                               values   = col.uncert) +
             scale_colour_manual(name   = element_blank(), 
                                 values = cols,
                                 breaks = c(gaugings.period,
                                            gaugings.other,
                                            MAPcurve,
                                            activ.stage),
                                 #labels = c(gaugings.period,
                                 #           gaugings.other,
                                 #           MAPcurve,
                                 #           activ.stage),
                                 guide = guide_legend(override.aes = list(
                                         linetype = c("blank", "blank","solid","solid"),
                                         shape = c(19, 19, NA, NA)))) +
          
             theme( plot.title          = element_text(hjust = 0.5)
                   ,plot.background     = element_rect(fill ="transparent", color = NA)
                   #,panel.grid.major   = element_line(size=0.4, linetype = "dashed")
                   ,panel.grid.minor    = element_blank()
                   ,panel.grid.major    = element_blank()
                   ,panel.background    = element_rect(fill ="transparent") 
                   #,axis.line          = element_line(colour = "black")
                   #,axis.text.x         = element_text(colour = c(rep('black', length(seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep))), 'green'))
                   ,axis.ticks          = element_line(colour = "black")
                   ,plot.margin         = unit(c(0.3,0.5,0.05,1.05),"cm") 
                   ,legend.position     = "none")
                   #,legend.position    = c(0.7, 0.2)
                   #,legend.box         = "vertical"
                   #,legend.key.size    = unit(0.5, "cm")
                   #,legend.position    = "bottom"
                   #,legend.direction   = "vertical"
                   #,legend.text        = element_text(size=12)
                   #,legend.spacing.x   = unit(0.05, 'cm')
                   #,legend.spacing.y   = unit(0.01, 'cm')
                   #,legend.key.height  = unit(0.5,"line")
                   #,legend.key.width   = unit(0.8,"line")
                   #,legend.key         = element_rect(colour = "transparent", fill = "transparent")
                   #,legend.background  = element_rect(colour = "transparent", fill = "transparent"))
             ggsave(plot = p, filename =paste0(workspace,"/RClog_it", iter,".png"), 
                    bg = "transparent", width = 6, height =5, dpi = 400)

        }    
             
        
  
  
  
  
  
        if (plot.RC.loglog == TRUE) {      
        #LOG LOG PLOT
        #---------------------------------------------------------------------------------------
        ploglog <- ggplot() + 
               coord_trans(xlim = c(Qminlog, Qmaxlog), ylim = c(Qminlog, Qmaxlog))+
               scale_y_log10(na.value=-10, breaks=ticks.log, labels=ticks.log) +
               scale_x_log10(na.value=-10, breaks=ticks.log, labels=ticks.log) +      
               #env.tot +
               theme_light(base_size = 10)+
               geom_ribbon(data = na.omit(maxpost.data), aes(x = x, ymin =(y-y2), ymax = (y+y2), fill = Ustruct.MAP ), alpha = 0.2) +
               geom_point(aes(x = h_Gaug, y = Q_Gaug, color =gaugings.other), size = 0.9) +
               geom_point(aes(x = hP, y = QP, color =gaugings.period), size = 0.9) +
               geom_errorbar(aes(x=h_Gaug, ymin =Q_Gaug-2*uQ_Gaug, ymax =Q_Gaug+2*uQ_Gaug, color =gaugings.other), width=.02, size = 0.2)+
               geom_errorbar(aes(x=hP, ymin =QP-2*uQP, ymax =QP+2*uQP, color =gaugings.period), width=.02, size = 0.2)+
               ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
               xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]"))) + 
               annotation_logticks(base = 10, sides = "l", scaled = TRUE,colour = "black", size = 0.3, linetype = 1)+
               geom_line(data = na.omit(maxpost.data), aes(x = x, y= y, color = MAPcurve),size = 0.8, na.rm = TRUE)+
               geom_vline(aes(xintercept = ktransit, color =activ.stage), size = 0.3, linetype ="dashed")+   
               scale_fill_manual(name=element_blank(), values= col.uncert) +
               scale_colour_manual(name = element_blank(), values=cols,
                                   breaks=c(gaugings.period,gaugings.other,MAPcurve,activ.stage),
                                   #labels = c(gaugings.period,gaugings.other,MAPcurve,activ.stage),
                                   guide = guide_legend(override.aes = list(
                                   linetype = c("blank", "blank","solid","solid"),
                                   shape = c(19, 19, NA, NA)))) +
               theme(plot.title = element_text(hjust = 0.5), 
                     plot.background = element_rect(fill ="transparent", color = NA),
                     #panel.grid.major=element_line(size=0.4, linetype = "dashed"), 
                     panel.grid.minor= element_blank(),
                     panel.grid.major = element_blank(),
                     panel.background = element_rect(fill ="transparent"), 
                     #axis.line = element_line(colour = "black"),
                     axis.ticks = element_line(colour = "black"),
                     plot.margin=unit(c(0.3,0.5,0.05,1.05),"cm"), 
                     #legend.position =c(0.7, 0.2), legend.box = "vertical", legend.key.size = unit(0.5, "cm"),
                     legend.position = "none")
               
                     ggsave(plot = ploglog, filename =paste(workspace,"/RCloglog_it", iter,".png", sep=""), 
                     bg = "transparent", width = 6, height =5, dpi = 400)             
        }
            
    #LEGEND:  
    #*************************************************************************************************       
    # ggsave(plot = g+scale_y_log10(breaks=c(0.01,0.1,1,10,100,500),labels=c(0.01,0.1,1,10,100,500)), 
    #        filename =paste(workspace,"/RClog_it", iter,".png", sep=""), 
    #        bg = "transparent", width = 6, height =5, dpi = 300)
    # glegend <- cowplot::get_legend(gg)         
    #*************************************************************************************************
    #figure <- ggarrange(g,p, glegend, labels = c("a) ", "b) ", "Legend "),  ncol = 1, nrow = 3)
    # figure <- plot_grid( g,p, glegend, labels = c("a) ", "b) ", "Legend "), hjust = -1, 
    #                     ncol = 1, nrow = 3, rel_heights = c(3, 3, 1.5))
    # ggsave(plot = figure, filename =paste(workspace,"/RC_both_it", iter,".png", sep=""), 
    #        bg = "transparent", width = 6, height =8, dpi = 300)

#-----------------------------------------------
list(g, p, mcmc, mcmc2)
#-----------------------------------------------
}




















###########################################################################################
Results_env.tot <- function(file.env.tot, file.X, file.model, colo) {
###########################################################################################
  #setwd(workspace)
  #source("BaMplots.R")
  model=ReadModel(file.model)
  env.tot = EnvelopLayer(file=file.env.tot, Xfile=file.X, color = colo, alpha = 0.3)
}







#################################################################################
Results_env.par <- function(file.env.par, file.X, file.model, colo) {
#################################################################################
  model=ReadModel(file.model)
  env.par = EnvelopLayer(file=file.env.par, Xfile=file.X , color = colo, alpha =0.5 )
  
}












#################################################################################
plot.mcmc <- function(workspace, iter) {
#################################################################################  
  setwd(workspace)
  file.summary = 'Results_Summary.txt'
  model=ReadModel('Config_model.txt')
  mcmc=MCMCplot(MCMCfile=paste0(workspace,"/Results_MCMC_Cooked.txt"), doLogPost=T , doDPar=T, type= "trace")
  ggsave(mcmc, filename =paste0(workspace,"/mcmc_P",iter,".png"), bg = "transparent", width = 12, height =7, dpi = 100)
  mcmc2=MCMCplot(MCMCfile=paste0(workspace,"/Results_MCMC_Cooked.txt"), doLogPost=T, doDPar=T, type='density', prior=model$par)
  ggsave(mcmc2, filename =paste0(workspace,"/mcmc2_P",iter,".png"), bg = "transparent", width = 12, height =7, dpi = 100)
}















########################################################################################
plot.mcmc.segment <- function(workspace, seg.iter, nS) {
######################################################################################
  setwd(workspace)
  file.summary    = 'Results_Summary.csv'
  model           = ReadModel('Config_Model.txt')
  vertical.length = length(model$par)*4  #Npar = Nseg*2
  mcmc = MCMCplot(doLogPost = T,
                  doPar     = T,
                  doDPar    = F, 
                  MCMCfile  = "Results_MCMC_Cooked.txt" , 
                  type      = "trace",  #="trace", # "histogram","density","scatterplot"
                  xlab      = '',
                  ylab      = '',
                  ncol      = 1, 
                  prior     = NULL,
                  burn      = 0, 
                  slim      = 1,
                  theme_size= 15)
  # ggsave(mcmc, filename =paste0(workspace,"/mcmc_it",seg.iter,".png"), 
  #        bg = "transparent", width = 12, height =vertical.length, dpi = 100)
  mcmc2 = MCMCplot(doLogPost = T,
                   doPar     = T,
                   doDPar    = F, 
                   MCMCfile  = "Results_MCMC_Cooked.txt" , 
                   type      = 'density', 
                   prior     = model$par,
                   xlab      = '',
                   ylab      = '',
                   ncol      = 1, 
                   burn      = 0, 
                   slim      = 1,
                   theme_size= 15)  
  # ggsave(mcmc2, filename =paste0(workspace,"/mcmc2_it",seg.iter,".png"), 
  #        bg = "transparent", width = 12, height =vertical.length, dpi = 100)
  pdf(paste0(workspace,"/mcmc_it",seg.iter, ".pdf"), 20, vertical.length, useDingbats=F)
  print(plot_grid(mcmc, mcmc2,
                  nrow=1, ncol = 2, 
                  labels = c("Trace plots", "Density plots"),
                  label_size = 24))
  dev.off()
  
          
          
          
          
  #----------------
  # Plot residuals:
  #----------------
  # g=Residualplot(nX=model$nX,nY=model$nY)
  # for(i in 1:length(g)){
  #     if(!is.null(g[[i]])){print(g[[i]])}
  # }
  
          
  #-------------------------------------------------------------------
  # CORRELATIONS BETWEEN PARAMETERS
  #-------------------------------------------------------------------
  if (nS >1) {
    pdf(paste0(workspace,"/correlation_matrix_it",seg.iter, ".pdf"), 14, 10, useDingbats=F)
    correlation.matrix = correlation.plot(MCMCfile = paste0(workspace,"/Results_MCMC_Cooked.txt"), 
                                          npar     = (model$nPar+2))
    dev.off()
    dev.set(dev.prev())
  }  
  
}









# FUNCTION NOT USED !!!  TO BE VERIFIED !!!
########################################################################################
Gelman.test <- function(x,n , m) {
########################################################################################
# Gelman function for the test of convergence of the MCMC:
#--------------------------------------------------------
# x = all chains combined
# n= number of iterations
# B/n = empirical between-chain variance
# W = the mean of the empirical variance within each chain
# var = the empirical variance from all chains combined
# m = number of chains
  
sigma2 =   (n-1)*W/n + B/n  #empirical variance from all chains combined
mu = mean(x)                #Sample mean of all chains combined
V = sigma2 + B/(n*m)        #Sample variance  of all chains combined
d = 2* v^2/var(V)           #degrees of freedom estimated by the method of moments

R = sqrt(((d+3)*V)/((d+1)*W))
return(R)
}

















################################################################################################################
Convergence.test <- function(dir.seg, npar, dir.plot) {
################################################################################################################
  mcmc.tot    <- read.table(paste0(dir.seg,"/Results_MCMC_Cooked.txt"), header = TRUE)
  len         <- length(mcmc.tot[,1])
  gel         =  NULL
  convergence =  TRUE
  # npar = ncontrols*(b,a,c) + (g1,g2)
  height.plot = npar*1
  
  #***************************************************************************************************************
  if (npar == 2){
    nroww = 1
    ncoll = 3
  } else {
    if (((npar+1) %% 5) ==0) {
      nroww = (npar+1)/5
      ncoll= 5
    } else { 
      nroww = ceil((npar+1)/5)
      ncoll = 5
    }
  }
  
  pdf(paste0(dir.plot,"/autocorrelation.pdf"), width = 10, height = 3*nroww, useDingbats=F)
  par(mfrow = c(nroww, ncoll))   #c(nrows, ncol)
  for (i in 1:npar) {
    par(mar = c(3, 3, 3, 3))  # Set the margin on all sides 
    autocorr.plot(x = mcmc.tot[,i],  auto.layout = FALSE)
    title(main = names(mcmc.tot)[i], cex.main = 2)
  }
  dev.off()
 # dev.set(dev.prev())
  #*************************************************************************************************************
  # jpeg(paste0(dir.plot,"/autocorrelation.jpg"), width = 10, height = height.plot, units = 'in', res = 400)
  # #if (((npar+1) %% 2) ==0) {ncoll = (npar+1)/2} else { ncoll = (npar+1)/2 + 0.5}
  # par(mfrow = c(ceil((npar+3)/3), 3))   #c(nrows, ncol)
  # for (i in 1:(ncol(mcmc.tot))) {
  #           par(mar = c(3, 3, 3, 3)) # Set the margin on all sides 
  #           autocorr.plot(x = mcmc.tot[,i], auto.layout = FALSE)
  #           title(main = names(mcmc.tot)[i], cex.main = 2)
  # }
  # dev.off()
  # dev.set(dev.prev())
  # 
  #*************************************************************************************************************
  file.name = paste0(dir.plot,"/convergence.csv")
  cat("Rc;", file = file.name, sep=" ")
  cat("Ru",  file = file.name, append = TRUE, sep="\n")
  #*************************************************************************************************************
  pdf(paste0(dir.plot,"/gelman.pdf"), width = 10, height = height.plot, useDingbats=F)
  par(mar=c(1,1,1,1))
  #if (((npar+1) %% 2) ==0) {ncoll = (npar+1)/2} else { ncoll = (npar+1)/2 + 0.5}
  par(mfrow = c(ceil((npar)/3), 3))   #c(nrows, ncol)
  for (i in 1:(npar)) {
    mcmc.list  = list(mcmc.tot[1:(len/4),i], 
                      mcmc.tot[(len/4+1):(len/2),i], 
                      mcmc.tot[(len/2+1):(len*3/4),i], 
                      mcmc.tot[(len*3/4+1):len,i])
    mcmc.list  = as.mcmc.list(lapply(mcmc.list, mcmc))
    gel[[i]]   = gelman.diag(x = mcmc.list, confidence=0.95, transform=FALSE)
    par(mar    = c(3, 3, 3, 3)) # Set the margin on all sides 
    gelman.plot(x = mcmc.list,auto.layout = FALSE) 
    title(main = names(mcmc.tot)[i], cex.main = 2)
    
    #acfplot(x = mcmc.list, lag.max=200)
    if (gel[[i]]$psrf[1] > 1.2) { 
          convergence = FALSE
    }
    cat(gel[[i]]$psrf[1], file = file.name, append = TRUE,sep="")
    cat(";", file = file.name, append = TRUE,sep=" ")
    cat(gel[[i]]$psrf[2], file = file.name, append = TRUE,sep="\n")
  }
  dev.off()
  #dev.set(dev.prev())
  #************************************************************************************************************
  Rc= "Reduction factor Rc"; Ru="Upper confidence limit Ru"; limit= "threshold";
  col.gelman = c("Reduction factor Rc"       = "blue", 
                 "Upper confidence limit Ru" = "red", 
                 "threshold"                 = "red")
  Rc.Ru.plot <- ggplot()+
    coord_cartesian(xlim = c(1, length(gel)))+
    #scale_y_continuous(breaks=seq(0,Qmax,Qmax/5))+
    scale_x_continuous(breaks=seq(1, length(gel), 1)) +
    geom_point(mapping = aes(x=seq(1,length(gel),1), unlist(gel, use.names = FALSE)[seq(1,length(gel)*2,2)] , 
                             col=Rc),shape=0, size = 4, stroke = 0.5) +
    geom_point(mapping = aes(x=seq(1,length(gel),1), unlist(gel, use.names = FALSE)[seq(2,length(gel)*2,2)], 
                             col=Ru), shape=4, size = 4, stroke = 0.5) +
    theme_bw(base_size = 15)+
    geom_hline(aes(yintercept = 1.2, col= limit), lwd =1, linetype = "dashed") +
    scale_colour_manual(name = "Legend", values=col.gelman ,
                        breaks=c(Rc, Ru, limit),
                        guide = guide_legend(override.aes = list(
                                linetype = c("blank", "blank", "dashed"),
                                shape = c(0, 4, NA)), 
                                title.position="top")) +
    xlab("Parameter index")+ 
    ylab("Gelman factor")+
    theme(text                  = element_text(size=15)
          ,plot.background      = element_rect(fill ="transparent", color = NA)
          ,panel.grid.major     = element_line(size=0.4, linetype = "dashed")
          ,panel.grid.minor     = element_blank()
          ,legend.position      = "bottom"
          ,legend.text          = element_text(size = 10) 
          ,legend.justification = c(0,1)
          ,legend.box           = "vertical"
          ,legend.box.just      = "left"
          ,legend.text.align    = 0)
  ggsave(Rc.Ru.plot, filename = paste0(dir.plot,"/Rc_Ru_gelman.png"), 
         bg = "transparent", device = "png", width = 6, height =4, dpi = 400)
  return(convergence)
}















################################################################################################################
Convergence.test.segment <- function(dir.seg, npar, dir.plot) {
################################################################################################################
  mcmc.tot    =  read.table(paste0(dir.seg,"/Results_MCMC_Cooked.txt"), header = TRUE)
  len         =  length(mcmc.tot[,1])
  gel         =  NULL
  convergence =  TRUE
  # npar = Nseg*2
  height.plot = npar*2
#***************************************************************************************************************
  if (npar == 2){
     nroww = 1
     ncoll = 3
  } else {
    if (((npar+1) %% 5) ==0) {
      nroww = (npar+1)/5
      ncoll= 5
    } else { 
      nroww = ceil((npar+1)/5)
      ncoll = 5
    }
  }
  
  pdf(paste0(dir.plot,"/autocorrelation.pdf"), width = 10, height = 3*nroww, useDingbats=F)
  par(mar=c(1,1,1,1))
  par(mfrow = c(nroww, ncoll))   #c(nrows, ncol)
  for (i in 1:(ncol(mcmc.tot))) {
         par(mar = c(3, 3, 3, 3))  # Set the margin on all sides 
         autocorr.plot(x = mcmc.tot[,i],  auto.layout = FALSE)
         title(main = names(mcmc.tot)[i], cex.main = 2)
  }
  dev.off()
  dev.set(dev.prev())
  
  #*************************************************************************************************************
  file.name = paste0(dir.plot,"/convergence.csv")
  cat("Rc;", file = file.name, sep=" ")
  cat("Ru",  file = file.name, append = TRUE, sep="\n")
  #*************************************************************************************************************
  pdf(paste0(dir.plot,"/gelman.pdf"), width = 10, height = 3*nroww, useDingbats=F)
  par(mar=c(1,1,1,1))
  par(mfrow = c(nroww, ncoll))   #c(nrows, ncol)
  for (i in 1:(npar+1)) {
    mcmc.list  = list(mcmc.tot[1:(len/4),i], 
                      mcmc.tot[(len/4+1):(len/2),i], 
                      mcmc.tot[(len/2+1):(len*3/4),i], 
                      mcmc.tot[(len*3/4+1):len,i])
    mcmc.list  = as.mcmc.list(lapply(mcmc.list, mcmc))
    gel[[i]]   = gelman.diag(x = mcmc.list, confidence=0.95, transform=FALSE)
    par(mar    = c(3, 3, 3, 3)) # Set the margin on all sides 
    gelman.plot(x = mcmc.list,auto.layout = FALSE) 
    title(main = names(mcmc.tot)[i], cex.main = 2)
    
    #acfplot(x = mcmc.list, lag.max=200)
    if (gel[[i]]$psrf[1] > 1.2) { 
      convergence = FALSE
    }
    cat(gel[[i]]$psrf[1], file = file.name, append = TRUE,sep="")
    cat(";", file = file.name, append = TRUE,sep=" ")
    cat(gel[[i]]$psrf[2], file = file.name, append = TRUE,sep="\n")
  }
  dev.off()
  dev.set(dev.prev())
  #************************************************************************************************************
  Rc= "Reduction factor Rc"
  Ru="Upper confidence limit Ru" 
  limit= "threshold"
  col.gelman = c("Reduction factor Rc"       = "blue", 
                 "Upper confidence limit Ru" = "red", 
                 "threshold"                 = "red")
  Rc.Ru.plot <- ggplot()+
    coord_cartesian(xlim = c(1, length(gel)))+
    #scale_y_continuous(breaks=seq(0,Qmax,Qmax/5))+
    scale_x_continuous(breaks=seq(1, length(gel), 1)) +
    geom_point(mapping = aes(x=seq(1,length(gel),1), unlist(gel, use.names = FALSE)[seq(1,length(gel)*2,2)] , 
                             col=Rc),shape=0, size = 4, stroke = 0.5) +
    geom_point(mapping = aes(x=seq(1,length(gel),1), unlist(gel, use.names = FALSE)[seq(2,length(gel)*2,2)], 
                             col=Ru), shape=4, size = 4, stroke = 0.5) +
    theme_bw(base_size = 15)+
    geom_hline(aes(yintercept = 1.2, col= limit), lwd =1, linetype = "dashed") +
    scale_colour_manual(name = "Legend", values=col.gelman ,
                        breaks=c(Rc, Ru, limit),
                        guide = guide_legend(override.aes = list(
                          linetype = c("blank", "blank", "dashed"),
                          shape = c(0, 4, NA)), 
                          title.position="top")) +
    xlab("Parameter index")+ 
    ylab("Gelman factor")+
    theme(text                  = element_text(size=15)
          ,plot.background      = element_rect(fill ="transparent", color = NA)
          ,panel.grid.major     = element_line(size=0.4, linetype = "dashed")
          ,panel.grid.minor     = element_blank()
          ,legend.position      = "bottom"
          ,legend.text          = element_text(size = 10) 
          ,legend.justification = c(0,1)
          ,legend.box           = "vertical"
          ,legend.box.just      = "left"
          ,legend.text.align    = 0)
  pdf(paste0(dir.plot,"/Rc_Ru_gelman.pdf"), 8, 5 , useDingbats=F)
  print(Rc.Ru.plot)
  dev.off()
  return(convergence)
}



















#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                       BaM PLOTS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ReadModel<-function(ModelFile){
  #------------------------------------------------------------
  # Read information on the fitted model in Config_Model.txt
  #------------------------------------------------------------
  # read model info
  k=0;ID=qscan(ModelFile,k)[1]
  k=k+1;nX=as.numeric(qscan(ModelFile,k)[1])
  k=k+1;nY=as.numeric(qscan(ModelFile,k)[1])
  k=k+1;nPar=as.numeric(qscan(ModelFile,k)[1])
  # parameter blocs
  par=list();m=0
  for(i in 1:nPar){
    k=k+1;name=qscan(ModelFile,k)[1]
    k=k+1;init=as.numeric(qscan(ModelFile,k)[1])
    k=k+1;dist=qscan(ModelFile,k)[1]
    k=k+1;foo=qscan(ModelFile,k,sep='?')
    if(dist=='Gaussian' | dist=='Uniform' | dist=='FlatPrior' | dist=='FlatPrior+' | dist=='LogNormal' | dist=='FlatPrior-'){
      ppar=as.numeric(strsplit(gsub(pattern=",",replacement=" ",foo),split=" ")[[1]][1:2])
    } else {
      ppar=c(NA,NA)
    }        
    if(dist!='FIX'){
      m=m+1;par[[m]]=list(name=name, init=init, prior=list(dist=dist, par=ppar))
    }
  }
  return(list(ID=ID,nX=nX,nY=nY,nPar=nPar,par=par))
}













#################################################################################################
MCMCplot<-function(MCMCfile,
                   doPar, #=T or F
                   doLogPost, # =T or F
                   doDPar,  # =T or F
                   type,  #="trace", # "histogram","density","scatterplot"
                   xlab, # =''
                   ylab, #=''
                   ncol, # number of columns in the plot
                   prior, # = NULL if traceplot type
                   burn,  # = 0 , fraction of mcmc to burn
                   slim,   # = 1, mcmc to slim
                   theme_size){
#################################################################################################
  # Plot MCMC samples
  #------------------------------------------------------------
  # read MCMC simulation, burn and slim if required
  X = read.table(MCMCfile, header=T)
  X = X[seq(burn+1, nrow(X), slim), ]
  # determine logpot column
  lpcol=which(names(X)=='LogPost')
  # determine which columns should be kept
  keep=c()
  if(doPar){keep=c(keep,1:(lpcol-1))}
  if(doLogPost){keep=c(keep,lpcol)}
  if(doDPar){if(ncol(X)>lpcol){keep=c(keep,(lpcol+1):ncol(X))}}
  X = X[,keep]
  if(type=="scatterplot"){
    m=ggpairs(X,diag=list(continuous='density'),
              lower=list(continuous='points'),upper=list(continuous='cor'),axisLabels='show')
  } else {
    X=cbind(X, indx = 1:nrow(X))
    Xm=melt(X, id.vars = 'indx')
    m=ggplot(Xm)
    m=m+switch(type,
               trace     = geom_line(aes(x = indx, y = value, colour = variable)),
               histogram = geom_histogram(aes(value, fill = variable, ..density..)),
               density   = geom_density(aes(value, fill = variable), colour=NA, alpha = 0.8))
    m=m+facet_wrap(~variable, scales='free', ncol=ncol)+
        xlab(xlab)+
        ylab(ylab)+
        theme(legend.position="none")


    # Plot priors if required
    if(doPar & (type=='density' | type=='histogram') & !is.null(prior)){
      #priorDF = data.frame(matrix(NA, length(prior), 100))
      y= matrix(NA, 100, length(prior))
      grid = matrix(NA, 100, length(prior)) 
      for(i in 1:length(prior)){
        grid[,i]=seq(min(X[,i]), max(X[,i]),length.out=100)
        if(prior[[i]]$prior$dist=='Gaussian'){
          y[,i]=dnorm(grid[,i],mean=prior[[i]]$prior$par[1],sd=prior[[i]]$prior$par[2])
        }
        if(prior[[i]]$prior$dist=='LogNormal'){
          y[,i]=dlnorm(grid[,i], meanlog = prior[[i]]$prior$par[1], sdlog =prior[[i]]$prior$par[2])
        }
        if(prior[[i]]$prior$dist=='Uniform'){
          y[,i]=dunif(grid[,i],min=prior[[i]]$prior$par[1],max=prior[[i]]$prior$par[2])
        }
        if(prior[[i]]$prior$dist=='FlatPrior'){
          y[,i]=dunif(grid[,i],min=-99999999,max=99999999)
        }
        if(prior[[i]]$prior$dist=='FlatPrior+'){
          y[,i]=dunif(grid[,i],min=0,max=99999999)
        }
        if(prior[[i]]$prior$dist=='FlatPrior-'){
          y[,i]= dunif(grid[,i],min=-99999999,max=0)
        }
        if(prior[[i]]$prior$dist=='Gaussian'  | prior[[i]]$prior$dist=='Uniform'   |
           prior[[i]]$prior$dist=='FlatPrior' | prior[[i]]$prior$dist=='FlatPrior+' | prior[[i]]$prior$dist=='FlatPrior-' |
           prior[[i]]$prior$dist=='LogNormal'){
           # priorDF = rbind(priorDF, data.frame(variable = prior[[i]]$name,
           #                                     x        = grid,
           #                                     y        = y,
           #                                     ind      = i))
        }
      }
      
      priorDF = data.frame(y)
      priorDF.grid = data.frame(grid)
      names(priorDF) = unlist(lapply(prior, '[[', 1))
      names(priorDF.grid) = unlist(lapply(prior, '[[', 1))
      priorDF.grid.bis    = cbind(priorDF.grid, indx = 1:nrow(priorDF.grid))
      priorDF.grid.bis.m  = melt(priorDF.grid.bis, id.vars = 'indx')
      priorDF.bis    = cbind(priorDF, indx = 1:nrow(priorDF))
      priorDF.bis.m  = cbind(melt(priorDF.bis, id.vars = 'indx'), grid = priorDF.grid.bis.m$value)

      m=ggplot(Xm)
      m=m+switch(type,
                 trace     = geom_line(aes(x = indx, y = value, colour = variable)),
                 histogram = geom_histogram( aes(value, fill = variable, ..density..)),
                 density   = geom_density(aes(value, fill = variable), colour=NA, alpha = 0.8)) +
         geom_line(data = priorDF.bis.m, aes(x=grid, y=value), colour="black")
      m=m+
        facet_wrap(~variable, scales='free', ncol= ncol)+
        xlab(xlab)+
        ylab(ylab)+
        theme(legend.position="none")
      
    }
  }
  m= m +
    theme_bw(base_size      = theme_size)+
    theme(legend.position   = "none"
          ,axis.text        = element_text(size=10)
          ,plot.margin      = unit(c(1, 0.5, 0.5, 0.5),"cm")
          ,panel.grid.minor = element_blank())
  return(m)
}






# #####################################################################################
# MCMCplot<-function(MCMCfile,
#                    doPar,
#                    doLogPost,
#                    doDPar,
#                    type, # "histogram","density","scatterplot"
#                    xlab, ylab,
#                    ncol, prior,
#                    burn, slim,
#                    theme_size){
# #####################################################################################
#   # Plot MCMC samples
#   #-------------------------------------------------
#   # read MCMC simulation, burn and slim if required
#   X=read.table(MCMCfile,header=T)
#   X=X[seq(burn+1,nrow(X),slim),]
#   # determine logpot column
#   lpcol=which(names(X)=='LogPost')
#   # determine which columns should be kept
#   keep=c()
#   if(doPar){keep=c(keep,1:(lpcol-1))}
#   if(doLogPost){keep=c(keep,lpcol)}
#   if(doDPar){if(ncol(X)>lpcol){keep=c(keep,(lpcol+1):ncol(X))}}
#   X=X[,keep]
#   if(type=="scatterplot"){
#     m=ggpairs(X,diag=list(continuous='density'),
#               lower=list(continuous='points'),upper=list(continuous='cor'),axisLabels='show')
#   } else {
#     X=cbind(X,indx=1:nrow(X))
#     Xm=melt(X,id.vars='indx')  
#     m=ggplot(Xm)
#     m=m+switch(type,
#                trace=geom_line(aes(x=indx,y=value,colour=variable)),
#                histogram=geom_histogram(aes(value,fill=variable,..density..)),
#                density=geom_density(aes(value,fill=variable),colour=NA))
#     m=m+facet_wrap(~variable, scales='free',ncol=ncol)+
#       xlab(xlab)+ylab(ylab)+
#       theme(legend.position="none")
#     
#     # Plot priors if required
#     if(doPar & (type=='density' | type=='histogram') & !is.null(prior)){ 
#       priorDF=data.frame()
#       for(i in 1:length(prior)){
#         grid=seq(min(X[,i]),max(X[,i]),length.out=100)
#         if(prior[[i]]$prior$dist=='Gaussian'){
#           y=dnorm(grid,mean=prior[[i]]$prior$par[1],sd=prior[[i]]$prior$par[2])
#         }
#         if(prior[[i]]$prior$dist=='Uniform'){
#           y=dunif(grid,min=prior[[i]]$prior$par[1],max=prior[[i]]$prior$par[2])
#         }
#         if(prior[[i]]$prior$dist=='LogNormal'){
#           y=dlnorm(grid,meanlog=prior[[i]]$prior$par[1],sdlog=prior[[i]]$prior$par[2])
#         }
#         if(prior[[i]]$prior$dist=='Gaussian' | prior[[i]]$prior$dist=='Uniform' | prior[[i]]$prior$dist=='LogNormal'){
#           priorDF=rbind(priorDF,data.frame(variable=prior[[i]]$name,x=grid,y=y))
#         }
#       }
#       m=m+geom_line(aes(x=x,y=y),
#                     data=priorDF,
#                     colour="black")+
#         theme_bw(base_size      = theme_size)
#     }   
#   }
#   return(m)
# }











##########################################################################################################
Residualplot<-function(ResFile='Results_Residuals.txt',
                       nX,nY,
                       xlab=c('Observed','Simulated','Simulated','Observed X'),
                       ylab=c('Simulated','Residuals','Standardized Residuals','Observed / simulated Y')){
##########################################################################################################
  # Perform residual analysis
  #------------------------------------------------------------
  R=read.table(ResFile,header=T)
  #------------------------------------------------------------
  # Yobs/Yunbiased vs. Ysim and Ysim vs. residuals
  DF=data.frame()
  for(i in 1:nY){
    DF=rbind(DF,
             data.frame(variable=i,
                        Yobs=R[[paste('Y',i,'_obs',sep='')]],
                        Yunbiased=R[[paste('Y',i,'_unbiased',sep='')]],
                        Ysim=R[[paste('Y',i,'_sim',sep='')]],
                        Residual=R[[paste('Y',i,'_res',sep='')]],
                        StandardizedResidual=R[[paste('Y',i,'_stdres',sep='')]]))
  }
  # Obs vs. sim plot
  m1=ggplot(DF)+
    geom_point(aes(x=Yobs,y=Ysim),colour='black',size=5)+
    facet_wrap(~variable,scales='free')+
    theme(legend.position="none")
  if( any(DF$Yobs!=DF$Yunbiased) ){
    m1=m1+geom_point(aes(x=Yunbiased,y=Ysim),size=3,colour='red')
    xl='Observed (black) / Unbiased (red)'
  } else {xl=xlab[1]}
  m1=m1+xlab(xl)+ylab(ylab[1])+geom_abline(intercept=0,slope=1)
  # sim vs. res plot
  m2=ggplot(DF)+
    geom_point(aes(x=Ysim,y=Residual),colour='black',size=5)+
    facet_wrap(~variable,scales='free')+
    theme(legend.position="none")+geom_abline(intercept=0,slope=0)+
    xlab(xlab[2])+ylab(ylab[2])
  # sim vs. standardized residuals plot
  m3=ggplot(DF)+
    geom_point(aes(x=Ysim,y=StandardizedResidual),colour='black',size=5)+
    facet_wrap(~variable,scales='free')+
    theme(legend.position="none")+geom_abline(intercept=0,slope=0)+
    xlab(xlab[3])+ylab(ylab[3])
  #------------------------------------------------------------
  # Xobs vs. Yobs/Ysim
  DF=data.frame()
  for(i in 1:nX){
    for(j in 1:nY){
      DF=rbind(DF,
               data.frame(Xvariable=i,Yvariable=j,
                          Xobs=R[[paste('X',i,'_obs',sep='')]],
                          Yobs=R[[paste('Y',j,'_obs',sep='')]],
                          Ysim=R[[paste('Y',j,'_sim',sep='')]]))
    }
  }
  m4=ggplot(DF)+
    geom_point(aes(x=Xobs,y=Yobs),colour='black',size=5)+
    geom_line(aes(x=Xobs,y=Ysim),colour='black')+
    facet_wrap(Yvariable~Xvariable,scales='free',ncol=nX)+
    theme(legend.position="none")+
    xlab(xlab[4])+ylab(ylab[4])
  #------------------------------------------------------------
  return(list(m1=m1,m2=m2,m3=m3,m4=m4))
  
}











#################################################################################################             
EnvelopLayer<-function(file, Xfile=NULL, color, alpha, sep=F){
################################################################################################# 
  # Return the ggplot layer from an envelop file 
  # WARNING: just a layer! need to call ggplot() before printing 
  #------------------------------------------------------------
  E=read.table(file,header=T)
  n=nrow(E)
  if(is.null(Xfile)){X=matrix(1:n,n,1)} else {X=read.table(Xfile,header=F)}
  w=sort(X[,1],index.return=T)
  if(!sep){ # uncertainty bands
    x=c(X[w$ix[1:n],1],X[w$ix[n:1],1])
    y=c(E[w$ix[1:n],2],E[w$ix[n:1],3])
    DF=data.frame(x=x,y=y)
    m=geom_polygon(aes(x=x,y=y,fill=color),data=DF,alpha=alpha, colour=NA)
  } else { # separated uncertainty bars
    x=X[w$ix[1:n],1]
    lo=E[w$ix[1:n],2]
    hi=E[w$ix[1:n],3]
    DF=data.frame(x=x,lo=lo,hi=hi)
    m=geom_crossbar(aes(x=x,y=lo,ymin=lo,ymax=hi),data=DF,fill=color,alpha=alpha,colour=NA)
  }
  return(m)
}










#################################################################################################
EnvelopLayerLog<-function(file,Xfile=NULL, color, alpha, sep=F){
#################################################################################################
  #------------------------------------------------------------
  # Return the ggplot layer from an envelop file 
  # WARNING: just a layer! need to call ggplot() before printing 
  #------------------------------------------------------------
  E=read.table(file,header=T)
  n=nrow(E)
  if(is.null(Xfile)){X=matrix(1:n,n,1)} else {X=read.table(Xfile,header=F)}
  w=sort(X[,1],index.return=T)
  if(!sep){ # uncertainty bands
    x=c(X[w$ix[1:n],1],X[w$ix[n:1],1])
    y=c(E[w$ix[1:n],2],E[w$ix[n:1],3])
    DF=data.frame(x=x,y=y)
    m=geom_polygon(aes(x=x,y=y),data=DF,fill=color,alpha=alpha,colour=NA,na.rm=T)
  } else { # separated uncertainty bars
    x=X[w$ix[1:n],1]
    lo=E[w$ix[1:n],2]
    hi=E[w$ix[1:n],3]
    DF=data.frame(x=x,lo=lo,hi=hi)
    m=geom_crossbar(aes(x=x,y=lo,ymin=lo,ymax=hi),data=DF,fill=color,alpha=alpha,colour=NA,na.rm=T)
  }
  return(m)
}







######################################################################
SpaghettiLayer<-function(file,Xfile=NULL,color='black',size=1,sep=F) {
######################################################################
  #------------------------------------------------------------
  # Return the ggplot layer from a spaghetti file 
  # WARNING: just a layer! need to call ggplot() before printing 
  #------------------------------------------------------------
  E=read.table(file,header=F)
  n=nrow(E);p=ncol(E)
  if(is.null(Xfile)){X=matrix(1:n,n,1)} else {X=read.table(Xfile,header=F)}
  for(i in 1:p){
    DF=data.frame(x=X[,1],y=E[,i])
    if(!sep) { # line
      mi=geom_line(aes(x=x,y=y),data=DF,colour=color,size=size)
    } else { # separated points
      mi=geom_point(aes(x=x,y=y),data=DF,colour=color,size=size)
    }
    if(i==1){m=mi} else {m=m+mi}
  }
  return(m)
}







#################################################################################################
plot.Qt <- function(t,Qmed,Q2,Q97) {
#################################################################################################
  df <- data.frame(t,Qmed,Q2,Q97)  
  ggplot() +
    #geom_point(data = df, aes(x = h_Gaug, y = log(Q_Gaug)), size = 2) +
    geom_line(data = df, aes(x = t, y = Qmed), color = "blue",size = 0.8) +
    #geom_line(data = CdT, aes(x = h_Gaug, y = log(Q_CM_97)), color = "red",size = 1) +
    #geom_line(data = CdT, aes(x = h_Gaug, y = log(Q_CM_2)), color = "red",size = 1) #+
    #geom_point(data = Qcrit, aes(x =c(k1,k2,k3) , y = c(Q1,Q2,Q3)), col="red",size = 3)+
    geom_ribbon(aes(x = t, ymin = Q2, ymax = Q97), alpha=0.5,fill = "blue") +
    xlim(0,200)+
    ylim(0,40)
}












#################################################################################################
plot.RC <- function(h,Qmed,Q2,Q97,h_Gaug,Q_Gaug,uQ_Gaug) {
#################################################################################################
  df <- data.frame(h,Qmed,Q2,Q97)
  dg <- data.frame(h_Gaug,Q_Gaug,uQ_Gaug)
  ggplot() +
    #geom_point(data = df, aes(x = h_Gaug, y = log(Q_Gaug)), size = 2) +
    geom_line(data = df, aes(x = h, y = Qmed), color = "blue",size = 0.8) +
    geom_point(data = dg, aes(x=h_Gaug, y=Q_Gaug)) +
    #geom_line(data = CdT, aes(x = h_Gaug, y = log(Q_CM_97)), color = "red",size = 1) +
    #geom_line(data = CdT, aes(x = h_Gaug, y = log(Q_CM_2)), color = "red",size = 1) #+
    #geom_point(data = Qcrit, aes(x =c(k1,k2,k3) , y = c(Q1,Q2,Q3)), col="red",size = 3)+
    geom_errorbar(data = dg, aes(x=h_Gaug, ymin=Q_Gaug-2*uQ_Gaug, ymax =Q_Gaug+2*uQ_Gaug), 
                  size = 0.7, width=15, col= "black") +
    geom_ribbon(data = df,aes(x = h, ymin = Q2, ymax = Q97), alpha=0.5,fill = "blue") +
    xlim(-1,2.5)+
    ylim(0,250)+
    theme_bw()
}














#################################################################################################
plot.log.RC <- function(h,Qmed,Q2,Q97,h_Gaug,Q_Gaug,uQ_Gaug) {
#################################################################################################
  df <- data.frame(h,Qmed,Q2,Q97)
  dg <- data.frame(h_Gaug,Q_Gaug,uQ_Gaug)
  
  ggplot() +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    #geom_point(data = df, aes(x = h_Gaug, y = log(Q_Gaug)), size = 2) +
    geom_line(data = df, aes(x = h, y = Qmed), color = "blue", size = 0.8) +
    geom_point(data = dg, aes(x=h_Gaug, y=Q_Gaug)) +
    geom_ribbon(data = df,aes(x = h, ymin = Q2, ymax =Q97),alpha=0.5,fill = "blue") +
    geom_errorbar(data= dg, aes(x=h_Gaug, ymin=Q_Gaug-2*uQ_Gaug, ymax =Q_Gaug+2*uQ_Gaug), 
                  size = 0.7, width=15, col= "black") +
    xlim(-1,2.5)+
    #ylim(0.1, 100)+
    theme_bw()
}



















#################################################################################################
qscan<-function(f,k,sep=' '){
#################################################################################################
  # quick scan - just a convenience wrapper
  return(scan(f,what='character',skip=k,nlines=1,quiet=TRUE,sep=sep))
}












#################################################################################################
Calendar2Days<-function(fileIN,
                        fileOUT,
                        format="%d/%m/%Y %H:%M:%S"){
  # read input file as text
  w=readLines(fileIN)
  # convert to date format
  z=as.POSIXct(w,format=format)
  # convert to seconds
  ww=unclass(z)
  # convert to days
  d=(ww-ww[1])/(60*60*24)
  # write result
  write(d,file=fileOUT,ncolumns=1)
}















##################################################################################################
correlation.plot = function(MCMCfile,
                            burn=0,
                            slim=1, 
                            npar) {
#################################################################################################
  #------------------------------------------------------------
  # Plot MCMC samples
  #------------------------------------------------------------
  # read MCMC simulation, burn and slim if required
  X=read.table(MCMCfile, header=T)
  X=X[seq(burn+1,nrow(X),slim), 1:npar]
  plt =  pairs.panels(X,
                      method = "pearson", # correlation method
                      hist.col = "#00AFBB",
                      density = TRUE,  # show density plots
                      ellipses = TRUE # show correlation ellipses
  )
  
  return(plt)
}





