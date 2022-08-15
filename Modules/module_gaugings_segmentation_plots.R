#################################################################################################################################
model.selection.plot <- function(criteria.df, 
                                 dir.seg.gaug, 
                                 seg.iter) {
#################################################################################################################################
  col =c("AIC" ="blue", 
         "BIC" ="gray40",
         "HQC" ="red",
         #"AICc (Hurvich and Tsai, 1989)" ="orange", 
         "DIC" ="green")
  AICcol ="AIC"; BICcol ="BIC"; HQCcol= "HQC";    DICcol= "DIC";  # AICccol= "AICc"; 

  
  bic.plot <- ggplot(criteria.df) + 
      theme_light(base_size = 15)+
      geom_line(aes(x = x, y = BIC,  color = BICcol), size = 0.1, na.rm = TRUE) +
      geom_line(aes(x = x, y = AIC,  color = AICcol), size = 0.1, na.rm = TRUE) +
      geom_line(aes(x = x, y = HQC,  color = HQCcol), size = 0.1, na.rm = TRUE) +
      #geom_line(aes(x = x, y = AICc, color =AICccol), size = 0.5) +
      geom_line(aes(x = x, y = DIC,  color =DICcol), size = 0.1, na.rm = TRUE) +
      scale_x_continuous(name = "Number of segments", expand = c(0,0), limits =c(1, length(criteria.df$BIC)))+
      #scale_y_continuous(name = "BIC, AIC, HQC, DIC" , expand = c(0,0),limits = c(min_grid, max_grid), breaks = seq(min_grid, max_grid, 50)) +
                         # breaks = scales::pretty_breaks(n = 2)) + 
      scale_y_continuous(name = "Model selection criteria" , expand = c(0,0)) +
      xlab("Number of segments") + 
      ylab("Model selection criteria") +
    
      geom_point(aes(x = x, y = BIC, color = BICcol), shape =0, size = 4, na.rm = TRUE) +
      geom_point(aes(x = x, y = AIC, color = AICcol), shape =2, size = 4, na.rm = TRUE) +
      geom_point(aes(x = x, y = HQC, color = HQCcol), shape =1, size = 4, na.rm = TRUE) +
      geom_point(aes(x = x, y = DIC, color = DICcol), shape =5, size = 4, na.rm = TRUE) +
      # geom_point(aes(x = x, y = AICc)), shape =6, size = 4, color = "orange", na.rm = TRUE) +  
      #
      geom_point(aes(x = BICmin, y = BIC[BICmin]), shape = 15, size = 4, color = "gray20", na.rm = TRUE) +
      geom_point(aes(x = AICmin, y = AIC[AICmin]), shape = 17, size = 4, color = "blue", na.rm = TRUE) +
      geom_point(aes(x = HQCmin, y = HQC[HQCmin]), shape = 19, size = 4, color = "red", na.rm = TRUE)+
      geom_point(aes(x = DICmin, y = DIC[DICmin]), shape = 23, size = 4, color = "green", fill="green", na.rm = TRUE)+
      # geom_point(aes(x = AICcmin, y = AICc[AICcmin]), shape = 25, size = 4, color = "orange", fill = "orange", na.rm = TRUE)+
      coord_cartesian(clip = 'off')+
      scale_colour_manual(name   = element_blank(),
                          values = col, 
                          breaks = c(AICcol, 
                                     BICcol, 
                                     HQCcol,
                                     DICcol),
                          guide = guide_legend(override.aes = list(
                          linetype = c("solid", "solid","solid", "solid"),
                          shape = c(2, 0, 1, 5)),
                          size  = c(4,4,4,4)) ) +
      theme(  plot.title        = element_text(hjust = 0.5) 
            ,plot.background   = element_rect(fill ="transparent", color = NA)
            ,panel.grid.major  = element_line(size=0.4, linetype = "dashed")
            ,panel.grid.minor  = element_blank()
            #,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            ,panel.background  = element_rect(fill ="transparent") 
            #,axis.line        = element_line(colour = "black")
            ,axis.ticks        = element_line(colour = "black")
            ,plot.margin       = unit(c(0.3,0.5,0.05,1.05),"cm")
            ,text              = element_text(size=10)
            ,legend.key        = element_rect(colour = "transparent", fill = "transparent")
            ,legend.background = element_rect(colour = "transparent", fill = "transparent")
            #,      axis.line = element_line(colour = "black"),
            #,legend.title = element_text(colour="blue", size=16, face="bold"),
            ,legend.position="bottom") 
    
  
  pdf(paste0(dir.seg.gaug,"/Criteria_it",seg.iter,".pdf"), 6, 4 , useDingbats=F)
  print(bic.plot)
  dev.off()
  
  return(bic.plot)
}


















#*************************************************************************************************************
segm_P_plot <- function(dir.seg.gaug, 
                        Shift.Q, 
                        Q10.ts, Q90.ts, ts.res, ts.res.before, ts.res.plus, 
                        Q10.mu.res, Q90.mu.res, mu.res, 
                        seg.iter, 
                        df.limni, 
                        ts.real,
                        resid.uncertaint) {
#*************************************************************************************************************
# Description:
# This function returns a plot with the results of the multi-changepoint segmentation of the residuals between 
# gaugings and rating curve for the iteration P:
#
  Utot.times   = "90% total uncertainty of change point times"
  MAPtimes     = "Change point times (MAP)"
  realtimes    = "Real rating shift time"
  Utot.mean    = "90% total uncertainty of segment mean"
  MAPmean      = "Segment mean (MAP)"
  resid        = "Gaugings residuals with 95% error bars"
  color.segm   = c("90% total uncertainty of change point times" = "blue",
                   "90% total uncertainty of segment mean"       = "red")
  col.segm     = c("Real rating shift time"                 = "red", 
                   "Change point times (MAP)"               = "blue", 
                   "Segment mean (MAP)"                     = "red", 
                   "Gaugings residuals with 95% error bars" = "black")
  color.segm_1 = c("90% total uncertainty of segment mean"  = "red")
  col.segm_1   = c("Segment mean (MAP)"                     = "red",
                   "Gaugings residuals with 95% error bars" = "blue")
  
  if (!is.null(ts.res[1])) {
    seg.plt <- ggplot()+
               geom_point(data = Shift.Q, aes(x = alpha_t, y = alpha, col=resid), 
               #color = Shift.Q$Q_Gaug), 
               shape = 1, size = 1.4)
    
    if (resid.uncertaint==TRUE) {
      seg.plt <-  seg.plt+
        scale_y_continuous(name   = bquote(.("Residual ") ~ .("[") ~ m^3*s^-1 ~ .("]")), # $r = Q_i - Q_{RC} (h_i) $"),
                           limits = c(min(Shift.Q$alpha - 2*Shift.Q$sigma.tot),
                                      max(Shift.Q$alpha + 2*Shift.Q$sigma.tot))) + 
        geom_errorbar(data = Shift.Q, 
                      aes(x=alpha_t, 
                          ymin= (alpha - 2*sigma.tot), 
                          ymax = (alpha + 2*sigma.tot), color = resid),
                      size = 0.1, width=0.01*(tail(Shift.Q$alpha_t,1)- Shift.Q$alpha_t[1]))
    } else {
      seg.plt <-  seg.plt+
        scale_y_continuous(name = bquote(.("Residual ") ~ .("[") ~ m^3*s^-1 ~ .("]"))) #$r = Q_i - Q_{RC} (h_i) $"))
    }
    seg.plt <-  seg.plt +
    coord_cartesian(clip = 'off') + 
    geom_rect( mapping= aes(xmin = Q10.ts,
                            xmax = Q90.ts,
                            ymin = -Inf, 
                            ymax = Inf, fill=Utot.times), 
               alpha=0.1 ) +
    # annotate("rect",xmin= Q10.ts, xmax=Q90.ts,
    #          ymin=-Inf, ymax=Inf, fill="blue", alpha=0.1) +
    geom_vline(aes(xintercept = ts.res, col=MAPtimes), lwd =0.5, linetype = "solid") +
    geom_vline(aes(xintercept = ts.real, col=realtimes), lwd =0.3, linetype ="longdash") +
    geom_rect(mapping = aes(xmin = ts.res.before, 
                            xmax = ts.res.plus, 
                            ymin = Q10.mu.res,
                            ymax = Q90.mu.res,
                            fill = Utot.mean), alpha=0.3) + 
    geom_segment(mapping=aes(x    = ts.res.before, 
                             y    = mu.res, 
                             xend = ts.res.plus, 
                             yend = mu.res, 
                             col  = MAPmean))+
    # annotate("rect",xmin= ts.res.before, xmax=ts.res.plus, ymin=(mu.res),
    #          ymax=(mu.res+max(Shift.Q$alpha + Shift.Q$sigma.tot)/100), fill="red", alpha=1) +
    # annotate("rect",xmin= ts.res.before, xmax=ts.res.plus, ymin=Q10.mu.res,
    #          ymax=Q90.mu.res, fill="red", alpha=0.3) +    
    scale_fill_manual(name     = element_blank(), 
                      values   = color.segm) +
    scale_colour_manual(name   = element_blank(), 
                        values = col.segm,
                        breaks = c(MAPtimes, MAPmean, resid),
                        labels = c(MAPtimes, MAPmean, resid),
                        guide  = guide_legend(override.aes = list(
                                 linetype = c("solid","longdash","solid", "blank"),
                                 shape = c(NA,NA, NA, 1)))) +
    theme_light(base_size = 15) +
    theme(text = element_text(size=10),
          plot.title = element_text(hjust = 0.5),
          #panel.grid.major=element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"),
          plot.margin=unit(c(0.3,0.5,0.05,1.05),"cm"),
          # theme(plot.title = element_text(hjust = 0.5),
          #       plot.background = element_rect(fill ="transparent", color = NA),
          #       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          #       panel.background = element_rect(fill ="transparent"), 
          #       axis.line = element_line(colour = "black"),
          #       axis.ticks = element_line(colour = "black"),
          #       plot.margin=unit(c(0.2,0.5,0.05,0.05),"cm"),
          legend.position ="none",
          #legend.position ="bottom", legend.box = "vertical", 
          # legend.key.size = unit(0.5, "cm"),
          # legend.text=element_text(size=8),
          # legend.spacing.x = unit(0.05, 'cm'),
          # legend.spacing.y = unit(0.05, 'cm'),
          # legend.key.height=unit(0.5,"line"),
          # legend.key.width=unit(0.8,"line"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.background = element_rect(colour = "transparent", fill = "transparent"))
    if (is.null(df.limni)==FALSE) {
         # seg.plt = seg.plt + scale_x_continuous(name = "time [days]", expand = c(0,0), 
         #                     limits =c(0,tail(df.limni$t_limni,1)))
          seg.plt = seg.plt + scale_x_continuous(name = "Time [days]", expand = c(0,0))
         
    } else {
         # seg.plt = seg.plt + scale_x_continuous(name = "time [days]", expand = c(0,0), 
         #                                     limits =c(0,tail(Shift.Q$alpha_t,1)))
         seg.plt = seg.plt + scale_x_continuous(name = "Time [days]", expand = c(0,0))
    }    
    #ggsave(seg.plt, filename =paste0(dir.seg.gaug,"/segm_it",seg.iter,".png"),width=4,height =3,dpi=300)
    pdf(paste0(dir.seg.gaug,"/segm_it",seg.iter,".pdf"), 5, 4 , useDingbats=F)
    print(seg.plt)
    dev.off()

    
    
    
    #****************
  } else {
    #****************
    seg.plt <- ggplot()+
    geom_point(data = Shift.Q, aes(x = alpha_t, y = alpha, col=resid), shape = 1, size = 1.4)
    if (resid.uncertaint==TRUE) {
          seg.plt <-  seg.plt +
          geom_errorbar(data = Shift.Q, aes(x=alpha_t, 
                                            ymin= (alpha - 2*sigma.tot), 
                                            ymax = (alpha + 2*sigma.tot), color = resid),
                        size = 0.1, width=0.01*(tail(Shift.Q$alpha_t,1)- Shift.Q$alpha_t[1]))
    }
    if (resid.uncertaint==TRUE) {
        seg.plt <-  seg.plt+
        scale_y_continuous(name   = bquote(.("Residual ") ~ .("[") ~ m^3*s^-1 ~ .("]")),    #$r = Q_i - Q_{RC} (h_i) $"),
                           limits = c(min(Shift.Q$alpha - 2*Shift.Q$sigma.tot),
                                      max(Shift.Q$alpha + 2*Shift.Q$sigma.tot))) + 
        geom_errorbar(data = Shift.Q, 
                      aes(x=alpha_t, 
                          ymin = (alpha - 2*sigma.tot), 
                          ymax = (alpha + 2*sigma.tot), color = resid),
                      size = 0.1, width=0.01*(tail(Shift.Q$alpha_t,1)- Shift.Q$alpha_t[1]))
    } else {
        seg.plt <-  seg.plt +
        scale_y_continuous(name = bquote(.("Residual ") ~ .("[") ~ m^3*s^-1 ~ .("]"))) # $r = Q_i - Q_{RC} (h_i) $"))
    }
    seg.plt <-  seg.plt +
    coord_cartesian(clip = 'off') + 
    geom_rect(data = Shift.Q, mapping = aes(xmin = alpha_t[1], 
                                            xmax = tail(alpha_t,1), 
                                            ymin = Q10.mu.res,
                                            ymax = Q90.mu.res, fill=Utot.mean), alpha=0.3) + 
    geom_segment(data = Shift.Q, mapping=aes(x =alpha_t[1], 
                                             y =mu.res, 
                                             xend = tail(alpha_t,1), 
                                             yend = mu.res, col =MAPmean)) +
    scale_fill_manual(name = element_blank(), values = color.segm_1) +
    scale_colour_manual(name     = element_blank(), 
                        values   = col.segm_1,
                        breaks   = c(MAPmean, resid),
                        labels   = c(MAPmean, resid),
                        guide    = guide_legend(override.aes = list(
                        linetype = c("solid", "blank"),
                        shape    = c(NA, 1)))) +
    theme_light(base_size = 15)+
    theme(text = element_text(size=14),
            #plot.title = element_text(hjust = 0.5),
            #panel.grid.major=element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent"),
            plot.margin=unit(c(0.3,0.5,0.05,1.05),"cm"),
            legend.position ="none",
            legend.key = element_rect(colour = "transparent", fill = "transparent"),
            legend.background = element_rect(colour = "transparent", fill = "transparent"))
    
    if (is.null(df.limni)==FALSE) {
      # seg.plt = seg.plt + scale_x_continuous(name = "time [days]", expand = c(0,0), 
      #                                        limits =c(0,tail(df.limni$t_limni,1)))
      seg.plt = seg.plt + scale_x_continuous(name = "Time [days]", expand = c(0,0))
    } else {
      # seg.plt = seg.plt + scale_x_continuous(name = "time [days]", expand = c(0,0), 
      #                                        limits =c(0,tail(Shift.Q$alpha_t,1)))
      seg.plt = seg.plt + scale_x_continuous(name = "Time [days]", expand = c(0,0))
    }    
    pdf(paste0(dir.seg.gaug,"/segm_it",seg.iter,".pdf"), 5, 4 , useDingbats=F)
    print(seg.plt)
    dev.off()
  } 
  
  
  return(seg.plt)
}

















#GaugingsShifts_P
#**************************************************************************************************
GaugingsShifts_P_plot <- function(dir.seg.gaug, 
                                  grid_RC.ylim, 
                                  grid_RC.xlim, 
                                  grid_RC.xstep, 
                                  grid_RC.ystep,
                                  ticks_RC.y.log,
                                  RC.x.labels, RC.y.labels, 
                                  gaugings, 
                                  df.RC, 
                                  seg.iter) {
#**************************************************************************************************
  # Linear plot:
  #****************
  gshifts <- ggplot(data = df.RC) + 
             theme_light(base_size = 10)+
             coord_cartesian(xlim = grid_RC.xlim,
                             ylim = grid_RC.ylim,
                             expand = c(0)) +
             scale_y_continuous(breaks=seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep))+
             scale_x_continuous(breaks=seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep)) + 
             geom_point(data = gaugings, aes(x = h, y = Q), size = 1.7, pch =21, color ="gray90", fill = "gray90") +
             geom_errorbar(data = df.RC, 
                  aes(x = hP, ymin =QP-2*uQP, ymax = QP +2*uQP), 
                  width=0.02*(grid_RC.xlim[2] - grid_RC.xlim[1]), size = 0.3,
                  color = "black") +
             geom_point(data = df.RC, aes(x = hP, y = QP), size = 1.7, pch =21, fill = df.RC$color) +
  
               # theme(plot.title = element_text(hjust = 0.5),
               #       plot.background = element_rect(fill ="transparent", color = NA),
               #       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               #       panel.background = element_rect(fill ="transparent"), 
               #       axis.line = element_line(colour = "black"),
               #       axis.ticks = element_line(colour = "black"),
             theme(text               = element_text(size=10)
                  ,plot.title         = element_text(hjust = 0.5)
                  ,panel.grid.major   = element_line(size=0.4, linetype = "dashed") 
                  ,panel.grid.minor   = element_blank()
                  ,plot.background    = element_rect(fill = "transparent", color = NA)
                  ,panel.background   = element_rect(fill = "transparent")
                  ,plot.margin        = unit(c(0.3,0.5,0.05,1.05),"cm")
                  ,axis.ticks         = element_line(colour = "black")) +
             ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
             xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]")))+
             coord_cartesian(clip = 'off')
   
             #ggsave(gshifts, filename =paste0(dir.seg.gaug,"/GaugingsShifts_it",seg.iter,".png"), width=6,height=5,dpi=300)
             pdf(paste0(dir.seg.gaug,"/GaugingsShifts_it",seg.iter,".pdf"), 5, 4 , useDingbats=F)
             print(gshifts)
             dev.off()
             
  
  
  #log plot:
  #********************************************************************************
  glogshifts <- ggplot(data = df.RC) + 
                coord_trans(xlim = grid_RC.xlim, ylim = grid_RC.ylim.log) +
                scale_y_log10(na.value = -10, 
                              breaks   = ticks_RC.y.log, 
                              labels   = ticks_RC.y.log) +
                scale_x_continuous(breaks = seq(grid_RC.xlim[1], 
                                                grid_RC.xlim[2], 
                                                grid_RC.xstep)) +
                geom_point(data = gaugings, aes(x = h, y = Q), size = 1.7, pch =21, color ="gray90", fill = "gray90") +
                geom_errorbar(data = df.RC, aes(x=hP, ymin =QP-2*uQP, ymax = QP+2*uQP), 
                              color = "black", width=0.02*(grid_RC.xlim[2] - grid_RC.xlim[1]), size = 0.3) +
                geom_point(data = df.RC, aes(x = hP, y = QP), size = 1.7, pch =21, fill = df.RC$color) +

                #scale_colour_manual(name = element_blank(), 
                #                    values=c("Gaugings of current period" = "black",
                #                             "Gaugings of other periods" ="gray",
                #                             "MAP rating curve" = "blue" ),
                #                    breaks=c("Gaugings of current period",
                #                             "Gaugings of other periods",
                #                             "MAP rating curve"),
                #                    guide = guide_legend(override.aes = list(
                #                            linetype = c("blank", "blank","solid"),
                #                            shape = c(19, 19, NA)))) +
                ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
                xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]"))) + 
                annotation_logticks(base = 10, sides = "l", scaled = TRUE,  colour = "black", size = 0.3, linetype = 1) +
                coord_cartesian(clip = 'off')+
                theme_light(base_size = 10)+
                theme(plot.title          = element_text(hjust = 0.5) 
                      ,plot.background    = element_rect(fill ="transparent", color = NA)
                      ,panel.grid.major   = element_line(size=0.4, linetype = "dashed")
                      ,panel.grid.minor   = element_blank()
                      #,panel.grid.major  = element_blank()
                      #,panel.grid.minor  = element_blank()
                      ,panel.background   = element_rect(fill ="transparent") 
                      #,axis.line         = element_line(colour = "black")
                      ,axis.ticks         = element_line(colour = "black")
                      ,plot.margin        = unit(c(0.3,0.5,0.05,1.05),"cm") 
                      ,legend.position    = "none")
                      #,legend.position   = c(0.7, 0.2)
                      #,legend.box        = "vertical"#
                      #,legend.key.size   = unit(0.5, "cm"),
                      #,legend.position   = "bottom"
                      #,legend.direction  = "vertical",
                      #,legend.text       = element_text(size=12),
                      #,legend.spacing.x  = unit(0.05, 'cm'),
                      #,legend.spacing.y  = unit(0.01, 'cm'),
                      #,legend.key.height = unit(0.5,"line"),
                      #,legend.key.width  = unit(0.8,"line"),
                      #,legend.key        = element_rect(colour = "transparent", fill = "transparent"),
                      #,legend.background = element_rect(colour = "transparent", fill = "transparent"))
                # ggsave(glogshifts, filename=paste0(dir.seg.gaug,"/GaugingsShiftsLog_it",seg.iter,".png"),
                #       width=6,height=5,dpi=300)
  
  pdf(paste0(dir.seg.gaug,"/GaugingsShiftsLog_it",seg.iter,".pdf"), 5, 4 , useDingbats=F)
  print(glogshifts)
  dev.off()
  
  
  
  
  
  
 
  # #log-log plot:
  # #***********************************************************************************************************
  # gloglogshifts <- ggplot(data = df.RC) + 
  #   coord_trans(limx = c(Qmin_grid.log, Qmax_grid.log), limy = c(Qmin_grid.log, Qmax_grid.log) , clip = "off")+
  #   scale_y_log10(na.value=-10, breaks=ticks.log, labels=ticks.log) +
  #   scale_x_log10(na.value=-10, breaks=ticks.log, labels=ticks.log) + 
  #   geom_point(data = gaugings, aes(x = h, y =Q), size = 1.2, color = "grey") +
  #   geom_point(aes(x = hP, y = QP), color = color, size = 1.2) +
  #   geom_errorbar(aes(x=hP, ymin =QP-2*uQP, ymax =QP+2*uQP), color = color, 
  #                 width=.01, size =0.1)+
  #   xlab(gaugings.labels[1]) + ylab(gaugings.labels[2]) +
  #   annotation_logticks(base = 10, sides = "l", scaled = TRUE,
  #                       colour = "black", size = 0.3, linetype = 1) +
  #   theme_light(base_size = 10)+
  #   theme(plot.title = element_text(hjust = 0.5), 
  #         plot.background = element_rect(fill ="transparent", color = NA),
  #         panel.grid.major=element_line(size=0.4, linetype = "dashed"),
  #         panel.grid.minor=element_blank(),
  #         #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         panel.background = element_rect(fill ="transparent"), 
  #         #axis.line = element_line(colour = "black"),
  #         axis.ticks = element_line(colour = "black"),
  #         plot.margin=unit(c(0.3,0.5,0.05,1.05),"cm"), 
  #         #legend.position =c(0.7, 0.2), legend.box = "vertical", legend.key.size = unit(0.5, "cm"),
  #         legend.position = "none")
  # 
  # ggsave(gloglogshifts, filename =paste0(dir.seg.gaug,"/GaugingsShiftsLogLog_it",seg.iter,".png"), width=6,height=5,dpi=300)
  # 
  #************************************************
  return(list(gshifts, glogshifts)) #gloglogshifts
}



















############################################################################################
# Stage h time series  with  shift times (before the times adjustment):
initial.ts.plot <- function(CdT.P, 
                            stage.record, 
                            tshift,
                            limni.labels, 
                            grid_limni.ylim,
                            dir.seg.gaug, 
                            seg.iter, 
                            t_Gaug, 
                            h_Gaug, 
                            mcmc.segment, 
                            nS,
                            df_peaks,
                            colors.period){
############################################################################################

  Utot.times     = "90% total uncertainty of change point times"
  MAPtimes       = "Change point times (MAP)"
  col.tflood     = "Time of flood"
  color.segm     = c("90% total uncertainty of change point times"="blue")
  col.segm       = c("Change point times (MAP)"="blue", 
                     "Time of flood"="red")
  
  #
  leg.gaugings           = "Gaugings"
  leg.stage.record       = "Stage record"
  leg.col.gaugings       = c( "Gaugings" = "black")
  leg.col.stage.gauging  = c("Stage record" = "gray70", 
                             "Gaugings" = "black")
  #
  leg.nearest.flood         = "Closest flood to tau MAP"
  leg.change.pointsMAP      = "Change points, tau MAP"
  leg.floods                = "Major flood peaks in the CI"
  leg.col.shifts            = c("Change points, tau MAP"      = "solid", 
                                "Closest flood to tau MAP"    = "dashed",
                                "Major flood peaks in the CI" = "dotted")
  
  
  #
  ts.res         = tshift$tau.MAP;
  tflood         = tshift$tflood;
  Q2.ts          = tshift$tau.q2;
  Q97.ts         = tshift$tau.q97;
  col.ts.distrib = colors.period[1:length(ts.res)] #rainbow((nS-1)) 
  X1             = mcmc.segment
  X2             = X1[,(nS+1):(2*nS-1)]
  
  if (nS ==2) {
    X = X2
  } else {
    X = data.frame(time    = X2[,1], 
                   ord     = rep(1, length(X2[,1])), 
                   colorrr = rep(col.ts.distrib[1], length(X2[,1])))
    for (orderr in 2:(nS-1)) {
         X = rbind(X,    
                   data.frame(time    = X2[,orderr], 
                              ord     = rep(orderr, length(X2[,1])),  
                              colorrr = rep(col.ts.distrib[orderr], length(X2[,1]))) )
    }
  }
  Xfake = X
  Xfake$time = -9999
  Xfake$colorrr = "gray60"
  #--------------------------------------------------------------------------------- plot 1
  # stage record with periods
  initial.ts.plot <- ggplot() 
        if (is.null(df.limni)==FALSE) {
          initial.ts.plot = initial.ts.plot + 
          geom_line(data = df.limni, aes(x = t_limni, y = h_limni, color = leg.stage.record), size = 0.2)+
          scale_x_continuous(name=element_blank(), expand = c(0,0), limits = c(0, tail(t_limni,1)))+
          scale_y_continuous(name=limni.labels[2], limits = grid_limni.ylim[1:2], expand = c(0,0))+
          coord_cartesian(clip = 'off')
        } else {
          initial.ts.plot = initial.ts.plot + 
          scale_x_continuous(name=element_blank(), expand = c(0,0), limits = c(0, tail(t_Gaug,1))) +
          scale_y_continuous(name=limni.labels[2], limits = c(min(h_Gaug), max(h_Gaug)), expand = c(0,0))+
          coord_cartesian(clip = 'off')
        }
  
        initial.ts.plot = initial.ts.plot +
        geom_point(aes(x = t_Gaug, y =h_Gaug, color = leg.gaugings), fill = "gray90", size = 2, pch =21) +
        geom_point(data = CdT.P, aes(x = tP, y = hP), fill = "black", size = 2, pch=21) +
        theme_classic(base_size = 15) +
        ylab(limni.labels[2]) +
        geom_rect(mapping= aes(xmin= Q2.ts, xmax=Q97.ts ,ymin=-Inf, ymax=Inf), fill=col.ts.distrib, alpha=0.2 ) +
        geom_vline(aes(xintercept = -9999, linetype = leg.change.pointsMAP ),  lwd =0.5, color = "gray60") +
        geom_vline(aes(xintercept = -9999, linetype = leg.nearest.flood), lwd =0.5, color = "gray60") +   
        geom_vline(aes(xintercept = -9999, linetype = leg.floods), lwd =0.5, color = "gray60") +
        geom_vline(xintercept  = ts.res, color = col.ts.distrib, lwd =1.5,  linetype = "solid") 
        
        if (!is.null(tflood)) {
          for (ss in 1:length(df_peaks)) {
            initial.ts.plot = initial.ts.plot +
              geom_vline(xintercept  = df_peaks[[ss]]$t_hmax,  color = col.ts.distrib[ss], lwd = 0.5, linetype = "dashed") +
              geom_point(aes(x= df_peaks[[ss]]$t_hmax, y=df_peaks[[ss]]$hmax), color = col.ts.distrib[ss], pch=1, size=3) + 
              geom_vline(xintercept  = tflood[ss],             color = col.ts.distrib[ss], lwd =2,    linetype = "dotted")+
              annotate("text", label = df_peaks[[ss]]$t_hmax,  x = df_peaks[[ss]]$t_hmax, y =grid_limni.ylim[2], size = 4, 
                               colour = col.ts.distrib[ss], parse = TRUE, vjust = 0, angle = 90)

          }
        }
        
        initial.ts.plot = initial.ts.plot +
        scale_colour_manual(  name     = element_blank(), 
                              values   = leg.col.stage.gauging,
                              breaks   = c(leg.stage.record, leg.gaugings),
                              labels   = c(leg.stage.record, leg.gaugings),
                              guide    = guide_legend(override.aes = list(
                                         linetype = c("solid", "blank"), shape = c(NA, 21),  fill = c(NA,"gray90")))) +
        scale_linetype_manual(name   = element_blank(),
                              values = leg.col.shifts,
                              breaks = c(leg.change.pointsMAP, leg.nearest.flood, leg.floods),
                              labels = c(leg.change.pointsMAP, leg.nearest.flood, leg.floods),
                              guide  = guide_legend(override.aes = list(
                              col =c("gray60", "gray60", "gray60")))) +
        #   
        theme(text               = element_text(size=10),
              plot.title         = element_text(hjust = 0.5),
              plot.background    = element_rect(fill = "transparent", color = NA),
              panel.background   = element_rect(fill = "transparent"),
              plot.margin        = unit(c(0.5,0.5,0, 0.5),"cm"),
              axis.title.y       = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0.3)),
              axis.text.x        = element_blank(),
              axis.line          = element_line(colour = "black", size = 0.4, linetype = "solid"),
              axis.ticks.x       = element_line(colour = "black", size = 0.4, linetype = "solid"),
              axis.ticks.y       = element_line(colour = "black", size = 0.4, linetype = "solid"),
              panel.grid.major   = element_blank(), 
              panel.grid.minor   = element_blank(),
              legend.position    = "none",
              legend.key         = element_rect(colour = "transparent", fill = "transparent"),
              legend.background  = element_rect(colour = "transparent", fill = "transparent"))
        


        
        
        
      #---------------------------------------------------------------------------- plot 2
      # Plotting the pdf of the shift times:
        init.ts.dens = ggplot()
        if (is.null(df.limni)==FALSE) {
          init.ts.dens= init.ts.dens + 
            scale_x_continuous(name=element_blank(), expand = c(0,0), limits = c(0, tail(t_limni,1)))
        } else {
          init.ts.dens= init.ts.dens + 
            scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits = c(0, tail(t_Gaug,1)))
        }
        
        init.ts.dens= init.ts.dens +  
          ylab("Scaled pdf")+
          theme_classic(base_size = 15)+ theme(text             = element_text(size=10),
                                               plot.title       = element_text(hjust = 0.5),
                                               plot.background  = element_rect(fill = "transparent", color = NA),
                                               panel.background = element_rect(fill = "transparent"),
                                               plot.margin      = unit(c(0, 0.5, 0.2, 0.85),"cm"),
                                               axis.title.y     = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0.5)),
                                               axis.line        = element_line(colour = "black",  size = 0.4, linetype = "solid"),
                                               axis.ticks.x     = element_line(colour = "black", size = 0.4, linetype = "solid"),
                                               axis.ticks.y     = element_line(colour = "black", size = 0.4, linetype = "solid"),
                                               panel.grid.major = element_blank(), 
                                               panel.grid.minor = element_blank(),
                                               legend.position  = "none",
                                               legend.key       = element_rect(colour = "transparent", fill = "transparent"),
                                               legend.background= element_rect(colour = "transparent", fill = "transparent"))
        if (nS == 2) {
          init.ts.dens= init.ts.dens + geom_density(aes(x= X, ..scaled..),
                                                    fill  = col.ts.distrib[1], 
                                                    colour = NA, alpha = 0.3)+
           scale_fill_manual(    name     = element_blank(),
                                  values =col.ts.distrib[1] )
        } else {
          init.ts.dens= init.ts.dens + geom_density(aes(x= X$time, ..scaled.., group =X$ord,
                                                    fill = X$colorrr), alpha=0.3)+
            scale_fill_manual(    name     = element_blank(),
                                  values = col.ts.distrib,
                                  breaks = col.ts.distrib)
        }

        # merge two plots:
        initial.ts.plot2 = plot_grid(initial.ts.plot, 
                                     init.ts.dens, 
                                     ncol = 1, nrow = 2, rel_heights = c(1, 0.5))
        
        pdf(paste0(dir.seg.gaug,"/Stage_record_segment_it",seg.iter,"_before_adjustment.pdf"), 16, 8 , useDingbats=F)
        print(initial.ts.plot2)
        dev.off()
  return(initial.ts.plot2)
}































##############################################################################################
# Stage h time series  with  shift times:
ts.plot <- function(dir.segment.res,    # directory
                    CdT.P,           # dataframe with gaugings with colors per periods
                    all.gaugings,    # dataframe with all gaugings 
                    df.limni,        # dataframe with stage record df(time, stage) 
                    df.mu,           #      
                    df.shift.times,  # dataframe with segment results
                    axis.labels,     # labels of x axis (e.g. time, stage)
                    grid_limni.ylim, # grid of y axis
                    iteration,       # iteration index
                    mcmc.segment,    # dataframe with all cooked mcmc
                    nS,
                    known.shift.times) {            # optimal number of segments
##############################################################################################
    j=0
    #plot:
    t.plot <- ggplot(data = CdT.P)
    if (is.null(df.limni)==FALSE) {
      t.plot= t.plot + 
        geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "black",size = 0.2)+
        scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    } else {
      t.plot= t.plot + 
        scale_x_continuous(name=element_blank(), expand = c(0,0), limits = c(0,tail(all.gaugings$t,1)))
    }
    
    #***********************************************************************************************
    if (!isFALSE(df.shift.times)) {
      pdf.ts = mcmc.segment[,(nS+1):(2*nS-1)]
      X1 = pdf.ts
      if (length(df.shift.times$tau.MAP) ==1) {
        X = X1
      } else {
        X = data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
        for (orderr in 1:(ncol(pdf.ts))) {
          X =rbind(X, data.frame(time= X1[,orderr], 
                                 ord = rep(orderr, length(X1[,1]))))
        }
      }
    } else {
        X1 = NULL
    }
    

    #************************************************************************************************
    t.plot = t.plot + 
             geom_point(data=all.gaugings, aes(x = t , y= h), size = 4, pch =21, color ="gray90", fill= "gray90")+
             geom_point(data=CdT.P, aes(x = tP , y= hP), size = 4, pch =21, fill= CdT.P$color)
    # if (plot.shift.times.on.limni == TRUE) {
    #   t.plot <- t.plot +
    #     geom_vline(aes(xintercept = data.annotate.gaug.adjust$t.adj), color = "red", lwd =0.3, linetype = "dashed")
    # }  
    t.plot = t.plot +
             scale_y_continuous(name   = axis.labels[2], expand = c(0,0), 
                                limits = c(grid_limni.ylim[1], grid_limni.ylim[2]), 
                                breaks = seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3]))+
             ylab(axis.labels[2]) +
             coord_cartesian(clip = 'off') +
             theme_bw(base_size=20) +
             theme( axis.text         = element_text(size=20)
                    ,axis.title       = element_text(size=30) #, face="bold")
                    ,panel.grid.major = element_blank()
                    ,panel.grid.minor = element_blank()
                    ,legend.text      = element_text(size=20)
                    ,legend.title     = element_text(size=30)
                    ,legend.key.size  = unit(1.5, "cm")
                    ,legend.position  = "none"
                    ,plot.margin      = unit(c(0.5,0.7,0, 0.5),"cm")
                    ,axis.title.y     = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))
                    ,axis.text.x      = element_blank()
                    ,axis.line.x      = element_blank())

    #*************************************************************************************************************
    t.plot2 <- ggplot() 
    if (!isFALSE(df.shift.times)) {
    if (length(df.shift.times$tau.MAP)==1) {
      t.plot2 = t.plot2 + 
        geom_density(aes(x  = X, ..scaled.. /2),
                     fill   = "blue", 
                     colour = NA,
                     alpha  = 0.3)
    } else {
      
      t.plot2 = t.plot2 +
        geom_density(aes(x  = X$time, ..scaled.. /2, group =X$ord),  
                     fill   = "blue", 
                     colour = NA, 
                     alpha  = 0.3)
    }
    t.plot2 = t.plot2 + 
      # annotate("rect", xmin= data.annotate.gaug$q2, xmax=data.annotate.gaug$q97, ymin=0, 
      #         ymax=0.5, fill="blue", alpha=0.1) +
      #geom_segment(data = data.annotate.gaug, aes(x = MAP, y = 0, yend =0.5, xend= MAP), size = 0.7, color ="blue")+
      geom_segment(data  = df.shift.times,
                   aes(x = treal, y = -0.5, yend =0, xend = treal), size = 0.7, color ="red")
    } else {
      if (is.null(df.limni)==FALSE) {
         t.plot2 = t.plot2 + 
         annotate("text",
               x = tail(df.limni$t_limni,1)/2,  
               y = 0.2, 
               label= "No shift detected", color = "gray", size=5) +
         annotate("text",
               x = tail(df.limni$t_limni,1)/2 , 
               y = - 0.2, 
               label = "No shift adjusted", color = "gray", size=5) 
      } else {
         t.plot2 = t.plot2 + 
         annotate("text",
               x = tail(CdT.P$tP,1)/2,  
               y = 0.2, 
               label= "No shift detected", color = "gray", size=6) +
         annotate("text",
                 x = tail(CdT.P$tP,1)/2 , 
                 y = - 0.2, 
                 label = "No shift adjusted", color = "gray", size=6) 
      }
    } 
      if (is.null(known.shift.times)){
          t.plot2 = t.plot2 + 
          scale_y_continuous(name = "pdf", expand = c(0,0), 
                               limits = c(-0.5,0.5),  breaks =c(-0.5, 0, 0.5),
                               labels = c(grid_limni.ylim[1], grid_limni.ylim[1], grid_limni.ylim[2])) +
          geom_hline(yintercept = c(0, 0.5), color="darkgray", linetype="dashed", size = 0.5)
      } else {
          t.plot2 = t.plot2 + 
          scale_y_continuous(name = "pdf", expand = c(0,0), 
                            limits = c(-1,0.5),  breaks =c(-1, -0.5, 0, 0.5),
                            labels = c(grid_limni.ylim[1], grid_limni.ylim[1], 
                                       grid_limni.ylim[2], grid_limni.ylim[2])) +   
          geom_hline(yintercept = c( -0.5,0, 0.5), color = "darkgray", linetype ="dashed", size = 0.5)+
          geom_point(aes(x = known.shift.times, y = -1), color = "black", size = 4, shape = 4, stroke = 2)
      }
      t.plot2 = t.plot2 + 
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      theme_bw(base_size=20)+
      theme(  axis.text        = element_text(size=20)
             ,axis.title       = element_text(size=30) #,face="bold")
             ,panel.grid.major = element_blank()
             ,panel.grid.minor = element_blank()
             ,legend.text      = element_text(size=20)
             ,legend.title     = element_text(size=30)
             ,legend.key.size  = unit(1.5, "cm")
             ,legend.position  = "none"
             ,plot.margin      = unit(c(0.5,0.7,0, 0.5),"cm")
             ,axis.title.y     = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), color= "white")
             ,axis.text.y      = element_text(color="white"))
    if (is.null(df.limni)==FALSE) {
        t.plot2 = t.plot2 + 
        scale_x_continuous(name = limni.labels[1], expand = c(0,0), limits = c(0,tail(df.limni$t_limni,1)))
    } else {
        t.plot2 = t.plot2 + 
        scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits = c(0,tail(all.gaugings$t,1)))
    }
    #***************************************************************************************************************
    t.plot3 = plot_grid( t.plot,
                         t.plot2, 
                         ncol = 1,
                         nrow = 2, 
                         rel_heights = c(1, 0.5))    # needs cowplot "package"
    
    pdf(paste0(dir.segment.res,"/Stage_record_segment_it",iteration,".pdf"), 16, 8 , useDingbats=F)
    print(t.plot3)
    dev.off()
    return(list(t.plot3, X1))
}




















###################################################################################################
#exemple of ts adjustment for the case study of Meyras (Figure for the paper):
exemple_ts.plot <- function(CdT.P, df.limni, Q10.ts, Q90.ts, ts.res, ts.real, station.name, 
                            limni.labels, grid_limni.ylim,  dir.seg.gaug, seg.iter, 
                            t_Gaug, h_Gaug, pdf.ts, time_period, df_peaks, save_name) {
###################################################################################################
  #col.ts.distrib = rainbow((nS-1)) 
  X1=pdf.ts
  if (length(ts.res)==1) {
      X = X1
  } else {
      X =data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
      for (orderr in 1:(ncol(pdf.ts))) {
        X =rbind(X, data.frame(time= X1[,orderr], ord = rep(orderr, length(X1[,1]))))
      }
  }
  #
  leg.gaugings           = "Gaugings"
  leg.stage.record       = "Stage record"
  leg.col.gaugings       = c( "Gaugings" = "black")
  leg.col.stage.gauging  = c("Stage record" = "gray10",  "Gaugings" = "black")
  leg.shift.times        = "Shift times, s"
  leg.change.points      = "change points, tau"
  leg.col.shifts         = c("change points, tau" = "solid", "Shift times, s" = "dotted")

  
  ts.exemple.plot <- ggplot()+
      geom_vline(xintercept = df_peaks$t_hmax, col= "red", lwd = 1) +
      geom_vline(aes(xintercept = ts.res, linetype = leg.change.points), col="blue", lwd =1.5) +
      #geom_segment(aes(x=ts.res, xend=ts.real, y = -1, yend =-0.8), col="red", lwd=2)+
      geom_vline(aes(xintercept = ts.real, linetype = leg.shift.times), col="red", lwd =1.5) +
      annotate("text", label = df_peaks$t_hmax , x = df_peaks$t_hmax, y =grid_limni.ylim[2], size = 4, colour = "red", parse = TRUE, vjust = 0, angle = 90)+
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni, color =leg.stage.record ), size = 1.2)+
      #scale_x_continuous(name=limni.labels[1], limits = c(1800, 2700) , expand = c(0,0)) +
      scale_x_continuous(name=limni.labels[1], limits = time_period, expand = c(0,0)) +
      scale_y_continuous(name=limni.labels[2], limits = c(grid_limni.ylim[1], grid_limni.ylim[2]), expand = c(0,0)) +
      xlab(limni.labels[1]) + ylab(limni.labels[2]) +
      #annotate("rect",xmin= (Q10.ts), xmax=(Q90.ts), ymin=-Inf, ymax=Inf, fill="blue", alpha=0.2) +
      #geom_vline(xintercept = ts.res, col= "blue", lwd = 1)+
      #geom_rect(mapping= aes(xmin= Q10.ts, xmax=Q90.ts ,ymin=-Inf, ymax=Inf), fill="blue", alpha=0.1) +
      coord_cartesian(clip = 'off') +
      #geom_point(aes(x=ts.real, y = c(-1)),  col="red", size = 2 , stroke =2,shape= 4)+
      # annotate("text", label = expression(t[s]^MAP), x = ts.res-2, y = -0.7, size = 3, colour = "blue", 
      #          parse = TRUE, vjust = 0, angle = 90)+
      # annotate("text", label = expression(t[flood]), x = ts.real+2, y = -0.7, size = 3, colour = "red", 
      #          parse = TRUE, angle = 90)+ 
      scale_linetype_manual(name   = element_blank(),
                            values = leg.col.shifts,
                            breaks = c(leg.change.points, leg.shift.times),
                            labels = c(unname(TeX("$ \\hat{\\tau_{1}}$")), unname(TeX("$t_{flood,1}$"))),
                            guide  = guide_legend(override.aes = list(col =c("blue", "red")))) +
      theme_bw(base_size=15)+
      theme(axis.text       = element_text(size=15),
          axis.title        = element_text(size=20),
          #,panel.grid.major= element_line(size=1.2)
          #,panel.grid.minor= element_line(size=0.8)
          plot.margin       = unit(c(1, 1.5, 1.2, 0.5),"cm"),
          axis.title.y      = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)),
          axis.line         = element_line(colour = "black",  size = 0.4, linetype = "solid"),
          axis.ticks.x      = element_line(colour = "black",  size = 0.4, linetype = "solid"),
          axis.ticks.y      = element_line(colour = "black",  size = 0.4, linetype = "solid"),
          axis.title.x      = element_blank(),
          panel.grid.major  = element_blank(), 
          panel.grid.minor  = element_blank(),
          legend.key.size   = unit(1, "cm"),
          legend.position   = "left",
          legend.box        = "horizontal",
          legend.direction  = "horizontal",
          legend.box.margin = margin(0, 0, 0, 50),
          #legend.spacing.x = unit(1.3, 'cm'),
          legend.text       = element_text(margin=margin(0,30,0,0.1), size=22),
          legend.justification = c(0,1),
          legend.key        = element_rect(colour = "transparent", fill = "transparent"),
          legend.background = element_rect(colour = "transparent", fill = "transparent"))
      
      if (!is.null(CdT.P)){
            ts.exemple.plot = ts.exemple.plot +
            #geom_point(aes(x = t_Gaug, y =h_Gaug), color = "gray60", size = 0.8) +
            geom_point(data = CdT.P, aes(x = tP, y = hP, color = leg.gaugings), pch = 21, size = 4, fill="blue")+
            scale_colour_manual(name     = element_blank(), 
                                values   = leg.col.stage.gauging,
                                breaks   = c(leg.stage.record, leg.gaugings),
                                labels   = c(leg.stage.record, leg.gaugings),
                                guide    = guide_legend(override.aes = list(
                                  linetype = c("solid", "blank"), shape = c(NA, 21),  fill = c(NA,"black"))))
      } else {
          ts.exemple.plot = ts.exemple.plot +
          scale_colour_manual(name     = element_blank(), 
                              values   = leg.col.stage.gauging[1],
                              breaks   = c(leg.stage.record),
                              labels   = c(leg.stage.record),
                              guide    = guide_legend(override.aes = list(
                              linetype = c("solid"),  fill = c(NA))))
      }
      


  
  #----------------------------------------------------------------------------------------
  Utau = "Posterior pdf of the change point times"
  leg.col.uncert = c("Posterior pdf of the change point times" = "blue")
  ts.exemple.plot2 <- ggplot() 
  ts.exemple.plot2 = ts.exemple.plot2 +
                     geom_density(aes(x= X$time, ..scaled.., group =1, fill= Utau),  colour=NA, alpha=0.3)
  ts.exemple.plot2 = ts.exemple.plot2 + 
    geom_segment( aes(x = ts.res,  y = 0, yend =1, xend= ts.res), size = 1.5, color ="blue")+
    geom_segment( aes(x = ts.real, y = 0, yend =1, xend= ts.real ), col="red", linetype="dotted", size =1.5) +
    scale_y_continuous(name  = "Scaled pdf", expand = c(0,0), limits = c(0,1)) + 
    scale_fill_manual(name   = element_blank(), 
                      breaks = Utau,
                      labels = unname(TeX("$Posterior \\; pdf \\; of \\; change \\; point \\; time \\; \\tau_{1} $")),
                      values = leg.col.uncert) +
    xlab(limni.labels[1]) + coord_cartesian(clip = 'off')+
    theme_bw(base_size=15) +
    theme(axis.text          = element_text(size=15)
          ,axis.title        = element_text(size=20)
          ,panel.grid.major  = element_blank()
          ,panel.grid.minor  = element_blank()
          ,plot.margin       = unit(c(0.5, 1.5, 0.5, 0.5),"cm")
          ,axis.title.y      = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))
          ,legend.key.size   = unit(1, "cm")
          #,legend.spacing.x = unit(0.2, "cm")
          ,legend.position   ="left"
          ,legend.box        = "horizontal"
          ,legend.direction  = "horizontal"
          ,legend.justification = c(0,1)
          ,legend.text       = element_text(size=22)
          ,legend.key        = element_rect(colour = "transparent", fill = "transparent")
          ,legend.background = element_rect(colour = "transparent", fill = "transparent"))+
    # scale_x_continuous(name=limni.labels[1], limits = c(1800, 2700) , expand = c(0,0))
    scale_x_continuous(name=limni.labels[1], limits = time_period, expand = c(0,0)) 

  
  
  
  ########################################################################################
  # Legend plot:
  Legend.title= ggdraw()+ draw_label("Legend",fontface ='bold',x = 0, hjust = 0,size =22)+
  theme(# add margin on the left of the drawing canvas,
          # so title is aligned with left edge of first plot
          plot.margin = margin(5, 0, 5, 0))
  glegend.A <- cowplot::get_legend( ts.exemple.plot )
  glegend.B <- cowplot::get_legend( ts.exemple.plot2 )
  glegends = plot_grid(glegend.A, glegend.B, ncol =2, nrow =1,  rel_widths = c(1, 0.75))
  Legends =  plot_grid(Legend.title, glegends,  ncol = 1,
                       # rel_heights values control vertical title margins
                       rel_heights = c(0.4, 1))
  ts.exemple.plot.without.legend = ts.exemple.plot + theme(legend.position="none")
  ts.exemple.plot2.without.legend = ts.exemple.plot2 + theme(legend.position="none")
  
  #save plot:
  ts.exemple.plot.3 = plot_grid( ts.exemple.plot.without.legend,
                                 ts.exemple.plot2.without.legend,
                                 #Legends,
                                 ncol = 1, nrow = 2, rel_heights = c(1, 0.5)) #, 0.3))
  
  # ggsave(ts.exemple.plot.3, filename = paste0(dir.seg.gaug,"/ts_exemple_paper",seg.iter,".png"), 
  #        device = "png", width = 16, height =8, dpi = 600, units = "in")
  # 
  pdf(paste0(dir.seg.gaug,"/", save_name, ".pdf"), 16, 8 , useDingbats=F)
  # pdf(paste0(dir.seg.gaug,"/Figure3.pdf"), 16, 8 , useDingbats=F)
  print(ts.exemple.plot.3)
  dev.off()
  return(ts.exemple.plot.3)
}




























#***********************************************************************************************************************
# Plotting final RC gaugings
final.RC.plot <- function(gaug, 
                          RCmaxpost, 
                          colorGauging, 
                          dir,
                          grid_RC.ylim, grid_RC.xlim, grid_RC.xstep, grid_RC.ystep,
                          ticks_RC.y.log,
                          RC.x.labels, 
                          RC.y.labels,
                          grid_RC.ylim.log) {
#***********************************************************************************************************************
  #-------------
  # Linear plot:
  #-------------
  #to do!
  
  #----------
  # Log plot:
  #----------
  p <- ggplot() + 
  coord_trans(xlim = grid_RC.xlim, 
              ylim = grid_RC.ylim.log) +
  scale_y_log10(na.value = -10,  
                breaks   = ticks_RC.y.log,   
                labels   = ticks_RC.y.log) +
  scale_x_continuous(breaks = seq(grid_RC.xlim[1], grid_RC.xlim[2],  grid_RC.xstep)) +  
  theme_light(base_size = 15)+
  ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
  xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]")))
  j=0
  for (i in 1:length(RCmaxpost)) {
  if (!is.null(RCmaxpost[[i]])) {
    j=j+1
    #p <- p + geom_ribbon(data = na.omit(maxpost.data[[i]]), aes(x = x, ymin = (y-y2) ,
    #                      ymax = (y+y2)), fill = colo[i], alpha = 0.3) 
    #p <- p + geom_line(data = na.omit(RCmaxpost[[i]]), aes(x = V1, y= V1.1), color = colo[i],size = 0.3, na.rm = TRUE)
    #p <- p + env[[i]]
    #p <- p + env.par[[i]]
    p <- p + geom_point(data = gaug[[i]],    aes(x = X.h. , y= X.Q.), size = 3, pch =21, fill = colorGauging[j])
    p <- p + geom_errorbar(data = gaug[[i]], aes(x = X.h.,  ymin =X.Q.-2*X.uQ. , ymax =X.Q. +2*X.uQ.), 
                           width=0.02*(grid_RC.xlim[2] - grid_RC.xlim[1]), size = 0.3, color = colorGauging[j])
  }
  }
  #coord_cartesian(xlim = c(-1, 2), ylim = c(0,6), expand = c(0))+
  p = p + 
  theme_bw(base_size=20)+
  theme( axis.text        = element_text(size=20)
        ,axis.title       = element_text(size=30, face="bold")
        ,panel.grid.major = element_line(size=1.2)
        ,legend.text      = element_text(size=20)
        ,legend.title     = element_text(size=30)
        ,legend.key.size  = unit(1.5, "cm")
        ,legend.position  = "none")+
  annotation_logticks(base = 10, sides = "l", scaled = TRUE, colour = "black", size = 0.5, linetype = 1)+
  scale_colour_manual(name="Error Bars",values=colorGauging) + theme(legend.position = c(0.8, 0.2))
  ggsave(p, filename =paste0(dir,"/RClog_shifts.png"),  device = "png", width = 16, height =8,
         dpi = 600, units = "in")
  #--------
  #log-log:
  #--------
  # ploglog = p + 
  #           scale_x_log10(na.value=-10, breaks=c(0.01, 0.1,1,10,100)) 
  #   
  # ggsave(ploglog, filename =paste0(dir.case_study,"/Results/segmentation_gaugings/RCloglog_shifts.png"), 
  #        bg = "transparent", width = 6, height =3.5, dpi = 400)
}






















#***********************************************************************************************************************
# Plotting final RC gaugings
final.RC.plot.2 <- function(gaug, 
                            dir,
                            grid_RC.ylim, grid_RC.xlim, grid_RC.xstep, grid_RC.ystep,
                            ticks_RC.y.log,
                            RC.x.labels, 
                            RC.y.labels,
                            grid_RC.ylim.log,
                            colors.period) {
  #***********************************************************************************************************************
  #----------
  # Log plot:
  #----------
  col.gaug = colors.period[1:tail(gaug$X.Period.,1)]
  
  p <- ggplot(data=gaug , aes(x=X.h., y=X.Q., fill = as.factor(X.Period.)))  + 
       coord_trans(xlim = grid_RC.xlim, ylim = grid_RC.ylim.log) +
       scale_y_log10(na.value = -10,  
                     breaks   = ticks_RC.y.log,   
                     labels   = ticks_RC.y.log) +
       scale_x_continuous(breaks = seq(grid_RC.xlim[1], grid_RC.xlim[2],  grid_RC.xstep)) +  
       theme_light(base_size = 15)+
       ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
       xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]")))+
       geom_errorbar(data=gaug, aes(x= X.h., ymin =X.Q.-2*X.uQ. , ymax =X.Q. +2*X.uQ.), 
                     width=0.02*(grid_RC.xlim[2] - grid_RC.xlim[1]), size = 0.3, color = gaug$color) +
       geom_point(size = 3, pch =21) + 
       #data = gaug, aes(x = X.h. , y= X.Q.), size = 3, pch =21, fill = gaug$color)+
       annotation_logticks(base = 10, sides = "l", scaled = TRUE, colour = "black", size = 0.5, linetype = 1)+
       scale_fill_manual(name  = "Periods",
                         values = col.gaug, 
                         breaks = unique(gaug$X.Period.),
                         labels = unique(gaug$X.Period.)) +  # c(paste0("P", seq(1,tail(gaug$X.Period.,1)))))+
       #scale_fill_brewer(palette="Spectral")+
       #labs(fill="X.Period.")+
       theme_bw(base_size=20) +
       theme( axis.text        = element_text(size=20)
             ,axis.title       = element_text(size=30, face="bold")
             ,panel.grid.major = element_blank() #element_line(size=1.2)
             ,panel.grid.minor = element_blank()
             ,legend.text      = element_text(size=10)
             ,legend.title     = element_text(size=20)
             ,legend.key.size  = unit(1, "cm")
             ,legend.position  = "right"
             ,legend.direction = "vertical")
  ggsave(p, filename = paste0(dir,"/RClog_shifts.png"),  device = "png", width = 16, height =8,
         dpi = 500, units = "in")
}




























#**********************************************************************************************************************
# Plotting final SHIFT TIMES in stage time series
t.final.plot <- function(t_limni,h_limni, t.q10, t.q90, RCmaxpost, gaug, data.annotate.off, data.annotate.gaug, colo) {
#********************************************************************************************************************** 
  color.ts = c("Gaugings segmentation"="blue", "Official segmentation"="black")
  ts.gaug = "Gaugings segmentation"; ts.offic = "Official segmentation"
  t.plot <- ggplot()
  if (exists("t_limni")==TRUE) { 
     if (is.null("t_limni")==FALSE) {
     t.plot <- t.plot + 
     geom_line(aes(x = t_limni, y = h_limni), color = "darkgrey", size = 0.2)
  }}
  
  for (i in 1:length(RCmaxpost)) {
  if (!is.null(RCmaxpost[[i]])) {
    t.plot <- t.plot + geom_point(data = gaug[[i]], aes(x = X.tP. , y= X.h.), size = 2, pch =21, fill =colo[i])
  }
  }
   
  if (exists("data.annotate.off")==TRUE)  {
     if (is.null("data.annotate.off")==FALSE) {
     t.plot <- t.plot +
     geom_point(data = data.annotate.off, aes(x = x, y = 1, color= ts.offic), size = 1.2, shape =4, stroke=0.5 )
     }
  }
  # geom_segment(data = data.annotate.off, aes(x = x, y = start, yend = finish, xend = x, 
  #                                            color = "Official segment."),size = 0.3) +
  t.plot <- t.plot + 
  geom_segment(data = data.annotate.gaug, aes(x = x, y = start, yend = finish, xend = x,
               color = ts.gaug), size = 0.4) +
  scale_x_continuous(name="time [days]", expand = c(0,0)) +
  #scale_y_continuous(name="Stage h [m]", expand = c(0,0), limits = c(1,6))+
  scale_y_continuous(name="Stage h [m]", expand = c(0,0), limits = c(0,40))+  
  xlab("Time [days]")+ ylab("Stage h [cm]")+
  annotate("rect", xmin= t.q10, xmax=t.q90, ymin=1, ymax=1.5, fill="blue", alpha=0.1) +
  coord_cartesian(clip = 'off')+
  
  scale_colour_manual(name = element_blank(), values=color.ts,
                      breaks=c(ts.gaug, ts.offic),
                      guide = guide_legend(override.aes = list(
                      linetype = c("solid", "blank"), shape = c(NA, 4)))) +
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(hjust = 0.5), 
          plot.background = element_rect(fill ="transparent", color = NA),
          #panel.grid.major=element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill ="transparent"), 
          #axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.margin=unit(c(0.3,0.5,0.05,1.05),"cm"), 
          #legend.position =c(0.7, 0.2), legend.box = "vertical", legend.key.size = unit(0.5, "cm"),
          #legend.position = "none")
          legend.position = "bottom") #legend.direction="vertical",
  # legend.text=element_text(size=12),
  # legend.spacing.x = unit(0.05, 'cm'),
  # legend.spacing.y = unit(0.01, 'cm'),
  # legend.key.height=unit(0.5,"line"),
  # legend.key.width=unit(0.8,"line"),
  # legend.key = element_rect(colour = "transparent", fill = "transparent"),
  # legend.background = element_rect(colour = "transparent", fill = "transparent"))
  
  ggsave(t.plot, filename =paste0(dir.case_study,"/Results/segmentation_gaugings/time_series.png"), bg = "transparent",
       width = 7, height =3.5, dpi = 400)
}
























# Plot stage h time series  with shift times ( + official times) with results of gaugings segmentation:
#######################################################################################################################
plot.time.shifts.gaugings <- function(dir, 
                                      g.s, 
                                      data.annotate.off, 
                                      data.annotate.gaug, 
                                      t_Gaug, h_Gaug, 
                                      df.limni, 
                                      limni.labels, 
                                      grid_limni.ylim,
                                      plot.shift.times.on.limni,
                                      pdf.ts,
                                      dates) {
#######################################################################################################################
    #plot 1:
    # filter the time series of stage record removing the long periods with missing data (putting a NA instead):
    if (!is.null(df.limni)){
       dt_limni = new_NA_limni= 0
       t_limni_filtered = df.limni$t_limni
       h_limni_filtered = df.limni$h_limni
    
       for (tt in 2:length(df.limni$t_limni)){
         dt_limni[tt] =  df.limni$t_limni[tt] - df.limni$t_limni[tt-1]
       }
       for (tt in 2:length(df.limni$t_limni)){
          if (dt_limni[tt] > 10*mean(dt_limni)){
              t_limni_filtered[(tt+ new_NA_limni): (length(t_limni_filtered)+1)] <- c(t_limni_filtered[tt+ new_NA_limni] -1, 
                                                                                      t_limni_filtered[(tt+ new_NA_limni) : length(t_limni_filtered)])
              h_limni_filtered[(tt+ new_NA_limni): (length(h_limni_filtered)+1)] <- c(NA,  h_limni_filtered[(tt+ new_NA_limni) : length(h_limni_filtered)])
              new_NA_limni = new_NA_limni +1 
          }
       } 
       df.limni_filtered = data.frame(t_limni = t_limni_filtered,  h_limni = h_limni_filtered)
    }
    
    
    
    # stage record plot:
    t.plot = ggplot()
    if (!is.null(df.limni)){
      t.plot= t.plot + 
        geom_line(data = df.limni_filtered, aes(x = t_limni, y = h_limni), color = "gray70",size = 0.2)+
        scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    } else {
      t.plot= t.plot + 
        scale_x_continuous(name=element_blank(), expand = c(0,0), limits = c(0,tail(g.s$X.tP.,1)))
    }
    t.plot <- t.plot + 
    geom_point(data = g.s, aes(x = X.tP. , y= X.h.), size = 4, pch =21, fill= g.s$color)
    
    if (plot.shift.times.on.limni == TRUE) {
       t.plot <- t.plot +
       geom_vline(aes(xintercept = data.annotate.gaug$t.adj), color = "red", lwd =0.3, linetype = "dashed")
    }  
    t.plot <- t.plot +
      scale_y_continuous(name   = limni.labels[2], expand = c(0,0), limits = c(grid_limni.ylim[1], grid_limni.ylim[2]), 
                         breaks = seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3]))+
      ylab(limni.labels[2])+
      coord_cartesian(clip = 'off')+
      theme_bw(base_size=10)+
      theme( axis.text        = element_text(size=20)
            ,axis.title       = element_text(size=30) #, face="bold")
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,legend.text      = element_text(size=20)
            ,legend.title     = element_text(size=30)
            ,legend.key.size  = unit(1.5, "cm")
            ,legend.position  = "none"
            ,plot.margin      = unit(c(0.5,0.7,0, 0.5),"cm")
            ,axis.title.y     = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))
            ,axis.text.x      = element_blank()
            ,axis.line.x      = element_blank())
    
    
    if (!is.null(data.annotate.gaug$t.adj)){
      if (!is.null(dates)){
        t.plot <- t.plot+
          geom_segment(aes(x     = data.annotate.gaug$t.adj, 
                           xend  = data.annotate.gaug$t.adj, 
                           y     = grid_limni.ylim[1],
                           yend  = grid_limni.ylim[1] + 0.10),
                       color = "red", 
                       size  = 0.5)
        t.plot <- t.plot +
                annotate("text", 
                         x     = data.annotate.gaug$t.adj, 
                         y     = grid_limni.ylim[1]+0.25,
                         label = substring(dates,7,10), color = "red", size=4)
      }
    }

    
    ############################################
    # shift times pdf:
    ############################################
    if (!is.null(data.annotate.gaug$MAP)) {
      X1=pdf.ts
      if (length(data.annotate.gaug$MAP) ==1) {
        X = X1
      } else {
        X =data.frame(time= X1[,1], ord = rep(1, length(X1[,1])))
        for (orderr in 1:(ncol(pdf.ts))) {
          X =rbind(X, data.frame(time= X1[,orderr], ord = rep(orderr, length(X1[,1]))))
        }
      }
    }
    # plot 2:
    t.plot2 <- ggplot() 
    if (!is.null(data.annotate.gaug$MAP)) {
    if (length(data.annotate.gaug$MAP)==1) {
      t.plot2 = t.plot2 + 
        geom_density(aes(x= X[,1], ..scaled.. /2),
                     fill="blue", 
                     colour=NA, alpha=0.3)
    } else {
      
      t.plot2 = t.plot2 +
        geom_density(aes(x= X$time, ..scaled.. /2, group =X$ord),  
                     fill= "blue", 
                     colour=NA, alpha=0.3)
    }
      t.plot2 = t.plot2 +
      geom_segment(data  = data.annotate.gaug, 
                   aes(x = t.adj, y = -0.5, yend = 0, xend=t.adj), 
                   size = 0.7, color ="red")
    } else {
      if (is.null(df.limni)==FALSE) {
        t.plot2 = t.plot2 + 
          annotate("text",
                   x = tail(df.limni$t_limni,1)/2,  
                   y = 0.2, 
                   label = "No shift detected", color = "gray", size=6) +
          annotate("text",
                   x = tail(df.limni$t_limni,1)/2 , 
                   y = - 0.2, 
                   label = "No shift adjusted", color = "gray", size=6) 
      } else {
        t.plot2 = t.plot2 + 
          annotate("text",
                   x    = tail(g.s$X.tP.,1)/2,  
                   y    = 0.2, 
                   label = "No shift detected", color = "gray", size=6) +
          annotate("text",
                   x     = tail(g.s$X.tP.,1)/2, 
                   y     = - 0.2, 
                   label = "No shift adjusted", color = "gray", size=6) 
      }
    }
    
    t.plot2 = t.plot2 + 
      # annotate("rect", xmin= data.annotate.gaug$q2, xmax=data.annotate.gaug$q97, ymin=0, 
      #         ymax=0.5, fill="blue", alpha=0.1) +
      #geom_segment(data = data.annotate.gaug, aes(x = MAP, y = 0, yend =0.5, xend= MAP), size = 0.7, color ="blue")
      scale_y_continuous(name = "pdf", expand = c(0,0), 
                         limits = c(-1,0.5),  breaks =c(-1, -0.5, 0, 0.5),
                         labels = c(grid_limni.ylim[1], grid_limni.ylim[1], grid_limni.ylim[2],grid_limni.ylim[2])) +
      xlab(limni.labels[1])+ 
      coord_cartesian(clip = 'off')+
      geom_hline(yintercept = c( -0.5,0, 0.5), color="darkgray", linetype="dashed", size = 0.5)+
      theme_bw(base_size=10)+
      theme( axis.text        = element_text(size=20)
            ,axis.title       = element_text(size=30) #,face="bold")
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,legend.text      = element_text(size=20)
            ,legend.title     = element_text(size=30)
            ,legend.key.size  = unit(1.5, "cm")
            ,legend.position  = "none"
            ,plot.margin      = unit(c(0.5,0.7,0, 0.5),"cm")
            ,axis.title.y     = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), color= "white")
            ,axis.text.y      = element_text(color="white"))
 
    if (is.null(df.limni)==FALSE) {
      t.plot2= t.plot2 + 
        scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
    } else {
      t.plot2= t.plot2 + scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits = c(0,tail(g.s$X.tP.,1)))
    }
    
    if (exists("data.annotate.off")==TRUE)  {
      if (is.null(data.annotate.off)==FALSE) {
        t.plot2 = t.plot2 +
          geom_point(data = data.annotate.off, aes(x = xeffect, y = -1), color= "black", size = 4, shape =4, stroke=2 )
        #geom_point(data = data.annotate.off, aes(x = xpotent, y = -1), color= "green", size = 4, shape =4, stroke=2 )
        
      } else {
        if (is.null(df.limni) == FALSE) {
          t.plot2 = t.plot2 + 
            annotate("text",
                     x    = tail(df.limni$t_limni,1)/2,  
                     y    = -0.85, 
                     label= "Official shift dates not available", color = "gray", size=6)
        } else {
          t.plot2 = t.plot2 + 
            annotate("text",
                     x    = tail(g.s$X.tP.,1)/2,  
                     y    = -0.85, 
                     label= "Official shift dates not available", color = "gray", size=6)
        }
      }
   }
  # Final plot:
  #*************************************************************************************************************
  t.plot3 = plot_grid( t.plot, t.plot2, ncol = 1, nrow = 2, rel_heights = c(1, 0.5)) #needs cowplot "package"
  pdf(paste0(dir,"/Segmented_stage_record.pdf"), 16, 8 , useDingbats = F)
  print(t.plot3)
  dev.off()
  #*************************************************************************************************************
}





























#######################################################################################################################
plot.time.shifts.gaugings.and.bt <- function(dir , nperiod, g.s, 
                                             data.annotate.off, data.annotate.off.bis,
                                             data.annotate.gaug, data.annotate.gaug.adjust, 
                                             colo, t_Gaug, h_Gaug, limni.labels, start.y.legend, grid_limni.ylim,
                                             dir.SPD.exe,  df.limni,  dir.SPD.results,
                                             shift.times, t.shift.for.b, h_G , 
                                             t_G, color_G, times.uncert, officialShiftsTime , ylimits) {
  #######################################################################################################################
  bt.MAP = as.numeric(read.table(paste(dir.SPD.exe,"/Results_Summary.txt",sep="")
                                 ,row.names=1,dec=".",sep="")[16,1:nperiod])
  bt.q10 = as.numeric(read.table(paste(dir.SPD.exe,"/Results_Summary.txt",sep="")
                                 ,row.names=1,dec=".",sep="")[7,1:nperiod])
  bt.q90 = as.numeric(read.table(paste(dir.SPD.exe,"/Results_Summary.txt",sep="")
                                 ,row.names=1,dec=".",sep="")[10,1:nperiod])
  shifts = t.shift.for.b 
  shifts.all = shift.times
  if (times.uncert ==TRUE) {   # /!\ this has to be changed !!!!!!!!!!!!!!
    # t.shifts = sort(c(shifts[,3]))
    t.shifts.before = c(0,shifts$treal)
    t.shifts.plus = c(shifts$treal, tail(df.limni$t_limni,1))
    # t.shifts.MAP = sort(shifts[,1])
    # t2.shifts = sort(c(shifts[,2]))
    # t97.shifts = sort(c(shifts[,4]))
    df.bt = data.frame(bt= bt.MAP, btq10 = bt.q10, btq90=bt.q90, t.shifts.before, t.shifts.plus)
  }else {
    t.shifts = shifts[-length(shifts)]
    t.shifts.before = c(0,t.shifts)
    t.shifts.plus = sort(c(t.shifts, tail(df.limni$t_limni,1)))
  }
  #plot:
  j=0
  t.plot <- ggplot()
  if (is.null(df.limni)==FALSE) {
    t.plot= t.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "black",size = 0.2)+
      scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
  } else {
    t.plot= t.plot + scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(t_Gaug,1)))
  }
  if (exists("data.annotate.off")==TRUE)  {
    if (is.null(data.annotate.off)==FALSE) {
      t.plot <- t.plot +
        geom_point(data = data.annotate.off, aes(x = x, y = start), color= "black", size = 4, shape =4, stroke=2 )
    }
  }
  if (exists("data.annotate.off.bis")==TRUE)  {
    if (is.null(data.annotate.off.bis)==FALSE) {
      t.plot <- t.plot +
        geom_point(data = data.annotate.off.bis, aes(x = x, y = start), color= "red", size = 4, shape =15)
    }
  }
  #
  if (times.uncert ==TRUE) {
    #geom_line(data = df.bt, aes(x=t.shifts.before, y =bt.MAP), color = "red", size = 1 )+
    t.plot = t.plot +
      #geom_rect(mapping= aes(xmin= shifts$t2, xmax=shifts$t97 ,ymin=-Inf, ymax=Inf), fill="blue", alpha=0.1) +
      #geom_vline(aes(xintercept = shifts$tMAP), color = "blue", lwd =0.5, linetype = "solid")+
      geom_vline(aes(xintercept = shifts$treal), color = "red", lwd =0.3, linetype = "dashed")+
      geom_segment(mapping= aes(x =t.shifts.before , y = bt.MAP, xend = t.shifts.plus, yend = bt.MAP),color = "red", size = 1) +
      geom_rect(mapping = aes(xmin= t.shifts.before, xmax=t.shifts.plus, ymin=bt.q10, ymax= bt.q90), fill="red", alpha=0.3) +
      geom_line(data = df.limni, aes(x = t_limni, y= h_limni), color="black", size =0.3)+
      geom_point(aes(x=t_G, y = h_G), color=color_G, size=3)+
      scale_y_continuous(name = expression("Stage h [m]"), limits =ylimits) + 
      theme_bw(base_size=20)+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
            ,legend.text=element_text(size=20),legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm"),legend.position="none"
            ,panel.grid = element_blank())
  } else {
    t.plot = t.plot +
      geom_vline(aes(xintercept = shifts.all), color = "blue", lwd =0.7, linetype = "solid") +
      geom_segment(mapping= aes(x =t.shifts.before , y = bt.MAP, xend = t.shifts.plus, yend = bt.MAP), 
                   color = "red", size = 1) +
      geom_rect(mapping = aes(xmin= t.shifts.before, xmax=t.shifts.plus, ymin=bt.q10,
                              ymax= bt.q90), fill="red", alpha=0.3) +
      geom_line(data = df.limni, aes(x = t_limni, y= h_limni), color="black", size =0.3)+
      scale_y_continuous(name = expression("Stage h [m]"), limits =ylimits) + 
      geom_point(aes(x=t_G, y = h_G), color=color_G, size=3)+
      theme_bw(base_size=20)+
      theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
            ,legend.text=element_text(size=20),legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm"),legend.position="none"
            ,panel.grid = element_blank())
  }
  #
  if (exists("data.annotate.gaug.adjust")==TRUE)  {
    if (is.null(data.annotate.gaug.adjust)==FALSE) {
      t.plot <- t.plot +
        geom_segment(data = data.annotate.gaug.adjust, aes(x = t.adj, xend = t.adj, y = start, yend = finish),
                     size = 0.4, color ="blue")
    }
  }
  t.plot <- t.plot + 
    annotate("rect", xmin= data.annotate.gaug$q2, xmax=data.annotate.gaug$q97, ymin=data.annotate.gaug$start, 
             ymax=data.annotate.gaug$finish, fill="blue", alpha=0.1) +
    geom_segment(data = data.annotate.gaug, aes(x = MAP, y = start, yend = finish, xend= MAP), size = 0.4, color ="blue")+
    scale_y_continuous(name=limni.labels[2], expand = c(0,0), limits = c(start.y.legend-1.5, grid_limni.ylim[2]), 
                       breaks=seq(grid_limni.ylim[1], grid_limni.ylim[2], grid_limni.ylim[3]))+
    xlab(limni.labels[1])+ 
    ylab(limni.labels[2])+
    coord_cartesian(clip = 'off')+
    geom_hline(yintercept = c(start.y.legend, start.y.legend-0.5, start.y.legend-1 ,start.y.legend - 1.5), color="darkgray", linetype="dashed", size = 0.5)+
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,panel.grid.major=element_blank(),panel.grid.minor=element_blank()
          ,legend.text=element_text(size=20),legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none"
          , plot.margin = margin(t=1, b=1, l=0.1, r=0, "cm"))
  #***********************************************************
  annot <- ggplot() +
    scale_x_continuous(name=element_blank(),  limits = c(0,2))+
    scale_y_continuous(name=element_blank(), expand = c(0,0), limits = c(start.y.legend-1.5, grid_limni.ylim[2]))+
    annotate(geom="text", x=0, y=start.y.legend-0.25, size= 7, label="1. Segmentation of gaugings",   color="blue", hjust = 0)+
    annotate(geom="text", x=0, y=start.y.legend-0.75, size= 7, label="2. After adjustment",   color="blue",  hjust = 0) +
    annotate(geom="text", x=0, y=start.y.legend-1.25, size= 7, label="Official dates of RC update",  color="black",  hjust = 0)+
    #annotate(geom="text", x=1, y=-0.75, size= 7, label="Different used methods:",   color="black")+
    #theme_bw(base_size = 10)+theme(panel.grid = element_blank())+
    geom_hline(yintercept = c(start.y.legend, start.y.legend-0.5,start.y.legend - 1 , start.y.legend-1.5, start.y.legend-2), color="gray", size = 0.5)+
    geom_text() +
    theme_void()
  tt.plot = plot_grid( t.plot, annot, ncol = 2, rel_widths = c(1, 0.32), align = 'h')   #needs cowplot "package"
  ggsave(tt.plot, filename =paste0(dir,"/time_series_with_bt.png"), 
         device = "png", width = 16, height =8, dpi = 400, units = "in")
}




























# Plot stage h time series  with shift times ( + official times) with results of gaugings segmentation:
########################################################################################################
plot.time.shifts.gaugings.4methods <- function(dir, 
                                               results.all.segmentations,
                                               color.periods, 
                                               df.limni,
                                               limni.labels, 
                                               grid.stage,
                                               limits.time,
                                               shift.dates,
                                               initial.time,
                                               index.winner.segmentation,
                                               index.other.segmentation,
                                               order.of.strategies,
                                               titles.legend
                                               ) {
########################################################################################################
  grid_stage  = define.grid(gaugings         = results.all.segmentations[[1]]$gaugings,
                            sequence.grid    = grid.stage,
                            df.limni         = df.limni)
  
  
  #------------------------------------------------
  #PLOT 1:   STAGE RECORD WITH GAUGINGS PER PERIOD
  #------------------------------------------------
  #Utot <- TeX("$Posterior \\; pdf \\; of \\; the \\; change \\; point \\; times \\tau $")
  leg.gaugings           = "Gaugings"
  leg.stage.record       = "Stage record" # (07/11/2001 - 29/10/2018)"
  leg.col.gaugings       = c( "Gaugings"    = "black")
  leg.col.stage.gauging  = c("Stage record" = "gray70", 
                             "Gaugings"     = "black")
  
  t.plot <- ggplot()
  if (is.null(df.limni)==FALSE) {
    t.plot= t.plot + 
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni, color = leg.stage.record), size = 0.3) +
      #scale_x_continuous(name=limni.labels[1], expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
      scale_x_continuous(name = limni.labels[1], expand = c(0,0), limits = limits.time)
  } else {
    t.plot= t.plot + 
      geom_line(data = df.limni, aes(x = -9999, y = 9999, color = leg.stage.record), size = 0.3) +
      scale_x_continuous(name = limni.labels[1], expand = c(0,0), 
                         limits =c(0,tail(results.all.segmentations[[index.winner.segmentation]]$gaugings$time, 1)))
  }
  
  t.plot = t.plot + 
    geom_point(data = results.all.segmentations[[index.winner.segmentation]]$gaugings.segment, 
               aes(x = X.tP. , y= X.h., color = leg.gaugings), size = 5, pch =21, 
               fill= results.all.segmentations[[index.winner.segmentation]]$gaugings.segment$color) +
    #geom_vline(aes(xintercept = data.annotate.gaug.adjust.1$t.adj), color = "red", lwd =0.3, linetype = "dashed")+
    scale_y_continuous(name=limni.labels[2], expand = c(0,0), limits = c(grid_stage$start.y.legend, 
                                                                         grid_stage$grid_limni.ylim.sim[2]), 
                       breaks=seq(grid_stage$grid_limni.ylim.sim[1], 
                                  grid_stage$grid_limni.ylim.sim[2], 
                                  grid_stage$grid_limni.ylim.sim[3]))+
    ylab(limni.labels[2])+
    coord_cartesian(clip = 'off')+
    scale_colour_manual(name   = element_blank(), 
                        values = leg.col.stage.gauging ,
                        breaks = c(leg.stage.record, leg.gaugings),
                        labels = c(leg.stage.record, leg.gaugings),
                        guide  = guide_legend(override.aes = list(
                                 linetype = c("solid", "blank"),
                                 shape    = c(NA, 21),
                                 fill     = c(NA, "black")))) +
    theme_light(base_size=20)+
    theme(  text = element_text(size=20),
            axis.text = element_text(size=15),
            plot.background = element_rect(fill ="transparent", color = NA),
            #panel.grid.major=element_line(size=0.4, linetype = "dashed"), 
            panel.grid.minor= element_blank(),
            panel.grid.major = element_blank(),
            panel.background = element_rect(fill ="transparent"), 
            #axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            plot.margin=unit(c(0.5, 0.5, 2, 0),"cm"),
            axis.title.y = element_text(margin = margin(t = 0, r = 19, b = 0, l = 0))
            ,legend.key.size = unit(1, "cm")
            ,legend.text = element_text(margin=margin(1, 0.3, 0,0), size=20)
            #,legend.spacing.y = unit(-0.5, "cm")
            ,legend.position ="left"
            ,legend.box = "horizontal"
            ,legend.direction = "horizontal"
            ,legend.box.margin=margin(0, 0, 0, 0)
            ,legend.spacing.x = unit(1, 'cm')
            ,legend.justification = c(0,1)
            ,legend.key = element_rect(colour = "transparent", fill = "transparent")
            ,legend.background = element_rect(colour = "transparent", fill = "transparent"))
  
  
  
  if (!is.null(results.all.segmentations[[index.winner.segmentation]]$data.annotate.segm.gaug$t.adj)) {
    if (shift.dates == TRUE) {
      ind = seq(1:length(results.all.segmentations[[index.winner.segmentation]]$data.annotate.segm.gaug$t.adj))
      t.plot <- t.plot +
        geom_segment(aes(x     = results.all.segmentations[[index.winner.segmentation]]$data.annotate.segm.gaug$t.adj[ind[lapply(ind, "%%", 2) != 0]],
                         xend  = results.all.segmentations[[index.winner.segmentation]]$data.annotate.segm.gaug$t.adj[ind[lapply(ind, "%%", 2) != 0]],
                         y     = grid_stage$grid_limni.ylim.sim[1], 
                         yend  = grid_stage$grid_limni.ylim.sim[1]  + 0.1), 
                     color = "red", 
                     size  = 1)+
        geom_segment(aes(x     = results.all.segmentations[[index.winner.segmentation]]$data.annotate.segm.gaug$t.adj[ind[lapply(ind, "%%", 2) == 0]],
                         xend  = results.all.segmentations[[index.winner.segmentation]]$data.annotate.segm.gaug$t.adj[ind[lapply(ind, "%%", 2) == 0]], 
                         y     = grid_stage$grid_limni.ylim.sim[1] + 0.1, 
                         yend  = grid_stage$grid_limni.ylim.sim[1]), 
                     color = "red", 
                     size  = 1)+
        annotate("text", 
                 x      = results.all.segmentations[[index.winner.segmentation]]$data.annotate.segm.gaug$t.adj[ind[lapply(ind, "%%", 2) != 0]],
                 y      = grid_stage$grid_limni.ylim.sim[1]  +  0.2, 
                 label  = results.all.segmentations[[index.winner.segmentation]]$real.dates.shifts[ind[lapply(ind, "%%", 2) != 0]], 
                 color  = "red", 
                 size   = 4) +
        geom_text(aes( 
          x      = results.all.segmentations[[index.winner.segmentation]]$data.annotate.segm.gaug$t.adj[ind[lapply(ind, "%%", 2) == 0]], 
          y      = grid_stage$grid_limni.ylim.sim[1] + 0.2, 
          label  = results.all.segmentations[[index.winner.segmentation]]$real.dates.shifts[ind[lapply(ind, "%%", 2) == 0]]), 
          color  = "red", 
          size   = 4)
    }
  }
  t.plot.withoutlegend     = t.plot  + theme(legend.position="none")
  ###########################################################################################
  
  
  
  
  
  
  # OTHER PLOTS BELOW:
  X        = list(); 
  X.new    = list(); 
  t.plot.2 = list(); 
  t.plot2.withoutlegend = list();
  
  #############################################################
  for (index.strategy in 1:length(results.all.segmentations)) {
  #############################################################
      if ( index.strategy != index.other.segmentation){
      ####################################################################################################  
          if (!is.null(results.all.segmentations[[index.strategy]]$data.annotate.segm.gaug$MAP)) {
               X[[index.strategy]] = results.all.segmentations[[index.strategy]]$pdf.detected.shift.times
               if (ncol(results.all.segmentations[[index.strategy]]$pdf.detected.shift.times)==1) {
                   X.new[[index.strategy]] = X[[index.strategy]]
               } else {
                   X.new[[index.strategy]] = data.frame(time = X[[index.strategy]][,1], 
                                                        ord  = rep(1, length(X[[index.strategy]][,1])))
                   for (orderr in 1:(ncol(results.all.segmentations[[index.strategy]]$pdf.detected.shift.times))) {
                        X.new[[index.strategy]] = rbind(X.new[[index.strategy]], 
                                                        data.frame(time = X[[index.strategy]][,orderr], 
                                                                   ord  = rep(orderr, length(X[[index.strategy]][,1]))))
                   }
               } 
               # assign(paste0("X",index.strategy), X[[index.strategy]])
               # assign(paste0("X.new",index.strategy), X.new[[index.strategy]])
               
          }


          #-------------------------------------------
          #PLOT SEGMENTATION RESULTS OF EACH STRATEGY
          #------------------------------------------
           Utau            = "Posterior pdf of tau"
           leg.col.uncert  = c("Posterior pdf of tau" = "blue")
           leg.shift.times = "Shift times, s"
           leg.col.shifts  = c("Shift times, s" = "solid")
               
          #if (exists("results.all.segmentations[[index.strategy]]$data.annotate.segm.gaug")) {
            t.plot.2[[index.strategy]]   = ggplot() 
            if (ncol(results.all.segmentations[[index.strategy]]$pdf.detected.shift.times)==1) {
              t.plot.2[[index.strategy]] = t.plot.2[[index.strategy]] +
                geom_density(data   = X[[index.strategy]],
                             aes(x  = V1, ..scaled.. ,
                             fill   = Utau),
                             colour = NA,
                             alpha  = 0.3)
            } else {
              
              t.plot.2[[index.strategy]] = t.plot.2[[index.strategy]] +
                geom_density(data       = X.new[[index.strategy]],
                             aes(x      = time, 
                                 group  = ord, 
                             ..scaled..,
                             fill   = Utau), 
                             colour = NA, 
                             alpha  = 0.3)
            }
            t.plot.2[[index.strategy]] = t.plot.2[[index.strategy]] + 
              # annotate("rect", xmin= data.annotate.gaug.1$q2, xmax=data.annotate.gaug.1$q97, 
              #          ymin=0, ymax=0.5, fill="blue", alpha=0.1) +
              # geom_segment(data = data.annotate.gaug.1, aes(x = MAP, y = 0, yend =1, xend= MAP), size = 0.8, color ="blue")+
              geom_segment(data = results.all.segmentations[[index.strategy]]$data.annotate.segm.gaug, 
                           aes(x = t.adj, y = -1, yend =0, xend=t.adj),
                           size = 0.9, color ="red") +
              geom_vline(aes(xintercept = -9999, linetype =leg.shift.times), col="red", lwd =0.9, show_legend = T) +
              scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
              xlab(limni.labels[1])+ 
              coord_cartesian(clip = 'off')+
              geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+     
              scale_fill_manual(name=element_blank(), 
                                breaks=Utau,
                                labels = unname(TeX("$Posterior \\; pdf \\; of \\; \\tau $")),
                                values = leg.col.uncert) +
              scale_linetype_manual(name=element_blank(),
                                    values=leg.col.shifts,
                                    breaks=leg.shift.times,
                                    guide = guide_legend(override.aes = list(
                                      col ="red"))) +
              theme_set(theme_gray(base_size = 20))+
              #theme_light(base_size = 20) +
              theme(axis.text             = element_text(size=20)
                    ,axis.title           = element_text(size=30, face="bold")
                    ,plot.title           = element_text(face = "bold", size = 22, color= "black", margin=margin(0,0,0,0))
                    ,panel.grid.major     = element_blank()
                    ,panel.grid.minor     = element_blank()
                    ,plot.margin          = unit(c(0,0.5,0.4,0.9),"cm")
                    ,axis.text.x          = element_blank()
                    ,axis.text.y          = element_blank()
                    ,axis.ticks.y         = element_blank()
                    ,axis.ticks.x         = element_blank()
                    ,legend.key.size      = unit(1, "cm")
                    ,legend.spacing.y     = unit(-0.5, "cm")
                    #,legend.spacing.x = unit(0.8, 'cm')
                    ,legend.text          = element_text(margin=margin(1,0.3,0,0), size=20)
                    ,legend.position      = "left"
                    ,legend.box           = "horizontal"
                    ,legend.direction     = "horizontal"
                    ,legend.justification = c(0,1)
                    ,legend.key           = element_rect(colour = "transparent", fill = "transparent")
                    ,legend.background    = element_rect(colour = "transparent", fill = "transparent")) +
              ggtitle(titles.legend[index.strategy])
            if (is.null(df.limni)==FALSE) {
              # t.plot2= t.plot2 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
              t.plot.2[[index.strategy]] = t.plot.2[[index.strategy]] +
                                          scale_x_continuous(name=element_blank(), expand = c(0,0), limits =limits.time)
            } else {
              t.plot.2[[index.strategy]] = t.plot.2[[index.strategy]] + 
                                           scale_x_continuous(name=element_blank(), 
                                                    expand = c(0,0), 
                                                    limits =c(0,tail(results.all.segmentations[[index.strategy]]$gaugings$time,1)))
            }
        
        
      ########
      } else {
      ########
           X.new[[index.strategy]] = list();
           X[[index.strategy]]     =list();
           # assign(paste0("X",index.strategy), X[[index.strategy]])
           # assign(paste0("X.new",index.strategy), X.new[[index.strategy]])
           #------------------------------------------------
           #PLOT: method with other method one single pass) Bin SEG
           #------------------------------------------------
           #if (exists("results.all.segmentations[[index.strategy]]$data.annotate.segm.gaug")) {
             t.plot.2[[index.strategy]] = ggplot() 
             t.plot.2[[index.strategy]] = t.plot.2[[index.strategy]] +  
                       annotate("rect", 
                                xmin  = results.all.segmentations[[index.strategy]]$data.annotate.segm.gaug$t2, 
                                xmax  = results.all.segmentations[[index.strategy]]$data.annotate.segm.gaug$t97 - 
                                        (tail(results.all.segmentations[[index.strategy]]$gaugings$time,1) -
                                           results.all.segmentations[[index.strategy]]$gaugings$time[1])/5000, 
                                ymin  = 0,
                                ymax  = 1, 
                                fill  = "blue", 
                                alpha = 0.3) +
               geom_segment(data = results.all.segmentations[[index.strategy]]$data.annotate.segm.gaug,
                            aes(x = t.adj, y = -1, yend =0, xend = t.adj), 
                            size = 0.9, color ="red")+
               scale_y_continuous(name = "", expand = c(0,0), limits = c(-1,1))+
               coord_cartesian(clip = 'off')+
               geom_hline(yintercept = c(0), color="darkgray", linetype="dashed", size = 0.3)+
               theme_set(theme_gray(base_size = 20))+
               #theme_light(base_size = 20) +
               #theme_classic()+
               theme(axis.text             = element_text(size=20)
                     ,axis.title           = element_text(size=30, face="bold")
                     ,plot.title           = element_text(face = "bold", size = 22, color= "black", margin=margin(0,0,0,0))
                     ,panel.grid.major     = element_blank()
                     ,panel.grid.minor     = element_blank()
                     ,plot.margin          = unit(c(0,0.5,0.4,0.9),"cm")
                     ,axis.text.x          = element_blank()
                     ,axis.text.y          = element_blank()
                     ,axis.ticks.y         = element_blank()
                     ,axis.ticks.x         = element_blank()
                     ,legend.key.size      = unit(1, "cm")
                     ,legend.spacing.y     = unit(-0.5, "cm")
                     #,legend.spacing.x = unit(0.8, 'cm')
                     ,legend.text          = element_text(margin=margin(1,0.3,0,0), size=20)
                     ,legend.position      = "left"
                     ,legend.box           = "horizontal"
                     ,legend.direction     = "horizontal"
                     ,legend.justification = c(0,1)
                     ,legend.key           = element_rect(colour = "transparent", fill = "transparent")
                     ,legend.background    = element_rect(colour = "transparent", fill = "transparent")) +
                ggtitle(titles.legend[index.strategy])
             if (is.null(df.limni)==FALSE) {
               #t.plot4 = t.plot4 + scale_x_continuous(name=element_blank(), expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1)))
               t.plot.2[[index.strategy]] = t.plot.2[[index.strategy]] + 
                                            scale_x_continuous(name = element_blank(), 
                                                              expand = c(0,0), limits =limits.time)
             } else {
               t.plot.2[[index.strategy]] = t.plot.2[[index.strategy]] + 
                                            scale_x_continuous(name=element_blank(), 
                                            expand = c(0,0), 
                                            limits =c(0,tail(results.all.segmentations[[index.strategy]]$gaugings$time,1)))
           }
      }
      t.plot2.withoutlegend[[index.strategy]]    = t.plot.2[[index.strategy]] + theme(legend.position="none")
    
  }
    

  
  
  #-----------------------------------
  #PLOT: OFFICIAL OR KNOWN SHIFT DATES
  #----------------------------------- 
  official.dates         = "Official dates of RC update"
  leg.col.official.dates = c("Official dates of RC update" = "black")
  t.plot.off = ggplot() 
  
  t.plot.off = t.plot.off + 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(0, 1))+
    coord_cartesian(clip = 'off')+
    #geom_hline(yintercept = c(0, 1), color="darkgray", linetype="dashed", size = 0.5)+
    theme_set(theme_grey(base_size = 20))+
    #theme_light(base_size = 20) +
    theme(axis.text             = element_text(size=20)
          ,axis.title           = element_text(size=30,face="bold")
          ,plot.title           = element_text(face = "bold", size = 22, color= "black", margin=margin(0,0,0,0))
          ,panel.grid.major     = element_blank()
          ,panel.grid.minor     = element_blank()
          ,legend.key.size      = unit(1, "cm")
          ,legend.spacing.y     = unit(-0.5, "cm")
          #,legend.spacing.x = unit(0., 'cm')
          ,legend.text          = element_text(margin=margin(1,0.3,0,0), size=20)
          ,legend.position      ="left"
          ,legend.box           = "horizontal"
          ,legend.direction     = "horizontal"
          ,legend.justification = c(0,1)
          ,legend.key           = element_rect(colour = "transparent", fill = "transparent")
          ,legend.background    = element_rect(colour = "transparent", fill = "transparent")
          ,plot.margin          = unit(c(0,0.5,0.5,0.9),"cm")
          ,axis.text.x          = element_blank()
          ,axis.text.y          = element_blank()
          ,axis.ticks.y         = element_blank()
          ,axis.ticks.x         = element_blank()) +
    ggtitle(tail(titles.legend,1))
  
  if (is.null(df.limni)==FALSE) {
      t.plot.off = t.plot.off + 
      scale_x_continuous(name = element_blank(), 
                         expand = c(0,0), limits =limits.time)
  } else {
      t.plot.off = t.plot.off + 
      scale_x_continuous(name=element_blank(), 
                         expand = c(0,0), 
                         limits =c(0,tail(results.all.segmentations[[index.strategy]]$gaugings$time,1)))
  }
  #if (exists("results.all.segmentations[[1]]$data.annotate.known")==TRUE)  {
  if (is.null(results.all.segmentations[[1]]$data.annotate.known)==FALSE) {
      t.plot.off = t.plot.off +
        geom_point(data = results.all.segmentations[[1]]$data.annotate.known, 
                   aes(x = xeffect, y = 0, 
                   color= official.dates), size = 4, shape =4, stroke=2 )+
        scale_color_manual(name = element_blank(), 
                           breaks = official.dates,
                           labels = official.dates,
                           values = leg.col.official.dates,
                           guide  = guide_legend(override.aes = list(
                                    linetype = c("blank"),
                                    shape = 4)))
      # geom_point(data = results.all.segmentations[[1]]$data.annotate.known,
      #            aes(x = xpotent, y = 0), color= "green", size = 4, shape =4, stroke=2 ) # add shift times that are not sure !!!)
  }
  t.plot.off.withoutlegend = t.plot.off    + theme(legend.position="none")
  t.plot2.withoutlegend[[length(results.all.segmentations)+1]] = t.plot.off.withoutlegend
  
  
  ####################################################################################
  # Legend plot:
  Legend.title <- ggdraw() + 
               draw_label("Legend",
               fontface = 'bold',
               x = 0, hjust = 0,
               size = 20) +
               theme(# add margin on the left of the drawing canvas,
                     # so title is aligned with left edge of first plot
               plot.margin = margin(15, 0, 15 , 0))
  glegend.A <- cowplot::get_legend( t.plot)
  glegend.B <- cowplot::get_legend( t.plot.2[[index.winner.segmentation]])
  glegend.C <- cowplot::get_legend( t.plot.off)
  glegends = plot_grid(glegend.A, 
                       glegend.B, 
                       glegend.C, 
                       ncol =3, nrow =1, rel_widths = c(0.8, 0.8, 0.6))
  Legends =  plot_grid(Legend.title, 
                       glegends,
                       ncol = 1,
                       # rel_heights values control vertical title margins
                       rel_heights = c(0.5, 1))
  ###########################################################################
  #save final plot:
  plot.2.all = plot_grid( plotlist    = t.plot2.withoutlegend, 
                          ncol        = 1, 
                          nrow        = length(t.plot2.withoutlegend), 
                          rel_heights = c(rep(1, length(t.plot2.withoutlegend)-1), 0.9)) 
  
  t.plot.FINAL = plot_grid( t.plot.withoutlegend,
                            plot.2.all,
                            #t.plot.off.withoutlegend,
                            Legends,
                            
                            ncol        = 1, 
                            nrow        = 3, 
                            rel_heights = c(1, 0.9, 0.2))

  # ggsave(t.plot.FINAL, filename =paste0(dir,"/time_series_4methods.png"), 
  #        device = "png", width = 16, height =14, dpi = 400, units = "in")
  pdf(paste0(dir,"/Figure5.pdf"), 16, 14 , useDingbats=F)
  print(t.plot.FINAL)
  dev.off()
}

























###############################################################################################################
legend.final <- function(station.name, dir.case_study, Q10.ts, Q90.ts, ts.res, ts.res.before, ts.res.plus,
                         Hmin_grid, Hmax_grid, seg.iter, t_Gaug, h_Gaug,
                         Shift.Q, BIC.df, df.limni, CdT.P,
                         Q10.mu.res, Q90.mu.res, mu.res,
                         plot.leg.Baratin ) {
###############################################################################################################
  # time series limni:
  #******************
  Utot.times = "90% total uncertainty of change point times"
  MAPtimes = "Change point times (MAP)"
  color.segm = c("90% total uncertainty of change point times"="black")
  col.segm = c("Change point times (MAP)"="black")
  
  if (!is.null(ts.res[1])) {
    ts.plt <- ggplot()+
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "darkgray",size = 0.2)+
      geom_point(aes(x = t_Gaug, y =h_Gaug), color = "gray60", size = 1) +
      geom_point(data = CdT.P, aes(x = tP, y = hP), color = CdT.P$color, size = 1) +
      theme_light(base_size = 15)  +
      scale_x_continuous(name="time [days]", expand = c(0,0)) +
      scale_y_continuous(name="Stage h [m]", limits =c(-1,4), expand =c(0,0))+ 
      xlab("Time [days]")+ ylab("Stage [m]") +
      geom_rect( mapping= aes(xmin= Q10.ts, xmax=Q90.ts ,ymin=-Inf, ymax=Inf, fill=Utot.times), 
                 alpha=0.1 ) +
      geom_vline(aes(xintercept = ts.res, col=MAPtimes), lwd =0.5, linetype = "longdash") +
      scale_fill_manual(name=element_blank(), values=color.segm) +
      scale_colour_manual(name = element_blank(), values=col.segm,
                          breaks=c(MAPtimes),
                          labels = c(MAPtimes),
                          guide = guide_legend(override.aes = list(
                                  linetype = c("longdash"),
                                  shape = c(NA)))) +
      theme(#text = element_text(size=10, family="LM Roman 10"),
            plot.title = element_text(hjust = 0.5),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent"),
            plot.margin=unit(c(0.3,0.5,0.05,1.05),"cm"),
            axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            #****** legend *****************
            legend.background = element_rect(colour = "transparent", fill = "transparent"),
            legend.justification=c(0,1),
            legend.margin=unit(0,"cm"),
            legend.box= "vertical",
            legend.box.just = "left",
            legend.key.size=unit(1,"lines"),
            legend.text.align=0,
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.text=element_text(size=10)
            )
            # legend.position ="none",                                                    
            # legend.key.size = unit(0.5, "cm"),
            # legend.text=element_text(size=8),
            # legend.spacing.x = unit(0.05, 'cm'),
            # legend.spacing.y = unit(0.05, 'cm'),
            # legend.key.height=unit(0.5,"line"),
            # legend.key.width=unit(0.8,"line")
  } else {
    ts.plt <- ggplot()+
      geom_line(data = df.limni, aes(x = t_limni, y = h_limni), color = "darkgrey",size = 0.2)+
      geom_point(data = CdT.P, aes(x = tP, y = hP), colour = "black",size = 1.5) +
      theme_light(base_size = 15) +  
      scale_x_continuous(name="time [day]",expand = c(0,0)) +
      scale_y_continuous(name="Stage h [m]",limits = c(Hmin_grid,Hmax_grid),expand = c(0,2)) +
      xlab("Time [day]")+ ylab("Stage h [m]")+
      theme(plot.background = element_rect(fill ="transparent", color = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill ="transparent"), 
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            plot.margin=unit(c(0.2,0.5,0.05,0.05),"cm"),
            #****************************************************************************
            legend.key = element_rect(colour = "transparent", fill = "transparent"),
            legend.background = element_rect(colour = "transparent", fill = "transparent"),
            legend.position ="bottom", 
            legend.box = "vertical",
            legend.direction="vertical",
            legend.text=element_text(size=10))
  }
  glegend.1 <- cowplot::get_legend(ts.plt)
  
  
  
  # time series residuals:
  #***********************
  Utot.mean = "90% total uncertainty of segment mean"
  MAPmean = "Segment mean (MAP)"
  resid ="Gaugings residuals with 95% error bars"
  #
  color.segm.2 = c("90% total uncertainty of change point times"= "black",
                 "90% total uncertainty of segment mean" = "black")
  col.segm.2 = c("Change point times (MAP)"= "black", 
               "Segment mean (MAP)"= "black",
               "Gaugings residuals with 95% error bars"= "black")
  
  if (!is.null(ts.res[1])) {
    seg.plt <- ggplot()+
      geom_point(data = Shift.Q, aes(x = Shift.Q$alpha_t, y = Shift.Q$alpha, col=resid), 
                 shape = 1, size = 1.4) +
      geom_errorbar(data = Shift.Q, aes(x=Shift.Q$alpha_t, 
                                        ymin= (Shift.Q$alpha - Shift.Q$sigma.tot), 
                                        ymax = (Shift.Q$alpha + Shift.Q$sigma.tot), color = resid),
                    size = 0.25, width=90) +
      theme_light(base_size = 15)+
      scale_x_continuous(name = "time [days]", expand = c(0,0), 
                         limits =c(Shift.Q$alpha_t[1], tail(df.limni$t_limni,1))) +
      scale_y_continuous(name = "residual",
                         #name = TeX("$\\epsilon = Q_G - Q_{RC} (h_G) [m]$"),
                         limits =c(min(Shift.Q$alpha - Shift.Q$sigma.tot),
                                   max(Shift.Q$alpha + Shift.Q$sigma.tot)))+ 
      geom_rect( mapping= aes(xmin= Q10.ts, xmax=Q90.ts ,ymin=-Inf, ymax=Inf, fill= Utot.times), 
                 alpha=0.1) +
      geom_vline(aes(xintercept = ts.res, col=MAPtimes), lwd =0.5, linetype = "longdash") +
      geom_rect(mapping = aes(xmin= ts.res.before, xmax=ts.res.plus, ymin=Q10.mu.res,
                              ymax=Q90.mu.res, fill=Utot.mean), alpha=0.3) + 
      geom_segment(mapping=aes(x =ts.res.before , y = mu.res, xend = ts.res.plus, yend = mu.res, col =MAPmean))+
      scale_fill_manual(name=element_blank(), values=color.segm.2) +
      scale_colour_manual(name=element_blank(), 
                          values=col.segm.2,
                          breaks=c(MAPtimes, MAPmean, resid),
                          labels=c(MAPtimes, MAPmean, resid),
                          guide=guide_legend(override.aes = list(
                          linetype=c("longdash","solid", "blank"),
                          shape=c(NA, NA, 1)))) +
      theme(
            text = element_text(size=10, family="LM Roman 10"),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent"),
            plot.margin=unit(c(0.3,0.5,0.05,1.05),"cm"),
            #*****************************************************************************
            legend.background = element_rect(colour = "transparent", fill = "transparent"),
            legend.justification=c(0,1),
            legend.margin=unit(0,"cm"),
            legend.box= "vertical",
            legend.box.just = "left",
            legend.key.size=unit(1,"lines"),
            legend.text.align=0,
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.text=element_text(size=10)
            )
            #legend.position ="none",
            # legend.key.size = unit(0.5, "cm"),
            # legend.text=element_text(size=8),
            # legend.spacing.x = unit(0.05, 'cm'),
            # legend.spacing.y = unit(0.05, 'cm'),
            # legend.key.height=unit(0.5,"line"),
            # legend.key.width=unit(0.8,"line")
    #*****
  } else {
    #*****
    seg.plt <- ggplot()+
      geom_point(data = Shift.Q, aes(x = Shift.Q$alpha_t, y = Shift.Q$alpha, color=Shift.Q$Q_Gaug), size = 1.4) +
      geom_errorbar(data = Shift.Q, aes(x=Shift.Q$alpha_t, ymin= (Shift.Q$alpha - Shift.Q$sigma.tot), 
                                        ymax = (Shift.Q$alpha + Shift.Q$sigma.tot), color=Shift.Q$Q_Gaug ), 
                    size = 0.3, width=40) +
      labs(color = "Discharge [m3.s-1]") +
      scale_color_gradient(low="black", high="green") +
      theme_light(base_size = 10) +
      scale_x_continuous(name = "time [days]", expand = c(0,0), 
                         limits =c(Shift.Q$alpha_t[1],tail(Shift.Q$alpha_t,1)+50)) +
      scale_y_continuous(name = TeX("$\\epsilon = Q_G - Q_{RC} (h_G) [m]$"),
                         limits =c(min(Shift.Q$alpha - Shift.Q$sigma.tot),
                                   max(Shift.Q$alpha + Shift.Q$sigma.tot)))+ 
      annotate("rect",xmin= Shift.Q$alpha_t[1], xmax= tail(Shift.Q$alpha_t,1), ymin=Q10.mu.res,
               ymax=Q90.mu.res, fill="red", alpha=0.3) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.background = element_rect(fill ="transparent", color = NA),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill ="transparent"), 
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            plot.margin=unit(c(0.2,0.5,0.05,0.05),"cm"),
            
            legend.background = element_rect(colour = "transparent", fill = "transparent"),
            legend.justification=c(0,1),
            legend.margin=unit(0,"cm"),
            legend.box= "vertical",
            legend.box.just = "left",
            legend.key.size=unit(1,"lines"),
            legend.text.align=0,
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.text=element_text(size=10)
            )
  } 
  
  glegend.2  <- cowplot::get_legend(seg.plt)
  
  
  
  #BIC  AIC
  #******************************************************
  grid = seq(-10000,10000, 50)
  min_grid <- NULL
  for (i in 1:length(grid)) {
    if((grid[i+1] > min(BIC.df$AIC, BIC.df$BIC)) & (grid[i] < min(BIC.df$AIC, BIC.df$BIC)) ) {
      min_grid = grid[i]
    }
  }
  max_grid <- NULL
  for (i in 1:length(grid)) {
    if((grid[i+1] >= max(BIC.df$AIC, BIC.df$BIC)) & (grid[i] <= max(BIC.df$AIC, BIC.df$BIC)) ) {
      max_grid = grid[i+1]
    }
  }
  col =c("AIC (Aikake, 1973)" ="lightblue", "BIC (Schwarz, 1978)" ="gray40")
  AICcol ="AIC (Aikake, 1973)"; BICcol ="BIC (Schwarz, 1978)";
  
  bic.plot <- ggplot(BIC.df) + 
      theme_light(base_size = 15)+
      geom_line(aes(x = x, y = BIC, color =BICcol), size = 0.5) +
      geom_line(aes(x = x, y = AIC, color =AICcol), size = 0.5) +
      scale_x_continuous(name = "Number of segments nS", expand = c(0,0), 
                         limits =c(1, length(BIC.df$BIC)))+
      scale_y_continuous(name = "BIC and AIC" , expand = c(0,0),
                         limits = c(min_grid, max_grid),
                         breaks = seq(min_grid, max_grid, 50)) +
      # breaks = scales::pretty_breaks(n = 2)) + 
      # theme_bw()+
      scale_colour_manual(name = element_blank(), values=col, breaks = c(AICcol, BICcol),
                          guide = guide_legend(override.aes = list(
                                  linetype = c("solid", "solid"),
                                  shape = c(0, 2)))) +
      xlab("Number of segments nS") + 
      ylab("BIC & AIC")+ coord_cartesian(clip = 'off')+
      geom_point(aes(x = x, y = BIC), shape =0, size = 4, color = "gray20") +
      geom_point(aes(x = x, y = AIC), shape =2, size = 4, color = "blue") +
      geom_point(aes(x = BICmin, y = BIC[BICmin]), shape = 15, size = 4, color = "gray20") +
      geom_point(aes(x = AICmin, y = AIC[AICmin]), shape = 17, size = 4, color = "blue") +
      theme(
           plot.title = element_text(hjust = 0.5), 
           plot.background = element_rect(fill ="transparent", color = NA),
           panel.grid.major=element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
           panel.background = element_rect(fill ="transparent"), 
           axis.ticks = element_line(colour = "black"),
           plot.margin=unit(c(0.3,0.5,0.05,1.05),"cm"),
           text = element_text(size=10, family="LM Roman 10"),
           
           #legend.key = element_rect(colour = "transparent", fill = "transparent"),
           legend.background = element_rect(colour = "transparent", fill = "transparent"),
           legend.justification=c(0,1),
           legend.margin=unit(0,"cm"),
           legend.box= "vertical",
           legend.box.just = "left",
           legend.key.size=unit(1,"lines"),
           legend.text.align=0,
           legend.key = element_blank(),
           legend.title = element_blank(),
           legend.text=element_text(size=10)
           )
  
    glegend.3  <- cowplot::get_legend(bic.plot)
    
    
    #******************************************************************************************************************
    legend.all <- plot_grid(plot.leg.Baratin,  glegend.1, glegend.2, glegend.3,
              labels = c("Legend ", " ", " ", " "),
              ncol = 4, nrow = 1) # rel_heights = c(3, 3, 3, 3))
    ggsave(legend.all, filename = paste0(dir.case_study,"/Results/segmentation_gaugings/it",seg.iter,"/all_legend.png"),
           bg = "transparent", width = 15, height =3, dpi = 400)
    #*******************************************************************************************************************
    
    return(list(glegend.1, glegend.2, glegend.3)) 
}



























#################################################################################################
segm_other_method_plot <- function(dir.seg.gaug, Shift.Q, Q10.ts, Q90.ts, ts.res, ts.res.before,
                                   ts.res.plus, Q10.mu.res, Q90.mu.res, 
                                   mu.res, seg.iter, df.limni, 
                                   ts.real, resid.uncertaint) {
##################################################################################################
  if (!is.null(ts.res[1])) {
      seg.plt <- ggplot()+
      geom_point(data = Shift.Q, aes(x = Shift.Q$alpha_t, y = Shift.Q$alpha, col=resid),
                 shape = 1, size = 1.4)
      if (resid.uncertaint==TRUE) {
        seg.plt <-  seg.plt+
        scale_y_continuous(name = TeX("Residual $r = Q_i - Q_{RC} (h_i) $"),
                             limits =c(min(Shift.Q$alpha - 2*Shift.Q$sigma.tot),
                                       max(Shift.Q$alpha + 2*Shift.Q$sigma.tot))) + 
        geom_errorbar(data = Shift.Q, 
                      aes(x=Shift.Q$alpha_t, 
                          ymin= (Shift.Q$alpha - 2*Shift.Q$sigma.tot), 
                          ymax = (Shift.Q$alpha + 2*Shift.Q$sigma.tot), color = resid),
                      size = 0.1, width=0.05*tail(Shift.Q$alpha_t,1))
      } else {
        seg.plt <-  seg.plt+
          scale_y_continuous(name = TeX("Residual $r = Q_i - Q_{RC} (h_i) $"))
      }
      seg.plt <-  seg.plt +
      coord_cartesian(clip = 'off') + 
      geom_rect( mapping= aes(xmin= Q10.ts, xmax=Q90.ts ,ymin=-Inf, ymax=Inf, fill=Utot.times), 
                 alpha=0.1 ) +
      geom_vline(aes(xintercept = ts.res, col=MAPtimes), lwd =0.5, linetype = "solid") +
      geom_vline(aes(xintercept = ts.real, col=realtimes), lwd =0.3, linetype ="longdash") +
      geom_rect(mapping = aes(xmin= ts.res.before, xmax=ts.res.plus, ymin=Q10.mu.res,
                              ymax=Q90.mu.res, fill=Utot.mean), alpha=0.3) + 
      geom_segment(mapping=aes(x =ts.res.before , y = mu.res, xend = ts.res.plus, yend = mu.res, col =MAPmean))+
      scale_fill_manual(name=element_blank(), values=color.segm) +
      scale_colour_manual(name = element_blank(), values=col.segm,
                          breaks=c(MAPtimes, MAPmean, resid),
                          labels = c(MAPtimes, MAPmean, resid),
                          guide = guide_legend(override.aes = list(
                            linetype = c("solid","longdash","solid", "blank"),
                            shape = c(NA,NA, NA, 1)))) +
      theme_light(base_size = 15) +
      theme(text = element_text(size=14),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent"),
            plot.margin=unit(c(0.3,0.5,0.05,1.05),"cm"),
            legend.position ="none",
            legend.key = element_rect(colour = "transparent", fill = "transparent"),
            legend.background = element_rect(colour = "transparent", fill = "transparent"))
    if (is.null(df.limni)==FALSE) {
      seg.plt = seg.plt + scale_x_continuous(name = "time [days]", expand = c(0,0))
    } else {
      seg.plt = seg.plt + scale_x_continuous(name = "time [days]", expand = c(0,0))
    }    
    ggsave(seg.plt, filename =paste0(dir.seg.gaug,"/segm_it",seg.iter,".png"),width=4,height =3,dpi=300)
  
    #****************
  } else {
    #****************
    seg.plt <- ggplot()+
      geom_point(data = Shift.Q, aes(x = Shift.Q$alpha_t, y = Shift.Q$alpha, col=resid),shape = 1, size = 1.4)
    if (resid.uncertaint==TRUE) {
      seg.plt <-  seg.plt+geom_errorbar(data = Shift.Q, aes(x=Shift.Q$alpha_t, 
                                                            ymin= (Shift.Q$alpha - 2*Shift.Q$sigma.tot), 
                                                            ymax = (Shift.Q$alpha + 2*Shift.Q$sigma.tot), color = resid),
                                        size = 0.1, width=0.05*tail(Shift.Q$alpha_t,1))
    }
    seg.plt <-  seg.plt+
      scale_y_continuous(name = TeX("Residual $r = Q_i - Q_{RC} (h_i)$"),
                         limits =c(min(Shift.Q$alpha - 2*Shift.Q$sigma.tot),
                                   max(Shift.Q$alpha + 2*Shift.Q$sigma.tot))) +
      coord_cartesian(clip = 'off') + 
      geom_rect(mapping = aes(xmin= Shift.Q$alpha_t[1], xmax=tail(Shift.Q$alpha_t,1), ymin=Q10.mu.res,
                              ymax=Q90.mu.res, fill=Utot.mean), alpha=0.3) + 
      geom_segment(mapping=aes(x =Shift.Q$alpha_t[1], y =mu.res, xend = tail(Shift.Q$alpha_t,1), yend = mu.res, col =MAPmean)) +
      
      scale_fill_manual(name=element_blank(), values=color.segm_1) +
      scale_colour_manual(name = element_blank(), values=col.segm_1,
                          breaks=c(MAPmean, resid),
                          labels = c(MAPmean, resid),
                          guide = guide_legend(override.aes = list(
                            linetype = c("solid", "blank"),
                            shape = c(NA, 1)))) +
      theme_light(base_size = 15)+
      theme(text = element_text(size=14),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent"),
            plot.margin=unit(c(0.3,0.5,0.05,1.05),"cm"),
            legend.position ="none",
            legend.key = element_rect(colour = "transparent", fill = "transparent"),
            legend.background = element_rect(colour = "transparent", fill = "transparent"))
    if (is.null(df.limni)==FALSE) {
        seg.plt = seg.plt + scale_x_continuous(name = "time [days]", expand = c(0,0))
    } else {
        seg.plt = seg.plt + scale_x_continuous(name = "time [days]", expand = c(0,0))
    }    
    ggsave(seg.plt, filename =paste0(dir.seg.gaug,"/segm_it",seg.iter,".png"),width=4,height =3,dpi=300)
  } 
  return(seg.plt)
}








# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Figure for paper Darienzo et al., 2020, with the main results of some iterations for the case study of Meyras
# ###############################################################################################################
# plot.figure.article = function(dir, iteration.indexes,
#                                df.limni,
#                                grid_RC.ylim, grid_RC.xlim, grid_RC.xstep, grid_RC.ystep,
#                                ticks_RC.y.log, RC.x.labels, RC.y.labels, u.m.Qgaug, u.m.Hgaug,
#                                ncontrol) {
# ###############################################################################################################  
#      # For each chosen iteration plot  3 subplots: 
#      # A) RC of reference, 
#           plot.A =list()
#           plot.A.legend =list()
#      # B) Criteria evolution with increasing nS
#           plot.B =list()
#           plot.B.legend =list()
#      # C) Results of the segmentation of residuals
#           plot.C =list()
#           plot.C.legend =list()
#           iteration.plots =list()
#           iteration.plots2 = list()
#   
#      for (i in 1:length(iteration.indexes)) {
#           # Let's study the iteration iteration.indexes[i] :
#           
#           ##########
#           # PLOT A :
#           ##########
#           dir.plot.A = paste0(dir, "/it",iteration.indexes[i],"/BaRatin")
#           setwd(dir.plot.A)
#           gaug.tot = read.table(paste0(dir.plot.A,"/df_gaug_tot",iteration.indexes[i],".txt"), 
#                                 sep ="\t", header =TRUE)
#           gaug.it = read.table(paste0(dir.plot.A,"/df_gaug_it",iteration.indexes[i],".txt"), 
#                                sep ="\t", header =TRUE)
#           file.X=paste0(dir.plot.A,"/Hgrid.txt")
#           file.env.tot=paste0(dir.plot.A,"/Qrc_TotalU.env")
#           file.spag.tot=paste0(dir.plot.A,"/Qrc_TotalU.spag")
#           file.env.par=paste0(dir.plot.A,"/Qrc_ParamU.env")
#           file.spag.par=paste0(dir.plot.A,"/Qrc_ParamU.spag")
#           file.spag.max=paste0(dir.plot.A,"/Qrc_Maxpost.spag")
#           file.summary = paste0(dir.plot.A,"/Results_Summary.txt")
#           model=ReadModel(ModelFile=paste0(dir.plot.A,"/Config_Model.txt"))
#           # legend:
#           MAPcurve <-  "MAP rating curve"
#           Ustruct.MAP <- "95% remnant uncertainty of the MAP RC"
#           Utot <- "95% total uncertainty"
#           gaugings.period <- "Gaugings of the current period"
#           gaugings.other <- "Gaugings of the other periods"
#           activ.stage <- "Control activation stages (MAP)"
#           col.uncert <- c( 
#                           "95% remnant uncertainty of the MAP RC" = "red") #"salmon") #"indianred1") 
#           cols <- c("Gaugings of the current period"="black", 
#                     "Gaugings of the other periods"="gray90",
#                     "MAP rating curve"="black")
#           #cols <- c("Current periods gaugings"= "black", "Other periods gaugings" = "gray", 
#           #          "MAP rating curve"= "blue") 
#           #al <- c("U_tot"=0.2, "U_par"=0.4, "U_MAP"=0.3)    #"U_par"="indianred1",
#           maxpost <- read.table(file = file.spag.max)
#           maxpost.par <-read.table(file = file.summary)
#           Hg <- read.table(file = file.X)
#           ktransit=NULL
#           for (n in 1:ncontrol) {
#                 ktransit[n] = maxpost.par[16, ncontrol*3+2+n]    
#           }
#           gamma1 <- maxpost.par[16,ncontrol*3+1];
#           str.err = 0
#           if (remnant.err.model == "Linear") {  
#                gamma2 <- maxpost.par[16,ncontrol*3+2]
#                str.err = gamma1+gamma2*maxpost
#           } else {
#                str.err = gamma1
#           }
#           env.tot = EnvelopLayer(file=file.env.tot, Xfile=file.X, color = Utot, alpha = 1)
#           #spag=SpaghettiLayer(file=file.spag, Xfile=file.X)
#           env.par = EnvelopLayer(file=file.env.par, Xfile=file.X , color = "U_par", alpha =1)
#           #env.max = EnvelopLayer(file=file.spag.max, Xfile=file.X , color = "black", alpha =1 )
#           maxpost.data <- data.frame(x = Hg[which(Hg >= ktransit[1]),1], 
#                                      y = maxpost[which(Hg >= ktransit[1]),1], 
#                                      y2 = str.err[which(Hg >= ktransit[1]),1])
#           new_line <- element_line(color = "grey", size = 0.1, linetype = 2)
#           loadfonts()
#           
#           plot.A[[i]] = ggplot() + 
#                         coord_trans(limx = grid_RC.xlim, limy = grid_RC.ylim) +
#                         ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
#                         xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]"))) + 
#                         scale_y_continuous(breaks=seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep))+
#                         scale_x_continuous(breaks=seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep)) +
#                         #annotation_logticks(base = 10, sides = "l", scaled = TRUE, colour = "black", size = 0.3, linetype = 1)+
#                         # scale_y_log10(na.value=-10, breaks=ticks_RC.y.log, labels=ticks_RC.y.log) +
#                         #               scale_x_continuous(breaks=seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep)) +
#                         geom_ribbon(data = na.omit(maxpost.data), aes(x = x, ymin =(y-y2) , 
#                                                                       ymax = (y+y2), 
#                                                                       fill = Ustruct.MAP ), alpha = 0.2)+
#                         geom_errorbar(data= gaug.tot, aes(x=h_Gaug, ymin =Q_Gaug-2*uQ_Gaug, 
#                                                            ymax =Q_Gaug+2*uQ_Gaug, color =gaugings.other), 
#                                                            width=.1, size = 0.2)+
#                         geom_point(data= gaug.tot,    aes(x = h_Gaug, y = Q_Gaug, 
#                                                           color =gaugings.other), size = 1.5) +
#                         geom_errorbar(data= gaug.it,  aes(x=hP, ymin =QP-2*uQP, 
#                                                           ymax =QP+2*uQP, 
#                                                           color =gaugings.period), 
#                                                       width=.1, size = 0.2)+
#                         geom_point(data= gaug.it,     aes(x = hP, y = QP, 
#                                                           color =gaugings.period), size = 1.5) +
#                         geom_line(data = na.omit(maxpost.data), 
#                                   aes(x = x, y= y, color = MAPcurve), 
#                                   size = 0.8, na.rm = TRUE)+
#                         #geom_vline(aes(xintercept = ktransit, color =activ.stage), size = 0.3, linetype ="dashed")+   
#                         
#                         scale_colour_manual(name = element_blank(), 
#                                             values=c(cols),
#                                             breaks=c(gaugings.period,
#                                                      gaugings.other,
#                                                      MAPcurve),
#                                             #labels = c(gaugings.period,gaugings.other,MAPcurve,activ.stage),
#                                             guide = guide_legend(override.aes = list(
#                                                     linetype = c("blank", "blank","solid"),
#                                                     shape = c(19, 19, NA)))) +
#                         scale_fill_manual(name=element_blank(), 
#                                           values= col.uncert) +
#                         theme_light(base_size = 20)+
#                         theme(  text = element_text(size=14),
#                                 axis.text = element_text(size=7),
#                                 plot.background = element_rect(fill ="transparent", color = NA),
#                                 #panel.grid.major=element_line(size=0.4, linetype = "dashed"), 
#                                 panel.grid.minor= element_blank(),
#                                 panel.grid.major = element_blank(),
#                                 panel.background = element_rect(fill ="transparent"), 
#                                 #axis.line = element_line(colour = "black"),
#                                 axis.ticks = element_line(colour = "black"),
#                                 plot.margin=unit(c(0.5, 0.2, 0, 0),"cm"),
#                                 legend.key.size = unit(1, "cm"),
#                                 legend.spacing.y = unit(-0.3, "cm"),
#                                 legend.position ="left",
#                                 legend.box = "vertical",
#                                 legend.direction = "vertical",
#                                 legend.justification = c(0,1),
#                                 legend.key = element_rect(colour = "transparent", fill = "transparent"),
#                                 legend.background = element_rect(colour = "transparent", fill = "transparent"))
#           
#           
#           
# 
#           
#           
#           
#           
#           ##########
#           # PLOT B :
#           ##########
#           dir.plot.B = paste0(dir, "/it",iteration.indexes[i],"/Segmentation")
#           setwd(dir.plot.B)
#           BIC.df = read.table(paste0(dir.plot.B,"/BIC.df_it", iteration.indexes[i] ,".txt"), sep ="\t", header =TRUE)
#           # grid = seq(-100000000, -10000000, 100)
#           # min_grid <- NULL
#           # for (gr in 1:length(grid)) {
#           #   if ((grid[gr+1] > min(BIC.df$AIC, BIC.df$BIC, BIC.df$HQC, BIC.df$DIC, na.rm =TRUE)) & #BIC.df$AICc)) &
#           #       (grid[gr] < min(BIC.df$AIC, BIC.df$BIC, BIC.df$HQC, BIC.df$DIC, na.rm =TRUE))) { #BIC.df$AICc, ))) {
#           #        min_grid = grid[gr]
#           #   }
#           # }
#           # max_grid <- NULL
#           # for (gr in 1:length(grid)) {
#           #   if((grid[gr+1] >= max(BIC.df$AIC, BIC.df$BIC, BIC.df$HQC,  BIC.df$DIC, na.rm =TRUE)) & #BIC.df$AICc,
#           #      (grid[gr] <= max(BIC.df$AIC, BIC.df$BIC, BIC.df$HQC, BIC.df$DIC, na.rm =TRUE)) ) { #BIC.df$AICc,
#           #       max_grid = grid[gr+1]
#           #   }
#           # }
#           col= c("AIC (Aikake, 1973)" ="blue", 
#                  "BIC (Schwarz, 1978)" ="gray40",
#                  "HQC (Hannan & Quinn, 1979)" ="red",
#                  #"AICc (Hurvich and Tsai, 1989)" ="orange", 
#                  "DIC (Gelman et al., 2004)" ="green")
#           colmin= c("min(AIC)" =17, 
#                  "min(BIC)" =15,
#                  "min(HQC)" =19,
#                  #"AICc (Hurvich and Tsai, 1989)" ="orange", 
#                  "min(DIC)" =23)
#           colmin2= c("min(AIC)" ="blue", 
#                      "min(BIC)" ="gray40",
#                      "min(HQC)" ="red",
#                      #"AICc (Hurvich and Tsai, 1989)" ="orange", 
#                      "min(DIC)" ="green")
#           #******************************************
#           AICcol ="AIC (Aikake, 1973)"; 
#           BICcol ="BIC (Schwarz, 1978)"; 
#           HQCcol= "HQC (Hannan & Quinn, 1979)";
#           # AICccol= "AICc (Hurvich and Tsai, 1989)"; 
#           DICcol= "DIC (Gelman et al., 2004)";
#           #*****************************************
#           AICcol2 ="min(AIC)"; 
#           BICcol2 ="min(BIC)"; 
#           HQCcol2 = "min(HQC)";
#           #AICccol= "AICc (Hurvich and Tsai, 1989)"; 
#           DICcol2= "min(DIC)";
#           
#           
#           plot.B[[i]] = ggplot(BIC.df) + 
#                         geom_line(aes(x = x, y = AIC,  color =AICcol), size = 0.5, na.rm = TRUE) +
#                         geom_line(aes(x = x, y = BIC,  color =BICcol), size = 0.5, na.rm = TRUE) +
#                         geom_line(aes(x = x, y = HQC,  color =HQCcol), size = 0.5, na.rm = TRUE) +
#                         #geom_line(aes(x = x, y = AICc, color =AICccol), size = 0.5) +
#                         geom_line(aes(x = x, y = DIC, color =DICcol), size = 0.5, na.rm = TRUE) +
#                         scale_x_continuous(name = "Number of segments K", expand = c(0,0), limits =c(1, length(BIC.df$BIC)))+
#                         scale_y_continuous(name = "Criterion value") +
#                                            # limits = c(min_grid, max_grid), 
#                                            #breaks = seq(min_grid, max_grid, 50))+
#                                            #breaks = scales::pretty_breaks(n = 2), 
#                                            #labels = scales::scientific ) + 
#                         coord_cartesian(clip = 'off')+
#                         # ylab("BIC, AIC, HQC, DIC") +
#                         geom_point(aes(x = x, y = AIC, color= AICcol), shape =2, size = 4, na.rm = TRUE) +
#                         geom_point(aes(x = x, y = BIC, color = BICcol ), shape =0, size = 4, na.rm = TRUE) +
#                         # geom_point(aes(x = x, y = AICc)), shape =6, size = 3, color = "orange", na.rm = TRUE) +    
#                         geom_point(aes(x = x, y = HQC, color= HQCcol), shape =1, size = 4, na.rm = TRUE) +
#                         geom_point(aes(x = x, y = DIC, color= DICcol), shape =5, size = 4, na.rm = TRUE) +
#                         geom_point(aes(x = AICmin, y = AIC[AICmin], color =AICcol2), fill="blue", shape = 17, size = 4,  na.rm = TRUE) +
#                         geom_point(aes(x = BICmin, y = BIC[BICmin], color =BICcol2), fill="gray40", shape = 15, size =4, na.rm = TRUE) +
#                         # geom_point(aes(x = AICcmin, y = AICc[AICcmin]), shape = 25, size = 5, color = "orange", fill = "orange", na.rm = TRUE)+
#                         geom_point(aes(x = HQCmin, y = HQC[HQCmin], color =HQCcol2), fill="red", shape = 19, size =4, na.rm = TRUE)+
#                         geom_point(aes(x = DICmin, y = DIC[DICmin], color =DICcol2), fill="green", shape = 23, size =4,  na.rm = TRUE) +
#                         scale_colour_manual(name=element_blank(),
#                                             values=c(col, colmin2), 
#                                 breaks = c(AICcol, 
#                                            BICcol, 
#                                            HQCcol,
#                                            DICcol,
#                                            AICcol2, 
#                                            BICcol2, 
#                                            HQCcol2,
#                                            DICcol2),
#                                 #labels = c(gaugings.period,gaugings.other,MAPcurve,activ.stage),
#                                 guide = guide_legend(override.aes = list(
#                                                      linetype = c("solid", "solid","solid", "solid",
#                                                                   "solid", "solid","solid", "solid"),
#                                                      shape = c(0, 2, 1, 5, 
#                                                                17,15,19,23))))  +
#                         # scale_shape_manual(name=element_blank(),
#                         #                   values=colmin2, 
#                         #                   breaks = c(AICcol2, 
#                         #                             #AICccol, 
#                         #                              BICcol2, 
#                         #                              HQCcol2,
#                         #                              DICcol2)
#                         #                  ) +
#                         xlab("Number of segments K") + 
#                         theme_light(base_size = 15)+
#                         theme(text = element_text(size=14),
#                               axis.text = element_text(size=7),
#                               plot.background = element_rect(fill ="transparent", color = NA),
#                               panel.grid.major=element_blank(),
#                               panel.grid.minor=element_blank(),
#                               #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                               panel.background = element_rect(fill ="transparent"), 
#                               #axis.line = element_line(colour = "black"),
#                               axis.ticks = element_line(colour = "black"),
#                               plot.margin=unit(c(0.5, 0.05, 0.15, 0.5),"cm"),
#                               #legend.key.size = unit(1, "cm"),
#                               #legend.key.size = unit(1, "cm"),
#                               legend.title = element_blank(),
#                               legend.key.size = unit(0.7, "cm"),
#                               legend.spacing.y = unit(-0.3, "cm"),
#                               legend.position ="left",
#                               legend.box = "vertical",
#                               legend.direction = "vertical",
#                               legend.justification = c(0,1),
#                               legend.key = element_rect(colour = "transparent", fill = "transparent"),
#                               legend.background = element_rect(colour = "transparent", fill = "transparent"))
#           
#           
#           
#           
#           
#           ##########
#           # PLOT C :
#           ##########
#           dir.plot.C = paste0(dir, "/it",iteration.indexes[i],"/Segmentation")
#           setwd(dir.plot.C)
#           #legend:
#           Utot.times = "95% total uncertainty of change point times"
#           MAPtimes ="Change point times (MAP)"
#           realtimes ="Adjusted shift time"
#           Utot.mean = "95% total uncertainty of segment mean"
#           MAPmean = "Segment mean (MAP)"
#           resid ="Residuals"
#           color.segm = c("95% total uncertainty of change point times"="blue")
#                         # "90% total uncertainty of segment mean" ="red")
#           col.segm.vert  = c("Adjusted shift time"  = "dotted", 
#                              "Change point times (MAP)"= "solid")
#           col.segm = c("Segment mean (MAP)"      = "red", 
#                        "Residuals"               = "black")          
#           #color.segm_1 = c("95% total uncertainty of segment mean" ="red")
#           col.segm_1 = c("Segment mean (MAP)"    = "red",
#                          "Residuals"             = "blue")
#           mu.df = read.table(paste0(dir.plot.C,"/df_mu_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
#           
#           iteration.name = c("ITERATION 0", "ITERATION 1.2", "ITERATION 1.3")
#           if (length(mu.df$mu.MAP) > 1) {
#               ts.plus.df = read.table(paste0(dir.plot.C,"/df.shift.times.plus_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
#               ts.df = read.table(paste0(dir.plot.C,"/df.shift.times_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
#               resid.df = read.table(paste0(dir.plot.C,"/df_residuals_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
#               tau.df = read.table(paste0(dir.plot.C,"/df_tau_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
#             
#               plot.C[[i]] = ggplot()+
#                             scale_y_continuous(name = TeX("Residual $r$"),
#                                                limits =c(min(resid.df$alpha - 2*resid.df$sigma.tot),
#                                                          max(resid.df$alpha + 2*resid.df$sigma.tot))) + 
#                             coord_cartesian(clip = 'off')+
#                             geom_errorbar(data = resid.df,    aes(x= alpha_t, 
#                                                                   ymin= (alpha - 2*sigma.tot), 
#                                                                   ymax = (alpha + 2*sigma.tot), 
#                                                                   color = resid),
#                                           size = 0.2, 
#                                           width=0.02*(tail(resid.df$alpha_t,1) - resid.df$alpha_t[1])) +
#                             geom_point(data = resid.df, aes(x = alpha_t, 
#                                                 y = alpha, 
#                                                 col=resid),   #color = Shift.Q$Q_Gaug), 
#                                        shape = 1, size = 0.9)+
#                             geom_rect(data=ts.df,  mapping= aes(xmin=  Q2.ts , 
#                                                                 xmax= Q97.ts  ,
#                                                                 ymin=-Inf, 
#                                                                 ymax=Inf, 
#                                                                 fill=Utot.times), 
#                                                                 alpha=0.1 ) +
#                             geom_vline(aes(xintercept = ts.df$ts.res, linetype =MAPtimes), col= "blue", lwd =0.8) +
#                             geom_vline(aes(xintercept = ts.df$ts.real, linetype =realtimes), col="red", lwd =0.5, show.legend = F) +
#                             # geom_rect(data=ts.plus.df, mapping = aes(xmin= ts.res.before, 
#                             #                                          xmax=ts.res.plus, 
#                             #                                          ymin=Q2.mu,
#                             #                                          ymax=Q97.mu, 
#                             #                                          fill=Utot.mean), alpha=0.3) + 
#                             geom_segment(data=ts.plus.df, mapping=aes(x =ts.res.before ,
#                                                                       y = mu.res,
#                                                                       xend = ts.res.plus,
#                                                                       yend = mu.res,
#                                                                       col =MAPmean), show.legend = F)+
#                 
#                             scale_linetype_manual(name=element_blank(), 
#                                                   values=col.segm.vert,
#                                                   breaks=c(realtimes, MAPtimes),
#                                                   guide = guide_legend(override.aes = list(
#                                                           col =c("red", "blue")))) +
#                 
#                             scale_fill_manual(name=element_blank(), values=color.segm) +
#                             scale_colour_manual(name = element_blank(), 
#                                                 values=col.segm,
#                                                 breaks=c(
#                                                          MAPmean, 
#                                                          resid),
#                                                 labels = c( 
#                                                            MAPmean, 
#                                                            resid),
#                                                 guide = guide_legend(override.aes = list(
#                                                        linetype = c("solid","blank"),
#                                                        shape = c(NA, 1)))) +
#                             theme_light(base_size = 15) +
#                             theme(text = element_text(size=14),
#                                   axis.text = element_text(size=7),
#                                   #panel.grid.major=element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
#                                   panel.grid.major = element_blank(), 
#                                   panel.grid.minor = element_blank(),
#                                   plot.background = element_rect(fill = "transparent", 
#                                                                  color = NA),
#                                   panel.background = element_rect(fill = "transparent"),
#                                   plot.margin=unit(c(0.5, 0.2, 0.15, 0.5),"cm"),
# 
#                                   # theme(plot.title = element_text(hjust = 0.5),
#                                   #       plot.background = element_rect(fill ="transparent", color = NA),
#                                   #       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                   #       panel.background = element_rect(fill ="transparent"),                                     #       axis.line = element_line(colour = "black"),
#                                   #       axis.ticks = element_line(colour = "black"),
#                                   #       plot.margin=unit(c(0.2,0.5,0.05,0.05),"cm"),
#                                   legend.key.size = unit(1, "cm"),
#                                   legend.spacing.y = unit(-0.3, "cm"),
#                                   legend.position ="left",
#                                   legend.box = "vertical",
#                                   legend.direction = "vertical",
#                                   legend.justification = c(0,1),
#                                   legend.key = element_rect(colour = "transparent", fill = "transparent"),
#                                   legend.background = element_rect(colour = "transparent", fill = "transparent"))
#                                   # legend.text=element_text(size=8),
#                                   # legend.spacing.x = unit(0.05, 'cm'),
#                                   # legend.spacing.y = unit(0.05, 'cm'),
#                                   # legend.key.height=unit(0.5,"line"),
#                                   # legend.key.width=unit(0.8,"line"),
#                             if (is.null(df.limni)==FALSE) {
#                                   # plot.C[[i]]  = plot.C[[i]]  + 
#                                   # scale_x_continuous(name = "time [days]", expand = c(0,0), 
#                                   #                   limits =c(0,tail(df.limni$t_limni,1)))
#                                   plot.C[[i]]  = plot.C[[i]]  + 
#                                   scale_x_continuous(name = "time [days]", expand = c(0,0))
#                             } else {
#                                   # plot.C[[i]]  = plot.C[[i]]  + 
#                                   # scale_x_continuous(name = "time [days]", expand = c(0,0), 
#                                   #                   limits =c(0,tail(Shift.Q$alpha_t,1)))
#                                   plot.C[[i]]  = plot.C[[i]]  + 
#                                   scale_x_continuous(name = "time [days]", expand = c(0,0))
#                             }    
#           } else {
#               resid.df = read.table(paste0(dir.plot.C,"/df_residuals_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
#               plot.C[[i]] = ggplot()+
#                             geom_point(data = resid.df, aes(x = alpha_t, 
#                                                             y = alpha, 
#                                                             col=resid),
#                                        shape = 1, size = 0.9)+
#                             geom_errorbar(data = resid.df, aes(x= alpha_t, 
#                                                               ymin= (alpha - 2*sigma.tot), 
#                                                               ymax= (alpha + 2*sigma.tot), 
#                                                               color = resid),
#                                           size = 0.2, 
#                                           width=0.03*(tail(resid.df$alpha_t,1) - resid.df$alpha_t[1])) +
#                             scale_y_continuous(name = TeX("Residual $r$"),
#                                                limits =c(min(resid.df$alpha - 2*resid.df$sigma.tot),
#                                                          max(resid.df$alpha + 2*resid.df$sigma.tot))) + 
#                             coord_cartesian(clip = 'off') + 
#                             # geom_rect(data=mu.df, mapping = aes(xmin= resid.df$alpha_t[1], 
#                             #                                     xmax=tail(resid.df$alpha_t,1), 
#                             #                                     ymin=mu.q2,
#                             #                                     ymax=mu.q97, 
#                             #                                     fill=Utot.mean), alpha=0.3) + 
#                             geom_segment(data=mu.df, mapping=aes(x =resid.df$alpha_t[1], 
#                                                                       y =mu.MAP, 
#                                                                       xend = tail(resid.df$alpha_t,1), 
#                                                                       yend = mu.MAP, 
#                                                                       col =MAPmean)) +
#                             scale_fill_manual(name=element_blank(), values=color.segm) +
#                             scale_colour_manual(name = element_blank(), 
#                                     values=col.segm,
#                                     breaks=c( MAPmean,  resid),
#                                     labels = c(MAPmean, resid),
#                                     guide = guide_legend(override.aes = list(
#                                                          linetype = c("solid","blank"),
#                                                          shape = c(NA, 1)))) +
#                             theme_light(base_size = 15)+
#                             theme(text = element_text(size=14),
#                                   axis.text = element_text(size=7),
#                                   #plot.title = element_text(hjust = 0.5),
#                                   #panel.grid.major=element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
#                                   panel.grid.major = element_blank(), 
#                                   panel.grid.minor = element_blank(),
#                                   plot.background = element_rect(fill = "transparent", color = NA),
#                                   panel.background = element_rect(fill = "transparent"),
#                                   plot.margin=unit(c(0.5,0.2,0.15, 0.5),"cm"),
#                                   legend.position ="none")
#                             if (is.null(df.limni)==FALSE) {
#                                #   plot.C[[i]]  =   plot.C[[i]]  + 
#                                # scale_x_continuous(name = "time [days]", expand = c(0,0), 
#                                # limits =c(0,tail(df.limni$t_limni,1)))
#                                plot.C[[i]]  =   plot.C[[i]]  +
#                                scale_x_continuous(name = "time [days]", expand = c(0,0))
#                             } else {
#                                # plot.C[[i]]  =   plot.C[[i]]  + 
#                                # scale_x_continuous(name = "time [days]", expand = c(0,0), 
#                                # limits =c(0,tail(Shift.Q$alpha_t,1)))
#                                plot.C[[i]]  =   plot.C[[i]]  + scale_x_continuous(name = "time [days]", expand = c(0,0))
#                             }    
#           }
#           
#           iteration.plots[[i]] =  plot_grid(plot.A[[i]] + theme(legend.position="none"), 
#                                             plot.B[[i]] + theme(legend.position="none"), 
#                                             plot.C[[i]] + theme(legend.position="none"),
#                                             #labels = c('A)', 'B)', 'C)'),    
#                                             label_size = 20, ncol =3, nrow =1)
#           #print(iteration.plots[[i]]) 
#           # now add the title
#           title <- ggdraw() + 
#                    draw_label(iteration.name[i],
#                                x = 0, hjust = 0,
#                               size = 15) +
#                    theme(# add margin on the left of the drawing canvas,
#                          # so title is aligned with left edge of first plot
#                          plot.margin = margin(0, 0, 0, 7))
#           iteration.plots2[[i]] =  plot_grid(title,
#                                              iteration.plots[[i]],  
#                                              ncol = 1,
#                                              # rel_heights values control vertical title margins
#                                              rel_heights = c(0.2, 1))
#      }    
#      ###################################################  
#       Legend.title <- ggdraw() + 
#                       draw_label("Legend",
#                                  fontface = 'bold',
#                                   x = 0, hjust = 0,
#                                   size = 15) +
#                       theme(# add margin on the left of the drawing canvas,
#                       # so title is aligned with left edge of first plot
#                       plot.margin = margin(0, 0, 0, 7))
#      glegend.A <- cowplot::get_legend(  plot.A[[1]] )
#      glegend.B <- cowplot::get_legend(  plot.B[[1]] )
#      glegend.C <- cowplot::get_legend(  plot.C[[1]] )
#      glegends =plot_grid(glegend.A, glegend.B, glegend.C, ncol =3, nrow =1)
#      Legends =  plot_grid(Legend.title, 
#                           glegends,
#                           ncol = 1,
#                           # rel_heights values control vertical title margins
#                           rel_heights = c(0.1, 1))
#           
#           # glegend.A <- get_legend(  plot.A[[1]] )
#           # glegend.B <- get_legend(  plot.B[[1]] )
#           # glegend.C <- get_legend(  plot.C[[1]] )
#           # blank_p <- plot_spacer() + theme_void()
#           # # combine legend 1 & 2
#           # leg12 <- plot_grid(glegend.A, glegend.B, glegend.C, blank_p,
#           #                    ncol = 4)
#      
#      
#      #points of diagram:
#      rectangles = data.frame(xmin=c(0, 6, 6, 6,       12, 12,   12, 12,     18, 18),
#                              ymin=c(-1, -1, 6, -8,    1, -3,    -6, -10,     -8.5, -11.5), 
#                              xmax=c(4, 10, 10, 10,      16, 16,    16, 16,       22,22), 
#                              ymax= c(1, 1, 8, -6 ,    3, -1,    -4,-8,       -6.5, -9.5))
#      arrows = data.frame(xmin=c(4, 4, 4,              10,10,      10,10,         16,16),
#                          ymin=c(0, 0, 0,              0,0,      -7,-7,       -9,-9), 
#                          xmax=c(6, 6, 6,              12,12,    12,12,       18,18), 
#                          ymax= c(0, 7, -7,            2, -2,    -5,-9,       -7.5,-10.5))
#      Textscheme =c("Iteration 0", 
#                    "Iteration 1.2", "Iteration 1.1", "Iteration 1.3", 
#                    "Iteration 2.2.1", "Iteration 2.2.2",
#                    "Iteration 2.3.1", "Iteration 2.3.2",
#                    "Iteration 3.3.2.1", "Iteration 3.3.2.2")
#      scheme = ggplot() +
#               ggtitle("A)") +
#               geom_rect(data=rectangles,  
#                         mapping= aes(xmin=  xmin , 
#                                             xmax= xmax  ,
#                                             ymin=ymin, 
#                                             ymax=ymax),
#                         color="gray40", fill="white",
#                         alpha=0.1 ) +
#               annotate("text", x=rectangles$xmin +2,  y=rectangles$ymin+1, 
#                                   label= Textscheme, color = "black", size=6) +
#       
#               geom_segment(data=arrows, aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
#                                      arrow = arrow(length = unit(0.03, "npc"))) +
#               theme_void()+
#               theme(plot.title = element_text(hjust = 0, size=20,face="bold"))
#      
#      title.a <- ggdraw() + 
#                 draw_label("B)", fontface = 'bold', x = 0, hjust = 0, size = 20) +
#                 theme(# add margin on the left of the drawing canvas,
#                       # so title is aligned with left edge of first plot
#                 plot.margin = margin(0, 0, 0, 7))
#      title.b <- ggdraw() + 
#                 draw_label("C)", fontface = 'bold', x = 0, hjust = 0, size = 20) +
#                 theme(# add margin on the left of the drawing canvas,
#                 # so title is aligned with left edge of first plot
#                 plot.margin = margin(0, 0, 0, 7))
#      title.c <- ggdraw() + 
#                 draw_label("D)", fontface = 'bold', x = 0, hjust = 0, size = 20) +
#                 theme(# add margin on the left of the drawing canvas,
#                 # so title is aligned with left edge of first plot
#                 plot.margin = margin(0, 0, 0, 7))
#      init.title =  plot_grid(title.a, title.b, title.c,   ncol = 3)
# 
#          # Final plot with all plots combined:
#      all.plots = plot_grid( scheme,
#                             init.title,
#                             iteration.plots2[[1]],
#                             iteration.plots2[[2]], 
#                             iteration.plots2[[3]],
#                             # iteration.plots2[[4]],
#                             Legends,
#                             label_size = 12, ncol =1, nrow=6,
#                             rel_heights = c(1, 0.1, 1, 1, 1, 1))
#       
#      ggsave(plot = all.plots, filename =paste0(dir,"/Figure_for_article.png"), 
#             bg = "transparent", width = 12, height =19, dpi = 400)
# }
#
# 



















































# Figure for paper Darienzo et al., 2020, with the main results of some iterations for the case study of Meyras
###############################################################################################################
plot.figure.article = function(dir,
                               iteration.indexes,
                               df.limni,
                               grid_RC.ylim,  grid_RC.xlim, 
                               grid_RC.xstep, grid_RC.ystep,
                               grid_RC.ylim.log,
                               ticks_RC.y.log, 
                               RC.x.labels, RC.y.labels, 
                               u.m.Qgaug, u.m.Hgaug,
                               ncontrol,
                               plot.Q.log) {
###############################################################################################################  

# For each chosen iteration plot  3 subplots: 
##*******************************************
  # A) RC of reference, 
  plot.A =list()
  plot.A.legend =list()
  # B) Criteria evolution with increasing nS
  plot.B =list()
  plot.B.legend =list()
  # C) Results of the segmentation of residuals
  plot.C =list()
  plot.C.legend =list()
  #
  iteration.plots =list()
  iteration.plots2 = list()
  limits.residual.axis =list()
  text.size.legend = 14
  
  
  
  
  
  
  
  for (i in 1:length(iteration.indexes)) {
  # Let's study the iteration iteration.indexes[i] :
    
    ##########
    # PLOT A :
    ##########
    dir.plot.A = paste0(dir, "/it",iteration.indexes[i],"/BaRatin")
    dir.plot.B = paste0(dir, "/it",iteration.indexes[i],"/Segmentation")
    setwd(dir.plot.A)
    gaug.tot = read.table(paste0(dir.plot.A,"/df_gaug_tot",iteration.indexes[i],".txt"), 
                          sep ="\t", header =TRUE)
    # gaug.it = read.table(paste0(dir.plot.A,"/df_gaug_it",iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
    gaug.it = read.table(paste0(dir.plot.B,"/df.gaug.periods.txt"), sep ="\t", header =TRUE)
    gaug.it$color = as.character(gaug.it$color)    
    
    
    #uncomment below just for the paper's figure:
    # if (i ==1) {
    #    gaug.it$color[which(gaug.it$color == "red")] = "red"   
    #    gaug.it$color[which(gaug.it$color == "chartreuse")] = "darkviolet"
    #    gaug.it$color[which(gaug.it$color == "blue")] = "chartreuse"  
    # 
    # } else if (i==2){
    #        gaug.it$color[which(gaug.it$color == "red")] = "chartreuse" 
    #        gaug.it$color[which(gaug.it$color == "blue")] = "orange"
    # } else if (i==3){
    #   gaug.it$color[which(gaug.it$color == "red")]   = "darkviolet"
    #   gaug.it$color[which(gaug.it$color == "blue")] = "cyan"
    # }
    
    if (i ==1) {
      gaug.it$color[which(gaug.it$color == "red")] = "red"   
      gaug.it$color[which(gaug.it$color == "green")] = "chartreuse"
      gaug.it$color[which(gaug.it$color == "blue")] = "darkviolet"  
      
    } else if (i==2){
      gaug.it$color[which(gaug.it$color == "red")] = "chartreuse" 
      gaug.it$color[which(gaug.it$color == "green")] = "orange"
    } else if (i==3){
      gaug.it$color[which(gaug.it$color == "red")]   = "darkviolet"
      gaug.it$color[which(gaug.it$color == "green")] = "cyan"
    }
    
    file.X        = paste0(dir.plot.A,"/Hgrid.txt")
    file.env.tot  = paste0(dir.plot.A,"/Qrc_TotalU.env")
    file.spag.tot = paste0(dir.plot.A,"/Qrc_TotalU.spag")
    file.env.par  = paste0(dir.plot.A,"/Qrc_ParamU.env")
    file.spag.par = paste0(dir.plot.A,"/Qrc_ParamU.spag")
    file.spag.max = paste0(dir.plot.A,"/Qrc_Maxpost.spag")
    file.summary  = paste0(dir.plot.A,"/Results_Summary.txt")
    model         = ReadModel(ModelFile=paste0(dir.plot.A,"/Config_Model.txt"))
    
    # legend:
    MAPcurve        <- "MAP rating curve"
    Ustruct.MAP     <- "95% rating curve uncertainty"
    Utot            <- "95% total uncertainty"
    gaugings.period <- "Gaugings of the current period \n with 95% uncertainty interval"
    gaugings.other  <- "Gaugings of the other periods"
    activ.stage     <- "Control activation stages (MAP)" 
    col.RC          <- c( "MAP rating curve" = "black")
    col.uncert      <- c( "95% rating curve uncertainty" = "red") #"salmon") #"indianred1") 
    col.gaug        <- c("Gaugings of the current period \n with 95% uncertainty interval" = 21) 
    #cols <- c("Current periods gaugings"= "black", "Other periods gaugings" = "gray", 
    #          "MAP rating curve"= "blue") 
    #al <- c("U_tot"=0.2, "U_par"=0.4, "U_MAP"=0.3)    #"U_par"="indianred1",
    maxpost     <- read.table(file = file.spag.max)
    maxpost.par <- read.table(file = file.summary)
    Hg          <- read.table(file = file.X)
    ktransit    = NULL
    for (n in 1:ncontrol) {
      ktransit[n] = maxpost.par[16, ncontrol*3+2+n]    
    }
    gamma1 <- maxpost.par[16,ncontrol*3+1];
    str.err = 0
    if (remnant.err.model == "Linear") {  
      gamma2 <- maxpost.par[16,ncontrol*3+2]
      str.err = gamma1+gamma2*maxpost
    } else {
      str.err = gamma1
    }
    env.tot = EnvelopLayer(file=file.env.tot, Xfile=file.X, color = Utot, alpha = 1)
    #spag=SpaghettiLayer(file=file.spag, Xfile=file.X)
    env.par = EnvelopLayer(file=file.env.par, Xfile=file.X , color = "U_par", alpha =1)
    #env.max = EnvelopLayer(file=file.spag.max, Xfile=file.X , color = "black", alpha =1 )
    maxpost.data <- data.frame(x = Hg[which(Hg >= ktransit[1]),1], 
                               y = maxpost[which(Hg >= ktransit[1]),1], 
                               y2 = str.err[which(Hg >= ktransit[1]),1])
    new_line     <- element_line(color = "grey", size = 0.1, linetype = 2)
    loadfonts()
    ######################

    plot.A[[i]] = ggplot()
      if (plot.Q.log == TRUE) {
          plot.A[[i]] =  plot.A[[i]]  +
          coord_trans(limx = grid_RC.xlim, 
                      limy = grid_RC.ylim.log) +
          scale_y_log10(na.value=-10,
                        breaks=ticks_RC.y.log,
                        labels=ticks_RC.y.log) +
          scale_x_continuous(breaks=seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep))
      } else  {
          coord_trans(limx = grid_RC.xlim, limy = grid_RC.ylim) +
          scale_y_continuous(breaks=seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep))+
          scale_x_continuous(breaks=seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep))
      }
      #########
      plot.A[[i]] = plot.A[[i]] +
      ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
      xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]"))) + 
      #annotation_logticks(base = 10, sides = "l", scaled = TRUE, colour = "black", size = 0.3, linetype = 1)+
      # scale_y_log10(na.value=-10, breaks=ticks_RC.y.log, labels=ticks_RC.y.log) +
      #               scale_x_continuous(breaks=seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep)) +
      geom_ribbon(data = na.omit(maxpost.data), aes(x = x, ymin =(y-y2) , 
                                                    ymax = (y+y2), 
                                                    fill = Ustruct.MAP ), alpha = 0.2)+
      # geom_errorbar(data= gaug.tot, aes(x=h_Gaug, ymin =Q_Gaug-2*uQ_Gaug,
      #                                   ymax =Q_Gaug+2*uQ_Gaug, color =gaugings.other),
      #               width=.08, size = 0.1) +
      # geom_point(data= gaug.tot,    aes(x = h_Gaug, y = Q_Gaug, 
      #                                   color =gaugings.other), size = 1) +
      geom_errorbar(data = gaug.it,  aes(x=hP, ymin = QP - 2*uQP,
                                        ymax =QP + 2*uQP),
                                        color = "black" , # gaug.it$color,
                    width=.08, size = 0.05)+
      geom_line(data = na.omit(maxpost.data), aes(x = x, y= y, color = MAPcurve), size = 0.8, na.rm = TRUE)+
      geom_point( aes(x = -999,  y = -999,   pch = gaugings.period), fill ="gray",     size = 2) +        
      geom_point(data= gaug.it,  aes(x = hP, y = QP), pch  =21 , fill = gaug.it$color, size = 2) +
      #geom_vline(aes(xintercept = ktransit, color =activ.stage), size = 0.3, linetype ="dashed")+   
      scale_colour_manual(name = element_blank(), 
                          values=c(col.RC),
                          breaks=c(MAPcurve),
                          guide = guide_legend(override.aes = list(
                                  linetype = c("solid"),
                                  shape = c(NA)))) +
      scale_fill_manual(name   = element_blank(), 
                        values = col.uncert) +
      scale_shape_manual(name  = element_blank(), 
                         values= col.gaug,
                         breaks= c(gaugings.period),
                         guide = guide_legend(override.aes = list(
                                              linetype = c("solid"),
                                              shape = c(21),
                                              fill  = "gray"))) +
      theme_light(base_size = 20)+
      theme(  text = element_text(size=14),
              axis.text = element_text(size=7),
              plot.background = element_rect(fill ="transparent", color = NA),
              #panel.grid.major=element_line(size=0.4, linetype = "dashed"), 
              panel.grid.minor= element_blank(),
              panel.grid.major = element_blank(),
              panel.background = element_rect(fill ="transparent"), 
              #axis.line = element_line(colour = "black"),
              axis.ticks = element_line(colour = "black"),
              plot.margin=unit(c(0.5, 0.2, 0, 0),"cm"),
              legend.key.size = unit(1, "cm"),
              legend.spacing.y = unit(-0.3, "cm"),
              legend.position ="left",
              legend.box = "vertical",
              legend.text = element_text(margin=margin(1,0,0,0.5), size=text.size.legend),
              legend.direction = "vertical",
              legend.justification = c(0,1),
              legend.key = element_rect(colour = "transparent", fill = "transparent"),
              legend.background = element_rect(colour = "transparent", fill = "transparent"))
    
    
    
    
    # if (i == length(iteration.indexes)){
    #   plot.A[[i]] =  plot.A[[i]] +
    #   theme(plot.margin=unit(c(0.5, 0.2, 1, 0),"cm"))
    # }
    
    
    
    
    
    
    
    
    ##########
    # PLOT B :
    ##########
    dir.plot.B = paste0(dir, "/it",iteration.indexes[i],"/Segmentation")
    setwd(dir.plot.B)
    BIC.df = read.table(paste0(dir.plot.B,"/BIC.df_it", iteration.indexes[i] ,".txt"), sep ="\t", header =TRUE)
    col= c("AIC" = "black", 
           "BIC" = "black",
           "HQC" = "black",
           "DIC" = "black")
    
    colmin= c("min(AIC)" =17,  #triangle
              "min(BIC)" =15,  # square
              "min(HQC)" =19,  # circle
              "min(DIC)" =23)  # romb
    
    colmin2= c("min(AIC)" ="black", 
               "min(BIC)" ="black",
               "min(HQC)" ="black",
               "min(DIC)" ="black")
    #******************************************
    AICcol ="AIC"; 
    BICcol ="BIC"; 
    HQCcol= "HQC";
    DICcol= "DIC";
    #*****************************************
    AICcol2 = "min(AIC)"; 
    BICcol2 = "min(BIC)"; 
    HQCcol2 = "min(HQC)";
    DICcol2 = "min(DIC)";
    
    
    plot.B[[i]] = ggplot(BIC.df) + 
      geom_line(aes(x = x, y = AIC),  color ="gray80", size = 0.2, na.rm = TRUE) +
      geom_line(aes(x = x, y = BIC),  color ="gray80", size = 0.2, na.rm = TRUE) +
      geom_line(aes(x = x, y = HQC),  color ="gray80", size = 0.2, na.rm = TRUE) +
      geom_line(aes(x = x, y = DIC),  color ="gray80", size = 0.2, na.rm = TRUE) +
      #geom_line(aes(x = x, y = AICc, color =AICccol), size = 0.5) +
      scale_x_continuous(name = "Number of segments K", expand = c(0,0), limits =c(1, length(BIC.df$BIC)))+
      scale_y_continuous(name = "Criterion value") +
      # limits = c(min_grid, max_grid), 
      # breaks = seq(min_grid, max_grid, 50))+
      # breaks = scales::pretty_breaks(n = 2), 
      # labels = scales::scientific ) + 
      coord_cartesian(clip = 'off')+
      geom_point(aes(x = x, y = AIC, color = AICcol),   shape =2, size = 3, na.rm = TRUE) +   # triangle
      geom_point(aes(x = x, y = BIC, color = BICcol ),  shape =0, size = 3, na.rm = TRUE) +   # square
      geom_point(aes(x = x, y = HQC, color = HQCcol),   shape =1, size = 3, na.rm = TRUE) +   # circle
      geom_point(aes(x = x, y = DIC, color = DICcol),   shape =5, size = 3, na.rm = TRUE) +   # romb
      #
      geom_point(aes(x = AICmin, y = AIC[AICmin], fill =AICcol2), shape = 17, size =3,  na.rm = TRUE) +
      geom_point(aes(x = BICmin, y = BIC[BICmin], fill =BICcol2), shape = 15, size =3,  na.rm = TRUE) +
      geom_point(aes(x = HQCmin, y = HQC[HQCmin], fill =HQCcol2), shape = 19, size =3,  na.rm = TRUE)+
      geom_point(aes(x = DICmin, y = DIC[DICmin], fill =DICcol2), shape = 23, size =3,  na.rm = TRUE) +
      #
      scale_colour_manual(name   = element_blank(),
                          values = c(col), 
                          breaks = c(AICcol, 
                                     BICcol, 
                                     HQCcol,
                                     DICcol),
                          #labels = c(gaugings.period,gaugings.other,MAPcurve,activ.stage),
                          guide = guide_legend(override.aes = list(
                            linetype = c("solid", "solid","solid", "solid"),
                            shape = c(2, 0, 1, 5))))  +
      scale_fill_manual(name = element_blank(),
                        values = colmin2,
                        breaks = c(AICcol2,
                                   BICcol2,
                                   HQCcol2,
                                   DICcol2),
                        guide = guide_legend(override.aes = list(
                                             shape = c(17,15,19,23)))) +
      #guides(fill=guide_legend(nrow=2,byrow=TRUE))+
      xlab("Number of segments K") + 
      theme_light(base_size = 15)+
      theme(text = element_text(size=14),
            axis.text = element_text(size=7),
            plot.background = element_rect(fill ="transparent", color = NA),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill ="transparent"), 
            axis.ticks = element_line(colour = "black"),
            plot.margin=unit(c(0.5, 0.05, 0.15, 0.5),"cm"),
            legend.title         = element_blank(),
            legend.text          = element_text(margin=margin(1, 0, 0 , 0.5), size=text.size.legend),
            legend.key.size      = unit(1, "cm"),
            legend.spacing.y     = unit(-0.3, "cm"),
            legend.position      = "left",
            legend.box           = "horizontal",
            legend.direction     = "vertical",
            legend.justification = c(0,1),
            legend.key           = element_rect(colour = "transparent", fill = "transparent"),
            legend.background    = element_rect(colour = "transparent", fill = "transparent"))
    # if (i == length(iteration.indexes)){
    #   plot.B[[i]] =  plot.B[[i]] +
    #     theme(plot.margin=unit(c(0.5, 0.2, 1, 0),"cm"))
    # }
    
    
    
    
    
    
    
    
    
    ##########
    # PLOT C :
    ##########
    dir.plot.C = paste0(dir, "/it",iteration.indexes[i],"/Segmentation")
    setwd(dir.plot.C)
    if (iteration.indexes[i] == 1) {
       limits.residual.axis[[i]] = c(-10,10) 
    } else if (iteration.indexes[i] == 5){
       limits.residual.axis[[i]] = c(-1,1) 
    } else if (iteration.indexes[i] == 8){
       limits.residual.axis[[i]] = c(-3,3) 
    }
    #legend:
    Utot.times = "95% uncertainty of change point times"
    MAPtimes   = "Change point times (MAP)"
    realtimes  = "Adjusted shift time"
    Utot.mean  = "95% uncertainty of segment means"
    MAPmean    = "Segment mean (MAP)"
    resid      = "Residuals with 95% uncertainty interval"
    #
    color.segm     = c("95% uncertainty of change point times"    = "blue",
                       "95% uncertainty of segment means"         = "gray70")
    col.segm.vert  = c("Adjusted shift time"                      = "dotted", 
                       "Change point times (MAP)"                 = "solid")
    col.segm       = c("Segment mean (MAP)"                       = "gray70", 
                       "Residuals with 95% uncertainty interval"  = "gray40") # "black") 
    col.segm_1     = c("Segment mean (MAP)"                       = "red",
                       "Residuals with 95% uncertainty interval"  = "blue")
    # "90% total uncertainty of segment mean" ="red")
    #color.segm_1 = c("95% total uncertainty of segment mean"  ="red")
    mu.df = read.table(paste0(dir.plot.C,"/df_mu_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
    iteration.name = c("ITERATION 0", "ITERATION 1.2", "ITERATION 1.3")

    if (length(mu.df$mu.MAP) > 1) {
      ts.plus.df = read.table(paste0(dir.plot.C,"/df.shift.times.plus_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
      ts.df = read.table(paste0(dir.plot.C,"/df.shift.times_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
      resid.df = read.table(paste0(dir.plot.C,"/df_residuals_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
      resid.df$color = gaug.it$color
      color.segments = unique(resid.df$color) 
      tau.df = read.table(paste0(dir.plot.C,"/df_tau_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
      
      error.max = 0
      error.min = 0
      for (iii in 1:length(resid.df$alpha)) {
           error.min[iii] = max((resid.df$alpha[iii] - 2*resid.df$sigma.tot[iii]), limits.residual.axis[[i]][1])
           error.max[iii] = min((resid.df$alpha[iii] + 2*resid.df$sigma.tot[iii]), limits.residual.axis[[i]][2])
      }
      resid.df$error.max = error.max
      resid.df$error.min = error.min
      #######
      plot.C[[i]] = ggplot()+
        scale_y_continuous(name = TeX("Residual $r$"),
                           # limits =c(min(resid.df$alpha - 2*resid.df$sigma.tot),
                           #           max(resid.df$alpha + 2*resid.df$sigma.tot))) + 
                           limits = limits.residual.axis[[i]],
                           expand = c(0,0))+
        coord_cartesian(clip = 'off')+
        geom_errorbar(data = resid.df,    aes(x     = alpha_t, 
                                              ymin  = error.min, 
                                              ymax  = error.max), 
                                              color = resid.df$color,
                      size = 0.001, 
                      width=0.02*(tail(resid.df$alpha_t,1) - resid.df$alpha_t[1])) +
        geom_point(aes(x = resid.df$alpha_t[1], y = -9999,  col=resid),   #color = Shift.Q$Q_Gaug), 
                   shape = 21, fill ="white", size = 1.2) +
        geom_point(data = resid.df, aes(x = alpha_t, 
                                        y = alpha), 
                                        col=resid.df$color,   #color = Shift.Q$Q_Gaug), 
                   shape = 21, fill ="white", size = 1.2) +
        geom_rect(data=ts.df,  mapping= aes(xmin=  Q2.ts , 
                                            xmax= Q97.ts  ,
                                            ymin=-Inf, 
                                            ymax=Inf, 
                                            fill=Utot.times),  alpha=0.1 ) +
        geom_vline(aes(xintercept = ts.df$ts.res, linetype  = MAPtimes),  col= "blue", lwd = 1) +
        geom_vline(aes(xintercept = ts.df$ts.real, linetype = realtimes), col="red",  lwd = 1) + #,  show.legend = F) +
        geom_rect(data=ts.plus.df, mapping = aes(xmin  = ts.res.before,
                                                 xmax  = ts.res.plus,
                                                 ymin  = Q2.mu,
                                                 ymax  = Q97.mu),
                                                 fill  = color.segments,
                                                 alpha = 0.2) +
        geom_rect(data=ts.plus.df, mapping = aes(xmin  = ts.res.before,
                                                 xmax  = ts.res.plus,
                                                 ymin  = -9999,
                                                 ymax  = -9998,
                                                 fill  = Utot.mean),
                                                 alpha = 0.2) +
        geom_segment(data=ts.plus.df, mapping=aes(x    = ts.res.before ,
                                                  y    = mu.res,
                                                  xend = ts.res.plus,
                                                  yend = mu.res),
                                                  col  = color.segments, 
                                                  #show.legend = T, 
                                                  size= 0.5)+
        geom_segment(data=ts.plus.df, mapping=aes(x    = ts.res.before ,
                                                  y    = -9999,
                                                  xend = ts.res.plus,
                                                  yend = -9999,
                                                  col  = MAPmean), 
                                                  size= 0.5)+
        
        scale_linetype_manual(name   = element_blank(), 
                              values = col.segm.vert,
                              breaks = c(realtimes, 
                                         MAPtimes),
                              guide  = guide_legend(override.aes = list(
                                                     col =c("red", "blue")))) +
        scale_fill_manual(name   = element_blank(), 
                          values = color.segm,
                          guide  = guide_legend(override.aes = list(
                                                linetype = c("solid", "solid"))))  +
        scale_colour_manual(name   = element_blank(), 
                            values = col.segm,
                            breaks = c(MAPmean, 
                                       resid),
                            labels = c(MAPmean, 
                                       resid),
                            guide  = guide_legend(override.aes = list(
                                                  linetype = c("solid","blank"),
                                                  shape = c(NA, 1)))) +
        theme_light(base_size = 15) +
        theme(text = element_text(size=14),
              axis.text = element_text(size=7),
              #panel.grid.major=element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.background  = element_rect(fill = "transparent", 
                                             color = NA),
              panel.background = element_rect(fill = "transparent"),
              plot.margin = unit(c(0.5, 0.2, 0.15, 0.5),"cm"),
              # axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
              
              # theme(plot.title = element_text(hjust = 0.5),
              #       plot.background = element_rect(fill ="transparent", color = NA),
              #       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              #       panel.background = element_rect(fill ="transparent"),                                     #       axis.line = element_line(colour = "black"),
              #       axis.ticks = element_line(colour = "black"),
              #       plot.margin=unit(c(0.2,0.5,0.05,0.05),"cm"),
              legend.key.size      = unit(1, "cm"),
              legend.text          = element_text(margin=margin(1.3, 0, 0, 0.5), size=text.size.legend),
              legend.spacing.y     = unit(-0.3, "cm"),
              legend.position      = "left",
              legend.box           = "vertical",
              legend.direction     = "vertical",
              legend.justification = c(0,1),
              legend.key           = element_rect(colour = "transparent", fill = "transparent"),
              legend.background    = element_rect(colour = "transparent", fill = "transparent"))
      # legend.text=element_text(size=8),
      # legend.spacing.x = unit(0.05, 'cm'),
      # legend.spacing.y = unit(0.05, 'cm'),
      # legend.key.height=unit(0.5,"line"),
      # legend.key.width=unit(0.8,"line"),
      if (is.null(df.limni)==FALSE) {
        # plot.C[[i]]  = plot.C[[i]]  + 
        # scale_x_continuous(name = "time [days]", expand = c(0,0), 
        #                   limits =c(0,tail(df.limni$t_limni,1)))
        plot.C[[i]]  = plot.C[[i]]  + 
          scale_x_continuous(name = "time [days]", expand = c(0,0))
      } else {
        # plot.C[[i]]  = plot.C[[i]]  + 
        # scale_x_continuous(name = "time [days]", expand = c(0,0), 
        #                   limits =c(0,tail(Shift.Q$alpha_t,1)))
        plot.C[[i]]  = plot.C[[i]]  + 
          scale_x_continuous(name = "time [days]", expand = c(0,0))
      }
    ##########
    } else {
    ##########
      resid.df = read.table(paste0(dir.plot.C,"/df_residuals_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
      plot.C[[i]] = ggplot()+
        geom_point(data = resid.df, aes(x = alpha_t, 
                                        y = alpha, 
                                        col=resid),
                   shape = 1, size = 0.9)+
        geom_errorbar(data = resid.df, aes(x     = alpha_t, 
                                           ymin  = (alpha - 2*sigma.tot), 
                                           ymax  = (alpha + 2*sigma.tot), 
                                           color = resid),
                      size = 0.2, 
                      width=0.03*(tail(resid.df$alpha_t,1) - resid.df$alpha_t[1])) +
        scale_y_continuous(name = TeX("Residual $r$"),
                           limits =c(min(resid.df$alpha - 2*resid.df$sigma.tot),
                                     max(resid.df$alpha + 2*resid.df$sigma.tot)),
                           expand = c(0,0)) + 
        coord_cartesian(clip = 'off') + 
        # geom_rect(data=mu.df, mapping = aes(xmin= resid.df$alpha_t[1], 
        #                                     xmax=tail(resid.df$alpha_t,1), 
        #                                     ymin=mu.q2,
        #                                     ymax=mu.q97, 
        #                                     fill=Utot.mean), alpha=0.3) + 
        geom_segment(data=mu.df, mapping=aes(x    =  resid.df$alpha_t[1], 
                                             y    =  mu.MAP, 
                                             xend =  tail(resid.df$alpha_t,1), 
                                             yend =  mu.MAP, 
                                             col  =  MAPmean),
                     size = 1.5) +
        scale_fill_manual(name=element_blank(), values=color.segm) +
        scale_colour_manual(name   = element_blank(), 
                            values = col.segm,
                            breaks = c( MAPmean,  resid),
                            labels = c(MAPmean, resid),
                            guide  = guide_legend(override.aes = list(
                                     linetype = c("solid","blank"),
                                     shape = c(NA, 1)))) +
        theme_light(base_size = 15)+
        theme(text = element_text(size=14),
              axis.text = element_text(size=7),
              #plot.title = element_text(hjust = 0.5),
              #panel.grid.major=element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "transparent", color = NA),
              panel.background = element_rect(fill = "transparent"),
              plot.margin=unit(c(0.5,0.2,0.15, 0.5),"cm"),
              legend.position ="none")
      if (is.null(df.limni)==FALSE) {
        #   plot.C[[i]]  =   plot.C[[i]]  + 
        # scale_x_continuous(name = "time [days]", expand = c(0,0), 
        # limits =c(0,tail(df.limni$t_limni,1)))
        plot.C[[i]]  =   plot.C[[i]]  +   scale_x_continuous(name = "time [days]", expand = c(0,0))
      } else {
        # plot.C[[i]]  =   plot.C[[i]]  + 
        # scale_x_continuous(name = "time [days]", expand = c(0,0), 
        # limits =c(0,tail(Shift.Q$alpha_t,1)))
        plot.C[[i]]  =   plot.C[[i]]  + scale_x_continuous(name = "time [days]", expand = c(0,0))
      }    
    }
    # if (i == length(iteration.indexes)){
    #   plot.C[[i]] =  plot.C[[i]] +
    #     theme(plot.margin=unit(c(0.5, 0.2, 2, 0.5),"cm"))
    # }
    

    
    
    #-------------------------------------------------------------------------------
    iteration.plots[[i]] =  plot_grid(plot.A[[i]] + theme(legend.position="none"), 
                                      plot.B[[i]] + theme(legend.position="none"), 
                                      plot.C[[i]] + theme(legend.position="none"),
                                      #labels = c('A)', 'B)', 'C)'),    
                                      label_size = 20, ncol =3, nrow =1)
    # now add the title
    title <- ggdraw() + 
             draw_label(iteration.name[i],
                 x = 0, hjust = 0,
                 size = 15) +
             theme(plot.margin = margin(0, 0, 0, 7))
    iteration.plots2[[i]] =  plot_grid(title,
                                       iteration.plots[[i]],  
                                       ncol = 1,
                                       # rel_heights values control vertical title margins
                                       rel_heights = c(0.2, 1))
  }    
  
  
  
  #########################################################################
  Legend.title <- ggdraw() + 
                  draw_label("Legend",
                              fontface = 'bold',
                              x = 0, hjust = 0,
                              size = 15) +
                  theme(# add margin on the left of the drawing canvas,
                        # so title is aligned with left edge of first plot
                        plot.margin = margin(10, 0, 5, 0))
  glegend.A <- cowplot::get_legend(  plot.A[[1]] )
  glegend.B <- cowplot::get_legend(  plot.B[[1]] )
  glegend.C <- cowplot::get_legend(  plot.C[[1]] )
  glegends =plot_grid(glegend.A, glegend.B, glegend.C, 
                      ncol =3, nrow =1)
  Legends =  plot_grid(Legend.title, 
                       glegends,
                       ncol = 1,
                       # rel_heights values control vertical title margins
                       rel_heights = c(0.3, 1))
  
  
  
  
  #########################################################################
  #points of diagram:
  rectangles = data.frame(xmin=c(0,    6, 6, 6,       12, 12,   12, 12,   12, 12,       18,   18 , 18 ),
                          
                          ymin=c(-1,   -1, 6, -8,     8, 4,    1, -3,    -6, -10,      -7, -10, -13), 
                        
                          xmax=c(4,    10, 10, 10,    16, 16,  16, 16,    16, 16,      22,  22 , 22 ), 
                          
                          ymax= c(1,    1, 8, -6 ,    10, 6 ,   3, -1,    -4,-8,        -5, -8, -11 ))
  
  
  arrows = data.frame(xmin=c(4, 4, 4,    10, 10,        10,10,      10,10,       16, 16, 16),
                      
                      ymin=c(0, 0, 0,     7, 7,         0,0,      -7,-7,         -9, -9, -9), 
                      
                      xmax=c(6, 6, 6,    12, 12,        12,12,    12,12,         18, 18, 18), 
                      
                      ymax= c(0, 7, -7,   9, 5,         2, -2,    -5,-9,        -6, -9, -12))
  Textscheme =c("Iteration 0", 
                "Iteration 1.2",     "Iteration 1.1",      "Iteration 1.3", 
                "Iteration 2.1.1",   "Iteration 2.1.2",
                "Iteration 2.2.1",   "Iteration 2.2.2",
                "Iteration 2.3.1",   "Iteration 2.3.2",
                "Iteration 3.3.2.1", "Iteration 3.3.2.2" ,"Iteration 3.3.2.3" )
  scheme = ggplot() +
    ggtitle("a)") +
    geom_rect(data=rectangles,  
              mapping= aes(xmin=  xmin , 
                           xmax= xmax  ,
                           ymin=ymin, 
                           ymax=ymax),
              color="gray50", fill="white",
              alpha=0.1, size =0.2) +
    annotate("text", x=rectangles$xmin +2,  y=rectangles$ymin+1, 
             label= Textscheme, color = "black", size=5) +
    geom_segment(data=arrows, aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
                 arrow = arrow(length = unit(0.03, "npc")), size=0.3, color="gray10", linetype ="dotted") +
    theme_void()+
    theme(plot.title = element_text(hjust = 0, size=20,face="bold"))
  
  title.a <- ggdraw() + 
    draw_label("b)", fontface = 'bold', x = 0, hjust = 0, size = 20) +
    theme(# add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))
  title.b <- ggdraw() + 
    draw_label("c)", fontface = 'bold', x = 0, hjust = 0, size = 20) +
    theme(# add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))
  title.c <- ggdraw() + 
    draw_label("d)", fontface = 'bold', x = 0, hjust = 0, size = 20) +
    theme(# add margin on the left of the drawing canvas,
          # so title is aligned with left edge of first plot
          plot.margin = margin(0, 0, 0, 7))
  init.title =  plot_grid(title.a, title.b, title.c,   ncol = 3)
  
  # Final plot with all plots combined:
  all.plots = plot_grid( scheme,
                         init.title,
                         iteration.plots2[[1]],
                         iteration.plots2[[2]], 
                         iteration.plots2[[3]],
                         # iteration.plots2[[4]],
                         Legends,
                         label_size = 12,
                         ncol =1, nrow=6,
                         rel_heights = c(1, 0.1, 1, 1, 1, 0.8))
  
  ggsave(plot = all.plots, filename =paste0(dir,"/Figure_for_article.png"), 
         bg = "transparent", width = 12, height =19, dpi = 400)
  pdf(paste0(dir,"/Figure4.pdf"), 12,18 ,useDingbats=F)
  print(all.plots)
  dev.off()
}

























# plot figure for each iteration of segmentation procedure
###############################################################################################################
plot.all.iterations = function(dir,
                               df.limni,
                               grid_RC.ylim,  grid_RC.xlim, 
                               grid_RC.xstep, grid_RC.ystep,
                               grid_RC.ylim.log,
                               ticks_RC.y.log, 
                               RC.x.labels, RC.y.labels, 
                               u.m.Qgaug, u.m.Hgaug,
                               ncontrol,
                               plot.Q.log) {
  ###############################################################################################################  
  
  # For each chosen iteration plot  3 subplots: 
  ##*******************************************
  # A) RC of reference, 
  plot.A =list()
  plot.A.legend =list()
  # B) Criteria evolution with increasing nS
  plot.B =list()
  plot.B.legend =list()
  # C) Results of the segmentation of residuals
  plot.C =list()
  plot.C.legend =list()
  iteration.plots =list()
  iteration.plots2 = list()
  
  text.size.legend = 14
  
  limits.residual.axis =list()
  
  
  
  
  
  for (i in 1:length(iteration.indexes)) {
    # Let's study the iteration iteration.indexes[i] :
    
    ##########
    # PLOT A :
    ##########
    dir.plot.A = paste0(dir, "/it",iteration.indexes[i],"/BaRatin")
    dir.plot.B = paste0(dir, "/it",iteration.indexes[i],"/Segmentation")
    setwd(dir.plot.A)
    gaug.tot = read.table(paste0(dir.plot.A,"/df_gaug_tot",iteration.indexes[i],".txt"), 
                          sep ="\t", header =TRUE)
    # gaug.it = read.table(paste0(dir.plot.A,"/df_gaug_it",iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
    gaug.it = read.table(paste0(dir.plot.B,"/df.gaug.periods.txt"), sep ="\t", header =TRUE)
    gaug.it$color = as.character(gaug.it$color)    
    
    
    #uncomment below just for the paper's figure:
    if (i ==1) {
      gaug.it$color[which(gaug.it$color == "red")] = "red"
      gaug.it$color[which(gaug.it$color == "green")] = "purple"
      gaug.it$color[which(gaug.it$color == "blue")] = "green"
    } else if (i==2){
      gaug.it$color[which(gaug.it$color == "red")] = "green"
      gaug.it$color[which(gaug.it$color == "green")] = "green"
    } else if (i==3){
      gaug.it$color[which(gaug.it$color == "red")]   = "purple"
      gaug.it$color[which(gaug.it$color == "green")] = "purple"
    }
    
    file.X=paste0(dir.plot.A,"/Hgrid.txt")
    file.env.tot=paste0(dir.plot.A,"/Qrc_TotalU.env")
    file.spag.tot=paste0(dir.plot.A,"/Qrc_TotalU.spag")
    file.env.par=paste0(dir.plot.A,"/Qrc_ParamU.env")
    file.spag.par=paste0(dir.plot.A,"/Qrc_ParamU.spag")
    file.spag.max=paste0(dir.plot.A,"/Qrc_Maxpost.spag")
    file.summary = paste0(dir.plot.A,"/Results_Summary.txt")
    model=ReadModel(ModelFile=paste0(dir.plot.A,"/Config_Model.txt"))
    
    # legend:
    MAPcurve <-  "MAP rating curve"
    Ustruct.MAP <- "95% remnant uncertainty of the MAP RC"
    Utot <- "95% total uncertainty"
    gaugings.period <- "Gaugings of the current period"
    gaugings.other <- "Gaugings of the other periods"
    activ.stage <- "Control activation stages (MAP)"
    col.uncert <- c( 
      "95% remnant uncertainty of the MAP RC" = "red") #"salmon") #"indianred1") 
    cols <- c("Gaugings of the current period"="black", 
              #"Gaugings of the other periods"="gray80",
              "MAP rating curve"="black")
    #cols <- c("Current periods gaugings"= "black", "Other periods gaugings" = "gray", 
    #          "MAP rating curve"= "blue") 
    #al <- c("U_tot"=0.2, "U_par"=0.4, "U_MAP"=0.3)    #"U_par"="indianred1",
    maxpost <- read.table(file = file.spag.max)
    maxpost.par <-read.table(file = file.summary)
    Hg <- read.table(file = file.X)
    ktransit=NULL
    for (n in 1:ncontrol) {
      ktransit[n] = maxpost.par[16, ncontrol*3+2+n]    
    }
    gamma1 <- maxpost.par[16,ncontrol*3+1];
    str.err = 0
    if (remnant.err.model == "Linear") {  
      gamma2 <- maxpost.par[16,ncontrol*3+2]
      str.err = gamma1+gamma2*maxpost
    } else {
      str.err = gamma1
    }
    env.tot = EnvelopLayer(file=file.env.tot, Xfile=file.X, color = Utot, alpha = 1)
    #spag=SpaghettiLayer(file=file.spag, Xfile=file.X)
    env.par = EnvelopLayer(file=file.env.par, Xfile=file.X , color = "U_par", alpha =1)
    #env.max = EnvelopLayer(file=file.spag.max, Xfile=file.X , color = "black", alpha =1 )
    maxpost.data <- data.frame(x = Hg[which(Hg >= ktransit[1]),1], 
                               y = maxpost[which(Hg >= ktransit[1]),1], 
                               y2 = str.err[which(Hg >= ktransit[1]),1])
    new_line <- element_line(color = "grey", size = 0.1, linetype = 2)
    loadfonts()
    ##########
    
    plot.A[[i]] = ggplot()
    if (plot.Q.log == TRUE) {
      plot.A[[i]] =  plot.A[[i]]  +
        coord_trans(limx = grid_RC.xlim, 
                    limy = grid_RC.ylim.log) +
        scale_y_log10(na.value=-10,
                      breaks=ticks_RC.y.log,
                      labels=ticks_RC.y.log) +
        scale_x_continuous(breaks=seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep))
    } else  {
      coord_trans(limx = grid_RC.xlim, limy = grid_RC.ylim) +
        scale_y_continuous(breaks=seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep))+
        scale_x_continuous(breaks=seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep))
    }
    
    plot.A[[i]] = plot.A[[i]] +
      ylab(bquote(.(RC.y.labels) ~ .("[") ~ m^3*s^-1 ~ .("]"))) +  
      xlab(bquote(.(RC.x.labels) ~ .("[") ~ m ~ .("]"))) + 
      #annotation_logticks(base = 10, sides = "l", scaled = TRUE, colour = "black", size = 0.3, linetype = 1)+
      # scale_y_log10(na.value=-10, breaks=ticks_RC.y.log, labels=ticks_RC.y.log) +
      #               scale_x_continuous(breaks=seq(grid_RC.xlim[1], grid_RC.xlim[2], grid_RC.xstep)) +
      geom_ribbon(data = na.omit(maxpost.data), aes(x = x, ymin =(y-y2) , 
                                                    ymax = (y+y2), 
                                                    fill = Ustruct.MAP ), alpha = 0.2)+
      # geom_errorbar(data= gaug.tot, aes(x=h_Gaug, ymin =Q_Gaug-2*uQ_Gaug, 
      #                                   ymax =Q_Gaug+2*uQ_Gaug, color =gaugings.other), 
      #               width=.1, size = 0.2)+
      # geom_point(data= gaug.tot,    aes(x = h_Gaug, y = Q_Gaug, 
      #                                   color =gaugings.other), size = 1) +
      geom_errorbar(data= gaug.it,  aes(x=hP, ymin =QP-2*uQP,
                                        ymax =QP+2*uQP),
                    color = "black" , # gaug.it$color,
                    width=.03, size = 0.1)+
      geom_line(data = na.omit(maxpost.data), 
                aes(x = x, y= y, color = MAPcurve), 
                size = 0.8, na.rm = TRUE)+
      geom_point( aes(x = -10, y = 10,   color = gaugings.period), size = 1) +        
      geom_point(data= gaug.it,  aes(x = hP, y = QP), pch  =21 , fill = gaug.it$color, size = 2) +
      
      #geom_vline(aes(xintercept = ktransit, color =activ.stage), size = 0.3, linetype ="dashed")+   
      scale_colour_manual(name = element_blank(), 
                          values=c(cols),
                          breaks=c(gaugings.period,
                                   #gaugings.other,
                                   MAPcurve),
                          #labels = c(gaugings.period,gaugings.other,MAPcurve,activ.stage),
                          guide = guide_legend(override.aes = list(
                            linetype = c("blank", "solid"),
                            shape = c(19, NA)))) +
      scale_fill_manual(name=element_blank(), 
                        values= col.uncert) +
      theme_light(base_size = 20)+
      theme(  text = element_text(size=14),
              axis.text = element_text(size=7),
              plot.background = element_rect(fill ="transparent", color = NA),
              #panel.grid.major=element_line(size=0.4, linetype = "dashed"), 
              panel.grid.minor= element_blank(),
              panel.grid.major = element_blank(),
              panel.background = element_rect(fill ="transparent"), 
              #axis.line = element_line(colour = "black"),
              axis.ticks = element_line(colour = "black"),
              plot.margin=unit(c(0.5, 0.2, 0, 0),"cm"),
              legend.key.size = unit(1, "cm"),
              legend.spacing.y = unit(-0.3, "cm"),
              legend.position ="left",
              legend.box = "vertical",
              legend.text = element_text(margin=margin(1,0,0,0.5), size=text.size.legend),
              legend.direction = "vertical",
              legend.justification = c(0,1),
              legend.key = element_rect(colour = "transparent", fill = "transparent"),
              legend.background = element_rect(colour = "transparent", fill = "transparent"))
    
    
    
    
    # if (i == length(iteration.indexes)){
    #   plot.A[[i]] =  plot.A[[i]] +
    #   theme(plot.margin=unit(c(0.5, 0.2, 1, 0),"cm"))
    # }
    
    
    
    
    
    
    
    
    ##########
    # PLOT B :
    ##########
    dir.plot.B = paste0(dir, "/it",iteration.indexes[i],"/Segmentation")
    setwd(dir.plot.B)
    BIC.df = read.table(paste0(dir.plot.B,"/BIC.df_it", iteration.indexes[i] ,".txt"), sep ="\t", header =TRUE)
    # grid = seq(-100000000, -10000000, 100)
    # min_grid <- NULL
    # for (gr in 1:length(grid)) {
    #   if ((grid[gr+1] > min(BIC.df$AIC, BIC.df$BIC, BIC.df$HQC, BIC.df$DIC, na.rm =TRUE)) & #BIC.df$AICc)) &
    #       (grid[gr] < min(BIC.df$AIC, BIC.df$BIC, BIC.df$HQC, BIC.df$DIC, na.rm =TRUE))) { #BIC.df$AICc, ))) {
    #        min_grid = grid[gr]
    #   }
    # }
    # max_grid <- NULL
    # for (gr in 1:length(grid)) {
    #   if((grid[gr+1] >= max(BIC.df$AIC, BIC.df$BIC, BIC.df$HQC,  BIC.df$DIC, na.rm =TRUE)) & #BIC.df$AICc,
    #      (grid[gr] <= max(BIC.df$AIC, BIC.df$BIC, BIC.df$HQC, BIC.df$DIC, na.rm =TRUE)) ) { #BIC.df$AICc,
    #       max_grid = grid[gr+1]
    #   }
    # }
    # col= c("AIC (Aikake, 1974)" ="blue", 
    #        "BIC (Schwarz, 1978)" ="gray40",
    #        "HQC (Hannan & Quinn, 1979)" ="red",
    #        #"AICc (Hurvich and Tsai, 1989)" ="orange", 
    #        "DIC (Gelman et al., 2004)" ="green")
    # colmin= c("min(AIC)" =17, 
    #           "min(BIC)" =15,
    #           "min(HQC)" =19,
    #           #"AICc (Hurvich and Tsai, 1989)" ="orange", 
    #           "min(DIC)" =23)
    # colmin2= c("min(AIC)" ="blue", 
    #            "min(BIC)" ="gray40",
    #            "min(HQC)" ="red",
    #            #"AICc (Hurvich and Tsai, 1989)" ="orange", 
    #            "min(DIC)" ="green")
    col= c("AIC (Akaike, 1974)" ="black", 
           "BIC (Schwarz, 1978)" ="black",
           "HQC (Hannan & Quinn, 1979)" ="black",
           #"AICc (Hurvich and Tsai, 1989)" ="orange", 
           "DIC (Gelman et al., 2004)" ="black")
    colmin= c("min(AIC)" =17, 
              "min(BIC)" =15,
              "min(HQC)" =19,
              #"AICc (Hurvich and Tsai, 1989)" ="orange", 
              "min(DIC)" =23)
    colmin2= c("min(AIC)" ="black", 
               "min(BIC)" ="black",
               "min(HQC)" ="black",
               #"AICc (Hurvich and Tsai, 1989)" ="orange", 
               "min(DIC)" ="black")
    #******************************************
    AICcol ="AIC (Akaike, 1974)"; 
    BICcol ="BIC (Schwarz, 1978)"; 
    HQCcol= "HQC (Hannan & Quinn, 1979)";
    # AICccol= "AICc (Hurvich and Tsai, 1989)"; 
    DICcol= "DIC (Gelman et al., 2004)";
    #*****************************************
    AICcol2 ="min(AIC)"; 
    BICcol2 ="min(BIC)"; 
    HQCcol2 = "min(HQC)";
    #AICccol= "AICc (Hurvich and Tsai, 1989)"; 
    DICcol2= "min(DIC)";
    
    
    plot.B[[i]] = ggplot(BIC.df) + 
      geom_line(aes(x = x, y = AIC),  color ="gray", size = 0.1, na.rm = TRUE) +
      geom_line(aes(x = x, y = BIC),  color ="gray", size = 0.1, na.rm = TRUE) +
      geom_line(aes(x = x, y = HQC),  color ="gray", size = 0.1, na.rm = TRUE) +
      geom_line(aes(x = x, y = DIC),  color ="gray", size = 0.1, na.rm = TRUE) +
      #geom_line(aes(x = x, y = AICc, color =AICccol), size = 0.5) +
      scale_x_continuous(name = "Number of segments K", expand = c(0,0), limits =c(1, length(BIC.df$BIC)))+
      scale_y_continuous(name = "Criterion value") +
      # limits = c(min_grid, max_grid), 
      #breaks = seq(min_grid, max_grid, 50))+
      #breaks = scales::pretty_breaks(n = 2), 
      #labels = scales::scientific ) + 
      coord_cartesian(clip = 'off')+
      # ylab("BIC, AIC, HQC, DIC") +
      geom_point(aes(x = x, y = AIC, color= AICcol), shape =2, size = 3, na.rm = TRUE) +
      geom_point(aes(x = x, y = BIC, color = BICcol ), shape =0, size = 3, na.rm = TRUE) +
      # geom_point(aes(x = x, y = AICc)), shape =6, size = 3, color = "orange", na.rm = TRUE) +    
      geom_point(aes(x = x, y = HQC, color= HQCcol), shape =1, size = 3, na.rm = TRUE) +
      geom_point(aes(x = x, y = DIC, color= DICcol), shape =5, size = 3, na.rm = TRUE) +
      geom_point(aes(x = AICmin, y = AIC[AICmin], color =AICcol2), fill="black", shape = 17, size = 3,  na.rm = TRUE) +
      geom_point(aes(x = BICmin, y = BIC[BICmin], color =BICcol2), fill="black", shape = 15, size =3, na.rm = TRUE) +
      # geom_point(aes(x = AICcmin, y = AICc[AICcmin]), shape = 25, size = 5, color = "orange", fill = "orange", na.rm = TRUE)+
      geom_point(aes(x = HQCmin, y = HQC[HQCmin], color =HQCcol2), fill="black", shape = 19, size =3, na.rm = TRUE)+
      geom_point(aes(x = DICmin, y = DIC[DICmin], color =DICcol2), fill="black", shape = 23, size =3,  na.rm = TRUE) +
      scale_colour_manual(name=element_blank(),
                          values=c(col, colmin2), 
                          breaks = c(AICcol, 
                                     BICcol, 
                                     HQCcol,
                                     DICcol,
                                     AICcol2, 
                                     BICcol2, 
                                     HQCcol2,
                                     DICcol2),
                          #labels = c(gaugings.period,gaugings.other,MAPcurve,activ.stage),
                          guide = guide_legend(override.aes = list(
                            linetype = c("solid", "solid","solid", "solid",
                                         "solid", "solid","solid", "solid"),
                            shape = c(0, 2, 1, 5, 
                                      17,15,19,23))))  +
      # scale_shape_manual(name=element_blank(),
      #                   values=colmin2, 
      #                   breaks = c(AICcol2, 
      #                             #AICccol, 
      #                              BICcol2, 
      #                              HQCcol2,
      #                              DICcol2)
      #                  ) +
      xlab("Number of segments K") + 
      theme_light(base_size = 15)+
      theme(text = element_text(size=14),
            axis.text = element_text(size=7),
            plot.background = element_rect(fill ="transparent", color = NA),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_rect(fill ="transparent"), 
            #axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            plot.margin=unit(c(0.5, 0.05, 0.15, 0.5),"cm"),
            #legend.key.size = unit(1, "cm"),
            #legend.key.size = unit(1, "cm"),
            legend.title = element_blank(),
            legend.text = element_text(margin=margin(1,0,0,0.5), size=text.size.legend),
            legend.key.size = unit(0.7, "cm"),
            legend.spacing.y = unit(-0.3, "cm"),
            legend.position ="left",
            legend.box = "vertical",
            legend.direction = "vertical",
            legend.justification = c(0,1),
            legend.key = element_rect(colour = "transparent", fill = "transparent"),
            legend.background = element_rect(colour = "transparent", fill = "transparent"))
    # if (i == length(iteration.indexes)){
    #   plot.B[[i]] =  plot.B[[i]] +
    #     theme(plot.margin=unit(c(0.5, 0.2, 1, 0),"cm"))
    # }
    
    
    
    
    
    ##########
    # PLOT C :
    ##########
    dir.plot.C = paste0(dir, "/it",iteration.indexes[i],"/Segmentation")
    setwd(dir.plot.C)
    if (iteration.indexes[i] == 1) {
      limits.residual.axis[[i]] = c(-10,10) 
    } else if (iteration.indexes[i] == 5){
      limits.residual.axis[[i]] = c(-1,1) 
    } else if (iteration.indexes[i] == 8){
      limits.residual.axis[[i]] = c(-3,3) 
    }
    #legend:
    Utot.times = "95% uncertainty of change point times"
    MAPtimes ="Change point times (MAP)"
    realtimes ="Adjusted shift time"
    Utot.mean = "95% total uncertainty of segment mean"
    MAPmean = "Segment mean (MAP)"
    resid ="Residuals with error bars at 95%"
    color.segm = c("95% uncertainty of change point times"="blue")
    # "90% total uncertainty of segment mean" ="red")
    col.segm.vert  = c("Adjusted shift time"     = "dotted", 
                       "Change point times (MAP)"= "solid")
    col.segm = c("Segment mean (MAP)"                 = "red", 
                 "Residuals with error bars at 95%"   = "black") # "gray50")          
    #color.segm_1 = c("95% total uncertainty of segment mean" ="red")
    col.segm_1 = c("Segment mean (MAP)"                  = "red",
                   "Residuals with error bars at 95%"    = "blue")
    mu.df = read.table(paste0(dir.plot.C,"/df_mu_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
    
    iteration.name = c("ITERATION 0", "ITERATION 1.2", "ITERATION 1.3")
    
    if (length(mu.df$mu.MAP) > 1) {
      ts.plus.df = read.table(paste0(dir.plot.C,"/df.shift.times.plus_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
      ts.df = read.table(paste0(dir.plot.C,"/df.shift.times_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
      resid.df = read.table(paste0(dir.plot.C,"/df_residuals_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
      resid.df$color = gaug.it$color
      tau.df = read.table(paste0(dir.plot.C,"/df_tau_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
      
      error.max = 0
      error.min = 0
      for (iii in 1:length(resid.df$alpha)) {
        error.min[iii] = max((resid.df$alpha[iii] - 2*resid.df$sigma.tot[iii]), limits.residual.axis[[i]][1])
        error.max[iii] = min((resid.df$alpha[iii] + 2*resid.df$sigma.tot[iii]), limits.residual.axis[[i]][2])
      }
      resid.df$error.max = error.max
      resid.df$error.min = error.min
      
      plot.C[[i]] = ggplot()+
        scale_y_continuous(name = TeX("Residual $r$"),
                           # limits =c(min(resid.df$alpha - 2*resid.df$sigma.tot),
                           #           max(resid.df$alpha + 2*resid.df$sigma.tot))) + 
                           limits = limits.residual.axis[[i]],
                           expand = c(0,0))+
        coord_cartesian(clip = 'off')+
        geom_errorbar(data = resid.df,    aes(x     = alpha_t, 
                                              ymin  = error.min, 
                                              ymax  = error.max), 
                      color = resid.df$color,
                      size = 0.001, 
                      width=0.02*(tail(resid.df$alpha_t,1) - resid.df$alpha_t[1])) +
        geom_point(aes(x = resid.df$alpha_t[1], y = -9999,  col=resid),   #color = Shift.Q$Q_Gaug), 
                   shape = 21, fill ="white", size = 1.2) +
        geom_point(data = resid.df, aes(x = alpha_t, 
                                        y = alpha), 
                   col=resid.df$color,   #color = Shift.Q$Q_Gaug), 
                   shape = 21, fill ="white", size = 1.2) +
        geom_rect(data=ts.df,  mapping= aes(xmin=  Q2.ts , 
                                            xmax= Q97.ts  ,
                                            ymin=-Inf, 
                                            ymax=Inf, 
                                            fill=Utot.times), 
                  alpha=0.1 ) +
        geom_vline(aes(xintercept = ts.df$ts.res, linetype =MAPtimes), col= "blue", lwd = 1) +
        geom_vline(aes(xintercept = ts.df$ts.real, linetype =realtimes), col="red", lwd = 0.8,  show.legend = F) +
        geom_rect(data=ts.plus.df, mapping = aes(xmin= ts.res.before,
                                                 xmax=ts.res.plus,
                                                 ymin=Q2.mu,
                                                 ymax=Q97.mu),
                  fill="red", 
                  alpha=0.2, 
                  show.legend=F) +
        geom_segment(data=ts.plus.df, mapping=aes(x    = ts.res.before ,
                                                  y    = mu.res,
                                                  xend = ts.res.plus,
                                                  yend = mu.res,
                                                  col  = MAPmean), 
                     show.legend = F, size= 0.7)+
        
        scale_linetype_manual(name=element_blank(), 
                              values=col.segm.vert,
                              breaks=c(realtimes, MAPtimes),
                              guide = guide_legend(override.aes = list(
                                col =c("red", "blue")))) +
        
        scale_fill_manual(name=element_blank(), 
                          values=color.segm,
                          guide = guide_legend(override.aes = list(
                            linetype = c("solid"))))  +
        scale_colour_manual(name = element_blank(), 
                            values=col.segm,
                            breaks=c( MAPmean,  resid),
                            labels = c(MAPmean, resid),
                            guide = guide_legend(override.aes = list(
                              linetype = c("solid","blank"),
                              shape = c(NA, 1)))) +
        theme_light(base_size = 15) +
        theme(text                 = element_text(size=14),
              axis.text            = element_text(size=7),
              #panel.grid.major    = element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
              panel.grid.major     = element_blank(), 
              panel.grid.minor     = element_blank(),
              plot.background      = element_rect(fill = "transparent",  color = NA),
              panel.background     = element_rect(fill = "transparent"),
              plot.margin          = unit(c(0.5, 0.2, 0.15, 0.5),"cm"),
              #       plot.title = element_text(hjust = 0.5),
              #       plot.background = element_rect(fill ="transparent", color = NA),
              #       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              #       panel.background = element_rect(fill ="transparent"),                                   
              #       axis.line = element_line(colour = "black"),
              #       axis.ticks = element_line(colour = "black"),
              #       plot.margin=unit(c(0.2,0.5,0.05,0.05),"cm"),
              legend.key.size      = unit(1, "cm"),
              legend.text          = element_text(margin=margin(1,0,0,0.5), size=text.size.legend),
              legend.spacing.y     = unit(-0.3, "cm"),
              legend.position      = "left",
              legend.box           = "vertical",
              legend.direction     = "vertical",
              legend.justification = c(0,1),
              legend.key           = element_rect(colour = "transparent", fill = "transparent"),
              legend.background    = element_rect(colour = "transparent", fill = "transparent"))

      if (is.null(df.limni)==FALSE) {
        # plot.C[[i]]  = plot.C[[i]]  + 
        # scale_x_continuous(name = "time [days]", expand = c(0,0), 
        #                   limits =c(0,tail(df.limni$t_limni,1)))
        plot.C[[i]]  = plot.C[[i]]  + 
          scale_x_continuous(name = "time [days]", expand = c(0,0))
      } else {
        # plot.C[[i]]  = plot.C[[i]]  + 
        # scale_x_continuous(name = "time [days]", expand = c(0,0), 
        #                   limits =c(0,tail(Shift.Q$alpha_t,1)))
        plot.C[[i]]  = plot.C[[i]]  + 
          scale_x_continuous(name = "time [days]", expand = c(0,0))
      }
      #########
    } else {
      resid.df = read.table(paste0(dir.plot.C,"/df_residuals_it", iteration.indexes[i],".txt"), sep ="\t", header =TRUE)
      plot.C[[i]] = ggplot()+
        geom_point(data = resid.df, aes(x = alpha_t, 
                                        y = alpha, 
                                        col=resid),
                   shape = 1, size = 0.9)+
        geom_errorbar(data = resid.df, aes(x= alpha_t, 
                                           ymin= (alpha - 2*sigma.tot), 
                                           ymax= (alpha + 2*sigma.tot), 
                                           color = resid),
                      size = 0.2, 
                      width=0.03*(tail(resid.df$alpha_t,1) - resid.df$alpha_t[1])) +
        scale_y_continuous(name = TeX("Residual $r$"),
                           limits =c(min(resid.df$alpha - 2*resid.df$sigma.tot),
                                     max(resid.df$alpha + 2*resid.df$sigma.tot)),
                           expand = c(0,0)) + 
        coord_cartesian(clip = 'off') + 
        # geom_rect(data=mu.df, mapping = aes(xmin= resid.df$alpha_t[1], 
        #                                     xmax=tail(resid.df$alpha_t,1), 
        #                                     ymin=mu.q2,
        #                                     ymax=mu.q97, 
        #                                     fill=Utot.mean), alpha=0.3) + 
        geom_segment(data=mu.df, mapping=aes(x =resid.df$alpha_t[1], 
                                             y =mu.MAP, 
                                             xend = tail(resid.df$alpha_t,1), 
                                             yend = mu.MAP, 
                                             col =  MAPmean)) +
        scale_fill_manual(name=element_blank(), values=color.segm) +
        scale_colour_manual(name = element_blank(), 
                            values=col.segm,
                            breaks=c( MAPmean,  resid),
                            labels = c(MAPmean, resid),
                            guide = guide_legend(override.aes = list(
                              linetype = c("solid","blank"),
                              shape = c(NA, 1)))) +
        theme_light(base_size = 15)+
        theme(text = element_text(size=14),
              axis.text = element_text(size=7),
              #plot.title = element_text(hjust = 0.5),
              #panel.grid.major=element_line(size=0.4, linetype = "dashed"), panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "transparent", color = NA),
              panel.background = element_rect(fill = "transparent"),
              plot.margin=unit(c(0.5,0.2,0.15, 0.5),"cm"),
              legend.position ="none")
      if (is.null(df.limni)==FALSE) {
        #   plot.C[[i]]  =   plot.C[[i]]  + 
        # scale_x_continuous(name = "time [days]", expand = c(0,0), 
        # limits =c(0,tail(df.limni$t_limni,1)))
        plot.C[[i]]  =   plot.C[[i]]  +   scale_x_continuous(name = "time [days]", expand = c(0,0))
      } else {
        # plot.C[[i]]  =   plot.C[[i]]  + 
        # scale_x_continuous(name = "time [days]", expand = c(0,0), 
        # limits =c(0,tail(Shift.Q$alpha_t,1)))
        plot.C[[i]]  =   plot.C[[i]]  + scale_x_continuous(name = "time [days]", expand = c(0,0))
      }    
    }
    # if (i == length(iteration.indexes)){
    #   plot.C[[i]] =  plot.C[[i]] +
    #     theme(plot.margin=unit(c(0.5, 0.2, 2, 0.5),"cm"))
    # }
    
    
    #-------------------------------------------------------------------------------
    iteration.plots[[i]] =  plot_grid(plot.A[[i]] + theme(legend.position="none"), 
                                      plot.B[[i]] + theme(legend.position="none"), 
                                      plot.C[[i]] + theme(legend.position="none"),
                                      #labels = c('A)', 'B)', 'C)'),    
                                      label_size = 20, ncol =3, nrow =1)
    # now add the title
    title <- ggdraw() + 
      draw_label(iteration.name[i],
                 x = 0, hjust = 0,
                 size = 15) +
      theme(plot.margin = margin(0, 0, 0, 7))
    iteration.plots2[[i]] =  plot_grid(title,
                                       iteration.plots[[i]],  
                                       ncol = 1,
                                       # rel_heights values control vertical title margins
                                       rel_heights = c(0.2, 1))
  }    
  
  
  
  #########################################################################
  Legend.title <- ggdraw() + 
    draw_label("Legend",
               fontface = 'bold',
               x = 0, hjust = 0,
               size = 15) +
    theme(# add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(10, 0, 10, 0))
  glegend.A <- cowplot::get_legend(  plot.A[[1]] )
  glegend.B <- cowplot::get_legend(  plot.B[[1]] )
  glegend.C <- cowplot::get_legend(  plot.C[[1]] )
  glegends =plot_grid(glegend.A, glegend.B, glegend.C, 
                      ncol =3, nrow =1)
  Legends =  plot_grid(Legend.title, 
                       glegends,
                       ncol = 1,
                       # rel_heights values control vertical title margins
                       rel_heights = c(0.4, 1))
  
  
  
  
  #########################################################################
  #points of diagram:
  rectangles = data.frame(xmin=c(0,    6, 6, 6,       12, 12,   12, 12,   12, 12,       18,   18 , 18 ),
                          
                          ymin=c(-1,   -1, 6, -8,     8, 4,    1, -3,    -6, -10,      -7, -10, -13), 
                          
                          xmax=c(4,    10, 10, 10,    16, 16,  16, 16,    16, 16,      22,  22 , 22 ), 
                          
                          ymax= c(1,    1, 8, -6 ,    10, 6 ,   3, -1,    -4,-8,        -5, -8, -11 ))
  
  
  arrows = data.frame(xmin=c(4, 4, 4,    10, 10,        10,10,      10,10,       16, 16, 16),
                      
                      ymin=c(0, 0, 0,     7, 7,         0,0,      -7,-7,         -9, -9, -9), 
                      
                      xmax=c(6, 6, 6,    12, 12,        12,12,    12,12,         18, 18, 18), 
                      
                      ymax= c(0, 7, -7,   9, 5,         2, -2,    -5,-9,        -6, -9, -12))
  Textscheme =c("Iteration 0", 
                "Iteration 1.2",     "Iteration 1.1",      "Iteration 1.3", 
                "Iteration 2.1.1",   "Iteration 2.1.2",
                "Iteration 2.2.1",   "Iteration 2.2.2",
                "Iteration 2.3.1",   "Iteration 2.3.2",
                "Iteration 3.3.2.1", "Iteration 3.3.2.2" ,"Iteration 3.3.2.3" )
  scheme = ggplot() +
    ggtitle("a)") +
    geom_rect(data=rectangles,  
              mapping= aes(xmin=  xmin , 
                           xmax= xmax  ,
                           ymin=ymin, 
                           ymax=ymax),
              color="gray50", fill="white",
              alpha=0.1, size =0.2) +
    annotate("text", x=rectangles$xmin +2,  y=rectangles$ymin+1, 
             label= Textscheme, color = "black", size=5) +
    geom_segment(data=arrows, aes(x = xmin, y = ymin, xend = xmax, yend = ymax),
                 arrow = arrow(length = unit(0.03, "npc")), size=0.3, color="gray10", linetype ="dotted") +
    theme_void()+
    theme(plot.title = element_text(hjust = 0, size=20,face="bold"))
  
  title.a <- ggdraw() + 
    draw_label("b)", fontface = 'bold', x = 0, hjust = 0, size = 20) +
    theme(# add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))
  title.b <- ggdraw() + 
    draw_label("c)", fontface = 'bold', x = 0, hjust = 0, size = 20) +
    theme(# add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))
  title.c <- ggdraw() + 
    draw_label("d)", fontface = 'bold', x = 0, hjust = 0, size = 20) +
    theme(# add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))
  init.title =  plot_grid(title.a, title.b, title.c,   ncol = 3)
  
  # Final plot with all plots combined:
  all.plots = plot_grid( scheme,
                         init.title,
                         iteration.plots2[[1]],
                         iteration.plots2[[2]], 
                         iteration.plots2[[3]],
                         # iteration.plots2[[4]],
                         Legends,
                         label_size = 12,
                         ncol =1, nrow=6,
                         rel_heights = c(1, 0.1, 1, 1, 1, 0.8))
  
  ggsave(plot = all.plots, filename =paste0(dir,"/Figure_for_article.png"), 
         bg = "transparent", width = 12, height =19, dpi = 400)
  pdf(paste0(dir,"/Figure4.pdf"), 12,18 ,useDingbats=F)
  print(all.plots)
  dev.off()
}


