







#####################################################################################################
bt_and_hcrit.plot = function(dir.sed.transp, 
                             ts, ts.before, ts.plus, 
                             bt2.df, bt1.df, hcr.df, 
                             df.limni, gaugings_SPD, 
                             ylimits , d50) {
#####################################################################################################
  #color.m = c("red", "blue", "green"); 
  #typelin = c("solid", "solid", "solid");
  col.b2 <- c("b2 (MAP and CI at 95%)");
  col.b1 <- c("b1 (MAP and CI at 95%)");
  col.hcr <- c("hcr (MAP and CI at 95%)");  
  col.gaug   <- c("Gaugings"); 
  shape.tshift <- c("Shift times (MAP)") 
  colssss <- c("b1 (MAP and CI at 95%)" = "blue",
               "b2 (MAP and CI at 95%)" = "red", 
               "hcr (MAP and CI at 95%)" =  "green",
               "Gaugings" = "black")
  colssss2 <- c("Shift times (MAP)" = "dashed")
  
  bt_and_hcrit.plot = ggplot()+
    geom_vline(aes(xintercept = ts$t.adj, linetype =shape.tshift), linetype = "dashed", col= "blue", lwd =0.6) +
    # geom_segment(mapping= aes(x =ts$t.adj , 
    #                           y = -Inf, 
    #                           xend = ts$t.adj, 
    #                           yend = bt2.df$mean[1:(length(bt2.df$mean)-1)],
    #                           linetype = shape.tshift),  
    #                           color = "blue", lwd =0.7) +
    geom_segment(mapping= aes(x     = ts.before , 
                              y     = bt2.df$mean, 
                              xend  = ts.plus, 
                              yend  = bt2.df$mean, 
                              color = col.b2), 
                 linetype ="solid", size = 1) +
    geom_segment(mapping= aes(x     = ts.before , 
                              y     = bt1.df$mean, 
                              xend  = ts.plus, 
                              yend  = bt1.df$mean, 
                              color = col.b1), 
                 linetype ="solid", size = 1) +
    geom_segment(mapping= aes(x     = ts.before , 
                              y     = hcr.df$mean, 
                              xend  = ts.plus,
                              yend  = hcr.df$mean, 
                              color = col.hcr), 
                 linetype ="solid", size = 1) +
    geom_rect(mapping = aes(xmin = ts.before, 
                            xmax = ts.plus, 
                            ymin = bt2.df$X2.5.,
                            ymax = bt2.df$X97.5.),
              fill="red", alpha=0.3) +
    geom_rect(mapping = aes(xmin = ts.before, 
                            xmax = ts.plus, 
                            ymin = bt1.df$X2.5.,
                            ymax = bt1.df$X97.5.),
              fill="blue", alpha=0.3) +
    geom_rect(mapping = aes(xmin = ts.before, 
                            xmax = ts.plus, 
                            ymin = hcr.df$mean - 2*hcr.df$stdev,
                            ymax = hcr.df$mean + 2*hcr.df$stdev), 
              fill="green", alpha=0.3)+
    geom_line(aes(x = df.limni$t_limni, y= df.limni$h_limni), color="gray70", size=0.2) +
    geom_point(aes(x=gaugings_SPD$t, y = gaugings_SPD$h, color= col.gaug), 
               fill=gaugings_SPD$color, pch=21, size=5)+
    annotate("text", x=0, y=hcr.df$mean[2] + 0.2,  
             label=paste0("Critical stage hcr with d =",d50," m"), color = "green", size=7, vjust =-1, hjust=0) +
    # annotate("text", x=0, y=bt2.df$X97.5.[2] - 0.2, 
    #          label="Riverbed mean elevation b", color = "red", size=5, hjust=0) +
    scale_y_continuous(name = expression("Stage h [m]"),
                       limits =c(ylimits[1], ylimits[2]),
                       breaks = seq(ylimits[1], ylimits[2], ylimits[3])) + 
    scale_x_continuous(name = expression("Time [days]"), expand = c(0,0))+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+
    theme(axis.text         = element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,legend.text      = element_text(size=18)
          ,legend.title     = element_text(size=30)
          ,legend.key.size  = unit(1.5, "cm")
          ,legend.position  = "bottom"
          ,legend.direction = "horizontal"
          #,legend.justification = "vertical"
          ,panel.grid       = element_blank()
          ,plot.margin      = unit(c(0.2, 0.8, 0.2, 0.2),"cm")) +
    scale_colour_manual(name = element_blank(),
                        values = colssss,
                        breaks = c(col.b1, col.b2, col.hcr , col.gaug),
                        guide = guide_legend(override.aes = list(
                                linetype = c("solid", "solid", "solid", "blank"),
                                 shape = c(NA, NA, NA, 21)))) +
    scale_linetype_manual(name=element_blank(), 
                          values=colssss2,
                          guide = guide_legend(override.aes = list( 
                                  col =c("blue"))))
  
  ggsave(bt_and_hcrit.plot, filename = paste0(dir.sed.transp,"/bt_and_hcrit.png"), 
         width = 18, height =8, dpi = 400)
  return(bt_and_hcrit.plot)
}





















#########################################################################################
bt_and_hcrit_zoom.plot = function(dir.sed.transp, 
                                  ts, ts.before, ts.plus, bt2.df, hcr.df, 
                                  df.limni, gaugings_SPD, ylimits, d50 ) {
#########################################################################################
  bt_and_hcrit_zoom.plot = ggplot() +
    # geom_rect(mapping= aes(xmin= ts$q2, xmax=ts$q97 ,
    #                        ymin=-Inf, ymax=Inf), 
    #           fill="blue",alpha=0.2) +
    #geom_vline(aes(xintercept = ts$tMAP), color = "blue", lwd =0.7, linetype = "dashed") +
    geom_vline(aes(xintercept = ts$t.adj), color = "blue", lwd =0.7, linetype = "solid")+ 
    geom_segment(mapping = aes(x =c(ts$t.adj[2]-100,  
                                    ts$t.adj[2]), 
                               y =c(bt2.df$mean[2], 
                                    bt2.df$mean[3]),
                               xend = c(ts$t.adj[2], ts$t.adj[2]+100), 
                               yend = c(bt2.df$mean[2], bt2.df$mean[3])), 
                 color = "red", size = 1) +
    geom_segment(mapping= aes(x =c(ts$t.adj[2]-100,  
                                   ts$t.adj[2]), 
                              y =c(hcr.df$mean[2], 
                                   hcr.df$mean[3]),
                              xend = c(ts$t.adj[2], ts$t.adj[2]+100), 
                              yend = c(hcr.df$mean[2], hcr.df$mean[3])), 
                 color = "green", size = 1) +
    #geom_rect(mapping= aes(xmin= ts_ST$t2[1], xmax=ts_ST$t97[1] ,ymin=-Inf, ymax=Inf), fill="blue",alpha=0.2) +
    geom_line(aes(x = df.limni$t_limni, y= df.limni$h_limni), color="black", size =0.3)+
    geom_point(aes(x=gaugings_SPD$t, y = gaugings_SPD$h), pch=21, fill=gaugings_SPD$color, size=5)+
    annotate("text", x=ts$t.adj[2] - 100, y=hcr.df$mean[2] + 0.1,  
             label=paste0("Critical stage hcr with d50=",d50," m"), color = "green", size=7, hjust=0) +
    annotate("text", x=ts$t.adj[2] - 100, y=bt2.df$X97.5.[2] - 0.1, 
             label="Riverbed mean elevation b", color = "red", size=7, hjust=0) +
    scale_y_continuous(name = expression("Stage h [m]"), limits =c(ylimits[1], ylimits[2] )) + 
    scale_x_continuous(name = expression("Time [days]"), limits =c(ts$t.adj[2] - 100, 
                                                                   ts$t.adj[2] + 100) , 
                       expand = c(0,0))+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+ 
    theme(axis.text=element_text(size=20)
          ,axis.title=element_text(size=30,face="bold")
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="none"
          ,panel.grid = element_blank()
          ,plot.margin=unit(c(0.2, 0.8, 0.2, 0.2),"cm"))
  ggsave(bt_and_hcrit_zoom.plot, filename = paste0(dir.sed.transp,"/bt_and_hcrit_zoom.png"), 
         width = 16, height =8, dpi = 400)
  return(bt_and_hcrit_zoom.plot)
}


























#########################################################################################
bt_hcrit_interpol_zoom.plot = function(dir.sed.transp, b, hcr, event.t, event.h, 
                                       ts, ts.before, ts.plus, bt2.df, hcr.df, 
                                       df.limni, gaugings_SPD, hlimits.ST, d50 ) {
#########################################################################################
  col.b2     <- c("b");
  col.hcr    <- c("hc");  
  col.gaug   <- c("Gaugings"); 
  col.event  <- c("stage of the morphogenic flood")
  shape.tshift <- c("Shift time") 
  colssss <- c("b" = "red", 
               "hc" =  "green",
               "Gaugings" = "black",
               "stage of the morphogenic flood"="black")
  colssss2 <- c("Shift time" = "dashed")
  #####################################
  bt_hcrit_interpol_zoom.plot = ggplot() +
    # geom_rect(mapping= aes(xmin= ts$q2, xmax=ts$q97 ,
    #                        ymin=-Inf, ymax=Inf), 
    #           fill="blue",alpha=0.2) +
    #geom_vline(aes(xintercept = ts$tMAP), color = "blue", lwd =0.7, linetype = "dashed") +
    geom_vline(aes(xintercept = ts$t.adj, linetype=shape.tshift), color = "blue", lwd =0.7) + 
    # geom_segment(mapping = aes(x =c(ts$t.adj[2]-100,  
    #                                 ts$t.adj[2]), 
    #                            y =c(bt2.df$mean[2], 
    #                                 bt2.df$mean[3]),
    #                            xend = c(ts$t.adj[2], ts$t.adj[2]+100), 
    #                            yend = c(bt2.df$mean[2], bt2.df$mean[3])), 
    #              color = "red", size = 1) +
    # geom_segment(mapping= aes(x =c(ts$t.adj[2]-100,  
    #                                ts$t.adj[2]), 
    #                           y =c(hcr.df$mean[2], 
    #                                hcr.df$mean[3]),
    #                            xend = c(ts$t.adj[2], ts$t.adj[2]+100), 
    #                           yend = c(hcr.df$mean[2], hcr.df$mean[3])), 
    #              color = "green", size = 1) +
    geom_line(aes(x=df.limni$t_limni, y = b, color=col.b2), size = 1) +
    geom_line(aes(x=df.limni$t_limni, y = hcr, color=col.hcr), size= 1) +
    #geom_rect(mapping= aes(xmin= ts_ST$t2[1], xmax=ts_ST$t97[1] ,ymin=-Inf, ymax=Inf), fill="blue",alpha=0.2) +
    geom_line(aes(x = df.limni$t_limni, y= df.limni$h_limni), color="gray30", size =0.4)+
    geom_point(aes(x=gaugings_SPD$t, y = gaugings_SPD$h, color = col.gaug), pch=21, fill=gaugings_SPD$color, size=5)+
    annotate("text", x=ts$t.adj[2] - 100, y=hcr.df$mean[2] + 0.2,  
             label=paste0("Triggering stage hc"), color = "green", size=9, hjust=0, vjust =-0.5) +
    annotate("text", x=ts$t.adj[2] - 100, y=bt2.df$X97.5.[2] - 0.2, 
             label="River bed elevation b", color = "red", size=9, hjust=0, vjust =-0.5) +
    geom_point(aes(x=event.t[[2]], y=event.h[[2]], color= col.event), pch= 19, fill="black", size=1.5)+
    scale_y_continuous(name = expression("Stage h [m]"), limits =c(hlimits.ST[1], hlimits.ST[2] )) + 
    scale_x_continuous(name = expression("Time [days]"), limits =c(ts$t.adj[2] - 100, 
                                                                   ts$t.adj[2] + 100) , expand = c(0,0))+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+ 
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold")
          ,legend.text=element_text(size=25)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="bottom"
          ,legend.direction = "horizontal"
          #,legend.justification = "vertical"
          ,panel.grid = element_blank()
          ,plot.margin=unit(c(0.2, 0.8, 0.2, 0.2),"cm"))+
    scale_colour_manual(name = element_blank(),
                        values = colssss,
                        breaks = c(col.b2, col.hcr , col.gaug, col.event),
                        guide = guide_legend(override.aes = list(
                          linetype = c("solid", "solid",  "blank", "blank"),
                          shape = c(NA, NA, 21, 19),
                          size =c(1, 1, 5, 1.5))))+
    scale_linetype_manual(name=element_blank(), 
                          values=colssss2,
                          guide = guide_legend(override.aes = list( 
                            col =c("blue"))))
  
  ggsave(bt_hcrit_interpol_zoom.plot, filename = paste0(dir.sed.transp,"/bt_and_hcrit_interpolat_zoom.png"), 
         width = 18, height =8, dpi = 400)
  return(bt_hcrit_interpol_zoom.plot)
}














####################################################################################################
y_limni.plot =function(dir.sed.transp, gaugings_SPD, df.limni, y_limni, y_Gaug, ycr, ylimits, d50) {
####################################################################################################
  col.ycr <- c("ycr (MAP)");  
  col.gaug   <- c("Gaugings"); 
  colssss <- c(
               "ycr (MAP)" =  "green",
               "Gaugings" = "black")
  y_limni.plot = ggplot() +
    geom_line( aes(x=df.limni$t_limni, y=y_limni),  color="gray60", size=0.3)+
    geom_point(aes(x=gaugings_SPD$t, y = y_Gaug, color =col.gaug), pch=21, fill=gaugings_SPD$color, size=5)+
    geom_hline(aes(yintercept = ycr, color=col.ycr), size = 1.5) +
    scale_y_continuous(name = expression("Water depth y [m]"),  limits =c(ylimits[1], ylimits[2] )) + 
    scale_x_continuous(name = expression("Time [days]"), expand = c(0,0))+
    coord_cartesian(clip = 'off')+
    #geom_point(aes(x=unlist(event.t), y=unlist(event.h)), color= "black", size=1)+
    annotate("text", x=0, y=ycr + 0.3,
             label=paste0("Critical ycr with d = ",d50, " m"), color = "green", size=7, hjust=0) +
    theme_bw(base_size=20)+ 
    theme(axis.text=element_text(size=20),axis.title=element_text(size=30,face="bold") 
          ,legend.text=element_text(size=20)
          ,legend.title=element_text(size=30)
          ,legend.key.size=unit(1.5, "cm")
          ,legend.position="bottom"
          ,legend.direction = "horizontal"
          ,panel.grid = element_blank()
          ,plot.margin=unit(c(0.2, 0.8, 0.2, 0.2),"cm"))+
    scale_colour_manual(name = element_blank(),
                        values = colssss,
                        breaks = c(col.ycr ,
                                   col.gaug),
                        guide = guide_legend(override.aes = list(
                          linetype = c("solid", "blank"),
                          shape = c(NA, 21),
                          size =c(1.5, 5))))
  ggsave(y_limni.plot, filename = paste0(dir.sed.transp,"/ylimni_d",d50,".png"), 
         width = 18, height =8, dpi = 400)
  return(y_limni.plot)
}





































##############################################################################################
qs.plot = function(dir.sed.transp.d50, 
                   df.limni, y_limni, 
                   qs, ycr, ts, 
                   t.event, y.event,
                   start.event.new, end.event.new, 
                   ylimits.st, d50,
                   gaugings_SPD, y_Gaug, 
                   qscum.tot) {
##############################################################################################
  df.event.effect=NULL; df.event.potent=NULL;
  #
  qs.plot = ggplot() + 
    geom_line( aes(x=df.limni$t_limni, y = y_limni), color="gray70", size=0.3)+
    geom_point(aes(x=gaugings_SPD$t, y = y_Gaug), pch=21, fill="gray60", size=3)+
    # geom_segment(aes(x = df.limni$t_limni, 
    #                  xend =df.limni$t_limni, 
    #                  y=0, 
    #                 yend = qs/1000), color="blue", size= 1, alpha = 0.5)+
    geom_hline(yintercept = ycr, color = "green", size=1.5)
  #**************************************************************************************
  # for (kkkk in 1:length(start.event.new)) {     # ALL POTENTIAL EVENTS
  #   df.event.potent[[kkkk]] = data.frame(t=df.limni$t_limni[start.event.new[[kkkk]] : end.event.new[[kkkk]]],
  #                                        y=y_limni[start.event.new[[kkkk]] : end.event.new[[kkkk]]])
  #   qs.plot = qs.plot + 
  #     geom_point(data= df.event.potent[[kkkk]], aes(x=t, y=y), color ="green", size = 0.4)
  # }
  #**************************************************************************************
  #geom_point(aes(x= ts$t.adj[sequence.ST], y = 0), color ="red", shape =4, size=4, stroke=2 )+
  qs.plot = qs.plot +
    coord_cartesian(clip = 'off')+
    annotate("text", x=50, y=ycr +0.2, 
             label=paste0("Critical water depth ycr with d = ",d50,"m")
             ,color = "green", size=5, hjust=0) +
    scale_y_continuous(name = expression("Water depth y [m]"), limits =c(ylimits.st[1], ylimits.st[2] ), 
                       expand = c(0,0)) + 
    scale_x_continuous(name = expression("Time [days]"), expand = c(0,0))+
    theme_bw(base_size=20)+ 
    theme( axis.text=element_text(size=20)
           , axis.title=element_text(size=30,face="bold")
           ,legend.text=element_text(size=20)
           ,legend.title=element_text(size=30)
           ,legend.key.size=unit(1.5, "cm")
           ,legend.position="none"
           ,panel.grid = element_blank()
           ,plot.margin=unit(c(0, 3.75, 0.2, 0),"cm")
           ,axis.text.y.right = element_text(size=10, color = "blue")
           ,axis.line.y.right = element_line(color = "blue") 
           ,axis.ticks.y.right = element_line(color = "blue")
           ,axis.title.y.right  = element_text(color = "blue", size=25,
                                               margin(t = 0, r = 0, b = 0, l = 0.5)))
  #**********************************************************************************
  # plot with continuous cumulative sed transport:
  qs.cum.plot = ggplot() + 
    #geom_line(aes(x=df.limni$t_limni, y=qs*100), color="blue", size=1.5) +
    geom_segment(aes(x=df.limni$t_limni, xend=df.limni$t_limni,
                     y= 0, yend = qs*200000), color="blue", size=1)+
    geom_line( aes(x=df.limni$t_limni, y = qscum.tot), color="red", size=1.5)+
    scale_y_continuous(name = "Volume of sediments Vs [m3]",
                       sec.axis = sec_axis(~ . /200000, name = "Sediment transport qs [m3/dt]"), 
                       expand = c(0,0)) + 
    scale_x_continuous(name = "Time [days]", 
                       expand = c(0,0))+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+ 
    theme( axis.text=element_text(size=20)
           , axis.title=element_text(size=25,face="bold")
           ,legend.text=element_text(size=20)
           ,legend.title=element_text(size=30)
           ,legend.key.size=unit(1.5, "cm")
           ,legend.position="none"
           ,panel.grid = element_blank()
           ,plot.margin=unit(c(0.8, 0, 0.2, 0.2),"cm")
           ,axis.text.y.right = element_text(size=10, color = "blue")
           ,axis.line.y.right = element_line(color = "blue") 
           ,axis.ticks.y.right = element_line(color = "blue")
           ,axis.title.y.right  = element_text(color = "blue", 
                                               size=25, margin(t = 0, r = 0, b = 0, l = 0)))
  
  
  
  
  
  #***********************************************************************************
  # annotation plot:
  qs.plot2 =  ggplot()+
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-4,0))+
    scale_x_continuous(name ="", expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1))) +
    coord_cartesian(clip = 'off')+
    #theme_bw(base_size = 20)+
    theme_set(theme_gray(base_size = 20))+
    theme( 
      axis.text=element_text(size=15)
      ,axis.title=element_text(size=20,face="bold")
      ,panel.grid.major=element_blank()
      ,panel.grid.minor=element_blank()
      ,legend.text=element_text(size=20)
      ,legend.title=element_text(size=30)
      ,legend.key.size=unit(1.5, "cm")
      ,legend.position="none"
      ,plot.margin=unit(c(0, 3.75, 0.2, 1.25),"cm")
      ,axis.text.x =element_blank()
      ,axis.text.y =element_blank()
      ,axis.ticks.y = element_blank()
      ,axis.ticks.x = element_blank())+
    annotate("text", x=10, y= -0.5, label=paste0("Shift times obtained from gaugings and stage-recessions: ")
             ,color = "violet", size=7, hjust = 0) +
    annotate("text", x=10, y= - 2.5, label=paste0("Selected reference morphogenic events:")
             ,color = "red", size=7, hjust=0) 
  
  #*******************************************************************
  # ALL EFFECTIVE EVENTS from step 1 of the retro analysis
  qs.plot2 = qs.plot2 + 
    geom_segment(aes(x    = ts$t.adj,
                     y    = rep(-2,length(ts$t.adj)), 
                     yend = rep(-1, length(ts$t.adj)), 
                     xend = ts$t.adj), 
                 size = 1, color="violet")
  #*******************************************************************
  for (jjjj in 1:length(t.event)) {      # SELECTED EFFECTIVE EVENTS from step 1 of the retro analysis
    df.event.effect[[jjjj]] = data.frame(t=t.event[[jjjj]], y=y.event[[jjjj]])
    qs.plot2 = qs.plot2 + 
      geom_segment(data= df.event.effect[[jjjj]] , 
                   aes(x=t[1], y =-4, yend =-3, xend=t[1]), 
                   size = 1, color="red")
  }
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  df.event.potent = NULL
  for (kkkk in 1:length(start.event.new)) {       #ALL POTENTIAL EVENTS
    df.event.potent[[kkkk]] = data.frame(t=df.limni$t_limni[start.event.new[[kkkk]] : end.event.new[[kkkk]]],
                                         y=y_limni[start.event.new[[kkkk]] : end.event.new[[kkkk]]])
    qs.plot = qs.plot +
              geom_point(data= df.event.potent[[kkkk]], aes(x=t, y=y), color ="blue", size = 1)
  }
  
  #**********************************************************
  qs.plot2 = qs.plot2 + 
    geom_hline(yintercept = c(-4,-3,-2,-1), color="darkgray", linetype="dashed", size = 0.3)
  #**********************************************************
  qs.plot3 = plot_grid( qs.plot, qs.cum.plot, qs.plot2, 
                        ncol = 1, nrow = 3, rel_heights = c(1,0.7, 0.25))
  #**********************************************************
  ggsave(qs.plot3, filename = paste0(dir.sed.transp.d50,"/qs_d",d50,".png") ,
         width = 16, height =16, dpi = 400)
  return(qs.plot)
}































########################################################################################################
qsc.t.plot = function(dir.sed.transp.d50, df.limni, y_limni, 
                      ylimits.st, hlimits.st,
                      t_Gaug, y_Gaug,
                      t.event, y.event, h.event,
                      
                      phi1, ycr1, hcr1, qsc_crit1, qs1,
                      start.event.new1, end.event.new1,
                      finalqs.cum1, 
                      
                      start.new1, end.new1,
                      
                      phi2, hcr2,
                      start.event.new2, end.event.new2,
  
                      phi3, hcr3,
                      start.event.new3, end.event.new3,
    
                      plot.event.index) {
########################################################################################################
  # Initialisation:
  df.event.effect1 =NULL; df.event.potent1=NULL; df.event.new1 =NULL;
  df.event.effect2 =NULL; df.event.potent2=NULL; df.event.new2 =NULL;
  df.event.effect3 =NULL; df.event.potent3=NULL; df.event.new3 =NULL;
  
  col.ycr.1 <- c("yc"); 
  col.hcr.1 <- c("hc (d/S0 = 12.5 m)"); 
  col.hcr.2 <- c("hc (d/S0 = 20 m)"); 
  col.hcr.3 <- c("hc (d/S0 = 30 m)"); 
  col.gaug.1   <- c("Gaugings"); 
  col.event.1 <- c("Data of the morphogenic flood (d/S0 = 12.5 m)")
  
  # colssss <- c(
  #   "yc"          = "green",
  #   "Gaugings"                      = "black",
  #   "Data of the morphogenic flood" = "green")
  
  colssss <- c(
               "hc (d/S0 = 12.5 m)"          = "green",
               "hc (d/S0 = 20 m)"            = "green",
               "hc (d/S0 = 30 m)"            = "green",
               "Gaugings"                    = "black")
               #"Data of the morphogenic flood (d/S0 = 12.5 m)" = "black")
  qs.cum.plot = ggplot() + 
    #geom_line( aes(x=df.limni$t_limni, y = y_limni), color="gray80", size=0.5)+
    #geom_point(aes(x=gaugings_SPD$t, y = y_Gaug, color = col.gaug), pch=21, fill="black", size=3)+    
    # annotate("text", x=50, y=ycr +0.2, 
    #          label=paste0("Critical water depth ycr with d = ",d50,"m")
    #          ,color = "green", size=7, hjust=0) +
    #geom_hline(aes(yintercept = ycr, color = col.ycr), size=1.5) +
    # scale_y_continuous(name = "Water depth y", 
    #                    limits =c(ylimits.st[1], ylimits.st[2]), 
    #                    expand = c(0,0)) + 
    geom_line( aes(x=df.limni$t_limni, y = df.limni$h_limni), color="gray50", size=0.5)+
    geom_point(aes(x=gaugings_SPD$t, y = h_Gaug, color = col.gaug), pch=21, fill="gray30", size=3)+
    geom_line( aes(x= df.limni$t_limni, y = hcr1, color = col.hcr.1), size =1.5, linetype ="solid") +
    geom_line( aes(x= df.limni$t_limni, y = hcr2, color = col.hcr.2), size =1.5, linetype ="dashed") +
    geom_line( aes(x= df.limni$t_limni, y = hcr3, color = col.hcr.3), size =1.5, linetype ="dotted") +
    coord_cartesian(clip = 'off')+
    # annotate("text", x=50, y=hcr[1] + 0.2, 
    #          label=paste0("Triggering stage hc(t)")
    #          ,color = "green", size=8, hjust=0) +
    scale_y_continuous(name = "Stage h [m]", 
                       limits =c(hlimits.st[1], hlimits.st[2]), 
                       expand = c(0,0)) + 
    scale_x_continuous(name = "Time (days)", expand = c(0,0))+
    theme_bw(base_size=20)+ 
    theme( axis.text=element_text(size=20)
           ,axis.text.x = element_blank()
           ,axis.title=element_text(size=25)
           ,legend.text=element_text(size=20)
           ,legend.title=element_text(size=30)
           ,legend.key.size=unit(1.5, "cm")
           ,legend.position="top"
           ,legend.direction ="horizontal"
           ,panel.grid = element_blank()
           ,plot.margin = unit(c(0, 4.5, 0.5, 0),"cm")
           ,axis.title.y  = element_text(margin(t = 0, r = 0, b = 0, l = 0))
           ,axis.title.x = element_blank())

  #**************************************************************************************
  # plot with continuous cumulative sed transport:
  #**************************************************************************************
  # col.qsc_crit <- c("Vsc* [m3]");  
  # col.cum.event <- c("Volume of sed. of the mrophogenic flood")
  # col.cum.tot   <- c( "Cumulative volume of sed. from t=0"); 
  # col.qst<- c( "Istantaneous sed. transport flux, qs(t)")
  # colssss.1 <- c(
  #   "Vsc* [m3]" =  "red",
  #    # "Volume of sed. of the mrophogenic flood" = "red",
  #   "Cumulative volume of sed. from t=0" = "red",
  #   "Istantaneous sed. transport flux, qs(t)" = "blue")
  # 
  qs.cum.plot1 = ggplot() + 
    #geom_line(aes(x=df.limni$t_limni, y=qs*100), color="blue", size=1.5) +
    geom_segment(aes(x=df.limni$t_limni, xend=df.limni$t_limni,
                     y= 0, yend = qs*200000), color= "blue", size=1)+
    geom_line( aes(x=df.limni$t_limni, y = qscum.tot), color= "red", size=1.5)+
    # geom_hline(aes(yintercept = ycr, color = col.qsc_crit), size=1.5) +
    # annotate("text", x=50, y=qsc_crit +0.2, 
    #          label=paste0("Critical volume of sediments Vsc* [m3]")
    #          ,color = "red", size=7, hjust=0) +
    scale_y_continuous(name = "Cumulative volume of sed. from t=0",
                       sec.axis = sec_axis(~ . /200000, name = "Sed. transport flux qs(t) [m3/s]"), 
                       expand = c(0,0)) + 
    scale_x_continuous(name = "", 
                       expand = c(0,0))+
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+ 
    theme( axis.text=element_text(size=20)
           ,axis.title=element_text(size=25)
           ,legend.text=element_text(size=20)
           ,legend.title=element_text(size=30)
           ,legend.key.size=unit(1.5, "cm")
           ,legend.position="top"
           ,panel.grid = element_blank()
           ,plot.margin=unit(c(1, 0, 1.5, 0.2),"cm")
           #,axis.title.x = element_blank()
           #,axis.text.x = element_blank()
           
           ,axis.ticks.y = element_line(color = "red")
           ,axis.title.y  = element_text(color = "red")
           ,axis.text.y = element_text(size=25, color="red")
           ,axis.line.y = element_line(color = "red") 
           
           ,axis.text.y.right = element_text(size=20, color = "blue")
           ,axis.line.y.right = element_line(color = "blue") 
           ,axis.ticks.y.right = element_line(color = "blue")
           ,axis.title.y.right  = element_text(color = "blue", size=25, margin(t = 0, r = 0, b = 0, l = 0)))
    # scale_colour_manual(name = element_blank(),
    #                     values = colssss.1,
    #                     breaks = c(col.qsc_crit , col.cum.event, col.cum.tot, col.qst),
    #                     guide = guide_legend(override.aes = list(
    #                       linetype = c("solid",  "solid", "solid"),
    #                       shape = c(NA, NA, NA),
    #                       size =c(1.5, 1.5, 1),
    #                       alpha= c(1,1,1))))
  ggsave(qs.cum.plot1, filename = paste0(dir.sed.transp.d50,"/plot_qs_cum_d",d50,".png"),
         width = 16, height =8, dpi = 400)
  #***********************************************************************************
  # annotation plot:
  qs.cum.plot2 =  ggplot()+
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-8,0))+
    scale_x_continuous(name ="Time [days]", expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1))) +
    coord_cartesian(clip = 'off')+
    #theme_bw(base_size = 20)+
    theme_set(theme_gray(base_size = 20))+
    theme( 
      axis.text=element_text(size=20)
      ,axis.title=element_text(size=25)
      ,panel.grid.major=element_blank()
      ,panel.grid.minor=element_blank()
      ,legend.text=element_text(size=20)
      ,legend.title=element_text(size=30)
      ,legend.key.size=unit(1.5, "cm")
      ,legend.position="none"
      ,plot.margin=unit(c(1, 4.5, 0.3, 2.5),"cm")
      #,axis.text.x =element_blank()
      ,axis.title.x = element_text(size=25, margin(t = 2, r = 0, b = 0, l = 0))
      ,axis.text.y =element_blank()
      ,axis.ticks.y = element_blank())+
      #,axis.ticks.x = element_blank())+
    annotate("text", x=10, y= -0.5, label=paste0("Reference shift times (from gaugings and stage-recessions): ")
             ,color = "red", size=6, hjust = 0) +
    # annotate("text", x=10, y= - 2.5, label=paste0("Selected reference morphogenic events:")
    #          ,color = "red", size=7, hjust=0) +
    # annotate("text", x=0, y= - 4.5, label="All potential morphogenic events (y > ycr):"
    #          ,color = "blue", size=7, hjust=0)
  annotate("text", x=0, y= - 2.5, label=paste0("All potential morphogenic events (d/S0 =",phi1,"m):")
           ,color = "blue", size=6, hjust=0) +
  annotate("text", x=0, y= - 4.5, label=paste0("All potential morphogenic events (d/S0 =",phi2,"m):")
           ,color = "blue", size=6, hjust=0)+
  annotate("text", x=0, y= - 6.5, label=paste0("All potential morphogenic events (d/S0 =",phi3,"m):")
           ,color = "blue", size=6, hjust=0)
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  

  # ALL EFFECTIVE EVENTS from step 1 of the retro analysis
  # qs.cum.plot2 = qs.cum.plot2 + 
  #   geom_segment(aes(x    = ts$t.adj,
  #                    y    = rep(-2,length(ts$t.adj)), 
  #                    yend = rep(-1, length(ts$t.adj)), 
  #                    xend = ts$t.adj), 
  #                size = 1, color="violet")
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (jjjj in 1:length(t.event)) {        # SELECTED EFFECTIVE EVENTS from step 1 of the retro analysis
    df.event.effect[[jjjj]] = data.frame(t=t.event[[jjjj]], 
                                         y=y.event[[jjjj]],
                                         h=h.event[[jjjj]])
    qs.cum.plot2 = qs.cum.plot2 + 
      geom_segment(data= df.event.effect[[jjjj]] , 
                   aes(x=t[1], 
                       y =-2, 
                       yend =-1, 
                       xend=t[1]), 
                   size = 1, color="red")
    qs.cum.plot = qs.cum.plot+
      geom_point(data= df.event.effect[[jjjj]] , 
                   aes(x=t[1], 
                       y =max(h)),
                   size = 10, color="red", pch=21, fill="NA")
  }
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (kkkk in 1:length(start.event.new1)) {       #ALL POTENTIAL EVENTS
    df.event.potent1[[kkkk]] = data.frame(t = df.limni$t_limni[start.event.new1[[kkkk]] : end.event.new1[[kkkk]]],
                                          y = y_limni[start.event.new1[[kkkk]] : end.event.new1[[kkkk]]],
                                          h = h_limni[start.event.new1[[kkkk]] : end.event.new1[[kkkk]]])
    # qs.cum.plot = qs.cum.plot +
    #   geom_point(data= df.event.potent[[kkkk]], aes(x=t, y=y, color= col.event), pch =19, fill ="black", size = 0.5)

    # qs.cum.plot = qs.cum.plot +
    #   geom_point(data= df.event.potent[[kkkk]], aes(x=t, y=h, color= col.event), pch =19, fill ="black", size = 0.5)
    # 
    # if (plot.event.index == TRUE) {
    #   qs.cum.plot = qs.cum.plot + 
    #     annotate("text", x= df.limni$t_limni[start.event.new[[kkkk]]], 
    #                      y= y_limni[start.event.new[[kkkk]]], 
    #                      label=kkkk , color = "blue", size=3, hjust=0)
    # }
    qs.cum.plot2 = qs.cum.plot2 + 
      geom_segment(data= df.event.potent1[[kkkk]], 
                   aes(x=t[1], y =-4, yend =-3, xend=t[1]), 
                   size = 1, color="blue")
  }
  
  for (kkkk in 1:length(start.event.new2)) {       #ALL POTENTIAL EVENTS
    df.event.potent2[[kkkk]] = data.frame(t = df.limni$t_limni[start.event.new2[[kkkk]] : end.event.new2[[kkkk]]],
                                          y = y_limni[start.event.new2[[kkkk]] : end.event.new2[[kkkk]]],
                                          h = h_limni[start.event.new2[[kkkk]] : end.event.new2[[kkkk]]])
    # qs.cum.plot = qs.cum.plot +
    #   geom_point(data= df.event.potent[[kkkk]], aes(x=t, y=y, color= col.event), pch =19, fill ="black", size = 0.5)
    
    # qs.cum.plot = qs.cum.plot +
    #   geom_point(data= df.event.potent[[kkkk]], aes(x=t, y=h, color= col.event), pch =19, fill ="black", size = 0.5)
    # 
    # if (plot.event.index == TRUE) {
    #   qs.cum.plot = qs.cum.plot + 
    #     annotate("text", x= df.limni$t_limni[start.event.new[[kkkk]]], 
    #                      y= y_limni[start.event.new[[kkkk]]], 
    #                      label=kkkk , color = "blue", size=3, hjust=0)
    # }
    qs.cum.plot2 = qs.cum.plot2 + 
      geom_segment(data= df.event.potent2[[kkkk]], 
                   aes(x=t[1], y =-6, yend =-5, xend=t[1]), 
                   size = 1, color="blue")
  }
  
  for (kkkk in 1:length(start.event.new3)) {       #ALL POTENTIAL EVENTS
    df.event.potent3[[kkkk]] = data.frame(t = df.limni$t_limni[start.event.new3[[kkkk]] : end.event.new3[[kkkk]]],
                                          y = y_limni[start.event.new3[[kkkk]] : end.event.new3[[kkkk]]],
                                          h = h_limni[start.event.new3[[kkkk]] : end.event.new3[[kkkk]]])
    # qs.cum.plot = qs.cum.plot +
    #   geom_point(data= df.event.potent[[kkkk]], aes(x=t, y=y, color= col.event), pch =19, fill ="black", size = 0.5)
    
    # qs.cum.plot = qs.cum.plot +
    #   geom_point(data= df.event.potent[[kkkk]], aes(x=t, y=h, color= col.event), pch =19, fill ="black", size = 0.5)
    # 
    # if (plot.event.index == TRUE) {
    #   qs.cum.plot = qs.cum.plot + 
    #     annotate("text", x= df.limni$t_limni[start.event.new[[kkkk]]], 
    #                      y= y_limni[start.event.new[[kkkk]]], 
    #                      label=kkkk , color = "blue", size=3, hjust=0)
    # }
    qs.cum.plot2 = qs.cum.plot2 + 
      geom_segment(data= df.event.potent3[[kkkk]], 
                   aes(x=t[1], y =-8, yend =-7, xend=t[1]), 
                   size = 1, color="blue")
  }
  # qs.cum.plot = qs.cum.plot +
  # scale_colour_manual(name = element_blank(),
  #                     values = colssss,
  #                     breaks = c(col.hcr , col.gaug, col.event),
  #                     #breaks = c(col.ycr , col.gaug, col.event),
  #                     guide = guide_legend(override.aes = list(
  #                       linetype = c("solid",  "blank", "blank"),
  #                       shape = c(NA, 21, 19),
  #                       size =c(1.5, 3, 0.5),
  #                       fill = c(NA, "gray", "green"))))
  qs.cum.plot = qs.cum.plot +
  scale_colour_manual(name = element_blank(),
                      values = colssss,
                      breaks = c(col.hcr.1 , col.hcr.2, col.hcr.3, col.gaug),
                      labels = c(paste0("hc (d/S0 = 12.5 m)"),
                                 paste0("hc (d/S0 = 17.5 m)"),
                                 paste0("hc (d/S0 = 32.5 m)"),
                                 "Gaugings"),
                      #breaks = c(col.ycr , col.gaug, col.event),
                      guide = guide_legend(override.aes = list(
                        linetype = c("solid",  "dashed", "dotted", "blank"),
                        shape = c(NA, NA, NA, 19),
                        size =c(1.5, 1.5, 1.5, 3),
                        fill = c(NA, NA, NA, "gray30"))))
  
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # qs.cum.plot2 = qs.cum.plot2 + 
  #   geom_hline(yintercept = c(-4,-3,-2,-1,0),
  #              color="darkgray", linetype="dashed", size = 0.3)
  qs.cum.plot2 = qs.cum.plot2 +
    geom_hline(yintercept = c(-8,-7,-6, -5,-4,-3,-2,-1,0),
               color="darkgray", linetype="dashed", size = 0.3)
  #*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # qs.cum.plot3 = plot_grid( qs.cum.plot,
  #                           qs.cum.plot1,
  #                           qs.cum.plot2, 
  #                           ncol = 1, nrow = 3, rel_heights = c(1, 0.5, 0.5))
  # 
  qs.cum.plot4 = plot_grid( qs.cum.plot,
                            qs.cum.plot2, 
                            ncol = 1, nrow = 2, rel_heights = c(1, 0.8))

  #save:
  # ggsave(qs.cum.plot3, filename = paste0(dir.sed.transp.d50,"/qs_cum_d",d50,".png"),
  #        width = 16, height =16, dpi = 400)
  
  ggsave(qs.cum.plot4, filename = paste0(dir.sed.transp.d50,"/detection_comparison.png"),
         width = 16, height =11, dpi = 400)

  # save pdf:
  # pdf(paste0(dir.sed.transp.d50,"/qs_cum_d",d50,".pdf"), 16, 16 ,useDingbats=F)
  # print(qs.cum.plot3)
  # dev.off()
  
  ###########################################
  return()
}



































#version with only one choice:
######################################################################################################################
qsc.t.plot.one.choice = function(dir.sed.transp.d50,
                                 df.limni.ST, 
                                 y_limni, 
                                 ylimits.st, hlimits.st, 
                                 gaugings_SPD,
                                 t.event1, y.event1, h.event1,
                                 phi1, hcr1, ycr1, qsc_crit1, qs1,
                                 
                                 start.event.index1, end.event.index1,
                                 start.event.new111, end.event.new111,
                                 finalqs.cum1, cum.qs1,
                                 start.new111, end.new111,
                                 
                                 qscum.tot1,
                                 plot.event.index,
                                 V.limits,
                                 V.for.text) {
  ######################################################################################################################
# Initialisation:
  df.event.effect  = NULL; 
  df.event.potent1 = NULL; 
  df.event.new1    = NULL;
  col.ycr.1        = c("yc"); 
  col.hcr.1        = c("hc"); 
  col.gaug.1       = c("Gaugings"); 
  col.event.1      = c("Data of the morphogenic flood")
  colssss          = c("hc"= "green",  "Gaugings" = "black")
  min.stage        = min(df.limni.ST$h_limni)
  
  qs.cum.plot = ggplot() + 
                  geom_line( aes(x=df.limni.ST$t_limni,   y = df.limni.ST$h_limni),     color="gray80", size=0.3)+
                  geom_point(aes(x=gaugings_SPD$t,        y = gaugings_SPD$h,           color = col.gaug.1), pch=21, 
                             fill="gray30", size=3)+
                  geom_line( aes(x= df.limni.ST$t_limni,  y = hcr1, color = col.hcr.1), size =1.5, linetype ="solid")+
                  coord_cartesian(clip = 'off') +
                  scale_y_continuous(name = "Stage h [m]", limits = c(min.stage, hlimits.st[2]), expand = c(0,0))+ 
                  scale_x_continuous(name = "Time (days)", limits = c(df.limni.ST$t_limni[1], tail(df.limni.ST$t_limni,1)), expand = c(0,0))+
                  theme_bw(base_size=20)+ 
                  theme( axis.text         = element_text(size=20)
                         ,axis.text.x      = element_blank()
                         ,axis.title       = element_text(size=25)
                         ,legend.text      = element_text(size=20)
                         ,legend.title     = element_text(size=30)
                         ,legend.key.size  = unit(1.5, "cm")
                         ,legend.position  ="top"
                         ,legend.direction ="horizontal"
                         ,panel.grid       = element_blank()
                         ,plot.margin      = unit(c(0, 2.5, 0.5, 0.5),"cm")
                         ,axis.title.y     = element_text(margin = margin(t = 0, r = 80, b = 0, l = 0))
                         ,axis.title.x     = element_blank())
   

  # plot the cumulative sed transport V(t):
  ##################################################################################################################################
  nevents.st           = length(start.event.index1)
  #start.event.new111   = start.event.new111[-1]
  #start.new111         = start.new111[-1]
  # plot the time series of V (cumulative sed transport of each event):
  qs.cum.plot1 = ggplot() +
                 geom_segment(aes(x     = df.limni.ST$t_limni[start.event.new111],
                                  xend  = df.limni.ST$t_limni[start.event.new111],
                                  y     = 0,
                                  yend  = finalqs.cum1), color="blue", size= 2, alpha = 1)+
                   annotate("rect",
                            xmin =  df.limni.ST$t_limni[start.event.index1],
                            xmax =  df.limni.ST$t_limni[end.event.index1],
                            ymin =  0,
                            ymax =  Inf,
                            fill = "pink", alpha =0.2)+
                   geom_vline(xintercept =  df.limni.ST$t_limni[start.event.index1], color = "red", linetype="dashed", size= 0.2)+
                   geom_vline(xintercept =  df.limni.ST$t_limni[end.event.index1],
                              color = "red", linetype="dashed", size = 0.2) +
                   annotate("text", 
                            x= (df.limni.ST$t_limni[start.event.index1] + 
                                  (df.limni.ST$t_limni[end.event.index1] - df.limni.ST$t_limni[start.event.index1])/2), 
                            y= V.for.text, 
                            label=seq(1,nevents.st, 1),   
                            color = "red", size=7) +
                   scale_y_continuous(name = TeX("$V \\; for \\; each \\; event \\; \\left[ m^3 \\right]$"),
                                      expand = c(0,0), limits= c(V.limits[1],V.limits[2])) +
                   scale_x_continuous(name = "", limits = c(df.limni.ST$t_limni[1], tail(df.limni.ST$t_limni,1)), expand = c(0,0))+       
                   coord_cartesian(clip = 'off')+
                   theme_bw(base_size=20)+
                   theme( axis.text=element_text(size=20)
                          ,axis.title=element_text(size=25)
                          ,legend.text=element_text(size=20)
                          ,legend.title=element_text(size=30)
                          ,legend.key.size=unit(1.5, "cm")
                          ,legend.position="top"
                          ,panel.grid    = element_blank()
                          ,plot.margin   = unit(c(0, 2.5, 0.5, 0.5),"cm")
                          ,axis.ticks.y  = element_line(color = "blue")
                          ,axis.title.y  = element_text(color = "blue", margin = margin(t = 0, r = 25, b = 0, l = 0))
                          ,axis.text.y   = element_text(size=16, color="blue")
                          ,axis.line.y   = element_line(color = "blue"))
  
  
  ######################################################################################################################
  # annotation plot:
  # will plots the detected shift times (the reference times and the potential ones).
  qs.cum.plot2 =  ggplot()+
                  scale_y_continuous(name = "", expand = c(0,0), limits = c(-4,0))+
                  scale_x_continuous(name = "Time (days)", limits = c(df.limni.ST$t_limni[1], tail(df.limni.ST$t_limni,1)), expand = c(0,0))+          
                  coord_cartesian(clip = 'off')+
                    #theme_bw(base_size = 20)+
                    theme_set(theme_gray(base_size = 20))+
                    theme( 
                      axis.text         = element_text(size=20)
                      ,axis.title       = element_text(size=25)
                      ,panel.grid.major = element_blank()
                      ,panel.grid.minor = element_blank()
                      ,legend.text      = element_text(size=20)
                      ,legend.title     = element_text(size=30)
                      ,legend.key.size  = unit(1.5, "cm")
                      ,legend.position  = "none"
                      ,plot.margin      = unit(c(1, 2.5, 0.5, 2.5),"cm")
                      ,axis.title.x    = element_text(size=25, margin = margin(t = 2, r = 10, b = 0, l = 0))
                      ,axis.text.y      = element_blank()
                      ,axis.ticks.y     = element_blank())+
                    annotate("text", x= 10, y= -0.5, 
                             label  = paste0("Reference shift times (e.g. from gaugings and stage-recessions): ")
                             ,color = "red", size=6, hjust = 0) +
                    annotate("text", x=0, y= - 2.5, label=paste0("All potential morphogenic events:"), color = "blue", size=6, hjust=0)
    

  ######################################################################################################################
  for (jjjj in 1:length(t.event1)) {        # SELECTED EFFECTIVE EVENTS from step 1 of the retro analysis
    df.event.effect[[jjjj]] = data.frame(t=t.event1[[jjjj]], 
                                         y=y.event1[[jjjj]],
                                         h=h.event1[[jjjj]],
                                         y=y.event1[[jjjj]])
    # update annotation plot with reference events:
    qs.cum.plot2 = qs.cum.plot2 + 
                   geom_segment(data= df.event.effect[[jjjj]],  aes(x=t[1], y =-2, yend =-1, xend=t[1]), size = 1, color="red")
    # update plot with reference events on stage record with circles:
    qs.cum.plot = qs.cum.plot+
                  geom_point(data= df.event.effect[[jjjj]], aes(x=t[1], y =max(h)),size = 10, color="red", pch=21, fill="NA")
  }
  

  # ADD ALL POTENTIAL EVENTS:
  ######################################################################################################################
  for (kkkk in 1:length(start.new111)) {        
    #print(paste0("Potent. Event n. ", kkkk))
    df.event.potent1[[kkkk]] = data.frame(t = df.limni.ST$t_limni[start.new111[[kkkk]] : end.new111[[kkkk]]],
                                          y = y_limni[start.new111[[kkkk]] : end.new111[[kkkk]]],
                                          h = h_limni[start.new111[[kkkk]] : end.new111[[kkkk]]])

    # update annotation plot with the potential shift events:
    qs.cum.plot2 = qs.cum.plot2 + 
                   geom_segment(data= df.event.potent1[[kkkk]],  aes(x=t[1], y =-4, yend =-3, xend=t[1]), 
                                size = 1, color="blue")
  }
  
  ######################################################################################################################
  # add the legend scale to plot 'qs.cum.plot' with the stage record:
  qs.cum.plot = qs.cum.plot +
                scale_colour_manual(name   = element_blank(),
                        values = colssss,
                        breaks = c(col.hcr.1, col.gaug.1),
                        labels = c("hc", "Gaugings"),
                        # breaks = c(col.ycr , col.gaug, col.event),
                        guide  = guide_legend(override.aes = list(
                                 linetype = c("solid", "blank"),
                                 shape = c(NA,  19), size =c(1.5,  3), fill = c(NA,  "gray50"))))
  
  ######################################################################################################################
  # Update annotation plot with separation lines:
  qs.cum.plot2 = qs.cum.plot2 +
                 geom_hline(yintercept = c(-4,-3,-2,-1, 0),
                 color="darkgray", linetype="dashed", size = 0.3)
  
  ######################################################################################################################
  # Combine plots 'qs.cum.plot' (stage record with reference events + gaugings) and  qsc (all events):
  qs.cum.plot1.bis = plot_grid( qs.cum.plot,
                                qs.cum.plot1, 
                                ncol = 1, nrow = 2,
                                rel_heights = c(1, 1))
  # save:
  pdf(paste0(dir.sed.transp.d50,"/plot_V_events.pdf"), 18, 8 ,useDingbats=F)
  print(qs.cum.plot1.bis)
  dev.off()
  
  
  
  
  
  # second plot:
  ######################################################################################################################
  qs.cum.plotbis = qs.cum.plot + theme(axis.title.y = element_text(size=25, margin = margin(t = 0, r = 45, b = 0, l = 0)))
  qs.cum.plot4 = plot_grid( qs.cum.plotbis,
                            qs.cum.plot2,
                            ncol = 1, nrow = 2, rel_heights = c(1, 0.6))
  pdf(paste0(dir.sed.transp.d50,"/ts_potential_all.pdf"), 18, 8 ,useDingbats=F)
  print(qs.cum.plot4)
  dev.off()
  ######################################################################################################################
}








































########################################################################################################
merge.qsc.t.plot = function(dir.sed.transp.d50, 
                            d50,
                            df.limni, y_limni, y.limits, 
                            t_Gaug, y_Gaug, 
                            ycr, 
                            merge.ts,
                            start.event.new, end.event.new,  # all new potential events
                            finalqs.cum,       # qsc of all new potential events
                            t.event, y.event,  # the events defined by step 1 of RA
                            merge.start.event.new, merge.end.event.new,  # merged new potential events
                            merge.finalqs.cum,    # qsc of merged new potential events
                            start.new, end.new, # final new potential events
                            qsc.new,
                            merge.qsc_crit ) { 
########################################################################################################
  # Initialisation:
  df.event.effect =NULL; df.event.potent=NULL; df.merge.event.new=NULL; df.event.new =NULL;
  #
  col.ycr <- c("ycr (MAP)");  
  col.gaug   <- c("Gaugings"); 
  col.event <- c("Data of the morphogenic flood")
  colssss <- c(
    "ycr (MAP)" =  "green",
    "Gaugings" = "black",
    "Data of the morphogenic flood" = "green")
  qs.cum.plot = ggplot() + 
    geom_line( aes(x=df.limni$t_limni, y = y_limni), color="gray80", size=0.5)+
    geom_point(aes(x=gaugings_SPD$t, y = y_Gaug, color = col.gaug), pch=21, fill="black", size=3)+
    geom_hline(aes(yintercept = ycr, color = col.ycr), size=1.5) +
    coord_cartesian(clip = 'off')+
    annotate("text", x=50, y=ycr +0.2, 
             label=paste0("Critical water depth ycr with d = ",d50,"m")
             ,color = "green", size=7, hjust=0) +
    scale_y_continuous(name = "Water depth y", 
                       limits =c(ylimits.st[1], ylimits.st[2] ), 
                       expand = c(0,0)) + 
    scale_x_continuous(name = "Time (days)", expand = c(0,0))+
    theme_bw(base_size=20)+ 
    theme( axis.text=element_text(size=20)
           ,axis.text.x = element_blank()
           ,axis.title=element_text(size=25)
           ,legend.text=element_text(size=20)
           ,legend.title=element_text(size=30)
           ,legend.key.size=unit(1.5, "cm")
           ,legend.position="top"
           ,legend.direction ="horizontal"
           ,panel.grid = element_blank()
           ,plot.margin = unit(c(0, 4.5, 0.5, 0),"cm")
           ,axis.title.y  = element_text(margin(t = 0, r = 0, b = 0, l = 0))
           ,axis.title.x = element_blank())
  
  #**********************************************************************************
  # plot with continuous cumulative sed transport:
  qs.cum.plot1 = ggplot() + 
    #geom_line(aes(x=df.limni$t_limni, y=qs*100), color="blue", size=1.5) +
    geom_segment(aes(x=df.limni$t_limni, xend=df.limni$t_limni,
                     y= 0, yend = qs*200000), color="blue", size=1)+
    geom_line( aes(x=df.limni$t_limni, y = qscum.tot), color="red", size=1.5)+
    scale_y_continuous(name = "Volume of sediments Vs (m3)",
                       sec.axis = sec_axis(~ . /200000, name = "Sediment transport qs (m3/dt)"), 
                       expand = c(0,0)) + 
    scale_x_continuous(name = "", 
                       expand = c(0,0))+
    geom_segment(aes(x = merge.ts, 
                     xend =merge.ts, 
                     y=0, 
                     yend = merge.finalqs.cum), 
                 color="red", size= 6, alpha = 0.4) +
    coord_cartesian(clip = 'off')+
    theme_bw(base_size=20)+ 
    theme( axis.text=element_text(size=20)
           ,axis.title=element_text(size=25)
           ,legend.text=element_text(size=20)
           ,legend.title=element_text(size=30)
           ,legend.key.size=unit(1.5, "cm")
           ,legend.position="none"
           ,panel.grid = element_blank()
           ,plot.margin=unit(c(0.8, 0, 2, 0.2),"cm")
           ,axis.title.x = element_blank()
           ,axis.text.x = element_blank()
           
           ,axis.ticks.y = element_line(color = "red")
           ,axis.title.y  = element_text(color = "red")
           ,axis.text.y = element_text(size=25, color="red")
           ,axis.line.y = element_line(color = "red") 
           
           ,axis.text.y.right = element_text(size=20, color = "blue")
           ,axis.line.y.right = element_line(color = "blue") 
           ,axis.ticks.y.right = element_line(color = "blue")
           ,axis.title.y.right  = element_text(color = "blue", size=25, margin(t = 0, r = 0, b = 0, l = 0)))
  #***********************************************************************************
  # annotation plot:
  qs.cum.plot2 =  ggplot()+ 
    scale_y_continuous(name = "", expand = c(0,0), limits = c(-8,0))+
    scale_x_continuous(name ="Time [days]", expand = c(0,0), limits =c(0,tail(df.limni$t_limni,1))) +
    coord_cartesian(clip = 'off')+
    #theme_bw(base_size = 20)+
    theme_set(theme_gray(base_size = 20))+
    theme( 
      axis.text=element_text(size=20)
      ,axis.title=element_text(size=25)
      ,panel.grid.major=element_blank()
      ,panel.grid.minor=element_blank()
      ,legend.text=element_text(size=20)
      ,legend.title=element_text(size=30)
      ,legend.key.size=unit(1.5, "cm")
      ,legend.position="none"
      ,plot.margin=unit(c(1, 4.5, 0.3, 2),"cm")
      #,axis.text.x =element_blank()
      ,axis.title.x = element_text(size=25, margin(t = 2, r = 0, b = 0, l = 0))
      ,axis.text.y =element_blank()
      ,axis.ticks.y = element_blank()) +
    #,axis.ticks.x = element_blank())+
  annotate("text", x=0, y= -0.5, label=" Shift times obtained from gaugings and stage-recessions: "
           ,color = "violet", size=5, hjust = 0) +
    annotate("text", x=0, y= - 2.5, label="Reference morphogenic events:"
             ,color = "red", size=5, hjust=0) +
    annotate("text", x=0, y= - 4.5, label="All potential morphogenic events (y > ycr):"
             ,color = "blue", size=5, hjust=0) +
    annotate("text", x=0, y= - 6.5, label="Merged potential morphogenic events:"
             ,color = "black", size=5, hjust=0) 
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # ALL EFFECTIVE EVENTS from step 1 of the retro analysis
  qs.cum.plot2 = qs.cum.plot2 + 
    geom_segment(aes(x    = ts$t.adj,
                     y    = rep(-2,length(ts$t.adj)), 
                     yend = rep(-1, length(ts$t.adj)), 
                     xend = ts$t.adj), 
                 size = 1, color="violet")
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (jjjj in 1:length(t.event)) {        # SELECTED EFFECTIVE EVENTS
    df.event.effect[[jjjj]] = data.frame(t=t.event[[jjjj]], y=y.event[[jjjj]])
    qs.cum.plot2 = qs.cum.plot2 + 
      geom_segment(data= df.event.effect[[jjjj]] , 
                   aes(x=t[1],
                       y =-4,
                       yend =-3, 
                       xend=t[1]), 
                   size = 1, color="red")
  }
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for (kkkk in 1:length(start.event.new)) {       #ALL POTENTIAL EVENTS
    df.event.potent[[kkkk]] = data.frame(t=df.limni$t_limni[start.event.new[[kkkk]] : end.event.new[[kkkk]]],
                                         y=y_limni[start.event.new[[kkkk]] : end.event.new[[kkkk]]])
    qs.cum.plot = qs.cum.plot +
      geom_point(data= df.event.potent[[kkkk]], aes(x=t, y=y, color= col.event), pch =19, fill ="green", size = 0.5)
    
    qs.cum.plot2 = qs.cum.plot2 + 
      geom_segment(data= df.event.potent[[kkkk]], 
                   aes(x=t[1], 
                       y =-6, 
                       yend =-5, 
                       xend=t[1]), 
                   size = 1, color="blue")
  }
  
  #######################################################################################################
  # add the merged new potential events
  for (llll in 1:length(merge.start.event.new)) { 
    df.merge.event.new[[llll]] = data.frame( t = df.limni$t_limni[merge.start.event.new[[llll]] : merge.end.event.new[[llll]]],
                                             y = y_limni[merge.start.event.new[[llll]] : merge.end.event.new[[llll]]])
    qs.cum.plot = qs.cum.plot +
      annotate("rect", 
               xmin =  df.merge.event.new[[llll]]$t[1], 
               xmax=   tail( df.merge.event.new[[llll]]$t,1), 
               ymin =  ycr, 
               ymax = Inf, fill = "black", alpha =0.1)
    qs.cum.plot2 = qs.cum.plot2 + 
      annotate("rect",
               xmin =  df.merge.event.new[[llll]]$t[1], 
               xmax=   tail( df.merge.event.new[[llll]]$t,1), 
               ymin =  -8, 
               ymax =  -7, 
               fill = "black", alpha = 0.1)
  }
  qs.cum.plot2 = qs.cum.plot2 + 
    geom_segment(
    aes(x=merge.ts,
        y =-8, 
        yend =-7,
        xend=merge.ts), 
    size = 1, color="black")
  #######################################################################################################
  qs.cum.plot = qs.cum.plot +
                scale_colour_manual(name = element_blank(),
                        values = colssss,
                        breaks = c(col.ycr , col.gaug, col.event),
                        guide = guide_legend(override.aes = list(
                          linetype = c("solid",  "blank", "blank"),
                          shape = c(NA, 21, 19),
                          size =c(1.5, 3, 0.5),
                          fill = c(NA, "black", "green"))))
  #######################################################################################################
  qs.cum.plot2 = qs.cum.plot2 + 
                 geom_hline(yintercept = c(-8,-7,-6,-5,-4,-3,-2,-1,0),
                            color="darkgray", 
                            linetype="dashed",
                            size = 0.3)
  ########################################################################################################
  qs.cum.plot3 = plot_grid( qs.cum.plot,
                            qs.cum.plot1,
                            qs.cum.plot2, 
                            ncol = 1, nrow = 3,
                            rel_heights = c(0.7, 0.5, 0.6))
  
  ########################################################################################################
  #save:
  ggsave(qs.cum.plot3, filename = paste0(dir.sed.transp.d50,"/qs_cum_merged_d",d50,".png"),
         width = 16, height =16, dpi = 400)
  # save pdf:
  # pdf(paste0(dir.sed.transp.d50,"/qs_cum_d",d50,".pdf"), 16, 16 ,useDingbats=F)
  # print(qs.cum.plot3)
  # dev.off()
  
  #######################################################################################################            
   return(list(qs.cum.plot3=qs.cum.plot3, 
               df.event.effect=df.event.effect, 
               df.event.potent=df.event.potent,
               df.merge.event.new=df.merge.event.new, 
               df.event.new=df.event.new))
}


















