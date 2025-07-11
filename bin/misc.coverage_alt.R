##alternate plotting format 

contigsize=LRcoverage %>% group_by(contig) %>% summarise(size=max(end), )

ggplot()+geom_rect(data=contigsize, aes(xmin=0, xmax=size, ymin=0.5, ymax=1.5), fill="grey")+
  geom_line(data=LRcoverage, aes(x=((start+end)/2), y=coverage_norm))+
  facet_grid(contig ~ ., switch="y")+
  theme(strip.text.y.left = element_text(angle = 0, hjust=1), axis.text.y=element_blank(),axis.ticks.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line.x = element_line(colour = "black"), strip.background =element_blank(), legend.position = c(0.75, 0.3))+
  ylim(0,3)
