
library(tidyverse)
library(data.table)

extend_meta <- fread('./other_sc_data/output/02-pbmc_concat_harm_anno_meta.csv',nThread = 16) 

ggplot(data = extend_meta, aes(x =  extend_meta$sample, fill =extend_meta$main_type))+
  geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(extend_meta$main_type))) + xlab('')+
  labs(fill="") +  scale_fill_manual(values=get_color(15,set_len = 11,set = 'Set1')) + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pie_df <- extend_meta$group %>% table() %>% data.frame() 
pie_df$. <- factor(pie_df$., levels = c('hc_child','hc','sle','sle_child','sle_flare','sle_treated','IFNbeta_stim'))
pie(x =pie_df$Freq ,labels = pie_df$Freq ,col = c('#377EB8','#4DAF4A','#A65628','#FF7F00','#984EA3','#E41A1C','#FFFF33'),
                                      lty=0,main = 'group',clockwise = T,)

legend("right", x=1.2,y=0.7,
       inset=-0.2, #图例为关键词时，inset = 分数 设置其相对位置(-向下，+向上)
       legend=extend_meta$group %>% table() %>% names(), #图例文字
       fill=c('#377EB8','#4DAF4A','#A65628','#FF7F00','#984EA3','#E41A1C','#FFFF33'), #图例填充颜色
       
       #box.lty=1, #图例最外大方框
       bty="n", #不要图例边框
       
       #title.adj=0, #图例标题的相对位置，0.5为默认，在中间。0最左，1为最右。
       #title="Type",#图例标题
       
       cex=1, #字体大小倍数
       #horiz = T,#横着显示，会覆盖掉ncol
       ncol=1, #图例显示为n列
       #byrow=T, #??
       
       xpd=T, #有这句话才能显示在图外
       x.intersp=0.5, #图例中文字离图片的水平距离
       #y.intersp=1, #图例中文字离图片的垂直距离
       
       text.width=0.3, #两个图例之间的距离
       #merge = TRUE,
       #title="Type", #图例标题
       border=NA #不要图例小方块描边
)
