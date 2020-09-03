# vim: set noexpandtab tabstop=2:

f=read.table(infile, header=T, sep='\t', quote="", check.names=F, comment.char='#', stringsAsFactors=F)
f$label=factor(f$label, unique(f$label))
pdf(NULL)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
pcount=ggplot(data=f, aes(x=label, y=counts)) +
geom_bar(stat='identity', show.legend=F, color='blue', width=0.5, fill='blue') +
scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
xlab('Breaks') +
ylab('Frequency') +
theme_bw() +
theme(
	plot.background = element_blank()
	, panel.grid.major = element_blank()
	, panel.grid.minor = element_blank()
	, panel.border = element_blank()
	, plot.title = element_text(hjust=0.5)
	, axis.line = element_line(color='black')
	, axis.text = element_text(color='black')
	, legend.position = 'top'
	, legend.title = element_blank()
	, axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=5)
	)
pdensity=ggplot(data=f, aes(x=label, y=density)) +
geom_bar(stat='identity', show.legend=F, color='blue', width=0.5, fill='blue') +
scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
xlab('Breaks') +
ylab('Density') +
theme_bw() +
theme(
	plot.background = element_blank()
	, panel.grid.major = element_blank()
	, panel.grid.minor = element_blank()
	, panel.border = element_blank()
	, plot.title = element_text(hjust=0.5)
	, axis.line = element_line(color='black')
	, axis.text = element_text(color='black')
	, legend.position = 'top'
	, legend.title = element_blank()
	, axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=5)
	)
ppi=ggplot(data=f, aes(x=label, y=pi)) +
geom_bar(stat='identity', show.legend=F, color='blue', width=0.5, fill='blue') +
scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
xlab('Breaks') +
ylab('Probability') +
theme_bw() +
theme(
	plot.background = element_blank()
	, panel.grid.major = element_blank()
	, panel.grid.minor = element_blank()
	, panel.border = element_blank()
	, plot.title = element_text(hjust=0.5)
	, axis.line = element_line(color='black')
	, axis.text = element_text(color='black')
	, legend.position = 'top'
	, legend.title = element_blank()
	, axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=5)
	)
ml=marrangeGrob(list(pcount, pdensity, ppi), nrow=3, ncol=1, top='')
ggsave(ml, file=outfile, width=width, height=height, useDingbats=F)
