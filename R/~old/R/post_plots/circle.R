# Libraries
library(ggraph)
library(igraph)

# create a data frame giving the hierarchical structure of your individuals
d1=data.frame(from="origin", to=paste("group", seq(1,8), sep=""))
d2=data.frame(from=rep(d1$to, each=1), to=paste("subgroup", seq(1,8), sep="_"))
hierarchy=rbind(d1, d2)

# create a vertices data.frame. One line per object of our hierarchy
vertices = data.frame(name = unique(c(as.character(hierarchy$from),
                      as.character(hierarchy$to))) )

# Create a graph object with the igraph library
mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )
# This is a network object, you visualize it as a network like shown in the network section!

all_leaves=paste("subgroup", seq(1,8), sep="_")

connect=as.data.frame(t(combn(all_leaves, 2)))
colnames(connect) <- c("from", "to")

# The connection object must refer to the ids of the leaves:
from = match( connect$from, vertices$name)
to = match( connect$to, vertices$name)

# plot
p=ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_conn_bundle(data = get_con(from = from, to = to),
                   colour="white", alpha=0.8, width=0.2, tension=0.5) +
  theme_void() +
  theme(legend.position = "none")

# just a blue uniform color. Note that the x*1.05 allows to make a space between the points and the connection ends
p <- p + geom_node_point(aes(filter = leaf, x = x*1.25, y=y*1.25),
                          colour="#fb9a99", size=2)

ggsave("inter_path_pink2.png", width=1, height=1, bg = "transparent")
