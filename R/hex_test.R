library(ggplot2)
library(hexSticker)




# Makes grid data and subplot for hex sticker -----------------------------

# object p must be available for plot codes to follow

# make random data for raster
d <- data.frame(expand.grid(x=1:10, y=1:10))
v <- runif(1:100)
d$v <- v

# make sublot
p <- ggplot() +
  geom_raster(data = d, aes(x=x, y=y, fill=v), alpha = 0.4) +
  scico::scale_fill_scico(palette = "vik") + #vik
  theme_void() + theme_transparent() +
  theme(legend.position = "none")

#--------------------------------------------------------------------------


# Make one plot -----------------------------------------------------------

# parameters - keep hex colours where possible
filename = XXXXXX # filename for plot including filepath (include ".png")
p_color = XXXXXX # colour for package name text
h_color = XXXXXX # colour of hex outline
h_fill = XXXXXX # colour of hex fill


hsticker <- sticker(p, package = "DBCAscatR", filename = filename,
              s_width = 1.4, s_height = 0.7, # these relate to subplot
              s_x = 1, s_y = 0.75, # these tweak location of subplot
              p_size = 22, p_color = p_color,
              h_size = 1, h_color= h_color, h_fill = h_fill)

# NOTE use plot command below to view result - if you don't it will NOT display
# as rendered
plot(hsticker)

#-------------------------------------------------------------------------


# A looped version to try a few things at once ----------------------------

# parameters
# colours to try
cols = c("#00560D", "#528D01", "#206E09", "#3BBE94", "#01A58A", "#01A58A")
# descriptors to add to filename
cnames = c("emerald", "gumdrop", "kelly", "lagoon", "mermaid", "cerulean")
# output filepath
fp = XXXXXX


for(i in seq_along(cols)){
  outfile <- paste0(fp, "/DBCAscatR_logo_", cnames[i], ".png")
  hsticker <- sticker(p, package="DBCAscatR", filename = outfile,
                s_width = 1.4, s_height = 0.7,
                s_x = 1, s_y = 0.75,
                p_size = 22, p_color = cols[i],
                h_size = 1, h_color= cols[i], h_fill = cols[i])
  plot(hsticker)
}

# can use below link to find different colours and tones
#  https://www.w3schools.com/colors/colors_picker.asp

# can add in other named colour vectors to the loop if text and borders etc are
# not the same
