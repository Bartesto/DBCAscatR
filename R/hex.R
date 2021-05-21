library(ggplot2)
library(hexSticker)
library(scico)
library(here)

# make random data for raster
d <- data.frame(expand.grid(x=1:10, y=1:10))
v <- runif(1:100)
d$v <- v

# make sublot
p2 <- ggplot() +
  geom_raster(data = d, aes(x=x, y=y, fill=v), alpha = 0.4) +
  scico::scale_fill_scico(palette = "vik") + #vik
  theme_void() + theme_transparent() +
  theme(legend.position = "none")


# output name
cols = c("#00560D", "#528D01", "#206E09", "#3BBE94", "#01A58A", "#01A58A")
cnames = c("emerald", "gumdrop", "kelly", "lagoon", "mermaid", "cerulean")




for(i in seq_along(cols)){
  outfile <- here("man", "figures", paste0("DBCAscatR_logo_", cnames[i], ".png"))
  t2 <- sticker(p2, package="DBCAscatR", filename = outfile,
                s_width = 1.4, s_height = 0.7,
                s_x = 1, s_y = 0.75,
                p_size = 22, p_color = "darkred",
                h_size = 1, h_color= "darkred", h_fill = cols[i])#"#003399" "snow2"
  plot(t2)
}

# find different tones at https://www.w3schools.com/colors/colors_picker.asp
# cerulean "#01A58A"
cshades = c("#01987f", "#017f6a", "#016555", "#004c3f")
cgrade = c("30", "25", "20", "15")

for(i in seq_along(cshades)){
  outfile <- here("man", "figures", paste0("DBCAscatR_logo_", cgrade[i], ".png"))
  t2 <- sticker(p2, package="DBCAscatR", filename = outfile,
                s_width = 1.4, s_height = 0.7,
                s_x = 1, s_y = 0.75,
                p_size = 22, p_color = cshades[i],
                h_size = 1, h_color= cshades[i], h_fill = "#01A58A")#"#003399" "snow2"
  plot(t2)
}





t2 <- sticker(p2, package="DBCAscatR", filename=outfile,
              s_width = 1.4, s_height = 0.7,
              s_x = 1, s_y = 0.75,
              p_size = 22, p_color = "lagoon",
              h_size = 1, h_color= "darkred", h_fill = "rosybrown3")#"#003399" "snow2"

# use this to view as viewer and png very different - plot call will view as
# png
plot(t2)
