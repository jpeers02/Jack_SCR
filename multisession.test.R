
source("~/SCR_func.R")



offsets <- expand.grid(dx = c(-700, 0, 700),
                       dy = c(-700, 0, 700))

traps_list <- lapply(1:nrow(offsets), function(i) {
  traps <- as.data.frame(traps.1)  
  traps$x <- traps[,1] + offsets$dx[i]
  traps$y <- traps[,2] + offsets$dy[i]
  traps <- traps[, c("x", "y")]     
  return(traps)
})

colors <- rainbow(9)
plot(NA, xlim = c(-900, 1200), ylim = c(-900, 1200),
     xlab = "X coordinate", ylab = "Y coordinate",
     main = "9 Trap Blocks Spaced by 700 metres")

for (i in seq_along(traps_list)) {
  traps <- traps_list[[i]]
  points(traps$x, traps$y, col = colors[i], pch = 19, cex = 1.5)
}

legend("topright", legend = paste("Block", 1:9),
       col = colors, pch = 19)

