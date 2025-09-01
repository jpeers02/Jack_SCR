
source("~/SCR_func.R")

offsets <- expand.grid(dx = c(-600, 0, 600),
                       dy = c(-600, 0, 600))

traps_list <- lapply(1:nrow(offsets), function(i) {
  traps <- traps.1
  traps[,1] <- traps[,1] + offsets$dx[i]
  traps[,2] <- traps[,2] + offsets$dy[i]
  data.frame(traps)
})

names(traps_list) <- paste0("Block", 1:9)

colors <- rainbow(9)
plot(NA, xlim = c(-600, 1000), ylim = c(-600, 1000),
     xlab = "X coordinate", ylab = "Y coordinate",
     main = "9 Trap Blocks Spaced by 600")

for (i in seq_along(traps_list)) {
  traps <- traps_list[[i]]
  points(traps$x, traps$y, col = colors[i], pch = 19, cex = 1.5)
}

legend("topright", legend = paste("Block", 1:9),
       col = colors, pch = 19)

traps_list

trial <- sim_scr_data_multi(10, 0.9, 75, c(-600, 1000), c(-600, 1000), traps_list =  traps_list, n_occasions = 10)


