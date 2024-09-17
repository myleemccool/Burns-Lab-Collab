#### Smiley face Graph
# Set up the plot window
plot(0, 0, type = "n", xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xlab = "", ylab = "", asp = 1)

# Draw the face circle
symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, bg = "yellow")

# Draw the eyes
points(-0.5, 0.5, pch = 16, cex = 2) # Left eye
points(0.5, 0.5, pch = 16, cex = 2)  # Right eye

# Draw the smile (using an arc)
curve(0.5 * (x^2), from = -0.5, to = 0.5, add = TRUE, col = "black", lwd = 2)

# Optionally, add a title or label
title(main = "Smiley Face!")