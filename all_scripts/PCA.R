library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(stringr)

pca_plot <- function(data,
                     factors,
                     color = "",
                     shape = "",
                     th = theme(),
                     text = F,
                     save = "") {
    pca <- t(log10(data + 1))
    pca <- prcomp(pca)
    df_out <- as.data.frame(pca$x) %>%
        rownames_to_column(var = "Sample") %>%
        right_join(factors) %>%
        column_to_rownames("Sample") %>%
        droplevels.data.frame()
    percentaget <- round(pca$sdev / sum(pca$sdev) * 100, 2)
    percentaget <- paste(colnames(df_out[, 1:8]), "(", paste(as.character(percentaget), "%", ")", sep = "")) # nolint

    if (color != "") {
        C <- df_out[, color]
    } else {
        C <- "black"
    }
    if (shape != "") {
        S <- df_out[, shape]
    } else {
        S <- 1
    }

    t <- ggplot(df_out, aes(x = PC1, y = PC2)) +
        geom_point(
            size = 4,
            show.legend = T
        ) +
        aes(color = as.factor(C), shape = as.factor(S)) +
        th +
        xlab(percentaget[1]) +
        ylab(percentaget[2]) +
        ggtitle("PCA on Log2(x+1) transformed TPM")

    if (color != "") {
        t <- t + guides(color = guide_legend(color))
    } else {
        t <- t + guides(color = "none")
    }
    if (shape != "") {
        t <- t + guides(shape = guide_legend(shape))
    } else {
        t <- t + guides(shape = "none")
    }
    if (text) {
        t <- t + geom_text(
            size = 2.5,
            nudge_y = 2,
            aes(label = rownames(df_out))
        )
    }

    print(t)

    if (save != "") {
        ggsave(save, t, device = "pdf", width = 10, height = 8)
    }
}
