## Code used to prepare the datasets of the package (not run at build time)

# PolSAR image: object `ImgMP_Sim` of SanFranciscoImage_150_150.Rdata
load("SanFranciscoImage_150_150.Rdata")
SanFrancisco150 <- ImgMP_Sim

# Ground truth: gt_SF_150x150.png (black = ocean, grey = vegetation,
# white = urban)
p <- png::readPNG("gt_SF_150x150.png")[, , 1]
SanFrancisco150_gt <- ifelse(p < 0.2, 1L, ifelse(p < 0.7, 2L, 3L))
storage.mode(SanFrancisco150_gt) <- "integer"

usethis::use_data(SanFrancisco150, SanFrancisco150_gt, overwrite = TRUE,
                  compress = "xz")
