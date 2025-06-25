library(dplyr)
`%>%` <- magrittr::`%>%`
## fertility matrix for 1974
fmat <- here::here("data")
ff <- utils::read.csv(paste0(fmat, "/fx_1974_on.txt"),skip = 2, sep = "")
ff <- as.data.frame(ff)
ff <- dplyr::filter(ff, Year == 1974)
ff <- ff[,c(2,3)]
ff
ff$Age <- as.numeric(gsub("[[:punct:]]", "", ff$Age))

lower_fert_age_1974 <- ff$Age[1]-1
upper_fert_age_1974 <- ff$Age[length(ff$Age)]+1

padding_low_1974 <- data.frame(Age = seq(0,lower_fert_age_1974),
                               ASFR = rep(0,12))
padding_high_1974 <- data.frame(Age = seq(upper_fert_age_1974,100),
                                ASFR = 0)

ff <- rbind(padding_low_1974,ff,padding_high_1974)
ff
f_mat <- matrix(0, nrow(ff), nrow(ff))
f_mat[1,] <- 0.5*ff$ASFR


#saveRDS(f_mat, "data/f_mat.Rds")
f_mat
## mortality matrix for 1974
mmat <- here::here("data")
mm <- utils::read.csv(paste0(mmat, "/qx_1974_on.txt"),skip = 2, sep = "")
mm <- as.data.frame(mm)
mm <- dplyr::filter(mm, Year == 1974)
mm <- dplyr::select(mm, c(Age,qx))
mm$Age <- as.numeric(gsub("[[:punct:]]", "", mm$Age))
mm <- dplyr::filter(mm, Age < 101)
mm$qx <- 1 - mm$qx
u_mat <- matrix(0, nrow(mm), nrow(mm))
diag(u_mat[-1,-ncol(u_mat)]) <- mm$qx[1:(length(mm$qx)-1)]
u_mat[101,101] <- mm$qx[length(mm$qx)]
#saveRDS(u_mat, "data/u_mat.Rds")
u_mat
### check matrix
f_mat
u_mat

ff

### convert to csv
df_fert_csv <- data.frame(Age = seq(0,100),
                          Fert = rep(0,101))
for(i in seq(12,55)){
  ffnew <-  ff %>% dplyr::filter(Age == i) %>% dplyr::select(ASFR) %>% unname()
  ffnew <- ffnew[[1]]
  df_fert_csv[(i+1),2] <- ffnew*0.5
}
df_fert_csv

f_mat[1,31]
pd_df <- df_fert_csv %>% dplyr::filter(Age == 14) %>% dplyr::select(Fert) %>% unname()
pd_df[[1]]

#write.csv(df_fert_csv, "data/1974_female_fert.csv")


mm$qx
df_mort_csv <- data.frame(Age = seq(0,100),
                          Mort = rep(0,101))
for(i in seq(0,100)){
  df_mort_csv[(i+1),2] <- mm$qx[(i+1)]
}
df_mort_csv

df_mort_csv[,2]
c(diag(u_mat[-1,-ncol(u_mat)]),u_mat[101,101])

#write.csv(df_mort_csv, "data/1974_female_mort.csv")

SD <- function(PM) {
  spectral_stuff <- eigen(PM)
  spectral_stuff <- Re(spectral_stuff$vectors[, which.max(abs(spectral_stuff$values))])
  # normalise...
  vec_lambda <- spectral_stuff/sum(spectral_stuff)
  return(vec_lambda)
}

pop_struct <- SD(u_mat+f_mat)
length(pop_struct)
df_ps <- data.frame(Age = seq(0,100),
                    ps = pop_struct)
pop_struct%>%sum()
#write.csv(pop_struct, "data/1974_female_ps.csv")
