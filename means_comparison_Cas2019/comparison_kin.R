
## Some comparisons of the mean distributions by age of kin, for fixed age of Focal.
## Note that I updated Caswell's model under my belief that pi should be defined strictly at
## the time step when Focal's q-th ancestor produced Focal's (q-1)-th, i.e., one step before (q-1) is in its first age-class
## This explains U %*% pi within the Caswell_2019.R file

### matrices of rates
f_mat <- readr::read_rds("data/f_mat.Rds")
u_mat <- readr::read_rds("data/u_mat.Rds")
f_int <- which(mothers_age(u_mat,f_mat)>0)%>%length()

fig_out <- "means_comparison_Cas2019/Figs/"

############################################## Ancestors ###############################################
m_cas <- mother_dist(u_mat, f_mat)
m_cas_50 <- m_cas[,51]
df_hal_mum_50 <- data.frame(age_kin = 0:100,
                            expectation = m_cas_50,
                            method = "Hal")
joe_mum_50 <- list()
for(i in 51:99){
  val <- unconditional_mother_pdf(50, i, u_mat, f_mat, 7)
  joe_mum_50[[(1+length(joe_mum_50))]] <- val
}
joe_mum_50 <- do.call("rbind", joe_mum_50) %>% as.data.frame()
df_joe_mum_50 <- joe_mum_50 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age,
                   expectation = mean,
                   method = "PFM")
## visual
mum50 <- rbind(df_joe_mum_50, df_hal_mum_50) %>%
  ggplot(aes(x = age_kin, y = expectation, color = method, shape = method)) +
  geom_point() + theme_bw() + theme(legend.position = "top")
mum50
ggsave(paste0(fig_out,"mother_F_50.png"), mum50)

## grandmother
gm_cas <- grand_mother_dist(u_mat, f_mat) ## Caswell 2019
gm_cas_0 <- gm_cas[,1]
df_hal_gran_0 <- data.frame(age_kin = 0:100,
                            expectation = gm_cas_0,
                            method = "Hal")
joe_gran_0 <- list() ## PMF kin
for(i in 30:85){
  val <- unconditional_grandmother_pdf(0, i, u_mat, f_mat, 7)
  df <- data.frame(number = seq(0,6),
                   prob = val$prob,
                   age = i)
  joe_gran_0[[(1+length(joe_gran_0))]] <- df
}
joe_gran_0 <- do.call("rbind", joe_gran_0)
df_joe_gran_0 <- joe_gran_0 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age,
                   expectation = mean,
                   method = "PFM")
gran0 <- rbind(df_joe_gran_0, df_hal_gran_0) %>%
  ggplot(aes(x = age_kin, y = expectation, color = method, shape = method)) +
  geom_point() + theme_bw() + theme(legend.position = "top")
gran0
ggsave(paste0(fig_out,"gran_F_0.png"), gran0)

############################################ Sisters ###########################
### Younger sisters (Focal at 50)
ys_cas <- younger_sis_dist(u_mat, f_mat) # Caswell 2019
ys_50_cas <- ys_cas[1:101,51]
df_hal_YS_50 <- data.frame(age_kin = 0:100,
                           expectation = ys_50_cas,
                           method = "Hal")
df_list_ys <- list() # PMF kin
for(i in 0:49){
  df <- ys_PMF(50, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_ys[[(1+length(df_list_ys))]] <- df1
}
df_list_ys_50 <- do.call("rbind", df_list_ys) %>% as.data.frame()
df_joe_YS_50 <- df_list_ys_50 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age,
                   expectation = mean,
                   method = "PFM")
library(ggplot2)
## visual
ys50 <- rbind(df_joe_YS_50, df_hal_YS_50) %>%
  ggplot(aes(x = age_kin, y = expectation, color = method, shape = method)) +
  geom_point() + theme_bw() + theme(legend.position = "top")
ys50
ggsave(paste0(fig_out,"younger_sis_F_50.png"), ys50)

### Older sisters (Focal at 50)
os_cas <- older_sis_dist(u_mat, f_mat) # Caswell 2019
os_50_cas <- os_cas[1:101,51]
df_hal_OS_50 <- data.frame(age_kin = 0:100,
                           expectation = os_50_cas,
                           method = "Hal")
df_list_os <- list() # PMF Kin
for(i in 51:(51+f_int+1)){
  df <- os_PMF(50, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_os[[(1+length(df_list_os))]] <- df1
}
df_list_os_50 <- do.call("rbind", df_list_os) %>% as.data.frame()
df_joe_OS_50 <- df_list_os_50 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age,
                   expectation = mean,
                   method = "PFM")
## visual
os50 <- rbind(df_joe_OS_50, df_hal_OS_50) %>%
  ggplot(aes(x = age_kin, y = expectation, color = method, shape = method)) +
  geom_point() + theme_bw() + theme(legend.position = "top")
os50
ggsave(paste0(fig_out,"older_sis_F_50.png"), os50)

################################################# Aunts ##########################################
### Younger aunts (Focal at 50)
ya_cas <- younger_aunts_dist(u_mat, f_mat) # Caswell 2019
ya_50_cas <- ya_cas[1:101,51]
df_hal_YA_50 <- data.frame(age_kin = 0:100,
                           expectation = ya_50_cas,
                           method = "Hal")
df_list_ya_50 <- list() # PMF Kin
for(i in 45:100){
  df <- ya_PMF(50, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_ya_50[[(1+length(df_list_ya_50))]] <- df1
}
df_list_ya_50 <- do.call("rbind", df_list_ya_50) %>% as.data.frame()
df_joe_YA_50 <- df_list_ya_50 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age,
                   expectation = mean,
                   method = "PFM")
### visual
ya50 <- rbind(df_joe_YA_50, df_hal_YA_50) %>%
  ggplot(aes(x = age_kin, y = expectation, color = method, shape = method)) +
  geom_point() + theme_bw() + theme(legend.position = "top")
ya50
ggsave(paste0(fig_out,"younger_aunt_F_50.png"), ya50)

### Older aunts (Focal at 50)
oa_cas <- older_aunts_dist(u_mat, f_mat)#
oa_50_cas <- oa_cas[1:101,51]
df_hal_OA_50 <- data.frame(age_kin = 0:100,
                           expectation = oa_50_cas,
                           method = "Hal")
df_list_oa_50 <- list()
for(i in 65:100){
  df <- oa_PMF(50, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_oa_50[[(1+length(df_list_oa_50))]] <- df1
}
df_list_oa_50 <- do.call("rbind", df_list_oa_50) %>% as.data.frame()
df_joe_OA_50 <- df_list_oa_50 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age,
                   expectation = mean,
                   method = "PFM")
### visual
oa50 <- rbind(df_joe_OA_50, df_hal_OA_50) %>%
  ggplot(aes(x = age_kin, y = expectation, color = method, shape = method)) +
  geom_point() + theme_bw() + theme(legend.position = "top")
oa50
ggsave(paste0(fig_out,"older_aunt_F_50.png"), oa50)

######################################## Cousins #####################################################
### Younger cousins (Focal at 20)
yc_cas <- younger_cousins_dist(u_mat, f_mat) # Caswell 2019
yc_20_cas <- yc_cas[1:101,21]
df_hal_YC_20 <- data.frame(age_kin = 0:100,
                           expectation = yc_20_cas,
                           method = "Hal")
df_list_cya_20 <- list() # PMF Kin
for(i in 0:60){
  df <- cya_PMF(20, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_cya_20[[(1+length(df_list_cya_20))]] <- df1
}
df_list_cya_20 <- do.call("rbind", df_list_cya_20) %>% as.data.frame()
df_joe_YC_20 <- df_list_cya_20 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age,
                   expectation = mean,
                   method = "PFM")
### visual
yc_20 <- rbind(df_joe_YC_20, df_hal_YC_20) %>%
  ggplot(aes(x = age_kin, y = expectation, color = method, shape = method)) +
  geom_point() + theme_bw() + theme(legend.position = "top")
ggsave(paste0(fig_out,"younger_cousin_F_20.png"), yc_20)
yc_20
### Older cousins (Focal at 20)
oc_cas <- older_cousins_dist(u_mat, f_mat) # Caswell 2019
oc_20_cas <- oc_cas[1:101,21]
df_hal_OC_20 <- data.frame(age_kin = 0:100,
                           expectation = oc_20_cas,
                           method = "Hal")
df_list_coa_20 <- list() # PMF Kin
for(i in 0:80){
  df <- coa_PMF(20, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_coa_20[[(1+length(df_list_coa_20))]] <- df1
}
df_list_coa_20 <- do.call("rbind", df_list_coa_20) %>% as.data.frame()
df_joe_OC_20 <- df_list_coa_20 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age,
                   expectation = mean,
                   method = "PFM")
### visual
oc_20 <- rbind(df_joe_OC_20, df_hal_OC_20) %>%
  ggplot(aes(x = age_kin, y = expectation, color = method, shape = method)) +
  geom_point() + theme_bw() + theme(legend.position = "top")
oc_20
ggsave(paste0(fig_out,"older_cousin_F_20.png"), oc_20)


############## Nieces


### Younger nieces (Focal at 20)
nys_cas <- younger_nieces_dist(u_mat, f_mat) # Caswell 2019
nys_40_cas <- nys_cas[1:101,41]
df_hal_NYS_40 <- data.frame(age_kin = 0:100,
                           expectation = nys_40_cas,
                           method = "Hal")
df_list_nys_40 <- list() # PMF Kin
for(i in 0:20){
  df <- nys_PMF(40, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_nys_40[[(1+length(df_list_nys_40))]] <- df1
}
df_list_nys_40 <- do.call("rbind", df_list_nys_40) %>% as.data.frame()
df_joe_NYS_40 <- df_list_nys_40 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age,
                   expectation = mean,
                   method = "PFM")
### visual
rbind(df_hal_NYS_40, df_joe_NYS_40) %>%
  ggplot(aes(x = age_kin, y = expectation, color = method, shape = method)) +
  geom_point() + theme_bw() + theme(legend.position = "top")


### Older nieces (Focal at 20)
nos_cas <- older_nieces_dist(u_mat, f_mat) # Caswell 2019
nos_40_cas <- nos_cas[1:101,41]
df_hal_NOS_40 <- data.frame(age_kin = 0:100,
                            expectation = nos_40_cas,
                            method = "Hal")
df_list_nos_40 <- list() # PMF Kin
for(i in 0:30){
  df <- nos_PMF(40, i, u_mat, f_mat, 7)
  df1 <- data.frame(number = seq(0,6))
  df1$prob <- df
  df1$age <- i
  df_list_nos_40[[(1+length(df_list_nos_40))]] <- df1
}
df_list_nos_40 <- do.call("rbind", df_list_nos_40) %>% as.data.frame()
df_joe_NOS_40 <- df_list_nos_40 %>%
  dplyr::mutate(X = number*prob,
                X2 = number^2*prob,
                X3 = number^3*prob) %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   s3 = sum(X3),
                   var = s2 - mean^2,
                   m3 = s3 + 2*mean^3 - 3*mean*s2,
                   skew = m3/((var^(0.5))^(3))) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(age_kin = age,
                   expectation = mean,
                   method = "PFM")
### visual
rbind(df_hal_NOS_40, df_joe_NOS_40) %>%
  ggplot(aes(x = age_kin, y = expectation, color = method, shape = method)) +
  geom_point() + theme_bw() + theme(legend.position = "top")


