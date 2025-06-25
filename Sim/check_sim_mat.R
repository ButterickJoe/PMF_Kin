library(readxl)
library(readr)
library(dplyr)
library(MASS)
library(tidyverse)


###########################


sis_df <- data.frame()
for(j in 1:575){

  data1 <- lapply(2:100, function(i){
    read_csv(paste0("C:/Users/jb4u23/OneDrive - University of Southampton/Butterick year 1/Python code/Kinship modelling/Simple Branching/Emulating_PDF_matrixmodel/Results_1/Time",i ,"replication", j ,".csv"),
             col_types = cols())})

  data11 <- do.call("rbind", data1) %>% as.data.frame() %>% dplyr::filter(Kin == "Younger siblings")
  data11$sim_no <- j
  sis_df <- rbind(sis_df, data11)


}



sis_df %>% head()
sis_df <- sis_df %>% dplyr::select(Year,`Kin ID`, Alive, `Age Kin`, `Age Focal`, `Focal ID`, sim_no)
colnames(sis_df) <- c("Year", "Kin_ID", "Alive", "Age_Kin", "Age_Focal", "Focal_ID", "sim_no")
sis_df %>% head()
sis_df1 <- sis_df %>%
  dplyr::group_by(Age_Focal, sim_no) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number, sim_no, Age_Focal) %>%
  dplyr::summarise(rel_number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Age_Focal) %>%
  dplyr::summarise(normaliser = sum(rel_number),
                   prob = rel_number/normaliser,
                   number = number) %>%
  dplyr::ungroup() %>%
  group_by(Age_Focal) %>%
  dplyr::summarise(X = number*prob,
                   X2 = number^2*prob,
                   expectation = sum(X),
                   variance = sum(X2) - expectation^2,
                   prob = prob,
                   cum_p = cumsum(prob),
                   number = number) %>%
  ungroup()

sis_df1 %>%
  ggplot(aes(x = Age_Focal)) + geom_line(aes(y = expectation)) + geom_point(aes(y = variance^(0.5)))


sis_df %>% filter(Year == 99) %>% dplyr::select(Age_Focal) %>% unique()




full_mother<- readRDS("C:/Temp/Branching/data/MUM_sim_py2.Rds")
full_mother <- full_mother %>% dplyr::filter(Age_Focal == 0)
fm <- full_mother %>% group_by(Focal_ID, sim_no, Age_Kin) %>%
  dplyr::summarise(number = sum(Alive)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number, Focal_ID, sim_no, Age_Kin) %>%
  dplyr::summarise(rel_number = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number, Age_Kin) %>%
  dplyr::summarise(occurance = sum(rel_number)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(number) %>%
  dplyr::summarise(total_occurance = sum(occurance),
                   prob_occurance = occurance/total_occurance,
                   Age_Kin = Age_Kin) %>%
  dplyr::ungroup()




fm %>% ggplot(aes(x = Age_Kin, y = prob_occurance)) + geom_point()

f_mat <- readr::read_rds("data/f_mat.Rds")
u_mat <- readr::read_rds("data/u_mat.Rds")

########### We plot the accumulated kin of Focal for each age of Focal

joe_mum_0 <- list()
for(i in 12:55){
  val <- unconditional_mother_pdf(0, i, u_mat, f_mat, 7)
  joe_mum_0[[(1+length(joe_mum_0))]] <- val
}
joe_mum_0 <- do.call("rbind", joe_mum_0) %>% as.data.frame()
joe_mum_0 %>% head()
joe_mum_0$X <- joe_mum_0$number*joe_mum_0$prob
joe_mum_0$X2 <- joe_mum_0$number^2*joe_mum_0$prob
joe_mum_0_means <- joe_mum_0 %>%
  dplyr::group_by(age) %>%
  dplyr::summarise(mean = sum(X),
                   s2 = sum(X2),
                   var = s2 - mean^2) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(method = "Age dist.",
                age_Focal = "At 20 years",
                kin = "Mothers")


fm %>% head()
m_IC <- joe_mum_0_means %>% dplyr::select(age, mean) %>%
  dplyr::transmute(Age = age, rhob1 = mean, method = "PMF")
s_IC <- fm %>% dplyr::select(Age_Kin, prob_occurance) %>%
  dplyr::transmute(Age = Age_Kin, rhob1 = prob_occurance, method = "MicroSsimulation")

fm$Age_Kin %>% unique()

mothers_plot <- rbind(m_IC, s_IC) %>% ggplot(aes(x = Age, y = rhob1, color = method)) + geom_point() +
  theme(legend.position = "top")

ggsave("age_mom_sim.png", mothers_plot)

sum(s_IC$rhob1)
saveRDS(s_IC, "data/alt_rhob1")


fff <- fert_dists(f_mat, 7)
for(i in 1:length(fff)){
  print(sum(fff[[i]]))
}
