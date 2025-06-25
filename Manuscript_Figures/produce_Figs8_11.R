library(ggplot2)
library(dplyr)

fig_out <- here::here("Manuscipt_Figures_code/Figures")

f_mat <- readr::read_rds(here::here("data","f_mat.Rds"))
u_mat <- readr::read_rds(here::here("data","u_mat.Rds"))

mothers_poss_age <- mothers_age(u_mat, f_mat)
actual_ages <- which(mothers_poss_age!=0)

### Younger Sisters #################################################
ysis_list_foc <- list()
for(i in actual_ages){
  for(foc in 1:20){
  for(j in 0:(foc-1)){
    pdf <- conditional_mothers_age_YS(foc, j, i, u_mat, f_mat, 7)
    temp_df <- data.frame(mo_age = i,
                          prob_YS = pdf,
                          number = seq(0,6),
                          age_YS = j,
                          age_Focal = foc)
    ysis_list_foc[[(1+length(ysis_list_foc))]] <- temp_df
  }
  }
}
ysis_list_foc <- do.call("rbind", ysis_list_foc)
YS_df_foc <- data.frame()
for(i in actual_ages){
  temp <- ysis_list_foc %>% dplyr::filter(mo_age == i)
  for(foc in temp$age_Focal%>%unique()){
    conv_list <- list()
    for(j in 0:(foc-1)){
      temp1 <- temp %>% dplyr::filter(age_YS == j, age_Focal == foc)
      val_con <- temp1 %>% dplyr::select(prob_YS) %>% unname()
      val_con <- val_con[[1]]
      conv_list[[(1+length(conv_list))]] <- val_con}
    val <- convoluion_nth((foc), conv_list)
    dat <- data.frame(number_YS = seq(0,6),
                    prob_YS = val,
                    age_mom_con = i,
                    age_Focal = foc)
    YS_df_foc <- rbind(YS_df_foc, dat)}
}
## chck sum to 1
YS_df_foc %>% dplyr::filter(age_Focal == 10) %>%
  dplyr::group_by(age_mom_con) %>%
  dplyr::summarise(check = sum(prob_YS)) %>%
  dplyr::ungroup()
YS_df_foc$number <- as.factor(YS_df_foc$number_YS)
YS_df_foc$age_Focal_2 <- paste0("Focal ", YS_df_foc$age_Focal)
YS_df_foc$age_Focal_3 = factor(YS_df_foc$age_Focal_2, levels=YS_df_foc$age_Focal_2%>%unique())

legend_title <- "# younger sisters"
ys_age_foc_mom <- YS_df_foc %>%
  ggplot(aes(x = age_mom_con, y = prob_YS, color = number)) +
  geom_line(aes(group = factor(number))) + facet_wrap(~age_Focal_3) +
  scale_color_viridis_d() + theme_bw() +
  theme(legend.position = "top") + xlab("Age of mother at Focal") +
  ylab("Prob focal has # younger sisters") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_color_manual(legend_title,
                     values = c("black",
                                "blue",
                                "purple",
                                "green",
                                "orange",
                                "red",
                                "pink"))
ys_age_foc_mom
#### Older sisters #########################################################################################
osis_list_foc <- list()
for(i in actual_ages){
  for(foc in 0:20){
    for(j in (foc+1):(foc+i-1)){
      pdf <- conditional_mothers_age_OS(foc, j, i, u_mat, f_mat, 7)
      temp_df <- data.frame(mo_age = i,
                            prob_OS = pdf,
                            number = seq(0,6),
                            age_OS = j,
                            age_Focal = foc)
      osis_list_foc[[(1+length(osis_list_foc))]] <- temp_df
    }
  }
}
osis_list_foc <- do.call("rbind", osis_list_foc)
OS_df_foc <- data.frame()
for(i in actual_ages){
  temp <- osis_list_foc %>% dplyr::filter(mo_age == i)
  for(foc in temp$age_Focal%>%unique()){
    conv_list <- list()
    for(j in (foc+1):(foc+i-1)){
      temp1 <- temp %>% dplyr::filter(age_OS == j, age_Focal == foc)
      val_con <- temp1 %>% dplyr::select(prob_OS) %>% unname()
      val_con <- val_con[[1]]
      conv_list[[(1+length(conv_list))]] <- val_con}
    val <- convoluion_nth(length(conv_list), conv_list)
    dat <- data.frame(number_OS = seq(0,6),
                      prob_OS = val,
                      age_mom_con = i,
                      age_Focal = foc)
    OS_df_foc <- rbind(OS_df_foc, dat)}
}
## check sum to 1
OS_df_foc %>% dplyr::filter(age_Focal == 10) %>%
  dplyr::group_by(age_mom_con) %>%
  dplyr::summarise(check = sum(prob_OS)) %>%
  dplyr::ungroup()
OS_df_foc$number <- as.factor(OS_df_foc$number_OS)
OS_df_foc$age_Focal_2 <- paste0("Focal ",OS_df_foc$age_Focal)
OS_df_foc$age_Focal_3 = factor(OS_df_foc$age_Focal_2,levels=OS_df_foc$age_Focal_2%>%unique())
legend_title <- "# older sisters"
os_age_foc_mom <- OS_df_foc %>%
  ggplot(aes(x = age_mom_con, y = prob_OS, color = number)) +
  geom_line(aes(group = factor(number))) + facet_wrap(~age_Focal_3) +
  scale_color_viridis_d() + theme_bw() +
  theme(legend.position = "top") + xlab("Age of mother at Focal") +
  ylab("Prob focal has # younger sisters") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_color_manual(legend_title,
                     values = c("black",
                                "blue",
                                "purple",
                                "green",
                                "orange",
                                "red",
                                "pink"))
os_age_foc_mom
#################### Combined Sisters ################################################################################

df_age_ys <- YS_df_foc %>%
  dplyr::transmute(number = as.factor(number_YS),
            prob = prob_YS,
            age_mom = age_mom_con,
            age_foc = age_Focal,
            age_Focal = age_Focal_3)
df_age_OS <- OS_df_foc %>%
  dplyr::transmute(number = as.factor(number_OS),
            prob = prob_OS,
            age_mom = age_mom_con,
            age_foc = age_Focal,
            age_Focal = age_Focal_3)
df_age_OS$kin <- "Older Sisters"
df_age_ys$kin <- "Younger Sisters"

#### combined sis
s_squared <- rbind(df_age_ys,df_age_OS)
conv_sis_df <- s_squared %>%
  dplyr::select(number,prob,age_mom,age_foc,kin) %>%
  tidyr::pivot_wider(names_from = kin, values_from = prob) %>%
  dplyr::transmute(number = number,
                   age_mom = age_mom,
                   age_foc = age_foc,
                   p_YS = `Younger Sisters`,
                   p_OS = `Older Sisters`,
                   age_Foc = paste0("Focal ", age_foc)) %>%
  dplyr::group_by(age_mom, age_foc, age_Foc) %>%
  dplyr::summarise(prob = convolve(p_YS, p_OS, conj = FALSE),
                   number = number) %>%
  dplyr::ungroup()
conv_sis_df$age_Foc2 <- factor(conv_sis_df$age_Foc, levels = conv_sis_df$age_Foc%>%unique())
c11 <- s_squared %>% dplyr::transmute(number = number,
                                 prob = prob,
                                 age_mom = age_mom,
                                 age_foc_num = age_foc,
                                 age_Foc = age_Focal,
                                 kin = kin)
c22 <- conv_sis_df %>% dplyr::transmute(number = number,
                                 prob = prob,
                                 age_mom = age_mom,
                                 age_foc_num = age_foc,
                                 age_Foc = age_Foc,
                                 kin = "Combined Sisters")

full_sis_df_prob <- rbind(c11,c22)
full_sis_df_prob$kin <- factor(full_sis_df_prob$kin, levels = c("Younger Sisters", "Older Sisters", "Combined Sisters"))
legend_title <- "Number of sisters"
## Figure 8 in manuscript
prob_plot <- full_sis_df_prob %>%
  dplyr::filter(age_foc_num %in% c(1, 2, 5, 10, 15, 20) ) %>%
  ggplot(aes(x = age_mom, y = prob, color = number)) +
  geom_line(aes(group = factor(number))) + facet_grid(age_Foc~kin) +
  theme_bw() +
  theme(legend.position = "top") + labs(x=expression('Age of mother at Focal, b'[1])) +
  theme(text = element_text(size = 14)) +
  ylab("Prob Focal has # sisters") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_color_manual(legend_title,
                     values = c("black",
                                "blue",
                                "purple",
                                "green",
                                "orange",
                                "red",
                                "pink"))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14)) + theme(strip.text = element_text(size = 14))
prob_plot
#ggsave(paste0(fig_out,"/conditioning_prob_sisters.png"), prob_plot, width = 9, height = 9)

################## The prob of no sisters ################
heat_df_YS <- YS_df_foc %>% dplyr::select(number_YS,prob_YS,age_mom_con,age_Focal) %>%
  tidyr::pivot_wider(names_from = number_YS, values_from = prob_YS) %>%
  dplyr::transmute(age_mom = age_mom_con,
                   age_Focal = age_Focal,
                   prob_no_sis = `0`,
                   prob_a_sis = `1`+`2`+`3`+`4`+`5`+`6`)
heat_df_OS <- OS_df_foc %>%
  dplyr::select(number_OS,prob_OS,age_mom_con,age_Focal) %>%
  tidyr::pivot_wider(names_from = number_OS, values_from = prob_OS) %>%
  dplyr::transmute(age_mom = age_mom_con,
                   age_Focal = age_Focal,
                   prob_no_sis = `0`,
                   prob_a_sis = `1`+`2`+`3`+`4`+`5`+`6`)
heat_df_YS$kin <- "Younger Sisters"
heat_df_OS$kin <- "Older Sisters"
heat_df_full <- s_squared %>% dplyr::select(number,prob,age_mom,age_foc,kin) %>%
  tidyr::pivot_wider(names_from = kin, values_from = prob) %>%
  dplyr::transmute(number = number,
                   age_mom = age_mom,
                   age_foc = age_foc,
                   p_YS = `Younger Sisters`,
                   p_OS = `Older Sisters`,
                   age_Foc = paste0("Focal age ", age_foc)) %>%
  dplyr::group_by(age_mom, age_foc, age_Foc) %>%
  dplyr::summarise(prob = convolve(p_YS, p_OS, conj = FALSE),
                   number = number) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = number, values_from = prob) %>%
  dplyr::transmute(age_mom = age_mom,
                   age_Focal = age_foc,
                   prob_no_sis = `0`,
                   prob_a_sis = `1`+`2`+`3`+`4`+`5`+`6`)
heat_df_full$kin <- "Combined Sisters"

legend_title_2 <- "Probability no sisters"
all_sis_heat <- rbind(heat_df_YS,heat_df_OS,heat_df_full)
all_sis_heat$kin2 <- factor(all_sis_heat$kin, levels = c("Younger Sisters","Older Sisters","Combined Sisters"))
## Figure 9 in manuscript
all_sis_heat_plot <- all_sis_heat %>%
  dplyr::filter(age_Focal != 0) %>%
  ggplot(aes(x = age_Focal, y = age_mom, fill = prob_no_sis)) +
  geom_tile() + facet_wrap(~kin2) +
  scale_fill_viridis_c(legend_title_2,
                       begin = 0, end = 1,
                       limits=c(0,1)) + theme_bw() +
  theme(legend.position = "top",
        legend.text = element_text(vjust = 1, angle = 45, size = 10),
        legend.title = element_text(vjust = 0.85, size = 14)) +
  labs(y = expression('Age of mother at Focal, b'[1])) + xlab("Age of Focal")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size=14)) + theme(strip.text = element_text(size = 14))
all_sis_heat_plot
#ggsave(paste0(fig_out,"/Heat_plot_age_mom.png"), all_sis_heat_plot, width = 9, height = 5)

##### sibsize
conv_sis_df$mum_age_f <- paste0("b[1]", "==" , conv_sis_df$age_mom)
to_add_f20 <- conv_sis_df %>% dplyr::filter(age_foc == 20, age_mom %in% seq(14,52,2)) %>%
  dplyr::mutate(X = prob*as.numeric(number),
         X2 = prob*as.numeric(number)^2)
to_add_f20n <- to_add_f20 %>%
  dplyr::group_by(age_mom) %>%
  dplyr::summarise(exp = sum(X),
            s2 = sum(X2),
            vari = s2 - exp^2,
            prob = prob,
            number = number,
            mum_age_f = mum_age_f) %>%
  dplyr::ungroup()
to_add_f20nn <- to_add_f20n %>% dplyr::mutate(prob = 0.05)
conv_sis_df %>%
  dplyr::filter(age_foc == 20, age_mom %in% seq(14,52,2)) %>%
  dplyr::group_by(age_mom) %>%
  summarise(check = sum(prob)) %>%
  dplyr::ungroup()

## Figure 10 in manuscript
mum_plot_sibsize <- conv_sis_df %>%
  dplyr::filter(age_foc == 20, age_mom %in% seq(14,52,2)) %>%
  ggplot(aes(x = number, y = prob)) +
  geom_bar(position = "dodge" , stat = "identity", alpha = 0.25) +
  facet_wrap(~mum_age_f , labeller = label_parsed) + theme_bw() +
  xlab("Sibsize (S)") + ylab("Probability")  +
  geom_errorbarh(data = to_add_f20nn, aes(xmin = ifelse(exp - sqrt(vari) < 0, 0, exp - sqrt(vari) ),
                                        xmax = exp + sqrt(vari) ),
                 alpha = 0.5, size = 0.5, height = 0.05, color = "red") +
  geom_vline(data = to_add_f20n, aes(xintercept = exp),
             color = "blue", size = 0.5, alpha = 0.5)+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size=14)) + theme(strip.text = element_text(size = 14))


mum_plot_sibsize
#ggsave(paste0(fig_out,"/Sibsize_age_mom.png"), mum_plot_sibsize, width = 9, height = 7)

########## Birth ranks prob
pij <- function(i, S, age_mother){
  test_df <- s_squared %>% dplyr::filter(age_foc == 20, age_mom == age_mother)
  j <- S - i
  pi <- test_df %>% dplyr::filter(kin == "Younger Sisters" , number == i) %>% dplyr::select(prob) %>% unname()
  pj <- test_df %>% dplyr::filter(kin == "Older Sisters" , number == j) %>% dplyr::select(prob) %>% unname()
  p <- pi[[1]]*pj[[1]]
  return(p)
}
birth_rank_df <- data.frame()
for(mother in actual_ages){
  for(sibsize in 0:6){
    for(ys in 0:sibsize){
      temp <- data.frame(age_mother = mother,
                         S = sibsize,
                         ys = ys,
                         os = sibsize - ys,
                         prob = pij(ys, sibsize, mother))
      birth_rank_df <- rbind(birth_rank_df, temp)
    }
  }
}
birth_rank_df$birth_rank <- ifelse(birth_rank_df$S == 0, "only child",
                                   paste0( birth_rank_df$os+1  ," in ", 1+birth_rank_df$S))
birth_rank_df$S_fac <- paste0("S==", birth_rank_df$S)
birth_rank_df$age_m_fac <- paste0("b[1] " ,"==", birth_rank_df$age_mother)
birth_rank_df$fecundity_mom <- mothers_poss_age[birth_rank_df$age_mother]
## check plot
birth_rank_df %>%
  dplyr::group_by(S,ys,os,birth_rank,S_fac) %>%
  dplyr::summarise(prob = sum(fecundity_mom*prob)) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = birth_rank, y = prob)) +
  geom_bar(position = "dodge" , stat = "identity") +
  facet_grid(~S_fac, scales = "free") +
  theme(axis.text.x = element_text(angle = 90))

weight_df <- birth_rank_df %>%
  dplyr::group_by(S,ys,os,birth_rank,S_fac) %>%
  dplyr::summarise(prob = sum(fecundity_mom*prob)) %>%
  dplyr::ungroup()

weight_df <- weight_df %>%
  dplyr::transmute(age_mother = "Weighted",
            S = S,
            ys = ys,
            os = os,
            prob = prob,
            birth_rank = birth_rank,
            S_fac = paste0("S==", S),
            age_m_fac = "Weighted",
            fecundity_mom = NA)
## Figure 11 in manuscript
birth_plot <- rbind(birth_rank_df, weight_df) %>%
  dplyr::filter(age_mother %in% c(15,20,25,30,35,40,45,"Weighted")) %>%
  ggplot(aes(x = birth_rank, y = prob)) + theme_bw() +
  facet_grid(age_m_fac~S_fac, scales = "free", labeller = label_parsed) +
  geom_bar(position = "stack", stat = "identity", alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90)) + xlab("Birth rank")+
  ylab("Probability")+
  theme(axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size=14),
        legend.title = element_text(size=14)) + theme(strip.text = element_text(size = 14))
birth_plot
#ggsave(paste0(fig_out,"/BirthRank_age_mom.png"), birth_plot, width = 9, height = 9)

