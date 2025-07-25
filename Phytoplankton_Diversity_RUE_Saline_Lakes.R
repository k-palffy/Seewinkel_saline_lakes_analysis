library(quantreg);library(mgcv);library(nlme)
library(vegan);library(FD);library(BiodiversityR)
library(CCA);library(CCP);library(MASS);library(AICcmodavg)
library(tidygam);library(drc);library(stringr);library(ggrepel)
library(tidyverse);library(investr);library(gridExtra);library(grid)
library(GGally);library(ggpubr);library(ggallin);library(ggfortify)
library(partykit);library(strucchange);library(ggparty);library(viridis)
library(scales);library(ggnewscale);library(lattice);library(ggforce)
library(readr)

# Loading data from figshare archive
# Download csv files from https://doi.org/10.6084/m9.figshare.29128589.v1
# or use the direct links below
url_tax16    <- "https://figshare.com/ndownloader/files/56554523" # Seewinkel_16S_OTU.csv
url_tax18    <- "https://figshare.com/ndownloader/files/56554520" # Seewinkel_18S_OTU.csv
url_traits16 <- "https://figshare.com/ndownloader/files/56554514" # Seewinkel_Cyano_OTU_trait_matrix.csv
url_traits18 <- "https://figshare.com/ndownloader/files/56554532" # Seewinkel_Euk_Phyto_OTU_trait_matrix.csv
url_abund16  <- "https://figshare.com/ndownloader/files/56554526" # Seewinkel_16S_OTU_abundance.csv
url_abund18  <- "https://figshare.com/ndownloader/files/56554529" # Seewinkel_18S_OTU_abundance.csv
url_env      <- "https://figshare.com/ndownloader/files/56554517" # Seewinkel_background_data.csv

tax16    <- as.data.frame(read_csv(url_tax16)) # 16S OTU taxonomic table 
tax18    <- as.data.frame(read_csv(url_tax18)) # 18S OTU taxonomic table
traits16 <- as.data.frame(read_csv(url_traits16))     # cyanobacterial OTU trait matrix
traits18 <- as.data.frame(read_csv(url_traits18)) # algal OTU trait matrix
abund16  <- as.data.frame(read_csv(url_abund16)) # 16S OTU abundance data
abund18  <- as.data.frame(read_csv(url_abund18)) # 18S OTU abundance data
env      <- as.data.frame(read_csv(url_env)) # environmental data

# FIGURE 1
# Functional group abundances (violin plot)
# Data processing
# Generating data frames for functional group abundances
# Functional trait and functional group assignments are detailed in Materials and methods.
abund16_t              <- t(abund16[,-c(1:4)])
abund16_fungroup       <- apply(abund16_t, 2, function(x) tapply(x, tax16$fungroup, sum))
abund16_fungroup       <- data.frame(abund16[,1:3],t(abund16_fungroup))
abund16_fungroup       <- rbind(abund16_fungroup, abund16_fungroup)
abund16_fungroup$group <- rep(c("Cyano","N.fixing_Cyano"), each = nrow(abund16_fungroup)/2)
abund16_fungroup$Cyano[(nrow(abund16_fungroup)/2+1):nrow(abund16_fungroup)] <- 0
abund16_fungroup$N.fixing_Cyano[1:nrow(abund16_fungroup)/2] <- 0
abund16_fungroup$abundance <- abund16_fungroup$Cyano + abund16_fungroup$N.fixing_Cyano
abund16_fungroup           <- abund16_fungroup[,c(1:3,6,7)]
abund16_fungroup$Year      <- as.character(abund16_fungroup$Year)

abund18_t           <- t(abund18[,-c(1:4)])
abund18_fungroup    <- apply(abund18_t, 2, function(x) tapply(x, tax18$fungroup, sum))
abund18_fungroup    <- data.frame(abund18[,1:3],t(abund18_fungroup))
group_number        <- ncol(abund18_fungroup)-3
abund18_fungroup_df <- as.data.frame(matrix(NA, nrow = nrow(abund18_fungroup)*group_number, ncol = 5))
colnames(abund18_fungroup_df) <- c(colnames(abund18_fungroup)[1:3],"group","abundance")
abund18_fungroup_df$Sample    <- rep(abund18_fungroup$Sample, times = group_number)
abund18_fungroup_df$Year    <- rep(abund18_fungroup$Year, times = group_number)
abund18_fungroup_df$Site    <- rep(abund18_fungroup$Site, times = group_number)
abund18_fungroup_df$group     <- rep(colnames(abund18_fungroup)[4:(3 + group_number)], each = nrow(abund18_fungroup))
abund18_fungroup_df$abundance <- c(abund18_fungroup$Attached.forms,
                                               abund18_fungroup$Autotr_flagellates,
                                               abund18_fungroup$Mixotrophs,
                                               abund18_fungroup$Non.mot_multicell_green,
                                               abund18_fungroup$Non.mot_unicell_green,
                                               abund18_fungroup$Other_Ochrophytes,
                                               abund18_fungroup$Picoeukaryotes,
                                               abund18_fungroup$Planktonic_Diatoms)
abund18_fungroup_df$Year      <- as.character(abund18_fungroup_df$Year)

# Fig 1.a
# plot for Cyanobacteria
colors_years <- c("2017" = "red", "2018" = "blue")
colors_years2 <- c("2017" = "pink", "2018" = "lightskyblue")
colors_years3 <- c("2017" = "red", "2018" = "dodgerblue2")

ggplot(aes(x = group, y = abundance, fill = forcats::fct_rev(Year)),
       data = abund16_fungroup) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 11)) +
  geom_violin(scale = "width", linewidth = 0.1) +
  geom_jitter(aes(color = forcats::fct_rev(Year)),
              shape = 16, position = position_jitterdodge(jitter.width = 0.4, dodge.width = 1),
              size = 1) +
  scale_fill_manual(values = colors_years2, breaks=c("2017","2018")) +
  guides(fill = guide_legend(title = "Year")) +
  scale_y_continuous(transform = "log1p", limits = c(-0.1,10000), breaks = c(0,10,100,1000,10000)) +
  coord_flip() +
  stat_summary(aes(group = forcats::fct_rev(Year)),
               geom = "point",
               fun = "mean",
               fill = "black", col = "black", size = 1.5, shape = 21,
               position = position_dodge(width = 0.9)
  ) +
  scale_color_manual(values = colors_years, breaks=c("2017","2018")) +
  scale_x_discrete(labels = c("N.fixing_Cyano" = expression(N[2]*"-fixing cyanobacteria"),
                              "Cyano" = "Cyanobacteria")) +
  xlab("Functional group") +
  ylab("OTU abundance") +
  guides(fill = guide_legend(title = "Year"), color = guide_legend(title = "Year"))

# Fig 1.b
# plot for algae
ggplot(aes(x = group, y = abundance, fill = forcats::fct_rev(Year)),
      data = abund18_fungroup_df) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 11)) +
  geom_violin(scale = "width", linewidth = 0.1) +
  geom_jitter(aes(color = forcats::fct_rev(Year)),
              shape = 16, position = position_jitterdodge(jitter.width = 0.4, dodge.width = 1),
              size = 1) +
  scale_fill_manual(values = colors_years2, breaks=c("2017","2018")) +
  guides(fill = guide_legend(title = "Year")) +
  scale_y_continuous(transform = "log1p", limits = c(-0.1,10000), breaks = c(0,10,100,1000,10000)) +
  coord_flip() +
  stat_summary(aes(group = forcats::fct_rev(Year)),
               geom = "point",
               fun = "mean",
               fill = "black", col = "black", size = 1.5, shape = 21,
               position = position_dodge(width = 0.9)
  ) +
  scale_color_manual(values = colors_years, breaks=c("2017","2018")) +
  scale_x_discrete(labels = c("Planktonic_Diatoms" = "Planktonic diatoms",
                              "Other_Ochrophytes" = "Other Ochrophya",
                              "Non.mot_unicell_green" = "Unicellular chlorophytes",
                              "Non.mot_multicell_green" = "Multicellular chlorophytes",
                              "Autotr_flagellates" = "Autotrophic flagellates",
                              "Attached.forms" = "Attached forms")) +
  xlab("Functional group") +
  ylab("OTU abundance") + 
  guides(fill = guide_legend(title = "Year"), color = guide_legend(title = "Year"))

# SUPPLEMENTARY FIGURE S1
# Scatterplot and correlation matrix
# of cyanobacteria, N2-fixing cyanobacteria and environmental variables
env$algal_free_TSS[env$TSS != 0] <- env$TSS[env$TSS != 0] - env$chl[env$TSS != 0]/10 # determining algal-free TSS

abund16_fungroup2      <- apply(abund16_t, 2, function(x) tapply(x, tax16$fungroup, sum))
abund16_fungroup2      <- data.frame(abund16[,1:3],t(abund16_fungroup2))
abund16_fungroup2      <- cbind(abund16_fungroup2,env[,c(4:6,10,11,14)])
abund16_fungroup2$Year <- as.character(abund16_fungroup2$Year)

abund16_fungroup_log         <- log10(abund16_fungroup2[,4:5]+1)
abund16_fungroup_log         <- cbind(abund16_fungroup2[,1:3], abund16_fungroup_log)
abund16_fungroup_log         <- cbind(abund16_fungroup_log,env[,c(4:6,10,11,14)])
abund16_fungroup_log$log_TP  <- log10(abund16_fungroup_log$TP)
abund16_fungroup_log$log_TN  <- log10(abund16_fungroup_log$TN)
abund16_fungroup_log$log_algal_free_TSS <- NA
abund16_fungroup_log$log_algal_free_TSS[env$TSS != 0] <- log10(abund16_fungroup_log$algal_free_TSS[env$TSS != 0])
abund16_fungroup_log         <- abund16_fungroup_log[,c(1:3,6:14,4,5)]
colnames(abund16_fungroup_log)[c(5,10:14)] <- c("Conductivity","log(TP)","log(TN)","log(af-TSS)",
                                                "log(cyano)","log(N-fix cyano)")

ggpairs(abund16_fungroup_log[,c(4:6,10:14)]) +
  theme_bw() +
  theme(panel.grid = element_blank())

# FIGURE 2.a
# Conductivity versus total cyanobacterial abundance scatterplot
ggplot(aes(x = cond, y = Cyano + N.fixing_Cyano, color = Year), data = abund16_fungroup2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 13, vjust = -0.4),
        axis.title.y = element_text(size = 13, vjust = 2),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.88,0.12),
        legend.background = element_rect(colour = "black", linewidth = 0.2)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors_years3) +
  scale_x_log10() +
  scale_y_continuous(transform = "pseudo_log", breaks = c(0, 100, 500, 1000, 3000)) +
  xlab("Conductivity (mS/cm)") +
  ylab("Total cyanobacterial OTU abundance") +
  guides(color = guide_legend(title = "Year"), fill = guide_legend(title = "Year"))

# FIGURE 2.b
# Canonical correspondence analysis (CCA) on eukaryotic algae
# log-transform variables with high variance
env$log_af_TSS <- log10(env$algal_free_TSS)
env$log_TP     <- log10(env$TP)
env$log_TN     <- log10(env$TN)

colnames(abund18_fungroup)[c(4:9,11)] <- c("Attached_Forms","Autotr_Flagellates","Mixotrophs","Multi_Chloro",
                                           "Non-mot_Uni_Chloro","Other_Ochrophytes","Planktonic_Diatoms")

cca_euk_log             <- cca(abund18_fungroup[,4:ncol(abund18_fungroup)]~., env[,c(4:6,15:17)],
                               na.action = na.exclude)
cca_euk_log_step        <- ordistep(cca_euk_log, scope = formula(cca_euk_log))
length(which(vif.cca(cca_euk_log_step) > 10)) > 0 # check superfluous variables
sign_var_log            <- which(anova.cca(cca_euk_log_step, by = "terms")$Pr < 0.05) # Filter non-significant variables
env.color               <- rep("steelblue1", length = ncol(env[,c(4:6,10,11,14)]))
env.color[sign_var_log] <- "blue"
temp.inertia_log        <- cca_euk_log_step$tot.chi
temp.eig_log            <- cca_euk_log_step$CCA$eig
cca.percent_log         <- c(temp.eig_log[1]/temp.inertia_log,temp.eig_log[2]/temp.inertia_log)

# Adjust graphical parameters in case plot window size affects layout 
# (e.g., legend position or title bounding box).
par(oma = c(0.4,0.3,0,0), mar = c(5,6,1,1), cex.axis = 1.35, cex.lab = 2, mgp = c(3.5, 1 ,0))
plot(cca_euk_log_step, choices = c(1, 2), display = c("sp"), scaling = "symmetric",
     type = "n",
     xlab = paste("CCA1 (",round(cca.percent_log[1]*100, digit = 2),"%)", sep = ""),
     ylab = paste("CCA2 (",round(cca.percent_log[2]*100, digit = 2),"%)", sep = ""),
     xlim =c(-1.5,2), ylim = c(-1.5,2.4),
     xaxt = "n", yaxt = "n")
axis(1, at = seq(-1.5, 2, by = 0.5))
axis(2, at = seq(-1.5, 2, by = 0.5))
points(cca_euk_log_step,
       display = "sites",
       choices = c(1, 2),
       scaling = "symmetric",
       col = "red", cex = 1.85,
       pch = 16)
points(cca_euk_log_step,
       display = "sites",
       choices = c(1, 2),
       scaling = "symmetric",
       select = c(23:44),
       col = "dodgerblue2", cex = 1.85,
       pch = 16)
points(cca_euk_log_step,
       display = "bp",
       choices = c(1, 2),
       scaling = "symmetric",
       col = env.color,
       lwd = 2,
       head.arrow = 0.1)
text(cca_euk_log_step,
     display = c("bp"),
     choices = c(1, 2),
     scaling = "symmetric",
     col = env.color, cex = 1.5,
     font = 2)
text(cca_euk_log_step,
     display = c("sp"),
     choices = c(1, 2),
     scaling = "symmetric",
     col = "black", cex = 1.4,
     font = 1)
legend(1.75, 2.4, legend = c("2017","2018"), title = "",
       col = c("red","dodgerblue2"), pt.cex = 1.85,
       pch = 16, text.width = 0.27, y.intersp = c(1.05,1.05), cex = 1.3)
text(x = 1.8, y = 2.3, labels = "Year", adj = 0, cex = 1.5)
rect(-0.75, 2.225, 1.25, 2.5, col = "white", border = "black", lwd = 1)
text(x = 0.25, y = 2.35, labels = "Autotrophic eukaryotes", cex = 1.7, font = 2)

# Determining Resource Use Efficiency (RUE) and Functional Diversity Indices
# in a new data frame
env_diversity     <- env
env_diversity$RUE <- env_diversity$chl / env_diversity$TN # RUE

# creating combined trait matrix
colnames(traits16)[colnames(traits16) == "otus"] <- "OTU"
traits16$mixotroph  <- 0
traits16$flagellate <- 0
traits16$silicate   <- "no"
traits16$pigment    <- "chl_a_phycob"

traits18$N_fixation <- 0

traits_combined           <- rbind(traits16[,c(1,8:15)], traits18[,c(1,13:15,19,11,12,18,17)])
rownames(traits_combined) <- traits_combined$OTU
traits_combined           <- traits_combined[,-1]

abund16_cyano           <- abund16
rownames(abund16_cyano) <- abund16_cyano$Sample
abund16_cyano           <- abund16_cyano[,colnames(abund16_cyano) %in% traits16$OTU]

abund18_alga            <- abund18
rownames(abund18_alga)  <- abund18_alga$Sample
abund18_alga            <- abund18_alga[,colnames(abund18_alga) %in% traits18$OTU]

abundance_combined_pres  <- cbind(abund16_cyano,abund18_alga)
abundance_combined_pres[abundance_combined_pres > 0] <- 1
zero_abs                 <- which(colSums(abundance_combined_pres) == 0)
abundance_combined_pres  <- abundance_combined_pres[,-zero_abs]
traits_combined          <- traits_combined[-zero_abs,]
traits_combined$form     <- as.factor(traits_combined$form) 
traits_combined$silicate <- as.factor(traits_combined$silicate)
traits_combined$pigment  <- as.factor(traits_combined$pigment)

rownames(traits16) <- traits16$OTU
rownames(traits18) <- traits18$OTU

traits16$form      <- as.factor(traits16$form)
traits18$form      <- as.factor(traits18$form)
traits18$pigment   <- as.factor(traits18$pigment)
traits18$silicate  <- as.factor(traits18$silicate)

zero_cyano   <- which(rowSums(abund16_cyano) == 0)
zero_sp      <- which(colSums(abund16_cyano) == 0)
zero_sp_alga <- which(colSums(abund18_alga) == 0)

# Functional diversity indices using 'FD' package
fundiv_alga  <- dbFD(traits18[-zero_sp_alga,c(11:15,17,18)], abund18_alga[,-zero_sp_alga],
                     asym.bin = 1,
                     calc.FRic = TRUE, calc.FGR = FALSE, calc.FDiv = FALSE,
                     calc.CWM = FALSE, corr = "cailliez")

fundiv_combined <- dbFD(traits_combined, abundance_combined_pres,
                        w.abun = FALSE, asym.bin = c(4,5),
                        calc.FRic = TRUE, calc.FGR = FALSE, calc.FDiv = FALSE,
                        calc.CWM = FALSE, corr = "cailliez")

# abundance-weighted indices for algae
env_diversity$FRic_alga <- fundiv_alga$FRic
env_diversity$FEve_alga <- fundiv_alga$FEve
env_diversity$FDis_alga <- fundiv_alga$FDis

# presence-based indices for algae and cyanobacteria
env_diversity$UTC       <- fundiv_combined$sing.sp
env_diversity$FRic_pr   <- fundiv_combined$FRic
env_diversity$FEve_pr   <- fundiv_combined$FEve
env_diversity$FDis_pr   <- fundiv_combined$FDis

# evenness of algal unique trait combinations (UTC)
identical_norow <- function(x,ya) {
  rownames(x)  <- NULL
  rownames(ya) <- NULL
  identical(x,ya)
}

identical_vector <- function(otu_check, df_check) {
  df_check  <- as.matrix(df_check)
  which(apply(df_check, MARGIN = 1, FUN = identical_norow, ya = otu_check))
}

alga_utcs <- unique(traits18[-zero_sp_alga,c(11:15,17,18)])
rownames(alga_utcs) <- NULL
alga_utcs$number    <- paste0("utc",c(1:nrow(alga_utcs)))

traits18_2 <- traits18[-zero_sp_alga,]
identity_vector <- apply(traits18_2[,c(11:15,17,18)],
                         MARGIN = 1,
                         FUN = identical_vector,
                         df_check = alga_utcs[,-ncol(alga_utcs)])

traits18_2$utc_category <- alga_utcs$number[identity_vector]

abund_alga_utcs   <- apply(abund18_alga[,-zero_sp_alga], 1, function(x) tapply(x, traits18_2$utc_category, sum))
abund_alga_utcs   <- t(abund_alga_utcs)

env_diversity$H_alga <- diversity(abund_alga_utcs, index = "shannon", base = exp(1))
env_diversity$J_alga <- env_diversity$H_alga/log(specnumber(abund_alga_utcs))

# FIGURE 3
# Canonical Correlation Analysis (CCOrA)
# on a set of environmental variables (Variable Group No. 1)
# and on diversity indices and RUE as a set of community-related variables  (Variable Group No. 2)
# Data preprocessing:
# Step 1: box-cox transformation for non-normally distributed variables
# checking distribution of all variables
for (i in colnames(env_diversity)[-(1:3)]) {
  print(densityplot(~env_diversity[[i]], main = i))
}
# density plots show non-normality in the values of
# cond, TN, algal_free_TSS, RUE, J_alga, FRic_alga, FRic_pr, FEve_pr
for(i in c("cond","TN","algal_free_TSS",
           "RUE","J_alga","FRic_alga","FRic_pr","FEve_pr")) {
  
  bc <- boxcox(lm(env_diversity[[i]] ~ 1), plotit = FALSE)
  lambda <- bc$x[which.max(bc$y)]
  if (abs(lambda) < 1e-8) {
    env_diversity[[paste0(i, "_bc")]] <- log(env_diversity[[i]])
  } else {
    env_diversity[[paste0(i, "_bc")]] <- (env_diversity[[i]]^lambda - 1) / lambda
  }
}

# Step 2: standardization
env_diversity_scaled <- scale(env_diversity[,-(1:3)])
env_diversity_scaled <- env_diversity_scaled[,c("Z","pH","FEve_alga","FDis_alga","UTC","FDis_pr",
                                                "cond_bc","TN_bc","algal_free_TSS_bc",
                                                "RUE_bc","J_alga_bc","FRic_alga_bc","FRic_pr_bc","FEve_pr_bc")]
env_diversity_scaled <- cbind(env_diversity[,1:3], env_diversity_scaled)
colnames(env_diversity_scaled)[10:ncol(env_diversity_scaled)] <-
  str_sub(colnames(env_diversity_scaled)[10:ncol(env_diversity_scaled)], end = -4,)

# CCorA
cc_total <- cc(env_diversity_scaled[,c("Z","pH","cond","TN","algal_free_TSS")],
               env_diversity_scaled[,c("FEve_alga","FDis_alga","UTC","FDis_pr","RUE",
                                       "J_alga","FRic_alga","FRic_pr","FEve_pr")])
# significance of canonical correlations
rho <- cc_total$cor # canonical correlations
n <- nrow(env_diversity_scaled)
p <- 5
q <- 9
p.asym(rho, n, p , q, tstat = "Wilks") # SUPPLEMENTARY TABLE S2

# Canonical loading plot (Fig. 3)
X <- env_diversity_scaled[,c("Z","pH","cond","TN","algal_free_TSS")]
Y <- env_diversity_scaled[,c("FEve_alga","FDis_alga","UTC","FDis_pr","RUE",
                             "J_alga","FRic_alga","FRic_pr","FEve_pr")]
U <- as.matrix(X) %*% cc_total$xcoef
V <- as.matrix(Y) %*% cc_total$ycoef
loadings_X <- cor(X, U, use = "pairwise.complete.obs")  # Canonical loadings for X variables
loadings_Y <- cor(Y, V)  # Canonical loadings for Y variables
loadings_df <- rbind(
  data.frame(var = rownames(loadings_X), x = loadings_X[,1], y = loadings_X[,2], group = "Environmental variables"),
  data.frame(var = rownames(loadings_Y), x = loadings_Y[,1], y = loadings_Y[,2], group = "Community-related variables")
)

# rescaled sample scores
U_scaled <- sweep(U[,1:2], 2, apply(U[,1:2], 2, function(x) max(abs(x), na.rm = TRUE)), FUN = "/")
U_scaled <- as.data.frame(U_scaled)
colnames(U_scaled) <- c("x", "y")
env_diversity$Year <- as.factor(env_diversity$Year)
U_scaled$Year      <- env_diversity$Year

loadings_df$var_label                   <- loadings_df$var
loadings_df$var_label[c(3,5:7,9,11:14)] <- c("Cond","afTSS","FEve[alga]","FDis[alga]",
                                             "FDis[pr]","J[alga]","FRic[alga]",
                                             "FRic[pr]","FEve[pr]")

ggplot() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 11)) +
  ggforce::geom_circle(aes(x0=0, y0=0, r=1), color="gray50", linetype="dashed") +
  ggforce::geom_circle(aes(x0=0, y0=0, r=0.5), color="gray50", linetype="dashed") +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  geom_point(data = U_scaled, aes(x = x, y = y, fill = Year), 
             shape = 21, color = "white", alpha = 0.35, size = 2.25) +
  geom_point(data = loadings_df, aes(x = x, y = y, color = group), size = 2, shape = 17) +
  geom_text_repel(data = loadings_df, aes(x = x, y = y, label = var_label, color = group),
                  size = 5, max.overlaps=Inf, show.legend = FALSE, parse = TRUE) +
  coord_fixed(xlim = c(-1.01, 1.01), ylim = c(-1.01, 1.01)) +
  labs(
    x = "Canonical loading 1",
    y = "Canonical loading 2"
  ) +
  scale_color_manual(values = c("Environmental variables" = "chocolate4", "Community-related variables" = "green4")) +
  scale_fill_manual(values = colors_years)

# SUPPLEMENTARY FIGURE S2
# scatterplot of sample scores
# showing the relationship between the two sets of canonical variates on the first dimension.
sitescores           <- data.frame(env_diversity$Year, cc_total$scores$xscores[,1], cc_total$scores$yscores[,1])
colnames(sitescores) <- c("Year","environmental_variables","community_variables")
sitescores$Year      <- as.character(sitescores$Year)
ggplot(sitescores, aes(x = environmental_variables, y = community_variables, color = Year)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        panel.grid = element_blank(),
        panel.background = element_rect(color = "white")) +
  geom_point(size = 3) +
  scale_color_manual(values = colors_years3) +
  xlab("Environmental variables") +
  ylab("Community-related variables")

# SUPPLEMENTARY FIGURE S3.a-c
# Scatterplot and correlation matrix
# of environmental and community-related variables
colnames(env_diversity)[c(5,15,17)] <- c("Cond","log(af-TSS)","log(TN)")
ggpairs(env_diversity[,c("Cond","log(af-TSS)","log(TN)","pH","RUE","FEve_pr","FDis_pr")]) +
  theme_bw() +
  theme(panel.grid = element_blank())

colnames(env_diversity)[c(19:21,27)] <- c("FRic_a","FEve_a","FDis_a","J_a")
ggpairs(env_diversity[,c("Cond","log(af-TSS)","log(TN)","pH","FRic_pr","FRic_a")]) +
  theme_bw() +
  theme(panel.grid = element_blank())
ggpairs(env_diversity[,c("Cond","log(af-TSS)","log(TN)","pH","FEve_a","FDis_a","J_a")]) +
  theme_bw() +
  theme(panel.grid = element_blank())

# REGRESSION MODELS
# Best fit model parameters and statistics in Table 2
colnames(env_diversity)[colnames(env_diversity) == "log(af-TSS)"] <- "log_algal_free_TSS"
colnames(env_diversity)[colnames(env_diversity) == "log(TN)"] <- "log_TN"

# (1) RUE as dependent variable
rue_gam    <- gam(RUE ~ s(log_algal_free_TSS), data = env_diversity)
rue_gamm_1 <- gamm(RUE ~ s(log_algal_free_TSS),
                 data = env_diversity,
                 weights = varFixed(~log_algal_free_TSS),
                 method = "ML")
rue_gamm_2 <- gamm(RUE ~ s(log_algal_free_TSS),
                   data = env_diversity,
                   weights = varPower(form = ~log_algal_free_TSS),
                   method = "ML")
rue_gamm_3 <- gamm(RUE ~ s(log_algal_free_TSS) + Cond,
                   data = env_diversity,
                   weights = varPower(form = ~log_algal_free_TSS),
                   method = "ML")
rue_gamm_4 <- gamm(RUE ~ s(log_algal_free_TSS) + s(Cond),
                   data = env_diversity,
                   weights = varPower(form = ~log_algal_free_TSS),
                   method = "ML")
rue_gamm_5 <- gamm(RUE ~ s(log_algal_free_TSS) + log_TN,
                   data = env_diversity,
                   weights = varPower(form = ~log_algal_free_TSS),
                   method = "ML")
rue_gamm_6 <- gamm(RUE ~ s(log_algal_free_TSS) + s(log_TN),
                   data = env_diversity,
                   weights = varPower(form = ~log_algal_free_TSS),
                   method = "ML")
rue_gamm_7 <- gamm(RUE ~ s(log_algal_free_TSS) + log_TN + Cond,
                   data = env_diversity,
                   weights = varPower(form = ~log_algal_free_TSS),
                   method = "ML")
rue_gamm_8 <- gamm(RUE ~ s(log_algal_free_TSS) + log_TN * Cond,
                   data = env_diversity,
                   weights = varPower(form = ~log_algal_free_TSS),
                   method = "ML")
rue_gamm_9 <- gamm(RUE ~ s(log_algal_free_TSS) + s(log_TN) + s(Cond),
                   data = env_diversity,
                   weights = varPower(form = ~log_algal_free_TSS),
                   method = "ML")

rue_lm_1    <- lm(RUE ~ log_algal_free_TSS, data = env_diversity)
rue_lm_2    <- lm(RUE ~ log_algal_free_TSS + Cond, data = env_diversity)
rue_lm_3    <- lm(RUE ~ log_algal_free_TSS * Cond, data = env_diversity)
rue_lm_4    <- lm(RUE ~ log_algal_free_TSS + Cond + TN, data = env_diversity)
rue_lm_5    <- lm(RUE ~ log_algal_free_TSS * Cond * TN, data = env_diversity)
rue_gls_1 <- gls(RUE ~ log_algal_free_TSS,
                 data = env_diversity,
                 weights = varFixed(~log_algal_free_TSS),
                 na.action = na.omit, method = "ML")
rue_gls_2 <- gls(RUE ~ log_algal_free_TSS,
                 data = env_diversity,
                 weights = varPower(form = ~log_algal_free_TSS),
                 na.action = na.omit, method = "ML")
rue_gls_3 <- gls(RUE ~ log_algal_free_TSS + Cond,
                 data = env_diversity,
                 weights = varFixed(~log_algal_free_TSS),
                 na.action = na.omit, method = "ML")
rue_gls_4 <- gls(RUE ~ log_algal_free_TSS * Cond,
                 data = env_diversity,
                 weights = varFixed(~log_algal_free_TSS),
                 na.action = na.omit, method = "ML")
rue_gls_5 <- gls(RUE ~ log_algal_free_TSS + Cond + log_TN,
                 data = env_diversity,
                 weights = varFixed(~log_algal_free_TSS),
                 na.action = na.omit, method = "ML")
rue_gls_6 <- gls(RUE ~ log_algal_free_TSS * Cond * log_TN,
                 data = env_diversity,
                 weights = varFixed(~log_algal_free_TSS),
                 na.action = na.omit, method = "ML")
rue_gls_7 <- gls(RUE ~ log_algal_free_TSS + Cond,
                 data = env_diversity,
                 weights = varPower(form = ~log_algal_free_TSS),
                 na.action = na.omit, method = "ML")
rue_gls_8 <- gls(RUE ~ log_algal_free_TSS * Cond,
                 data = env_diversity,
                 weights = varPower(form = ~log_algal_free_TSS),
                 na.action = na.omit, method = "ML")
rue_gls_9 <- gls(RUE ~ log_algal_free_TSS + Cond + log_TN,
                 data = env_diversity,
                 weights = varPower(form=~log_algal_free_TSS),
                 na.action = na.omit, method = "ML")
rue_gls_10 <- gls(RUE ~ log_algal_free_TSS * Cond * log_TN,
                  data = env_diversity,
                  weights = varPower(form=~log_algal_free_TSS),
                  na.action = na.omit, method = "ML")
AIC(rue_gam,
    rue_gamm_1$lme, rue_gamm_2$lme, rue_gamm_3$lme, rue_gamm_4$lme, rue_gamm_5$lme,
    rue_gamm_6$lme, rue_gamm_7$lme, rue_gamm_8$lme, rue_gamm_9$lme,
    rue_lm_1, rue_lm_2, rue_lm_3, rue_lm_4, rue_lm_5,
    rue_gls_1, rue_gls_2, rue_gls_3, rue_gls_4, rue_gls_5,
    rue_gls_6, rue_gls_7, rue_gls_8, rue_gls_9, rue_gls_10)
summary(rue_gls_7) # model showing the lowest AIC score
# Root mean squared error calculation
pred <- predict(rue_gls_7)[1]
rmse <- sqrt(mean((env_diversity$RUE - pred)^2))
print(rmse)

# (2) FEve as dependent variable
feve_lm_1    <- lm(FEve_pr ~ log_algal_free_TSS, data = env_diversity)
feve_lm_2    <- lm(FEve_pr ~ log_algal_free_TSS + Cond, data = env_diversity)
feve_lm_3    <- lm(FEve_pr ~ log_algal_free_TSS * Cond, data = env_diversity)
feve_lm_4    <- lm(FEve_pr ~ log_algal_free_TSS + Cond + TN, data = env_diversity)
feve_lm_5    <- lm(FEve_pr ~ log_algal_free_TSS * Cond * TN, data = env_diversity)
feve_gls_1 <- gls(FEve_pr ~ log_algal_free_TSS, data = env_diversity,
                   weights = varFixed(~log_algal_free_TSS),
                   na.action = na.omit, method = "ML")
feve_gls_2 <- gls(FEve_pr ~ log_algal_free_TSS,
                  data = env_diversity,
                  weights = varPower(form = ~log_algal_free_TSS),
                  na.action = na.omit, method = "ML")
feve_gls_3 <- gls(FEve_pr ~ log_algal_free_TSS + Cond,
                  data = env_diversity,
                  weights = varFixed(~log_algal_free_TSS),
                  na.action = na.omit, method = "ML")
feve_gls_4 <- gls(FEve_pr ~ log_algal_free_TSS * Cond,
                  data = env_diversity,
                  weights = varFixed(~log_algal_free_TSS),
                  na.action = na.omit, method = "ML")
feve_gls_5 <- gls(FEve_pr ~ log_algal_free_TSS + Cond + log_TN,
                  data = env_diversity,
                  weights = varFixed(~log_algal_free_TSS),
                  na.action = na.omit, method = "ML")
feve_gls_6 <- gls(FEve_pr ~ log_algal_free_TSS * Cond * log_TN,
                  data = env_diversity,
                  weights = varFixed(~log_algal_free_TSS),
                  na.action = na.omit, method = "ML")
feve_gls_7 <- gls(FEve_pr ~ log_algal_free_TSS + Cond,
                  data = env_diversity,
                  weights = varPower(form = ~log_algal_free_TSS),
                  na.action = na.omit, method = "ML")
feve_gls_8 <- gls(FEve_pr ~ log_algal_free_TSS * Cond,
                  data = env_diversity,
                  weights = varPower(form = ~log_algal_free_TSS),
                  na.action = na.omit, method = "ML")
feve_gls_9 <- gls(FEve_pr ~ log_algal_free_TSS + Cond + log_TN,
                  data = env_diversity,
                  weights = varPower(form = ~log_algal_free_TSS),
                  na.action = na.omit, method = "ML")
feve_gls_10 <- gls(FEve_pr ~ log_algal_free_TSS * Cond * log_TN,
                  data = env_diversity,
                  weights = varPower(form = ~log_algal_free_TSS),
                  na.action = na.omit, method = "ML")

AIC(feve_lm_1, feve_lm_2, feve_lm_3, feve_lm_4, feve_lm_5,
    feve_gls_1, feve_gls_2, feve_gls_3, feve_gls_4, feve_gls_5, feve_gls_6,
    feve_gls_7, feve_gls_8, feve_gls_9, feve_gls_10)
summary(feve_gls_1) # model showing the lowest AIC score
# Root mean squared error calculation
pred <- predict(feve_gls_1)[1]
rmse <- sqrt(mean((env_diversity$FEve_pr - pred)^2))
print(rmse)

# (3) FDis as dependent variable
fdis_lm_1  <- lm(FDis_pr ~ log_TN + Cond, data = env_diversity)
fdis_lm_2  <- lm(FDis_pr ~ log_TN * Cond, data = env_diversity)
fdis_gls_1 <- gls(FDis_pr ~ log_TN * Cond, data = env_diversity,
                  weights = varFixed(~log_TN),
                  na.action = na.omit, method = "ML")
fdis_gls_2 <- gls(FDis_pr ~ log_TN * Cond, data = env_diversity,
                  weights = varFixed(~Cond),
                  na.action = na.omit, method = "ML")
fdis_gls_3 <- gls(FDis_pr ~ log_TN * Cond, data = env_diversity,
                  weights = varPower(form = ~log_TN),
                  na.action = na.omit, method = "ML")
fdis_gls_4 <- gls(FDis_pr ~ log_TN * Cond, data = env_diversity,
                  weights = varPower(form = ~Cond),
                  na.action = na.omit, method = "ML")
fdis_gls_5 <- gls(FDis_pr ~ log_TN * Cond, data = env_diversity,
                  weights = varComb(varFixed(~log_TN), varPower(form = ~Cond)),
                  na.action = na.omit, method = "ML")

fdis_gamm_1 <- gamm(FDis_pr ~ log_TN + s(Cond), data = env_diversity,
                    weights = varFixed(~log_TN),
                    na.action = na.omit, method = "ML")
fdis_gamm_2 <- gamm(FDis_pr ~ log_TN + s(Cond), data = env_diversity,
                    weights = varPower(form = ~log_TN),
                    na.action = na.omit, method = "ML")

AIC(fdis_lm_1, fdis_lm_2,
    fdis_gls_1, fdis_gls_2, fdis_gls_3, fdis_gls_4, fdis_gls_5,
    fdis_gamm_1$lme, fdis_gamm_2$lme)
summary(fdis_gls_1)  # model showing the lowest AIC score
# Root mean squared error calculation
pred <- predict(fdis_gls_1)[1]
rmse <- sqrt(mean((env_diversity$FEve_pr - pred)^2))
print(rmse)

# FIGURE 4
# FIGURE 4.a
# RUE vs. algal-free TSS and conductivity
# For median conductivity (main model fit, solid line)
x_seq <- seq(min(env_diversity$log_algal_free_TSS, na.rm=TRUE),
             max(env_diversity$log_algal_free_TSS, na.rm=TRUE), length.out=100)
cond_median <- median(env_diversity$Cond, na.rm=TRUE)
newdat <- data.frame(log_algal_free_TSS = x_seq, Cond = cond_median)

preds <- predictSE.gls(rue_gls_7, newdat, se.fit=TRUE)
newdat$fit <- preds$fit
newdat$se.fit <- preds$se.fit
alpha <- 0.05
z <- qnorm(1 - alpha/2)
newdat$lower <- newdat$fit - z * newdat$se.fit
newdat$upper <- newdat$fit + z * newdat$se.fit

ggplot(aes(x = log_algal_free_TSS, y = RUE, color = Cond), data = env_diversity) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 13, vjust = -0.4),
        axis.title.y = element_text(size = 13, vjust = 2),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.15,0.8),
        legend.background = element_rect(colour = "black", linewidth = 0.2)) +
  geom_point(size = 3) +
  scale_color_viridis_c(end = 0.9, option = "G", direction = -1) +
  xlab("log10(Algal-free TSS) (mg/L)") +
  ylab("RUE") +
  guides(color = guide_colorbar(title = "Conductivity \n(mS/cm)")) +
  geom_ribbon(data = newdat, aes(x = log_algal_free_TSS, ymin = lower, ymax = upper),
              inherit.aes = FALSE, fill = "grey50", alpha = 0.2) +
  geom_line(data = newdat, aes(x = log_algal_free_TSS, y = fit),
            inherit.aes = FALSE, color = "black", linewidth = 0.5)

# FIGURE 4.b
# FEve_pr vs. algal-free TSS
newdat <- data.frame(log_algal_free_TSS = x_seq)
preds <- predictSE.gls(feve_gls_1, newdat, se.fit=TRUE)
newdat$fit <- preds$fit
newdat$se.fit <- preds$se.fit
alpha <- 0.05
z <- qnorm(1 - alpha/2)
newdat$lower <- newdat$fit - z * newdat$se.fit
newdat$upper <- newdat$fit + z * newdat$se.fit

ggplot(aes(x = log_algal_free_TSS, y = FEve_pr), data = env_diversity) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 13, vjust = -0.4),
        axis.title.y = element_text(size = 13, vjust = 2),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "inside",
        legend.position.inside = c(0.15,0.8),
        legend.background = element_rect(colour = "black", linewidth = 0.2)) +
  geom_point(size = 3) +
  scale_color_viridis_c(end = 0.9, option = "G", direction = -1) +
  xlab("log10(Algal-free TSS) (mg/L)") +
  ylab(expression(FEve[pr])) +
  guides(color = guide_colorbar(title = "Conductivity (mS/cm)")) +
  geom_ribbon(data = newdat, aes(x = log_algal_free_TSS, ymin = lower, ymax = upper),
              inherit.aes = FALSE, fill = "grey50", alpha = 0.2) +
  geom_line(data = newdat, aes(x = log_algal_free_TSS, y = fit),
            inherit.aes = FALSE, color = "black", linewidth = 0.5)

# FIGURE 4.c
# FDis_pr vs. TN and conductivity
# For median conductivity (main model fit, solid line)
x_seq <- seq(min(env_diversity$log_TN, na.rm=TRUE),
             max(env_diversity$log_TN, na.rm=TRUE), length.out=100)
newdat <- data.frame(log_TN = x_seq, Cond = cond_median)

preds <- predictSE.gls(fdis_gls_1, newdat, se.fit=TRUE)
newdat$fit <- preds$fit
newdat$se.fit <- preds$se.fit
alpha <- 0.05
z <- qnorm(1 - alpha/2)
newdat$lower <- newdat$fit - z * newdat$se.fit
newdat$upper <- newdat$fit + z * newdat$se.fit

ggplot(aes(x = log_TN, y = FDis_pr, color = Cond), data = env_diversity) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 13, vjust = -0.4),
        axis.title.y = element_text(size = 13, vjust = 2),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.85,0.15),
        legend.background = element_rect(colour = "black", linewidth = 0.2)) +
  geom_point(size = 3) +
  scale_color_viridis_c(end = 0.9, option = "G", direction = -1) +
  xlab(expression("log10(TN) ("*mu*"g/L)")) +
  ylab(expression(FDis[pr])) +
  guides(color = guide_colorbar(title = "Conductivity \n(mS/cm)")) +
  geom_ribbon(data = newdat, aes(x = log_TN, ymin = lower, ymax = upper),
              inherit.aes = FALSE, fill = "grey50", alpha = 0.2) +
  geom_line(data = newdat, aes(x = log_TN, y = fit),
            inherit.aes = FALSE, color = "black", linewidth = 0.5)

# FIGURE 5
# Quantile regression
# Definition of Bragg equation used for nonlinear (unimodal) relationships
# source: https://www.statforbiology.com/2020/stat_nls_selfstarting/
bragg.3.fun <- function(X, b, d, e){
  d * exp(- b * (X - e)^2)
}

bragg.3.init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["X"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  
  d <- max(y)
  e <- x[which.max(y)]
  
  pseudoY <- log( y / d )
  pseudoX <- (x - e)^2
  coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
  b <- - coefs[1]
  start <- c(b, d, e)
  names(start) <- mCall[c("b", "d", "e")]
  start
}

# self-starting function for initial model parameter values
NLS.bragg.3 <- selfStart(bragg.3.fun, bragg.3.init, parameters=c("b", "d", "e"))

# Quantile regression on FDis_pr as independent variable
min_FDis <- min(env_diversity$FDis_pr)
max_FDis <- max(env_diversity$FDis_pr)

nlrq_rue_fdis       <- nlrq(RUE ~ NLS.bragg.3(FDis_pr, b, d, e), data = env_diversity, tau = 0.8)
pred_nlrq_rue_fdis  <- predict(nlrq_rue_fdis,
                               newdata = list(FDis_pr = seq(min_FDis, max_FDis, 0.001)))
pred_nlrq_rue_fdis  <- cbind(pred_nlrq_rue_fdis, FDis_pr = seq(min_FDis, max_FDis, 0.001))

poly_df <- rbind(data.frame(FDis_pr = min_FDis, pred_nlrq_rue_fdis = 0),
                 pred_nlrq_rue_fdis,
                 data.frame(FDis_pr = max_FDis, pred_nlrq_rue_fdis = 0))

env_diversity$Year <- as.factor(env_diversity$Year)

ggplot(aes(x = FDis_pr, y = RUE, color = Year), data = env_diversity) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 12, vjust = -0.4),
        axis.title.y = element_text(size = 12, vjust = 2),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  geom_polygon(aes(x = FDis_pr, y = pred_nlrq_rue_fdis), color = NA,
               fill = "palegreen", alpha = 0.5, data = poly_df) +
  geom_point(size = 2) +
  scale_color_manual(values = colors_years3) +
  geom_path(aes(x = FDis_pr, y = pred_nlrq_rue_fdis),
            color = "grey25", linewidth = 0.5, alpha = 0.75, linetype = "dashed",
            data = pred_nlrq_rue_fdis) +
  ylab("Resource use efficiency") +
  xlab(expression("FDis"["pr"])) +
  guides(color = guide_legend(title = "Year"))

# TABLE 3
# Summary statistics
summary(nlrq_rue_fdis)

# SUPPLEMENTARY FIGURE S4
# Quantile regression on FEve_pr as independent variable
min_FEve <- min(env_diversity$FEve_pr)
max_FEve <- max(env_diversity$FEve_pr)

nlrq_rue_FEve       <- nlrq(RUE ~ NLS.bragg.3(FEve_pr, b, d, e), data = env_diversity, tau = 0.8)
pred_nlrq_rue_FEve  <- predict(nlrq_rue_FEve,
                               newdata = list(FEve_pr = seq(min_FEve, max_FEve, 0.001)))
pred_nlrq_rue_FEve  <- cbind(pred_nlrq_rue_FEve, FEve_pr = seq(min_FEve, max_FEve, 0.001))

poly_df <- rbind(data.frame(FEve_pr = min_FEve, pred_nlrq_rue_FEve = 0),
                 pred_nlrq_rue_FEve,
                 data.frame(FEve_pr = max_FEve, pred_nlrq_rue_FEve = 0))

ggplot(aes(x = FEve_pr, y = RUE, color = Year), data = env_diversity) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 12, vjust = -0.4),
        axis.title.y = element_text(size = 12, vjust = 2),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  geom_polygon(aes(x = FEve_pr, y = pred_nlrq_rue_FEve), color = NA,
               fill = "palegreen", alpha = 0.5, data = poly_df) +
  geom_point(size = 2) +
  scale_color_manual(values = colors_years3) +
  geom_path(aes(x = FEve_pr, y = pred_nlrq_rue_FEve),
            color = "grey25", linewidth = 0.5, alpha = 0.75, linetype = "dashed",
            data = pred_nlrq_rue_FEve) +
  ylab("Resource use efficiency") +
  xlab(expression("FEve"["pr"])) +
  guides(color = guide_legend(title = "Year"))

# SUPPLEMENTARY TABLE S3
# Summary statistics
summary(nlrq_rue_FEve)
