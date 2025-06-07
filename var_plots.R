library(ggplot2)
library(dplyr)
library(stringr)

############################################################################################

# Definir o diretório de trabalho
setwd("C:/path/")


dados_vcf <- read.csv("file.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)

#187547917


# Verificar as primeiras linhas do dataframe combinado
head(dados_vcf)

# Verificar a estrutura do dataframe combinado
str(dados_vcf)


mapeamento <- c(
  "827" = "rs28358569",
  "1095" = "rs267606618",
  "1189" = "rs28358571",
  "1243" = "rs28358572",
  "1555" = "rs267606617",
  "10398" = "rs2853826",
  "11778" = "rs199476112",
  "14484" = "rs199476104"
)


dados_vcf$Variant_Name <- mapeamento[as.character(dados_vcf$Position)]
dados_vcf <- dados_vcf[dados_vcf$Variant_Name != "" & !is.na(dados_vcf$Variant_Name), ]
###################################################################################################

data_haplog <- "sabe_af_haplogroup.csv"
data_haplog <- read.csv(data_haplog, header=TRUE, sep=";", stringsAsFactors=FALSE)
data_haplog$SampleID <- as.character(data_haplog$SampleID)
data_haplog$Haplogroup <- as.character(data_haplog$Haplogroup)
data_haplog$Cluster <- as.character(data_haplog$Cluster)

data_autoss <- "sabe_aut.csv"
data_autoss <- read.csv(data_autoss, header=TRUE, sep=";", stringsAsFactors=FALSE)
data_autoss$SampleID <- as.character(data_autoss$SampleID)
data_autoss$AUT <- as.character(data_autoss$AUT)
###################################################################################################

dados_vcf$SampleID <- as.character(dados_vcf$SampleID)
dados_init <- left_join(dados_vcf, data_haplog %>% select(SampleID, Haplogroup, Cluster), by = "SampleID")
dados <- left_join(dados_init, data_autoss %>% select(SampleID, AUT), by = "SampleID")

dados$origin <- ifelse(nchar(dados$SampleID) == 9, "SABE", "BIPMed")

dados <- dados %>%
  mutate(SampleID = str_replace(SampleID, "_p.*$", ""))

###################################################################################################

color_palette <- c("A" = "#FF9999",
                   "B" = "#FF9999", 
                   "C" = "#FF9999", 
                   "D" = "#FF9999", 
                   "G" = "#80EF80", 
                   "H" = "#FFFFB5",
                   "HV" = "#FFFFB5",
                   "I" = "#FFFFB5", 
                   "J" = "#FFFFB5", 
                   "K" = "#FFFFB5", 
                   "L0" = "#B3EBF2", 
                   "L1" = "#B3EBF2",
                   "L2" ="#B3EBF2", 
                   "L3" = "#B3EBF2", 
                   "L4" = "#B3EBF2", 
                   "M" = "#80EF80", 
                   "N" = "#80EF80", 
                   "R" = "#FFFFB5", 
                   "T" = "#FFFFB5",
                   "U" = "#FFFFB5",
                   "V" = "#FFFFB5",
                   "X" = "#FFFFB5",
                   "W" = "#FFFFB5")


ggplot(data = data_haplog, aes(x = Cluster, fill = Cluster)) +
  geom_bar(color = "black") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +  # Updated here
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Cluster", y = "Number of samples", title = "Haplogroups Cluster Distribution - AF") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(breaks = seq(0, 300, by = 50), expand = c(0, 0), limits = c(0, 300))

data_haplog <- data_haplog %>%
  mutate(Cluster_g = case_when(
    Cluster %in% c("A", "B", "C", "D") ~ "AMR",
    Cluster %in% c("G", "M", "N") ~ "ASI",
    Cluster %in% c("L0", "L1", "L2", "L3", "L4") ~ "AFR",
    Cluster %in% c("H", "HV", "I", "J", "K", "R", "T", "U", "V", "W", "X") ~ "EUR",
  ))


###################################################################################################

# Filtrar os dados para a variante específica
variant_rs28358569 <- dados %>% filter(Variant_Name == "rs28358569")

# Contar o número de amostras por Haplogroup e Cluster
count_rs28358569 <- variant_rs28358569 %>%
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

ggplot(count_rs28358569, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs28358569 (m.827A>G)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 70, by=10), expand=c(0, 0), limits=c(0, 70))

###################################################################################################

# Filtrar os dados para a variante específica
variant_rs28358571 <- dados %>% filter(Variant_Name == "rs28358571")

# Contar o número de amostras por Haplogroup e Cluster
count_rs28358571 <- variant_rs28358571 %>% 
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

ggplot(count_rs28358571, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs28358571 (m.1189T>C)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 5, by=1), expand=c(0, 0), limits=c(0, 5))

###################################################################################################

# Filtrar os dados para a variante específica
variant_rs28358572 <- dados %>% filter(Variant_Name == "rs28358572")

# Contar o número de amostras por Haplogroup e Cluster
count_rs28358572 <- variant_rs28358572 %>% 
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

ggplot(count_rs28358572, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs28358572 (m.1243T>C)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 5, by=1), expand=c(0, 0), limits=c(0, 5))

###################################################################################################

# Filtrar os dados para a variante específica
variant_rs2853826 <- dados %>% filter(Variant_Name == "rs2853826")

# Contar o número de amostras por Haplogroup e Cluster
count_rs2853826 <- variant_rs2853826 %>% 
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

ggplot(count_rs2853826, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G)") +
  theme(axis.text.x = element_blank()) +  # Remove os nomes dos haplogrupos
  scale_y_continuous(breaks=seq(0, 50, by=10), expand=c(0, 0), limits=c(0, 50))

ggplot(count_rs2853826, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 50, by=10), expand=c(0, 0), limits=c(0, 50))



##################################################################################################

filtered_A <- count_rs2853826 %>%
  filter(Cluster == "A")

# Criar o gráfico com os dados filtrados
ggplot(filtered_A, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster A") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 5, by=1), expand=c(0, 0), limits=c(0, 5))

##################################################################################################

filtered_C <- count_rs2853826 %>%
  filter(Cluster == "C")

# Criar o gráfico com os dados filtrados
ggplot(filtered_C, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster C") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 50, by=10), expand=c(0, 0), limits=c(0, 50))

##################################################################################################

filtered_D <- count_rs2853826 %>%
  filter(Cluster == "D")

# Criar o gráfico com os dados filtrados
ggplot(filtered_D, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster D") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 30, by=10), expand=c(0, 0), limits=c(0, 30))

##################################################################################################

filtered_G <- count_rs2853826 %>%
  filter(Cluster == "G")

# Criar o gráfico com os dados filtrados
ggplot(filtered_G, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster G") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 3, by=1), expand=c(0, 0), limits=c(0, 3))

##################################################################################################

filtered_I <- count_rs2853826 %>%
  filter(Cluster == "I")

# Criar o gráfico com os dados filtrados
ggplot(filtered_I, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster I") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 3, by=1), expand=c(0, 0), limits=c(0, 3))

##################################################################################################

filtered_J <- count_rs2853826 %>%
  filter(Cluster == "J")

# Criar o gráfico com os dados filtrados
ggplot(filtered_J, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster J") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 5, by=1), expand=c(0, 0), limits=c(0, 5))

##################################################################################################

filtered_K <- count_rs2853826 %>%
  filter(Cluster == "K")

# Criar o gráfico com os dados filtrados
ggplot(filtered_K, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster K") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 5, by=1), expand=c(0, 0), limits=c(0, 5))

##################################################################################################

filtered_L0 <- count_rs2853826 %>%
  filter(Cluster == "L0")

# Criar o gráfico com os dados filtrados
ggplot(filtered_L0, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster L0") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 10, by=1), expand=c(0, 0), limits=c(0, 10))

##################################################################################################

filtered_L1 <- count_rs2853826 %>%
  filter(Cluster == "L1")

# Criar o gráfico com os dados filtrados
ggplot(filtered_L1, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster L1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 20, by=2), expand=c(0, 0), limits=c(0, 20))

##################################################################################################

filtered_L2 <- count_rs2853826 %>%
  filter(Cluster == "L2")

# Criar o gráfico com os dados filtrados
ggplot(filtered_L2, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster L2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 15, by=1), expand=c(0, 0), limits=c(0, 15))

##################################################################################################

filtered_L3 <- count_rs2853826 %>%
  filter(Cluster == "L3")

# Criar o gráfico com os dados filtrados
ggplot(filtered_L3, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster L3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 10, by=1), expand=c(0, 0), limits=c(0, 10))

##################################################################################################

filtered_L4 <- count_rs2853826 %>%
  filter(Cluster == "L4")

# Criar o gráfico com os dados filtrados
ggplot(filtered_L4, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster L4") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 3, by=1), expand=c(0, 0), limits=c(0, 3))

##################################################################################################

filtered_M <- count_rs2853826 %>%
  filter(Cluster == "M")

# Criar o gráfico com os dados filtrados
ggplot(filtered_M, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster M") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 3, by=1), expand=c(0, 0), limits=c(0, 3))

##################################################################################################

filtered_N <- count_rs2853826 %>%
  filter(Cluster == "N")

# Criar o gráfico com os dados filtrados
ggplot(filtered_N, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G) - Cluster N") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 3, by=1), expand=c(0, 0), limits=c(0, 3))

##################################################################################################

# Filtrar os dados para a variante específica
variant_rs267606617 <- dados %>% filter(Variant_Name == "rs267606617")

# Contar o número de amostras por Haplogroup e Cluster
count_rs267606617 <- variant_rs267606617 %>%
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

ggplot(count_rs267606617, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs267606617 (m.1555A>G)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 3, by=1), expand=c(0, 0), limits=c(0, 3))

##################################################################################################

# Filtrar os dados para a variante específica
variant_rs28665937 <- dados %>% filter(Variant_Name == "rs28665937")

# Contar o número de amostras por Haplogroup e Cluster
count_rs28665937 <- variant_rs28665937 %>%
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

ggplot(count_rs28665937, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs28665937 (m.6413T>C)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 5, by=1), expand=c(0, 0), limits=c(0, 5))


##################################################################################################

# Filtrar os dados para a variante específica
variant_rs199476112 <- dados %>% filter(Variant_Name == "rs199476112")

# Contar o número de amostras por Haplogroup e Cluster
count_rs199476112 <- variant_rs199476112 %>%
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

ggplot(count_rs199476112, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs199476112 (m.11778G>A)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 3, by=1), expand=c(0, 0), limits=c(0, 3))

##################################################################################################

# Filtrar os dados para a variante específica
variant_rs1603223785 <- dados %>% filter(Variant_Name == "rs1603223785")

# Contar o número de amostras por Haplogroup e Cluster
count_rs1603223785 <- variant_rs1603223785 %>%
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

ggplot(count_rs1603223785, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs1603223785 (m.12530A>G)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 3, by=1), expand=c(0, 0), limits=c(0, 3))

##################################################################################################

# Filtrar os dados para a variante específica
variant_rs199476104 <- dados %>% filter(Variant_Name == "rs199476104")

# Contar o número de amostras por Haplogroup e Cluster
count_rs199476104 <- variant_rs199476104 %>%
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

ggplot(count_rs199476104, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs199476104 (m.14484T>C)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 3, by=1), expand=c(0, 0), limits=c(0, 3))

##################################################################################################

# Filtrar os dados para a variante específica
variant_rs267606618 <- dados %>% filter(Variant_Name == "rs267606618")

# Contar o número de amostras por Haplogroup e Cluster
count_rs267606618 <- variant_rs267606618 %>%
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

ggplot(count_rs267606618, aes(x = Haplogroup, y = Samples, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Haplogroup", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs267606618 (m.1095T>C)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 3, by=1), expand=c(0, 0), limits=c(0, 3))

##############################################################


# Carregar o pacote dplyr
library(dplyr)

# Vamos supor que o dataframe seja 'data_haplog', que já está carregado com seus dados.

# Contar o número de haplogrupos únicos por cluster
result <- data_haplog %>%
  group_by(Cluster) %>%
  summarise(count_haplogroups = n())

# Exibir o resultado
print(result)

# Contar haplogrupos únicos por cluster
result_unique <- data_haplog %>%
  group_by(Cluster) %>%
  summarise(count_haplogroups = n_distinct(Haplogroup))

# Exibir o resultado
print(result_unique)












###################################################################################################

# Filtrar os dados para a variante específica
variant_rs2853826 <- dados %>% filter(Variant_Name == "rs2853826")

# Contar o número de amostras por Haplogroup e Cluster
count_rs2853826 <- variant_rs2853826 %>% 
  group_by(Haplogroup, Cluster) %>% 
  summarise(Samples = n()) %>% 
  ungroup()

new_dataset <- count_rs2853826 %>%
  group_by(Cluster) %>%
  summarise(Ninds = sum(Samples)) %>%
  ungroup()

ggplot(new_dataset, aes(x = Cluster, y = Ninds, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = color_palette, name = "Cluster", na.value = "lightgrey") +
  labs(x = "Cluster", y = "Number of samples") +
  theme_bw() +
  labs(title = "rs2853826 (m.10398A>G)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks=seq(0, 150, by=10), expand=c(0, 0), limits=c(0, 150))
###################################################################################################
