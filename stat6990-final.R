# The following program holds the analysis for \Evaluating the effect of New
# Zealands Earthquake Commission Act using natural disasterdamage claims\ by Barnabas Amlalo,
# Rory McDermott, and Tripp Modzelewski. Last updated 04/29/2026.
# 
# This code does all data manipulation, runs all tests, and produces all figures for the final paper.

library(CASdatasets)
library(dplyr)
library(readxl)
library(ggplot2)
library(MatchIt)
library(tidyr)
library(stringr)
library(lubridate)
library(car)
library(CausalImpact)


###Load in all data
nzcathist <- read_excel("nzcathist.xlsx")
auscathist <- read_excel("auscathist.xlsx")

### Create outcome variable that is cost accrued, inflation-adjusted to 2026 USD as of 4/20/2026
nzcathist$USDtoday <- nzcathist$NormCost2014*1000/827.003
auscathist$USDtoday <- auscathist$NormCost2014*1000/995.645

### combining datasets with country indicator
nzcathist$country <- "NZ"
auscathist$country <- "AUS"
cols <- c("Year", "Quarter", "FirstDay", "Event", "Type", "Location", "USDtoday", "country", "exposure")
data <- rbind(nzcathist[cols], auscathist[cols])

#drop non-natural disaster events
data <- data %>% filter(Type != "Other") %>% filter(Type != "Power outage") 

data$treated <- ifelse(data$country == "NZ", 1, 0) #Is obs in treated population (in NZ)
data$log_exposure <- log1p(data$exposure) #reduce influence of highly populated areas
data$Type <- as.factor(data$Type)

data <- data %>% drop_na(log_exposure, Year, Type) %>% filter(exposure>0)



#### Begin work toward DiD and matching

#balance check, Looking for SMD < 0.25 for continuous effects. Type effects are all 0 because of exact matching. best caliper size ~ 0.3
for(i in c(0.2, 0.3)){
matched <- matchit(treated ~ log_exposure + Year, data = data, exact= ~ Type, replace = FALSE, caliper = i)
print(summary(matched))
}
matched.data <- match.data(matched)
plot(matched, type = "jitter")
plot(matched, type = "qq")

# post-event indicator
matched.data <- matched.data %>% mutate(post = Year >= 1994)
data <- data %>% mutate(post = Year >= 1994)

#Figure 2
ggplot(data, aes(x = Year, y = log1p(USDtoday), colour = country)) +
  geom_point() +
  geom_smooth(
    data = subset(data, post == FALSE),
    se = TRUE
  ) +
  geom_smooth(
    data = subset(data, post == TRUE),
    se = TRUE
  ) + scale_color_manual(values = c("AUS" = "red", "NZ" = "blue")) +
  geom_vline(xintercept = 1994) + theme_classic() + ggtitle("Event cost and year by country with loess curve fit") + ylab("ln(1 + cost) (cost in millions, 2026 USD)")

# figure 3
ggplot(matched.data, aes(x = Year, y = log1p(USDtoday), colour = country)) +
  geom_point() +
  geom_smooth(
    data = subset(matched.data, post == FALSE),
    se = TRUE
  ) +
  geom_smooth(
    data = subset(matched.data, post == TRUE),
    se = TRUE
  ) + scale_color_manual(values = c("AUS" = "red", "NZ" = "blue")) +
  geom_vline(xintercept = 1994) + theme_classic() + ggtitle("Event cost and year by country with loess curve fit") + ylab("ln(1 + cost) (cost in millions, 2026 USD)")


# DiD on matched dataset
did <- lm(
  log1p(USDtoday) ~ treated*post + Year + Type + log_exposure, data = matched.data, weights = weights
)
summary(did)
# DiD on full dataset
did <- lm(
  log1p(USDtoday) ~ treated*post + Year + Type + log_exposure, data = data  
)
summary(did)


#parallel trends check
pre_data <- matched.data %>% filter(post == FALSE)

summary(lm(log1p(USDtoday) ~ treated*Year + Type + log_exposure, data = pre_data))

##############################################################################

####### Interrupted Time Series Analyses

####### Analysis 1: Segmented Regression

####### SR Part 1: Data Manipulation

nz_event_dummy <- nzcathist$Year > 1994
nz_event_dummy <- as.numeric(nz_event_dummy)

nz_log_exposure <- log1p(nzcathist$exposure)

nzcathist_enriched <- nzcathist
nzcathist_enriched$dummy <- nz_event_dummy
nzcathist_enriched$log_exposure <- nz_log_exposure

head(nzcathist_enriched[, -c(1,2,3)]) 

####### SR Part 2: Model Building

es_lm <- lm( log1p(USDtoday) ~ Year + dummy + Year:dummy + Type + log_exposure, data = nzcathist_enriched )
es_lm
summary(es_lm)

####### SR Part 2: Assumptions


es_dwtest <- durbinWatsonTest(es_lm)
es_dwtest

ggplot(nzcathist_enriched, aes( y = log1p(USDtoday), x = Year) )+
  geom_smooth() + geom_vline(xintercept = 1994) + geom_point()

####### Analysis 2: Bayesian Structural Time Series

####### BSTS Part 1: Data Manipulation

nzcathist_enriched_before <- nzcathist_enriched[nzcathist_enriched$dummy == 0,]
nzcathist_enriched_after <- nzcathist_enriched[nzcathist_enriched$dummy == 1,]

X_before = nzcathist_enriched_before[, c("Type", "log_exposure", "Year")]
X_after = nzcathist_enriched_after[, c("Type", "log_exposure", "Year")]

n_before = nrow(nzcathist_enriched_before)
n_after = nrow(nzcathist_enriched_after)
n = nrow(nzcathist_enriched)

####### BSTS Part 2: Model Building


type_dummies <- model.matrix(~ nzcathist_enriched$Type - 1)

nzcathist_final <- nzcathist_enriched[,c("USDtoday", "log_exposure", "Year")]
nzcathist_final$USDtoday <- log1p(nzcathist_final$USDtoday)
nzcathist_final <- cbind(nzcathist_final, type_dummies)
head(nzcathist_final)

flipped_matrix <- apply(nzcathist_final, 2, rev)
flipped_matrix

model_2 <- CausalImpact(flipped_matrix, c(1,n_before), c( (1+n_before), n_after ) )

summary(model_2)

plot(model_2)


##############################################################################
# 1. DESCRIPTIVE STATISTICS TABLES

# Summary stats for Exposure (population at risk)
exposure_summary <- data %>%
  group_by(country) %>%
  summarise(
    mean_exposure = mean(exposure, na.rm = TRUE),
    median_exposure = median(exposure, na.rm = TRUE),
    sd_exposure = sd(exposure, na.rm = TRUE),
    min_exposure = min(exposure, na.rm = TRUE),
    max_exposure = max(exposure, na.rm = TRUE)
  )
print(exposure_summary)

# Summary stats for NormCost (claims severity in USD)
normcost_summary <- data %>%
  group_by(country) %>%
  summarise(
    mean_cost = mean(USDtoday, na.rm = TRUE),
    median_cost = median(USDtoday, na.rm = TRUE),
    sd_cost = sd(USDtoday, na.rm = TRUE),
    total_cost = sum(USDtoday, na.rm = TRUE)
  )

print(normcost_summary)

# Summary stats by Disaster Type
type_summary <- data %>%
  group_by(country, Type) %>%
  summarise(
    n_events = n(),
    mean_cost = mean(USDtoday, na.rm = TRUE),
    total_cost = sum(USDtoday, na.rm = TRUE),
    mean_exposure = mean(exposure, na.rm = TRUE)
  )

print(type_summary)



# 2. EXPOSURE VISUAL (Population at risk)

# Boxplot of Exposure by Country
p1 <- ggplot(data, aes(x = country, y = exposure, fill = country)) +
  geom_boxplot() +  # Log scale because exposure varies widely
  labs(
    title = "Population Exposure to Natural Disasters",
    subtitle = "NZ vs Australia",
    x = "Country",
    y = "Exposure (population at risk)"  ) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("NZ" = "darkblue", "AUS" = "darkred"))

print(p1)



# 3. CLAIMS SEVERITY (NormCost) VISUAL


# Boxplot of USDtoday by Country
p2 <- ggplot(data, aes(x = country, y = log1p(USDtoday), fill = country)) +
  geom_boxplot() +
  labs(
    title = "Claims Severity (USD)",
    subtitle = "NZ vs Australia (log scale)",
    x = "Country",
    y = "ln(1 + cost) (cost in millions)",
    caption = "Note: Inflation-adjusted to 2026 USD"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("NZ" = "darkblue", "AUS" = "darkred"))

print(p2)

# Histogram of claims severity distribution
p3 <- ggplot(data, aes(x = USDtoday, fill = country)) +
  geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
  scale_x_log10() +
  labs(
    title = "Distribution of Claims Severity",
    subtitle = "NZ (blue) vs Australia (red)",
    x = "Claim Amount (USD, log scale)",
    y = "Number of Events"
  ) +
  theme_classic() +
  scale_fill_manual(values = c("NZ" = "darkblue", "AUS" = "darkred"))

print(p3)


# 4. DISASTER TYPE VISUAL


# Barplot: Number of events by Type and Country
p4 <- ggplot(type_summary, aes(x = Type, y = n_events, fill = country)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Number of Disaster Events by Type",
    subtitle = "NZ vs Australia",
    x = "Disaster Type",
    y = "Number of Events"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("NZ" = "darkblue", "AUS" = "darkred"))

print(p4)

# Barplot: Total claims cost by Type and Country
p5 <- ggplot(type_summary, aes(x = Type, y = total_cost, fill = country)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Total Claims Cost by Disaster Type",
    subtitle = "NZ vs Australia (in millions USD)",
    x = "Disaster Type",
    y = "Total Cost (Millions USD)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("NZ" = "darkblue", "AUS" = "darkred"))

print(p5)


# 5. EARTHQUAKE-SPECIFIC COMPARISON


# Filter for Earthquake events only
quake_data <- data %>% filter(Type == "Earthquake")

# Barplot: Number of earthquakes by country
p6 <- ggplot(quake_data, aes(x = country, fill = country)) +
  geom_bar() +
  labs(
    title = "Number of Earthquake Events",
    subtitle = "NZ vs Australia",
    x = "Country",
    y = "Number of Earthquakes"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("NZ" = "darkblue", "AUS" = "darkred"))

print(p6)

# Barplot: Total claims cost from earthquakes by country
quake_cost_summary <- quake_data %>%
  group_by(country) %>%
  summarise(
    n_quakes = n(),
    total_cost_millions = sum(USDtoday, na.rm = TRUE) ,
    avg_cost_millions = mean(USDtoday, na.rm = TRUE) 
  )


print(quake_cost_summary)

p7 <- ggplot(quake_cost_summary, aes(x = country, y = total_cost_millions, fill = country)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Total Earthquake Claims Cost",
    subtitle = "NZ vs Australia (in millions USD)",
    x = "Country",
    y = "Total Cost (Millions USD)"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("NZ" = "darkblue", "AUS" = "darkred"))

print(p7)
