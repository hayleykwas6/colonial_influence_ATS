# 🌍 Antarctic Treaty System Participation Study

## 🧊 Overview
This repository contains the data and reproducible code for the research presented.

Our study examines global patterns of Antarctic Treaty System (ATS) participation, 
identifying three distinct clusters of State engagement and demonstrating that colonial 
history is the strongest correlate of participation among the variables tested.

## 📊 Key Findings
- 🗺️ Three distinct clusters of ATS participation identified across 204 States
- 🏛️ Colonial status explained the most variance of ATS participation
- 🌐 NATO membership, BRICS membership, continent, and State age are also significant explainers
- 🪖 Military spending (% of budget) significantly correlates with ATS particiption, but did not signficantly drive clustering
- 💰 Other wealth and spending patterns show no significant association with participation

## 🗂️ Repository Structure
```r
ATS_participation
├── 📄 ATS_Data.csv       # State-level data on ATS participation and predictor variables
├── 📄 ATS_Code.R         # Complete analysis script
└── 📄 README.md          # This file
```

## 🔧 Prerequisites
**Required Software**
- R (≥ 4.5.2)

**Required R Packages**
```r
install.packages(c("vegan", "car", "dplyr", "ggplot2", 
                   "factoextra", "cluster", "patchwork",
                   "countrycode", "maps", "mice"))
```

## 🚀 Getting Started
**Step 1: Clone the Repository**
```bash
git clone [repository URL]
cd ATS_participation
```

**Step 2: Update the Data Path**

In `ATS_code.R`, update the file path to match your local directory:
```r
df <- read.csv("your/path/here/ATSData.csv")
```

**Step 3: Run the Analysis**

Open `ATS_code.R` in R and run sequentially. The script is organized into 
five sections:
1. Load and prep data
2. Create and analyze Gower distance matrix (NMDS + PERMANOVA)
3. Conduct cluster analysis
4. Identify drivers of clusters
5. Visualize global clusters

## 📋 Data Description
`ATSData.csv` contains 204 countries with the following variables:

| Variable | Description |
|---|---|
| `AT` | Antarctic Treaty status (CP, NonCP, Neither) |
| `CCAMLR` | CCAMLR membership (Member, Acceding, Neither) |
| `Colonization.Status` | Colonial history (Colonizer, Colonized, Both, Neither) |
| `In.NATO` | NATO membership (Yes, No, Partner) |
| `In.BRICS` | BRICS membership (Yes, No) |
| `Continent` | Continent |
| `IndepEstYr` | Year of independence or establishment |
| `Age` | State age (years) |
| `GDP` | Total GDP (USD) |
| `GDPperCap` | GDP per capita (USD) |
| `Edu_Total` | Total education spending (USD) |
| `Edu_Per` | Education spending as % of GDP |
| `Mil_Total` | Total military spending (USD) |
| `Mil_Per` | Military spending as % of GDP |
| `NatResRents` | Natural resource rents as % of GDP |

Colonial status is defined as holding colonies during or after 1870, capturing 
the period of peak colonial expansion and its aftermath.

## 🔬 Methods Summary
- Gower dissimilarity matrix constructed from AT and CCAMLR status
- NMDS ordination (k=2)
- Hierarchical clustering using Ward's method (k=3)
- PERMANOVA with 9,999 permutations (vegan::adonis2)
- Collinearity assessed via VIF (all values < 5.5)

## 📚 Associated Publication
**Citation:**
[Authors]. ([Year]). [Title]. [Journal].

## 👥 Authors & Contact
**Authors:** Hayley Kwasniewski, Dr. Peyton Thomas, Dr. Vasco Chavez-Molina, Rosie Sanchez, and Dr. Cassandra Brooks

📧 Corresponding author contact: hayley.kwasniewski@colorado.edu

🏛️ Affiliations: University of Colorado Boulder, Institute for Arctic and Alpine Research (INSTAAR)

---
*Last updated: April 2026*
