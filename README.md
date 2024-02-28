# BINP29 Case Study: The origin of human malaria [EDIT]
Herein, we go back a few years when the origin of _Plasmodium falciparum_ was still not clear. Answers to questions and our workflow/scripts will be added to this `README`.

## Programs used for this analysis
All programs installed succesfully through `mamba` and their versions
```
#Name           Version         Channel
bbmap           39.06           bioconda
seqkit          2.7.0           bioconda

```

## Setting the scene
Malaria is a disease caused by parasites within the phylum **Apicomplexa**, with the most prevalent and severe species affecting humans being  **_Plasmodium falciparum_**. However, it is an ongoing debate on the evolution of _P. falciparum_. The main point of contention is whether it is more closely related to other mammalian parasites or if it originated from a bird malaria parasite that changed its host. The genomes of many _Plasmodium_ species have been sequenced, including, notably, the first genome sequence of an avian malaria parasite, **_Haemoproteus tartakovskyi_** (Ht). An extensive phylogenetic analysis may provide answers to this controversy.

### Species and Genomic Information [EDIT]
meta sample info

| Abbr.  | Species                | Host     | Genome Size | Genes | Genomic GC |
|----|----------------------------|----------|-------------|-------|------------|
| Pb | Plasmodium berghei         | rodents  |             |       |            |
| Pc | Plasmodium cynomolgi       | macaques |             |       |            |
| Pf | Plasmodium falciparum      | humans   |             |       |            |
| Pk | Plasmodium knowlesi        | lemurs   |             |       |            |
| Pv | Plasmodium vivax           | humans   |             |       |            |
| Py | Plasmodium yoelii          | rodents  |             |       |            |
| Ht | Haemoproteus tartakovskyi  | birds    |             |       |            |
| Tg | Toxoplasma gondii          | humans   |             |       |            |


### Do you think that in a phylogenetic tree the parasites that use similar hosts will group together?
Typically, yes, but it's not always straightforward. Especially when considering humans as hosts, it's unlikely they'll cluster neatly. Humans originated in Africa and have since spread globally. Given the potential for host-switching and the diversity of species radiations, it's not a given that Apicomplexa species would necessarily cluster based on their use of similar hosts.
Birds are not geographically set, whereas rodents and primates are mote constrained to their biogeographical setting. Knowing that malaria occurs in Africa and Asia, and that different primates exist there. We might see some interesting results.


## Setting up work directory
### Raw assembly of _H. tartakovskyi_
```
mkdir 00_RAW
ln -s /resources/binp29/Data/malaria/Haemoproteus_tartakovskyi.raw.genome.gz 00_RAW/
```
### Database information
```
mkdir DB
ln -s /resources/binp29/Data/malaria/SwissProt.fasta DB
ln -s /resources/binp29/Data/malaria/taxonomy.dat DB
ln -s /resources/binp29/Data/malaria/uniprot_sprot.dat DB
```
### Existing gene predictions
Other students have already run gene prediction on the _Plasmodium spp._ genomes and _Toxoplasma gondii_ using `gmes_petap.pl` from the GeneMark program (version unknown).

```
# Create directory for gene predictions and soft link them.
mkdir 01_GENEMARK
ln -s /tmp/Prediction/* 01_GENEMARK/
ln -s /resources/binp29/Data/malaria/Tg.gff.gz 01_GENEMARK/
```
### Basic structure
```
.
├── 00_RAW
│   └── Haemoproteus_tartakovskyi.raw.genome.gz -> /resources/binp29/Data/malaria/Haemoproteus_tartakovskyi.raw.genome.gz
├── 01_GENEMARK
│   ├── fixed_Pk.gtf -> /tmp/Prediction/fixed_Pk.gtf
│   ├── genemark.Pb.gtf -> /tmp/Prediction/genemark.Pb.gtf
│   ├── genemark.Pc.gtf -> /tmp/Prediction/genemark.Pc.gtf
│   ├── genemark.Pf.gtf -> /tmp/Prediction/genemark.Pf.gtf
│   ├── genemark.Pk.gtf -> /tmp/Prediction/genemark.Pk.gtf
│   ├── genemark.Pv.gtf -> /tmp/Prediction/genemark.Pv.gtf
│   ├── genemark.Py.gtf -> /tmp/Prediction/genemark.Py.gtf
│   └── Tg.gff.gz -> /resources/binp29/Data/malaria/Tg.gff.gz
├── DB
│   ├── SwissProt.fasta -> /resources/binp29/Data/malaria/SwissProt.fasta
│   ├── taxonomy.dat -> /resources/binp29/Data/malaria/taxonomy.dat
│   └── uniprot_sprot.dat -> /resources/binp29/Data/malaria/uniprot_sprot.dat
└── README.md
```
## Host decontamination of the novel _Haemoproteus tartakovskyi_ assembly.
Both shotgun and paired end sequencing was used to sequence the _H. tartakovskyi_ genome (454 pyroseq). However, as avian blood cells have nuclei there is an abundance of host contamination that has to be delt with. To decontaminate we will use **GC-content** and **contig length** thresholds as the parasite and host GC content is expected to differ. Any contig less than 3000 bp will be removed.

Summary statistics of the **raw_assembly** using `stats.sh`from `bbmap`. We see that we have at least 2750 contigs after hard filtering on 3000 bp length.
```
# Summary statistics
stats.sh 00_RAW/Haemoproteus_tartakovskyi.raw.genome.gz


A       C       G       T       N       IUPAC   Other   GC      GC_stdev
0.3482  0.1407  0.1333  0.3778  0.0000  0.0000  0.0000  0.2740  0.1071

Main genome scaffold total:             15048
Main genome contig total:               15048
Main genome scaffold sequence total:    27.427 MB
Main genome contig sequence total:      27.427 MB       0.000% gap
Main genome scaffold N/L50:             1371/5.219 KB
Main genome contig N/L50:               1371/5.219 KB
Main genome scaffold N/L90:             7572/721
Main genome contig N/L90:               7572/721
Max scaffold length:                    64.494 KB
Max contig length:                      64.494 KB
Number of scaffolds > 50 KB:            1
% main genome in scaffolds > 50 KB:     0.24%


Minimum         Number          Number          Total           Total           Scaffold
Scaffold        of              of              Scaffold        Contig          Contig  
Length          Scaffolds       Contigs         Length          Length          Coverage
--------        --------------  --------------  --------------  --------------  --------
    All                 15,048          15,048      27,426,784      27,426,710   100.00%
     50                 15,048          15,048      27,426,784      27,426,710   100.00%
    100                 15,048          15,048      27,426,784      27,426,710   100.00%
    250                 12,367          12,367      26,974,066      26,973,993   100.00%
    500                  9,767           9,767      26,034,309      26,034,249   100.00%
   1 KB                  5,463           5,463      22,904,063      22,904,026   100.00%
 2.5 KB                  2,750           2,750      18,741,094      18,741,069   100.00%
   5 KB                  1,442           1,442      14,078,731      14,078,713   100.00%
  10 KB                    466             466       7,209,139       7,209,130   100.00%
  25 KB                     33              33       1,033,424       1,033,424   100.00%
  50 KB                      1               1          64,494          64,494   100.00%
```

### GC-distribution
Instead of re-inventing the wheel we will use [seqkit](https://github.com/shenwei356/seqkit) which is a _"ultrafast toolkit for FASTA/Q file manipulation"_. With this tool we can get summary statistics for each contig as well as to apply our filters.
```
mkdir 02_DE-HOST

# Generating a .tsv file with GC% and length for each contig.
seqkit fx2tab --only-id --length --gc --threads 24 --out-file 02_DE-HOST/raw_contig_stats.tsv 00_RAW/Haemoproteus_tartakovskyi.raw.genome.gz

# The contig sequence was included as well, so lets remove the second column
awk '{print $1 "\t" $3 "\t" $4}' 02_DE-HOST/raw_contig_stats.tsv > temp.tsv
mv temp.tsv 02_DE-HOST/raw_contig_stats.tsv

head 02_DE-HOST/raw_contig_stats.tsv 
contig00001     64494   27.77
contig00002     42559   25.90
contig00003     40672   26.60
contig00004     35838   28.12
```
Lets now plot it in `R` using `tidyverse` library.
```
library(tidyverse)

# Read in the data
raw <- read.table("02_DE-HOST/raw_contig_stats.tsv", sep = "\t", header = F)

# Add column names
colnames(raw) <- c("contig", "length", "gc")

# Set custome theme
fig <- theme(
  panel.grid.major = element_blank(),  
  panel.grid.minor = element_blank(),  
  axis.ticks = element_line(size = 1, color = "black"),  
  axis.text.x = element_text(face = "bold", color = "black"),  
  axis.text.y = element_text(face = "bold", color = "black"),  
  axis.title.x = element_text(face = "bold", size = 12),  
  axis.title.y = element_text(face = "bold", size = 12),  
  panel.background = element_blank(),  
  panel.border = element_rect(colour = "black", fill=NA), 
  legend.position = "bottom",
  legend.title = element_blank(),  
  legend.background = element_blank(),  
  legend.key = element_blank(),  
  legend.text = element_text(face = "bold")    
)

# Plot density plot of GC distribution using ggplot
gc_dist <- raw %>% 
  ggplot(aes(x = gc)) +
  geom_density(fill = "steelblue1", colour = "black", alpha = 0.5) +
  labs(title = "Density plot of GC distribution",
       x = "GC content",
       y = "Density") +
  geom_vline(xintercept = 31.5, linetype = "dashed", color = "black") +
  annotate("text", x = 10, y = 0.04, label = "H. tartakovskyi", color = "black", 
           size = 4) +
  annotate("text", x = 45, y = 0.03, label = "Avian host", color = "black", 
           size = 4) +
  theme_bw() +
  fig

# Save the plot as a PNG file
ggsave("figures/gc_dist.png", plot = gc_dist, width = 8, height = 6, units = "in", 
       dpi = 300)
```
![GC% Distribution](figures/gc_dist.png)

We expect **Ht** to have a lesse GC% than the avian host. Hence, we set the GC threshold as 30 and keep everything that is below.

### Filtering
We had 15048 contigs originally and after filtering we have 2222 contigs with an GC-content of 25.68 % (compared to 27.40 % we had before).
```
# Get the contig ids that meet thresholds.
awk '$2 >= 3000 && $3 < 30 {print $1}' 02_DE-HOST/raw_contig_stats.tsv > 02_DE-HOST/keep_ids.txt

# Generate new fasta using said contig ids.
seqkit grep --threads 24 -f keep_ids.txt 00_RAW/Haemoproteus_tartakovskyi.raw.genome.gz > 02_DE-HOST/contigs-deHOST.fasta
```
Here we see the summary statistics of our decontaminated assembly.
```
> stats.sh 02_DE-HOST/contigs-deHOST.fasta 
 
A       C       G       T       N       IUPAC   Other   GC      GC_stdev
0.3569  0.1321  0.1247  0.3864  0.0000  0.0000  0.0000  0.2568  0.0219

Main genome scaffold total:             2222
Main genome contig total:               2222
Main genome scaffold sequence total:    16.825 MB
Main genome contig sequence total:      16.825 MB       0.000% gap
Main genome scaffold N/L50:             604/8.676 KB
Main genome contig N/L50:               604/8.676 KB
Main genome scaffold N/L90:             1740/3.989 KB
Main genome contig N/L90:               1740/3.989 KB
Max scaffold length:                    64.494 KB
Max contig length:                      64.494 KB
Number of scaffolds > 50 KB:            1
% main genome in scaffolds > 50 KB:     0.38%


Minimum         Number          Number          Total           Total           Scaffold
Scaffold        of              of              Scaffold        Contig          Contig  
Length          Scaffolds       Contigs         Length          Length          Coverage
--------        --------------  --------------  --------------  --------------  --------
    All                  2,222           2,222      16,824,648      16,824,630   100.00%
 2.5 KB                  2,222           2,222      16,824,648      16,824,630   100.00%
   5 KB                  1,374           1,374      13,501,854      13,501,838   100.00%
  10 KB                    451             451       6,996,687       6,996,680   100.00%
  25 KB                     33              33       1,033,424       1,033,424   100.00%
  50 KB                      1               1          64,494          64,494   100.00%
```

