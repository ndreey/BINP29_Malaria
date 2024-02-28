# BINP29 Case Study: The origin of human malaria [EDIT]
Herein, we go back a few years when the origin of _Plasmodium falciparum_ was still not clear. Answers to questions and our workflow/scripts will be added to this `README`.

## Setting the scene
Malaria is a disease caused by parasites within the phylum **Apicomplexa**, with the most prevalent and severe species affecting humans being  **_Plasmodium falciparum_**. However, it is an ongoing debate on the evolution of _P. falciparum_. The main point of contention is whether it is more closely related to other mammalian parasites or if it originated from a bird malaria parasite that changed its host. The genomes of many _Plasmodium_ species have been sequenced, including, notably, the first genome sequence of an avian malaria parasite, **_Haemoproteus tartakovskyi_**. An extensive phylogenetic analysis may provide answers to this controversy.

### Species and Genomic Information [EDIT]
meta sample info

| #  | Species                    | Host     | Genome Size | Genes | Genomic GC |
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


## Existing gene predictions
Other students have already run gene prediction on the _Plasmodium spp._ genomes and _Toxoplasma gondii_ using `gmes_petap.pl` from the GeneMark program (version unknown).

```
# Create directory for gene predictions and soft link them.
mkdir 00_GENEMARK
ln -s /tmp/Prediction/* 00_GENEMARK/
ln -s 
```

## Processing the novel _Haemoproteus tartakovskyi_ assembly
Both shotgun and paired end sequencing was used to sequence the _H. tartakovskyi_ genome (454 pyroseq). However, as avian blood cells have nuclei there is an abundance of host contamination. 