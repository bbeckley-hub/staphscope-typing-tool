#!/usr/bin/env python3
"""
StaphScope Comprehensive Lineage Database - UPDATED & CORRECTED
Complete MRSA/MSSA/LA-MRSA lineage reference for prediction engine
Author: Brown Beckley <brownbeckley94@gmail.com>
Scientific Review: Updated with literature-validated corrections
Date: 2025
Send a quick mail for any issues or further explanations.
Affiliation: University of Ghana Medical School-Department of Medical Bioichemistry
"""

# =============================================================================
# COMPREHENSIVE STAPHYLOCOCCUS AUREUS LINEAGE DATABASE 
# =============================================================================

LINEAGE_DATABASE = {
    # =========================================================================
    # HEALTHCARE-ASSOCIATED MRSA (HA-MRSA) LINEAGES
    # =========================================================================
    
    "CC5": {
        "primary_name": "New York/Japan Clone",
        "type": "Healthcare-associated MRSA",
        "subtypes": ["USA100", "New York/Japan", "Pediatric", "EMRSA-3"],
        "sequence_types": [5, 105, 125, 225, 228, 231, 235, 236, 240, 245, 685, 1034, 1057, 1089],
        "common_spa_types": ["t002", "t045", "t067", "t106", "t242", "t311", "t548", "t688"],
        "sccmec_types": ["II", "IV"],
        "geographic_distribution": {
            "regions": ["Global", "North America", "Europe", "Asia", "South America"],
            "prevalence": "High",
            "notes": "Dominant healthcare clone in many regions, SCCmec II predominant"
        },
        "virulence_profile": {
            "toxins": ["sea", "sek", "seq", "tst"],
            "adhesins": ["fnbA", "fnbB", "clfA", "clfB"],
            "immune_evasion": ["scn", "chp", "sak"],
            "enzymes": ["aur", "spl", "geh"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Erythromycin", "Clindamycin", "Aminoglycosides"],
            "common_genes": ["mecA", "erm(A)", "erm(C)", "aac(6')-aph(2'')"],
            "typical_patterns": ["Multidrug-resistant", "Often VISA/VRSA precursors"]
        },
        "clinical_significance": "Major healthcare-associated pandemic clone causing bacteremia, pneumonia, surgical site infections",
        "outbreak_potential": "HIGH",
        "risk_level": "HIGH",
        "key_references": ["PMID: 22527128", "PMID: 22617140", "PMID: 19680247", "PMID: 27992523"]
    },
    
    "CC6": {
        "primary_name": "UK/Irish Hospital Clone",
        "type": "Healthcare-associated MRSA",
        "subtypes": ["UK", "Irish", "Early Hospital", "EMRSA-4"],
        "sequence_types": [6, 7, 258, 259, 1221],
        "common_spa_types": ["t304", "t701", "t843", "t955"],
        "sccmec_types": ["IV"],
        "geographic_distribution": {
            "regions": ["UK", "Ireland", "Europe"],
            "prevalence": "Medium",
            "notes": "Early hospital clone in UK and Ireland"
        },
        "virulence_profile": {
            "toxins": ["sea", "sek", "seq"],
            "adhesins": ["fnbA", "clfA", "clfB"],
            "immune_evasion": ["scn", "chp"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Erythromycin"],
            "common_genes": ["mecA", "erm(C)"],
            "typical_patterns": ["Healthcare-associated resistance"]
        },
        "clinical_significance": "Early hospital-associated clone in UK and Ireland",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["PMID: 22189119", "PMID: 24395241", "PMID: 27992523"]
    },
    
    "CC7": {
        "primary_name": "Pediatric/Neonatal Clone",
        "type": "Healthcare-associated/Community-associated MRSA",
        "subtypes": ["Pediatric", "Neonatal", "USA Pediatric"],
        "sequence_types": [7, 6, 72, 254, 1223],
        "common_spa_types": ["t091", "t148", "t657", "t1234"],
        "sccmec_types": ["IV", "V"],
        "geographic_distribution": {
            "regions": ["North America", "Europe", "Global"],
            "prevalence": "Medium",
            "notes": "Associated with pediatric and neonatal infections"
        },
        "virulence_profile": {
            "toxins": ["lukS-PV", "lukF-PV", "sea"],
            "adhesins": ["fnbA", "clfA"],
            "immune_evasion": ["scn", "sak"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Often susceptible to non-beta-lactams", "Pediatric resistance patterns"]
        },
        "clinical_significance": "Pediatric and neonatal infections, both community and healthcare-associated",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["PMID: 23055880", "PMID: 23177801", "PMID: 19178545"]
    },
    
    "CC8": {
        "primary_name": "USA300/Pandemic CA-MRSA",
        "type": "Healthcare-associated/Community-associated MRSA", 
        "subtypes": ["USA300", "USA500", "EMRSA-2", "EMRSA-6", "Berlin", "South German"],
        "sequence_types": [8, 72, 241, 247, 250, 254, 260, 1139, 1223],
        "common_spa_types": ["t008", "t024", "t064", "t121", "t334", "t622"],
        "sccmec_types": ["IV", "IVa", "V"],
        "geographic_distribution": {
            "regions": ["Global", "North America", "Europe"],
            "prevalence": "Very High", 
            "notes": "USA300 is major community clone with global spread"
        },
        "virulence_profile": {
            "toxins": ["lukS-PV", "lukF-PV", "sea", "sek", "seq"],
            "adhesins": ["fnbA", "clfA", "sdrC"],
            "immune_evasion": ["scn", "chp", "sak"],
            "other": ["ACME", "PSM-mec"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Fluoroquinolones", "Erythromycin"],
            "common_genes": ["mecA", "norA", "erm(C)"],
            "typical_patterns": ["Community-associated resistance", "Fluoroquinolone resistance common"]
        },
        "clinical_significance": "Includes both HA-MRSA and CA-MRSA, USA300 is hypervirulent community clone",
        "outbreak_potential": "VERY HIGH",
        "risk_level": "HIGH",
        "key_references": ["PMID: 26484389", "PMID: 20610826", "PMID: 27992523"]
    },
    
    "CC12": {
        "primary_name": "European MSSA/MRSA Clone",
        "type": "Healthcare-associated MRSA",
        "subtypes": ["European", "MSSA-background", "UK"],
        "sequence_types": [12, 13, 218, 1226],
        "common_spa_types": ["t160", "t189", "t230", "t344", "t878"],
        "sccmec_types": ["IV"],
        "geographic_distribution": {
            "regions": ["Europe", "UK", "Global"],
            "prevalence": "Medium",
            "notes": "Common MSSA background that acquired SCCmec"
        },
        "virulence_profile": {
            "toxins": ["sea", "sek", "seq"],
            "adhesins": ["fnbA", "clfA", "sdrC"],
            "immune_evasion": ["scn", "chp", "sak"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Erythromycin"],
            "common_genes": ["mecA", "erm(C)"],
            "typical_patterns": ["Variable resistance patterns"]
        },
        "clinical_significance": "Healthcare-associated lineage from common MSSA background",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["PMID: 21521706", "PMID: 17956927"]
    },
    
    "CC15": {
        "primary_name": "European MRSA",
        "type": "Healthcare-associated",
        "subtypes": ["European"],
        "sequence_types": [15, 22, 508, 582, 583, 1047, 1136],
        "common_spa_types": ["t084", "t085", "t163", "t346", "t586", "t078"],
        "sccmec_types": ["IV"],
        "geographic_distribution": {
            "regions": ["Global"],
            "prevalence": "Medium",
            "notes": "Healthcare-associated lineage"
        },
        "virulence_profile": {
            "toxins": ["Variable"],
            "adhesins": ["fnbA", "clfA"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Variable"]
        },
        "clinical_significance": "Healthcare-associated lineage",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["PMID: 23675030", "PMID: 35056595", "PMID: 27992523"]
    },
    
    "CC22": {
        "primary_name": "EMRSA-15/Barnim Clone",
        "type": "Healthcare-associated MRSA",
        "subtypes": ["EMRSA-15", "Middle Eastern", "Gaza", "Barnim"],
        "sequence_types": [22, 217, 221, 224, 228, 230, 233, 1057, 1137],
        "common_spa_types": ["t022", "t032", "t223", "t309", "t852", "t899"],
        "sccmec_types": ["IV", "IVh"],
        "geographic_distribution": {
            "regions": ["Global", "UK", "Europe", "Middle East"],
            "prevalence": "High",
            "notes": "Major epidemic clone in UK and European hospitals"
        },
        "virulence_profile": {
            "toxins": ["eta", "etb", "tst"],
            "adhesins": ["fnbA", "clfA", "clfB"],
            "immune_evasion": ["scn", "chp"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Erythromycin", "Gentamicin", "Mupirocin"],
            "common_genes": ["mecA", "erm(C)", "aac(6')-aph(2'')", "mupA"],
            "typical_patterns": ["Multidrug-resistant", "Often gentamicin-resistant"]
        },
        "clinical_significance": "Major epidemic clone in UK and Europe, causes nosocomial infections",
        "outbreak_potential": "HIGH",
        "risk_level": "HIGH",
        "key_references": ["PMID: 35056595", "PMID: 11310446", "PMID: 23675030", "https://doi.org/10.1016/J.JIPH.2010.09.004", "DOI:10.4084/MJHID.2021.050"]
    },
    
    "CC25": {
        "primary_name": "European MRSA",
        "type": "Healthcare-associated",
        "subtypes": ["European"],
        "sequence_types": [25, 1172, 1173, 1131],
        "common_spa_types": ["t078", "t186", "t280", "t324"],
        "sccmec_types": ["IV"],
        "geographic_distribution": {
            "regions": ["Europe"],
            "prevalence": "Low",
            "notes": "Healthcare-associated lineage"
        },
        "virulence_profile": {
            "toxins": ["Variable"],
            "adhesins": ["fnbA", "clfA"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Variable"]
        },
        "clinical_significance": "Healthcare-associated lineage",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["", "PMID: 27605711", "PMID: 27992523"]
    },
    
    "CC30": {
        "primary_name": "EMRSA-16/Southwest Pacific Clone",
        "type": "Healthcare-associated/Community-associated MRSA",
        "subtypes": ["EMRSA-16", "Southwest Pacific", "Iberian", "Chilean"],
        "sequence_types": [30, 34, 36, 37, 38, 39, 834, 1092],
        "common_spa_types": ["t018", "t019", "t021", "t024", "t037", "t138", "t318", "t937"],
        "sccmec_types": ["II", "IV"],
        "geographic_distribution": {
            "regions": ["Global", "Europe", "South America", "Oceania"],
            "prevalence": "High",
            "notes": "Major healthcare clone with community transmission"
        },
        "virulence_profile": {
            "toxins": ["tst", "sea", "sec", "sel", "lukS-PV", "lukF-PV"],
            "adhesins": ["fnbA", "clfA"],
            "immune_evasion": ["scn", "sak"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Erythromycin", "Clindamycin", "Aminoglycosides"],
            "common_genes": ["mecA", "erm(A)", "erm(C)"],
            "typical_patterns": ["Variable resistance patterns"]
        },
        "clinical_significance": "Major healthcare clone, also community-associated in some regions",
        "outbreak_potential": "HIGH",
        "risk_level": "HIGH",
        "key_references": ["PMID: 22586109", "PMID: 23284024", "PMID: 27992523"]
    },
    
    "CC45": {
        "primary_name": "USA600/Berlin Clone",
        "type": "Healthcare-associated MRSA",
        "subtypes": ["USA600", "Berlin", "Hanover", "Canadian"],
        "sequence_types": [45, 47, 49, 50, 51, 52, 508, 1093, 1134],
        "common_spa_types": ["t015", "t026", "t038", "t040", "t067", "t071", "t186", "t231", "t269", "t742"],
        "sccmec_types": ["IV", "V"],
        "geographic_distribution": {
            "regions": ["Global", "Europe", "North America"],
            "prevalence": "Medium",
            "notes": "Healthcare-associated lineage with moderate spread"
        },
        "virulence_profile": {
            "toxins": ["sea", "sek", "seq"],
            "adhesins": ["fnbA", "clfA"],
            "enzymes": ["aur", "spl"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Erythromycin", "Tetracycline"],
            "common_genes": ["mecA", "erm(C)", "tet(K)"],
            "typical_patterns": ["Healthcare-associated resistance patterns"]
        },
        "clinical_significance": "Healthcare-associated lineage with moderate virulence",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["PMID: 27605711", "PMID: 23887918", "PMID: 32822004"]
    },
    
    "CC51": {
        "primary_name": "Asian Healthcare MRSA",
        "type": "Healthcare-associated",
        "subtypes": ["Asian"],
        "sequence_types": [51, 149, 150],
        "common_spa_types": ["t037", "t045", "t067"],
        "sccmec_types": ["III"],
        "geographic_distribution": {
            "regions": ["Global"],
            "prevalence": "Medium",
            "notes": "Healthcare-associated lineage"
        },
        "virulence_profile": {
            "toxins": ["Variable"],
            "adhesins": ["fnbA", "clfA"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Multi-drug resistant"]
        },
        "clinical_significance": "Healthcare-associated lineage",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["PMID: 19178545"]
    },
    
    "CC239": {
        "primary_name": "Brazilian/Hungarian Clone",
        "type": "Healthcare-associated MRSA",
        "subtypes": ["Brazilian", "Hungarian", "Portuguese"],
        "sequence_types": [239, 241, 247],
        "common_spa_types": ["t037", "t030", "t045"],
        "sccmec_types": ["III"],
        "geographic_distribution": {
            "regions": ["South America", "Europe", "Asia"],
            "prevalence": "High",
            "notes": "Major epidemic clone, recombinant lineage (CC8/CC30 hybrid)"
        },
        "virulence_profile": {
            "toxins": ["tst", "sea", "sec"],
            "adhesins": ["fnbA", "clfA"],
            "immune_evasion": ["scn", "chp"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Aminoglycosides", "Macrolides"],
            "common_genes": ["mecA", "aac(6')-aph(2'')", "erm(A)"],
            "typical_patterns": ["Multidrug-resistant"]
        },
        "clinical_significance": "Major epidemic healthcare clone, recombinant lineage",
        "outbreak_potential": "HIGH",
        "risk_level": "HIGH",
        "key_references": ["https://doi.org/10.1128/jcm.43.10.5069-5073.2005", "PMID: 38990431", "PMID: 36504833", "https://doi.org/10.1089/fpd.2013.171"]
    },
    
    # =========================================================================
    # COMMUNITY-ASSOCIATED MRSA (CA-MRSA) LINEAGES  
    # =========================================================================
    
    "CC1": {
        "primary_name": "USA400/Western Samoan Clone",
        "type": "Community-associated MRSA",
        "subtypes": ["USA400", "WA-MRSA", "Western Samoan", "MW2"],
        "sequence_types": [1, 122, 124, 128, 129, 1229, 1247],
        "common_spa_types": ["t127", "t128", "t136", "t159", "t174", "t178", "t227", "t321", "t384", "t573"],
        "sccmec_types": ["IV", "IVa", "V"],
        "geographic_distribution": {
            "regions": ["Global", "USA", "Australia", "Europe", "Asia"],
            "prevalence": "Medium",
            "notes": "Early CA-MRSA clone with pediatric association"
        },
        "virulence_profile": {
            "toxins": ["lukS-PV", "lukF-PV", "sek", "seq"],
            "adhesins": ["fnbA", "clfA"],
            "immune_evasion": ["scn", "chp"],
            "enzymes": ["aur", "spl"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Often susceptible to other antibiotics", "Community resistance pattern"]
        },
        "clinical_significance": "Early CA-MRSA clone, causes skin infections, necrotizing pneumonia",
        "outbreak_potential": "HIGH",
        "risk_level": "HIGH",
        "key_references": ["PMID: 34223815", "https://doi.org/10.1186/s40001-024-02076-z", "https://doi.org/10.3389/fmicb.2019.00139"]
    },
    
    "CC59": {
        "primary_name": "Taiwan/Asian-Pacific Clone",
        "type": "Community-associated MRSA", 
        "subtypes": ["Taiwan", "Asian-Pacific", "Korean"],
        "sequence_types": [59, 338, 339, 340, 341, 375, 954, 1123],
        "common_spa_types": ["t437", "t441", "t454", "t632", "t1635"],
        "sccmec_types": ["IV", "V"],
        "geographic_distribution": {
            "regions": ["Asia-Pacific", "Taiwan", "China", "Korea", "USA", "Australia"],
            "prevalence": "High",
            "notes": "Major CA-MRSA clone in Asia with increasing global spread"
        },
        "virulence_profile": {
            "toxins": ["lukS-PV", "lukF-PV", "eta", "etb"],
            "adhesins": ["fnbA", "clfA"],
            "enzymes": ["aur"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Erythromycin", "Clindamycin"],
            "common_genes": ["mecA", "erm(B)", "erm(C)"],
            "typical_patterns": ["Often erythromycin-resistant"]
        },
        "clinical_significance": "Major CA-MRSA clone in Asia, causes skin infections, necrotizing fasciitis",
        "outbreak_potential": "HIGH",
        "risk_level": "HIGH",
        "key_references": ["PMID: 32544568", "PMID: 20211891", "PMID: 23966409", "https://doi.org/10.1371/journal.pone.0070602", " DOI https://doi.org/10.2147/IDR.S284781 "]
    },
    
    "CC80": {
        "primary_name": "European CA-MRSA",
        "type": "Community-associated MRSA",
        "subtypes": ["European", "Middle Eastern", "Mediterranean"],
        "sequence_types": [80, 1154, 1239, 1250],
        "common_spa_types": ["t044", "t131", "t376", "t860"],
        "sccmec_types": ["IV"],
        "geographic_distribution": {
            "regions": ["Europe", "Middle East", "North Africa"],
            "prevalence": "Variable",
            "notes": "Major CA-MRSA clone in Europe and Middle East"
        },
        "virulence_profile": {
            "toxins": ["lukS-PV", "lukF-PV", "eta", "etb"],
            "adhesins": ["fnbA", "clfA"],
            "enzymes": ["aur"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Often susceptible to other antibiotics", "Often fusidic acid resistant"]
        },
        "clinical_significance": "Major CA-MRSA clone in Europe and Middle East",
        "outbreak_potential": "HIGH", 
        "risk_level": "HIGH",
        "key_references": ["PMID: 25078407", "https://doi.org/10.1016/j.cmi.2017.06.022", "https://doi.org/10.1128/jcm.01381-07"]
    },
    
    "CC88": {
        "primary_name": "African/Mediterranean CA-MRSA",
        "type": "Community-associated MRSA",
        "subtypes": ["African", "Middle Eastern", "Malaysian"],
        "sequence_types": [88, 78, 1002, 1187, 1249],
        "common_spa_types": ["t786", "t325", "t861", "t1324", "t1869"],
        "sccmec_types": ["IV", "V"],
        "geographic_distribution": {
            "regions": ["North Africa", "Middle East", "Southeast Asia", "Mediterranean"],
            "prevalence": "Medium",
            "notes": "Emerging CA-MRSA clone in specific regions"
        },
        "virulence_profile": {
            "toxins": ["lukS-PV", "lukF-PV"],
            "adhesins": ["fnbA", "clfA"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Variable resistance patterns"]
        },
        "clinical_significance": "Emerging CA-MRSA clone in developing regions",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["PMID: 22745670", "https://doi.org/10.3390/microorganisms12010017",  "PMID: 27992523"]
    },
    
    "CC93": {
        "primary_name": "Queensland/Australian Clone",
        "type": "Community-associated MRSA",
        "subtypes": ["Queensland", "Australian"],
        "sequence_types": [93, 1179, 1180],
        "common_spa_types": ["t202", "t334", "t573", "t1340"],
        "sccmec_types": ["IV", "V"],
        "geographic_distribution": {
            "regions": ["Australia", "New Zealand"],
            "prevalence": "Medium",
            "notes": "CA-MRSA clone primarily in Australia and New Zealand"
        },
        "virulence_profile": {
            "toxins": ["lukS-PV", "lukF-PV"],
            "adhesins": ["fnbA", "clfA"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Often susceptible to other antibiotics"]
        },
        "clinical_significance": "CA-MRSA clone in Australia and New Zealand",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": [" DOI: 10.1111/apm.12181", "https://doi.org/10.1111/apm.12181", "https://doi.org/10.1016/j.cmi.2016.11.002"]
    },
    
    "CC97": {
        "primary_name": "Bovine/Human Associated",
        "type": "Animal-associated/Human",
        "subtypes": ["Bovine", "European"],
        "sequence_types": [97, 71, 1265, 1303, 1394, 1715],
        "common_spa_types": ["t173", "t267", "t359", "t911", "t1344", "t1736", "t2247"],
        "sccmec_types": ["IV", "V"],
        "geographic_distribution": {
            "regions": ["Europe", "Global"],
            "prevalence": "Medium",
            "notes": "Originally bovine-associated, now human infections"
        },
        "virulence_profile": {
            "toxins": ["sea", "sek", "seq"],
            "adhesins": ["fnbA", "clfA"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Variable"]
        },
        "clinical_significance": "Originally bovine-associated, now human infections",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["https://doi.org/10.1016/j.vetmic.2022.109374", "https://doi.org/10.3390/foods9040439", "PMID: 26590279", "PMID: 27992523"]
    },
    
    "CC121": {
        "primary_name": "European CA-MRSA",
        "type": "Community-associated MRSA",
        "subtypes": ["European CA-MRSA"],
        "sequence_types": [121, 954, 1125, 1529, 1530],
        "common_spa_types": ["t159", "t314", "t342", "t525", "t645", "t852"],
        "sccmec_types": ["IV", "V"],
        "geographic_distribution": {
            "regions": ["Europe", "Global"],
            "prevalence": "Medium",
            "notes": "Emerging CA-MRSA clone"
        },
        "virulence_profile": {
            "toxins": ["lukS-PV", "lukF-PV"],
            "adhesins": ["fnbA", "clfA"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Often susceptible to other antibiotics"]
        },
        "clinical_significance": "Emerging CA-MRSA clone",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["https://doi.org/10.4269/ajtmh.16-0746"]
    },
    
    "CC772": {
        "primary_name": "Bengal Bay Clone", 
        "type": "Community-associated MRSA",
        "subtypes": ["Bengal Bay", "Indian Subcontinent"],
        "sequence_types": [772, 2884],
        "common_spa_types": ["t657", "t1081"],
        "sccmec_types": ["IV", "V"],
        "geographic_distribution": {
            "regions": ["Indian Subcontinent", "UK", "Middle East"],
            "prevalence": "Medium",
            "notes": "PVL-positive, often associated with travel"
        },
        "virulence_profile": {
            "toxins": ["lukS-PV", "lukF-PV", "eta", "etb"],
            "clinical": "Severe skin and soft tissue infections"
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Variable resistance patterns"]
        },
        "clinical_significance": "CA-MRSA clone from Indian subcontinent with global travel-associated spread",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["PMID: 33122734", "PMID: 29238933", " https://doi.org/10.3389/fmicb.2019.02505"]
    },
    
    "CC1153": {
        "primary_name": "Arabian Gulf/North African Clone",
        "type": "Community-associated MRSA", 
        "subtypes": ["Arabian Gulf", "North African", "fusC-positive"],
        "sequence_types": [1153, 1526, 1527],
        "common_spa_types": ["t657", "t1081", "t1869"],
        "sccmec_types": ["V", "chimeric"],
        "geographic_distribution": {
            "regions": ["Arabian Gulf", "North Africa", "Middle East"],
            "prevalence": "Emerging",
            "notes": "PVL-positive with fusidic acid resistance (fusC)"
        },
        "virulence_profile": {
            "toxins": ["lukS-PV", "lukF-PV"],
            "other": ["fusC positive"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Fusidic Acid"],
            "common_genes": ["mecA", "fusC"],
            "typical_patterns": ["Often fusidic acid resistant"]
        },
        "clinical_significance": "Emerging CA-MRSA clone in Arabian Gulf and North Africa",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM", 
        "key_references": ["PMID: 30951513", "PMID: 36671279"]
    },
    
    # =========================================================================
    # LIVESTOCK-ASSOCIATED MRSA (LA-MRSA) LINEAGES
    # =========================================================================
    
    "CC9": {
        "primary_name": "Asian LA-MRSA",
        "type": "Livestock-associated MRSA",
        "subtypes": ["Asian LA-MRSA", "Chinese"],
        "sequence_types": [9, 834, 835, 1187, 1347],
        "common_spa_types": ["t899", "t1430", "t1939", "t337"],
        "sccmec_types": ["III", "IV", "V"],
        "geographic_distribution": {
            "regions": ["Asia", "China", "Korea", "Increasing global"],
            "prevalence": "Medium",
            "notes": "Livestock-associated in Asia, emerging human infections"
        },
        "virulence_profile": {
            "toxins": ["Limited virulence factors"],
            "adhesins": ["Livestock adaptation"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecA"],
            "typical_patterns": ["Multi-drug resistant"]
        },
        "clinical_significance": "Livestock-associated in Asia, emerging human infections",
        "outbreak_potential": "MEDIUM",
        "risk_level": "MEDIUM",
        "key_references": ["PMID: 36363707", "https://doi.org/10.3389/fmicb.2012.00103"]
    },
    
    "CC101": {
        "primary_name": "CC101 (reported in animals and humans)",
        "type": "Occasional livestock-associated / human isolates",
        "subtypes": [],
        "sequence_types": ["ST101 reported in some studies; other STs possible"],
        "common_spa_types": ["spa types reported but variable; see PubMLST for details"],
        "sccmec_types": ["variable / not consistently reported"],
        "geographic_distribution": {
            "regions": ["Reported in livestock/food-animal studies (various countries)"],
            "prevalence": "Low",
            "notes": "Reported in animal surveillance studies; not a dominant LA-MRSA lineage (see refs)."
        },
        "virulence_profile": {
            "toxins": ["variable / not consistently reported"],
            "adhesins": ["variable"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["variable; MRSA and MSSA reports exist"],
            "common_genes": ["mecA may be present in MRSA isolates but not universal"],
            "typical_patterns": ["Low prevalence; patterns vary by study"]
        },
        "clinical_significance": "Occasional livestock-associated isolates; low prevalence in surveys",
        "outbreak_potential": "LOW",
        "risk_level": "LOW",
        "key_references": ["PMID: 31920996", "PMID: 19178545"]
    },
    
    "CC130": {
        "primary_name": "mecC-MRSA/Livestock-associated",
        "type": "Livestock-associated MRSA",
        "subtypes": ["mecC-MRSA", "Bovine", "Wildlife", "European"],
        "sequence_types": [130, 425, 1945, 1947, 1948],
        "common_spa_types": ["t843", "t1736", "t3256", "t3570"],
        "sccmec_types": ["XI"],
        "geographic_distribution": {
            "regions": ["Europe", "UK", "Global"],
            "prevalence": "Low but emerging",
            "notes": "mecC-positive MRSA, zoonotic transmission from livestock and wildlife"
        },
        "virulence_profile": {
            "toxins": ["Limited virulence factors"],
            "adhesins": ["Livestock adaptation"],
            "other": ["mecC gene", "Lower virulence in human models"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin"],
            "common_genes": ["mecC"],
            "typical_patterns": ["Often susceptible to other antibiotics", "False negative in standard mecA PCR"]
        },
        "clinical_significance": "mecC-positive MRSA with zoonotic transmission, lower human virulence",
        "outbreak_potential": "LOW",
        "risk_level": "LOW",
        "key_references": ["https://doi.org/10.3168/jds.2013-7378", "PMID: 33841383", "PMID: 27992523", " https://doi.org/10.1007/s00248-019-01328-4"]
    },
    
    "CC398": {
        "primary_name": "LA-MRSA/Pig-associated",
        "type": "Livestock-associated MRSA",
        "subtypes": ["Pig-associated", "Cattle-associated", "Poultry-associated"],
        "sequence_types": [398, 813, 929, 1232, 1234, 1580],
        "common_spa_types": ["t011", "t034", "t108", "t1250", "t1451", "t1456", "t1457", "t2383", "t571", "t899", "t1939"],
        "sccmec_types": ["IV", "V", "III"],
        "geographic_distribution": {
            "regions": ["Global", "Europe", "North America", "Asia"],
            "prevalence": "High in livestock",
            "notes": "Dominant livestock-associated clone, occupational risk for farmers/vets"
        },
        "virulence_profile": {
            "toxins": ["Limited human virulence factors"],
            "adhesins": ["Livestock adaptation"],
            "immune_evasion": ["Often lacks human immune evasion cluster"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin", "Tetracycline", "Lincosamides"],
            "common_genes": ["mecA", "tet(M)", "tet(K)", "lincosamide resistance"],
            "typical_patterns": ["Livestock-associated resistance patterns", "Often tetracycline-resistant"]
        },
        "clinical_significance": "Dominant livestock-associated clone, occupational infections, can cause human infections",
        "outbreak_potential": "MEDIUM in occupational settings",
        "risk_level": "MEDIUM",
        "key_references": ["PMID: 19701815", "PMID: 27605711", "PMID: 27992523","PMID: 32822004", "https://doi.org/10.1371/journal.pone.0010990"]
    },
    
    "CC425": {
        "primary_name": "mecC-MRSA/Wildlife Clone",
        "type": "Livestock-associated/Wildlife MRSA",
        "subtypes": ["mecC-MRSA", "Wildlife", "Hedgehog", "European"],
        "sequence_types": [425, 130, 1945, 1946],
        "common_spa_types": ["t843", "t2346", "t3256", "t3570"],
        "sccmec_types": ["XI"],
        "geographic_distribution": {
            "regions": ["Europe", "UK", "Scandinavia"],
            "prevalence": "Low but emerging",
            "notes": "mecC-MRSA associated with wildlife, particularly hedgehogs"
        },
        "virulence_profile": {
            "toxins": ["Limited virulence factors"],
            "adhesins": ["Wildlife adaptation"],
            "other": ["mecC gene", "Wildlife adaptation genes"]
        },
        "resistance_profile": {
            "antimicrobial_classes": ["Methicillin (mecC)"],
            "common_genes": ["mecC"],
            "typical_patterns": ["mecC-positive, often susceptible to other antibiotics"]
        },
        "clinical_significance": "mecC-MRSA associated with wildlife, potential zoonotic transmission",
        "outbreak_potential": "LOW",
        "risk_level": "LOW",
        "key_references": ["PMID: 39650146", "https://doi.org/10.1186/1751-0147-55-6", "https://doi.org/10.1038/s41586-021-04265-w"]
    }
}

# =============================================================================
# SPECIALIZED LINEAGE PROFILES - UPDATED
# =============================================================================

SPECIALIZED_LINEAGES = {
    "PVL_POSITIVE": {
        "CC8-USA300": {
            "st": 8,
            "sccmec": "IV",
            "spa": "t008",
            "virulence": ["lukS-PV", "lukF-PV", "sea", "sek", "seq"],
            "clinical": "Necrotizing pneumonia, skin infections, osteomyelitis",
            "risk": "VERY HIGH"
        },
        "CC1-USA400": {
            "st": 1, 
            "sccmec": "IV",
            "spa": "t128",
            "virulence": ["lukS-PV", "lukF-PV", "sek", "seq"],
            "clinical": "Necrotizing pneumonia, severe skin infections",
            "risk": "HIGH"
        },
        "CC80-European": {
            "st": 80,
            "sccmec": "IV", 
            "spa": "t044",
            "virulence": ["lukS-PV", "lukF-PV", "eta", "etb"],
            "clinical": "Severe skin infections, necrotizing fasciitis",
            "risk": "HIGH"
        },
        "CC59-Taiwan": {
            "st": 59,
            "sccmec": "V",
            "spa": "t437", 
            "virulence": ["lukS-PV", "lukF-PV", "eta", "etb"],
            "clinical": "Skin infections, necrotizing fasciitis",
            "risk": "HIGH"
        },
        "CC1153-Arabian": {
            "st": 1153,
            "sccmec": "V",
            "spa": "t657",
            "virulence": ["lukS-PV", "lukF-PV"],
            "clinical": "Skin and soft tissue infections",
            "risk": "MEDIUM"
        }
    },
    
    "TST_POSITIVE": {
        "CC30-EMRSA16": {
            "st": 30,
            "sccmec": "II",
            "spa": "t019",
            "virulence": ["tst", "sea", "sec"],
            "clinical": "Toxic shock syndrome, menstrual TSS",
            "risk": "HIGH"
        },
        "CC5-USA100": {
            "st": 5,
            "sccmec": "II", 
            "spa": "t002",
            "virulence": ["tst", "sea", "sek", "seq"],
            "clinical": "Toxic shock syndrome, healthcare-associated",
            "risk": "HIGH"
        },
        "CC239-Brazilian": {
            "st": 239,
            "sccmec": "III",
            "spa": "t037",
            "virulence": ["tst", "sea", "sec"],
            "clinical": "Toxic shock syndrome, healthcare-associated",
            "risk": "HIGH"
        }
    },
    
    "EXFOLIATIVE_TOXIN": {
        "CC121-ETA": {
            "st": 121,
            "sccmec": "IV",
            "virulence": ["eta"],
            "clinical": "Staphylococcal scalded skin syndrome",
            "risk": "MEDIUM"
        },
        "CC80-ETB": {
            "st": 80,
            "sccmec": "IV",
            "virulence": ["etb"], 
            "clinical": "Bullous impetigo, scalded skin syndrome",
            "risk": "MEDIUM"
        },
        "CC772-BengalBay": {
            "st": 772,
            "sccmec": "IV",
            "virulence": ["eta", "etb"],
            "clinical": "Bullous impetigo, scalded skin syndrome",
            "risk": "MEDIUM"
        }
    },
    
    "FUSIDIC_ACID_RESISTANT": {
        "CC1153-Arabian": {
            "st": 1153,
            "sccmec": "V", 
            "spa": "t657",
            "resistance": ["fusC"],
            "clinical": "Skin infections in Middle East/North Africa",
            "risk": "MEDIUM"
        },
        "CC80-European": {
            "st": 80,
            "sccmec": "IV",
            "spa": "t044", 
            "resistance": ["fusC"],
            "clinical": "Community skin infections in Europe",
            "risk": "MEDIUM"
        }
    },

    "MEC_C_POSITIVE": {
        "CC130-LA-MRSA": {
            "st": 130,
            "sccmec": "XI",
            "spa": "t843", 
            "virulence": ["mecC"],
            "clinical": "Zoonotic infections from livestock, often missed by routine tests",
            "risk": "MEDIUM"
        },
        "CC425-Wildlife": {
            "st": 425,
            "sccmec": "XI",
            "spa": "t843",
            "virulence": ["mecC"],
            "clinical": "Wildlife-associated, potential zoonotic transmission",
            "risk": "LOW"
        }
    },

    "PEDIATRIC_ASSOCIATED": {
        "CC7-Pediatric": {
            "st": 7,
            "sccmec": "IV",
            "spa": "t091",
            "virulence": ["lukS-PV", "lukF-PV"],
            "clinical": "Neonatal and pediatric infections",
            "risk": "MEDIUM"
        },
        "CC1-USA400": {
            "st": 1,
            "sccmec": "IV", 
            "spa": "t128",
            "virulence": ["lukS-PV", "lukF-PV"],
            "clinical": "Pediatric necrotizing pneumonia",
            "risk": "HIGH"
        }
    }
}

# =============================================================================
# GEOGRAPHIC DISTRIBUTION PATTERNS - UPDATED
# =============================================================================

GEOGRAPHIC_LINEAGES = {
    "NORTH_AMERICA": {
        "USA300": {"cc": "CC8", "prevalence": "High", "setting": "Community/Hospital"},
        "USA100": {"cc": "CC5", "prevalence": "High", "setting": "Hospital"},
        "USA400": {"cc": "CC1", "prevalence": "Medium", "setting": "Community"},
        "USA500": {"cc": "CC8", "prevalence": "Low", "setting": "Hospital"},
        "USA600": {"cc": "CC45", "prevalence": "Low", "setting": "Hospital"},
        "USA800": {"cc": "CC8", "prevalence": "Low", "setting": "Hospital"}
    },
    
    "EUROPE": {
        "EMRSA-15": {"cc": "CC22", "prevalence": "High", "setting": "Hospital"},
        "EMRSA-16": {"cc": "CC30", "prevalence": "High", "setting": "Hospital"},
        "Berlin": {"cc": "CC45", "prevalence": "Medium", "setting": "Hospital"},
        "Irish-1": {"cc": "CC8", "prevalence": "Medium", "setting": "Hospital"},
        "European CA-MRSA": {"cc": "CC80", "prevalence": "Medium", "setting": "Community"},
        "Barnim": {"cc": "CC22", "prevalence": "Medium", "setting": "Hospital"}
    },
    
    "ASIA": {
        "Taiwan": {"cc": "CC59", "prevalence": "High", "setting": "Community"},
        "Japanese": {"cc": "CC5", "prevalence": "High", "setting": "Hospital"},
        "Korean": {"cc": "CC5", "prevalence": "Medium", "setting": "Hospital"},
        "Chinese": {"cc": "CC9", "prevalence": "Medium", "setting": "Livestock/Human"},
        "Bengal Bay": {"cc": "CC772", "prevalence": "Medium", "setting": "Community"},
        "Vietnamese": {"cc": "CC59", "prevalence": "Medium", "setting": "Community"}
    },
    
    "SOUTH_AMERICA": {
        "Chilean": {"cc": "CC5", "prevalence": "High", "setting": "Hospital"},
        "Brazilian": {"cc": "CC239", "prevalence": "High", "setting": "Hospital"},
        "Argentinian": {"cc": "CC30", "prevalence": "Medium", "setting": "Hospital"}
    },
    
    "AFRICA": {
        "African CA-MRSA": {"cc": "CC88", "prevalence": "Medium", "setting": "Community"},
        "North African": {"cc": "CC80", "prevalence": "Medium", "setting": "Community"},
        "Arabian Gulf": {"cc": "CC1153", "prevalence": "Emerging", "setting": "Community"}
    },
    
    "OCEANIA": {
        "Southwest Pacific": {"cc": "CC30", "prevalence": "High", "setting": "Hospital/Community"},
        "Queensland": {"cc": "CC93", "prevalence": "Medium", "setting": "Community"},
        "Western Australian": {"cc": "CC1", "prevalence": "Medium", "setting": "Community"}
    },
    
    "MIDDLE_EAST": {
        "Arabian Gulf": {"cc": "CC1153", "prevalence": "High", "setting": "Community"},
        "Middle Eastern": {"cc": "CC80", "prevalence": "Medium", "setting": "Community"},
        "Gaza": {"cc": "CC22", "prevalence": "Medium", "setting": "Hospital"}
    },

    "SCANDINAVIA": {
        "mecC-MRSA": {"cc": "CC130", "prevalence": "Emerging", "setting": "Livestock/Wildlife"},
        "Baltic CA-MRSA": {"cc": "CC80", "prevalence": "Medium", "setting": "Community"}
    },

    "SOUTHEAST_ASIA": {
        "Taiwan Clone": {"cc": "CC59", "prevalence": "High", "setting": "Community"},
        "Vietnamese": {"cc": "CC59", "prevalence": "Medium", "setting": "Community"},
        "Malaysian": {"cc": "CC88", "prevalence": "Medium", "setting": "Community"}
    }
}

# =============================================================================
# LINEAGE PREDICTION RULES - ENHANCED
# =============================================================================

PREDICTION_RULES = {
    "high_confidence": {
        "criteria": ["exact_st_match", "spa_type_match", "sccmec_match", "virulence_match"],
        "min_matches": 3,
        "confidence": "HIGH"
    },
    "medium_confidence": {
        "criteria": ["exact_st_match", "spa_type_match", "sccmec_match", "virulence_match"],
        "min_matches": 2,
        "confidence": "MEDIUM" 
    },
    "low_confidence": {
        "criteria": ["st_in_cc_range", "geographic_consistent", "resistance_consistent"],
        "min_matches": 2,
        "confidence": "LOW"
    },
    "hybrid_lineages": {
        "CC239": {
            "parental_ccs": ["CC8", "CC30"],
            "characteristics": ["Recombinant", "Healthcare-associated", "Multidrug-resistant"],
            "confidence_note": "Hybrid lineage requires specific markers"
        }
    }
}

# =============================================================================
# SCIENTIFIC VALIDATION NOTES
# =============================================================================

SCIENTIFIC_NOTES = {
    "spa_types": "Partial list - spa types vary geographically and temporally. Refer to PubMLST for current distributions.",
    "sequence_types": "Core STs shown - each CC contains multiple related STs via MLST allele profiles.",
    "sccmec_types": "Common types shown - variants and novel types exist. Some lineages show SCCmec flexibility.",
    "geographic": "Distributions change over time due to global spread and travel.",
    "virulence": "Gene presence varies within lineages - not all strains carry all listed virulence factors.",
    "validation": "Database validated against recent literature and PubMLST data.",
    "last_updated": "2025"
}

# =============================================================================
# LINEAGE DATABASE CLASS
# =============================================================================

class LineageDatabase:
    """Comprehensive S. aureus lineage database for StaphScope - UPDATED"""
    
    def __init__(self):
        self.lineages = LINEAGE_DATABASE
        self.specialized = SPECIALIZED_LINEAGES
        self.geographic = GEOGRAPHIC_LINEAGES
        self.rules = PREDICTION_RULES
        self.notes = SCIENTIFIC_NOTES
    
    def get_all_lineages(self) -> dict:
        """Return complete lineage database"""
        return {
            "lineages": self.lineages,
            "specialized_lineages": self.specialized,
            "geographic_distribution": self.geographic,
            "prediction_rules": self.rules,
            "scientific_notes": self.notes
        }
    
    def get_lineage_by_cc(self, clonal_complex: str) -> dict:
        """Get lineage data by clonal complex"""
        return self.lineages.get(clonal_complex, {})
    
    def get_lineages_by_type(self, lineage_type: str) -> dict:
        """Get all lineages of a specific type (HA-MRSA, CA-MRSA, LA-MRSA)"""
        return {cc: data for cc, data in self.lineages.items() if data.get("type") == lineage_type}
    
    def predict_lineage_from_st(self, sequence_type: int) -> tuple:
        """Predict clonal complex from sequence type"""
        for cc, data in self.lineages.items():
            if sequence_type in data["sequence_types"]:
                return cc, data["primary_name"]
        return "Unknown", "Unknown lineage"
    
    def predict_lineage_from_spa(self, spa_type: str) -> list:
        """Predict possible lineages from spa type"""
        matches = []
        for cc, data in self.lineages.items():
            if spa_type in data["common_spa_types"]:
                matches.append((cc, data["primary_name"]))
        return matches if matches else [("Unknown", "Unknown lineage")]
    
    def get_specialized_profiles(self, profile_type: str) -> dict:
        """Get specialized lineage profiles (PVL_POSITIVE, TST_POSITIVE, etc.)"""
        return self.specialized.get(profile_type, {})
    
    def get_geographic_distribution(self, region: str) -> dict:
        """Get lineage distribution for a specific geographic region"""
        return self.geographic.get(region, {})
    
    def save_database(self, output_path: str):
        """Save database to JSON file"""
        import json
        with open(output_path, 'w') as f:
            json.dump(self.get_all_lineages(), f, indent=2)
    
    def load_database(self, input_path: str):
        """Load database from JSON file"""
        import json
        with open(input_path, 'r') as f:
            data = json.load(f)
            self.lineages = data.get("lineages", {})
            self.specialized = data.get("specialized_lineages", {})
            self.geographic = data.get("geographic_distribution", {})
            self.rules = data.get("prediction_rules", {})
            self.notes = data.get("scientific_notes", {})

# =============================================================================
# VALIDATION AND UTILITY FUNCTIONS
# =============================================================================

def validate_lineage_database():
    """
    Validate the integrity of the lineage database
    """
    issues = []
    
    for cc, data in LINEAGE_DATABASE.items():
        # Check required fields
        required_fields = ["primary_name", "type", "sequence_types", "common_spa_types", 
                          "sccmec_types", "clinical_significance"]
        for field in required_fields:
            if field not in data:
                issues.append(f"Missing field '{field}' in {cc}")
        
        # Check geographic distribution structure
        if "geographic_distribution" not in data:
            issues.append(f"Missing geographic_distribution in {cc}")
        else:
            geo_fields = ["regions", "prevalence", "notes"]
            for field in geo_fields:
                if field not in data["geographic_distribution"]:
                    issues.append(f"Missing geographic_distribution.{field} in {cc}")
    
    return issues

def export_lineage_summary():
    """
    Export a summary of all lineages for quick reference
    """
    summary = {}
    for cc, data in LINEAGE_DATABASE.items():
        summary[cc] = {
            "primary_name": data["primary_name"],
            "type": data["type"],
            "risk_level": data.get("risk_level", "UNKNOWN"),
            "outbreak_potential": data.get("outbreak_potential", "UNKNOWN")
        }
    return summary

# =============================================================================
# MAIN EXECUTION AND TESTING
# =============================================================================

if __name__ == "__main__":
    # Validate database integrity
    print("=== VALIDATING LINEAGE DATABASE ===")
    issues = validate_lineage_database()
    if issues:
        print(f"Found {len(issues)} issues:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("✓ Database validation passed")
    
    # Create database instance
    db = LineageDatabase()
    
    # Test lineage prediction
    print("\n=== TESTING LINEAGE PREDICTION ===")
    test_cases = [
        {"st": 8, "spa": "t008"},
        {"st": 5, "spa": "t002"}, 
        {"st": 30, "spa": "t021"},
        {"st": 398, "spa": "t011"}
    ]
    
    for test in test_cases:
        cc, name = db.predict_lineage_from_st(test["st"])
        spa_matches = db.predict_lineage_from_spa(test["spa"])
        print(f"\nST{test['st']} + spa{test['spa']}:")
        print(f"  ST Prediction: {cc} ({name})")
        print(f"  SPA Matches: {spa_matches}")
    
    # Export summary
    print(f"\n=== DATABASE SUMMARY ===")
    summary = export_lineage_summary()
    print(f"Total lineages: {len(summary)}")
    
    # Count by type
    type_counts = {}
    for cc, data in summary.items():
        lineage_type = data["type"]
        type_counts[lineage_type] = type_counts.get(lineage_type, 0) + 1
    
    print("\nLineages by type:")
    for lineage_type, count in type_counts.items():
        print(f"  {lineage_type}: {count}")
    
    # Save database
    db.save_database("staphscope_lineage_database.json")
    print(f"\n💾 Database saved to: staphscope_lineage_database.json")
    print("✅ StaphScope Lineage Database UPDATED & CORRECTED!")
