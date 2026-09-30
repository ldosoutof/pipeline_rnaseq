# Intégration CIBERSORTx (déconvolution cellulaire, LM22)

Estime les proportions de 22 types cellulaires immunitaires (signature LM22) par
échantillon, à partir de la matrice TPM du run. Exécuté via **Apptainer** (rootless).

## 1. Prérequis (installation unique)

### a) Apptainer (à faire installer par l'admin — nécessite root)
```bash
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update && sudo apt install -y apptainer
which apptainer      # vérifier
```

### b) Compte + token CIBERSORTx (sur cibersortx.stanford.edu)
1. Créer un compte sur https://cibersortx.stanford.edu
2. Menu **Download** : récupérer le **token**, et télécharger le fichier de
   signature **LM22.txt**.
3. Stocker le token dans un fichier hors dépôt, en lecture restreinte :
   ```bash
   mkdir -p ~/secrets
   echo -n 'VOTRE_TOKEN' > ~/secrets/cibersortx_token
   chmod 600 ~/secrets/cibersortx_token
   ```

### c) Construire l'image Apptainer depuis l'image Docker officielle
```bash
apptainer build cibersortx_fractions.sif docker://cibersortx/fractions
# placer le .sif dans un emplacement stable, ex :
#   /appli/containers/cibersortx_fractions.sif
```

### d) Placer LM22.txt dans un emplacement stable
```bash
# ex : /dataref/cibersortx/LM22.txt
```

## 2. Configuration (config.yml, section cibersortx)

```yaml
cibersortx:
  enabled:    true
  sif:        "/appli/containers/cibersortx_fractions.sif"
  lm22:       "/dataref/cibersortx/LM22.txt"
  email:      "votre.email@chu-nantes.fr"
  token_file: "/home/ldosoutoferreira/secrets/cibersortx_token"
  perm:       100
  qn:         "FALSE"    # RNA-seq : FALSE (TRUE réservé au microarray)
```

## 3. Activation dans le pipeline

- Copier `scripts/ensg_to_symbol_cibersortx.py` et `rules/07_cibersortx.smk`
  dans le pipeline.
- Ajouter l'include de la règle dans le Snakefile principal (`pipeline.smk`) :
  `include: "../rules/07_cibersortx.smk"`
- Ajouter la sortie aux cibles finales (`rule all`) **si** `cibersortx.enabled` :
  `rules.cibersortx_fractions.output.results`

## 4. Sortie

`pipeline_v0/cibersortx/CIBERSORTx_Results.txt` : fractions cellulaires par
échantillon (une ligne par échantillon, une colonne par type cellulaire LM22,
+ P-value, Correlation, RMSE).

## Points de vigilance

- **Conversion ENSG -> symbole** : la matrice TPM du pipeline est en identifiants
  Ensembl ; LM22 est en symboles HGNC. Le script `ensg_to_symbol_cibersortx.py`
  convertit via le GTF du pipeline et **somme** les TPM des ENSG partageant un
  même symbole (standard, sans perte d'information).
- **Mode absolu non supporté** en conteneur (`--absolute TRUE` renvoie des zéros).
  On reste en mode relatif (fractions).
- **QN FALSE** pour du RNA-seq (la quantile normalization est destinée au microarray).
- **Pertinence biologique (diagnostic)** : LM22 est conçu pour tissu
  sanguin/immunitaire. Valider qu'elle correspond à vos échantillons avant tout
  usage diagnostique.
- **Dépendance externe** : le token appelle le serveur Stanford à chaque exécution.
  À évaluer pour la reproductibilité en contexte ISO 15189.

## Références
- Newman et al., *Nature Biotechnology* 2019 (CIBERSORTx).
- Documentation : https://cibersortx.stanford.edu
