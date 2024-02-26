source('scripts/load_libs.R')

#' yshah Sunday, Apr 10, 2022 04:13:09 PM
#' Setup phenotypes
eid = fread('scripts/temp/ukb26867/eid.csv', colClasses = 'character')

birth_year = fread('scripts/temp/ukb26867/birth_year.csv')
colnames(birth_year) = 'birth_year'
birth_year$birth_year = as.numeric(birth_year$birth_year)
bmi = fread('scripts/temp/ukb26867/bmi.csv', nThread = 20)
bmi = bmi %>%
    mutate(bmi = apply(bmi, 1, sum, na.rm = TRUE)) %>%
    dplyr::select(bmi)
pca = fread('scripts/temp/ukb26867/pca.csv', nThread = 20)
colnames(pca) = paste0('PC', 1:ncol(pca))
sex = fread('scripts/temp/ukb26867/sex.csv')
colnames(sex) = 'sex'
sex$sex = factor(ifelse(sex$sex == 1, 'Male', 'Female'), levels = c('Female', 'Male'))
smoking = fread('scripts/temp/ukb26867/smoking_status.csv', nThread = 20)
smoking = smoking[,1]
colnames(smoking) = 'smoking_status'
smoking = smoking %>%
    mutate(smoking_status = case_when(
               smoking_status == 0 ~ 'Never',
               smoking_status == 1 ~ 'Previous',
               smoking_status == 2 ~ 'Current',
               smoking_status == -3 ~ 'No Response',
               is.na(smoking_status) ~ 'Missing'
           ))
alcohol = fread('scripts/temp/ukb26867/alcohol_status.csv', nThread = 20)
alcohol = alcohol[,1]
colnames(alcohol) = 'alcohol_status'
alcohol = alcohol %>%
    mutate(alcohol_status = case_when(
               alcohol_status == 0 ~ 'Never',
               alcohol_status == 1 ~ 'Previous',
               alcohol_status == 2 ~ 'Current',
               alcohol_status == -3 ~ 'No Response',
               is.na(alcohol_status) ~ 'Missing'
           ))
father.history = fread('scripts/temp/ukb26867/father_history.csv', nThread = 20)
father.history = mclapply(1:nrow(father.history), function(rw){
    dt = father.history[rw, ]
    ct.na = rowSums(is.na(dt))
    if (ct.na == ncol(father.history)) return(FALSE)
    ct.pca = rowSums(dt == 13, na.rm = TRUE)
    if (ct.pca > 0) return(TRUE)
    return(FALSE)
}, mc.cores = 30) %>% unlist
sibling.history = fread('scripts/temp/ukb26867/sibling_history.csv', nThread = 20)
sibling.history = mclapply(1:nrow(sibling.history), function(rw){
    dt = sibling.history[rw, ]
    ct.na = rowSums(is.na(dt))
    if (ct.na == ncol(sibling.history)) return(FALSE)
    ct.pca = rowSums(dt == 13, na.rm = TRUE)
    if (ct.pca > 0) return(TRUE)
    return(FALSE)
}, mc.cores = 30) %>% unlist
primary.death = fread('scripts/temp/ukb26867/death_cause_primary_icd10.csv', nThread = 20)
primary.death = mclapply(1:nrow(primary.death), function(rw){
    dt = primary.death[rw, ]
    ct.na = rowSums(is.na(dt))
    if (ct.na == ncol(primary.death)) return(FALSE)
    ct.pca = rowSums(dt == 'C61', na.rm = TRUE)
    if (ct.pca > 0) return(TRUE)
    return(FALSE)
}, mc.cores = 30) %>% unlist
secondary.death = fread('scripts/temp/ukb26867/death_cause_secondary_icd10.csv', nThread = 20)
secondary.death = mclapply(1:nrow(secondary.death), function(rw){
    dt = secondary.death[rw, ]
    ct.na = rowSums(is.na(dt))
    if (ct.na == ncol(secondary.death)) return(FALSE)
    ct.pca = rowSums(dt == 'C61', na.rm = TRUE)
    if (ct.pca > 0) return(TRUE)
    return(FALSE)
}, mc.cores = 30) %>% unlist
age_at_assessment = fread('scripts/temp/ukb26867/age_at_assessment.csv')[,1]
colnames(age_at_assessment) = 'age_at_assessment'


## Combine fields for phenotype table
phenotypes = cbind(eid, sex, birth_year, bmi, smoking, alcohol, pca, father_history = father.history,
                   sibling_history = sibling.history, primary_pca_death = primary.death,
                   secondary_pca_death = secondary.death, age_at_assessment = age_at_assessment$age_at_assessment) %>%
    mutate(family_history = father_history | sibling_history, death_pca = primary_pca_death | secondary_pca_death)

## Ancestry
afr = fread('~/data/UKBiobank/custom_qc/fam_files/african_fam')$V1
eur = fread('~/data/UKBiobank/custom_qc/fam_files/euro_fam')$V1
eas = fread('~/data/UKBiobank/custom_qc/fam_files/asian_fam')$V1
phenotypes = phenotypes %>%
    mutate(ancestry = case_when(
               eid %in% afr ~ 'AFR',
               eid %in% eur ~ 'EUR',
               eid %in% eas ~ 'EAS',
               !eid %in% c(afr, eur, eas) ~ 'Other'
           ))

## Prostate cancer status
self.rep = fread('scripts/temp/ukb26867/cancer_self_report.csv')
self.rep = mclapply(1:nrow(self.rep), function(rw){
    dt = self.rep[rw, ]
    ct.na = rowSums(is.na(dt))
    if (ct.na == ncol(self.rep)) return(FALSE)
    ct.pca = rowSums(dt == 1044, na.rm = TRUE)
    if (ct.pca > 0) return(TRUE)
    return(FALSE)
}, mc.cores = 30) %>% unlist()

hesin = fread('~/data/UKBiobank/hesin_diagnosis/hesin_diag.txt', nThread = 30)
hesin.pca = hesin %>% filter(diag_icd10 == 'C61' | diag_icd9 == 185 | diag_icd9_nb == 185) %>% pull(eid) %>% unique
length(hesin.pca)
phenotypes$type = 'Control'
phenotypes$type[self.rep] = 'Case'
phenotypes$type[phenotypes$eid %in% hesin.pca] = 'Case'

## Other cancer status
hesin.other.cancer = hesin %>%
    filter(startsWith(diag_icd10, 'C1') | startsWith(diag_icd10, 'C0') | startsWith(diag_icd10, 'C2') |
           startsWith(diag_icd10, 'C3') | startsWith(diag_icd10, 'C4') | startsWith(diag_icd10, 'C5') |
           startsWith(diag_icd10, 'C6') | startsWith(diag_icd10, 'C7') | startsWith(diag_icd10, 'C8') |
           startsWith(diag_icd10, 'C9')) %>% filter(diag_icd10 != 'C61') %>%
    pull(eid) %>% unique
phenotypes$other_cancer = FALSE
phenotypes$other_cancer[phenotypes$eid %in% hesin.other.cancer] = TRUE

## Remove bad samples
bad.het = read.table('~/data/UKBiobank/custom_qc/kick/hetrozy.txt')$V1
bad.sex = read.table('~/data/UKBiobank/custom_qc/kick/sexchromo.txt')$V1
bad.rel = read.table('~/data/UKBiobank/custom_qc/kick/excessrel.txt')$V1
bad.samp = c(bad.het, bad.sex, bad.rel) %>% unique
phenotypes = phenotypes %>% filter(!eid %in% bad.samp)

## Keep samples present in scott's eid list
genet.present = fread('~/data/UKBiobank/calls/ukbb.10.fam')
phenotypes = phenotypes %>% filter(eid %in% genet.present$V2)

## Annotate whether participant is alive
alive = fread('~/data/UKBiobank/custom_qc/fam_files/qc_alive_fam')
phenotypes = phenotypes %>% mutate(alive = ifelse(eid %in% alive$V2, TRUE, FALSE))

## Save raw phenotypes
setkey(phenotypes, 'eid')
saveRDS(phenotypes, 'db/phenotypes/phenotypes.rds')







############################################################################################################
############################################################################################################
############                                   Test/Train split                                 ############
############                                                                                    ############
############################################################################################################
############################################################################################################

## stratified split!
source('scripts/load_libs.R')
phenotypes = readRDS('db/phenotypes/phenotypes.rds') %>%
    filter(sex == 'Male', ancestry != 'Other')

split_ratio = 2/3
pheno1 = phenotypes %>% mutate(split_criteria = paste0(ancestry, '_', type))
set.seed(2022)
tt_split = initial_split(pheno1, prop = split_ratio, strata = 'split_criteria')
training_pheno = training(tt_split)
testing_pheno = testing(tt_split)
pheno_combined = rbind(training_pheno %>% mutate(train_test = 'Training'),
                       testing_pheno %>% mutate(train_test = 'Testing'))
setkey(pheno_combined, 'eid')

pheno_combined$train_test = factor(pheno_combined$train_test, levels = c('Training', 'Testing'))
pheno_combined$type = factor(pheno_combined$type, levels = c('Control', 'Case'))

saveRDS(pheno_combined, 'db/phenotypes/pheno_male_split.rds')

############################################################################################################
############################################################################################################
############                                   Supp table stats                                 ############
############                                                                                    ############
############################################################################################################
############################################################################################################
source('scripts/load_libs.R')
library(tableone)
phenotypes = readRDS('db/phenotypes/pheno_male_split.rds') %>%
    filter(sex == 'Male', ancestry != 'Other')

var_label(phenotypes) = list(ancestry = 'Ancestry', birth_year = 'Birth Year',
                             smoking_status = 'Smoking Status', alcohol_status = 'Alcohol Status',
                             type = 'Status')

#' yshah Monday, Jan 23, 2023 11:13:07 AM
#' Case vs Control
tab = CreateTableOne(vars = c('ancestry', 'birth_year', 'bmi', 'smoking_status', 'alcohol_status'),
                     strata = 'type',
                     data = phenotypes,
                     test = TRUE,
                     testNormal = t.test)
tab = print(tab, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.table(tab, paste0(basesave, 'tables/Supplementary_Table_1.tsv'), sep = '\t', quote = F)

#' yshah Monday, Jan 23, 2023 11:13:17 AM
#' Training vs testing
tab = CreateTableOne(vars = c('type', 'ancestry', 'birth_year', 'bmi', 'smoking_status', 'alcohol_status'),
                     strata = 'train_test',
                     data = phenotypes,
                     test = TRUE,
                     testNormal = t.test)
tab = print(tab, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.table(tab, paste0(basesave, 'tables/Supplementary_Table_2.tsv'), sep = '\t', quote = F)

############################################################################################################
############################################################################################################
############                                   Match phenotypes                                 ############
############                                                                                    ############
############################################################################################################
############################################################################################################

#' yshah Sunday, Apr 10, 2022 05:08:39 PM
#' Match case control for males
source('scripts/load_libs.R')
phenotypes = readRDS('db/phenotypes/pheno_male_split.rds') %>%
    filter(sex == 'Male', ancestry != 'Other')


match_data = function(dat, basefile = 'training', match_ratio = 1, match_method = 'nearest'){
    mkdir(paste0(basesave, 'UKBB_matching/'))
    message(paste0(Sys.time(), ' Starting matching'))
    t0 = Sys.time()
    mod_match = matchit(type ~ age_at_assessment + family_history + bmi + ancestry + smoking_status + alcohol_status,
                        data = dat, method = match_method, ratio = match_ratio)
    t1 = Sys.time()
    message(paste0('Finished matching in ', (t1-t0)/60 , ' minutes'))
    filename = paste0(basesave, 'UKBB_matching/', basefile, '_', match_method, '_ratio_', match_ratio, '.png')
    png(filename, width = 5, height = 4.5, unit = 'in', res = 300)
    print(love.plot(mod_match) + scale_color_aaas())
    dev.off()
    return(match.data(mod_match))
}

mkdir('db/phenotypes/matched/')

saveRDS(match_data(phenotypes %>% filter(train_test == 'Training'),
                   basefile = 'training',
                   match_method = 'nearest', match_ratio = 1),
        'db/phenotypes/matched/training_nearest_1.rds')
saveRDS(match_data(phenotypes %>% filter(train_test == 'Training'),
                   basefile = 'training',
                   match_method = 'nearest', match_ratio = 2),
        'db/phenotypes/matched/training_nearest_2.rds')
saveRDS(match_data(phenotypes %>% filter(train_test == 'Training'),
                   basefile = 'training',
                   match_method = 'nearest', match_ratio = 3),
        'db/phenotypes/matched/training_nearest_3.rds')
saveRDS(match_data(phenotypes %>% filter(train_test == 'Training'),
                   basefile = 'training',
                   match_method = 'nearest', match_ratio = 4),
        'db/phenotypes/matched/training_nearest_4.rds')

saveRDS(match_data(phenotypes %>% filter(train_test == 'Testing'),
                   basefile = 'testing',
                   match_method = 'nearest', match_ratio = 1),
        'db/phenotypes/matched/testing_nearest_1.rds')
saveRDS(match_data(phenotypes %>% filter(train_test == 'Testing'),
                   basefile = 'testing',
                   match_method = 'nearest', match_ratio = 2),
        'db/phenotypes/matched/testing_nearest_2.rds')
saveRDS(match_data(phenotypes %>% filter(train_test == 'Testing'),
                   basefile = 'testing',
                   match_method = 'nearest', match_ratio = 3),
        'db/phenotypes/matched/testing_nearest_3.rds')
saveRDS(match_data(phenotypes %>% filter(train_test == 'Testing'),
                   basefile = 'testing',
                   match_method = 'nearest', match_ratio = 4),
        'db/phenotypes/matched/testing_nearest_4.rds')



############################################################################################################
############################################################################################################
############                                 Functions for scoring                              ############
############                                                                                    ############
############################################################################################################
############################################################################################################


## Setup functions to train dataset
make_model = function(training.data, dependant.var, main.predictor, covars, nrepeats = 3, nfolds = 3,
                      saveModel = FALSE, savePath = ''){
    print(paste0(Sys.time(), ' Training ', main.predictor, '...'))
    form = as.formula(paste0(dependant.var, ' ~ ', main.predictor, ' + ', paste0(covars, collapse = ' + ')))
    form0 = as.formula(paste0(dependant.var, ' ~ ', paste0(covars, collapse = ' + ')))
    train_control = trainControl(method = 'repeatedcv', number = nfolds, repeats = nrepeats,
                                 savePredictions = 'all', classProbs = TRUE,
                                 summaryFunction = twoClassSummary)
    mod = train(form = form, data = training.data, method = 'glm', family = 'binomial', trControl = train_control, metric = 'ROC')
    mod0 = train(form = form0, data = training.data, method = 'glm', family = 'binomial', trControl = train_control, metric = 'ROC')
    if (saveModel & savePath != ''){
        message(paste0('Saving model: ', main.predictor))
        mod.save = stripGlmLR(mod$finalModel)
        saveRDS(mod.save, paste0(savePath, '/', main.predictor, '.rds'))
    }
    return(list(fullModel = mod, nullModel = mod0))
}

prediction_aucs = function(model, testing.data, dependant.var = 'type', main.predictor, metric = 'closest.topleft'){
    mod = model
    ## preds = predict(mod, testing.data, type = 'prob')
    ## mod.eval = MLeval::evalm(data.frame(preds, testing.data %>% pull(all_of(dependant.var))), showplots = FALSE, silent = TRUE)
    roc.obj = roc(testing.data %>% pull(all_of(dependant.var)) ~ predict(mod, testing.data, type = 'prob')$Case, quiet = TRUE)
    aucs = ci.auc(roc.obj) %>% as.numeric
    metrics = coords(roc.obj, 'best', ret = c('sens', 'spec', 'ppv', 'npv'), as.matrix = TRUE, best.method = metric)
    dt = data.table(PGS = main.predictor,
                    AUROC = aucs[2],
                    AUROC.lower = aucs[1],
                    AUROC.upper = aucs[3])
    dt = cbind(dt, metrics)
    setDT(dt)
    return(dt)
}

get_nagelkerke = function(mod, mod0, main.predictor){
    n = fmsb::NagelkerkeR2(mod$finalModel)$R2
    n0 = fmsb::NagelkerkeR2(mod0$finalModel)$R2
    dt = data.table(PGS = main.predictor, nagelkerke_full = n, nagelkerke_null = n0)
    return(dt)
}


############################################################################################################
############################################################################################################
############                                Plot the phenotypes                                 ############
############                                                                                    ############
############################################################################################################
############################################################################################################


#' yshah Sunday, Apr 10, 2022 04:50:59 PM
#' Plot phenotypes
source('scripts/load_libs.R')
phenotypes = readRDS('db/phenotypes/pheno_male_split.rds') %>%
    mutate(ancestry = ifelse(ancestry == 'EAS', 'ASN', ancestry)) %>%
    filter(sex == 'Male', ancestry != 'Other')

plot_pcs = function(dat, x = 'PC1', y = 'PC2'){
    dat %>% filter(sex == 'Male') %>%
        arrange(desc(ancestry)) %>%
        ggscatter(x = x, y = y, color = 'ancestry', alpha = 1/10, size = 1/2) + theme_few(14) +
        scale_color_manual(values = cols[c('AFR', 'ASN', 'EUR')]) +
        labs(color = 'Population') +
        guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
}

plt1 = plot_pcs(phenotypes)
plt2 = plot_pcs(phenotypes, 'PC3', 'PC4')
plt3 = plot_pcs(phenotypes, 'PC5', 'PC6')
plt4 = plot_pcs(phenotypes, 'PC7', 'PC8')

png(paste0(basesave, 'PCA_UKBB_all_males.png'), width = 6, height = 6, unit = 'in', res = 300)
ggarrange(plotlist = list(plt1, plt2, plt3, plt4), align = 'hv', nrow = 2, ncol = 2, common.legend = TRUE)
dev.off()

## Map current annotations to Prive et al
PC_UKBB <- bigreadr::fread2("~/data/UKBiobank/phenotypes/ukb26867.csv.gz", select = paste0("22009-0.", 1:16), nThread = 30)
all_centers <- read.csv("https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv", stringsAsFactors = FALSE)
all_sq_dist <- apply(all_centers[-1], 1, function(one_center) {
      rowSums(sweep(PC_UKBB, 2, one_center, '-')^2)
})
thr_sq_dist <- max(dist(all_centers[-1])^2) * 0.002 / 0.16
group <- apply(all_sq_dist, 1, function(x) {
    grp <- NA
    ind <- which.min(x)
    if (isTRUE(x[ind] < thr_sq_dist)) {
        grp <- all_centers$Ancestry[ind]
        ## We used a more stringent cutoff for the Ashkenazi group
        if (grp == "Ashkenazi" && x[ind] > 12.5^2) grp <- NA
    }
    grp
})
eid = fread('scripts/temp/ukb26867/eid.csv', colClasses = 'character')
prive = data.table(eid = eid$eid, Prive = group) %>%
    inner_join(phenotypes %>% dplyr::select(eid, ancestry, sex))
prive[is.na(prive)] = 'Other'

png(paste0(basesave, 'population_classification_vs_prive_hm.png'), width = 5, height = 5, unit = 'in', res = 300)
dt = prive %>% dplyr::count(ancestry, Prive) %>% group_by(Prive) %>%
    dplyr::summarise(ancestry, frac = 100*n/sum(n)) %>%
    pivot_wider(id_cols = Prive, names_from = ancestry, values_from = frac, values_fill = 0) %>%
    column_to_rownames('Prive')
Heatmap(dt, col = viridis::viridis(200),
        name = 'Percent of\nPrive',
        row_title = 'Prive et al., 2022 (16 PCs)',
        column_title = 'Our Populations (40 PCs)')
dev.off()


############################################################################################################
############################################################################################################
############                                   Evaluate PGScatalog                              ############
############                                                                                    ############
############################################################################################################
############################################################################################################

source('scripts/load_libs.R')
phenotypes = readRDS('db/phenotypes/pheno_male_split.rds') %>%
    filter(sex == 'Male', ancestry != 'Other')

pgs.dl = 'files/pgscatalog/downloaded/'
fl = list.files(pgs.dl)

## Process scores - save only rsid effect allele and effect weight cols. Include header
pgs_ss = lapply(fl, function(this.fl){
    if (file.ready(paste0('files/pgscatalog/processed/', gsub('.{3}$', '', basename(this.fl))))) return(NA)
    message(this.fl)
    dat = fread(paste0(pgs.dl, this.fl)) %>%
        dplyr::select(rsID, effect_allele, effect_weight)
    colnames(dat)[3] = gsub('.{7}$', '', basename(this.fl))
    fwrite(dat, paste0('files/pgscatalog/processed/', gsub('.{3}$', '', basename(this.fl))),
           quote = F, row.names = F, sep = '\t', col.names = F)
})

pairs = data.table(pair = gsub('.{4}$', '', list.files('files/pgscatalog/processed/')),
                   pgs_file = paste0('files/pgscatalog/processed/', list.files('files/pgscatalog/processed/')) %>% normalizePath(),
                   cores = 20, mem = 30)
setDT(pairs); setkey(pairs, 'pair')

pgs.jb = Job('tasks/plink_score.task', pairs, rootdir = 'Flow/pgscatalog/', mem = 30, cores = 20)
run(pgs.jb['ready'], mc.cores = 1)
Flow::update(pgs.jb)

scores = outputs(pgs.jb)
scores = mclapply(scores$pair, function(this.pair){
    dat = fread(scores[this.pair, scores], colClasses = 'character')
    dat = data.table(FID = dat$FID, IID = dat$IID, SCORE = dat$total_sum, PGS = this.pair)
    dat$SCORE = as.numeric(dat$SCORE)
    return(dat)
}, mc.cores = 10) %>% rbindlist

scores = scores %>%
    pivot_wider(names_from = PGS, values_from = SCORE, id_cols = FID) %>%
    dplyr::rename(eid = FID)

nfolds = 5
nrepeats = 5

training.dat = phenotypes %>% filter(train_test == 'Training') %>% merge(scores, by = 'eid')
testing.dat = phenotypes %>% filter(train_test == 'Testing') %>% merge(scores, by = 'eid')

pgs_models = mclapply(names(scores[,2:ncol(scores)]), function(pgs){
    mod = make_model(training.data = training.dat, dependant.var = 'type', main.predictor = pgs,
                     covars = paste0('PC', 1:10), nrepeats = nrepeats, nfolds = nfolds)
    mod1 = mod$fullModel
    roc.obj = roc(testing.dat$type ~ predict(mod1, testing.dat, type = 'prob')$Case, quiet = TRUE)
    curve.df = data.table(PGS = pgs,
                          sensitivity = roc.obj$sensitivities,
                          specificity = roc.obj$specificities,
                          AUC = paste0(round(as.numeric(ci.auc(roc.obj))[2],2 ), ' [',
                                       round(as.numeric(ci.auc(roc.obj))[1], 2), '-',
                                       round(as.numeric(ci.auc(roc.obj))[3], 2), ']'))
    anc.aucs = lapply(c('AFR', 'EAS', 'EUR'), function(this.anc){
        anc.test = testing.dat %>% filter(ancestry == this.anc)
        roc.obj = roc(anc.test$type ~ predict(mod1, anc.test, type = 'prob')$Case, quiet = TRUE)
        anc.auc = data.table(t(as.numeric(ci.auc(roc.obj)))) %>% mutate(PGS = pgs, ancestry = this.anc)
        colnames(anc.auc)[1:3] = c('AUC.low', 'AUC', 'AUC.up')
        return(anc.auc)
    }) %>% rbindlist()
    return(list(curve.df, anc.aucs, mod1))
}, mc.cores = 6)

overall.auc = pgs_models %>% purrr::map(1) %>% rbindlist
ancestry.auc = pgs_models %>% purrr::map(2) %>% rbindlist %>% mutate(ancestry = ifelse(ancestry == 'EAS', 'ASN', ancestry))
models = pgs_models %>% purrr::map(3)

mkdir(paste0(basesave, 'PGScatalog_eval/'))

png(paste0(basesave, 'PGScatalog_eval/total_pop_auc_curve.png'), width = 8, height = 6, unit = 'in', res = 300)
overall.auc %>%
    mutate(lab = paste0(PGS, ' (', AUC, ')')) %>%
    ggplot(aes(x = (1-specificity), y = sensitivity)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    facet_wrap(~lab) +
    labs(x = '1 - Specificity', y = 'Sensitivity') +
    theme_bw(14)
dev.off()

png(paste0(basesave, 'PGScatalog_eval/ancestry_evaluation.png'), width = 8, height = 6, unit = 'in', res = 300)
ancestry.auc %>%
    ggplot(aes(x = AUC, xmin = AUC.low, xmax = AUC.up, y = PGS, fill = ancestry)) +
    geom_col(position = position_dodge()) + geom_errorbar(width = 0.3, position = position_dodge(0.9)) +
    scale_fill_manual(values = cols[c('AFR', 'ASN', 'EUR')]) +
    labs(x = 'AUROC', y = NULL, fill = 'Evaluated\nPopulation') +
    theme_few(14) +
    coord_flip() +
    theme(legend.position = 'right')
dev.off()

mean_cl_normal(ancestry.auc %>% filter(ancestry == 'EUR') %>% pull(AUC)) 

png(paste0(basesave, 'PGScatalog_eval/ancestry_average_aucs.png'), width = 3, height = 4, unit = 'in', res = 300)
ancestry.auc %>%
    ggplot(aes(x = ancestry, y = AUC)) +
    stat_summary(fun = mean, geom = 'point', size = 3, shape = 15) +
    stat_summary(fun.data = mean_cl_normal, geom = 'errorbar', width = 0.3) +
    geom_beeswarm(aes(color = ancestry)) +
    ylim(c(0.52, 0.75)) +
    geom_line(aes(group = PGS), color = 'gray70', lty = 2) +
    stat_compare_means(method = 't.test', paired = TRUE,
                       comparison = list(c('AFR', 'ASN'), c('AFR', 'EUR')),
                       label.y = c(0.73, 0.7)) +
    labs(x = NULL) +
    scale_color_manual(values = cols[c('AFR', 'EAS', 'EUR')]) +
    theme_few(14) +
    theme(legend.position = 'none')
dev.off()

auc_pvals = lapply(1:6, function(i){
    mod = models[i][[1]]
    afr = testing.dat %>% filter(ancestry == 'AFR')
    eas = testing.dat %>% filter(ancestry == 'EAS')
    eur = testing.dat %>% filter(ancestry == 'EUR')
    roc.afr = roc(afr$type ~ predict(mod, afr, type = 'prob')$Case, quiet = TRUE)
    roc.eur = roc(eur$type ~ predict(mod, eur, type = 'prob')$Case, quiet = TRUE)
    roc.eas = roc(eas$type ~ predict(mod, eas, type = 'prob')$Case, quiet = TRUE)
    afr = roc.test(roc.afr, roc.eur)$p.value
    eas = roc.test(roc.eas, roc.eur)$p.value
    dt = data.table(PGS = names(scores)[i+1],
                    AFR = afr, EAS = eas)
}) %>% rbindlist()

auc_pvals


############################################################################################################
############################################################################################################
############                               Correlate created scores                             ############
############                                                                                    ############
############################################################################################################
############################################################################################################

source('scripts/load_libs.R')
scores = readRDS('db/combined_scores.rds')
scores[,2:ncol(scores)] = -scores[,2:ncol(scores)]
phenotypes = readRDS('db/phenotypes/pheno_male_split.rds') %>%
    filter(sex == 'Male', ancestry != 'Other')

cor.score = cor(scores[, 2:ncol(scores)])
mean_cl_normal(as.numeric(cor.score[upper.tri(cor.score, diag = F)]))

ha0 = data.frame(Derivation = str_split(colnames(cor.score), '\\.') %>% purrr::map(1) %>% unlist,
                Method = str_split(colnames(cor.score), '\\.') %>% purrr::map(3) %>% unlist,
                row.names = colnames(cor.score)) %>%
    mutate(Derivation = case_when(
               Derivation == 'african' ~ 'AFR',
               Derivation == 'european' ~ 'EUR',
               Derivation == 'eastasian' ~ 'ASN',
               Derivation == 'total' ~ 'Total'
           ))
ha = HeatmapAnnotation(Derivation = ha0$Derivation, Method = ha0$Method,
                       col = list(Derivation = cols[c('AFR', 'ASN', 'EUR', 'Total')]),
                       annotation_name_gp= gpar(fontsize = 10))
hm = Heatmap(cor.score, top_annotation = ha, show_row_names = F, show_column_names = F,
             name = 'Pearson', col = circlize::colorRamp2(breaks = c(-1, 0, 1), colors = c('blue', 'white', 'red')))
hc = fastcluster::hclust(dist(cor.score))
my_clusters = dendextend::cutree(hc, k = 3)

ht = draw(Heatmap(cor.score, top_annotation = ha, show_row_names = F, show_column_names = F,
                  name = 'Pearson', col = circlize::colorRamp2(breaks = c(-1, 0, 1), colors = c('blue', 'white', 'red')),
                  split = my_clusters, column_split = my_clusters))

dev.off()
png(paste0(basesave, 'score_pearson_correlation.png'), width = 5.5, height = 5, unit = 'in', res = 600)
print(ht)
dev.off()

tab = as.data.frame(my_clusters) %>%
    rownames_to_column('score') %>%
    data.table() %>%
    mutate(ancestry = str_split(score, '\\.') %>% purrr::map(1) %>% unlist(),
           anc = ifelse(ancestry == 'total', 'european', ancestry),
           method = str_split(score, '\\.') %>% purrr::map(3) %>% unlist(),
           combo = paste0(anc, method))
chisq.test(table(tab$my_clusters, tab$anc))
chisq.test(table(tab$my_clusters, tab$ancestry))
tab = tab[method != 'xpass']
chisq.test(table(tab$my_clusters, tab$method))

tb = tab %>% filter(anc %in% c('african', 'eastasian'))
chisq.test(table(tb$my_clusters, tb$method))

rm(cor.score)
gc()

## Repeat for cases onle
scores = scores %>% filter(eid %in% phenotypes[type == 'Case', eid])
cor.score = cor(scores[, 2:ncol(scores)])
mean_cl_normal(as.numeric(cor.score[upper.tri(cor.score, diag = F)]))

ha0 = data.frame(Derivation = str_split(colnames(cor.score), '\\.') %>% purrr::map(1) %>% unlist,
                Method = str_split(colnames(cor.score), '\\.') %>% purrr::map(3) %>% unlist,
                row.names = colnames(cor.score)) %>%
    mutate(Derivation = case_when(
               Derivation == 'african' ~ 'AFR',
               Derivation == 'european' ~ 'EUR',
               Derivation == 'eastasian' ~ 'ASN',
               Derivation == 'total' ~ 'Total'
           ))
ha = HeatmapAnnotation(Derivation = ha0$Derivation, Method = ha0$Method,
                       col = list(Derivation = cols[c('AFR', 'ASN', 'EUR', 'Total')]),
                       annotation_name_gp= gpar(fontsize = 10))
hm = Heatmap(cor.score, top_annotation = ha, show_row_names = F, show_column_names = F,
             name = 'Pearson', col = circlize::colorRamp2(breaks = c(-1, 0, 1), colors = c('blue', 'white', 'red')))
ht = draw(Heatmap(cor.score, top_annotation = ha, show_row_names = F, show_column_names = F,
                  cluster_rows = function(x) fastcluster::hclust(dist(x)),
                  cluster_columns = function(x) fastcluster::hclust(dist(x)),
                  name = 'Pearson', col = circlize::colorRamp2(breaks = c(-1, 0, 1), colors = c('blue', 'white', 'red'))))
dg = column_dend(ht)

dev.off()
png(paste0(basesave, 'score_pearson_correlation_only_cases.png'), width = 5, height = 4.5, unit = 'in', res = 600)
print(hm)
dev.off()

rm(cor.score)
gc()


############################################################################################################
############################################################################################################
############                               Evaluate created scores                              ############
############                                                                                    ############
############################################################################################################
############################################################################################################

source('scripts/load_libs.R')
scores = readRDS('db/combined_scores.rds')
scores[,2:ncol(scores)] = -scores[,2:ncol(scores)]
phenotypes = readRDS('db/phenotypes/pheno_male_split.rds') %>%
    filter(sex == 'Male', ancestry != 'Other')


## Put functions together for SPORE
spore_models = function(training.data, testing.data, scores, dependant.var = 'type', main.predictor, covars, nrepeats = 5, nfolds = 5,
                        saveModel = FALSE, savePath = '', metric = 'closest.topleft'){
    ## Setup data
    scores = scores %>% dplyr::select(eid, all_of(main.predictor))
    training.data = merge(training.data, scores %>% mutate(eid = as.character(eid)), by = 'eid')
    testing.data = merge(testing.data, scores %>% mutate(eid = as.character(eid)), by = 'eid')
    ## Train model
    models = make_model(training.data = training.data, dependant.var = dependant.var, main.predictor = main.predictor,
                        covars = paste0('PC', 1:10), saveModel = saveModel, savePath = savePath)
    ## Pseudo R2
    nagelkerke = get_nagelkerke(models$fullModel, models$nullModel, main.predictor)
    ## Evaluation AUCs
    evals = lapply(c('AFR', 'EAS', 'EUR'), function(anc){
        pred = prediction_aucs(model = models$fullModel, testing.data = testing.data %>% filter(ancestry == anc),
                               main.predictor = main.predictor) %>%
            mutate(Evaluation = anc)
    }) %>% rbindlist()
    evals = rbind(evals,
                  prediction_aucs(model = models$fullModel, testing.data = testing.data, main.predictor = main.predictor) %>%
                  mutate(Evaluation = 'Total'))
    ## Format results
    nagelkerke = nagelkerke %>%
        mutate(Derivation = str_split(PGS, '\\.') %>% purrr::map(1) %>% unlist,
           Method = str_split(PGS, '\\.') %>% purrr::map(3) %>% unlist,
           Derivation = case_when(
               Derivation == 'african' ~ 'AFR',
               Derivation == 'european' ~ 'EUR',
               Derivation == 'eastasian' ~ 'EAS',
               Derivation == 'total' ~ 'Total'
           )) %>%
        dplyr::select(PGS, Method, Derivation, nagelkerke_full, nagelkerke_null)
    evals = evals %>%
        mutate(Derivation = str_split(PGS, '\\.') %>% purrr::map(1) %>% unlist,
               Method = str_split(PGS, '\\.') %>% purrr::map(3) %>% unlist,
               Derivation = case_when(
                   Derivation == 'african' ~ 'AFR',
                   Derivation == 'european' ~ 'EUR',
                   Derivation == 'eastasian' ~ 'EAS',
                   Derivation == 'total' ~ 'Total'
               )) %>%
        dplyr::select(PGS, Method, Derivation, Evaluation, everything())
    if (saveModel & savePath != ''){
        model_paths = data.table(PGS = main.predictor,
                                 model_path = normalizePath(paste0(savePath, '/', main.predictor, '.rds')))
        return(list(evals, nagelkerke, model_paths))
    }
    return(list(evals, nagelkerke))
}

## Get SPORE models
pgs.names = colnames(scores)[-1]

mkdir('db/models/SPORE/matching_1/')
match1.spore = mclapply(pgs.names, function(this.pgs){
    spore_models(training.data = readRDS('db/phenotypes/matched/training_nearest_1.rds'),
                 testing.data = readRDS('db/phenotypes/matched/testing_nearest_1.rds'),
                 scores = scores, dependant.var = 'type', main.predictor = this.pgs,
                 covars = paste0('PC', 1:10), saveModel = TRUE,
                 savePath = normalizePath('db/models/SPORE/matching_1/'))
}, mc.cores = 25)

mkdir('db/models/SPORE/matching_2/')
match2.spore = mclapply(pgs.names, function(this.pgs){
    spore_models(training.data = readRDS('db/phenotypes/matched/training_nearest_2.rds'),
                 testing.data = readRDS('db/phenotypes/matched/testing_nearest_2.rds'),
                 scores = scores, dependant.var = 'type', main.predictor = this.pgs,
                 covars = paste0('PC', 1:10), saveModel = TRUE,
                 savePath = normalizePath('db/models/SPORE/matching_2/'))
}, mc.cores = 25)

mkdir('db/models/SPORE/matching_3/')
match3.spore = mclapply(pgs.names, function(this.pgs){
    spore_models(training.data = readRDS('db/phenotypes/matched/training_nearest_3.rds'),
                 testing.data = readRDS('db/phenotypes/matched/testing_nearest_3.rds'),
                 scores = scores, dependant.var = 'type', main.predictor = this.pgs,
                 covars = paste0('PC', 1:10), saveModel = TRUE,
                 savePath = normalizePath('db/models/SPORE/matching_3/'))
}, mc.cores = 25)

mkdir('db/models/SPORE/no_matching/')
nomatch.spore = mclapply(pgs.names, function(this.pgs){
    spore_models(training.data = phenotypes %>% filter(train_test == 'Training'),
                 testing.data = phenotypes %>% filter(train_test == 'Testing'),
                 scores = scores, dependant.var = 'type', main.predictor = this.pgs,
                 covars = paste0('PC', 1:10), saveModel = TRUE,
                 savePath = 'db/models/SPORE/no_matching/')
}, mc.cores = 25)

spore.models = rbind(
    match1.spore %>% purrr::map(3) %>% rbindlist %>% mutate(Matching = '1:1'),
    match2.spore %>% purrr::map(3) %>% rbindlist %>% mutate(Matching = '1:2'),
    match3.spore %>% purrr::map(3) %>% rbindlist %>% mutate(Matching = '1:3'),
    nomatch.spore %>% purrr::map(3) %>% rbindlist %>% mutate(Matching = 'None')
)

spore.r2 = rbind(
    match1.spore %>% purrr::map(2) %>% rbindlist %>% mutate(Matching = '1:1'),
    match2.spore %>% purrr::map(2) %>% rbindlist %>% mutate(Matching = '1:2'),
    match3.spore %>% purrr::map(2) %>% rbindlist %>% mutate(Matching = '1:3'),
    nomatch.spore %>% purrr::map(2) %>% rbindlist %>% mutate(Matching = 'None')
)

spore.auc = rbind(
    match1.spore %>% purrr::map(1) %>% rbindlist %>% mutate(Matching = '1:1'),
    match2.spore %>% purrr::map(1) %>% rbindlist %>% mutate(Matching = '1:2'),
    match3.spore %>% purrr::map(1) %>% rbindlist %>% mutate(Matching = '1:3'),
    nomatch.spore %>% purrr::map(1) %>% rbindlist %>% mutate(Matching = 'None')
)

## fwrite(spore.models, 'db/models/SPORE/model_paths.txt', row.names = F, quote = F, sep = '\t')
## fwrite(spore.r2, 'db/models/SPORE/model_r2.txt', row.names = F, quote = F, sep = '\t')
## fwrite(spore.auc, 'db/models/SPORE/model_auc.txt', row.names = F, quote = F, sep = '\t')

spore.models = fread('db/models/SPORE/model_paths.txt')
spore.r2 = fread('db/models/SPORE/model_r2.txt')
spore.auc = fread('db/models/SPORE/model_auc.txt')

## Plot the effect of matching
png(paste0(basesave, 'UKBB_matching/matching_auc.png'), width = 11, height = 6, unit = 'in', res = 300)
spore.auc[, .(PGS, Derivation, Evaluation, AUROC, sensitivity, specificity, Matching)] %>%
    dplyr::rename(Sensitivity = sensitivity, Specificity = specificity) %>%
    melt(c('PGS', 'Derivation', 'Evaluation', 'Matching')) %>%
    mutate(Derivation = ifelse(Derivation == 'EAS', 'ASN', Derivation),
           Evaluation = ifelse(Evaluation == 'EAS', 'ASN', Evaluation),
           Matching = factor(Matching, levels = c('1:1', '1:2', '1:3', 'None'))) %>%
    ggplot(aes(x = Matching, y = value)) +
    geom_jitter(alpha = 1/3, aes(color = Derivation)) +
    geom_boxplot(alpha = 2/3, outlier.shape = NA) +
    labs(x = 'Matching', y = 'Value', color = 'Derivation') +
    facet_grid(vars(variable), vars(Evaluation)) +
    stat_compare_means(label.y.npc = 0.8, label.x.npc = 0.1) +
    scale_color_manual(values = cols, limits = force) +
    theme_few(14) +
    guides(color = guide_legend(override.aes = list(size = 6, alpha = 1)))
dev.off()

my.match = 'None'
aucs = spore.auc %>%
    filter(Matching == my.match) %>%
    mutate(Derivation = ifelse(Derivation == 'EAS', 'ASN', Derivation),
           Evaluation = ifelse(Evaluation == 'EAS', 'ASN', Evaluation))
r2s = spore.r2 %>%
    filter(Matching == my.match) %>%
    mutate(Derivation = ifelse(Derivation == 'EAS', 'ASN', Derivation))
mkdir(paste0(basesave, 'model_eval/'))

png(paste0(basesave, 'model_eval/nagelkerke.png'), width = 5, height = 3.5, unit = 'in',res = 300)
r2s %>%
    ggplot(aes(x = Derivation, y = nagelkerke_full)) +
    geom_quasirandom(aes(color = Method)) +
    geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 2/3) +
    labs(x = NULL, color = NULL, y = expression(R[Nagelkerke]^2)) +
    theme_few(14) +
    scale_color_manual(values = cols, limits = force) +
    guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

## Relative R2
clump.r2 = r2s %>% filter(Method == 'clump') %>% group_by(Derivation) %>% dplyr::summarize(clump.mean = mean(nagelkerke_full))
png(paste0(basesave, 'model_eval/nagelkerke_relative.png'), width = 6, height = 3.5, unit = 'in', res = 300)
r2s %>% merge(clump.r2) %>% mutate(relative_r2 = nagelkerke_full/clump.mean) %>%
    filter(Derivation != 'Total') %>%
    ggplot(aes(x = Derivation, y = relative_r2, fill = Method)) +
    stat_summary(fun = mean, geom = 'bar', position = position_dodge(0.9)) +
    geom_point(show.legend = F, position = position_dodge(0.9)) +
    scale_fill_manual(values = cols, limits = force) +
    theme_few(14) + labs(x = NULL, fill = NULL, y = expression(Relative~R[Nagelkerke]^2))
dev.off()

od = aucs %>% group_by(Method) %>% dplyr::summarize(med = median(AUROC)) %>% arrange(-med) %>% pull(Method)
png(paste0(basesave, 'model_eval/auc_by_method_only.png'), width = 6, height = 3.5, unit = 'in', res = 300)
aucs %>%
    mutate(Method = factor(Method, levels = od)) %>%
    ggplot(aes(x = Method, y = AUROC)) +
    geom_quasirandom(size = 1/2, aes(color = Evaluation)) +
    geom_boxplot(outlier.shape = NA, alpha = 2/3) +
    labs(x = NULL, color = 'Evaluation') +
    scale_color_manual(values = cols, limits = force) +
    stat_compare_means(method = 'wilcox',
                       comparisons = list(c('prscsx', 'prscs'),
                                          c('clump', 'prscsx'),
                                          c('clump', 'prscs')),
                       label = 'p.signif') +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_few(14)
dev.off()

## Heatmap for best model
conti.samps = data.frame(Cases = c(10368, 8611, 85554, 107247),
                         Controls = c(10986, 18809, 91972, 127006),
                         row.names = c('AFR', 'ASN', 'EUR', 'Total'))
conti.samps$Cases = log10(conti.samps$Cases); conti.samps$Controls = log10(conti.samps$Controls)
anno.conti = rowAnnotation(dbGaP = anno_barplot(conti.samps, border = FALSE,
                                                gp = gpar(fill = cols[c('Case', 'Control')] %>% unname,
                                                          col = cols[c('Case', 'Control')])))
testing.samps = data.frame(Cases = c(101, 49, 4231, 4381),
                           Controls = c(1142, 2183, 65709, 69034),
                           row.names = c(c('AFR', 'ASN', 'EUR', 'Total')))
testing.samps$Cases = log10(testing.samps$Cases); testing.samps$Controls = log10(testing.samps$Controls)
anno.testing = HeatmapAnnotation(UKBB = anno_barplot(testing.samps, border = TRUE,
                                                       gp = gpar(fill = cols[c('Case', 'Control')],
                                                                 col = cols[c('Case', 'Control')])))

hm.auc = aucs %>%
    mutate(comb = paste0(Derivation, '_', Evaluation)) %>%
    group_by(comb) %>%
    arrange(-AUROC) %>%
    dplyr::slice(1:1) %>%
    data.table() %>%
    pivot_wider(id_cols = Derivation, names_from = Evaluation, values_from = AUROC) %>%
    column_to_rownames('Derivation')
hm.method = aucs %>%
    mutate(comb = paste0(Derivation, '_', Evaluation)) %>%
    group_by(comb) %>%
    arrange(-AUROC) %>%
    dplyr::slice(1:1) %>%
    data.table() %>%
    pivot_wider(id_cols = Derivation, names_from = Evaluation, values_from = Method) %>%
    column_to_rownames('Derivation')

lgd_list = list(
    Legend(labels = c('Case', 'Control'), title = 'Status',
           legend_gp = gpar(fill = cols[c('Case', 'Control')]))
)

hm = Heatmap(hm.auc, name = 'AUROC', col = viridis::viridis(100),
             right_annotation = anno.conti, top_annotation = anno.testing,
             cell_fun = function(j, i, x, y, w, h, col) {
                 grid.text(hm.method[i, j], x, y)
             },
             row_title = 'Derivation', row_title_side = 'left',
             column_title = 'Evaluation', column_title_side = 'bottom')

png(paste0(basesave, 'model_eval/auc_heatmap.png'), width = 5, height = 5, unit = 'in', res = 300)
draw(hm)
dev.off()
ppng(draw(lgd_list[[1]]), width = 5, height = 5, unit = 'in', res = 300)

#' yshah Sunday, Jul 10, 2022 09:48:57 PM
#' Which scores are missing?
scores = readRDS('db/combined_scores.rds')
nms = names(scores)
nms = nms[-1]
dt = data.table(score = nms) %>%
    mutate(pop = str_split(score, '\\.') %>% purrr::map(1) %>% unlist,
           param = str_split(score, '\\.') %>% purrr::map(2) %>% unlist,
           method = str_split(score, '\\.') %>% purrr::map(3) %>% unlist,
           present = TRUE)

dt = dt %>%
    pivot_wider(id_cols = c('param', 'method'), names_from = pop, values_from = present, values_fill = FALSE)
dt$present = rowSums(dt[,3:6])

#' yshah Tuesday, May 03, 2022 09:35:37 PM
#' Look at best scores in EAS and AFR
source('scripts/load_libs.R')
scores = readRDS('db/combined_scores.rds')
scores[,2:ncol(scores)] = -scores[,2:ncol(scores)]
phenotypes = readRDS('db/phenotypes/pheno_male_split.rds') %>%
    filter(sex == 'Male', ancestry != 'Other')

covars = c(paste0('PC', 1:10), 'age_at_assessment', 'family_history', 'bmi')
this.match = 'None'
this.anc = 'EAS'
this.pgs = spore.auc %>% filter(Evaluation == this.anc, Matching == this.match) %>% arrange(-AUROC) %>% pull(PGS) %>% head(1)

mod = make_model(training.data = phenotypes %>% filter(train_test == 'Training') %>% merge(scores, by = 'eid'),
           dependant.var = 'type', main.predictor = this.pgs,
           covars = covars)
mod.eas = mod$fullModel$finalModel
res = phenotypes %>% filter(train_test == 'Testing') %>% merge(scores, by = 'eid') %>%
    mutate(family_historyTRUE = ifelse(family_history, 1, 0))
res2 = res %>% filter(ancestry == this.anc)
auc.eas = roc(res2$type ~ predict(mod$fullModel, res2, type = 'prob')$Case, quiet = TRUE)
mod.eas = mod$fullModel
imp.eas = data.frame(varImp(mod.eas)$importance) %>%
    rownames_to_column('feature') %>%
    dplyr::rename(EAS = Overall) %>%
    mutate(feature = ifelse(feature == this.pgs, 'PRS', feature))
res$logOR = predict(mod$fullModel$finalModel, res)
res.eas = res %>% dplyr::select(ancestry, type, logOR) %>% mutate(score = 'Best ASN', logOR = scale(logOR))

this.anc = 'AFR'
this.pgs = spore.auc %>% filter(Evaluation == this.anc, Matching == this.match) %>% arrange(-AUROC) %>% pull(PGS) %>% head(1)
mod = make_model(training.data = phenotypes %>% filter(train_test == 'Training') %>% merge(scores, by = 'eid'),
                 dependant.var = 'type', main.predictor = this.pgs,
                 covars = covars)
mad.afr = mod$fullModel$finalModel
res = phenotypes %>% filter(train_test == 'Testing') %>% merge(scores, by = 'eid') %>%
    mutate(family_historyTRUE = ifelse(family_history, 1, 0))
res2 = res %>% filter(ancestry == this.anc)
auc.afr = roc(res2$type ~ predict(mod$fullModel, res2, type = 'prob')$Case, quiet = TRUE)
mod.afr = mod$fullModel
imp.afr = data.frame(varImp(mod.afr)$importance) %>%
    rownames_to_column('feature') %>%
    dplyr::rename(AFR = Overall) %>%
    mutate(feature = ifelse(feature == this.pgs, 'PRS', feature))
res$logOR = predict(mod$fullModel$finalModel, res)
res.afr = res %>% dplyr::select(ancestry, type, logOR) %>% mutate(score = 'Best AFR', logOR = scale(logOR))

this.anc = 'EUR'
this.pgs = spore.auc %>% filter(Evaluation == this.anc, Matching == this.match) %>% arrange(-AUROC) %>% pull(PGS) %>% head(1)
mod = make_model(training.data = phenotypes %>% filter(train_test == 'Training') %>% merge(scores, by = 'eid'),
                 dependant.var = 'type', main.predictor = this.pgs,
                 covars = covars)
mod.eur = mod$fullModel$finalModel
res = phenotypes %>% filter(train_test == 'Testing') %>% merge(scores, by = 'eid') %>%
    mutate(family_historyTRUE = ifelse(family_history, 1, 0))
res2 = res %>% filter(ancestry == this.anc)
mod.eur = mod$fullModel
imp.eur = data.frame(varImp(mod.eur)$importance) %>%
    rownames_to_column('feature') %>%
    dplyr::rename(EUR = Overall) %>%
    mutate(feature = ifelse(feature == this.pgs, 'PRS', feature))
auc.eur = roc(res2$type ~ predict(mod$fullModel$finalModel, res2), quiet = TRUE)
res$logOR = predict(mod$fullModel$finalModel, res)
res.eur = res %>% dplyr::select(ancestry, type, logOR) %>% mutate(score = 'Best EUR', logOR = scale(logOR))

res = rbind(res.afr, res.eas, res.eur) %>%
    mutate(ancestry = ifelse(ancestry == 'EAS', 'ASN', ancestry))
png(paste0(basesave, 'model_eval/combined_genetic_clinical.png'), width = 6, height = 5, unit = 'in', res = 300)
ggdensity(res, x = 'logOR', fill = 'type') +
    facet_grid(vars(score), vars(ancestry)) +
    labs(x = 'Scaled Risk', y = 'Density', fill = NULL) +
    theme(legend.position = 'right') +
    scale_fill_manual(values = cols, limits = force) +
    theme_few(14)
dev.off()

## Age distribution of phenotypes
plt = phenotypes %>%
    mutate(ancestry = ifelse(ancestry == 'EAS', 'ASN', ancestry)) %>%
    ggplot(aes(x = type, y = age_at_assessment)) +
    geom_violin() + ggforce::geom_sina(aes(color = type)) + geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 2/3) +
    labs(x = NULL, y = 'Age at assessment', color = NULL) +
    facet_grid(vars(train_test), vars(ancestry)) +
    theme_few(14) +
    ylim(30, 80) +
    stat_compare_means(method = 't.test', comparisons = list(c('Control', 'Case'))) +
    scale_color_manual(values = cols, limits = force) +
    guides(color = guide_legend(override.aes = list(size = 5)))
png(paste0(basesave, 'age_distribution_ukbb.png'), width = 8, height = 5, unit = 'in', res = 300)
plt
dev.off()

## Get AUC for each added covariate
covar.list = list(c('PGS'), c('age_at_assessment', 'PGS'), c('age_at_assessment', 'PGS', 'family_history'),
                  c('age_at_assessment', 'PGS', 'family_history', 'bmi'), c('age_at_assessment', 'family_history', 'bmi'))
combined.data = merge(phenotypes, scores)
testing.data = combined.data %>% filter(train_test == 'Testing')
training.data = combined.data %>% filter(train_test == 'Training')

covariate_aucs = mclapply(c('AFR', 'EAS', 'EUR'), function(this.anc){
    this.pgs = spore.auc %>% filter(Evaluation == this.anc, Matching == 'None') %>% arrange(-AUROC) %>% pull(PGS) %>% head(1)
    evals = mclapply(covar.list, function(covar.l){
        covar.l = unlist(covar.l)
        covar.form = covar.l
        if (length(covar.l) > 1) covar.form = paste0(covar.l, collapse = '+')
        form = as.formula(paste0('type ~ ', paste0('PC', 1:10, collapse = '+'), '+', covar.form))
        set.seed(22)
        train_control = trainControl(method = 'repeatedcv', number = 5, repeats = 5,
                                     savePredictions = 'all', classProbs = TRUE,
                                     summaryFunction = twoClassSummary)
        tr.dat = training.data %>% dplyr::rename(PGS = all_of(this.pgs))
        ts.dat = testing.data %>% dplyr::rename(PGS = all_of(this.pgs))
        set.seed(22)
        mod = train(form = form, data = tr.dat, method = 'glm', family = 'binomial', trControl = train_control, metric = 'ROC')
        ts.dat = ts.dat %>% filter(ancestry == this.anc)
        ts.dat$preds = predict(mod, ts.dat, type = 'prob')$Case
        roc.obj = roc(ts.dat$type ~ ts.dat$preds)
        aucs = data.table(t(as.numeric(ci.auc(roc.obj))))
        colnames(aucs) = c('AUC.lower', 'AUC', 'AUC.upper')
        covar.l[covar.l == 'age_at_assessment'] = 'Age'
        covar.l[covar.l == 'family_history'] = 'Family history'
        covar.l[covar.l == 'bmi'] = 'BMI'
        aucs = aucs %>% mutate(PGS = this.pgs, Derivation = paste0('Best ', this.anc), Evaluation = this.anc,
                               covars = paste0(covar.l, collapse = ', ') )
        return(aucs)
    }, mc.cores = 5) %>% rbindlist
    return(evals)
}, mc.cores = 3)

covariate_aucs = rbindlist(covariate_aucs)

png(paste0(basesave, 'model_eval/auc_by_covariates.png'), width = 7, height = 4, unit = 'in', res = 300)
covariate_aucs %>%
    mutate(Derivation = ifelse(Derivation == 'Best EAS', 'Best ASN', Derivation),
           Evaluation = ifelse(Evaluation == 'EAS', 'ASN', Evaluation)) %>%
    filter(covars != 'Age, Family history, BMI') %>%
    mutate(covars = factor(covars, levels = c('PGS', 'Age, PGS', 'Age, PGS, Family history', 'Age, PGS, Family history, BMI') %>% rev)) %>%
    ggplot(aes(x = AUC, xmin = AUC.lower, xmax = AUC.upper, y = covars, label = round(AUC, 2))) +
    geom_col(position = position_dodge(0.9)) + geom_errorbar(width = 0.3, position = position_dodge(0.9)) +
    geom_text(color = 'white', position = position_dodge(), aes(x = 0.15)) +
    facet_wrap(~Derivation, nrow = 1) +
    labs(x = 'AUROC', y = NULL) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
    theme_few(14)
dev.off()

#' yshah Tuesday, May 17, 2022 01:21:53 PM
#' What is the AUC for models that have the best R2
best.r2 = r2s %>% group_by(Derivation) %>% arrange(-nagelkerke_full) %>% dplyr::slice(1:1) %>% pull(PGS)
aucs1 = aucs %>% filter(PGS %in% best.r2)
aucs2 = aucs %>% group_by(Evaluation) %>% arrange(-AUROC) %>% dplyr::slice(1:1) %>% pull(PGS)
aucs2 = aucs %>% filter(PGS %in% aucs2)


lapply(c('AFR', 'EAS', 'EUR'), function(this.anc){
    best.r2.pgs = r2s %>%
        mutate(Derivation = ifelse(Derivation == 'ASN', 'EAS', Derivation)) %>%
        filter(Derivation == this.anc) %>% arrange(-nagelkerke_full) %>% dplyr::slice(1:1) %>% pull(PGS)
    best.auc.pgs = aucs %>%
        mutate(Evaluation = ifelse(Evaluation == 'ASN', 'EAS', Evaluation)) %>%
        filter(Evaluation == this.anc) %>% arrange(-AUROC) %>% dplyr::slice(1:1) %>% pull(PGS)
    r2.mod = readRDS(spore.models %>% filter(Matching == 'None', PGS == best.r2.pgs) %>% pull(model_path))
    auc.mod = readRDS(spore.models %>% filter(Matching == 'None', PGS == best.auc.pgs) %>% pull(model_path))
    ts.dat = phenotypes %>% filter(train_test == 'Testing', ancestry == this.anc) %>% merge(scores)
    ts.dat$r2 = predict(r2.mod, ts.dat)
    ts.dat$auc = predict(auc.mod, ts.dat)
    r2.roc = roc(ts.dat$type ~ ts.dat$r2)
    auc.roc = roc(ts.dat$type ~ ts.dat$auc)
    ts.stat = roc.test(r2.roc, auc.roc)$p.value
    dt = rbind(data.table(pop = this.anc,
                          PGS = best.r2.pgs,
                          type = 'Best R2',
                          sensitivity = r2.roc$sensitivities,
                          specificity = r2.roc$specificities),
               data.table(pop = this.anc,
                          PGS = best.auc.pgs,
                          type = 'Best AUROC',
                          sensitivity = auc.roc$sensitivities,
                          specificity = auc.roc$specificities)
               ) %>%
        mutate(pop = ifelse(pop == 'EAS', 'ASN', pop))
    this.anc = ifelse(this.anc == 'EAS', 'ASN', this.anc)
    plt = dt %>%
        mutate(type = paste0(type, ' (', PGS, ')')) %>%
        ggplot(aes(x = 1-specificity, y = sensitivity, color = type)) +
        geom_line(size = 1.5, alpha = 2/3) + theme_few(14) +
        geom_abline(intercept = 0, lty = 2, slope = 1) +
        labs(x = '1-Specificity', y = 'Sensitivity', color = NULL, title = paste0(this.anc, ', p = ', round(ts.stat, 4))) +
        theme(legend.position = 'top') +
        guides(color=guide_legend(nrow=2,byrow=TRUE))
    system(paste0('mkdir -p ', basesave, 'model_eval/best_auc_vs_best_r2/'))
    png(paste0(basesave, 'model_eval/best_auc_vs_best_r2/', this.anc, '.png'), width = 4.5, height = 4.5, unit = 'in', res = 300)
    print(plt)
    dev.off()
})
