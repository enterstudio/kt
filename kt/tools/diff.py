import os
import time
import sys

from random import randint
from random import gauss
from argparse import ArgumentParser

import numpy as np
import pandas as pd

import h5py

from scipy.stats import mannwhitneyu
from pyensembl import EnsemblRelease


def get_foldchange(scores_tpm, num_query):
    c = 1.00
    mean_query = np.mean(scores_tpm[:num_query])
    mean_subject = np.mean(scores_tpm[num_query:])
    log2fc = np.log2(mean_query + c) - np.log2(mean_subject + c)

    return (log2fc, mean_query, mean_subject)


def get_median(data):
    if len(data) % 2 == 0:
        return (data[int(len(data)*0.5) - 1] + data[int(len(data)*0.5)]) / 2
    return data[len(data)/2]


def trans_to_gene(values, ids_by_genes, transcripts):
    trans_by_gene = []
    tmp = []
    for ids in ids_by_genes:
        tmp.append(sum(values[ids]))
        trans_by_gene.append(','.join(transcripts[ids]))

    return (tmp, trans_by_gene)


def auc_u_test(scores_tpm, num_query, num_subject):
    # Don't need to sort when we use mannwhitneyu in auc_u_test
    (u_value, p_value) = mannwhitneyu(scores_tpm[:num_query], scores_tpm[num_query:], alternative="two-sided")
    auc = u_value / (num_query * num_subject)

    return(auc, p_value)


def cut_version(id_str):
    index_version = id_str.find('.')
    if index_version != -1:
        id_str = id_str[0:index_version]

    return (id_str)


def create_bs_sample(scores_tpm, num_query, num_subject):
    scores_tpm_bs = np.copy(scores_tpm)

    for i in xrange(num_query):
        i_bs = randint(0, num_query-1)
        scores_tpm_bs[i] = scores_tpm[i_bs]

    for i in xrange(num_subject):
        i_bs = randint(0, num_subject-1)
        scores_tpm_bs[num_query + i] = scores_tpm[num_query + i_bs]

    return (scores_tpm_bs)


def all_auc(file_out, list_ids, design, num_query, num_subject, kallisto,
            num_kal_bs, num_sample_bs, min_exp, min_l2fc, ensembl,
            trans_by_gene=None):
    print(time.strftime('%X') + ": computing AUC and FC on: " + str(len(list_ids)) + " ids")

    # If need it, you can add sample
    # df = pd.DataFrame({'class':design[args.CLASS].values,'score':0}) #,'sample':design["sample"]})

    with open(file_out, 'w') as w_csv:
        if args.GENE:
            w_csv.write("gene_name\tauc\tlog2_fc\tpvalue\ttpm_query\ttpm_subject\tauc_ci_lo\tauc_ci_hi\tfc_ci_lo\tfc_ci_hi\tpv_ci_lo\tpv_ci_hi\treverse\tgene_id\tlist_trans\n")
        else:
            w_csv.write("trans_id\tgene_name\tauc\tlog2_fc\tpvalue\ttpm_query\ttpm_subject\tauc_ci_lo\tauc_ci_hi\tfc_ci_lo\tfc_ci_hi\tpv_ci_lo\tpv_ci_hi\treverse\tgene_id\n")

        cpt_id = 0
        cpt_not_exp = 0

        # Trick to avoid generate n_sample random number for each id.
        # It's done one time and are good for each id.
        # all_rand_gauss = [gauss(0, 0.0000001) for sample in design["sample"]]
        mod_value = 10000
        if args.GENE:
            mod_value = 1000
        for select_id in list_ids:
            if cpt_id % mod_value == 0:
                print(time.strftime('%X') + ": " + str(cpt_id))

            aucs_tpm = []
            pvalues_tpm = []
            foldchanges = []
            tpm_query = []
            tpm_subject = []
            cpt_bs_not_exp = 0
            for cpt_kal_bs in xrange(0, num_kal_bs):
                sum_score = 0
                cpt_sample = 0
                scores_tpm = []
                for sample in design["sample"]:
                    score_count = kallisto[sample]["tpm"][cpt_kal_bs][cpt_id]
                    score_tpm = (score_count / kallisto[sample]["total_mass"][cpt_kal_bs]) * 1e6
                    # sum_score += score_tpm

                    # Trick to avoid identical values for mannwhitneyu
                    # all_rand_gauss[cpt_sample])
                    scores_tpm.append(score_tpm + gauss(0, 0.0000001))
                    cpt_sample += 1

                sum_score = sum(scores_tpm)
                (l2fc, mean_query, mean_subject) = get_foldchange(scores_tpm, num_query)
                abs_l2fc = abs(l2fc)
                if sum_score < min_exp or abs_l2fc < min_l2fc:
                    cpt_bs_not_exp += 1
                    if cpt_bs_not_exp > num_kal_bs/2:
                        cpt_not_exp += 1
                        break

                if num_sample_bs == 0:
                    (auc, pvalue) = auc_u_test(scores_tpm, num_query, num_subject)
                    aucs_tpm.append(auc)
                    pvalues_tpm.append(pvalue)
                    foldchanges.append(l2fc)
                    tpm_query.append(mean_query)
                    tpm_subject.append(mean_subject)
                else:
                    for cpt_bs_sample in xrange(0, num_sample_bs):
                        df_bs = create_bs_sample(scores_tpm, num_query, num_subject)
                        (auc, pvalue) = auc_u_test(df_bs, num_query, num_subject)
                        (l2fc, mean_query, mean_subject) = get_foldchange(df_bs, num_query)
                        aucs_tpm.append(auc)
                        pvalues_tpm.append(pvalue)
                        foldchanges.append(l2fc)
                        tpm_query.append(mean_query)
                        tpm_subject.append(mean_subject)

            if sum_score >= min_exp and abs_l2fc >= min_l2fc:
                aucs_tpm = sorted(aucs_tpm)
                auc_ci_lo = aucs_tpm[int(len(aucs_tpm)*0.025)]
                auc_ci_hi = aucs_tpm[int(len(aucs_tpm)*0.975)]
                median_auc = get_median(aucs_tpm)

                pvalues_tpm = sorted(pvalues_tpm)
                pv_ci_lo = pvalues_tpm[int(len(pvalues_tpm)*0.025)]
                pv_ci_hi = pvalues_tpm[int(len(pvalues_tpm)*0.975)]
                median_pvalue = get_median(pvalues_tpm)

                foldchanges = sorted(foldchanges)
                fc_ci_lo = foldchanges[int(len(foldchanges)*0.025)]
                fc_ci_hi = foldchanges[int(len(foldchanges)*0.975)]
                median_fc = get_median(foldchanges)

                tpm_query = sorted(tpm_query)
                tpm_query = get_median(tpm_query)

                tpm_subject = sorted(tpm_subject)
                tpm_subject = get_median(tpm_subject)

                reverse = 1
                if median_auc < 0.5:
                    median_auc = 1.0 - median_auc
                    tmp = auc_ci_lo
                    auc_ci_lo = 1.0 - auc_ci_hi
                    auc_ci_hi = 1.0 - tmp
                    reverse = -1

                if args.GENE:
                    gene_name = select_id
                    try:
                        gene_id = ensembl.gene_ids_of_gene_name(gene_name)[0]
                    except ValueError as e:
                        gene_name = "NA"
                        gene_id = "NA"

                else:
                    select_id = cut_version(select_id)
                    try:
                        gene_name = ensembl.gene_name_of_transcript_id(select_id)
                        gene_id = ensembl.gene_ids_of_gene_name(gene_name)[0]
                    except ValueError as e:
                        gene_name = "NA"
                        gene_id = "NA"

                line = ""
                if args.GENE:
                    line = gene_name + "\t"
                else:
                    line = select_id + "\t" + gene_name + "\t"

                line = line + str(median_auc) + "\t" + str(median_fc) + "\t" + str(median_pvalue) + "\t"
                line = line + str(tpm_query) + "\t" + str(tpm_subject) + "\t"
                line = line + str(auc_ci_lo) + "\t" + str(auc_ci_hi) + "\t"
                line = line + str(fc_ci_lo) + "\t" + str(fc_ci_hi) + "\t"
                line = line + str(pv_ci_lo) + "\t" + str(pv_ci_hi) + "\t"
                line = line + str(reverse) + "\t" + gene_id

                if args.GENE:
                    line = line + "\t" + trans_by_gene[cpt_id]

                line = line + "\n"

                w_csv.write(line)

            cpt_id += 1

        print("Number of trans/gene filtred (<" + str(min_exp) + " TPM or <" + str(min_l2fc) + "abs(L2FC)): " + str(cpt_not_exp))


def main_diff(args, argparser):

    if not os.path.exists(args.PATH_OUT):
        os.makedirs(args.PATH_OUT)

    file_out = ""
    if args.GENE:
        file_out = args.PATH_OUT + "/AUC_genes_prediction.txt"
    else:
        file_out = args.PATH_OUT + "/AUC_trans_prediction.txt"

    ensembl = EnsemblRelease(args.ENSEMBL)

    design = pd.read_csv(args.DESIGN, sep="\t")
    design[args.CLASS] = [1 if condition == args.QUERY else 0 for condition in design[args.CLASS]]
    design = design.sort_values(by=args.CLASS, ascending=False)

    kall_folder_design = False
    if ("kallisto" in list(design.columns.values)):
        kall_folder_design = True
    else:
        if not os.path.exists(args.FOLDER):
            sys.stderr.write('You need to kallisto column or -f option in our command line\n')
            sys.exit()

    id_sample1 = np.where(design[args.CLASS] == 1)[0]
    id_sample2 = np.where(design[args.CLASS] == 0)[0]
    num_query = len(id_sample1)
    num_subject = len(id_sample2)

    kallisto = dict()
    trans_by_gene = []

    file_name = ""
    if (kall_folder_design):
        file_name = os.path.join(str(design["kallisto"][0]), "abundance.h5")
    else:
        file_name = os.path.join(args.FOLDER,
                                 str(design["sample"][0]),
                                 "abundance.h5")

    f = h5py.File(file_name, 'r')
    if args.BS_KAL == -1:
        num_kal_bs = f["aux/num_bootstrap"][0]
    else:
        if args.BS_KAL <= f["aux/num_bootstrap"][0]:
            num_kal_bs = args.BS_KAL
        else:
            sys.stderr.write("Bootstrap in kallisto output file is lower than: " + str(args.BS_KAL) + "\n")

    transcripts = f["aux/ids"].value
    if args.GENE:
        print(time.strftime('%X') + ": transcripts to genes")
        genes = np.array([])

        cpt_trans = 0
        cpt_dont_found = 0
        for trans in transcripts:
            cpt_trans += 1
            trans = cut_version(trans)
            try:
                genes = np.append(genes,
                                  ensembl.gene_name_of_transcript_id(trans))
            except ValueError as e:
                if cpt_dont_found < 10:
                    sys.stderr.write("WARNING! Gene name not found for: " + trans + "\n")
                    if cpt_dont_found == 9:
                        sys.stderr.write("...\n")
                cpt_dont_found += 1
                genes = np.append(genes, "unknow_" + str(cpt_dont_found))

        if(cpt_dont_found > 0):
            sys.stderr.write("WARNING! " + str(cpt_dont_found) + "/" + str(cpt_trans) + " gene name not found!\n" )

        uniq_genes = np.unique(genes)

        print(time.strftime('%X') + ": found genes ids")
        ids_genes = []
        i_esr1 = []
        for gene in uniq_genes:
            i = np.where(genes == gene)[0]
            ids_genes.append(i)

    transcripts_len = f["aux/eff_lengths"].value

    print(time.strftime('%X') + ": loading h5")
    for index, row in design.iterrows():
        sample = row['sample']
        print(time.strftime('%X') + ": \tfor sample: " + sample)

        file_name = ""
        if (kall_folder_design):
            file_name = os.path.join(str(row['kallisto']), "abundance.h5")
        else:
            file_name = os.path.join(args.FOLDER, str(sample), "abundance.h5")

        file_h5 = os.path.join(file_name)
        f = h5py.File(file_h5, 'r')

        tpm = dict()
        total_mass = dict()
        for cpt_bs in xrange(0, num_kal_bs):
            tmp = f["bootstrap/bs" + str(cpt_bs)].value / transcripts_len
            if args.GENE:
                tmp, trans_by_gene = trans_to_gene(tmp, ids_genes, transcripts)
            tpm[cpt_bs] = tmp
            total_mass[cpt_bs] = np.sum(tmp)

        kallisto[sample] = dict()
        kallisto[sample]["tpm"] = tpm
        kallisto[sample]["total_mass"] = total_mass

    min_l2fc = args.L2FC
    min_exp = args.EXP
    num_sample_bs = args.BS_SAMPLE

    if args.GENE:
        all_auc(file_out, uniq_genes, design, num_query, num_subject, kallisto,
                num_kal_bs, num_sample_bs, min_exp, min_l2fc, ensembl,
                trans_by_gene)
    else:
        all_auc(file_out, transcripts, design, num_query, num_subject,
                kallisto, num_kal_bs, num_sample_bs, min_exp, min_l2fc,
                ensembl)
