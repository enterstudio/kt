from .common import *


def get_argparser_diff(parser):
    parser.add_argument("-d",
                        dest="DESIGN",
                        help="Tabulated file who discribe your design, with the path of each kallisto folder (or -f)",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-f",
                        dest="FOLDER",
                        help="if all kallisto folder are in the same path folder, you can specify this folder her and not in the design file for each sample",
                        type=str,
                        default="")
    parser.add_argument("-o",
                        dest="PATH_OUT",
                        help="Path to the directory output")
    parser.add_argument("-e",
                        dest="EXP",
                        help="filter trans/genes with sum(TPM) < (Default:0.1)",
                        type=float,
                        default=0.1)
    parser.add_argument("-l",
                        dest="L2FC",
                        help="filter trans/genes with log2FC < (Default:0.2)",
                        type=float,
                        default=0.2)
    parser.add_argument("-bs-kal",
                        dest="BS_KAL",
                        help="To set number of kallisto bootstrap used (Default: equal to the number of bs in the kallisto output)",
                        type=int,
                        default=-1)
    parser.add_argument("-bs-sample",
                        dest="BS_SAMPLE",
                        help="To set number of bootstrap on sample (Default: 10)",
                        type=int,
                        default=10)
    parser.add_argument("--gene",
                        dest="GENE",
                        help="Calculate AUCs on genes (vs Transcripts)",
                        action='store_true')
    parser.add_argument("--ensembl",
                        dest="ENSEMBL",
                        help="Ensemble version (Default:79)",
                        type=int,
                        default=79)
    parser.add_argument("--class",
                        dest="CLASS",
                        help="Name column of class in the design file (Default: class)",
                        type=str,
                        default="class")
    parser.add_argument("--query",
                        dest="QUERY",
                        help="Query value in the class column (Default: relapse)",
                        type=str,
                        default="relapse")

    parser.set_defaults(GENE=False)
