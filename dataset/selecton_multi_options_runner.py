#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
import os
import itertools

SELECTON_WEB_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton/programs/selecton/selecton"
SELECTON_WEB_WRT_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton-wrt/programs/selecton-wrt/selecton"
SELECTON_WEB_WRTM_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton-wrtm/programs/selecton-wrtm/selecton"
SELECTON_WEB_WRTM_REV_PATH= "/sternadi/home/volume1/taliakustin/selecton-web/selecton-wrtm-reversible/programs/selecton/selecton"
SELECTON_WEB_APOBEC_CONTEXT_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton-apobec_context_all_branches/programs/selecton-apobec_context_all_branches/selecton"

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("--idx", dest="idx", help="index number")
    parser.add_option("-a", "--aln", dest="aln", help="alignment file")
    parser.add_option("-t", "--tree", dest="tree", help="tree file")
    parser.add_option("-o", "--output", dest="output", help="output directory")
    parser.add_option("-p", "--protein", dest="protein", help="protein name")
    (options, args) = parser.parse_args()


    idx = int(options.idx)
    aln = options.aln
    tree = options.tree
    output = options.output
    protein = options.protein


    alpha = [0.6, 1, 1.7]
    beta = [0.7, 1, 2]
    kappa = [2, 3, 4]
    m1 = [0.5, 1, 2]
    m2 = [0.5, 1, 2]
    omega = [1.2, 2, 3.3]
    prob = [0.01, 0.05, 0.2]


    dic = {}

    for i,option in enumerate(itertools.product(alpha, beta, kappa, m1, m2,  omega, prob)):
        dic[i+1] = option


    a, x, k, m1, m2, w, p = dic[idx]

    output_GA = output + "/apobec_context_model" \
                + "_a%s_b%s_k%s_m1-%s_m2-%s_w%s_p%s/" % (str(a), str(x), str(k),
                                                      str(m1), str(m2), str(w),
                                                      str(p)) \
                + protein

    if not os.path.isdir("/home/el"):
        os.system("mkdir %s" %output_GA.split(protein)[0])

    output_path = "%s_output.txt" % output_GA
    if (os.path.exists(output_path) and os.stat(output_path).st_size == 0) or not os.path.exists(output_path):
        print("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                                    "-c %s_color.txt -t %s_tree.txt -y %s -z %s -w %f -p %f -a %f -x %f -k %f -j 40"
                                                    % (SELECTON_WEB_APOBEC_CONTEXT_PATH, aln, tree, output_GA, output_GA, output_GA
                                                       , output_GA, output_GA, output_GA, m1,  m2, w, 1 - p, a, x, k))

        os.system("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                                    "-c %s_color.txt -t %s_tree.txt -y %s -z %s -w %f -p %f -a %f -x %f -k %f -j 40"
                                                    % (SELECTON_WEB_APOBEC_CONTEXT_PATH, aln, tree, output_GA, output_GA, output_GA
                                                       , output_GA, output_GA, output_GA, m1, m2, w, 1 - p, a, x, k))

    else:
        print("output file already exists - didn't run selecton")



if __name__ == "__main__":
    main()

