from cyvcf2 import VCF
import numpy as np
import pandas as pd
import subprocess
import gzip
from uuid import uuid4
import sys
import argparse
import os

def main():

    how_high_image = \
    "    __  __                 __    _       __  ___  \n" \
    "   / / / /___ _      __   / /_  (_)___ _/ /_/__ \ \n" \
    "  / /_/ / __ \ | /| / /  / __ \/ / __ `/ __ \/ _/ \n" \
    " / __  / /_/ / |/ |/ /  / / / / / /_/ / / / /_/   \n" \
    "/_/ /_/\____/|__/|__/  /_/ /_/_/\__, /_/ /_(_)    \n" \
    "                               /____/             \n" \


    help = "\nHOW HIGH? Returns the allele frequency of high impact variants for a defined feature type"

    parser = argparse.ArgumentParser(description=how_high_image + help, prog='tool', formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog,max_help_position=40))

    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    #required
    required.add_argument('-v', '--vcf', type=str, help='Vcf file (indexed and annotated with SNPEff)', required=True)
    required.add_argument('-a', '--annotation', type=str, help='Annotation gff', required=True)
    required.add_argument('-w', '--window_size', type=int, help='Window size', required=False)
    required.add_argument('-p', '--populations', type=str, help='Populations file.\nCsv file with id and population in each row: [id,pop]', required=True)
    required.add_argument('-o', '--output', type=str, help='Output name (gzipped if ends with .gz)', required=True)

    # optional
    optional.add_argument('-m', '--min_high', type=int, help='Mininum number of high impact allel per site per population\n[Default = 1, i.e. singletons are discarded]', default=1, required=False)
    optional.add_argument('-t', '--type', choices = ("exon", "intron", "gene", "CDS", ), default = "exon", help='Feature type. Choose from: gene, CDS, exon and intron [Default = exon]', required=False)

    args = parser.parse_args()

    # check parameters
    print("[HOW HIGH?] Checking command line arguments")

    if os.path.exists(args.vcf) is not True:
        raise Exception(f"[HOW HIGH?] ERROR: The specified VCF {args.vcf} does not exist")

    if not os.path.exists(args.vcf + ".csi") | os.path.exists(args.vcf + ".tbi"):
        raise Exception('[HOW HIGH?] ERROR: The vcf is not indexed. The vcf can be indexed with "bcftools index" or "tabix -p vcf"')

    if os.access(os.path.dirname("./" + args.output), os.W_OK) is not True:
        raise Exception(f"[HOW HIGH?] ERROR: Cannot write output to {args.output}")

    if os.path.exists(args.populations) is not True:
        raise Exception(f"[HOW HIGH?] ERROR: The specified populations file {args.populations} does not exist")

    _stdout = subprocess.check_output(f"bcftools query -l {args.vcf}", shell=True, text=True)
    samples = [sample_id for sample_id in _stdout.split("\n") if sample_id]

    pop_df = pd.read_csv(args.populations, names=["id", "pop"])

    if not np.all(np.isin(pop_df.id, samples)):
        raise Exception(f"[HOW HIGH?] ERROR: The following sample in the population file could not be found in the VCF:\n{', '.join(pop_df.id[~np.isin(pop_df.id, vcf_header.samples)])}")

    pops = {}
    unique_pops = np.unique(pop_df["pop"]).astype(str)
    for pop in unique_pops:
        ids = pop_df[pop_df['pop'] == pop]['id']
        pops[pop] = dict()
        pops[pop]['idx'] = np.where(np.isin(samples, ids))[0]
        pops[pop]['n_callable'] = 0
        pops[pop]['n_high'] = 0

    tmp_bed = f"tmp_mask_{uuid4().hex}.bed.gz"
    tmp_vcf = f"tmp_filtered_{uuid4().hex}.vcf.gz"

    print(f"[HOW HIGH?] Subsetting {args.type}s ")

    run_command(f"""grep "{args.type}" {args.annotation} | cut -f 1,4,5 | bedtools sort | bedtools merge -i - | bgzip > {tmp_bed}""")
    run_command(f"""bcftools filter -R {tmp_bed} {args.vcf} -Oz -o {tmp_vcf}""")
    run_command(f"""bcftools index {tmp_vcf}""")

    data = subprocess.check_output(f"bcftools index -s {tmp_vcf}", shell=True, text=True)
    chromosomes = pd.DataFrame([x.split('\t') for x in data[:-1].split('\n')], columns=["id", "chr_size", "vcf_line"])

    print("[HOW HIGH?] Extracting high impact allele frequency in window")

    output_file = gzip.open(args.output, "wt") if args.output.endswith(".gz") else open(args.output, "wt")

    vcf = VCF(tmp_vcf)

    for iteration, chr in chromosomes.iterrows():

        windows = create_windows(args.window_size, chr.chr_size)

        for i, window in enumerate(windows):

            for variant in vcf(f"{chr.id}:{window[0]}-{window[1]}"):

                genotypes = np.array(variant.genotypes)

                for pop in pops:
                    pops[pop]['n_callable'] += np.sum(np.all(genotypes[pops[pop]['idx'], 0:2] != -1, axis=1))


                if variant.is_snp:

                    annotation = variant.INFO.get('ANN')

                    if "HIGH" in annotation:

                        # get alleles
                        alleles = np.array(variant.ALT)
                        alleles = np.insert(alleles, 0, variant.REF)

                        # get alleles impacts and effects
                        effects = [[] for y in range(len(alleles))]
                        effect_impacts = np.zeros(len(alleles), dtype=int)
                        for annotation_string in annotation.split(","):
                            ann = annotation_string.split("|")
                            if ann[2] == "HIGH":
                                effect_impacts[int(np.where(alleles == ann[0])[0])] += 1
                                effects[int(np.where(alleles == ann[0])[0])].append(ann[1])

                        # get number of high impact per pop
                        for pop in pops:

                            count = np.array([np.sum(genotypes[pops[pop]['idx'], 0:2] == i) for i in np.arange(len(alleles))])

                            #n_high = np.sum(count[np.where(effect_impacts > 0)[0]])

                            n_high = 0
                            for idx in np.where(effect_impacts > 0)[0]:
                                if count[idx] > args.min_high:
                                    n_high += count[idx]

                            pops[pop]['n_high'] += n_high

            for pop in pops:
                output_file.write(f"{chr.id}\t{window[0]}\t{window[1]}\t{pops[pop]['n_high']}\t{pops[pop]['n_callable']}\t{pop}\n")
                pops[pop]['n_callable'] = 0
                pops[pop]['n_high'] = 0

    output_file.close()

    run_command(f"rm {tmp_vcf} && rm {tmp_vcf}.csi && rm {tmp_bed}")

    print("[HOW HIGH?] Ran successfuly!")


def create_windows(window_size, chr_len):
    """Create windows"""

    window_starts = np.array([*range(0, int(chr_len), window_size)])
    window_ends = np.array([*range(0 + window_size, int(chr_len) + window_size, window_size)])
    window_ends[-1] = int(chr_len)
    windows = np.stack((window_starts, window_ends)).T
    return windows

def run_command(args):
    try:
        call = subprocess.run(args, text=True, capture_output=True, shell=True)
        call.check_returncode()
    except subprocess.CalledProcessError as cpe:
        if call.stdout:
            sys.exit('[HOW HIGH] ERROR: The command following command failed with error code %r:\n[X] => %s\n[X] (STDOUT): %r\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stdout, cpe.stderr))
        sys.exit('[HOW HIGH] ERROR: The command following command failed with error code %r:\n[X] => %s\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stderr.rstrip("\n")))


if __name__ == "__main__":

   main()
