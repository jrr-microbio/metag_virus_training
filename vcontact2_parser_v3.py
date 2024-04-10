"""

The goal of this parser will be to report what types of genomes (column "genome_category") are present in each cluster (column "VC").  Note: We had previously done this where we parsed based on the name "this_study" and "input_db" to give us an output that reported three categories: 1) novel genera clustering with input_db ("this study" clustering with "input db" 2) known taxonomy ("this study" clustering with "known taxonomy" + anything else, and 3) Novel genera ("this study" clustering to "this study" only).

What is different this time is that it will report the "mix" of all possible combinations of categories based on the "genome category".

Example 1: A VC with "genome_category" that includes refseq_genome, river genome, and wastewater genome will be labeled as "Refseq, River, and Wastewater Genome Cluster".
Example 2: A VC with "genome_category" that includes erpe genome, river genome, wastewater genome will be labeled as "Erpe, River, and Wastewater Cluster".

The output should just be the same exact input table with an additional column labeled as "auto-typed VC category" that has the information on what genomes are present in each cluster. The reason why this is different is because previously I wasn't too concerned about a vMAG if it was not clustering to my Erpe database - but now I'm trying to populate the venn diagram below - and I need info on the other groups as well. This will also be incredibly useful for Kayla and GROWdb - as these datasets are massive and need some manual parsing if we want to get it into a similar structure. I have manually categorized with the help of your script the Erpe-related genomes already, so I can spot-check the script based on that.

"""

import re
import pandas as pd
import numpy as np
import click


def check_for_other(df):
    for i in df['Genome'].unique():
        j = i.replace('this_study_', '').replace('input_db', '')
        if len(re.findall('[A-z,0-9]*_', j)) > 1:
            return True
    return False


def check_all_genome_categories(df:pd.DataFrame):
    all_genome_categories = list(set(df['genome_category']))
    all_genome_categories = [str(i) for i in all_genome_categories]
    all_genome_categories.sort()
    if len(all_genome_categories) > 1:
        all_genome_categories[-1] = 'and ' + all_genome_categories[-1] + ' Cluster'
    return pd.Series(', '.join(all_genome_categories))


@click.command()
@click.option('-i', '--input', help="A vContact2 file, to parse",
              required=True)
@click.option('-o', '--output', help="output name, optional",
              required=False)
def parse(input, output:str = None):
    """Parse the command line, and read input"""
    data = pd.read_csv(input)
    data_gc = data.groupby('VC').apply(check_all_genome_categories).rename(columns={0: "auto-typed VC category"})
    data = data.merge(data_gc, left_on='VC', right_index=True)
    if output is None:
        output = input.replace('.csv', "_auto_typed_v3.csv")
    data.to_csv(output, index=False)

if __name__ == '__main__':
    parse()



"""
import os
os.system("python3 vcontact2_parser_v2.py -i vContact2_new_output_parser.csv")
exit
"""

def test_parser_2():
    """
    Example 1: A VC with "genome_category" that includes refseq_genome, river genome, and wastewater genome will be labeled as "Refseq, River, and Wastewater Genome Cluster".
    Example 2: A VC with "genome_category" that includes erpe genome, river genome, wastewater genome will be labeled as "Erpe, River, and Wastewater Cluster".
    """
    data = pd.DataFrame({'VC': ['ex1', 'ex1', 'ex1',
                                'ex2', 'ex2', 'ex2', ],
                  'genome_category': ['refseq_genome', 'river_genome', 'wastewater',
                                      'erpe_genome', 'river genome', 'wastewater_genome']})
    expected = pd.DataFrame({'VC': ['ex1', 'ex1', 'ex1',
                                    'ex2', 'ex2', 'ex2', ],
                  'genome_category': ['refseq_genome', 'river_genome', 'wastewater',
                                      'erpe_genome', 'river genome', 'wastewater_genome'],
                  'auto-typed VC category': [
                      'refseq_genome, river_genome, and wastewater Cluster',
                      'refseq_genome, river_genome, and wastewater Cluster',
                      'refseq_genome, river_genome, and wastewater Cluster',
                      'erpe_genome, river genome, and wastewater_genome Cluster',
                      'erpe_genome, river genome, and wastewater_genome Cluster',
                      'erpe_genome, river genome, and wastewater_genome Cluster',
                  ]})
    csv_name = 'temp_csv.csv'
    csv_out = csv_name.replace('.csv', "_auto_typed_v2.csv")
    data.to_csv(csv_name)
    parse(csv_name)
    result = pd.read_csv(csv_out)
    assert result.equals(expected)
