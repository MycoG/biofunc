import pandas as pd

def load_bed(path:str, names:list=None) -> pd.DataFrame:
    """
    Load BED file as pandas DataFrame. Automatically looks for a "{path}.header" file in the same directory to use as column labels

    :param path: path to BED file
    :type path: str
    :param names: names to be used for column labels
    :type names: list
    """
    try : 
        with open(path+".header", 'r') as f:
            column_labels = f.read().strip().split("\n")
    except :
        column_labels = names

    df = pd.read_csv(path, names=column_labels, sep='\t', index_col=False)
    return df


def save_bed(df:pd.DataFrame, out:str, header=True) -> None:
    """
    save pandas dataframe as a BED file and header file

    :param df: dataframe
    :type df: pd.DataFrame
    :param out: output filename
    :type out: str
    """
    df.to_csv(out, sep="\t", index=False, header=False)
    if header:
        with open(out+".header", 'w') as f:
            f.write("\n".join(list(df.columns)) + "\n")
    return